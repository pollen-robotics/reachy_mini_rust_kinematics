#![allow(dead_code)]
use nalgebra::{Matrix3, Vector3};

/// Create rotation matrix from euler angles (xyz intrinsic = Z * Y * X matrix order)
/// This matches scipy's R.from_euler('xyz', angles)
pub fn rotation_from_euler_xyz(x: f64, y: f64, z: f64) -> Matrix3<f64> {
    let cx = x.cos();
    let sx = x.sin();
    let cy = y.cos();
    let sy = y.sin();
    let cz = z.cos();
    let sz = z.sin();

    // Intrinsic xyz = Rz * Ry * Rx (matrix multiplication order)
    // Result[i,j] = sum over k of Rz[i,k] * (Ry * Rx)[k,j]
    Matrix3::new(
        cy * cz,
        cz * sx * sy - cx * sz,
        cx * cz * sy + sx * sz,
        cy * sz,
        cx * cz + sx * sy * sz,
        cx * sy * sz - cz * sx,
        -sy,
        cy * sx,
        cx * cy,
    )
}

/// Extract euler angles (XYZ order) from rotation matrix
pub fn euler_from_rotation_xyz(r: &Matrix3<f64>) -> [f64; 3] {
    let sy = r[(0, 2)];

    if sy.abs() < 0.99999 {
        let x = (-r[(1, 2)]).atan2(r[(2, 2)]);
        let y = sy.asin();
        let z = (-r[(0, 1)]).atan2(r[(0, 0)]);
        [x, y, z]
    } else {
        // Gimbal lock
        let x = r[(2, 1)].atan2(r[(1, 1)]);
        let y = if sy > 0.0 {
            std::f64::consts::FRAC_PI_2
        } else {
            -std::f64::consts::FRAC_PI_2
        };
        let z = 0.0;
        [x, y, z]
    }
}

/// Align vectors: find rotation that aligns 'from' to 'to'
/// Similar to scipy.spatial.transform.Rotation.align_vectors
pub fn align_vectors(from: &Vector3<f64>, to: &Vector3<f64>) -> Matrix3<f64> {
    let from_n = from.normalize();
    let to_n = to.normalize();

    let dot = from_n.dot(&to_n);

    // If vectors are nearly parallel
    if dot > 0.99999 {
        return Matrix3::identity();
    }

    // If vectors are nearly opposite
    if dot < -0.99999 {
        // Find a perpendicular axis
        let mut perp = Vector3::new(1.0, 0.0, 0.0).cross(&from_n);
        if perp.norm() < 0.001 {
            perp = Vector3::new(0.0, 1.0, 0.0).cross(&from_n);
        }
        let axis = perp.normalize();
        // Rotate 180 degrees around perpendicular axis
        let k = Matrix3::new(
            0.0, -axis.z, axis.y, axis.z, 0.0, -axis.x, -axis.y, axis.x, 0.0,
        );
        return Matrix3::identity() + 2.0 * k * k;
    }

    // General case: Rodrigues' rotation formula
    let cross = from_n.cross(&to_n);
    let s = cross.norm();
    let c = dot;

    let k = Matrix3::new(
        0.0, -cross.z, cross.y, cross.z, 0.0, -cross.x, -cross.y, cross.x, 0.0,
    );

    Matrix3::identity() + k + k * k * ((1.0 - c) / (s * s))
}

/// Create rotation matrix from ZYZ euler angles
/// ZYZ convention: first rotate by alpha around Z, then by beta around new Y, then by gamma around new Z
pub fn rotation_from_euler_zyz(alpha: f64, beta: f64, gamma: f64) -> Matrix3<f64> {
    let ca = alpha.cos();
    let sa = alpha.sin();
    let cb = beta.cos();
    let sb = beta.sin();
    let cg = gamma.cos();
    let sg = gamma.sin();

    Matrix3::new(
        ca * cg - cb * sa * sg,
        -ca * sg - cb * cg * sa,
        sb * sa,
        cg * sa + ca * cb * sg,
        ca * cb * cg - sa * sg,
        -ca * sb,
        sb * sg,
        cg * sb,
        cb,
    )
}
/// Extract ZYZ euler angles from rotation matrix
/// Returns [alpha, beta, gamma] where rotation is Rz(alpha) * Ry(beta) * Rz(gamma)
pub fn euler_from_rotation_zyz(r: &Matrix3<f64>) -> [f64; 3] {
    // r33 = cos(beta)
    let cb = r[(2, 2)].clamp(-1.0, 1.0);
    let beta = cb.acos();

    // Use a sin(beta) test for gimbal lock
    let sb = beta.sin();

    if sb.abs() > 1e-8 {
        // alpha from r13 = sin(alpha) sin(beta), r23 = -cos(alpha) sin(beta)
        let alpha = r[(0, 2)].atan2(-r[(1, 2)]);
        // gamma from r31 = sin(beta) sin(gamma), r32 = sin(beta) cos(gamma)
        let gamma = r[(2, 0)].atan2(r[(2, 1)]);
        [alpha, beta, gamma]
    } else {
        // Gimbal lock: beta ≈ 0 or beta ≈ π
        // Only alpha+gamma (or alpha-gamma) is observable depending on convention.
        // A common choice: set alpha = 0 and put everything into gamma.
        let alpha = 0.0;

        // When beta ≈ 0, R ≈ Rz(alpha+gamma)
        // Angle = atan2(r21, r11) = atan2(r[(1,0)], r[(0,0)])
        // (also works reasonably for beta ≈ π with some convention choices)
        let gamma = r[(1, 0)].atan2(r[(0, 0)]);
        [alpha, beta, gamma]
    }
}
