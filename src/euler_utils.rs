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
pub  fn align_vectors(from: &Vector3<f64>, to: &Vector3<f64>) -> Matrix3<f64> {
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