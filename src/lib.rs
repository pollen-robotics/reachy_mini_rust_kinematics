use nalgebra::{DVector, Matrix3, Matrix3x6, Matrix4, MatrixXx6, Vector3};
use pyo3::prelude::*;
use pyo3_stub_gen::{
    define_stub_info_gatherer,
    derive::{gen_stub_pyclass, gen_stub_pymethods},
};
use serde::Deserialize;

const HEAD_Z_OFFSET: f64 = 0.177;

#[gen_stub_pyclass]
#[pyclass(frozen)]
struct ReachyMiniRustKinematics {
    inner: std::sync::Mutex<Kinematics>,
}

#[gen_stub_pymethods]
#[pymethods]
impl ReachyMiniRustKinematics {
    #[new]
    fn new(motor_arm_length: f64, rod_length: f64) -> Self {
        Self {
            inner: std::sync::Mutex::new(Kinematics::new(motor_arm_length, rod_length)),
        }
    }

    fn add_branch(&self, branch_platform: [f64; 3], t_world_motor: [[f64; 4]; 4], solution: f64) {
        let branch_platform: Vector3<f64> =
            Vector3::new(branch_platform[0], branch_platform[1], branch_platform[2]);

        let t_world_motor: Matrix4<f64> = Matrix4::new(
            t_world_motor[0][0],
            t_world_motor[0][1],
            t_world_motor[0][2],
            t_world_motor[0][3],
            t_world_motor[1][0],
            t_world_motor[1][1],
            t_world_motor[1][2],
            t_world_motor[1][3],
            t_world_motor[2][0],
            t_world_motor[2][1],
            t_world_motor[2][2],
            t_world_motor[2][3],
            t_world_motor[3][0],
            t_world_motor[3][1],
            t_world_motor[3][2],
            t_world_motor[3][3],
        );
        self.inner
            .lock()
            .unwrap()
            .add_branch(branch_platform, t_world_motor, solution);
    }

    #[pyo3(signature = (t_world_platform, body_yaw=None))]
    fn inverse_kinematics(
        &self,
        t_world_platform: [[f64; 4]; 4],
        body_yaw: Option<f64>,
    ) -> Vec<f64> {
        let t_world_platform = Matrix4::new(
            t_world_platform[0][0],
            t_world_platform[0][1],
            t_world_platform[0][2],
            t_world_platform[0][3],
            t_world_platform[1][0],
            t_world_platform[1][1],
            t_world_platform[1][2],
            t_world_platform[1][3],
            t_world_platform[2][0],
            t_world_platform[2][1],
            t_world_platform[2][2],
            t_world_platform[2][3],
            t_world_platform[3][0],
            t_world_platform[3][1],
            t_world_platform[3][2],
            t_world_platform[3][3],
        );
        self.inner
            .lock()
            .unwrap()
            .inverse_kinematics(t_world_platform, body_yaw)
    }

    fn reset_forward_kinematics(&self, t_world_platform: [[f64; 4]; 4]) {
        let t_world_platform = Matrix4::new(
            t_world_platform[0][0],
            t_world_platform[0][1],
            t_world_platform[0][2],
            t_world_platform[0][3],
            t_world_platform[1][0],
            t_world_platform[1][1],
            t_world_platform[1][2],
            t_world_platform[1][3],
            t_world_platform[2][0],
            t_world_platform[2][1],
            t_world_platform[2][2],
            t_world_platform[2][3],
            t_world_platform[3][0],
            t_world_platform[3][1],
            t_world_platform[3][2],
            t_world_platform[3][3],
        );
        self.inner
            .lock()
            .unwrap()
            .reset_forward_kinematics(t_world_platform);
    }

    #[pyo3(signature = (joint_angles, body_yaw=None))]
    fn forward_kinematics(&self, joint_angles: [f64; 6], body_yaw: Option<f64>) -> [[f64; 4]; 4] {
        let t = self
            .inner
            .lock()
            .unwrap()
            .forward_kinematics(joint_angles.to_vec(), body_yaw);
        [
            [t[(0, 0)], t[(0, 1)], t[(0, 2)], t[(0, 3)]],
            [t[(1, 0)], t[(1, 1)], t[(1, 2)], t[(1, 3)]],
            [t[(2, 0)], t[(2, 1)], t[(2, 2)], t[(2, 3)]],
            [t[(3, 0)], t[(3, 1)], t[(3, 2)], t[(3, 3)]],
        ]
    }

    #[pyo3(signature = (t_world_platform, body_yaw=None, max_relative_yaw=None, max_body_yaw=None))]
    fn inverse_kinematics_safe(
        &self,
        t_world_platform: [[f64; 4]; 4],
        body_yaw: Option<f64>,
        max_relative_yaw: Option<f64>,
        max_body_yaw: Option<f64>,
    ) -> Vec<f64> {
        let t_world_platform = Matrix4::new(
            t_world_platform[0][0],
            t_world_platform[0][1],
            t_world_platform[0][2],
            t_world_platform[0][3],
            t_world_platform[1][0],
            t_world_platform[1][1],
            t_world_platform[1][2],
            t_world_platform[1][3],
            t_world_platform[2][0],
            t_world_platform[2][1],
            t_world_platform[2][2],
            t_world_platform[2][3],
            t_world_platform[3][0],
            t_world_platform[3][1],
            t_world_platform[3][2],
            t_world_platform[3][3],
        );
        self.inner.lock().unwrap().inverse_kinematics_safe(
            t_world_platform,
            body_yaw,
            max_relative_yaw,
            max_body_yaw,
        )
    }
}

struct Branch {
    branch_platform: Vector3<f64>,
    t_world_motor: Matrix4<f64>,
    solution: f64,
    jacobian: Matrix3x6<f64>,
}

pub struct Kinematics {
    motor_arm_length: f64,
    rod_length: f64,
    t_world_platform: Matrix4<f64>,
    line_search_maximum_iterations: usize,
    branches: Vec<Branch>,
    body_yaw: f64,
}

impl Kinematics {
    pub fn new(motor_arm_length: f64, rod_length: f64) -> Self {
        let t_world_platform = Matrix4::identity();
        let line_search_maximum_iterations = 16;

        let branches = Vec::new();
        Self {
            motor_arm_length,
            rod_length,
            t_world_platform,
            line_search_maximum_iterations,
            branches,
            body_yaw: 0.0,
        }
    }

    pub fn add_branch(
        &mut self,
        branch_platform: Vector3<f64>,
        t_world_motor: Matrix4<f64>,
        solution: f64,
    ) {
        // Building a 3x6 jacobian relating platform velocity to branch anchor point
        // linear velocity Linear velocity is kept as identity and angular velocity is
        // using Varignon's formula w x p, which Is anti-symmetric -p x w and used in
        // matrix form [-p]

        let mut jacobian: Matrix3x6<f64> = Matrix3x6::zeros();
        let mut slice = jacobian.view_mut((0, 0), (3, 3));
        slice += Matrix3::identity();
        let p = -branch_platform;
        let mut slice = jacobian.view_mut((0, 3), (3, 3));
        slice[(0, 1)] = -p.z;
        slice[(0, 2)] = p.y;
        slice[(1, 0)] = p.z;
        slice[(1, 2)] = -p.x;
        slice[(2, 0)] = -p.y;
        slice[(2, 1)] = p.x;

        self.branches.push(Branch {
            branch_platform,
            t_world_motor,
            solution,
            jacobian,
        });
    }

    fn wrap_angle(angle: f64) -> f64 {
        angle
            - (2.0 * std::f64::consts::PI)
                * ((angle + std::f64::consts::PI) * (1.0 / (2.0 * std::f64::consts::PI))).floor()
    }

    pub fn inverse_kinematics_safe(
        &mut self,
        t_world_platform: Matrix4<f64>,
        body_yaw: Option<f64>,
        max_relative_yaw: Option<f64>,
        max_body_yaw: Option<f64>,
    ) -> Vec<f64> {
        let mut joint_angles: Vec<f64> = vec![0.0; self.branches.len() + 1];
        let mut body_yaw_target = 0.0;
        // if body yaw is specified, rotate the platform accordingly
        if body_yaw.is_some() {
            body_yaw_target = body_yaw.unwrap();
            // first verify if the body yaw is within the allowed limits
            // relative yaw is the yaw difference between the current platform yaw and body yaw
            // it should stays within +/- max_relative_yaw
            if let Some(max_rel_yaw) = max_relative_yaw {
                let current_yaw = t_world_platform[(0, 1)].atan2(t_world_platform[(0, 0)]);
                let relative_yaw = body_yaw_target - current_yaw;
                body_yaw_target = current_yaw + relative_yaw.clamp(-max_rel_yaw, max_rel_yaw);
                body_yaw_target = -body_yaw_target;
            }
            // then clamp the body yaw within +/- max_body_yaw
            // this is physically limited by the mechanical design
            if let Some(max_body_yaw) = max_body_yaw {
                body_yaw_target = body_yaw_target.clamp(-max_body_yaw, max_body_yaw);
            }
        }

        // construct the joint angles vector
        joint_angles[0] = body_yaw_target;
        joint_angles[1..]
            .copy_from_slice(&self.inverse_kinematics(t_world_platform, Some(body_yaw_target)));
        joint_angles
    }

    #[allow(non_snake_case)]
    pub fn inverse_kinematics(
        &mut self,
        t_world_platform: Matrix4<f64>,
        body_yaw: Option<f64>,
    ) -> Vec<f64> {
        let mut joint_angles: Vec<f64> = vec![0.0; self.branches.len()];
        let rs = self.motor_arm_length;
        let rp = self.rod_length;

        let mut t_world_platform_target = t_world_platform;

        // if body yaw is specified, rotate the platform accordingly
        if body_yaw.is_some() {
            let yaw = body_yaw.unwrap();
            let rotation = nalgebra::Rotation3::from_axis_angle(
                &nalgebra::Unit::new_normalize(Vector3::z()),
                -yaw,
            );
            let t_yaw = rotation.to_homogeneous();
            t_world_platform_target = t_yaw * t_world_platform;
        }

        for (k, branch) in self.branches.iter().enumerate() {
            let t_world_motor_inv = branch.t_world_motor.try_inverse().unwrap();
            let branch_motor = t_world_motor_inv
                * t_world_platform_target
                * Matrix4::new(
                    1.0,
                    0.0,
                    0.0,
                    branch.branch_platform.x,
                    0.0,
                    1.0,
                    0.0,
                    branch.branch_platform.y,
                    0.0,
                    0.0,
                    1.0,
                    branch.branch_platform.z,
                    0.0,
                    0.0,
                    0.0,
                    1.0,
                );
            let px = branch_motor[(0, 3)];
            let py = branch_motor[(1, 3)];
            let pz = branch_motor[(2, 3)];

            let x = px.powi(2) + 2.0 * px * rs + py.powi(2) + pz.powi(2) - rp.powi(2) + rs.powi(2);
            let y = 2.0 * py * rs
                + branch.solution
                    * (-(px.powi(4))
                        - 2.0 * px.powi(2) * py.powi(2)
                        - 2.0 * px.powi(2) * pz.powi(2)
                        + 2.0 * px.powi(2) * rp.powi(2)
                        + 2.0 * px.powi(2) * rs.powi(2)
                        - py.powi(4)
                        - 2.0 * py.powi(2) * pz.powi(2)
                        + 2.0 * py.powi(2) * rp.powi(2)
                        + 2.0 * py.powi(2) * rs.powi(2)
                        - pz.powi(4)
                        + 2.0 * pz.powi(2) * rp.powi(2)
                        - 2.0 * pz.powi(2) * rs.powi(2)
                        - rp.powi(4)
                        + 2.0 * rp.powi(2) * rs.powi(2)
                        - rs.powi(4))
                    .sqrt();

            joint_angles[k] = Self::wrap_angle(2.0 * y.atan2(x));
        }
        joint_angles
    }

    pub fn reset_forward_kinematics(&mut self, t_world_platform: Matrix4<f64>) {
        self.t_world_platform = t_world_platform;
    }

    #[allow(non_snake_case)]
    pub fn forward_kinematics(
        &mut self,
        joint_angles: Vec<f64>,
        body_yaw: Option<f64>,
    ) -> Matrix4<f64> {
        if self.branches.len() != 6 {
            panic!("Forward kinematics requires exactly 6 joint angles");
        }

        let mut J = MatrixXx6::<f64>::zeros(6);
        let mut errors = DVector::<f64>::zeros(6);
        let mut arms_motor: Vec<Vector3<f64>> = Vec::new();

        for k in 0..self.branches.len() {
            let branch = &self.branches[k];

            // Computing the position of motor arm in the motor frame
            let arm_motor = self.motor_arm_length
                * Vector3::new(joint_angles[k].cos(), joint_angles[k].sin(), 0.0);
            arms_motor.push(arm_motor);

            // Expressing the tip of motor arm in the platform frame
            // Convert arm_motor to homogeneous coordinates for multiplication
            let arm_motor_hom = arm_motor.push(1.0);
            let arm_platform_hom =
                self.t_world_platform.try_inverse().unwrap() * branch.t_world_motor * arm_motor_hom;
            let arm_platform = arm_platform_hom.fixed_rows::<3>(0).into_owned();

            // Computing the current distance
            let current_distance = (arm_platform - branch.branch_platform).norm();

            // Computing the arm-to-branch vector in platform frame
            let arm_branch_platform: Vector3<f64> = branch.branch_platform - arm_platform;

            // Computing the jacobian of the distance
            let mut slice = J.view_mut((k, 0), (1, 6));
            slice += arm_branch_platform.transpose() * branch.jacobian;
            errors[k] = self.rod_length - current_distance;
        }

        // If the error is sufficiently high, performs a line-search along the direction given by the jacobian inverse
        if errors.norm() > 1e-6 {
            let mut V = J.pseudo_inverse(1e-6).unwrap() * errors.clone();
            for _i in 0..self.line_search_maximum_iterations {
                let mut T: Matrix4<f64> = Matrix4::identity();
                T[(0, 3)] = V[0];
                T[(1, 3)] = V[1];
                T[(2, 3)] = V[2];

                let norm = V.fixed_rows::<3>(3).norm();
                if norm.abs() > 1e-6 {
                    let tail = V.fixed_rows::<3>(3).normalize();
                    let axis = nalgebra::Unit::new_normalize(tail);
                    let rotation = nalgebra::Rotation3::from_axis_angle(&axis, norm);
                    let linear = rotation.matrix();
                    let mut slice = T.view_mut((0, 0), (3, 3));
                    slice.copy_from(linear);
                }
                let t_world_platform2 = self.t_world_platform * T;

                let mut new_errors = DVector::<f64>::zeros(self.branches.len());
                for k in 0..self.branches.len() {
                    let branch = &self.branches[k];

                    let arm_motor_hom = arms_motor[k].push(1.0);
                    let arm_platform_hom = t_world_platform2.try_inverse().unwrap()
                        * branch.t_world_motor
                        * arm_motor_hom;
                    let arm_platform = arm_platform_hom.fixed_rows::<3>(0).into_owned();
                    let current_distance = (arm_platform - branch.branch_platform).norm();

                    new_errors[k] = self.rod_length - current_distance;
                }

                if new_errors.norm() < errors.norm() {
                    self.t_world_platform = t_world_platform2;
                    break;
                } else {
                    for j in 0..V.len() {
                        V[j] *= 0.5;
                    }
                }
            }
        }

        // prepare the retun value by applying body yaw if specified
        let mut t_world_platform = self.t_world_platform;

        // rotate the body around Z if body_yaw is specified
        if let Some(yaw) = body_yaw {
            // remove the z offset
            t_world_platform[(2, 3)] -= HEAD_Z_OFFSET;
            // rotate
            let rotation = nalgebra::Rotation3::from_axis_angle(
                &nalgebra::Unit::new_normalize(Vector3::z()),
                yaw,
            );
            let t_yaw = rotation.to_homogeneous();
            t_world_platform = t_yaw * t_world_platform;
            // re-apply the z offset
            t_world_platform[(2, 3)] += HEAD_Z_OFFSET;
        }

        t_world_platform
    }
}

#[pyo3::pymodule]
fn reachy_mini_rust_kinematics(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<ReachyMiniRustKinematics>()?;
    Ok(())
}

define_stub_info_gatherer!(stub_info);

#[cfg(test)]
mod tests {
    use std::fs;

    use super::*;

    #[allow(non_snake_case)]
    #[derive(Deserialize)]
    struct Motor {
        branch_position: Vec<f64>,
        T_motor_world: Vec<Vec<f64>>,
        solution: f64,
    }

    fn initialize_kinematics() -> Kinematics {
        let mut kinematics = Kinematics::new(0.038, 0.09);
        let data = fs::read_to_string("motors.json").expect("Unable to read file");
        let motors: Vec<Motor> = serde_json::from_str(&data).expect("Unable to parse JSON");
        for motor in motors {
            let branch_position = nalgebra::Vector3::new(
                motor.branch_position[0],
                motor.branch_position[1],
                motor.branch_position[2],
            );
            let T_motor_world = nalgebra::Matrix4::new(
                motor.T_motor_world[0][0],
                motor.T_motor_world[0][1],
                motor.T_motor_world[0][2],
                motor.T_motor_world[0][3],
                motor.T_motor_world[1][0],
                motor.T_motor_world[1][1],
                motor.T_motor_world[1][2],
                motor.T_motor_world[1][3],
                motor.T_motor_world[2][0],
                motor.T_motor_world[2][1],
                motor.T_motor_world[2][2],
                motor.T_motor_world[2][3],
                motor.T_motor_world[3][0],
                motor.T_motor_world[3][1],
                motor.T_motor_world[3][2],
                motor.T_motor_world[3][3],
            );
            let solution = if motor.solution != 0.0 { 1.0 } else { -1.0 };
            kinematics.add_branch(
                branch_position,
                T_motor_world.try_inverse().unwrap(),
                solution,
            );
        }

        // Test inverse kinematics
        let t_world_platform =
            nalgebra::Matrix4::new_translation(&nalgebra::Vector3::new(0.0, 0.0, HEAD_Z_OFFSET));

        kinematics.reset_forward_kinematics(t_world_platform);
        kinematics
    }

    #[test]
    fn test_inverse_kinematics() {
        let mut kinematics = initialize_kinematics();
        let t_world_platform =
            nalgebra::Matrix4::new_translation(&nalgebra::Vector3::new(0.0, 0.0, HEAD_Z_OFFSET));
        let r = kinematics.inverse_kinematics(t_world_platform, None);
        let expected_res = [
            0.5469084013213722,
            -0.6911929467384811,
            0.6290593106168814,
            -0.6290625053607944,
            0.6911968541984359,
            -0.5469156644896231,
        ];
        assert!(
            r.iter()
                .zip(expected_res.iter())
                .all(|(a, b)| (a - b).abs() < 1e-6)
        );
    }

    #[test]
    fn test_forward_kinematics() {
        let mut kinematics = initialize_kinematics();
        let joints = vec![0.3, 0.0, 0.0, 0.0, 0.0, 0.0];
        let mut t = kinematics.forward_kinematics(joints, None);
        t[(2, 3)] -= HEAD_Z_OFFSET;
        let t_flat = t.as_slice().to_vec();
        let expected_res = [
            [
                0.9500544312762279,
                0.2784767970164817,
                0.14088027234444309,
                0.0,
            ],
            [
                -0.303810126556134,
                0.928527235235177,
                0.21339301869663957,
                0.0,
            ],
            [
                -0.07138616542684614,
                -0.24553583638638127,
                0.9667554853403684,
                0.0,
            ],
            [
                -0.032726798589555045,
                0.01043105987061096,
                -0.04973316892312707,
                1.0,
            ],
        ];
        let expected_flat: Vec<f64> = expected_res
            .iter()
            .flat_map(|row| row.iter())
            .copied()
            .collect();
        assert!(
            t_flat
                .iter()
                .zip(expected_flat.iter())
                .all(|(a, b)| (a - b).abs() < 1e-6)
        );
    }

    // test ik + fk consistency
    #[test]
    fn test_ik_fk_consistency() {
        let mut kinematics = initialize_kinematics();
        let t_world_platform =
            nalgebra::Matrix4::new_translation(&nalgebra::Vector3::new(0.0, 0.0, HEAD_Z_OFFSET));
        let r = kinematics.inverse_kinematics(t_world_platform, None);
        kinematics.reset_forward_kinematics(t_world_platform);
        let mut t = kinematics.forward_kinematics(r.clone(), None);
        for _ in 0..100 {
            t = kinematics.forward_kinematics(r.clone(), None);
        }
        let t_flat = t.as_slice().to_vec();
        let expected_res = t_world_platform.as_slice().to_vec();
        assert!(
            t_flat
                .iter()
                .zip(expected_res.iter())
                .all(|(a, b)| (a - b).abs() < 1e-6)
        );
    }
    // test ik + fk consistency with body yaw
    #[test]
    fn test_ik_fk_consistency_body_yaw() {
        let body_yaw = 0.1;
        let mut kinematics = initialize_kinematics();
        let t_world_platform =
            nalgebra::Matrix4::new_translation(&nalgebra::Vector3::new(0.0, 0.0, HEAD_Z_OFFSET));
        let r = kinematics.inverse_kinematics(t_world_platform, Some(body_yaw));
        kinematics.reset_forward_kinematics(t_world_platform);
        let mut t = kinematics.forward_kinematics(r.clone(), Some(body_yaw));
        for _ in 0..100 {
            t = kinematics.forward_kinematics(r.clone(), Some(body_yaw));
        }
        let t_flat = t.as_slice().to_vec();
        let expected_res = t_world_platform.as_slice().to_vec();
        assert!(
            t_flat
                .iter()
                .zip(expected_res.iter())
                .all(|(a, b)| (a - b).abs() < 1e-4)
        );
    }
}
