use nalgebra::{DVector, Matrix3, Matrix3x6, Matrix4, MatrixXx6, Vector3};
use pyo3::prelude::*;
use pyo3_stub_gen::{
    define_stub_info_gatherer,
    derive::{gen_stub_pyclass, gen_stub_pymethods},
};

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

    fn inverse_kinematics(&self, t_world_platform: [[f64; 4]; 4]) -> Vec<f64> {
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
            .inverse_kinematics(t_world_platform)
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

    fn forward_kinematics(&self, joint_angles: [f64; 6]) -> [[f64; 4]; 4] {
        let t = self
            .inner
            .lock()
            .unwrap()
            .forward_kinematics(joint_angles.to_vec());
        [
            [t[(0, 0)], t[(0, 1)], t[(0, 2)], t[(0, 3)]],
            [t[(1, 0)], t[(1, 1)], t[(1, 2)], t[(1, 3)]],
            [t[(2, 0)], t[(2, 1)], t[(2, 2)], t[(2, 3)]],
            [t[(3, 0)], t[(3, 1)], t[(3, 2)], t[(3, 3)]],
        ]
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

    #[allow(non_snake_case)]
    pub fn inverse_kinematics(&mut self, t_world_platform: Matrix4<f64>) -> Vec<f64> {
        let mut joint_angles: Vec<f64> = vec![0.0; self.branches.len()];
        let rs = self.motor_arm_length;
        let rp = self.rod_length;

        for (k, branch) in self.branches.iter().enumerate() {
            let t_world_motor_inv = branch.t_world_motor.try_inverse().unwrap();
            let branch_motor = t_world_motor_inv
                * t_world_platform
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
    pub fn forward_kinematics(&mut self, joint_angles: Vec<f64>) -> Matrix4<f64> {
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

        self.t_world_platform
    }
}

#[pyo3::pymodule]
fn reachy_mini_rust_kinematics(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<ReachyMiniRustKinematics>()?;
    Ok(())
}

define_stub_info_gatherer!(stub_info);
