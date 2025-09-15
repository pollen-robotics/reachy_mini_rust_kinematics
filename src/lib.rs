use core::slice;

use nalgebra::{DVector, Matrix3, Matrix3x6, Matrix4, MatrixXx6, Vector3};

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

    pub fn inverse_kinematics(&mut self, t_world_platform: Matrix4<f64>) -> Vec<f64> {
        // TODO octuple check this against cpp if something is acting weird
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

    // reference cpp implementation
    // Eigen::Affine3d Kinematics::forward_kinematics(Eigen::VectorXd joint_angles) {
    //     if (branches.size() != 6) {
    //         throw std::runtime_error("Forward kinematics requires exactly 6 branches");
    //     }

    //     Eigen::MatrixXd J = Eigen::MatrixXd::Zero(6, 6);
    //     Eigen::VectorXd errors = Eigen::VectorXd::Zero(6);
    //     std::vector<Eigen::Vector3d> arms_motor;

    //     for (int k = 0; k < branches.size(); k++) {
    //         Branch &branch = branches[k];

    //         // Computing the position of motor arm in the motor frame
    //         Eigen::Vector3d arm_motor =
    //             motor_arm_length *
    //             Eigen::Vector3d(cos(joint_angles[k]), sin(joint_angles[k]), 0);
    //         arms_motor.push_back(arm_motor);

    //         // Expressing the tip of motor arm in the platform frame
    //         Eigen::Vector3d arm_platform =
    //             T_world_platform.inverse() * branch.T_world_motor * arm_motor;

    //         // Computing the current distance
    //         double current_distance = (arm_platform - branch.branch_platform).norm();

    //         // Computing the arm-to-branch vector in platform frame
    //         Eigen::Vector3d armBranch_platform = branch.branch_platform - arm_platform;

    //         // Computing the jacobian of the distance
    //         J.block(k, 0, 1, 6) = armBranch_platform.transpose() * branch.jacobian;
    //         errors(k) = rod_length - current_distance;
    //     }

    //     // If the error is sufficiently high, performs a line-search along the
    //     // direction given by the jacobian inverse
    //     if (errors.norm() > 1e-6) {
    //         Eigen::VectorXd V = J.inverse() * errors;

    //         for (int i = 0; i < line_search_maximum_iterations; i++) {
    //         Eigen::Affine3d T = Eigen::Affine3d::Identity();
    //         T.translation() = V.head(3);

    //         double norm = V.tail(3).norm();
    //         if (fabs(norm) > 1e-6) {
    //             T.linear() =
    //                 Eigen::AngleAxisd(norm, V.tail(3).normalized()).toRotationMatrix();
    //         }
    //         Eigen::Affine3d T_world_platform2 = T_world_platform * T;

    //         Eigen::VectorXd new_errors(6);
    //         for (int k = 0; k < branches.size(); k++) {
    //             Branch &branch = branches[k];

    //             Eigen::Vector3d arm_platform =
    //                 T_world_platform2.inverse() * branch.T_world_motor * arms_motor[k];
    //             double current_distance =
    //                 (arm_platform - branch.branch_platform).norm();

    //             new_errors(k) = rod_length - current_distance;
    //         }

    //         if (new_errors.norm() < errors.norm()) {
    //             T_world_platform = T_world_platform2;
    //             break;
    //         } else {
    //             V = V * 0.5;
    //         }
    //         }
    //     }

    //     return T_world_platform;
    //     }

    pub fn forward_kinematics(&mut self, joint_angles: Vec<f64>) -> Matrix4<f64> {
        // replicate cpp implementation here

        if joint_angles.len() != self.branches.len() {
            panic!("Forward kinematics requires exactly 6 joint angles");
        }

        let mut J = MatrixXx6::<f64>::zeros(self.branches.len());
        let mut errors = DVector::<f64>::zeros(self.branches.len());
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
            let arm_branch_platform = branch.branch_platform - arm_platform;

            // Computing the jacobian of the distance
            // J.block(k, 0, 1, 6) = armBranch_platform.transpose() * branch.jacobian;
            // errors(k) = rod_length - current_distance;
            let mut slice = J.view_mut((k, 0), (1, 6));
            slice += arm_branch_platform.transpose() * branch.jacobian;
            errors[k] = self.rod_length - current_distance;
        }

        // If the error is sufficiently high, performs a line-search along the direction given by the jacobian inverse
        if errors.norm() > 1e-6 {
            let mut V = J.pseudo_inverse(1e-6).unwrap() * errors.clone();
            for i in 0..self.line_search_maximum_iterations {
                let mut T = Matrix4::identity();
                T[(0, 3)] = V[0];
                T[(1, 3)] = V[1];
                T[(2, 3)] = V[2];

                let norm = (V[3].powi(2) + V[4].powi(2) + V[5].powi(2)).sqrt();
                if norm.abs() > 1e-6 {
                    let axis = Vector3::new(V[3], V[4], V[5]) / norm;
                    let cos_theta = norm.cos();
                    let sin_theta = norm.sin();
                    let one_minus_cos = 1.0 - cos_theta;

                    let rotation = Matrix3::new(
                        cos_theta + axis.x * axis.x * one_minus_cos,
                        axis.x * axis.y * one_minus_cos - axis.z * sin_theta,
                        axis.x * axis.z * one_minus_cos + axis.y * sin_theta,
                        axis.y * axis.x * one_minus_cos + axis.z * sin_theta,
                        cos_theta + axis.y * axis.y * one_minus_cos,
                        axis.y * axis.z * one_minus_cos - axis.x * sin_theta,
                        axis.z * axis.x * one_minus_cos - axis.y * sin_theta,
                        axis.z * axis.y * one_minus_cos + axis.x * sin_theta,
                        cos_theta + axis.z * axis.z * one_minus_cos,
                    );
                    T.fixed_slice_mut::<3, 3>(0, 0).copy_from(&rotation);
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
                    // divide V by 2
                    for j in 0..V.len() {
                        V[j] *= 0.5;
                    }
                }
            }
        }

        self.t_world_platform
    }
}
