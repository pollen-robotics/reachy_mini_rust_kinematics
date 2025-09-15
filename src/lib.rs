use std::collections::linked_list;

use nalgebra::{Matrix3, Matrix3x6, Matrix4, Vector3};

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

        // === cpp code ===
        // Eigen::MatrixXd jacobian = Eigen::MatrixXd::Zero(3, 6);
        // jacobian.block(0, 0, 3, 3) = Eigen::Matrix3d::Identity();

        // Eigen::Vector3d p = -branch_platform;
        // jacobian.block(0, 3, 3, 3) << 0, -p.z(), p.y(), p.z(), 0, -p.x(), -p.y(),
        //     p.x(), 0;
        // branches.push_back({branch_platform, T_world_motor, solution, jacobian});
        // =================
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

            // rust : y.atan2(x);
            // cpp : atan2(y, x)

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

        // cpp code below
        // Eigen::VectorXd Kinematics::inverse_kinematics(Eigen::Affine3d T_world_platform) {
        //     Eigen::VectorXd joint_angles(branches.size());

        //     double rs = motor_arm_length;
        //     double rp = rod_length;

        //     for (int k = 0; k < branches.size(); k++) {
        //         Branch &branch = branches[k];

        //         Eigen::Vector3d branch_motor = branch.T_world_motor.inverse() *
        //                                     T_world_platform * branch.branch_platform;
        //         double px = branch_motor.x();
        //         double py = branch_motor.y();
        //         double pz = branch_motor.z();

        //         joint_angles[k] =
        //             2 *
        //             atan2(
        //                 (2 * py * rs +
        //                 branch.solution *
        //                     sqrt(
        //                         -(pow(px, 4)) - 2 * pow(px, 2) * pow(py, 2) -
        //                         2 * pow(px, 2) * pow(pz, 2) + 2 * pow(px, 2) * pow(rp, 2) +
        //                         2 * pow(px, 2) * pow(rs, 2) - pow(py, 4) -
        //                         2 * pow(py, 2) * pow(pz, 2) + 2 * pow(py, 2) * pow(rp, 2) +
        //                         2 * pow(py, 2) * pow(rs, 2) - pow(pz, 4) +
        //                         2 * pow(pz, 2) * pow(rp, 2) - 2 * pow(pz, 2) * pow(rs, 2) -
        //                         pow(rp, 4) + 2 * pow(rp, 2) * pow(rs, 2) - pow(rs, 4))),
        //                 (pow(px, 2) + 2 * px * rs + pow(py, 2) + pow(pz, 2) - pow(rp, 2) +
        //                 pow(rs, 2)));

        //         joint_angles[k] = wrap_angle(joint_angles[k]);
        //     }

        //     return joint_angles;
        // }
    }
}
