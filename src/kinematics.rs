use nalgebra::{DVector, Matrix3, Matrix3x6, Matrix4, MatrixXx6, Vector3};
use serde::Deserialize;
use std::fs;

use super::euler_utils::{
    align_vectors, euler_from_rotation_xyz, euler_from_rotation_zyz, rotation_from_euler_xyz,
    rotation_from_euler_zyz,
};

pub const HEAD_Z_OFFSET: f64 = 0.177;
pub const STEWARD_ROD_LENGTH: f64 = 0.09;
pub const MOTOR_ARM_LENGTH: f64 = 0.04;

struct Branch {
    branch_platform: Vector3<f64>,
    t_world_motor: Matrix4<f64>,
    solution: f64,
    jacobian: Matrix3x6<f64>,
    limits: Option<(f64, f64)>,
}

#[allow(non_snake_case)]
#[derive(Deserialize)]
struct Motor {
    branch_position: Vec<f64>,
    T_motor_world: Vec<Vec<f64>>,
    solution: f64,
    limits: Vec<f64>,
}

#[allow(non_snake_case)]
#[derive(Deserialize)]
struct PassiveKinematics {
    #[serde(deserialize_with = "deserialize_matrix4")]
    t_xl330_in_platform_frame: Matrix4<f64>,
    orientation_offset_in_servo_arm_frame: Vec<[f64; 3]>,
    stewart_rod_direction_in_passive_frame: Vec<[f64; 3]>,
}

fn deserialize_matrix4<'de, D>(deserializer: D) -> Result<Matrix4<f64>, D::Error>
where
    D: serde::Deserializer<'de>,
{
    let array: [[f64; 4]; 4] = serde::Deserialize::deserialize(deserializer)?;
    Ok(Matrix4::new(
        array[0][0],
        array[0][1],
        array[0][2],
        array[0][3],
        array[1][0],
        array[1][1],
        array[1][2],
        array[1][3],
        array[2][0],
        array[2][1],
        array[2][2],
        array[2][3],
        array[3][0],
        array[3][1],
        array[3][2],
        array[3][3],
    ))
}

pub struct Kinematics {
    motor_arm_length: f64,
    rod_length: f64,
    head_z_offset: f64,
    t_world_platform: Matrix4<f64>,
    line_search_maximum_iterations: usize,
    branches: Vec<Branch>,
    passives: PassiveKinematics,
}

impl Kinematics {
    pub fn new(motor_arm_length: f64, rod_length: f64, head_z_offset: f64) -> Self {
        let t_world_platform = Matrix4::identity();
        let line_search_maximum_iterations = 16;
        let passives = PassiveKinematics {
            t_xl330_in_platform_frame: Matrix4::identity(),
            orientation_offset_in_servo_arm_frame: Vec::new(),
            stewart_rod_direction_in_passive_frame: Vec::new(),
        };
        let head_z_offset = head_z_offset;
        let branches = Vec::new();
        Self {
            motor_arm_length,
            rod_length,
            head_z_offset,
            t_world_platform,
            line_search_maximum_iterations,
            branches,
            passives,
        }
    }

    pub fn add_branch(
        &mut self,
        branch_platform: Vector3<f64>,
        t_world_motor: Matrix4<f64>,
        solution: f64,
        limits: Option<(f64, f64)>,
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
            limits,
        });
    }

    pub fn init_passive_kinematics(
        &mut self,
        t_xl330_in_platform_frame: Matrix4<f64>,
        orientation_offset_in_servo_arm_frame: Vec<[f64; 3]>,
        stewart_rod_direction_in_passive_frame: Vec<[f64; 3]>,
    ) {
        self.passives = PassiveKinematics {
            t_xl330_in_platform_frame,
            orientation_offset_in_servo_arm_frame,
            stewart_rod_direction_in_passive_frame,
        };
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
        max_lean_angle: Option<f64>,
    ) -> Vec<f64> {
        let mut joint_angles: Vec<f64> = vec![0.0; self.branches.len() + 1];
        let mut body_yaw_target = 0.0;

        // if body yaw is specified, rotate the platform accordingly
        if body_yaw.is_some() {
            body_yaw_target = -body_yaw.unwrap();
            // first verify if the body yaw is within the allowed limits
            // relative yaw is the yaw difference between the current platform yaw and body yaw
            // it should stay within +/- max_relative_yaw
            if let Some(max_rel_yaw) = max_relative_yaw {
                let z_pos = t_world_platform[(2, 3)] - self.head_z_offset;
                let mut max_rel_yaw_adapt = max_rel_yaw;
                // reduce progressively max_relative_yaw to 0 at - 4cm
                if (z_pos < 0.0) && (z_pos >= -0.04) {
                    max_rel_yaw_adapt = (0.04 + z_pos) / 0.04 * max_rel_yaw;
                } else if z_pos < -0.04 {
                    max_rel_yaw_adapt = 0.0;
                }

                let current_yaw = t_world_platform[(0, 1)].atan2(t_world_platform[(0, 0)]);
                let relative_yaw = body_yaw_target - current_yaw;
                body_yaw_target =
                    current_yaw + relative_yaw.clamp(-max_rel_yaw_adapt, max_rel_yaw_adapt);
            }
            // then clamp the body yaw within +/- max_body_yaw
            // this is physically limited by the mechanical design
            if let Some(max_body_yaw) = max_body_yaw {
                body_yaw_target = body_yaw_target.clamp(-max_body_yaw, max_body_yaw);
            }
            body_yaw_target = -body_yaw_target;
        }

        // Extract rotation from platform transform and convert to ZYZ Euler angles
        let mut t_world_platform_clamped = t_world_platform;
        if let Some(max_angle) = max_lean_angle {
            // Extract 3x3 rotation matrix from 4x4 transform
            let rotation = t_world_platform.fixed_view::<3, 3>(0, 0).into_owned();
            // Convert to ZYZ Euler angles
            let mut euler_angles = euler_from_rotation_zyz(&rotation);

            // Clamp the middle angle (beta) within [-max_angle, max_angle]
            euler_angles[1] = euler_angles[1].clamp(-max_angle, max_angle);

            // Convert back to rotation matrix
            let clamped_rotation =
                rotation_from_euler_zyz(euler_angles[0], euler_angles[1], euler_angles[2]);
            // Update the transform with the clamped rotation
            for i in 0..3 {
                for j in 0..3 {
                    t_world_platform_clamped[(i, j)] = clamped_rotation[(i, j)];
                }
            }
        }

        // construct the joint angles vector
        joint_angles[0] = body_yaw_target;
        joint_angles[1..].copy_from_slice(
            &self.inverse_kinematics(t_world_platform_clamped, Some(body_yaw_target)),
        );

        // clamp each joint angle within its limits if specified
        for (i, branch) in self.branches.iter().enumerate() {
            if let Some((min_limit, max_limit)) = branch.limits {
                joint_angles[i + 1] = joint_angles[i + 1].clamp(min_limit, max_limit);
            }
        }
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

        let mut found_solution: bool = false;
        let mut tries_left = 5;

        while !found_solution && tries_left > 0 {
            tries_left -= 1;
            found_solution = true;

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
                let arm_platform_hom = self.t_world_platform.try_inverse().unwrap()
                    * branch.t_world_motor
                    * arm_motor_hom;
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

            // if the head is lower than 7cm below the initial position,
            // or if the head orientation is too extreme (looking down or rotated more than 80 degrees)
            // or if the head yaw is more than 100 degrees from the body yaw
            // RETRY with a small random orientation offset
            let eulers = euler_from_rotation_xyz(
                &self.t_world_platform.fixed_view::<3, 3>(0, 0).into_owned(),
            );
            let b_yaw = body_yaw.unwrap_or(0.0);
            if self.t_world_platform[(2, 3)] < 0.07
                || eulers[0].abs() > 1.4
                || eulers[1].abs() > 1.4
                || (eulers[2] - b_yaw).abs() > 2.5
            {
                // Use retry counter for variation
                let t = eulers[0] * 12563.618033988749895; // golden ratio for better distribution
                let x = 0.2 * (t * 32.1).sin() - 0.1;
                let y = 0.2 * (t * 43.7).cos() - 0.1;
                let z = 0.2 * (t * 15.3).sin() - 0.1 + b_yaw;
                let mut T_new =
                    Matrix4::new_translation(&Vector3::new(0.0, 0.0, self.head_z_offset));
                if self.t_world_platform[(2, 3)] < 0.07 {
                    T_new[(2, 3)] -= 0.02; // a heuristic to help convergence when too low
                }
                let mut slice = T_new.view_mut((0, 0), (3, 3));
                slice.copy_from(&rotation_from_euler_xyz(x, y, z).fixed_view::<3, 3>(0, 0));
                self.reset_forward_kinematics(T_new);
                eprintln!(
                    "Retrying forward kinematics with orientation offset ({:.3}, {:.3}, {:.3}), tries left: {}",
                    x, y, z, tries_left
                );
                found_solution = false;
            }
        }

        // if (tries_left == 0) && !found_solution {
        //     panic!("ERRORL forward kinematics did not converge after maximum retries.");
        // }

        // prepare the retun value by applying body yaw if specified
        let mut t_world_platform = self.t_world_platform;

        // rotate the body around Z if body_yaw is specified
        if let Some(yaw) = body_yaw {
            // remove the z offset
            t_world_platform[(2, 3)] -= self.head_z_offset;
            // rotate
            let rotation = nalgebra::Rotation3::from_axis_angle(
                &nalgebra::Unit::new_normalize(Vector3::z()),
                yaw,
            );
            let t_yaw = rotation.to_homogeneous();
            t_world_platform = t_yaw * t_world_platform;
            // re-apply the z offset
            t_world_platform[(2, 3)] += self.head_z_offset;
        }

        t_world_platform
    }

    /// Create a Kinematics instance from a JSON configuration file
    ///
    /// # Arguments
    /// * `path` - Path to the JSON configuration file
    ///
    /// # Returns
    /// Result containing the initialized Kinematics instance or an error message
    pub fn from_json_file(path: &str) -> Result<Self, String> {
        let data = fs::read_to_string(path).map_err(|e| format!("Unable to read file: {}", e))?;
        Self::from_json_string(&data)
    }

    /// Create a Kinematics instance from a JSON string
    ///
    /// # Arguments
    /// * `json_data` - JSON string containing kinematics configuration
    ///
    /// # Returns
    /// Result containing the initialized Kinematics instance or an error message
    pub fn from_json_string(json_data: &str) -> Result<Self, String> {
        let data_deserialized: serde_json::Value =
            serde_json::from_str(json_data).map_err(|e| format!("Unable to parse JSON: {}", e))?;

        let rod_length = data_deserialized["rod_length"]
            .as_f64()
            .ok_or("Unable to parse rod_length")?;
        let motor_arm_length = data_deserialized["motor_arm_length"]
            .as_f64()
            .ok_or("Unable to parse motor_arm_length")?;

        let mut head_z_offset = HEAD_Z_OFFSET;
        if let Some(z_offset) = data_deserialized["head_z_offset"].as_f64() {
            head_z_offset = z_offset;
        }

        let mut kinematics = Kinematics::new(motor_arm_length, rod_length, head_z_offset);

        let motors: Vec<Motor> = serde_json::from_str(&data_deserialized["motors"].to_string())
            .map_err(|e| format!("Unable to parse motors: {}", e))?;

        for motor in motors {
            let branch_position = Vector3::new(
                motor.branch_position[0],
                motor.branch_position[1],
                motor.branch_position[2],
            );
            let t_motor_world = Matrix4::new(
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

            let mut limits = None;
            if motor.limits.len() == 2 {
                limits = Some((motor.limits[0], motor.limits[1]));
            }

            kinematics.add_branch(
                branch_position,
                t_motor_world.try_inverse().unwrap(),
                solution,
                limits,
            );
        }

        let passives: PassiveKinematics =
            serde_json::from_str(&data_deserialized["passive_joint_kinematics"].to_string())
                .map_err(|e| format!("Unable to parse passives: {}", e))?;

        kinematics.init_passive_kinematics(
            passives.t_xl330_in_platform_frame,
            passives.orientation_offset_in_servo_arm_frame,
            passives.stewart_rod_direction_in_passive_frame,
        );

        let t_world_platform =
            Matrix4::new_translation(&Vector3::new(0.0, 0.0, kinematics.head_z_offset));
        kinematics.reset_forward_kinematics(t_world_platform);

        Ok(kinematics)
    }

    /// Calculate passive joint angles from head joints and head pose
    ///
    /// # Arguments
    /// * `head_joints` - Array of 7 floats: [yaw_body, stewart_1, ..., stewart_6]
    /// * `head_pose` - 4x4 transformation matrix as 16 floats (row-major)
    ///
    /// # Returns
    /// Array of 21 floats: passive joint angles [p1_x, p1_y, p1_z, ..., p7_x, p7_y, p7_z]
    #[allow(non_snake_case)]
    pub fn calculate_passive_joints(
        &mut self,
        head_joints: &[f64],
        head_pose: &[[f64; 4]; 4],
    ) -> Vec<f64> {
        if self.passives.orientation_offset_in_servo_arm_frame.len() < 6
            || self.passives.stewart_rod_direction_in_passive_frame.len() < 6
        {
            return vec![0.0; 21];
        }

        if head_joints.len() < 7 || head_pose.len() < 4 || head_pose[0].len() < 4 {
            return vec![0.0; 21];
        }

        let body_yaw = head_joints[0];

        // Build head pose matrix from row-major input
        let mut pose = Matrix4::new(
            head_pose[0][0],
            head_pose[0][1],
            head_pose[0][2],
            head_pose[0][3],
            head_pose[1][0],
            head_pose[1][1],
            head_pose[1][2],
            head_pose[1][3],
            head_pose[2][0],
            head_pose[2][1],
            head_pose[2][2],
            head_pose[2][3],
            head_pose[3][0],
            head_pose[3][1],
            head_pose[3][2],
            head_pose[3][3],
        );

        // Add head Z offset
        pose[(2, 3)] += self.head_z_offset;

        // Inverse rotation: rotate pose around Z by -body_yaw
        let cos_yaw = body_yaw.cos();
        let sin_yaw = body_yaw.sin();
        let r_z_inv = Matrix4::new(
            cos_yaw, sin_yaw, 0.0, 0.0, -sin_yaw, cos_yaw, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
            0.0, 1.0,
        );
        pose = r_z_inv * pose;

        // Pre-compute passive correction rotations
        let passive_corrections: Vec<Matrix3<f64>> = self
            .passives
            .orientation_offset_in_servo_arm_frame
            .iter()
            .map(|offset| rotation_from_euler_xyz(offset[0], offset[1], offset[2]))
            .collect();

        let mut passive_joints = vec![0.0; 21];
        let mut last_r_servo_branch = Matrix3::identity();
        let mut last_r_world_servo = Matrix3::identity();

        // T_motor_servo_arm: translation by motor_arm_length along X
        let t_motor_servo_arm = Vector3::new(self.motor_arm_length, 0.0, 0.0);

        // For each of the 6 stewart motors
        for (i, motor) in self.branches.iter().enumerate() {
            let stewart_joint = head_joints[i + 1];

            // Extract pose rotation and translation
            let pose_rot = pose.fixed_view::<3, 3>(0, 0).into_owned();
            let pose_trans = Vector3::new(pose[(0, 3)], pose[(1, 3)], pose[(2, 3)]);

            // Calculate branch position on platform in world frame
            let branch_pos = motor.branch_platform;
            let branch_pos_world = pose_rot * branch_pos + pose_trans;

            // Compute servo rotation (rotating around Z axis)
            let cos_z = stewart_joint.cos();
            let sin_z = stewart_joint.sin();
            let r_servo = Matrix3::new(cos_z, -sin_z, 0.0, sin_z, cos_z, 0.0, 0.0, 0.0, 1.0);

            // T_world_motor from motor data
            let t_world_motor = motor.t_world_motor;
            let t_world_motor_rot = t_world_motor.fixed_view::<3, 3>(0, 0).into_owned();
            let t_world_motor_trans = Vector3::new(
                t_world_motor[(0, 3)],
                t_world_motor[(1, 3)],
                t_world_motor[(2, 3)],
            );

            // Compute world servo arm position
            let servo_pos_local = r_servo * t_motor_servo_arm;
            let p_world_servo_arm = t_world_motor_rot * servo_pos_local + t_world_motor_trans;

            // Apply passive correction to orientation
            let r_world_servo = t_world_motor_rot * r_servo * passive_corrections[i];

            // Vector from servo arm to branch in world frame
            let vec_servo_to_branch = branch_pos_world - p_world_servo_arm;

            // Transform to servo frame (use transpose for inverse of rotation)
            let vec_servo_to_branch_in_servo = r_world_servo.transpose() * vec_servo_to_branch;

            // Rod direction in passive frame
            let rod_dir = Vector3::new(
                self.passives.stewart_rod_direction_in_passive_frame[i][0],
                self.passives.stewart_rod_direction_in_passive_frame[i][1],
                self.passives.stewart_rod_direction_in_passive_frame[i][2],
            );

            // Normalize and get straight line direction
            let norm_vec = vec_servo_to_branch_in_servo.norm();
            let straight_line_dir = vec_servo_to_branch_in_servo / norm_vec;

            // Align rod direction to actual direction
            let r_servo_branch = align_vectors(&rod_dir, &straight_line_dir);
            let euler = euler_from_rotation_xyz(&r_servo_branch);

            passive_joints[i * 3] = euler[0];
            passive_joints[i * 3 + 1] = euler[1];
            passive_joints[i * 3 + 2] = euler[2];

            // Save for 7th passive joint calculation
            if i == 5 {
                last_r_servo_branch = r_servo_branch;
                last_r_world_servo = r_world_servo;
            }
        }

        // 7th passive joint (XL330 on the head)
        // Head XL330 target orientation
        let t_head_xl330_rot = Matrix3::new(
            self.passives.t_xl330_in_platform_frame[(0, 0)],
            self.passives.t_xl330_in_platform_frame[(0, 1)],
            self.passives.t_xl330_in_platform_frame[(0, 2)],
            self.passives.t_xl330_in_platform_frame[(1, 0)],
            self.passives.t_xl330_in_platform_frame[(1, 1)],
            self.passives.t_xl330_in_platform_frame[(1, 2)],
            self.passives.t_xl330_in_platform_frame[(2, 0)],
            self.passives.t_xl330_in_platform_frame[(2, 1)],
            self.passives.t_xl330_in_platform_frame[(2, 2)],
        );
        let pose_rot = pose.fixed_view::<3, 3>(0, 0).into_owned();
        let r_head_xl330 = pose_rot * t_head_xl330_rot;

        // Current rod orientation with correction for 7th passive joint
        let r_rod_current = last_r_world_servo * last_r_servo_branch * passive_corrections[6];

        // Compute relative rotation
        let r_dof = r_rod_current.transpose() * r_head_xl330;
        let euler_7 = euler_from_rotation_xyz(&r_dof);

        passive_joints[18] = euler_7[0];
        passive_joints[19] = euler_7[1];
        passive_joints[20] = euler_7[2];

        passive_joints
    }
}

pub fn create_solver() -> Kinematics {
    Kinematics::from_json_file("kinematics_data.json")
        .expect("Failed to create solver from kinematics_data.json")
}

#[cfg(test)]
mod tests {

    use super::*;

    fn initialize_kinematics() -> Kinematics {
        Kinematics::from_json_file("kinematics_data.json")
            .expect("Failed to initialize kinematics for tests")
    }

    #[test]
    fn test_inverse_kinematics() {
        let mut kinematics = initialize_kinematics();
        let t_world_platform = nalgebra::Matrix4::new_translation(&nalgebra::Vector3::new(
            0.0,
            0.0,
            kinematics.head_z_offset,
        ));
        let r = kinematics.inverse_kinematics(t_world_platform, None);
        println!("IK result: {:?}", r);
        let expected_res = [
            0.6265471361608503,
            -0.6265425506844711,
            0.6265402507582252,
            -0.6265361796605958,
            0.6265371961549008,
            -0.62653923545576,
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

        let mut t = kinematics.forward_kinematics(joints.clone(), None);
        let t_flat = t.as_slice().to_vec();
        let expected_res = [
            [
                0.9434358487250847,
                0.30211528027029916,
                0.13658388180007697,
                0.0,
            ],
            [
                -0.3274858015533049,
                0.9134576760121603,
                0.2415535632431293,
                0.0,
            ],
            [
                -0.051786572790330436,
                -0.2726195729614102,
                0.9607271825637964,
                0.0,
            ],
            [
                -0.01413499352355952,
                0.01341068963871718,
                0.11999593761525008,
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
        let t_world_platform = nalgebra::Matrix4::new_translation(&nalgebra::Vector3::new(
            0.0,
            0.0,
            kinematics.head_z_offset,
        ));
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
        let t_world_platform = nalgebra::Matrix4::new_translation(&nalgebra::Vector3::new(
            0.0,
            0.0,
            kinematics.head_z_offset,
        ));
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

    #[test]
    fn test_identity_pose_zero_joints() {
        let mut kinematics = initialize_kinematics();
        // Test: Identity pose, zero joints
        let head_joints = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let head_pose = [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        let expected = [
            0.0022508907,
            0.0362949623,
            -0.1238610683,
            -0.0222426253,
            0.0013675279,
            -0.1273488284,
            -0.0036008297,
            -0.0641988484,
            -0.1120216899,
            0.0018793787,
            -0.0298951753,
            0.1255567074,
            -0.0021551464,
            -0.0346164750,
            -0.1243428060,
            0.0018360718,
            0.0291668900,
            -0.1257263345,
            0.0018226962,
            0.0291985444,
            -0.1257131448,
        ];

        let result = kinematics.calculate_passive_joints(&head_joints, &head_pose);
        println!("{:?}", result);
        assert_eq!(result.len(), 21);

        let tolerance = 0.01; // Allow 1% error
        for i in 0..21 {
            let diff = (result[i] - expected[i]).abs();
            assert!(
                diff < tolerance,
                "Mismatch at index {}: got {}, expected {}, diff {}",
                i,
                result[i],
                expected[i],
                diff
            );
        }
    }

    #[test]
    fn test_small_body_yaw() {
        let mut kinematics = initialize_kinematics();
        // Test: Small body yaw
        let head_joints = [0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let head_pose = [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        let expected = [
            0.0023094851,
            0.0309104488,
            -0.1491418088,
            -0.0265536010,
            -0.0035773668,
            -0.1030629683,
            -0.0044785419,
            -0.0648270895,
            -0.1379017245,
            0.0017013496,
            -0.0337621624,
            0.1006896894,
            -0.0021646104,
            -0.0288928516,
            -0.1495473876,
            0.0016750546,
            0.0331768126,
            -0.1008825400,
            0.0920552079,
            0.0746590292,
            -0.0940957704,
        ];

        let result = kinematics.calculate_passive_joints(&head_joints, &head_pose);
        assert_eq!(result.len(), 21);

        let tolerance = 0.01;
        for i in 0..21 {
            let diff = (result[i] - expected[i]).abs();
            assert!(
                diff < tolerance,
                "Mismatch at index {}: got {}, expected {}, diff {}",
                i,
                result[i],
                expected[i],
                diff
            );
        }
    }

    #[test]
    fn test_all_stewart_joints() {
        let mut kinematics = initialize_kinematics();
        // Test: All stewart joints at 0.5
        let head_joints = [0.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5];
        let head_pose = [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        let expected = [
            0.0201470224,
            0.0664757285,
            -0.5883623150,
            -0.0050969762,
            -0.0349257327,
            0.2740303711,
            -0.0565607056,
            -0.1953238381,
            -0.5621706414,
            -0.0002505518,
            -0.0018002749,
            -0.2765717423,
            -0.0178861002,
            -0.0589442498,
            -0.5890751964,
            -0.0004703285,
            0.0033795988,
            0.2765574117,
            0.0420138661,
            0.0441513789,
            -0.2210345269,
        ];

        let result = kinematics.calculate_passive_joints(&head_joints, &head_pose);
        assert_eq!(result.len(), 21);

        let tolerance = 0.01;
        for i in 0..21 {
            let diff = (result[i] - expected[i]).abs();
            assert!(
                diff < tolerance,
                "Mismatch at index {}: got {}, expected {}, diff {}",
                i,
                result[i],
                expected[i],
                diff
            );
        }
    }
}
