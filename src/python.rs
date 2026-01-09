use super::kinematics::{HEAD_Z_OFFSET, Kinematics, MOTOR_ARM_LENGTH, STEWARD_ROD_LENGTH};
use nalgebra::{Matrix4, Vector3};

// Python bindings (enabled with "python" feature, which is default)
#[cfg(feature = "python")]
use pyo3::prelude::*;
#[cfg(feature = "python")]
use pyo3_stub_gen::{
    define_stub_info_gatherer,
    derive::{gen_stub_pyclass, gen_stub_pymethods},
};

#[cfg(feature = "python")]
#[gen_stub_pyclass]
#[pyclass(frozen)]
struct ReachyMiniRustKinematics {
    inner: std::sync::Mutex<Kinematics>,
}
#[cfg(feature = "python")]
#[gen_stub_pymethods]
#[pymethods]
impl ReachyMiniRustKinematics {
    #[new]
    #[pyo3(signature = (motor_arm_length=None, rod_length=None, head_z_offset=None, json_file_path=None))]
    fn new(
        motor_arm_length: Option<f64>,
        rod_length: Option<f64>,
        head_z_offset: Option<f64>,
        json_file_path: Option<String>,
    ) -> Self {
        if let Some(json_file_path) = json_file_path {
            if json_file_path != "" {
                let kinematics = Kinematics::from_json_file(&json_file_path).unwrap();
                return Self {
                    inner: std::sync::Mutex::new(kinematics),
                };
            }
        }
        let motor_arm_length = motor_arm_length.unwrap_or(MOTOR_ARM_LENGTH);
        let rod_length = rod_length.unwrap_or(STEWARD_ROD_LENGTH);
        let head_z_offset = head_z_offset.unwrap_or(HEAD_Z_OFFSET);
        Self {
            inner: std::sync::Mutex::new(Kinematics::new(
                motor_arm_length,
                rod_length,
                head_z_offset,
            )),
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
            .add_branch(branch_platform, t_world_motor, solution, None);
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

    #[pyo3(signature = (t_world_platform, body_yaw=None, max_relative_yaw=None, max_body_yaw=None, max_lean_angle=Some(0.78)))]
    fn inverse_kinematics_safe(
        &self,
        t_world_platform: [[f64; 4]; 4],
        body_yaw: Option<f64>,
        max_relative_yaw: Option<f64>,
        max_body_yaw: Option<f64>,
        max_lean_angle: Option<f64>,
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
            max_lean_angle,
        )
    }

    fn init_passive_kinematics(
        &self,
        t_xl330_in_platform_frame: [[f64; 4]; 4],
        orientation_offset_in_servo_arm_frame: Vec<[f64; 3]>,
        stewart_rod_direction_in_passive_frame: Vec<[f64; 3]>,
    ) {
        self.inner.lock().unwrap().init_passive_kinematics(
            Matrix4::new(
                t_xl330_in_platform_frame[0][0],
                t_xl330_in_platform_frame[0][1],
                t_xl330_in_platform_frame[0][2],
                t_xl330_in_platform_frame[0][3],
                t_xl330_in_platform_frame[1][0],
                t_xl330_in_platform_frame[1][1],
                t_xl330_in_platform_frame[1][2],
                t_xl330_in_platform_frame[1][3],
                t_xl330_in_platform_frame[2][0],
                t_xl330_in_platform_frame[2][1],
                t_xl330_in_platform_frame[2][2],
                t_xl330_in_platform_frame[2][3],
                t_xl330_in_platform_frame[3][0],
                t_xl330_in_platform_frame[3][1],
                t_xl330_in_platform_frame[3][2],
                t_xl330_in_platform_frame[3][3],
            ),
            orientation_offset_in_servo_arm_frame,
            stewart_rod_direction_in_passive_frame,
        );
    }

    fn calculate_passive_joints(
        &self,
        head_joints: [f64; 7],
        head_pose: [[f64; 4]; 4],
    ) -> [f64; 21] {
        self.inner
            .lock()
            .unwrap()
            .calculate_passive_joints(&head_joints, &head_pose)
            .try_into()
            .unwrap()
    }
}

// ============================================================================
// Python Bindings
// ============================================================================
#[cfg(feature = "python")]
#[pyo3::pymodule]
fn reachy_mini_rust_kinematics(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<ReachyMiniRustKinematics>()?;
    Ok(())
}

#[cfg(feature = "python")]
define_stub_info_gatherer!(stub_info);
