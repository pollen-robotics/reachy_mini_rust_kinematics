
// WASM bindings (enabled with "wasm" feature)
#[cfg(target_arch = "wasm32")]
use wasm_bindgen::prelude::*;

use super::kinematics::{Kinematics};
use nalgebra::{Matrix4};


// ============================================================================
// WASM Bindings
// ============================================================================

#[cfg(target_arch = "wasm32")]
#[wasm_bindgen]
pub struct WasmKinematics {
    inner: Kinematics,
}

#[cfg(target_arch = "wasm32")]
#[wasm_bindgen]
impl WasmKinematics {
    /// Create a new kinematics solver from JSON string
    /// 
    /// # Arguments
    /// * `json_data` - JSON configuration as a string (fetch from file in JS)
    /// 
    /// # Returns
    /// WasmKinematics instance or error
    #[wasm_bindgen(constructor)]
    pub fn new(json_data: &str) -> Result<WasmKinematics, JsValue> {
        let kinematics = Kinematics::from_json_string(json_data)
            .map_err(|e| JsValue::from_str(&format!("Failed to load kinematics: {}", e)))?;

        Ok(WasmKinematics { inner: kinematics })
    }

    /// Inverse kinematics: calculate joint angles from platform pose
    /// 
    /// # Arguments
    /// * `t_world_platform` - 4x4 transformation matrix as flat array (16 floats, row-major)
    /// * `body_yaw` - Optional body yaw angle in radians (pass NaN for None)
    /// 
    /// # Returns
    /// Array of 6 joint angles in radians
    #[wasm_bindgen(js_name = inverseKinematics)]
    pub fn inverse_kinematics(&mut self, t_world_platform: &[f64], body_yaw: f64) -> Vec<f64> {
        if t_world_platform.len() != 16 {
            return vec![0.0; 6];
        }

        let t_world_platform = Matrix4::new(
            t_world_platform[0],
            t_world_platform[1],
            t_world_platform[2],
            t_world_platform[3],
            t_world_platform[4],
            t_world_platform[5],
            t_world_platform[6],
            t_world_platform[7],
            t_world_platform[8],
            t_world_platform[9],
            t_world_platform[10],
            t_world_platform[11],
            t_world_platform[12],
            t_world_platform[13],
            t_world_platform[14],
            t_world_platform[15],
        );

        let body_yaw_opt = if body_yaw.is_nan() {
            None
        } else {
            Some(body_yaw)
        };
        
        self.inner.inverse_kinematics(t_world_platform, body_yaw_opt)
    }

    /// Forward kinematics: calculate platform pose from joint angles
    /// 
    /// # Arguments
    /// * `joint_angles` - Array of 6 joint angles in radians
    /// * `body_yaw` - Optional body yaw angle in radians (pass NaN for None)
    /// 
    /// # Returns
    /// 4x4 transformation matrix as flat array (16 floats, row-major)
    #[wasm_bindgen(js_name = forwardKinematics)]
    pub fn forward_kinematics(&mut self, joint_angles: &[f64], body_yaw: f64) -> Vec<f64> {
        if joint_angles.len() != 6 {
            return vec![1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0];
        }

        let body_yaw_opt = if body_yaw.is_nan() {
            None
        } else {
            Some(body_yaw)
        };

        let t = self.inner.forward_kinematics(joint_angles.to_vec(), body_yaw_opt);

        vec![
            t[(0, 0)], t[(0, 1)], t[(0, 2)], t[(0, 3)],
            t[(1, 0)], t[(1, 1)], t[(1, 2)], t[(1, 3)],
            t[(2, 0)], t[(2, 1)], t[(2, 2)], t[(2, 3)],
            t[(3, 0)], t[(3, 1)], t[(3, 2)], t[(3, 3)],
        ]
    }

    /// Calculate passive joint angles from head joints and head pose
    /// 
    /// # Arguments
    /// * `head_joints` - Array of 7 floats: [yaw_body, stewart_1, ..., stewart_6]
    /// * `head_pose` - 4x4 transformation matrix as flat array (16 floats, row-major)
    /// 
    /// # Returns
    /// Array of 21 floats: passive joint angles [p1_x, p1_y, p1_z, ..., p7_x, p7_y, p7_z]
    #[wasm_bindgen(js_name = calculatePassiveJoints)]
    pub fn calculate_passive_joints(&mut self, head_joints: &[f64], head_pose: &[f64]) -> Vec<f64> {
        if head_joints.len() != 7 || head_pose.len() != 16 {
            return vec![0.0; 21];
        }

        // Convert flat head_pose array to [[f64; 4]; 4]
        let head_pose_matrix = [
            [head_pose[0], head_pose[1], head_pose[2], head_pose[3]],
            [head_pose[4], head_pose[5], head_pose[6], head_pose[7]],
            [head_pose[8], head_pose[9], head_pose[10], head_pose[11]],
            [head_pose[12], head_pose[13], head_pose[14], head_pose[15]],
        ];

        self.inner.calculate_passive_joints(head_joints, &head_pose_matrix)
    }

    /// Reset forward kinematics to a specific platform pose
    /// 
    /// # Arguments
    /// * `t_world_platform` - 4x4 transformation matrix as flat array (16 floats, row-major)
    #[wasm_bindgen(js_name = resetForwardKinematics)]
    pub fn reset_forward_kinematics(&mut self, t_world_platform: &[f64]) {
        if t_world_platform.len() != 16 {
            return;
        }

        let t = Matrix4::new(
            t_world_platform[0],
            t_world_platform[1],
            t_world_platform[2],
            t_world_platform[3],
            t_world_platform[4],
            t_world_platform[5],
            t_world_platform[6],
            t_world_platform[7],
            t_world_platform[8],
            t_world_platform[9],
            t_world_platform[10],
            t_world_platform[11],
            t_world_platform[12],
            t_world_platform[13],
            t_world_platform[14],
            t_world_platform[15],
        );

        self.inner.reset_forward_kinematics(t);
    }

    /// Safe inverse kinematics with limits
    /// 
    /// # Arguments
    /// * `t_world_platform` - 4x4 transformation matrix as flat array (16 floats, row-major)
    /// * `body_yaw` - Optional body yaw angle in radians (pass NaN for None)
    /// * `max_relative_yaw` - Optional max relative yaw limit (pass NaN for None)
    /// * `max_body_yaw` - Optional max body yaw limit (pass NaN for None)
    /// 
    /// # Returns
    /// Array of 7 values: [body_yaw_target, stewart_1, ..., stewart_6]
    #[wasm_bindgen(js_name = inverseKinematicsSafe)]
    pub fn inverse_kinematics_safe(
        &mut self,
        t_world_platform: &[f64],
        body_yaw: f64,
        max_relative_yaw: f64,
        max_body_yaw: f64,
        max_lean_angle: f64
    ) -> Vec<f64> {
        if t_world_platform.len() != 16 {
            return vec![0.0; 7];
        }

        let t = Matrix4::new(
            t_world_platform[0],
            t_world_platform[1],
            t_world_platform[2],
            t_world_platform[3],
            t_world_platform[4],
            t_world_platform[5],
            t_world_platform[6],
            t_world_platform[7],
            t_world_platform[8],
            t_world_platform[9],
            t_world_platform[10],
            t_world_platform[11],
            t_world_platform[12],
            t_world_platform[13],
            t_world_platform[14],
            t_world_platform[15],
        );

        let body_yaw_opt = if body_yaw.is_nan() { None } else { Some(body_yaw) };
        let max_rel_yaw_opt = if max_relative_yaw.is_nan() { None } else { Some(max_relative_yaw) };
        let max_body_yaw_opt = if max_body_yaw.is_nan() { None } else { Some(max_body_yaw) };
        let max_lean_angle_opt = if max_lean_angle.is_nan() { None } else { Some(max_lean_angle) };

        self.inner.inverse_kinematics_safe(t, body_yaw_opt, max_rel_yaw_opt, max_body_yaw_opt, max_lean_angle_opt)
    }
}

/// Initialize WASM module
#[cfg(target_arch = "wasm32")]
#[wasm_bindgen(start)]
pub fn wasm_init() {
    // Set panic hook for better error messages in the browser console
    #[cfg(feature = "console_error_panic_hook")]
    console_error_panic_hook::set_once();
}
