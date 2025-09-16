from reachy_mini_rust_kinematics import ReachyMiniRustKinematics


import time
import json
import numpy as np

kin = ReachyMiniRustKinematics(0.038, 0.09)
head_z_offset = 0.177

with open("motors.json", "r") as f:
    motors = json.load(f)

for motor in motors:
    kin.add_branch(
        np.array(motor["branch_position"]),
        np.linalg.inv(motor["T_motor_world"]),
        1 if motor["solution"] else -1,
    )

T_world_platform = np.eye(4)
T_world_platform[:3, 3][2] = 0.177
# T_world_platform = tf.translation_matrix((0, 0, 0.177))
r = kin.inverse_kinematics(T_world_platform)
print("Inverse kinematics : ", r)

kin.reset_forward_kinematics(T_world_platform)
joints = np.array([0.3, 0.0, 0.0, 0.0, 0.0, 0.0])


T = np.array(kin.forward_kinematics(joints))
print(T)
T[2, 3] -= head_z_offset
print("Forward kinematics\n", T)
