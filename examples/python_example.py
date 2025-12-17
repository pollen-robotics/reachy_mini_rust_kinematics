#!/usr/bin/env python3
"""
Simple Python example using the reachy_mini_rust_kinematics library.

This example demonstrates:
- Creating a kinematics solver from JSON configuration
- Running inverse kinematics (IK)
- Running forward kinematics (FK)
- Calculating passive joint angles
"""

import numpy as np
from reachy_mini_rust_kinematics import ReachyMiniRustKinematics


def main():
    print("=== Reachy Mini Rust Kinematics - Python Example ===\n")

    # Create kinematics solver (manually configure)
    # In real use, you'd load these from kinematics_data.json
    motor_arm_length = 0.04
    rod_length = 0.08
    
    solver = ReachyMiniRustKinematics(json_file_path="kinematics_data.json")  # Placeholder, configure manually below
    
    # Add branches (example values - in real use, load from JSON)
    # These would typically come from kinematics_data.json
    print("Note: This example shows API usage.")
    print("For full functionality, use the pre-configured solver from kinematics_data.json\n")

    # Example 1: Identity pose at head height
    print("Example 1: Inverse Kinematics")
    print("-" * 50)
    
    # Create identity transformation at head height (0.177m)
    t_world_platform = np.array([
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.177],
        [0.0, 0.0, 0.0, 1.0]
    ])
    
    print(f"Target platform pose:\n{t_world_platform}\n")
    
    # Note: This will work once branches are added
    joint_angles = solver.inverse_kinematics(t_world_platform)
    print(f"Joint angles (rad): {joint_angles}\n")

    # Example 2: Forward Kinematics
    print("\nExample 2: Forward Kinematics")
    print("-" * 50)
    
    # Example joint angles (all zeros)
    joint_angles = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    print(f"Input joint angles: {joint_angles}")
    
    # Reset FK state to initial pose
    solver.reset_forward_kinematics(t_world_platform)
    
    # Calculate platform pose
    pose = solver.forward_kinematics(joint_angles)
    print(f"Resulting platform pose:\n{pose}\n")

    # Example 3: Safe IK with limits
    print("\nExample 3: Safe Inverse Kinematics with Limits")
    print("-" * 50)
    
    body_yaw = 0.1  # 0.1 radians
    max_relative_yaw = 0.5
    max_body_yaw = 0.8
    
    print(f"Body yaw: {body_yaw} rad")
    print(f"Max relative yaw: {max_relative_yaw} rad")
    print(f"Max body yaw: {max_body_yaw} rad")
    
    result = solver.inverse_kinematics_safe(
        t_world_platform,
        body_yaw,
        max_relative_yaw,
        max_body_yaw
    )
    print(f"Safe IK result (7 values): {result}\n")

    # Example 4: Passive Joints Calculation
    print("\nExample 4: Calculate Passive Joints")
    print("-" * 50)
    
    # Head joints: [yaw_body, stewart_1, ..., stewart_6]
    head_joints = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    
    # Head pose (identity matrix)
    head_pose = np.eye(4)
    
    print(f"Head joints: {head_joints}")
    print(f"Head pose:\n{head_pose}")
    
    # Note: Requires passive kinematics initialization
    passive_joints = solver.calculate_passive_joints(head_joints, head_pose)
    print(f"\nPassive joints (21 values - 7 joints × 3 DOF):")
    for i in range(7):
        joint_num = i + 1
        x, y, z = passive_joints[i*3:(i+1)*3]
        print(f"  Passive joint {joint_num}: x={x:.6f}, y={y:.6f}, z={z:.6f}")

    print("\n" + "=" * 50)
    print("Example completed!")
    print("\nTo use with full configuration:")
    print("1. Ensure kinematics_data.json is in the working directory")
    print("2. Load and parse the JSON to configure all branches")
    print("3. Initialize passive kinematics from the JSON data")


if __name__ == "__main__":
    main()
