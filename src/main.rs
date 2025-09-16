use reachy_mini_rust_kinematics::Kinematics;
use serde::Deserialize;
use std::fs;

#[allow(non_snake_case)]
#[derive(Deserialize)]
struct Motor {
    branch_position: Vec<f64>,
    T_motor_world: Vec<Vec<f64>>,
    solution: f64,
}
#[allow(non_snake_case)]
fn main() {
    println!("Hello, world!");
    let data = fs::read_to_string("motors.json").expect("Unable to read file");
    let motors: Vec<Motor> = serde_json::from_str(&data).expect("JSON was not well-formatted");
    let mut kinematics = Kinematics::new(0.038, 0.09);

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

    let head_z_offset = 0.177;

    // Test inverse kinematics
    let t_world_platform =
        nalgebra::Matrix4::new_translation(&nalgebra::Vector3::new(0.0, 0.0, head_z_offset));

    kinematics.reset_forward_kinematics(t_world_platform);
    let r = kinematics.inverse_kinematics(t_world_platform);
    println!("Inverse kinematics {:?}", r);

    // Test forward kinematics
    kinematics.reset_forward_kinematics(t_world_platform);
    let joints = vec![0.3, 0.0, 0.0, 0.0, 0.0, 0.0];
    let mut t = kinematics.forward_kinematics(joints);

    // remove head_z_offset
    t[(2, 3)] -= head_z_offset;
    println!("Forward kinematics:");
    for i in 0..4 {
        for j in 0..4 {
            print!("{:>10.6} ", t[(i, j)]);
        }
        println!();
    }
}
