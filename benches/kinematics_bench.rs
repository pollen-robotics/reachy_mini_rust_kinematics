use criterion::{black_box, criterion_group, criterion_main, Criterion};
use nalgebra::{Matrix4, UnitQuaternion, Vector3};
use reachy_mini_rust_kinematics::Kinematics;
use serde::Deserialize;
use std::fs;

const HEAD_Z_OFFSET: f64 = 0.177;

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
        let branch_position =
            Vector3::new(motor.branch_position[0], motor.branch_position[1], motor.branch_position[2]);
        #[rustfmt::skip]
        let t_motor_world = Matrix4::new(
            motor.T_motor_world[0][0], motor.T_motor_world[0][1], motor.T_motor_world[0][2], motor.T_motor_world[0][3],
            motor.T_motor_world[1][0], motor.T_motor_world[1][1], motor.T_motor_world[1][2], motor.T_motor_world[1][3],
            motor.T_motor_world[2][0], motor.T_motor_world[2][1], motor.T_motor_world[2][2], motor.T_motor_world[2][3],
            motor.T_motor_world[3][0], motor.T_motor_world[3][1], motor.T_motor_world[3][2], motor.T_motor_world[3][3],
        );
        let solution = if motor.solution != 0.0 { 1.0 } else { -1.0 };
        kinematics.add_branch(branch_position, t_motor_world.try_inverse().unwrap(), solution);
    }

    let t_world_platform = Matrix4::new_translation(&Vector3::new(0.0, 0.0, HEAD_Z_OFFSET));
    kinematics.reset_forward_kinematics(t_world_platform);
    kinematics
}

fn bench_inverse_kinematics(c: &mut Criterion) {
    let mut kinematics = initialize_kinematics();
    let t_world_platform = Matrix4::new_translation(&Vector3::new(0.0, 0.0, HEAD_Z_OFFSET));

    c.bench_function("ik_home_position", |b| {
        b.iter(|| kinematics.inverse_kinematics(black_box(t_world_platform), None))
    });
}

fn bench_inverse_kinematics_with_body_yaw(c: &mut Criterion) {
    let mut kinematics = initialize_kinematics();
    let t_world_platform = Matrix4::new_translation(&Vector3::new(0.0, 0.0, HEAD_Z_OFFSET));

    c.bench_function("ik_with_body_yaw", |b| {
        b.iter(|| kinematics.inverse_kinematics(black_box(t_world_platform), Some(0.3)))
    });
}

fn bench_inverse_kinematics_tilted(c: &mut Criterion) {
    let mut kinematics = initialize_kinematics();

    let rotation = UnitQuaternion::from_euler_angles(0.1, 0.15, 0.0);
    let mut t_world_platform = rotation.to_homogeneous();
    t_world_platform[(2, 3)] = HEAD_Z_OFFSET;

    c.bench_function("ik_tilted_pose", |b| {
        b.iter(|| kinematics.inverse_kinematics(black_box(t_world_platform), None))
    });
}

fn bench_forward_kinematics(c: &mut Criterion) {
    let mut kinematics = initialize_kinematics();

    // Get valid joint angles from IK at home position
    let t_world_platform = Matrix4::new_translation(&Vector3::new(0.0, 0.0, HEAD_Z_OFFSET));
    let joints = kinematics.inverse_kinematics(t_world_platform, None);

    c.bench_function("fk_from_home_joints", |b| {
        b.iter(|| kinematics.forward_kinematics(black_box(joints.clone()), None))
    });
}

fn bench_forward_kinematics_with_body_yaw(c: &mut Criterion) {
    let mut kinematics = initialize_kinematics();

    let t_world_platform = Matrix4::new_translation(&Vector3::new(0.0, 0.0, HEAD_Z_OFFSET));
    let joints = kinematics.inverse_kinematics(t_world_platform, Some(0.3));

    c.bench_function("fk_with_body_yaw", |b| {
        b.iter(|| kinematics.forward_kinematics(black_box(joints.clone()), Some(0.3)))
    });
}

fn bench_forward_kinematics_perturbed(c: &mut Criterion) {
    let mut kinematics = initialize_kinematics();

    let joints = vec![0.3, 0.0, 0.0, 0.0, 0.0, 0.0];

    c.bench_function("fk_perturbed_joints", |b| {
        b.iter(|| kinematics.forward_kinematics(black_box(joints.clone()), None))
    });
}

fn bench_ik_fk_roundtrip(c: &mut Criterion) {
    let mut kinematics = initialize_kinematics();
    let t_world_platform = Matrix4::new_translation(&Vector3::new(0.0, 0.0, HEAD_Z_OFFSET));

    c.bench_function("ik_fk_roundtrip", |b| {
        b.iter(|| {
            let joints = kinematics.inverse_kinematics(black_box(t_world_platform), None);
            kinematics.forward_kinematics(black_box(joints), None)
        })
    });
}

fn bench_inverse_kinematics_safe(c: &mut Criterion) {
    let mut kinematics = initialize_kinematics();
    let t_world_platform = Matrix4::new_translation(&Vector3::new(0.0, 0.0, HEAD_Z_OFFSET));

    c.bench_function("ik_safe", |b| {
        b.iter(|| {
            kinematics.inverse_kinematics_safe(
                black_box(t_world_platform),
                Some(0.3),
                Some(0.5),
                Some(1.0),
            )
        })
    });
}

criterion_group!(
    benches,
    bench_inverse_kinematics,
    bench_inverse_kinematics_with_body_yaw,
    bench_inverse_kinematics_tilted,
    bench_inverse_kinematics_safe,
    bench_forward_kinematics,
    bench_forward_kinematics_with_body_yaw,
    bench_forward_kinematics_perturbed,
    bench_ik_fk_roundtrip,
);
criterion_main!(benches);
