# Building and Using WASM Module

This document explains how to build and use the WebAssembly (WASM) version of the Reachy Mini kinematics library.

## Prerequisites

Install `wasm-pack`:
```bash
cargo install wasm-pack
```

## Building for WASM

Build the WASM module with the wasm feature enabled:

```bash
wasm-pack build --target web --no-default-features --features wasm
```

Or for Node.js:
```bash
wasm-pack build --target nodejs --no-default-features --features wasm
```

This will generate files in the `pkg/` directory.

## Usage in JavaScript/TypeScript

### Loading the Module (Web)

```javascript
import init, { WasmKinematics } from './pkg/reachy_mini_rust_kinematics.js';

// Load kinematics_data.json
const response = await fetch('kinematics_data.json');
const jsonData = await response.text();

// Initialize WASM module
await init();

// Create kinematics solver
const kinematics = new WasmKinematics(jsonData);
```

### Inverse Kinematics

```javascript
// Create identity transform at head height
const t_world_platform = new Float64Array([
    1, 0, 0, 0,
    0, 1, 0, 0,
    0, 0, 1, 0.177,
    0, 0, 0, 1
]);

// Calculate joint angles (pass NaN for no body yaw)
const jointAngles = kinematics.inverseKinematics(t_world_platform, NaN);
console.log('Joint angles:', jointAngles);
```

### Forward Kinematics

```javascript
// Joint angles from IK or elsewhere
const jointAngles = new Float64Array([0.3, 0.0, 0.0, 0.0, 0.0, 0.0]);

// Calculate platform pose (pass NaN for no body yaw)
const pose = kinematics.forwardKinematics(jointAngles, NaN);
console.log('Platform pose:', pose);
```

### Safe Inverse Kinematics with Limits

```javascript
const t_world_platform = new Float64Array([
    1, 0, 0, 0,
    0, 1, 0, 0,
    0, 0, 1, 0.177,
    0, 0, 0, 1
]);

// With body yaw and limits
const bodyYaw = 0.1;
const maxRelativeYaw = 0.5;
const maxBodyYaw = 0.8;

const result = kinematics.inverseKinematicsSafe(
    t_world_platform, 
    bodyYaw, 
    maxRelativeYaw, 
    maxBodyYaw
);
// Returns: [body_yaw_target, stewart_1, ..., stewart_6]
console.log('Safe IK result:', result);
```

### Calculate Passive Joints

```javascript
// Head joints: [yaw_body, stewart_1, ..., stewart_6]
const headJoints = new Float64Array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);

// Head pose (4x4 matrix, row-major)
const headPose = new Float64Array([
    1, 0, 0, 0,
    0, 1, 0, 0,
    0, 0, 1, 0,
    0, 0, 0, 1
]);

// Calculate 21 passive joint angles (7 joints × 3 DOF)
const passiveJoints = kinematics.calculatePassiveJoints(headJoints, headPose);
console.log('Passive joints:', passiveJoints);
```

### Reset Forward Kinematics

```javascript
const initialPose = new Float64Array([
    1, 0, 0, 0,
    0, 1, 0, 0,
    0, 0, 1, 0.177,
    0, 0, 0, 1
]);

kinematics.resetForwardKinematics(initialPose);
```

## TypeScript Types

For TypeScript users, the generated `.d.ts` file in the `pkg/` directory provides type definitions.

Example:
```typescript
import init, { WasmKinematics } from './pkg/reachy_mini_rust_kinematics.js';

const kinematics: WasmKinematics = new WasmKinematics(jsonData);
const jointAngles: Float64Array = kinematics.inverseKinematics(pose, NaN);
```

## Notes

- All angles are in **radians**
- Transformation matrices are **4×4** in **row-major order** (16 floats)
- Use `NaN` for optional parameters when you don't want to specify them
- The `kinematics_data.json` file must be embedded at initialization time
