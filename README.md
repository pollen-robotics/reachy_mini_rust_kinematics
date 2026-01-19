# Reachy Mini Rust Kinematics 

Translation of https://github.com/pollen-robotics/reachy_mini_cpp_kinematics

Analytical Inverse Kinematics, Numerical Forward Kinematics

## To install locally 
```bash
pip install maturin
```

## To build the wheel
```bash
pip install -e . --verbose
```

See the examples in the `examples/` directory.

## To install the wheel

```bash
cd `target/wheels`
pip install reachy_mini_rust_kinematics...
```

## WASM Build Instructions

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

See [WASM_BUILD.md](WASM_BUILD.md)
Also see an example in [examples](./examples/wasm_example.html)


```js
import init, { WasmKinematics } from './pkg/reachy_mini_rust_kinematics.js'; //or './reachy-mini-js/reachy_mini_rust_kinematics.js'

// Initialize WASM module
await init();
console.log('WASM module initialized');

// Load kinematics configuration from JSON file
const response = await fetch('./kinematics_data.json');
if (!response.ok) {
    throw new Error(`Failed to fetch kinematics_data.json: ${response.statusText}`);
}
const jsonData = await response.text();
console.log('Kinematics data loaded');

// Create kinematics solver with JSON string
kinematics = new WasmKinematics(jsonData);
console.log('Kinematics module initialized');

// Example: Solve IK for a given end-effector pose
const t_world_platform = [1.0, 0.0, 0.0, 0.5,
                            0.0, 1.0, 0.0, 0.0,
                            0.0, 0.0, 1.0, 0.5,
                            0.0, 0.0, 0.0, 1.0];
const body_yaw = 0.0;
const ikSolutions = kinematics.inverseKinematics(t_world_platform, body_yaw);
console.log('IK Solutions:', ikSolutions);
```
