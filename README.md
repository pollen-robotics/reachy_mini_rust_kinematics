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