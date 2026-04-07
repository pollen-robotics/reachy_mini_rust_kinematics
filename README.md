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

## To install the wheel

```bash
cd `target/wheels`
pip install reachy_mini_rust_kinematics...
```

## Cross-compiling for Raspberry Pi Zero (from Ubuntu)

The Pi Zero uses an ARMv6 processor. To cross-compile the pure Rust library (without Python bindings):

### 1. Install the toolchain

```bash
sudo apt install gcc-arm-linux-gnueabihf
rustup target add arm-unknown-linux-gnueabihf
```

### 2. Build

```bash
cargo build --target arm-unknown-linux-gnueabihf --no-default-features --release
```

The `--no-default-features` flag disables the `python` feature, which avoids the pyo3 dependency (requires a host Python).

### 3. Run benchmarks on the Pi

Cross-compile the benchmark binary without running it:

```bash
cargo bench --target arm-unknown-linux-gnueabihf --no-default-features --no-run
```

Then copy the binary and `motors.json` to the Pi and run:

```bash
scp target/arm-unknown-linux-gnueabihf/release/deps/kinematics_bench-* pi@<pi-ip>:~/
scp motors.json pi@<pi-ip>:~/
ssh pi@<pi-ip> ./kinematics_bench-* --bench
```

