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

## Cross-compiling for Raspberry Pi (from Ubuntu)

The `--no-default-features` flag disables the `python` feature, which avoids the pyo3 dependency (requires a host Python).

### Pi Zero (32-bit, ARMv6)

```bash
sudo apt install gcc-arm-linux-gnueabihf
rustup target add arm-unknown-linux-gnueabihf
cargo build --target arm-unknown-linux-gnueabihf --no-default-features --release
```

### Pi Zero 2 (64-bit, AArch64)

```bash
sudo apt install gcc-aarch64-linux-gnu
rustup target add aarch64-unknown-linux-gnu
cargo build --target aarch64-unknown-linux-gnu --no-default-features --release
```

### Running benchmarks on the Pi

Cross-compile the benchmark binary without running it:

```bash
# For Pi Zero:
cargo bench --target arm-unknown-linux-gnueabihf --no-default-features --no-run
# For Pi Zero 2 (64-bit):
cargo bench --target aarch64-unknown-linux-gnu --no-default-features --no-run
```

Then copy the binary and `motors.json` to the Pi and run:

```bash
scp target/<target>/release/deps/kinematics_bench-* pi@<pi-ip>:~/
scp motors.json pi@<pi-ip>:~/
ssh pi@<pi-ip> ./kinematics_bench-* --bench
```

