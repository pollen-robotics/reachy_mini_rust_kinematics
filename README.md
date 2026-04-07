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

## Cross-compiling for Raspberry Pi Zero 2 (from Ubuntu)

The `--no-default-features` flag disables the `python` feature, which avoids the pyo3 dependency (requires a host Python).

```bash
sudo apt install gcc-aarch64-linux-gnu
rustup target add aarch64-unknown-linux-gnu
cargo build --target aarch64-unknown-linux-gnu --no-default-features --release
```

### Running benchmarks on the Pi

Cross-compile the benchmark binary without running it:

```bash
cargo bench --target aarch64-unknown-linux-gnu --no-default-features --no-run
```

Then copy the binary and `motors.json` to the Pi and run:

```bash
scp target/aarch64-unknown-linux-gnu/release/deps/kinematics_bench-* pi@<pi-ip>:~/
scp motors.json pi@<pi-ip>:~/
ssh pi@<pi-ip> ./kinematics_bench-* --bench
```

