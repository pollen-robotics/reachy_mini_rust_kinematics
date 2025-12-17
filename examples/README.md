# Examples

This directory contains example code demonstrating how to use the Reachy Mini Rust Kinematics library in different environments.

## Python Example

**File:** `python_example.py`

Demonstrates using the library from Python with PyO3 bindings.

### Running:

```bash
# Make sure the library is installed
cd ..
maturin develop --release

# Run the example
python examples/python_example.py
```

### What it shows:
- Creating a kinematics solver
- Running inverse kinematics (IK)
- Running forward kinematics (FK)
- Safe IK with joint limits
- Calculating passive joint angles

## JavaScript/WASM Example

**File:** `wasm_example.html`

Interactive web page demonstrating the library compiled to WebAssembly.

### Setup:

1. **Build the WASM module:**
   ```bash
   cd ..
   wasm-pack build --target web --no-default-features --features wasm
   ```

2. **Copy the generated pkg directory:**
   ```bash
   cp -r pkg examples/
   ```

3. **Serve the example:**
   ```bash
   cd examples
   python -m http.server 8000
   ```

4. **Open in browser:**
   Navigate to `http://localhost:8000/wasm_example.html`

### What it shows:
- Loading and initializing WASM module
- Interactive IK/FK calculations
- Safe IK with constraints
- Passive joints calculation
- Real-time updates in the browser

## Requirements

### Python Example
- Python 3.7+
- numpy
- maturin (for building)
- The compiled reachy_mini_rust_kinematics library

### WASM Example
- wasm-pack
- Modern web browser with WASM support
- HTTP server (e.g., Python's http.server)
- kinematics_data.json configuration file

## Notes

- Both examples expect `kinematics_data.json` to be available
- The Python example shows API usage patterns even when data isn't loaded
- The WASM example requires the pkg/ directory from wasm-pack build
- For production use, properly configure all branches from the JSON data
