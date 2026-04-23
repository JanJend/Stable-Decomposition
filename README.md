
# Multiparameter Persistence Pruning

A C++ implementation for computing the pruning of a multiparameter persistence module.

## Overview

This tool reads a presentation matrix from a file, computes the pruning, and optionally saves the result as a minimal presentation matrix.

## Building
The following builds two binaries `pruning` and `tests` inside the directory `build/`.
```bash
mkdir build && cd build
cmake ..
make
```
To build with debug symbols, replace `cmake ..` with `CMAKE_BUILD_TYPE=Debug cmake ..`.

## Usage

```bash
./build/pruning <file_path> [options]
```

### Options

- `--delta <value>` - Set delta threshold value (default: 1% of the range of all degrees)
- `--no-output` - Skip saving the output file
- `--no-timers` - Disable timing output (doesn't work yet)

### Examples

```bash
# Basic usage with default delta
./matrix_pruning input.scc

# Custom delta value
./matrix_pruning input.scc --delta 0.05

# Process without saving output
./matrix_pruning input.scc --no-output --no-timers

# Combine options
./matrix_pruning input.scc --delta 0.01 --no-timers
```

## Input Format

The program expects input files in `.scc` - sparse chain complex - format (link to paper).

## Output

Output files are automatically named with the pattern: `<input_name>_pru<delta><extension>`

Example: `torus3_largestcomp.scc` → `torus3_largestcomp_pru0.0200.scc`

## Default Behavior

If no arguments are provided, the program uses `tests/torus3_largestcomp.scc` as the default input file.

## Version

## Version

Current version: **0.2**

## Authors

- [Havard Bjerkevik]
- [Jan Jendrysiak](https://github.com/jendjan)
- [Fabian Lenzen](https://gitlab.com/flenzen)

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

## Citation

If you use this software in your research, please cite:

```bibtex
@software{multiparameter_pruning,
  author = {Havard Bjerkevik and Jan Jendrysiak and Fabian Lenzen},
  title = {Multiparameter Persistence Pruning},
  year = {2025},
  url = {https://github.com/jendjan/Stable-Decomposition},
  license = {GPL-3.0}
}
