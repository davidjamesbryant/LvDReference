# LvDReference
This is a reference implementation for the LvD algorithm.

The algorithm is describe in
https://arxiv.org/abs/2601.19064 and in subsequent versions

This implementation is used to illustrate the algorithm and to make the simulation benchmarks in the paper.



## Repository Structure

This repository contains reference implementations of the LvD likelihood algorithm and associated tools, developed as part of ongoing research.


## Building

Requires CMake >= 3.20 and a C++17 compiler.

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

This builds all executables (`RunSimulation`, `TreeHeightSimulation`, `TreeHeights`, `LvDLikelihoods`, `BranchLengthMCMC`) into `build/`.

To build a single target:

```
cmake --build build --target TreeHeightSimulation -j
```

To use the reference (unoptimized, easier-to-validate) tensor calculations, add `-DUSE_REFERENCE_TENSOR=ON` to the first `cmake` command.
