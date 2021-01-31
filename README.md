[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Build Status](https://travis-ci.com/muhos/ParaFROST.svg?token=YXUywHfBSpqMqyUKnyT4&branch=master)](https://travis-ci.com/muhos/ParaFROST)
# ParaFROST
ParaFROST stands for Parallel Formal Reasoning of Satisfiability. It is a parallel SAT solver with GPU-accelerated inprocessing capable of harnessing NIVIDA CUDA-enabled GPUs in applying modern inprocessing tecnhiques in parallel. The CDCL search is built from scratch with various optimisations based on CaDiCaL heuristics (see our paper: ). The inprocessing engine extends our previous work in SIGmA simplifier with new data structures, parallel garbage collection and more.

# Install
Run `make -C gpu` to build the gpu solver or `make -C cpu` to build a CPU-only version. Add `assert=1` argument with the make command to enable assertions or `debug=1` to collect debugging information for both the CPU and GPU solvers. 
The `parafrost` binary and the library `libparafrost.a` will be created by default in build subdirectory. 
Make sure you have a CUDA-capable GPU before running the GPU solver.

The solver has a complete artefact for performance evaluation and comparisons with CaDiCaL solver. More information can be found in: https://gears.win.tue.nl/software/parafrost

# Run
The solver can be used via the command `parafrost [<option> ...][<infile>.<cnf>][<option> ...]`. For more options, type `parafrost -h` or `parafrost --helpmore`.
