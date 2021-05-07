[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Build Status](https://travis-ci.com/muhos/ParaFROST.svg?token=YXUywHfBSpqMqyUKnyT4&branch=master)](https://travis-ci.com/muhos/ParaFROST)
# ParaFROST
ParaFROST stands for Parallel Formal Reasoning On SaTisfiability. It is a parallel SAT solver with GPU-accelerated inprocessing capable of harnessing NIVIDA CUDA-enabled GPUs in applying modern inprocessing tecnhiques in parallel. The CDCL search is built from scratch with various optimisations based on CaDiCaL heuristics (see our paper in [TACAS'21](https://gears.win.tue.nl/papers/parafrost_gpu.pdf)). The inprocessing engine extends our previous work in SIGmA simplifier with new data structures, parallel garbage collection and more.

# Install
## GPU solver
To install the GPU solver, make sure you have a CUDA-capable GPU with pre-installed NVIDIA driver and CUDA toolkit.

For installing the driver + CUDA, run the following commands:<br>

`wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/cuda-ubuntu1804.pin`<br>
`sudo mv cuda-ubuntu1804.pin /etc/apt/preferences.d/cuda-repository-pin-600`<br>
`sudo apt-key adv --fetch-keys https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/7fa2af80.pub`<br>
`sudo add-apt-repository "deb https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/ /"`<br>
`sudo apt-get update`<br>
`sudo apt-get -y install cuda`<br>

Now the GPU solver is ready to install by running the makefile via the command `make -C gpu`. 
The `parafrost` binary and the library `libparafrost.a` will be created by default in the solver local directory.<br>

## CPU solver
To build a CPU-only version of the solver, run `make -C cpu`.<br>

## Debug and Testing
Add `assert=1` argument with the make command to enable assertions or `debug=1` to collect debugging information for both the CPU and GPU solvers.<br>

The solver has a complete artifact for performance evaluation and comparisons with CaDiCaL solver.<br>
More information can be found in: https://gears.win.tue.nl/software/parafrost <br>

# Run
The solver can be used via the command `parafrost [<option> ...][<infile>.<cnf>][<option> ...]`.<br>
For more options, type `parafrost -h` or `parafrost --helpmore`.
