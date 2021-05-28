[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Build Status](https://travis-ci.com/muhos/ParaFROST.svg?token=YXUywHfBSpqMqyUKnyT4&branch=master)](https://travis-ci.com/muhos/ParaFROST)
# ParaFROST
ParaFROST stands for Parallel Formal Reasoning On SaTisfiability. 
It is a parallel SAT solver with GPU-accelerated inprocessing capable of harnessing NVIDIA CUDA-enabled GPUs in applying modern inprocessing tecnhiques in parallel. 
The CDCL search is built from scratch with various optimisations based on CaDiCaL heuristics (see our paper in [TACAS'21](https://gears.win.tue.nl/papers/parafrost_gpu.pdf)).
The inprocessing engine extends our previous work in SIGmA simplifier with new data structures, parallel garbage collection and more.

# Install

To install either the CPU or the GPU solvers, use the `install.sh` script which has the following usage:

`
 usage: install.sh [ <option> ... ]

 where '<option>' is one of the following

       -h or --help          print this usage summary
       -c or --cpu           install CPU solver
       -g or --gpu           install GPU solver (if CUDA exists)
       -w or --wall          compile with '-Wall' flag
       -d or --debug         compile with debugging inf|ormation
       -t or --assert        enable only code assertions
       -p or --pedantic      compile with '-pedantic' flag
       -l or --logging       enable logging (needed for verbosity level > 2)
       -s or --statistics    enable costly statistics (may impact runtime)
       -a or --all           enable all above flags except 'assert'
       --clean=<target>      remove old installation of <cpu | gpu | all> solvers
       --standard=<n>        compile with <11 | 14 | 17 > c++ standard
       --extra="flags"       pass extra "flags" to the compiler(s)
`

## GPU solver
To build the GPU solver, make sure you have a CUDA-capable GPU with pre-installed NVIDIA driver and CUDA toolkit.

For installing the driver + CUDA, run the following commands:<br>

`wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/cuda-ubuntu1804.pin`<br>
`sudo mv cuda-ubuntu1804.pin /etc/apt/preferences.d/cuda-repository-pin-600`<br>
`sudo apt-key adv --fetch-keys https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/7fa2af80.pub`<br>
`sudo add-apt-repository "deb https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/ /"`<br>
`sudo apt-get update`<br>
`sudo apt-get -y install cuda`<br>

Now the GPU solver is ready to install by running the install script via the command `./install -g'. 
The `parafrost` binary and the library `libparafrost.a` will be created by default in the solver local directory.<br>

## CPU solver
To build a CPU-only version of the solver, run `./install -c`.<br>

## Debug and Testing
Add `-t` argument with the install command to enable assertions or `-d` to collect debugging information for both the CPU and GPU solvers.<br>

The solver has a complete artifact for performance evaluation and comparisons with CaDiCaL solver.<br>
More information can be found in: https://gears.win.tue.nl/software/parafrost <br>

# Usage
The solver can be used via the command `parafrost [<option> ...][<infile>.<cnf>][<option> ...]`.<br>
For more options, type `parafrost -h` or `parafrost --helpmore`.