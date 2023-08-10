[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Build Status](https://github.com/muhos/ParaFROST/actions/workflows/test-build.yml/badge.svg)](https://github.com/muhos/ParaFROST/actions/workflows/test-build.yml)
# ParaFROST
ParaFROST stands for Parallel Formal Reasoning about SaTisfiability. 
It is a parallel SAT solver with GPU-accelerated inprocessing capable of harnessing NVIDIA CUDA-enabled GPUs in applying modern inprocessing tecnhiques in parallel. 
The CDCL search is built from scratch with various optimisations based on CaDiCaL heuristics (see our paper in [FMSD'23](https://link.springer.com/article/10.1007/s10703-023-00432-z)).
The inprocessing engine extends our previous work in SIGmA simplifier with new data structures, parallel garbage collection and more.

# Install

To install either the CPU or the GPU solvers, use the `install.sh` script which has the following usage:


&nbsp; usage: `install.sh [ <option> ... ]`<br>
&nbsp; where `<option>` is one of the following

       -h or --help          print this usage summary
       -n or --less          print less verbose messages
       -q or --quiet         be quiet (make steps still be saved in the log)
       -c or --cpu           install CPU solver
       -g or --gpu           install GPU solver (if CUDA exists)
       -w or --wall          compile with '-Wall' flag
       -d or --debug         compile with debugging information
       -t or --assert        enable only code assertions
       -p or --pedantic      compile with '-pedantic' flag
       -l or --logging       enable logging (needed for verbosity level > 2)
       -s or --statistics    enable costly statistics (may impact runtime)
       -a or --all           enable all above flags except 'assert'
	   --ncolors             disable colors in all solver outputs
       --clean=<target>      remove old installation of <cpu | gpu | all> solvers
       --standard=<n>        compile with <11 | 14 | 17> c++ standard
       --gextra="flags"      pass extra "flags" to the GPU compiler (nvcc)
	   --cextra="flags"      pass extra "flags" to the CPU compiler (g++)


## GPU solver
To build the GPU solver, make sure you have a CUDA-capable GPU with pre-installed NVIDIA driver and CUDA toolkit.

For installing CUDA v12, run the following commands on a Ubuntu 22.04:<br>

`wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/cuda-keyring_1.0-1_all.deb`<br>
`sudo dpkg -i cuda-keyring_1.0-1_all.deb`<br>
`sudo apt-get update`<br>
`sudo apt-get -y install cuda`<br>

The source code is also platform-compatible with Windows and WSL2. To install CUDA on those platforms, follow the
installation guide in https://docs.nvidia.com/cuda/.

Now the GPU solver is ready to install by running the install script via the command `./install.sh -g`. 
The `parafrost` binary and the library `libparafrost.a` will be created by default in the build directory.<br>

## CPU solver
To build a CPU-only version of the solver, run `./install.sh -c`.<br>

## Debug and Testing
Add `-t` argument with the install command to enable assertions or `-d` to collect debugging information for both the CPU and GPU solvers.<br>

The solver has a complete artifact for performance evaluation and comparisons with CaDiCaL solver.<br>
More information can be found in: https://gears.win.tue.nl/software/parafrost <br>

# Usage
The solver can be used via the command `parafrost [<option> ...][<infile>.<cnf>][<option> ...]`.<br>
For more options, type `parafrost -h` or `parafrost --helpmore`.

# Incremental Solving
ParaFROST supports incremental solving to `add`/`remove` variables or clauses incrementally while solving with assumptions if needed. A fully configurable interface to integrate ParaFROST with CBMC model checker is created here (https://github.com/muhos/gpu4bmc). A similar interface can be created to work with ParaFROST in any SAT-based bounded model checker.
