[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Build Status](https://app.travis-ci.com/muhos/ParaFROST.svg?branch=master)](https://app.travis-ci.com/muhos/ParaFROST)
# ParaFROST
ParaFROST stands for Parallel Formal Reasoning On SaTisfiability. 
It is a parallel SAT solver with GPU-accelerated inprocessing capable of harnessing NVIDIA CUDA-enabled GPUs in applying modern inprocessing tecnhiques in parallel. 
The CDCL search is built from scratch with various optimisations based on CaDiCaL heuristics (see our paper in [TACAS'21](https://gears.win.tue.nl/papers/parafrost_gpu.pdf)).
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
       --clean=<target>      remove old installation of <cpu | gpu | all> solvers
       --standard=<n>        compile with <11 | 14 | 17 > c++ standard
       --extra="flags"       pass extra "flags" to the compiler(s)


## GPU solver
To build the GPU solver, make sure you have a CUDA-capable GPU with pre-installed NVIDIA driver and CUDA toolkit.

For installing CUDA v11.7, run the following commands on a Ubuntu 20.04:<br>

`wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/cuda-ubuntu2004.pin`<br>
`sudo mv cuda-ubuntu2004.pin /etc/apt/preferences.d/cuda-repository-pin-600`<br>
`sudo apt-key adv --fetch-keys https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/3bf863cc.pub`<br>
`sudo add-apt-repository "deb https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/ /"`<br>
`sudo apt-get update`<br>
`sudo apt-get -y install cuda`<br>

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

# References
If you are using ParaFROST, please cite the following papers:

```
@inproceedings{Parafrost-OsamaWB-tacas21,
  author    = {Muhammad Osama and
               Anton Wijs and
               Armin Biere},
  editor    = {Jan Friso Groote and
               Kim Guldstrand Larsen},
  title     = {{SAT} Solving with {GPU} Accelerated Inprocessing},
  booktitle = {Tools and Algorithms for the Construction and Analysis of Systems
               - 27th International Conference, {TACAS} 2021, Held as Part of the
               European Joint Conferences on Theory and Practice of Software, {ETAPS}
               2021, Luxembourg City, Luxembourg, March 27 - April 1, 2021, Proceedings,
               Part {I}},
  series    = {Lecture Notes in Computer Science},
  volume    = {12651},
  pages     = {133--151},
  publisher = {Springer},
  year      = {2021},
  url       = {https://doi.org/10.1007/978-3-030-72016-2\_8},
  doi       = {10.1007/978-3-030-72016-2\_8}
}
```
```
@inproceedings{ParalelSimp-OsamaW-tacas19,
  author    = {Muhammad Osama and
               Anton Wijs},
  editor    = {Tom{\'{a}}s Vojnar and
               Lijun Zhang},
  title     = {Parallel {SAT} Simplification on {GPU} Architectures},
  booktitle = {Tools and Algorithms for the Construction and Analysis of Systems
               - 25th International Conference, {TACAS} 2019, Held as Part of the
               European Joint Conferences on Theory and Practice of Software, {ETAPS}
               2019, Prague, Czech Republic, April 6-11, 2019, Proceedings, Part
               {I}},
  series    = {Lecture Notes in Computer Science},
  volume    = {11427},
  pages     = {21--40},
  publisher = {Springer},
  year      = {2019},
  url       = {https://doi.org/10.1007/978-3-030-17462-0\_2},
  doi       = {10.1007/978-3-030-17462-0\_2}
}
```
