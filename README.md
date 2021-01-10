# ParaFROST
ParaFROST stands for Parallel Formal Reasoning of Satisfiability. It is a parallel SAT solver with GPU-accelerated inprocessing capable of harnessing NIVIDA CUDA-enabled GPUs in applying modern inprocessing tecnhiques in parallel. The CDCL search is built from scratch with various optimisations based on CaDiCaL heuristics (see our paper: ). The inprocessing engine extends our previous work in with new data structures, parallel garbage collection and more.

Run "make -C gpu" to build the gpu solver or "make -C cpu" to build a CPU-only version.

Run the "artefact" to test and compare the performance of ParaFROST (GPU) with its sequential counterpart and CaDiCaL.
