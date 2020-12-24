# ParaFROST
ParaFROST stands for Parallel Formal Reasoning of Satisfiability. It is parallel SAT solver with GPU-accelerated inprocessing capable of harnessing NIVIDA CUDA-enabled GPUs in applying modern inprocessing tecnhiques in parallel. The CDCL search is built from scratch based on CaDiCaL heuristics (see our paper: ). The inprocessing engine extends our previous work in __, __ with parallel garbage cllection and more.

Run "make -C gpu" to build the gpu solver or "make -C cpu" to build a CPU-only version.

Run the "artefact" to generate the reported graphs in __.
