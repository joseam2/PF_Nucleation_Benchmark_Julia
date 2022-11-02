# PF_Nucleation_Benchmark_Julia

## Introduction and Motivation
This repository represents a phase field benchmark described in PF Hub: https://pages.nist.gov/pfhub/benchmarks/benchmark8.ipynb/

The paper that was written describing this repository is published in the IMMI journal and is available here: https://arxiv.org/abs/2210.16866

## Methods
This benchmark was solved and explored using Julia code with CUDA parallelization and finite difference solvers.
The code and data are stored in this repository.

## Code and Data
curv_nuc.jl is the pf code used to develop all the raw data found in this repository.
Each part has a folder that corresponds to the data created for that part.
In parts 2 and 3 there are additional folders where the data for different values of (inital radius / critical radius) are stored.
There is also an additional folder labelled phi_domains which has the data for the full phi array at the time labelled in the title of the file.

## Results
[Part 1 data results](part1_volfrac.png)
[Part 2 data results](2dp2.png)
[Part 3 data results](2dp3.png)
