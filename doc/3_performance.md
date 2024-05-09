Performance                         {#perf_page}
===========

[TOC]



# Vectorization
<!-- raw data 
COMPHY-D-19 Figure 8 양식과 동일하게  -->

- All the computations were executed on the Nurion manycore cluster at the Korea Institute of Science and Technology Information (KISTI). The Nurion consists of 8,305 Cray CS500 nodes interconnected by the Intel Omni-Path Architecture. Each node has a 68-Core Intel Xeon Phi 7250 processor with 16GB of high bandwidth on a chip and 96 GB of main memory. Each simulation is repeated for 10 time steps to average any small fluctuations in the execution time. An executable program was built using an Intel compiler (version 19.0.1.144) with the flags of full optimization (-O3) and automatic vectorization according to the hardware architecute (-xMIC-AVX512).

# Parallel benchmark
## Open MP
<!-- graph - raw data -->

## MPI
<!-- Figure 6 -->
<!-- 설명/테스트 장비 -->
<div class="darkmode_inverted_image">\image html perf_parallel_cpu_mpi.png width=80%</div>

The strong scalability of PaScaL_TDMA_single and PaScaL_TDMA_many for the tridiagonal systems of the three-dimensional array. The results for the fixed grid size of \\(512^3\\) and \\(2,048^3\\) are plotted in blue dashed and red solid lines, respectively. The measured performance of PaScaL_TDMA_single and PaScaL_TDMA_many are represented as square and triangular symbols, respectively. (a) Computation time, (b) communication time, (c) total execution time, and (d) speedup curve starting from 8 cores.

## 3D Heat Transformer Benchmark 
<div class="darkmode_inverted_image">\image html perf_parallel_cpu_heat.png width=80%</div>

The scalability of PaScaL_TDMA_many for solving the three-dimensional heat equation. The dashed lines represent the ideal speedup. (a) Strong scalability and (b) weak scalability.

## GPU benchmark 

We evaluated the computational performance and energy efficiency of the GPU implementation of PaScaL_TDMA 2.0 on the NEURON cluster at the Korea Institute of Science and Technology Information (KISTI). The cluster consisted of two AMD EPYC 7543 processors (hosts) and eight NVLnk-connected NVIDIA A100 GPUs (devices) per compute node. The results were compared with those obtained on the NURION cluster at KISTI, which features an Intel Xeon Phi 7750 Knight Landing (KNL) processor per compute node. Intel OneAPI 22.2 and NVIDIA HPC SDK 22.7 were used to compile PaScaL_TDMA 2.0 on the NURION and NEURON clusters, respectively. 
<!-- 논문 그림 복사e -->
<!-- 논문 해당 내용 1-2줄 -->
<div class="darkmode_inverted_image">
\image html perf_parallel_gpu.png width=60%
</div>
(a) Strong scalability cureves with wall-clock time results as a fuction of the number of CPUs/GPUs and (b) energy consumption results obtained with two KNL CPUs and two A100 GPUs, respectivley. The sonsumed energies by the host CPUs in the GPU cluster, AMD EPYC, are also plotted along with the results measured with A100 GPUs (For interpretation of the colors in the figure.


<div class="section_buttons">

| Previous          |                              Next |
|:------------------|----------------------------------:|
| [Installation](install_page.html) | [Theory](theory_page.html) |
</div>