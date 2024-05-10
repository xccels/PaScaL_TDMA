Performance                         {#perf_page}
===========

[TOC]

# Vectorization

@pre 
- Tested in [CPU environment](perf_page.html#experimental-setup) proposed by [Kim et al. (2021)](reference_page.html)

<!-- Vectorization vs. no vectorization -->
<div class="darkmode_inverted_image">\image html perf_vectorization.png width=50%</div>

<!-- grid가 증가할 수록, AVX Enabled의 효율이 증가한다 ~ 수치적인 값(00% 감소 등..) -->

The figure presents that the efficiency of AVX Enabled increases as the grid dimension increases.

# CPU Parallel benchmark

@pre 
- Tested in [CPU environment](perf_page.html#experimental-setup) proposed by [Kim et al. (2021)](reference_page.html)

## MPI and OpenMP parallelization

<!-- # MPI and OpenMP parallelization -->
<div class="darkmode_inverted_image">\image html perf_parallel_openmp.png width=50%</div>

<!-- MPI를 사용할 경우 OpenMP보다 효율이 좋다. core(2^0, 2^6)에 따른 Wall-clock time 변화 양상 -->


## MPI w/ OpenMP hybrid parallelization
<div class="darkmode_inverted_image">\image html perf_parallel_mpi_openmp_hybrid.png width=45%</div>

<!-- Total number of cores가 64개인 경우, theread 수에 따른 효율을 나타낸 그래프. 하드웨어의 제한조건 상, threads 배열이 1, 8, 64개인 경우 비효율적임을 확인할 수 있음. -->



## MPI: PaScaL_TDMA_single vs. PaScaL_TDMA_many
<!-- ## MPI -->
<!-- Figure 6 -->
<!-- 설명/테스트 장비 -->
<div class="darkmode_inverted_image">\image html perf_parallel_cpu_mpi.png width=80%</div>

The strong scalability of `PaScaL_TDMA_single` and `PaScaL_TDMA_many` for the tridiagonal systems of the three-dimensional array. The results for the fixed grid size of \\(512^3\\) and \\(2,048^3\\) are plotted in blue dashed and red solid lines, respectively. The measured performance of PaScaL_TDMA_single and PaScaL_TDMA_many are represented as square and triangular symbols, respectively. (a) Computation time, (b) communication time, (c) total execution time, and (d) speedup curve starting from 8 cores.

- Details of the benchmark are presented by [Kim et al. (2021)](reference_page.html)



## 3D Heat Transformer Benchmark 
<div class="darkmode_inverted_image">\image html perf_parallel_cpu_heat.png width=80%</div>
<!-- CPU Setting -->
The scalability of PaScaL_TDMA_many for solving the three-dimensional heat equation. The dashed lines represent the ideal speedup. (a) Strong scalability and (b) weak scalability.

- We apply the numerical method developed by Pan et al. to solve the heat equation without considering the convection term.
- This results for up to 262,144 cores were obtained. The analysis was based on the total execution time, including not only the time for solving the tridiagonal systems but also the processing time for constructing the coefficient matrices and right-hand sides of the tridiagonal systems.

- Details of the benchmark are presented by [Kim et al. (2021)](reference_page.html)

# GPU Parallel benchmark 

@pre 
- Tested in [GPU environment](perf_page.html#experimental-setup) proposed by [Yang et al. (2023)](reference_page.html)

<div class="darkmode_inverted_image">
\image html perf_parallel_gpu.png width=60%
</div>

(a) Strong scalability cureves with wall-clock time results as a fuction of the number of CPUs/GPUs and (b) energy consumption results obtained with two KNL CPUs and two A100 GPUs, respectivley. 


- Details of the benchmark are presented by [Yang et al. (2023)](reference_page.html)


# Experimental Setup {#experimental-setup}           
 <div class="tabbed">

- <b class="tab-title">CPU</b>
    @note
    **[CPU Environment] proposed by [Kim et al. (2021)](reference_page.html)** <br>
    All the computations were executed on the Nurion manycore cluster at the Korea Institute of Science and Technology Information (KISTI). The Nurion consists of 8,305 Cray CS500 nodes interconnected by the Intel Omni-Path Architecture. Each node has a 68-Core Intel Xeon Phi 7250 processor with 16GB of high bandwidth on a chip and 96 GB of main memory. Each simulation is repeated for 10 time steps to average any small fluctuations in the execution time. An executable program was built using an Intel compiler (version 19.0.1.144) with the flags of full optimization (-O3) and automatic vectorization according to the hardware architecute (-xMIC-AVX512).<br>
    
- <b class="tab-title">GPU</b>
    @note
    **[GPU Environment] proposed by [Yang et al. (2023)](reference_page.html)** <br>
    We evaluated the computational performance and energy efficiency of the GPU implementation of PaScaL_TDMA 2.0 on the NEURON cluster at the Korea Institute of Science and Technology Information (KISTI). The cluster consisted of two AMD EPYC 7543 processors (hosts) and eight NVLnk-connected NVIDIA A100 GPUs (devices) per compute node. The results were compared with those obtained on the NURION cluster at KISTI, which features an Intel Xeon Phi 7750 Knight Landing (KNL) processor per compute node. Intel OneAPI 22.2 and NVIDIA HPC SDK 22.7 were used to compile PaScaL_TDMA 2.0 on the NURION and NEURON clusters, respectively. 
</div>



<div class="section_buttons">

| Previous          |                              Next |
|:------------------|----------------------------------:|
| [Installation](install_page.html) | [Theory](theory_page.html) |
</div>