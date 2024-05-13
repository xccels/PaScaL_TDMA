Performance                         {#perf_page}
===========

[TOC]

# Vectorization


<!-- Vectorization vs. no vectorization -->
<div class="darkmode_inverted_image">\image html perf_vectorization.png width=50%</div>

<!-- 매니코어 프로세서 레벨의 벡터화 적용: 벡터화 적용을 위해 PaScaL_TDMA 라이브러리를 대상으로 2차원 배열변수의 접근 순서 벡터와에 유리하게 변경하여, 의존성이 없는 배열 인덱싱을 연속적으로 수행할 수 있게 배열 재설계 -->
<!-- 벡터화 적용 결과 512^3 격자에서 wall-colock 42% 감소 및 72%의 성능 향상 달성  -->
- A multicore processor-level vectorization application reduces wall-clock time to process a \\(512^3\\) grid by 42% and improves performance by 72%.

- Tested in [CPU environment](perf_page.html#experimental-setup) proposed by [Kim et al. (2021)](reference_page.html)


# CPU Parallel benchmark

@pre 
- Tested in [CPU environment](perf_page.html#experimental-setup) proposed by [Kim et al. (2021)](reference_page.html)

## MPI and OpenMP parallelization

<!-- # MPI and OpenMP parallelization -->
<div class="darkmode_inverted_image">\image html perf_parallel_openmp.png width=50%</div>

<!-- 격자수가 큰 경우(512^3) OpenMP/MPI 병렬 확장성이 모두 우수한 반면, 격자수가 작은 경우(64^3), 병렬성능의 저하가 관찰됨. 전반적으로 MPI 병렬화가 OpenMP 병렬화에 비해 우수함.-->

- The scalability of parallelization, as exemplified by OpenMP and MPI, is optimal when the number of grids is high (\\(512^3\\)). However, when the grid number is low (\\(64^3\\)), the parallel performance is compromised.


## MPI w/ OpenMP hybrid parallelization
<div class="darkmode_inverted_image">\image html perf_parallel_mpi_openmp_hybrid.png width=45%</div>


<!-- 64^3 격자에 대해 총 코어수 64개로 고정하고, MPI+OpenMP hybrid병렬 성능을 측정. 적절한 OpenMP 스레드에 대해 hybrid 병렬화가 우수하며, 스레드 수를 2개로 두었을 때, MPI 병렬화에 비해 34%의 성능 향상을 보임. -->

- The parallel performance of the MPI with OpenMP hybrid approach is done in \\(64^3\\) grids and 64 total cores. 
- Hybrid parallelization represents an effective approach for appropriate OpenMP threads, and when the number of threads is set to two, it shows a 34% performance improvement compared to MPI parallelization.



## MPI: PaScaL_TDMA_single vs. PaScaL_TDMA_many
<!-- ## MPI -->
<!-- Figure 6 -->
<!-- 설명/테스트 장비 -->
<div class="darkmode_inverted_image">\image html perf_parallel_cpu_mpi.png width=80%</div>


- The strong scalability of `PaScaL_TDMA_single` and `PaScaL_TDMA_many` for the tridiagonal systems of the three-dimensional array. 
- The results for the fixed grid size of \\(512^3\\) and \\(2,048^3\\) are plotted in blue dashed and red solid lines, respectively. 
- The measured performance of PaScaL_TDMA_single and PaScaL_TDMA_many are represented as square and triangular symbols, respectively. 
- (a) Computation time, (b) communication time, (c) total execution time, and (d) speedup curve starting from 8 cores.
- Details of the benchmark are presented by [Kim et al. (2021)](reference_page.html)


## 3D Heat Transformer Benchmark 
<div class="darkmode_inverted_image">\image html perf_parallel_cpu_heat.png width=80%</div>
<!-- CPU Setting -->
- The scalability of PaScaL_TDMA_many for solving the three-dimensional heat equation. The dashed lines represent the ideal speedup. (a) Strong scalability and (b) weak scalability. 
- The analysis was based on the total execution time, including not only the time for solving the tridiagonal systems but also the processing time for constructing the coefficient matrices and right-hand sides of the tridiagonal systems.
- Details of the benchmark are presented by [Kim et al. (2021)](reference_page.html)

# GPU Parallel benchmark 

@pre 
- Tested in [GPU environment](perf_page.html#experimental-setup) proposed by [Yang et al. (2023)](reference_page.html)

<div class="darkmode_inverted_image">
\image html perf_parallel_gpu.png width=60%
</div>

- (a) Strong scalability cureves with wall-clock time results as a fuction of the number of CPUs/GPUs and (b) energy consumption results obtained with two KNL CPUs and two A100 GPUs, respectivley. 
- Details of the benchmark are presented by [Yang et al. (2023)](reference_page.html)


# Experimental Setup {#experimental-setup}           
 <div class="tabbed">

- <b class="tab-title">CPU</b>
    @note
    **[CPU Environment]** <br>
    All the computations were executed on **the Nurion manycore cluster** at the Korea Institute of Science and Technology Information (KISTI). The Nurion consists of **8,305 Cray CS500 nodes** interconnected by the Intel Omni-Path Architecture. Each node has a **68-Core Intel Xeon Phi 7250 processor with 16GB** of high bandwidth on a chip and 96 GB of main memory. Each simulation is repeated for 10 time steps to average any small fluctuations in the execution time. An executable program was built using an Intel compiler (version 19.0.1.144) with the flags of full optimization (-O3) and automatic vectorization according to the hardware architecute (-xMIC-AVX512).<br>
    This environment proposed by [Kim et al. (2021)](reference_page.html)
    
- <b class="tab-title">GPU</b>
    @note
    **[GPU Environment]** <br>
    We evaluated the computational performance and energy efficiency of the GPU implementation of PaScaL_TDMA 2.0 on **the NEURON cluster** at the Korea Institute of Science and Technology Information (KISTI). The cluster consisted of **two AMD EPYC 7543 processors (hosts) and eight NVLnk-connected NVIDIA A100 GPUs (devices) per compute node.** The results were compared with those obtained on the NURION cluster at KISTI, which features an Intel Xeon Phi 7750 Knight Landing (KNL) processor per compute node. Intel OneAPI 22.2 and NVIDIA HPC SDK 22.7 were used to compile PaScaL_TDMA 2.0 on the NURION and NEURON clusters, respectively. <br>
    This environment proposed by [Yang et al. (2023)](reference_page.html)
</div>



<div class="section_buttons">

| Previous          |                              Next |
|:------------------|----------------------------------:|
| [Installation](install_page.html) | [Theory](theory_page.html) |
</div>