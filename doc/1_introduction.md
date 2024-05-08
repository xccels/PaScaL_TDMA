Introduction                       {#intro_page}
============

[TOC]
# Overview
Parallel and Scalable Library for TriDiagonal Matrix Algorithm

PaScaL_TDMA provides an efficient and scalable computational procedure to solve many tridiagonal systems in multi-dimensional partial differential equations. The modified Thomas algorithm proposed by [Laszlo et al.(2016)](reference_page.html) and the newly designed communication scheme have been used to reduce the communication overhead in solving many tridiagonal systems.

This library is for both single and many tridiagonal systems of equations. The main algorithm for a tridiagonal matrix consists of the following five steps: 

- (1) Transform the partitioned submatrices in the tridiagonal systems into modified submatrices:
        Each computing core transforms the partitioned submatrices in the tridiagonal systems of equations into the modified forms by applying modified Thomas algorithm.
- (2) Construct reduced tridiagonal systems from the modified submatrices:
        The reduced tridiagonal systems are constructed by collecting the first and last rows of the modified submatrices from each core using MPI_Ialltoallw.
- (3) Solve the reduced tridiagonal systems:
        The reduced tridiagonal systems constructed in Step 2 are solved by applying the Thomas algorithm.
- (4) Distribute the solutions of the reduced tridiagonal systems:
        The solutions of the reduced tridiagonal systems in Step 3 are distributed to each core using MPI_Ialltoallw.
        This communication is an exact inverse of the communication in Step 2.
- (5) Update the other unknowns in the modified tridiagonal systems:
        The remaining unknowns in the modified submatrices in Step 1 are solved in each computing core with the solutions obtained in Step 3 and Step 4.
    
Step 1 and Step 5 are similar to the method proposed by [Laszlo et al.(2016)](reference_page.html) which uses parallel cyclic reduction (PCR) algorithm to build and solve the reduced tridiagonal systems. Instead of using the PCR, we develop an all-to-all communication scheme using the MPI_Ialltoall function after the modified Thomas algorithm is executed. The number of coefficients for the reduced tridiagonal systems are greatly reduced, so we can avoid the communication bandwidth problem, which is a main bottleneck for all-to-all communications. Our algorithm is also distinguished from the work of [Mattor et al. (1995)](reference_page.html) which assembles the undetermined coefficients of the temporary solutions in a single processor using MPI_Gather, where load imbalances are serious.


# Features

## Single and many tridiagonal system
### Single tridiagonal system
<div class="darkmode_inverted_image">\image html intro_2.png width=50%</div>
- A communication method using MPI_Gather for a single diagonal system.
- \\(\texttt{PaScaL_TDMA_plan_single_create}\\): A subroutine to create a type variable of \\(\texttt{plan_single}\\), for a single tridiagonal system

- \\(\texttt{PaScaL_TDMA_plan_single_destroy}\\): A subroutine to dellocate the \\(\texttt{plan_single}\\) and all associated data

- \\(\texttt{PaScaL_TDMA_plan_single_solve}\\): A subroutine for solving a single tridiagonal system with the \\(\texttt{MPI_gather}\\) and \\(\texttt{MPI_scatter}\\) fuctions, as shown in upper figure

### Many tridiagonal system
<div class="darkmode_inverted_image">\image html intro_1.png width=90%</div>
- The present method using MPI_Alltoall for many tridiagonal systems. 
- \\(\texttt{PaScaL_TDMA_plan_many_create}\\):
  A subroutine to create a type variable of \\(\texttt{plan_many}\\), for a many tridiagonal systems. The definition of the derived data type depends on the version of MPI library.

- \\(\texttt{PaScaL_TDMA_plan_many_destroy}\\): A subroutine to dellocate the \\(\texttt{plan_many}\\) and all associated data

- \\(\texttt{PaScaL_TDMA_plan_many_solve}\\): A subroutine for solving a many tridiagonal systems with the \\(\texttt{MPI_Alltoallw}\\) fuction, as shown in upper figure
  

## MPI parallelization
<div class="darkmode_inverted_image">\image html intro_3.png width=90%</div>
- Distributed memory systems is required for large scale problem
- The parallel computation of implicit method in distributed memory systems requires massive data communications
- This leads to a major bottleneck in a large scale cluster system with distributed memories

## GPU implementation
In PaScaL_TDMA 2.0, multi-GPU acceleration is implemented using NVIDIA CUDA. CUDA-related features are as follows:
- (1) Incorporation of CUDA kernels into the loop structures of the existing algorithm, that are modified to exploit more GPU threads.
- (2) Utilization of shared memory using pipeline copy of variables in device memory to reduce the amount of device memory access.
- (3) CUDA-aware MPI communication for rapid communication with the support of hardward
- (4) Use of 3D-array for practical applications and accordingly the use of CUDA threads more than in a single dimension with the 3-D array.
- (5) Depreciation on functions for single tridiagonal matrix as they are rarely used for three-dimensional problems.

# Authors

- Ki-Ha Kim (k-kiha@yonsei.ac.kr), School of Mathematics and Computing (Computational Science and Engineering), Yonsei University (v1.0, v2.0)
- Mingyu Yang (yang926@yonsei.ac.kr), School of Mathematics and Computing (Computational Science and Engineering), Yonsei University (v2.0)
- Ji-Hoon Kang (jhkang@kisti.re.kr), Korea Institute of Science and Technology Information (v1.0, v2.0)
- Oh-Kyoung Kwon (okkwon@kisti.re.kr), Korea Institute of Science and Technology Information (v2.0)
- Jung-Il Choi (jic@yonsei.ac.kr), School of Mathematics and Computing (Computational Science and Engineering), Yonsei University (v1.0, v2.0)

# Versions
<div class="tabbed">

- <b class="tab-title">1.0</b>
    - MPI parallelization proposed by [Kim et al. (2021)](reference_page.html)
    
- <b class="tab-title">2.0</b>
    - GPU implementation proposed by [Yang et al. (2023)](reference_page.html)

</div>


# Citation

Please use the following bibtex, when you refer to this project.

    @article{ykkkc2023,
        title = "PaScaL_TDMA 2.0: A multi-GPU-based library for solving massive tridiagonal systems",
        author = "Yang, Mingyu Yang and Kang, Ji-Hoon Kang and Kim, Ki-Ha Kim and Kwon, Oh-Kyoung and Choi, Jung-Il",
        journal = "Computer Physics Communications",
        volume = "290",
        pages = "108779",
        year = "2023",
        issn = "0010-4655",
        doi = "https://doi.org/10.1016/j.cpc.2023.108785"
    }
    @article{kkpc2020,
        title = "PaScaL_TDMA: A library of parallel and scalable solvers for massive tridiagonal system",
        author = "Kim, Ki-Ha and Kang, Ji-Hoon and Pan, Xiaomin and Choi, Jung-Il",
        journal = "Computer Physics Communications",
        volume = "260",
        pages = "107722",
        year = "2021",
        issn = "0010-4655",
        doi = "https://doi.org/10.1016/j.cpc.2020.107722"
    }

    @misc{PaScaL_TDMA2019,
        title  = "Parallel and Scalable Library for TriDiagonal Matrix Algorithm",
        author = "Kim, Ki-Ha and Kang, Ji-Hoon and Choi, Jung-Il",
        url    = "https://github.com/MPMC-Lab/PaScaL_TDMA",
        year   = "2019"
    }

## Credits

Thanks for writing this library and inspiring feedback on GitHub!

Special thanks to all the contributors:
<br><br>
<a href="https://github.com/xccels/PaScaL_TDMA/graphs/contributors">
    <img src="https://contrib.rocks/image?repo=xccels/PaScaL_TDMA" />
</a>


<div class="section_buttons">

|                        Read Next |
|---------------------------------:|
| [Installation](install_page.html) |

</div>