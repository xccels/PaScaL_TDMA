Installation                        {#install_page}
============

[TOC]
# Downloads
The repository can be cloned as follows:

```
git clone https://github.com/MPMC-Lab/PaScaL_TDMA.git
```
Alternatively, the source files can be downloaded through github menu 'Download ZIP'.




# Compile
## Prerequisites
Prerequisites to compile PaScaL_TDMA are as follows:
* MPI (OpenMPI, IntelMPI, MPICH, .. and so on)
* Fortran compiler (NVIDIA HPC SDK, Intel fortran, GNU fortran, .. and so on )
  * ([`nvfortran`](https://developer.nvidia.com/hpc-sdk-downloads) for GPU runs, NVIDIA HPC SDK 21.1 or higher)

### Option1: NVIDIA HPC SDK + Open MPI
[![](https://img.shields.io/badge/HPC_SDK-21.1_or_higher-nvidia.svg?logo=nvidia)](https://developer.nvidia.com/hpc-sdk-downloads)
  ![](https://img.shields.io/badge/Tested-HPC_SDK_24.3_with_CUDA_12.3-silver.svg?logo=cachet)
* The NVIDIA HPC SDK includes a pre-compiled version of OpenMPI.
* Change `24.3` to `Your SDK version`
  
  ```shell
  cd <your-path>/PaScaL_TDMA
  export TDMA_PATH=$(pwd)
  export NVARCH=`uname -s`_`uname -m`
  export PATH=/opt/nvidia/hpc_sdk/$NVARCH/24.3/compilers/bin:$PATH
  export cuda_version=$(nvcc --version | grep "release" | sed 's/.*release //' | sed 's/,//' | cut -d' ' -f1); export cuda_version
  export MPI=/opt/nvidia/hpc_sdk/$NVARCH/24.3/comm_libs/$cuda_version/openmpi4/latest
  export PATH=$PATH:/opt/nvidia/hpc_sdk/$NVARCH/24.3/compilers/bin
  export opt=module
  ```
  
* Add shortcuts (*Optional*) 
  
  ```shell
  # Also, change `24.3` to `Your SDK version`
  echo 'export NVARCH=`uname -s`_`uname -m`' >> ~/.bashrc
  echo 'alias  nvmpirun="/opt/nvidia/hpc_sdk/$NVARCH/24.3/comm_libs/mpi/bin/mpirun"' >> ~/.bashrc

  source ~/.bashrc
  ```

* Build PaScaL_TDMA 

  * Activate `FC = $(MPI)/bin/mpifort` in `Makefile.inc`

  ```shell
  make all
  ```

* Test examples
  ```shell 
  cd examples
  nvmpirun -np 2 ex1_single.out
  ```

<!-- TODO: Make shell script for configuration -->
_____
### Option2: Intel Fortran + Open MPI
[![](https://img.shields.io/badge/Intel_HPC_Toolkit-Fortran_Compiler-blue.svg?logo=intel)](https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html)
![](https://img.shields.io/badge/Tested-oneAPI_2024.1_+_OpenMPI--5.0.2-silver.svg?logo=cachet)

* Install `Intel oneAPI Base Toolkit` and `Intel HPC Toolkit`
* Install `Open MPI`

  ```shell
  cd ~/PaScaL_TDMA
  export TDMA_PATH=$(pwd)

  # Or use 'source /opt/intel/oneapi/<version>/oneapi-vars.sh'
  source /opt/intel/oneapi/setvars.sh
  export MPI=<your-path>/openmpi
  export opt=J

  ```
* Build PaScaL_TDMA 
  
  * Activate `FC = $(MPI)/bin/mpifort` in `Makefile.inc`

  
  ```shell
  make all
  ```

* Test examples
  ```shell 
  cd examples
  mpirun -np 2 ./ex1_single.out
  ```

_____
### Option3: Intel Fortran + Intel MPI
[![](https://img.shields.io/badge/Intel_HPC_Toolkit-Fortran_Compiler-blue.svg?logo=intel)](https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html)
![](https://img.shields.io/badge/Tested-oneAPI_2024.1-silver.svg?logo=cachet)


* Install `Intel oneAPI Base Toolkit` and `Intel HPC Toolkit`


  ```shell
  cd ~/PaScaL_TDMA
  export TDMA_PATH=$(pwd)

  # Or use 'source /opt/intel/oneapi/<version>/oneapi-vars.sh'
  source /opt/intel/oneapi/setvars.sh
  export opt=module

  ```
* Build PaScaL_TDMA 
  
  * Activate `FC = mpiifort` in `Makefile.inc`

  
  ```shell
  make all
  ```

* Test examples
  ```shell 
  cd examples
  mpirun -np 2 ./ex1_single.out
  ```

___
### Compile and build
* Build PaScaL_TDMA
    ```
	make lib
	```
* Build an example problem after build PaScaL_TDMA

    ```
	make example
	```

* Build all

  ```
	make all
	```

* Clean build
  ```shell
  export TDMA_PATH=$(pwd)
  make clean
  ```


# Mores on compile option
The `Makefile` in root directory is to compile the source code, and is expected to work for most systems. The 'Makefile.inc' file in the root directory can be used to change the compiler (and MPI wrapper) and a few pre-defined compile options depending on compiler, execution environment and et al.

# Running the example
After building the example file, an executable binary, `*.out`, is built in the `run` folder. The `PARA_INPUT.inp` file in the `run` folder is a pre-defined input file, and the `*.out` can be executed as follows:

```
mpirun -np 8 ./a.out ./PARA_INPUT.inp
```
# GPU power monitoring
In the `tool` folder, there is a Python script `gpu_power_monitor.py` that can be used to monitor and print real-time GPU power usage with timestamps. To use this script, you will need to install the `pynvml` library.

# Folder structure
* `src` : source files of PaScaL_TDMA 2.0.
* `example` : source files of an example problem for 3D heat-transfer equation.
* `include` : header files are created after building
* `lib` : a static library of PaScaL_TDMA 2.0 is are created after building
* `doc` : documentation
* `run` : an executable binary file for the example problem is created after building.
* `tool` : contains useful scripts and tools.

<div class="section_buttons">

| Previous          |                              Next |
|:------------------|----------------------------------:|
| [Introduction](index.html) | [Performance](perf_page.html) |
</div>