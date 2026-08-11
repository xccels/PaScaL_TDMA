# =============================================================================
# PaScaL_TDMA - Python C-API wrapper library
#
#   Target: lib/libpascal_tdma_capi.so
#           Shared library bundling the CUDA solver (PaScaL_TDMA_cuda) and the
#           bind(C) wrapper, called from Python via ctypes with CuPy arrays.
#
#   Compiler and flags are reused from 01_Fortran/Makefile.inc (FC = mpifort).
#
#   Usage:
#     make                 # build lib/libpascal_tdma_capi.so
#     make GPU_CC=cc90     # change target GPU architectures (V100=cc70, A100=cc80, H100=cc90)
#     make clean
# =============================================================================
TDMA_PATH ?= $(CURDIR)/../01_Fortran
include $(TDMA_PATH)/Makefile.inc

# nvfortran 25.x removed -Mcuda: replace it with -cuda and select GPU targets.
# Default is a V100+A100+H100 fat binary. The architecture of the GPU that runs
# the code must be included in the target list: with nprocs >= 2 the solver
# executes kernel paths that fail with SIGILL when only PTX-JITed from another
# target.
GPU_CC ?= cc70,cc80,cc90
FLAG   := $(subst -Mcuda,-cuda,$(FLAG)) -gpu=$(GPU_CC)

SRC    = $(TDMA_PATH)/src
OBJDIR = ./obj
LIBDIR = ./lib
TARGET = $(LIBDIR)/libpascal_tdma_capi.so

# Module dependency order: para_range -> tdmas_cuda -> PaScaL_TDMA_cuda -> capi wrapper
OBJS = $(OBJDIR)/para_range.o \
       $(OBJDIR)/tdmas_cuda.o \
       $(OBJDIR)/pascal_tdma_cuda.o \
       $(OBJDIR)/pascal_tdma_capi.o

all: $(TARGET)

$(TARGET): $(OBJS)
	mkdir -p $(LIBDIR)
	$(FC) $(FLAG) -shared $(OBJS) -o $(TARGET)

# Library sources (01_Fortran/src)
$(OBJDIR)/%.o: $(SRC)/%.f90
	mkdir -p $(OBJDIR)
	$(FC) $(FLAG) -c $< -o $@ -module $(OBJDIR) -I$(OBJDIR)

# Wrapper source (02_python/src)
$(OBJDIR)/%.o: src/%.f90
	mkdir -p $(OBJDIR)
	$(FC) $(FLAG) -c $< -o $@ -module $(OBJDIR) -I$(OBJDIR)

clean:
	rm -rf $(OBJDIR) $(TARGET) *.mod *__genmod.f90
