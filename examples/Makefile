include ../Makefile.inc

INC = $(TDMA_PATH)/include
LIB = $(TDMA_PATH)/lib -lpascal_tdma

EXE1= ex1_single.out
EXE2= ex2_many_2D.out
EXE3= ex3_many_thread_3D.out
EXE4= ex4_many_rhs_3D.out

topology = mpi_topology_2D.o 

main1 = tdm_single.o
main2 = tdm_many_2D.o
main3 = tdm_many_thread_3D.o
main4 = tdm_many_rhs_3D.o

OBJS1   = $(topology) $(main1)
OBJS2   = $(topology) $(main2)
OBJS3   = $(topology) $(main3)
OBJS4   = $(topology) $(main4)

all:$(OBJS1) $(OBJS2) $(OBJS3) $(OBJS4)
	$(FC) $(FLAG) $(OBJS1) -o $(EXE1) -I$(INC) -L$(LIB)
	$(FC) $(FLAG) $(OBJS2) -o $(EXE2) -I$(INC) -L$(LIB)
	$(FC) $(FLAG) $(OBJS3) -o $(EXE3) -I$(INC) -L$(LIB)
	$(FC) $(FLAG) $(OBJS4) -o $(EXE4) -I$(INC) -L$(LIB)
	cd convection_3D;make all

%.o: %.f90
	$(FC) $(FLAG) -c $^ -o $@  -I$(INC)

%(EXE1)	: %(OBJS1)
%(EXE2)	: %(OBJS2)
%(EXE3)	: %(OBJS3)
%(EXE4)	: %(OBJS4)

clean:
	rm -f *.mod *.o *__genmod.f90 $(EXE1) $(EXE2) $(EXE3) $(EXE4)
	cd convection_3D;make clean

