include $(TDMA_PATH)/Makefile.inc

LIB = libtimer.a

OBJDIR = ./obj

TIMER = $(OBJDIR)/timer.o #$(OBJDIR)/tdmas_cuda.o

OBJS   = $(TIMER)

all: $(LIB) inc

$(LIB): $(OBJS)
	ar qc $@ $(OBJS)
	cp $@ $(TDMA_PATH)/lib

inc: $(LIB)
	cp $(OBJDIR)/*.mod $(TDMA_PATH)/include

$(OBJDIR)/%.o: %.f90
	$(FC) $(FLAG) -c $^ -o $@  -$(opt) $(OBJDIR) -I$(OBJDIR)
$(OBJDIR)/%.o:$(OBJDIR)/%.mod
# $(FC) $(FLAG) -c $^ -o $@  -module $(OBJDIR) -I$(OBJDIR)

%(PTDMA)	 : %(RANGE) %(TDMAS)

clean:
	rm -rf ../lib ../include *__genmod.f90 $(LIB)
	rm -rf $(OBJDIR)

