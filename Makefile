# Compiler and complilation flag are included in Makefile.inc
include Makefile.inc

lib:
	mkdir -p include;mkdir -p lib
	cd timer; mkdir -p obj; make all
	cd src; mkdir -p obj; make all

example:
	cd examples; make all

all:
	mkdir -p include; mkdir -p lib
	cd timer; mkdir -p obj; make all
	cd src; mkdir -p obj; make all
	cd examples; make all

clean:
	cd timer; make clean
	cd src; make clean
	cd examples; make clean
