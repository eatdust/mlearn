ext=$(shell uname | cut -c1-3)

ifeq ($(ext),Lin)
FC=gfortran
#FCFLAGS= -O3 -fopenmp -ffast-math -march=haswell  -fprefetch-loop-arrays -funroll-loops
FCFLAGS=-O -DDOUBLEFANN
LFLAGS= -ldoublefann 
#FC=ifort
#FCFLAGS = -g -openmp -mkl -assume byterecl -free
#LFLAGS = -lfann
#CPP=cpp -P -traditional-cpp 
endif 

OBJSCOM = prec.o ioml.o

OBJSFANN = fnnprec.o fann.o fnn.o

fnnmain.$(ext): $(OBJSCOM) $(OBJSFANN) fnnmain.o
	$(FC) $(FCFLAGS) $(OBJSCOM) $(OBJSFANN) fnnmain.o -o $@ $(LFLAGS)

testfnn.$(ext): $(OBJSCOM) $(OBJSFANN) testfnn.o
	$(FC) $(FCFLAGS) $(OBJSCOM) $(OBJSFANN) testfnn.o -o $@ $(LFLAGS)

fnd.o: fann.o

%.o: %.f90
	$(FC) $(FCFLAGS) -c $*.f90

%.o: %.F90
	$(FC) $(FCFLAGS) -c $*.F90

%.o: %.f95
	$(FC) $(FCFLAGS) -c $*.f95

%.o: %.F95
	$(FC) $(FCFLAGS) -c $*.F95

%.o: ../common/%.f90
	$(FC) $(FCFLAGS) -c $<

%.o: %.c
	$(CC) $(CFLAGS) -c $<

ifeq ($(FC),ifort)
%.ifort: %.F03
	$(CPP) $< -o $@
%.ifort: %.f03
	$(CPP) $< -o $@
%.o: %.ifort
	$(FC) $(FCFLAGS) -c -free -Tf$< 

else
%.o: %.F03
	$(FC) $(FCFLAGS) -c $<

%.o: %.f03
	$(FC) $(FCFLAGS) -c $<

endif

clean:
	rm -f *.ifort *.o *.mod *.Lin

