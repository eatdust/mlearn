FC=gfortran
FCFLAGS=-g
#-DDOUBLEFANN -DFIXEDFANN
LFLAGS= -lfann 
#-ldoublefann -lfixedfann

example: fann.o example.o
	$(FC) $(FCFLAGS) fann.o example.o -o $@ $(LFLAGS)

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

%.o: %.f03
	$(FC) $(FCFLAGS) -c $<

%.o: %.F03
	$(FC) $(FCFLAGS) -c $<

clean:
	rm -f *.o *.mod example

