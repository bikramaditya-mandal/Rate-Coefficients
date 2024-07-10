FC=gfortran
OBJ=engines.f90 Generate_Rates.f90
compute_rate: $(OBJ)
	$(FC) -o $@ $(OBJ)

clean:
	rm -f *.o
	rm -f *.mod
	rm -f compute_rate
	rm -f *.out
	rm -f MQCT_TACS.dat


