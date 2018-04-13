FC = nagfor
FFLAGS = -O2
INC = -I.

fd1d_heat_explicit.exe : fd1d_heat_explicit.o
	$(FC) $(FFLAGS) $< -L. -lfmw -o $@

Types_mod.o : Types_mod.f90
	$(FC) $(FFLAGS)  -c $<
RHS_mod.o : RHS_mod.f90 Types_mod.o 
	$(FC) $(FFLAGS) $(INC) -c $<
IO_mod.o : IO_mod.f90 Types_mod.o
	$(FC) $(FFLAGS) $(INC) -c $<
CFL_mod.o : CFL_mod.f90 Types_mod.o
	$(FC) $(FFLAGS) $(INC) -c $<
Solver_mod.o : Solver_mod.f90 RHS_mod.o Types_mod.o
	$(FC) $(FFLAGS) $(INC) -c $<
fd1d_heat_explicit.o : fd1d_heat_explicit.f90 Types_mod.o Solver_mod.o CFL_mod.o RHS_mod.o IO_mod.o
	$(FC) $(FFLAGS) $(INC) -c $<
	ar rcs libfmw.a CFL_mod.o IO_mod.o RHS_mod.o Solver_mod.o Types_mod.o

.PHONY: clean
clean: 
	rm -f *.mod *.o *.exe

