FC = nagfor
FFLAGS = -O2
INC = -I. -I/usr/local/netcdf-4.6.1/include -I/usr/local/plplot-5.13.0/lib/fortran/modules/plplot
LIB = -L/usr/local/netcdf-4.6.1/lib -lnetcdff -lnetcdf -L/usr/local/plplot-5.13.0/lib -lplplotfortran -lplplot

fd1d_heat_explicit.exe : fd1d_heat_explicit.o libfmw.a
	$(FC) $< -L. -lfmw  $(LIB) -o $@

libfmw.a : CFL_mod.o IO_mod.o RHS_mod.o Solver_mod.o Types_mod.o 
	ar rcs $@ $^ 
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
fd1d_heat_explicit.o : fd1d_heat_explicit.f90 libfmw.a
	$(FC) $(FFLAGS) $(INC) -c $<

.PHONY: clean
clean: 
	rm -f *.mod *.o *.exe

