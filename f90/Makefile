
# set compile environment
# ----------------------------------------------------------------------
FC = "gfortran-9"

# determine compile options from environment set above
# ----------------------------------------------------------------------

ifeq ($(FC),ifort)
	COMP_SWITCH := -Difort
	#BUG := -g -O0 -traceback -fpe:0 -debug all -C -fp-model source -assume minus0
	BUG := -g -O2 -traceback -fpe:0 -debug all -xSSE4.2 -ip -fp-model source -align -assume minus0

	MKL_INC = /opt/intel/Compiler/11.1/073/mkl/include
	MKL_LIB = /opt/intel/Compiler/11.1/073/mkl/lib/em64t
  LAPACK = -L$(MKL_LIB) -I$(MKL_INC) -lmkl_lapack95_lp64 $(MKL_LIB)/libmkl_intel_lp64.a $(MKL_LIB)/libmkl_sequential.a $(MKL_LIB)/libmkl_core.a $(MKL_LIB)/libmkl_sequential.a $(MKL_LIB)/libmkl_core.a $(MKL_LIB)/libmkl_sequential.a $(MKL_LIB)/libmkl_core.a -lpthread

endif
ifeq ($(FC),"gfortran-9")

  # mac debug:
  # lldb <executable>
  # >> process launch -i <file>
  # where <file> is std/carrot pipe in

	COMP_SWITCH :=
	BUG := -g -O0 -m64 -fbounds-check -ffpe-trap=invalid,overflow
	#BUG := -gdwarf-2 -O2 -m64 -funsafe-math-optimizations
	#BUG := -gdwarf-2 -O2 -m64 -ffpe-trap=invalid,overflow

	# for gfortran on 10.9 macs
	#LAPACK = -L/System/Library/Frameworks/Accelerate.framework/Versions/Current/Frameworks/vecLib.framework/Versions/Current/ -lBLAS -lLAPACK
	#LAPACK = -framework Accelerate

	#
	MKL_INC = /usr/local/opt/lapack/include
	MKL_LIB = /usr/local/opt/lapack/lib
  LAPACK = -L$(MKL_LIB) -I$(MKL_INC) -lblas -llapack

endif


SIA2_OBJ := sia2_constants.o kei_hacks.o sia2_parameters.o sia2_state.o \
    sia2_types.o sia2_grid.o sia2_flux_heat.o sia2_desalination.o


OBJ1D  := sia2_constants.o kei_hacks.o sia2_parameters.o sia2_state.o \
  sia2_types.o sia2_grid.o sia2_flux_heat.o sia2_desalination.o\
  kei_ocncommon.o \
  kei_icecommon.o kei_ecocommon.o \
  kei_ice.o kei_eco.o kei_init.o kei_fluxes.o \
  kei_subs1D.o kei_atm.o kei_ocn.o kei_kpp.o

F901D  := kei_parameters.f90 kei_common.f90 kei_link.f90


FFLAGS :=  $(BUG)


# compile directives
# ----------------------------------------------------------------------


kei: $(OBJ1D)
	f2py -c --f90exec="gfortran-9" --f90flags="$(FFLAGS) $(LAPACK)" $(OBJ1D) $(F901D) -m kei

# compile routines
kei_fluxes.o: kei_hacks.o
	$(FC) -c $(FFLAGS) kei_fluxes.f90
kei_ocn.o	:
	$(FC) -c $(FFLAGS) kei_ocn.f90
kei_kpp.o   :
	$(FC) -c $(FFLAGS) kei_kpp.f90
kei_subs1D.o:
	$(FC) -c $(FFLAGS) kei_subs1D.f90
kei_atm.o	: kei_atm.f90
	$(FC) -c $(FFLAGS) kei_atm.f90
kei_init.o	: kei_ecocommon.o kei_parameters.o kei_common.o kei_icecommon.o
	$(FC) -c $(FFLAGS) kei_init.f90

# compile modules
kei_link.o	: $(OBJ1D) kei_parameters.o kei_common.o
	$(FC) -c $(FFLAGS) kei_link.f90 $(LAPACK)
kei_ice.o	: kei_ice.f90 kei_hacks.o
	$(FC) -c $(FFLAGS) kei_ice.f90
kei_eco.o	: kei_ecocommon.o kei_parameters.o kei_common.o kei_hacks.o
	$(FC) -c $(FFLAGS) kei_eco.f90
kei_ecocommon.o	:
	$(FC) -c $(FFLAGS) kei_ecocommon.f90
kei_icecommon.o	: kei_icecommon.f90 $(SIA2_OBJ)
	$(FC) -c $(FFLAGS) kei_icecommon.f90
kei_ocncommon.o	:
	$(FC) -c $(FFLAGS) kei_ocncommon.f90
kei_common.o	: kei_icecommon.o
	$(FC) -c $(FFLAGS) kei_common.f90
kei_parameters.o	:
	$(FC) -c $(FFLAGS) kei_parameters.f90

sia2_flux_heat.o	: sia2_constants.o sia2_parameters.o sia2_types.o
	$(FC) -c $(FFLAGS) sia2_flux_heat.f90
sia2_grid.o	: sia2_constants.o sia2_parameters.o sia2_types.o sia2_desalination.o
	$(FC) -c $(FFLAGS) sia2_grid.f90
sia2_desalination.o	: sia2_constants.o sia2_parameters.o sia2_state.o sia2_types.o
	$(FC) -c $(FFLAGS) sia2_desalination.f90
sia2_state.o	: sia2_parameters.o sia2_types.o
	$(FC) -c $(FFLAGS) sia2_state.f90
sia2_types.o	: sia2_constants.o sia2_parameters.o
	$(FC) -c $(FFLAGS) sia2_types.f90
sia2_parameters.o	: sia2_constants.o
	$(FC) -c $(FFLAGS) sia2_parameters.f90
kei_hacks.o	: kei_parameters.o kei_common.o kei_subs1D.o kei_ecocommon.o
	$(FC) -c $(FFLAGS) kei_hacks.f90
sia2_constants.o	:
	$(FC) -c $(FFLAGS) sia2_constants.f90



#	-------
#	CLEANUP
#	-------
clean :
	rm -f *.o
	rm -f *.mod
	rm -f *.so
	rm -f -r KEI.run.dSYM

# end
