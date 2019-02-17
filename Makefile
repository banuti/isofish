FC =     mpif90
#FC = gfortran-4.9
.SUFFIXES:
.SUFFIXES: .f .f90 .F90 .o

# FC =     mpiifort
# LDFLAGS = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core \
          ${MKLROOT}/../compiler/lib/intel64/libiomp5.so \
          ${MKLROOT}/../compiler/lib/intel64/libimf.so
# INCLUDE = -I${MKLROOT}/include

# for mpif90
FFLAGS = -fno-align-commons -O3


EXE = isofish

# ----- Libs sources
SOURCES = flowprops.f90 \
					modoutput.f90 \
					modanalysis.f90 \
					boundaries.f90 \
					modconvect.f90 \
					moddiffuse.f90 \
					modmass.f90 \
					modconfinement.f90 \
					poisson.f \
					isofish.f90


OBJECT = $(SOURCES:.f90=.o)

$(EXE): $(OBJECT)
	$(FC) $(FFLAGS) -o $@ $^


#
#
# # ----- General rules
# $(MKL_OBJ): $(MKL)
# 	$(FC) $(FFLAGS) -c $<
%.o: %.f90
	$(FC) $(FFLAGS) -c $<

# ----- Utility targets
.PHONY: clean print
clean:
	@\rm -f *.o *.mod *.x *.out
print-%  :
	@echo $* = $($*)
