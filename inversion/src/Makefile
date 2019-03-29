CMD = DSurfTomoAA_ITE
FC = gfortran
FFLAGS  = -O3 -ffixed-line-length-none -ffloat-store\
           -W  -fbounds-check -m64 -mcmodel=medium
F90SRCS = lsmrDataModule.f90 lsmrblasInterface.f90\
          lsmrblas.f90 lsmrModule.f90 delsph.f90\
	  aprod.f90 gaussian.f90 MainLS_iter.f90 GaussianLS.f90\
	  CalSigamNorm.f90 CalAzimTraveltime.f90 TikhRegul.f90
FSRCS =  surfdisp96.f
FMP = CalSurfG.f90 CalSurfGAniso.f90 rpathsAzim.f90 FwdAzimuthalAniMap.f90
OBJS =  $(F90SRCS:%.f90=%.o) $(FSRCS:%.f=%.o) $(FMP:%.f90=%.o)
all:$(CMD)
$(CMD):$(OBJS)
	$(FC) -fopenmp $^ -o $@
	echo "Make successfully!"
%.o: %.f90
	$(FC) -fopenmp $(FFLAGS) -c $(@F:.o=.f90)  -o $@
%.o: %.f
	$(FC) $(FFLAGS) -c $(@F:.o=.f) -o $@
#%.o: %.f90
#	$(FC) $(FFLAGS) -c $(@F:.o=.f90) -o $@
clean:
	rm *.o *.mod $(CMD)
