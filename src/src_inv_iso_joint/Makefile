CMD = DSurfJt
FC = gfortran
FFLAGS  = -O3 -ffixed-line-length-none -ffloat-store\
           -W  -fbounds-check -m64 -mcmodel=medium
F90SRCS = lsmrDataModule.f90 lsmrblasInterface.f90\
          lsmrblas.f90 lsmrModule.f90 delsph.f90\
	  aprod.f90 gaussian.f90 Main_Jt.f90\
	  CalSigamNorm.f90 CalAzimTraveltime.f90 TikhRegul.f90
FSRCS =  surfdisp96.f tregn96_subroutine.f
FMP =   rpathsAzim.f90 FwdAzimuthalAniMap.f90
OBJS = CalSurfG.o $(F90SRCS:%.f90=%.o) $(FSRCS:%.f=%.o) $(FMP:%.f90=%.o) CalSurfGAniso_Joint.o depthkernelTI.o
all:$(CMD)
$(CMD):$(OBJS)
	$(FC) -fopenmp $^ -o $@
	echo "Make successfully!"
	mv $(CMD) ../../bin/

CalSurfG.o: CalSurfG.f90
	$(FC) -fopenmp $(FFLAGS) -c $<  -o $@		
%.o: %.f90
	$(FC) $(FFLAGS) -c $(@F:.o=.f90)  -o $@
%.o: %.f
	$(FC) $(FFLAGS) -c $(@F:.o=.f) -o $@
CalSurfGAniso_Joint.o: CalSurfGAniso_Joint.f90
	$(FC) -fopenmp $(FFLAGS) -c $<  -o $@	
		
depthkernelTI.0: depthkernelTI.f90
	$(FC) $(FFLAGS) -c $<  -o $@

clean:
	rm *.o *.mod ../../bin/$(CMD)
