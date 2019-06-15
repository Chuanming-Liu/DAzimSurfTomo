CMD = SurfAAForward
FC = gfortran
FFLAGS  = -O3 -ffixed-line-length-none -ffloat-store\
           -W  -fbounds-check -m64 -mcmodel=medium
F90SRCS = delsph.f90 gaussian.f90 MainForward.f90\
	    FwdAzimuthalAniMap.f90   	   
FSRCS =  surfdisp96.f tregn96_subroutine.f
FMP = rpathsAzim.f90 
OBJS =  CalSurfG.o $(F90SRCS:%.f90=%.o) $(FSRCS:%.f=%.o) $(FMP:%.f90=%.o) FwdTraveltimeCPS.o depthkernelTI.o
all:$(CMD)
$(CMD):$(OBJS)
	$(FC) -fopenmp  $^ -o $@
	echo "Make SurfAzimForward successfully!"
	mv ./SurfAAForward ../../bin
CalSurfG.o: CalSurfG.f90
	$(FC) -fopenmp $(FFLAGS) -c $<  -o $@	

%.o: %.f90
	$(FC) $(FFLAGS) -c $(@F:.o=.f90) -o $@
%.o: %.f
	$(FC) $(FFLAGS) -c $(@F:.o=.f) -o $@
FwdTraveltimeCPS.o: FwdTraveltimeCPS.f90
	$(FC) -fopenmp $(FFLAGS) -c $<  -o $@
		
depthkernelTI.o: depthkernelTI.f90
	$(FC) $(FFLAGS) -c $<  -o $@

	
clean:
	rm *.o *.mod ../../bin/$(CMD)
