# makefile for creating G matrix
MOD_DIR	= MODULE/
O_DIR= ../      
FC	= gfortran
# -C:link -c: compile
FFLAGS	= -ffixed-line-length-132 -O4 -c -g
PROG	= $(O_DIR)cbasic_module.o $(O_DIR)nrtype.o $(O_DIR)nrutil.o $(O_DIR)sort_module.o $(O_DIR)interpol_module.o $(O_DIR)ModeNorm_module.o buildG_dcdm1D_module_TI.o


OBJS	= $(MOD_DIR)cbasic_module.f90 $(MOD_DIR)nrtype.f90 $(MOD_DIR)nrutil.f90 $(MOD_DIR)sort_module.f90 $(MOD_DIR)interpol_module.f90 $(MOD_DIR)ModeNorm_module.f90 buildG_dcdm1D_module_TI.f90 



all :	Mineos 

Mineos: $(OBJS) 
	$(FC) $(FFLAGS) $(OBJS) 


clean: 
	rm -f *.mod *.o

