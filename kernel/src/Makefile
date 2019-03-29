#
#	Make file for program Main_DirectInv_Aniso_kernel
#

# Ways to use this make file 
#    > `make all' 
#    > `make clean'
#    

#---- start user specific make  (Edit below) ----------------------------------
USER_DIR  = .
USER_DIR_F = $(USER_DIR)/DK_src
USER_DIR_MOD=$(USER_DIR)/Mineos_Module
USER_DIR_RD=$(USER_DIR)/Reader

FC	= gfortran
FFLAGS_DK = -ffixed-line-length-132 -O3 -c -I../Mineos_Module  -mcmodel=medium -m64 -W -fbounds-check -ffloat-store
FFLAGS = -ffixed-line-length-132 -O3  -mcmodel=medium -m64 -W -fbounds-check -ffloat-store




USER_OBJ_M =$(USER_DIR_MOD)/buildG_dcdm1D_module_TI.o\
	      $(USER_DIR_MOD)/cbasic_module.o\
            $(USER_DIR_MOD)/interpol_module.o\
            $(USER_DIR_MOD)/ModeNorm_module.o\
            $(USER_DIR_MOD)/nrtype.o\
            $(USER_DIR_MOD)/nrutil.o\
            $(USER_DIR_MOD)/sort_module.o

USER_DK_F = CalSenKernel.f90 Main_DirectInv_Aniso_kernel.f90 mineos_model_output.f90
USER_DK_OBJ = $(USER_DIR_F)/CalSenKernel.o $(USER_DIR_F)/Main_DirectInv_Aniso_kernel.o $(USER_DIR_F)/mineos_model_output.o
USER_MAIN = DKernel

USER_RD_F =  $(USER_DIR_RD)/Read_DKernel.f90
USER_READ = DKReader
all : user

user :
	cd $(USER_DIR_MOD); make clean; make
	cd ../
	cd $(USER_DIR_F); rm -f *.o;$(FC) $(FFLAGS_DK) $(USER_DK_F)
	cd ../;	rm -f *.o *.mod
	$(FC) $(FFLAGS) -o $(USER_MAIN)  $(USER_DK_OBJ) $(USER_OBJ_M)
	$(FC) $(FFLAGS) -o $(USER_READ)  $(USER_RD_F)
	mv $(USER_MAIN) $(USER_READ)  ../bin
	echo "Make DK successfully!"
clean: 
	rm -f *.o *.mod
	rm -f $(USER_DIR_F)/*.o
	rm -f $(USER_DIR_MOD)/*.o *.mod
	rm -f ../bin/$(USER_MAIN)
	rm -f ../bin/$(USER_READ)


