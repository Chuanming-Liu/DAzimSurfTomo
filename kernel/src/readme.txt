Makefile read me
Step 1.
Compile minos_bran_reall_kern_TI_dcdL.f --> mineos_kern_dcdL_TI 
cd ./Mineos_kern
make -f makefile_mineos_kern
   ** make sure the absolute path of mineos_kern_dcdL_TI is indicated properly in ./source_dc_dALCF/Mineos_Module/buildG_dcdm1D_module.f90,
   ** otherwise change the path
Step 2.
cd ..
make clean
make all
