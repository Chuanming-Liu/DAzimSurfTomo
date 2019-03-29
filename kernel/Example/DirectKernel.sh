#!/bin/bash

exe=../bin
DK=$exe/DKernel
DKReader=$exe/DKReader

# path
path=`pwd`
#
# infold=$path/infile
infold=$path/infile

# output=$path/output
output=$path/outputRV

kernel_fold=$output/kernel
Gramdd_fold=$output/Gramdd
mode_fold=$output/Sph_mode
mdl_fold=$output/pro_mdl

Lon_fold=$path/LonKernel

if [ ! -d $Lon_fold ]; then
   mkdir $Lon_fold
fi

if [ ! -d $output ]; then
   mkdir $output
fi

if [ ! -d $kernel_fold ]; then
   mkdir $kernel_fold
fi

if [ ! -d $Gramdd_fold ]; then
   mkdir $Gramdd_fold
fi

if [ ! -d $out_fold ]; then
   mkdir $out_fold
fi

if [ ! -d $mode_fold ]; then
  mkdir $mode_fold
 fi

if [ ! -d $mdl_fold ]; then
   mkdir $mdl_fold
fi
# 1.-------------------------------------------------
echo "_______________RUN DK LON CYCLE_____________"
#MOD			                        c: input iso-Vs model file
#2 2 9                            c: nx ny nz (grid number in lat lon and depth direction)
#1	                              c: ny_clc 1:ny
#3                                c: 'minthk' = layerthick / sublayernum, can be 1, 2 ,3
#2                                c: kmaxRc (followed by periods)
#40 50
#PREMQL6ic_12742u_intp_2km.card   c: reference Mineos format model
nx=16
ny=16
nz=4
minthk=2
kmaxRc=36
mdl=$infold/MOD
refmdl=$infold/PREMQL6ic_12742u_intp_2km.card

#ny=1
for ((i=1;i<=$ny;i=i+1));do
echo "ny_cycle=" $i
(echo "$mdl
$nx $ny $nz
$i
$minthk
$kmaxRc
5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40
$refmdl" | nohup $DK >log_DK_$i 2>&1 ) &

done
# 2.---------------------------------------------
# zip
sw=0
#last=Lon_*"$ny"_SenG_Vs_Gsc.dat

while [[ $sw -ne $ny ]];do
 sw=`find ./ -name "Lon*.dat" | wc -l `
 echo $sw
 sleep 1m
done

mv *.Gdused $Gramdd_fold
mv *.fre $mode_fold
mv *.card $mdl_fold
mv *dL_dA* $kernel_fold
rm *.split
Lonfile=Lonlist.txt
ls Lon*.dat > $Lonfile

# 3.-----------------------------------------
# Read_Dkernel
# input list
#2 2 9                         c: nx ny nz (grid number in lat lon and depth direction)
#2                             c: kmaxRc (followed by periods)
#Lonlist.txt                   c: Lonfile
#Sen_dcdL_dcdA.dat             c: output kernel file

(echo "$nx $ny $nz
$kmaxRc
$Lonfile" | nohup $DKReader >log_Read 2>&1 ) &

sw=0
while [[ $sw -ne 1 ]];do
 sw=`find ./ -name "Sen_dcdL_dA_dC_dF.dat" | wc -l `
 sleep 10m
done

mv *Lon*.dat log_DK* $Lon_fold


echo "DK finished!"
