#!/bin/bash

#SBATCH -J AmzKernel
#SBATCH -o AmzKernel_%j.out
#SBATCH -e AmzKernel_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --partition=shas
#SBATCH --qos=normal
#SBATCH --time=02:30:00
#SBATCH --mem=MaxMemPerNode
module load intel
# module load fftw
#. ~/.my.bashforSEED2COR

path=`pwd`
cppexe=$path/DirectKernel.sh
$cppexe
