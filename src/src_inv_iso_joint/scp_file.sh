#!/bin/bash
inpath=`pwd`
path=/home/liu/programs/DSurf_package
sshpass -p lccmm1234 scp -r $inpath liu@marin.colorado.edu:$path
