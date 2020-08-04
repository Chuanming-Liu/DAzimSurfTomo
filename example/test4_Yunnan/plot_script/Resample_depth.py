#!/usr/bin/env python
import numpy as np 
import pdb

input               = 'Gc_Gs_model.inv'
output              = 'Gc_Gs_model_reSmp.inv'
data                = np.loadtxt(input, dtype=np.float64)

nx                  = 38 -2; ny = 42 - 2; nz = 18 - 1;
ind                 = 0
lineArr             = [] 
for zz in np.arange(0, nz):
    for yy in np.arange(0, ny, 2):
        for xx in np.arange(0, nx, 2):
            index   = xx + yy * nx + zz * nx * ny
            lineArr.append(index)
dataRe              = data[lineArr, :]
fmt                 = '%9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f'
np.savetxt(output, dataRe, delimiter=" ", fmt=fmt)