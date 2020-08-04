#!/usr/bin/env python
import numpy as np
import pdb

input                   = 'period_Azm_tomo.inv'
output                  = 'Period_Azm_reSmp.inv'

data                    = np.loadtxt(input, dtype=np.float64)

nx                      = 38 -2
ny                      = 42 - 2 
nt                      = 36

# add percentage column
newdata                 = np.zeros((data.shape[0], data.shape[1]+1));
newdata[:, 0:4]         = data[:, 0:4]
newdata[:, 5:]          = data[:, 4:]

num                     = nx * ny

for zz in np.arange(nt):
    rows                = np.arange(num*zz, num*(zz+1))
    meanp               = np.mean(data[rows, 3])
    newdata[rows, 3]    = (data[rows, 3]-meanp)/meanp*100
    newdata[rows, 4]    = meanp   


lineArr                 = [] 
for zz in np.arange(0, nt):
    for yy in np.arange(0, ny, 2):
        for xx in np.arange(0, nx, 2):
            index       = xx + yy * nx + zz * nx * ny
            lineArr.append(index)

dataRe                  = newdata[lineArr, :]
fmt                     = '%9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f'
np.savetxt(output, dataRe, delimiter=" ", fmt=fmt)
