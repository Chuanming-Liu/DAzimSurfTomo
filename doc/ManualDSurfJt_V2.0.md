## User's Manual for DSurfJt (V2.0)

Chuanming Liu (chuanmingliu@foxmail.com) and Huajian Yao (hjyao@ustc.edu.cn)

### 1. Description

DSurfJt is a Rayleigh wave inversion program which can directly invert Rayleigh wave dispersion data to 3-D depth-dependent Vsv velocity and azimuthal anisotropy, which does not need the conventional intermediate step of tomography. This method was developed at the University of Science and Technology of China.  The inversion frame and isotropic inversion part are based on the [DSurfTomo](https://github.com/HongjianFang/DSurfTomo) (Fang et al. 2015) with the same initial model and input data format. The fast marching method (Rawlinson et al. 2004) is used to compute Rayleigh wave traveltime and ray paths at each period. DSurfJt includes two inversion mode for isotropic inversion and joint inversion of Vsv and azimuthal anisotropy.  Please refer to our following paper for the details.


### 2. Installation
This program has been tested successfully on Debian and MacOS platforms with gfortran. 

Please check the Makefile in ./src for installation.

### 3. Input file for inversion (para.in)

```
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c INPUT PARAMETERS
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
surfphase_forward_RV3th.dat	     	 c: traveltime data file 
17 17 4                              c: nx ny nz (grid number in lat, lon and depth direction)
26.5  101.25                         c: goxd gozd (upper left point,[lat,lon])
0.25 0.25                            c: dvxd dvzd (grid interval in lat and lon direction)
2                                    c: sablayers (2~5)
2.912 4.142                          c: minimum and maximum Vsv
1000                                 c: max(sources, receivers)
0.2                                  c: sparsity fraction
5				    				 c: maxmum of interation for joint inversion 
F                                    c: iso-mode (T: isotropic inversion; F: joint inversion)
cccccccc control parameters
950	           		                 c: weight for Vs
8.5		           	                 c: weight for Gc, Gs
0		           	                 c: damp
cccccccc periods
36                                   c: kmaxRc number of periods (followed by periods)
5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40
```
Note:

1. `nx, ny, nz` means the grid number of whole model (MOD) including the boundary grids, the grid number of the model, which updates in the inversion, is`(nx-2)*(ny-2)*(nz-1)`.

2. `goxd, gozd ` indicate the origin point (excluding the boundary points) in latitude (from north to south) and longitude (from west to east). Please make sure the inverted region full includes all the sources and receivers. 

 The inverted region: `Lat: goxd-(nx-3)*dvxd ~ goxd; Lon: gozd ~ gozd+(ny-3)*dvzd`

3. `sublayers` represents how many sublayers used to transform knot grids to layers in order to calculate the depth kernel. 
 
4. `sparsity fraction` parameter means how sparsity the sensitivity matrix is, 2-10 percent will be enough for most cases.

5. `weight` is the balancing parameter between data fitting term and smoothing regularization term. 

6. `damp`is the input parameter for LSQR, it controls the amplitude of the inverted parameter.


### 4. Data 

The format of dispersion data ,same as DSurfTomo, is as followed: 

```
# 25.148500 121.511100 1 2 0
25.158529 121.476890 0.7990
25.133539 121.499190 1.0420 
# 25.158529 121.476890 1 2 0 
25.119850 121.473190 0.6460 
# 25.128920 121.417420 1 2 0 
25.119850 121.473190 0.9430 
# 25.119850 121.473190 1 2 0 
25.090361 121.462250 0.8280 
25.083694 121.435220 1.0870 
25.133539 121.499190 1.3910
```

Lines begin with '#' represent the sources, followed by source latitude, source longitude, period index (integer), wave type and velocity type.

Each source is then followed by the receiver data: the first two columns are the latitude and longitude of the receivers, the third column is phase or group velocity (surface wave dispersion measurements). Period index (integer): index of the period vector that is listed in the parameter file para.in.

Wave type (integer): 2 for Rayleigh wave and 1 for Love wave 

Velocity type (integer): 0 for phase velocity and 1 for group velocity.

For DSurfJt fouces on the inversion of azimuthal anisotropy, the input data is limit to the Rayleigh wave phase velocity.

### 5. Initial Model (MOD)

The file name of the initial model must be 'MOD', the content looks like:

```
  0.0 10.0 35.0 60.0
  3.200  3.200  3.200  3.200  3.200  3.200  3.200  3.200  3.200  3.200  3.200  3.200  3.200  
  3.200  3.200  3.200  3.200  3.200  3.200  3.200  3.200  3.200  3.200  3.200  3.200  3.200  
  3.200  3.200  3.200  3.200  3.200  3.200  3.200  3.200  3.200  3.200  3.200  3.200  3.200  
  3.200  3.200  3.200  3.200  3.200  3.200  3.200  3.200  3.200  3.200  3.200  3.200  3.200  
  3.200  3.200  3.200  3.200  3.200  3.200  3.200  3.200  3.200  3.200  3.200  3.200  3.200  
  3.200  3.200  3.200  3.200  3.200  3.200  3.200  3.200  3.200  3.200  3.200  3.200  3.200  
```

The first line is the depths (km) of grid points in the vertical direction.

Then followed by shear velocity values. The order is altitude first, then longitude, followed by depth. Each row represents shear velocity values at different latitude at a single longitude and a certain depth, then followed by next longitude, then depth.

In case of a 3D initial velocity model, note we have boundary values included in this file. The original point will be in the upper northwest corner.


### 6. Output files

#### isotropic mode

1.`DSurfTomo.dat` (grid model)

```
col 1: longitude (degree)
col 2: latitude (degree)
col 3: depth (km)
col 4: Vsv (km/s)
```

2.`MOD_Ref`: same as `MOD`, which can be used as initial model for joint inversion.

#### joint inversion mode (or anisotropic inversion mode) 

1.`Gc_Gs_model.inv` (layered model)
 
```
col 1: longitude (degree)
col 2: latitude (degree)
col 3: depth (lower interface) (km)
col 4: Vsv (layer averaged Vsv) (km/s)
col 5: fast direction (from north) (degree)
col 6: amplitude 
col 7: Gc/L (%)
col 8: Gs/L (%)
```
2.`period_Azm_tomo.inv` (layered model)
 
```
col 1: longitude (degree)
col 2: latitude (degree)
col 3: period (s)
col 4: Rayleigh wave phase velocity (km/s)
col 5: fast direction (from north) (degree)
col 6: relative amplitude (/a0_c)
col 7: amplitude (sqrt(a1**2+a2**2))
col 8: a1_cos
col 9: a2_sin
```

### References

Liu, C., Yao, H., Yang, H. ....


Fang, H., Yao, H., Zhang, H., Huang, Y. C., & van der Hilst, R. D., 2015. Direct inversion of surface wave dispersion for three-dimensional shallow crustal structure based on ray tracing: methodology and application. Geophysical Journal International, 201(3), 1251-1263.

Rawlinson, N. & Sambridge, M., 2004. Wave front evolution in strongly heterogeneous layered media using the fast marching method, Geophys. J. Int., 156(3), 631–647


  