## DAzimSurfTomo

DAzimSurfTomo is a package of direct inversion of surface wave for 3-D isotropic Vsv and azimuthal anisotropy without conventional tomography. Please refer to [Liu et al. (2019)](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2018JB016920) for the details of the method. The [fast marching method](http://rses.anu.edu.au/~nick/waves.html) (Rawlinson et al., 2004) is used to compute period-dependent surface wave traveltime and ray paths. The forward computation of surface wave is based on the Thomson-Haskell method (the [code of Herrman](http://www.eas.slu.edu/eqc/eqccps.html)) (Herrmann, 2013). The inversion frame is similar to [DSurfTomo](https://github.com/HongjianFang/DSurfTomo) (Fang et al. 2015) for isotropic Vs inversion with the same initial model and input data format.  

Please check the [manual](https://github.com/Chuanming-Liu/DAzimSurfTomo/blob/master/doc/Manual_DAzimSurfTomo_V2.0.md) in ./doc for usage.

V1.0: [Mineos](https://geodynamics.org/cig/software/mineos/) is used in the calculation of frequency-dependent phase velocities. Only azimuthal anisotropy is inverted. (Aug, 2017)

V2.0: Both isotropic Vsv perturbation and azimuthal anisotropy are inverted. The transfer matrix method (Herrmann, 2013) is used to calculate frequency-dependent phase velocities. (Jun, 2019)

References:

Liu, C., Yao, H., Yang, H., Shen, W., Fang, H., Hu, S., Qiao, L., 2019. Direct inversion for three-dimensional shear wavespeed azimuthal anisotropy based on surface-wave ray tracing: methodology and application to Yunnan, southwest China. Journal of Geophysics Research: Solid Earth. 124(11), 11394-11413.

Fang, H., Yao, H., Zhang, H., Huang, Y. C., & van der Hilst, R. D., 2015. Direct inversion of surface wave dispersion for three-dimensional shallow crustal structure based on ray tracing: methodology and application. Geophysical Journal International, 201(3), 1251-1263.

Rawlinson, N. & Sambridge, M., 2004. Wave front evolution in strongly heterogeneous layered media using the fast marching method, Geophys. J. Int., 156(3), 631–647.

Herrmann, R. B., 2013. Computer programs in seismology: An evolving tool for instruction and research. Seism Research Letters, 84(6),1081–1088.


