Modification history:

CalSurfG_azimuth.f90
About modefication
1. the initial code: CalSurfG.f90 is from SurfTomo15_new/srcsmooth supplied by
Fang in about May/2016
2. The CalSurfG.f90 seems edited in model add interval which differs from the
original code with readme file.


Change 1.
about subroutine raypaths
note danger
1.  subroutine azdist varibles precision is double precision, rpath is
real(kind=i10)
2. for a ray path, calculate to which point: for original cycle fdm sum to the
front point of nrp (the source ponit),and so is the azimuthal code
but I still do not why one cycle to calculate two point and end with nrp-1
point which  means j+2 is source ponit.
两点可能会产生错误，1. fdm 最后的计算公式j 和j+1 ,我直接在后面乘以cos2pi
sin2pis
2. 和源程序一样，每次算两个点，算到j=nrp-3, 抛去了最后的震源点

3. 第一次程序有问题,先把所有raypath方位角都算了,其实不用,所以p1.p2不是问题
   azdist()input and output all in degree, so, wrong
   Fang code all latitude is in colatitude; subroutine adizst use latitude in degree; cos use in rad
   New modification: 1. add fmcd=0; 2. cal rgpsi 3.fdmc,fdms.
   Promble: should /Ck^2 but use Cj^2 and Cp^2 Fang means they seem in close value.

change 2. grid model-> layer model
   nz grid inversion point -> nz-1 inversion lay
   for last point nz, which sensitivity kernel stands for half space, only use nz-1 in the isotropic model-----means still use nz-1 inversion grid in depth. for mineos doesn't contain half space, ref model. 

change 3. one inversion---> iteration inversion ----> one inversion 
for now, just keep the iteration inversion 

Change 4. April 4. 2017--Test6
  GaussianLS.f90
 twoXYpow2=4.*(LcorrXY**2)--->twoXYpow2=2.*(LcorrXY**2)--Line 129,71
 twoZZpow2=4.*(LcorrZZ**2)---> twoZZpow2=2.*(LcorrZZ**2)--Line 153,93
   MainLS.f90 about the Norm2(res)-
   add line 552 ---resSigma(i)=resbst(i)*1/sigmaT(i): resNorm will get small
Change 5. April. 7 2017--Test 6
Add different sigma Gsc
        sigmaGcs(1)=sigmaGc
        sigmaGcs(2)=sigmaGs


Change 6. April. 13 (Test7)
Add different sigma Gc,Gs at different layer.
same layer sigmaGc=sigmaGc
different layer same grid sigmaGc=(sigmaGc(k1)+sigmaGs(k2))/2

Change 7. fix bug.
MainLS.f90
Line569:  resbst(i)=obst(i)-fwdT(i)

!-----------------------------source_inversion_test8 (test8)----------------------------!
Change 8. add damping term  
misfit function: ||W(Gm-d)||2^2+alpha*||Lm||2^2+damp*||Em||2^2

Change 9.  1. fix bug of isotropci dispersion phase velocity------> 
 	     2. fix bug of CalResidualNorm--fix the wrong ResNorm2 -> CalResidualNorm
	     3. fix bug of result traveltime 
!-----------------------------source_inversion_iter (test8)----------------------------!
Change 10. Apirl, 24.
           1. change the one step inversion of dVs/Vs, Gc/L, Gs/L into the iteration inversion staggered inversion.
Change 11. Apirl,27 fix bug.
           subroutine: DampingE:
           col(nar+1)=para2+maxvp Line 35,43
           But: Iteration, problem at the step isotropic inversion.
!-----------------------------source_inversion_iter_diff (test9---test7)----------------------------!
Change 12. May. 2017
          1. separate the inversion of dVs from Gc/L Gs/L, still staggered inversion.---
	    2. add TikhonovRegul.f90----for the guassian regularization works bad.

!-----------------------------source_inversion_iter_diff (Big debug)----------------------------!
Change 12. July 5. 2017
  rpathsAzim.f90
                  rdc1=vi(m)*wi(l)/vel**2
                  rdc2=vio(m)*wio(l)/velo**2
                  rd1= -(rdc1+rdc2)*dinc/2.0
                  rd2=fdm(ivzt-2+l,ivxt-2+m)
                  fdm(ivzt-2+l,ivxt-2+m)=rd1+rd2
! anisotropy part cos--Liu
                  rd1= -(rdc1*cos(2.0*rgpsi)+rdc2*cos(2.0*rgpsi))*dinc/2.0
                  rd2=fdm(ivzt-2+l,ivxt-2+m)
                  fdmc(ivzt-2+l,ivxt-2+m)=rd1+rd2
! anisotropy part sin
                  rd1= -(rdc1*sin(2.0*rgpsi)+rdc2*sin(2.0*rgpsi))*dinc/2.0
                  rd2=fdm(ivzt-2+l,ivxt-2+m)
                  fdms(ivzt-2+l,ivxt-2+m)=rd1+rd2

2019-05-30
use Tregn.f90 from cps to calculate the sensivity kernel
potenial problem: cannot use -fopenmp for depthkernelTI.f90 
but will not affect the time cost

2019-06-03
improve iso-inversion part. Unknown small difference in 6-th iso-inversion, which may be related to gfortran version (works on gcc version 4.8.5 summit/marine, not on gcc version 6.3.0)

2019-06-04 
Makefile: %.o: %.f90
	$(FC) $(FFLAGS) -c $(@F:.o=.f90)  -o $@
will cause LSQR crashed in Ani-inversion.
with or without, will make almost same result for Ios-inversion.

2019-06-04
Modify the inversion into joint inversion for isotropic part and anisotropic part simultaneously.




