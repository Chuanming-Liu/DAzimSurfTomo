c----------------------------------------------------------------------c
c                                                                      c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                 c
c      VOLUME III                                                      c
c                                                                      c
c      PROGRAM: SLEGN96                                                c
c                                                                      c
c      COPYRIGHT 2010                                                  c
c      R. B. Herrmann                                                  c
c      Department of Earth and Atmospheric Sciences                    c
c      Saint Louis University                                          c
c      221 North Grand Boulevard                                       c
c      St. Louis, Missouri 63103                                       c
c      U. S. A.                                                        c
c                                                                      c
c----------------------------------------------------------------------c
c
c
c     This program calculates the group velocity and partial
c     derivatives of Rayleigh waves for any plane multi-layered TI
c     model.  The propagator-matrix, instead of numerical-
c     integration method is used.
c     than Harkrider formalisms are concerned.
c
c     R. B. Herrmann, St. Louis
c     University, Oct. 10, 1981.  Modified for use in surface
c     wave inversion, with addition of spherical earth flattening
c     transformation and numerical calculation of group velocity
c     partial derivatives by David R. Russell, St. Louis
c     University, Jan. 1984.
c
c     Rewrite of theory to agree more with hspec96 wavenumber
c     integration code
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c Revision history:
c       07 AUG 2002 - make string lengths 120 characters from 80
c       13 OCT 2006 - verbose output of energy integrals
c       26 SEP 2008 - fixed undefined LOT in subroutine up
c       14 JUN 2009 - reformulate in general terms as per book
c       01 AUG 2010 - change code to agree with book  for
c                     migration to TI.
c                     explicit none
c       20 MAR 2012 - corrected compiler warnings
c       05 JUN 2019 - modified by Lili Feng and Chuanming Liu
c                     to calculate Rayleigh wave phase velocity
c                     sensitivity kernels.
c       10 JUN 2019 - add subroutine: sphere_tdisp96
c                     assimilate other subroutines
c                     Chuanming Liu.
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        subroutine tregn96(nl_in, d_in, TA_in, TC_in, TF_in, TL_in,
     1      TN_in, TRho_in, qai_in, qbi_in, etapi_in, etasi_in,
     1      frefpi_in, frefsi_in,  Nt_in, t_in, cp_in,
     1      dcdah_out, dcdbv_out, dcdn_out)
C     1      d_out, TA_out, TC_out, TF_out, TL_out, TN_out, TRho_out,
C     1      u_out, ur_out, tur_out, uz_out, tuz_out,
C     1      dcdh_out,dcdav_out,dcdah_out,
C     1      dcdbv_out,dcdbh_out,dcdn_out,dcdr_out)

        implicit none
c-----
c       common blocks
c-----
        integer NL
        parameter(NL=200)
c-----
c       LIN - unit for FORTRAN read from terminal
c       LOT - unit for FORTRAN write to terminal
c       LER - unit for FORTRAN error output to terminal
c       NL  - number of layers in model
c-----
        common/timodel/d(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qai(NL),qbi(NL),etapi(NL),etasi(NL),
     3      frefpi(NL), frefsi(NL)
        real d,TA,TC,TN,TL,TF,TRho,qai,qbi,etapi,etasi,frefpi,frefsi

        common/depref/refdep
        real refdep

        common/timod/  zd(NL),zta(NL),ztc(NL),ztf(NL),
     1      ztl(NL),ztn(NL),zrho(NL),zqai(NL),zqbi(NL),
     2      zetap(NL),zetas(NL),zfrefp(NL),zfrefs(NL)
        real*8 zd, zta, ztc, ztf, ztl, ztn, zrho,
     1      zqai, zqbi, zetap, zetas,
     1      zfrefp, zfrefs
        common/pari/ mmax
        integer mmax

        common/eigfun/ ur(NL),uz(NL),tz(NL),tr(NL),uu0(4),
     *    dcdah(NL),dcdav(NL),dcdbh(NL),dcdbv(NL),
     2    dcdn(NL),dcdr(NL),dcdh(NL)
        real*8 ur, uz, tz, tr, uu0,
     1     dcdah, dcdav, dcdbh, dcdbv, dcdn, dcdr, dcdh

        real*4 sdcdah(NL), sdcdav(NL), sdcdbh(NL), sdcdbv(NL),
     1     sdcdn(NL), sdcdh(NL), sdcdr(NL)
        real*4 spur(NL), sptr(NL), spuz(NL), sptz(NL)

        common/wateri/iwat(NL)
        integer iwat
c-----
c       model characterization
c-----
        common/modspec/allfluid
        logical allfluid

        common/sumi/   sumi0,sumi1,sumi2,sumi3,flagr,are,ugr
        real*8 sumi0, sumi1, sumi2, sumi3, flagr, are, ugr

        common/sphereL/vtp(NL),dtp(NL),rtp(NL)
        real*4 vtp,dtp,rtp

c-----
c       function prototype
c-----
        integer lgstr

c-----
c       local variables
c-----
        integer LER, LIN, LOT
        parameter(LER=0,LIN=5,LOT=6)
        integer ipar(20)
        integer nipar(20)
        integer i,j,k,l
        integer iunit,iiso,iflsph,idimen,icnvel,ierr
        integer MAXMOD
        parameter(MAXMOD=2000)
        integer lss, lsso, lrr, lrro, n1, n2, npts, ifirst
        integer ifunc, mode, nsph, nper, nmodes, mmaxot

        logical dotmp
        logical ext
        logical dogam
        logical dderiv
        logical verbose
        logical nwlyrs, nwlyrr
        logical dolove, dorayl

        character mname*120
        character hsfile*120, hrfile*120, title*120
        character*12 fname(2)

        real*4 fpar(10)
        real*4 s1, dephs, dephr, faclov, facray, hr, hs, dt
        real*4 deplw, depthr, depths, depup
        real*4 ohs, ohr
        real*4 rare, rtr, rtz, rur0, rur, ruz
        real*4 sare, sd2ur, sd2uz, sduz, sur, sur0, suz, sdur
        real*4 sum, sumgr, sumgv, sumkr
        real*4 wvnsrc, wvnrec
        real*4 twopi

        real*8 durdz, duzdz, d2urdz, d2uzdz
        real*8 rurdz, ruzdz
        real*8 c, omega, wvno, gammar, csph, usph
        real*8 cp(MAXMOD), t

c   Input arrays, added by LF & LCM
        integer NP
        parameter(NP=60)
        real*4 hs_in, hr_in, ohr_in, ohs_in, refdep_in
        integer nl_in, iflsph_in
        real*4 d_in(NL),TA_in(NL),TC_in(NL),TF_in(NL)
        real*4 TL_in(NL),TN_in(NL),TRho_in(NL)
        real*4 qai_in(NL),qbi_in(NL),etapi_in(NL)
        real*4 etasi_in(NL),frefpi_in(NL),frefsi_in(NL)
        integer Nt_in, idisp
        integer*4 mode_in(NP)
        real*4 t_in(NP),cp_in(NP)

        logical dogam_in
        ! Output arrays, added by LF & LCM
        ! in LF's version, output kernel in real*8,
        ! But output is from sdcdah.. which is in real*4
        real*8 u_out(NP)
        real*4 ur_out(NP,NL), tur_out(NP,NL), uz_out(NP,NL)
        real*4 tuz_out(NP,NL)
        real*4 dcdh_out(NP,NL),dcdav_out(NP,NL)
        real*4 dcdah_out(NP,NL),dcdn_out(NP,NL)
        real*4 dcdbv_out(NP,NL),dcdbh_out(NP,NL)
        real*4 dcdr_out(NP,NL)

        real*4 d_out(NL),TA_out(NL),TC_out(NL),TF_out(NL)
        real*4 TL_out(NL),TN_out(NL),TRho_out(NL)
c-----
c       machine dependent initialization
c-----
c        call mchdep() by LCM.
C        implicit none
c-----
c       parse command line information
c-----
        ! by LCM initialization
        hs_in = 0.0
        hr_in = 0.0
        ohr_in = 0.0
        ohs_in = 0.0
        refdep_in = 0.0
        ! turn on attenuation
        dogam_in = .true.
        ! 0 flat, 1: spherical
        iflsph_in = 1

cc        call gcmdln(hsfile,hrfile,hs,hr,dotmp,dogam,dderiv,
cc     1      nipar,verbose)
        verbose = .false.
        dogam   = dogam_in
        dderiv  = .true.
        nipar(4) = 1
        nipar(5) = 1
        nipar(6) = 1
        nipar(7) = 1
        nipar(8) = 1
        nipar(9) = 1
        nipar(10) = 1
        nipar(11) = 1

c-----
c       wired in file names     - output of this program
c                           is always in the
c               binary file tregn96.egn or tregn96.der
c                               - input from tdisp96 is always
c               binary file tdisp96.ray, DELETED!
c-----

c-----
c       get control parameters from tdisp96.dat, DELETED !
c-----
        ohs = ohs_in
        ohr = ohr_in
        refdep  = refdep_in
c-----
c       get the earth model, DELETED !
c-----
        twopi=6.283185307179586d+00
        iflsph = iflsph_in
        mmax = nl_in
        idisp=1
C   assign model arrays, by LF
        do 339 i=1,mmax
            d(i) = d_in(i)
            TA(i) = TA_in(i)
            TC(i) = TC_in(i)
            TF(i) = TF_in(i)
            TN(i) = TN_in(i)
            TL(i) = TL_in(i)
            Trho(i) = TRho_in(i)
            qai(i) = qai_in(i)
            qbi(i) = qbi_in(i)
            etapi(i) = etapi_in(i)
            etasi(i) = etasi_in(i)
            frefpi(i) = frefpi_in(i)
            frefsi(i) = frefsi_in(i)
            if(TN(i).le.1.0e-4*TA(i))then
                iwat(i) = 1
            else
                iwat(i) = 0
            endif
 339   continue
c   - End of get model, by LF

C        WRITE(6,*)"tregn96 check1"
C        do i=1,mmax
C         WRITE(6,*)i,d(i),TA(i),TC(i),TL(i),TN(i),TF(i),
C     1      TRho(i)
C        enddo
C        write(6,*)"refdep:", refdep


c-----
c       check for water only
c-----
        allfluid = .true.
        do 1200 i=1,mmax
            if(tn(i) .gt. 0.0001*ta(i))then
                allfluid = .false.
            endif
 1200   continue
c-----
c       get the Q information into the program
c       since this is not carried through in the tdisp96 output
c-----
        do 1234 i=1,mmax
            if(.not.dogam)then
                qai(i) = 0.0
                qbi(i) = 0.0
            endif
            if(qai(i).gt.1.0)qai(i)=1.0/qai(i)
            if(qbi(i).gt.1.0)qbi(i)=1.0/qbi(i)
            zqai(i) = qai(i)
            zqbi(i) = qbi(i)
            if(frefpi(i).le.0.0)frefpi(i) = 1.0
            if(frefsi(i).le.0.0)frefsi(i) = 1.0
            zfrefp(i) = frefpi(i)
            zfrefs(i) = frefsi(i)
            zetap(i) = etapi(i)
            zetas(i) = etasi(i)
 1234   continue
        nsph = iflsph
c-----
c       get source depth, and sphericity control
c-----
        if(hs.lt.-1.0E+20)then
            depths = ohs
        else
            depths = hs
        endif
        if(hr.lt.-1.0E+20)then
            depthr = ohr
        else
            depthr = hr
        endif
        depthr = depthr + refdep
        depths = depths + refdep
c-----
c       if there is a spherical model, map the source
c            and receiver depth
c       in the spherical model to the equivalent depth in the flat model
c-----
        if(nsph.gt.0)then
                dephs = 6371.0 *alog( 6371.0/(6371.0 - (depths-refdep)))
                dephr = 6371.0 *alog( 6371.0/(6371.0 - (depthr-refdep)))
        else
                dephs = depths
                dephr = depthr
        endif
c----
c       see if the file tdisp96.ray exists, DELETED
c-----

c-----
c       open output file of tdisp96 and of tregn96, DELETED !
c-----

c-----
c       obtain the earth model: note if the original was spherical,
c       this will be the transformed model, DELETED !
c-----
c----   Added by LCM. gtsmdt deleted, which get flat model from tdisp96.ray
c----   To perform spherical model to flat model transformation
c----   Use the sphere subroutine.
        if (nsph .gt. 0)then
            call sphere_tdisp96()
        endif
        do 340 i=1,mmax
            d_out(i)=d(i)
            TA_out(i)=TA(i)
            TC_out(i)=TC(i)
            TF_out(i)=TF(i)
            TN_out(i)=TN(i)
            TL_out(i)=TL(i)
            TRho_out(i)=Trho(i)
  340   continue
c-----
c       define modulus of rigidity, also get current transformed model
c       parameters
c       set fpar(1) = refdep
c       set ipar(1)   = 1 if medium is spherical
c       set ipar(2)   = 1 if source is in fluid
c       set ipar(3)   = 1 if receiver is in fluid
c       set ipar(4)   = 1 if eigenfunctions are output with -DER flag
c       set ipar(5)   = 1 if dc/dh are output with -DER flag
c       set ipar(6)   = 1 if dc/dav are output with -DER flag
c       set ipar(7)   = 1 if dc/dbv are output with -DER flag
c       set ipar(8)   = 1 if dc/dr are output with -DER flag
c       set ipar(9)   = 1 if dc/dah are output with -DER flag
c       set ipar(10)  = 1 if dc/dn are output with -DER flag
c       set ipar(11)  = 1 if dc/dbh are output with -DER flag
c-----
        deplw = 0.0
        depup = 0.0
        ipar(2) = 0
        ipar(3) = 0
        do 185 i=1,mmax
            zd(i)   = d(i)
            zta(i)  = ta(i)
            ztc(i)  = tc(i)
            ztl(i)  = tl(i)
            ztn(i)  = tn(i)
            ztf(i)  = tf(i)
            zrho(i) = trho(i)
            depup = deplw + d(i)
            if(tn(i) .lt. 0.0001*ta(i))then
                iwat(i) = 1
                if(depths.ge.deplw .and. depths.lt.depup)then
                    ipar(2) = 1
                endif
                if(depthr.ge.deplw .and. depthr.lt.depup)then
                    ipar(3) = 1
                endif
            else
                iwat(i) = 0
            endif
            deplw = depup
c           WRITE(6,*)i,zd(i),zta(i),ztc(i),ztl(i),ztn(i),ztf(i),
c     1          zrho(i)
  185   continue
        do 186 i=4,11
            ipar(i) = nipar(i)
  186   continue
c    putmdt, DELETED !

c-----
c               If a spherical model is used, reconstruct the thickness
c               of the original model to get the mapping from
c               spherical to flat. We need this for correct
c               spherical partial derivatives.
c
c               THIS MUST BE DONE BEFORE SOURCE LAYER INSERTION
c-----
        if(nsph.gt.0)then
                call bldsph()
        endif
c-----
c       split a layer at the source depth
c       the integer ls gives the interface at which the eigenfunctions
c       are to be defined
c-----
        call insert(dephs,nwlyrs,lsso)
        call insert(dephr,nwlyrr,lrro)
        call srclyr(dephs,lss)
        call srclyr(dephr,lrr)
        !-----debug
        !       output the new model with source layer
        !-----
        !       WRITE(6,*)"tregn96 check1"
        !       write(6,*)'lss,lrr:',lss,lrr
        !        do  i=1,mmax
        !            WrItE(6,*)i,iwat(i),zd(i),zta(i),ztc(i),
        !     1    ztf(i),ztl(i),ztn(i),zrho(i)
        !        enddo
        !-----
        !       end debug
        !-----
        twopi=2.*3.141592654
        ifirst = 0

        ! Main cycle!
  400   continue
        !       gtshed, COMMENTED !
        !        call gtshed(1,ifunc,mode,t,ierr)
        !       ifunc, mode, t, ierr

        !       Currently only support fundamental mode!
        !       mode_in(idisp)
        ifunc = 2
        mode  = 1
C        t     = t_in(idisp)
        if (idisp.gt.Nt_in)then
            ierr=200
            go to 700
        else
            ierr=0
        endif
        t     = t_in(idisp) ! cp from 441
        idisp = idisp+1
        if(ierr.ne.0)go to 700
        s1=t
c-----
c       DEBUG OUTPUT
c-----
c        write(6,*) ifunc,mode,s1
c-----
c       END DEBUG OUTPUT
c-----
c       puthed, COMMENTED !
c        call puthed(2,ifunc,mode,s1)

        omega=twopi/t
        if(ifunc.lt.0) go to 700
        if(mode.le.0) go to 400
c COMMENTED
c        read(1) (cp(k),k=1,mode)

        do 600 k=1,mode
                c=cp_in(idisp-1)
c-----
c       main part.
c-----
            wvno=omega/c
            call svfunc(omega,wvno)
C            WRITE(6,*)'mode=',k,'wvno,omega,c=',wvno,omega,c
            call energy(omega,wvno,mmax)
c-----
c       the gamma routine will use the spherical model, but the
c       frequency dependence and Q of the original model
c-----
C        WRITE(6,*)'c(noq)=',c
            if(dogam)then
                call gammap(omega,wvno,gammar)
                c = omega/wvno
C        WRITE(6,*)'c(q)=',c
            else
                gammar = 0.0d+00
            endif
c------
c     also check for possible conversion errors in IEEE
c     conversion from double precision to single precision
c-----
            if(dabs(uu0(1)).lt.1.0d-36)uu0(1)=0.0d+00
            if(dabs(uu0(2)).lt.1.0d-36)uu0(2)=0.0d+00
            if(dabs(uu0(3)).lt.1.0d-36)uu0(3)=0.0d+00
            if(dabs(uu0(4)).lt.1.0d-36)uu0(4)=0.0d+00
            if(dabs(c).lt.1.0d-36)c=0.0d+00
            if(dabs(ugr).lt.1.0d-36)ugr=0.0d+00
            if(dabs(are).lt.1.0d-36)are=0.0d+00
            if(dabs(sumi0).lt.1.0d-36)sumi0=0.0d+00
            if(dabs(sumi1).lt.1.0d-36)sumi1=0.0d+00
            if(dabs(sumi2).lt.1.0d-36)sumi2=0.0d+00
            if(dabs(flagr).lt.1.0d-36)flagr=0.0d+00
            if(dabs(gammar).lt.1.0d-36)gammar=0.0d+00

c-----
c       output necessary eigenfunction values for
c       source excitation
c-----
            mmaxot = mmax
            if(nsph.gt.0)then
c-----
c           sphericity correction for partial derivatives
c           of the original model
c-----
                call sprayl(omega,c,mmaxot,csph,usph,ugr)
                wvno = omega / csph
                u_out(idisp-1) = usph
            else
                u_out(idisp-1) = ugr
            endif
c            WRITE(6,*)'c(flat)=',c,' u(flat)=',ugr, idisp, Nt_in
            u_out(idisp-1) = ugr

C STOP HERE
c-----
c      get the derivatives of the eigenfunctions required
c      for source excitation from the definition of stress. For
c      completeness get the second derivative from the first
c      derivatives and the equation of motion for the medium
c-----
            duzdz = ( tz(lss) + wvno*ztf(lss)*ur(lss))/ztc(lss)
            if(iwat(lss).eq.1)then
                durdz  = wvno*uz(lss)
                d2urdz = wvno*duzdz
            else
                durdz = -wvno*uz(lss) + tr(lss)/ztl(lss)
                d2urdz = (-wvno*duzdz*(ztf(lss)+ztl(lss))
     1              - ur(lss)*(zrho(lss)*omega*omega -
     2              wvno*wvno*zta(lss)))/ztl(lss)
            endif

            d2uzdz = ( - uz(lss) *( zrho(lss)*omega*omega -
     1          ztl(lss)*wvno*wvno) +
     2          wvno*durdz*(ztl(lss)+ztf(lss)))/ ztc(lss)

            if(dabs(duzdz).lt.1.0d-36)duzdz=0.0d+00
            if(dabs(durdz).lt.1.0d-36)durdz=0.0d+00
            if(dabs(d2uzdz).lt.1.0d-36)d2uzdz=0.0d+00
            if(dabs(d2urdz).lt.1.0d-36)d2urdz=0.0d+00

            ruzdz = ( tz(lrr) + wvno*ztf(lrr)*ur(lss))/ztc(lrr)
            if(iwat(lrr).eq.1)then
                rurdz = wvno*uz(lss)
            else
                rurdz = -wvno*uz(lss) + tr(lss)/ztl(lrr)
            endif
            if(dabs(ruzdz).lt.1.0d-36)ruzdz=0.0d+00
            if(dabs(rurdz).lt.1.0d-36)rurdz=0.0d+00

            sur = sngl(ur(lss))
            sdur = sngl(durdz)
            sd2ur = sngl(d2urdz)
            suz = sngl(uz(lss))
            sduz = sngl(duzdz)
            sd2uz = sngl(d2uzdz)
            sare = sngl(are)
            wvnsrc = sngl(wvno)
            sur0 = sngl(uu0(1))

            if( verbose ) then
                 WRITE(LOT,2)t,c,ugr,gammar,
     1                sumi0,sumi1,sumi2,sumi3,
     2                flagr,are,sur0
    2    format(' T=',e15.7,'  C=',e15.7,'  U=',e15.7,'  G=',e15.7/
     1          'I0=',e15.7,' I1=',e15.7,' I2=',e15.7,' I3=',e15.7/
     2          ' L=',e15.7,' AR=',e15.7,'  E=',f15.7)
            endif

            rur = sngl(ur(lrr))
            rtr = sngl(tr(lrr))
            ruz = sngl(uz(lrr))
            rtz = sngl(tz(lrr))
            rare = sngl(are)
            wvnrec = sngl(wvno)
            rur0 = sngl(uu0(1))

            sumkr = 0.0
            sumgr = 0.0
            sumgv = 0.0
            if(nsph.gt.0)then
                wvno = omega/ csph
                ugr = usph
            endif
        if(dderiv)then
c----C
c               if a layer was inserted, get the partial
c               derivative for the original layer.
c           The sequence is cannot be changed.
c           Originally the model grows because the source layer
c           is added and then the receiver layer
c           This we must first strip the receiver and lastly the
c           source
c----
c       initialize
c-----
C            WRITE(6,*)'nwlyrr,lrro:',nwlyrr,lrro
C            WRITE(6,*)'nwlyrs,lsso:',nwlyrs,lsso
            if(nwlyrr)then
                call collap(lrro+1,mmaxot)
            endif
            if(nwlyrs)then
                call collap(lsso+1,mmaxot)
            endif
            call chksiz(ur,spur,mmaxot)
            call chksiz(tr,sptr,mmaxot)
            call chksiz(uz,spuz,mmaxot)
            call chksiz(tz,sptz,mmaxot)
            call chksiz(dcdh,sdcdh,mmaxot)
            call chksiz(dcdah,sdcdah,mmaxot)
            call chksiz(dcdav,sdcdav,mmaxot)
            call chksiz(dcdbh,sdcdbh,mmaxot)
            call chksiz(dcdbv,sdcdbv,mmaxot)
            call chksiz(dcdn ,sdcdn ,mmaxot)
            call chksiz(dcdr,sdcdr,mmaxot)
c-----
c           up to this point the dcdh are changes to phase velocity if
c           if the layer boundary changes. Here we change this to mean
c           the dc/dh for a change in layer thickness
c
c           A layer becomes thicker if the base increases and the top
c           decreases its position. The dcdh to this point indicates
c           the effect of moving a boundary down. Now we convert to
c           the effect of changing a layer thickness.
c-----
            do 505 i=1,mmaxot-1
                sum = 0.0
                do 506 j=i+1,mmaxot
                    sum = sum + sdcdh(j)
  506           continue
                sdcdh(i) = sum
  505       continue
            sdcdh(mmaxot) = 0.0

c subroutine putdrt(lun,lorr,wvno,u,gamma,
c     1      sur,sdur,sd2ur,suz,sduz,sd2uz,sare,wvnsrc,sur0,
c     2      rur,rtr,ruz,rtz,rare,wvnrec,rur0,
c     3      sumkr,sumgr,sumgv,mmax,dcdh,
c     4      dcdav,dcdah, dcdbv, dcdbh, dcdn,dcdr,
c     5      ur,tur,uz,tuz,ipar)

C   - OUTPUT eigenfunctions/kernels, COMMENTED!
c            call putdrt(2,6,sngl(wvno),sngl(ugr),
c     1          sngl(gammar),
c     1          sur,sdur,sd2ur,suz,sduz,sd2uz,sare,wvnsrc,sur0,
c     2          rur,rtr,ruz,rtz,rare,wvnrec,rur0,
c     3          sumkr,sumgr,sumgv,mmaxot,
c     4          sdcdh,sdcdav,sdcdah,sdcdbv,sdcdbh,sdcdn,sdcdr,
c     5          spur,sptr,spuz,sptz,ipar)

        do 559 i=1,mmaxot
            ur_out(idisp-1, i)=spur(i)
            tur_out(idisp-1, i)=sptr(i)
            uz_out(idisp-1, i)=spuz(i)
            tuz_out(idisp-1, i)=sptz(i)
            dcdh_out(idisp-1, i)=sdcdh(i)
            dcdav_out(idisp-1, i)=sdcdav(i)
            dcdah_out(idisp-1, i)=sdcdah(i)
            dcdbv_out(idisp-1, i)=sdcdbv(i)
            dcdbh_out(idisp-1, i)=sdcdbh(i)
            dcdn_out(idisp-1, i)=sdcdn(i)
            dcdr_out(idisp-1, i)=sdcdr(i)
  559       continue

        else
C  7 sur = 0
C  8 sdur NaN
C  9 suz 0
C  10 sduz NaN
C  11 rur 0
C  12 rtr 0
C  13 ruz 0
C  14 rtz   0

C   - OUTPUT eigenfunctions, COMMENTED!
c            call putegn(2,2,1,sngl(wvno),sngl(ugr),
c     1          sngl(gammar),
c     1          sur,sdur,suz,sduz,sare,wvnsrc,sur0,
c     2          rur,rtr,ruz,rtz,rare,wvnrec,rur0,
c     3          sumkr,sumgr,sumgv)
        endif
c-----
c       DEBUG OUTPUT
c-----
C           write(6,*)  'tregn96 check3'
C           write(6,*) wvno,c,ugr,are
C           write(6,*) uu(ls), dut
cc-----
c       END DEBUG OUTPUT
c-----
  600   continue
        go to 400
  700   continue
c-----
c       close input file from tdisp96 and output file of this program
c-----
        do 900 i=1,2
                close(i,status='keep')
  900   continue

 9999   continue
        end

c   End of main function
c------------------------------------------------------------------------------c

c-----  add sphere function
c-----
c       This is subroutine tdosph() in VOLVI/src/tspec96
c       modified to work for the Love Rayleigh problem
c
c       Transform spherical earth to flat earth
c
c       Schwab, F. A., and L. Knopoff (1972).
c           Fast surface wave and free
c       mode computations, in
c           Methods in Computational Physics, Volume 11,
c       Seismology: Surface Waves and Earth Oscillations,
c           B. A. Bolt (ed),
c       Academic Press, New York
c
c       Love Wave Equations  44, 45 , 41 pp 112-113
c       Rayleigh Wave Equations 102, 108, 109 pp 142, 144
c
c
c       We will treat all as P-SV for the heck of it
c       This requires more work
c-----
c       mmax    I*4 number of layers
c       TA     R   A
c       TC     R   C
c       TF     R   F
c       TL     R   L
c       TN     R   N
c                  note  density not required
c       TD     R   layer thickness
c       v() R   array of velocities
c       h() R   array of layer thicknesses
c       ipsvsh  I       1 - get P time
c                       2 - get SV time
c                       3 - get SH time
c       refdep R   Reference depth for the model specification
c
c       Note we need the constants here.  Since the velocities
c       must increase with depth, e.g., vf = vs (a/r)
c       and that density  varies
c       as rhof = rhos (a/r)^-P, [not the TI surface wave code has not yet
c        been written], then using the model that m = rho beta^2, we have
c
c       TA = rho VA^2,
c       TAf = rhof * VAf^2 = rhos (a/r)^-P VAs^2 (a/r)^2
c           = (a/r)^2-P TAs
c-----
        subroutine sphere_tdisp96()
        integer ifunc
        integer NL
        parameter (NL=200)

        common/timodel/td(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),qai(NL),qbi(NL),etapi(NL),etasi(NL),
     2      frefpi(NL), frefsi(NL)

        real td, TA, TC, TN, TL, TF, TRho
        real qai, qbi, etapi, etasi, frefpi, frefsi
        common/depref/refdep
        real refdep

        common/pari/ mmax
        integer mmax

        double precision z0,z1,r0,r1,ar,tmp
        ! double precision va,vc,cl,vn

        real radius
        radius = 6371.

        ar=radius
        r0=ar + refdep
        td(mmax)=1.0
        ! For Rayleigh wave ifunc=2
        ifunc=2

        do 10 i=1,mmax
            r1=r0-dble(td(i))
            z0=ar*dlog(ar/r0)
            z1=ar*dlog(ar/r1)
            td(i)=z1-z0
            !-----
            !        attempt 7 15 2007 - use standard rule but at mid layer depth as per DGH
            !-----
            TMP=(ar+ar)/(r0+r1)
            if(ifunc.eq.2)then
                !-----
                !                SV
                !-----
                rhosph = trho(i)
                trho(i) = rhosph * tmp**(-2.275)

                ta(i)=ta(i)*tmp**(-0.2750)
                tc(i)=tc(i)*tmp**(-0.2750)
                !  tf(i)=tf(i)*tmp**(-0.2750)
                elsph = tl(i)
                tl(i)  =elsph*tmp**(-0.2750)
                ensph = tn(i)
                tn(i)=ensph*tmp**(-0.2750)
            else if(ifunc.eq.1)then
                rhosph    = trho(i)
                !  vl = sqrt(tl(i)/rhosph)*tmp
                !  vn = sqrt(tn(i)/rhosph)*tmp
                trho(i) = rhosph * tmp**(-5)
                elsph = tl(i)
                tl(i)=elsph*tmp**(-3.0)
                ensph = tn(i)
                tn(i)=ensph*tmp**(-3.0)
                ! WRITE(6,*)i,tl(i),vl*vl*trho(i),tn(i),vn*vn*trho(i)
            endif
            r0 = r1
  10    continue
        td(mmax)=0.0
        return
        end
c End of sphere_dp
c------------------------------------------------------------------------------c



        subroutine insert(dph,newlyr,ls)
        implicit none
c-----
c       procdure arguments
c-----
        real*4 dph
        logical newlyr
        integer ls
c-----
c       common blocks
c-----
        integer NL
        parameter(NL=200)
        common/timod/  zd(NL),zta(NL),ztc(NL),ztf(NL),
     1      ztl(NL),ztn(NL),zrho(NL),zqai(NL),zqbi(NL),
     2      zetap(NL),zetas(NL),zfrefp(NL),zfrefs(NL)
        real*8 zd, zta, ztc, ztf, ztl, ztn, zrho,
     1      zqai, zqbi, zetap, zetas,
     1      zfrefp, zfrefs
        common/pari/ mmax
        integer mmax
        common/wateri/iwat(NL)
        integer iwat
c-----
c       local variables
c------
        integer LER, LIN, LOT
        parameter (LER=0, LIN=5, LOT=6)
        real*8 dep, dp, dphh, hsave
        integer m, m1

c-----
c       Insert a depth point into the model by splitting a layer
c       so that the point appears at the top boundary of the layer
c       dph = depth of interest
c       newlyr  - L .true. layer added
c               .false. no layer added to get source eigenfunction
c-----
c-----
c       determine the layer in which the depth dph lies.
c       if necessary, adjust  layer thickness at the base
c-----
c       Here determine layer corresponding to specific depth dph
c       If the bottom layer is not thick enough, extend it
c
c       dep     - depth to bottom of layer
c       dphh    - height of specific depth above bottom of the layer
c-----
        dep = 0.0
        dp = 0.0
        dphh=-1.0
        ls = 1
        do 100 m =1, mmax
                dp = dp + zd(m)
                dphh = dp - dph
                if(m.eq.mmax)then
                        if(zd(mmax).le.0.0d+00 .or. dphh.lt.0.0)then
                                zd(mmax) = (dph - dp)
                        endif
                endif
                dep = dep + zd(m)
                dphh = dep - dph
                ls = m
                if(dphh.ge.0.0) go to 101
  100   continue
  101   continue
c-----
c       In the current model, the depth point is in the ls layer
c       with a distance dphh to the bottom of the layer
c
c       Do not create unnecessary layers, e.g.,
c            at surface and internally
c       However do put in a zero thickness layer
c            at the base if necessary
c-----
        if(dph .eq. 0.0)then
            newlyr = .false.
                return
        else if(dphh .eq. 0.0 .and. ls.ne.mmax)then
            ls = ls + 1
            newlyr = .false.
                return
        else
            newlyr = .true.
c-----
c               adjust layering
c-----
                 do 102 m = mmax,ls,-1
                       m1=m+1
                        zd(m1) = zd(m)
                        zta(m1) = zta(m)
                        ztc(m1) = ztc(m)
                        ztf(m1) = ztf(m)
                        ztl(m1) = ztl(m)
                        ztn(m1) = ztn(m)
                        zrho(m1) = zrho(m)
                        zqai(m1) = zqai(m)
                        zqbi(m1) = zqbi(m)
                        zfrefp(m1) = zfrefp(m)
                        zfrefs(m1) = zfrefs(m)
                        zetap(m1) = zetap(m)
                        zetas(m1) = zetas(m)
                        iwat(m1) = iwat(m)
  102           continue
                hsave=zd(ls)
                zd(ls) = hsave - dphh
                zd(ls+1) = dphh
                mmax = mmax + 1
        endif
        return
        end

        subroutine srclyr(depth,lmax)
        implicit none
        real depth
        integer lmax
        integer LER, LIN, LOT
        parameter (LER=0, LIN=5, LOT=6)
        integer NL
        parameter(NL=200)
        common/timod/  zd(NL),zta(NL),ztc(NL),ztf(NL),
     1      ztl(NL),ztn(NL),zrho(NL),zqai(NL),zqbi(NL),
     2      zetap(NL),zetas(NL),zfrefp(NL),zfrefs(NL)
        real*8 zd, zta, ztc, ztf, ztl, ztn, zrho,
     1      zqai, zqbi, zetap, zetas,
     1      zfrefp, zfrefs
        common/pari/ mmax
        integer mmax
        common/wateri/iwat(NL)
        integer iwat

        real dep
        integer i
c-----
c       Find source/receiver boundary. It is assumed that
c       it will lie upon a boundary
c
c       lmax = source layer
c       depth = source depth
c-----
        dep = 0.0
        do 100 i=1,mmax
            if(abs(depth - dep).le.0.001*zd(i))then
                lmax = i
                return
            endif
            dep = dep + zd(i)
  100   continue
        return
        end

        subroutine gcmdln(hsfile,hrfile,hs,hr,dotmp,dogam,dderiv,
     1      nipar,verbose)
c-----
c       parse the command line arguments
c-----
c       hsfile  C*120   - name of source depth file
c       hrfile  C*120   - name of receiver depth file
c       hs  R*4 source depth (single one specified
c       hr  R*4 receiver depth (single one specified
c       dotmp   L   - .true. use file tdisp96.lov
c       dogam   L   - .true. incorporate Q
c       dderiv  L   - .true. output depth dependent values
c       nipar   I*4 - array o integer controls
c           set nipar(4)   = 1 if eigenfunctions are output with -DER flag
c           set nipar(5)   = 1 if dc/dh are output with -DER flag
c           set nipar(6)   = 1 if dc/dav are output with -DER flag
c           set nipar(7)   = 1 if dc/dbv are output with -DER flag
c           set nipar(8)   = 1 if dc/dr are output with -DER flag
c           set nipar(9)   = 1 if dc/dah are output with -DER flag
c           set nipar(10)  = 1 if dc/dn are output with -DER flag
c           set nipar(11)  = 1 if dc/dbh are output with -DER flag
c       verbose L   - .true. output information on energy integrals
c-----
        implicit none
c-----
c       procedure arguments
c-----
        character hrfile*120, hsfile*120
        real hs, hr
        logical dotmp, dogam, dderiv
        integer nipar(20)
        logical verbose

c-----
c       procedure prototypes
c-----
        integer mnmarg
c-----
c       local arguments
c-----
        character name*40
        integer i, nmarg

        hrfile = ' '
        hsfile = ' '
        hs = -1.0E+21
        hr = -1.0E+21
        dotmp = .false.
        dogam = .true.
        dderiv = .false.
        nipar(4) = 0
        nipar(5) = 0
        nipar(6) = 0
        nipar(7) = 0
        nipar(8) = 0
        nipar(9) = 0
        nipar(10) = 0
        nipar(11) = 0
        verbose = .false.
        nmarg = mnmarg()
        i = 0
 1000   continue
            i = i + 1
            if(i.gt.nmarg)go to 2000
                call mgtarg(i,name)
                if(name(1:3).eq.'-HR')then
                    i = i + 1
                    call mgtarg(i,name)
                    read(name,'(bn,f20.0)')hr
                else if(name(1:3).eq.'-HS')then
                    i = i + 1
                    call mgtarg(i,name)
                    read(name,'(bn,f20.0)')hs
C               else if(name(1:4).eq.'-FHR')then
C                   i = i + 1
C                   call mgtarg(i,hrfile)
C               else if(name(1:4).eq.'-FHS')then
C                   i = i + 1
C                   call mgtarg(i,hsfile)
                else if(name(1:2) .eq. '-T')then
                    dotmp = .true.
                else if(name(1:4) .eq. '-NOQ')then
                    dogam = .false.
                else if(name(1:4) .eq. '-DER')then
                    nipar(4) = 1
                    nipar(5) = 1
                    nipar(6) = 1
                    nipar(7) = 1
                    nipar(8) = 1
                    nipar(9) = 1
                    nipar(10) = 1
                    nipar(11) = 1
                    dderiv = .true.
                else if(name(1:3).eq.'-DE')then
                    dderiv = .true.
                    nipar(4) = 1
                else if(name(1:3).eq.'-DH')then
                    dderiv = .true.
                    nipar(5) = 1
                else if(name(1:3).eq.'-DB')then
                    dderiv = .true.
                    nipar(7) = 1
                else if(name(1:3).eq.'-DR')then
                    dderiv = .true.
                    nipar(8) = 1
                else if(name(1:2).eq.'-V')then
                    verbose = .true.
                else if(name(1:2) .eq. '-?')then
                    call usage(' ')
                else if(name(1:2) .eq. '-h')then
                    call usage(' ')
                endif
        go to 1000
 2000   continue
        return
        end

        subroutine usage(ostr)
        implicit none
c-----
c       procedure arguments
c-----
        character ostr*(*)
c-----
c       local variables
c-----
        integer LER
        parameter (LER=0)
        integer lostr
c-----
c       function prototype
c-----
        integer lgstr

        if(ostr .ne. ' ')then
            lostr = lgstr(ostr)
            write(LER,*)ostr(1:lostr)
        endif
        write(LER,*)'Usage: tregn96 ',
C    1      ' -FHR recdep -FHS srcdep -HS hs -HR hr ' ,
     1      ' -HS hs -HR hr ' ,
     2       ' [-NOQ] [-T] [-DER -DE -DH -DB -DR ] [-?] [-h]'
C       write(LER,*)
C     1 '-FHS srcdep (overrides -HS )  Name of source depth  file'
C       write(LER,*)
C     1 '-FHR recdep (overrides -HR )  Name of receiver depth  file'
        write(LER,*)
     1  '-HS hs      (default 0.0 )  Source depth '
        write(LER,*)
     1  '-HR hr      (default 0.0 )  Receiver depth'
        write(LER,*)
     1  '-NOQ        (default Q used) perfectly elastic'
        write(LER,*)
     1  '-T          (default false) use tdisp96.lov not disp96.lov'
        write(LER,*)
     1  '-DER        (default false) output depth dependent values'
        write(LER,*)
     1  '-DE         (default false) output eigenfunctions(depth)'
        write(LER,*)
     1  '-DH         (default false) output DC/DH(depth)'
        write(LER,*)
     1  '-DB         (default false) output DC/DB(depth)'
        write(LER,*)
     1  '-DR         (default false) output DC/DR(depth)'
        write(LER,*)
     1  '-V          (default false) list energy integrals'
        write(LER,*)
     1  '-?          (default none )  this help message '
        write(LER,*)
     1  '-h          (default none )  this help message '
        stop
        end

        subroutine collap(ls,mmaxot)
        implicit none
c-----
c       routine arguments
c-----
        integer ls, mmaxot
c-----
c       common blocks
c-----
        integer NL
        parameter(NL=200)
        common/eigfun/ ur(NL),uz(NL),tz(NL),tr(NL),uu0(4),
     *    dcdah(NL),dcdav(NL),dcdbh(NL),dcdbv(NL),
     2    dcdn(NL),dcdr(NL),dcdh(NL)
        real*8 ur, uz, tz, tr, uu0,
     1     dcdah, dcdav, dcdbh, dcdbv, dcdn, dcdr, dcdh
c-----
c       local arguments
c-----
        integer i

        do 501 i = ls-1,mmaxot
            if(i .eq. ls -1)then
                dcdah(i) = dcdah(i) + dcdah(i+1)
                dcdav(i) = dcdav(i) + dcdav(i+1)
                dcdbh(i) = dcdbh(i) + dcdbh(i+1)
                dcdbv(i) = dcdbv(i) + dcdbv(i+1)
                dcdr(i) = dcdr(i) + dcdr(i+1)
                dcdn(i) = dcdn(i) + dcdn(i+1)
            endif
            if(i.gt.ls)then
                dcdah(i-1) = dcdah(i)
                dcdav(i-1) = dcdav(i)
                dcdbh(i-1) = dcdbh(i)
                dcdbv(i-1) = dcdbv(i)
                dcdh(i-1) = dcdh(i)
                dcdr(i-1) = dcdr(i)
                dcdn(i-1) = dcdn(i)
                ur(i-1) = ur(i)
                tr(i-1) = tr(i)
                uz(i-1) = uz(i)
                tz(i-1) = tz(i)
            endif
  501   continue
        mmaxot = mmaxot - 1
        return
        end

        subroutine chksiz(dp,sp,mmaxot)
c-----
c       correctly convert double precision to single precision
c-----
        real*8 dp(mmaxot)
        real*4 sp(mmaxot)
            do 610 i=1,mmaxot
                if(dabs(dp(i)).lt.1.0d-36)then
                    sp(i) = 0.0
                else
                    sp(i) = sngl(dp(i))
                endif
  610       continue
        return
        end

        subroutine sprayl(om,c,mmax,csph,usph,ugr)
c-----
c       Transform spherical earth to flat earth
c       and relate the corresponding flat earth dispersion to spherical
c
c       Schwab, F. A., and L. Knopoff (1972).
c          Fast surface wave and free
c       mode computations, in
c          Methods in Computational Physics, Volume 11,
c       Seismology: Surface Waves and Earth Oscillations,
c          B. A. Bolt (ed),
c       Academic Press, New York
c
c       Rayleigh Wave Equations 111, 114 p 144
c
c       Partial with respect to parameter uses the relation
c       For phase velocity, for example,
c
c       dcs/dps = dcs/dpf * dpf/dps, c is phase velocity, p is
c       parameter, and f is flat model
c
c       om      R*8     angular frequency
c       c       R*8     phase velocity
c       mmax    I*4     number of layers
c-----
        implicit none
c-----
c       procedure arguments
c-----
        real*8 om, c, csph,usph,ugr
        integer mmax
c-----
c       common blocks
c-----
        integer NL
        parameter(NL=200)
        common/eigfun/ ur(NL),uz(NL),tz(NL),tr(NL),uu0(4),
     *    dcdah(NL),dcdav(NL),dcdbh(NL),dcdbv(NL),
     2    dcdn(NL),dcdr(NL),dcdh(NL)
        real*8 ur, uz, tz, tr, uu0,
     1     dcdah, dcdav, dcdbh, dcdbv, dcdn, dcdr, dcdh
        common/sphereL/vtp(NL),dtp(NL),rtp(NL)
        real*4 vtp,dtp,rtp
c-----
c       local arguments
c-----
        real*8 ar, tm
        integer i

        ar = 6370.0d0
        tm=sqrt(1.+(c/(2.*ar*om))**2)
            do 20 i=1,mmax
                dcdah(i)=dcdah(i)*  vtp(i)/(tm**3)
                dcdav(i)=dcdav(i)*  vtp(i)/(tm**3)
                dcdbh(i)=dcdbh(i)*  vtp(i)/(tm**3)
                dcdbv(i)=dcdbv(i)*  vtp(i)/(tm**3)
C  check dcdn mapping
C                dcdh(i)=dcdh(i)*  dtp(i)/(tm**3)
                dcdr(i)=dcdr(i)*  rtp(i)/(tm**3)
   20       continue
c       write(6,*)'c flat=',c,' csph=',c/tm
        usph = ugr*tm
        csph = c/tm
        return
        end

        subroutine bldsph()
c-----
c       Transform spherical earth to flat earth
c
c       Schwab, F. A., and L. Knopoff (1972).
c          Fast surface wave and free
c       mode computations, in
c          Methods in Computational Physics, Volume 11,
c       Seismology: Surface Waves and Earth Oscillations,
c          B. A. Bolt (ed),
c       Academic Press, New York
c
c       Love Wave Equations  44, 45 , 41 pp 112-113
c       Rayleigh Wave Equations 102, 108, 109 pp 142, 144
c
c     Revised 28 DEC 2007 to use mid-point, assume linear variation in
c     slowness instead of using average velocity for the layer
c     Use the Biswas (1972:PAGEOPH 96, 61-74, 1972) density mapping
c
c       This is for Rayleigh Waves
c
c-----
        implicit none
        integer NL
        parameter (NL=200)
c-----
c       common blocks
c-----
        common/sphereL/vtp(NL),dtp(NL),rtp(NL)
        real*4 vtp,dtp,rtp
        common/timod/  zd(NL),zta(NL),ztc(NL),ztf(NL),
     1      ztl(NL),ztn(NL),zrho(NL),zqai(NL),zqbi(NL),
     2      zetap(NL),zetas(NL),zfrefp(NL),zfrefs(NL)
        real*8 zd, zta, ztc, ztf, ztl, ztn, zrho,
     1      zqai, zqbi, zetap, zetas,
     1      zfrefp, zfrefs
        common/pari/ mmax
        integer mmax
c-----
c       local arguments
c-----
        real*8 ar, dr, r0, r1, z0, z1
        real*8 tmp
        integer i
c-----
c       vtp is the factor used to convert spherical
c          velocity to flat velocity
c       rtp is the factor used to convert spherical
c          density  to flat density
c       dtp is the factor used to convert spherical
c          boundary to flat boundary
c-----
c-----
c       duplicate computations of srwvds(IV)
c-----
        ar=6370.0d0
        dr=0.0d0
        r0=ar
        z0 = 0.0d+00
        zd(mmax)=1.0
        do 10 i=1,mmax
            r1 = r0 * dexp(-zd(i)/ar)
            if(i.lt.mmax)then
                z1 = z0 + zd(i)
            else
                z1 = z0 + 1.0d+00
            endif
            TMP=(ar+ar)/(r0+r1)
            vtp(i) = TMP
            rtp(i) = TMP**(-2.275)
            dtp(i) = sngl(ar/r0)
            r0 = r1
            z0 = z1
   10   continue
C        write(6,*)'vtp:',(vtp(i),i=1,mmax)
C        write(6,*)'rtp:',(rtp(i),i=1,mmax)
C        write(6,*)'dtp:',(dtp(i),i=1,mmax)
c-----
c       at this point the model information is no longer used and
c       will be overwritten
c-----
        return
        end

        function ffunc(nub,dm)
        implicit none
        complex*16 ffunc
        complex*16 nub
        real*8 dm
        complex*16 exqq
        complex*16 argcd
c-----
c       get the f function
c-----
        if(cdabs(nub).lt. 1.0d-08)then
            ffunc = dm
        else
            argcd = nub*dm
            if(dreal(argcd).lt.40.0)then
                  exqq = cdexp(-2.0d+00*argcd)
            else
                  exqq=0.0d+00
            endif
            ffunc = (1.d+00-exqq)/(2.d+00*nub)
        endif
        return
        end

        function gfunc(nub,dm)
        implicit none
        complex*16 gfunc
        complex*16 nub
        real*8 dm
        complex*16 argcd
        argcd = nub*dm
        if(dreal(argcd).lt.75)then
             gfunc = cdexp(-argcd)*dm
        else
             gfunc = dcmplx(0.0d+00, 0.0d+00)
        endif
        return
        end

        function h1func(nua,nub,dm)
        implicit none
        complex*16 h1func
        complex*16 nua, nub
        real*8 dm
        complex*16 argcd
        complex*16 exqq
        if(cdabs(nub+nua).lt. 1.0d-08)then
            h1func = dm
        else
            argcd = (nua+nub)*dm
            if(dreal(argcd).lt.40.0)then
                  exqq = cdexp(-argcd)
            else
                  exqq=0.0d+00
            endif
            h1func = (1.d+00-exqq)/(nub+nua)
        endif
        return
        end

        function h2func(nua,nub,dm)
        implicit none
        complex*16 h2func
        complex*16 nua, nub
        real*8 dm
        complex*16 argcd
        complex*16 exqq, exqp
        if(cdabs(nub-nua).lt. 1.0d-08)then
c-----
c           this should never occur for surface waves
c-----
            h2func =  dm
        else
            argcd = nua*dm
            if(dreal(argcd).lt.40.0)then
                  exqp = cdexp(-argcd)
            else
                  exqp=0.0d+00
            endif
            argcd = nub*dm
            if(dreal(argcd).lt.40.0)then
                  exqq = cdexp(-argcd)
            else
                  exqq=0.0d+00
            endif
            h2func = (exqq - exqp)/(nua-nub)
        endif
        return
        end

        subroutine cnormc(ee,ex,nmat)
c-----
c       This routine is an important step to control over- or
c       underflow.
c       The Haskell or Dunkin vectors are normalized before
c       the layer matrix stacking.
c       Note that some precision will be lost during normalization.
c
c       This is use with compound and has complex arguments
c-----
        implicit none
        complex*16 ee(nmat)
        real*8 ex
        integer nmat

        real*8 t1
        complex*16 ztmp
        integer i
        ex = 0.0d+00
        t1 = 0.0d+00
        do  i = 1,nmat
          if(cdabs(ee(i)).gt.t1) then
               t1 = cdabs(ee(i))
          endif
        enddo
        if(t1.lt.1.d-40) t1=1.d+00
        do i =1,nmat
          ztmp=ee(i)
          ztmp=ztmp/t1
          ee(i)=ztmp
        enddo
c-----
c       store the normalization factor in exponential form.
c-----
        ex=dlog(t1)
        return
        end

        subroutine rnormc(ee,ex,nmat)
c-----
c       This routine is an important step to control over- or
c       underflow.
c       The Haskell or Dunkin vectors are normalized before
c       the layer matrix stacking.
c       Note that some precision will be lost during normalization.
c
c       This is used with hska which returns real values
c-----
        implicit none
        real*8 ee(nmat)
        real*8 ex
        integer nmat

        real*8 t1
        real*8 ztmp
        integer i
        ex = 0.0d+00
        t1 = 0.0d+00
        do  i = 1,nmat
          if(dabs(ee(i)).gt.t1) then
               t1 = dabs(ee(i))
          endif
        enddo
        if(t1.lt.1.d-40) t1=1.d+00
        do i =1,nmat
          ztmp=ee(i)
          ztmp=ztmp/t1
          ee(i)=ztmp
        enddo
c-----
c       store the normalization factor in exponential form.
c-----
        ex=dlog(t1)
        return
        end

c-----
c      this section of code is specific to the isotropic problem
c      the TI code will use the same subroutine calls. The
c      conversion of isotropic to TI involves replacement of
c      subroutines and/or pieces of code, especially the common
c      block describing the model parameters
c-----

      subroutine svfunc(omega,wvno)
c
c     This combines the Haskell vector from sub down and
c     Dunkin vector from sub up to form the eigenfunctions.
c
      implicit none
c-----
c     command arguments
c-----
        real*8 omega, wvno
c-----
c     common blocks
c-----
        integer NL
        parameter (NL=200)
        common/timod/  zd(NL),zta(NL),ztc(NL),ztf(NL),
     1      ztl(NL),ztn(NL),zrho(NL),zqai(NL),zqbi(NL),
     2      zetap(NL),zetas(NL),zfrefp(NL),zfrefs(NL)
        real*8 zd, zta, ztc, ztf, ztl, ztn, zrho,
     1      zqai, zqbi, zetap, zetas,
     1      zfrefp, zfrefs
        common/pari/ mmax
        integer mmax
        common/wateri/iwat(NL)
        integer iwat
        common/eigfun/ ur(NL),uz(NL),tz(NL),tr(NL),uu0(4),
     *    dcdah(NL),dcdav(NL),dcdbh(NL),dcdbv(NL),
     2    dcdn(NL),dcdr(NL),dcdh(NL)
        real*8 ur, uz, tz, tr, uu0,
     1     dcdah, dcdav, dcdbh, dcdbv, dcdn, dcdr, dcdh
        common/dunk/   cd(NL,5),exe(NL),exa(NL)
        complex*16 cd
        real*8 exe, exa
        common/hask/   vv(NL,4)
        real*8 vv
c-----
c       model characterization
c-----
        common/modspec/allfluid
        logical allfluid

c-----
c       unknown common check for clean up
c-----
        common/hwat/wh(NL,2),hex(NL)
        real*8 wh, hex
c-----
c       internal variables
c-----
        real*8 fr
        integer i
        real*8 uu1, uu2, uu3, uu4, ext, f1213, fact
        real*8 cd1, cd2, cd3, cd4, cd5, cd6
        real*8 tz1, tz2, tz3, tz4
        integer jwat
c-----
c       get compound matrix from bottom to top
c-----
        call up(omega,wvno,fr)
        call down(omega,wvno)


c-----
c       get propagator from top to bottom
c-----
        f1213 = -cd(1,2)
        ur(1) = cd(1,3)/cd(1,2)
        uz(1) = 1.0d+00
        tz(1) = 0.0d+00
        tr(1) = 0.0d+00
        uu0(1) = ur(1)
        uu0(2) = 1.0
        uu0(3) = fr
        uu0(4) = fr
C       WRITE(6,*)'ur(1):',ur(1)
c------
c      following
c                -1|kl    -1
c       Bk =    X  |     Z
c                  |34    4l
c             -----------------
c                 -1|24
c                R  |
c                   |34
c-----
c       Compound Matrix (6x6)   Here (5x5)
c
c       X|12/12 ==               cd(1)
c       X|12/13 ==               cd(2)
c       X|12/14 ==               cd(3)
c       X|12/23 ==              -cd(3)
c       X|12/24 ==               cd(4)
c       X|12/34 ==               cd(5)
c       X|12/21 =  -X|12/12
c       X|12/ij =  -X|12/ji
c
c       where
c
c        |12
c       X|  = X|12/ij == X(1,i)*X(2,j) - X(1,j)*X(2,i)
c        |ij

c   and using the shorthand
c                  6x6 matrix                    a5x5 reduced matrix
c       True [11 12 13 14 15 16 ] = Reduced [ 11 12 13 -13 14 15 ]
c     we have
c-----
c   In long hand with full notation and using the complex matrix short hand, e.g.,
C   that the denominator is -R12 and recalling R11 == 0
c
c
c    |    |         |                               | |      |
c    |    |         |           12      12      12  | |   -1 |
c    | B1 |         |   0      X|     -X|      X|   | |  Z   |
c    |    |         |           34      24      23  | |   41 |
c    |    |         |                               | |      |
c    |    |         |   12              12      12  | |   -1 |
c    | B2 |         | -X|       0      X|     -X|   | |  Z   |
c    |    |     1   |   34              14      13  | |   42 |
c    |    | = _____ |                               | |      |
c    |    |      12 |   12      12              12  | |   -1 |
c    | B3 |   -R|   |  X|     -X|       0      X|   | |  Z   |
c    |    |      13 |   24      14              12  | |   43 |
c    |    |         |                               | |      |
c    |    |         |   12      12      12          | |   -1 |
c    | B4 |         | -X|      X|     -X|       0   | |  Z   |
c    |    |         |   23      13      12          | |   44 |
c
c---------
c    a property of the propagator matrix is that
c    the elements of Z^-1 are related to the elements of Z as
c
c                   |  Z     Z    -Z    -Z     |
c                   |   44    34    24    14   |
c                   |                          |
c                   |  Z     Z    -Z    -Z     |
c        -1         |   43    33    23    13   |
c      Z      =     |                          |
c                   | -Z    -Z     Z     Z     |
c                   |   42    32    22    12   |
c                   |                          |
c                   | -Z    -Z     Z     Z     |
c                   |   41    31    21    11   |
c-----
c    Here the vv are the first column of the  Haskell down, and the
c    cd are the first row if the X up. At the surface R=X
c------
        ext = 0.0
        do 200 i=2,mmax
         cd1=  cd(i,1)
         cd2=  cd(i,2)
         cd3=  cd(i,3)
         cd4= -cd(i,3)
         cd5=  cd(i,4)
         cd6=  cd(i,5)
         tz1 = -vv(i,4)
         tz2 = -vv(i,3)
         tz3 =  vv(i,2)
         tz4 =  vv(i,1)

            uu1 =
     1                  tz2*cd6 - tz3*cd5 + tz4*cd4
            uu2 =
     1        -tz1*cd6          + tz3*cd3 - tz4*cd2
            uu3 =
     1         tz1*cd5 - tz2*cd3          + tz4*cd1
            uu4 =
     1        -tz1*cd4 + tz2*cd2 - tz3*cd1

            ext=exa(i)+exe(i) -exe(1)
            if(ext.gt.-80.0 .and. ext.lt.80.0 ) then
                fact=dexp(ext)

                ur(i)=uu1*fact/f1213
                uz(i)=uu2*fact/f1213
                tz(i)=uu3*fact/f1213
                tr(i)=uu4*fact/f1213
C        WRITE(6,*)'i,B:',i,ur(i),uz(i),tz(i),tr(i)
            else
                ur(i) = 0.0
                uz(i) = 0.0
                tz(i) = 0.0
                tr(i) = 0.0
            endif
  200   continue
c-----
c       correction for fluid layers on top if not
c       all fluid
c-----
        if(.not. allfluid)then
                jwat = 0
                do 300 i=1,mmax
                    if(iwat(i).gt.0)then
                            jwat = i
                    else
                            go to 301
                    endif
  300   continue
  301   continue
                if(jwat.gt.0)then
                        do i=1,jwat
                            ur(i) = 0.0
                            tr(i) = 0.0
                        enddo
                endif
        endif

CRBHc-----
CRBHc       Continue KLUDGE if top layers are water
CRBHc       The problem here was that even though the
CRBHc          phase velocity is determined,
CRBHc       the computed eigenfunctions do not match perfectly,
CRBHc       and we do not match the free surface requirement that Tz=0 at
CRBHc       top of fluid. So RBH used the fluid propagator from top
CRBHc       down to the solid-fluid interface, the Haskell from
CRBHc       base up to the same interface, and then adjust the amplitudes
CRBHc       to match. The purpose is to get a better estimate of the
CRBHc       eigenfunction in the fluid.
CRBHc
CRBHc       Ultimately use Kennett
CRBHc-----
CRBHc       first get the largest exponent
CRBHc-----
CRBH        if(iwat(1).eq.1)then
CRBHc-----
CRBHc           OK first layer is fluid
CRBHc-----
CRBH            hmax = 0.0
CRBH            jwat=0
CRBH            do 300 i=1,mmax-1
CRBH                if(iwat(i).eq.1)then
CRBH                    jwat = i
CRBH                    hmax = hex(i+1)
CRBH                endif
CRBH  300       continue
CRBHc-----
CRBHc           now get correction factor
CRBHc-----
CRBH            corrd = sqrt(uz(jwat+1)**2 + tz(jwat+1)**2 )
CRBH            corrw = sqrt(wh(jwat+1,1)**2 + wh(jwat+1,2)**2 )
CRBH            if(corrd.eq.0.0)then
CRBHc-----
CRBHc       give precedence to surface values
CRBHc-----
CRBH                fac = 0.0
CRBH            else
CRBH                fac = corrw/corrd
CRBH            endif
CRBH            fac = fac * exp(hmax - hmax)
CRBH            fac = fac * 2.0
CRBHc-----
CRBHc           eventually normalize
CRBHc-----
CRBH        IF(.not. ALLFLUID)then
CRBH            do 400 i=1,mmax
CRBH                if(i.le.jwat)then
CRBH                    efac = exp(hex(i) - hmax)
CRBH                    efac = efac * 2.0
CRBH                    ur(i) = 0.0
CRBH                    uz(i) = wh(i,1)*efac
CRBH                    tz(i) = wh(i,2)*efac
CRBH                    tr(i) = 0.0
CRBH                else
CRBH                    ur(i) = ur(i) * fac
CRBH                    uz(i) = uz(i) * fac
CRBH                    tr(i) = tr(i) * fac
CRBH                    tz(i) = tz(i) * fac
CRBH                endif
CRBH  400       continue
CRBH        ENDIF
CRBH        endif

        return
        end

        subroutine up(omega,wvno,fr)
c-----
c       This finds the values of the Dunkin vectors at
c       each layer boundaries from bottom layer upward.
c
c       Note that these are the reduced vectors, e.g.,
c       The true compound (G An-1 ... A1 H) is
c       True [11 12 13 14 15 16 ] = Reduced [ 11 12 13 -13 14 15 ]
c-----
        implicit none
c-----
c       command line arguments
c-----
        real*8 omega, wvno, fr
c-----
c       common blocks
c-----
        integer NL
        parameter (NL=200)
        common/timod/  zdd(NL),zta(NL),ztc(NL),ztf(NL),
     1      ztl(NL),ztn(NL),zrho(NL),zqai(NL),zqbi(NL),
     2      zetap(NL),zetas(NL),zfrefp(NL),zfrefs(NL)
        real*8 zdd, zta, ztc, ztf, ztl, ztn, zrho,
     1      zqai, zqbi, zetap, zetas,
     1      zfrefp, zfrefs
        common/pari/ mmax
        integer mmax
        common/wateri/iwat(NL)
        integer iwat
        common/eigfun/ ur(NL),uz(NL),tz(NL),tr(NL),uu0(4),
     *    dcdah(NL),dcdav(NL),dcdbh(NL),dcdbv(NL),
     2    dcdn(NL),dcdr(NL),dcdh(NL)
        real*8 ur, uz, tz, tr, uu0,
     1     dcdah, dcdav, dcdbh, dcdbv, dcdn, dcdr, dcdh
        common/dunk/   cd(NL,5),exe(NL),exa(NL)
        complex*16 cd
        real*8 exe, exa
c-----
c       model characterization
c-----
        common/modspec/allfluid
        logical allfluid
c-----
c       unknown common check for clean up
c-----
        common/save/   dd(5,5),aa(4,4),ex1,ex2
        real*8 dd, aa, ex1, ex2
        common/aamatx/ ww(NL),xx(NL),yy(NL),zz(NL),
     1      cospp(NL),cosqq(NL)
        real*8 ww, xx, yy, zz, cospp, cosqq
        common/engerw/ wra,wd,wba
        real*8 wra,wd,wba
c-----
c       internal arguments
c-----

c-----
c       The save labeled common is used to pass the
c       Dunkin 5x5 and Haskell 4x4 matrices to the main program
c       exp(ex1) is the scaling for the Dunkin-Thrower compound matrix
c       exp(ex2) is the scaling for the Haskell matrices
c-----

        complex*16 rp, rsv, p, q
        real*8 wvno2, om2
        complex*16 gbr(2,5)
        complex*16 cr, ee(5)
        real*8  exn
        complex*16 ca(5,5)
        real*8 pex, svex

        complex*16 cosp,rsinp,sinpr,cossv,rsinsv,sinsvr
        real*8 exsum
        complex*16 NP, NSV
        COMPLEX*16 X11, X21, X31, X41
        COMPLEX*16 X12, X22, X32, X42
        complex*16 Za,Zb,Zc,Zd,Ze,Zf
c-----
c       norms

        integer i,j,m,mmm1,nmat
c-----
c       initialize base vector
c-----
c-----
c       set up starting values for bottom halfspace
c-----
        wvno2=wvno*wvno
        om2 = omega*omega
        call evalg(0,mmax,mmax-1,gbr,1,
     1      wvno,omega,om2,wvno2)
c-----
c       note for surface waves, the ra, rb are real for the halfspace as long
c       as the phase velocity is less than the respective halfspace velocity
c-----

            cd(mmax,1)= dreal(gbr(1,1))
            cd(mmax,2)= dreal(gbr(1,2))
            cd(mmax,3)= dreal(gbr(1,3))
            cd(mmax,4)= dreal(gbr(1,4))
            cd(mmax,5)= dreal(gbr(1,5))
            exe(mmax) = 0.0d+00
C         WRITE(6,*)'gbr(1,1):',gbr(1,1)
C         WRITE(6,*)'gbr(1,2):',gbr(1,2)
C         WRITE(6,*)'gbr(1,3):',gbr(1,3)
C         WRITE(6,*)'gbr(1,4):',gbr(1,4)
C         WRITE(6,*)'gbr(1,5):',gbr(1,5)
c------
c       matrix multiplication from bottom layer upward
c------
        mmm1=mmax-1
        exsum = 0.0
        do 500 m = mmm1,1,-1
           call  gettiegn(Za,Zb,Zc,Zd,Ze,Zf,om2,wvno2,rp, rsv,
     1      x11,x21,x31,x41,x12,x22,x32,x42,NP, NSV,m,omega,wvno,
     2      dcmplx(1.0d+00,0.0d+00),dcmplx(1.0d+00,0.0d+00))
                p = rp  * zdd(m)
                q = rsv * zdd(m)
                call varsv(p,q,rp, rsv,
     1              cosp, cossv, rsinp, rsinsv,
     1              sinpr, sinsvr, pex,svex,iwat(m),zdd(m))
C       WRITE(6,*)'up   :rp,rsv,pex,svex:',rp,rsv,pex,svex
                call dnka_tregn(CA,cosp,rsinp,sinpr,cossv,rsinsv,
     1              sinsvr, NP,NSV,
     1              x11,x21,x31,x41,x12,x22,x32,x42,
     1              sngl(ZRho(m)),iwat(m),pex+svex,om2)

C                call dnka(CA,cosp,rsinp,sinpr,cossv,rsinsv,sinsvr,
C     1              NP,NSV,
C     1              x11,x21,x31,x41,x12,x22,x32,x42,
C     1              sngl(ZRho(m)),iwat(m),pex+svex,om2)
C



            nmat = 5
            do 200 i=1,nmat
                cr=0.0d+00
                do 100 j=1,nmat
                    cr=cr+cd(m+1,j)*ca(j,i)
  100           continue
                ee(i)=cr
  200       continue
            exn= 0.0d+00
            call cnormc(ee,exn,nmat)
            exsum = exsum + pex + svex + exn
            exe(m) = exsum
            do 300 i = 1,nmat
                cd(m,i)=ee(i)
  300       continue
  500   continue
c-----
c       define period equation
c-----
        fr=cd(1,1)
        return
        end

        subroutine dnka_tregn(CA,cosp,rsinp,sinpr,cossv,rsinsv,sinsvr,
     1      NP,NSV, x11,x21,x31,x41,x12,x22,x32,x42,
     1      TRho,iwat,ex,om2)
        implicit none
c   change name from dnka to dnka_tregn
c-----
c       command line variables
c-----
        complex*16 CA(5,5)
        complex*16 cosp,rsinp,sinpr,cossv,rsinsv,sinsvr
        complex*16 NP,NSV
        complex*16 x11,x21,x31,x41,x12,x22,x32,x42
        real*4 Trho
        integer iwat
        real*8 ex, om2

c-----
c       internal variables
c-----
        complex*16 tca(6,6)
        complex*16 x(4,2)
        complex*16 a1, a2
        complex*16 c1,ls1,s1l,c2,ls2,s2l
        real*8 dfac
        integer i,j

c-----
c       introduce conciseness to reduce operations
c-----
c-----
c        A11     A12     A13    -A13     A15     A16
c        A21     A22     A23    -A23     A25     A15
c        A31     A32     A33    1-A33   -A23    -A13
c       -A31    -A32    1-A33    A33     A23     A13
c        A51     A52    -A32     A32     A22     A12
c        A61     A51    -A31     A31     A21     A11
c-----
c       this will be multipled on the left by the G matrix
c
c       [ G11   G12 G13 -G13    G15 G16 ]
c
c-----
c       or on the right by
c
c       [ H11   H21 H31 -H31    H51 H61  ] ^T
c-----
c       the number of multiplications can be reduced from
c            36 to 25 if we define a new matrices
c       related to the original matrices by
c-----
c         A11     A12     A13         A15     A16
c         A21     A22     A23         A25     A15
c       2 A31   2 A32    2 A33 -1   -2 A23  -2 A13
c         A51     A52    -A32         A22     A12
c         A61     A51    -A31         A21     A11
c-----
c
c       [ G11   G12  G13    G15 G16  ]
c       [ H11   H21 2 H31   H51 H61  ] ^T
c
c-----
c       this means that some of the original definitions of the
c           Aij elements must be changed for the
c       definition of the modified 5x5 compount A matrix
c
c       old 6x6                 new 5x5
c       A11 = 1 - A33 + A22     1 - (1/2)(new A33 + 1) + new A2)
c       A53 = -A32              A43 = - (1/2) new A32
c       A63 = -A31              A53 = - (1/2) new A31
c-----
c       To recover the needed elements, we note that the
c           old G14 = -old G14 = new G13
c-----

        if(ex.gt.35.0d+00)then
            dfac = 0.0d+00
        else
            dfac = dexp(-ex)
        endif
        if(iwat.eq.1)then
c-----
c       fluid layer
c-----
            do 100 j=1,5
                do 101 i=1,5
                    ca(i,j) = 0.0d+00
  101           continue
  100       continue
            ca(3,3) = dfac
            ca(1,1) = cosp
            ca(5,5) = cosp
            ca(1,2) = -rsinp/(Trho*om2)
            ca(2,1) = - Trho*sinpr*om2
            ca(2,2) = cosp
            ca(4,4) = cosp
            ca(4,5) = ca(1,2)
            ca(5,4) = ca(2,1)
        else
            a1 = 0.5/np
            a2 = 0.5/nsv
            c1  = 2.*a1*cosp
            ls1 = 2.*a1*rsinp
            s1l = 2.*a1*sinpr
            c2  = 2.*a2*cossv
            ls2 = 2.*a2*rsinsv
            s2l = 2.*a2*sinsvr
            x(1,1) = x11
            x(2,1) = x21
            x(3,1) = x31
            x(4,1) = x41
            x(1,2) = x12
            x(2,2) = x22
            x(3,2) = x32
            x(4,2) = x42
C           WRITE(6,*)'a1,a2:',a1,a2
C           WRITE(6,*)'c1,ls1,s1l:',c1,ls1,s1l
C           WRITE(6,*)'c2,ls2,s2l:',c2,ls2,s2l


c-------------------------------------------------------------
C      if(i.eq.1 .and. j.eq.1)then
C      zout =
C    1  (
C    1     (+1)*x(1,1)*x(4,1)* (-1)*x(2,1)*x(3,1)* c1 * c1   +
C    1     (+1)*x(1,1)*x(4,1)* (-1)*x(2,2)*x(3,2)* c1 * c2   +
C    1     (+1)*x(1,2)*x(4,2)* (-1)*x(2,1)*x(3,1)* c2 * c1   +
C    1     (+1)*x(1,2)*x(4,2)* (-1)*x(2,2)*x(3,2)* c2 * c2
C    1   )  - (
C    1     (-1)*x(1,1)*x(3,1)* (+1)*x(2,1)*x(4,1)* s1l* ls1  +
C    1     (-1)*x(1,1)*x(3,1)* (+1)*x(2,2)*x(4,2)* s1l* s2l  +
C    1     (-1)*x(1,2)*x(3,2)* (+1)*x(2,1)*x(4,1)* ls2* ls1  +
C    1     (-1)*x(1,2)*x(3,2)* (+1)*x(2,2)*x(4,2)* ls2* s2l
C    1   )
       tca(1,1) =
     1  (
c    1     (+1)*x(1,1)*x(4,1)* (-1)*x(2,1)*x(3,1)* c1 * c1   +
c    1     (+1)*x(1,1)*x(3,1)* (+1)*x(2,1)*x(4,1)* s1l* ls1  +
     1     (-1)*x(1,1)*x(3,1)* (+1)*x(2,1)*x(4,1)*dfac*a1*a1*4 +
     1     (+1)*x(1,1)*x(4,1)* (-1)*x(2,2)*x(3,2)* c1 * c2   +
     1     (+1)*x(1,2)*x(4,2)* (-1)*x(2,1)*x(3,1)* c2 * c1   +
     1     (+1)*x(1,1)*x(3,1)* (+1)*x(2,2)*x(4,2)* s1l* s2l  +
     1     (+1)*x(1,2)*x(3,2)* (+1)*x(2,1)*x(4,1)* ls2* ls1  +
c    1     (+1)*x(1,2)*x(4,2)* (-1)*x(2,2)*x(3,2)* c2 * c2   +
c    1     (+1)*x(1,2)*x(3,2)* (+1)*x(2,2)*x(4,2)* ls2* s2l
     1     (-1)*x(1,2)*x(3,2)* (+1)*x(2,2)*x(4,2)*dfac*a2*a2*4
     1   )
c-------------------------------------------------------------
C      else if(i.eq.1 .and. j.eq.2)then
C      zout =
C    1  (
C    1     (+1)*x(1,1)*x(4,1)* (-1)*x(2,1)*x(2,1)* c1 * ls1  +
C    1     (+1)*x(1,1)*x(4,1)* (-1)*x(2,2)*x(2,2)* c1 * s2l  +
C    1     (+1)*x(1,2)*x(4,2)* (-1)*x(2,1)*x(2,1)* c2 * ls1  +
C    1     (+1)*x(1,2)*x(4,2)* (-1)*x(2,2)*x(2,2)* c2 * s2l
C    1   )  - (
C    1     (-1)*x(1,1)*x(2,1)* (+1)*x(2,1)*x(4,1)* c1 * ls1  +
C    1     (-1)*x(1,1)*x(2,1)* (+1)*x(2,2)*x(4,2)* c1 * s2l  +
C    1     (-1)*x(1,2)*x(2,2)* (+1)*x(2,1)*x(4,1)* c2 * ls1  +
C    1     (-1)*x(1,2)*x(2,2)* (+1)*x(2,2)*x(4,2)* c2 * s2l
C    1   )
       tca(1,2) =
     1  (
     1     (+1)*x(1,1)*x(4,1)* (-1)*x(2,2)*x(2,2)* c1 * s2l  +
     1     (+1)*x(1,2)*x(4,2)* (-1)*x(2,1)*x(2,1)* c2 * ls1
     1   )  - (
     1     (-1)*x(1,1)*x(2,1)* (+1)*x(2,2)*x(4,2)* c1 * s2l  +
     1     (-1)*x(1,2)*x(2,2)* (+1)*x(2,1)*x(4,1)* c2 * ls1
     1   )
c-------------------------------------------------------------
C      else if(i.eq.1 .and. j.eq.3)then
C      zout =
C    1  (
C    1     (+1)*x(1,1)*x(4,1)* (+1)*x(2,1)*x(1,1)* c1 * c1   +
C    1     (+1)*x(1,1)*x(4,1)* (+1)*x(1,2)*x(2,2)* c1 * c2   +
C    1     (+1)*x(1,2)*x(4,2)* (+1)*x(2,1)*x(1,1)* c2 * c1   +
C    1     (+1)*x(1,2)*x(4,2)* (+1)*x(1,2)*x(2,2)* c2 * c2
C    1   )  - (
C    1     (+1)*x(1,1)*x(1,1)* (+1)*x(2,1)*x(4,1)* s1l* ls1  +
C    1     (+1)*x(1,1)*x(1,1)* (+1)*x(2,2)*x(4,2)* s1l* s2l  +
C    1     (+1)*x(1,2)*x(1,2)* (+1)*x(2,1)*x(4,1)* ls2* ls1  +
C    1     (+1)*x(1,2)*x(1,2)* (+1)*x(2,2)*x(4,2)* ls2* s2l
C    1   )
       tca(1,3) =
     1  (
     1     (+1)*x(1,1)*x(4,1)* (+1)*x(2,1)*x(1,1)*4.*a1*a1*dfac  +
     1     (+1)*x(1,1)*x(4,1)* (+1)*x(1,2)*x(2,2)* c1 * c2   +
     1     (+1)*x(1,2)*x(4,2)* (+1)*x(2,1)*x(1,1)* c2 * c1   +
     1     (-1)*x(1,1)*x(1,1)* (+1)*x(2,2)*x(4,2)* s1l* s2l  +
     1     (-1)*x(1,2)*x(1,2)* (+1)*x(2,1)*x(4,1)* ls2* ls1  +
     1     (+1)*x(1,2)*x(4,2)* (+1)*x(1,2)*x(2,2)*4.*a2*a2*dfac
     1   )
c-------------------------------------------------------------
C      else if(i.eq.1 .and. j.eq.4)then
C      zout =
C    1  (
C    1     (-1)*x(1,1)*x(3,1)* (-1)*x(2,1)*x(2,1)* s1l* ls1  +
C    1     (-1)*x(1,1)*x(3,1)* (-1)*x(2,2)*x(2,2)* s1l* s2l  +
C    1     (-1)*x(1,2)*x(3,2)* (-1)*x(2,1)*x(2,1)* ls2* ls1  +
C    1     (-1)*x(1,2)*x(3,2)* (-1)*x(2,2)*x(2,2)* ls2* s2l
C    1   )  - (
C    1     (-1)*x(1,1)*x(2,1)* (-1)*x(2,1)*x(3,1)* c1 * c1   +
C    1     (-1)*x(1,1)*x(2,1)* (-1)*x(2,2)*x(3,2)* c1 * c2   +
C    1     (-1)*x(1,2)*x(2,2)* (-1)*x(2,1)*x(3,1)* c2 * c1   +
C    1     (-1)*x(1,2)*x(2,2)* (-1)*x(2,2)*x(3,2)* c2 * c2
C    1   )
       tca(1,4) =
     1  (
     1     (-1)*x(1,2)*x(3,2)* (-1)*x(2,1)*x(2,1)* ls2* ls1  +
     1     (+1)*x(1,1)*x(2,1)* (-1)*x(2,1)*x(3,1)*4.*a1*a1*dfac +
     1     (-1)*x(1,1)*x(3,1)* (-1)*x(2,2)*x(2,2)* s1l* s2l  +
     1     (+1)*x(1,1)*x(2,1)* (-1)*x(2,2)*x(3,2)* c1 * c2   +
     1     (+1)*x(1,2)*x(2,2)* (-1)*x(2,1)*x(3,1)* c2 * c1   +
     1     (+1)*x(1,2)*x(2,2)* (-1)*x(2,2)*x(3,2)*4.*a2*a2*dfac
     1   )
c-------------------------------------------------------------
C      else if(i.eq.1 .and. j.eq.5)then
C      zout =
C    1  (
C    1     (-1)*x(1,1)*x(3,1)* (+1)*x(2,1)*x(1,1)* s1l* c1   +
C    1     (-1)*x(1,1)*x(3,1)* (+1)*x(1,2)*x(2,2)* s1l* c2   +
C    1     (-1)*x(1,2)*x(3,2)* (+1)*x(2,1)*x(1,1)* ls2* c1   +
C    1     (-1)*x(1,2)*x(3,2)* (+1)*x(1,2)*x(2,2)* ls2* c2
C    1   )  - (
C    1     (+1)*x(1,1)*x(1,1)* (-1)*x(2,1)*x(3,1)* s1l* c1   +
C    1     (+1)*x(1,1)*x(1,1)* (-1)*x(2,2)*x(3,2)* s1l* c2   +
C    1     (+1)*x(1,2)*x(1,2)* (-1)*x(2,1)*x(3,1)* ls2* c1   +
C    1     (+1)*x(1,2)*x(1,2)* (-1)*x(2,2)*x(3,2)* ls2* c2
C    1   )
       tca(1,5) =
     1  (
     1     (-1)*x(1,1)*x(3,1)* (+1)*x(1,2)*x(2,2)* s1l* c2   +
     1     (-1)*x(1,2)*x(3,2)* (+1)*x(2,1)*x(1,1)* ls2* c1
     1   )  - (
     1     (+1)*x(1,1)*x(1,1)* (-1)*x(2,2)*x(3,2)* s1l* c2   +
     1     (+1)*x(1,2)*x(1,2)* (-1)*x(2,1)*x(3,1)* ls2* c1
     1   )
c-------------------------------------------------------------
C      else if(i.eq.1 .and. j.eq.6)then
C      zout =
C    1  (
C    1     (-1)*x(1,1)*x(2,1)* (+1)*x(2,1)*x(1,1)* c1 * c1   +
C    1     (-1)*x(1,1)*x(2,1)* (+1)*x(1,2)*x(2,2)* c1 * c2   +
C    1     (-1)*x(1,2)*x(2,2)* (+1)*x(2,1)*x(1,1)* c2 * c1   +
C    1     (-1)*x(1,2)*x(2,2)* (+1)*x(1,2)*x(2,2)* c2 * c2
C    1   )  - (
C    1     (+1)*x(1,1)*x(1,1)* (-1)*x(2,1)*x(2,1)* s1l* ls1  +
C    1     (+1)*x(1,1)*x(1,1)* (-1)*x(2,2)*x(2,2)* s1l* s2l  +
C    1     (+1)*x(1,2)*x(1,2)* (-1)*x(2,1)*x(2,1)* ls2* ls1  +
C    1     (+1)*x(1,2)*x(1,2)* (-1)*x(2,2)*x(2,2)* ls2* s2l
C    1   )
       tca(1,6) =
     1  (
     1     (-1)*x(1,1)*x(2,1)* (+1)*x(1,2)*x(2,2)* c1 * c2   +
     1     (-1)*x(1,2)*x(2,2)* (+1)*x(2,1)*x(1,1)* c2 * c1   +
     1     (-1)*x(1,1)*x(1,1)* (-1)*x(2,2)*x(2,2)* s1l* s2l  +
     1     (-1)*x(1,2)*x(1,2)* (-1)*x(2,1)*x(2,1)* ls2* ls1  +
     1     (-1)*x(1,1)*x(2,1)* (+1)*x(2,1)*x(1,1)*4*a1*a1*dfac +
     1     (-1)*x(1,2)*x(2,2)* (+1)*x(1,2)*x(2,2)*4*a2*a2*dfac
     1   )
c-------------------------------------------------------------
C      else if(i.eq.2 .and. j.eq.1)then
C      zout =
C    1  (
C    1     (+1)*x(1,1)*x(4,1)* (-1)*x(3,1)*x(3,1)* c1 * s1l  +
C    1     (+1)*x(1,1)*x(4,1)* (-1)*x(3,2)*x(3,2)* c1 * ls2  +
C    1     (+1)*x(1,2)*x(4,2)* (-1)*x(3,1)*x(3,1)* c2 * s1l  +
C    1     (+1)*x(1,2)*x(4,2)* (-1)*x(3,2)*x(3,2)* c2 * ls2
C    1   )  - (
C    1     (-1)*x(1,1)*x(3,1)* (+1)*x(3,1)*x(4,1)* s1l* c1   +
C    1     (-1)*x(1,1)*x(3,1)* (+1)*x(3,2)*x(4,2)* s1l* c2   +
C    1     (-1)*x(1,2)*x(3,2)* (+1)*x(3,1)*x(4,1)* ls2* c1   +
C    1     (-1)*x(1,2)*x(3,2)* (+1)*x(3,2)*x(4,2)* ls2* c2
C    1   )
       tca(2,1) =
     1  (
     1     (+1)*x(1,1)*x(4,1)* (-1)*x(3,2)*x(3,2)* c1 * ls2  +
     1     (+1)*x(1,2)*x(4,2)* (-1)*x(3,1)*x(3,1)* c2 * s1l
     1   )  - (
     1     (-1)*x(1,1)*x(3,1)* (+1)*x(3,2)*x(4,2)* s1l* c2   +
     1     (-1)*x(1,2)*x(3,2)* (+1)*x(3,1)*x(4,1)* ls2* c1
     1   )
c-------------------------------------------------------------
C      else if(i.eq.2 .and. j.eq.2)then
C      zout =
C    1  (
C    1     (+1)*x(1,1)*x(4,1)* (-1)*x(3,1)*x(2,1)* c1 * c1   +
C    1     (+1)*x(1,1)*x(4,1)* (-1)*x(3,2)*x(2,2)* c1 * c2   +
C    1     (+1)*x(1,2)*x(4,2)* (-1)*x(3,1)*x(2,1)* c2 * c1   +
C    1     (+1)*x(1,2)*x(4,2)* (-1)*x(3,2)*x(2,2)* c2 * c2
C    1   )  - (
C    1     (-1)*x(1,1)*x(2,1)* (+1)*x(3,1)*x(4,1)* c1 * c1   +
C    1     (-1)*x(1,1)*x(2,1)* (+1)*x(3,2)*x(4,2)* c1 * c2   +
C    1     (-1)*x(1,2)*x(2,2)* (+1)*x(3,1)*x(4,1)* c2 * c1   +
C    1     (-1)*x(1,2)*x(2,2)* (+1)*x(3,2)*x(4,2)* c2 * c2
C    1   )
       tca(2,2) =
     1  (
     1     (+1)*x(1,1)*x(4,1)* (-1)*x(3,2)*x(2,2)* c1 * c2   +
     1     (+1)*x(1,1)*x(2,1)* (+1)*x(3,2)*x(4,2)* c1 * c2   +
     1     (+1)*x(1,2)*x(4,2)* (-1)*x(3,1)*x(2,1)* c2 * c1   +
     1     (+1)*x(1,2)*x(2,2)* (+1)*x(3,1)*x(4,1)* c2 * c1
     1   )
c-------------------------------------------------------------
C      else if(i.eq.2 .and. j.eq.3)then
C      zout =
C    1  (
C    1     (+1)*x(1,1)*x(4,1)* (+1)*x(3,1)*x(1,1)* c1 * s1l  +
C    1     (+1)*x(1,1)*x(4,1)* (+1)*x(3,2)*x(1,2)* c1 * ls2  +
C    1     (+1)*x(1,2)*x(4,2)* (+1)*x(3,1)*x(1,1)* c2 * s1l  +
C    1     (+1)*x(1,2)*x(4,2)* (+1)*x(3,2)*x(1,2)* c2 * ls2
C    1   )  - (
C    1     (+1)*x(1,1)*x(1,1)* (+1)*x(3,1)*x(4,1)* s1l* c1   +
C    1     (+1)*x(1,1)*x(1,1)* (+1)*x(3,2)*x(4,2)* s1l* c2   +
C    1     (+1)*x(1,2)*x(1,2)* (+1)*x(3,1)*x(4,1)* ls2* c1   +
C    1     (+1)*x(1,2)*x(1,2)* (+1)*x(3,2)*x(4,2)* ls2* c2
C    1   )
       tca(2,3) =
     1  (
     1     (+1)*x(1,1)*x(4,1)* (+1)*x(3,2)*x(1,2)* c1 * ls2  +
     1     (+1)*x(1,2)*x(4,2)* (+1)*x(3,1)*x(1,1)* c2 * s1l
     1   )  - (
     1     (+1)*x(1,1)*x(1,1)* (+1)*x(3,2)*x(4,2)* s1l* c2   +
     1     (+1)*x(1,2)*x(1,2)* (+1)*x(3,1)*x(4,1)* ls2* c1
     1   )
c-------------------------------------------------------------
C      else if(i.eq.2 .and. j.eq.4)then
C      zout =
C    1  (
C    1     (-1)*x(1,1)*x(3,1)* (-1)*x(3,1)*x(2,1)* s1l* c1   +
C    1     (-1)*x(1,1)*x(3,1)* (-1)*x(3,2)*x(2,2)* s1l* c2   +
C    1     (-1)*x(1,2)*x(3,2)* (-1)*x(3,1)*x(2,1)* ls2* c1   +
C    1     (-1)*x(1,2)*x(3,2)* (-1)*x(3,2)*x(2,2)* ls2* c2
C    1   )  - (
C    1     (-1)*x(1,1)*x(2,1)* (-1)*x(3,1)*x(3,1)* c1 * s1l  +
C    1     (-1)*x(1,1)*x(2,1)* (-1)*x(3,2)*x(3,2)* c1 * ls2  +
C    1     (-1)*x(1,2)*x(2,2)* (-1)*x(3,1)*x(3,1)* c2 * s1l  +
C    1     (-1)*x(1,2)*x(2,2)* (-1)*x(3,2)*x(3,2)* c2 * ls2
C    1   )
       tca(2,4) =
     1  (
     1     (-1)*x(1,1)*x(3,1)* (-1)*x(3,2)*x(2,2)* s1l* c2   +
     1     (-1)*x(1,2)*x(3,2)* (-1)*x(3,1)*x(2,1)* ls2* c1
     1   )  - (
     1     (-1)*x(1,1)*x(2,1)* (-1)*x(3,2)*x(3,2)* c1 * ls2  +
     1     (-1)*x(1,2)*x(2,2)* (-1)*x(3,1)*x(3,1)* c2 * s1l
     1   )
c-------------------------------------------------------------
C      else if(i.eq.2 .and. j.eq.5)then
C      zout =
C    1  (
C    1     (-1)*x(1,1)*x(3,1)* (+1)*x(3,1)*x(1,1)* s1l* s1l  +
C    1     (-1)*x(1,1)*x(3,1)* (+1)*x(3,2)*x(1,2)* s1l* ls2  +
C    1     (-1)*x(1,2)*x(3,2)* (+1)*x(3,1)*x(1,1)* ls2* s1l  +
C    1     (-1)*x(1,2)*x(3,2)* (+1)*x(3,2)*x(1,2)* ls2* ls2
C    1   )  - (
C    1     (+1)*x(1,1)*x(1,1)* (-1)*x(3,1)*x(3,1)* s1l* s1l  +
C    1     (+1)*x(1,1)*x(1,1)* (-1)*x(3,2)*x(3,2)* s1l* ls2  +
C    1     (+1)*x(1,2)*x(1,2)* (-1)*x(3,1)*x(3,1)* ls2* s1l  +
C    1     (+1)*x(1,2)*x(1,2)* (-1)*x(3,2)*x(3,2)* ls2* ls2
C    1   )
       tca(2,5) =
     1  (
     1     (-1)*x(1,1)*x(3,1)* (+1)*x(3,2)*x(1,2)* s1l* ls2  +
     1     (-1)*x(1,2)*x(3,2)* (+1)*x(3,1)*x(1,1)* ls2* s1l
     1   )  - (
     1     (+1)*x(1,1)*x(1,1)* (-1)*x(3,2)*x(3,2)* s1l* ls2  +
     1     (+1)*x(1,2)*x(1,2)* (-1)*x(3,1)*x(3,1)* ls2* s1l
     1   )
c-------------------------------------------------------------
C      else if(i.eq.2 .and. j.eq.6)then
C      zout =
C    1  (
C    1     (-1)*x(1,1)*x(2,1)* (+1)*x(3,1)*x(1,1)* c1 * s1l  +
C    1     (-1)*x(1,1)*x(2,1)* (+1)*x(3,2)*x(1,2)* c1 * ls2  +
C    1     (-1)*x(1,2)*x(2,2)* (+1)*x(3,1)*x(1,1)* c2 * s1l  +
C    1     (-1)*x(1,2)*x(2,2)* (+1)*x(3,2)*x(1,2)* c2 * ls2
C    1   )  - (
C    1     (+1)*x(1,1)*x(1,1)* (-1)*x(3,1)*x(2,1)* s1l* c1   +
C    1     (+1)*x(1,1)*x(1,1)* (-1)*x(3,2)*x(2,2)* s1l* c2   +
C    1     (+1)*x(1,2)*x(1,2)* (-1)*x(3,1)*x(2,1)* ls2* c1   +
C    1     (+1)*x(1,2)*x(1,2)* (-1)*x(3,2)*x(2,2)* ls2* c2
C    1   )
       tca(2,6) =
     1  (
     1     (-1)*x(1,1)*x(2,1)* (+1)*x(3,2)*x(1,2)* c1 * ls2  +
     1     (-1)*x(1,2)*x(2,2)* (+1)*x(3,1)*x(1,1)* c2 * s1l
     1   )  - (
     1     (+1)*x(1,1)*x(1,1)* (-1)*x(3,2)*x(2,2)* s1l* c2   +
     1     (+1)*x(1,2)*x(1,2)* (-1)*x(3,1)*x(2,1)* ls2* c1
     1   )
c-------------------------------------------------------------
C      else if(i.eq.3 .and. j.eq.1)then
C      zout =
C    1  (
C    1     (+1)*x(1,1)*x(4,1)* (-1)*x(4,1)*x(3,1)* c1 * c1   +
C    1     (+1)*x(1,1)*x(4,1)* (-1)*x(4,2)*x(3,2)* c1 * c2   +
C    1     (+1)*x(1,2)*x(4,2)* (-1)*x(4,1)*x(3,1)* c2 * c1   +
C    1     (+1)*x(1,2)*x(4,2)* (-1)*x(4,2)*x(3,2)* c2 * c2
C    1   )  - (
C    1     (-1)*x(1,1)*x(3,1)* (+1)*x(4,1)*x(4,1)* s1l* ls1  +
C    1     (-1)*x(1,1)*x(3,1)* (+1)*x(4,2)*x(4,2)* s1l* s2l  +
C    1     (-1)*x(1,2)*x(3,2)* (+1)*x(4,1)*x(4,1)* ls2* ls1  +
C    1     (-1)*x(1,2)*x(3,2)* (+1)*x(4,2)*x(4,2)* ls2* s2l
C    1   )
       tca(3,1) =
     1  (
     1     (+1)*x(1,1)*x(4,1)* (-1)*x(4,2)*x(3,2)* c1 * c2   +
     1     (+1)*x(1,2)*x(4,2)* (-1)*x(4,1)*x(3,1)* c2 * c1   +
     1     (+1)*x(1,2)*x(3,2)* (+1)*x(4,1)*x(4,1)* ls2* ls1  +
     1     (+1)*x(1,1)*x(3,1)* (+1)*x(4,2)*x(4,2)* s1l* s2l  +
     1     (+1)*x(1,1)*x(4,1)* (-1)*x(4,1)*x(3,1)*4*a1*a1*dfac +
     1     (+1)*x(1,2)*x(4,2)* (-1)*x(4,2)*x(3,2)*4*a2*a2*dfac
     1   )
c-------------------------------------------------------------
C      else if(i.eq.3 .and. j.eq.2)then
C      zout =
C    1  (
C    1     (+1)*x(1,1)*x(4,1)* (-1)*x(4,1)*x(2,1)* c1 * ls1  +
C    1     (+1)*x(1,1)*x(4,1)* (-1)*x(4,2)*x(2,2)* c1 * s2l  +
C    1     (+1)*x(1,2)*x(4,2)* (-1)*x(4,1)*x(2,1)* c2 * ls1  +
C    1     (+1)*x(1,2)*x(4,2)* (-1)*x(4,2)*x(2,2)* c2 * s2l
C    1   )  - (
C    1     (-1)*x(1,1)*x(2,1)* (+1)*x(4,1)*x(4,1)* c1 * ls1  +
C    1     (-1)*x(1,1)*x(2,1)* (+1)*x(4,2)*x(4,2)* c1 * s2l  +
C    1     (-1)*x(1,2)*x(2,2)* (+1)*x(4,1)*x(4,1)* c2 * ls1  +
C    1     (-1)*x(1,2)*x(2,2)* (+1)*x(4,2)*x(4,2)* c2 * s2l
C    1   )
       tca(3,2) =
     1  (
     1     (+1)*x(1,1)*x(4,1)* (-1)*x(4,2)*x(2,2)* c1 * s2l  +
     1     (+1)*x(1,2)*x(4,2)* (-1)*x(4,1)*x(2,1)* c2 * ls1
     1   )  - (
     1     (-1)*x(1,1)*x(2,1)* (+1)*x(4,2)*x(4,2)* c1 * s2l  +
     1     (-1)*x(1,2)*x(2,2)* (+1)*x(4,1)*x(4,1)* c2 * ls1
     1   )
c-------------------------------------------------------------
C      else if(i.eq.3 .and. j.eq.3)then
C      zout =
C    1  (
C    1     (+1)*x(1,1)*x(4,1)* (+1)*x(4,1)*x(1,1)* c1 * c1   +
C    1     (+1)*x(1,1)*x(4,1)* (+1)*x(4,2)*x(1,2)* c1 * c2   +
C    1     (+1)*x(1,2)*x(4,2)* (+1)*x(4,1)*x(1,1)* c2 * c1   +
C    1     (+1)*x(1,2)*x(4,2)* (+1)*x(4,2)*x(1,2)* c2 * c2
C    1   )  - (
C    1     (+1)*x(1,1)*x(1,1)* (+1)*x(4,1)*x(4,1)* s1l* ls1  +
C    1     (+1)*x(1,1)*x(1,1)* (+1)*x(4,2)*x(4,2)* s1l* s2l  +
C    1     (+1)*x(1,2)*x(1,2)* (+1)*x(4,1)*x(4,1)* ls2* ls1  +
C    1     (+1)*x(1,2)*x(1,2)* (+1)*x(4,2)*x(4,2)* ls2* s2l
C    1   )
       tca(3,3) =
     1  (
     1     (+1)*x(1,1)*x(4,1)* (+1)*x(4,2)*x(1,2)* c1 * c2   +
     1     (+1)*x(1,2)*x(4,2)* (+1)*x(4,1)*x(1,1)* c2 * c1   +
     1     (-1)*x(1,1)*x(1,1)* (+1)*x(4,2)*x(4,2)* s1l* s2l  +
     1     (-1)*x(1,2)*x(1,2)* (+1)*x(4,1)*x(4,1)* ls2* ls1  +
     1     (+1)*x(1,1)*x(4,1)* (+1)*x(4,1)*x(1,1)*4*a1*a1*dfac +
     1     (+1)*x(1,2)*x(4,2)* (+1)*x(4,2)*x(1,2)*4*a2*a2*dfac
     1   )
c-------------------------------------------------------------
C      else if(i.eq.3 .and. j.eq.4)then
C      zout =
C    1  (
C    1     (-1)*x(1,1)*x(3,1)* (-1)*x(4,1)*x(2,1)* s1l* ls1  +
C    1     (-1)*x(1,1)*x(3,1)* (-1)*x(4,2)*x(2,2)* s1l* s2l  +
C    1     (-1)*x(1,2)*x(3,2)* (-1)*x(4,1)*x(2,1)* ls2* ls1  +
C    1     (-1)*x(1,2)*x(3,2)* (-1)*x(4,2)*x(2,2)* ls2* s2l
C    1   )  - (
C    1     (-1)*x(1,1)*x(2,1)* (-1)*x(4,1)*x(3,1)* c1 * c1   +
C    1     (-1)*x(1,1)*x(2,1)* (-1)*x(4,2)*x(3,2)* c1 * c2   +
C    1     (-1)*x(1,2)*x(2,2)* (-1)*x(4,1)*x(3,1)* c2 * c1   +
C    1     (-1)*x(1,2)*x(2,2)* (-1)*x(4,2)*x(3,2)* c2 * c2
C    1   )
       tca(3,4) =
     1  (
     1     (-1)*x(1,1)*x(3,1)* (-1)*x(4,2)*x(2,2)* s1l* s2l  +
     1     (-1)*x(1,2)*x(3,2)* (-1)*x(4,1)*x(2,1)* ls2* ls1  +
     1     (+1)*x(1,1)*x(2,1)* (-1)*x(4,2)*x(3,2)* c1 * c2   +
     1     (+1)*x(1,2)*x(2,2)* (-1)*x(4,1)*x(3,1)* c2 * c1   +
     1     (+1)*x(1,1)*x(2,1)* (-1)*x(4,1)*x(3,1)*4*a1*a1*dfac +
     1     (+1)*x(1,2)*x(2,2)* (-1)*x(4,2)*x(3,2)*4*a2*a2*dfac
     1   )
c-------------------------------------------------------------
C      else if(i.eq.3 .and. j.eq.5)then
C      zout =
C    1  (
C    1     (-1)*x(1,1)*x(3,1)* (+1)*x(4,1)*x(1,1)* s1l* c1   +
C    1     (-1)*x(1,1)*x(3,1)* (+1)*x(4,2)*x(1,2)* s1l* c2   +
C    1     (-1)*x(1,2)*x(3,2)* (+1)*x(4,1)*x(1,1)* ls2* c1   +
C    1     (-1)*x(1,2)*x(3,2)* (+1)*x(4,2)*x(1,2)* ls2* c2
C    1   )  - (
C    1     (+1)*x(1,1)*x(1,1)* (-1)*x(4,1)*x(3,1)* s1l* c1   +
C    1     (+1)*x(1,1)*x(1,1)* (-1)*x(4,2)*x(3,2)* s1l* c2   +
C    1     (+1)*x(1,2)*x(1,2)* (-1)*x(4,1)*x(3,1)* ls2* c1   +
C    1     (+1)*x(1,2)*x(1,2)* (-1)*x(4,2)*x(3,2)* ls2* c2
C    1   )
       tca(3,5) =
     1  (
     1     (-1)*x(1,1)*x(3,1)* (+1)*x(4,2)*x(1,2)* s1l* c2   +
     1     (-1)*x(1,2)*x(3,2)* (+1)*x(4,1)*x(1,1)* ls2* c1
     1   )  - (
     1     (+1)*x(1,1)*x(1,1)* (-1)*x(4,2)*x(3,2)* s1l* c2   +
     1     (+1)*x(1,2)*x(1,2)* (-1)*x(4,1)*x(3,1)* ls2* c1
     1   )
c-------------------------------------------------------------
C      else if(i.eq.3 .and. j.eq.6)then
C      zout =
C    1  (
C    1     (-1)*x(1,1)*x(2,1)* (+1)*x(4,1)*x(1,1)* c1 * c1   +
C    1     (-1)*x(1,1)*x(2,1)* (+1)*x(4,2)*x(1,2)* c1 * c2   +
C    1     (-1)*x(1,2)*x(2,2)* (+1)*x(4,1)*x(1,1)* c2 * c1   +
C    1     (-1)*x(1,2)*x(2,2)* (+1)*x(4,2)*x(1,2)* c2 * c2
C    1   )  - (
C    1     (+1)*x(1,1)*x(1,1)* (-1)*x(4,1)*x(2,1)* s1l* ls1  +
C    1     (+1)*x(1,1)*x(1,1)* (-1)*x(4,2)*x(2,2)* s1l* s2l  +
C    1     (+1)*x(1,2)*x(1,2)* (-1)*x(4,1)*x(2,1)* ls2* ls1  +
C    1     (+1)*x(1,2)*x(1,2)* (-1)*x(4,2)*x(2,2)* ls2* s2l
C    1   )
       tca(3,6) =
     1  (
     1     (-1)*x(1,1)*x(2,1)* (+1)*x(4,2)*x(1,2)* c1 * c2   +
     1     (-1)*x(1,2)*x(2,2)* (+1)*x(4,1)*x(1,1)* c2 * c1   +
     1     (-1)*x(1,1)*x(1,1)* (-1)*x(4,2)*x(2,2)* s1l* s2l  +
     1     (-1)*x(1,2)*x(1,2)* (-1)*x(4,1)*x(2,1)* ls2* ls1  +
     1     (-1)*x(1,1)*x(2,1)* (+1)*x(4,1)*x(1,1)*4*a1*a1*dfac +
     1     (-1)*x(1,2)*x(2,2)* (+1)*x(4,2)*x(1,2)*4*a2*a2*dfac
     1   )
c-------------------------------------------------------------
C      else if(i.eq.4 .and. j.eq.1)then
C      zout =
C    1  (
C    1     (+1)*x(2,1)*x(4,1)* (-1)*x(3,1)*x(3,1)* ls1* s1l  +
C    1     (+1)*x(2,1)*x(4,1)* (-1)*x(3,2)*x(3,2)* ls1* ls2  +
C    1     (+1)*x(2,2)*x(4,2)* (-1)*x(3,1)*x(3,1)* s2l* s1l  +
C    1     (+1)*x(2,2)*x(4,2)* (-1)*x(3,2)*x(3,2)* s2l* ls2
C    1   )  - (
C    1     (-1)*x(2,1)*x(3,1)* (+1)*x(3,1)*x(4,1)* c1 * c1   +
C    1     (-1)*x(2,1)*x(3,1)* (+1)*x(3,2)*x(4,2)* c1 * c2   +
C    1     (-1)*x(2,2)*x(3,2)* (+1)*x(3,1)*x(4,1)* c2 * c1   +
C    1     (-1)*x(2,2)*x(3,2)* (+1)*x(3,2)*x(4,2)* c2 * c2
C    1   )
       tca(4,1) =
     1  (
     1     (+1)*x(2,1)*x(4,1)* (-1)*x(3,2)*x(3,2)* ls1* ls2  +
     1     (+1)*x(2,2)*x(4,2)* (-1)*x(3,1)*x(3,1)* s2l* s1l  +
     1     (+1)*x(2,1)*x(3,1)* (+1)*x(3,2)*x(4,2)* c1 * c2   +
     1     (+1)*x(2,2)*x(3,2)* (+1)*x(3,1)*x(4,1)* c2 * c1   +
     1     (+1)*x(2,1)*x(3,1)* (+1)*x(3,1)*x(4,1)*4*a1*a1*dfac +
     1     (+1)*x(2,2)*x(3,2)* (+1)*x(3,2)*x(4,2)*4*a2*a2*dfac
     1   )
c-------------------------------------------------------------
C      else if(i.eq.4 .and. j.eq.2)then
C      zout =
C    1  (
C    1     (+1)*x(2,1)*x(4,1)* (-1)*x(3,1)*x(2,1)* ls1* c1   +
C    1     (+1)*x(2,1)*x(4,1)* (-1)*x(3,2)*x(2,2)* ls1* c2   +
C    1     (+1)*x(2,2)*x(4,2)* (-1)*x(3,1)*x(2,1)* s2l* c1   +
C    1     (+1)*x(2,2)*x(4,2)* (-1)*x(3,2)*x(2,2)* s2l* c2
C    1   )  - (
C    1     (-1)*x(2,1)*x(2,1)* (+1)*x(3,1)*x(4,1)* ls1* c1   +
C    1     (-1)*x(2,1)*x(2,1)* (+1)*x(3,2)*x(4,2)* ls1* c2   +
C    1     (-1)*x(2,2)*x(2,2)* (+1)*x(3,1)*x(4,1)* s2l* c1   +
C    1     (-1)*x(2,2)*x(2,2)* (+1)*x(3,2)*x(4,2)* s2l* c2
C    1   )
       tca(4,2) =
     1  (
     1     (+1)*x(2,1)*x(4,1)* (-1)*x(3,2)*x(2,2)* ls1* c2   +
     1     (+1)*x(2,2)*x(4,2)* (-1)*x(3,1)*x(2,1)* s2l* c1
     1   )  - (
     1     (-1)*x(2,1)*x(2,1)* (+1)*x(3,2)*x(4,2)* ls1* c2   +
     1     (-1)*x(2,2)*x(2,2)* (+1)*x(3,1)*x(4,1)* s2l* c1
     1   )
c-------------------------------------------------------------
C      else if(i.eq.4 .and. j.eq.3)then
C      zout =
C    1  (
C    1     (+1)*x(2,1)*x(4,1)* (+1)*x(3,1)*x(1,1)* ls1* s1l  +
C    1     (+1)*x(2,1)*x(4,1)* (+1)*x(3,2)*x(1,2)* ls1* ls2  +
C    1     (+1)*x(2,2)*x(4,2)* (+1)*x(3,1)*x(1,1)* s2l* s1l  +
C    1     (+1)*x(2,2)*x(4,2)* (+1)*x(3,2)*x(1,2)* s2l* ls2
C    1   )  - (
C    1     (+1)*x(2,1)*x(1,1)* (+1)*x(3,1)*x(4,1)* c1 * c1   +
C    1     (+1)*x(2,1)*x(1,1)* (+1)*x(3,2)*x(4,2)* c1 * c2   +
C    1     (+1)*x(1,2)*x(2,2)* (+1)*x(3,1)*x(4,1)* c2 * c1   +
C    1     (+1)*x(1,2)*x(2,2)* (+1)*x(3,2)*x(4,2)* c2 * c2
C    1   )
       tca(4,3) =
     1  (
     1     (+1)*x(2,1)*x(4,1)* (+1)*x(3,2)*x(1,2)* ls1* ls2  +
     1     (+1)*x(2,2)*x(4,2)* (+1)*x(3,1)*x(1,1)* s2l* s1l  +
     1     (-1)*x(2,1)*x(1,1)* (+1)*x(3,2)*x(4,2)* c1 * c2   +
     1     (-1)*x(1,2)*x(2,2)* (+1)*x(3,1)*x(4,1)* c2 * c1   +
     1     (-1)*x(2,1)*x(1,1)* (+1)*x(3,1)*x(4,1)*4*a1*a1*dfac +
     1     (-1)*x(1,2)*x(2,2)* (+1)*x(3,2)*x(4,2)*4*a2*a2*dfac
     1   )
c-------------------------------------------------------------
C      else if(i.eq.4 .and. j.eq.4)then
C      zout =
C    1  (
C    1     (-1)*x(2,1)*x(3,1)* (-1)*x(3,1)*x(2,1)* c1 * c1   +
C    1     (-1)*x(2,1)*x(3,1)* (-1)*x(3,2)*x(2,2)* c1 * c2   +
C    1     (-1)*x(2,2)*x(3,2)* (-1)*x(3,1)*x(2,1)* c2 * c1   +
C    1     (-1)*x(2,2)*x(3,2)* (-1)*x(3,2)*x(2,2)* c2 * c2
C    1   )  - (
C    1     (-1)*x(2,1)*x(2,1)* (-1)*x(3,1)*x(3,1)* ls1* s1l  +
C    1     (-1)*x(2,1)*x(2,1)* (-1)*x(3,2)*x(3,2)* ls1* ls2  +
C    1     (-1)*x(2,2)*x(2,2)* (-1)*x(3,1)*x(3,1)* s2l* s1l  +
C    1     (-1)*x(2,2)*x(2,2)* (-1)*x(3,2)*x(3,2)* s2l* ls2
C    1   )
       tca(4,4) =
     1  (
     1     (-1)*x(2,1)*x(3,1)* (-1)*x(3,2)*x(2,2)* c1 * c2   +
     1     (-1)*x(2,2)*x(3,2)* (-1)*x(3,1)*x(2,1)* c2 * c1   +
     1     (+1)*x(2,1)*x(2,1)* (-1)*x(3,2)*x(3,2)* ls1* ls2  +
     1     (+1)*x(2,2)*x(2,2)* (-1)*x(3,1)*x(3,1)* s2l* s1l  +
     1     (-1)*x(2,1)*x(3,1)* (-1)*x(3,1)*x(2,1)*4*a1*a1*dfac +
     1     (-1)*x(2,2)*x(3,2)* (-1)*x(3,2)*x(2,2)*4*a2*a2*dfac
     1   )
c-------------------------------------------------------------
C      else if(i.eq.4 .and. j.eq.5)then
C      zout =
C    1  (
C    1     (-1)*x(2,1)*x(3,1)* (+1)*x(3,1)*x(1,1)* c1 * s1l  +
C    1     (-1)*x(2,1)*x(3,1)* (+1)*x(3,2)*x(1,2)* c1 * ls2  +
C    1     (-1)*x(2,2)*x(3,2)* (+1)*x(3,1)*x(1,1)* c2 * s1l  +
C    1     (-1)*x(2,2)*x(3,2)* (+1)*x(3,2)*x(1,2)* c2 * ls2
C    1   )  - (
C    1     (+1)*x(2,1)*x(1,1)* (-1)*x(3,1)*x(3,1)* c1 * s1l  +
C    1     (+1)*x(2,1)*x(1,1)* (-1)*x(3,2)*x(3,2)* c1 * ls2  +
C    1     (+1)*x(1,2)*x(2,2)* (-1)*x(3,1)*x(3,1)* c2 * s1l  +
C    1     (+1)*x(1,2)*x(2,2)* (-1)*x(3,2)*x(3,2)* c2 * ls2
C    1   )
       tca(4,5) =
     1  (
     1     (-1)*x(2,1)*x(3,1)* (+1)*x(3,2)*x(1,2)* c1 * ls2  +
     1     (-1)*x(2,2)*x(3,2)* (+1)*x(3,1)*x(1,1)* c2 * s1l
     1   )  - (
     1     (+1)*x(2,1)*x(1,1)* (-1)*x(3,2)*x(3,2)* c1 * ls2  +
     1     (+1)*x(1,2)*x(2,2)* (-1)*x(3,1)*x(3,1)* c2 * s1l
     1   )
c-------------------------------------------------------------
C      else if(i.eq.4 .and. j.eq.6)then
C      zout =
C    1  (
C    1     (-1)*x(2,1)*x(2,1)* (+1)*x(3,1)*x(1,1)* ls1* s1l  +
C    1     (-1)*x(2,1)*x(2,1)* (+1)*x(3,2)*x(1,2)* ls1* ls2  +
C    1     (-1)*x(2,2)*x(2,2)* (+1)*x(3,1)*x(1,1)* s2l* s1l  +
C    1     (-1)*x(2,2)*x(2,2)* (+1)*x(3,2)*x(1,2)* s2l* ls2
C    1   )  - (
C    1     (+1)*x(2,1)*x(1,1)* (-1)*x(3,1)*x(2,1)* c1 * c1   +
C    1     (+1)*x(2,1)*x(1,1)* (-1)*x(3,2)*x(2,2)* c1 * c2   +
C    1     (+1)*x(1,2)*x(2,2)* (-1)*x(3,1)*x(2,1)* c2 * c1   +
C    1     (+1)*x(1,2)*x(2,2)* (-1)*x(3,2)*x(2,2)* c2 * c2
C    1   )
       tca(4,6) =
     1  (
     1     (-1)*x(2,1)*x(2,1)* (+1)*x(3,2)*x(1,2)* ls1* ls2  +
     1     (-1)*x(2,2)*x(2,2)* (+1)*x(3,1)*x(1,1)* s2l* s1l  +
     1     (-1)*x(2,1)*x(1,1)* (-1)*x(3,2)*x(2,2)* c1 * c2   +
     1     (-1)*x(1,2)*x(2,2)* (-1)*x(3,1)*x(2,1)* c2 * c1   +
     1     (-1)*x(2,1)*x(1,1)* (-1)*x(3,1)*x(2,1)*4*a1*a1*dfac +
     1     (-1)*x(1,2)*x(2,2)* (-1)*x(3,2)*x(2,2)*4*a2*a2*dfac
     1   )
c-------------------------------------------------------------
C      else if(i.eq.5 .and. j.eq.1)then
C      zout =
C    1  (
C    1     (+1)*x(2,1)*x(4,1)* (-1)*x(4,1)*x(3,1)* ls1* c1   +
C    1     (+1)*x(2,1)*x(4,1)* (-1)*x(4,2)*x(3,2)* ls1* c2   +
C    1     (+1)*x(2,2)*x(4,2)* (-1)*x(4,1)*x(3,1)* s2l* c1   +
C    1     (+1)*x(2,2)*x(4,2)* (-1)*x(4,2)*x(3,2)* s2l* c2
C    1   )  - (
C    1     (-1)*x(2,1)*x(3,1)* (+1)*x(4,1)*x(4,1)* c1 * ls1  +
C    1     (-1)*x(2,1)*x(3,1)* (+1)*x(4,2)*x(4,2)* c1 * s2l  +
C    1     (-1)*x(2,2)*x(3,2)* (+1)*x(4,1)*x(4,1)* c2 * ls1  +
C    1     (-1)*x(2,2)*x(3,2)* (+1)*x(4,2)*x(4,2)* c2 * s2l
C    1   )
       tca(5,1) =
     1  (
     1     (+1)*x(2,1)*x(4,1)* (-1)*x(4,2)*x(3,2)* ls1* c2   +
     1     (+1)*x(2,2)*x(4,2)* (-1)*x(4,1)*x(3,1)* s2l* c1
     1   )  - (
     1     (-1)*x(2,1)*x(3,1)* (+1)*x(4,2)*x(4,2)* c1 * s2l  +
     1     (-1)*x(2,2)*x(3,2)* (+1)*x(4,1)*x(4,1)* c2 * ls1
     1   )
c-------------------------------------------------------------
C      else if(i.eq.5 .and. j.eq.2)then
C      zout =
C    1  (
C    1     (+1)*x(2,1)*x(4,1)* (-1)*x(4,1)*x(2,1)* ls1* ls1  +
C    1     (+1)*x(2,1)*x(4,1)* (-1)*x(4,2)*x(2,2)* ls1* s2l  +
C    1     (+1)*x(2,2)*x(4,2)* (-1)*x(4,1)*x(2,1)* s2l* ls1  +
C    1     (+1)*x(2,2)*x(4,2)* (-1)*x(4,2)*x(2,2)* s2l* s2l
C    1   )  - (
C    1     (-1)*x(2,1)*x(2,1)* (+1)*x(4,1)*x(4,1)* ls1* ls1  +
C    1     (-1)*x(2,1)*x(2,1)* (+1)*x(4,2)*x(4,2)* ls1* s2l  +
C    1     (-1)*x(2,2)*x(2,2)* (+1)*x(4,1)*x(4,1)* s2l* ls1  +
C    1     (-1)*x(2,2)*x(2,2)* (+1)*x(4,2)*x(4,2)* s2l* s2l
C    1   )
       tca(5,2) =
     1  (
     1     (+1)*x(2,1)*x(4,1)* (-1)*x(4,2)*x(2,2)* ls1* s2l  +
     1     (+1)*x(2,2)*x(4,2)* (-1)*x(4,1)*x(2,1)* s2l* ls1  +
     1     (+1)*x(2,1)*x(2,1)* (+1)*x(4,2)*x(4,2)* ls1* s2l  +
     1     (+1)*x(2,2)*x(2,2)* (+1)*x(4,1)*x(4,1)* s2l* ls1
     1   )
c-------------------------------------------------------------
C      else if(i.eq.5 .and. j.eq.3)then
C      zout =
C    1  (
C    1     (+1)*x(2,1)*x(4,1)* (+1)*x(4,1)*x(1,1)* ls1* c1   +
C    1     (+1)*x(2,1)*x(4,1)* (+1)*x(4,2)*x(1,2)* ls1* c2   +
C    1     (+1)*x(2,2)*x(4,2)* (+1)*x(4,1)*x(1,1)* s2l* c1   +
C    1     (+1)*x(2,2)*x(4,2)* (+1)*x(4,2)*x(1,2)* s2l* c2
C    1   )  - (
C    1     (+1)*x(2,1)*x(1,1)* (+1)*x(4,1)*x(4,1)* c1 * ls1  +
C    1     (+1)*x(2,1)*x(1,1)* (+1)*x(4,2)*x(4,2)* c1 * s2l  +
C    1     (+1)*x(1,2)*x(2,2)* (+1)*x(4,1)*x(4,1)* c2 * ls1  +
C    1     (+1)*x(1,2)*x(2,2)* (+1)*x(4,2)*x(4,2)* c2 * s2l
C    1   )
       tca(5,3) =
     1  (
     1     (+1)*x(2,1)*x(4,1)* (+1)*x(4,2)*x(1,2)* ls1* c2   +
     1     (+1)*x(2,2)*x(4,2)* (+1)*x(4,1)*x(1,1)* s2l* c1
     1   )  - (
     1     (+1)*x(2,1)*x(1,1)* (+1)*x(4,2)*x(4,2)* c1 * s2l  +
     1     (+1)*x(1,2)*x(2,2)* (+1)*x(4,1)*x(4,1)* c2 * ls1
     1   )
c-------------------------------------------------------------
C      else if(i.eq.5 .and. j.eq.4)then
C      zout =
C    1  (
C    1     (-1)*x(2,1)*x(3,1)* (-1)*x(4,1)*x(2,1)* c1 * ls1  +
C    1     (-1)*x(2,1)*x(3,1)* (-1)*x(4,2)*x(2,2)* c1 * s2l  +
C    1     (-1)*x(2,2)*x(3,2)* (-1)*x(4,1)*x(2,1)* c2 * ls1  +
C    1     (-1)*x(2,2)*x(3,2)* (-1)*x(4,2)*x(2,2)* c2 * s2l
C    1   )  - (
C    1     (-1)*x(2,1)*x(2,1)* (-1)*x(4,1)*x(3,1)* ls1* c1   +
C    1     (-1)*x(2,1)*x(2,1)* (-1)*x(4,2)*x(3,2)* ls1* c2   +
C    1     (-1)*x(2,2)*x(2,2)* (-1)*x(4,1)*x(3,1)* s2l* c1   +
C    1     (-1)*x(2,2)*x(2,2)* (-1)*x(4,2)*x(3,2)* s2l* c2
C    1   )
       tca(5,4) =
     1  (
     1     (-1)*x(2,1)*x(3,1)* (-1)*x(4,2)*x(2,2)* c1 * s2l  +
     1     (-1)*x(2,2)*x(3,2)* (-1)*x(4,1)*x(2,1)* c2 * ls1
     1   )  - (
     1     (-1)*x(2,1)*x(2,1)* (-1)*x(4,2)*x(3,2)* ls1* c2   +
     1     (-1)*x(2,2)*x(2,2)* (-1)*x(4,1)*x(3,1)* s2l* c1
     1   )
c-------------------------------------------------------------
C      else if(i.eq.5 .and. j.eq.5)then
C      zout =
C    1  (
C    1     (-1)*x(2,1)*x(3,1)* (+1)*x(4,1)*x(1,1)* c1 * c1   +
C    1     (-1)*x(2,1)*x(3,1)* (+1)*x(4,2)*x(1,2)* c1 * c2   +
C    1     (-1)*x(2,2)*x(3,2)* (+1)*x(4,1)*x(1,1)* c2 * c1   +
C    1     (-1)*x(2,2)*x(3,2)* (+1)*x(4,2)*x(1,2)* c2 * c2
C    1   )  - (
C    1     (+1)*x(2,1)*x(1,1)* (-1)*x(4,1)*x(3,1)* c1 * c1   +
C    1     (+1)*x(2,1)*x(1,1)* (-1)*x(4,2)*x(3,2)* c1 * c2   +
C    1     (+1)*x(1,2)*x(2,2)* (-1)*x(4,1)*x(3,1)* c2 * c1   +
C    1     (+1)*x(1,2)*x(2,2)* (-1)*x(4,2)*x(3,2)* c2 * c2
C    1   )
       tca(5,5) =
     1  (
     1     (-1)*x(2,1)*x(3,1)* (+1)*x(4,2)*x(1,2)* c1 * c2   +
     1     (-1)*x(2,2)*x(3,2)* (+1)*x(4,1)*x(1,1)* c2 * c1   +
     1     (-1)*x(2,1)*x(1,1)* (-1)*x(4,2)*x(3,2)* c1 * c2   +
     1     (-1)*x(1,2)*x(2,2)* (-1)*x(4,1)*x(3,1)* c2 * c1
     1   )
c-------------------------------------------------------------
C      else if(i.eq.5 .and. j.eq.6)then
C      zout =
C    1  (
C    1     (-1)*x(2,1)*x(2,1)* (+1)*x(4,1)*x(1,1)* ls1* c1   +
C    1     (-1)*x(2,1)*x(2,1)* (+1)*x(4,2)*x(1,2)* ls1* c2   +
C    1     (-1)*x(2,2)*x(2,2)* (+1)*x(4,1)*x(1,1)* s2l* c1   +
C    1     (-1)*x(2,2)*x(2,2)* (+1)*x(4,2)*x(1,2)* s2l* c2
C    1   )  - (
C    1     (+1)*x(2,1)*x(1,1)* (-1)*x(4,1)*x(2,1)* c1 * ls1  +
C    1     (+1)*x(2,1)*x(1,1)* (-1)*x(4,2)*x(2,2)* c1 * s2l  +
C    1     (+1)*x(1,2)*x(2,2)* (-1)*x(4,1)*x(2,1)* c2 * ls1  +
C    1     (+1)*x(1,2)*x(2,2)* (-1)*x(4,2)*x(2,2)* c2 * s2l
C    1   )
       tca(5,6) =
     1  (
     1     (-1)*x(2,1)*x(2,1)* (+1)*x(4,2)*x(1,2)* ls1* c2   +
     1     (-1)*x(2,2)*x(2,2)* (+1)*x(4,1)*x(1,1)* s2l* c1
     1   )  - (
     1     (+1)*x(2,1)*x(1,1)* (-1)*x(4,2)*x(2,2)* c1 * s2l  +
     1     (+1)*x(1,2)*x(2,2)* (-1)*x(4,1)*x(2,1)* c2 * ls1
     1   )
c-------------------------------------------------------------
C      else if(i.eq.6 .and. j.eq.1)then
C      zout =
C    1  (
C    1     (+1)*x(3,1)*x(4,1)* (-1)*x(4,1)*x(3,1)* c1 * c1   +
C    1     (+1)*x(3,1)*x(4,1)* (-1)*x(4,2)*x(3,2)* c1 * c2   +
C    1     (+1)*x(3,2)*x(4,2)* (-1)*x(4,1)*x(3,1)* c2 * c1   +
C    1     (+1)*x(3,2)*x(4,2)* (-1)*x(4,2)*x(3,2)* c2 * c2
C    1   )  - (
C    1     (-1)*x(3,1)*x(3,1)* (+1)*x(4,1)*x(4,1)* s1l* ls1  +
C    1     (-1)*x(3,1)*x(3,1)* (+1)*x(4,2)*x(4,2)* s1l* s2l  +
C    1     (-1)*x(3,2)*x(3,2)* (+1)*x(4,1)*x(4,1)* ls2* ls1  +
C    1     (-1)*x(3,2)*x(3,2)* (+1)*x(4,2)*x(4,2)* ls2* s2l
C    1   )
       tca(6,1) =
     1  (
     1     (+1)*x(3,1)*x(4,1)* (-1)*x(4,2)*x(3,2)* c1 * c2   +
     1     (+1)*x(3,2)*x(4,2)* (-1)*x(4,1)*x(3,1)* c2 * c1   +
     1     (+1)*x(3,1)*x(3,1)* (+1)*x(4,2)*x(4,2)* s1l* s2l  +
     1     (+1)*x(3,2)*x(3,2)* (+1)*x(4,1)*x(4,1)* ls2* ls1  +
     1     (+1)*x(3,1)*x(4,1)* (-1)*x(4,1)*x(3,1)*4*a1*a1*dfac +
     1     (+1)*x(3,2)*x(4,2)* (-1)*x(4,2)*x(3,2)*4*a2*a2*dfac
     1   )
c-------------------------------------------------------------
C      else if(i.eq.6 .and. j.eq.2)then
C      zout =
C    1  (
C    1     (+1)*x(3,1)*x(4,1)* (-1)*x(4,1)*x(2,1)* c1 * ls1  +
C    1     (+1)*x(3,1)*x(4,1)* (-1)*x(4,2)*x(2,2)* c1 * s2l  +
C    1     (+1)*x(3,2)*x(4,2)* (-1)*x(4,1)*x(2,1)* c2 * ls1  +
C    1     (+1)*x(3,2)*x(4,2)* (-1)*x(4,2)*x(2,2)* c2 * s2l
C    1   )  - (
C    1     (-1)*x(3,1)*x(2,1)* (+1)*x(4,1)*x(4,1)* c1 * ls1  +
C    1     (-1)*x(3,1)*x(2,1)* (+1)*x(4,2)*x(4,2)* c1 * s2l  +
C    1     (-1)*x(3,2)*x(2,2)* (+1)*x(4,1)*x(4,1)* c2 * ls1  +
C    1     (-1)*x(3,2)*x(2,2)* (+1)*x(4,2)*x(4,2)* c2 * s2l
C    1   )
       tca(6,2) =
     1  (
     1     (+1)*x(3,1)*x(4,1)* (-1)*x(4,2)*x(2,2)* c1 * s2l  +
     1     (+1)*x(3,2)*x(4,2)* (-1)*x(4,1)*x(2,1)* c2 * ls1
     1   )  - (
     1     (-1)*x(3,1)*x(2,1)* (+1)*x(4,2)*x(4,2)* c1 * s2l  +
     1     (-1)*x(3,2)*x(2,2)* (+1)*x(4,1)*x(4,1)* c2 * ls1
     1   )
c-------------------------------------------------------------
c      else if(i.eq.6 .and. j.eq.3)then
C      zout =
C    1  (
C    1     (+1)*x(3,1)*x(4,1)* (+1)*x(4,1)*x(1,1)* c1 * c1   +
C    1     (+1)*x(3,1)*x(4,1)* (+1)*x(4,2)*x(1,2)* c1 * c2   +
C    1     (+1)*x(3,2)*x(4,2)* (+1)*x(4,1)*x(1,1)* c2 * c1   +
C    1     (+1)*x(3,2)*x(4,2)* (+1)*x(4,2)*x(1,2)* c2 * c2
C    1   )  - (
C    1     (+1)*x(3,1)*x(1,1)* (+1)*x(4,1)*x(4,1)* s1l* ls1  +
C    1     (+1)*x(3,1)*x(1,1)* (+1)*x(4,2)*x(4,2)* s1l* s2l  +
C    1     (+1)*x(3,2)*x(1,2)* (+1)*x(4,1)*x(4,1)* ls2* ls1  +
C    1     (+1)*x(3,2)*x(1,2)* (+1)*x(4,2)*x(4,2)* ls2* s2l
C    1   )
       tca(6,3) =
     1  (
     1     (+1)*x(3,1)*x(4,1)* (+1)*x(4,2)*x(1,2)* c1 * c2   +
     1     (+1)*x(3,2)*x(4,2)* (+1)*x(4,1)*x(1,1)* c2 * c1   +
     1     (-1)*x(3,1)*x(1,1)* (+1)*x(4,2)*x(4,2)* s1l* s2l  +
     1     (-1)*x(3,2)*x(1,2)* (+1)*x(4,1)*x(4,1)* ls2* ls1  +
     1     (+1)*x(3,1)*x(4,1)* (+1)*x(4,1)*x(1,1)*4*a1*a1*dfac +
     1     (+1)*x(3,2)*x(1,2)* (+1)*x(4,2)*x(4,2)*4*a2*a2*dfac
     1   )
c-------------------------------------------------------------
C      else if(i.eq.6 .and. j.eq.4)then
C      zout =
C    1  (
C    1     (-1)*x(3,1)*x(3,1)* (-1)*x(4,1)*x(2,1)* s1l* ls1  +
C    1     (-1)*x(3,1)*x(3,1)* (-1)*x(4,2)*x(2,2)* s1l* s2l  +
C    1     (-1)*x(3,2)*x(3,2)* (-1)*x(4,1)*x(2,1)* ls2* ls1  +
C    1     (-1)*x(3,2)*x(3,2)* (-1)*x(4,2)*x(2,2)* ls2* s2l
C    1   )  - (
C    1     (-1)*x(3,1)*x(2,1)* (-1)*x(4,1)*x(3,1)* c1 * c1   +
C    1     (-1)*x(3,1)*x(2,1)* (-1)*x(4,2)*x(3,2)* c1 * c2   +
C    1     (-1)*x(3,2)*x(2,2)* (-1)*x(4,1)*x(3,1)* c2 * c1   +
C    1     (-1)*x(3,2)*x(2,2)* (-1)*x(4,2)*x(3,2)* c2 * c2
C    1   )
       tca(6,4) =
     1  (
     1     (-1)*x(3,1)*x(3,1)* (-1)*x(4,2)*x(2,2)* s1l* s2l  +
     1     (-1)*x(3,2)*x(3,2)* (-1)*x(4,1)*x(2,1)* ls2* ls1  +
     1     (+1)*x(3,1)*x(2,1)* (-1)*x(4,2)*x(3,2)* c1 * c2   +
     1     (+1)*x(3,2)*x(2,2)* (-1)*x(4,1)*x(3,1)* c2 * c1   +
     1     (+1)*x(3,1)*x(2,1)* (-1)*x(4,1)*x(3,1)*4*a1*a1*dfac +
     1     (+1)*x(3,2)*x(3,2)* (-1)*x(4,2)*x(2,2)*4*a2*a2*dfac
     1   )
c-------------------------------------------------------------
C      else if(i.eq.6 .and. j.eq.5)then
C      zout =
C    1  (
C    1     (-1)*x(3,1)*x(3,1)* (+1)*x(4,1)*x(1,1)* s1l* c1   +
C    1     (-1)*x(3,1)*x(3,1)* (+1)*x(4,2)*x(1,2)* s1l* c2   +
C    1     (-1)*x(3,2)*x(3,2)* (+1)*x(4,1)*x(1,1)* ls2* c1   +
C    1     (-1)*x(3,2)*x(3,2)* (+1)*x(4,2)*x(1,2)* ls2* c2
C    1   )  - (
C    1     (+1)*x(3,1)*x(1,1)* (-1)*x(4,1)*x(3,1)* s1l* c1   +
C    1     (+1)*x(3,1)*x(1,1)* (-1)*x(4,2)*x(3,2)* s1l* c2   +
C    1     (+1)*x(3,2)*x(1,2)* (-1)*x(4,1)*x(3,1)* ls2* c1   +
C    1     (+1)*x(3,2)*x(1,2)* (-1)*x(4,2)*x(3,2)* ls2* c2
C    1   )
       tca(6,5) =
     1  (
     1     (-1)*x(3,1)*x(3,1)* (+1)*x(4,2)*x(1,2)* s1l* c2   +
     1     (-1)*x(3,2)*x(3,2)* (+1)*x(4,1)*x(1,1)* ls2* c1
     1   )  - (
     1     (+1)*x(3,1)*x(1,1)* (-1)*x(4,2)*x(3,2)* s1l* c2   +
     1     (+1)*x(3,2)*x(1,2)* (-1)*x(4,1)*x(3,1)* ls2* c1
     1   )
c-------------------------------------------------------------
C      else if(i.eq.6 .and. j.eq.6)then
C      zout =
C    1  (
C    1     (-1)*x(3,1)*x(2,1)* (+1)*x(4,1)*x(1,1)* c1 * c1   +
C    1     (-1)*x(3,1)*x(2,1)* (+1)*x(4,2)*x(1,2)* c1 * c2   +
C    1     (-1)*x(3,2)*x(2,2)* (+1)*x(4,1)*x(1,1)* c2 * c1   +
C    1     (-1)*x(3,2)*x(2,2)* (+1)*x(4,2)*x(1,2)* c2 * c2
C    1   )  - (
C    1     (+1)*x(3,1)*x(1,1)* (-1)*x(4,1)*x(2,1)* s1l* ls1  +
C    1     (+1)*x(3,1)*x(1,1)* (-1)*x(4,2)*x(2,2)* s1l* s2l  +
C    1     (+1)*x(3,2)*x(1,2)* (-1)*x(4,1)*x(2,1)* ls2* ls1  +
C    1     (+1)*x(3,2)*x(1,2)* (-1)*x(4,2)*x(2,2)* ls2* s2l
C    1   )
       tca(6,6) =
     1  (
     1     (-1)*x(3,1)*x(2,1)* (+1)*x(4,2)*x(1,2)* c1 * c2   +
     1     (-1)*x(3,2)*x(2,2)* (+1)*x(4,1)*x(1,1)* c2 * c1   +
     1     (-1)*x(3,2)*x(1,2)* (-1)*x(4,1)*x(2,1)* ls2* ls1  +
     1     (-1)*x(3,1)*x(1,1)* (-1)*x(4,2)*x(2,2)* s1l* s2l  +
     1     (-1)*x(3,1)*x(2,1)* (+1)*x(4,1)*x(1,1)*4*a1*a1*dfac +
     1     (-1)*x(3,2)*x(2,2)* (+1)*x(4,2)*x(1,2)*4*a2*a2*dfac
     1   )
C       endif
c-------------------------------------------------------------
C       do i=1,6
C          do j=1,6
C            WRITE(6,*)'tca(',i,',',j,')=',i,j,tca(i,j)
C          enddo
C       enddo
c-----
c       note theoreticalla CA34 = 1 - CA33
c       however we factor out the term exp(ex) for numerical
c       reason.  Thus the expression here
c       is CA34 = dfac - CA33
c-----
        tca(1,4) = - tca(1,3)
        tca(2,4) = - tca(2,3)
        tca(2,6) =   tca(1,5)
        tca(3,4) = dfac - tca(3,3)
        tca(3,5) = - tca(2,3)
        tca(3,6) = - tca(1,3)
        tca(4,1) = - tca(3,1)
        tca(4,2) = - tca(3,2)
        tca(4,3) = dfac - tca(3,3)
        tca(4,4) =   tca(3,3)
        tca(4,5) =   tca(2,3)
        tca(4,6) =   tca(1,3)
        tca(5,3) = - tca(3,2)
        tca(5,4) =   tca(3,2)
        tca(5,5) =   tca(2,2)
        tca(5,6) =   tca(1,2)
        tca(6,2) =   tca(5,1)
        tca(6,3) = - tca(3,1)
        tca(6,4) =   tca(3,1)
        tca(6,5) =   tca(2,1)
        tca(6,6) =   tca(1,1)
C       do i=1,6
C          do j=1,6
C             CA(i,j) = tca(i,j)
C            WRITE(6,*)'CA (',i,',',j,')=',i,j,CA(i,j),CA(i,j)/tca(i,j)
C          enddo
C       enddo
        CA(1,1) = tca(1,1)
        CA(1,2) = tca(1,2)
        CA(1,3) = tca(1,3)
        CA(1,4) = tca(1,5)
        CA(1,5) = tca(1,6)
        CA(2,1) = tca(2,1)
        CA(2,2) = tca(2,2)
        CA(2,3) = tca(2,3)
        CA(2,4) = tca(2,5)
        CA(2,5) = tca(1,5)
        CA(3,1) = 2*tca(3,1)
        CA(3,2) = 2*tca(3,2)
        CA(3,3) = 2*tca(3,3) -dfac
        CA(3,4) = -2*tca(2,3)
        CA(3,5) = -2*tca(1,3)
        CA(4,1) = tca(5,1)
        CA(4,2) = tca(5,2)
        CA(4,3) = - tca(3,2)
        CA(4,4) = tca(2,2)
        CA(4,5) = tca(1,2)
        CA(5,1) = tca(6,1)
        CA(5,2) = tca(5,1)
        CA(5,3) = - tca(3,1)
        CA(5,4) = tca(2,1)
        CA(5,5) = tca(1,1)
        endif
        return
        end

        subroutine evalg(jbdry,m,m1,gbr,indx,
     1      wvno,om,om2,wvno2)
        implicit none
        integer jbdry, m, m1, indx
        complex*16 gbr(2,5)
        real*8 wvno,om,wvno2,om2
        integer NL
        parameter(NL=200)
        common/timod/d(NL),TA(NL),TC(NL),TF(NL),TL(NL),TN(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL),
     3      frefp(NL), frefs(NL)
        real*8 d,TA,TC,TN,TL,TF,TRho,qa,
     1       qb,etap,etas,frefp,frefs
        common/modlly/mmax
        integer mmax
        common/modspec/allfluid
        logical allfluid

        common/emat/e,einv,ra,rb
        complex*16 ra,rb
        complex*16 e(4,4), einv(4,4)


        complex*16 cg(6)
        complex*16 g(4,4)
        COMPLEX*16 zA,zB,zC,zD,zE,zF
        COMPLEX*16 rp, rsv
        COMPLEX*16 NP, NSV
        COMPLEX*16 X11, X21, X31, X41
        COMPLEX*16 X12, X22, X32, X42
        integer iwat

c       complex*16 e(4,4), einv(4,4)

c       integer i,j,k
c       complex*16 zsum


C        WRITE(6,*)'allfluid:',allfluid
c-----
c       set up halfspace conditions
c-----
            if(TL(m).eq. 0.0 .or.TN(m).eq.0.0)then
                iwat = 1
            else
                iwat = 0
            endif
c-----
c       HALFSPACE
c-----
            call gettiegn(za,zb,zc,zd,ze,zf,om2,wvno2,rp, rsv,
     1          x11,x21,x31,x41,x12,x22,x32,x42,NP, NSV,m,om,wvno,
     2          dcmplx(1.0d+00,0.0d+00),dcmplx(1.0d+00,0.0d+00))
            ra = rp
            rb = rsv
C       WRITE(6,*)'evalg:ra,rb:',ra,rb
            if(iwat.eq.0)then
c-----
c               ELASTIC HALFSPACE
c       This is E^-1
c-----
        G(1,1) =  x41*rp /( 2.*rp*NP)
        G(2,1) =  x42    /( 2.*rsv*NSV)
        G(3,1) = -x41*rp /(-2.*rp*NP)
        G(4,1) =  x42    /(-2.*rsv*NSV)
        G(1,2) = -x31    /( 2.*rp*NP)
        G(2,2) = -x32*rsv/( 2.*rsv*NSV)
        G(3,2) = -x31    /(-2.*rp*NP)
        G(4,2) =  x32*rsv/(-2.*rsv*NSV)
        G(1,3) = -x21*rp /( 2.*rp*NP)
        G(2,3) = -x22    /( 2.*rsv*NSV)
        G(3,3) =  x21*rp /(-2.*rp*NP)
        G(4,3) = -x22    /(-2.*rsv*NSV)
        G(1,4) =  x11    /( 2.*rp*NP)
        G(2,4) =  x12*rsv/( 2.*rsv*NSV)
        G(3,4) =  x11    /(-2.*rp*NP)
        G(4,4) = -x12*rsv/(-2.*rsv*NSV)
c-----
c       This is E
c-----
           E(1,1) =  x11
           E(2,1) =  x21 * rp
           E(3,1) =  x31
           E(4,1) =  x41 * rp

           E(1,2) =  x12 * rsv
           E(2,2) =  x22
           E(3,2) =  x32 * rsv
           E(4,2) =  x42

           E(1,3) =  x11
           E(2,3) = -x21 * rp
           E(3,3) =  x31
           E(4,3) = -x41 * rp

           E(1,4) = -x12 * rsv
           E(2,4) =  x22
           E(3,4) = -x32 * rsv
           E(4,4) =  x42

        EINV(1,1) = G(1,1)
        EINV(2,1) = G(2,1)
        EINV(3,1) = G(3,1)
        EINV(4,1) = G(4,1)
        EINV(1,2) = G(1,2)
        EINV(2,2) = G(2,2)
        EINV(3,2) = G(3,2)
        EINV(4,2) = G(4,2)
        EINV(1,3) = G(1,3)
        EINV(2,3) = G(2,3)
        EINV(3,3) = G(3,3)
        EINV(4,3) = G(4,3)
        EINV(1,4) = G(1,4)
        EINV(2,4) = G(2,4)
        EINV(3,4) = G(3,4)
        EINV(4,4) = G(4,4)
c
C          do j=1,4
C              do i = 1,4
C                   write(6,*)'E   (',I ,',' , J, ')=',E(i,J)
C             enddo
C          enddo
C          do j=1,4
C              do i = 1,4
C                   write(6,*)'EINV(',I ,',' , J, ')=',G(i,J)
C             enddo
C          enddo
C          do i=1,4
C              do j = 1,4
C                zsum = dcmplx(0.0d+00,0.0d+00)
C                do k=1,4
C                zsum = zsum + E(i,k)*g(k,j)
C                enddo
C                write(6,*)'E INV(',I ,',' , J, ')=',ZSUM
C             enddo
C          enddo

c CG(1) = G 12 12 = 11 22 - 12 21
c CG(2) = G 12 13 = 11 23 - 13 21
c CG(3) = G 12 14 = 11 24 - 14 21
c CG(4) = G 12 23 = 12 23 - 13 22
c CG(5) = G 12 24 = 12 24 - 14 22
c CG(6) = G 12 34 = 13 24 - 14 23
                CG(1) = G(1,1)*G(2,2) - G(1,2)*G(2,1)
                CG(2) = G(1,1)*G(2,3) - G(1,3)*G(2,1)
                CG(3) = G(1,1)*G(2,4) - G(1,4)*G(2,1)
                CG(4) = G(1,2)*G(2,3) - G(1,3)*G(2,2)
                CG(5) = G(1,2)*G(2,4) - G(1,4)*G(2,2)
                CG(6) = G(1,3)*G(2,4) - G(1,4)*G(2,3)
                gbr(indx,1) = CG(1)
                gbr(indx,2) = CG(2)
                gbr(indx,3) = CG(3)
                gbr(indx,4) = CG(5)
                gbr(indx,5) = CG(6)

            else if(iwat.eq.1)then
c-----
c               FLUID HALFSPACE
c-----
                if(allfluid)then
                    gbr(indx,1) = dble(TRho(m))*om2
                    gbr(indx,2) = -rp
                    gbr(indx,3) = dcmplx(0.0d+00,0.0d+00)
                    gbr(indx,4) = dcmplx(0.0d+00,0.0d+00)
                    gbr(indx,5) = dcmplx(0.0d+00,0.0d+00)
                else
                    gbr(indx,1) = dcmplx(0.0d+00,0.0d+00)
                    gbr(indx,2) = dcmplx(0.0d+00,0.0d+00)
                    gbr(indx,3) = dcmplx(0.0d+00,0.0d+00)
                    gbr(indx,4) = -dble(TRho(m))*om2
                    gbr(indx,5) = rp
                endif
            endif
        return
        end

        subroutine gettiegn(a,b,c,d,e,f,omega2,wvno2,rp, rsv,
     1      x11,x21,x31,x41,x12,x22,x32,x42,NP, NSV,m,omg,wvn,
     2      atna,atnb)
        implicit none
        COMPLEX*16 A,B,C,D,E,F
        real*8 wvno2, omega2
        COMPLEX*16 rp, rsv
c-----
c       get the modified eigen vectors x11 and x31 have
c           lambda1 (rp ) factored out
c               modified eigen vectors x22 and x42 have
c           lambda2 (rsv) factored out
c-----
        COMPLEX*16 X11, X21, X31, X41
        COMPLEX*16 X12, X22, X32, X42
c-----
c       norms
c-----
        COMPLEX*16 NP, NSV
        integer m
        real*8 omg, wvn
        complex*16 atna, atnb

        COMPLEX*16 xka2, xkb2

        integer NL
        parameter (NL=200)
        common/timod/td(NL),TA(NL),TC(NL),TF(NL),TL(NL),TN(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL),
     3      frefp(NL), frefs(NL)
        real*8 td,TA,TC,TN,TL,TF,TRho,qa,
     1       qb,etap,etas,frefp,frefs
c-----
c       internal variables
c-----
        COMPLEX*16 L2(2)
        COMPLEX*16 bb, cc
        COMPLEX*16 CDSQRT
        COMPLEX*16 SRTBB4AC
        COMPLEX*16 ddef, aabc

        COMPLEX*16 ZFAC

c-----
c       first test to see if a fluid layer - if it is fluid, the
c       eigenfunctions are specially computed and we need only the
c       rp
c-----
C       WRITE(6,*)'TI:',m,TA(m),TC(m),TL(m),TN(m),TRHO(m)
        if(TL(m).eq.0.0 .or. TN(m).eq.0.0)then
            rp = cdsqrt(dcmplx(wvno2 -omega2*TRho(m)/(TA(m)), 0.0d+00))
            rsv = dcmplx(0.0d+000, 0.0d+00)
            return
        endif


        a = wvn * TF(m) / (TC(m))
        b = 1.0/(TC(m))
        c = - TRho(m)*omg*omg + wvn*wvn *
     1      (TA(m) -TF(m)*TF(m)/(TC(m)))
        d = - wvn
        e = 1.0/(TL(m))
        f = - TRho(m)*omg*omg

c-----
c       do algebra first to avoid numerical problems
c-----
        ddef = wvn*wvn - TRho(m)*omg*omg/(TL(m))
        aabc = wvn*wvn*TA(m)/TC(m) - TRho(m)*omg*omg/(TC(m))

c-----
c       Do the QUASI P and SV - WE MUST BE CAREFUL HERE CONCERNING
c       BRANCH CUTS OF THE SQUARE ROOT
c-----
c-----
c       The characteristic equation to be solved is
c
c       L^4 - L^2[ 2 ad +ec +fb ] + [ (d^2+ef)(a^2+bc)] = 0
c-----
        bb = 2.0d+00 * a*d + e*c +f*b
        cc = ddef * aabc
c----
c       ensure that the SQRT(bb*bb - 4.0D+00*cc) is in the
c       I and II quadrants
c-----

        SRTBB4AC = CDSQRT(bb*bb - 4.0D+00*cc)
        IF(DIMAG(SRTBB4AC) .lt.0.0D+00)THEN
            SRTBB4AC = - SRTBB4AC
        ENDIF
c-----
c       recall solutions of ax^2 + bx + c =0 are
c            a(x-x1)(x-x2)
c       thus if we know x1, x2 = c/x1
c
c-----
c       Determine L^2 with care for roundoff
c-----

        IF(DREAL(BB) .LT.0.0D+00 .AND. DREAL(SRTBB4AC).LT.0.0D+00)THEN
            L2(2) = ( bb - SRTBB4AC) / 2.0d+00
            if(cdabs(L2(2)).gt.0.0)then
                 L2(1) = cc/L2(2)
            else
                 L2(1) = ( bb + SRTBB4AC) / 2.0d+00
            endif
        ELSE
            L2(1) = ( bb + SRTBB4AC) / 2.0d+00
            if(cdabs(L2(1)).gt.0.0)then
                 L2(2) = cc/L2(1)
            else
                 L2(2) = ( bb - SRTBB4AC) / 2.0d+00
            endif
        ENDIF
c-----
c       Use the Lambda^2 values to form
c       xka^2 == k^2 - L(1)^2
c       xkb^2 == k^2 - L(2)^2
c       Associate the smallest xka, xkb with the P!
c-----
        xka2 = wvno2 - L2(1)
        xkb2 = wvno2 - L2(2)
        if(cdabs(xkb2) .lt. cdabs(xka2))THEN
                ZFAC = L2(1)
                L2(1) = L2(2)
                L2(2) = ZFAC
        endif
        rp  = CDSQRT(L2(1))
        rsv = CDSQRT(L2(2))
c-----
c        safety check to so that Real part > 0
c-----
        if(dreal(rp).lt.0.0d+00)then
               rp = - rp
        endif
        if(dreal(rsv).lt.0.0d+00)then
               rsv = - rsv
        endif

c-----
c       get the norms - note that the true norm will be
c            2  NP amd 2 L(2) NSV
c       The factorization permits us to use the sin nz/n or n sin nz
c-----
C        NP  = (  L2(1)*(-2*a*b*d + 2*a*a*e + b*c*e - b*b*f)
C     1      + (a*a+b*c)*(2*b*d*d - 2*a*d*e + b*e*f - c*e*e) )
C        NSV = (- L2(2)*(2*b*d*d - 2*a*d*e - c*e*e + b*e*f)
C     1      + (d*d+e*f)*(2*a*b*d - 2*a*a*e + b*b*f - b*c*e) )
c-----
c       get the modified eigen vectors x11 and x31 have
c           lambda1 (rp ) factored out
c               modified eigen vectors x22 and x42 have
c           lambda2 (rsv) factored out
c-----
        x12 =              (b*d - a*e)
        x22 =  b*L2(2) - e*(b*c + a*a)
        x32 =    L2(2) -   (a*d + c*e)
        x42 = -a*L2(2) + d*(b*c + a*a)

        x11 = -e*L2(1) + b*(d*d + e*f)
        x21 =             ( b*d - a*e)
        x31 =  d*L2(1) - a*(d*d + e*f)
        x41 =  - ( L2(1) -  a*d - b*f)

c-----
c       TEST
c       Force the eigenfunctions to be as given in 7.4.4
c       note this will not work if wnv = 0
c-----
         if(wvn.ne.0.0)then
               zfac = wvn / x11
               x11  = x11 *zfac
               x21  = x21 *zfac
               x31  = x31 *zfac
               x41  = x41 *zfac

               zfac = wvn / x22
               x12  = x12 * zfac
               x22  = x22 * zfac
               x32  = x32 * zfac
               x42  = x42 * zfac
         endif
c-----
c       REDEFINE HERE USING THE adjusted eigenvectors
c       Note that TRUE NP  = 2 * np  * rp
c       Note that TRUE NSV = 2 * nsv * rsv     , where
c       but also note that since the 11 31 22 42 13 33 24 44 have
c       an eigenvalue factored out the normalization using these
c       elements will not use the rp and rsv
c-----
        np   = x11*x41 - x21*x31
        nsv  = x12*x42 - x22*x32

        return
        end

        subroutine varsv(p,q, rp, rsv,
     1      cosp, cosq, rsinp, rsinq,
     1      sinpr, sinqr, pex,svex,iwat,dm)
c-----
c       p = rp  * h
c       q = rsv * h
c       rp  vertical wave number for P
c       rsv vertical wave number for SV
c       cosp=cosh(p)  rsinp =rp *sinh(p)  sinpr = sinh(p)/rp
c       cosq=cosh(q)  rsinsv=rsv*sinh(p)  sinpq = sinh(p)/rsv
c         The sin rd/r = d sin(rd)/(rd)
c              so check the size o rd
c              if(cdabs(p) .lt.1.0e-4)
c                       sinpr = dm
c              else
c                       sinpr = sinh(p)/rp
c       cosq=cosh(q)  rsinsv=rsv*sinh(p)
c              if(cdabs(q) .lt.1.0e-4)
c                       sinqr = dm
c              else
c                       sinqr = sinh(q)/rsv
c-----
        implicit none
        COMPLEX*16 p, q
        COMPLEX*16 rp, rsv
        complex*16 cosp, cosq
        complex*16 rsinp, rsinq
        complex*16 sinpr, sinqr
        REAL *8 pex,svex
        integer iwat
        real*8 dm

        REAL*8 pr, pi, qr, qi
        COMPLEX*16 epp, epm, eqp, eqm
        COMPLEX*16 sinp, sinq

        REAL*8 PFAC, SVFAC

        pex  = 0.0d+00
        svex = 0.0d+00
        pr = dreal(p)
        pi = dimag(p)
        qr = dreal(q)
        qi = dimag(q)
        pex   = pr
        if(iwat.eq.1)then
c-----
c       fluid layer
c-----
            epp = dcmplx(dcos(pi), dsin(pi))/2.0
            epm = dconjg(epp)
            if(pr.lt.15.) then
                pfac=dexp(-2.*pr)
            else
                pfac  = 0.0d+00
            endif
            cosp = epp + pfac*epm
            sinp = epp - pfac*epm
            rsinp = rp *sinp
            if(dabs(pr) .lt. 1.0e-5 .and. cdabs(rp).lt.1.0e-5)then
                 sinpr = dm
            else
                 sinpr = (sinp/rp)
            endif
            cosq  = 1.0d+00
            rsinq = 0.0d+00
            sinqr = 0.0d+00
        else
c-----
c       elastic layer
c-----
            svex = qr
            epp = dcmplx(dcos(pi), dsin(pi))/2.0
            epm = dconjg(epp)
            eqp = dcmplx(dcos(qi), dsin(qi))/2.0
            eqm = dconjg(eqp)
            if(pr.lt.15.) then
                pfac=dexp(-2.*pr)
            else
                pfac  = 0.0d+00
            endif
            cosp = (epp + pfac*epm)
            sinp = epp - pfac*epm
            rsinp = (rp *sinp)
            if(dabs(pr) .lt. 1.0e-5 .and. cdabs(rp).lt.1.0e-5)then
                 sinpr = dm
            else
                 sinpr = (sinp/rp)
            endif
C           COSP  =COSP*DEXP(PEX)
C           SINPR=SINPR*DEXP(PEX)
C           RSINP=RSINP*DEXP(PEX)

            if(qr.lt.15.) then
                svfac=dexp(-2.*qr)
            else
                svfac  = 0.0d+00
            endif
            cosq = (eqp + svfac*eqm)
            sinq = eqp - svfac*eqm
            rsinq = (rsv*sinq)
            if(dabs(qr) .lt. 1.0e-5 .and. cdabs(rsv).lt.1.0e-5)then
                 sinqr = dm
            else
                 sinqr = (sinq/rsv)
            endif
C           COSQ =COSQ*DEXP(SVEX)
C           SINQR=SINQR*DEXP(SVEX)
C           RSINQ=RSINQ*DEXP(SVEX)

        endif
        return
        end

        subroutine hska(AA,tcosp,trsinp,tsinpr,tcossv,trsinsv,tsinsvr,
     1      NP,NSV,
     1      X11, X21, X31, X41,X12, X22, X32, X42,
     2      TRho,iwat,ex,om2)
c-----
c       Changes
c
c       01 May 2002  - defined cosp = tcosp/NP to
c             reduce nmber of complex divisions
c-----
        implicit none
        real*8 AA(4,4)
        complex*16 tcosp , tcossv
        complex*16 trsinp, trsinsv
        complex*16 tsinpr, tsinsvr
        COMPLEX*16 NP, NSV
        COMPLEX*16 X11, X21, X31, X41
        COMPLEX*16 X12, X22, X32, X42
        real TRho
        integer iwat
        real*8 ex, dfac
        real*8  om2
        complex*16 zrho

c-----
c       introduce shorthand to reduce work
c-----
        COMPLEX*16 cosp , cossv
        COMPLEX*16 rsinp, rsinsv
        COMPLEX*16 sinpr, sinsvr

        integer i, j
        zrho = dcmplx(dble(TRho),0.0d+00)
        if(iwat.eq.1)then
c-----
c       fluid layer
c-----
            do 100 j=1,4
                do 101 i=1,4
                    AA(i,j) = 0.0d+00
  101           continue
  100       continue
            if(ex.gt.35.0d+00)then
                dfac = 0.0d+00
            else
                dfac = dexp(-ex)
            endif
            AA(1,1) = dfac
            AA(4,4) = dfac
            AA(2,2) = tcosp
            AA(3,3) = tcosp
            AA(2,3) = -trsinp/(zrho*om2)
            AA(3,2) = - zrho*om2*tsinpr
        else
c-----
c       elastic layer
c-----
            cosp   = tcosp/NP
            sinpr  = tsinpr/NP
            rsinp  = trsinp/NP
            cossv  = tcossv/NSV
            sinsvr = tsinsvr/NSV
            rsinsv = trsinsv/NSV

            AA(1,1) = dreal(  x11*x41*cosp  + x12*x42*cossv  )
            AA(1,2) = dreal(- x11*x31*sinpr - x12*x32*rsinsv )
            AA(1,3) = dreal(- x11*x21*cosp  - x12*x22*cossv  )
            AA(1,4) = dreal(  x11*x11*sinpr + x12*x12*rsinsv )
            AA(2,1) = dreal(  x21*x41*rsinp + x22*x42*sinsvr )
            AA(2,2) = dreal(- x21*x31*cosp  - x22*x32*cossv  )
            AA(2,3) = dreal(- x21*x21*rsinp - x22*x22*sinsvr )
            AA(3,1) = dreal(  x31*x41*cosp  + x32*x42*cossv  )
            AA(3,2) = dreal(- x31*x31*sinpr - x32*x32*rsinsv )
            AA(4,1) = dreal(  x41*x41*rsinp + x42*x42*sinsvr )
            AA(2,4) = - AA(1,3)
            AA(3,3) =   AA(2,2)
            AA(3,4) = - AA(1,2)
            AA(4,2) = - AA(3,1)
            AA(4,3) = - AA(2,1)
            AA(4,4) =   AA(1,1)
        endif
        return
        end

        subroutine down(omega,wvno)
c-----
c       This finds the values of the Haskell vectors at
c       each layer boundaries from top layer downward.
c-----
        implicit none
c-----
c     command arguments
c-----
        real*8 omega, wvno
c-----
c     common blocks
c-----
        integer NL
        parameter (NL=200)
        common/timod/ tzd(NL),zta(NL),ztc(NL),ztf(NL),
     1      ztl(NL),ztn(NL),zrho(NL),zqai(NL),zqbi(NL),
     2      zetap(NL),zetas(NL),zfrefp(NL),zfrefs(NL)
        real*8 tzd, zta, ztc, ztf, ztl, ztn, zrho,
     1      zqai, zqbi, zetap, zetas,
     1      zfrefp, zfrefs
        common/pari/ mmax
        integer mmax
        common/wateri/iwat(NL)
        integer iwat
        common/eigfun/ ur(NL),uz(NL),tz(NL),tr(NL),uu0(4),
     *    dcdah(NL),dcdav(NL),dcdbh(NL),dcdbv(NL),
     2    dcdn(NL),dcdr(NL),dcdh(NL)
        real*8 ur, uz, tz, tr, uu0,
     1     dcdah, dcdav, dcdbh, dcdbv, dcdn, dcdr, dcdh
        common/dunk/   cd(NL,5),exe(NL),exa(NL)
        complex*16 cd
        real*8 exe, exa
        common/hask/   vv(NL,4)
        real*8 vv
c-----
c       model characterization
c-----
        common/modspec/allfluid
        logical allfluid
c-----
c      internal variables
c-----
       real*8 wvno2, om2
       integer m
       real*8 aa(4,4)
       integer i,j
       real*8 pex, svex
       real*8 ex2
       complex*16 rp, rsv, p, q
       complex*16 cosp,rsinp,sinpr,cossv,rsinsv,sinsvr
       real*8 cc, aa0(4)
       real*8 exsum
       complex*16 NP, NSV
        COMPLEX*16 X11, X21, X31, X41
        COMPLEX*16 X12, X22, X32, X42
        real*8 dfac, cpex
        complex*16 Za,Zb,Zc,Zd,Ze,Zf

       om2 = omega*omega
       wvno2 = wvno*wvno
c-----
c      propagate down to get the first column of the
c      Am .. A2 A1
c-----
c-----
c      initialize the top surface for the
c      first column of the Haskell propagator
c-----
      do i=1,4
         if(i.eq.1)then
             vv(1,i) = 1.0d+00
         else
             vv(1,i) = 0.0d+00
         endif
      enddo
      exa(1) = 0.0

c------
c       matrix multiplication from top layer downward
c------
       exsum  = 0.0
       do 500 m= 1, mmax -1
           call  gettiegn(Za,Zb,Zc,Zd,Ze,Zf,om2,wvno2,rp, rsv,
     1      x11,x21,x31,x41,x12,x22,x32,x42,NP, NSV,m,omega,wvno,
     2      dcmplx(1.0d+00,0.0d+00),dcmplx(1.0d+00,0.0d+00))
                p = rp  * tzd(m)
                q = rsv * tzd(m)
                call varsv(p,q,rp, rsv,
     1              cosp, cossv, rsinp, rsinsv,
     1              sinpr, sinsvr, pex,svex,iwat(m),tzd(m))
C       WRITE(6,*)'down:rp,rsv,pex,svex:',rp,rsv,pex,svex
c-----
c       For isotropic, Re p > Re sv do adjsut only the
c           SV part of the Haskell matrix
c       For TI this is not always true, so carefully adjust
c-----
            if(pex .gt. svex)then
c-----
c               PEX > SVEX, factor out PEX
c-----
                if((pex-svex).gt. 40.0d+00)then
                    dfac = 0.0d+00
                else
                    dfac = dexp(-(pex-svex))
                endif
                cpex = pex
                call hska(AA,cosp,rsinp,sinpr,
     1              dfac*cossv,dfac*rsinsv,dfac*sinsvr,NP,NSV,
     1              X11, X21, X31, X41,X12, X22, X32, X42,
     2              sngl(Zrho(m)),iwat(m),pex,om2)
        else
c-----
c               SVEX > PEX, factor out SVEX
c-----
                if((svex-pex).gt. 40.0d+00)then
                    dfac = 0.0d+00
                else
                    dfac = dexp(-(svex-pex))
                endif
                cpex = svex
                call hska(AA,dfac*cosp,dfac*rsinp,dfac*sinpr,
     1              cossv,rsinsv,sinsvr,NP,NSV,
     1              X11, X21, X31, X41,X12, X22, X32, X42,
     2              sngl(Zrho(m)),iwat(m),pex,om2)
        endif

            do 300 i=1,4
                cc=0.0d+00
                do 200 j=1,4
                    cc=cc+aa(i,j)*vv(m,j)
  200           continue
                aa0(i)=cc
  300       continue
            ex2 = 0.0
            call rnormc(aa0,ex2,4)
            exsum = exsum + cpex + ex2
            exa(m+1)= exsum

            do 400 i=1,4
                vv(m+1,i)=aa0(i)
  400       continue
C           WRITE(6,*)'m,aa0:',m,aa0
  500   continue
c-----
c       vv is  column 1 of the Haskell propagator at the top of layer m
c-----
        return
        end

        subroutine gammap(omega,wvno,gammar)
c-----
c       This routine finds the attenuation gamma value.
c
c-----
        implicit none
c-----
c       routine arguments
c-----
        real*8 omega, wvno, gammar
c-----
c       common blocks
c-----
        integer NL
        parameter (NL=200)
        common/timod/  zd(NL),zta(NL),ztc(NL),ztf(NL),
     1      ztl(NL),ztn(NL),zrho(NL),zqai(NL),zqbi(NL),
     2      zetap(NL),zetas(NL),zfrefp(NL),zfrefs(NL)
        real*8 zd, zta, ztc, ztf, ztl, ztn, zrho,
     1      zqai, zqbi, zetap, zetas,
     1      zfrefp, zfrefs
        common/pari/ mmax
        integer mmax
        common/wateri/iwat(NL)
        integer iwat
        common/eigfun/ ur(NL),uz(NL),tz(NL),tr(NL),uu0(4),
     *    dcdah(NL),dcdav(NL),dcdbh(NL),dcdbv(NL),
     2    dcdn(NL),dcdr(NL),dcdh(NL)
        real*8 ur, uz, tz, tr, uu0,
     1     dcdah, dcdav, dcdbh, dcdbv, dcdn, dcdr, dcdh
c-----
c       internal variables
c-----
        real*8 dc, pi, x, omgref, c
        integer i
        real*8 a12, a14, a21, a23
        real*8 ah,av,bh,bv,eta,rho
        real*8 TA,TC,TF,TL,TN


        gammar=0.0
        dc = 0.0
        pi = 3.141592653589493d+00
        do 100 i=1,mmax
           call getmat(i,wvno,omega,a12, a14, a21, a23,
     1      ah,av,bh,bv,eta,rho,
     2      TA,TC,TF,TL,TN,iwat(i))
            if(iwat(i).eq.  0)then
c WHOA
                x=dcdbh(i)*bh*zqbi(i) + dcdbv(i)*bv*zqbi(i)
                gammar = gammar + x
                omgref=2.0*pi*zfrefs(i)
                dc = dc + dlog(omega/omgref)*x/pi
            endif
            x=dcdav(i)*av*zqai(i) + dcdah(i)*ah*zqai(i)
            gammar = gammar + x
            omgref=2.0*pi*zfrefp(i)
            dc = dc + dlog(omega/omgref)*x/pi
  100   continue
        c=omega/wvno
        gammar=0.5*wvno*gammar/c
        c=c+dc
        wvno=omega/c
        return
        end

        subroutine energy(om,wvno,mmax)
c-----
c       determine the energy integrals and also the
c       phase velocity partials
c-----
        implicit none
c-----
c       subroutine arguments
c-----
        real*8 om, wvno
        integer mmax
c-----
c       common blocks
c-----
        integer NL
        parameter (NL=200)
        common/sumi/   sumi0,sumi1,sumi2,sumi3,flagr,are,ugr
        real*8 sumi0, sumi1, sumi2, sumi3, flagr, are, ugr
        common/eigfun/ ur(NL),uz(NL),tz(NL),tr(NL),uu0(4),
     *    dcdah(NL),dcdav(NL),dcdbh(NL),dcdbv(NL),
     2    dcdn(NL),dcdr(NL),dcdh(NL)
        real*8 ur, uz, tz, tr, uu0,
     1     dcdah, dcdav, dcdbh, dcdbv, dcdn, dcdr, dcdh
        common/wateri/iwat(NL)
        integer iwat

c-----
c       external function call prototypes
c-----
        real*8 intijr
c-----
c       internal variables
c-----
        real*8 om2, wvno2
        real*8 TA, TC, TF, TL, TN
        real*8 INT11, INT13, INT22, INT44, INT33, INT24
        real*8 URUR, UZUZ, DURDUR, DUZDUZ, URDUZ, UZDUR
        real*8 c
        real*8 fac
        real*8 facah, facav, facbh, facbv,  facr
        real*8 facn

        integer m
c-----
c       coefficients of the ODE
c-----
        real*8 a12, a14, a21, a23
        real *8 ah, av, bh, bv, eta, rho

        integer TYPELYR
c-----
c       TYPELYR = -1  layer represents upper halfspace
c                  0  layer is true intenal layer
c                 +1  layer represents lower halfspace
c-----

c----
c       initialize
c----

        sumi0 = 0.0
        sumi1 = 0.0
        sumi2 = 0.0
        sumi3 = 0.0
        c = om/wvno
        om2 = om*om
        wvno2 = wvno*wvno

        do m=1,mmax
           call getmat(m,wvno,om,a12, a14, a21, a23,
     1      ah,av,bh,bv,eta,rho,
     2      TA,TC,TF,TL,TN,iwat(m))
C          WRITE(6,*)m,wvno,om,a12, a14, a21, a23,
C    1      ah,av,bh,bv,eta,rho
c------
c       get the integrals for layers over a halfspace
c       this is here is we ever adopt the code to the coal
c       seam problem
c-----
           if(m.eq.mmax)then
                   typelyr = 1
           else
                   typelyr = 0
           endif
           INT11 = intijr(1,1,m,typelyr,om,om2,wvno,wvno2)
           INT13 = intijr(1,3,m,typelyr,om,om2,wvno,wvno2)
           INT22 = intijr(2,2,m,typelyr,om,om2,wvno,wvno2)
           INT24 = intijr(2,4,m,typelyr,om,om2,wvno,wvno2)
           INT33 = intijr(3,3,m,typelyr,om,om2,wvno,wvno2)
           INT44 = intijr(4,4,m,typelyr,om,om2,wvno,wvno2)
C          WRITE(6,*)'m:',m
C          WRITE(6,*)'INT11:',INT11
C          WRITE(6,*)'INT13:',INT13
C          WRITE(6,*)'INT22:',INT22
C          WRITE(6,*)'INT24:',INT24
C          WRITE(6,*)'INT33:',INT33
C          WRITE(6,*)'INT44:',INT44


           if(iwat(m).eq.1)then
c-----
c       fluid - not these are for the 4x4 formulation
c-----
             URUR   = INT22*(wvno/(rho*om2))**2
             UZUZ   = INT11
             URDUZ  =  - (wvno/(rho*om2))*a12*INT22
             DUZDUZ = a12*a12*INT22
             sumi0  = sumi0 + rho*(URUR + UZUZ)
             sumi1  = sumi1 + TA*URUR
             sumi2  = sumi2 - TF*URDUZ
             sumi3  = sumi3 + TC*DUZDUZ
             facah = rho*ah*(URUR -2.*eta*URDUZ/wvno)
             facav = rho*av*DUZDUZ / wvno2
             dcdah(m) = facah
             dcdav(m) = facav
             facr = -0.5*c*c*(URUR + UZUZ)
             dcdr(m) = 0.5*(av*facav + ah*facah )/rho + facr
C            WRITE(6,*)'URUR  :',URUR
C            WRITE(6,*)'UZUZ  :',UZUZ
C            WRITE(6,*)'URDUZ :',URDUZ
C            WRITE(6,*)'DUZDUZ:',DUZDUZ
c-----
c            define the ur in the fliud from the Tz - this will be at
c            the top of the layer
c-----
             ur(m) = - wvno*tz(m)/(rho*om2)

           else
c-----
c       solid
c-----
             URUR   = INT11
             UZUZ   = INT22
             DURDUR = a12*a12*INT22 + 2.*a12*a14*INT24 + a14*a14*INT44
             DUZDUZ = a21*a21*INT11 + 2.*a21*a23*INT13 + a23*a23*INT33
             URDUZ  = a21*INT11 + a23*INT13
             UZDUR  = a12*INT22 + a14*INT24
             sumi0  = sumi0 + rho*(URUR + UZUZ)
             sumi1  = sumi1 + TL*UZUZ + TA*URUR
             sumi2  = sumi2 + TL*UZDUR - TF*URDUZ
             sumi3  = sumi3 + TL*DURDUR + TC*DUZDUZ
C            WRITE(6,*)'URUR  :',URUR
C            WRITE(6,*)'UZUZ  :',UZUZ
C            WRITE(6,*)'URDUZ :',URDUZ
C            WRITE(6,*)'UZDUR :',UZDUR
C            WRITE(6,*)'DUZDUZ:',DUZDUZ
C            WRITE(6,*)'DURDUR:',DURDUR

c-----
c            partial derivatives of phase velocity with
c            respect to medium parameters. Note that these are
c            later divided by (U sumi0) when these are finalized
c
c            note that the code distinguishes between
c            ah and av, bh and bv and uses eta. These are
c            actually TI parameters, and for isotropic media
c            ah = av, bh = bv aned eta = 1. The code is written this
c            way for easier conversion to a TI case
c-----
             facah = rho*ah*(URUR -2.*eta*URDUZ/wvno)
             facav = rho*av*DUZDUZ / wvno2
             facbh = 0.0
             facbv = rho*bv*(UZUZ + 2.*UZDUR/wvno + DURDUR/wvno2 +
     1           4.*eta*URDUZ/wvno )

             facn = - TF*URDUZ/(wvno*eta)
             dcdah(m) = facah
             dcdav(m) = facav
             dcdbh(m) = facbh
             dcdbv(m) = facbv
             dcdn(m)  = facn
C       WRITE(6,*)'ah,av,bh,bv,eta:',ah,av,bh,bv,eta
C       WRITE(6,*)'facah:',facah
C       WRITE(6,*)'facav:',facav
C       WRITE(6,*)'facbh:',facbh
C       WRITE(6,*)'facbv:',facbv
C       WRITE(6,*)'m,dcda(m),dcdb(m):',m,dcda(m),dcdb(m)
C       WRITE(6,*)'m,dcdah(m),dcdav(m):',m,dcdah(m),dcdav(m)

             facr = -0.5*c*c*(URUR + UZUZ)
c-----
c      this is correct for TI
             dcdr(m) = 0.5*(av*facav + ah*facah + bv*facbv)/rho + facr

c-----
c       partial with layer thickness needs the value at the layer, e.g.,
c       dcdh(1) is top surface
c-----
           endif
        enddo
c-----
c       determine final parameters
c-----

        flagr=om2*sumi0-wvno2*sumi1-2.d+00*wvno*sumi2-sumi3
        ugr=(wvno*sumi1+sumi2)/(om*sumi0)
        are=wvno/(2.d+00*om*ugr*sumi0)
        fac = are*c/wvno2
C       WRITE(6,*)'FLAGR,UGR,ARE:',FLAGR,UGR,ARE
c-----
c       use the final factors to apply the 1/(U I0) factor
c----
        do m=1,mmax
           dcdah(m) = dcdah(m) /(ugr*sumi0)
           dcdav(m) = dcdav(m) /(ugr*sumi0)
           dcdbh(m) = dcdbh(m) /(ugr*sumi0)
           dcdbv(m) = dcdbv(m) /(ugr*sumi0)
           dcdr(m)  = dcdr(m)  /(ugr*sumi0)
           dcdn(m)  = dcdn(m)  /(ugr*sumi0)
C       WRITE(6,*)m,dcdah(m),dcdav(m),dcdbh(m),dcdbv(m),dcdr(m),
C    1      dcdn(m)
        enddo
c-----
c       determine the dcdh
c-----
        call getdcdh(om2,wvno,wvno2,fac)
        return
        end

        subroutine getmat(m,wvno,om,a12, a14, a21, a23,
     1         ah,av,bh,bv,eta,rho,
     2         TA,TC,TF,TL,TN,iwat)

c-----
c       get matrix elements of the ODE
c       We do not require all eight non-zero elements
c       also return the model characterization parameters
c-----
        implicit none
c-----
c       procedure arguments
c-----
        real*8 wvno,om,a12, a14, a21, a23
        real*8 ah,av,bh,bv,eta,rho
        real*8 TA,TC,TF,TL,TN
        integer m, iwat
c-----
c       common blocks
c-----

        integer NL
        parameter (NL=200)
        common/timod/  zd(NL),zta(NL),ztc(NL),ztf(NL),
     1      ztl(NL),ztn(NL),zrho(NL),zqai(NL),zqbi(NL),
     2      zetap(NL),zetas(NL),zfrefp(NL),zfrefs(NL)
        real*8 zd, zta, ztc, ztf, ztl, ztn, zrho,
     1      zqai, zqbi, zetap, zetas,
     1      zfrefp, zfrefs
        common/pari/ mmax
        integer mmax

        if(iwat.eq.1)then
c-----
c            non-gravitating fluid - using 2x2 ODE
c-----
             ah= sqrt(ZTA(m)/zrho(m))
             av= sqrt(ZTC(m)/zrho(m))
             bh= 0.0
             bv= 0.0
             rho=zrho(m)
             eta = 1.0

             TL = 0.0
             TN = 0.0
             TC = ZTC(m)
             TA = ZTA(m)
             TF = ZTF(m)

             a12 = - ( wvno*wvno - om*om/(ah*ah))/(rho*om*om)
        else
c-----
c            elastic - using 4x4 ODE
c-----
             ah= sqrt(ZTA(m)/zrho(m))
             av= sqrt(ZTC(m)/zrho(m))
             bh= sqrt(ZTN(m)/zrho(m))
             bv= sqrt(ZTL(m)/zrho(m))
             rho=zrho(m)


             TL = ZTL(m)
             TN = ZTN(m)
             TC = ZTC(m)
             TA = ZTA(m)
             TF = ZTF(m)
             eta = TF/(TA-2.*TL)

             a12 = -wvno
             a14 = 1.0/TL
             a21 = wvno * TF/TC
             a23 = 1.0/TC
        endif

        return
        end

        function intijr(i,j,m,typelyr,om,om2,wvno,wvno2)
        implicit none
c-----
c       procedure arguments
c-----
        real*8 intijr
        integer i,j,m
        real*8 om, om2, wvno,wvno2
        integer TYPELYR
c-----
c       TYPELYR = -1  layer represents upper halfspace
c                  0  layer is true internal layer
c                 +1  layer represents lower halfspace
c-----
c NOTE DO NOT PERMIT RA RB to be ZERO add a small number TAKEN CARE OF
c IN FUNC CALLS
c beware of E matrix for fluid - internally I use a 4x4
c but hte book uses a 2x2 in the 8-10 problem
c need typelyr, wvno, om, om2, wvno2, m, mmax, do not need medium
c-----
c       Potential coefficients
c       kmpu, kmsu are the upward coefficients at bottom of layer
c       km1pd, km1sd are the downward coefficients at top of layer
c-----
        integer NL
        parameter (NL=200)
        common/timod/  zd(NL),zta(NL),ztc(NL),ztf(NL),
     1      ztl(NL),ztn(NL),zrho(NL),zqai(NL),zqbi(NL),
     2      zetap(NL),zetas(NL),zfrefp(NL),zfrefs(NL)
        real*8 zd, zta, ztc, ztf, ztl, ztn, zrho,
     1      zqai, zqbi, zetap, zetas,
     1      zfrefp, zfrefs
        common/pari/ mmax
        integer mmax
        common/wateri/iwat(NL)
        integer iwat
        common/eigfun/ ur(NL),uz(NL),tz(NL),tr(NL),uu0(4),
     *    dcdah(NL),dcdav(NL),dcdbh(NL),dcdbv(NL),
     2    dcdn(NL),dcdr(NL),dcdh(NL)
        real*8 ur, uz, tz, tr, uu0,
     1     dcdah, dcdav, dcdbh, dcdbv, dcdn, dcdr, dcdh
        common/emat/e,einv,ra,rb
        complex*16 e(4,4), einv(4,4)
        complex*16 ra, rb

c-----
c       external function definitions
c-----
        complex*16 ffunc, gfunc, h1func, h2func
c-----
c       internal variables
c-----
        complex*16 cintijr, FA, GA, FB, GB, H1, H2
        complex*16 gbr(2,5)
        complex*16 kmpu, kmsu, km1pd, km1sd


c-----
c       call evalg to get the E and EINV matrices
c       evalg knows about water - we ignore everything else
c-----
        call evalg(0,m,m-1,gbr,1,
     1      wvno,om,om2,wvno2)
c-----
c       for an elastic solid
c
c                       T        PU  SvU  PD  SvD T
c       [Ur, Uz, Tr, Tz]  = E [ K  ,K   ,K  ,K   ]
c and
c          PU  SvU  PD  SvD T   -1                 T
c       [ K  ,K   ,K  ,K   ] = E   [Ur, Uz, Tr, Tz]
c
c       The text used the notation K m-1 to represent the potentials
c       at the top of a layer and K m those at the bottom of the layer
c       Because of Fortran indexing UZ(1) is the Z displacement at the
c       top of layer 1 and UZ(2) is the value at the bottom.  This is
c       reason for the slight incondistency in presentation below
c-----



        if(iwat(m).eq.1)then
c-----
c       fluid
c-----
            if(typelyr .lt. 0)then
              kmpu = einv(1,1)*uz(m+1)
     1             + einv(1,2)*tz(m+1)
              cintijr = e(i,1)*e(j,1)*kmpu*kmpu/(2.0*ra)
            else if(typelyr .eq. 0)then
              km1pd = einv(2,1)*uz(m)
     1              + einv(2,2)*tz(m)
              kmpu = einv(1,1)*uz(m+1)
     1             + einv(1,2)*tz(m+1)
              FA=ffunc(ra,zd(m))
              GA=gfunc(ra,zd(m))
              cintijr = e(i,1)*e(j,1)*kmpu*kmpu*FA
     1          +(e(i,1)*e(j,2)+e(i,2)*e(j,1))*kmpu*km1pd*GA
     2          + e(i,2)*e(j,2)*km1pd*km1pd*FA
            else
              km1pd = einv(2,1)*uz(m)
     2              + einv(2,2)*tz(m)
              cintijr = e(i,2)*e(j,2)*km1pd*km1pd/(2.0*ra)
            endif
        else
c-----
c       solid
c-----
            if(typelyr .lt. 0)then
c-----
c             potential coefficients for upper halfspace
c             based on the displacement, stress at top
c             Basically m = 1
c-----
              kmpu = einv(1,1)*ur(m) + einv(1,2)*uz(m)
     2             + einv(1,3)*tz(m) + einv(1,4)*tr(m)
              kmsu = einv(2,1)*ur(m) + einv(2,2)*uz(m)
     2             + einv(2,3)*tz(m) + einv(2,4)*tr(m)
              cintijr = e(i,1)*e(j,1)*kmpu*kmpu/(2.0*ra)
     1         + (e(i,1)*e(j,2)+e(i,2)*e(j,1))*kmpu*kmsu/(ra+rb)
     2         + e(i,2)*e(j,2)*kmsu*kmsu/(2.0*rb)
            else if(typelyr .eq. 0)then
c-----
c             downward potentials coefficients at top of layer m
c-----
              km1pd = einv(3,1)*ur(m) + einv(3,2)*uz(m)
     1              + einv(3,3)*tz(m) + einv(3,4)*tr(m)
              km1sd = einv(4,1)*ur(m) + einv(4,2)*uz(m)
     1              + einv(4,3)*tz(m) + einv(4,4)*tr(m)
c-----
c             upward potentials coefficients at bottom of layer m
c-----
              kmpu = einv(1,1)*ur(m+1) + einv(1,2)*uz(m+1)
     1             + einv(1,3)*tz(m+1) + einv(1,4)*tr(m+1)
              kmsu = einv(2,1)*ur(m+1) + einv(2,2)*uz(m+1)
     1             + einv(2,3)*tz(m+1) + einv(2,4)*tr(m+1)
              FA=ffunc(ra,zd(m))
              GA=gfunc(ra,zd(m))
              FB=ffunc(rb,zd(m))
              GB=gfunc(rb,zd(m))
              H1=h1func(ra,rb,zd(m))
              H2=h2func(ra,rb,zd(m))
              cintijr = e(i,1)*e(j,1)*kmpu*kmpu*FA
     1          + e(i,3)*e(j,3)*km1pd*km1pd*FA
     1          + e(i,2)*e(j,2)*kmsu*kmsu*FB
     1          + e(i,4)*e(j,4)*km1sd*km1sd*FB
     4          + H1*((e(i,1)*e(j,2)+e(i,2)*e(j,1))*kmpu*kmsu +
     5             (e(i,3)*e(j,4)+e(i,4)*e(j,3))*km1pd*km1sd)
     6          + H2*((e(i,1)*e(j,4)+e(i,4)*e(j,1))*kmpu*km1sd +
     7             (e(i,2)*e(j,3)+e(i,3)*e(j,2))*km1pd*kmsu)
     8          + GA*(e(i,1)*e(j,3)+e(i,3)*e(j,1))*kmpu*km1pd
     9          + GB*(e(i,2)*e(j,4)+e(i,4)*e(j,2))*kmsu*km1sd

            else
c-----
c             downward potential coefficients for lower halfspace
c             based on the displacement, stress at bottom
c             Basically m = mmax
c-----
              km1pd = einv(3,1)*ur(m)+einv(3,2)*uz(m)
     2            +einv(3,3)*tz(m)+einv(3,4)*tr(m)
              km1sd = einv(4,1)*ur(m)+einv(4,2)*uz(m)
     2            +einv(4,3)*tz(m)+einv(4,4)*tr(m)
              cintijr = e(i,3)*e(j,3)*km1pd*km1pd/(2.0*ra)
     1         +(e(i,3)*e(j,4)+e(i,4)*e(j,3))*km1pd*km1sd/(ra+rb)
     2         +e(i,4)*e(j,4)*km1sd*km1sd/(2.0*rb)

            endif
        endif

        intijr = dreal(cintijr)
        return
        end

        subroutine getdcdh(om2,wvno,wvno2,fac)
        implicit none
c-----
c       procedure parameters
c-----
        real*8 om2,wvno,wvno2,fac
c-----
c       common blocks
c-----
        integer NL
        parameter (NL=200)
        common/timod/  zd(NL),zta(NL),ztc(NL),ztf(NL),
     1      ztl(NL),ztn(NL),zrho(NL),zqai(NL),zqbi(NL),
     2      zetap(NL),zetas(NL),zfrefp(NL),zfrefs(NL)
        real*8 zd, zta, ztc, ztf, ztl, ztn, zrho,
     1      zqai, zqbi, zetap, zetas,
     1      zfrefp, zfrefs
        common/pari/ mmax
        integer mmax
        common/wateri/iwat(NL)
        integer iwat
        common/eigfun/ ur(NL),uz(NL),tz(NL),tr(NL),uu0(4),
     *    dcdah(NL),dcdav(NL),dcdbh(NL),dcdbv(NL),
     2    dcdn(NL),dcdr(NL),dcdh(NL)
        real*8 ur, uz, tz, tr, uu0,
     1     dcdah, dcdav, dcdbh, dcdbv, dcdn, dcdr, dcdh
c-----
c       internal variables
c-----
        real*8 tur, tuz, ttz, ttr
        real*8 dfac, gfac1, gfac2, gfac3, gfac4, gfac5, gfac6
c  DEVELOPMENT
        real*8 da, dc, dl, drho

        real*8 duzdzp, daur2
        real*8 duzdzm, drur2
        real*8 URB, DURDZM, DURDZP


        integer m
c-----

        do  m=1,mmax

            tuz = uz(m)
            ttz = tz(m)
            ttr = tr(m)
            if(iwat(m).eq.1)then
                tur = -wvno*ttz/(zrho(m)*om2)
            else
                tur = ur(m)
            endif
c-----
c       this assumes that the top is a halfspace
c-----
            if(m.eq.1)then
                drho = zrho(1) - 0.0
                da  = zta(1) - 0.0
                dc  = ztc(1) - 0.0
                dl  = ztl(1) - 0.0
                duzdzp = ( tz(m) + wvno*ztf(m)*tur)/ztc(m)
                if(iwat(m) .eq.1)then
                    durdzp = wvno*tuz
                else
                    durdzp = -wvno*tuz + ttr/ztl(m)
                endif
                drur2 = tur*tur*drho
                daur2 = tur*tur*da

                gfac1 =  om2*drho*tuz**2
                gfac2 =  om2*drur2
                gfac3 = -wvno2*dl*tuz**2
                gfac4 = -wvno2*daur2
                gfac5 =  (ztc(m)*duzdzp**2)
                gfac6 =  (ztl(m)*durdzp**2 )
            else
                drho = zrho(m) - zrho(m-1)
                da  = zta(m) - zta(m-1)
                dc  = ztc(m) - ztc(m-1)
                dl  = ztl(m) - ztl(m-1)
c---MAS
                duzdzp = (ttz + wvno*ztf(m)  *tur)/ztc(m)
                if(iwat(m).eq.1)then
                    durdzp = wvno*tuz
                else
                    durdzp = -wvno*tuz + ttr/ztl(m)
                endif
                if(iwat(m-1).eq.1 )then
                    durdzm = wvno*tuz
                else
                    durdzm = -wvno*tuz + ttr/ztl(m-1)
                endif
c-----
c       attempt to fix for water layer, since Ur is not continuous
c       across fluid - solid interface or fluid-fluid interface
c-----
                if(iwat(m-1).eq.1 .and. iwat(m).eq.0 )then
                    URB = -wvno*tz(m)/(zrho(m-1)*om2)
                    drur2 = tur*tur*zrho(m)-URB*URB*zrho(m-1)
                    daur2 = tur*tur*zta(m) -URB*URB*zta(m-1)
                    duzdzm = (ttz + wvno*ztf(m-1)*URB)/ztc(m-1)
                else if(iwat(m-1).eq.1 .and. iwat(m).eq.1 )then
                    URB = -wvno*tz(m)/(zrho(m-1)*om2)
                    drur2 = tur*tur*zrho(m)- URB*URB*zrho(m-1)
                    daur2 = tur*tur*zta(m)  - URB*URB*zta(m-1)
                    duzdzm = (ttz + wvno*ztf(m-1)*URB)/ztc(m-1)
                else
                    drur2 = tur*tur*drho
                    daur2 = tur*tur*da
                    duzdzm = (ttz + wvno*ztf(m-1)*ur(m))/ztc(m-1)
                endif


c---NO MAS
                gfac1 =  om2*drho*tuz**2
                gfac2 =  om2*drur2
                gfac3 = -wvno2*dl*tuz**2
                gfac4 = -wvno2*daur2
                gfac5 =  (ztc(m)*duzdzp**2 - ztc(m-1)*duzdzm**2)
                gfac6 =  (ztl(m)*durdzp**2 - ztl(m-1)*durdzm**2)
            endif
            dfac = fac * (
     1          gfac1 + gfac2 + gfac3 + gfac4
     2          + gfac5 + gfac6 )
            if(dabs(dfac).lt.1.0d-38)then
                dfac = 0.0d+00
            endif
            dcdh(m) = dfac
        enddo
        return
        end


!------------------------------------------------------------------------------!
! Added by LCM.

        function lgstr(str)
c-----
c       function to find the length of a string
c       this will only be used with file system path names
c       thus the first blank
c       indicates the end of the string
c-----
        implicit none
        character*(*) str
        integer*4 lgstr
        integer n, i
        n = len(str)
        lgstr = 1
        do 1000 i=n,1,-1
            lgstr = i
            if(str(i:i).ne.' ')goto 100
 1000   continue
  100   continue
        return
        end
!------------------------------------------------------------------------------!
        subroutine mgtarg(i,name)
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME V                                                       c
c                                                                     c
c      SUBROUTINE: MGTARG                                             c
c                                                                     c
c      COPYRIGHT 1996 R. B. Herrmann                                  c
c                                                                     c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
c-----
c       return the i'th command line argument
c
c       This version works on SUN, IBM RS6000
c-----
        implicit none
        integer i
        character name*(*)
            call getarg(i,name)
        return
        end
!------------------------------------------------------------------------------!
        function mnmarg()
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME V                                                       c
c                                                                     c
c      SUBROUTINE: MNMARG                                             c
c                                                                     c
c      COPYRIGHT 1996 R. B. Herrmann                                  c
c                                                                     c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
        implicit none
        integer iargc
        integer mnmarg
            mnmarg = iargc()
        return
        end
