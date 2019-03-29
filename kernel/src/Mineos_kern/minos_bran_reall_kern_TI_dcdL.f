c
c MINEOS version 1.0 by Guy Masters, John Woodhouse, and Freeman Gilbert
c
c This program is free software; you can redistribute it and/or modify
c it under the terms of the GNU General Public License as published by
c the Free Software Foundation; either version 2 of the License, or
c (at your option) any later version.
c
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c
c You should have received a copy of the GNU General Public License
c along with this program; if not, write to the Free Software
c Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
c
c****************************************************************************
c  minos_bran program produces the solution of the seismic normal mode
c  eigenvalue-eigenfunction problem.
c
c****************************************************************************
c  this program uses one input and two output files
c  the input file contains the model
c  the output files are 1) a model listing + a summary of mode properties
c  and 2) a file for the eigenfunctions
c  the program then asks for some control info described below
c
c structure of model file
c  card 1  :   title (up to 256 chars long)
c  card 2  :   ifanis,tref,ifdeck              (unformatted)
c           ifanis=1 for anisotropic model, 0 for isotropic
c           tref=ref period(secs) of model for dispersion correction.
c           if tref is .le. 0. no correction is made.
c           ifdeck=1 for card deck model, 0 for polynomial model.
c             *** card deck model ***
c  card 3  :   n,nic,noc                       (unformatted)
c           n=no of levels,nic=index of solid side of icb,noc=index of
c           fluid side of mcb.note that n must be .le. 223.
c  card 4-n:  r,rho,vpv,vsv,qkappa,qshear,vph,vsh,eta (f8.0,3f9.2,2f9.1,2f9.2,
c                                                           f9.5)
c           if isotropic vp=vpv,vs=vsv and vph,vsh,eta are unspecified.
c           if the q model is not specified then no disp. correction is made.
c           s.i. units are used,e.g. radius in metres and velocities in m/s.
c             *** polynomial model ***
c  card 3  :   nreg,nic,noc,rx                 (unformatted)
c           nreg is the number of regions in the model,nic and noc as before
c           and rx is the normalising radius for the polynomials.
c           rx is given in kms and is usually 6371.
c  card 4  :   nlay,r1,r2                      (unformatted)
c           nlay is the number of levels to be used in the region extending
c           from radius r1 to r2 (in kms).
c  card 5-9(iso) or card 5-12(ani) : coefs     (5f9.5)
c           5 sets of coefficients are required for each region of an isotropic
c           model and 8 sets for an anisotropic model.each polynomial can  be
c           up to a quartic in r/rx ( 5 coefficients) and are ordered thusly:
c           rho,vpv,vsv,qkappa,qshear,vph,vsh,eta. the coeffs are given in the
c           usual mixed seismological units (rho in g/cc,vel in km/s etc.)
c           conversion to s.i. is done by the program.
c           cards 4-9(iso) or 4-12(ani) are repeated for each region of the
c           model
c
c control parameters  (from screen)
c   eps,wgrav
c           eps controls the accuracy of the runge-kutta integration scheme.
c           the relative accuracy of an eigen frequency will be 2-3 x eps.
c           it also controls the precision with which a root is found and
c           the minimum relative separation of two roots with the same angular
c           order.it is safe to set eps=1.d-7.
c           wgrav is the frequency in millihertz above which gravitational terms
c           are neglected-this gives about a factor of 3 increase in speed.
c   jcom
c           jcom=1 radial modes, 2 toroidal modes, 3 spheroidal modes,
c           4 inner core toroidal modes.
c   lmin,lmax,wmin,wmax,nmin,nmax
c           lmin - lmax defines the range of angular orders to be computed.
c           if jcom=1 this is ignored. wmin - wmax defines the frequency range
c           to be computed (in millihertz)
c           nmin-nmax specifies the branch numbers to be computed -- nmin=0
c           is the fundamental mode
c
c  model listing
c    this is an ascii file which lists the model and mode properties
c    i.e. phase velocity in km/s,frequency,period,group velocity in km/s,
c    q and a parameter which is the ratio of kinetic to potential energy
c    minus 1 which should be small ( of order eps )if the eigenfunction is
c    accurate and if there are enough radial knots in the model to allow
c    quodratures to be done accurately. (you will probably see some degradation
c    in this parameter for strongly exponential modes such as stoneley modes ).
c
c  eigenfunction file
c****** file name "none" will suppress calculation of eigenfunctions
c    this is a fixed record length binary file with an entry for each mode
c    written as
c      write(ioeig) (abuf(i),i=1,nvec)
c    where nvec is 6+6*npts for spheroidal radial modes and 6+2*npts for toroidal
c    and modes. the first five words of abuf are n,l,frequecy,q and
c    group velocity. the rest of abuf contains w(1..npts) and wp(1..npts)
c    for toroidal modes and u(1..npts),up(1..npts),v(1..npts),vp(1..npts),
c    p(1..npts) and pp(1..npts) for spheroidal modes. these are as in
c    woodhouse and dahlen 1978 except that w,wp,v and vp must be divided
c    by sqrt(l*(l+1)). the normalisation is such that
c       frequency**2 times integral(rho*w*w*r*r) is 1 for toroidal modes
c       frequency**2 times integral(rho*(u*u+l*(l+1)*v*v)*r*r) is 1 for
c    spheroidal modes
c-------------------------------------------------------------------------
c  comment of H.Y.:
c  (1)
c     if eigenfunction s is expressed in the form,
c     s  = (U r^ + V/k surface_gradient - W/k surface_curl ) Y_lm
c      where Y_lm is real spherical harmonics
c
c    then spheroidal and toroidal mode are normalized in such way that
c      wrad_norm**2 * integral_[0,1] (rho_norm*(U**2+V**2)* r**2) dr = 1
c      wrad_norm**2 * integral_[0,1] (rho_norm*(W**2)     * r**2) dr = 1
c
c-------------------------------------------------------------------------
c
c    the model has been normalised such that a density of 5515 mg/m**3 = 1 ;
c    pi*g=1 where g is the grav constant (6.6723e-11) and the radius of the
c    earth (rn=6371000 m) is 1. these normalizations result in
c      acceleration normalisation = pi*g*rhobar*rn
c      velocity normalisation     = rn*sqrt(pi*g*rhobar)= vn
c      frequency normalization    = vn/rn
c------------------------------------------------------------------------
c  H.Y. Yang, June 2013
c    modified SUB model, modout
c    correct variable l+0.5 --> sfl3 for phase velocity
c    ****Note that the stored model for toroidal modes are different from my
c         original code***
c------------------------------------------------------------------------
c  H.Y. Yang, Aug 2013
c    Using cubic spline method to get P of radial modes instead of using
c    trapzoidal method
c------------------------------------------------------------------------
c  H.Y. Yang, Mar 2015
c  -- real l------------
c  change integer angular order l, lmin, lmax to real(KIND=r8b)
c  change input argument l for sub bsf to be real(KIND=r8b), instead of integer
c  add input argument dl for sub wtable
c
c  NOTE: avoid to calculate mode with angular degree l>0 and l<1
c
c  -- kernel ----
c  output kernel
c  add idfkern in argument of sub wtable, rotspl, model
c  add ichoix in argument of sub wtable, rotspl
c  allow to calculate TI kernel
c
c**************************************************************************
c
      PROGRAM mineos_reall_kern
      INTEGER, PARAMETER :: idfkern=20
      CHARACTER(LEN=300) :: fmodel, frename, funname, fkern
      INTEGER :: nmodes
      INTEGER :: nn0(2), jcom0
      INTEGER :: ichoix
      REAL*8  :: ll0(2), dl
      REAL*8 :: eps0, wgrav0, ww0(2)

      WRITE(*,*) 'Input model file:'
      READ(*,'(A)') fmodel
      WRITE(*,*) 'Output file name for modal summaries (ASCII):'
      READ(*,'(A)') frename
      WRITE(*,*) 'Output file name for modal eigenfunctions (binary):'
      READ(*,'(A)') funname
      WRITE(*,*) 'Enter eps and wgrav:'
      read(*,*) eps0,wgrav0
      WRITE(*,*) 'Enter jcom (1=rad;2=tor;3=sph;4=ictor)'
      read(*,*) jcom0
      print *,'Enter lmin,lmax,wmin,wmax,nmin,nmax,dl'
      read(*,*) ll0,ww0,nn0,dl
      WRITE(*,*) 'Output filename for sensitivity kernel (binary):'
      READ(*,*) fkern
      WRITE(*,*) 'Enter output kernel type (dw/dv=1, dc/dv=2, dc/dA=3):'
      READ(*,*) ichoix

c      eps0=max(eps0,1.d-12)
      ww0(1)=max(ww0(1),0.01d0)

      open(7,file=fmodel,status='old',form='formatted',iostat=iret,ACTION='READ')
      open(8,file=frename,form='formatted',iostat=iret,ACTION='WRITE')
      open(9,file=funname,form='unformatted',ACCESS='SEQUENTIAL',
     +        iostat=iret,ACTION='WRITE')
      WRITE(9) 'VFUN1.02'
      WRITE(9) jcom0,REAL(ww0(1)),REAL(ww0(2)),FLOOR(ll0(1)),CEILING(ll0(2)),
     +         REAL(wgrav0)

      OPEN(idfkern,FILE=fkern,FORM='UNFORMATTED',ACCESS='SEQUENTIAL',ACTION='WRITE')
!      WRITE(idfkern) 'KERN1.02'
      WRITE(idfkern) jcom0, REAL(ww0(1)), REAL(ww0(2)), FLOOR(ll0(1)),CEILING(ll0(2)),
     +               REAL(wgrav0)

      call model(7,8,9,idfkern)
      close(7)

      nmodes=0
      ifreq=1
      if(funname(1:4).eq.'none') ifreq=0
      call wtable(8,9,ifreq,nmodes,jcom0,ll0,ww0,nn0,eps0,wgrav0,
     +            dl,idfkern,ichoix)
      IF (nmodes.GT.0) THEN
        close(8)
        close(9)
        CLOSE(idfkern)
        WRITE(*,'(/I6,A)') nmodes,' modes in total.'
      ELSE
        CLOSE(8,STATUS='DELETE')
        CLOSE(9,STATUS='DELETE')
        CLOSE(idfkern,STATUS='DELETE')
      END IF
      stop
      END PROGRAM mineos_reall_kern

      subroutine wtable(iout,ioeig,ifreq,ncall,jcom0,ll0,ww0,nn0,eps0,wgrav0,
     +                  dl,idfkern,ichoix)
c*** makes up table of frequencies ***
      implicit real*8(a-h,o-z)

      INTEGER :: jcom0, nn0(2), idfkern,ichoix
      REAL*8  :: ll0(2),dl
      REAL*8  :: ww0(2),eps0, wgrav0

      REAL*8 :: l, lmin, lmax

      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      common/shanks/b(46),c(10),dx,step(8),stepf,maxo,in
      common/mtab/we(2),de(2),ke(2),wtry,bm
      dimension wt(2)
      data inss/5/
      common/counter/icount,ibad

      cmhz=pi/500.d0
      stepf=1.d0

      eps=eps0
      wgrav=wgrav0
      jcom=jcom0
      lmin=ll0(1)
      lmax=ll0(2)
      wmin=ww0(1)
      wmax=ww0(2)
      normin=nn0(1)
      normax=nn0(2)

      eps1=eps
      eps2=eps
      wgrav=wgrav*cmhz
c      write(iout,100) eps,eps1,wgrav
c  100 format(/,'integration precision =',g12.4,'  root precision =',
c     +   g12.4,'  gravity cut off =',g12.4,' rad/s',///,6x,'mode',
c     +   8x,'phs vel',7x,'w(mhz)',10x,'t(secs)',6x,'grp vel(km/s)',
c     +   8x,'q',13x,'raylquo',/)
      WRITE(*,100) eps, eps1, wgrav
100   format(/,5X,'integration precision =',g12.4,'  root precision =',
     +   g12.4,'  gravity cut off =',g12.4,' rad/s',//)
      call steps(eps)
      if(jcom.lt.1.or.jcom.gt.4) jcom=3

      wmin=wmin*cmhz
      wmax=wmax*cmhz

c            To avoid modes with angular degree >0 & <1 : H.Y., Jan 2016
c                     Debug
      IF (lmin.LE.1.d0 .AND. dl.LT.1.d0) lmin=1.d0-dl
c
      if(jcom.eq.1) then
        lmin=0.d0
        lmax=0.d0
      end if
      normin=max(normin,0)
      normax=max(normax,normin)
	ncall = 0
      do 50 nor=normin,normax
      wt(1)=wmin
      wt(2)=wmax
      wtry=0.d0
      bm=0.d0
      icount=0
      ibad=0

      nlmax=((lmax-lmin)/dl)+1
      DO 10 il=1,nlmax
      l=lmin+(il-1)*dl
      if(wt(1).ge.wmax) goto 50
      nord=nor
      nup=nor
      ndn=nor-1
      if(l.eq.1) then
        nup=ndn
        ndn=ndn-1
      end if
      knsw=1
      maxo=inss
      fl=l
      fl1=fl+1.d0
      fl2=fl+fl1
      fl3=fl*fl1
      sfl3=sqrt(fl3)
      we(1)=wt(1)
      we(2)=wt(2)
      if(wtry.ne.0.d0) we(2)=wtry
      call detqn(we(1),ke(1),de(1),0)
      if(ke(1).gt.ndn) goto 10
      call detqn(we(2),ke(2),de(2),0)
      if(ke(2).lt.nup) then
         we(2)=wt(2)
         call detqn(we(2),ke(2),de(2),0)
         if(ke(2).lt.nup) goto 50
      end if
c*** bracket this mode ***
      wx=0.5d0*(we(1)+we(2))
      ktry=0
   15 if(ke(1).eq.ndn.and.ke(2).eq.nup) goto 40
      ktry=ktry+1
      if(ktry.gt.50) goto 10
      call detqn(wx,kx,dx,0)
      if(kx.le.ndn) then
        we(1)=wx
        ke(1)=kx
        de(1)=dx
      else
        we(2)=wx
        ke(2)=kx
        de(2)=dx
      end if
      wx=0.5d0*(we(1)+we(2))
      goto 15
c*** find roots ***
   40 knsw=0
      maxo=8
c		added argument for number of times called (mhr)
	ncall = ncall + 1
      call rotspl(eps1,wt,iout,ioeig,ifreq,ncall,idfkern,ichoix)
   10 continue

      WRITE(*,900) nor,wt(1)/cmhz,wt(2)/cmhz, icount, ibad
900   FORMAT(' radial order n =', I5,'  (',
     +         F7.2,' - ', F7.2,' mHz)  # =', I6, ', bad =', I4 )

   50 continue
      return
      end

      subroutine rotspl(eps1,wt,iout,ioeig,ifreq,ncall,idfkern,ichoix)
c*** find roots by spline interpolation ***
      implicit real*8(a-h,o-z)
      REAL*8 :: l
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      common/mtab/we(2),de(2),ke(2),wtry,bm
      dimension x(20),det(20),qx(3,20),wrk(60),wt(*),kchar(4)
      character*2 kchar
      data tol/1.d-9/,itmax/15/,kchar/' S',' T',' S',' C'/
      if(de(1)*de(2).gt.0.d0) return
      nord=ke(2)
      if(l.eq.1) nord=nord+1
      det(1)=de(1)
      det(2)=de(2)
      x(1)=we(1)
      x(2)=we(2)
    5 c=x(1)
      if(dabs(det(1)).gt.dabs(det(2))) c=x(2)
   10 ntry=2
      grad=(x(1)-x(2))/(det(1)-det(2))
      b=x(1)-det(1)*grad
   15 t=dabs(b*eps1)
      if(dabs(b-c).lt.t) goto 65
      call detqn(b,knt,fb,0)
      ind=1
      do 20 m=2,ntry
      ind=ind+1
   20 if(b.lt.x(m)) goto 25
   25 ntry=ntry+1
      j2=ntry
   30 j1=j2-1
      x(j2)=x(j1)
      det(j2)=det(j1)
      if(j1.eq.ind) goto 35
      j2=j2-1
      goto 30
   35 x(ind)=b
      det(ind)=fb
      idn=0
      do 40 m=2,ntry
      idn=idn+1
      iup=idn+1
   40 if(det(idn)*det(iup).le.0.d0) goto 45
   45 ind=iup
      if(dabs(det(idn)).lt.dabs(det(iup))) ind=idn
      c=x(ind)
      if(ntry.ge.itmax) goto 60
      call dsplin(ntry,x,det,qx,wrk)
      del=-det(ind)/qx(1,ind)
   50 delx=-det(ind)/(qx(1,ind)+del*qx(2,ind))
      if(dabs(delx-del).lt.tol) goto 55
      if(del*delx.lt.0.d0) goto 60
      del=delx
      goto 50
   55 b=c+delx
      if(b.ge.x(idn).and.b.le.x(iup)) goto 15
   60 x(1)=x(idn)
      x(2)=x(iup)
      det(1)=det(idn)
      det(2)=det(iup)
      goto 10
c
c
c                             write out frequencies ***
c
   65 call detqn(b,knt,fb,ifreq)
      tcom=2.d0*pi/b
      wmhz=1000.d0/tcom
c      cvel=b*rn/(l+0.5)/1000.
      cvel=b*(rn/1000.d0)/SQRT((l+1.d0)*l)
      if(ifreq.eq.1) then
        wdiff=(b-wray*wn)/b
        gcom=vn*cg/1000.d0
        qmod=0.d0
        if(qinv.gt.0.d0) qmod=1.d0/qinv
c        write(iout,200) nord,kchar(jcom),l,cvel,wmhz,tcom,gcom,qmod,
c     +      wdiff
c  200   format(i5,a2,i5,6g16.7)
         WRITE(iout,200) nord,kchar(jcom),l,b,wmhz,tcom,gcom,qmod,wdiff
200      format(i5,a2,E14.6,6g16.7)
c         WRITE(*,200) nord,kchar(jcom),l,b,wmhz,tcom,gcom,qmod,wdiff

        call modout(b,qmod,gcom,ioeig,ncall,wdiff)

c
c                           output kernel
c
        CALL kernel_velrhopTI(ichoix,idfkern,b,qmod,gcom,wdiff)

      else
c                              we estimate group velocity
c                                 using previous frequencies
c
        cg=5000./vn
        if(bm.ne.0.d0) cg=(b-bm)/wn
        gcom=vn*cg/1000.d0
c		added print of gcom. qmod & wdiff not computed w/o egnfcns. (mhr)
c       write(iout,200) nord,kchar(jcom),l,cvel,wmhz,tcom,gcom
        WRITE(iout,201) nord,kchar(jcom),l,b,wmhz,tcom,gcom
201     format(I5,A2,E14.6,4g16.7)
        WRITE(*,201) nord,kchar(jcom),l,b,wmhz,tcom,gcom
      end if
      wtry=b+2.d0*cg*wn
      wt(1)=b
      bm=b
      return
      end

      subroutine modout(wcom,qmod,gcom,ioeig,ncall,wdiff)
c
c     H.Y. Yang, June 2013
c     different from HY original code
c      nrec: # of stored radial knots for model
c     nvec: header + stored eigenfunctions as a function of radius
c           # of radial knots = nrec        for spheroidal and radial modes
c                             = nrec-nrecoc for toroidal modes
c                             = nrecic      for inner-core toroidal modes
c
c
        implicit real*8(a-h,o-z)
        parameter (mk=30000, nhead=6)
        integer krec(mk)
        dimension abuf(6*mk+nhead)

        REAL*8 :: l
        common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     1    fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
        common/rindx/nic,noc,nsl,nicp1,nocp1,nslp1,n
        common/eifx/a(14,mk),ndum(mk)
        common/eigsav/krec,nrec,nrecic,nrecoc,nkrec
        common/counter/icount,ibad

        common r(mk),fmu(mk),flam(mk),qshear(mk),
     1    qkappa(mk),xa2(mk),xlam(mk),rho(mk),
     2    qro(3,mk),g(mk),qg(3,mk),fcon(mk),
     3    fspl(3,mk),lcon(mk),lspl(3,mk),ncon(mk),
     4    nspl(3,mk),ccon(mk),cspl(3,mk),acon(mk),
     5    aspl(3,mk)

c        REAL*8 pU(mk),SO(mk),dS,sumall
c        REAL rrec(mk)
        REAL*8 pU(mk), qpU(3,mk), wrk(3,mk), sumall(mk),sumtop,sum1
        INTEGER nr, nvec

        if(jcom.eq.1) then
          nr=nrec
          nvec=nhead+6*nr
        elseif(jcom.eq.2) then
          nr=nrec-nrecoc
          nvec=nhead+2*nr
        elseif(jcom.eq.3) then
          nr=nrec
          nvec=nhead+6*nr
        elseif(jcom.eq.4) then
          nr=nrecic
          nvec=nhead+2*nr
        else
          STOP 'jcom is out of range, STOP'
        endif
        nkrec=nr

        if(jcom.eq.2) then
          do i=1,noc
            a(1,i)=0.d0
            a(2,i)=0.d0
          end do
        elseif(jcom.eq.4) then
          do i=nic+1,n
            a(1,i)=0.d0
            a(2,i)=0.d0
          end do
        else
        endif

        abuf(1)=dble(nord)
        abuf(2)=dble(l)
        abuf(3)=dble(wcom)
        abuf(4)=dble(qmod)
        abuf(5)=dble(gcom)
        abuf(6)=dble(wdiff)
        jst=nhead
        j=0
        if(jcom.eq.1) then
c          SO=0.d0
c          DO i=1,n
c            IF(krec(i).EQ.1) THEN
c              j=j+1
c              pU(j)=rho(i)*a(1,i)
c              rrec(j)=r(i)
c           Using trapzoidal method to get integrand
c              dS = 0.d0
c              IF ( j>=2 ) THEN
c                j1=j-1
c                dS=(pU(j)+pU(j1))*(rrec(j)-rrec(j1))/2.d0
c                SO(j)=SO(j1)+dS
c              ENDIF
c            ENDIF
c          END DO
c          sumall=SO(j)
c          factor=4.d0
c          j=0
c          do i=1,n
c            if(krec(i).eq.1) then
c              j=j+1
c              abuf(j+jst)=a(1,i)
c              abuf(j+nr+jst)=a(2,i)
c              abuf(j+2*nr+jst)=0.d0
c              abuf(j+3*nr+jst)=0.d0
c              abuf(j+4*nr+jst)=factor*(sumall-SO(j))
c              abuf(j+5*nr+jst)=-factor*pU(j)
c            endif
c          end do
c
          factor=4.d0
          DO i=1,n
            pU(i)=rho(i)*a(1,i)
          END DO
          CALL drspln(1,n,r,pU,qpU,wrk)
          CALL intgrl(r,1,n,pU,qpU(1,1:n),qpU(2,1:n),qpU(3,1:n),sumall(1:n))
          sumtop=sumall(n)
          j=0
          DO i=1,n
            IF (krec(i).EQ.1) THEN
              sum1=sumtop-sumall(i)
              j=j+1
              abuf(j+jst)=a(1,i)
              abuf(j+nr+jst)=a(2,i)
              abuf(j+2*nr+jst)=0.d0
              abuf(j+3*nr+jst)=0.d0
              abuf(j+4*nr+jst)=factor*sum1
              abuf(j+5*nr+jst)=-factor*pU(i)
            END IF
          END DO

        elseif(jcom.eq.2) then
          do i=noc+1,n
            if(krec(i).eq.1) then
              j=j+1
              abuf(j+jst)=a(1,i)
              abuf(j+nr+jst)=a(2,i)
            endif
          end do
        elseif(jcom.eq.3) then
          do i=1,n
            if(krec(i).eq.1) then
              j=j+1
              abuf(j+jst)=a(1,i)
              abuf(j+nr+jst)=a(2,i)
              abuf(j+2*nr+jst)=a(3,i)
              abuf(j+3*nr+jst)=a(4,i)
              abuf(j+4*nr+jst)=a(5,i)
              abuf(j+5*nr+jst)=a(6,i)
            endif
          end do
        elseif(jcom.eq.4) then
          do i=1,nic
            if(krec(i).eq.1) then
              j=j+1
              abuf(j+jst)=a(1,i)
              abuf(j+nr+jst)=a(2,i)
            endif
          end do
        endif
c       WRITING EIGENFUNCTIONSS

        write(ioeig) (abuf(i),i=1,nvec)
        call flush(ioeig)
        icount=icount+1
        IF (ABS(wdiff)>=1.d-5) ibad=ibad+1

        return
      end

      subroutine model(iin,iout,ioeig,idfkern)
c
c  modified by H.Y. Yang, June 2013
c
c  Note that knots below the cmb is also stored for toroidal modes. (!!!)
c
c
      implicit real*8(a-h,o-z)
      parameter (mk=30000, nhead=6)
      save
      CHARACTER*128 :: ktitle
      integer*4 ititle(20)
      real*8 lcon,ncon,lspl,nspl

      integer krec(mk)
      real*4 rrec(mk),rhorec(mk),vpvrec(mk),vsvrec(mk),vphrec(mk)
      real*4 vshrec(mk),qkapparec(mk),qshearrec(mk),etarec(mk)
      REAL*8 :: l

      common r(mk),fmu(mk),flam(mk),qshear(mk),qkappa(mk),
     + xa2(mk),xlam(mk),rho(mk),qro(3,mk),g(mk),qg(3,mk),
     + fcon(mk),fspl(3,mk),lcon(mk),lspl(3,mk),ncon(mk),
     + nspl(3,mk),ccon(mk),cspl(3,mk),acon(mk),aspl(3,mk)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      common/eifx/vpv(mk),vph(mk),vsv(mk),vsh(mk),eta(mk),wrk(mk*10)
      common/rindx/nic,noc,nsl,nicp1,nocp1,nslp1,n
      common/eigsav/krec,nrec,nrecic,nrecoc,nkrec

      data bigg,tau/6.6723d-11,1.d3/,rhobar/5515.d0/
      pi=3.14159265358979d0

      krec=0

      read(iin,'(A)') ktitle
      read(iin,*) ifanis,tref,ifdeck
      if(ifdeck.eq.0) go to 1000

c*** card deck model ***
      read(iin,*) n,nic,noc, lrec
      read(iin,*) (r(i),rho(i),vpv(i),vsv(i),
     +     qkappa(i),qshear(i),vph(i),vsh(i),eta(i),i=1,n)

       irstart=1
        if(lrec.eq.0) then
c                                            store whole radial knots, but
          do i=irstart,n
            krec(i)=1
          end do
        elseif(lrec.eq.1) then
c                                            store knots with integer radius,
c
          i1=irstart
          do j=0,NINT(r(n)/tau)
            rcur=real(j)*1.e3
            xdiff=1.e20
            do i=i1,n
              if(xdiff.gt.abs(rcur-real(r(i)))) then
                xdiff=abs(rcur-real(r(i)))
                ki=i
              endif
            end do
            krec(ki)=1
            i1=ki
          end do
          do i=irstart,n-1
            if(r(i).eq.r(i+1)) then
              krec(i)=1
              krec(i+1)=1
            endif
          end do
          krec(irstart)=1
          krec(n)=1
        endif
      go to 2000

c*** polynomial model ***
 1000 read(iin,*) nreg,nic,noc,rx
      rx=rx*tau
      n=0
      knt=0
      jj=8 !5
      if(ifanis.ne.0) jj=8
      print*,ifanis,nreg,nic,noc,rx
      do nn=1,nreg
        read(iin,*) nlay,r1,r2
        r1=r1*tau
        r2=r2*tau
        dr=(r2-r1)/float(nlay-1)
        do i=1,nlay
          n=n+1
          r(n)=r1+dr*float(i-1)
        enddo
       print*,nlay,r1,r2
        do j=1,jj
          read(iin,*) (wrk(i),i=1,5)
          print*,wrk(1:5)
  110     format(5f9.5)
          if(ifanis.eq.2) then
            do i=1,nlay
              ind=knt+i
              rt=r(ind)/rx
              val=wrk(1)+rt*(wrk(2)+rt*(wrk(3)+rt*(wrk(4)+rt*wrk(5))))
              if(j.eq.1) rho(ind)=val*tau
              if(j.eq.2) vph(ind)=val*tau
              if(j.eq.3) vsv(ind)=val*tau
              if(j.eq.4 .and. val.gt.0.) qshear(ind)=1./val
              if(j.eq.5 .and. val.gt.0.) qkappa(ind)=1./val
              if(j.eq.6) vpv(ind)=sqrt(val)*vph(ind)
              if(j.eq.7) vsh(ind)=sqrt(val)*vsv(ind)
              if(j.eq.8) eta(ind)=val
            enddo
          else
            do i=1,nlay
              ind=knt+i
              rt=r(ind)/rx
              val=wrk(1)+rt*(wrk(2)+rt*(wrk(3)+rt*(wrk(4)+rt*wrk(5))))
              if(j.eq.1) rho(ind)=val*tau
              if(j.eq.2) vpv(ind)=val*tau
              if(j.eq.3) vsv(ind)=val*tau
              if(j.eq.4 .and. val.gt.0.) qshear(ind)=1./val
              if(j.eq.5 .and. val.gt.0.) qkappa(ind)=1./val
              if(j.eq.6) vph(ind)=val*tau
              if(j.eq.7) vsh(ind)=val*tau
              if(j.eq.8) eta(ind)=val
            enddo
          end if
        enddo
        knt=knt+nlay
      enddo

 2000 if(ifanis.eq.0) then
      do i=1,n
        vph(i)=vpv(i)
        vsh(i)=vsv(i)
        eta(i)=1.d0
      enddo
      end if

      if(iout.ge.0) then
c*** write out model ***

        j=0
        nrecic=0
        nrecoc=0
        do i=1,n
          if(krec(i).eq.1) then
            j=j+1
            rrec(j)=real(r(i))
            rhorec(j)=real(rho(i))
            vpvrec(j)=real(vpv(i))
            vsvrec(j)=real(vsv(i))
            vphrec(j)=real(vph(i))
            vshrec(j)=real(vsh(i))
            qkapparec(j)=real(qkappa(i))
            qshearrec(j)=real(qshear(i))
            etarec(j)=real(eta(i))
            if(i.eq.nic) nrecic=j
            if(i.eq.noc) nrecoc=j
          endif
        end do
        nrec=j

c                                   write model into binary file
        WRITE(*,*) 'ifanis 1:', ifanis
        write(ioeig) nrec,nrecic,nrecoc,ifanis,real(tref),(rrec(i),i=1,nrec),
     1    (rhorec(i),i=1,nrec),(vpvrec(i),i=1,nrec),
     2    (vsvrec(i),i=1,nrec),(vphrec(i),i=1,nrec),
     3    (vshrec(i),i=1,nrec),(qkapparec(i),i=1,nrec),
     4    (qshearrec(i),i=1,nrec),(etarec(i),i=1,nrec)
        call flush(ioeig)

        write(idfkern) nrec,nrecic,nrecoc,ifanis,real(tref),(rrec(i),i=1,nrec),
     1    (rhorec(i),i=1,nrec),(vpvrec(i),i=1,nrec),
     2    (vsvrec(i),i=1,nrec),(vphrec(i),i=1,nrec),
     3    (vshrec(i),i=1,nrec),(qkapparec(i),i=1,nrec),
     4    (qshearrec(i),i=1,nrec),(etarec(i),i=1,nrec)
        call flush(idfkern)

      end if

c*** normalise and spline ***
      rn=r(n)
      gn=pi*bigg*rhobar*rn
      vn2=gn*rn
      vn=sqrt(vn2)
      wn=vn/rn
      do i=1,n
        r(i)=r(i)/rn
        if(i.gt.1)then
          if(dabs(r(i)-r(i-1)).lt.1.d-7) r(i)=r(i-1)
        endif
        if(qshear(i).gt.0.d0) qshear(i)=1.d0/qshear(i)
        if(qkappa(i).gt.0.d0) qkappa(i)=1.d0/qkappa(i)
        rho(i)=rho(i)/rhobar
        acon(i)=rho(i)*vph(i)*vph(i)/vn2
        ccon(i)=rho(i)*vpv(i)*vpv(i)/vn2
        lcon(i)=rho(i)*vsv(i)*vsv(i)/vn2
        ncon(i)=rho(i)*vsh(i)*vsh(i)/vn2
        fcon(i)=eta(i)*(acon(i)-2.d0*lcon(i))
        fmu(i)=(acon(i)+ccon(i)-2.d0*fcon(i)+5.d0*ncon(i)
     +   +6.d0*lcon(i))/15.d0
        flam(i)=(4.d0*(acon(i)+fcon(i)-ncon(i))+ccon(i))/9.d0
     +    -2.d0*fmu(i)/3.d0
        rat=4.d0*fmu(i)/(3.d0*(flam(i)+2.d0*fmu(i)))
        xlam(i)=((1.d0-rat)*qkappa(i)-.5d0*rat*qshear(i))/(1.d0
     +    -1.5d0*rat)
        xa2(i)=(1.d0-rat)*qkappa(i)+rat*qshear(i)
      enddo
      call drspln(1,n,r,rho,qro,wrk)
c      OPEN(UNIT=20,FILE='gspl.out')
c      WRITE(UNIT=20,FMT='(3E15.7)') (qro(1:3,i),i=1,n)
c      CLOSE(UNIT=20)
      call grav(g,rho,qro,r,n)
c      OPEN(UNIT=21,FILE='grav.out')
c      WRITE(UNIT=21,FMT='(E15.7)') (g(i),i=1,n)
c      CLOSE(UNIT=21)
      call drspln(1,n,r,g,qg,wrk)
      call drspln(1,n,r,fcon,fspl,wrk)
      call drspln(1,n,r,lcon,lspl,wrk)
      if(ifanis.ne.0) then
        call drspln(1,n,r,acon,aspl,wrk)
        call drspln(1,n,r,ccon,cspl,wrk)
        call drspln(1,n,r,ncon,nspl,wrk)
      end if
      nsl=n+1
   65 nsl=nsl-1
      if(vsv(nsl).le.0.d0) go to 65
      nicp1=nic+1
      nocp1=noc+1
      nslp1=nsl+1
      tref=0.5d0*tref/pi

c
c                         write gravity (m/s^2) at surface on screen
c      WRITE(ioeig) g(n)*gn
      WRITE(*,'(/,6X,A,2(E23.16,A/))') 'gravity at surface =', g(n)*gn ,' m/s^2',
     #  g(n),' non-dimensional'
      WRITE(*,*) 'rho, C, L=', rho(n),ccon(n),lcon(n)
c

      return
      end

      subroutine detqn(wdim,knt,det,ifeif)
c**** supevises the integration of the equations,it returns the value
c**** of the secular determinant as det and the count of zero crossings.
      implicit real*8(a-h,o-z)
      parameter (mk=30000)
      save
      real*8 lcon,ncon,lspl,nspl
      REAL*8 :: l
      common r(mk),fmu(mk),flam(mk),qshear(mk),qkappa(mk),
     + xa2(mk),xlam(mk),rho(mk),qro(3,mk),g(mk),qg(3,mk),
     + fcon(mk),fspl(3,mk),lcon(mk),lspl(3,mk),ncon(mk),
     + nspl(3,mk),ccon(mk),cspl(3,mk),acon(mk),aspl(3,mk)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      common/eifx/a(14,mk),dum(mk)
      common/rindx/nic,noc,nsl,nicp1,nocp1,nslp1,n
      dimension ass(14),vf(mk),zi(4)
      iback=0
      w=wdim/wn
      wsq=w*w
      iexp=0
      kount=0
      kg=0
      fct=0.d0
      if(tref.gt.0.d0) fct=2.d0*dlog(tref*wdim)/pi
      goto (2,3,1,3),jcom
    1 if(wdim.le.wgrav) kg=1
      nvefm=2+kg*3
      nvesm=5+kg*9
      call sdepth(wdim,ls)
      if(ls.le.2) then
        r10=4.5d-4*(fl+.5d0)/wdim
        if(r10.lt.r(2)) then
          r(1)=r10
          g(1)=rho(1)*r(1)*1.333333333333333d0
          ls=1
        end if
      end if
      if(ls.le.nic) then
c*** propagate through inner core ***
      call spsm(ls,nvesm,ass)
        call sprpmn(ls,nic,ass,vf,nvesm,iexp)
        r(1)=0.d0
        g(1)=0.d0
        call sfbm(ass,kg,iback)
      end if
      if(ls.le.noc) then
c*** propagate through outer core ***
        is=max(ls,nicp1)
        if(is.eq.ls) call fpsm(ls,nvefm,ass)
        call fprpmn(is,noc,ass,vf,nvefm,iexp)
        call fsbm(ass,kg,iback)
      end if
      is=max(ls,nocp1)
      if(is.eq.ls) call spsm(ls,nvesm,ass)
c*** propagate through mantle ***
      call sprpmn(is,nsl,ass,vf,nvesm,iexp)
      if(nsl.eq.n) then
        dnorm=a(1,nsl)*a(1,nsl)
        do i=2,nvesm
          dnorm=dnorm+a(i,nsl)*a(i,nsl)
        enddo
        det=a(5,nsl)/sqrt(dnorm)
      else
        call sfbm(ass,kg,iback)
c*** propagate through ocean ***
        call fprpmn(nslp1,n,ass,vf,nvefm,iexp)
        if(kg.eq.0) then
           det=a(2,n)/sqrt(a(1,n)*a(1,n)+a(2,n)*a(2,n))
        else
          det=a(5,n)/sqrt(a(1,n)**2+a(2,n)**2+a(3,n)**2+
     +   a(4,n)**2+a(5,n)**2)
        end if
      end if
      if(ls.gt.noc) det=-det
      if(knsw.eq.1) then
        if(ls.gt.noc) kount=kount-2
        irem=mod(kount,2)
        if(irem.eq.0.and.det.lt.0.d0) kount=kount+1
        if(irem.ne.0.and.det.gt.0.d0) kount=kount+1
        knt=kount
      end if
      if(ifeif.eq.0) return
c*** this does eigenfunction calculation for spheroidal modes ***
      iback=1
      jexp=0
      nbakf=1+kg*3
      nbaks=4+kg*10
      do i=1,nbaks
        ass(i)=0.d0
      enddo
      if(n.ne.nsl) then
        if(kg.eq.0) then
          ass(1)=dsign(1.d0,a(1,n))
        else
          asi1=a(3,n)*a(3,n)
          asi2=a(4,n)*a(4,n)
          if(asi2.le.asi1) ass(1)=dsign(1.d0,a(3,n))
          if(asi2.gt.asi1) ass(2)=dsign(1.d0,a(2,n))
        end if
        call fprpmn(n,nslp1,ass,vf,nbakf,jexp)
        call fsbm(ass,kg,iback)
      else
        asi1=a(3,n)*a(3,n)
        asi2=a(4,n)*a(4,n)
        if(kg.ne.0) then
          asi1=asi1+a(12,n)*a(12,n)
          asi2=asi2+a(11,n)*a(11,n)
        end if
        if(asi2.le.asi1) ass(1)=dsign(1.d0,a(3,n))
        if(asi2.gt.asi1) ass(2)=dsign(1.d0,a(2,n))
      end if
      nto=max(ls,nocp1)
      call sprpmn(nsl,nto,ass,vf,nbaks,jexp)
      if(nto.eq.ls) goto 90
      call sfbm(ass,kg,iback)
      nto=max(ls,nicp1)
      call fprpmn(noc,nto,ass,vf,nbakf,jexp)
      if(nto.eq.ls) goto 90
      call fsbm(ass,kg,iback)
      nto=max(ls,2)
      call sprpmn(nic,nto,ass,vf,nbaks,jexp)
   90 if(dabs(det).gt.1.d-4.and.ls.le.noc) call remedy(ls)
      call eifout(ls)
      return
c*** radial modes ***
    2 ls=2
      call rps(ls,ass)
      call rprop(ls,n,ass)
      det=a(2,n)/dsqrt(a(1,n)*a(1,n)+a(2,n)*a(2,n))
      knt=kount-1
      if(ifeif.eq.0) return
      a(1,1)=0.d0
      a(2,1)=0.d0
      do i=ls,n
        ff=fcon(i)*(1.d0+xlam(i)*fct)
        cc=ccon(i)*(1.d0+xa2(i)*fct)
        a(2,i)=(a(2,i)-2.d0*ff*a(1,i)/r(i))/cc
      enddo
      zi(1)=0.d0
      zi(2)=0.d0
      zi(3)=0.d0
      do i=ls,n
        im=i-1
        if(r(i).ne.r(im)) call gauslv(r(im),r(i),im,zi,3)
      enddo
      rnrm=1.d0/(w*dsqrt(zi(1)))
      cg=0.d0
      qinv=zi(2)/(wsq*zi(1))
      wray=dsqrt(zi(3)/zi(1))
      do i=2,n
        do j=1,2
          a(j,i)=a(j,i)*rnrm
        enddo
      enddo
      return
c*** toroidal modes ***
    3 if(jcom.eq.2) then
        nb=nocp1
        n2=nsl
c---------------------------
c TEST
c      WRITE(*,*) 'before eig:', a(1:2,n-2)
c----------------------

        ass(1)=1.d0
        ass(2)=0.d0
      else
        nb=2
        a(1,1)=0.d0
        a(2,1)=0.d0
        n2=nic
      end if
      q=0.d0
      ls=nb
      IF (r(ls)==0.d0) ls=ls+1
      call startl(ls,n2,fmu,ls,q)
      if(ls.ne.nocp1) call tps(ls,ass)
      call tprop(ls,n2,ass)
      det=a(2,n2)/dsqrt(a(1,n2)*a(1,n2)+a(2,n2)*a(2,n2))
      if(knsw.eq.1) then
        knt=kount-2
        if(jcom.ne.4) then
          if(l.ne.1) then
            knt=knt+1
            irem=mod(knt,2)
            if(irem.eq.0.and.det.gt.0.d0) knt=knt+1
            if(irem.ne.0.and.det.lt.0.d0) knt=knt+1
          end if
        end if
      end if
      if(ifeif.eq.0) return
      do i=ls,n2
        a(2,i)=a(1,i)/r(i)+a(2,i)/(lcon(i)*(1.d0+qshear(i)*fct))
      enddo
      if(ls.ne.nb) then
        do i=nb,ls-1
          a(1,i)=0.d0
          a(2,i)=0.d0
        enddo
      end if
      do i=1,4
        zi(i)=0.d0
      enddo
      do i=ls,n2
        im=i-1
        if(r(i).ne.r(im)) call gauslv(r(im),r(i),im,zi,4)
      enddo
      rnrm=1.d0/(w*dsqrt(zi(1)))
c      cg=(fl+0.5d0)*zi(2)/(w*zi(1))
      cg=sfl3*zi(2)/(w*zi(1))
      qinv=zi(3)/(wsq*zi(1))
      wray=dsqrt(zi(4)/zi(1))
c-----------------------
c TEST
c      WRITE(*,*) 'D=1km:', a(1:2,n-2), 'nsl-1:', nsl-1, a(1:2,nsl-1)
c----------------------

      do i=ls,n2
        do j=1,2
          a(j,i)=a(j,i)*rnrm
        enddo
      enddo
      return
      end

      subroutine sprpmn(jf,jl,f,h,nvesm,iexp)
c*** propagate a minor vector in a solid region from level jf to jl ***
      implicit real*8(a-h,o-z)
      parameter (mk=30000)
      save
      real*8 lcon,ncon,lspl,nspl
      REAL*8 :: l
      common r(mk),fmu(mk),flam(mk),qshear(mk),qkappa(mk),
     + xa2(mk),xlam(mk),rho(mk),qro(3,mk),g(mk),qg(3,mk),
     + fcon(mk),fspl(3,mk),lcon(mk),lspl(3,mk),ncon(mk),
     + nspl(3,mk),ccon(mk),cspl(3,mk),acon(mk),aspl(3,mk)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      common/eifx/ar(14,mk),inorm(mk),jjj(mk)
      common/shanks/b(46),c(10),dx,step(8),stepf,maxo,in
      dimension f(*),h(nvesm,*),s(14),fp(14),rne(6)
      data econst/1048576.d0/
      maxo1=maxo-1
      jud=1
      if(jl.lt.jf) jud=-1
      y=r(jf)
      i=jf
      go to 45
   10 x=y
      y=r(i)
      if(y.eq.x) goto 45
      iq=min(i,i-jud)
      qff=1.d0+xlam(iq)*fct
      qll=1.d0+qshear(iq)*fct
      qaa=1.d0+xa2(iq)*fct
      xi=g(i)/y
      vpsq=(flam(i)+2.d0*fmu(i))/rho(i)
      vssq=fmu(i)/rho(i)
      alfsq=(wsq+4.d0*rho(i)+xi)/vpsq
      betasq=wsq/vssq
      delsq=sqrt((betasq-alfsq)**2+4.d0*fl3*xi*xi/(vssq*vpsq))
      fksq=.5d0*(alfsq+betasq+delsq)-fl3/(x*x)
      qt=sqrt(abs(fksq))+sqrt(abs(fksq-delsq))+2.d0/r(iq)
      q=(qt+float(kg)*sfl3/x)/stepf
      del=float(jud)*step(maxo)/q
      dxs=0.d0
   15   y=x+del
        if(float(jud)*(y-r(i)).gt.0.d0) y=r(i)
        dx=y-x
        if(dx.ne.dxs) call baylis(q,maxo1)
        dxs=dx
        do j=1,nvesm
          s(j)=f(j)
        enddo
        do ni=1,in
          z=x+c(ni)
          call derms(iq,z,f,h(1,ni),0,qff,qll,qaa)
          call rkdot(f,s,h,nvesm,ni)
        enddo
        if(knsw.eq.1) then
          call derms(iq,y,f,fp,1,qff,qll,qaa)
          call zknt(s,h,f,fp,x,y,1)
        end if
        x=y
        if(y.ne.r(i)) goto 15
   45 size=abs(f(1))
      do j=2,nvesm
        size=max(size,dabs(f(j)))
      enddo
   55 if(size.lt.1024.d0) goto 65
      do j=1,nvesm
        f(j)=f(j)/econst
      enddo
      size=size/econst
      iexp=iexp+20
      goto 55
   65 if(iback.ne.0) then
        inorm(i)=inorm(i)+iexp
        if(kg.ne.0) then
          t1=f(4)+f(8)
          t2=t1+f(4)
          t1=t1+f(8)
          t3=f(8)-f(4)
          rne(1)=ar(6,i)*f(10)-ar(14,i)*f(9)+ar(13,i)*t3
     1          -ar(1,i)*f(7)-ar(7,i)*f(6)+ar(8,i)*f(5)
     2          +ar(12,i)*f(3)-ar(2,i)*f(2)+ar(3,i)*f(1)
          rne(2)=ar(6,i)*f(13)+ar(14,i)*t2+ar(13,i)*f(12)
     1          -ar(1,i)*f(11)-ar(9,i)*f(6)-ar(7,i)*f(5)
     2          +ar(11,i)*f(3)-ar(4,i)*f(2)-ar(2,i)*f(1)
          rne(3)=ar(6,i)*f(14)-ar(7,i)*t1-ar(8,i)*f(12)
     1          +ar(13,i)*f(11)-ar(9,i)*f(9)+ar(14,i)*f(7)
     2          +ar(10,i)*f(3)+ar(11,i)*f(2)+ar(12,i)*f(1)
          rne(4)=ar(14,i)*f(14)+ar(7,i)*f(13)+ar(12,i)*f(12)
     1          -ar(2,i)*f(11)-ar(9,i)*f(10)-ar(11,i)*t3
     2          +ar(4,i)*f(7)+ar(10,i)*f(5)+ar(5,i)*f(1)
          rne(5)=ar(13,i)*f(14)+ar(8,i)*f(13)-ar(12,i)*t2
     1          -ar(3,i)*f(11)+ar(7,i)*f(10)-ar(11,i)*f(9)
     2          -ar(2,i)*f(7)+ar(10,i)*f(6)+ar(5,i)*f(2)
          rne(6)=ar(1,i)*f(14)+ar(13,i)*f(13)-ar(2,i)*t1
     1          -ar(3,i)*f(12)+ar(14,i)*f(10)-ar(4,i)*f(9)
     2          -ar(11,i)*f(6)-ar(12,i)*f(5)+ar(5,i)*f(3)
        else
          rne(1)=-ar(1,i)*f(3)+ar(2,i)*f(2)-ar(3,i)*f(1)
          rne(2)=-ar(1,i)*f(4)+ar(4,i)*f(2)+ar(2,i)*f(1)
          rne(3)=-ar(2,i)*f(4)+ar(4,i)*f(3)-ar(5,i)*f(1)
          rne(4)=-ar(3,i)*f(4)-ar(2,i)*f(3)-ar(5,i)*f(2)
        end if
        do jj=1,6
          ar(jj,i)=rne(jj)
        enddo
      else
        inorm(i)=iexp
        do j=1,nvesm
          ar(j,i)=f(j)
        enddo
      end if
      if(i.eq.jl) return
      i=i+jud
      go to 10
      end

      subroutine fprpmn(jf,jl,f,h,nvefm,iexp)
c*** propagate the minor vector in a fluid region from level jf to jl ***
      implicit real*8(a-h,o-z)
      parameter (mk=30000)
      save
      real*8 lcon,ncon,lspl,nspl
      REAL*8 :: l
      common r(mk),fmu(mk),flam(mk),qshear(mk),qkappa(mk),
     + xa2(mk),xlam(mk),rho(mk),qro(3,mk),g(mk),qg(3,mk),
     + fcon(mk),fspl(3,mk),lcon(mk),lspl(3,mk),ncon(mk),
     + nspl(3,mk),ccon(mk),cspl(3,mk),acon(mk),aspl(3,mk)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      common/eifx/ar(14,mk),inorm(mk),jjj(mk)
      common/shanks/b(46),c(10),dx,step(8),stepf,maxo,in
      dimension f(*),h(nvefm,1),s(5),fp(5)
      data econst/1048576.d0/
      if(nvefm.eq.1) then
        do i=jl,jf
          inorm(i)=inorm(i)+iexp
          do j=1,2
            ar(j,i)=ar(j,i)*f(1)
          enddo
        enddo
        return
      end if
      maxo1=maxo-1
      jud=1
      if(jl.lt.jf) jud=-1
      y=r(jf)
      i=jf
      go to 45
   10 x=y
      y=r(i)
      if(y.eq.x) goto 45
      iq=min(i,i-jud)
      qff=1.d0+xlam(iq)*fct
      xi=g(i)/y
      alfsq=(wsq+4.d0*rho(i)+xi-fl3*xi*xi/wsq)*rho(i)/flam(i)
      q=(sqrt(abs(alfsq-fl3/(x*x)))+1.d0/r(iq)+float(kg)*sfl3/x)/stepf
      del=float(jud)*step(maxo)/q
      dxs=0.d0
   15   y=x+del
        if(float(jud)*(y-r(i)).gt.0.d0) y=r(i)
        dx=y-x
        if(dx.ne.dxs) call baylis(q,maxo1)
        dxs=dx
        do j=1,nvefm
          s(j)=f(j)
        enddo
        do ni=1,in
          z=x+c(ni)
          call dermf(iq,z,f,h(1,ni),0,qff)
          call rkdot(f,s,h,nvefm,ni)
        enddo
        if(knsw.eq.1) then
          call dermf(iq,y,f,fp,1,qff)
          call zknt(s,h,f,fp,x,y,0)
        end if
        x=y
        if(y.ne.r(i)) go to 15
   45 size=dabs(f(1))
      do j=2,nvefm
        size=dmax1(size,dabs(f(j)))
      enddo
   55 if(size.lt.1024.d0) goto 65
      do j=1,nvefm
        f(j)=f(j)/econst
      enddo
      size=size/econst
      iexp=iexp+20
      goto 55
   65 if(iback.ne.0) then
        inorm(i)=inorm(i)+iexp
        rne2   =-ar(1,i)*f(4)+ar(4,i)*f(2)+ar(2,i)*f(1)
        ar(1,i)=-ar(1,i)*f(3)+ar(2,i)*f(2)-ar(3,i)*f(1)
        rne3   =-ar(2,i)*f(4)+ar(4,i)*f(3)-ar(5,i)*f(1)
        ar(4,i)=-ar(3,i)*f(4)-ar(2,i)*f(3)-ar(5,i)*f(2)
        ar(2,i)=rne2
        ar(3,i)=rne3
      else
        inorm(i)=iexp
        do j=1,nvefm
          ar(j,i)=f(j)
        enddo
      end if
      if(i.eq.jl) return
      i=i+jud
      go to 10
      end

      subroutine derms(iq,z,f,fp,iknt,qff,qll,qaa)
c*** calculates minor vector derivative (fp) in a solid ***
      implicit real*8(a-h,o-z)
      parameter (mk=30000)
      save
      real*8 nn,ll,lcon,ncon,lspl,nspl
      REAL*8 :: l
      common r(mk),fmu(mk),flam(mk),qshear(mk),qkappa(mk),
     + xa2(mk),xlam(mk),rho(mk),qro(3,mk),g(mk),qg(3,mk),
     + fcon(mk),fspl(3,mk),lcon(mk),lspl(3,mk),ncon(mk),
     + nspl(3,mk),ccon(mk),cspl(3,mk),acon(mk),aspl(3,mk)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      dimension f(*),fp(*)
      if(iknt.ne.0) then
        if(kg.eq.0) then
          fp(3)=s22*f(1)-2.d0*t12*f(2)+b33*f(3)+c11*f(5)
          fp(4)=-s11*f(1)+2.d0*t21*f(2)-b33*f(4)-c22*f(5)
          fp(5)=-2.d0*s12*f(2)+s11*f(3)-s22*f(4)-b11*f(5)
        else
          fp(3)=s22*f(1)+b32*f(2)+b33*f(3)+c11*f(5)+b313*f(13)
          fp(4)=-s11*f(1)+b42*f(2)+b44*f(4)-c22*f(5)+b414*f(14)
          fp(5)=b52*f(2)+s11*f(3)-s22*f(4)+b55*f(5)-b313*f(11)
     +      +b414*f(12)
        end if
        return
      end if
      t=z-r(iq)
      if(t.eq.0.d0) then
        ro=rho(iq)
        gr=g(iq)
        ff=fcon(iq)*qff
        ll=lcon(iq)*qll
        nn=ncon(iq)*qll
        cc=ccon(iq)*qaa
        aa=acon(iq)*qaa
      else
        ro=rho(iq)+t*(qro(1,iq)+t*(qro(2,iq)+t*qro(3,iq)))
        gr=g(iq)+t*(qg(1,iq)+t*(qg(2,iq)+t*qg(3,iq)))
        ff=(fcon(iq)+t*(fspl(1,iq)+t*(fspl(2,iq)+t*fspl(3,iq))))*qff
        ll=(lcon(iq)+t*(lspl(1,iq)+t*(lspl(2,iq)+t*lspl(3,iq))))*qll
        if(ifanis.eq.0) then
          nn=ll
          cc=ff+ll+ll
          aa=cc
        else
          nn=(ncon(iq)+t*(nspl(1,iq)+t*(nspl(2,iq)+t*nspl(3,iq))))*qll
          cc=(ccon(iq)+t*(cspl(1,iq)+t*(cspl(2,iq)+t*cspl(3,iq))))*qaa
          aa=(acon(iq)+t*(aspl(1,iq)+t*(aspl(2,iq)+t*aspl(3,iq))))*qaa
        end if
      end if
      zr=1.d0/z
      sfl3z=sfl3*zr
      rogr=ro*gr
      c11=1.d0/cc
      c22=1.d0/ll
      dmg=aa-nn-ff*ff*c11
      zdmg=zr*dmg
      t11=-2.d0*ff*zr*c11+zr
      t12=sfl3z*ff*c11
      t21=-sfl3z
      t22=zr+zr
      s22=-ro*wsq
      s11=s22+4.d0*zr*(zdmg-rogr)
      s22=s22+zr*zr*(fl3*(dmg+nn)-nn-nn)
      s12=sfl3z*(rogr-zdmg-zdmg)
      if(kg.eq.0) then
        s11=s11+4.d0*ro*ro
        if(iback.eq.0) then
          b11=t11+t22
          b33=t11-t22
          fp(1)=b11*f(1)+c22*f(3)-c11*f(4)
          fp(2)=s12*f(1)-t21*f(3)+t12*f(4)
          fp(3)=s22*f(1)-2.d0*t12*f(2)+b33*f(3)+c11*f(5)
          fp(4)=-s11*f(1)+2.d0*t21*f(2)-b33*f(4)-c22*f(5)
          fp(5)=-2.d0*s12*f(2)+s11*f(3)-s22*f(4)-b11*f(5)
        else
          fp(1)=t22*f(1)-t21*f(2)-c22*f(3)
          fp(2)=-t12*f(1)+t11*f(2)-c11*f(4)
          fp(3)=-s22*f(1)+s12*f(2)-t22*f(3)+t12*f(4)
          fp(4)=s12*f(1)-s11*f(2)+t21*f(3)-t11*f(4)
        endif
      else
        t31=-4.d0*ro
        t33=-fl*zr
        s13=-fl1*zr*ro
        s23=ro*sfl3z
        if(iback.eq.0) then
          b11=t11+t22-t33
          b33=t11-t22-t33
          b44=t22-t11-t33
          b55=-t11-t22-t33
          b32=-t12-t12
          b42=t21+t21
          b52=-s12-s12
          b313=-s23-s23
          b414=s13+s13
          b914=t31+t31
          fp(1)=b11*f(1)+c22*f(3)-c11*f(4)
          fp(2)=s12*f(1)-t33*f(2)-t21*f(3)+t12*f(4)-s13*f(13)-s23*f(14)
          fp(3)=s22*f(1)+b32*f(2)+b33*f(3)+c11*f(5)+b313*f(13)
          fp(4)=-s11*f(1)+b42*f(2)+b44*f(4)-c22*f(5)+b414*f(14)
          fp(5)=b52*f(2)+s11*f(3)-s22*f(4)+b55*f(5)-b313*f(11)
     +        +b414*f(12)
          fp(6)=4.d0*f(1)-b55*f(6)+c22*f(8)-c11*f(9)
          fp(7)=4.d0*f(2)+s12*f(6)+t33*f(7)-t21*f(8)+t12*f(9)-t31*f(13)
          fp(8)=4.d0*f(3)+s22*f(6)+b32*f(7)-b44*f(8)+c11*f(10)
          fp(9)=4.d0*f(4)-s11*f(6)+b42*f(7)-b33*f(9)-c22*f(10)
     +         +b914*f(14)
          fp(10)=4.d0*f(5)+b52*f(7)+s11*f(8)-s22*f(9)-b11*f(10)
     +          +b914*f(12)
          fp(11)=-t31*f(2)+s13*f(7)+s23*f(9)-t11*f(11)+t21*f(12)
     +      -s11*f(13)+s12*f(14)
          fp(12)=t31*f(3)+s23*f(7)-s13*f(8)+t12*f(11)-t22*f(12)
     +      +s12*f(13)-s22*f(14)
          fp(13)=s23*f(6)-c11*f(11)+t11*f(13)-t12*f(14)
          fp(14)=-t31*f(1)+s13*f(6)-c22*f(12)-t21*f(13)+t22*f(14)
        else
          b11=t22+t33
          b22=t11+t33
          b33=t11+t22
          b55=t22-t33
          b66=t11-t33
          b99=t11-t22
          t4=f(4)+f(8)
          t5=t4+f(8)
          t4=t4+f(4)
          fp(1)=b11*f(1)-t21*f(2)-t31*f(3)-4.d0*f(5)+c22*f(7)
          fp(2)=-t12*f(1)+b22*f(2)-4.d0*f(6)+c11*f(11)
          fp(3)=b33*f(3)-c22*f(9)+c11*f(12)
          fp(4)=-s23*f(1)+s13*f(2)+t31*f(6)
          fp(5)=s13*f(3)+b55*f(5)-t21*f(6)-c22*f(10)
          fp(6)=s23*f(3)-t12*f(5)+b66*f(6)-c11*f(13)
          fp(7)=s22*f(1)-s12*f(2)-b55*f(7)+t31*f(9)+4.d0*f(10)+t12*f(11)
          fp(8)=s23*f(1)-s12*f(3)-t21*f(9)+t12*f(12)
          fp(9)=s23*f(2)-s22*f(3)-t12*t5+b99*f(9)-c11*f(14)
          fp(10)=s23*(f(4)-f(8))-s22*f(5)+s12*f(6)+s13*f(9)-b11*f(10)
     1      +t12*f(13)
          fp(11)=-s12*f(1)+s11*f(2)-t4*t31+t21*f(7)-b66*f(11)+4.d0*f(13)
          fp(12)=-s13*f(1)+s11*f(3)+t21*t5-t31*f(5)-b99*f(12)+c22*f(14)
          fp(13)=-t4*s13+s12*f(5)-s11*f(6)+t21*f(10)-s23*f(12)-b22*f(13)
          fp(14)=s12*t5-s13*f(7)-s11*f(9)+t31*f(10)-s23*f(11)+s22*f(12)
     1      -b33*f(14)
        end if
      end if
      return
      end

      subroutine dermf(iq,z,f,fp,iknt,qff)
c*** calculates minor vector derivative (fp) in a fluid ***
      implicit real*8(a-h,o-z)
      parameter (mk=30000)
      save
      real*8 lcon,ncon,lspl,nspl
      REAL*8 :: l
      common r(mk),fmu(mk),flam(mk),qshear(mk),qkappa(mk),
     + xa2(mk),xlam(mk),rho(mk),qro(3,mk),g(mk),qg(3,mk),
     + fcon(mk),fspl(3,mk),lcon(mk),lspl(3,mk),ncon(mk),
     + nspl(3,mk),ccon(mk),cspl(3,mk),acon(mk),aspl(3,mk)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      dimension f(*),fp(*)
      if(iknt.ne.0) then
        if(kg.eq.0) then
          fp(1)=t11*f(1)+c11*f(2)
          fp(2)=(s11-t21*ro)*f(1)-t11*f(2)
        else
          fp(3)=s22*f(1)-(t12+t12)*f(2)+b33*f(3)+c11*f(5)
          fp(4)=-s11*f(1)+(t21+t21)*f(2)-b33*f(4)-4.d0*f(5)
          fp(5)=-(s12+s12)*f(2)+s11*f(3)-s22*f(4)-b11*f(5)
        end if
        return
      end if
      t=z-r(iq)
      if(t.eq.0.d0) then
        ro=rho(iq)
        flu=fcon(iq)*qff
        gr=g(iq)
      else
        ro=rho(iq)+t*(qro(1,iq)+t*(qro(2,iq)+t*qro(3,iq)))
        flu=(fcon(iq)+t*(fspl(1,iq)+t*(fspl(2,iq)+t*fspl(3,iq))))*qff
        gr=g(iq)+t*(qg(1,iq)+t*(qg(2,iq)+t*qg(3,iq)))
      end if
      t21=-4.d0*ro
      zr=1.d0/z
      t12=fl3*zr*zr/wsq
      t11=gr*t12-zr
      s11=ro*(gr*gr*t12-wsq)+t21*gr*zr
      c11=-t12/ro+1.d0/flu
      if(kg.eq.0) then
        fp(1)=t11*f(1)+c11*f(2)
        fp(2)=(s11-t21*ro)*f(1)-t11*f(2)
      else
        t22=-fl*zr
        s22=ro*t12
        b11=t11+t22
        s12=ro*b11
        if(iback.eq.0) then
          b33=t11-t22
          fp(1)=b11*f(1)+4.d0*f(3)-c11*f(4)
          fp(2)=s12*f(1)-t21*f(3)+t12*f(4)
          fp(3)=s22*f(1)-(t12+t12)*f(2)+b33*f(3)+c11*f(5)
          fp(4)=-s11*f(1)+(t21+t21)*f(2)-b33*f(4)-4.d0*f(5)
          fp(5)=-(s12+s12)*f(2)+s11*f(3)-s22*f(4)-b11*f(5)
        else
          fp(1)=t22*f(1)-t21*f(2)-4.d0*f(3)
          fp(2)=-t12*f(1)+t11*f(2)-c11*f(4)
          fp(3)=-s22*f(1)+s12*f(2)-t22*f(3)+t12*f(4)
          fp(4)=s12*f(1)-s11*f(2)+t21*f(3)-t11*f(4)
        end if
      end if
      return
      end

      subroutine eifout(lsmin)
c*** massages spheroidal mode eigenfunctions before output ***
      implicit real*8(a-h,o-z)
      parameter (mk=30000)
      save
      real*8 ll,lcon,ncon,lspl,nspl
      REAL*8 :: l
      common r(mk),fmu(mk),flam(mk),qshear(mk),qkappa(mk),
     + xa2(mk),xlam(mk),rho(mk),qro(3,mk),g(mk),qg(3,mk),
     + fcon(mk),fspl(3,mk),lcon(mk),lspl(3,mk),ncon(mk),
     + nspl(3,mk),ccon(mk),cspl(3,mk),acon(mk),aspl(3,mk)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      common/eifx/a(14,mk),inorm(mk),jjj(mk)
      common/rindx/nic,noc,nsl,nicp1,nocp1,nslp1,n
      dimension zi(4)
      i1=min0(nic,max0(2,lsmin))
      i2=nic
    5 if(i1.eq.i2) goto 20
      do iq=i1,i2
        ff=fcon(iq)*(1.d0+xlam(iq)*fct)
        ll=lcon(iq)*(1.d0+qshear(iq)*fct)
        zr=1.d0/r(iq)
        sfl3z=sfl3*zr
        d=1.d0/(ccon(iq)*(1.d0+xa2(iq)*fct))
        v=a(2,iq)
        if(kg.eq.0) then
          a(2,iq)=(zr-2.d0*ff*d*zr)*a(1,iq)+sfl3z*ff*d*v+d*a(3,iq)
          a(4,iq)=-sfl3z*a(1,iq)+(zr+zr)*v+a(4,iq)/ll
          a(5,iq)=0.d0
          a(6,iq)=0.d0
        else
          a(2,iq)=(zr-2.d0*ff*d*zr)*a(1,iq)+sfl3z*ff*d*v+d*a(4,iq)
          a(4,iq)=-sfl3z*a(1,iq)+(zr+zr)*v+a(5,iq)/ll
          a(5,iq)=a(3,iq)
          a(6,iq)=4.d0*(a(6,iq)-rho(iq)*a(1,iq))-fl*zr*a(5,iq)
        end if
        a(3,iq)=v
      enddo
   20 if(i2.eq.nsl) goto 25
      i1=min(nsl,max(lsmin,nocp1))
      i2=nsl
      goto 5
   25 i1=min(noc,max(lsmin,nicp1))
      i2=noc
   30 if(i1.eq.i2) goto 50
      do iq=i1,i2
        zr=1.d0/r(iq)
        sfl3z=sfl3*zr
        ffi=1.d0/(flam(iq)*(1.d0+xlam(iq)*fct))
        if(kg.eq.0) then
          p=a(2,iq)
          a(5,iq)=0.d0
          a(6,iq)=0.d0
        else
          p=a(3,iq)
          a(5,iq)=a(2,iq)
          a(6,iq)=4.d0*(a(4,iq)-rho(iq)*a(1,iq))-fl*zr*a(5,iq)
        end if
        a(3,iq)=sfl3z*(g(iq)*a(1,iq)-p/rho(iq)+a(5,iq))/wsq
        a(2,iq)=sfl3z*a(3,iq)-a(1,iq)*zr+p*ffi
        a(4,iq)=sfl3z*(a(1,iq)+p*(qro(1,iq)/(rho(iq)**2)+g(iq)*ffi)/wsq)
      enddo
   50 if(n.eq.nsl.or.i2.eq.n) goto 55
      i1=nslp1
      i2=n
      goto 30
   55 imax=0
      do iq=lsmin,n
        imax=max(inorm(iq),imax)
      enddo
      do iq=lsmin,n
        iexp=inorm(iq)-imax
        al=0.d0
        if(iexp.ge.-80) al=2.d0**iexp
        do j=1,6
          a(j,iq)=a(j,iq)*al
        enddo
      enddo
      lsm1=max(1,lsmin-1)
      do i=1,lsm1
        do j=1,6
          a(j,i)=0.d0
        enddo
      enddo
      if(l.gt.1.or.lsmin.gt.2) goto 75
      a(2,1)=1.5d0*a(1,2)/r(2)-.5d0*a(2,2)
      a(4,1)=1.5d0*a(3,2)/r(2)-.5d0*a(4,2)
   75 do j=1,4
        zi(j)=0.d0
      enddo
      i1=max0(lsmin,2)
      do iq=i1,n
        ip=iq-1
        if(r(iq).ne.r(ip)) call gauslv(r(ip),r(iq),ip,zi,4)
      enddo
      cg=zi(2)/(w*zi(1))
      wray=dsqrt(2.d0*zi(4)/zi(1))
      qinv=2.d0*zi(3)/(wsq*zi(1))
      rnorm=1.d0/(w*dsqrt(zi(1)))
      do iq=i1,n
        zr=1.d0/r(iq)
        a(1,iq)=a(1,iq)*zr
        a(2,iq)=(a(2,iq)-a(1,iq))*zr
        a(3,iq)=a(3,iq)*zr
        a(4,iq)=(a(4,iq)-a(3,iq))*zr
        a(5,iq)=a(5,iq)*zr
        a(6,iq)=(a(6,iq)-a(5,iq))*zr
        a(1,iq)=a(1,iq)*rnorm
        a(2,iq)=a(2,iq)*rnorm
        a(3,iq)=a(3,iq)*rnorm
        a(4,iq)=a(4,iq)*rnorm
        a(5,iq)=a(5,iq)*rnorm
        a(6,iq)=a(6,iq)*rnorm
      enddo
      if(lsmin.gt.2.or.l.gt.2) return
      if(l.eq.1) then
        a(1,1)=a(1,2)-.5d0*a(2,2)*r(2)
        a(2,1)=0.d0
        a(3,1)=a(3,2)-.5d0*a(4,2)*r(2)
        a(4,1)=0.d0
        a(6,1)=1.5d0*a(5,2)/r(2)-.5d0*a(6,2)
      else
        a(2,1)=1.5d0*a(1,2)/r(2)-.5d0*a(2,2)
        a(4,1)=1.5d0*a(3,2)/r(2)-.5d0*a(4,2)
      end if
      return
      end

      subroutine gauslv(r1,r2,iq,fint,nint)
c*** fifth order gauss-legendre integration ***
      implicit real*8(a-h,o-z)
      save
      dimension fint(*),vals(4),vals1(4),sum(4),w(2),x(2)
      data w,x/.478628670499366d0,.236926885056189d0,
     +         .538469310105683d0,.906179845938664d0/
      y1=.5d0*(r2+r1)
      y2=.5d0*(r2-r1)
      call intgds(y1,iq,vals)
      do j=1,nint
        sum(j)=.568888888888889d0*vals(j)
      enddo
      do i=1,2
        t1=x(i)*y2
        call intgds(y1+t1,iq,vals)
        call intgds(y1-t1,iq,vals1)
        do j=1,nint
          sum(j)=sum(j)+w(i)*(vals(j)+vals1(j))
        enddo
      enddo
      do j=1,nint
        fint(j)=fint(j)+y2*sum(j)
      enddo
      return
      end

      subroutine sdepth(wdim,ls)
c*** finds starting level,ls, for a given l and w ***
      implicit real*8(a-h,o-z)
      parameter (mk=30000)
      save
      real*8 lcon,ncon,lspl,nspl
      REAL*8 :: l
      common r(mk),fmu(mk),flam(mk),qshear(mk),qkappa(mk),
     + xa2(mk),xlam(mk),rho(mk),qro(3,mk),g(mk),qg(3,mk),
     + fcon(mk),fspl(3,mk),lcon(mk),lspl(3,mk),ncon(mk),
     + nspl(3,mk),ccon(mk),cspl(3,mk),acon(mk),aspl(3,mk)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      common/rindx/nic,noc,nsl,nicp1,nocp1,nslp1,n
      data aw,bw,dw/-2.d-3,2.25d-3,1.28d-3/
      q=0.d0
      w=wdim/wn
      wsoc=aw+dw*fl
      if(wdim.gt.wsoc) goto 10
      call startl(nocp1,nsl,fmu,ls,q)
      if(ls.eq.nsl) ls=ls-1
      if(ls.gt.nocp1) return
   10 wsic=aw+bw*fl
      if(wdim.gt.wsic) goto 20
      call startl(nicp1,noc,flam,ls,q)
      if(ls.eq.noc) ls=ls-1
      if(ls.gt.nicp1) return
   20 call startl(2,nic,fmu,ls,q)
      if(ls.eq.nic) ls=ls-1
      return
      end

      subroutine sfbm(ass,kg,iback)
c*** convert minor vector at a solid/fluid boundary ***
      implicit real*8(a-h,o-z)
      save
      dimension ass(14),as(14)
      do j=1,14
        as(j)=ass(j)
        ass(j)=0.d0
      enddo
      if(iback.ne.1) then
        if(kg.eq.0) then
          ass(1)=as(3)
          ass(2)=as(5)
        else
          ass(1)=as(8)
          ass(2)=-as(12)
          ass(3)=as(3)
          ass(4)=-as(10)
          ass(5)=as(5)
        end if
      else
        if(kg.eq.0) then
          ass(1)=-as(3)
        else
          ass(1)=as(7)
          ass(2)=-as(9)
          ass(3)=-as(10)
          ass(4)=-as(14)
        end if
      end if
      return
      end

      subroutine fsbm(ass,kg,iback)
c*** convert minor vector at a fluid/solid boundary ***
      implicit real*8(a-h,o-z)
      save
      dimension ass(14),as(14)
      do j=1,14
        as(j)=ass(j)
        ass(j)=0.d0
      end do
      if(iback.ne.1) then
        if(kg.eq.0) then
          ass(1)=as(1)
          ass(4)=-as(2)
        else
          ass(6)=as(1)
          ass(14)=as(2)
          ass(1)=as(3)
          ass(9)=as(4)
          ass(4)=-as(5)
        end if
      else
        if(kg.eq.0) then
          ass(1)=-as(1)
        else
          ass(1)=-as(1)
          ass(3)=-as(2)
          ass(5)=-as(3)
          ass(12)=as(4)
        end if
      end if
      return
      end

      subroutine zknt(s,sp,f,fp,x,y,ifsol)
c*** given minor vector and derivs,constructs mode count ***
      implicit real*8(a-h,o-z)
      save
      REAL*8 :: l
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      dimension s(*),sp(*),f(*),fp(*),xs(4),val(4)
      if(ifsol.eq.0.and.kg.eq.0) goto 5
      y1=s(5)
      y2=f(5)
      y1p=sp(5)
      y2p=fp(5)
      t1=s(3)-s(4)
      t2=f(3)-f(4)
      t1p=sp(3)-sp(4)
      t2p=fp(3)-fp(4)
      goto 10
    5 y1=s(2)
      y2=f(2)
      y1p=sp(2)
      y2p=fp(2)
      t1=s(1)
      t2=f(1)
      t1p=sp(1)
      t2p=fp(1)
   10 h=y-x
      ns=0
      if(kount.ne.0) goto 15
      a1=y2-y1
      a2=0.d0
      a3=0.d0
      a22=0.d0
      a33=0.d0
      goto 50
   15 a1=h*y1p
      a2=-h*(2.d0*y1p+y2p)+3.d0*(y2-y1)
      a3=h*(y1p+y2p)-2.d0*(y2-y1)
      a33=3.d0*a3
      a22=2.d0*a2
      if(a3.ne.0.d0) goto 20
      if(a2.eq.0.d0) goto 50
      xs(2)=-a1/a22
      if(xs(2).ge.0.d0.and.xs(2).le.1.d0) ns=1
      goto 50
   20 disc=a2*a2-a1*a33
      if(disc) 50,25,30
   25 xs(2)=-a2/a33
      if(xs(2).ge.0.d0.and.xs(2).le.1.d0) ns=1
      goto 50
   30 disc=dsqrt(disc)
      tr1=(-a2+disc)/a33
      tr2=(-a2-disc)/a33
      if(dabs(a33).gt.dabs(a1)) goto 35
      fac=a1/a33
      tr1=fac/tr1
      tr2=fac/tr2
   35 if(tr1.lt.0.d0.or.tr1.gt.1.d0) goto 40
      xs(2)=tr1
      ns=1
   40 if(tr2.lt.0.d0.or.tr2.gt.1.d0) goto 50
      ns=ns+1
      xs(ns+1)=tr2
      if(ns.lt.2) goto 50
      if(tr2.ge.tr1) goto 50
      xs(2)=tr2
      xs(3)=tr1
   50 val(1)=y1
      xs(1)=0.d0
      ns2=ns+2
      val(ns2)=y2
      xs(ns2)=1.d0
      if(ns.eq.0) goto 60
      ns1=ns+1
      do 55 j=2,ns1
      t=xs(j)
   55 val(j)=y1+t*(a1+t*(a2+t*a3))
   60 ift=0
      do 100 j=2,ns2
      if(val(j-1)*val(j).gt.0.d0) goto 100
      if(val(j-1).ne.0.d0) goto 65
      tes=t1*a1
      goto 90
   65 rt1=0.5d0*(xs(j-1)+xs(j))
      rt=rt1
      do 70 i=1,5
      v=y1+rt*(a1+rt*(a2+rt*a3))
      vp=a1+rt*(a22+rt*a33)
      add=-v/vp
      rt=rt+add
      if(dabs(add).lt.1.d-5) goto 75
      if(dabs(rt-rt1).le..5d0) goto 70
      rt=rt1
      goto 75
   70 continue
   75 if(ift.ne.0) goto 85
      if(kount.ne.0) goto 80
      b1=t2-t1
      b2=0.d0
      b3=0.d0
      goto 85
   80 b1=h*t1p
      b2=-h*(2.d0*t1p+t2p)+3.d0*(t2-t1)
      b3=h*(t1p+t2p)-2.d0*(t2-t1)
      ift=1
   85 tes=t1+rt*(b1+rt*(b2+rt*b3))
      vp=a1+rt*(a22+rt*a33)
      tes=tes*vp
   90 if(tes.lt.0.d0) kount=1+kount
      if(tes.gt.0.d0) kount=kount-1
  100 continue
      return
      end

      subroutine baylis(q,maxo1)
c    baylis returns the coefficients for rks integration.
c    see e. baylis shanks(1966 a. m. s.) and references therein for the
c    coefficients. the eight runge-kutta-shanks formulae are (1-1) (2-2)
c    (3-3) (4-4) (5-5) (6-6) (7-7) (8-10). for orders greater than 4 the
c    formulae are approximate rather than exact so incurring less roundoff.
      implicit real*8(a-h,o-z)
      save
      common/shanks/b(46),c(10),dx,step(8),stepf,maxo,i
      ds=q*dabs(dx)
      do 10 j=1,maxo1
      if(ds.gt.step(j)) go to 10
      i=j
      go to 15
   10 continue
      i=maxo
   15 c(1)=0.d0
      go to (1,2,3,4,5,6,7,8),i
    1 b(1)=dx
      return
    2 c(2)=dx
      b(1)=dx
      b(2)=.5d0*dx
      b(3)=1.d0
      return
    3 c(2)=.5d0*dx
      c(3)=dx
      b(1)=c(2)
      b(2)=-dx
      b(3)=-2.d0
      b(4)=.16666666666667d0*dx
      b(5)=4.d0
      b(6)=1.d0
      return
    4 c(2)=.01d0*dx
      c(3)=.6d0*dx
      c(4)=dx
      b(1)=c(2)
      b( 2)=-.17461224489790d+02*dx
      b( 3)=-.10343618513324d+01
      b( 4)= .59691275167780d+02*dx
      b( 5)=-.10140620414448d+01
      b( 6)= .30814908546230d-01
      b( 7)=-.25555555555556d+01*dx
      b( 8)=-.11165449632656d+01
      b( 9)=-.22568165070006d+00
      b(10)=-.49077733860351d-01
      return
    5 c( 2)= 1.1111111111111d-04*dx
      c( 3)= 3.0d-01*dx
      c( 4)= 7.5d-01*dx
      c( 5)= dx
      b( 1)=c(2)
      b( 2)=-.40470000000000d+03*dx
      b( 3)=-.10007412898443d+01
      b( 4)= .25301250000000d+04*dx
      b( 5)=-.10004446420631d+01
      b( 6)= .74107010523195d-03
      b( 7)=-.11494333333333d+05*dx
      b( 8)=-.10004929965491d+01
      b( 9)= .52629261224803d-03
      b(10)=-.12029545422812d-03
      b(11)= .92592592592593d-01*dx
      b(12)= .00000000000000d+00
      b(13)= .47619047619048d+01
      b(14)= .42666666666667d+01
      b(15)= .77142857142857d+00
      return
    6 c(2)=3.3333333333333d-03*dx
      c(3)=.2d0*dx
      c(4)=.6d0*dx
      c(5)=9.3333333333333d-01*dx
      c(6)=dx
      b( 1)=c(2)
      b( 2)=-.58000000000000d+01*dx
      b( 3)=-.10344827586207d+01
      b( 4)= .64600000000000d+02*dx
      b( 5)=-.10216718266254d+01
      b( 6)= .30959752321982d-01
      b( 7)=-.62975802469136d+03*dx
      b( 8)=-.10226149961576d+01
      b( 9)= .24906685695466d-01
      b(10)=-.37737402568887d-02
      b(11)=-.54275714285714d+04*dx
      b(12)=-.10225567867765d+01
      b(13)= .25375487829097d-01
      b(14)=-.31321559234596d-02
      b(15)= .12921040478749d-03
      b(16)= .53571428571429d-01*dx
      b(17)= .00000000000000d+00
      b(18)= .61868686868687d+01
      b(19)= .77777777777778d+01
      b(20)= .40909090909091d+01
      b(21)=-.38888888888889d+00
      return
    7 c(2)=5.2083333333333d-03*dx
      c(3)=1.6666666666667d-01*dx
      c(4)=.5d0*dx
      c(5)=dx
      c(6)=8.3333333333333d-01*dx
      c(7)=dx
      b( 1)=c(2)
      b( 2)=-.25000000000000d+01*dx
      b( 3)=-.10666666666667d+01
      b( 4)= .26166666666667d+02*dx
      b( 5)=-.10421204027121d+01
      b( 6)= .61228682966918d-01
      b( 7)=-.64500000000000d+03*dx
      b( 8)=-.10450612653163d+01
      b( 9)= .51262815703925d-01
      b(10)=-.77519379844961d-02
      b(11)=-.93549382716049d+02*dx
      b(12)=-.10450293206756d+01
      b(13)= .48394546673620d-01
      b(14)=-.11877268228307d-01
      b(15)=-.39590894094358d-03
      b(16)= .35111904761905d+03*dx
      b(17)=-.10446476812124d+01
      b(18)= .52479782656724d-01
      b(19)=-.71200922221468d-02
      b(20)=-.61029361904114d-03
      b(21)= .27463212856852d-02
      b(22)= .46666666666667d-01*dx
      b(23)= .57857142857143d+01
      b(24)= .78571428571429d+01
      b(25)= .00000000000000d+00
      b(26)= b(23)
      b(27)= .10000000000000d+01
      return
    8 c(2)=.14814814814815d0*dx
      c(3)=.22222222222222d0*dx
      c(4)=.33333333333333d0*dx
      c(5)= .5d0*dx
      c(6)=.66666666666667d0*dx
      c(7)=.16666666666667d0*dx
      c(8)=dx
      c(9)=.83333333333333d0*dx
      c(10)=dx
      b( 1)=c(2)
      b( 2)= .55555555555556d-01*dx
      b( 3)= .30000000000000d+01
      b( 4)= .83333333333333d-01*dx
      b( 5)= .00000000000000d+00
      b( 6)= .30000000000000d+01
      b( 7)= .12500000000000d+00*dx
      b( 8)= .00000000000000d+00
      b( 9)= .00000000000000d+00
      b(10)= .30000000000000d+01
      b(11)= .24074074074074d+00*dx
      b(12)= .00000000000000d+00
      b(13)=-.20769230769231d+01
      b(14)= .32307692307692d+01
      b(15)= .61538461538461d+00
      b(16)= .90046296296295d-01*dx
      b(17)= .00000000000000d+00
      b(18)=-.13881748071980d+00
      b(19)= .24832904884319d+01
      b(20)=-.21182519280206d+01
      b(21)= .62467866323908d+00
      b(22)=-.11550000000000d+02*dx
      b(23)=-.35064935064935d+00
      b(24)= .50389610389610d+01
      b(25)=-.28398268398268d+01
      b(26)= .52813852813853d+00
      b(27)=-.34632034632035d+01
      b(28)=-.44097222222222d+00*dx
      b(29)=-.14173228346457d+00
      b(30)= .53385826771654d+01
      b(31)=-.35905511811023d+01
      b(32)= .70866141732284d-01
      b(33)=-.45354330708661d+01
      b(34)=-.31496062992126d-01
      b(35)= .18060975609756d+01*dx
      b(36)=-.54692775151925d-01
      b(37)= .47967589466576d+01
      b(38)=-.22795408507765d+01
      b(39)= .48615800135044d-01
      b(40)=-.34031060094530d+01
      b(41)=-.40513166779204d-01
      b(42)= .48615800135044d+00
      b(43)= .48809523809524d-01*dx
      b(44)= .65853658536585d+00
      b(45)= .66341463414634d+01
      b(46)= .52682926829268d+01
      i=10
      return
      end

      subroutine grav(g,rho,qro,r,n)
c*** given rho and spline coeffs,computes gravity ***
      implicit real*8(a-h,o-z)
      save
      dimension g(1),rho(1),qro(3,1),r(1)
      g(1)=0.d0
      do i=2,n
        im1=i-1
        del=r(i)-r(im1)
        rn2=r(im1)*r(im1)
        trn=2.d0*r(im1)
        c1=rho(im1)*rn2
        c2=(qro(1,im1)*rn2+trn*rho(im1))*0.5d0
        c3=(qro(2,im1)*rn2+trn*qro(1,im1)+rho(im1))/3.d0
        c4=(qro(3,im1)*rn2+trn*qro(2,im1)+qro(1,im1))*.25d0
        c5=(trn*qro(3,im1)+qro(2,im1))*0.2d0
        g(i)=(g(im1)*rn2+4.d0*del*(c1+del*(c2+del*(c3+del*(c4
     +   +del*(c5+del*qro(3,im1)/6.d0))))))/(r(i)*r(i))
      enddo
      return
      end

      subroutine startl(jf,jl,v,ls,q)
c*** finds start level between jf and jl using velocityv and ang. ord. l.
c*** upon entry q is the value of the exponent at r(jf) or at the turning
c*** point(q=0) depending on previous calls to startl. upon exit q is the
c*** value of the exponent at the starting level ls.
      implicit real*8(a-h,o-z)
      parameter (mk=30000)
      save
      real*8 lcon,ncon,lspl,nspl
      REAL*8 :: l
      common r(mk),fmu(mk),flam(mk),qshear(mk),qkappa(mk),
     + xa2(mk),xlam(mk),rho(mk),qro(3,mk),g(mk),qg(3,mk),
     + fcon(mk),fspl(3,mk),lcon(mk),lspl(3,mk),ncon(mk),
     + nspl(3,mk),ccon(mk),cspl(3,mk),acon(mk),aspl(3,mk)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      common/rindx/nic,noc,nsl,nicp1,nocp1,nslp1,n
      dimension rrlog(mk),p(mk),v(*)
      data ifirst/1/
      if(ifirst.eq.1) then
        ifirst=0
        vertno=-dlog(eps)
        do i=3,n
          rrlog(i)=.5d0*dlog(r(i)/r(i-1))
        enddo
      end if
      do j=jf,jl
        pp=fl3-wsq*r(j)*r(j)*rho(j)/v(j)
        if(pp.le.0.d0) goto 15
        p(j)=dsqrt(pp)
      enddo
   15 p(j)=0.d0
   20 k=j
      j=j-1
      if(j.le.jf) go to 25
      q=q+rrlog(k)*(p(j)+p(k))
      if(q.lt.vertno) go to 20
      ls=j
      return
   25 ls=jf
      return
      end

      subroutine steps(eps)
c*** computes 8 dimensionless step sizes for rks integration
      implicit real*8(a-h,o-z)
      save
      common/shanks/b(46),c(10),dx,step(8),stepf,maxo,in
      ps=dlog(eps)
      fac=1.d0
      do n=1,8
        fn=n+1
        fac=fac*fn
        x=(dlog(fac)+ps)/fn
        x=exp(x)
        s=x
        do i=1,n
          s=x*dexp(-s/fn)
        enddo
        step(n)=s
      enddo
      return
      end

      subroutine drspln(i1,i2,x,y,q,f)
      implicit real*8(a-h,o-z)
      save
c   rspln computes cubic spline interpolation coefficients
c   for y(x) between grid points i1 and i2 saving them in q.  the
c   interpolation is continuous with continuous first and second
c   derivitives.  it agrees exactly with y at grid points and with the
c   three point first derivitives at both end points (i1 and i2).
c   x must be monotonic but if two successive values of x are equal
c   a discontinuity is assumed and seperate interpolation is done on
c   each strictly monotonic segment.  the arrays must be dimensioned at
c   least - x(i2), y(i2), q(3,i2), and f(3,i2).  f is working storage
c   for rspln.
c                                                     -rpb
      dimension x(1),y(1),q(3,1),f(3,1),yy(3)
      equivalence (yy(1),y0)
      data yy/3*0.d0/
      j1=i1+1
      y0=0.d0
c   bail out if there are less than two points total.
      if(i2-i1)13,17,8
 8    a0=x(j1-1)
c   search for discontinuities.
      do 3 i=j1,i2
      b0=a0
      a0=x(i)
      if(a0-b0)3,4,3
 3    continue
 17   j1=j1-1
      j2=i2-2
      go to 5
 4    j1=j1-1
      j2=i-3
c   see if there are enough points to interpolate (at least three).
 5    if(j2+1-j1)9,10,11
c   only two points.  use linear interpolation.
 10   j2=j2+2
      y0=(y(j2)-y(j1))/(x(j2)-x(j1))
      do 15 j=1,3
      q(j,j1)=yy(j)
 15   q(j,j2)=yy(j)
      go to 12
c   more than two points.  do spline interpolation.
 11   a0=0.d0
      h=x(j1+1)-x(j1)
      h2=x(j1+2)-x(j1)
      y0=h*h2*(h2-h)
      h=h*h
      h2=h2*h2
c   calculate derivitive at near end.
      b0=(y(j1)*(h-h2)+y(j1+1)*h2-y(j1+2)*h)/y0
      b1=b0
c   explicitly reduce banded matrix to an upper banded matrix.
      do 1 i=j1,j2
      h=x(i+1)-x(i)
      y0=y(i+1)-y(i)
      h2=h*h
      ha=h-a0
      h2a=h-2.d0*a0
      h3a=2.d0*h-3.*a0
      h2b=h2*b0
      !print*,'hello',i
      q(1,i)=h2/ha
      q(2,i)=-ha/(h2a*h2)
      q(3,i)=-h*h2a/h3a
      f(1,i)=(y0-h*b0)/(h*ha)
      f(2,i)=(h2b-y0*(2.d0*h-a0))/(h*h2*h2a)
      f(3,i)=-(h2b-3.d0*y0*ha)/(h*h3a)
      a0=q(3,i)
 1    b0=f(3,i)
c   take care of last two rows.
      i=j2+1
      h=x(i+1)-x(i)
      y0=y(i+1)-y(i)
      h2=h*h
      ha=h-a0
      h2a=h*ha
      h2b=h2*b0-y0*(2.d0*h-a0)
      q(1,i)=h2/ha
      f(1,i)=(y0-h*b0)/h2a
      ha=x(j2)-x(i+1)
      y0=-h*ha*(ha+h)
      ha=ha*ha
c   calculate derivitive at far end.
      y0=(y(i+1)*(h2-ha)+y(i)*ha-y(j2)*h2)/y0
      q(3,i)=(y0*h2a+h2b)/(h*h2*(h-2.d0*a0))
      q(2,i)=f(1,i)-q(1,i)*q(3,i)
c   solve upper banded matrix by reverse iteration.
      do 2 j=j1,j2
      k=i-1
      q(1,i)=f(3,k)-q(3,k)*q(2,i)
      q(3,k)=f(2,k)-q(2,k)*q(1,i)
      q(2,k)=f(1,k)-q(1,k)*q(3,k)
 2    i=k
      q(1,i)=b1
c   fill in the last point with a linear extrapolation.
 9    j2=j2+2
      do 14 j=1,3
 14   q(j,j2)=yy(j)
c   see if this discontinuity is the last.
 12   if(j2-i2)6,13,13
c   no.  go back for more.
 6    j1=j2+2
      if(j1-i2)8,8,7
c   there is only one point left after the latest discontinuity.
 7    do 16 j=1,3
 16   q(j,i2)=yy(j)
c   fini.
 13   return
      end

      subroutine dsplin(n,x,y,q,f)
      implicit real*8(a-h,o-z)
      save
      dimension x(*),y(*),q(3,1),f(3,1),yy(3)
      equivalence (yy(1),y0)
      data yy/3*0.d0/
      a0=0.d0
      j2=n-2
      h=x(2)-x(1)
      h2=x(3)-x(1)
      y0=h*h2*(h2-h)
      h=h*h
      h2=h2*h2
      b0=(y(1)*(h-h2)+y(2)*h2-y(3)*h)/y0
      b1=b0
      do 5 i=1,j2
      h=x(i+1)-x(i)
      y0=y(i+1)-y(i)
      h2=h*h
      ha=h-a0
      h2a=h-2.d0*a0
      h3a=2.d0*h-3.*a0
      h2b=h2*b0
      q(1,i)=h2/ha
      q(2,i)=-ha/(h2a*h2)
      q(3,i)=-h*h2a/h3a
      f(1,i)=(y0-h*b0)/(h*ha)
      f(2,i)=(h2b-y0*(2.d0*h-a0))/(h*h2*h2a)
      f(3,i)=-(h2b-3.d0*y0*ha)/(h*h3a)
      a0=q(3,i)
    5 b0=f(3,i)
      i=j2+1
      h=x(i+1)-x(i)
      y0=y(i+1)-y(i)
      h2=h*h
      ha=h-a0
      h2a=h*ha
      h2b=h2*b0-y0*(2.d0*h-a0)
      q(1,i)=h2/ha
      f(1,i)=(y0-h*b0)/h2a
      ha=x(j2)-x(i+1)
      y0=-h*ha*(ha+h)
      ha=ha*ha
      y0=(y(i+1)*(h2-ha)+y(i)*ha-y(j2)*h2)/y0
      q(3,i)=(y0*h2a+h2b)/(h*h2*(h-2.d0*a0))
      q(2,i)=f(1,i)-q(1,i)*q(3,i)
      do 10 j=1,j2
      k=i-1
      q(1,i)=f(3,k)-q(3,k)*q(2,i)
      q(3,k)=f(2,k)-q(2,k)*q(1,i)
      q(2,k)=f(1,k)-q(1,k)*q(3,k)
   10 i=k
      q(1,i)=b1
      do 15 j=1,3
   15 q(j,n)=yy(j)
      return
      end

      subroutine rprop(jf,jl,f)
c*** propagates soln ,f, for radial modes from jf to jl ***
      implicit real*8(a-h,o-z)
      parameter (mk=30000)
      save
      real*8 lcon,ncon,lspl,nspl,nn
      REAL*8 :: l
      common r(mk),fmu(mk),flam(mk),qshear(mk),qkappa(mk),
     + xa2(mk),xlam(mk),rho(mk),qro(3,mk),g(mk),qg(3,mk),
     + fcon(mk),fspl(3,mk),lcon(mk),lspl(3,mk),ncon(mk),
     + nspl(3,mk),ccon(mk),cspl(3,mk),acon(mk),aspl(3,mk)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      common/eifx/a(14,mk),dum(mk)
      common/shanks/b(46),c(10),dx,step(8),stepf,maxo,in
      dimension h(2,10),s(2),f(2)
      maxo1=maxo-1
      y=r(jf)
      vy=sqrt((flam(jf)+2.d0*fmu(jf))/rho(jf))
      i=jf
   10 a(1,i)=f(1)
      a(2,i)=f(2)
      if(i.eq.jl) return
      iq=i
      i=i+1
      x=y
      y=r(i)
      if(y.eq.x) goto 10
      qff=1.d0+xlam(iq)*fct
      qll=1.d0+qshear(iq)*fct
      qaa=1.d0+xa2(iq)*fct
      vx=vy
      vy=sqrt((flam(i)+2.d0*fmu(i))/rho(i))
      q=max(w/vx+1.d0/x,w/vy+1.d0/y)
      del=step(maxo)/q
      dxs=0.d0
   15   y=x+del
        if(y.gt.r(i)) y=r(i)
        dx=y-x
        if(dx.ne.dxs) call baylis(q,maxo1)
        dxs=dx
        s(1)=f(1)
        s(2)=f(2)
        do ni=1,in
          z=x+c(ni)
          t=z-r(iq)
          ro=rho(iq)+t*(qro(1,iq)+t*(qro(2,iq)+t*qro(3,iq)))
          gr=g(iq)+t*(qg(1,iq)+t*(qg(2,iq)+t*qg(3,iq)))
          ff=(fcon(iq)+t*(fspl(1,iq)+t*(fspl(2,iq)+t*fspl(3,iq))))*qff
          if(ifanis.eq.0) then
           nn=(lcon(iq)+t*(lspl(1,iq)+t*(lspl(2,iq)+t*lspl(3,iq))))*qll
           cc=ff+nn+nn
           aa=cc
          else
           nn=(ncon(iq)+t*(nspl(1,iq)+t*(nspl(2,iq)+t*nspl(3,iq))))*qll
           cc=(ccon(iq)+t*(cspl(1,iq)+t*(cspl(2,iq)+t*cspl(3,iq))))*qaa
           aa=(acon(iq)+t*(aspl(1,iq)+t*(aspl(2,iq)+t*aspl(3,iq))))*qaa
          end if
          z=1.d0/z
          a21=-ro*wsq+4.d0*z*(z*(aa-nn-ff*ff/cc)-ro*gr)
          h(1,ni)=(f(2)-2.d0*ff*z*f(1))/cc
          h(2,ni)=a21*f(1)+2.d0*z*f(2)*(ff/cc-1.d0)
          call rkdot(f,s,h,2,ni)
        enddo
        if(knsw.eq.1) then
          if(s(2)*f(2).le.0.d0) then
            if(f(2).eq.0.d0) then
              tes=-s(2)*s(1)
            else
              tes=f(2)*s(1)-s(2)*f(1)
            endif
            if(tes.lt.0.d0) kount=kount+1
            if(tes.gt.0.d0) kount=kount-1
          end if
        end if
        x=y
        if(y.ne.r(i)) go to 15
      go to 10
      end

      subroutine tprop(jf,jl,f)
c*** propagates f from jf to jl - toroidal modes ***
      implicit real*8(a-h,o-z)
      parameter (mk=30000)
      save
      real*8 lcon,ncon,lspl,nspl,nn,ll
      REAL*8 :: l
      common r(mk),fmu(mk),flam(mk),qshear(mk),qkappa(mk),
     + xa2(mk),xlam(mk),rho(mk),qro(3,mk),g(mk),qg(3,mk),
     + fcon(mk),fspl(3,mk),lcon(mk),lspl(3,mk),ncon(mk),
     + nspl(3,mk),ccon(mk),cspl(3,mk),acon(mk),aspl(3,mk)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      common/eifx/a(14,mk),dum(mk)
      common/shanks/b(46),c(10),dx,step(8),stepf,maxo,in
      dimension h(2,10),s(2),f(2)
      fl3m2=fl3-2.d0
      maxo1=maxo-1
      y=r(jf)
      qy=1.d0/y+sqrt(abs(rho(jf)*wsq/fmu(jf)-fl3/(y*y)))
      i=jf
   10 a(1,i)=f(1)
      a(2,i)=f(2)
      if(i.eq.jl) return
      iq=i
      i=i+1
      x=y
      y=r(i)
      if(y.eq.x) goto 10
      qll=1.d0+qshear(iq)*fct
      qx=qy
      qy=1.d0/y+sqrt(abs(rho(i)*wsq/fmu(i)-fl3/(y*y)))
      q=max(qx,qy)
      del=step(maxo)/q
      dxs=0.d0
   15   y=x+del
        if(y.gt.r(i)) y=r(i)
        dx=y-x
        if(dx.ne.dxs) call baylis(q,maxo1)
        dxs=dx
        s(1)=f(1)
        s(2)=f(2)
        do ni=1,in
          z=x+c(ni)
          t=z-r(iq)
          ro=rho(iq)+t*(qro(1,iq)+t*(qro(2,iq)+t*qro(3,iq)))
          ll=(lcon(iq)+t*(lspl(1,iq)+t*(lspl(2,iq)+t*lspl(3,iq))))*qll
          nn=ll
          if(ifanis.ne.0) nn=(ncon(iq)+
     +        t*(nspl(1,iq)+t*(nspl(2,iq)+t*nspl(3,iq))))*qll
          z=1.d0/z
          h(1,ni)=z*f(1)+f(2)/ll
          h(2,ni)=(nn*fl3m2*z*z-ro*wsq)*f(1)-3.d0*z*f(2)
          call rkdot(f,s,h,2,ni)
        enddo
        if(knsw.eq.1) then
          if(s(2)*f(2).le.0.d0) then
            if(f(2).eq.0.d0) then
              tes=-s(2)*s(1)
            else
              tes=f(2)*s(1)-s(2)*f(1)
            endif
            if(tes.lt.0.d0) kount=kount+1
            if(tes.gt.0.d0) kount=kount-1
          end if
        end if
        x=y
        if(y.ne.r(i)) goto 15
      go to 10
      end

      subroutine rkdot(f,s,h,nvec,ni)
c*** performs dot product with rks coefficients ***
      implicit real*8(a-h,o-z)
      save
      common/shanks/b(46),c(10),dx,step(8),stepf,maxo,in
      dimension s(*),f(*),h(nvec,*)
      goto (1,2,3,4,5,6,7,8,9,10),ni
    1 do 21 j=1,nvec
   21 f(j)=s(j)+b(1)*h(j,1)
      return
    2 do 22 j=1,nvec
   22 f(j)=s(j)+b(2)*(h(j,1)+b(3)*h(j,2))
      return
    3 do 23 j=1,nvec
   23 f(j)=s(j)+b(4)*(h(j,1)+b(5)*h(j,2)+b(6)*h(j,3))
      return
    4 do 24 j=1,nvec
   24 f(j)=s(j)+b(7)*(h(j,1)+b(8)*h(j,2)+b(9)*h(j,3)+b(10)*h(j,4))
      return
    5 do 25 j=1,nvec
   25 f(j)=s(j)+b(11)*(h(j,1)+b(12)*h(j,2)+b(13)*h(j,3)+b(14)*h(j,4)+
     +b(15)*h(j,5))
      return
    6 do 26 j=1,nvec
   26 f(j)=s(j)+b(16)*(h(j,1)+b(17)*h(j,2)+b(18)*h(j,3)+b(19)*h(j,4)+
     +b(20)*h(j,5)+b(21)*h(j,6))
      return
    7 do 27 j=1,nvec
   27 f(j)=s(j)+b(22)*(h(j,1)+b(23)*h(j,3)+b(24)*h(j,4)+b(25)*h(j,5)+
     +b(26)*h(j,6)+b(27)*h(j,7))
      return
    8 do 28 j=1,nvec
   28 f(j)=s(j)+b(28)*(h(j,1)+b(29)*h(j,3)+b(30)*h(j,4)+b(31)*h(j,5)+
     +b(32)*h(j,6)+b(33)*h(j,7)+b(34)*h(j,8))
      return
    9 do 29 j=1,nvec
   29 f(j)=s(j)+b(35)*(h(j,1)+b(36)*h(j,3)+b(37)*h(j,4)+b(38)*h(j,5)+
     +b(39)*h(j,6)+b(40)*h(j,7)+b(41)*h(j,8)+b(42)*h(j,9))
      return
   10 do 30 j=1,nvec
   30 f(j)=s(j)+b(43)*(h(j,1)+h(j,10)+b(44)*(h(j,4)+h(j,6))+
     +b(45)*h(j,5)+b(46)*(h(j,7)+h(j,9)))
      return
      end

      subroutine intgds(rr,iq,vals)
c*** interpolates integrands for normalisation,cg,q etc..for use with gauslv.
      implicit real*8(a-h,o-z)
      parameter (mk=30000)
      save
      real*8 lcon,ncon,lspl,nspl,nn,ll
      REAL*8 :: l
      common r(mk),fmu(mk),flam(mk),qshear(mk),qkappa(mk),
     + xa2(mk),xlam(mk),rho(mk),qro(3,mk),g(mk),qg(3,mk),
     + fcon(mk),fspl(3,mk),lcon(mk),lspl(3,mk),ncon(mk),
     + nspl(3,mk),ccon(mk),cspl(3,mk),acon(mk),aspl(3,mk)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      common/eifx/ar(14,mk),dum(mk)
      dimension q(3),qp(3),vals(*)
      data d1,d2,d3,d4,d5,d6,d7/.111111111111111d0,
     + 0.066666666666667d0,0.666666666666667d0,1.333333333333333d0,
     + 2.666666666666667d0,3.333333333333333d0,5.333333333333333d0/
      t=rr-r(iq)
      hn=1.d0/(r(iq+1)-r(iq))
      hsq=hn*hn
      qff=1.d0+xlam(iq)*fct
      qll=1.d0+qshear(iq)*fct
      iq1=iq+1
      ifun=3
      if(jcom.ne.3) ifun=1
      do 10 i=1,ifun
      i2=2*i
      i1=i2-1
      a=((ar(i2,iq)+ar(i2,iq1))+2.d0*hn*(ar(i1,iq)-ar(i1,iq1)))*hsq
      b=-(2.d0*ar(i2,iq)+ar(i2,iq1))*hn-3.d0*(ar(i1,iq)-ar(i1,iq1))*hsq
      q(i)=(ar(i1,iq)+t*(ar(i2,iq)+t*(b+t*a)))/rr
   10 qp(i)=ar(i2,iq)+t*(2.d0*b+t*3.d0*a)
      rro=(rho(iq)+t*(qro(1,iq)+t*(qro(2,iq)+t*qro(3,iq))))*rr
      gr=g(iq)+t*(qg(1,iq)+t*(qg(2,iq)+t*qg(3,iq)))
      ff=(fcon(iq)+t*(fspl(1,iq)+t*(fspl(2,iq)+t*fspl(3,iq))))*qff
      ll=(lcon(iq)+t*(lspl(1,iq)+t*(lspl(2,iq)+t*lspl(3,iq))))*qll
      if(ifanis.ne.0) goto 15
      nn=ll
      cc=ff+ll+ll
      aa=cc
      goto 20
   15 qaa=1.d0+xa2(iq)*fct
      nn=(ncon(iq)+t*(nspl(1,iq)+t*(nspl(2,iq)+t*nspl(3,iq))))*qll
      cc=(ccon(iq)+t*(cspl(1,iq)+t*(cspl(2,iq)+t*cspl(3,iq))))*qaa
      aa=(acon(iq)+t*(aspl(1,iq)+t*(aspl(2,iq)+t*aspl(3,iq))))*qaa
   20 qrka=d1*(4.d0*(aa+ff-nn)+cc)
     1     *(qkappa(iq)+t*hn*(qkappa(iq1)-qkappa(iq)))
      qrmu=d2*(aa+cc-2.d0*ff+5.d0*nn+6.d0*ll)
     1     *(qshear(iq)+t*hn*(qshear(iq1)-qshear(iq)))
      if(jcom.ne.3) goto 25
      q1sq=q(1)*q(1)
      q2sq=q(2)*q(2)
      vals(1)=rr*rro*(q1sq+q2sq)
c      fac=(fl+.5d0)/sfl3
      fac=sfl3/sfl3
      vals(2)=(sfl3*(ll*q1sq+aa*q2sq)+q(2)*((rro*gr+2.d0*(nn-aa-ll)+ff)
     +   *q(1)+rro*q(3)-ff*qp(1))+ll*qp(2)*q(1))*fac
     +   +.25d0*q(3)*(qp(3)+fl*q(3))
      t2=qrka+d7*qrmu
      t3=qrka+d4*qrmu
      t4=qrka+d6*qrmu
      t5=qrka-d5*qrmu
      t6=qrka-d3*qrmu
      vals(3)=.5d0*((fl3*qrmu+t2)*q1sq+(2.d0*qrmu+fl3*t3)*q2sq)
     1 -q(1)*sfl3*t4*q(2)+q(1)*(t5*qp(1)+sfl3*qrmu*qp(2))+q(2)*(-2.d0*
     2 qrmu*qp(2)-sfl3*t6*qp(1))+.5d0*(t3*qp(1)*qp(1)+qrmu*qp(2)*qp(2))
      vals(4)=.5d0*((fl3*ll+4.d0*(rro*(rro-gr)+aa-nn-ff)+cc)*q1sq+
     +(4.d0*ll-nn-nn+fl3*aa)*q2sq +fl*fl*.25d0*q(3)*q(3)+cc*qp(1)*qp(1)+
     +ll*qp(2)*qp(2)+.25d0*qp(3)*qp(3))+q(3)*(rro*sfl3*q(2)+fl*.25d0*qp
     +(3))+q(1)*(sfl3*(rro*gr+2.d0*(nn-aa-ll)+ff)*q(2)+rro*(qp(3)-q(3))+
     +(ff+ff-cc)*qp(1)+sfl3*ll*qp(2))-q(2)*(sfl3*ff*qp(1)+(ll+ll)*qp(2))
      return
   25 q(1)=q(1)*rr
      vals(1)=rr*rro*q(1)*q(1)
      if(jcom.eq.1) goto 30
      vals(2)=nn*q(1)*q(1)
      t1=(rr*qp(1)-q(1))**2
      t2=(fl3-2.d0)*q(1)*q(1)
      vals(3)=(t1+t2)*qrmu
      vals(4)=t1*ll+t2*nn
      return
   30 t1=(rr*qp(1)+2.d0*q(1))**2
      t2=d4*(rr*qp(1)-q(1))**2
      vals(2)=t1*qrka+t2*qrmu
      vals(3)=rr*qp(1)*(cc*rr*qp(1)+4.d0*ff*q(1))+4.d0*q(1)*q(1)
     +    *(aa-nn-rro*gr)
      return
      end

      subroutine fpsm(ls,nvefm,ass)
c*** spheroidal mode start solution in a fluid region using sph. bessel fns.
      implicit real*8(a-h,o-z)
      parameter (mk=30000)
      save
      real*8 lcon,ncon,lspl,nspl
      REAL*8 :: l
      common r(mk),fmu(mk),flam(mk),qshear(mk),qkappa(mk),
     + xa2(mk),xlam(mk),rho(mk),qro(3,mk),g(mk),qg(3,mk),
     + fcon(mk),fspl(3,mk),lcon(mk),lspl(3,mk),ncon(mk),
     + nspl(3,mk),ccon(mk),cspl(3,mk),acon(mk),aspl(3,mk)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      dimension ass(*)
      x=r(ls)
      fla=flam(ls)*(1.d0+xlam(ls)*fct)
      vpsq=fla/rho(ls)
      xi=g(ls)/x
      qsq=(wsq+float(kg)*4.d0*rho(ls)+xi-fl3*xi*xi/wsq)/vpsq
      zsq=qsq*x*x
      call bfs(l,zsq,eps,fp)
      if(kg.ne.0) then
        u=(fl-fp)/qsq
        c1=fl*g(ls)-wsq*x
        c2=fl2*c1*0.25d0/x-rho(ls)*fl
        ass(1)=-x*fl*vpsq-c1*u
        ass(2)=-x*fl*fla
        ass(3)=-fl*fl2*vpsq*0.25d0-u*c2
        ass(4)=x*fla*c1
        ass(5)=-x*fla*c2
      else
        ass(1)=-(fl3*xi/wsq+fp)/qsq
        ass(2)=x*fla
      end if
      sum=ass(1)*ass(1)
      do i=2,nvefm
        sum=sum+ass(i)*ass(i)
      enddo
      sum=1.d0/dsqrt(sum)
      if(ass(nvefm).lt.0.d0) sum=-sum
      do i=1,nvefm
        ass(i)=ass(i)*sum
      enddo
      return
      end

      subroutine spsm(ls,nvesm,ass)
c*** spheroidal mode start solution in a solid region using sph. bessel fns.
      implicit real*8(a-h,o-z)
      parameter (mk=30000)
      save
      real*8 lcon,ncon,lspl,nspl
      REAL*8 :: l
      common r(mk),fmu(mk),flam(mk),qshear(mk),qkappa(mk),
     + xa2(mk),xlam(mk),rho(mk),qro(3,mk),g(mk),qg(3,mk),
     + fcon(mk),fspl(3,mk),lcon(mk),lspl(3,mk),ncon(mk),
     + nspl(3,mk),ccon(mk),cspl(3,mk),acon(mk),aspl(3,mk)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      dimension a(6,2),e(15),ass(*)
      x=r(ls)
      ro=rho(ls)
      fu=fmu(ls)*(1.d0+qshear(ls)*fct)
      flu=flam(ls)*(1.d0+xlam(ls)*fct)+2.d0*fu
      vssq=fu/ro
      vpsq=flu/ro
      zeta=4.d0*ro
      xi=g(ls)/x
      alfsq=(wsq+float(kg)*zeta+xi)/vpsq
      betasq=wsq/vssq
      delsq=dsqrt((betasq-alfsq)**2+4.d0*fl3*xi*xi/(vpsq*vssq))
      fksq=.5d0*(alfsq+betasq+delsq)
      qsq=fksq-delsq
      do k=1,2
        if(k.eq.1) then
          zsq=qsq*x*x
          b=xi/(vssq*(betasq-qsq))
        else
          zsq=fksq*x*x
          b=-flu/(fu*fl3*b)
        end if
        call bfs(l,zsq,eps,fp)
        a(1,k)=fl3*b+fp
        a(2,k)=1.d0+b+b*fp
        a(3,k)=-zsq
        a(4,k)=b*a(3,k)
        a(5,k)=1.d0
        a(6,k)=fl1-fl3*b
      enddo
      jj=3+2*kg
      kk=jj+1
      ll=0
      do i=1,jj
        i1=i+1
        do j=i1,kk
          ll=ll+1
          e(ll)=a(i,1)*a(j,2)-a(j,1)*a(i,2)
        enddo
      enddo
      if(kg.eq.0) then
        ass(1)=x*x*e(1)
        ass(2)=fu*x*sfl3*(2.d0*e(1)-e(5))
        ass(3)=fu*x*(e(3)-2.d0*e(1))
        ass(4)=x*(flu*e(4)+4.d0*fu*e(1))
        ass(5)=fu*(flu*(e(6)+2.d0*e(4))+4.d0*fu*(fl3*(e(5)-e(1))
     +     -e(3)+2.d0*e(1)))
      else
        c0=wsq-xi*fl
        c1=ro*fl+0.25d0*fl2*c0
        c2=2.d0*fu/x
        c3=c2*(fl-1.d0)
        ass(6)=x*x*(c0*e(1)-zeta*(fl*e(8)-e(4)))
        ass(14)=flu*(fl*e(6)-e(2))
        ass(13)=fu*sfl3*(fl*e(7)-e(3))
        ass(1)=x*(c1*e(1)-ro*(fl*e(9)-e(5)))
        ass(7)=x*flu*(c0*e(2)-zeta*fl*e(11))/sfl3+c2*sfl3*ass(6)
        ass(8)=x*fu*(c0*e(3)-zeta*fl*e(13))-c2*ass(6)
        ass(12)=(flu*fl*e(10)+2.d0*(ass(14)+sfl3*ass(13)))*fu/x
        ass(2)=flu*(c1*e(2)-ro*fl*e(12))/sfl3+c2*sfl3*ass(1)
        ass(3)=fu*(c1*e(3)-ro*fl*e(14))-c2*ass(1)
        ass(9)=(x*c0*ass(14)+sfl3*ass(7)-c3*fl*ass(6))/fl
        ass(11)=(sfl3*ass(12)+c3*(sfl3*ass(14)-fl*ass(13)))/fl
        ass(4)=(c1*ass(14)+sfl3*ass(2)-c3*fl*ass(1))/fl
        ass(10)=(x*c0*ass(11)-c3*(sfl3*ass(9)+fl*ass(7)))/sfl3
        ass(5)=(c1*ass(11)-c3*(sfl3*ass(4)+fl*ass(2)))/sfl3
      end if
      sum=ass(1)*ass(1)
      do i=2,nvesm
        sum=sum+ass(i)*ass(i)
      enddo
      sum=1.d0/dsqrt(sum)
      if(ass(5).lt.0.d0) sum=-sum
      do i=1,nvesm
        ass(i)=ass(i)*sum
      enddo
      return
      end

      subroutine rps(i,a)
c*** radial mode start soln using sph bessel fns.
      implicit real*8(a-h,o-z)
      parameter (mk=30000)
      save
      real*8 lcon,ncon,lspl,nspl
      REAL*8 :: l
      common r(mk),fmu(mk),flam(mk),qshear(mk),qkappa(mk),
     + xa2(mk),xlam(mk),rho(mk),qro(3,mk),g(mk),qg(3,mk),
     + fcon(mk),fspl(3,mk),lcon(mk),lspl(3,mk),ncon(mk),
     + nspl(3,mk),ccon(mk),cspl(3,mk),acon(mk),aspl(3,mk)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      dimension a(2)
      fla=flam(i)*(1.d0+xlam(i)*fct)
      sig=fla+2.d0*fmu(i)*(1.d0+qshear(i)*fct)
      zsq=r(i)*r(i)*rho(i)*(wsq+4.d0*g(i)/r(i))/sig
      call bfs(1.d0,zsq,eps,fp)
      a(1)=r(i)
      a(2)=sig*fp+2.d0*fla
      return
      end

      subroutine tps(i,a)
c*** toroidal mode start soln using sph bessel fns.
      implicit real*8(a-h,o-z)
      parameter (mk=30000)
      save
      real*8 lcon,ncon,lspl,nspl
      REAL*8 :: l
      common r(mk),fmu(mk),flam(mk),qshear(mk),qkappa(mk),
     + xa2(mk),xlam(mk),rho(mk),qro(3,mk),g(mk),qg(3,mk),
     + fcon(mk),fspl(3,mk),lcon(mk),lspl(3,mk),ncon(mk),
     + nspl(3,mk),ccon(mk),cspl(3,mk),acon(mk),aspl(3,mk)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      dimension a(2)
      fu=fmu(i)*(1.d0+qshear(i)*fct)
      zsq=r(i)*r(i)*wsq*rho(i)/fu
      call bfs(l,zsq,eps,fp)
      a(1)=r(i)
      a(2)=fu*(fp-1.d0)
      return
      end

      subroutine bfs(l,xsq,eps,fp)
c  this routine calculates spherical bessel function of the ist kind.
c  fp is equivalent to (r*dj/dr)/j
c  where r is radius and j is the sbf of order l and argument x=k*r
c  the technique employs the continued fraction approach
c  described in w. lentz's article in applied optics, vol.15, #3, 1976
c
c  H.Y. Yang, Mar 2015: change l and lp1 to real(KIND=r8b)
c
      implicit real*8(a-h,o-z)
      save
      REAL*8 :: l, lp1
      real*8 numer,nu
c*** positive argument uses continued fraction
      if(xsq.gt.0.d0) then
        x=sqrt(xsq)
        lp1=l+1.d0
        rx=2.0d0/x
        nu=lp1-0.5d0
        rj=nu*rx
        rx=-rx
        denom=(nu+1.d0)*rx
        numer=denom+1.0d0/rj
        rj=rj*numer/denom
        nm1=1
    5     nm1=nm1+1
          rx=-rx
          a3=(nu+nm1)*rx
          denom=a3+1.d0/denom
          numer=a3+1.d0/numer
          ratio=numer/denom
          rj=rj*ratio
          if(abs(abs(ratio)-1.d0).gt.eps) goto 5
        fp=rj*x-lp1
      else
c  series solution
        f=1.d0
        fp=l
        a=1.d0
        b=l+l+1.d0
        c=2.d0
        d=l+2.d0
   15     a=-a*xsq/(c*(b+c))
          f=f+a
          fp=fp+a*d
          if(abs(a*d).lt.eps) goto 20
          c=c+2.d0
          d=d+2.d0
          goto 15
   20   fp=fp/f
      end if
      return
      end

      subroutine remedy(ls)
c    obtains the eigenfunction of an awkward spheroidal mode by
c    integrating to the icb or the mcb.
      implicit real*8(a-h,o-z)
      parameter (mk=30000)
      save
      REAL*8 :: l
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      common/eifx/ar(14,mk),inorm(mk),jjj(mk)
      common/arem/a(6,3,mk)
      common/rindx/nic,noc,nsl,nicp1,nocp1,nslp1,n
      dimension af(4,2),as(6,3),afr(4)
c MB commented two next lines
cx    print 900,ls
cx900 format('in remedy with start level : ',i4)
      iexp=0
      do k=1,2
        do j=1,4
          af(j,k)=0.d0
        enddo
      enddo
      af(1,1)=1.d0
      if(kg.eq.1) af(2,2)=1.d0
      if(nsl.ne.n) then
        do i=nslp1,n
          do k=1,3
            do j=1,6
              a(j,k,i)=0.d0
            enddo
          enddo
        enddo
        call fprop(n,nslp1,af,iexp)
      end if
      call fsbdry(af,as,kg)
      do k=1,3
        do j=1,6
          a(j,k,nsl)=as(j,k)
        enddo
      enddo
      if(n.ne.nsl) call ortho(n,nsl,as,kg)
c*** now we are at the bottom of the ocean
      call sprop(n,nsl,nocp1,as,iexp)
      call sfbdry(n,nocp1,as,af,kg)
      imtch=noc
      do i=1,4
        afr(i)=ar(i,noc)
      enddo
c*** now we are at the cmb -- test to see if we might have to go to the icb
      if(ls.le.nic) then
        icomp=0
        call match(n,noc,kg,af,afr,icomp)
        if(icomp.eq.0) return
        call fprop(noc,nicp1,af,iexp)
        imtch=nic
        do i=1,4
          afr(i)=ar(i,nicp1)
        enddo
      end if
      icomp=-1
      call match(n,imtch,kg,af,afr,icomp)
      return
      end

      subroutine fprop(jf,jl,f,iexp)
c    fprop propagates the fundamental matrix f from jf to jl (a fluid region)
      implicit real*8(a-h,o-z)
      parameter (mk=30000)
      save
      real*8 lcon,ncon,lspl,nspl
      REAL*8 :: l
      common r(mk),fmu(mk),flam(mk),qshear(mk),qkappa(mk),
     + xa2(mk),xlam(mk),rho(mk),qro(3,mk),g(mk),qg(3,mk),
     + fcon(mk),fspl(3,mk),lcon(mk),lspl(3,mk),ncon(mk),
     + nspl(3,mk),ccon(mk),cspl(3,mk),acon(mk),aspl(3,mk)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      common/eifx/ar(14,mk),inorm(mk),jjj(mk)
      common/arem/a(6,3,mk)
      common/rindx/nic,noc,nsl,nicp1,nocp1,nslp1,n
      common/shanks/b(46),c(10),dx,step(8),sdum,idum,in
      dimension f(4,2),s(4,2),h(4,2,10)
      data econst/1048576.d0/
      kk=kg+1
      jj=2*kk
      jud=1
      if(jl.lt.jf) jud=-1
      y=r(jf)
      i=jf
      go to 80
   10 x=y
      y=r(i)
      if(y.eq.x) goto 80
      iq=min(i,i-jud)
      qff=1.d0+xlam(iq)*fct
      xi=g(i)/y
      alfsq=(wsq+4.d0*rho(i)+xi-fl3*xi*xi/wsq)*rho(i)/flam(i)
      q=max(sfl3/x,sqrt(abs(alfsq-fl3/(x*x)))+1.d0/r(iq))
      del=jud*step(8)/q
      dxs=0.d0
   15 y=x+del
      if(float(jud)*(y-r(i)).gt.0.d0) y=r(i)
      dx=y-x
      if(dx.ne.dxs) call baylis(q,7)
      dxs=dx
      do k=1,kk
        do j=1,jj
          s(j,k)=f(j,k)
        enddo
      enddo
      d=fl3/wsq
      do 40 ni=1,in
      z=x+c(ni)
      t=z-r(iq)
      zr=1.d0/z
      ro=rho(iq)+t*(qro(1,iq)+t*(qro(2,iq)+t*qro(3,iq)))
      flu=(fcon(iq)+t*(fspl(1,iq)+t*(fspl(2,iq)+t*fspl(3,iq))))*qff
      gr=(g(iq)+t*(qg(1,iq)+t*(qg(2,iq)+t*qg(3,iq))))*zr
      t21=-4.d0*ro
      t12=d*zr*zr
      t11=(gr*d-1.d0)*zr
      s11=-ro*(wsq+4.d0*gr-gr*gr*d)
      c11=-t12/ro+1.d0/flu
      if(kg.eq.0) then
        s11=s11-t21*ro
      else
        t22=-fl*zr
        s22=ro*t12
        s12=ro*(t11+t22)
      end if
      do 70 k=1,kk
      if(kg.eq.0) then
        h(1,k,ni)=t11*f(1,k)+c11*f(2,k)
        h(2,k,ni)=s11*f(1,k)-t11*f(2,k)
      else
        h(1,k,ni)=t11*f(1,k)+t12*f(2,k)+c11*f(3,k)
        h(2,k,ni)=t21*f(1,k)+t22*f(2,k)+4.d0*f(4,k)
        h(3,k,ni)=s11*f(1,k)+s12*f(2,k)-t11*f(3,k)-t21*f(4,k)
        h(4,k,ni)=s12*f(1,k)+s22*f(2,k)-t12*f(3,k)-t22*f(4,k)
      end if
      do 70 j=1,jj
      go to (701,702,703,704,705,706,707,708,709,710),ni
  701 f(j,k)=s(j,k)+b(1)*h(j,k,1)
      go to 70
  702 f(j,k)=s(j,k)+b(2)*(h(j,k,1)+b(3)*h(j,k,2))
      go to 70
  703 f(j,k)=s(j,k)+b(4)*(h(j,k,1)+b(5)*h(j,k,2)+b(6)*h(j,k,3))
      go to 70
  704 f(j,k)=s(j,k)+b(7)*(h(j,k,1)+b(8)*h(j,k,2)+b(9)*h(j,k,3)+
     +b(10)*h(j,k,4))
      go to 70
  705 f(j,k)=s(j,k)+b(11)*(h(j,k,1)+b(12)*h(j,k,2)+b(13)*h(j,k,3)+
     +b(14)*h(j,k,4)+b(15)*h(j,k,5))
      go to 70
  706 f(j,k)=s(j,k)+b(16)*(h(j,k,1)+b(17)*h(j,k,2)+b(18)*h(j,k,3)+
     +b(19)*h(j,k,4)+b(20)*h(j,k,5)+b(21)*h(j,k,6))
      go to 70
  707 f(j,k)=s(j,k)+b(22)*(h(j,k,1)+b(23)*h(j,k,3)+b(24)*h(j,k,4)+
     +b(25)*h(j,k,5)+b(26)*h(j,k,6)+b(27)*h(j,k,7))
      go to 70
  708 f(j,k)=s(j,k)+b(28)*(h(j,k,1)+b(29)*h(j,k,3)+b(30)*h(j,k,4)+
     +b(31)*h(j,k,5)+b(32)*h(j,k,6)+b(33)*h(j,k,7)+b(34)*h(j,k,8))
      go to 70
  709 f(j,k)=s(j,k)+b(35)*(h(j,k,1)+b(36)*h(j,k,3)+b(37)*h(j,k,4)+
     +b(38)*h(j,k,5)+b(39)*h(j,k,6)+b(40)*h(j,k,7)+b(41)*h(j,k,8)+
     +b(42)*h(j,k,9))
      go to 70
  710 f(j,k)=s(j,k)+b(43)*(h(j,k,1)+h(j,k,10)+b(45)*h(j,k,5)+
     +b(44)*(h(j,k,4)+h(j,k,6))+b(46)*(h(j,k,7)+h(j,k,9)))
   70 continue
   40 continue
      x=y
      if(y.ne.r(i)) go to 15
   80 size=0.d0
      do k=1,kk
        do j=1,jj
          size=max(size,abs(f(j,k)))
        enddo
      enddo
   82 if(size.lt.1024.d0) goto 84
      do k=1,kk
        do j=1,jj
          f(j,k)=f(j,k)/econst
        enddo
      enddo
      size=size/econst
      iexp=iexp+20
      goto 82
   84 inorm(i)=iexp
      do k=1,kk
        do j=1,jj
          a(j,k,i)=f(j,k)
        enddo
      enddo
      if(i.eq.jl) return
      i=i+jud
      go to 10
      end

      subroutine sprop(li,jf,jl,f,iexp)
c    sprop propagates the fundamental matrix f from jf to jl (a solid region)
c    if jorth=1 the columns of f are orthogonalized at each level
c    except in regions of oscillatory p and s.
      implicit real*8(a-h,o-z)
      parameter (mk=30000)
      save
      real*8 lcon,ncon,lspl,nspl,nn,ll
      common r(mk),fmu(mk),flam(mk),qshear(mk),qkappa(mk),
     + xa2(mk),xlam(mk),rho(mk),qro(3,mk),g(mk),qg(3,mk),
     + fcon(mk),fspl(3,mk),lcon(mk),lspl(3,mk),ncon(mk),
     + nspl(3,mk),ccon(mk),cspl(3,mk),acon(mk),aspl(3,mk)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      common/eifx/ar(14,mk),inorm(mk),jjj(mk)
      common/arem/a(6,3,mk)
      common/rindx/nic,noc,nsl,nicp1,nocp1,nslp1,n
      common/shanks/b(46),c(10),dx,step(8),sdum,idum,in
      dimension f(6,3),s(6,3),h(6,3,10)
      data econst/1048576.d0/
      kk=kg+2
      jj=2*kk
      jud=1
      if(jl.lt.jf) jud=-1
      y=r(jf)
      i=jf
      go to 80
   10 x=y
      y=r(i)
      if(x.eq.y) goto 80
      iq=min0(i,i-jud)
      qff=1.d0+xlam(iq)*fct
      qll=1.d0+qshear(iq)*fct
      qaa=1.d0+xa2(iq)*fct
      xi=g(i)/y
      vpsq=(flam(i)+2.d0*fmu(i))/rho(i)
      vssq=fmu(i)/rho(i)
      alfsq=(wsq+4.d0*rho(i)+xi)/vpsq
      betasq=wsq/vssq
      delsq=dsqrt((betasq-alfsq)**2+4.d0*fl3*xi*xi/(vssq*vpsq))
      fksq=.5d0*(alfsq+betasq+delsq)
      al=fl3/(x*x)
      jorth=1
      aq=fksq-delsq-al
      if(aq.gt.0.d0) jorth=0
      qs=sqrt(abs(fksq-al))+1.d0/r(iq)
      qf=sqrt(abs(aq))+1.d0/r(iq)
      q=max(sfl3/x,qs,qf)
      del=jud*step(8)/q
      dxs=0.d0
   15 y=x+del
      if(float(jud)*(y-r(i)).gt.0.d0) y=r(i)
      dx=y-x
      if(dx.ne.dxs) call baylis(q,7)
      dxs=dx
      do k=1,kk
        do j=1,jj
          s(j,k)=f(j,k)
        enddo
      enddo
      do 50 ni=1,in
      z=x+c(ni)
      t=z-r(iq)
      ro=rho(iq)+t*(qro(1,iq)+t*(qro(2,iq)+t*qro(3,iq)))
      gr=g(iq)+t*(qg(1,iq)+t*(qg(2,iq)+t*qg(3,iq)))
      ff=(fcon(iq)+t*(fspl(1,iq)+t*(fspl(2,iq)+t*fspl(3,iq))))*qff
      ll=(lcon(iq)+t*(lspl(1,iq)+t*(lspl(2,iq)+t*lspl(3,iq))))*qll
      if(ifanis.eq.0) then
        nn=ll
        cc=ff+ll+ll
        aa=cc
      else
        nn=(ncon(iq)+t*(nspl(1,iq)+t*(nspl(2,iq)+t*nspl(3,iq))))*qll
        cc=(ccon(iq)+t*(cspl(1,iq)+t*(cspl(2,iq)+t*cspl(3,iq))))*qaa
        aa=(acon(iq)+t*(aspl(1,iq)+t*(aspl(2,iq)+t*aspl(3,iq))))*qaa
      end if
      zr=1.d0/z
      sfl3z=sfl3*zr
      rogr=ro*gr
      c11=1.d0/cc
      c22=1.d0/ll
      dmg=aa-nn-ff*ff*c11
      zdmg=zr*dmg
      t11=-2.d0*ff*zr*c11+zr
      t12=sfl3z*ff*c11
      t21=-sfl3z
      t22=zr+zr
      s22=-ro*wsq
      s11=s22+4.d0*zr*(zdmg-rogr)
      s22=s22+zr*zr*(fl3*(dmg+nn)-nn-nn)
      s12=sfl3z*(rogr-zdmg-zdmg)
      if(kg.eq.0)  then
        s11=s11+4.d0*ro*ro
      else
        t31=-4.d0*ro
        t33=-fl*zr
        s13=-fl1*zr*ro
        s23=ro*sfl3z
      end if
      do 70 k=1,kk
      if(kg.eq.0) then
        h(1,k,ni)=t11*f(1,k)+t12*f(2,k)+c11*f(3,k)
        h(2,k,ni)=t21*f(1,k)+t22*f(2,k)+c22*f(4,k)
        h(3,k,ni)=s11*f(1,k)+s12*f(2,k)-t11*f(3,k)-t21*f(4,k)
        h(4,k,ni)=s12*f(1,k)+s22*f(2,k)-t12*f(3,k)-t22*f(4,k)
      else
        h(1,k,ni)=t11*f(1,k)+t12*f(2,k)+c11*f(4,k)
        h(2,k,ni)=t21*f(1,k)+t22*f(2,k)+c22*f(5,k)
        h(3,k,ni)=t31*f(1,k)+t33*f(3,k)+4.d0*f(6,k)
        h(4,k,ni)=s11*f(1,k)+s12*f(2,k)+s13*f(3,k)-t11*f(4,k)-t21*f(5,k)
     +    -t31*f(6,k)
        h(5,k,ni)=s12*f(1,k)+s22*f(2,k)+s23*f(3,k)-t12*f(4,k)-t22*f(5,k)
        h(6,k,ni)=s13*f(1,k)+s23*f(2,k)-t33*f(6,k)
      end if
      do 70 j=1,jj
      go to (701,702,703,704,705,706,707,708,709,710),ni
  701 f(j,k)=s(j,k)+b(1)*h(j,k,1)
      go to 70
  702 f(j,k)=s(j,k)+b(2)*(h(j,k,1)+b(3)*h(j,k,2))
      go to 70
  703 f(j,k)=s(j,k)+b(4)*(h(j,k,1)+b(5)*h(j,k,2)+b(6)*h(j,k,3))
      go to 70
  704 f(j,k)=s(j,k)+b(7)*(h(j,k,1)+b(8)*h(j,k,2)+b(9)*h(j,k,3)+
     +b(10)*h(j,k,4))
      go to 70
  705 f(j,k)=s(j,k)+b(11)*(h(j,k,1)+b(12)*h(j,k,2)+b(13)*h(j,k,3)+
     +b(14)*h(j,k,4)+b(15)*h(j,k,5))
      go to 70
  706 f(j,k)=s(j,k)+b(16)*(h(j,k,1)+b(17)*h(j,k,2)+b(18)*h(j,k,3)+
     +b(19)*h(j,k,4)+b(20)*h(j,k,5)+b(21)*h(j,k,6))
      go to 70
  707 f(j,k)=s(j,k)+b(22)*(h(j,k,1)+b(23)*h(j,k,3)+b(24)*h(j,k,4)+
     +b(25)*h(j,k,5)+b(26)*h(j,k,6)+b(27)*h(j,k,7))
      go to 70
  708 f(j,k)=s(j,k)+b(28)*(h(j,k,1)+b(29)*h(j,k,3)+b(30)*h(j,k,4)+
     +b(31)*h(j,k,5)+b(32)*h(j,k,6)+b(33)*h(j,k,7)+b(34)*h(j,k,8))
      go to 70
  709 f(j,k)=s(j,k)+b(35)*(h(j,k,1)+b(36)*h(j,k,3)+b(37)*h(j,k,4)+
     +b(38)*h(j,k,5)+b(39)*h(j,k,6)+b(40)*h(j,k,7)+b(41)*h(j,k,8)+
     +b(42)*h(j,k,9))
      go to 70
  710 f(j,k)=s(j,k)+b(43)*(h(j,k,1)+h(j,k,10)+b(45)*h(j,k,5)+
     +b(44)*(h(j,k,4)+h(j,k,6))+b(46)*(h(j,k,7)+h(j,k,9)))
   70 continue
   50 continue
      x=y
      if(y.ne.r(i)) go to 15
   80 size=0.d0
      do k=1,kk
        do j=1,jj
          size=max(size,abs(f(j,k)))
        enddo
      enddo
   82 if(size.lt.1024.d0) goto 84
      do k=1,kk
        do j=1,jj
          f(j,k)=f(j,k)/econst
        enddo
      enddo
      size=size/econst
      iexp=iexp+20
      goto 82
   84 inorm(i)=iexp
      do k=1,kk
        do j=1,jj
          a(j,k,i)=f(j,k)
        enddo
      enddo
      if(jorth.eq.1) call ortho(li,i,f,kg)
      if(i.eq.jl) return
      i=i+jud
      go to 10
      end

      subroutine sfbdry(jf,jl,as,af,kg)
c*** The tangential traction scalar is forced to vanish at the solid
c*** side of a s/f boundary(level jl). a(j,3,i) is elliminated for
c*** i=jf...jl and af is loaded from a at level jl.
      implicit real*8(a-h,o-z)
      parameter (mk=30000)
      save
      common/arem/a(6,3,mk)
      dimension as(6,*),af(4,*)
      n1=min(jf,jl)
      n2=max(jf,jl)
      if(kg.eq.0) then
        if(abs(as(4,2)).gt.abs(as(4,1))) then
          i1=1
          i2=2
        else
          i1=2
          i2=1
        end if
        rat=-as(4,i1)/as(4,i2)
        do i=n1,n2
          do j=1,4
            a(j,1,i)=a(j,i1,i)+rat*a(j,i2,i)
          enddo
        enddo
        af(1,1)=a(1,1,jl)
        af(2,1)=a(3,1,jl)
      else
        ab53=abs(as(5,3))
        do k=1,2
          if(ab53.gt.dabs(as(5,k))) then
            i1=k
            i2=3
          else
            i1=3
            i2=k
          end if
          rat=-as(5,i1)/as(5,i2)
          do i=n1,n2
            do j=1,6
              a(j,k,i)=a(j,i1,i)+rat*a(j,i2,i)
            enddo
          enddo
          af(1,k)=a(1,k,jl)
          af(2,k)=a(3,k,jl)
          af(3,k)=a(4,k,jl)
          af(4,k)=a(6,k,jl)
        enddo
      end if
      return
      end

      subroutine fsbdry(af,as,kg)
c    fsbdry creates solid fundamental matrix as from fluid fundamental matrix
c    af. It is presumed that fsbdry is used to cross a f/s boundary.
      implicit real*8(a-h,o-z)
      save
      dimension af(4,*),as(6,*)
      do i=1,3
        do j=1,6
          as(j,i)=0.d0
        enddo
      enddo
      if(kg.eq.0) then
        as(1,1)=af(1,1)
        as(3,1)=af(2,1)
        as(2,2)=1.d0
      else
        do k=1,2
          as(1,k)=af(1,k)
          as(3,k)=af(2,k)
          as(4,k)=af(3,k)
          as(6,k)=af(4,k)
        enddo
        as(2,3)=1.d0
      end if
      return
      end

      subroutine match(n,j,kg,af,afr,icomp)
c*** matches boundary conditions during downward integration --
c*** first try is at cmb, second try is at icb.
      implicit real*8(a-h,o-z)
      parameter (mk=30000)
      save
      common/eifx/ar(14,mk),inorm(mk),jjj(mk)
      common/arem/a(6,3,mk)
      dimension af(4,*),afr(*),afi(4)
      k=j+2
      rms=0.d0
      fnor=0.d0
      if(kg.eq.0) then
        c=(af(1,1)*afr(1)+af(2,1)*afr(2))/(af(1,1)**2+af(2,1)**2)
        do i=1,2
          afi(i)=af(i,1)*c
          rms=rms+(afi(i)-afr(i))**2
          fnor=fnor+afr(i)*afr(i)
        enddo
        rms=sqrt(rms/fnor)
        if(icomp.ge.0) then
          if(rms.ge.1.d-3) then
            icomp=1
            return
          end if
        end if
        idiff=inorm(j)-inorm(j+1)
        inorm(j+1)=inorm(j)
        do i=k,n
          inorm(i)=inorm(i)+idiff
          do jj=1,4
            ar(jj,i)=c*a(jj,1,i)
          enddo
        enddo
      else
        a2=(af(3,1)*afr(1)-af(1,1)*afr(3))/(af(1,2)*af(3,1)-af(1,1)
     +   *af(3,2))
        a1=(af(3,2)*afr(1)-af(1,2)*afr(3))/(af(1,1)*af(3,2)-af(3,1)
     +   *af(1,2))
        do i=1,4
          afi(i)=a1*af(i,1)+a2*af(i,2)
          rms=rms+(afi(i)-afr(i))**2
          fnor=fnor+afr(i)*afr(i)
        enddo
        rms=sqrt(rms/fnor)
        if(icomp.ge.0) then
          if(rms.ge.1.d-3) then
            icomp=1
            return
          end if
        end if
        idiff=inorm(j)-inorm(j+1)
        inorm(j+1)=inorm(j)
        do i=k,n
          inorm(i)=inorm(i)+idiff
          do jj=1,6
            ar(jj,i)=a1*a(jj,1,i)+a2*a(jj,2,i)
          enddo
        enddo
      end if
      return
      end

      subroutine ortho(li,lc,b,kg)
c   Finds the orthogonal matrix v such that the columns of b*v are orthogonal,
c   the array a is replaced by a*v for levels li - lc. Array b is replaced
c   by b*v and is then ready for entry to sprop at level lc.This is intended
c   to diminish the onset of degeneracy caused by rapid exponential growth
c   in the mantle for modes with deeply turning S and shallowly turning P.
      implicit real*8(a-h,o-z)
      parameter (mk=30000)
      save
      common/arem/a(6,3,mk)
      dimension b(6,*),as(6,3)
      i1=min(lc,li)
      i2=max(lc,li)
      nc=kg+2
      nr=2*nc
      call svd(b,nr,nc)
      do i=i1,i2
        do j=1,nc
          do k=1,nr
            as(k,j)=0.d0
            do l=1,nc
              as(k,j)=as(k,j)+a(k,l,i)*b(l,j)
            enddo
          enddo
        enddo
        do j=1,nc
          do k=1,nr
            a(k,j,i)=as(k,j)
          enddo
        enddo
      enddo
      do j=1,nc
        do k=1,nr
          b(k,j)=a(k,j,lc)
        enddo
      enddo
      return
      end

      subroutine svd(a,mrow,ncol)
c    section i chapter 10 wilkenson and reinsch (1971 ,springer).
c    the matrix a is overwritten with v(ncol,ncol), the right side orthogonal
c    matrix in the svd decomposition. for use only in eos subs as ,to reduce
c    branching points, i have used the fact that ncol is lt mrow.
      implicit real*8(a-h,o-z)
      save
      !dimension a(6,*),e(3),q(3)
      dimension a(6,*),e(4),q(3)
      eps=1.5d-14
      tol=1.d-293
      g=0.d0
      x=0.d0
      do 60 i=1,ncol
      l=i+1
      e(i)=g
      s=0.d0
      do j=i,mrow
        s=s+a(j,i)*a(j,i)
      enddo
      if(s.le.tol) then
        q(i)=0.d0
        if(l.gt.ncol) goto 60
      else
        q(i)=dsign(sqrt(s),-a(i,i))
        h=a(i,i)*q(i)-s
        a(i,i)=a(i,i)-q(i)
        if(l.gt.ncol) go to 60
        do j=l,ncol
          s=0.d0
          do k=i,mrow
            s=s+a(k,i)*a(k,j)
          enddo
          f=s/h
          do k=i,mrow
            a(k,j)=a(k,j)+f*a(k,i)
          enddo
        enddo
      end if
      s=0.d0
      do j=l,ncol
        s=s+a(i,j)*a(i,j)
      enddo
      if(s.lt.tol) then
        g=0.d0
      else
        g=dsign(dsqrt(s),-a(i,l))
        h=a(i,l)*g-s
        a(i,l)=a(i,l)-g
        do j=l,ncol
          e(j)=a(i,j)/h
        enddo
        do j=l,mrow
          s=0.d0
          do k=l,ncol
            s=s+a(j,k)*a(i,k)
          enddo
          do k=l,ncol
            a(j,k)=a(j,k)+s*e(k)
          enddo
        enddo
      end if
   60 x=dmax1(dabs(q(i))+dabs(e(i)),x)
      goto 100
   75 if(g.ne.0.d0) then
        h=a(i,l)*g
        do j=l,ncol
          a(j,i)=a(i,j)/h
        enddo
        do j=l,ncol
          s=0.d0
          do k=l,ncol
            s=s+a(i,k)*a(k,j)
          enddo
          do k=l,ncol
            a(k,j)=a(k,j)+s*a(k,i)
          enddo
        enddo
      end if
      do j=l,ncol
        a(i,j)=0.d0
        a(j,i)=0.d0
      enddo
  100 a(i,i)=1.d0
      !print*,'ncol,i:',ncol,i
      g=e(i)
      l=i
      i=i-1
      if(i.ge.1)go to 75
      ep=eps*x
      k=ncol
  105 l=k
  110 if(dabs(e(l)).le.ep)go to 125
      if(dabs(q(l-1)).le.ep) go to 115
      l=l-1
      if(l.ge.1)go to 110
  115 c=0.d0
      s=1.d0
      do i=l,k
        f=s*e(i)
        e(i)=c*e(i)
        if(abs(f).le.ep)go to 125
        g=q(i)
        h=sqrt(f*f+g*g)
        c=g/h
        s=-f/h
        q(i)=h
      enddo
  125 z=q(k)
      if(l.ne.k) then
        x=q(l)
        y=q(k-1)
        g=e(k-1)
        h=e(k)
        f=((y-z)*(y+z)+(g-h)*(g+h))/(2.d0*h*y)
        g=dsqrt(f*f+1.d0)
        f=((x-z)*(x+z)+h*(y/(f+dsign(g,f))-h))/x
        c=1.d0
        s=1.d0
        lp1=l+1
        do i=lp1,k
          g=e(i)
          y=q(i)
          h=s*g
          g=c*g
          z=dsqrt(f*f+h*h)
          im1=i-1
          e(im1)=z
          c=f/z
          s=h/z
          f=s*g+c*x
          g=c*g-s*x
          h=s*y
          y=c*y
          do j=1,ncol
            x=a(j,im1)
            z=a(j,i)
            a(j,im1)=c*x+s*z
            a(j,i)=c*z-s*x
          enddo
          z=sqrt(f*f+h*h)
          q(im1)=z
          c=f/z
          s=h/z
          f=s*y+c*g
          x=c*y-s*g
        enddo
        e(l)=0.d0
        e(k)=f
        q(k)=x
        go to 105
      end if
      if(z.lt.0.d0) then
        q(k)=-z
        do j=1,ncol
          a(j,k)=-a(j,k)
        enddo
      end if
      k=k-1
      if(k.ge.1)go to 105
      return
      end

      SUBROUTINE intgrl(r,nir,ner,f,s1,s2,s3,sum0)
c         Compute the integral of f[i] from i=nir to i=nir:ner for
c         radii values as in model stored in f
c         sum0: result (array)
! modified by Chuanming, for forrtl warning, array temporary
! add         parameter (mk=30000)
        IMPLICIT NONE
        INTEGER mk
        PARAMETER (mk=30000)
*        REAL*8 sum0(*)
*        REAL*8 r(*), f(*), s1(*), s2(*), s3(*)
        REAL*8 sum0(mk)
        REAL*8 r(mk), f(mk), s1(mk), s2(mk), s3(mk)
        INTEGER nir, ner

        REAL*8 third,fourth,fifth,sixth
        REAL*8 rji,dsum
        INTEGER nir1, i, j

        third=1.d0/3.d0
        fourth=1.d0/4.d0
        fifth=1.d0/5.d0
        sixth=1.d0/6.d0

        sum0(1)=0.d0
        nir1 = nir + 1
        DO i=nir1,ner
          j = i-1
          rji = r(i) - r(j)
          dsum=rji*(f(j)+rji*(.50*s1(j)+rji*(third*s2(j)+rji*
     +              fourth*s3(j))))+2.0*r(j)*rji*rji*(.50*f(j)+rji*(third*s1(j)+rji*
     +              (fourth*s2(j)+rji*fifth*s3(j))))+rji*rji*rji*(third*f(j)+rji*
     +            (fourth*s1(j)+rji*(fifth*s2(j)+rji*sixth*s3(j))))
          sum0(i)=sum0(j)+dsum
        END DO
      END SUBROUTINE intgrl

      SUBROUTINE kernel_velrhop(ichoice,fkernid,wcom,qmod,gcom,wdiff)
c----------------------------------------------------------------------------
c  <Input>
c     ichoice: =1, dw/d(alpha,beta,rhop,d)
c              =2, dc/d(alpha,beta,rhop,d)
c              =3, dw/d(alpha,betap,sqrt(mu),d)
c              =4, dc/d(alpha,betap,sqrt(mu),d)
c     save K\alpha, K\beta, K\rhop and K\d in a(11:14,mk)
c     save K\alpha, K\betap, K\sqrt(mu), and K\d in a(11:14,mk)
c
c H.Y. Yang, Mar 2015
c allow to calculate kernels for TI model
c
c H.Y. Yang, Jan 2012
c Ouput kernels are normalized according D.&T., rather than the rules of MINEOS, but
c are still non-dimensional.
c debug for Love wave

c
c H.Y. Yang, Nov 2011
c Parameters such as w, wsq, cg, rho, flam, fmu are dimensionaless.
c
c----------------------------------------------------------------------------
        IMPLICIT REAL*8(a-h,o-z)
        PARAMETER(mk=30000,nhead=6)
        REAL*8 lcon,ncon,lspl,nspl

        INTEGER :: ichoice, fkernid
        REAL*8 :: wcom, qmod, gcom

        INTEGER :: krec(mk), nrec
        REAL*8 :: l
        COMMON/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback

        COMMON/eifx/a(14,mk),ndum(mk)
        COMMON/rindx/nic,noc,nsl,nicp1,nocp1,nslp1,n
        COMMON r(mk),fmu(mk),flam(mk),qshear(mk),
     +qkappa(mk),xa2(mk),xlam(mk),rho(mk),
     +qro(3,mk),g(mk),qg(3,mk),fcon(mk),
     +fspl(3,mk),lcon(mk),lspl(3,mk),ncon(mk),
     +nspl(3,mk),ccon(mk),cspl(3,mk),acon(mk),
     +aspl(3,mk)
        COMMON/eigsav/krec,nrec,nrecic,nrecoc,nkrec

        DOUBLE PRECISION, ALLOCATABLE :: kern(:)
        INTEGER :: i, i1, istart, j,Ier
        DOUBLE PRECISION :: rr2, rr,rri, twoukv, qint(mk),qy, dr, dq,
     +     twowi, fac, fac2, fkappa, xx, zz, alpha, beta
        LOGICAL :: isset1, isset2
        DOUBLE PRECISION :: tmpk(3), wnorm
        INTEGER :: nvec

        isset1=(ichoice == 1) .OR. (ichoice==2)
        isset2=(ichoice == 3) .OR. (ichoice==4)
        IF (.NOT.(isset1 .OR. isset2)) STOP 'please check ichoice in sub kernel_velrhop, Terminate!!'
        IF (ifanis /= 0) STOP 'model has to be isotropic, Terminated!'
        fac=4.d0/3.d0
        fac2=2.d0/3.d0
        twowi=0.5d0/w

        IF (jcom == 3) THEN
c--------------------------------------------------------
c   % Spheroidal mode %
c-------------------------------------------------------
          qy=0.d0
          DO i = n, 1, -1
            rr2=r(i)**2
            rr=r(i)
            rri=1.d0/rr
            IF (rr==0.d0) rri=rr
            twoukv=2.d0*a(1,i)-sfl3*a(3,i)
            xx=rr*a(4,i)-a(3,i)+sfl3*a(1,i)
            qint(i)=rho(i)*a(1,i)*twoukv*rri
c           ONLY use trapzoid integration
            IF (i < n) THEN
              i1=i+1
              dr=r(i1)-r(i)
              dq=(qint(i)+qint(i1))*dr*0.5d0
              qy=qy+dq
            END IF
            fkappa=flam(i)+fac2*fmu(i)

c--------------------------------------------------------
c           a(11:14,:)=[eqn(9.13)-eqn(9.16)]/2w in D.&T.,1998
c--------------------------------------------------------
            a(11,i)=twowi*(rr*a(2,i)+twoukv)**2
            a(12,i)=twowi*(((2.d0*rr*a(2,i)-twoukv)**2)/3.d0+
     +                     xx**2+
     +                     (fl3-2.0d0)*a(3,i)**2)
            a(13,i)=twowi*((-wsq*rr2*(a(1,i)**2+a(3,i)**2))+
     +                      8.d0*rho(i)*rr2*a(1,i)**2+
     +                      2.d0*(rr2*a(1,i)*a(6,i)+sfl3*rr*a(3,i)*a(5,i))+
     +                     -2.d0*g(i)*rr*a(1,i)*twoukv+
     +                     -8.d0*rr2*qy)
            a(14,i)=-(fkappa*a(11,i)+fmu(i)*a(12,i)+rho(i)*a(13,i))+
     +               (fkappa*a(2,i)*(rr*a(2,i)+twoukv)+
     +                fac2*fmu(i)*a(2,i)*(2.d0*rr*a(2,i)-twoukv)+
     +                fmu(i)*a(4,i)*xx)*rr/w
            alpha=SQRT((flam(i)+2.d0*fmu(i))/rho(i))
            beta =SQRT(fmu(i)/rho(i))
            IF (isset1) THEN
c--------------------------------------------------------
c           a(11:13,:)=[eqn(9.21)-eqn(9.23)] in D.&T.,1998
c--------------------------------------------------------
              a(13,i)=a(13,i)+
     +                a(12,i)*beta**2+
     +                a(11,i)*(alpha**2-fac*beta**2)
              a(12,i)=2.d0*rho(i)*beta*(a(12,i)-fac*a(11,i))
              a(11,i)=2.d0*rho(i)*alpha*a(11,i)

            ELSE IF (isset2) THEN
              tmpk(1:3)=a(11:13,i)
              a(11,i)=2.d0*rho(i)*alpha*tmpk(1)
              IF (beta==0.d0) THEN
                a(12,i)= 0.d0
                a(13,i)= 0.d0   !\ I am not sure
              ELSE
                a(12,i)= -2.d0*rho(i)/beta*(tmpk(3)+alpha**2*tmpk(1))
                a(13,i)=  2.0*SQRT(rho(i))*beta*(
     +                      ((alpha/beta)**2-fac)*tmpk(1)+
     +                                            tmpk(2)+
     +                          1.d0/beta**2     *tmpk(3)  )
              END IF
            END IF
          END DO

        ELSE IF (jcom == 2) THEN
c-----------------------------------------------------------
c   % Toroidal mode %
c-----------------------------------------------------------
          DO i = n, nocp1 , -1
            rr=r(i)
            rr2=rr**2
            rri=1.d0/rr
            zz=rr*a(2,i)-a(1,i)
            a(11,i)=0.d0
            a(12,i)=twowi*(zz**2+(fl3-2.d0)*a(1,i)**2)
            a(13,i)=twowi*(-wsq*rr2*a(1,i)**2)
            a(14,i)=-(fkappa*a(11,i)+fmu(i)*a(12,i)+rho(i)*a(13,i))+
     +                fmu(i)*a(2,i)*zz*rr/w

            alpha=SQRT((flam(i)+2.d0*fmu(i))/rho(i))
            beta =SQRT(fmu(i)/rho(i))
            IF (isset1) THEN
              a(13,i)=a(13,i)+
     +                a(12,i)*beta**2+
     +                a(11,i)*(alpha**2-fac*beta**2)
              a(12,i)=2.d0*rho(i)*beta*(a(12,i)-fac*a(11,i))
              a(11,i)=2.d0*rho(i)*alpha*a(11,i)

            ELSE IF (isset2) THEN
              tmpk(1:3)=a(11:13,i)
              a(11,i)=2.d0*rho(i)*alpha*tmpk(1)
              IF (beta==0.d0) THEN
                a(12,i)= 0.d0
                a(13,i)= 0.d0   !\ I am not sure
              ELSE
                a(12,i)= -2.d0*rho(i)/beta*(tmpk(3)+alpha**2*tmpk(1))
                a(13,i)=  2.0*SQRT(rho(i))*beta*(
     +                      ((alpha/beta)**2-fac)*tmpk(1)+
     +                                            tmpk(2)+
     +                          1.d0/beta**2     *tmpk(3)  )
              END IF
            END IF
          END DO

        ELSE
          STOP 'jcom out of range in kernel_velrhop, Terminated!'
        END IF

        IF ((ichoice==2) .OR. (ichoice==4)) THEN
c-------------------------------------------------------------
c         a(11:14,:)=eqn(11.94-95) in D.& T., 1998
c------------------------------------------------------------
          a(11:14,1:n)=w/(cg*fl3)*a(11:14,1:n)
        END IF

c------------------------------------------------------------
c   % write out sensitivity kernel %
c-----------------------------------------------------------
        nvec=nhead+4*nkrec
        IF ( .NOT. ALLOCATED(kern)) THEN
          ALLOCATE(kern(nvec),STAT=Ier)
          IF (Ier/=0) STOP 'allocate error of kern in sub kernel_velrhop'
        END IF
        kern=0.d0
        kern(1)=DBLE(nord)
        kern(2)=DBLE(l)
        kern(3)=DBLE(wcom)
        kern(4)=DBLE(qmod)
        kern(5)=DBLE(gcom)
        kern(6)=DBLE(wdiff)

        IF(jcom.EQ.2) THEN
          istart=nocp1
        ELSE IF (jcom.EQ.3) THEN
          istart=1
        ENDIF
        j=0
        DO i = istart, n
          IF (krec(i).EQ.1) THEN
            j=j+1
            kern(j+nhead)=a(11,i)
            kern(j+nkrec+nhead)=a(12,i)
            kern(j+2*nkrec+nhead)=a(13,i)
            kern(j+3*nkrec+nhead)=a(14,i)
          END IF
        END DO
c    --------------------------------------------------------------------------------
c      normalized kernels by the rules of D.&.T.
c      since kernel is composed of U^2-like stuff, so
c        kern should be multipled by wnorm^2
c    -----------------------------------------------------------------------------------
        wnorm=kern(3)/wn
        kern(nhead+1:nvec)=wnorm**2*kern(nhead+1:nvec)

        WRITE(UNIT=fkernid) (kern(i),i=1,nvec)
        CALL FLUSH(fkernid)

      END SUBROUTINE kernel_velrhop


      SUBROUTINE kernel_velrhopTI(ichoice,fkernid,wcom,qmod,gcom,wdiff)
c----------------------------------------------------------------------------
c calculate absolute sensitivity kernel dw/dm
c         dw = int [dw/dm]dm dr
c         dw can be replaced by dc (see $ichoice)
c
c  <Input>
c     ichoice: =1, dw/d(alphab,betah,rhop,d,eta,alphah,betav)
c              =2, dc/d(alphab,betah,rhop,d,eta,alphah,betav)
c		   =3, dc/d(A,C,F,L,N,rhop,d)
c     save in kern(nhead+nrec*7)
c
c H.Y. Yang, Mar 2015
c -allow to calculate kernels for TI model
c -Using cubic spline integration to calculate element of Krho, instead
c   of trapzoidal method (have already checked).
c -using drspln and intgrl
c
c H.Y. Yang, Jan 2012
c Ouput kernels are normalized according D.&T., rather than the rules of MINEOS, but
c are still non-dimensional.
c debug for Love wave
c
c
c H.Y. Yang, Nov 2011
c Parameters such as w, wsq, cg (group velocity), rho, fcon, lcon, ncon, acon, ccon are dimensionaless.
c and physical dispersion is not accounted in elastic parameters (e.g. ccon)
c a(1:6,:)=(U,dU,V,dV,P,dP) or a(1:2,:)=(W,dW) whose normalization follows the rules in D.&T. (1998),
c  rather than MINEOS
c
c Note fl3 : l*(l+1)=k^2
c      sfl3: k
c      cg : group velocity
c      w  : eigenfreq
c      wsq: w^2
c
c----------------------------------------------------------------------------
        IMPLICIT REAL*8(a-h,o-z)
        PARAMETER (mk=30000, nhead=6)
        REAL*8 lcon,ncon,lspl,nspl

        INTEGER :: ichoice, fkernid
c                               wcom: w in rad/s, gcom: gv in km/s
        REAL*8 :: wcom, qmod, gcom

        INTEGER :: krec(mk), nrec
        REAL*8 :: l
        COMMON/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback

        COMMON/eifx/a(14,mk),ndum(mk)
        COMMON/rindx/nic,noc,nsl,nicp1,nocp1,nslp1,n
        COMMON r(mk),fmu(mk),flam(mk),qshear(mk),
     +qkappa(mk),xa2(mk),xlam(mk),rho(mk),
     +qro(3,mk),g(mk),qg(3,mk),fcon(mk),
     +fspl(3,mk),lcon(mk),lspl(3,mk),ncon(mk),
     +nspl(3,mk),ccon(mk),cspl(3,mk),acon(mk),
     +aspl(3,mk)
        COMMON/eigsav/krec,nrec,nrecic,nrecoc,nkrec

        DOUBLE PRECISION, ALLOCATABLE :: kern(:)
        INTEGER :: i, i1, istart, j,Ier
        DOUBLE PRECISION :: rr2, rr,rri, twoukv,
     +     twowi, fac, fac2, xx, zz, alphav, alphah, betaa, betah, eta
        LOGICAL :: isset1, isset2
        DOUBLE PRECISION :: tmpk(7), wnorm
c       DOUBLE PRECISION :: qint(mk),qy, dr, dq
        DOUBLE PRECISION :: qint(mk),qqint(3,mk), wrk(3,mk), sumy(mk),qy
        isset1=(ichoice == 1) .OR. (ichoice==2) .OR. (ichoice==3)
        isset2=(ichoice == 0) .OR. (ichoice==4)
        IF (.NOT.(isset1)) STOP 'please check ichoice in sub kernel_velrhop, Terminate!!'
        IF (ifanis /= 1 .AND. ifanis /=0) STOP 'ifanis is out of range, Terminated!'
        fac=4.d0/3.d0
        fac2=2.d0/3.d0
        twowi=0.5d0/w

        IF (jcom == 3) THEN
c--------------------------------------------------------
c   % Spheroidal mode %
c-------------------------------------------------------
c!---------
c          qy=0.d0
c!---------
          DO i = n, 1, -1
            rri=1.d0/r(i)
            IF (r(i)==0.d0) rri=0.d0
            twoukv=2.d0*a(1,i)-sfl3*a(3,i)
            qint(i)=rho(i)*a(1,i)*twoukv*rri
c!----------------
c                 ONLY use trapzoid integration to get the integral of
c                 rho(r')*U(r')*twoukv(r')/r' over r to a
c            IF (i < n) THEN
c              i1=i+1
c              dr=r(i1)-r(i)
c              dq=(qint(i)+qint(i1))*dr*0.5d0
c              qy=qy+dq
c            END IF
c!-------------------
          END DO
c
c                          Cubic spline integration
c
          CALL drspln(1,n,r,qint,qqint,wrk)
!                         Note the integral starts from 0 to r
          CALL intgrl(r,1,n,qint,qqint(1,1:n),qqint(2,1:n),qqint(3,1:n),
     +                            sumy(1:n))
!                            but what we need is from r to Ra -->
!                            arrange sumy: qy
!

          DO i = n, 1, -1
            qy=sumy(n)-sumy(i)
            rr2=r(i)**2
            rr=r(i)
            rri=1.d0/rr
            IF (rr==0.d0) rri=rr
            twoukv=2.d0*a(1,i)-sfl3*a(3,i)
            xx=rr*a(4,i)-a(3,i)+sfl3*a(1,i)

c                             elastic parameters
            alphav=SQRT(ccon(i)/rho(i))
            alphah=SQRT(acon(i)/rho(i))
            betah =SQRT(ncon(i)/rho(i))
            betav =SQRT(lcon(i)/rho(i))
            eta   =fcon(i)/(acon(i)-2.d0*lcon(i))
c--------------------------------------------------------
c           a(11:14,:)=2w*[Kc, Kn, Krho, Kd]
c                     =eqs.[(9.31), (9.33), (9.15), (9.35)] in D.&T.(1998)
c          &
c           a(8:10,:) =2w*[Kf, Ka, Kl]
c                     =eqs.[(9.34), (9.32), (9.31)]
c--------------------------------------------------------
            a(11,i)=(rr*a(2,i))**2                   !2w*Kc (9.31) - related to Vpv
            a(9, i)= twoukv**2                       !2w*Ka (9.31) -            Vph
            a(12,i)=-twoukv**2+(fl3-2.0d0)*a(3,i)**2 !2w*Kn (9.33) -            Vsh
            a(10,i)= xx**2                           !2w*Kl (9.32) -            Vsv
            a( 8,i)= 2.d0*rr*a(2,i)*twoukv           !2w*Kf (9.34)

                                                     !2w*Krho (9.15)
            a(13,i)=(-wsq*rr2*(a(1,i)**2+a(3,i)**2))+
     +                      8.d0*rho(i)*rr2*a(1,i)**2+
     +                      2.d0*(rr2*a(1,i)*a(6,i)+sfl3*rr*a(3,i)*a(5,i))+
     +                     -2.d0*g(i)*rr*a(1,i)*twoukv+
     +                     -8.d0*rr2*qy
                                                     !2w*Kd (9.35)
            a(14,i)= a(11,i)*ccon(i)
     +              -a(9,i) *acon(i)
     +              -a(10,i)*lcon(i)
     +              -a(12,i)*ncon(i)
     +              -a(13,i)*rho(i)
     +              +(2.d0*rr*a(4,i)*xx)*lcon(i)

            IF ((ichoice==1) .OR. (ichoice==2)) THEN
c--------------------------------------------------------
c           a(11:14,:)=[Kav, Kbh, Krho', Kd] = (9.40, 9.41, 9.43,[9.35/2w])
c           a(8:10, :)=[Keta, Kah, Kbv] = (9.42, 9.40, 9.41) in D.&T.,1998
c--------------------------------------------------------
              tmpk(1:7)=a(8:14,i) !tmpk=2w*(Kf, Ka, Kl, Kc, Kn, Krho, Kd)

              a(11,i)=rho(i)/w*alphav*tmpk(4)   !Kav   (9.40)
              a(12,i)=rho(i)/w*betah *tmpk(5)   !Kbh   (9.41)
                                                !Krho' (9.43)
              a(13,i)=(alphav**2*tmpk(4)+
     +                 alphah**2*(tmpk(2)+     eta*tmpk(1))+
     +                 betav**2 *(tmpk(3)-2.d0*eta*tmpk(1))+
     +                 betah**2 * tmpk(5)+
     +                            tmpk(6))*twowi
                                                !Kd    (9.35)/2w
              a(14,i)=tmpk(7)*twowi
                                                !Keta  (9.42)
              a(8, i)=rho(i)*(alphah**2-2.d0*betav**2)*tmpk(1)*twowi
                                                !Kah   (9.40)
              a(9, i)=rho(i)/w*alphah*(tmpk(2)+eta*tmpk(1))
                                                !Kbv   (9.41)
              a(10,i)=rho(i)/w*betav*(tmpk(3)-2.d0*eta*tmpk(1))

		ELSE IF (ichoice==3) THEN
!----------------------------------------------------------------------
!             dc/d(A,C,F,L,N,rhop,d)
!----------------------------------------------------------------------
	         a(11,i)=a(11,i)/w*0.5d0                !K_C
	         a(9, i)=a(9, i)/w*0.5d0                !K_A
		   a(12,i)=a(12,i)/w*0.5d0			!K_N
		   a(10,i)=a(10,i)/w*0.5d0			!K_L
		   a( 8,i)=a( 8,i)/w*0.5d0		      !K_F
		   a(13,i)=a(13,i)/w*0.5d0		      !K_rho
		   a(14,i)=a(14,i)/w*0.5d0		      !K_d

            ELSE IF (isset2) THEN
              STOP 'ichoice is not allowd, Terminated!'
            END IF
          END DO

        ELSE IF (jcom == 2) THEN
c-----------------------------------------------------------
c   % Toroidal mode %
c-----------------------------------------------------------
          DO i = n, nocp1 , -1
            rr=r(i)
            rr2=rr**2
            rri=1.d0/rr
            zz=rr*a(2,i)-a(1,i)  ! rdW-W
c--------------------------------------------------------
c           tmpk=2w*(Kn, Kl, Krho, Kd)
c               =eqs.[(9.33), (9.32), (9.15), (9.35)] in D.&T.(1998)
c           Note that tmpk for jcom=2 is different from that for jcom=3
c--------------------------------------------------------
            tmpk(1)=(fl3-2.d0)*a(1,i)**2        !2w*Kn    (9.33)
            tmpk(2)=zz**2                       !2w*Kl    (9.32)
            tmpk(3)=-wsq*rr2*a(1,i)**2          !2w*Krho  (9.15)
                                                !2w*Kd    (9.35)
            tmpk(4)=-tmpk(1)*ncon(i)
     +              -tmpk(2)*lcon(i)
     +              -tmpk(3)*rho(i)
     +              +2.d0*rr*a(2,i)*zz*lcon(i)

c            IF (i==n-2) THEN
c              WRITE(*,*) 'D=1km :', r(i)*rn, a(1,i),a(2,i), tmpk(3)
c            END IF

c                             elastic parameters
            alphav=SQRT(ccon(i)/rho(i))
            alphah=SQRT(acon(i)/rho(i))
            betah =SQRT(ncon(i)/rho(i))
            betav =SQRT(lcon(i)/rho(i))
            eta   =fcon(i)/(acon(i)-2.d0*lcon(i))

            IF (isset1) THEN
c--------------------------------------------------------
c           a(11:14,:)=[Kav=0, Kbh, Krho', Kd] = (9.40, 9.41, 9.43,[9.35/2w])
c           a(8:10, :)=[Keta=0, Kah=0, Kbv] = (9.42, 9.40, 9.41) in D.&T.,1998
c--------------------------------------------------------
              a(8,i)=0.d0                    !Keta (9.42)
              a(9,i)=0.d0                    !Kah  (9.40)
              a(11,i)=0.d0                   !Kav  (9.40)
              a(10,i)=rho(i)/w*betav*tmpk(2) !Kbv  (9.41)
              a(12,i)=rho(i)/w*betah*tmpk(1) !Kbh  (9.41)
                                             !Krho'(9.43)
              a(13,i)=(betav**2*tmpk(2)+
     +                 betah**2*tmpk(1)+
     +                 tmpk(3))*twowi
              a(14,i)=tmpk(4)*twowi          !Kd   (9.35)/2w
            ELSE IF (ichoice==3) THEN
!----------------------------------------------------------------------
!           ichoice=3  output dc/d(A,C,F,L,N)
!----------------------------------------------------------------------
	        a(11,i)=0.d0                   !K_C
	        a(9, i)=0.d0                   !K_A
		  a(12,i)=tmpk(1)/w*0.5d0	   !K_N
		  a(10,i)=tmpk(2)/w*0.5d0	   !K_L
		  a( 8,i)=0.d0   		         !K_F
		  a(13,i)=tmpk(3)/w*0.5d0	   !K_rho
		  a(14,i)=tmpk(4)/w*0.5d0	   !K_d
            ELSE IF (isset2) THEN ! change
              STOP 'ichoice is not allowd, Terminated!'
            END IF
          END DO

        ELSE
          STOP 'jcom out of range in kernel_velrhopTI, Terminated!'
        END IF


        IF ((ichoice==2).OR.(ichoice==3)) THEN
c                                   !until now K=dw/(dm,dd), where m and d are
c                                   !velocity and discontinuity, respectively
c                                   !Now we are looking for K'=dc/(dm,dd),
c                                   !so K'=c^2/(U*w)*K [eqs.11.94-11.95], where
c                                   ! c: phase velocity = w/sfl3 (in prog)
c                                   ! U: group velocity = cg (in prog)
          a(8:14,1:n)=w/(cg*fl3)*a(8:14,1:n)
        END IF

c------------------------------------------------------------
c   % write out sensitivity kernel %
c-----------------------------------------------------------
        nvec=nhead+7*nkrec
c        WRITE(*,*) 'nkrec, nvec=', nkrec, nvec
        IF ( .NOT. ALLOCATED(kern)) THEN
         ALLOCATE(kern(nvec),STAT=Ier)
         IF (Ier/=0) STOP 'allocate error of kern in sub kernel_velrhop'
        END IF
        kern=0.d0
        kern(1)=DBLE(nord)
        kern(2)=DBLE(l)
        kern(3)=DBLE(wcom)
        kern(4)=DBLE(qmod)
        kern(5)=DBLE(gcom)
        kern(6)=DBLE(wdiff)
*        kern(7)=DBLE(ichoice)
        IF(jcom.EQ.2) THEN
          istart=nocp1
        ELSE IF (jcom.EQ.3) THEN
          istart=1
        ENDIF
        j=0
        DO i = istart, n
          IF (krec(i).EQ.1) THEN
            j=j+1                                        ! Spheriodal/Toriodal   ! dc/d(A,C,F,N,L)
            kern(j+nhead)=a(11,i)                        ! dc/dav -                dc/dC
            kern(j+nkrec+nhead)=a(9,i)                   ! dc/dah                  dc/dA
            kern(j+2*nkrec+nhead)=a(10,i)                ! dc/dbv                  dc/dL
            kern(j+3*nkrec+nhead)=a(12,i)                ! dc/dbh                  dc/dN
            kern(j+4*nkrec+nhead)=a( 8,i)                ! dc/deta                 dc/dF
            kern(j+5*nkrec+nhead)=a(13,i)                ! dc/drho'                dc/drho
            kern(j+6*nkrec+nhead)=a(14,i)                ! dc/dd                   dc/dd
          END IF
        END DO

c    --------------------------------------------------------------------------------
c      normalized kernels by the rules of D.&.T.
c      since kernel is composed of U^2-like stuff, so
c        kern should be multipled by wnorm^2
c    -----------------------------------------------------------------------------------
        wnorm=kern(3)/wn
        kern(nhead+1:nvec)=wnorm**2*kern(nhead+1:nvec)

        WRITE(UNIT=fkernid) (kern(i),i=1,nvec)
        CALL FLUSH(fkernid)

      END SUBROUTINE kernel_velrhopTI

