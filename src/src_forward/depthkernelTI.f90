! use tregn96 from cps to calculate the dcdL and dcdA
subroutine depthkernelTI(nx, ny, nz, vel, pvRc, iwave, igr, kmaxRc, tRc, depz, minthk, Lsen_Gsc)
    use omp_lib
    implicit none

    integer nx, ny, nz
    real vel(nx, ny, nz)
    integer iwave,igr
    real minthk
    real depz(nz)
    integer kmaxRc
    real*8 tRc(kmaxRc)
    ! output
    real*8 pvRc(nx*ny, kmaxRc)
    ! parameter list
    real vpz(nz),vsz(nz),rhoz(nz)
    integer mmax,iflsph,mode,rmax
    integer ii,jj,k,i,j,jjj
    integer, parameter:: NL=200
    integer, parameter:: NP=60

    real*8 cgRc(NP)
    real rdep(NL),rvp(NL),rvs(NL),rrho(NL),rthk(NL)
    ! for tregn96
    real t_in(kmaxRc), cp_in(kmaxRc)
    real TA_in(NL), TC_in(NL), TF_in(NL)
    real TL_in(NL), TN_in(NL), TRho_in(NL)
    real qp(NL), qs(NL), etap(NL)
    real etas(NL), frefp(NL), frefs(NL)

    real*4 dcdah(NP,NL),dcdn(NP,NL)
    real*4 dcdbv(NP,NL)
    real*4 dcR_dL, dcR_dA
    real*4 Lsen_Gsc(nx*ny, kmaxRc, nz-1)
    integer nsublay(NL), post

    mmax=nz
    iflsph=1
    mode=1
    pvRc=0.0
    Lsen_Gsc=0.0
    write(6, *) ' depth kernel parallel:'
    !$omp parallel &
    !$omp default(private) &
    !$omp shared(nx,ny,nz,vel,minthk,mmax,depz,kmaxRc) &
    !$omp shared(tRc,pvRc,iflsph,iwave,mode,igr,Lsen_Gsc)
    !$omp do
    do jj=1,ny
        do ii=1,nx
            post=ii+(jj-1)*nx
            vsz(1:nz)=vel(ii,jj,1:nz)
            ! some other emperical relationship maybe better,
            do k=1,nz
                vpz(k)=0.9409 + 2.0947*vsz(k) - 0.8206*vsz(k)**2+ &
                0.2683*vsz(k)**3 - 0.0251*vsz(k)**4
                rhoz(k)=1.6612*vpz(k) - 0.4721*vpz(k)**2 + &
                0.0671*vpz(k)**3 - 0.0043*vpz(k)**4 + &
                0.000106*vpz(k)**5
            enddo
            ! change from refineGrid2LayerMdl into refineLayerMdl
            ! call refineGrid2LayerMdl(minthk, mmax, depz, vpz, vsz, rhoz, rmax, rdep, &
            ! rvp, rvs, rrho, rthk)
            call refineLayerMdl(minthk, mmax, depz, vpz, vsz, rhoz, rmax, rdep, &
            rvp, rvs, rrho, rthk, nsublay)

            call surfdisp96(rthk, rvp, rvs, rrho, rmax, iflsph, iwave, mode, igr, kmaxRc, &
            tRc, cgRc)
            pvRc(ii+(jj-1)*nx,1:kmaxRc)=cgRc(1:kmaxRc)
            !print*,cgRc(1:kmaxRc)
            !------------------------------------------------------------------!
            do i = 1, rmax
                TA_in(i)=rrho(i)*rvp(i)**2
                TC_in(i)=TA_in(i)
                TL_in(i)=rrho(i)*rvs(i)**2
                TN_in(i)=TL_in(i)
                TF_in(i)=1.0*(TA_in(i) - 2 * TL_in(i))
                TRho_in(i)=rrho(i)
            enddo
            qp(1:rmax)=150.0
            qs(1:rmax)=50.0
            etap(1:rmax)=0.00
            etas(1:rmax)=0.00
            frefp(1:rmax)=1.00
            frefs(1:rmax)=1.00

            cp_in(1:kmaxRc)=sngl(cgRc(1:kmaxRc))
            t_in(1:kmaxRc)=sngl(tRc(1:kmaxRc))

            ! ! write(6, *)'tregn96'
            call tregn96(rmax, rthk, TA_in, TC_in, TF_in, TL_in, TN_in, TRho_in, &
            qp, qs, etap, etas, frefp, frefs,  &
            kmaxRc, t_in, cp_in(1:kmaxRc),&
            dcdah, dcdbv, dcdn)
            !
            ! ! write(*,*)"nsublay:", nsublay(1:nz)
            do i=1,kmaxRc  ! period
                k=0
                do j=1,nz-1                ! inversion layer
                    do jjj=1,nsublay(j)    ! refined layer k-th in jth inversion layer
                        k=k+1
                        dcR_dA = 0.5/(rrho(k)*rvp(k))*dcdah(i, k) - TF_in(k)/((TA_in(k)-2.0*TL_in(k))**2)*dcdn(i,k)
                        dcR_dL = 0.5/(rrho(k)*rvs(k))*dcdbv(i, k) + 2.0*TF_in(k)/((TA_in(k)-2.0*TL_in(k))**2)*dcdn(i,k)
                        Lsen_Gsc(post,i,j)=Lsen_Gsc(post,i,j)+dcR_dA*TA_in(k)+dcR_dL*TL_in(k)
                    enddo
                enddo
            enddo

        enddo
    enddo
    !$omp end do
    !$omp end parallel
end subroutine

! do i=1, kmaxRc
! 	do j=1, rmax
! 		dcR_dA(i, j)= 0.5/(rrho(j)*rvp(j))*dcdah(i, j) - TF_in(j)/((TA_in(j)-2.0*TL_in(j))**2)*dcdn(i,j)
! 		dcR_dL(i, j)= 0.5/(rrho(j)*rvs(j))*dcdbv(i, j) + 2.0*TF_in(j)/((TA_in(j)-2.0*TL_in(j))**2)*dcdn(i,j)
! 	enddo
! enddo

! subroutine test96(rmax, rthk, TA_in, TC_in, TF_in, TL_in, TN_in, TRho_in, &
!             qp, qs, etap, etas, frefp, frefs,  &
!             kmaxRc, t_in, cp_in,&
!             dcdah, dcdbv, dcdn)
!
!     implicit none
!     integer rmax
!     integer, parameter:: NL=200
!     real rthk(NL)
!     real t_in(kmaxRc), cp_in(kmaxRc)
!     real TA_in(NL), TC_in(NL), TF_in(NL)
!     real TL_in(NL), TN_in(NL), TRho_in(NL)
!     real qp(NL), qs(NL), etap(NL)
!     real etas(NL), frefp(NL), frefs(NL)
!
!
! end subroutine
