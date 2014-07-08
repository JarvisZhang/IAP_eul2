#include <misc.h>
#include <params.h>

module restart_dynamics
!--------------------------------------------------------------------
! Modified by ZhangHe, 2007.5.25
! Update : ZhangHe, 2008.6.5
!          ZhangHe, 2008.6.21
!          ZhangHe, 2008.7.21, write the IAP dynamical fields from 
!                              (beglatdyn --> endlatdyn) to (beglat --> endlat) 
!--------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid   
                                         !jjr
! =====================================================
   use IAP_grid,     only: NX, NL, NZ,IB,IE,period       
! =====================================================
   use prognostics
   use IAP_prog
   use constituents, only: pnats, pcnst, ppcnst
   use ppgrid, only: pcols, pver
   use binary_io     ! for wrtout_r8 & readin_r8
   use scanslt,      only: scanslt_alloc, lammp, phimp, sigmp, qfcst
   use massfix,      only: alpha, hw1, hw2, hw3
#if ( defined BFB_CAM_SCAM_IOP )
   use iop,          only: alphasav, dqfx3savm1, divq3dsav, divt3dsav, t3sav, u3sav, &
                           v3sav, t2sav, q3sav, pssav, tssav, fixmassav, betasav,    &
                           init_iop_fields
#endif

  use mod_comm, only: mp_send3d, mp_recv3d
   use abortutils, only: endrun
   use constituents, only: cnst_need_pdeldry

   implicit none

CONTAINS

   subroutine write_restart_dynamics (nrg)

#include <comqfl.h>

!
! Input arguments
!
      integer, intent(in) :: nrg     ! Unit number
!
! Local workspace
!
      integer :: ioerr   ! error status
      integer :: k,i, j,m, jd,num
! =================== zhh 2008.7.21 ===================
      real(r8), allocatable :: tmp2d(:,:)
      real(r8), allocatable :: tmp3d(:,:,:)
      real(r8), allocatable :: temp3d(:,:,:)
      real(r8), allocatable :: tmp3dz(:,:,:)
! =================== zhh 2008.7.21 ===================
!--------------------------------------------------------------------
!
! allocate arrays
!jjr      allocate ( tmp2d(NX,   beglat:endlat) )
      allocate ( tmp2d(IB:IE,   beglat:endlat) )
      allocate ( temp3d(plon,beglev:endlev,beglat:endlat) )
      allocate ( tmp3d(IB:IE,beglev:endlev,beglat:endlat) )
      allocate ( tmp3dz(IB:IE,beglev:endlevp,beglat:endlat) )
!jjr      allocate ( tmp3d(NX,NL,beglat:endlat) )
!jjr      allocate ( tmp3dz(NX,NZ,beglat:endlat) )
!
     num=plon*plat
     call wrtout(nrg, strip2d, phis, num, 2)

! Write fields for cam physics
!
!jjr      call wrtout_r8 (nrg,phis  ,plon )
      num = plnlv*plat
  call wrtout(nrg, strip3dxzy, omga, num, 3)

!jjr      call wrtout_r8 (nrg,omga  ,plnlv)
!
! Write fields u3,v3,t3,q3,ps at time n3m2 & n3m1
!
      call wrtout(nrg, strip3dxzy, u3(:,:,:,n3m1), num, 3)
      call wrtout(nrg, strip3dxzy, v3(:,:,:,n3m1), num, 3)
      call wrtout(nrg, strip3dxzy, t3(:,:,:,n3m1), num, 3)
!      call wrtout_r8 (nrg,u3(1,1,beglat,n3m1)  ,plnlv)
!      call wrtout_r8 (nrg,v3(1,1,beglat,n3m1)  ,plnlv)
!      call wrtout_r8 (nrg,t3(1,1,beglat,n3m1)  ,plnlv)
      num = plon*plat
      call wrtout(nrg, strip2d, ps(:,:,n3m1), num, 2)

!      call wrtout_r8 (nrg,ps(1,beglat,n3m1)  ,plon)
       num = plnlv*plat
      call wrtout(nrg, strip3dxzy, u3(:,:,:,n3m2), num, 3)
      call wrtout(nrg, strip3dxzy, v3(:,:,:,n3m2), num, 3)
      call wrtout(nrg, strip3dxzy, t3(:,:,:,n3m2), num, 3)

!      call wrtout_r8 (nrg,u3(1,1,beglat,n3m2)  ,plnlv)
!      call wrtout_r8 (nrg,v3(1,1,beglat,n3m2)  ,plnlv)
!      call wrtout_r8 (nrg,t3(1,1,beglat,n3m2)  ,plnlv)
      num = plon*plat
      call wrtout(nrg, strip2d, ps(:,:,n3m2), num, 2)
        num = plnlv*plat
!      call wrtout_r8 (nrg,ps(1,beglat,n3m2)  ,plon)
         do m=1,pcnst+pnats
             do j=beglat,endlat
                do k=beglev,endlev
                     do i=1,plon
                    temp3d(i,k,j)=q3(i,k,m,j,n3m2)
                     end do
                 end do 

              end do
      call wrtout(nrg, strip3dxzy, temp3d, num, 3)

!      call wrtout_r8 (nrg,q3(1,1,1,beglat,n3m2),plnlv*(pcnst+pnats))
        end do
         do m=1,pcnst+pnats
             do j=beglat,endlat
                do k=beglev,endlev
                     do i=1,plon
                    temp3d(i,k,j)=q3(i,k,m,j,n3m1)
                     end do
                 end do
              end do
      call wrtout(nrg, strip3dxzy, temp3d, num, 3)

!      call wrtout_r8 (nrg,q3(1,1,1,beglat,n3m1),plnlv*(pcnst+pnats))
      enddo
      if (cnst_need_pdeldry) then
   call wrtout(nrg, strip3dxzy, pdeld(:,:,:,n3), num, 3)
   call wrtout(nrg, strip3dxzy, pdeld(:,:,:,n3m1), num, 3)
   call wrtout(nrg, strip3dxzy, pdeld(:,:,:,n3m2), num, 3)

!         call wrtout_r8 (nrg,pdeld(1,1,beglat,n3)  ,plnlv)
!         call wrtout_r8 (nrg,pdeld(1,1,beglat,n3m1)  ,plnlv)
!         call wrtout_r8 (nrg,pdeld(1,1,beglat,n3m2)  ,plnlv)
      endif !
!
! Write slt arrays (trajectory mid-point coordinates and 
! slt forcast of moisture and constituents
!
      num = plnlv*plat
    call wrtout(nrg, strip3dxzy, lammp, num, 3)
    call wrtout(nrg, strip3dxzy, phimp, num, 3)
    call wrtout(nrg, strip3dxzy, sigmp, num, 3)
      do m=1,pcnst
             do j=beglat,endlat
                do k=beglev,endlev
                     do i=1,plon
                    temp3d(i,k,j)=qfcst(i,k,m,j)
                     end do
                 end do
              end do
      call wrtout(nrg, strip3dxzy, temp3d, num, 3)

!      call wrtout_r8 (nrg,lammp,plnlv)
!      call wrtout_r8 (nrg,phimp,plnlv)
!      call wrtout_r8 (nrg,sigmp,plnlv)
!      call wrtout_r8 (nrg,qfcst,plnlv*pcnst)
!
      end do
! Write dynamical fields for restart
!
! for 2-D
! =================== zhh 2008.7.21 ===================
      do j = beglat, endlat
		 tmp2d(IB:IE,j) = P(IB:IE,plat+1-j)   
      end do
     num=plon*plat
     call wrtout(nrg, strip2d, tmp2d, num, 2)


!jjr      call wrtout_r8 (nrg, tmp2d(1,beglat), NX)
!
      do j = beglat, endlat
		 tmp2d(IB:IE,j) = PT(IB:IE,plat+1-j)   
      end do
     call wrtout(nrg, strip2d, tmp2d, num, 2)

!      call wrtout_r8 (nrg, tmp2d(1,beglat), NX)
!
      do j = beglat, endlat
		 tmp2d(IB:IE,j) = Pstar1(IB:IE,plat+1-j)   
      end do
     call wrtout(nrg, strip2d, tmp2d, num, 2)

!      call wrtout_r8 (nrg, tmp2d(1,beglat), NX)
!
      do j = beglat, endlat
		 tmp2d(IB:IE,j) = Psa(IB:IE,plat+1-j)   
      end do
     call wrtout(nrg, strip2d, tmp2d, num, 2)

!      call wrtout_r8 (nrg, tmp2d(1,beglat), NX)
!
      do j = beglat, endlat
		 tmp2d(IB:IE,j) = PS2(IB:IE,plat+1-j)   
      end do
     call wrtout(nrg, strip2d, tmp2d, num, 2)

!      call wrtout_r8 (nrg, tmp2d(1,beglat), NX)
!
      do j = beglat, endlat
		 tmp2d(IB:IE,j) = GHS(IB:IE,plat+1-j)   
      end do
!      call wrtout_r8 (nrg, tmp2d(1,beglat), NX)
     call wrtout(nrg, strip2d, tmp2d, num, 2)

!
! for 3-D
      do j = beglat, endlat
		 tmp3d(IB:IE,beglev:endlev,j) = U(IB:IE,beglev:endlev,plat+1-j)   
      end do
      num = plnlv*plat
  call wrtout(nrg, strip3dxzy, tmp3d, num, 3)

!      call wrtout_r8 (nrg, tmp3d(1,1,beglat), NX*NL)
!
      do j = beglat, endlat
		 tmp3d(IB:IE,beglev:endlev,j) = V(IB:IE,beglev:endlev,plat+1-j)   
      end do
  call wrtout(nrg, strip3dxzy, tmp3d, num, 3)

!      call wrtout_r8 (nrg, tmp3d(1,1,beglat), NX*NL)
!
      do j = beglat, endlat
		 tmp3d(IB:IE,:,j) = T(IB:IE,:,plat+1-j)   
      end do
  call wrtout(nrg, strip3dxzy, tmp3d, num, 3)

!      call wrtout_r8 (nrg, tmp3d(1,1,beglat), NX*NL)
!
      do j = beglat, endlat
		 tmp3d(IB:IE,:,j) = Q(IB:IE,:,plat+1-j)   
      end do
  call wrtout(nrg, strip3dxzy, tmp3d, num, 3)

!      call wrtout_r8 (nrg, tmp3d(1,1,beglat), NX*NL)
!
      do j = beglat, endlat
		 tmp3d(IB:IE,:,j) = UT(IB:IE,:,plat+1-j)   
      end do
  call wrtout(nrg, strip3dxzy, tmp3d, num, 3)

!      call wrtout_r8 (nrg, tmp3d(1,1,beglat), NX*NL)
!
      do j = beglat, endlat
		 tmp3d(IB:IE,:,j) = VT(IB:IE,:,plat+1-j)   
      end do
  call wrtout(nrg, strip3dxzy, tmp3d, num, 3)

!      call wrtout_r8 (nrg, tmp3d(1,1,beglat), NX*NL)
!
      do j = beglat, endlat
		 tmp3d(IB:IE,:,j) = TT(IB:IE,:,plat+1-j)   
      end do
  call wrtout(nrg, strip3dxzy, tmp3d, num, 3)

!      call wrtout_r8 (nrg, tmp3d(1,1,beglat), NX*NL)
!
      num=plevp*plon*plat   
      do j = beglat, endlat
	 tmp3dz(IB:IE,beglev:endlevp,j) = WS(IB:IE,beglev:endlevp,plat+1-j)   
      end do
     call wrtout(nrg, strip3dxzyp, tmp3dz, num, 3)

!      call wrtout_r8 (nrg, tmp3dz(1,1,beglat), NX*NZ)
! =================== zhh 2008.7.21 ===================
!
! Write global integrals
!
      if (masterproc) then
         write(nrg, iostat=ioerr) tmass0, fixmas, hw1,    hw2,  &
                                  hw3, alpha
         if (ioerr /= 0 ) then
            write (6,*) 'WRITE ioerror ',ioerr,' on i/o unit = ',nrg
            call endrun ('WRITE_RESTART_DYNAMICS')
         end if
      end if

#if ( defined BFB_CAM_SCAM_IOP )
!
! Write scam values
!
!jjr no need     call wrtout_r8 (nrg,alphasav(1,beglat),pcnst)
!     call wrtout_r8 (nrg,dqfx3savm1(1,1,1,beglat),plnlv*pcnst)       
!     call wrtout_r8 (nrg,divq3dsav(1,1,1,beglat),plnlv*ppcnst)
!     call wrtout_r8 (nrg,divt3dsav(1,1,beglat),plnlv)       
!     call wrtout_r8 (nrg,t3sav(1,1,beglat),plnlv)       
!     call wrtout_r8 (nrg,u3sav(1,1,beglat),plnlv)
!     call wrtout_r8 (nrg,v3sav(1,1,beglat),plnlv)
!     call wrtout_r8 (nrg,t2sav(1,1,beglat),plnlv)
!     call wrtout_r8 (nrg,q3sav(1,1,1,beglat),plnlv*ppcnst)
!     call wrtout_r8 (nrg,pssav(1,beglat),plon)
!     call wrtout_r8 (nrg,tssav(1,beglat),plon)
!     call wrtout_r8 (nrg,fixmassav(beglat),1)
!     call wrtout_r8 (nrg,betasav(beglat),1)
#endif
! ============================ for test =================================
      do j = beglat, endlat
         if ( j == 3.and.beglev.le.7.and.endlev.ge.7 ) then
            jd = plat+1-j
            write(6,*) '--------------- At sub. write_restart_dynamics --------------------'
!            write(6,*) 'UT(14,7,180) =', UT(14,7,180), ' T(14,7,180) =', T(14,7,180)
!            write(6,*) 'VT(14,7,180) =', VT(14,7,180), 'TT(14,7,180)=',TT(14,7,180)
!            write(6,*) 'Q(14,7,180) =', Q(14,7,180), ' ps(14,2,n3m1) =', ps(14,2,n3m1)
!            write(6,*) 'p(14,180) =', P(14,180), ' pT(14,180) =', pT(14,180)
!            write(6,*) 'v3(14,7,2,n3m1) =', v3(14,7,2,n3m1)
!            write(6,*) 'q3(14,7,2,2,n3m1) =', q3(14,7,2,2,n3m1)

            write(6,*) 'jd =', jd
            write(6,*) 'UT(4,7,jd) =', UT(4,7,jd), ' T(4,7,jd) =', T(4,7,jd)
            write(6,*) 'VT(4,7,jd) =', VT(4,7,jd), 'TT(4,7,jd)=',TT(4,7,jd)
            write(6,*) 'Q(4,7,jd) =', Q(4,7,jd), ' ps(4,3,n3m1) =', ps(4,3,n3m1)
            write(6,*) 'v3(4,7,3,n3m1) =', v3(4,7,3,n3m1)
            write(6,*) 'q3(4,7,2,3,n3m1) =', q3(4,7,2,3,n3m1)
            write(6,*) '-----------------------------------'

         end if
!         write(6,*) 'VT(4,1,',jd,') =', VT(4,1,jd)
!         write(6,*) 'PT(4,',jd,') =', PT(4,jd)
      end do
! =============================== zhh ===================================
! deallocate
      deallocate(temp3d)
      deallocate(tmp2d)
      deallocate(tmp3d)
      deallocate(tmp3dz)
!
      return
   end subroutine write_restart_dynamics
   subroutine wrtout(iu, decomp, arr, lenarr, ndim)
!-----------------------------------------------------------------------
! Wrapper routine to write restart binary file
!-----------------------------------------------------------------------
      use shr_kind_mod, only: r8 => shr_kind_r8
      use pmgrid
      use decompmodule, only: decomptype
#if ( defined SPMD )
      use spmd_dyn, only: comm_y
      use parutilitiesmodule, only: commglobal, pargather
#endif
!------------------------------Arguments--------------------------------
      integer iu                 ! Logical unit
      type (decomptype):: decomp ! Decomposition descriptor
      integer lenarr             ! Global length of array
#if defined( SPMD )
      real(r8) arr(*)            ! Array to be gathered
#else
      real(r8) arr(lenarr)       ! Array (SMP-only)
#endif
      integer ndim               ! dimensionality (2 or 3) of array
!---------------------------Local variables-----------------------------
      integer ioerr              ! errorcode
#if ( defined SPMD )
      real(r8), allocatable :: bufres(:)
#endif
!-----------------------------------------------------------------------
#if ( defined SPMD )
      if ( masterproc ) then
         allocate( bufres(lenarr) )
      else
         allocate( bufres(1) )
      endif
      if (ndim .eq. 2 .and. twod_decomp .eq. 1) then
         if (myid_z .eq. 0) call pargather( comm_y, 0, arr, decomp, bufres )
      else
         call pargather( commglobal, 0, arr, decomp, bufres )
      endif

! WS 01.01.03: This code is OK
      if (masterproc) then
         write (iu,iostat=ioerr) bufres
         if (ioerr /= 0 ) then
            write (6,*) 'WRTOUT ioerror ',ioerr,' on i/o unit = ',iu
            call endrun ('WRTOUT')
         end if
      endif
      deallocate( bufres )
#else
      write (iu,iostat=ioerr) arr
      if (ioerr /= 0 ) then
         write (6,*) 'wrt ioerror ',ioerr,' on i/o unit = ',iu
         call endrun ('WRTOUT')
      end if
#endif
      return
   end subroutine wrtout
   subroutine lrreadin(iu, decomp, arr, lenarr, ndim)
!-----------------------------------------------------------------------
! Wrapper routine to read real variable from restart binary file
!-----------------------------------------------------------------------
      use shr_kind_mod, only: r8 => shr_kind_r8
      use pmgrid
      use decompmodule, only : decomptype
#if ( defined SPMD )
      use spmd_dyn, only: comm_y, comm_z
      use parutilitiesmodule, only : commglobal, parscatter, parcollective, BCSTOP
#endif
!------------------------------Arguments--------------------------------
      integer iu                 ! Logical unit
      type (decomptype):: decomp ! Decomposition descriptor
      integer lenarr             ! Global size of array
#if defined( SPMD )
      real(r8) arr(*)            ! Array to be gathered
#else
      real(r8) arr(lenarr)       ! Array (SMP-only)
#endif
      integer ndim               ! dimensionality (2 or 3) of array
!---------------------------Local variables-----------------------------
      integer ioerr              ! errorcode
#ifdef SPMD
      real(r8), allocatable :: bufres(:)
#endif
!-----------------------------------------------------------------------
#ifdef SPMD
      if (masterproc) then
         allocate (bufres(lenarr))
         read (iu,iostat=ioerr) bufres
         if (ioerr /= 0 ) then
            write (6,*) 'LRREADIN ioerror ',ioerr,' on i/o unit = ',iu
            call endrun
         end if
      else
         allocate (bufres(1))
      endif
      if (ndim .eq. 2 .and. twod_decomp .eq. 1) then
         if (myid_z .eq. 0) call parscatter( comm_y, 0, bufres, decomp, arr )
         call parcollective( comm_z, BCSTOP, plon*(endlat-beglat+1), arr )
      else
         call parscatter( commglobal, 0, bufres, decomp, arr )
      endif
      deallocate (bufres)
#else
      read (iu,iostat=ioerr) arr
      if (ioerr /= 0 ) then
         write (6,*) 'LRREADIN ioerror ',ioerr,' on i/o unit = ',iu
         call endrun
      end if
#endif
      return
   end subroutine lrreadin



!#######################################################################

   subroutine read_restart_dynamics (nrg)
!jjr
      use dynamics_vars, only: dynamics_init
      use time_manager, only: get_step_size
!jjr
#if ( defined SPMD )
      use mpishorthand
#endif

#include <comqfl.h>
#include <comctl.h>

!
! Input arguments
!
      integer, intent(in) :: nrg     ! Unit number
!
! Local workspace
!
      integer :: ioerr   ! error status
      integer :: k,i,j,m,jd,num        !zhh 2007.12.20
      real(r8):: dtime
! =================== zhh 2008.7.21 ===================
      real(r8), allocatable :: tmp2d(:,:)
      real(r8), allocatable :: tmp3d(:,:,:)
      real(r8), allocatable :: temp3d(:,:,:)
      real(r8), allocatable :: tmp3dz(:,:,:)
      integer:: src3,dest3
! =================== zhh 2008.7.21 ===================
!--------------------------------------------------------------------
!
! allocate arrays
      allocate ( tmp2d(IB:IE,   beglat:endlat) )
      allocate ( tmp3d(IB:IE,beglev:endlev,beglat:endlat) )
      allocate ( temp3d(plon,beglev:endlev,beglat:endlat) )
      allocate ( tmp3dz(IB:IE,beglev:endlevp,beglat:endlat) )
!
! Read fields for cam physics
      dtime = get_step_size()
      call dynamics_init( dtime, iord, jord, nsplit,  &
                          plon, plat, plev, ppcnst,   &
                          beglonxy, endlonxy,         &
                          beglatxy, endlatxy,         &
                          beglat,   endlat,           &
                          beglev,   endlev )
      call initialize_prognostics
      call initialize_IAPprog     !zhh 2008.6.10
      num = plon*plat
      call lrreadin(nrg, strip2d, phis, num, 2)

      num = plnlv*plat
      call lrreadin(nrg, strip3dxzy, omga, num, 3)
!      call readin_r8 (nrg,phis  ,plon )
!      call readin_r8 (nrg,omga  ,plnlv)
!
! Read fields u3,v3,t3,q3,ps at time indices n3m2 & n3m1
!
    call lrreadin(nrg, strip3dxzy, u3(:,:,:,n3m1), num, 3)
    call lrreadin(nrg, strip3dxzy, v3(:,:,:,n3m1), num, 3)
    call lrreadin(nrg, strip3dxzy, t3(:,:,:,n3m1), num, 3)

!      call readin_r8 (nrg,u3(1,1,beglat,n3m1)  ,plnlv)
!      call readin_r8 (nrg,v3(1,1,beglat,n3m1)  ,plnlv)
!      call readin_r8 (nrg,t3(1,1,beglat,n3m1)  ,plnlv)
      num = plon*plat
      call lrreadin(nrg, strip2d, ps(:,:,n3m1), num, 2)

!      call readin_r8 (nrg,ps(1,beglat,n3m1)  ,plon)
       num = plnlv*plat
    call lrreadin(nrg, strip3dxzy, u3(:,:,:,n3m2), num, 3)
    call lrreadin(nrg, strip3dxzy, v3(:,:,:,n3m2), num, 3)
    call lrreadin(nrg, strip3dxzy, t3(:,:,:,n3m2), num, 3)
 
!      call readin_r8 (nrg,u3(1,1,beglat,n3m2)  ,plnlv)
!      call readin_r8 (nrg,v3(1,1,beglat,n3m2)  ,plnlv)
!      call readin_r8 (nrg,t3(1,1,beglat,n3m2)  ,plnlv)
      num = plon*plat
      call lrreadin(nrg, strip2d, ps(:,:,n3m2), num, 2)

!      call readin_r8 (nrg,ps(1,beglat,n3m2)  ,plon)
       num = plnlv*plat
         do m=1,pcnst+pnats
    call lrreadin(nrg, strip3dxzy, temp3d, num, 3)
             do j=beglat,endlat
                do k=beglev,endlev
                     do i=1,plon
                    q3(i,k,m,j,n3m2)=temp3d(i,k,j)
                     end do
                 end do
              end do
           end do
         do m=1,pcnst+pnats
    call lrreadin(nrg, strip3dxzy, temp3d, num, 3)
             do j=beglat,endlat
                do k=beglev,endlev
                     do i=1,plon
                    q3(i,k,m,j,n3m1)=temp3d(i,k,j)
                     end do
                 end do
              end do
           end do
!      call readin_r8 (nrg,q3(1,1,1,beglat,n3m2),plnlv*(pcnst+pnats))
!      call readin_r8 (nrg,q3(1,1,1,beglat,n3m1),plnlv*(pcnst+pnats))
      if (cnst_need_pdeldry) then
    call lrreadin(nrg, strip3dxzy, pdeld(:,:,:,n3), num, 3)
    call lrreadin(nrg, strip3dxzy, pdeld(:,:,:,n3m1), num, 3)
    call lrreadin(nrg, strip3dxzy, pdeld(:,:,:,n3m2), num, 3)

!         call readin_r8 (nrg,pdeld(1,1,beglat,n3)  ,plnlv)
!         call readin_r8 (nrg,pdeld(1,1,beglat,n3m1)  ,plnlv)
!         call readin_r8 (nrg,pdeld(1,1,beglat,n3m2)  ,plnlv)
      endif !
!
! Write slt arrays (trajectory mid-point coordinates and 
! slt forcast of moisture and constituents
!
      call scanslt_alloc()
    call lrreadin(nrg, strip3dxzy,lammp , num, 3)
    call lrreadin(nrg, strip3dxzy,phimp , num, 3)
    call lrreadin(nrg, strip3dxzy,sigmp ,num, 3)

!      call readin_r8 (nrg,lammp,plnlv)
!      call readin_r8 (nrg,phimp,plnlv)
!      call readin_r8 (nrg,sigmp,plnlv)
       do m=1,pcnst
    call lrreadin(nrg, strip3dxzy,temp3d , num, 3)
             do j=beglat,endlat
                do k=beglev,endlev
                     do i=1,plon
                    qfcst(i,k,m,j)=temp3d(i,k,j)
                     end do
                 end do
              end do
       end do    
!      call readin_r8 (nrg,qfcst,plnlv*pcnst)
!
! Read dynamical fields
! for 2-D
      num = plon*plat
      call lrreadin(nrg, strip2d, tmp2d, num, 2)

!      call readin_r8 (nrg, tmp2d(1,beglat), NX)
      do j = beglat, endlat
		 P(IB:IE,plat+1-j) = tmp2d(IB:IE,j)   
       call  period(P(1,plat+1-j))
      end do
      call lrreadin(nrg, strip2d, tmp2d, num, 2)
!
!      call readin_r8 (nrg, tmp2d(1,beglat), NX)
      do j = beglat, endlat
		 PT(IB:IE,plat+1-j) = tmp2d(IB:IE,j)   
       call  period(PT(1,plat+1-j))
      end do
!

      call lrreadin(nrg, strip2d, tmp2d, num, 2)
!      call readin_r8 (nrg, tmp2d(1,beglat), NX)
      do j = beglat, endlat
		 Pstar1(IB:IE,plat+1-j) = tmp2d(IB:IE,j)   
           call  period(Pstar1(1,plat+1-j))
      end do
!
      call lrreadin(nrg, strip2d, tmp2d, num, 2)
!      call readin_r8 (nrg, tmp2d(1,beglat), NX)
      do j = beglat, endlat
		 Psa(IB:IE,plat+1-j) = tmp2d(IB:IE,j)   
          call  period(Psa(1,plat+1-j))
      end do
!
      call lrreadin(nrg, strip2d, tmp2d, num, 2)
!      call readin_r8 (nrg, tmp2d(1,beglat), NX)
      do j = beglat, endlat
		 PS2(IB:IE,plat+1-j) = tmp2d(IB:IE,j)   
           call  period(PS2(1,plat+1-j))
      end do
!
      call lrreadin(nrg, strip2d, tmp2d, num, 2)
!      call readin_r8 (nrg, tmp2d(1,beglat), NX)
      do j = beglat, endlat
		 GHS(IB:IE,plat+1-j) = tmp2d(IB:IE,j)   
            call  period(GHS(1,plat+1-j))

      end do
!

! for 3-D
      num = plnlv*plat
      call lrreadin(nrg, strip3dxzy, tmp3d, num, 3)

!      call readin_r8 (nrg, tmp3d(1,1,beglat), NX*NL)
      do j = beglat, endlat
		 U(IB:IE,:,plat+1-j) = tmp3d(IB:IE,:,j)   
            do k=beglev,endlev
          call period(U(1,k,plat+1-j))
            end do
      end do
!
      call lrreadin(nrg, strip3dxzy, tmp3d, num, 3)
!      call readin_r8 (nrg, tmp3d(1,1,beglat), NX*NL)
      do j = beglat, endlat
		 V(IB:IE,:,plat+1-j) = tmp3d(IB:IE,:,j)   
           do k=beglev,endlev
          call period(V(1,k,plat+1-j))
            end do

      end do
!
      call lrreadin(nrg, strip3dxzy, tmp3d, num, 3)
!      call readin_r8 (nrg, tmp3d(1,1,beglat), NX*NL)
      do j = beglat, endlat
		 T(IB:IE,:,plat+1-j) = tmp3d(IB:IE,:,j)   
                  do k=beglev,endlev
          call period(T(1,k,plat+1-j))
            end do

      end do
!
      call lrreadin(nrg, strip3dxzy, tmp3d, num, 3)
!      call readin_r8 (nrg, tmp3d(1,1,beglat), NX*NL)
      do j = beglat, endlat
		 Q(IB:IE,:,plat+1-j) = tmp3d(IB:IE,:,j)   
                   do k=beglev,endlev
          call period(Q(1,k,plat+1-j))
            end do

      end do
!
      call lrreadin(nrg, strip3dxzy, tmp3d, num, 3)
!      call readin_r8 (nrg, tmp3d(1,1,beglat), NX*NL)
      do j = beglat, endlat
		 UT(IB:IE,:,plat+1-j) = tmp3d(IB:IE,:,j)   
           do k=beglev,endlev
          call period(UT(1,k,plat+1-j))
            end do

      end do
!
      call lrreadin(nrg, strip3dxzy, tmp3d, num, 3)
!      call readin_r8 (nrg, tmp3d(1,1,beglat), NX*NL)
      do j = beglat, endlat
		 VT(IB:IE,:,plat+1-j) = tmp3d(IB:IE,:,j)   
            do k=beglev,endlev
          call period(VT(1,k,plat+1-j))
            end do

      end do
!
      call lrreadin(nrg, strip3dxzy, tmp3d, num, 3)
!      call readin_r8 (nrg, tmp3d(1,1,beglat), NX*NL)
      do j = beglat, endlat
		 TT(IB:IE,:,plat+1-j) = tmp3d(IB:IE,:,j)   
            do k=beglev,endlev
          call period(TT(1,k,plat+1-j))
            end do

      end do
!
         num=plevp*plon*plat
      call lrreadin(nrg, strip3dxzyp, tmp3dz, num, 3)
!      call readin_r8 (nrg, tmp3dz(1,1,beglat), NX*NZ)
      do j = beglat, endlat
	 WS(IB:IE,beglev:endlevp,plat+1-j) = tmp3dz(IB:IE,beglev:endlevp,j)   
                 do k=beglev,endlevp
          call period(WS(1,k,plat+1-j))
            end do

      end do
!test jjr
! Read global integrals
!
!down
 !test jjr
!      q3=0.
    if(npr_z.gt.1) then
    src3=iam+npr_y
    dest3=iam-npr_y
    if(myid_z.eq.0) dest3=-1
    if(myid_z.eq.(npr_z-1)) src3=-1
    call mp_send3d(dest3 , src3, NX,NZ,NY,                     &
                      1, NX,beglev,endlev+1,beglatdynex,endlatdynex,        &
                      1, NX,beglev,beglev, beglatdynex,endlatdynex,WS )
    call mp_recv3d(src3, NX,NZ,NY,                     &
                      1, NX,beglev,endlev+1,beglatdynex,endlatdynex,        &
                      1, NX,endlev+1,endlev+1,beglatdynex,endlatdynex,WS )
   endif

       if(iam.eq.32) print*,'test WS1',ws(4,beglev,beglatdyn+1)
      if (masterproc) then
         read (nrg, iostat=ioerr) tmass0, fixmas, hw1,    hw2,  &
                                  hw3, alpha
         if (ioerr /= 0 ) then
            write (6,*) 'WRITE ioerror ',ioerr,' on i/o unit = ',nrg
            call endrun ('READ_RESTART_DYNAMICS')
         end if
         print*,'test WS',ws(4,endlev+1,beglatdyn+1),npr_z
      end if
!
! ============================ for test =================================
      do j = beglat, endlat
		 if ( j == 3.and.beglev.le.7.and.endlev.ge.7 ) then
            jd = plat+1-j
            write(6,*) '--------------- At sub. read_restart_dynamics --------------------'
            write(6,*) 'jd =', jd
            write(6,*) 'UT(4,7,jd) =', UT(4,7,jd), ' T(4,7,jd) =', T(4,7,jd)
            write(6,*) 'VT(4,7,jd) =', VT(4,7,jd), 'TT(4,7,jd)=',TT(4,7,jd)
            write(6,*) 'Q(4,7,jd) =', Q(4,7,jd), ' ps(4,3,n3m1) =', ps(4,3,n3m1)
            write(6,*) 'v3(4,7,3,n3m1) =', v3(4,7,3,n3m1)
            write(6,*) 'q3(4,7,2,3,n3m1) =', q3(4,7,2,3,n3m1)
            write(6,*) '--------------------------------------------------------------'
         end if
!         write(6,*) 'VT(4,1,',jd,') =', VT(4,1,jd)
!         write(6,*) 'PT(4,',jd,') =', PT(4,jd)
      end do
! =============================== zhh ===================================
!
#if ( defined SPMD )
!===================== delete by zhh ===============================
   call mpibcast (tmass0,1         ,mpir8  ,0,mpicom)      
   call mpibcast (fixmas,1         ,mpir8  ,0,mpicom)
   call mpibcast (hw1   ,pcnst     ,mpir8  ,0,mpicom)
   call mpibcast (hw2   ,pcnst     ,mpir8  ,0,mpicom)
   call mpibcast (hw3   ,pcnst     ,mpir8  ,0,mpicom)   
   call mpibcast (alpha ,pcnst     ,mpir8  ,0,mpicom)
!===================== delete by zhh ===============================
#endif
#if ( defined BFB_CAM_SCAM_IOP )
!
! Read scam values
!
!jjr     call init_iop_fields(ps, t3, u3, v3, q3, nocopy=.true. )

!     call readin_r8 (nrg,alphasav(1,beglat),pcnst)
!     call readin_r8 (nrg,dqfx3savm1(1,1,1,beglat),plnlv*pcnst)       
!     call readin_r8 (nrg,divq3dsav(1,1,1,beglat),plnlv*ppcnst)
!     call readin_r8 (nrg,divt3dsav(1,1,beglat),plnlv)       
!     call readin_r8 (nrg,t3sav(1,1,beglat),plnlv)       
!     call readin_r8 (nrg,u3sav(1,1,beglat),plnlv)
!     call readin_r8 (nrg,v3sav(1,1,beglat),plnlv)
!     call readin_r8 (nrg,t2sav(1,1,beglat),plnlv)
!     call readin_r8 (nrg,q3sav(1,1,1,beglat),plnlv*ppcnst)
!     call readin_r8 (nrg,pssav(1,beglat),plon)
!     call readin_r8 (nrg,tssav(1,beglat),plon)
!     call readin_r8 (nrg,fixmassav(beglat),1)
!     call readin_r8 (nrg,betasav(beglat),1)
#endif
!
! deallocate
      deallocate(tmp2d)
      deallocate(tmp3d)
      deallocate(tmp3dz)
      deallocate(temp3d)
!
      return

   end subroutine read_restart_dynamics
   subroutine lrreadin4(iu, decomp, arr, lenarr, ndim)
!-----------------------------------------------------------------------
! Wrapper routine to read real*4 variable from restart binary file
!-----------------------------------------------------------------------
      use shr_kind_mod, only: r8 => shr_kind_r8, r4 => shr_kind_r4
      use pmgrid
      use decompmodule, only : decomptype
#if ( defined SPMD )
      use spmd_dyn, only: comm_y, comm_z
      use parutilitiesmodule, only : commglobal, parscatterreal4, parcollective1dreal4, BCSTOP
#endif
!------------------------------Arguments--------------------------------
      integer iu                 ! Logical unit
      type (decomptype):: decomp ! Decomposition descriptor
      integer lenarr             ! Global size of array
#if defined( SPMD )
      real(r4) arr(*)            ! Array to be gathered
#else
      real(r4) arr(lenarr)       ! Array (SMP-only)
#endif
      integer ndim               ! dimensionality (2 or 3) of array
!---------------------------Local variables-----------------------------
      integer ioerr              ! errorcode
#ifdef SPMD
      real(r4), allocatable :: bufres(:)
#endif
!-----------------------------------------------------------------------
#ifdef SPMD
      if (masterproc) then
         allocate (bufres(lenarr))
         read (iu,iostat=ioerr) bufres
         if (ioerr /= 0 ) then
            write (6,*) 'LRREADIN4 ioerror ',ioerr,' on i/o unit = ',iu
            call endrun
         end if
      else
         allocate (bufres(1))
      endif
      if (ndim .eq. 2 .and. twod_decomp .eq. 1) then
         if (myid_z .eq. 0) call parscatterreal4( comm_y, 0, bufres, decomp, arr )
         call parcollective1dreal4( comm_z, BCSTOP, plon*(endlat-beglat+1), arr )
      else
         call parscatterreal4( commglobal, 0, bufres, decomp, arr )
      endif
      deallocate (bufres)
#else
      read (iu,iostat=ioerr) arr
      if (ioerr /= 0 ) then
         write (6,*) 'LRREADIN4 ioerror ',ioerr,' on i/o unit = ',iu
         call endrun
      end if
#endif
      return
   end subroutine lrreadin4

   subroutine lrreadini(iu, decomp, arr, lenarr, ndim)
!-----------------------------------------------------------------------
! Wrapper routine to read integer variable from restart binary file
!-----------------------------------------------------------------------
      use shr_kind_mod, only: r8 => shr_kind_r8
      use pmgrid
      use decompmodule, only : decomptype
#if ( defined SPMD )
      use spmd_dyn, only: comm_y, comm_z
      use parutilitiesmodule, only : commglobal, parscatter, parcollective, BCSTOP
#endif
!------------------------------Arguments--------------------------------
      integer iu                 ! Logical unit
      type (decomptype):: decomp ! Decomposition descriptor
      integer lenarr             ! Global size of array
#if defined( SPMD )
      integer arr(*)             ! Array to be gathered
#else
      integer arr(lenarr)        ! Array (SMP-only)
#endif
      integer ndim               ! dimensionality (2 or 3) of array
!---------------------------Local variables-----------------------------
      integer ioerr              ! errorcode
#ifdef SPMD
      integer, allocatable :: bufres(:)
#endif
!-----------------------------------------------------------------------
#ifdef SPMD
      if (masterproc) then
         allocate (bufres(lenarr))
         read (iu,iostat=ioerr) bufres
         if (ioerr /= 0 ) then
            write (6,*) 'LRREADINI ioerror ',ioerr,' on i/o unit = ',iu
            call endrun
         end if
      else
         allocate (bufres(1))
      endif
      if (ndim .eq. 2 .and. twod_decomp .eq. 1) then
         if (myid_z .eq. 0) call parscatter( comm_y, 0, bufres, decomp, arr )
         call parcollective( comm_z, BCSTOP, plon*(endlat-beglat+1), arr )
      else
         call parscatter( commglobal, 0, bufres, decomp, arr )
      endif
      deallocate (bufres)
#else
      read (iu,iostat=ioerr) arr
      if (ioerr /= 0 ) then
         write (6,*) 'LRREADINI ioerror ',ioerr,' on i/o unit = ',iu
         call endrun
      end if
#endif
      return
   end subroutine lrreadini


end module restart_dynamics
