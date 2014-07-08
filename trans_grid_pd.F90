#include <misc.h>
#include <params.h>

subroutine trans_grid_pd( fu, fv, t2, qminus, beglat, endlat, pcnst )
!----------------------------------------------------------------------------------
! Purpose: Transform tendencies from physics grid to dynamical grid.
! Author : ZhangHe
! Completed: 2007.4.11
! Update : ZhangHe, 2007.10.11
!          ZhangHe, 2007.12.21
!          WuJianping, 2008.5, parallel version
!          ZhangHe, 2008.6.11, available for both serial & parallel
!----------------------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use IAP_grid,  only: NX, NY, NL, EX, IB, IE, period
   use pmgrid,    only: plon, plat, plev, beglatdyn, endlatdyn, iam,beglev,endlev,beglatdynex,endlatdynex
   use IAP_prog,  only: Q, Qliq, Qice
   use tendency,  only: SU, SV, ST
   use Dyn_const, only: DXPVN, DXPVS, PMTOP   
   use mathconst, only: ZERO, HALF
#if (defined SPMD)
      use mod_comm, only: mp_send3d, mp_recv3d
      use pmgrid,     only : npr_y,iam
#endif

   implicit none
!----------------------------------------Arguments-----------------------------------------------
   integer,  intent(in) :: beglat
   integer,  intent(in) :: endlat
   integer,  intent(in) :: pcnst
   real(r8), intent(in) :: fu (plon,beglev:endlev,beglat:endlat)
   real(r8), intent(in) :: fv (plon,beglev:endlev,beglat:endlat)
   real(r8), intent(in) :: t2 (plon,beglev:endlev,beglat:endlat)
   real(r8), intent(in) :: qminus(plon,beglev:endlev,pcnst,beglat:endlat)
!-------------------------------------Local workspace--------------------------------------------
   real(r8), allocatable :: futmp(:,:,:)    ! u wind tendency
   real(r8), allocatable :: fvtmp(:,:,:)    ! v wind tendency
   integer  :: I, Jp, Jd, K, J,dest,src
   
!-----------------------------------------------------------------------------
!
   allocate(futmp(NX,beglev:endlev,beglatdynex:endlatdynex))
   allocate(fvtmp(NX,beglev:endlev,beglatdynex:endlatdynex))
!
   DO Jp = beglat, endlat
      Jd = plat + 1 - Jp
      DO K = beglev,endlev 
         DO I = 1, plon
            futmp(I+EX,K,Jd) = fu(I,K,Jp)
            fvtmp(I+EX,K,Jd) =-fv(I,K,Jp)
            ST   (I+EX,K,Jd) = t2(I,K,Jp)
            Q    (I+EX,K,Jd) = qminus(I,K,1,Jp)
            Qliq (I+EX,K,Jd) = qminus(I,K,2,Jp)
            Qice (I+EX,K,Jd) = qminus(I,K,3,Jp)
!           print*,'tst q3',Q(I+ex,k,jd),i,k,jd,iam
         END DO
         call period( futmp(1,K,Jd) )
         call period( fvtmp(1,K,Jd) )
         call period( ST   (1,K,Jd) )
         call period( Q    (1,K,Jd) )
         call period( Qliq (1,K,Jd) )
         call period( Qice (1,K,Jd) )
      END DO
   END DO       
!
!      TRANSFORM (DU,DV) FROM P-GRID TO (U,V)-GRID
!     
#if (defined SPMD)
!JJR
      src = iam-1
      dest  = iam+1
      if ( mod(iam,npr_y) == 0 ) src = -1
      if ( mod(iam+1,npr_y) == 0 ) dest = -1
      call mp_send3d( dest, src, NX,  NL,NY,                     &
                      1, NX,beglev,endlev, beglatdynex, endlatdynex,       &
                      1, NX,beglev,endlev, beglatdyn, beglatdyn,  fvtmp )
      call mp_recv3d( src, NX,  NL,NY,                             &
                      1, NX,beglev,endlev, beglatdynex, endlatdynex,       &
                      1, NX,beglev,endlev, endlatdynex, endlatdynex,fvtmp )

! Latitude distribution for plat=16 on 3 processors.             !wjp 2007.05
!  procId:      iam=0     |      iam=1     |      iam=2          !wjp 2007.05
!  CAM:    01 02 03 04 05 | 06 07 08 09 10 | 11 12 13 14 15 16   !wjp 2007.05
!  IAP:    16 15 14 13 12 | 11 10 09 08 07 | 06 05 04 03 02 01   !wjp 2007.05
!                left     |     current    |      right          !wjp 2007.05
!JJr   call mpi_move_right(fvtmp(1,1,beglatdyn), fvtmp(1,1,endlatdyn+1), NX*NL)
#endif
!
   DO J = beglatdyn, endlatdyn
      if ( J>1 .and. J<plat ) then
         DO K = beglev,endlev 
            DO I = IB, IE
               SU(I,K,J) = ( futmp(I,K,J) + futmp(I-1,K,J) ) * HALF
               SV(I,K,J) = fvtmp(I,K,J)*DXPVN(J) + fvtmp(I,K,J+1)*DXPVS(J)
            END DO
            call period( SU(1,K,J) )
            call period( SV(1,K,J) )
         END DO

      else if ( J==1 ) then
         DO K = beglev,endlev 
            DO I = IB, IE
               SU(I,K,J) = ZERO
               SV(I,K,J) = fvtmp(I,K,J)*DXPVN(J) + fvtmp(I,K,J+1)*DXPVS(J)
            END DO
            call period( SU(1,K,J) )
            call period( SV(1,K,J) )
         END DO

      else  ! J = plat       
         DO K = beglev,endlev 
            DO I = 1, NX
               SU(I,K,J) = ZERO
               SV(I,K,J) = ZERO
            END DO
         END DO
      end if
   END DO
!
!=============================== zhh ============================================
!!   if (masterproc) then
!!      write (6,*) 'At sub. trans_grid_pd'
!!      write (6,*) 'fv(1,26,115)=', fv(1,26,115), 'fv(1,26,114)=', fv(1,26,114)
!!      write (6,*) 'gfv(1,26,115)=', gfv(1,26,115), 'gfv(1,26,114)=', gfv(1,26,114)
!!      write (6,*) 'fvtmp(4,15,26)=', fvtmp(4,15,26), 'fvtmp(4,16,26)=', fvtmp(4,16,26)
!!      write (6,*) 'sv(4,15,26)=', sv(4,15,26), 'DXPVN(15)=', DXPVN(15), 'DXPVS(15)=', DXPVS(15)
!2007.8.13
!!      write (6,*) 'ST(4,1,26)=', ST(4,1,26), ' gt2(1,26,128)=', gt2(1,26,128)
!      do k = 1, NL
!         write (6,*) 'ST(4,1,',k,')=', ST(4,1,k), ' gt2(1,',k,',128)=', gt2(1,k,128)
!      end do
!!   end if
!!   stop    !!
!============================ 2007.7.20 =========================================
   deallocate(futmp)
   deallocate(fvtmp)

   return
end subroutine trans_grid_pd
