#include <misc.h>
#include <params.h>

subroutine init_trans_pd
!----------------------------------------------------------------------------------
! Purpose: Transform initial data from physics grid to dynamical grid.
! Author : ZhangHe
! Completed: 2007.5.7
! Update : ZhangHe, 2007.10.11, transform direction of V-wind from north(v3) to 
!                               south(V)
!          ZhangHe, 2007.12.21, add Qliq & Qice
!          ZhangHe, 2008.4.8, add calculation of P
!          WuJianping, 2008.05, parallel version
!----------------------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use IAP_grid,    only: NX, NL, IB, IE, EX, period,NY
   use pmgrid
!,      only: plon, plat, plev, beglat, endlat, beglatdyn, endlatdyn, masterproc, iam
   use IAP_prog,    only: P, U, V, T, Q, Qliq, Qice, PS2, GHS
   use prognostics, only: u3, v3, t3, q3, ps, phis, n3m2
   use Dyn_const,   only: DXPVN, DXPVS, PMTOP, STFRAM   
   use mathconst,   only: ZERO, HALF
   use constituents, only: pnats, pcnst, ppcnst
   use flexib,      only: PHALFC
#if (defined SPMD || defined COUP_CSM)
   use mpishorthand, only: mpir8, mpicom
   use spmd_dyn,     only: npes
   use mod_comm, only: mp_send3d, mp_recv3d

#endif

   implicit none
!-------------------------------------Local workspace--------------------------------------------
   real(r8), allocatable :: utmp(:,:,:)    ! u wind 
   real(r8), allocatable :: vtmp(:,:,:)    ! v wind 
   integer  :: I, J, Jp, Jd, K, m,src,dest
!-----------------------------------------------------------------------------
!
   allocate(utmp(NX,beglev:endlev,beglatdynex:endlatdynex))
   allocate(vtmp(NX,beglev:endlev,beglatdynex:endlatdynex))
!
! for 2D variables
   DO Jp = beglat, endlat
      Jd = plat + 1 - Jp
      DO I = 1, plon
         PS2(I+EX,Jd) = ps(I,Jp,n3m2) / 100.0E0   !zhh 2007.7.19
         GHS(I+EX,Jd) = phis(I,Jp)
      END DO
!--------- zhh debug 2013-02-21 ----------
      if (jp==10) then
!         print*, 'n3m2 =', n3m2, '  plat =', plat
!         print*, 'jd =', jd, '  EX =', EX
         print*, 'ps(10,10,1) =', ps(10,10,1)
         print*, 'ps2(13,352) =', ps2(13,352)
      end if
      if (jp==352) then
         print*, 'ps(7,352,1) =', ps(7,352,1)
         print*, 'ps2(10,10) =', ps2(10,10)
      end if
!--------- zhh debug 2013-02-21 ----------
      call period( PS2(1,Jd) )
      call period( GHS(1,Jd) )
   END DO  
   DO J = beglatdyn, endlatdyn    
      DO I = 1, NX
         P(I,J) = PS2(I,J) - PMTOP         !zhh 2008.4.8
      END DO
   END DO       
!
!   for 3D variables
   DO Jp = beglat, endlat
      Jd = plat + 1 - Jp
      DO K = beglev,endlev 
         DO I = 1, plon
            utmp(I+EX,K,Jd) = u3(I,K,Jp,n3m2)
            vtmp(I+EX,K,Jd) = -v3(I,K,Jp,n3m2)
            T   (I+EX,K,Jd) = t3(I,K,Jp,n3m2)
            Q   (I+EX,K,Jd) = q3(I,K,1,Jp,n3m2)
!=============================== zhh =====================================
            Qliq(I+EX,K,Jd) = q3(I,K,2,Jp,n3m2)
            Qice(I+EX,K,Jd) = q3(I,K,3,Jp,n3m2)
!        print*,'tst q3',Q(I+ex,k,jd),i,k,jd,iam
!============================ 2007.12.21 =================================
         END DO
         call period( utmp(1,K,Jd) )
         call period( vtmp(1,K,Jd) )
         call period( T   (1,K,Jd) )
         call period( Q   (1,K,Jd) )
!=============================== zhh =====================================
         call period( Qliq(1,K,Jd) )
         call period( Qice(1,K,Jd) )
!============================ 2007.12.21 =================================
      END DO
   END DO       
!
!      TRANSFORM (U,V) FROM P-GRID TO (U,V)-GRID
!     
#if (defined SPMD)
! Latitude distribution for plat=16 on 3 processors.             !wjp 2007.05
!  procId:      iam=0     |      iam=1     |      iam=2          !wjp 2007.05
!  CAM:    01 02 03 04 05 | 06 07 08 09 10 | 11 12 13 14 15 16   !wjp 2007.05
!  IAP:    16 15 14 13 12 | 11 10 09 08 07 | 06 05 04 03 02 01   !wjp 2007.05
!                right    |     current    |       left          !wjp 2007.05
!   call mpi_move_right(vtmp(1,1,beglatdyn), vtmp(1,1,endlatdyn+1), NX*NL)
      src = iam-1
      dest  = iam+1
      if ( mod(iam,npr_y) == 0 ) src = -1
      if ( mod(iam+1,npr_y) == 0 ) dest = -1
      call mp_send3d( dest, src, NX, NL, NY,                      &
                      1, NX,beglev,endlev, beglatdynex, endlatdynex,       &
                      1, NX,beglev,endlev, beglatdyn, beglatdyn,  vtmp )
      call mp_recv3d( src, NX, NL, NY,                             &
                      1, NX, beglev,endlev,beglatdynex, endlatdynex,       &
                      1, NX, beglev,endlev,endlatdynex, endlatdynex, vtmp )

#endif
!
   DO J = beglatdyn, endlatdyn
      if ( J>1 .and. J<plat ) then
         DO K = beglev, endlev
            DO I = IB, IE
               U(I,K,J) = ( utmp(I,K,J) + utmp(I-1,K,J) ) * HALF
               V(I,K,J) = vtmp(I,K,J)*DXPVN(J) + vtmp(I,K,J+1)*DXPVS(J)
            END DO
            call period( U(1,K,J) )
            call period( V(1,K,J) )
         END DO

      else if ( J==1 ) then
         DO K = beglev,endlev 
            DO I = IB, IE
               U(I,K,J) = ZERO
               V(I,K,J) = vtmp(I,K,J)*DXPVN(J) + vtmp(I,K,J+1)*DXPVS(J)
            END DO
            call period( U(1,K,J) )
            call period( V(1,K,J) )
         END DO

      else  ! J = plat       
         DO K = beglev,endlev
            DO I = 1, NX
               U(I,K,J) = ZERO
               V(I,K,J) = ZERO
            END DO
         END DO
      end if
   END DO
!
   deallocate(utmp)
   deallocate(vtmp)

   return
end subroutine init_trans_pd
