#include <misc.h>
#include <params.h>

subroutine init_trans_IAP
!----------------------------------------------------------------------------------
! Purpose: Transform P, U, V, T, Q to PT, UT, VT, TT, QT at first time step.
! Author : ZhangHe
! Completed: 2007.5.8
! Update : ZhangHe, 2007.12.21
!----------------------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use IAP_grid,   only: NX, NY, NL
   use pmgrid,     only: beglatdyn, endlatdyn, masterproc,beglev,endlev
   use stdatm,     only: PSB, P00, TB
   use Trans_coef, only: PTU, PTV, PTT, TRANSC 
   use IAP_prog,   only: U, V, T, Q, Qliq, Qice, P, PS2, Psa, Pstar1,     &
                         PT, UT, VT, TT, QT, QTliq, QTice  
   use Dyn_const,  only: PMTOP   

   implicit none
!------------------------------Local workspace--------------------------------
   real(r8) :: sqP0
   integer  :: I, J, K
!-----------------------------------------------------------------------------
!
   sqP0 = sqrt(P00)
! 
   DO J = beglatdyn, endlatdyn
      DO I = 1, NX
         Psa(I,J) = PS2(I,J) - PSB(I,J)
         P  (I,J) = PS2(I,J) - PMTOP
         Pstar1(I,J) = P(I,J)
         PT(I,J)     = sqrt( Pstar1(I,J) ) / sqP0
      END DO
!--------- zhh debug 2013-02-21 ----------
      if (j==10) then
         print*, 'ps2(10,10) =', ps2(10,10), ' PSB(10,10) =', PSB(10,10)
         print*, 'PMTOP =', PMTOP, ' sqP0 =', sqP0
         print*, 'PT(10,10) =', PT(10,10)
      end if
!--------- zhh debug 2013-02-21 ----------
   END DO       
!
!------------------- CALCULATE TRANSFORMATION COEFFECIENTS -------------------------
   CALL TRANSC   
!-----------------------------------------------------------------------------------
!
   DO J = beglatdyn, endlatdyn
      DO K = beglev,endlev 
         DO I = 1, NX
            UT(I,K,J) = U(I,K,J) * PTU(I,J)
            VT(I,K,J) = V(I,K,J) * PTV(I,J)
            TT(I,K,J) = ( T(I,K,J)-TB(I,K,J) ) * PTT(I,J)
            QT(I,K,J) = Q(I,K,J) * P(I,J)
            QTliq(I,K,J) = Qliq(I,K,J) * P(I,J)
            QTice(I,K,J) = Qice(I,K,J) * P(I,J)
         END DO
      END DO
   END DO
!=============================== zhh ============================================
!   if (masterproc) then
!wjp 2008.05.05      write (6,*) 'At sub. init_trans_IAP'
!wjp 2008.05.05      write (6,*) 'v(4,15,26)=', v(4,15,26), 'VT(4,15,26)=', VT(4,15,26)
!!      write (6,*) 'stop sub. init_trans_IAP for test -zhh'
!!      stop
!   end if
!============================ 2007.7.19 =========================================
! ============================= for test ====================================
!wjp 2008.05.05   write(*,*) '===================================================================='
!wjp 2008.05.05   write(*,*) 'U(105,81,1)  =', U(105,81,1), ' V(105,81,1)  =', V(105,81,1)
!wjp 2008.05.05   write(*,*) 'U(105,81,26) =', U(105,81,26), ' V(105,81,26) =', V(105,81,26)
!wjp 2008.05.05   write(*,*) 'T(5,88,1) =', T(5,88,1), ' T(5,88,26) =', T(5,88,26)
!wjp 2008.05.05   write(*,*) 'PLY(105,81,1)  =', PLY(105,81,1), ' PLY(105,2,1)  =', PLY(105,2,1)
!wjp 2008.05.05   write(*,*) 'PLY(105,81,27) =', PLY(105,81,27), ' PLY(105,2,27) =', PLY(105,2,27)
!wjp 2008.05.05   write(*,*) 'Pstar1(10,2)=', Pstar1(10,2)
!wjp 2008.05.05   write(*,*) 'P(10,2)     =', P(10,2)
!wjp 2008.05.05   write(*,*) 'Pstar1(9,90)=', Pstar1(9,90)
!wjp 2008.05.05   write(*,*) 'P(9,90)     =', P(9,90)
!wjp 2008.05.05   write(*,*) '===================================================================='
!   stop 'init_trans_IAP-test'
! ============================ zhh 2007.11.3 =============================

   return
end subroutine init_trans_IAP
