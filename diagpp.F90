#include <misc.h>
#include <params.h>

SUBROUTINE DIAGPP
!------------------------------------------------------------------------------------------------
! Purpose: COMPUTE PRESSURE AT MODEL & INTERFACE SIGMA LAYERS
! Original version : DIAGPP.f (IAP 9L)
! Reconstructed & supplemented : ZhangHe, add the calculation of PIN
! Completed : 2005.9.1
!------------------------------------------------------------------------------------------------
   use IAP_grid,   only : NX, NY, NL, NZ
   use Dyn_const,  only : PMTOP, SIG, SIGL
   use IAP_prog,   only : P, PLY, PIN
   use pmgrid,     only : beglatdyn, endlatdyn,beglev,endlev
! ========================= test ==========================
   use pmgrid,     only: masterproc, iam
! ===================== zhh 2008.6.12 =====================

   implicit none
!----------------------------------Local workspace-----------------------------------------------
   integer  :: I,J,K             ! loop index
!------------------------------------------------------------------------------------------------
!
! ========================= test ==========================
!!   if (masterproc) print*, 'at beginning of sub. DIAGPP'
! ===================== zhh 2008.6.12 =====================
!---------------------------- compute model layer pressure ------------------------------
   DO J = beglatdyn, endlatdyn
      DO K = 1 ,NL
         DO I = 1 ,NX
            PLY(I,K,J) = P(I,J)*SIGL(K) + PMTOP
         END DO
      END DO
      DO I = 1 ,NX
         PLY(I,NZ,J)   = P(I,J)         + PMTOP
      END DO
   END DO
! ========================= test ==========================
!!   if (masterproc) print*, 'PLY(1,1,65) =', PLY(1,1,65)
!!   if (masterproc) print*, 'P(1,65) =', P(1,65)
! ===================== zhh 2008.6.12 =====================
!---------------------------- compute interface layer pressure ------------------------------
   do J = beglatdyn, endlatdyn
!      do K = 2 ,NL      
      do K = 1 ,NL      
         do I = 1 ,NX
            PIN(I,K,J) = P(I,J)*SIG(K)  + PMTOP
         end do
      end do
	  
      do I = 1 ,NX
        
        PIN(I,1,J)   = PMTOP
       PIN(I,NZ,J)   = PMTOP + P(I,J)
      end do
   end do
      
   RETURN
END
