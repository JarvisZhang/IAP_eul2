#include <misc.h>
#include <params.h>

SUBROUTINE NLITTI
!------------------------------------------------------------------------------------------------
! Purpose: NONLINEAR ITERATIVE TIME INTEGRATION & splitting method
! Original version : NLITTI.f (IAP 21L)
! Reconstructed & add splitting method in : ZhangHe
! Completed : 2006.2.20
! Update    : 2006.12.10, Zhanghe, delete CALL pstartend2 when integration for advection term  
!             2007.04.23, ZhangHe, 1) delete dumb parameter 'SETUV' & 'NCTDCB'
!                                  2) 'ptuvtend1' ==> 'tend_lin'; 'ptuvtend2' ==> 'tend_adv' 
!             2007.05.07, ZhangHe, 1) 'pstartend' ==> 'tend_pstar' ;
!                                  2) 'nliter2'   ==> 'nliter_uvtp' ;
!                                  3) 'nliter3'   ==> 'nliter_uvt'
!             2007.12.22, ZhangHe, change calling QPDATA from lin timestep to adv timestep
!             2008.04.10, ZhangHe, add calculation of UVW0, update calling nliter_uvtp
!                                  move statement of Istar to module flexib
!             2008.04.23, ZhangHe, delete calling QPDATA & calculation of UVW0
!             2008.06.05, WuJianping, for parallel version
!------------------------------------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use IAP_grid,  only : NX, NY, NL,ex
   use IAP_prog,  only : Psa, UT, VT, TT, Pstar1, U, V, T, P, Q,PT
   use flexib,    only : DTDY, DTlin, DTadv, Ndt, IBCFFT, Istar,testnum
   use mathconst, only : HALF
   use pmgrid,    only : masterproc, beglatdyn,endlatdyn,iam,beglev,endlev,plev,plon,plat     !zhh 2007.7.20
! ============================ zhh 2007.12.14 =============================
   use Trans_coef, only: trans_antiIAP, trans_IAP
   use smoother,   only: SMOTHP, SMOOTH
! ============================ zhh 2007.12.14 =============================

   implicit none
!----------------------------------Local workspace-----------------------------------------------
   real(r8) :: PsaL(NX,beglatdyn:endlatdyn)
   real(r8) :: PstL(NX,beglatdyn:endlatdyn)
   real(r8) :: UL(NX,beglev:endlev,beglatdyn:endlatdyn)
   real(r8) :: VL(NX,beglev:endlev,beglatdyn:endlatdyn)
   real(r8) :: TL(NX,beglev:endlev,beglatdyn:endlatdyn)
   real(r8) :: time_begin, time_end
   integer  :: I, J, K, n, ITERNM,tmp1    ! loop index
!------------------------------------------------------------------------------------------------

![1]  integration for linear(adaption) term  
!    print*,'test ndt',ndt
   do n = 1, Ndt
      DO J = beglatdyn, endlatdyn
         DO I = 1, NX
            PsaL(I,J)    = Psa(I,J)
            PstL(I,J)    = Pstar1(I,J)
         END DO
      END DO
      DO J = beglatdyn, endlatdyn
         DO K = beglev,endlev 
            DO I = 1, NX
               UL(I,K,J) = UT(I,K,J)
               VL(I,K,J) = VT(I,K,J)
               TL(I,K,J) = TT(I,K,J)
            END DO
         END DO
      END DO
!
!-----------------------------------------------
      CALL tend_pstar( Istar )
!!      if (masterproc) print*, '----- 1st calling tend_pstar ----, n =', n
!      print*, '----- 1st calling tend_pstar ----, n =', n, iam
!-----------------------------------------------
      CALL tend_lin       ! 2007.4.23
!      if (masterproc) print*, '----- 1st calling tend_lin ----, n =', n
!-----------------------------------------------
      CALL nliter_uvtp(UL, VL, TL, PsaL, PstL, DTlin, Istar)
! ============================================

      DO ITERNM = 1,1
!------------------------------------------------------------------------------------------------
         CALL tend_pstar( Istar )
!         if (masterproc) print*, '----- 2nd calling tend_pstar ----, n =', n
! ============================================      

         CALL tend_lin
!         if (masterproc) print*, '----- 2nd calling tend_lin ----, n =', n
!       tmp1 = 0
!         do while(tmp1.eq.0)
!          call sleep(2)
!       enddo

! ======================================================================
         CALL nliter_uvtp(UL, VL, TL, PsaL, PstL, DTlin, Istar)

!  if(ndt.eq.1)      print*,'test60',q(i+ex,k,j),i,k,j,iam

! ======================================================================
!------------------------------------------------------------------------------------------------
         DO J = beglatdyn, endlatdyn
            DO I = 1, NX
               Psa(I,J)  = (Psa(I,J) + PsaL(I,J)) * HALF
               Pstar1(I,J)  = (Pstar1(I,J) + PstL(I,J)) * HALF     
            END DO
         END DO
         DO J = beglatdyn, endlatdyn
            DO K = beglev,endlev 
               DO I = 1, NX
                  UT(I,K,J) = (UT(I,K,J) + UL(I,K,J)) * HALF
                  VT(I,K,J) = (VT(I,K,J) + VL(I,K,J)) * HALF
                  TT(I,K,J) = (TT(I,K,J) + TL(I,K,J)) * HALF
               END DO
            END DO
         END DO

!------------------------------------------------------------------------------------------------
! ============================================
         CALL tend_pstar( Istar )
!         if (masterproc) print*, '----- 3rd calling tend_pstar ----, n =', n
! ============================================       
    
         CALL tend_lin
!         if (masterproc) print*, '----- 3rd calling tend_lin ----, n =', n

! ====================================================================
         CALL nliter_uvtp(UL, VL, TL, PsaL, PstL, DTlin, Istar)
! ====================================================================
!------------------------------------------------------------------------------------------------
      END DO
   end do    ! for n = 1, Nt

!********************************************************************************
![2]  integration for advection term
   do n = 1, 1
      DO J = beglatdyn, endlatdyn
         DO K = beglev,endlev 
            DO I = 1, NX
               UL(I,K,J) = UT(I,K,J)
               VL(I,K,J) = VT(I,K,J)
               TL(I,K,J) = TT(I,K,J)
            END DO
         END DO
      END DO
!
! ============================================

      CALL tend_adv
!      if (masterproc) print*, '----- 1st calling tend_adv ----, n =', n

      CALL nliter_uvt(UL, VL, TL, DTadv)
! ============================================

!

!  if(ndt.eq.1)      print*,'test60',q(i+ex,k,j),i,k,j,iam

      DO ITERNM = 1,1
!------------------------------------------------------------------------------------------------
! ============================================
         CALL tend_adv
!         if (masterproc) print*, '----- 2nd calling tend_adv ----, n =', n
         CALL nliter_uvt(UL, VL, TL, DTadv)
! ============================================
!------------------------------------------------------------------------------------------------
         DO J = beglatdyn, endlatdyn
            DO K = beglev, endlev
               DO I = 1, NX
                  UT(I,K,J) = (UT(I,K,J) + UL(I,K,J)) * HALF
                  VT(I,K,J) = (VT(I,K,J) + VL(I,K,J)) * HALF
                  TT(I,K,J) = (TT(I,K,J) + TL(I,K,J)) * HALF
               END DO
            END DO
         END DO
!------------------------------------------------------------------------------------------------
! ============================================
         CALL tend_adv
!         if (masterproc) print*, '----- 3rd calling tend_adv ----, n =', n
         CALL nliter_uvt(UL, VL, TL, DTadv)
! ============================================
!------------------------------------------------------------------------------------------------
      END DO
! ============================ zhh 2007.12.15 =============================
       call t_startf('other')
      call trans_antiIAP

      CALL SMOTHP  ! SMOOTHING P & PT BY 2-ORDER SHAPIRO SMOOTHER(LON & LAT)

      CALL SMOOTH

      call trans_IAP
         call t_stopf('other')

! ============================ zhh 2007.12.15 =============================
   end do
   RETURN
END
