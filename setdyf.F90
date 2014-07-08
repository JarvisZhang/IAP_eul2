#include <misc.h>
#include <params.h>

!=================================================================================
SUBROUTINE SETDYF
!---------------------------------------------------------------------------------
! Purpose: SET CONSTANTS USED IN DYNAMIC FRAME INTEGRATION
! Original version: ENTRY SETDYF at DYFRAM.f (IAP 9L)
! Reconstructed & revised : ZhangHe
! Completed : 2005.10.13
! Update    : 2006.10.10, Zhanghe, revise the calling of sub. CONPDA     
!           : 2006.12.12, Zhanghe, revise the calling of sub. STFRAM     
!           : 2007.05.03, Zhanghe,
!           : 2007.08.25, Zhanghe, add the calling of sub. setfle
!           : 2007.11.11, Zhanghe, move calling of STFRAM to sub. stepon 
!                                  move calling of setfle to sub. STFRAM
!---------------------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use flexib,     only: DTDY, PHALFC, EPS, IDSOSS, setfle   ! DFS0, DTHDFS 
   use mathconst,  only: ONE
   use stdatm,     only: PS0, SETMSA     
   use Dyn_const,  only: PMTOP, DLAT, DLON, STFRAM, CONPDA         
   use Engy_const, only: STENGC
   use vapor,      only: IQHADV 
   use smoother,   only: STSMOT
   use sm9h,       only: STSM9C
   use Filt,       only: CONFIL
   use hdif,       only: STDFSC

   IMPLICIT NONE
!----------------------------------Local workspace----------------------------------------------
   real(r8) :: FIRST                ! initial parameter to control WATER VAPOR prediction
!------------------------------------------------------------------------------------------------

!  FOR THE MODEL STANDARD ATMOSPHERE
!--------------- SET CONST OF THE MODEL STANDARD ATMOSPHERE ----------------------
   CALL SETMSA( EPS )
  
!---------------------------------------------------------------------------------
!------------------ FOR THE ENERGY BUGET CALCULATION ----------------------------
   CALL STENGC    

!---------------------------------------------------------------------------------
!
!--------------- FOR THE ADVECTION INTEGRATION OF Q-H2O --------------------------
   FIRST = -ONE
   IF ( IQHADV.EQ.+2 ) THEN
!-----------------------------------------------------------------------------------------
      CALL CONPDA( FIRST )  ! at module Dyn_const 
!-----------------------------------------------------------------------------------------
   ENDIF

!----------------------------- FOR THE SMOOTHERS/FILTERS ---------------------------------
   CALL STSMOT   

!-----------------------------------------------------------------------------------------
!--------------------- SET CONSTANTS USED IN SUB.SM9HAS ----------------------------------
   CALL STSM9C(1)
!-----------------------------------------------------------------------------------------
!--------------------- SET CONSTANTS USED IN SUB.FILTER ----------------------------------
   CALL CONFIL( IDSOSS,DLAT,DLON )  

!-----------------------------------------------------------------------------------------
!------------- Set constants for computation of the horizontal diffusion ----------------- 
   CALL STDFSC       
!-----------------------------------------------------------------------------------------
   RETURN
END

