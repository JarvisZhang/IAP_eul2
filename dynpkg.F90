#include <misc.h>
#include <params.h>

subroutine dynpkg(adv_state,t2,fu,fv,qminus,flx_net,t3,u3,v3,q3,ps,     &
                  omga,phis,etamid,cwava,detam,dtime,pcnst,n3,npt )
!----------------------------------------------------------------------- 
! Purpose: Driving routines for dynamics and horizontal diffusion.
! Author : ZhangHe
! Completed: 2007.4.28
! Update:  ZhangHe, 2007.12.17, 'q3(plon,plev,plat,npt)' ==> q3(plon,plev,ppcnst,plat,npt)
!          ZhangHe, 2008.4.20, add calling sub. mass_engy 
!          ZhangHe, 2008.4.27, add calling sub. trans_af_mass 
!          ZhangHe, 2010.8.10, add NSLT 
!          ZhangHe, 2011.05.04, do not call HDIFUS and mass_engy if adiabatic run
!-----------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid,     only: plon, plev, plat, beglat, endlat,beglev,endlev,beglatdyn,iam,endlatdyn
   use flexib,     only: DTadv, DTIMEQ, NENG, NFST
   use IAP_prog,   only: U, V, T, Q, PLY, PT, UT, VT, TT, QT, GHS, P, Pstar1, PS2, psa,  &
                         Qliq, QTliq
   use IAP_grid,only:IE1,EX
   use tendency,   only: NADDSS, SU, SV, ST
   use hdif,       only: HDIFUS
   use vapor,      only: AVNEGQ
   use constituents, only: ppcnst                    !zhh 2007.12.17
   use scanslt,      only: advection_state           !zhh 2008.4.21
   use Trans_coef,   only: trans_af_mass             !zhh 2008.4.27
! ======================= added by zhh 2008.9.10 =========================
   use prognostics,  only: n3m1
! ======================= added by zhh 2008.9.10 =========================

   implicit none
#include <comctl.h>   ! adiabatic
!------------------------------Arguments--------------------------------
   integer,  intent(in) :: pcnst    ! number of advected constituents (including water vapor)
   integer,  intent(in) :: n3       ! current time level
   integer,  intent(in) :: npt      ! number of time levels in the dycore
   real(r8), intent(in) :: fu (plon,beglev:endlev,beglat:endlat)   ! tendency of u-wind
   real(r8), intent(in) :: fv (plon,beglev:endlev,beglat:endlat)   ! tendency of v-wind
   real(r8), intent(in) :: t2 (plon,beglev:endlev,beglat:endlat)   ! tendency of temperature
   real(r8), intent(inout) :: qminus(plon,beglev:endlev,pcnst,beglat:endlat) ! constituents
   real(r8), intent(in) :: flx_net(plon,beglat:endlat)  ! net flux from physics
   real(r8), intent(in) :: dtime                        ! timestep of stepon    
!
   real(r8), intent(inout) :: t3(plon,beglev:endlev,beglat:endlat,npt)  ! temperature
   real(r8), intent(inout) :: u3(plon,beglev:endlev,beglat:endlat,npt)  ! u-wind component
   real(r8), intent(inout) :: v3(plon,beglev:endlev,beglat:endlat,npt)  ! v-wind component
   real(r8), intent(inout) :: q3(plon,beglev:endlev,ppcnst,beglat:endlat,npt)  ! specific humidity
   real(r8), intent(inout) :: omga(plon,beglev:endlev,beglat:endlat)    ! p-surface vertical velocity
   real(r8), intent(inout) :: ps(plon,beglat:endlat,npt)       ! surface pressure
   real(r8), intent(inout) :: phis(plon,beglat:endlat)         ! surface geopotential
! ======================= added by zhh 2008.1.2 ====================
   type(advection_state), intent(inout) :: adv_state        ! Advection state data
   real(r8), intent(in) :: etamid(plev)     ! vertical coords at midpoints 
   real(r8), intent(inout) :: cwava(plat)   ! weight applied to global integrals
   real(r8), intent(inout) :: detam(plev)   ! intervals between vert full levs.
! ======================= added by zhh 2008.1.2 ====================

!---------------------------Local workspace-----------------------------
   integer :: NSEQ     ! dynamical computation times of every dtime
   integer :: NSLT     ! slt computation times of every DTadv
   integer :: j,k,i
! ======================= added by zhh 2008.9.10 =========================
   real(r8) :: dqphy(plon,beglev:endlev,beglat:endlat)     ! q tendency due to physics 
! ======================= added by zhh 2008.9.10 =========================
!------------------------------------------------------------------------------------------------
!
! ======================= added by zhh 2008.9.10 =========================
! compute q tendency due to physics
   do j=beglat, endlat
      do k=beglev,endlev
         do i=1,plon
!            dqphy(i,k,j) = (q3(i,k,1,j,n3) - q3(i,k,1,j,n3m1))/dtime
            dqphy(i,k,j) = (qminus(i,k,1,j) - q3(i,k,1,j,n3m1))/dtime
         end do
      end do
   end do
! ======================= added by zhh 2008.9.10 =========================

   NSEQ = dtime / DTadv + 0.01D0
   NSLT = DTadv / DTIMEQ + 0.01D0
! ============================================================================
   call trans_grid_pd( fu, fv, t2, qminus, beglat, endlat, pcnst ) !end jjr mpp 520
!!   print*, '--------- at dynpkg: after trans_grid_pd -------------'
!=============================================================================
!

! print*,'test jjr3',fu(50,1,92),fv(50,1,92),t2(50,1,92),qminus(50,4,1,92)
   IF ( NADDSS == 0 ) THEN
! ============ do first half-step horizontal diffusion ===================
      if(.not. adiabatic) CALL HDIFUS( dtime, PLY, U, V, T )  !end jjr mpp
! ========================================================================
   ENDIF
!
!!         print*,'test jjr in dynpkg1',iam

!!   print*, '--------- at dynpkg: before dyfram -------------'
   call t_startf('DYFRAM')
! ============== perform the dynamic integration cycle ==================
   call  DYFRAM( NENG,NFST,NSEQ,NSLT, adv_state, detam, etamid, cwava,   &
                 t3, u3, v3, q3, ps, omga, phis, n3, npt )    
! =======================================================================
!!   print*, '--------- at dynpkg: after dyfram -------------'

   call t_stopf('DYFRAM')
!!       print*,'test jjr in dynpkg2',iam


!
! =======================================================================
!  do other half-step horizontal diffusion ( NADDSS == 0 )
!  do first half-step horizontal diffusion ( NADDSS /= 0 )
   if(.not. adiabatic) CALL HDIFUS( dtime, PLY, U, V, T)
! =======================================================================
!

   IF ( NADDSS /= 0 ) THEN  
! ============ do other half-step horizontal diffusion ==================
      if(.not. adiabatic) CALL HDIFUS( dtime, PLY, U, V, T )
! =======================================================================
   ENDIF




!=============================== 2007.12.20 ===============================
   call trans_grid_dp( 1, t3, u3, v3, q3, omga, ps, phis, n3, npt )
!!   print*,'test jjr in dynpkg3',iam


!========================================================================
! mass and energy correction
   if ((.not. ideal_phys) .and. (.not. adiabatic)) call mass_engy (dtime, cwava, etamid, flx_net, fu, fv, t2, dqphy)  !zhh 2008.9.11
!
!!   print*,'test jjr in dynpkg4',iam


! Convert ps, t3, q3 to P, T, Q
   if ((.not. ideal_phys) .and. (.not. adiabatic)) call trans_af_mass( ps, t3, q3, ppcnst, n3, npt )  !zhh 2008.4.27
!

! ============================ for test =================================
   do j = beglatdyn, endlatdyn
      if ( j == 65.and.endlev==26 ) then
         write(6,*) '--------------- At the end of sub. dynpkg --------------------'
         write(6,*) 'UT(4,26,65) =', UT(4,26,65), ' T(4,26,65) =', T(4,26,65)
         write(6,*) 'SU(4,26,65) =', SU(4,26,65), ' SV(4,26,65) =', SV(4,26,65)
!wjp 2011.04.16         write(6,*) 'Q(4,26,65) =', Q(4,26,65), ' ps(4,64,n3) =', ps(4,64,n3)
!!         write(6,*) 'v3(4,26,64,n3) =', v3(4,26,64,n3)
!wjp 2011.04.16         write(6,*) 'q3(4,26,2,64,n3) =', q3(4,26,2,64,n3)
        write(6,*) '--------------------------------------------------------------'
      end if
   end do
! =============================== zhh ===================================
!
   return

end subroutine dynpkg

