#include <misc.h>
#include <params.h>

module prognostics

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Prognostic variables held in-core for convenient access.
! q3 is specific humidity (water vapor) and other constituents.
! pcnst is advected constituents, pnats is non-advected.
! 
! Author: G. Grant
! 
! Revised by: ZhangHe, 2007.4.28
!-----------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid!jjr
   use infnan,       only: inf
   use constituents, only: pcnst, pnats, cnst_need_pdeldry,ppcnst


   implicit none

   private

!!   public ps, u3, v3, t3, q3, qminus, vort, div, dpsl, dpsm, dps, omga, phis, hadv, pdeld
   public ps, u3, v3, t3, q3, qminus, omga, phis, hadv, pdeld
   public n3, n3m1, n3m2, ptimelevels
   public initialize_prognostics
   public shift_time_indices

   integer, parameter :: ptimelevels = 3  ! number of time levels in the dycore
   integer :: n3   = 3                    ! current time level
   integer :: n3m1 = 2                    ! previous time level
   integer :: n3m2 = 1                    ! time level before previous
!-----------------------------------------------------------------------
! for CAM:  three time level  
!           n3m2     n3m1      n3
!             |________|________|     f(n3) = f(n3m2) + df * 2dt
!            n-2      n-1       n
!
! for IAP:  two time level (time n3m2 is useless) 
!           n3m2     n3m1      n3
!             |________|________|     f(n3) = f(n3m1) + df * dt
!            n-2      n-1       n
!-----------------------------------------------------------------------
   real(r8), allocatable :: ps(:,:,:)       ! surface pressure
   real(r8), allocatable :: u3(:,:,:,:)     ! u-wind component
   real(r8), allocatable :: v3(:,:,:,:)     ! v-wind component
   real(r8), allocatable :: t3(:,:,:,:)     ! temperature
   real(r8), allocatable :: pdeld(:,:,:,:)  ! layer thickness dry (Pa)
   real(r8), allocatable :: q3(:,:,:,:,:)   ! specific humidity
   real(r8), allocatable :: qminus(:,:,:,:) ! constituents
   real(r8), allocatable :: hadv  (:,:,:,:) ! horizontal advection tendency

!!   real(r8), allocatable :: vort(:,:,:,:)   ! vorticity
!!   real(r8), allocatable :: div(:,:,:,:)    ! divergence

!!   real(r8), allocatable :: dpsl(:,:)       ! longitudinal pressure gradient
!!   real(r8), allocatable :: dpsm(:,:)       ! meridional pressure gradient
!!   real(r8), allocatable :: dps(:,:)        ! pressure gradient
   real(r8), allocatable :: phis(:,:)       ! surface geopotential
   real(r8), allocatable :: omga(:,:,:)     ! vertical velocity

!================================================================================================
CONTAINS
!================================================================================================

!================================================================================================
   subroutine initialize_prognostics
!-----------------------------------------------------------------------------------------------
! Purpose:  Allocate and initialize the prognostic arrays.
!jjr divide klevel
!-----------------------------------------------------------------------------------------------
      allocate (ps    (plon                 ,beglat:endlat    ,ptimelevels))
      allocate (u3    (plon,beglev:endlev            ,beglat:endlat,ptimelevels))
      allocate (v3    (plon,beglev:endlev            ,beglat:endlat,ptimelevels))
      allocate (t3    (plon,beglev:endlev            ,beglat:endlat,ptimelevels))
      allocate (q3    (plon,beglev:endlev,pcnst+pnats,beglat:endlat,ptimelevels))
      allocate (qminus(plon,beglev:endlev,pcnst      ,beglat:endlat  ))
      allocate (hadv  (plon,beglev:endlev,pcnst      ,beglat:endlat  ))

!!      allocate (vort  (plon,plev,beglat:endlat,ptimelevels))   
!!      allocate (div   (plon,plev,beglat:endlat,ptimelevels))    

!!      allocate (dpsl  (plon,beglat:endlat))        
!!      allocate (dpsm  (plon,beglat:endlat))        
!!      allocate (dps   (plon,beglat:endlat))         
      allocate (phis  (plon,beglat:endlat))        
      allocate (omga  (plon,beglev:endlev,beglat:endlat))    
!jjr      if (cnst_need_pdeldry) then
         allocate (pdeld (plon,beglev:endlev,beglat:endlat,ptimelevels))
!jjr      else
!jjr        allocate (pdeld (1,1,1:1,ptimelevels))
!jjr      endif

      ps(:,:,:)       = inf
      u3(:,:,:,:)     = inf
      v3(:,:,:,:)     = inf
      t3(:,:,:,:)     = inf
      pdeld(:,:,:,:)  = inf
      q3(:,:,:,:,:)   = inf
      qminus(:,:,:,:) = inf
      hadv  (:,:,:,:) = inf

!!      vort(:,:,:,:)   = inf
!!      div (:,:,:,:)   = inf

!!      dpsl  (:,:) = inf
!!      dpsm  (:,:) = inf
!!      dps   (:,:) = inf
      phis  (:,:) = inf
      omga  (:,:,:) = inf

      return
   end subroutine initialize_prognostics

!================================================================================================
   subroutine shift_time_indices
!-----------------------------------------------------------------------------------------------
! Purpose: 
! Shift the indices that keep track of which index stores
! the relative times (current time, previous, time before previous etc).
!-----------------------------------------------------------------------------------------------
      integer :: itmp

      itmp = n3m2

!-----------------------------------
      n3m2 = n3m1
      n3m1 = n3
!zhh 2007.12.17      n3m2 = n3     
!-----------------------------------
      n3   = itmp
      return        !zhh 2007.4.28
   end subroutine shift_time_indices

end module prognostics
