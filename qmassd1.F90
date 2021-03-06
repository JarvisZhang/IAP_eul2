#include <misc.h>
#include <params.h>

subroutine qmassd1(cwava   ,etamid  ,w       ,q1      ,q2      , &
                  pdel    ,hwn     ,nlon    )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute comtribution of current latitude to global integral of
! q2*|q2 - q1|*eta
! This is a measure of the difference between the fields before and
! after the SLT "forecast" weighted by the approximate mass of the tracer.
! It is used in the "fixer" which enforces conservation in constituent
! fields transport via SLT. 
! 
! Method: 
! Reference Rasch and Williamson, 1991, Rasch, Boville and Brasseur, 1995
! 
! Author: J. Olson
! 
!-----------------------------------------------------------------------
!
! $Id: qmassd.F90,v 1.1.2.3 2004/09/23 17:20:10 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid,       only: plev, plon,beglev,endlev
  use constituents, only: pcnst

  implicit none

!------------------------------Arguments--------------------------------
  integer , intent(in)  :: nlon                 ! longitude dimension  
  real(r8), intent(in)  :: cwava                ! normalization factor
  real(r8), intent(in)  :: etamid(plev)         ! vertical coords at midpoints 
  real(r8), intent(in)  :: w                    ! gaussian weight this latitude
  real(r8), intent(in)  :: q1(plon,beglev:endlev)        ! constituents (pre -SLT)
  real(r8), intent(in)  :: q2(plon,beglev:endlev)        ! constituents (post-SLT)
  real(r8), intent(in)  :: pdel(plon,beglev:endlev)      ! pressure diff between interfaces
  real(r8), intent(inout) :: hwn(pcnst)         ! accumulator for global integrals
!
! cwava l/(g*plon)
! w     Gaussian weight.
! q1    Untransported q-field.
! q2    Transported   q-field.
! pdel  array of pressure differences between layer interfaces (used for mass weighting)
! hwn   Mass averaged constituent in units of kg/m**2.
!-----------------------------------------------------------------------

!---------------------------Local variables-----------------------------
  integer i,k               ! longitude and level indices
  real(r8) hwava            ! accumulator
!-----------------------------------------------------------------------
!
  hwava = 0.0
  do k=beglev,endlev
     do i=1,nlon
        hwava = hwava + (q2(i,k)* etamid(k)*(abs(q1(i,k) - q2(i,k))))*pdel(i,k)
     end do
  end do
!
! The 0.5 factor arises because gaussian weights sum to 2
!
  hwn(1) = hwn(1) + cwava*w*hwava*0.5

  return
end subroutine qmassd1

