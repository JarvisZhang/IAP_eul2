#include <misc.h>
#include <params.h>
module commap

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid, only: plev, plat, plon
   use pspect, only: pmmax, pnmax

!   private

!!   public bps, sq, rsq, slat, w, cs, href, ecref, clat, clon, latdeg, bm1
   public bps, slat, w, cs, href, ecref, clat, clon, latdeg, bm1
   public tau, londeg, t0, xm

   real(r8) :: bps(plev)          ! coefficient for ln(ps) term in divergence eqn
!================== delete by zhh (2007.5.25) ===========================
!!   real(r8) :: sq(pnmax)          ! n(n+1)/a^2  (del^2 response function)
!!   real(r8) :: rsq(pnmax)         ! a^2/(n(n+1))
!============================== delete by zhh ===========================
#if ( defined SCAM )
   real(r8) :: slat(plat)         ! |sine latitude| (hemisphere)
#else
   real(r8) :: slat(plat/2)       ! |sine latitude| (hemisphere)
#endif
   real(r8) :: w(plat)            ! gaussian weights (hemisphere)
#if ( defined SCAM )
   real(r8) :: cs(plat)           ! cosine squared latitude (hemisphere)
#else
   real(r8) :: cs(plat/2)         ! cosine squared latitude (hemisphere)
#endif
   real(r8) :: href(plev,plev)    ! reference hydrostatic equation matrix
   real(r8) :: ecref(plev,plev)   ! reference energy conversion matrix
   real(r8) :: clat(plat)         ! model latitudes (radians)
   real(r8) :: clon(plon,plat)    ! model longitudes (radians)
   real(r8) :: latdeg(plat)       ! model latitudes (degrees)
   real(r8) :: bm1(plev,plev,pnmax)     ! transpose of right eigenvectors of semi-implicit matrix
   real(r8) :: tau(plev,plev )    ! matrix for reference d term in thermodynamic eqn
   real(r8) :: londeg(plon,plat)  ! model longitudes (degrees)
   real(r8) :: t0(plev)           ! Reference temperature for t-prime computations
   real(r8) :: xm(pmmax)          ! m (longitudinal wave number)
   real(r8) w_staggered(plat-1)      ! integration weights for the staggered wind arrays
 ! JJR
   real(r8) clat_staggered (plat-1)  ! model latitudes on staggered grid (radians)

end module commap
