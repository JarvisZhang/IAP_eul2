#include <misc.h>
#include <params.h>

subroutine herxin(pf      ,pkcnst  ,fb      ,fxl     ,fxr     , &
                  x       ,xdp     ,idp     ,jdp     ,fint    , &
                  nlon    ,nlonex  )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! 
! Method: 
! For each departure point in the latitude slice being forecast,
! interpolate (using equally spaced Hermite cubic formulas) to its
! x value at each latitude required for later interpolation in the y
! direction.
! 
! Author: 
! Original version:  J. Olson
! Standardized:      J. Rosinski, June 1992
! Reviewed:          D. Williamson, P. Rasch, August 1992
! Reviewed:          D. Williamson, P. Rasch, March 1996
!
!-----------------------------------------------------------------------
!
! $Id: herxin.F90,v 1.1.44.5 2005/02/06 19:40:49 rosinski Exp $
! $Author: rosinski $
!
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid,       only: plev, plon,beglev,endlev
   use scanslt,      only: plond, beglatex, endlatex, platd, nxpt
   use rgrid,        only: fullgrid
   use abortutils, only: endrun
!-----------------------------------------------------------------------
   implicit none
!------------------------------Parameters-------------------------------
#include <parslt.h>
!------------------------------Arguments--------------------------------
!
! Input arguments
!
   integer, intent(in) :: pf                   ! dimension (number of fields)
   integer, intent(in) :: pkcnst               ! dimension,=p3d
!
   real(r8), intent(in) :: fb (plond,beglev:endlev,pkcnst,beglatex:endlatex) ! field
   real(r8), intent(in) :: fxl(plond,beglev:endlev,pf,beglatex:endlatex)     ! left  x derivative
   real(r8), intent(in) :: fxr(plond,beglev:endlev,pf,beglatex:endlatex)     ! right x derivative
   real(r8), intent(in) :: x(plond,platd)      ! longitudinal grid coordinates
   real(r8), intent(in) :: xdp(plon,beglev:endlev)      ! departure point coordinates
!
   integer, intent(in) :: idp(plon,beglev:endlev,4)     ! longitude index of dep pt.
   integer, intent(in) :: jdp(plon,beglev:endlev)       ! latitude  index of dep pt.
   integer, intent(in) :: nlon
   integer, intent(in) :: nlonex(platd)
!
! Output arguments
!
   real(r8), intent(out) :: fint(plon,beglev:endlev,ppdy,pf) ! x-interpolants
!
!-----------------------------------------------------------------------
!
!  pf      Number of fields being interpolated.
!  pkcnst  Dimensioning construct for 3-D arrays.
!  fb      extended array of data to be interpolated.
!  fxl     x derivatives at the left edge of each interval containing 
!          the departure point
!  fxr     x derivatives at the right edge of each interval containing 
!          the departure point
!  x       Equally spaced x grid values in extended arrays.
!  xdp     xdp(i,k) is the x-coordinate (extended grid) of the
!          departure point that corresponds to global grid point (i,k)
!          in the latitude slice being forecasted.
!  idp     idp(i,k) is the index of the x-interval (extended grid) that
!          contains the departure point corresponding to global grid
!          point (i,k) in the latitude slice being forecasted.
!          Note that
!                x(idp(i,k)) .le. xdp(i,k) .lt. x(idp(i,k)+1) .
!  jdp     jdp(i,k) is the index of the y-interval (extended grid) that
!          contains the departure point corresponding to global grid
!          point (i,k) in the latitude slice being forecasted.
!          Suppose yb contains the y-coordinates of the extended array
!          and ydp(i,k) is the y-coordinate of the departure point
!          corresponding to grid point (i,k).  Then,
!                yb(jdp(i,k)) .le. ydp(i,k) .lt. yb(jdp(i,k)+1) .
!  fint    (fint(i,k,j,n),j=1,ppdy) contains the x interpolants at each
!          latitude needed for the y derivative estimates at the
!          endpoints of the interval that contains the departure point
!          for grid point (i,k).  The last index of fint allows for
!          interpolation of multiple fields.
!
!---------------------------Local workspace-----------------------------
!
   integer i,j,k,m           ! indices
!
   real(r8) dx (platd)           ! x-increment
   real(r8) rdx(platd)           ! 1./dx
   real(r8) xl (plon,beglev:endlev)       ! |
   real(r8) xr (plon,beglev:endlev)       ! |
   real(r8) hl (plon,beglev:endlev)       ! | --interpolation coeffs
   real(r8) hr (plon,beglev:endlev)       ! |
   real(r8) dhl(plon,beglev:endlev)       ! |
   real(r8) dhr(plon,beglev:endlev)       ! |

   integer n

!
!-----------------------------------------------------------------------
!
   if(ppdy .ne. 4) then
      call endrun ('HERXIN:Fatal error: ppdy must be set to 4')
   end if
!
   if (fullgrid) then
      dx (1) = x(nxpt+2,1) - x(nxpt+1,1)
      rdx(1) = 1./dx(1)
      do k=beglev,endlev
         do i=1,nlon
            xl (i,k) = ( x(idp(i,k,1)+1,1) - xdp(i,k) )*rdx(1)
         enddo
      enddo
      do k=beglev,endlev
         do i=1,nlon
            xr (i,k) = 1. - xl(i,k)
            hl (i,k) = ( 3.0 - 2.0*xl(i,k) )*xl(i,k)**2
            hr (i,k) = ( 3.0 - 2.0*xr(i,k) )*xr(i,k)**2
            dhl(i,k) = -dx(1)*( xl(i,k) - 1. )*xl(i,k)**2
            dhr(i,k) =  dx(1)*( xr(i,k) - 1. )*xr(i,k)**2
         end do
      end do
!
! x interpolation at each latitude needed for y interpolation.
! Once for each field.
! 
      do m = 1,pf
         do n=1,4
            do k = beglev,endlev
               do i = 1,nlon
                  fint(i,k,n,m) = &
                       fb (idp(i,k,1)  ,k,m,jdp(i,k)+(n-2))*hl (i,k) + &
                       fb (idp(i,k,1)+1,k,m,jdp(i,k)+(n-2))*hr (i,k) + &
                       fxl(idp(i,k,1)  ,k,m,jdp(i,k)+(n-2))*dhl(i,k) + &
                       fxr(idp(i,k,1)  ,k,m,jdp(i,k)+(n-2))*dhr(i,k)      
               enddo
            enddo
         enddo
      enddo
!
   else
!
      do j = 1,platd
         dx (j) = x(nxpt+2,j) - x(nxpt+1,j)
         rdx(j) = 1./dx(j)
      end do
!
      do k=beglev,endlev
         do i=1,nlon
            xl (i,k) = ( x(idp(i,k,1)+1,jdp(i,k)-1) - xdp(i,k) )*  &
               rdx(jdp(i,k)-1)
            xr (i,k) = 1. - xl(i,k)
            hl (i,k) = ( 3.0 - 2.0*xl(i,k) )*xl(i,k)**2
            hr (i,k) = ( 3.0 - 2.0*xr(i,k) )*xr(i,k)**2
            dhl(i,k) = -dx(jdp(i,k)-1)*( xl(i,k) - 1. )*xl(i,k)**2
            dhr(i,k) =  dx(jdp(i,k)-1)*( xr(i,k) - 1. )*xr(i,k)**2
         end do
      end do
!
! x interpolation at each latitude needed for y interpolation.
! Once for each field.
! 
      do m = 1,pf
         do k = beglev,endlev
            do i = 1,nlon
               fint(i,k,1,m) = &
                  fb (idp(i,k,1)  ,k,m,jdp(i,k)-1)*hl (i,k) + &
                  fb (idp(i,k,1)+1,k,m,jdp(i,k)-1)*hr (i,k) + &
                  fxl(idp(i,k,1)  ,k,m,jdp(i,k)-1)*dhl(i,k) + &
                  fxr(idp(i,k,1)  ,k,m,jdp(i,k)-1)*dhr(i,k)
            end do
         end do
      end do

      do k=beglev,endlev
         do i=1,nlon
            xl (i,k) = ( x(idp(i,k,2)+1,jdp(i,k)) - xdp(i,k) )* rdx(jdp(i,k))
            xr (i,k) = 1. - xl(i,k)
            hl (i,k) = ( 3.0 - 2.0*xl(i,k) )*xl(i,k)**2
            hr (i,k) = ( 3.0 - 2.0*xr(i,k) )*xr(i,k)**2
            dhl(i,k) = -dx(jdp(i,k))*( xl(i,k) - 1. )*xl(i,k)**2
            dhr(i,k) =  dx(jdp(i,k))*( xr(i,k) - 1. )*xr(i,k)**2
         end do
      end do
!
! x interpolation at each latitude needed for y interpolation.
! Once for each field.
! 
      do m = 1,pf
         do k = beglev,endlev
            do i = 1,nlon
               fint(i,k,2,m) = &
                  fb (idp(i,k,2)  ,k,m,jdp(i,k)  )*hl (i,k) + &
                  fb (idp(i,k,2)+1,k,m,jdp(i,k)  )*hr (i,k) + &
                  fxl(idp(i,k,2)  ,k,m,jdp(i,k)  )*dhl(i,k) + &
                  fxr(idp(i,k,2)  ,k,m,jdp(i,k)  )*dhr(i,k)
            end do
         end do
      end do

      do k=beglev,endlev
         do i=1,nlon
            xl (i,k) = ( x(idp(i,k,3)+1,jdp(i,k)+1) - xdp(i,k) )* rdx(jdp(i,k)+1)
            xr (i,k) = 1. - xl(i,k)
            hl (i,k) = ( 3.0 - 2.0*xl(i,k) )*xl(i,k)**2
            hr (i,k) = ( 3.0 - 2.0*xr(i,k) )*xr(i,k)**2
            dhl(i,k) = -dx(jdp(i,k)+1)*( xl(i,k) - 1. )*xl(i,k)**2
            dhr(i,k) =  dx(jdp(i,k)+1)*( xr(i,k) - 1. )*xr(i,k)**2
         end do
      end do
!
! x interpolation at each latitude needed for y interpolation.
! Once for each field.
! 
      do m = 1,pf
         do k = beglev,endlev
            do i = 1,nlon
               fint(i,k,3,m) = &
                  fb (idp(i,k,3)  ,k,m,jdp(i,k)+1)*hl (i,k) + &
                  fb (idp(i,k,3)+1,k,m,jdp(i,k)+1)*hr (i,k) + &
                  fxl(idp(i,k,3)  ,k,m,jdp(i,k)+1)*dhl(i,k) + &
                  fxr(idp(i,k,3)  ,k,m,jdp(i,k)+1)*dhr(i,k)
            end do
         end do
      end do
!
      do k=beglev,endlev
         do i=1,nlon
            xl (i,k) = ( x(idp(i,k,4)+1,jdp(i,k)+2) - xdp(i,k) )*rdx(jdp(i,k)+2)
            xr (i,k) = 1. - xl(i,k)
            hl (i,k) = ( 3.0 - 2.0*xl(i,k) )*xl(i,k)**2
            hr (i,k) = ( 3.0 - 2.0*xr(i,k) )*xr(i,k)**2
            dhl(i,k) = -dx(jdp(i,k)+2)*( xl(i,k) - 1. )*xl(i,k)**2
            dhr(i,k) =  dx(jdp(i,k)+2)*( xr(i,k) - 1. )*xr(i,k)**2
         end do
      end do
!
! x interpolation at each latitude needed for y interpolation.
! Once for each field.
! 
      do m = 1,pf
         do k = beglev,endlev
            do i = 1,nlon
               fint(i,k,4,m) = &
                  fb (idp(i,k,4)  ,k,m,jdp(i,k)+2)*hl (i,k) + &
                  fb (idp(i,k,4)+1,k,m,jdp(i,k)+2)*hr (i,k) + &
                  fxl(idp(i,k,4)  ,k,m,jdp(i,k)+2)*dhl(i,k) + &
                  fxr(idp(i,k,4)  ,k,m,jdp(i,k)+2)*dhr(i,k)
            end do
         end do
      end do
   end if
!
   return
end subroutine herxin
