#include <misc.h>
#include <params.h>
#if ( defined SCAM )
#include <max.h>
#endif

subroutine initcom

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Initialize Model commons, including COMCON, COMHYB, COMMAP, COMSPE,
! and COMTRCNM
! 
! Method: 
! 
! Author: 
! Original version:  CCM1
! Standardized:      L. Bath, Jun 1992
!                    L. Buja, Feb 1996
!
!-----------------------------------------------------------------------
!
! $Id: initcom.F90,v 1.16.2.11 2005/02/06 19:40:50 rosinski Exp $
! $Author: rosinski $
!
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid,       only: plat, plev, plon, plevp, dyngrid_set, &
                           masterproc, iam
   use pspect
   use comspe
   use rgrid,     only: nlon, beglatpair, wnummax, nmmax, fullgrid
   use scanslt,   only: nlonex, platd, j1
   use gauaw_mod, only: gauaw
   use commap
!!   use dynconst, only: rearth, ra, dynconsti
   use physconst, only: rair, ra              !! zhh 2007.5.9
   use time_manager, only: get_step_size
   use abortutils, only: endrun
#if ( defined SCAM )
   use scamMod, only :columnlat,columnlon
#endif
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
#include <comctl.h>
!-----------------------------------------------------------------------
#include <comfft.h>
!-----------------------------------------------------------------------
#include <comhyb.h>
!-----------------------------------------------------------------------

!
! Local workspace
!
   real(r8) zsi(plat)      ! sine of latitudes
   real(r8) zw(plat)       ! Gaussian weights
   real(r8) zra2           ! ra squared
   real(r8) zalp(2*pspt)   ! Legendre function array
   real(r8) zdalp(2*pspt)  ! Derivative array
   real(r8) zslat          ! sin of lat  and cosine of colatitude

   integer i           ! longitude index
   integer j           ! Latitude index
   integer k           ! Level index
   integer kk          ! Level index
   integer kkk         ! Level index
   integer m,lm,mr,lmr ! Indices for legendre array
   integer n           ! Index for legendre array
   integer nkk         ! Print control variables
   integer ik1         ! Print index temporary variable
   integer ik2         ! Print index temporary variable
   integer itmp        ! Dimension of polynomial arrays temporary.
   integer iter        ! Iteration index
   real(r8)    zdt         ! Time step for settau

   logical lprint      ! Debug print flag
   integer irow        ! Latitude pair index
   integer lat         ! Latitude index

   real(r8) xlat           ! Latitude (radians)
   real(r8) pi             ! Mathematical pi (3.14...)
   real(r8) dtime          ! timestep size [seconds]
   real(r8) sum,dp
!
!-----------------------------------------------------------------------
!!   call dynconsti
!
   lprint = masterproc .and. .FALSE.

   dtime = get_step_size()

!!   call hdinti  (rearth  ,dtime   )
!
! Initialize commap.  Set hybrid level dependent arrays
!=========================================
   call hycoef
!=========================================

#if ( !defined SCAM )
!
! NMAX dependent arrays
!
   if (pmmax.gt.plon/2) then
      call endrun ('INITCOM:mmax=ptrm+1 .gt. plon/2')
   end if
#endif

!================= delete by zhh (2007.5.25) ===========================
!!   zra2 = ra*ra
!!   do j=2,pnmax
!!      sq(j)  = j*(j-1)*zra2
!!      rsq(j) = 1./sq(j)
!!   end do
!!   sq(1)  = 0.
!!   rsq(1) = 0.
!============================== delete by zhh ===========================
!
! MMAX dependent arrays
!
   do j=1,pmmax
      xm(j) = j-1
   end do
#if ( !defined SCAM )
!
! Gaussian latitude dependent arrays
!
! ========================== zhh ==========================
   pi = 4.0*atan(1.0)
!
   call gauaw(zsi     ,zw      ,plat    )
   do irow=1,plat/2
      slat(irow) = zsi(irow)
!!      w(irow)              = zw(irow)
!!      w(plat - irow + 1)   = zw(irow)
      cs(irow)  = 1. - zsi(irow)*zsi(irow)
      xlat = asin(slat(irow))
      clat(irow) = -xlat
      clat(plat - irow + 1) = xlat
   end do
   do lat=1,plat
      latdeg(lat) = clat(lat)*45./atan(1._r8)
   end do
!
!!   print*, 'for CAM: latdeg(1) =', latdeg(1), 'latdeg(64) =', latdeg(64)
!!   print*, 'latdeg(128) =', latdeg(128)
!
   do irow=1,plat/2
      xlat = pi/2.0 - (irow-1)*pi / dble(plat-1)
      clat(irow) = -xlat
      clat(plat - irow + 1) = xlat
      slat(irow) = sin(xlat)
      cs(irow)  = 1. - slat(irow)*slat(irow)
   end do
! ======================= 2008.1.3 ========================
   do lat=1,plat
      latdeg(lat) = clat(lat)*45./atan(1._r8)
   end do
!
!!   print*, 'for IAP: latdeg(1) =', latdeg(1), 'latdeg(64) =', latdeg(64)
!!   print*, 'latdeg(128) =', latdeg(128)
#endif
   do lat = 1, plat-1
      clat_staggered(lat) = (clat(lat) + clat(lat+1)) / 2.0
   end do

! Weights are defined as cos(phi)*(delta-phi)
! For a sanity check, the sum of w across all lats should be 2, or 1 across
! half of the latitudes.

   do lat = 2, plat-1
      w(lat) = sin(clat_staggered(lat)) - sin(clat_staggered(lat-1))
   end do
   w(1) = sin(clat_staggered(1)) + 1.
   w(plat) = w(1)

   sum = 0.
   do lat=1,plat
      if (iam .eq. 0) write (6,*) 'initcom: lat, clat, w ', lat, clat(lat), w(lat)
      sum = sum + w(lat)
   end do

   if (abs(sum - 2._r8) > 1.e-8) then
      write(6,*) 'INITCOM 1: weights do not sum to 2. sum=',sum
      call endrun
   end if

   dp = pi / float(plat-1)
   do lat = 1, plat-1
      w_staggered(lat) = sin(clat(lat+1)) - sin(clat(lat))
   end do
   sum = 0.
   do lat=1,plat-1
      sum = sum + w_staggered(lat)
   end do

   if (abs(sum - 2._r8) > 1.e-8) then
      write(6,*) 'INITCOM 2: weights do not sum to 2. sum=',sum
      call endrun
   end if

!
! Integration matrices of hydrostatic equation(href) and conversion
! term(a).  href computed as in ccm0 but isothermal bottom ecref
! calculated to conserve energy
!
   do k=1,plev
      do kk=1,plev
         href(kk,k) = 0.
         ecref(kk,k) = 0.
      end do
   end do
!
! Mean atmosphere energy conversion term is consistent with continiuty
! Eq.  In ecref, 1st index = column; 2nd index = row of matrix.
! Mean atmosphere energy conversion term is energy conserving
!
   do k=1,plev
      ecref(k,k) = 0.5/hypm(k) * hypd(k)
      do kk=1,k-1
         ecref(kk,k) = 1./hypm(k) * hypd(kk)
      end do
   end do
!
! Reference hydrostatic integration matrix consistent with conversion
! term for energy conservation.  In href, 1st index = column; 
! 2nd index = row of matrix.
!
   do k = 1,plev
      do kk = k,plev
         href(kk,k) = ecref(k,kk)*hypd(kk)/hypd(k)
      end do
   end do
!
! Print statements
!
   if (lprint) then
      nkk = plev/13
      if (mod(plev,13).ne.0) nkk = nkk + 1
      write(6,*)' '
      write(6,*)'INITCOM: Hydrostatic matrix href'
      do kk=1,nkk
         ik1 = 1 + (kk-1)*13
         ik2 = min0( ik1+12, plev )
         write(6,9920) (k,k=ik1,ik2)
         do kkk=1,plev
            write(6,9910) kkk,(href(kkk,k),k=ik1,ik2)
         end do
      end do
      write(6,*)' '
      write(6,*)'INITCOM: Thermodynamic matrix ecref'
      do kk=1,nkk
         ik1 = 1 + (kk-1)*13
         ik2 = min0( ik1+12, plev )
         write(6,9920) (k,k=ik1,ik2)
         do kkk=1,plev
            write(6,9910) kkk,(ecref(kkk,k),k=ik1,ik2)
         end do
      end do
   end if
!
! Multiply href by r
!
   do k=1,plev
      do kk=1,plev
         href(kk,k) = href(kk,k)*rair
      end do
   end do
!
! Compute truncation parameters
!
   if (masterproc) then
      write(6,9950) ptrm,ptrn,ptrk
   end if
!
! Compute semi-implicit timestep constants (COMSPE)
!
   zdt = dtime
   if (.not.nlres) zdt = 0.5*zdt
!
!!   call settau(zdt)
!
! Determine whether full or reduced grid
!
   fullgrid = .true.
   do j=1,plat
      if (masterproc) then
         write(6,*)'nlon(',j,')=',nlon(j),' wnummax(',j,')=',wnummax(j)
      end if
      if (nlon(j).lt.plon) fullgrid = .false.
   end do

#if ( defined SCAM )
   do j=1,plat
      slat(j) = 1.0_r8 * sin(4.0_r8*atan(1.0_r8)*columnLat/180._r8)
      w(j)   = 2.0_r8/plat
      cs(j)  = 10._r8 - slat(j)*slat(j)
!         itmp = 2*pspt - 1
!         call phcs  (zalp    ,zdalp   ,itmp    ,zslat    )
!         call reordp(j       ,itmp    ,zalp    ,zdalp   )
   end do

!
! Latitude array (list of latitudes in radians)
!
   xlat = asin(slat(1))
   clat(1) = xlat

   clat(1)=columnLat*atan(1._r8)/45.
   latdeg(1) = clat(1)*45._r8/atan(1._r8)
   clon(1,1)   = 4.0_r8*atan(1._r8)*mod((columnLon+360._r8),360._r8)/180._r8
   londeg(1,1) = mod((columnLon+360._r8),360._r8)
!
! SCAM not yet able to handle reduced grid.
!
   if (.not. fullgrid) then
      call endrun ('INITCOM: SCAM not yet configured to handle reduced grid')
   end if
#else
!
! Compute constants related to Legendre transforms
! Compute and reorder ALP and DALP
!
!wjp 2007.06   allocate( alp  (pspt,plat/2) )
!wjp 2007.06   allocate( dalp (pspt,plat/2) )
!wjp 2007.06   do j=1,plat/2
!wjp 2007.06      zslat = slat(j)
!wjp 2007.06      itmp = 2*pspt - 1
!wjp 2007.06      call phcs  (zalp    ,zdalp   ,itmp    ,zslat    )
!wjp 2007.06      call reordp(j       ,itmp    ,zalp    ,zdalp   )
!wjp 2007.06   end do
!
! Copy and save local ALP and DALP
!
!wjp 2007.06  allocate( lalp  (lpspt,plat/2) )
!wjp 2007.06   allocate( ldalp (lpspt,plat/2) )
!wjp 2007.06   do j=1,plat/2
!wjp 2007.06      do lm=1,numm(iam)
!wjp 2007.06         m = locm(lm,iam)
!wjp 2007.06         mr = nstart(m)
!wjp 2007.06         lmr = lnstart(lm)
!wjp 2007.06         do n=1,nlen(m)
!wjp 2007.06            lalp(lmr+n,j) = alp(mr+n,j)
!wjp 2007.06            ldalp(lmr+n,j) = dalp(mr+n,j)
!wjp 2007.06         end do
!wjp 2007.06      end do
!wjp 2007.06   end do
!
! Mirror latitudes south of south pole
!
   lat = 1
   do j=j1-2,1,-1
      nlonex(j) = nlon(lat)
      lat = lat + 1
   end do
   nlonex(j1-1) = nlon(1)     ! south pole
!
! Real latitudes
!
   j = j1
   do lat=1,plat
      nlonex(j) = nlon(lat)
      j = j + 1
   end do
   nlonex(j1+plat) = nlon(plat)  ! north pole
!
! Mirror latitudes north of north pole
!
   lat = plat
   do j=j1+plat+1,platd
      nlonex(j) = nlon(lat)
      lat = lat - 1
   end do
!
! Longitude array
!
!!   pi = 4.0*atan(1.0)
   do lat=1,plat
      do i=1,nlon(lat)
         londeg(i,lat) = (i-1)*360./nlon(lat)
         clon(i,lat)   = (i-1)*2.0*pi/nlon(lat)
      end do
   end do

   do j=1,plat/2
      nmmax(j) = wnummax(j) + 1
   end do
   do m=1,pmmax
      do irow=1,plat/2
         if (nmmax(irow) .ge. m) then
            beglatpair(m) = irow
            goto 10
         end if
      end do
      call endrun ('INITCOM: Should not ever get here')
10    continue
   end do
!
! Set up trigonometric tables for fft
!
!wjp 2010.04.13   do j=1,plat
!wjp 2010.04.13      call set99(trig(1,j),ifax(1,j),nlon(j))
!wjp 2010.04.13   end do
!
! Set flag indicating dynamics grid is now defined.
! NOTE: this ASSUMES initcom is called after spmdinit.  The setting of nlon done here completes
! the definition of the dynamics grid.
!
   dyngrid_set = .true.
#endif

   return

9910 format( 1x,i3,13f9.5)
9920 format(/,      13i9)
9950 format(/,'     Truncation Parameters',/,'     NTRM = ',i4,/, &  !wjp 2007.06
     '     NTRN = ',i4,/,'     NTRK = ',i4,/)                        !wjp 2007.06

end subroutine initcom
