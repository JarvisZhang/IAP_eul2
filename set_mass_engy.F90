#include <misc.h>
#include <params.h>

subroutine set_mass_engy ( ztodt,lat, nlon, cwava, &
                       qfcst   ,qminus  ,etamid  , &
                       ps      ,u3      ,v3      , &
                       t3      ,flx_net ,phis    , &
                       psm1    ,u3m1    ,v3m1    , &
                       t3m1    ,hw2l    ,hw3l    , &
                       hwxl    ,engy1lat, engy2lat, difftlat  )
!
!-----------------------------------------------------------------------
! Purpose: Set coefficient used for mass and energy fixer
! Author: ZhangHe
! Completed: 2008.4.17
!-----------------------------------------------------------------------
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid,       only: plon, plev, plat,beglev,endlev,endlevp,iam
   use commap,       only: w     
   use constituents, only: pcnst, pnats
   use Dyn_const,    only: PMTOP
! ======================= zhh 2008.6.16 =========================
   use pmgrid,       only: masterproc,myid_z
! ======================= zhh 2008.6.16 =========================

   implicit none
!----------------------------------------------------------------------------------
#include <comqfl.h> ! tmass, tmass0, tmassf, qmassf ,fixmas ,qmass1, qmass2, pdela
!----------------------------------------------------------------------------------

!----------------------------------Arguments------------------------------------
! input:
   integer, intent(in) :: lat                 ! latitude index
   integer, intent(in) :: nlon                ! number of longitudes

   real(r8), intent(in) :: ztodt              ! timestep 
   real(r8), intent(in) :: cwava              ! normalization factor (1/g*plon)
   real(r8), intent(in) :: qfcst(plon,beglev:endlev,pcnst)  ! temporary q at time n+1
   real(r8), intent(in) :: qminus(plon,beglev:endlev,pcnst) ! q at time n
   real(r8), intent(in) :: etamid(plev)       ! vertical coords at midpoints

   real(r8), intent(in) :: ps(plon)        ! surface pressure at time n+1
   real(r8), intent(in) :: u3(plon,beglev:endlev)   ! u-wind at time n+1
   real(r8), intent(in) :: v3(plon,beglev:endlev)   ! v-wind at time n+1 
   real(r8), intent(in) :: t3(plon,beglev:endlev)   ! temperature at time n+1
   real(r8), intent(in) :: phis(plon)      ! surface geopotential at time n+1
   real(r8), intent(in) :: flx_net(plon)   ! net flux from physics at time n
! 
   real(r8), intent(in) :: psm1(plon)      ! surface pressure at time n
   real(r8), intent(in) :: u3m1(plon,beglev:endlev) ! u-wind at time n 
   real(r8), intent(in) :: v3m1(plon,beglev:endlev) ! v-wind at time n 
   real(r8), intent(in) :: t3m1(plon,beglev:endlev) ! temperature at time n 
!
! output:
   real(r8), intent(out) :: hw2l(pcnst)   ! component of slt global mass integrals (q)
   real(r8), intent(out) :: hw3l(pcnst)   ! component of slt global mass integrals (q)
   real(r8), intent(out) :: hwxl(pcnst,4) ! component of slt global mass integrals (others)

   real(r8), intent(out) :: engy1lat     ! total energy integral per latitude at time n-1
   real(r8), intent(out) :: engy2lat     ! total energy integral per latitude at time n
   real(r8), intent(out) :: difftlat     ! delta-temperature integral per latitude

!---------------------------Local workspace-----------------------------
   real(r8) :: pint(plon,beglev:endlev+1)    ! Pressure at model interfaces
   real(r8) :: pmid(plon,beglev:endlev)      ! Pressure at model levels
   real(r8) :: pdel(plon,beglev:endlev)      ! Layer thickness (pint(k+1) - pint(k))
!
   integer  :: i,m,k               ! indices
   real(r8) :: sum, tmp          ! accumulator
   real(r8) :: hcwavaw           ! 
   real(r8) :: dotprod           ! dot product
!-----------------------------------------------------------------------
!
! Accumulate mass integrals
!
   sum = 0.
   do i=1,nlon
      sum = sum + ( ps(i)-PMTOP )
   end do
  
   tmass(lat) = w(lat)*cwava*sum
! 
! Compute integrals
!     
   call plevs00(nlon,plon,plev,psm1,pint,pmid,pdel)
   call engy_te1 (cwava,w(lat),t3m1,u3m1,v3m1,phis,pdel, tmp  ,nlon)
   engy1lat = tmp 
!
! Include top/bottom flux integral to energy integral
!
   call flxint  (w(lat) ,flx_net ,tmp  ,nlon )
  if(myid_z.eq.0)  engy1lat = engy1lat + tmp *ztodt
!
! Calculate SLT moisture, constituent, energy, and temperature integrals
! 
   call plevs00(nlon,plon,plev,ps,pint,pmid,pdel)
!
   hcwavaw   = 0.5*cwava*w(lat)
   engy2lat = 0.
   difftlat = 0.
   do m=1,pcnst
      hw2l(m) = 0.
      hw3l(m) = 0.
      hwxl(m,1) = 0.
      hwxl(m,2) = 0.
      hwxl(m,3) = 0.
      hwxl(m,4) = 0.
      do k=beglev,endlev
         dotprod = 0.
         do i=1,nlon
!------------ zhh test 2013-10-08 -------------
!            print*, 'i,k,m =', i,k,m
!            print*, 'qfcst(i,k,m) =',qfcst(i,k,m)
!            print*, 'pdel(i,k) =',pdel(i,k)
!            print*, 'dotprod =', dotprod
!------------ zhh test 2013-10-08 -------------
            dotprod = dotprod + qfcst(i,k,m)*pdel(i,k)
         end do
         hw2l(m) = hw2l(m) + hcwavaw*dotprod
      end do
   end do


!
! ======================= zhh 2008.6.16 =========================
!   if ( lat == 64 ) then
!      print*, '------------ before calling engy_te ---------------------'
!      print*, 'w(lat) =', w(lat), 'engy2lat =', engy2lat
!      print*, 't3(1,26) =', t3(1,26), 'u3(1,26) =', u3(1,26)
!      print*, 'v3(1,26) =', v3(1,26), 'phis(1) =', phis(1)
!      print*, 'pdel(1,26) =', pdel(1,26), 'cwava =', cwava
!      print*, 'nlon =', nlon
!   end if
! ======================= zhh 2008.6.16 =========================
   call engy_te1  (cwava ,w(lat) ,t3  ,u3  ,v3 ,phis    ,pdel, engy2lat ,nlon)
! ======================= zhh 2008.6.16 =========================
!   if ( lat == 64 ) then
!      print*, '------------ after calling engy_te ---------------------'
!      print*, 'engy2lat =', engy2lat
!   end if
! ======================= zhh 2008.6.16 =========================
   call engy_tdif1(cwava ,w(lat) ,t3  ,t3m1             ,pdel, difftlat ,nlon)
!
   call qmassd1 (cwava, etamid, w(lat), qminus, qfcst, &
                pdel,  hw3l, nlon)

   if (pcnst.gt.1) then
      call xqmass (cwava, etamid, w(lat), qminus, qfcst, &
                   qminus, qfcst, pdel  , hwxl  , nlon      )
   end if
!
  return
end subroutine set_mass_engy

