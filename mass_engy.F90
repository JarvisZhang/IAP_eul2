#include <misc.h>
#include <params.h>

!=================================================================================
!!subroutine mass_engy (ztodt, cwava, etamid, flx_net)
subroutine mass_engy (ztodt, cwava, etamid, flx_net, fu, fv, t2, dqphy)  !zhh 2008.9.11
!---------------------------------------------------------------------------------
! Purpose: Driving routines for mass and energy correction
! Original version: scan2.F90 (cam3.1)
! Modified : ZhangHe
! Completed : 2008.4.18
! Update:  Zhanghe, 2008.6.8, add sub. realloc5 for parallel version
!---------------------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid,       only: plon, plat, plev, plevp, beglat, endlat, masterproc, iam&
                    ,beglev,endlev,beglatdyn,npr_z
   use prognostics,  only: n3, n3m1, ps, u3, v3, q3, t3, qminus,   &
                           phis, omga, hadv, pdeld
   use rgrid,        only: nlon
   use scanslt,      only: hw1lat, engy1lat, qfcst
   use constituents, only: cnst_need_pdeldry, pcnst, pnats
   use physconst,    only: cpair, rga
   use massfix,      only: hw1,hw2,hw3,alpha
   use IAP_grid,only: ie1,ex
#ifdef SPMD
   use mpishorthand, only: mpicom, mpir8
 use parutilitiesmodule, only : sumop, parcollective
#endif

   implicit none
!-----------------------------------------------------------------------
!------------------------------Commons----------------------------------
#include <comqfl.h> ! tmass, tmass0, tmassf, qmassf ,fixmas ,qmass1, qmass2, pdela
!-----------------------------------------------------------------------
#include <comctl.h> ! adiabatic, ideal_phys
!-----------------------------------------------------------------------
!
!--------------------------Input arguments------------------------------
!
   real(r8), intent(in) :: ztodt                ! timestep
   real(r8), intent(in) :: cwava(plat)          ! weight applied to global integrals
   real(r8), intent(in) :: etamid(plev)         ! vertical coords at midpoints 
   real(r8), intent(in) :: flx_net(plon,beglat:endlat) ! net flux from physics
! ========================= zhh 2008.9.11 ===============================
   real(r8), intent(in) :: fu(plon,beglev:endlev,beglat:endlat)  ! u-wind tendency due to physics 
   real(r8), intent(in) :: fv(plon,beglev:endlev,beglat:endlat)  ! v-wind tendency due to physics 
   real(r8), intent(in) :: t2(plon,beglev:endlev,beglat:endlat)  ! temperature tendency due to physics 
   real(r8), intent(in) :: dqphy(plon,beglev:endlev,beglat:endlat)  ! q tendency due to physics 
! ========================= zhh 2008.9.11 ===============================
!
!---------------------------Local workspace-----------------------------
!
   real(r8) :: engy1         ! component of global energy integral (for time step n)
   real(r8) :: engy2         ! component of global energy integral (for time step n+1)
   real(r8) :: difft         ! component of global delta-temp integral ( (n+1) - n )
   real(r8) :: beta          ! 
   real(r8) :: hwx(pcnst,4)  ! component of constituent global mass integral
   real(r8) :: engy2lat(plat)     ! lat contribution to total energy integral
   real(r8) :: difftlat(plat)     ! lat contribution to delta-temperature integral
   real(r8) :: hw2l(pcnst,plat)   !  | latitudinal contributions to the
   real(r8) :: hw3l(pcnst,plat)   ! <  components of global mass integrals
   real(r8) :: hwxl(pcnst,4,plat) !  |
   real(r8):: temp(5*pcnst+2)
!
   real(r8) :: residual          ! residual energy integral
   real(r8) :: tmp               ! accumulator
   integer  :: lat, m, n, jdry,tt,i         ! indices
!
!-----------------------------------------------------------------------
!
! Set coefficient used for mass and energy fixer
!
   engy1  = 0.
   engy2  = 0.
   difft  = 0.
   do m=1,pcnst
      hw1(m)=0.
      hw2(m) = 0.
      hw3(m) = 0.
      do n=1,4
         hwx(m,n) = 0.
      end do
   end do
      tmassf = 0.
       do m=1,5*pcnst+2
        temp(m)=0.
       end do

   do lat=beglat,endlat
!        do m=beglev,endlev
!           do n=1,pcnst
!             do i=1,plon
!          print*,'test qminus', qminus(i,m,n,lat),qfcst(i,m,n,lat),i,m,n,lat
!           end do
!        end do 
!       end do
      call set_mass_engy (ztodt, lat, nlon(lat), cwava(lat),  &
                   qfcst(:,:,:,lat), qminus(:,:,:,lat), etamid, &
                   ps(:,lat,n3), u3(:,:,lat,n3), v3(:,:,lat,n3), &
                   t3(:,:,lat,n3), flx_net(1,lat), phis(1,lat), & 
                   ps(:,lat,n3m1), u3(:,:,lat,n3m1), v3(:,:,lat,n3m1), &
                   t3(:,:,lat,n3m1), hw2l(1,lat), hw3l(1,lat), &
                   hwxl(:,:,lat), engy1lat(lat), engy2lat(lat), difftlat(lat) )
!      engy1  = engy1  + engy1lat(lat)
      temp(1)=temp(1)+engy1lat(lat)
!      engy2  = engy2  + engy2lat(lat)
      temp(2)=temp(2)+ engy2lat(lat)
!      difft  = difft  + difftlat(lat)
      temp(3)=temp(3)+ difftlat(lat)
!      hw1(1) = hw1(1) + hw1lat(1,lat)
      temp(4)=temp(4)+hw1lat(1,lat)
!      hw2(1) = hw2(1) + hw2l(1,lat)
      temp(5)=temp(5)+hw2l(1,lat)
!      hw3(1) = hw3(1) + hw3l(1,lat)
      temp(6)=temp(6)+hw3l(1,lat)
!      tmassf = tmassf + tmass(lat)
      temp(7)=temp(7)+ tmass(lat)
        tt=7
     do m = 2,pcnst
        tt=tt+1
!jjr      do lat=1,plat
!         hw1(m) = hw1(m) + hw1lat(m,lat)
       temp(tt)=temp(tt)+ hw1lat(m,lat)

      do n = 1,4
         tt=tt+1
!jjr         do lat=1,plat
!            hwx(m,n) = hwx(m,n) + hwxl(m,n,lat)
       temp(tt)=temp(tt)+hwxl(m,n,lat)
!jjr         end do
      end do
     end do

   end do
       
#ifdef SPMD
#ifdef TIMING_BARRIERS
   call t_startf ('sync_realloc5')
   call mpibarrier (mpicom)
   call t_stopf ('sync_realloc5')
#endif
        call parcollective(mpicom,sumop,5*pcnst+2,temp)
!jjr   call t_startf('realloc5')
!jjr   call realloc5 (hw2l   ,hw3l   ,tmass   ,hw1lat  , &
!jjr                  AP4_1x1IAP4_1x1hwxl  ,engy1lat,engy2lat, difftlat   )
!jjr   call t_stopf('realloc5')
#endif
!jjr
!           print*,'test engy 3',temp(n),n


      engy1  = temp(1) 
      engy2  = temp(2)
      difft  = temp(3)
      hw1(1) = temp(4)
      hw2(1) = temp(5)
      hw3(1) = temp(6)
      tmassf = temp(7)/npr_z
        tt=7
   do m = 2,pcnst
        tt=tt+1
         hw1(m) =temp(tt)

      do n = 1,4
         tt=tt+1
            hwx(m,n) = temp(tt)
      end do
  end do
!jjr
! Initialize moisture, mass, energy, and temperature integrals
!
!jjr   hw1(1) = 0.
!   engy1  = 0.
!  engy2  = 0.
!   difft  = 0.
!   do m=1,pcnst
!      hw2(m) = 0.
!      hw3(m) = 0.
!      do n=1,4
!         hwx(m,n) = 0.
!      end do
!   end do
!
! Sum water and energy integrals over latitudes
!JJR mpi_allreduce
!
!   do lat=1,plat
!      engy1  = engy1  + engy1lat(lat)
!      engy2  = engy2  + engy2lat(lat)
!      difft  = difft  + difftlat(lat)
!      hw1(1) = hw1(1) + hw1lat(1,lat)
!      hw2(1) = hw2(1) + hw2l(1,lat)
!      hw3(1) = hw3(1) + hw3l(1,lat)
!   end do
!
! Accumulate and normalize global integrals for mass fixer (dry mass of
! atmosphere is held constant).
!
!   tmassf = 0.
!   do lat = 1, plat
!      tmassf = tmassf + tmass(lat)
!   end do
!
! Compute atmospheric mass fixer coefficient
!
tmassf = tmassf * 0.5
   qmassf = hw1(1)
   if (adiabatic .or. ideal_phys) then
      fixmas = tmass0/tmassf
   else
      fixmas = (tmass0 + qmassf)/tmassf
   end if
!
! Compute alpha for water ONLY
!
   hw2(1)    = fixmas*hw2(1)
   hw3(1)    = fixmas*hw3(1)
   if(hw3(1) .ne. 0.) then
      alpha(1)  = ( hw1(1) - hw2(1) )/hw3(1)
   else
      alpha(1)  = 1.
   endif
!
! Compute alpha for non-water constituents
!
   do m = 2,pcnst
!jjr      hw1(m) = 0.
!jjr      do lat=1,plat
!jjr         hw1(m) = hw1(m) + hw1lat(m,lat)
!      end do
!      do n = 1,4
!         do lat=1,plat
!            hwx(m,n) = hwx(m,n) + hwxl(m,n,lat)
!         end do
!      end do
      hw2(m) = hwx(m,1) - alpha(1)*hwx(m,2)
      hw3(m) = hwx(m,3) - alpha(1)*hwx(m,4)
      hw2(m) = fixmas*hw2(m)
      hw3(m) = fixmas*hw3(m)
      if(hw3(m) .ne. 0.) then
         alpha(m)  = ( hw1(m) - hw2(m) )/hw3(m)
      else
         alpha(m)  = 1.
      end if
   end do
!
! Compute beta for energy
!
   engy2    = fixmas*engy2
   difft    = fixmas*difft
   residual = (engy2 - engy1)/ztodt
   if(difft .ne. 0.) then
     beta = -residual*ztodt/(cpair*difft)
   else
     beta = 0.
   endif

!=========================== zhh 2008.4.23 =========================
!   if (masterproc) then
!      print*, 'fixmas =', fixmas
!      print*, 'engy2 =', engy2
!      print*, 'engy1 =', engy1
!      print*, 'beta =', beta
!   end if
!=========================== zhh 2008.4.23 =========================

   call t_startf ('mass_engy_fix')
!
!wjp 2011.04.16 !$OMP PARALLEL DO PRIVATE (LAT)
!


   do lat=beglat,endlat
      if ( cnst_need_pdeldry ) then
         jdry = lat
      else
         jdry = 1
      endif
      call mass_engy_fix (ztodt,     lat,   nlon(lat),        u3(:,:,lat,n3m1),    &
                          u3(:,:,lat,n3),   v3(:,:,lat,n3m1), v3(:,:,lat,n3),      &
                          t3(:,:,lat,n3m1), t3(:,:,lat,n3),   q3(:,:,:,lat,n3m1),  &
                          q3(:,:,:,lat,n3), ps(1,lat,n3m1),   ps(1,lat,n3),        &
                          alpha,            etamid,           qfcst(:,:,:,lat),    &
                          qminus(:,:,:,lat), beta,            hadv(:,:,:,lat) ,    &
! =================================== zhh 2008.9.11 ==================================
                          fu(:,:,lat),      fv(:,:,lat),      t2(:,:,lat),         &
                          dqphy(:,:,lat),   omga(:,:,lat),    pdeld(:,:,jdry,n3) )
! =================================== zhh 2008.9.11 ==================================

!-------------- zhh 2013-10-09 --------------
!      print*, 'lat =', lat
!      call mass_engy_fix (ztodt,     lat,   nlon(lat),    u3(1,beglev:endlev,lat,n3m1),    &
!              u3(1,beglev:endlev,lat,n3), v3(1,beglev:endlev,lat,n3m1), v3(1,beglev:endlev,lat,n3),      &
!              t3(1,beglev:endlev,lat,n3m1), t3(1,beglev:endlev,lat,n3), q3(1,beglev:endlev,1:pcnst+pnats,lat,n3m1),  &
!              q3(1,beglev:endlev,1:pcnst+pnats,lat,n3), ps(1,lat,n3m1),   ps(1,lat,n3),        &
!              alpha,            etamid,           qfcst(:,:,:,lat),    &
!              qminus(1,beglev:endlev,1:pcnst,lat), beta, hadv(1,beglev:endlev,1:pcnst,lat) ,    &
!              fu(1,beglev:endlev,lat), fv(1,beglev:endlev,lat), t2(1,beglev:endlev,lat),         &
!              dqphy(1,beglev:endlev,lat), omga(1,beglev:endlev,lat), pdeld(1,beglev:endlev,jdry,n3) )
!-------------- zhh 2013-10-09 --------------
   end do

   call t_stopf ('mass_engy_fix')
!
   return
end subroutine mass_engy
!
#ifdef SPMD1
subroutine realloc5 (hw2l   ,hw3l   ,tmass   ,hw1lat  , &
                     hwxl  ,engy1lat,engy2lat, difftlat      )
!-----------------------------------------------------------------------
!
! Purpose: Reallocation routine for slt variables.
!
! Method: MPI_Allgatherv (or point-to-point implementation)
! 
! Author:  J. Rosinski
! Standardized:      J. Rosinski, Oct 1995
!                    J. Truesdale, Feb. 1996
! Modified: P. Worley, December 2003, October 2004
!           ZhangHe, 2008.6.8
!
! $Id: scan2.F90,v 1.10.2.13 2005/03/02 20:52:29 bundy Exp $
! $Author: bundy $
!
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid, only: iam, numlats, plat, beglat, endlat
   use mpishorthand, only: mpicom, mpir8
   use spmd_dyn
   use constituents, only: pcnst
!---------------------------------Parameters----------------------------------
   integer, parameter :: msgtag  = 5000
!---------------------------------Commons-------------------------------------
#include <comsta.h>
!-----------------------------------------------------------------------
!
! Input arguments
!
!!   real(r8), intent(inout) :: hw2al(pcnst,plat)
   real(r8), intent(inout) :: hw2l(pcnst,plat)
!!   real(r8), intent(inout) :: hw3al(pcnst,plat)
   real(r8), intent(inout) :: hw3l(pcnst,plat)
   real(r8), intent(inout) :: tmass (plat)
   real(r8), intent(inout) :: hw1lat(pcnst,plat)
!!   real(r8), intent(inout) :: hwxal(pcnst,4,plat)
   real(r8), intent(inout) :: hwxl(pcnst,4,plat)
!                                                ! -
   real(r8), intent(inout)   :: engy1lat (plat)  ! lat contribution to total energy (n)
!!   real(r8), intent(inout)   :: engy2alat(plat)  ! lat contribution to total energy (n+1)
   real(r8), intent(inout)   :: engy2lat(plat)  ! lat contribution to total energy (n+1)
!!   real(r8), intent(inout)   :: difftalat(plat)  ! lat contribution to delta-T integral
   real(r8), intent(inout)   :: difftlat(plat)  ! lat contribution to delta-T integral
!
!---------------------------Local workspace-----------------------------
!
   integer procid
   integer bufpos
   integer procj
   integer step, i, j, m, jstrt
   integer beglat_p, endlat_p, numlats_p, jstrt_p
!=================== zhh 2008.6.8 ===================
   integer npc1    ! npc1 = 3 (hw1lat, hw2l, hw3l)
   integer npc4    ! npc4 = 1 (hwxl)
   integer npc     ! npc = npc1 + npc4*4
   integer nme     ! nme = 4 (tmass, engy1lat, engy2lat, difftlat)
   integer ntl     ! ntl = npc*pcnst + nme
!=================== zhh 2008.6.8 ===================
!
   logical, save :: first = .true.
   integer, save :: sndcnt
   integer, allocatable, save :: sndcnts(:), sdispls(:)
   integer, allocatable, save :: rcvcnts(:), rdispls(:)
   integer, allocatable, save :: pdispls(:)
!-----------------------------------------------------------------------
   if (first) then
!=================== zhh 2008.6.8 ===================
      npc1 = 3
      npc4 = 1
      npc = npc1 + npc4*4
      nme = 4
      ntl = npc*pcnst + nme
!=================== zhh 2008.6.8 ===================
! Compute send/recv/put counts and displacements
      allocate(sndcnts(0:npes-1))
      allocate(sdispls(0:npes-1))
      allocate(rcvcnts(0:npes-1))
      allocate(rdispls(0:npes-1))
      allocate(pdispls(0:npes-1))
!
! Compute send count
!!      sndcnt = (pcnst*(5 + 2*4) + 6)*numlats
      sndcnt = ntl*numlats     !zhh 2008.6.8
      sndcnts(:) = 0
      do step=1,allgather_steps
         procid = allgather_proc(step)
         sndcnts(procid) = sndcnt
      enddo
!   
      sdispls(0) = 0
      do procid=1,npes-1
        sdispls(procid) = 0
      enddo
!
! Compute recv counts and displacements
      rcvcnts(:) = 0
      do step=1,allgather_steps
         procid = allgather_proc(step)
!!         rcvcnts(procid) = (pcnst*(5 + 2*4) + 6)*nlat_p(procid)
         rcvcnts(procid) = ntl*nlat_p(procid)  !zhh 2008.6.8
      enddo
!!      rcvcnts(iam) = (pcnst*(5 + 2*4) + 6)*numlats
      rcvcnts(iam) = ntl*numlats     !zhh 2008.6.8
!   
      rdispls(0) = 0
      do procid=1,npes-1
        rdispls(procid) = rdispls(procid-1) + rcvcnts(procid-1)
      enddo
!
      call mpialltoallint(rdispls, 1, pdispls, 1, mpicom)
!
      first = .false.
   endif
!
! Fill send buffer
   jstrt = beglat - 1
   bufpos = 0
! tmass
   do j=1,numlats
      buf1(bufpos+j) = tmass(jstrt+j)
   enddo
   bufpos = bufpos + numlats
! engy1lat
   do j=1,numlats
      buf1(bufpos+j) = engy1lat(jstrt+j)
   enddo
   bufpos = bufpos + numlats
! engy2alat
!!   do j=1,numlats
!!      buf1(bufpos+j) = engy2alat(jstrt+j)
!!   enddo
!!   bufpos = bufpos + numlats
!
! engy2lat
   do j=1,numlats
      buf1(bufpos+j) = engy2lat(jstrt+j)
   enddo
   bufpos = bufpos + numlats
! difftalat
!!   do j=1,numlats
!!      buf1(bufpos+j) = difftalat(jstrt+j)
!!   enddo
!!   bufpos = bufpos + numlats
!
! difftlat
   do j=1,numlats
      buf1(bufpos+j) = difftlat(jstrt+j)
   enddo
   bufpos = bufpos + numlats
!hw1lat
   do j=beglat,endlat
      do m=1,pcnst
         buf1(bufpos+m) = hw1lat(m,j)
      enddo
      bufpos = bufpos + pcnst
   enddo
!hw2al
!!   do j=beglat,endlat
!!      do m=1,pcnst
!!         buf1(bufpos+m) = hw2al(m,j)
!!      enddo
!!      bufpos = bufpos + pcnst
!!   enddo
!
!hw2l
   do j=beglat,endlat
      do m=1,pcnst
         buf1(bufpos+m) = hw2l(m,j)
      enddo
      bufpos = bufpos + pcnst
   enddo
!hw3al
!!   do j=beglat,endlat
!!      do m=1,pcnst
!!         buf1(bufpos+m) = hw3al(m,j)
!!      enddo
!!      bufpos = bufpos + pcnst
!!   enddo
!
!hw3l
   do j=beglat,endlat
      do m=1,pcnst
         buf1(bufpos+m) = hw3l(m,j)
      enddo
      bufpos = bufpos + pcnst
   enddo
!hwxal
!!   do j=beglat,endlat
!!      do i=1,4
!!         do m=1,pcnst
!!            buf1(bufpos+m) = hwxal(m,i,j)
!!         enddo
!!         bufpos = bufpos + pcnst
!!      enddo
!!   enddo
!
!hwxl
   do j=beglat,endlat
      do i=1,4
         do m=1,pcnst
            buf1(bufpos+m) = hwxl(m,i,j)
         enddo
         bufpos = bufpos + pcnst
      enddo
   enddo
!
! Gather the data
!
   if (dyn_allgather .eq. 0) then
      call mpiallgatherv(buf1, sndcnt, mpir8, &
                         buf2, rcvcnts, rdispls, mpir8, &
                         mpicom)
   else
      call altalltoallv(dyn_allgather, iam, npes, &
                        allgather_steps, allgather_proc, &
                        buf1, spmdbuf_siz, sndcnts, sdispls, mpir8, &
                        buf2, spmdbuf_siz, rcvcnts, rdispls, mpir8, &
                        msgtag, pdispls, mpir8, buf2win, mpicom)
   endif
!
! Copy out of message buffers
!
!wjp 2011.04.16 !$OMP PARALLEL DO PRIVATE (STEP, PROCID, BEGLAT_P, ENDLAT_P, NUMLATS_P, BUFPOS, JSTRT_P, I, J, M)
!CSD$ PARALLEL DO PRIVATE (STEP, PROCID, BEGLAT_P, ENDLAT_P, NUMLATS_P, BUFPOS, JSTRT_P, I, J, M)
   do step=1,allgather_steps
      procid = allgather_proc(step)
      beglat_p = cut(1,procid)
      endlat_p = cut(2,procid)
      numlats_p = nlat_p(procid)
      bufpos = rdispls(procid)
! tmass
      jstrt_p  = beglat_p - 1
      do j=1,numlats_p
         tmass(jstrt_p+j) = buf2(bufpos+j)
      enddo
      bufpos = bufpos + numlats_p
! engy1lat
      jstrt_p  = beglat_p - 1
      do j=1,numlats_p
         engy1lat(jstrt_p+j) = buf2(bufpos+j)
      enddo
      bufpos = bufpos + numlats_p
! engy2alat
!!      jstrt_p  = beglat_p - 1
!!      do j=1,numlats_p
!!         engy2alat(jstrt_p+j) = buf2(bufpos+j)
!!      enddo
!!      bufpos = bufpos + numlats_p
! engy2lat
      jstrt_p  = beglat_p - 1
      do j=1,numlats_p
         engy2lat(jstrt_p+j) = buf2(bufpos+j)  
      enddo
      bufpos = bufpos + numlats_p
! difftalat
!!      jstrt_p  = beglat_p - 1
!!      do j=1,numlats_p
!!         difftalat(jstrt_p+j) = buf2(bufpos+j)
!!      enddo
!!      bufpos = bufpos + numlats_p
! difftlat
      jstrt_p  = beglat_p - 1
      do j=1,numlats_p
         difftlat(jstrt_p+j) = buf2(bufpos+j)
      enddo
      bufpos = bufpos + numlats_p
! hw1lat
      do j=beglat_p,endlat_p
         do m=1,pcnst
            hw1lat(m,j) = buf2(bufpos+m)
         enddo
         bufpos = bufpos + pcnst
      enddo
! hw2al
!!      do j=beglat_p,endlat_p
!!         do m=1,pcnst
!!            hw2al(m,j) = buf2(bufpos+m)
!!         enddo
!!         bufpos = bufpos + pcnst
!!      enddo
! hw2l
      do j=beglat_p,endlat_p
         do m=1,pcnst
            hw2l(m,j) = buf2(bufpos+m)
         enddo
         bufpos = bufpos + pcnst
      enddo
! hw3al
!!      do j=beglat_p,endlat_p
!!         do m=1,pcnst
!!            hw3al(m,j) = buf2(bufpos+m)
!!         enddo
!!         bufpos = bufpos + pcnst
!!      enddo
! hw3l
      do j=beglat_p,endlat_p
         do m=1,pcnst
            hw3l(m,j) = buf2(bufpos+m)
         enddo
         bufpos = bufpos + pcnst
      enddo
! hwxal
!!      do j=beglat_p,endlat_p
!!         do i=1,4
!!            do m=1,pcnst
!!               hwxal(m,i,j) = buf2(bufpos+m)
!!            enddo
!!            bufpos = bufpos + pcnst
!!         enddo
!!      enddo
! hwxl
      do j=beglat_p,endlat_p
         do i=1,4
            do m=1,pcnst
               hwxl(m,i,j) = buf2(bufpos+m)
            enddo
            bufpos = bufpos + pcnst
         enddo
      enddo
!
   end do
!CSD$ END PARALLEL DO
!
   return
end subroutine realloc5
#endif
