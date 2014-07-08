#include <misc.h>
#include <params.h>
#if ( defined SCAM )
#include <max.h>
#endif
subroutine stepon
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Loop over time, calling driving routines for physics, dynamics, 
! transport
! 
! Method: 
! 
! Author: 
! Original version:  CCM1
! Standardized:      J. Rosinski, June 1992
! Reviewed:          B. Boville, D. Williamson, August 1992
! Reviewed:          B. Boville, D. Williamson, April 1996
! Restructured:      J. Truesdale, May 1999
! Modified:          ZhangHe, 2007.5.25
! Update:            ZhangHe, 2008.4.18, move calling of sub. STFRAM to 
!                                        sub. inital & sub. read_restart
!
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use history, only: wshist, wrapup
   use pmgrid,  only: plev, plat, plevp, plon, beglat, endlat, &
                      masterproc, twod_decomp,beglev,endlev,&
                      beglonxy,endlonxy,beglatxy,endlatxy,iam
   use scanslt, only: advection_state, scanslt_initial, scanslt_final, qfcst
   use rgrid,   only: nlon
!!   use prognostics, only: ps, u3, v3, t3, q3, qminus, div, dpsl, dpsm, omga, phis,  &
   use prognostics, only: ps, u3, v3, t3, q3, qminus, omga, phis, pdeld,    &
                          n3, n3m2, n3m1, ptimelevels, shift_time_indices   !zhh
   use restart, only: write_restart
   use constituents, only: pcnst,pnats
#if (defined COUP_CSM)
   use ccsm_msg, only: csmstop, ccsmfin
#endif
#if (defined SPMD)
   use spmd_dyn,      only: ijk_yz_to_xy,& !buf1, buf1win, buf2, buf2win, &
                            ikj_yz_to_xy,ijk_xy_to_yz   !  spmdbuf_siz,
    use mpishorthand, only : mpicom
    use parutilitiesmodule, only : sumop, parcollective
    use mod_comm, only: mp_sendirr, mp_recvirr

#endif


   use ppgrid,         only: begchunk, endchunk
   use physics_types,  only: physics_state, physics_tend
   use phys_buffer,    only: pbuf
   use dp_coupling,    only: d_p_coupling, p_d_coupling
   use commap,         only: clat
   use physconst, only: gravit
   use time_manager, only: advance_timestep, get_step_size, get_nstep, &
                           is_first_step, is_first_restart_step, &
                           is_last_step, is_end_curr_day, get_curr_date
   use times, only: time_dyn, time_phy, time_all, times_set, times_out
#if ( defined BFB_CAM_SCAM_IOP )
   use iop, only: init_iop_fields,fixmassav,betasav,alphasav,divt3dsav,divq3dsav,dqfx3savm1
#endif
#if ( defined SCAM )
   use scamMod, only :use_iop,use_pert_frc,doiopupdate,switch
   use iop, only: init_iop_fields,fixmassav,betasav,alphasav,divt3dsav,divq3dsav,dqfx3savm1
#endif
! ========================== zhh =============================
   use Dyn_const,  only: SIGL   !zhh 2008.1.2         
! ====================== 2007.11.11 ==========================
!
   implicit none
#if (defined SPMD)
#include <mpif.h>
#endif

!-----------------------------------------------------------------------
#include <comctl.h>
!-----------------------------------------------------------------------
#include <comhyb.h>
!-----------------------------------------------------------------------
#include <comlun.h>
!-----------------------------------------------------------------------
#if ( defined SCAM )
#include <comfrc.h>
!-----------------------------------------------------------------------
#endif
!
   type(physics_state), allocatable :: phys_state(:)
   type(physics_tend ), allocatable :: phys_tend(:)

   real(r8) gw(plat)                     ! Gaussian weights
   real(r8) detam(plev)                  ! intervals between vert full levs.
   real(r8) cwava(plat)                  ! weight applied to global integrals
   real(r8) etamid(plev)                 ! vertical coords at midpoints 
   real(r8), allocatable :: t2(:,:,:)    ! temp tendency
   real(r8), allocatable :: fu(:,:,:)    ! u wind tendency
   real(r8), allocatable :: fv(:,:,:)    ! v wind tendency
   real(r8) flx_net(plon,beglat:endlat)  ! net flux from physics
   real(r8) coslat(plon)
   real(r8) rcoslat(plon)
!!   real(r8) rpmid(plon,plev)
!!   real(r8) pdel(plon,plev)     ! pressure layer thickness
!!   real(r8) pint(plon,plevp)    ! pressure at model interfaces
!!   real(r8) pmid(plon,plev)     ! pressure at model levels
   real(r8) dtime               ! timestep size
   real(r8) ztodt               ! timestep size
   real(r8) :: wcstart, wcend   ! wallclock timestamp at start, end of timestep
   real(r8) :: usrstart, usrend ! user timestamp at start, end of timestep
   real(r8) :: sysstart, sysend ! sys timestamp at start, end of timestep
   type(advection_state) :: adv_state ! Advection state data
!
   integer i,k,lat,j,l,l1           ! longitude,level,latitude indices
   integer jd1,jd2              ! latitude index for pdeld
   integer iter
   integer :: yr, mon, day      ! year, month, and day components of a date
   integer :: ncsec             ! current time of day [seconds]
! ========================== zhh 2008.4.7 =============================
   real(r8)  :: time_begin, time_end   ! begin & end of the CPU time
   real(r8)  :: time_begin2, time_end2 ! begin & end of the CPU time
! ========================== zhh 2008.4.7 =============================
!-----------------------------------------------------------------------
! The following arrays are for secondary 2D x-y decomposition
!-----------------------------------------------------------------------
     real(r8), allocatable :: dummy3(:,:,:)

   real(r8), allocatable :: phisxy(:,:)    ! Surface geopotential
   real(r8), allocatable ::   psxy(:,:,:)    ! Surface pressure
   real(r8), allocatable :: omgaxy(:,:,:)  ! vertical pressure velocity
   real(r8), allocatable ::  u3xy(:,:,:,:)  ! Staggered grid winds, latitude
   real(r8), allocatable ::  v3xy(:,:,:,:)  ! Satggered grid winds, longitude
   real(r8), allocatable :: pdeldxy(:,:,:,:)  ! delta pressure
   real(r8), allocatable ::   t3xy(:,:,:,:)  ! virtual potential temperature
   real(r8), allocatable ::   q3xy(:,:,:,:,:)! Moisture and constituents
!   real(r8), allocatable ::  qminusxy(:,:,:,:)  ! finite-volume mean pk
!   real(r8), allocatable ::   hadvxy(:,:,:,:)  ! pe**cappa
!   real(r8), allocatable ::   pexy(:,:,:)  ! edge pressure
!   real(r8), allocatable :: pilnxy(:,:,:)  ! ln(pe)
!   real(r8), allocatable :: dudtxy(:,:,:)
!   real(r8), allocatable :: dvdtxy(:,:,:)
   real(r8), allocatable :: dummy3xy(:,:,:)
   integer:: tmp1

!-----------------------------------------------------------------------

! Other local variables


#if ( defined SCAM )
   integer dummy
   real(r8) omegapert
   real(r8) pertval              ! perturbation valie
   real(r8) drand48
   external drand48
#endif
!
! Externals
!
   logical, external :: rstwr  ! whether or not to write restart files
!
!-----------------------------------------------------------------------
!
   call t_startf ('stepon_startup')
   dtime = get_step_size()
!
! Define eta coordinates: Used for calculation etadot vertical velocity 
! for slt.
!
   do k=1,plev
! ========================== zhh =========================
!!      etamid(k) = hyam(k) + hybm(k)
      etamid(k) = SIGL(k)
! ======================= 2008.1.2 =======================
   end do
   call scanslt_initial( adv_state, etamid, gravit, gw, detam, cwava )
!
! Initial guess for trajectory midpoints in spherical coords.
! nstep = 0:  use arrival points as initial guess for trajectory midpoints.
! nstep > 0:  use calculated trajectory midpoints from previous time 
! step as first guess.
! NOTE:  reduce number of iters necessary for convergence after nstep = 1.
!
   if (is_first_step()) then
!
! Initialize for dynamics
!
      omga(:,:,:)=0.0

#if (defined SPMD)
      time_begin2=mpi_wtime()!wjp 2011.04.16
      call times_set         !wjp 2011.04.16
  
#endif  
      call init_trans_pd     !! jjr endmpp

      call initial_dyn       ! zhh 2007.8.23! jjr 0514

	  call init_trans_IAP    !! added by zhh  !jjr 0514

      CALL DIAGHI( 0 )       !zhh 2013-02-21

      if(masterproc) write(6,*) 'this is first timestep, omga is set to 0.0'
! ========================= delete by zhh ======================================
!!      do lat=beglat,endlat
!!#if ( !defined SCAM )
!!         do i=1,nlon(lat)
!!            coslat(i) = cos(clat(lat))
!!            rcoslat(i) = 1./coslat(i)
!!         end do
!!#endif
!     
! Set current time pressure arrays for model levels etc.
!
!!         call plevs0(nlon(lat), plon, plev, ps(1,lat,n3), pint, pmid, pdel)
!
!!         do k=1,plev
!!            do i=1,nlon(lat)
!!               rpmid(i,k) = 1./pmid(i,k)
!!            end do
!!         end do

!!#if ( ! defined SCAM )
!
! Calculate vertical pressure velocity (omga = dp/dt)
!
!!         call omcalc (rcoslat, div(1,1,lat,n3), u3(1,1,lat,n3), v3(1,1,lat,n3), dpsl(1,lat), &
!!                      dpsm(1,lat), pmid, pdel, rpmid   ,pint(1,plevp), &
!!                      omga(1,1,lat), nlon(lat))
!!#else
         
!!         omga(1,:,lat)=wfld(:)
!!#endif
!!      end do
! ========================= delete by zhh ======================================
   else if (is_first_restart_step()) then
      
#if (defined SPMD)
      time_begin2=mpi_wtime()!wjp 2011.04.16
      call times_set         !wjp 2011.04.16
#endif    
      call initial_dyn       !zhh 2007.5.25
	  if(masterproc) write(6,*) ' this is first restart timestep '

   end if

   allocate(phys_state(begchunk:endchunk))
   allocate(phys_tend(begchunk:endchunk))
   allocate(t2(plon,plev,beglat:endlat))
   allocate(fu(plon,plev,beglat:endlat))
   allocate(fv(plon,plev,beglat:endlat))
!----------------------------------------------------------
! Allocate variables for secondary 2D xy decomposition
!----------------------------------------------------------
  allocate (dummy3(plon,beglat:endlat,beglev:endlev))


   allocate (phisxy(beglonxy:endlonxy      , beglatxy:endlatxy))
   allocate (  psxy(beglonxy:endlonxy      , beglatxy:endlatxy, ptimelevels))
   allocate (omgaxy(beglonxy:endlonxy, plev, beglatxy:endlatxy))
   allocate ( u3xy(beglonxy:endlonxy,plev, beglatxy:endlatxy, ptimelevels))  ! now unghosted, was N1
   allocate ( v3xy(beglonxy:endlonxy,plev, beglatxy:endlatxy, ptimelevels))
   allocate (pdeldxy(beglonxy:endlonxy,plev, beglatxy:endlatxy, ptimelevels))
   allocate (  t3xy(beglonxy:endlonxy,plev, beglatxy:endlatxy, ptimelevels))
   allocate (  q3xy(beglonxy:endlonxy,plev,pcnst+pnats, beglatxy:endlatxy, ptimelevels))
!   allocate ( qminusxy(beglonxy:endlonxy,plev,pcnst, beglatxy:endlatxy))
!   allocate (  hadvxy(beglonxy:endlonxy, plev,pcnst,beglatxy:endlatxy))
!   allocate (  pexy(beglonxy:endlonxy, plevp, beglatxy:endlatxy))
!   allocate (pilnxy(beglonxy:endlonxy, plevp, beglatxy:endlatxy))
!   allocate (dudtxy(beglonxy:endlonxy, plev, beglatxy:endlatxy))
!   allocate (dvdtxy(beglonxy:endlonxy, plev, beglatxy:endlatxy))
   allocate (dummy3xy(beglonxy:endlonxy,beglatxy:endlatxy,plev))

!-----------------------------------------------------------------------

!
! Beginning of basic time step loop
!
   call t_stopf ('stepon_startup')

! Begin time loop.
!   if (twod_decomp .eq. 1) then
!JJR
!   else
!      do j = beglat,endlat
!         do i = 1,plon
!            phisxy(i,j) = phis(i,j)
!         enddo
!      enddo
!   endif


   do     !!!
      call t_startf('stepon_st')
      if (twod_decomp .eq. 1) then
#if defined( SPMD )
      do k = beglev,endlev
         do j = beglat,endlat
            do i = 1,plon
               dummy3(i,j,k) = phis(i,j)
            enddo
         enddo
      enddo
      call mp_sendirr(dummy3, ijk_yz_to_xy%SendDesc,              &
                      ijk_yz_to_xy%RecvDesc, dummy3xy)
      call mp_recvirr(dummy3xy, ijk_yz_to_xy%RecvDesc )
      do j = beglatxy,endlatxy
         do i = beglonxy,endlonxy
            phisxy(i,j) = dummy3xy(i,j,1)
         enddo
      enddo

      call mp_sendirr(omga, ikj_yz_to_xy%SendDesc,              &
                      ikj_yz_to_xy%RecvDesc, omgaxy)

      call mp_recvirr(omgaxy, ikj_yz_to_xy%RecvDesc )
!JJR      do l=1,ptimelevels
      call mp_sendirr(u3(:,:,:,n3m1), ikj_yz_to_xy%SendDesc,              &
                      ikj_yz_to_xy%RecvDesc, u3xy(:,:,:,n3m1))
     call mp_recvirr(u3xy(:,:,:,n3m1), ikj_yz_to_xy%RecvDesc )

      call mp_sendirr(v3(:,:,:,n3m1), ikj_yz_to_xy%SendDesc,              &
                      ikj_yz_to_xy%RecvDesc,v3xy(:,:,:,n3m1))
     call mp_recvirr(v3xy(:,:,:,n3m1), ikj_yz_to_xy%RecvDesc )

      call mp_sendirr(t3(:,:,:,n3m1), ikj_yz_to_xy%SendDesc,              &
                      ikj_yz_to_xy%RecvDesc, t3xy(:,:,:,n3m1))
      call mp_recvirr(t3xy(:,:,:,n3m1), ikj_yz_to_xy%RecvDesc )
      call mp_sendirr(pdeld(:,:,:,n3m1), ikj_yz_to_xy%SendDesc,              &
                      ikj_yz_to_xy%RecvDesc, pdeldxy(:,:,:,n3m1))
      call mp_recvirr(pdeldxy(:,:,:,n3m1), ikj_yz_to_xy%RecvDesc )

      do k = beglev,endlev
         do j = beglat,endlat
            do i = 1,plon
               dummy3(i,j,k) = ps(i,j,n3m1)
            enddo
         enddo
      enddo
     call mp_sendirr(dummy3, ijk_yz_to_xy%SendDesc,              &
                      ijk_yz_to_xy%RecvDesc, dummy3xy)

     call mp_recvirr(dummy3xy, ijk_yz_to_xy%RecvDesc )
      do j = beglatxy,endlatxy
         do i = beglonxy,endlonxy
            psxy(i,j,n3m1) = dummy3xy(i,j,1)
         enddo
      end do
       do l1=1,pcnst+pnats
        do k = beglev,endlev
         do j = beglat,endlat
            do i =1,plon
               dummy3(i,j,k) = q3(i,k,l1,j,n3m1)
            enddo
          enddo
        enddo
       call mp_sendirr(dummy3, ijk_yz_to_xy%SendDesc,              &
                      ijk_yz_to_xy%RecvDesc, dummy3xy)
       call mp_recvirr(dummy3xy, ijk_xy_to_yz%RecvDesc )
      do k=1,plev
       do j = beglatxy,endlatxy
         do i = beglonxy,endlonxy
            q3xy(i,k,l1,j,n3m1) = dummy3xy(i,j,k)
         enddo
      end do
      enddo

    enddo !jjr end l1


!JJR    end do
#endif
     else !jjr need define _xy
    psxy(:,:,n3m1)=ps(:,:,n3m1)
    t3xy(:,:,:,n3m1)=t3(:,:,:,n3m1)
    u3xy(:,:,:,n3m1)=u3(:,:,:,n3m1)
    v3xy(:,:,:,n3m1)= v3(:,:,:,n3m1)
    q3xy(:,:,:,:,n3m1)=q3(:,:,:,:,n3m1)
    omgaxy=omga
    phisxy=phis
    pdeldxy(:,:,:,n3m1)=pdeld(:,:,:,n3m1)
     end if

!
      if (masterproc .and. print_step_cost) then
         call t_stampf (wcstart, usrstart, sysstart)
      end if

!!      ztodt = 2.0*dtime
      ztodt = dtime         !!
!
!----------------------------------------------------------
! PHYSPKG  Call the Physics package
!----------------------------------------------------------
!
!!          do j=beglatxy,endlatxy   
!!            do i=beglonxy,endlonxy
!!             if ( j == 3.and.i==4 ) then
!!            write(6,*) '--------------- test jjr in stepon --------------------'
!!            write(6,*) 'v3(4,7,3,n3m1) =', v3(4,7,3,n3m1)
!!            write(6,*) 'q3(4,7,2,3,n3m1) =', q3(4,7,2,3,n3m1)
!!             endif
!!         enddo
!!         enddo
      call t_stopf('stepon_st')
!!      if (masterproc) print*, '--------------- before d_p_coupling ------------------'
      call t_startf('d_p_coupling')
      call d_p_coupling (psxy(:,:,n3m1), t3xy(:,:,:,n3m1), u3xy(:,:,:,n3m1), &
                         v3xy(:,:,:,n3m1), q3xy(:,:,:,:,n3m1), &
                         omgaxy, phisxy, phys_state, phys_tend, pbuf, pdeldxy(:,:,:,n3m1))

      call t_stopf('d_p_coupling')
!!      if (masterproc) print*, '--------------- after d_p_coupling ------------------'

      call t_startf('phys_driver')
! =========================== 2008.4.7 ===========================
! ========================= zhh 2008.9.11 ===============================
!!      do lat = beglat, endlat
!!         if ( lat==87 ) then
!!            print*, '--------------------- before physics ---------------------'
!!            print*, 'q3(62,8,1,87,n3m1) =', q3(62,8,1,87,n3m1)
!!            print*, 'q3(62,8,1,87,n3  ) =', q3(62,8,1,87,n3)
!!         end if
!!      end do
! ========================= zhh 2008.9.11 ===============================
#if (defined SPMD)
      time_begin=mpi_wtime()
#endif
      if (ideal_phys) then
         call phys_idealized(phys_state, phys_tend, ztodt, etamid)
      else if (adiabatic) then
         call phys_adiabatic(phys_state, phys_tend)
      else
!!         if (masterproc) print*, '--------------- before physpkg ------------------'
         call physpkg(phys_state, gw, ztodt, phys_tend, pbuf)
         if (masterproc) print*, '--------------- after physpkg ------------------'
      end if
#if (defined SPMD)
      time_end=mpi_wtime()
#endif
      time_phy=time_phy+time_end-time_begin
! =========================== 2008.4.7 ===========================
      call t_stopf('phys_driver')
   
      call t_startf('p_d_coupling')
      call p_d_coupling (phys_state, phys_tend, t2, fu, fv, flx_net, &
                         qminus(:,:,:,:),  q3(:,:,:,:,n3))

      call t_stopf('p_d_coupling') !end mpp jjr 2012.5.15
!!      if (masterproc) print*, '--------------- after p_d_coupling ------------------'
! ========================= zhh 2008.9.11 ===============================
!!      do lat = beglat, endlat
!!         if ( lat==87 ) then
!!            print*, '--------------------- after physics ---------------------'
!!            print*, 'q3(62,8,1,87,n3m1) =', q3(62,8,1,87,n3m1)
!!            print*, 'q3(62,8,1,87,n3  ) =', q3(62,8,1,87,n3)
!!            print*, 'qminus(62,8,1,87 ) =', qminus(62,8,1,87)
!!         end if
!!      end do
! ========================= zhh 2008.9.11 ===============================
!
!----------------------------------------------------------
! DYNPKG Call the Dynamics Package
!----------------------------------------------------------
!============================ zhh ===============================
! =========================== 2008.4.7 ===========================
#if (defined SPMD)
      time_begin=mpi_wtime()
#endif
      call t_startf('dynpkg')
!jjr t3,v3,q3,ps,omga,phis
!====================== revised by zhh 08.04.28 ==============================
!!           if(masterproc) write(6,*) 'test before dynpkg'


      call dynpkg(adv_state,t2,fu,fv,qminus,flx_net,t3,u3,v3,q3,ps,     &
                  omga,phis,etamid,cwava,detam,dtime,pcnst,n3,ptimelevels )
               if(masterproc) write(6,*) 'test after dynpkg'

!=============================================================================
      call t_stopf('dynpkg')
#if (defined SPMD)     
      time_end=mpi_wtime()
#endif
      time_dyn=time_dyn+time_end-time_begin
! ========================= zhh 2008.9.11 ===============================
!!      do lat = beglat, endlat
!!         if ( lat==87 ) then
!!            print*, '--------------------- after dynamics ---------------------'
!!            print*, 'q3(62,8,1,87,n3m1) =', q3(62,8,1,87,n3m1)
!!            print*, 'q3(62,8,1,87,n3  ) =', q3(62,8,1,87,n3)
!!         end if
!!      end do
! ========================= zhh 2008.9.11 ===============================
! =========================== 2008.4.7 ===========================
!!
!! Dump state variables to IC file
!!
      call diag_dynvar_ic (phis, ps(1,beglat,n3), t3(1,beglev,beglat,n3), u3(1,beglev,beglat,n3), &
                           v3(1,beglev,beglat,n3), q3(1,beglev,1,beglat,n3) )
!
! Shift time level pointers
!
!=============================================================================
      call shift_time_indices ()
!=============================================================================

      call t_startf('stepon_st')
      if (is_first_step() .or. is_first_restart_step()) then
         call print_memusage ('stepon after dynpkg')
      end if

! Set end of run flag.

#if ( ! defined COUP_CSM )
      if (is_last_step()) nlend = .true.
#else
      if (csmstop) then
         if ( masterproc ) write(6,*)'atm: Stopping at the end of this day'
         if (is_end_curr_day()) nlend = .true.
      end if
#endif

#if ( ! defined SCAM )
!
!------------------------------------------------------------------------------
! History and restart logic: Write and/or dispose history tapes if required
!------------------------------------------------------------------------------
!
      call t_startf ('wshist')
      call wshist ()
      call t_stopf ('wshist')
!
! Write restart file
!
      if (rstwr() .and. nrefrq /= 0) then
         call t_startf ('write_restart')
         call write_restart
         call t_stopf ('write_restart')
      end if
!
! Dispose necessary files
!
      call t_startf ('wrapup')
      call wrapup
      call t_stopf ('wrapup')

      if (masterproc .and. print_step_cost) then
         call t_stampf (wcend, usrend, sysend)
         write(6,'(a,3f8.3,a)')'Prv timestep wallclock, usr, sys=', &
                               wcend-wcstart, usrend-usrstart, sysend-sysstart, ' seconds'
      end if
!================================================================
! Advance timestep before returning to top of loop
!
      call advance_timestep()
      call get_curr_date(yr, mon, day, ncsec)
!================================================================
      if (is_last_step()) then 
         if(masterproc) write(6,*) 'At the end loop, the run time is:', yr, mon, day, ncsec
      end if
!================================================================

#else                 
! else if running SCAM then end timestep and dealloc
      nlend = .true.
#endif
      call t_stopf('stepon_st')
!
! Check for end of run
!
      if (nlend) then
!!         if(masterproc)  write(6,*) 'before calling scanslt_final' !zhh
         call scanslt_final( adv_state )
!!         if(masterproc)  write(6,*) 'success calling scanslt_final' !zhh
         call print_memusage ('End stepon')
!
         deallocate(phys_state)
         deallocate(phys_tend)
         deallocate(t2)
         deallocate(fu)
         deallocate(fv)
#ifdef COUP_CSM
         call ccsmfin
#endif
         return
      end if

#if (defined SPMD)
      time_end2=mpi_wtime()
      time_all=time_end2-time_begin2
      call times_out
#endif
   end do  ! End of timestep loop

end subroutine stepon
