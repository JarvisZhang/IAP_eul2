#include <misc.h>
#include <params.h>
!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  inital --- Define initial conditions for first run of case
!
! !INTERFACE:
subroutine inital

! !USES:
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid
   use inidat      , only: read_inidat
   use buffer
   use radae, only: initialize_radbuffer
   use pspect
   use prognostics
   use comsrf
   use dynamics_vars, only : dynamics_init
   use chem_surfvals, only: chem_surfvals_init
   use phys_grid, only: phys_grid_init
   use constituents, only: ppcnst, cnst_need_pdeldry
   use phys_buffer,  only: pbuf_allocate
   use time_manager, only: timemgr_init, get_step_size
   use filenames, only: ncdata, bnd_topo
#if (defined COUP_CSM)
   use ccsm_msg, only: initialize_ccsm_msg
#endif
   use ioFileMod
   use rgrid,        only: nlon
 ! ========================== zhh =============================
   use Dyn_const,    only: STFRAM
   use IAP_prog,     only: initialize_IAPprog  !zhh 2008.6.10
 ! ====================== 2007.04.18 ==========================

#if (defined SPMD)
!jjr    use spmd_dyn,     only: spmdbuf
#endif
    use scanslt,      only: scanslt_alloc
#if (defined BFB_CAM_SCAM_IOP )
    use history,         only: initialize_iop_history
#endif


!-----------------------------------------------------------------------
   implicit none
!------------------------------Parameters-------------------------------
#include <comctl.h>
!-----------------------------------------------------------------------
#include <comlun.h>
!-----------------------------------------------------------------------
#include <comqfl.h>
!-----------------------------------------------------------------------
   include 'netcdf.inc'

! !DESCRIPTION:
!
!   Define initial conditions for first run of case
! 
! !REVISION HISTORY:
!
!   92.06.01      Bath          Creation from CCM1
!   96.03.01      Acker         Modifications
!   96.04.01      Boville       Reviewed 
!   01.06.17      Sawyer        Added call to dynamics_init
!   01.07.12      Sawyer        Added arguments to dynamics_init
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   integer k                  ! indices
   character(len=256) locfn   ! local filename
   real(r8) :: dtime          ! timestep size
    integer n,i,lat            ! index
    real(r8) pdel(plon,plev)     ! pressure layer thickness
    real(r8) pint(plon,plevp)    ! pressure at model interfaces
    real(r8) pmid(plon,plev)     ! pressure at model levels

!
!-----------------------------------------------------------------------
!
! Obtain initial and topography datasets
!
   if (masterproc) then
      call getfil(ncdata, locfn)
      call wrap_open(locfn, NF_NOWRITE, ncid_ini)

! Backward compatibility: look for topography data on initial file if topo file name not provided.
      if (trim(bnd_topo) /= 'bnd_topo') then

      call wrap_open(locfn, NF_NOWRITE, ncid_topo)
         call getfil(bnd_topo, locfn)

      else
         ncid_topo = ncid_ini
      end if
   end if
!
! Check for consistent settings on initial dataset
!
   call readinitial(ncid_ini)

! Initialize time manager.
 

   call timemgr_init()
    call chem_surfvals_init()
   dtime = get_step_size()
!jjr
   call dynamics_init( dtime, iord, jord, nsplit, &
                       plon, plat, plev, ppcnst,  &
                       beglonxy, endlonxy,        &
                       beglatxy, endlatxy,        &
                       beglat,   endlat,          &
                       beglev,   endlev )

!jjr
! Initialize prognostics variables 
!
   call initialize_prognostics
   call initialize_IAPprog     !zhh 2008.6.10
   call scanslt_alloc()
!

        

! Initialize commons
!
   call initcom

!

! Define physics data structures
!
   call phys_grid_init


#if (defined COUP_CSM)
!
! Initialize ccsm arrays (must be done after phys_grid_init where
! begchunk and endchunk are defined
!

   call initialize_ccsm_msg

#endif
#if (defined SPMD)
! Allocate communication buffers for
! collective communications in realloc
! routines and in dp_coupling
!jjr   call spmdbuf ()
#endif

!
! Initialize buffer, comsrf, and radbuffer variables 
! (which must occur after the call to phys_grid_init)
!
   call pbuf_allocate('global')
   call initialize_buffer

   call initialize_comsrf

   call initialize_radbuffer

 
! Initialize ghg surface values before default initial distributions 
! are set in inidat.
! waccm_mozart requires after phys_grid_init
!jjr   call chem_surfvals_init   
!
! Read in initial data
#if (defined BFB_CAM_SCAM_IOP )
   call initialize_iop_history
#endif
!


!  Set constants used in dynamical calculations depend on model resolution
   call STFRAM     !zhh 2008.6.11

!jjr end

   call read_inidat


! Close the topographic dataset
! Backward compatibility: don't close if ncid_topo = ncid_ini
   if (masterproc .and. (ncid_topo /= ncid_ini)) call wrap_close(ncid_topo)
! If dry-type tracers are present, initialize pdeld
! First, set current time pressure arrays for model levels etc. to get pdel
!
      do lat=beglat,endlat
         call plevs0(nlon(lat), plon, plev, ps(1,lat,1), pint, pmid, pdel)
         if (  cnst_need_pdeldry ) then
            do k=1,plev !jjr
               do i=1,nlon(lat)
                  pdeld(i,k,lat,1) = pdel(i,k)*(1.-q3(i,k,1,lat,1))
               end do !i
            end do !k
         endif !cnst_need_pdeldry
      end do !lat
!
! Make all time levels of prognostics contain identical data.
! Fields to be convectively adjusted only *require* n3 time
! level since copy gets done in linems.
!
   do n=2,ptimelevels
      ps(:,:,n)     = ps(:,:,1)
      u3(:,:,:,n)   = u3(:,:,:,1)
      v3(:,:,:,n)   = v3(:,:,:,1)
      t3(:,:,:,n)   = t3(:,:,:,1)
      q3(1:plon,:,:,:,n) = q3(1:plon,:,:,:,1)
!!zhh      vort(:,:,:,n) = vort(:,:,:,1)
!!zhh      div(:,:,:,n)  = div(:,:,:,1)
      if (  cnst_need_pdeldry )  pdeld(1:plon,:,:,n) = pdeld(1:plon,:,:,1)
   end do

   call print_memusage ('post-inidat')

   return
!EOC
end subroutine inital
!----------------------------------------------------------------------- 
