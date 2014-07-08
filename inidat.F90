#include <misc.h>
#include <params.h>
#if ( defined SCAM )
#include <max.h>
#endif

module inidat

!----------------------------------------------------------------------- 
! 
! Purpose: Read initial dataset and process fields as appropriate
!
! Method: Read and process one field at a time
! 
! Author: J. Olson  May 2004
! 
!-----------------------------------------------------------------------

   use pmgrid
   use rgrid
   use prognostics
   use ncdio_atm
   use shr_kind_mod, only: r8 => shr_kind_r8
   use abortutils  , only: endrun
#if ( defined SPMD )
   use mpishorthand
#endif
#if ( defined SCAM )
   use getnetcdfdata
   use scamMod, only :initlonidx, initlatidx, initTimeIdx, modelfile, analysisfile
#endif

   implicit none

#if ( defined SCAM )
#include <runtype.h>
#endif

#include <comlun.h>

PRIVATE
   include 'netcdf.inc'
!
! Public interfaces
!
   public :: read_inidat
!
! Private module data
!
   integer ixcldice,ixcldliq  ! indices into q3 array for cloud liq and cloud ice
   real(r8), allocatable :: ps_tmp  (:,:  )
   real(r8), allocatable :: phis_tmp(:,:  )
   real(r8), allocatable :: q3_tmp  (:,:,:)
   real(r8), allocatable :: cld_liq (:,:,:)
   real(r8), allocatable :: cld_ice (:,:,:)
   real(r8), allocatable :: t3_tmp  (:,:,:)
   real(r8), allocatable :: arr3d_a (:,:,:)
   real(r8), allocatable :: arr3d_b (:,:,:)

   logical readvar            ! inquiry flag:  true => variable exists on netCDF file

   logical :: fill_ends          ! For SCAM, flag if ends are filled or not
   integer :: STATUS             ! For SCAM, Status of NetCDF operation
   logical :: have_surfdat       ! If have the surface dataset or not

contains

   subroutine read_inidat(datatype)
!
!-----------------------------------------------------------------------
!
! Purpose:
! Read initial dataset and spectrally truncate as appropriate.
!
!-----------------------------------------------------------------------
!
! $Id: inidat.F90,v 1.20.4.26 2005/03/08 17:06:41 eaton Exp $
! $Author: eaton $
!
!-----------------------------------------------------------------------
!
    use buffer
    use comsrf    ! tsice,  
    use phys_grid,    only: scatter_field_to_chunk, gather_chunk_to_field
    use phys_buffer,  only: pbuf, pbuf_times, pbuf_get_fld_idx
    use constituents, only: ppcnst, cnst_name, cnst_read_iv, cnst_get_ind

#include <comctl.h>
!
! Arguments
!
   integer, optional, intent(in) :: datatype !  in:  model or analysis dataset
   integer type
!
!---------------------------Local workspace-----------------------------
!
    integer i,c,m,n,j                         ! indices
    integer ncol

    real(r8), pointer, dimension(:,:,:,:) :: cldptr
    real(r8), pointer, dimension(:,:    ) :: arr2d_tmp
    real(r8), pointer, dimension(:,:    ) :: arr2d
    character*16 fieldname                  ! field name

    character*16 :: subname='READ_INIDAT'   ! subroutine name
!
!-----------------------------------------------------------------------
!     May 2004 revision described below (Olson)
!-----------------------------------------------------------------------
!
! This routine reads and processes fields one at a time to minimize 
! memory usage.
!
!   State fields (including PHIS) are read into a global array on 
!     masterproc, processed, and scattered to all processors on the
!     appropriate grid 
!
!   Physics fields are read in and scattered immediately to all
!     processors on the physics grid.
!
!-----------------------------------------------------------------------
!
#if (defined SCAM )
   if (.not.present (datatype)) then
     write(6,*)'must specify datatype when calling this routine as part of SCAM'
     call endrun
   else
     type=datatype
   endif

   call read_inidat_scam(type)
#endif
!
!-------------------------------------
! Allocate memory for temporary arrays
!-------------------------------------
!
! Note if not masterproc still might need to allocate array for spmd case
! since each processor calls MPI_scatter 
!
    allocate ( ps_tmp  (plon,plat     ) )
    allocate ( phis_tmp(plon,plat     ) )
    allocate ( q3_tmp  (plon,plev,plat) )
    allocate ( cld_liq (plon,plev,plat) )
    allocate ( cld_ice (plon,plev,plat) )
    allocate ( t3_tmp  (plon,plev,plat) )
!
!---------------------
! Read required fields
!---------------------

!
!-----------
! 3-D fields
!-----------
!jjr
    allocate ( arr3d_a (plon,plev,plat) )
    allocate ( arr3d_b (plon,plat,plev) )
!jjr
    if(masterproc) then

       fieldname = 'U'
       call infld(fieldname, ncid_ini, 'lon', 'lev', 'lat', 1, plon, 1, plev, 1, plat, &
                                                     arr3d_a, readvar, grid_map='global')
       if(.not. readvar) call endrun()
    end if
      call process_inidat('U')
    if(masterproc) then
     
        fieldname = 'V'
       call infld(fieldname, ncid_ini, 'lon', 'lev', 'lat', 1, plon, 1, plev, 1, plat, &
                                                     arr3d_a, readvar, grid_map='global')
       if(.not. readvar) call endrun()

     end if
     call process_inidat('V')
!!
    if (masterproc) then
       write (6,*) 'Success calling process_inidat(UV) '
    endif
!!    
    if(masterproc) then

       fieldname = 'T'
       call infld(fieldname, ncid_ini, 'lon', 'lev', 'lat', 1, plon, 1, plev, 1, plat, &
                                                      t3_tmp, readvar, grid_map='global')
       if(.not. readvar) call endrun()

    end if
    call process_inidat('T')
!
! Get indices for cloud liquid and cloud ice so that these
! fields can be saved off for later use
!
    call cnst_get_ind('CLDLIQ', ixcldliq)
    call cnst_get_ind('CLDICE', ixcldice)
!
! Constituents (read and process one at a time)
!
    do m = 1,ppcnst

       if(masterproc) then
          readvar   = .false.
          fieldname = cnst_name(m)
          if(cnst_read_iv(m)) &
             call infld(fieldname, ncid_ini, 'lon', 'lev', 'lat', 1, plon, 1, plev, 1, plat, &
                                                           arr3d_a, readvar, grid_map='global')
       end if
       call process_inidat('CONSTS', m_cnst=m)

    end do

    deallocate ( arr3d_a  )
    deallocate ( arr3d_b  )
!
!-----------
! 2-D fields
!-----------
!
    if(masterproc) then

       fieldname = 'PHIS'
       readvar   = .false.
       if (ideal_phys .or. aqua_planet) then
          phis_tmp(:,:) = 0.
       else
! ============================== zhh 2007.11.8 ================================
          call infld(fieldname, ncid_topo, 'lon', 'lat', 1, plon, 1, plat, &
                                       phis_tmp, readvar, grid_map='global')
! ============================== zhh 2007.11.8 ================================
       end if

    end if
    call process_inidat('PHIS')

    if(masterproc) then

       fieldname = 'PS'
       call infld(fieldname, ncid_ini, 'lon', 'lat', 1, plon, 1, plat, &
                                    ps_tmp  , readvar, grid_map='global')
       if(.not. readvar) call endrun()

    end if
    call process_inidat('PS')
    print*,'test in inidat 1',iam,begchunk,endchunk,ideal_phys,aqua_planet
!------------ zhh debug 2013-02-21 ----------
    do j = beglat, endlat
       if (j==10) then
          print*, 'u3(10,26,10,1) =', u3(10,26,10,1), ' v3(10,26,10,1) =', v3(10,26,10,1)
          print*, 't3(10,26,10,1) =', t3(10,26,10,1), ' ps(10,10,1) =', ps(10,10,1)
       end if
    end do
!------------ zhh debug 2013-02-21 ----------

    fieldname = 'SGH'
    if (ideal_phys .or. aqua_planet) then
       sgh(:,:) = 0.
    else
       call infld(fieldname, ncid_topo, 'lon', 'lat', 1, pcols, begchunk, endchunk, &
                                                      sgh, readvar, grid_map='phys')
       if(.not. readvar) call endrun()
    end if

    fieldname = 'SGH30'
    if (ideal_phys .or. aqua_planet) then
       sgh30(:,:) = 0.
    else
       call infld(fieldname, ncid_topo, 'lon', 'lat', 1, pcols, begchunk, endchunk, &
                                                      sgh30, readvar, grid_map='phys')
       if(.not. readvar) then
          if (masterproc) write (6,*) trim(subname), ' Warning: SGH30 not found on topo dataset.'
          if (masterproc) write (6,*) '        The field SGH30 will be filled using data from SGH.'
          do i = begchunk,endchunk
             ncol = get_ncols_p(i)
             sgh30(:ncol,i) = sgh(:ncol,i)
          end do
       end if
    end if

    fieldname = 'LANDM_COSLAT'
    if (aqua_planet) then
       landm(:,:) = 0.
    else
       call infld(fieldname, ncid_topo, 'lon', 'lat', 1, pcols, begchunk, endchunk, &
                                                    landm, readvar, grid_map='phys')
       if (.not. readvar) &
          call endrun(trim(subname)//' ERROR: LANDM_COSLAT not found on topo dataset.')
    end if

#if ( ! defined COUP_CSM )

    fieldname = 'TSICE'
    call infld(fieldname, ncid_ini, 'lon', 'lat', 1, pcols, begchunk, endchunk, &
                                              tsice   , readvar, grid_map='phys')
    if(.not. readvar) call endrun()

    fieldname = 'SNOWHICE'
    call infld(fieldname, ncid_ini, 'lon', 'lat', 1, pcols, begchunk, endchunk, &
                                              snowhice, readvar, grid_map='phys')
    if(.not. readvar) call endrun()
!
! read TS1, TS2, TS3, TS4 in "plevmx" loop
!
    allocate ( arr2d(1:pcols,begchunk:endchunk) )
    do m = 1,plevmx
       fieldname = tsnam(m)
       call infld(fieldname, ncid_ini, 'lon', 'lat', 1, pcols, begchunk, endchunk, &
                                                    arr2d, readvar, grid_map='phys')
       if(.not. readvar) call endrun()
       do i = begchunk,endchunk
          surface_state2d(i)%tssub(:,m) = arr2d(:,i)
       end do
    end do
    deallocate ( arr2d )

    fieldname = 'LANDFRAC'
    if (aqua_planet) then
       landfrac(:,:) = 0.
    else
       call infld(fieldname, ncid_topo, 'lon', 'lat', 1, pcols, begchunk, endchunk, &
                                                 landfrac, readvar, grid_map='phys')
       if(.not. readvar) then
          fieldname = 'FLAND'
          call infld(fieldname, ncid_ini, 'lon', 'lat', 1, pcols, begchunk, endchunk, &
                                                    landfrac, readvar, grid_map='phys')
          if(.not. readvar) call endrun()
       end if
    end if

#if ( defined COUP_SOM )

    fieldname = 'SICTHK'
    call infld(fieldname, ncid_ini, 'lon', 'lat', 1, pcols, begchunk, endchunk, &
                                              sicthk  , readvar, grid_map='phys')
    if(.not. readvar) call endrun()

    fieldname = 'ICEFRAC'
    call infld(fieldname, ncid_ini, 'lon', 'lat', 1, pcols, begchunk, endchunk, &
                                              icefrac , readvar, grid_map='phys')
    if(.not. readvar) call endrun()

    fieldname = 'TSOCN'
    call infld(fieldname, ncid_ini, 'lon', 'lat', 1, pcols, begchunk, endchunk, &
                                              tsocn   , readvar, grid_map='phys')
    if(.not. readvar) call endrun()

! define an initial ocean fraction and non-land ice fraction
! The 1st "where" stmt used to be done in update_srf_fractions (dev45)
!jjr change pcols to ncol
    do c = begchunk,endchunk
       ncol = get_ncols_p(c) !jjr
       where (icefrac(:pcols,c) + landfrac(:pcols,c) > 1.0)
          icefrac(:pcols,c) = 1. - landfrac(:pcols,c)
       end where

       where (landfrac(:pcols,c) < 1.)
          aice(:pcols,c) = icefrac(:pcols,c)/(1. - landfrac(:pcols,c))
       elsewhere
          aice(:pcols,c) = 0.
       end where
       ocnfrac(:pcols,c) = 1. - landfrac(:pcols,c) - icefrac(:pcols,c)
    end do

    write(6,*)'INIDAT: ocnfrac=',ocnfrac(1,begchunk)
!
! Master needs global landfrac
!
    call gather_chunk_to_field(1,1,1,plon,landfrac,landfrac_field)
!   write(6,*)'INIDAT iam=',iam,' landfrac=',landfrac
!   write(6,*)'INIDAT iam=',iam,' landfrac_field=',landfrac_field
!
!JR Could read in Focn from initial dataset if available
!
    Focn(:,:) = 0.
#else
    Focn(:,:) = inf
    frzmlt(:,:) = 0.  ! needs to be 0, otherwise test in tstm always true
    tsocn(:,:) = inf
#endif
#endif
!
! Integrals of mass, moisture and geopotential height
! (fix mass of moisture as well)
!
    call global_int

    deallocate ( ps_tmp   )
    deallocate ( phis_tmp )
!
! --------------------------------------------------------------------
! Read optional fields (if not found, will be arbitrarily initialized)
! --------------------------------------------------------------------
!

!
! 2-D fields
!
    fieldname = 'PBLH'
    call infld(fieldname, ncid_ini, 'lon', 'lat', 1, pcols, begchunk, endchunk, &
                                                 pblht, readvar, grid_map='phys')
    if(.not. readvar) then
       pblht(:,:) = 0.
       if (masterproc) write(6,*) trim(fieldname), ' initialized to 0.'
    end if

    fieldname = 'TPERT'
    call infld(fieldname, ncid_ini, 'lon', 'lat', 1, pcols, begchunk, endchunk, &
                                                 tpert, readvar, grid_map='phys')
    if(.not. readvar) then
       tpert(:,:) = 0.
       if (masterproc) write(6,*) trim(fieldname), ' initialized to 0.'
    end if

    fieldname = 'QPERT'
    qpert(:,:,:) = 0.
    allocate ( arr2d(1:pcols,begchunk:endchunk) )
    call infld(fieldname, ncid_ini, 'lon', 'lat', 1, pcols, begchunk, endchunk, &
                                                 arr2d, readvar, grid_map='phys')
    if(.not. readvar) then
       if (masterproc) write(6,*) trim(fieldname), ' initialized to 0.'
       arr2d(:,:) = 0.
    end if
    do i = begchunk,endchunk
       qpert(:,1,:) = arr2d(:,:)
    end do
    deallocate ( arr2d )

    fieldname = 'TBOT'
    allocate ( arr2d(1:pcols,begchunk:endchunk) )
    call infld(fieldname, ncid_ini, 'lon', 'lat', 1, pcols, begchunk, endchunk, &
                                                 arr2d, readvar, grid_map='phys')
    if(.not. readvar) then
       allocate ( arr2d_tmp(1:plon,1:plat) )
       if (masterproc) write(6,*) trim(fieldname), ' initialized with lowest level of T'
       if (masterproc) arr2d_tmp(:plon,:) = t3_tmp(:plon,plev,:)
       call scatter_field_to_chunk(1,1,1,plon,arr2d_tmp,arr2d)
       deallocate ( arr2d_tmp )
    end if
    do i = begchunk,endchunk
       surface_state2d(i)%tbot(:) = arr2d(:,i)
    end do
    deallocate ( arr2d )

    fieldname = 'TSICERAD'
    call infld(fieldname, ncid_ini, 'lon', 'lat', 1, pcols, begchunk, endchunk, &
                                             tsice_rad, readvar, grid_map='phys')
    if(.not. readvar) then
       tsice_rad(:,:) = tsice(:,:)
       if (masterproc) write(6,*) trim(fieldname), ' initialized with TSICE.'
    end if
!
! 3-D fields
!    
    fieldname = 'CLOUD'
    m = pbuf_get_fld_idx('CLD')
    cldptr => pbuf(m)%fld_ptr(1,1:pcols,1:pver,begchunk:endchunk,1:pbuf_times)
    call infld(fieldname, ncid_ini, 'lon', 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
                                                       cldptr(:,:,:,1), readvar, grid_map='phys')

    if(.not. readvar) then
       cldptr(:,:,:,1) = 0.
       if (masterproc) write(6,*) trim(fieldname), ' initialized to 0.'
    end if
    do n = 2, pbuf_times
       cldptr(:,:,:,n) = cldptr(:,:,:,1)
    end do

    fieldname = 'QCWAT'
    m = pbuf_get_fld_idx(fieldname)
    cldptr => pbuf(m)%fld_ptr(1,1:pcols,1:pver,begchunk:endchunk,1:pbuf_times)
    call infld(fieldname, ncid_ini, 'lon', 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
                                                       cldptr(:,:,:,1), readvar, grid_map='phys')
    if(.not. readvar) then
       allocate ( arr3d_a(1:plon,1:plev,1:plat) )
       if (masterproc) write(6,*) trim(fieldname), ' initialized with Q'
       if (masterproc) arr3d_a(:plon,:,:) = q3_tmp(:plon,:,:)
       call scatter_field_to_chunk(1,plev,1,plon,arr3d_a,cldptr(:,:,:,1))
       deallocate ( arr3d_a )
    end if
    do n = 2, pbuf_times
       cldptr(:,:,:,n) = cldptr(:,:,:,1)
    end do

    fieldname = 'LCWAT'
    m = pbuf_get_fld_idx(fieldname)
    cldptr => pbuf(m)%fld_ptr(1,1:pcols,1:pver,begchunk:endchunk,1:pbuf_times)
    call infld(fieldname, ncid_ini, 'lon', 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
                                                       cldptr(:,:,:,1), readvar, grid_map='phys')
    if(.not. readvar) then
       allocate ( arr3d_a(1:plon,1:plev,1:plat) )
       if (masterproc) write(6,*) trim(fieldname), ' initialized with CLDICE + CLDLIQ'
       if (masterproc) arr3d_a(:plon,:,:) = cld_ice(:plon,:,:) + cld_liq(:plon,:,:)
       call scatter_field_to_chunk(1,plev,1,plon,arr3d_a,cldptr(:,:,:,1))
       deallocate ( arr3d_a )
    end if
    do n = 2, pbuf_times
       cldptr(:,:,:,n) = cldptr(:,:,:,1)
    end do

    fieldname = 'TCWAT'
    m = pbuf_get_fld_idx(fieldname)
    cldptr => pbuf(m)%fld_ptr(1,1:pcols,1:pver,begchunk:endchunk,1:pbuf_times)
    call infld(fieldname, ncid_ini, 'lon', 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
                                                       cldptr(:,:,:,1), readvar, grid_map='phys')
    if(.not. readvar) then
       allocate ( arr3d_a(1:plon,1:plev,1:plat) )
       if (masterproc) write(6,*) trim(fieldname), ' initialized with T'
       if (masterproc) arr3d_a(:plon,:,:) = t3_tmp(:plon,:,:)
       call scatter_field_to_chunk(1,plev,1,plon,arr3d_a,cldptr(:,:,:,1))
       deallocate ( arr3d_a )
    end if
    do n = 2, pbuf_times
       cldptr(:,:,:,n) = cldptr(:,:,:,1)
    end do

    deallocate ( q3_tmp  )
    deallocate ( cld_liq )
    deallocate ( cld_ice )
    deallocate ( t3_tmp  )

#if ( defined SCAM )
      STATUS = NF_CLOSE( NCID_ini )
#endif

    return

  end subroutine read_inidat

#if (defined SCAM )
!*********************************************************************C
   subroutine read_inidat_scam(datatype)
!-----------------------------------------------------------------------
!
! Purpose:
! Read initial level information for single column version of cam
!
!-----------------------------------------------------------------------
!
! $Id: inidat.F90,v 1.20.4.26 2005/03/08 17:06:41 eaton Exp $
! $Author: eaton $
!
!-----------------------------------------------------------------------
   use time_manager, only: get_curr_date

   character*(MAX_PATH_LEN) initdata
   real(r8) time( MAX_RECORDS )
   integer, intent(in) :: datatype       !  in:  model or analysis dataset
   integer max_recs
   parameter(max_recs=MAX_RECORDS)
   integer i, delta
   integer date( MAX_RECORDS ), datesec( MAX_RECORDS )
   integer time_dimID, time_varID,tmpid
   integer lev_dimID
   integer date_varID, datesec_varID
   integer ntime
   integer curr_date
   integer :: ncsec                             ! current time of day [seconds]
   integer :: yr, mon, day                      ! year, month, and day component
   integer iexcept, ieee_handler, myhandler
   logical use_nf_real     ! used for 4 byte mode 

!
!     Trap ieee exceptions on SUN for debugging purposes
!
#if ( defined sun )
   iexcept = ieee_handler( 'set', 'common', myhandler )
   if ( iexcept .ne. 0 ) write( 6,* ) 'ieee trapping not supported here'
#endif

#if USE_4BYTE_REAL
   use_nf_real = .true.
#else
   use_nf_real = .false.
#endif
   have_surfdat = .false.
!
!     Open initial dataset
!
      if ( datatype .eq. MODEL ) then
        initdata = modelfile
        fill_ends = .TRUE.
      endif

      if ( datatype .eq. ANAL ) then
        initdata = analysisfile
        fill_ends = .FALSE.
      endif

      STATUS = NF_OPEN( initdata, NF_NOWRITE, ncid_ini )
      if( STATUS .NE. NF_NOERR )  then
        write( 6,* )'ERROR - readinitdata.F:','Cant open model data netcdf file', &
             initdata
        STATUS = NF_CLOSE( NCID_ini )
        return
      endif

      ! CAM now reads the topography fields from a separate dataset.  The SCAM GUI 
      ! will need to be updated to read that file.  For now the only option is
      ! backwards compatibility, i.e., the topo fields must be on the initial 
      ! dataset.
      ncid_topo = ncid_ini

!     -----------------------------------------------------------------
!     Get dimension IDs and sizes
!     -----------------------------------------------------------------

   STATUS =  NF_INQ_DIMID ( NCID_INI, 'time', time_dimID )
   if ( STATUS .NE. NF_NOERR ) then
      write( 6,* )'ERROR - readinitdata.F:', 'Cant get variable dim for time'
      STATUS = NF_CLOSE( NCID_INI )
      return
   endif

   STATUS = NF_INQ_DIMLEN( NCID_INI, time_dimID, ntime )

   if ( ntime .GT. MAX_RECORDS ) then
      write( *, * ) 'ERROR: readinitdata.F:',initdata, 'has too many records', &
         'It should contain at most', max_recs ,' records'
      STATUS = NF_CLOSE( NCID_INI )
      return
   endif

   STATUS =  NF_INQ_DIMID ( NCID_INI, 'lev', lev_dimID )
   if ( STATUS .NE. NF_NOERR ) then
      write( 6,* )'ERROR - readinitdata.F:', 'Cant get variable dim for lev'
      STATUS = NF_CLOSE( NCID_INI )
      return
   endif

!     -----------------------------------------------------------------
!     Get  time variables. Variable called 'time' is current day of year
!     -----------------------------------------------------------------

   STATUS = NF_INQ_VARID( NCID_INI, 'time', time_varID )
   if ( STATUS .NE. NF_NOERR ) then
      write( 6,* )'ERROR - readinitdata.F:','Cant get variable ID for time'
      STATUS = NF_CLOSE( NCID_INI )
      return
   endif

   if (use_nf_real) then
      STATUS = NF_GET_VAR_REAL( NCID_INI, time_varID, time )
   else
      STATUS = NF_GET_VAR_DOUBLE( NCID_INI, time_varID, time )
   endif

   if ( STATUS .NE. NF_NOERR ) then
      write( 6,* )'ERROR - readinitdata:','Cant get variable time'
      STATUS = NF_CLOSE( NCID_INI )
      return
   endif

   STATUS = NF_INQ_VARID( NCID_INI, 'date', date_varID )
   if ( STATUS .NE. NF_NOERR ) then
      write( 6,* )'ERROR - readinitdata.F:','Cant get variable ID for date'
      STATUS = NF_CLOSE( NCID_INI )
      return
   endif
   STATUS = NF_GET_VAR_INT( NCID_INI, date_varID, date )
   if ( STATUS .NE. NF_NOERR ) then
      write( 6,* )'ERROR - readinitdata:','Cant get variable date'
      STATUS = NF_CLOSE( NCID_INI )
      return
   endif

   STATUS =  NF_INQ_VARID( NCID_INI, 'datesec', datesec_varID )
   if ( STATUS .NE. NF_NOERR )  then
      write( 6,* )'ERROR - readinitdata.F:','Cant get variable ID for datasec'
      STATUS = NF_CLOSE( NCID_INI )
      return
   endif

   STATUS = NF_GET_VAR_INT( NCID_INI, datesec_varID, datesec )
   if ( STATUS .NE. NF_NOERR ) then
      write( 6,* )'ERROR - readinitdata:','Cant get variable datesec'
      STATUS = NF_CLOSE( NCID_INI )
      return
   endif
!
!     Find the index into the date array that corresponds to the date that
!     is closest to mcdate. 
!
   initTimeIdx = 1

   call get_curr_date(yr,mon,day,ncsec)
   curr_date = yr*10000 + mon*100 + day

   delta = abs( curr_date - mod( date( 1 )/100, 100 ))
   do i = 1, ntime
      if ( delta .GT. abs( curr_date - mod( date( i )/100, 100 )))  then
         initTimeIdx = i
         delta = abs( curr_date - mod( date( i )/100, 100 ))
      endif
   end do

  return

   end subroutine read_inidat_scam
#endif

!*********************************************************************
!*********************************************************************

  subroutine process_inidat(fieldname, m_cnst)
!
!-----------------------------------------------------------------------
!
! Purpose:
! Post-process input fields
!
!-----------------------------------------------------------------------
!
! $Id: inidat.F90,v 1.30.4.32 2005/03/08 17:06:48 eaton Exp $
! $Author: eaton $
!
!-----------------------------------------------------------------------
!
    use history,      only: fillvalue
    use constituents, only: cnst_name, qmin
    use chemistry   , only: chem_implements_cnst, chem_init_cnst
    use aerosol_intr, only: aerosol_implements_cnst, aerosol_init_cnst
    use tracers     , only: tracers_implements_cnst, tracers_init_cnst
    use stratiform,      only: stratiform_implements_cnst, stratiform_init_cnst
#if ( defined SPMD )
    use spmd_dyn          , only: comm_y, comm_z
    use parutilitiesmodule, only: parcollective2d, BCSTOP
#endif

#include <comlun.h>
#include <perturb.h>

!
! Input arguments
!
    character(len=*), intent(in)           :: fieldname ! fields to be processed
    integer         , intent(in), optional :: m_cnst    ! constituent index
!
!---------------------------Local workspace-----------------------------
!
    integer  i,j,k,lat                     ! grid and constituent indices
    real(r8) pertval                       ! perturbation value
    integer  varid                         ! netCDF variable id
    integer  ret, attlen                   ! netcdf return values
    logical  phis_hires                    ! true => PHIS came from hi res topo
    real(r8), allocatable :: uv_local (:,:,:)
   character*256 text
    character*80 trunits                   ! tracer untis

    real(r8), pointer, dimension(:,:,:) :: q_tmp

    character*16 :: subname='PROCESS_INIDAT' ! subroutine name

    select case (fieldname)

!----------
! Process U
!----------

    case ('U')

!
! SJL: initialize j=1 because later on arr3d_b will be copied to u3s using f90 array syntax
!
       if(masterproc) then
!
! Transpose array
!

!$omp parallel do private(i, j, k)
          do k = 1, plev
             do j = 1, plat
                do i = 1, plon
                   arr3d_b(i,j,k) = arr3d_a(i,k,j)
                enddo
             enddo
          enddo
       end if
!
! Scatter u3
!

#if ( defined SPMD )
       allocate( uv_local(plon,beglat:endlat,beglev:endlev) )
       call scatter( arr3d_b, strip3dxyz, uv_local, mpicom )

!$omp parallel do private(i,j,k)
       do k = beglev,endlev
          do j = beglat,endlat
             do i = 1,plon
                u3(i,k,j,1) = uv_local(i,j,k) !jjr
             enddo
          enddo
       enddo
       deallocate( uv_local )
#else

!$omp parallel do private(i, j, k)
       do k = beglev,endlev
          do j = beglat,endlat
             do i = 1,plon
                u3(i,k,j,1) = arr3d_b(i,j,k)
             enddo
          enddo
       enddo
#endif

!----------
! Process V
!----------

    case ('V')
       if(masterproc) then
!
! Transpose array
!

!$omp parallel do private(i, j, k)
          do k = 1, plev
             do j = 1, plat
                do i = 1, plon
                   arr3d_b(i,j,k) = arr3d_a(i,k,j)
                enddo
             enddo
          enddo
       end if
!
! Scatter v3
!

#if ( defined SPMD )
       allocate( uv_local(plon,beglat:endlat,beglev:endlev) )
       call scatter( arr3d_b, strip3dxyz, uv_local, mpicom )
!$omp parallel do private(i,j,k)
      do k = beglev,endlev
         do j=beglat,endlat
            do i = 1,plon
               v3(i,k,j,1) = uv_local(i,j,k)
            enddo
         enddo
      enddo
#else
!$omp parallel do private(i, j, k)
      do k = beglev,endlev
         do j = beglat,endlat
            do i = 1,plon
               v3(i,k,j,1) = arr3d_b(i,j,k)
            enddo
         enddo
      enddo
#endif

!----------
! Process T
!----------

    case ('T')

!
! Add random perturbation to temperature if required
!
       if(masterproc) then

          if (pertlim.ne.0.0) then
             write(6,*)trim(subname), ':  Adding random perturbation bounded by +/-', &
                       pertlim,' to initial temperature field'
             do lat = 1,plat
                do k = 1,plev
                   do i = 1,nlon(lat)
                      call random_number (pertval)
                      pertval = 2.*pertlim*(0.5 - pertval)
                      t3_tmp(i,k,lat) = t3_tmp(i,k,lat)*(1. + pertval)
                   end do
                end do
             end do
          end if
!
! Average T at the poles.
!

                                                                !$omp parallel do private(k)
!jjr          do k = 1, plev
!jjr             call xpavg(t3_tmp(1,k,   1), plon)
!jjr             call xpavg(t3_tmp(1,k,plat), plon)
!jjr          enddo

       end if   ! end of if-masterproc

! Scatter t3

#if ( defined SPMD )
      call scatter( t3_tmp, strip3dxzy, t3(1,beglev,beglat,1), mpicom )
#else
!$omp parallel do private(i, j, k)
      do j = beglat,endlat
         do k = beglev,endlev
            do i = 1,plon
               t3(i,k,j,1) = t3_tmp(i,k,j)
            enddo
         enddo
      enddo
#endif

!---------------------
! Process Constituents
!---------------------

    case ('CONSTS')

       if (.not. present(m_cnst)) then
          call endrun('  '//trim(subname)//' Error:  m_cnst needs to be present in the'// &
                      ' argument list')
       end if

       if(masterproc) then
!
! Check that all tracer units are in mass mixing ratios
!
          if(readvar) then
             call wrap_inq_varid    (NCID_INI,cnst_name(m_cnst), varid)
             call wrap_get_att_text (NCID_INI,varid,'units',trunits)
             if (trunits(1:5) .ne. 'KG/KG' .and. trunits(1:5) .ne. 'kg/kg') then
                call endrun('  '//trim(subname)//' Error:  Units for tracer ' &
                            //trim(cnst_name(m_cnst))//' must be in KG/KG')
             end if
                                                                                                                              !
! Constituents not read from initial file are initialized by the package that implements them.
!
          else
             if(m_cnst == 1) then
                call endrun('  '//trim(subname)//' Error:  Q must be on Initial File')
             end if

             allocate ( q_tmp(plon,plev,plat) )
             q_tmp(:,:,:) = 0.
             write(6,*) 'Warning:  Not reading ',cnst_name(m_cnst), ' from IC file.'
             if (stratiform_implements_cnst(cnst_name(m_cnst))) then
                call stratiform_init_cnst(cnst_name(m_cnst), q_tmp)
                write(6,*) '          ', cnst_name(m_cnst), ' initialized by "stratiform_init_cnst"'
             else if (chem_implements_cnst(cnst_name(m_cnst))) then
                call chem_init_cnst(cnst_name(m_cnst), q_tmp)
                write(6,*) '          ', cnst_name(m_cnst), ' initialized by "chem_init_cnst"'
             else if (tracers_implements_cnst(cnst_name(m_cnst))) then
                call tracers_init_cnst(cnst_name(m_cnst), q_tmp)
                write(6,*) '          ', cnst_name(m_cnst), ' initialized by "tracers_init_cnst"'
             else if (aerosol_implements_cnst(cnst_name(m_cnst))) then
                call aerosol_init_cnst(cnst_name(m_cnst), q_tmp)
                write(6,*) '          ', cnst_name(m_cnst), ' initialized by "aerosol_init_cnst"'
             else
                write(6,*) '          ', cnst_name(m_cnst), ' set to 0.'
             end if

             arr3d_a(:plon,:,:) = q_tmp(:,:,:)

             deallocate ( q_tmp )
          end if

!$omp parallel do private(lat)
          do lat = 1,plat
             call qneg3(trim(subname), lat   ,nlon(lat),plon   ,plev    , &
                  m_cnst, m_cnst, qmin(m_cnst) ,arr3d_a(1,1,lat))
          end do

!$omp parallel do private(k)
!jjr          do k = 1, plev
!jjr             call xpavg(arr3d_a(:,k,   1), plon)
!jjr             call xpavg(arr3d_a(:,k,plat), plon)
!jjr          enddo
!
! if "Q", "CLDLIQ", or "CLDICE", save off for later use
                                                                                                  !
          if(m_cnst == 1       ) q3_tmp (:plon,:,:) = arr3d_a(:plon,:,:)
          if(m_cnst == ixcldliq) cld_liq(:plon,:,:) = arr3d_a(:plon,:,:)
          if(m_cnst == ixcldice) cld_ice(:plon,:,:) = arr3d_a(:plon,:,:)

!$omp parallel do private(j,k)
          do k = 1, plev
             do j = 1, plat
                arr3d_b(:,j,k) = arr3d_a(:,k,j)
             enddo
          enddo

       end if   ! end of if-masterproc

       call scatter_q_field_to_block(arr3d_b, m_cnst) !jjr

!-----------
! Process PS
!-----------

    case ('PS')

!
! Average PS at the poles.
!
       if(masterproc) then
!jjr          call xpavg(ps_tmp(1,   1), plon)
!jjr          call xpavg(ps_tmp(1,plat), plon)
       end if

#if ( defined SPMD )
       if (myid_z .eq. 0) then
          call scatter( ps_tmp, strip2d, ps(1,beglat,1), comm_y )!jjr
       endif
       if (twod_decomp .eq. 1) then
          call parcollective2d( comm_z, BCSTOP, plon, endlat-beglat+1, ps(1,beglat,1) )
       endif
#else
       ps(:plon,:,1) = ps_tmp(:plon,:)
#endif

!-------------
! Process PHIS
!-------------
    case ('PHIS')

       if(masterproc) then
!
! Check for presence of 'from_hires' attribute to decide whether to filter
!
          if(readvar) then
             call wrap_inq_varid (ncid_topo, 'PHIS', varid)
             ret = nf_inq_attlen (ncid_topo, varid, 'from_hires', attlen)
             if (ret.eq.NF_NOERR .and. attlen.gt.256) then
                call endrun('  '//trim(subname)//' Error:  from_hires attribute length is too long')
             end if
             ret = nf_get_att_text (ncid_topo, varid, 'from_hires', text)
             if (ret.eq.NF_NOERR .and. text(1:4).eq.'true') then
                phis_hires = .true.
!                write(6,*) trim(subname), ': Will filter input PHIS: attribute from_hires is true'
             else
                phis_hires = .false.
!                write(6,*)trim(subname), ': Will not filter input PHIS: attribute ', &
!                          'from_hires is either false or not present'
             end if
!
! Average PHIS at the poles.
!
!jjr             call xpavg(phis_tmp(1,   1), plon)
!jjr             call xpavg(phis_tmp(1,plat), plon)

          end if

       end if    ! end of if-masterproc

#if ( defined SPMD )
      if (myid_z .eq. 0) then
         call scatter( phis_tmp, strip2d, phis, comm_y )
      endif
      if (twod_decomp .eq. 1) then
         call parcollective2d( comm_z, BCSTOP, plon, endlat-beglat+1, phis )
      endif
#else
      phis(:,:) = phis_tmp(:,:)
#endif

    end select
    return

  end subroutine process_inidat

                                                                             
  subroutine process_inidat1(fieldname, m_cnst)
!
!-----------------------------------------------------------------------
!
! Purpose:
! Post-process input fields
!
! Revised by: ZhangHe, 2007.5.7, 
!-----------------------------------------------------------------------
!
! $Id: inidat.F90,v 1.20.4.26 2005/03/08 17:06:41 eaton Exp $
! $Author: eaton $
!
!-----------------------------------------------------------------------
!
    use commap
    use comspe
!!    use spetru
    use constituents, only: cnst_name, qmin
    use chemistry   , only: chem_implements_cnst, chem_init_cnst
    use aerosol_intr, only: aerosol_implements_cnst, aerosol_init_cnst
    use tracers     , only: tracers_implements_cnst, tracers_init_cnst
    use stratiform,      only: stratiform_implements_cnst, stratiform_init_cnst
#if ( defined SPMD )
    use spmd_dyn, only: npes, compute_gsfactors, comm_y, comm_z
#endif

#include <perturb.h>

!
! Input arguments
!
    character(len=*), intent(in)           :: fieldname ! fields to be processed
    integer         , intent(in), optional :: m_cnst    ! constituent index
!
!---------------------------Local workspace-----------------------------
!
    integer i,j,k,n,lat,irow               ! grid and constituent indices
    real(r8) pertval                       ! perturbation value
    integer  varid                         ! netCDF variable id
    integer  ret, attlen                   ! netcdf return values
    logical  phis_hires                    ! true => PHIS came from hi res topo
    character*256 text
    character*80 trunits                   ! tracer untis

    real(r8), pointer, dimension(:,:,:) :: q_tmp
    real(r8), pointer, dimension(:,:,:) :: tmp3d_a, tmp3d_b, tmp3d_extend
    real(r8), pointer, dimension(:,:  ) :: tmp2d_a, tmp2d_b

!================ delete by zhh ==================
!!#if ( defined BFB_CAM_SCAM_IOP )
!!    real(r8), allocatable :: ps_sav(:,:)
!!    real(r8), allocatable :: u3_sav(:,:,:)
!!    real(r8), allocatable :: v3_sav(:,:,:)
!!#endif
!=================== 2007.5.7 ====================

#if ( defined SPMD )
    integer :: numperlat                   ! number of values per latitude band
    integer :: numsend(0:npes-1)           ! number of items to be sent
    integer :: numrecv                     ! number of items to be received
    integer :: displs(0:npes-1)            ! displacement array
#endif

    character*16 :: subname='PROCESS_INIDAT' ! subroutine name

!===================================================================
!
    select case (fieldname)

!------------
! Process U/V
!------------

    case ('UV')

!
       allocate ( tmp3d_a(plon,plev,plat) )
       allocate ( tmp3d_b(plon,plev,plat) )

       if(masterproc) then
!================ delete by zhh ==================
!
! Spectral truncation
!
!!#if ( defined SCAM )
!!          tmp3d_a(:,:,:) = 0.
!!          tmp3d_b(:,:,:) = 0.
!!#else
!!#if ( defined BFB_CAM_SCAM_IOP )
!!          allocate ( u3_sav (plon,plev,plat) )             
!!          allocate ( v3_sav (plon,plev,plat) )             
!!          u3_sav(:plon,:plev,:plat) = arr3d_a(:plon,:plev,:plat)
!!          v3_sav(:plon,:plev,:plat) = arr3d_b(:plon,:plev,:plat)
!!          call spetru_uv(u3_sav  ,v3_sav  ,vort=tmp3d_a, div=tmp3d_b)
!!          deallocate ( u3_sav )
!!          deallocate ( v3_sav )
!!#else
!!zhh          call spetru_uv(arr3d_a ,arr3d_b ,vort=tmp3d_a, div=tmp3d_b)
!!#endif
!!#endif
!=================== 2007.5.7 ====================

       end if

#if ( defined SPMD )
       numperlat = plnlv    ! plnlv = plon*plev
       call compute_gsfactors (numperlat, numrecv, numsend, displs)
       call mpiscatterv (arr3d_a ,numsend, displs, mpir8,u3  (1,1,beglat,1)  ,numrecv, mpir8,0,mpicom)
       call mpiscatterv (arr3d_b ,numsend, displs, mpir8,v3  (1,1,beglat,1)  ,numrecv, mpir8,0,mpicom)
!
!!       call mpiscatterv (tmp3d_a ,numsend, displs, mpir8,vort(1,1,beglat,1) ,numrecv, mpir8,0,mpicom)
!!       call mpiscatterv (tmp3d_b ,numsend, displs, mpir8,div (1,1,beglat,1) ,numrecv, mpir8,0,mpicom)
#else
       u3    (:,:,:,1) = arr3d_a(:plon,:plev,:plat)
       v3    (:,:,:,1) = arr3d_b(:plon,:plev,:plat)
!
!!       vort  (:,:,:,1) = tmp3d_a(:,:,:)
!!       div   (:,:,:,1) = tmp3d_b(:,:,:)
#endif
!
       deallocate ( tmp3d_a )
       deallocate ( tmp3d_b )

!----------
! Process T
!----------

    case ('T')
!
! Add random perturbation to temperature if required
!
       if(masterproc) then

          if (pertlim.ne.0.0) then
             write(6,*)trim(subname), ':  Adding random perturbation bounded by +/-', &
                       pertlim,' to initial temperature field'
             do lat = 1,plat
                do k = 1,plev
                   do i = 1,nlon(lat)
                      call random_number (pertval)
                      pertval = 2.*pertlim*(0.5 - pertval)
                      t3_tmp(i,k,lat) = t3_tmp(i,k,lat)*(1. + pertval)
                   end do
                end do
             end do
          end if
!
! Spectral truncation
!

!================ delete by zhh ==================
!!#if ( ! defined BFB_CAM_SCAM_IOP ) && ( ! defined SCAM )
!!          call spetru_3d_scalar(t3_tmp)
!!#endif
!=================== 2007.5.7 ====================

       end if   ! end of if-masterproc

#if ( defined SPMD )
       numperlat = plnlv
       call compute_gsfactors (numperlat, numrecv, numsend, displs)
       call mpiscatterv (t3_tmp  ,numsend, displs, mpir8,t3(1,1,beglat,1) ,numrecv, mpir8,0,mpicom)
#else
       t3    (:,:,:,1) = t3_tmp(:plon,:plev,:plat)
#endif

!---------------------
! Process Constituents
!---------------------

    case ('CONSTS')

!
       if (.not. present(m_cnst)) then
          call endrun('  '//trim(subname)//' Error:  m_cnst needs to be present in the'// &
                      ' argument list')
       end if

       allocate ( tmp3d_extend(plon,plev,beglat:endlat) )

       if(masterproc) then
!
! Check that all tracer units are in mass mixing ratios
!
          if(readvar) then
             call wrap_inq_varid    (NCID_INI,cnst_name(m_cnst), varid)
             call wrap_get_att_text (NCID_INI,varid,'units',trunits)
             if (trunits(1:5) .ne. 'KG/KG' .and. trunits(1:5) .ne. 'kg/kg') then
                call endrun('  '//trim(subname)//' Error:  Units for tracer ' &
                            //trim(cnst_name(m_cnst))//' must be in KG/KG')
             end if
!
! Constituents not read from initial file are initialized by the package that implements them.
!
          else
             if(m_cnst == 1) then
                call endrun('  '//trim(subname)//' Error:  Q must be on Initial File')
             end if

             allocate ( q_tmp(plon,plev,plat) )
             q_tmp(:,:,:) = 0.
             write(6,*) 'Warning:  Not reading ',cnst_name(m_cnst), ' from IC file.'
             if (stratiform_implements_cnst(cnst_name(m_cnst))) then
                call stratiform_init_cnst(cnst_name(m_cnst), q_tmp)
                write(6,*) '          ', cnst_name(m_cnst), ' initialized by "stratiform_init_cnst"'
             else if (chem_implements_cnst(cnst_name(m_cnst))) then
                call chem_init_cnst(cnst_name(m_cnst), q_tmp)
                write(6,*) '          ', cnst_name(m_cnst), ' initialized by "chem_init_cnst"'
             else if (tracers_implements_cnst(cnst_name(m_cnst))) then
                call tracers_init_cnst(cnst_name(m_cnst), q_tmp)
                write(6,*) '          ', cnst_name(m_cnst), ' initialized by "tracers_init_cnst"'
             else if (aerosol_implements_cnst(cnst_name(m_cnst))) then
                call aerosol_init_cnst(cnst_name(m_cnst), q_tmp)
                write(6,*) '          ', cnst_name(m_cnst), ' initialized by "aerosol_init_cnst"'
             else
                write(6,*) '          ', cnst_name(m_cnst), ' set to 0.'
             end if

             arr3d_a(:plon,:,:) = q_tmp(:,:,:)

             deallocate ( q_tmp )
          end if

!wjp 2011.04.16 !$OMP parallel do private(lat)
          do lat = 1,plat
             call qneg3(trim(subname), lat   ,nlon(lat),plon   ,plev    , &
                  m_cnst, m_cnst, qmin(m_cnst) ,arr3d_a(1,1,lat))
          end do
!
! if "Q", "CLDLIQ", or "CLDICE", save off for later use
!
          if(m_cnst == 1       ) q3_tmp (:plon,:,:) = arr3d_a(:plon,:,:)
          if(m_cnst == ixcldliq) cld_liq(:plon,:,:) = arr3d_a(:plon,:,:)
          if(m_cnst == ixcldice) cld_ice(:plon,:,:) = arr3d_a(:plon,:,:)

       end if   ! end of if-masterproc

#if ( defined SPMD )
       numperlat = plnlv
       call compute_gsfactors (numperlat, numrecv, numsend, displs)
       call mpiscatterv (arr3d_a  , numsend, displs, mpir8, tmp3d_extend(1,1,beglat) ,numrecv, mpir8,0,mpicom)
       q3(:,:,m_cnst,:,1) = tmp3d_extend(:,:,beglat:endlat)
#else
       q3(:,:plev,m_cnst,:,1) = arr3d_a(:plon,:plev,:plat)
#endif
       deallocate ( tmp3d_extend )

!-----------
! Process PS
!-----------

    case ('PS')

!
!!       allocate ( tmp2d_a(plon,plat) )
!!       allocate ( tmp2d_b(plon,plat) )

       if(masterproc) then
!================ delete by zhh ==================
!
! Spectral truncation
!

!!#if ( defined SCAM )
!!          tmp2d_a(:,:) = 0.
!!          tmp2d_b(:,:) = 0.
!!#else
!!#if ( defined BFB_CAM_SCAM_IOP )
!!          allocate ( ps_sav(plon,plat) )                   
!!          ps_sav(:plon,:plat)=ps_tmp(:plon,:plat)
!!          call spetru_ps(ps_sav, tmp2d_a, tmp2d_b)
!!          deallocate ( ps_sav )
!!#else
!!          call spetru_ps(ps_tmp, tmp2d_a, tmp2d_b)
!!#endif
!!#endif
!=================== 2007.5.7 ====================

       end if
#if ( defined SPMD )
       numperlat = plon
       call compute_gsfactors (numperlat, numrecv, numsend, displs)
       call mpiscatterv (ps_tmp,  numsend, displs, mpir8,ps(1,beglat,1), numrecv, mpir8,0,mpicom)
!!       call mpiscatterv (tmp2d_a ,numsend, displs, mpir8,dpsl(1,beglat) ,numrecv, mpir8,0,mpicom)
!!       call mpiscatterv (tmp2d_b ,numsend, displs, mpir8,dpsm(1,beglat) ,numrecv, mpir8,0,mpicom)
#else
!!       dpsl(:,:) = tmp2d_a(:,:)
!!       dpsm(:,:) = tmp2d_b(:,:)
         ps(:plon,:,1) = ps_tmp(:plon,:)
#endif

!!       deallocate ( tmp2d_a )
!!       deallocate ( tmp2d_b )

!-------------
! Process PHIS
!-------------

    case ('PHIS')

!
       if(masterproc) then
!
! Check for presence of 'from_hires' attribute to decide whether to filter
!
          if(readvar) then
             call wrap_inq_varid (ncid_topo, 'PHIS', varid)
             ret = nf_inq_attlen (ncid_topo, varid, 'from_hires', attlen)
             if (ret.eq.NF_NOERR .and. attlen.gt.256) then
                call endrun('  '//trim(subname)//' Error:  from_hires attribute length is too long')
             end if
             ret = nf_get_att_text (ncid_topo, varid, 'from_hires', text)
             if (ret.eq.NF_NOERR .and. text(1:4).eq.'true') then
                phis_hires = .true.
                write(6,*) trim(subname), ': Will filter input PHIS: attribute from_hires is true'
             else
                phis_hires = .false.
                write(6,*)trim(subname), ': Will not filter input PHIS: attribute ', &
                          'from_hires is either false or not present'
             end if

          end if
!================ delete by zhh ==================
!
! Spectral truncation
!
!!#if ( ! defined BFB_CAM_SCAM_IOP ) && ( ! defined SCAM )
!!          call spetru_phis  (phis_tmp, phis_hires)
!!#endif
!=================== 2007.5.7 ====================

       end if   ! end of if-masterproc
          
#if ( defined SPMD )
       numperlat = plon
       call compute_gsfactors (numperlat, numrecv, numsend, displs)
       call mpiscatterv (phis_tmp  ,numsend, displs, mpir8,phis  (1,beglat) ,numrecv, mpir8,0,mpicom)
#else
!wjp 2011.04.16 !$OMP parallel do private(lat)
       do lat = 1,plat
          phis(:nlon(lat),lat) = phis_tmp(:nlon(lat),lat)
       end do
#endif

    end select

!
    return

  end subroutine process_inidat1

!*********************************************************************

  subroutine global_int
!
!-----------------------------------------------------------------------
!
! Purpose:
! Compute global integrals of mass, moisture and geopotential height
! and fix mass of atmosphere
!
!-----------------------------------------------------------------------
!
! $Id: inidat.F90,v 1.20.4.26 2005/03/08 17:06:41 eaton Exp $
! $Author: eaton $
!
!-----------------------------------------------------------------------
!
    use commap       ! w
    use physconst,    only: gravit
#if ( defined SPMD )
    use parutilitiesmodule, only: parcollective2d, BCSTOP
    use mpishorthand
    use spmd_dyn, only: npes, compute_gsfactors,comm_y,comm_z
#endif

#include <comctl.h>
#include <comhyb.h>
#include <comqfl.h>

!
!---------------------------Local workspace-----------------------------
!
    integer i,k,lat,ihem,irow  ! grid indices
    real(r8) pdelb(plon,plev)  ! pressure diff between interfaces
                               ! using "B" part of hybrid grid only
    real(r8) pssum             ! surface pressure sum
    real(r8) dotproda          ! dot product
    real(r8) dotprodb          ! dot product
    real(r8) zgssum            ! partial sums of phis
    real(r8) hyad (plev)       ! del (A)
    real(r8) tmassf_tmp        ! Global mass integral
    real(r8) qmass1_tmp        ! Partial Global moisture mass integral
    real(r8) qmass2_tmp        ! Partial Global moisture mass integral
    real(r8) qmassf_tmp        ! Global moisture mass integral
    real(r8) zgsint_tmp        ! Geopotential integral

#if ( defined SPMD )
    integer :: numperlat         ! number of values per latitude band
    integer :: numsend(0:npes-1) ! number of items to be sent
    integer :: numrecv           ! number of items to be received
    integer :: displs(0:npes-1)  ! displacement array
#endif
!
!-----------------------------------------------------------------------
!
    if(masterproc) then
!        
! Initialize mass and moisture integrals for summation
! in a third calculation loop (assures bit-for-bit compare
! with non-random history tape).
!
       tmassf_tmp = 0.
       qmass1_tmp = 0.
       qmass2_tmp = 0.
       zgsint_tmp = 0.
!
! Compute pdel from "A" portion of hybrid vertical grid for later use in global integrals
!
       do k = 1,plev
          hyad(k) = hyai(k+1) - hyai(k)
       end do
       do k = 1,plev
          do i = 1,plon
             pdela(i,k) = hyad(k)*ps0
          end do
       end do
!
! Compute integrals of mass, moisture, and geopotential height
!
#if ( defined SCAM )
       do irow = 1,plat
#else
       do irow = 1,plat/2
#endif
          do ihem = 1,2
             if (ihem.eq.1) then
                lat = irow
             else
                lat = plat - irow + 1
             end if
!              
! Accumulate average mass of atmosphere
!
             call pdelb0 (ps_tmp(1,lat),pdelb   ,nlon(lat))
             pssum  = 0.
             do i = 1,nlon(lat)
                pssum  = pssum  + ps_tmp  (i,lat)
             end do
             tmassf_tmp = tmassf_tmp + w(irow)*pssum/nlon(lat)

             zgssum = 0.
             do i = 1,nlon(lat)
                zgssum = zgssum + phis_tmp(i,lat)
             end do
             zgsint_tmp = zgsint_tmp + w(irow)*zgssum/nlon(lat)
!
! Calculate global integrals needed for water vapor adjustment
!
             do k = 1,plev
                dotproda = 0.
                dotprodb = 0.
                do i = 1,nlon(lat)
                   dotproda = dotproda + q3_tmp(i,k,lat)*pdela(i,k)
                   dotprodb = dotprodb + q3_tmp(i,k,lat)*pdelb(i,k)
                end do
                qmass1_tmp = qmass1_tmp + w(irow)*dotproda/nlon(lat)
                qmass2_tmp = qmass2_tmp + w(irow)*dotprodb/nlon(lat)
             end do
          end do
       end do                  ! end of latitude loop
!
! Normalize average mass, height
!
       tmassf_tmp = tmassf_tmp*.5/gravit
       qmass1_tmp = qmass1_tmp*.5/gravit
       qmass2_tmp = qmass2_tmp*.5/gravit
       zgsint_tmp = zgsint_tmp*.5/gravit
       qmassf_tmp = qmass1_tmp + qmass2_tmp
!
! Globally avgd sfc. partial pressure of dry air (i.e. global dry mass):
!
       tmass0 = 98222./gravit
       if (ideal_phys ) tmass0 =  100000./gravit
       if (aqua_planet) tmass0 = (101325._r8-245._r8)/gravit
       write(6,800) tmassf_tmp,tmass0,qmassf_tmp
       write(6,810) zgsint_tmp
800    format(/72('*')//'INIDAT: Mass of initial data before correction = ' &
              ,1p,e20.10,/,' Dry mass will be held at = ',e20.10,/, &
              ' Mass of moisture after removal of negatives = ',e20.10) 
810    format('INIDAT: Globally averaged geopotential ', &
              'height = ',f16.10,' meters'//72('*')/)

#if (! defined SCAM )
!
! Compute and apply an initial mass fix factor which preserves horizontal
! gradients of ln(ps).
!
       if (adiabatic .or. ideal_phys) then
          fixmas = tmass0/tmassf_tmp
       else
          fixmas = (tmass0 + qmass1_tmp)/(tmassf_tmp - qmass2_tmp)
       end if
       do lat = 1,plat
          do i = 1,nlon(lat)
             ps_tmp(i,lat) = ps_tmp(i,lat)*fixmas
          end do
       end do
#endif 
!
! Global integerals
!
       tmassf = tmassf_tmp
       qmass1 = qmass1_tmp
       qmass2 = qmass2_tmp
       qmassf = qmassf_tmp
       zgsint = zgsint_tmp

    end if   ! end of if-masterproc

#if ( defined SPMD )
    call mpibcast (tmass0,1,mpir8,0,mpicom)
    call mpibcast (tmassf,1,mpir8,0,mpicom)
    call mpibcast (qmass1,1,mpir8,0,mpicom)
    call mpibcast (qmass2,1,mpir8,0,mpicom)
    call mpibcast (qmassf,1,mpir8,0,mpicom)
    call mpibcast (zgsint,1,mpir8,0,mpicom)
#endif

#if ( defined SPMD )
       if (myid_z .eq. 0) then
          call scatter( ps_tmp, strip2d, ps(1,beglat,1), comm_y )!jjr
       endif
       if (twod_decomp .eq. 1) then
          call parcollective2d( comm_z, BCSTOP, plon, endlat-beglat+1, ps(1,beglat,1) )
       endif

!jjr
!jjr    numperlat = plon
!jjr    call compute_gsfactors (numperlat, numrecv, numsend, displs)
!jjr    call mpiscatterv (ps_tmp    ,numsend, displs, mpir8,ps    (1,beglat,1) ,numrecv, mpir8,0,mpicom)
#else
!wjp 2011.04.16 !$OMP parallel do private(lat)
    do lat = 1,plat
       ps(:nlon(lat),lat,1) = ps_tmp(:nlon(lat),lat)
    end do
#endif
    return

  end subroutine global_int
!-----------------------------------------------------------------------
!BOP
! !IROUTINE: scatter_q_field_to_block --- scatter a 3D constituent array to prognostic array q3
!
! !INTERFACE:
   subroutine scatter_q_field_to_block(xyz, cnst_idx)

! !USES:
      use shr_kind_mod, only: r8 => shr_kind_r8
#if ( defined SPMD )
      use mpishorthand, only: mpicom
#endif

! !INPUT PARAMETERS:
      real(r8), dimension(plon,plat,plev), intent(in) :: &
         xyz        ! 3D constituent field
      integer, intent(in) :: &
         cnst_idx   ! constituent index in prognostic array q3

! !DESCRIPTION:
!
! Scatter a 3D constituent array from the master processor to the
! q3 array in the prognostics module.
!
! !REVISION HISTORY:
!
!   02.07.31   Eaton      Initial version
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer j, k
      real(r8), allocatable :: xyz_local(:,:,:)

!-----------------------------------------------------------------------

#if ( defined SPMD )

      allocate( xyz_local(plon,beglat:endlat,beglev:endlev) )
      call scatter( xyz, strip3dxyz, xyz_local, mpicom )
      do k = beglev,endlev
         do j = beglat,endlat
            q3(:,k,cnst_idx,j,1) = xyz_local(:,j,k)
         enddo
      enddo
      deallocate( xyz_local )

#else

      do j = beglat,endlat
         do k = beglev,endlev
            q3(:,k,cnst_idx,j,1) = xyz(:,j,k)
         enddo
      enddo

#endif

      return

!EOC
   end subroutine scatter_q_field_to_block


end module inidat
