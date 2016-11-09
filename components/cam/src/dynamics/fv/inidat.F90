module inidat
!----------------------------------------------------------------------- 
! 
! Purpose: Read initial dataset and process fields as appropriate
!
! Method: Read and process one field at a time
! 
! Author: J. Olson  May 2004
!
! !HISTORY:
!   2004.05.01   Olson      Creation
!   2005.07.12   Sawyer     Revised to use grid information
!   2005.11.10   Sawyer     Now using dyn_import/export_t containers
!   2006.04.13   Sawyer     Removed dependency on prognostics
!   2006.06.28   Sawyer     Changed T3 from IKJ to IJK, q3 to tracer
!   2009.04.21   Edwards    Changed to parallel io
!   2010 Nov     A. Gettelman and C. Craig  put micro/macro physics into separate routines
!
!
!-----------------------------------------------------------------------
   use shr_kind_mod,       only: r8 => shr_kind_r8
   use spmd_utils,         only: masterproc
   use pmgrid,             only: plev
   use ncdio_atm,          only: infld
   use dyn_internal_state, only: get_dyn_state_grid
   use dynamics_vars,      only: T_FVDYCORE_GRID
   use dyn_comp,           only: dyn_import_t
   use cam_abortutils,         only: endrun
   use phys_grid,          only: get_ncols_p
   use cam_control_mod,    only: ideal_phys, aqua_planet, moist_physics
   use fv_control_mod,     only: tmass0
   use cam_logfile,        only: iulog
   use pio, only : file_desc_t, io_desc_t, pio_double, pio_freedecomp, pio_setdebuglevel, &
                   pio_noerr, pio_inq_varid, pio_inq_attlen, pio_get_att



#if ( defined SPMD )
   use mpishorthand
#endif

   implicit none

PRIVATE
   include 'netcdf.inc'

   real(r8), allocatable :: ps_tmp  (:,:  )
   real(r8), allocatable :: phis_tmp(:,:  )

   logical readvar            ! inquiry flag:  true => variable exists on netCDF file

   real(r8), parameter ::  D0_0                    =  0.0_r8
   real(r8), parameter ::  D0_5                    =  0.5_r8
   real(r8), parameter ::  D1_0                    =  1.0_r8
   real(r8), parameter ::  D2_0                    =  2.0_r8
   real(r8), parameter ::  D1E5                    =  1.0e5_r8
   real(r8), parameter ::  DRY_AIR_SLP             =  98222._r8 

   public :: read_inidat

contains



  subroutine global_int(grid, dyn_in)
!
!-----------------------------------------------------------------------
!
! Purpose:
! Compute global integral of geopotential height
!
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
!-----------------------------------------------------------------------
!
    use commap,       only: w
    use physconst,    only: gravit
    use fv_control_mod, only: zgsint    
    use spmd_utils, only : iam
#if ( defined SPMD )
    use spmd_dyn,       only: mpicom_xy, npes_xy
    use mod_comm,       only: mp_sendirr, mp_recvirr
#endif
    type (T_FVDYCORE_GRID), intent(in) :: grid
    type (dyn_import_t), target, intent(inout) :: dyn_in
    real (r8), pointer :: phisxy(:,:)
    real (r8), allocatable :: phisglob(:,:)

!
!---------------------------Local workspace-----------------------------
!
    integer i,lat              ! grid indices
    real(r8) zgssum            ! partial sums of phis
    real(r8) zgsint_tmp        ! Geopotential integral
    
    integer :: im              ! From grid for convenience
    integer :: jm              ! From grid for convenience
!
!-----------------------------------------------------------------------
!
    phisxy => dyn_in%phis
    im = grid%im
    jm = grid%jm
    allocate (phisglob(im,jm))
#if ( defined SPMD )
    if (iam .lt. npes_xy) then
       call mp_sendirr( mpicom_xy, grid%g_2dxy_r8%SendDesc,                          &
                        grid%g_2dxy_r8%RecvDesc, phisxy, phisglob,                   &
                        modc=grid%modc_gather )
       call mp_recvirr( mpicom_xy, grid%g_2dxy_r8%SendDesc,                          &
                        grid%g_2dxy_r8%RecvDesc, phisxy, phisglob,                   &
                        modc=grid%modc_gather )
    endif
#else
    phisglob(:,:) = phisxy(:,:)
#endif

    if(masterproc) then

       zgsint_tmp = D0_0
!              
! Accumulate average geopotential
!
       do lat = 1, jm
          zgssum = D0_0
          do i = 1, im
             zgssum = zgssum + phisglob(i,lat)
          end do
          zgsint_tmp = zgsint_tmp + w(lat)*zgssum/im
       end do
!
! Normalize average height
!
       zgsint_tmp = zgsint_tmp*D0_5/gravit
       zgsint     = zgsint_tmp
!
! Globally avgd sfc. partial pressure of dry air (i.e. global dry mass):
! (appears to have no relevance to FV dycore but left in for now)
!
       tmass0 = DRY_AIR_SLP/gravit
       if (ideal_phys) tmass0 = D1E5/gravit
       write(iulog,800) zgsint
800    format(/72('*')//'INIDAT:  Globally averaged geopotential height = ',f16.10,' meters'//72('*')/)

    end if    ! end of if-masterproc

    deallocate (phisglob)

    return

  end subroutine global_int


   subroutine read_inidat( ncid_ini, ncid_topo, dyn_in )
!
!-----------------------------------------------------------------------
!
! Purpose:
! Read initial dataset and spectrally truncate as appropriate.
!
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
!-----------------------------------------------------------------------
!
    use constituents,     only: cnst_name, cnst_read_iv, cnst_get_ind
    use cam_history_support,      only: fillvalue
!
! Arguments
!
   type(file_desc_t),        intent(inout)    :: ncid_ini
   type(file_desc_t),        intent(inout)    :: ncid_topo
   type (dyn_import_t), target, intent(inout) :: dyn_in   ! dynamics import
!
!---------------------------Local workspace-----------------------------
!
    integer i,m,n                           ! indices
    integer ncol
    integer :: im, jm, km, ntotq            ! From grid for convenience
    integer :: ierror                       ! For allocation errors

    type (T_FVDYCORE_GRID), pointer       :: grid
    character*16 fieldname                  ! field name
    integer :: ifirstxy, ilastxy, jfirstxy, jlastxy, k
     
    character*16 :: subname='READ_INIDAT'   ! subroutine name

!
!-----------------------------------------------------------------------
!     May 2004 revision described below (Olson)
!-----------------------------------------------------------------------
!
! This routine reads and processes fields one at a time to minimize 
! memory usage.
!
!   State fields (including PHIS) are read in on the
!     appropriate grid and processed
!
!
!-----------------------------------------------------------------------
!
    
    grid   => get_dyn_state_grid()
    im     =  grid%im
    jm     =  grid%jm
    km     =  grid%km
    ntotq  =  grid%ntotq
    ifirstxy =  grid%ifirstxy
    ilastxy  =  grid%ilastxy
    jfirstxy =  grid%jfirstxy
    jlastxy  =  grid%jlastxy


!---------------------
! Read required fields
!---------------------
	     
!
!-----------
! 2-D fields
!-----------
!

    fieldname = 'PS'
    call infld(fieldname, ncid_ini, 'lon', 'lat', ifirstxy, ilastxy, jfirstxy, jlastxy, &
         dyn_in%ps  , readvar, gridname='fv_centers')
    if(.not. readvar) call endrun(trim(subname)//' ERROR: PS not found on initial dataset.')

    call process_inidat(ncid_ini, ncid_topo, grid, dyn_in, 'PS')

    fieldname = 'PHIS'
    readvar   = .false.
    if (ideal_phys .or. aqua_planet) then
       dyn_in%phis(:,:) = D0_0
    else
       allocate(phis_tmp(im, jm))
       call infld(fieldname, ncid_topo, 'lon', 'lat', ifirstxy, ilastxy, jfirstxy, jlastxy, &
            dyn_in%phis, readvar, gridname='fv_centers')
       if (.not. readvar) &
            call endrun(trim(subname)//' ERROR: PHIS not found on topo dataset.')
       call process_inidat(ncid_ini, ncid_topo, grid, dyn_in, 'PHIS')
    end if


!
!-----------
! 3-D fields
!-----------
!

    fieldname = 'US'
    dyn_in%u3s = fillvalue

    call infld(fieldname, ncid_ini, 'lon', 'slat', 'lev',  ifirstxy, ilastxy, jfirstxy, jlastxy, 1, km, &
        dyn_in%u3s, readvar, gridname='fv_u_stagger')

    if(.not. readvar) call endrun()

    fieldname = 'VS'
    call infld(fieldname, ncid_ini, 'slon', 'lat', 'lev',  ifirstxy, ilastxy, jfirstxy, jlastxy, 1, km, &
         dyn_in%v3s, readvar, gridname='fv_v_stagger')
    if(.not. readvar) call endrun()

    fieldname = 'T'

    call infld(fieldname, ncid_ini, 'lon', 'lat', 'lev', ifirstxy, ilastxy, jfirstxy, jlastxy, 1, km, &
       dyn_in%t3, readvar, gridname='fv_centers')
    call process_inidat(ncid_ini, ncid_topo, grid, dyn_in, 'T')
        

    ! Constituents (read and process one at a time)

    do m = 1,ntotq
       readvar   = .false.
       fieldname = cnst_name(m)
       if(cnst_read_iv(m)) then
          call infld(fieldname, ncid_ini, 'lon', 'lat', 'lev', ifirstxy, ilastxy, jfirstxy, jlastxy, 1, km, &
               dyn_in%tracer(:,:,:,m), readvar, gridname='fv_centers')
       end if
	
       call process_inidat(ncid_ini, ncid_topo, grid, dyn_in, 'CONSTS', m_cnst=m)

    end do
!
! Global integral of geopotential height (diagnostic for now)
!
    if (.NOT. (ideal_phys .or. aqua_planet)) then
       call global_int(grid, dyn_in)
    end if

    return

  end subroutine read_inidat




!*********************************************************************

  subroutine process_inidat(ncid_ini, ncid_topo, grid, dyn_in, fieldname, m_cnst)
!
!-----------------------------------------------------------------------
!
! Purpose:
! Post-process input fields
!
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
!-----------------------------------------------------------------------
!
    use constituents, only: cnst_name, qmin
    use chemistry   , only: chem_implements_cnst, chem_init_cnst
    use carma_intr,   only: carma_implements_cnst, carma_init_cnst
    use tracers     , only: tracers_implements_cnst, tracers_init_cnst
    use aoa_tracers , only: aoa_tracers_implements_cnst, aoa_tracers_init_cnst
    use clubb_intr,   only: clubb_implements_cnst, clubb_init_cnst
    use stratiform,   only: stratiform_implements_cnst, stratiform_init_cnst
    use microp_driver, only: microp_driver_implements_cnst, microp_driver_init_cnst
    use phys_control,  only: phys_getopts
    use co2_cycle   , only: co2_implements_cnst, co2_init_cnst
    use unicon_cam,   only: unicon_implements_cnst, unicon_init_cnst
#if ( defined SPMD )
    use mod_comm,           only: mp_sendirr, mp_recvirr
    use spmd_dyn, only: npes_xy, mpicom_xy
#endif
    use spmd_utils, only: iam
    use cam_control_mod, only: pertlim
    use dyn_grid, only : get_block_gcol_d, get_horiz_grid_dim_d
!
! Input arguments
!
    type(file_desc_t),           intent(inout)    :: ncid_ini
    type(file_desc_t),           intent(inout)    :: ncid_topo
    type (T_FVDYCORE_GRID), target, intent(inout) :: grid   ! dynamics state grid
    type (dyn_import_t),     target, intent(inout) :: dyn_in      ! dynamics import
    character(len=*),                intent(in)    :: fieldname   ! fields to be processed
    integer,               optional, intent(in)    :: m_cnst      ! constituent index

!
!---------------------------Local workspace-----------------------------
!
    integer :: i,j,k                     ! grid and constituent indices
    integer :: nglon, nglat, rndm_seed_sz
    integer, allocatable :: rndm_seed(:)
    real(r8) :: pertval                       ! perturbation value
    integer ::   varid                         ! netCDF variable id
    integer  ret, attlen                   ! netcdf return values
    real(r8), allocatable :: uv_local (:,:,:)
    character*256 text
    character*256 trunits                  ! tracer untis

    real(r8), pointer :: phisxy(:,:), psxy(:,:), t3xy(:,:,:)
    real(r8), pointer :: u3sxy(:,:,:), v3sxy(:,:,:)

    integer :: im, jm, km, jfirst, jlast
    integer :: ifirstxy, ilastxy, jfirstxy, jlastxy

    integer :: ierror                       ! For allocation errors
    real(r8) :: xsum(grid%km)               ! temp array for parallel sums

    character*16 :: subname='PROCESS_INIDAT' ! subroutine name
    real(r8) :: polesum(1)
    real(r8), pointer :: tracer(:,:,:,:), q3tmp(:,:)
    integer :: lpes_xy
    integer, allocatable :: gcid(:)
    integer :: blksiz, ib


#if ( defined SPMD )
    lpes_xy=npes_xy
#else
    lpes_xy=iam+1
#endif

    u3sxy  => dyn_in%u3s
    v3sxy  => dyn_in%v3s

    psxy   => dyn_in%ps
    phisxy => dyn_in%phis
    t3xy   => dyn_in%t3

    im     =  grid%im
    jm     =  grid%jm
    km     =  grid%km
    jfirst =  grid%jfirst
    jlast  =  grid%jlast

    ifirstxy =  grid%ifirstxy
    ilastxy  =  grid%ilastxy
    jfirstxy =  grid%jfirstxy
    jlastxy  =  grid%jlastxy


    select case (fieldname)

!----------
! Process U
!----------

    case ('U')


!----------
! Process V
!----------

    case ('V')

!----------
! Process T
!----------

    case ('T')

       if(iam>=lpes_xy) return

       if (pertlim .ne. D0_0) then

          ! Add random perturbation to temperature if required

          if(masterproc) then
             write(iulog,*)trim(subname), ':  Adding random perturbation bounded by +/-', &
                  pertlim,' to initial temperature field'
          end if

          call get_horiz_grid_dim_d(nglon, nglat)
          call random_seed(size=rndm_seed_sz)
          allocate(rndm_seed(rndm_seed_sz))

          do j = jfirstxy, jlastxy
             do i = ifirstxy, ilastxy
                ! seed random_number generator based on global column index
                rndm_seed = i + (j-1)*nglon
                call random_seed(put=rndm_seed)
                do k = 1, km
                   call random_number(pertval)
                   pertval = D2_0*pertlim*(D0_5 - pertval)
                   t3xy(i,j,k) = t3xy(i,j,k)*(D1_0 + pertval)
                end do
             end do
          end do
             
          deallocate(rndm_seed)
          
       end if
!
! Average T at the poles.
!
      if ( jfirstxy == 1 ) then
         call par_xsum(grid, t3xy(:,1,:), km, xsum )
         do k=1, km
            do i=ifirstxy, ilastxy
               t3xy(i,1,k) = xsum(k) / real(im,r8)
            enddo
         enddo
      endif
      if ( jlastxy == jm ) then
         call par_xsum(grid, t3xy(:,jm,:), km, xsum )
         do k=1, km
            do i=ifirstxy, ilastxy
               t3xy(i,jm,k) = xsum(k) / real(im,r8)
            enddo
         enddo
      endif


!---------------------
! Process Constituents
!---------------------

    case ('CONSTS')

       if (.not. present(m_cnst)) then
          call endrun('  '//trim(subname)//' Error:  m_cnst needs to be present in the'// &
                      ' argument list')
       end if
       tracer => dyn_in%tracer	

!
! Check that all tracer units are in mass mixing ratios
!
       if(readvar) then
          ret = pio_inq_varid    (NCID_INI,cnst_name(m_cnst), varid)
          ret = pio_get_att (NCID_INI,varid,'units',trunits)
          if (trunits(1:5) .ne. 'KG/KG' .and. trunits(1:5) .ne. 'kg/kg') then
             call endrun('  '//trim(subname)//' Error:  Units for tracer ' &
                  //trim(cnst_name(m_cnst))//' must be in KG/KG')
          end if
!
! Constituents not read from initial file are initialized by the package that implements them.
!
       else
          if(iam>=lpes_xy) return
          if(m_cnst == 1 .and. moist_physics) then
             call endrun('  '//trim(subname)//' Error:  Q must be on Initial File')
          end if
          if(masterproc) write(iulog,*) 'Warning:  Not reading ',cnst_name(m_cnst), ' from IC file.'

          blksiz = (jlastxy-jfirstxy+1)*(ilastxy-ifirstxy+1)
          allocate( gcid(blksiz),stat=ierror )
          allocate( q3tmp(blksiz,plev),stat=ierror )
          q3tmp(:,:) = D0_0
          tracer(:,:,:,m_cnst) = D0_0	
          ib=0

	  call get_block_gcol_d(iam+1,blksiz,gcid)

          if (microp_driver_implements_cnst(cnst_name(m_cnst))) then
             call microp_driver_init_cnst(cnst_name(m_cnst),q3tmp , gcid)
             if(masterproc) write(iulog,*) '          ', cnst_name(m_cnst), ' initialized by "microp_driver_init_cnst"'
          else if (clubb_implements_cnst(cnst_name(m_cnst))) then
             call clubb_init_cnst(cnst_name(m_cnst),q3tmp , gcid)
             if(masterproc) write(iulog,*) '          ', cnst_name(m_cnst), ' initialized by "clubb_init_cnst"'
          else if (stratiform_implements_cnst(cnst_name(m_cnst))) then
             call stratiform_init_cnst(cnst_name(m_cnst),q3tmp , gcid)
             if(masterproc) write(iulog,*) '          ', cnst_name(m_cnst), ' initialized by "stratiform_init_cnst"'
          else if (chem_implements_cnst(cnst_name(m_cnst))) then
             call chem_init_cnst(cnst_name(m_cnst),q3tmp, gcid)
             if(masterproc) write(iulog,*) '          ', cnst_name(m_cnst), ' initialized by "chem_init_cnst"'
          else if (tracers_implements_cnst(cnst_name(m_cnst))) then
             call tracers_init_cnst(cnst_name(m_cnst),q3tmp, gcid)
             if(masterproc) write(iulog,*) '          ', cnst_name(m_cnst), ' initialized by "tracers_init_cnst"'
          else if (aoa_tracers_implements_cnst(cnst_name(m_cnst))) then
             call aoa_tracers_init_cnst(cnst_name(m_cnst),q3tmp, gcid)
             if(masterproc) write(iulog,*) '          ', cnst_name(m_cnst), ' initialized by "aoa_tracers_init_cnst"'
          else if (carma_implements_cnst(cnst_name(m_cnst))) then
             call carma_init_cnst(cnst_name(m_cnst),q3tmp, gcid)
             if(masterproc) write(iulog,*) '          ', cnst_name(m_cnst), ' initialized by "carma_init_cnst"'
          else if (co2_implements_cnst(cnst_name(m_cnst))) then
             call co2_init_cnst(cnst_name(m_cnst),q3tmp, gcid)
             if(masterproc) write(iulog,*) '          ', cnst_name(m_cnst), ' initialized by "co2_init_cnst"'
          else if (unicon_implements_cnst(cnst_name(m_cnst))) then
             call unicon_init_cnst(cnst_name(m_cnst),q3tmp, gcid)
             if(masterproc) write(iulog,*) '          ', cnst_name(m_cnst), ' initialized by "unicon_init_cnst"'
          else
             if(masterproc) write(iulog,*) '          ', cnst_name(m_cnst), ' set to 0.'
          end if

          do k=1,km
             ib=0
             do j=jfirstxy,jlastxy
                do i=ifirstxy,ilastxy                
                   ib=ib+1
                   tracer(i,j,k,m_cnst) = q3tmp(ib,k)
                end do
             end do
          end do
          deallocate(q3tmp)
          deallocate(gcid)
       end if   ! end of if-readvar

       do k=1,km
          do j=jfirstxy,jlastxy
             do i=ifirstxy,ilastxy                
                tracer(i,j,k,m_cnst) = max(tracer(i,j,k,m_cnst), qmin(m_cnst))
             end do
          end do
       end do

      if(iam>=lpes_xy) return
!
! Compute polar average
!
      if ( jfirstxy == 1 ) then
         call par_xsum(grid, dyn_in%tracer(:,1,:,m_cnst), km, xsum )
         do k=1, km
            do i=ifirstxy, ilastxy
               dyn_in%tracer(i,1,k,m_cnst) = xsum(k) / real(im,r8)
            enddo
         enddo
      endif
      if ( jlastxy == jm ) then
         call par_xsum(grid, dyn_in%tracer(:,jm,:,m_cnst), km, xsum )
         do k=1, km
            do i=ifirstxy, ilastxy
               dyn_in%tracer(i,jm,k,m_cnst) = xsum(k) / real(im,r8)
            enddo
         enddo
      endif
      


!-----------
! Process PS
!-----------

    case ('PS')

!
! Average PS at the poles.
!
      if ( jfirstxy == 1 ) then
         if (size(psxy,2) > 0) then
            call par_xsum(grid, psxy(:,1:1), 1, xsum(1:1) )
            do i=ifirstxy, ilastxy
               psxy(i,1) = xsum(1) / real(im,r8)
            end do
         end if
      endif
      if ( jlastxy == jm ) then
         call par_xsum(grid, psxy(:,jm:jm), 1, xsum(1:1) )
         do i=ifirstxy, ilastxy
            psxy(i,jm) = xsum(1) / real(im,r8)
         enddo
      endif



!-------------
! Process PHIS
!-------------

    case ('PHIS')

      ! Average PHIS at the poles.
      if ( jfirstxy == 1 ) then
         if (size(phisxy,2) > 0) then
            call par_xsum(grid, phisxy(:,1:1), 1, xsum(1:1) )
            do i=ifirstxy, ilastxy
               phisxy(i,1) = xsum(1) / real(im,r8)
            end do
         end if
      endif
      if ( jlastxy == jm ) then
         call par_xsum(grid, phisxy(:,jm:jm), 1, xsum(1:1) )
         do i=ifirstxy, ilastxy
            phisxy(i,jm) = xsum(1) / real(im,r8)
         enddo
      endif

    end select

    return

  end subroutine process_inidat

!*********************************************************************


end module inidat

