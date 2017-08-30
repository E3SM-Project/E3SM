module physpkg
  !-----------------------------------------------------------------------
  ! Purpose:
  !
  ! Provides the interface to CAM physics package
  !
  ! Revision history:
  ! Aug  2005,  E. B. Kluzek,  Creation of module from physpkg subroutine
  ! 2005-10-17  B. Eaton       Add contents of inti.F90 to phys_init().  Add
  !                            initialization of grid info in phys_state.
  ! Nov 2010    A. Gettelman   Put micro/macro physics into separate routines
  ! July 2015   B. Singh       Added code for unified convective transport
  !-----------------------------------------------------------------------

  use shr_kind_mod,     only: r8 => shr_kind_r8, i8 => shr_kind_i8, longchar=>SHR_KIND_CL
  use physconst,        only: latvap, latice, rh2o
  use physics_types,    only: physics_state, physics_tend,&
       physics_ptend, physics_tend_init,&
       physics_type_alloc, physics_ptend_dealloc,&
       physics_state_alloc, physics_state_dealloc, physics_tend_alloc, physics_tend_dealloc
  use ppgrid,           only: begchunk, endchunk, pcols, pver, pverp, psubcols
  use constituents,     only: pcnst, cnst_name, cnst_get_ind
  use glb_verif_smry,   only: tp_stat_smry, global_smry_init

  use cam_logfile,     only: iulog
  
  use cam_abortutils,      only : endrun !BSINGH

  implicit none
  private
  
  save

  ! Public methods
  public phys_init   ! Public initialization method
  public phys_final  ! Public finalization method
  !
  ! Private module data
  !
  ! Physics package options
  !======================================================================= 
contains



  !======================================================================= 

!!$subroutine phys_inidat()
!!$    use cam_abortutils, only : endrun
!!$
!!$    use cam_initfiles,       only: initial_file_get_id, topo_file_get_id
!!$    use ncdio_atm,           only: infld
!!$    integer :: lchnk, m, n, i, k, ncol
!!$    type(file_desc_t), pointer :: fh_ini, fh_topo
!!$    character(len=8) :: fieldname
!!$    real(r8), pointer :: cldptr(:,:,:,:), convptr_3d(:,:,:,:)
!!$    real(r8), pointer :: tptr(:,:), tptr3d(:,:,:), tptr3d_2(:,:,:)
!!$    real(r8), pointer :: qpert(:,:)
!!$
!!$    character*11 :: subname='phys_inidat' ! subroutine name
!!$    logical :: found=.false., found2=.false.
!!$    integer :: ierr
!!$    character(len=8) :: dim1name, dim2name
!!$    integer :: ixcldice, ixcldliq
!!$    integer                   :: grid_id  ! grid ID for data mapping
!!$    nullify(tptr,tptr3d,tptr3d_2,cldptr,convptr_3d)
!!$
!!$    fh_ini=>initial_file_get_id()
!!$
!!$    !   dynamics variables are handled in dyn_init - here we read variables needed for physics 
!!$    !   but not dynamics
!!$
!!$    !
!!$    ! 3-D fields
!!$    !
!!$
!!$    allocate(tptr3d(pcols,pver,begchunk:endchunk))
!!$
!!$    fieldname='CLOUD'
!!$    m = pbuf_get_index('CLD')
!!$    call infld(fieldname, fh_ini, dim1name, 'lev', dim2name, 1, pcols, 1, pver, begchunk, endchunk, &
!!$         tptr3d, found, gridname='physgrid')
!!$    if(found) then
!!$       do n = 1, dyn_time_lvls
!!$          call pbuf_set_field(pbuf2d, m, tptr3d, (/1,1,n/),(/pcols,pver,1/))
!!$       end do
!!$    else
!!$       call pbuf_set_field(pbuf2d, m, 0._r8)
!!$       if (masterproc) write(iulog,*) trim(fieldname), ' initialized to 0.'
!!$    end if
!!$
!!$
!!$end subroutine phys_inidat


subroutine phys_init( phys_state, phys_tend, chunk_smry, domain_smry, nstep)

    !----------------------------------------------------------------------- 
    ! 
    ! Initialization of physics package.
    ! 
    !-----------------------------------------------------------------------

    use physconst,          only: rair, cpair, gravit, stebol, tmelt, &
                                  latvap, latice, rh2o, rhoh2o, pstd, zvir, &
                                  karman, rhodair, physconst_init 
    use ascii_io, only : read_state_ascii
    use ppgrid,             only: pver, begchunk, endchunk
!!$    use check_energy,       only: check_energy_timestep_init
!!$    use physics_buffer,     only: physics_buffer_desc, pbuf_initialize, pbuf_get_index
!!$    use ref_pres,           only: pref_edge, pref_mid
!!$    use tracers,            only: tracers_init
!!    use clubb_intr,         only: clubb_ini_cam
!!$
!!$
!!$
    ! Input/output arguments
    type(physics_state), pointer :: phys_state(:)     ! (beginchunk:endchunk)
    type(physics_tend ), pointer :: phys_tend (:)     ! (beginchunk:endchunk)
    type(tp_stat_smry),  pointer :: chunk_smry (:,:)  ! (nstatfld,beginchunk:endchunk)
    type(tp_stat_smry),  pointer :: domain_smry(:)    ! (nstatfld)
    integer, intent(IN) :: nstep

    ! local variables
    integer :: lchnk
    integer :: ichnk

    character(len=longchar) :: exe_name, namelist_filename
    character(len=longchar) :: ic_filepath

    namelist /test_nl/ ic_filepath
    !-----------------------------------------------------------------------

    call physics_type_alloc( phys_state, phys_tend, begchunk, endchunk, pcols )
    call global_smry_init( chunk_smry, domain_smry, begchunk, endchunk )
    call physconst_init()

    ! Parse command line argument to locate namelist file
    call getarg(0, exe_name)
    call getarg(1, namelist_filename)

    ! Read namelist to get path to IC files 
    write(iulog,*) 'Looking for test_nl in ', trim(namelist_filename)

    open(unit=10,file=trim(namelist_filename))
    read(10,nml=test_nl)
    close(10)

    write(iulog,*) 'File path for initial conditions is ' , trim(ic_filepath)
  
    ! Read in the initial conditions 
    do ichnk=begchunk,endchunk
       call read_state_ascii(trim(ic_filepath), nstep, phys_state(ichnk))
    end do

end subroutine phys_init

  !

subroutine phys_final( phys_state, phys_tend )
!!$    use physics_buffer, only : physics_buffer_desc, pbuf_deallocate
!!$    use chemistry, only : chem_final
!!$    use carma_intr, only : carma_final
!!$    use wv_saturation, only : wv_sat_final
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Finalization of physics package
    ! 
    !-----------------------------------------------------------------------
    ! Input/output arguments
    type(physics_state), pointer :: phys_state(:)
    type(physics_tend ), pointer :: phys_tend(:)

    deallocate(phys_state)
    deallocate(phys_tend)
end subroutine phys_final


end module physpkg
