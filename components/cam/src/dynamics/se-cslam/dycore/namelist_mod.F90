module namelist_mod
  !-----------------
  use cam_logfile,    only: iulog
  !-----------------
  use params_mod,     only: recursive, sfcurve
  !-----------------
  use shr_string_mod, only: shr_string_toUpper
  use shr_kind_mod,   only: r8=>shr_kind_r8
  !-----------------
  use control_mod,    only:   &
       partmethod,            & ! Mesh partitioning method (METIS)
       multilevel,            &
       numnodes,              &
       tasknum,               & ! used dg model in AIX machine
       remapfreq,             & ! number of steps per remapping call
       remap_type,            & ! selected remapping option
       statefreq,             & ! number of steps per printstate call
       runtype,               &
       cubed_sphere_map,      &
       prescribed_wind,       &
       limiter_option,        &
       nu,                    &
       nu_s,                  &
       nu_q,                  &
       nu_div,                &
       nu_top,                &
       hypervis_scaling,      & ! use tensor HV instead of scalar coefficient
       disable_diagnostics,   & ! Use to disable diagnostics for timing reasons
       hypervis_power,        &
       columnpackage,         &
       tracer_transport_type, &
       TRACERTRANSPORT_CONSISTENT_SE_FVM

  !-----------------
  use thread_mod, only : omp_get_max_threads, max_num_threads, horz_num_threads, vert_num_threads, tracer_num_threads
  !-----------------
  use dimensions_mod, only : ne, np, npdg, nnodes, nmpi_per_node, npart, qsize, qsize_d, set_mesh_dimensions
  !-----------------
  !-----------------
  use cam_abortutils, only: endrun
  use parallel_mod,   only: parallel_t, partitionfornodes, useframes
  !-----------------


  use interpolate_mod, only : set_interp_parameter, get_interp_parameter

!=============================================================================!
  implicit none
  private
!
! This module should contain no global data and should only be 'use'd to
!    call one of the public interfaces below
!
  public :: homme_set_defaults
  public :: homme_postprocess_namelist

 contains

  ! ============================================
  ! homme_set_defaults:
  !
  !  Set default values for namelist variables
  !
  ! ============================================
  subroutine homme_set_defaults()
    npart               = 1
    useframes           = 0
    multilevel          = 1
    numnodes            = -1
    runtype             = 0
    statefreq           = 1
    remapfreq           = 240
    remap_type          = "parabolic"
    tasknum             =-1
    columnpackage       = "none"
    nu_top              = 0
    ne                  = 0
    disable_diagnostics = .false.

  end subroutine homme_set_defaults

  subroutine homme_postprocess_namelist(mesh_file, par)
    use mesh_mod,        only: MeshOpen

    ! Dummy arguments
    character(len=*),  intent(in) :: mesh_file
    type (parallel_t), intent(in) :: par

    ! Local variable
    real(kind=r8) :: dt_max
    character(len=*), parameter :: subname = 'HOMME_POSTPROCESS_NAMELIST: '

    if(par%masterproc) then
      write(iulog, *) subname, 'omp_get_max_threads() = ', max_num_threads
    end if

    if((vert_num_threads > 1) .and. (limiter_option .ne. 8)) then
       if(par%masterproc) then
         write(iulog, *) subname, 'WARNING: vertical threading on supported for limiter_option != 8 '
       end if
       vert_num_threads = 1
    endif

    if (ne /= 0) then
      if (mesh_file /= "none" .and. mesh_file /= "/dev/null") then
        if (par%masterproc) then
          write(iulog, *) subname, "mesh_file:", trim(mesh_file),        &
               " and ne:",ne," are both sepcified in the input file."
          write(iulog,*) "            Specify one or the other, but not both."
        end if
        call endrun(subname//"Do not specify ne if using a mesh file input.")
      end if
    end if
    if (par%masterproc) then
      write(iulog,*) subname, "Mesh File:", trim(mesh_file)
    end if
    if (ne == 0) then
      if (par%masterproc) then
        write (iulog,*) subname, "Opening Mesh File:", trim(mesh_file)
      end if
      call set_mesh_dimensions()
      call MeshOpen(mesh_file, par)
    end if

    ! set map
    if (cubed_sphere_map < 0) then
      if (ne == 0) then
        cubed_sphere_map = 2  ! element_local for var-res grids
      else
        cubed_sphere_map = 0  ! default is equi-angle gnomonic
      end if
    end if

    if ((cubed_sphere_map /= 0) .AND.                                         &
        tracer_transport_type .eq. TRACERTRANSPORT_CONSISTENT_SE_FVM) then
      if (par%masterproc) then
        write(iulog, *) subname, 'fvm transport and require equi-angle gnomonic cube sphere mapping.'
        write(iulog, *) '         Set cubed_sphere_map = 0 or comment it out all together.                          '
      end if
      call endrun(subname//"ERROR: fvm transport and cubed_sphere_map>0")
    end if
    if (par%masterproc) then
      write (iulog,*) subname, "Reference element projection: cubed_sphere_map=",cubed_sphere_map
    end if

    !logic around different hyperviscosity options
    if (hypervis_power /= 0) then
      if (hypervis_scaling /= 0) then
        if (par%masterproc) then
          write(iulog, *) subname, 'Both hypervis_power and hypervis_scaling are nonzero.'
          write(iulog, *) '        (1) Set hypervis_power=1, hypervis_scaling=0 for HV based on an element area.'
          write(iulog, *) '        (2) Set hypervis_power=0 and hypervis_scaling=1 for HV based on a tensor.'
          write(iulog, *) '        (3) Set hypervis_power=0 and hypervis_scaling=0 for constant HV.'
        end if
        call endrun(subname//"ERROR: hypervis_power>0 and hypervis_scaling>0")
      end if
    end if

    if((prescribed_wind /= 0) .and. (prescribed_wind /= 1))then
      call endrun(subname//'prescribed_wind should be either 0 or 1')
    end if

    ! some default diffusion coefficiets
    if (nu_s < 0) then
      nu_s = nu
    end if
    if (nu_q < 0) then
      nu_q = nu
    end if
    if (nu_div < 0) then
      nu_div = nu
    end if

    if (multilevel <= 0) then
      nmpi_per_node = 1
    end if

    nnodes = npart / nmpi_per_node

    if((numnodes > 0) .and. (multilevel == 1)) then
      nnodes = numnodes
      nmpi_per_node = npart/nnodes
    end if

    ! ====================================================================
    !  Do not perform node level partitioning if you are only on one node
    ! ====================================================================
    if((nnodes .eq. 1) .and. PartitionForNodes) then
      PartitionForNodes = .FALSE.
    end if

  end subroutine homme_postprocess_namelist
end module namelist_mod
