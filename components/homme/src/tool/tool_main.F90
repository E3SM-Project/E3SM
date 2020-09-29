#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

! This program is an offline tool. See the documentation in the following
! subroutines, implemented below; each is a tool:
!   * topo_convert
!   * topo_pgn_to_smoothed
!
! See also homme/test/tool/README
!
program tool_main
  use prim_driver_mod,  only: prim_init1, prim_init2, prim_finalize, deriv1, smooth_topo_datasets
  use hybvcoord_mod,    only: hvcoord_t, hvcoord_init
  use parallel_mod,     only: parallel_t, initmp, syncmp, haltmp, abortmp
  use hybrid_mod,       only: hybrid_t, hybrid_create
  use dimensions_mod,   only: nelemd,ne,np
  use domain_mod,       only: domain1d_t
  use element_mod,      only: element_t
  use common_io_mod,    only: output_dir, infilenames, output_varnames2
  use time_mod,         only: timelevel_t
  use control_mod,      only: vfile_mid, vfile_int, theta_hydrostatic_mode, topology, test_case,&
       vanalytic, theta_hydrostatic_mode, ftype, smooth_phis_numcycle
  use common_io_mod,    only: tool
  use kinds,            only: iulog
  use interpolate_driver_mod, only : pio_read_phis, interpolate_driver
  use native_mapping, only: create_native_mapping_files
  use gll_subcell_grid, only: write_gll_subcell_grid
  use prim_movie_mod,   only : prim_movie_output, prim_movie_finish,prim_movie_init




  implicit none

  type (element_t),  pointer  :: elem(:)
  type (hybrid_t)             :: hybrid         ! parallel structure for shared memory/distributed memory
  type (parallel_t)           :: par            ! parallel structure for distributed memory programming
  type (domain1d_t), pointer  :: dom_mt(:)
  type (hvcoord_t)            :: hvcoord        ! hybrid vertical coordinate struct
  type (TimeLevel_t)          :: tl             ! Main time level struct
  integer :: ithr, nets, nete, ierr
  
  par = initmp()

  call set_namelist_defaults()

  call prim_init1(elem, par, dom_mt, tl)
  theta_hydrostatic_mode=.true.  ! disable some NH tests

  ! Set up fake threading; this offline tool doesn't thread.
  ithr = 0
  hybrid = hybrid_create(par,ithr,1)
  nets = 1
  nete = nelemd

  hvcoord = hvcoord_init(vfile_mid, vfile_int, .true., hybrid%masterthread, ierr)
  if (ierr /= 0) call haltmp("error in hvcoord_init")

  call prim_init2(elem, hybrid, nets, nete, tl, hvcoord)

  select case(tool)
     case('topo_convert')
        call topo_convert(par, elem)
     case ('topo_pgn_to_smoothed')
        call topo_pgn_to_smoothed(par, elem)
     case ('topo_gll_to_smoothed')
        call topo_gll_to_smoothed(elem, hybrid, tl)
     case ('interpolate_tool')
        ! Interpolate a netcdf file from one grid to another
        call interpolate_tool(elem, hybrid)
     case ('gll_mapping_file')
        ! output HOMME's internal interpolation operator as a mapping file
        call gll_mapping_file(elem, hybrid, tl)
     case ('grid_template_tool')
        call grid_template_tool(elem, hybrid, tl)
     case ('gll_subcell_grid')
        call write_gll_subcell_grid(par, elem)
     case('none')
        if (par%masterproc) then
           write(iulog,*) 'homme_tool was given tool="none"; exiting without doing anything'
        end if
     case default
        if (par%masterproc) then
           write(iulog,*) 'homme_tool was given tool=', trim(tool), ', which is not a tool'
        end if
     end select

  call syncmp(hybrid%par)  ! wait for I/O tasks to finish     
  call prim_finalize()

  call haltmp("exiting homme_tool")

contains

  subroutine set_namelist_defaults()
    ! Set some values that are irrelevant to the tool but make the namelist
    ! processor succeed.

    use time_mod, only: tstep
    use dimensions_mod, only: qsize_d, qsize

    tstep = 1
    topology = 'cube'
    test_case = 'dcmip2012_test1_1'
    vanalytic = 1
    qsize = qsize_d
  end subroutine set_namelist_defaults

  subroutine topo_convert(par, elem)
    ! Convert a pure-GLL topography file to GLL-physgrid topography file. This
    ! is useful if you have an existing v1 GLL topo file and want to run an
    ! equivalent physgrid-based simulation without having to run the topography
    ! tool chain from scratch.
    !
    ! Namelist example:
    !
    ! &ctl_nl
    ! ne = 30
    ! /
    ! &vert_nl
    ! /
    ! &analysis_nl
    ! tool = 'topo_convert'
    ! infilenames = 'USGS-gtopo30_ne30np4_16xdel2-PFC-consistentSGH.nc', 'USGS-gtopo30_ne30np4pg2_16xdel2-PFC-consistentSGH_converted'
    ! /

    use gllfvremap_util_mod, only: gfr_convert_topo

    type (parallel_t), intent(in) :: par
    type (element_t), intent(inout) :: elem(:)

    if (min(len(trim(infilenames(1))), len(trim(infilenames(2)))) == 0) then
       call abortmp('homme_tool: topo_convert requires infilenames 1 and 2 to be defined')
    end if
    call gfr_convert_topo(par, elem, 2, infilenames(1), infilenames(2))
  end subroutine topo_convert

  subroutine topo_pgn_to_smoothed(par, elem)
    ! Use a pure-physgrid non-smoothed topography file to create a smoothed
    ! GLL-physgrid topography file. This is part of the physgrid topography full
    ! tool chain.
    !   The input is an unsmoothed pg4 topography file output by
    ! cube_to_target. The output is a GLL-physgrid file with smoothed topography
    ! and GLL-physgrid consistent PHIS_d and PHIS fields, respectively. This
    ! file is then input for a second run of cube_to_target, this time to
    ! compute SGH, SGH30, LANDFRAC, and LANDM_COSLAT on physgrid.
    !
    !
    ! Namelist example:
    ! NOTE:  smooth_phis_numcycle, smooth_phis_nudt and hypervis_order settings
    ! are grid dependent.  See test/tool/toposmooth_gll.nl for recommended settings
    ! for other grids   
    !
    ! &ctl_nl
    ! ne = 30
    ! smooth_phis_numcycle = 16
    ! smooth_phis_nudt = 28e7
    ! hypervis_scaling = 0
    ! hypervis_order = 2
    ! se_ftype = 2 ! actually output NPHYS; overloaded use of ftype
    ! /
    ! &vert_nl
    ! /
    ! &analysis_nl
    ! tool = 'topo_pgn_to_smoothed'
    ! infilenames = 'ne30pg4_c2t_topo.nc', 'USGS-gtopo30_ne30np4pg2_16xdel2'
    ! /

    use gllfvremap_util_mod, only: gfr_pgn_to_smoothed_topo

    type (parallel_t), intent(in) :: par
    type (element_t), intent(inout) :: elem(:)

    integer :: output_nphys, stat

    if (min(len(trim(infilenames(1))), len(trim(infilenames(2)))) == 0) then
       call abortmp('homme_tool: topo_pgn_to_smoothed requires infilenames 1 and 2 to be defined')
    end if
    output_nphys = ftype
    stat = gfr_pgn_to_smoothed_topo(par, elem, output_nphys, infilenames(1), infilenames(2))
    if (stat /= 0) then
       call abortmp('homme_tool: gfr_pgn_to_smoothed_topo returned an error code')
    end if
  end subroutine topo_pgn_to_smoothed


  subroutine topo_gll_to_smoothed(elem,hybrid,tl)
    use hybrid_mod, only : hybrid_t
    use viscosity_mod, only : smooth_phis

    implicit none
    ! Read, smooth and output a GLL topo field
    ! input specified in infilenames(1)
    ! output: phis_smoothed1.nc
    !
    ! To convert the output into a standalone topo file:
    ! % ncks -O -v geos,lat,lon phis_smoothed1.nc phis_smoothed2.nc
    ! % ncrename -v geos,PHIS  phis_smoothed2.nc
    !
    ! Namelist example: see test/tool/toposmooth_gll.nl
    !
    type (hybrid_t)      , intent(in) :: hybrid
    type (element_t)     , intent(inout), target :: elem(:)
    type (TimeLevel_t)   , intent(in)     :: tl     

    character*8 :: varname

    if (smooth_phis_numcycle>0) then
       if (len(trim(infilenames(1)))==0 ) then
          call abortmp('homme_tool: topo_gll_to_smoothed requires infilenames 1 to be defined')
       end if
       varname='PHIS'
       call pio_read_phis(elem,hybrid%par,'PHIS')
       call smooth_topo_datasets(elem, hybrid, 1, nelemd)
       test_case = 'phis-smoothed'
    else
       ! in this case, we just output PHIS from the initial condition
       ! used to generate test data
       test_case = 'phis-baroclinic'
    endif

    ! call output
    call prim_movie_init( elem, par, hvcoord, tl )
    call prim_movie_output(elem, tl, hvcoord, par)
    call prim_movie_finish
  end subroutine
    

  
  subroutine interpolate_tool(elem,hybrid)
    use hybrid_mod, only : hybrid_t

    implicit none
    ! Interpolate a native grid E3SM atmosphere history file to a lat/lon history file
    ! notes:
    !   vertical coordinates have to match those used in E3SM
    !   homme_tool has to be compiled with PLEV=72  (matching E3SM file)
    !   infilesnames(1) = name of E3SM history file to interpolate
    !   output will use the sane name but with an interp1.nc extension
    !
    ! Namelist example: see test/tool/interpolate.nl

    type (hybrid_t)      , intent(in) :: hybrid
    type (element_t)     , intent(inout), target :: elem(:)


    ! debug
    call interpolate_driver(elem, hybrid)

  end subroutine


  subroutine gll_mapping_file(elem,hybrid,tl)
    use hybrid_mod, only : hybrid_t

    implicit none
    ! create mapping files from HOMME's internal interpolation algorithm
    ! notes:
    !   mapping file output:  map_ne30np4_to_dstfile.nc
    !   The SE coordinates and areas are written to "coords1.nc" and need
    !   to be converted to radians and the copied into xc_a, yc_a, and area_a
    !
    ! Namelist example: see test/tool/mappingfiles.nl

    type (hybrid_t)      , intent(in) :: hybrid
    type (element_t)     , intent(inout), target :: elem(:)
    type (TimeLevel_t)   , intent(in)     :: tl     

    vanalytic=0  ! read vcoord files given in namelist
    call create_native_mapping_files( hybrid, elem,'bilin')
    !call create_native_mapping_files( hybrid, elem,'native')

    test_case = 'coords'
    call prim_movie_init( elem, hybrid%par, hvcoord, tl )
    call prim_movie_output(elem, tl, hvcoord, par)
    call prim_movie_finish

  end subroutine
    
  subroutine grid_template_tool(elem,hybrid,tl)
    use hybrid_mod, only : hybrid_t
    use viscosity_mod, only : smooth_phis

    implicit none
    ! Read, smooth and output a GLL topo field
    ! input specified in infilenames(1)
    ! output: phis_smoothed1.nc
    !
    ! To convert the output into a standalone topo file:
    ! % ncks -O -v geos,lat,lon phis_smoothed1.nc phis_smoothed2.nc
    ! % ncrename -v geos,PHIS  phis_smoothed2.nc
    !
    ! Namelist example: see test/tool/toposmooth_gll.nl
    !
    type (hybrid_t)      , intent(in) :: hybrid
    type (element_t)     , intent(inout), target :: elem(:)
    type (TimeLevel_t)   , intent(in)     :: tl     
    character*80 :: nestr,npstr

    ! call output
    write(nestr,*) ne
    write(npstr,*) np
    test_case = 'ne' // trim(adjustl(nestr)) // 'np' // trim(adjustl(npstr)) // '_tmp'
    call prim_movie_init( elem, par, hvcoord, tl )
    call prim_movie_output(elem, tl, hvcoord, par)
    call prim_movie_finish
  end subroutine
    

  

  
end program tool_main
