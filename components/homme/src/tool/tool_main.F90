#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

! This program is an offline tool. See the documentation in the following
! subroutines, implemented below; each is a tool:
!   * topo_convert

program tool_main
  use prim_driver_mod,  only: prim_init1, prim_init2, prim_finalize
  use hybvcoord_mod,    only: hvcoord_t, hvcoord_init
  use parallel_mod,     only: parallel_t, initmp, syncmp, haltmp, abortmp
  use hybrid_mod,       only: hybrid_t, hybrid_create
  use dimensions_mod,   only: nelemd
  use domain_mod,       only: domain1d_t
  use element_mod,      only: element_t
  use common_io_mod,    only: output_dir, infilenames
  use time_mod,         only: timelevel_t
  use control_mod,      only: vfile_mid, vfile_int
  use common_io_mod,    only: tool
  use kinds,            only: iulog

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
     case('none')
        if (par%masterproc) then
           write(iulog,*) 'homme_tool was given tool="none"; exiting without doing anything'
        end if
     case default
        if (par%masterproc) then
           write(iulog,*) 'homme_tool was given tool=', trim(tool), ', which is not a tool'
        end if
     end select

  call prim_finalize()

  call haltmp("exiting homme_tool")

contains

  subroutine set_namelist_defaults()
    ! Set some values that are irrelevant to the tool but make the namelist
    ! processor succeed.

    use time_mod, only: tstep
    use control_mod, only: topology, test_case, vanalytic
    use dimensions_mod, only: qsize_d, qsize

    tstep = 1
    topology = 'cube'
    test_case = 'dcmip2012_test1_1'
    vanalytic = 1
    qsize = qsize_d
  end subroutine set_namelist_defaults

  subroutine topo_convert(par, elem)
    ! Convert pure-GLL topography file to GLL-physgrid topography file. Namelist
    ! example:
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
       call abortmp('homme_tool: gfr_convert_topo requires infilenames 1 and 2 to be defined')
    end if
    call gfr_convert_topo(par, elem, 2, infilenames(1), infilenames(2))
  end subroutine topo_convert
  
end program tool_main
