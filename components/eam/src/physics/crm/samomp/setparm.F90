module setparm_mod
  use task_util_mod
  implicit none

contains

  subroutine setparm

    ! initialize parameters:

    use vars
    !use micro_params
    use params
    use microphysics, only: micro_setparm
    use sgs, only: sgs_setparm

    implicit none

    integer icondavg, ierr

    dt    = CRM_DT
    dx    = CRM_DX
    dy    = CRM_DY
    
    rank  = 0   ! in MMF model, rank = 0
    !------------------------------------
    !  Set parameters


    ! Allow only special cases for separate output:

    output_sep = output_sep.and.RUN3D
    if(output_sep)  save2Dsep = .true.

    if(RUN2D) dy=dx

    if(RUN2D.and.YES3D.eq.1) then
      print*,'Error: 2D run and YES3D is set to 1. Exitting...'
      call task_abort()
    endif
    if(RUN3D.and.YES3D.eq.0) then
      print*,'Error: 3D run and YES3D is set to 0. Exitting...'
      call task_abort()
    endif

    if(ny.eq.1) dy=dx
    dtn = dt

    notopened2D = .true.
    notopened3D = .true.

    !        call zero_instr_diag() ! initialize instruments output
    call sgs_setparm() ! read in SGS options from prm file.
    call micro_setparm() ! read in microphysical options from prm file.

    if(navgmom_x.lt.0.or.navgmom_y.lt.0) then
      nstatmom        = 1
      nstatmomstart    = 99999999
      nstatmomend      = 999999999
    end if

    masterproc = rank.eq.0

  end subroutine setparm
end module setparm_mod
