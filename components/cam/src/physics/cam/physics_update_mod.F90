module physics_update_mod 
  !============================================================================
  !Author: Balwinder Singh
  !
  !Purpose: An interface for physics update call:
  !For perturbation growth test, we need to output certain variables after every
  !physics_update call (when the state variables are updated). Placing outfld
  !directly in the existing physics_update call resulted in hard to resolve
  !circular dependencies. The alternative was to create this interface which calls 
  !physics_update and outfld subroutines one after another.

  !The circular dependencies, if we add outfld call to physics_types module, were:
  ! [NOTE: '<-' can be read as "depends upon" as the example below can be read as 
  !"physics_type depends upon cam_history" (due to the outfld call in physics_type)]

  !1. physics_type<-cam_history<-subcol_utils<-physics_types
  !2. physics_types<-cam_history<-chem_surfvals<-mo_flbc<-phys_gmean<-physics_type

  !This module helps breaking these circular dependencies but it requires changes 
  !to other parts of the code where we need to replace the "use" statement of
  !physics_update call with physics_update_mod module
  !============================================================================

  use spmd_utils,    only: masterproc
  use cam_abortutils,only: endrun
  use cam_history,   only: outfld, fieldname_len  
  use shr_kind_mod,  only: r8 => shr_kind_r8
  use ppgrid,        only: pcols, pver, begchunk
  use physics_types, only: physics_update_main, physics_ptend, physics_state, physics_tend

  implicit none
  private
  public  :: physics_update_init, physics_update, get_var
  
  save

  character(len = 25), parameter :: fname = 'pergro_ptend_names.txt'
  character(len = fieldname_len) :: plist(70)
  logical :: pergro_test_active
  integer :: unitn, pid

  !Following arrays and variables are declared so that we can add all variables in a loop to the history files for pergro test.
  !For adding any new variables, we need to do the following:
  !
  !1. Add variable name to 'hist_vars' (and increment nvars_prtrb_hist variable accordingly), if a variable is part of the 
  !   constituent array ("q" array), add the _exact_ name as in cnst_add call(e.g.  NUMLIQ, CLDICE etc.)
  !2. If the variable is not present in the constituent array,add a "case" statement for that variable in the "select case" 
  !   construct in get_var function in this module

  integer, public, parameter :: nvars_prtrb_hist = 11
  character(len=6), public, parameter :: hist_vars(nvars_prtrb_hist) = ['s     ', 't     ', 'Q     ', 'v     ', &
       'CLDLIQ', 'NUMLIQ', 'CLDICE', 'NUMICE', 'num_a1','num_a2','num_a3']
  
contains 

  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  subroutine physics_update_init()
    !Purpose: Initialize variables for physics update interface and pergro test

    use phys_control, only: phys_getopts
    use units,        only: getunit

    integer :: stat

    call phys_getopts(pergro_test_active_out = pergro_test_active)

    if (pergro_test_active) then
       !open file for writing ptend names
       pid = huge(1)
       if(masterproc) then
          pid = 0   !initialize pid for masterproc only
          unitn = getunit()
          open( unitn,file=fname,status='replace', form='formatted', action='write', position='append' )
          write(unitn,*,iostat=stat)'topphysbc' ! topphysbc record variables values at top of tphysbc
          if( stat > 0 ) then
             call endrun('Error writing topphysbc string  in '//fname//' file')
          endif
       endif
    endif

  end subroutine physics_update_init

  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------

  subroutine physics_update(state, ptend, dt, tend)
    !purpose: This subroutine calls physics_update_main (old physics_update)
    !and also output variables for pergro test

    use time_manager,  only: is_first_step

    
    !Arguments
    type(physics_ptend), intent(inout)  :: ptend   ! Parameterization tendencies
    type(physics_state), intent(inout)  :: state   ! Physics state variables
    real(r8),            intent(in)     :: dt      ! time step
    
    !optional arguments
    type(physics_tend ), intent(inout), optional  :: tend  ! Physics tendencies over timestep
    
    !Local vars
    character(len = fieldname_len)   :: pname, varname, vsuffix
    
    logical                          :: outfld_active, add_pname
    integer                          :: lchnk, stat, ip, ihist
    
    lchnk = state%lchnk
    
    !IMPORTANT:Store ptend%name as it will be modified in physics_update_main call
    pname = ptend%name
    
    !if nothing is to be updated in physics_update_main, DO NOT output (using outfld calls) 
    !PERGRO variables ("t_...", "s_..." etc.) below
    !Note: The following logical flag is required as sometimes "pname" is an empty string
    !      due to some stub routine calls (e.g. "iondrag_calc") where ptend is an 
    !      intent-out. This causes issues as intent-out will cause ptend%name to become undefined

    outfld_active = .true. !decides whether to call outfld calls or not
    if (.not. (any(ptend%lq(:)) .or. ptend%ls .or. ptend%lu .or. ptend%lv)) outfld_active = .false.
    
    !call the old physics update call
    call physics_update_main (state, ptend, dt, tend)
    
    if (pergro_test_active .and. outfld_active) then
       
       !write text file to be used for the post processing
       if(masterproc .and. lchnk == begchunk .and. is_first_step()) then
          !Here we write a text file to use for the post processing. We list
          !all the ptend names in this file. We do not want duplicates as it will
          !confuse the post processing script, therefore we skip ptend names when
          !they already exist in the plist

          !Find if this pname already exist in plist
          add_pname = .true.!decides whether to add pname to plist or not
          do ip = 1 , pid
             if (trim(adjustl(pname)) == trim(adjustl(plist(ip)))) then
                !already exists, do NOT add in plist
                add_pname = .false.
                exit
             endif
          enddo
          if (add_pname) then
             write(unitn,*,iostat=stat) pname
             if( stat > 0 ) then
                call endrun('Error writing '//pname// 'in '//fname//' file')
             endif
             !increment index pid and add pname to the list
             pid = pid + 1
             plist(pid) = trim(adjustl(pname))
          endif
       endif
          
       !call outfld
       do ihist = 1 , nvars_prtrb_hist
          vsuffix  = trim(adjustl(hist_vars(ihist)))
          varname  = trim(adjustl(vsuffix))//'_'//trim(adjustl(pname)) ! form variable name
          !find the prognostic variable associated with this hist_vars(ihist) via "get_var" function
          call outfld( trim(adjustl(varname)), get_var(state,vsuffix), pcols, lchnk )          
       enddo
    endif
  end subroutine physics_update
  
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------

  function get_var(state,hist_var) result (prg_var)
    !Purpose: Find which state variable to output based on the hist_var string

    use constituents, only: cnst_get_ind
    
    character(len=fieldname_len), intent(in)  :: hist_var
    type(physics_state),          intent(in)  :: state  
    real(r8)                                  :: prg_var(pcols,pver)

    !local vars
    integer :: idx

    !see if the hist_var exists in constituent array
    call cnst_get_ind(trim(adjustl(hist_var)), idx, abort=.false.)
    
    if (idx .ne. -1 ) then ! idx == -1  means, variable doesn't exists in the constituent array

       prg_var(1:pcols,1:pver) = state%q(1:pcols,1:pver,idx)

    else !variable doesn't exists in the constituent array

       select case (trim(adjustl(hist_var)))
       case('s')
          prg_var(1:pcols,1:pver) = state%s(1:pcols,1:pver)
       case('t')
          prg_var(1:pcols,1:pver) = state%t(1:pcols,1:pver)
       case('v')
          prg_var(1:pcols,1:pver) = state%v(1:pcols,1:pver)
       case default
          call endrun('physics_update_mod.F90 - func get_var, unrecognized variable: '// trim(adjustl(hist_var)))
       end select
    endif
  end function get_var
     
end module physics_update_mod
