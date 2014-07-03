

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 MODULE overflows

!testing out some changes here

!BOP
! !MODULE: overflows
! !DESCRIPTION:
!  This module contains data types and routines for computing 
!  parameterized overflows. Overflows are sub-grid scale flows
!  along topography thought to be important for bottom water
!  formation.
!
! !REVISION HISTORY:
!  SVN:
!  

! !USES:

   use POP_KindsMod
   use POP_ErrorMod
   use POP_BlocksMod
   use POP_CommMod
   use POP_ConfigMod
   use POP_DistributionMod
   use POP_DomainSizeMod
   use POP_FieldMod
   use POP_GridHorzMod
   use POP_HaloMod
   use POP_RedistributeMod

   use POP_SolversMod


   use blocks
   use broadcast
   use communicate
   use constants
   use domain
   use exit_mod
   use global_reductions
   use grid
   use io_types
   use kinds_mod
   use overflow_type
   use ovf_utils
   use prognostic
   use time_management
   use registry
   use state_mod

   use timers, only: timer_start, timer_stop, get_timer

   !*** ccsm
   use gather_scatter
   use shr_sys_mod

   implicit none
   private
   save

  include 'mpif.h'

! !PUBLIC MEMBER FUNCTIONS:

   public :: init_overflows1,           &   ! initial.F90
             init_overflows2,           &   ! initial.F90
             init_overflows_kmt,        &
             init_overflows_mask,       &
             init_overflows3,           &   ! initial.F90
             init_overflows4,           &   ! initial.F90
             init_overflows5,           &   ! initial.F90
             ovf_write_restart,         &   ! step_mod.F90
             ovf_read_restart,          &   ! init_overflows1
             ovf_read_broadcast,        &   
             ovf_advt,                  &   ! advection.F90
             ovf_wtkb_check,            &   ! advection.F90
             ovf_UV_check,              &  
             ovf_Utlda,                 &   ! baroclinic.F90
             ovf_driver,                &   ! step_mod.F90
             ovf_reg_avgs,              &
             ovf_transports,            &
             ovf_loc_prd,               &
             ovf_W,                     &
             ovf_UV,                    &
             ovf_rhs_brtrpc_momentum,   &   ! barotropic.F90
             ovf_brtrpc_renorm,         &   ! barotropic.F90
             ovf_rhs_brtrpc_continuity, &   ! barotropic.F90
             ovf_solvers_9pt,           &   ! barotropic.F90
             ovf_HU,                    &
             ovf_UV_solution                ! step_mod.F90

! !PUBLIC DATA MEMBERS:

!-----------------------------------------------------------------------
!     list of nomenclature definitions
!-----------------------------------------------------------------------
! 
! ovf    = overflow
! inf    = inflow (refering to inflow region)
! src    = source (either region or grid box)
! ent    = entrainment (either region or grid box)
! prd    = product (either region or grid box)
! reg    = region (for averaging density and tracers over region)
! adj    = adjacent (for averaging density and tracers over adjacent boxes)
! num    = number (usually refers to actual number used based on input)
! no.    = number (usually refers to actual number used based on input)
! locs   = locations (i.e. grid boxes)
! orient = orientation (1,2,3 or 4; refers to grid box sidewall)
! params = parameters
! ssb    = shelf-slope break- shelf/slope transition to abyssal depth
! 
!-----------------------------------------------------------------------
!     define overflow types and parameters
!-----------------------------------------------------------------------


   integer (POP_i4) ::     &
      errorCode

!EOC
!-----------------------------------------------------------------------
!
!  controls for frequency and output of diagnostics
!
!-----------------------------------------------------------------------

   integer (int_kind) ::  &
      ovf_diag_unit        ! i/o unit for overflow output diagnostics file

   character (char_len) ::  &
      ccsm_diag_date

   logical (log_kind) ::    &
      lccsm = .false.
   
   character (10) ::  &
      cdate            ! character date

!-----------------------------------------------------------------------
!
!  ovf timers
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
        timer_ovf         ! timer for ovf





!***********************************************************************

 contains

!***********************************************************************
!EOP
! !IROUTINE: init_overflows1
! !INTERFACE:

 subroutine init_overflows1

! !DESCRIPTION:
!  This routine is the first of four which together initialize the overflow 
!  parameterization. It reads the namelist and overflow_infile (text file 
!  containing ovf info). See info file comments for description of text file 
!  format. This routine also computes prd region limits based on prd input, 
!  writes out to stdout and overflows_diag_outfile, and then broadcasts ovf info 
!  to all processors.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------
   save

   namelist /overflows_nml/  overflows_on, overflows_interactive,      &
                             overflows_infile, overflows_diag_outfile,       &
                             overflows_restart_type, overflows_restfile

   integer (i4) ::          &
      index,                &! overflow index
      nu,                   &! unit for overflow input file
      nml_error,            &! namelist i/o error flag
      ovf_error,            &! ovf i/o error flag
      num_req,              &! number requested for error message
      imin,                 &! i index for checking input order
      jmin,                 &! j index for checking input order
      kmin,                 &! k index for checking input number of levels
      ornt                   ! orientation for checking constancy

   character (88) :: line ! temporary for line of text input/output

   integer   (int_kind)  :: &
      n,                    &! ovf loop index
      m,                    &! sub-ovf loop index
      nn,                   &! tracer loop index 
      mp,                   &! sub-ovf sub-loop index
      i,j,                  &! horizontal loop indices
      k,                    &! vertical loop index
      iblock,               &! local block address
      ib,ie,jb,je,          &! local domain index boundaries
      di,dj                  ! orientation adjustments in i,j

   type (block) ::          &
      this_block             ! block information for current block

!-----------------------------------------------------------------------
!  read overflow namelist
!-----------------------------------------------------------------------

   if (my_task == master_task) then
      write(stdout,blank_fmt)
      write(stdout,ndelim_fmt)
      write(stdout,blank_fmt)
      write(stdout,'(a11)') ' Overflows:'
      write(stdout,blank_fmt)
      call shr_sys_flush(stdout)
   endif

   overflows_on           = .false.
   overflows_interactive  = .false.
   overflows_infile       = 'unknown_ovf_infile'
   overflows_diag_outfile = 'unknown_ovf_outfile'
   overflows_restart_type = 'ccsm_startup'
   overflows_restfile     = 'unknown_ovf_restfile'

   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old',iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif
      do while (nml_error > 0)
         read(nml_in, nml=overflows_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   endif

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
      call exit_POP(sigAbort,'ERROR reading overflows_nml')
   endif

   if (my_task == master_task) then
      write(stdout,blank_fmt)
      write(stdout,'(a33)') ' overflows_nml namelist settings:'
      write(stdout,blank_fmt)
      write(stdout,overflows_nml)
      write(stdout,blank_fmt)
      call shr_sys_flush(stdout)
   endif

   call broadcast_scalar(overflows_on, master_task)
   call broadcast_scalar(overflows_interactive, master_task)
   call broadcast_scalar(overflows_infile, master_task)
   call broadcast_scalar(overflows_diag_outfile, master_task)
   call broadcast_scalar(overflows_restart_type, master_task)
   call broadcast_scalar(overflows_restfile, master_task)

!-----------------------------------------------------------------------
!  if overflows off, exit
!-----------------------------------------------------------------------

   if( .not. overflows_on ) return



! don't initialize the groups yet
  ovf_groups%init = .false.


!-----------------------------------------------------------------------
!
!  determine if this is a ccsm coupled run
!-----------------------------------------------------------------------

   lccsm                    = registry_match('lccsm')

!-----------------------------------------------------------------------
!  overflows on; read overflows info file if ccsm_startup; otherwise 
!  read restart data
!-----------------------------------------------------------------------

   if( overflows_restart_type == 'ccsm_startup' ) then

   ovf_error = 0
   call get_unit(nu)

!-----------------------------------------------------------------------
!  master task section
!-----------------------------------------------------------------------

   if (my_task == master_task) then

      open(nu, file=overflows_infile, status='old',iostat=ovf_error)

      write(stdout,2345) ovf_error
      2345 format(' after open nu   ovf_error=',i5)
      write(stdout,'(a41)') 'reading overflows_infile: contents echoed'
      call shr_sys_flush(stdout)

      do m=1,40
         read(nu,'(a88)') line
         write(stdout,'(a88)') line
      end do

      read(nu,*) num_ovf
      write(stdout,*) num_ovf
      call shr_sys_flush(stdout)
      if( num_ovf <= 0 .or. num_ovf > max_ovf ) then
         ovf_error = 1
         num_req   = num_ovf
         goto 10
      end if

      do n=1,num_ovf
        ovf(n)%interactive = overflows_interactive
        read(nu,*) index,ovf(n)%name
        write(stdout,*) index,ovf(n)%name

        read(nu,*) ovf(n)%ovf_params%lat
        read(nu,*) ovf(n)%ovf_params%width
        read(nu,*) ovf(n)%ovf_params%source_thick
        read(nu,*) ovf(n)%ovf_params%distnc_str_ssb
        read(nu,*) ovf(n)%ovf_params%bottom_slope
        read(nu,*) ovf(n)%ovf_params%bottom_drag

        write(stdout,*) ovf(n)%ovf_params%lat
        write(stdout,*) ovf(n)%ovf_params%width
        write(stdout,*) ovf(n)%ovf_params%source_thick
        write(stdout,*) ovf(n)%ovf_params%distnc_str_ssb
        write(stdout,*) ovf(n)%ovf_params%bottom_slope
        write(stdout,*) ovf(n)%ovf_params%bottom_drag
        call shr_sys_flush(stdout)

! kmt changes if any
        read(nu,*) ovf(n)%num_kmt
        write(stdout,*) ovf(n)%num_kmt
        if( ovf(n)%num_kmt < 0 .or. ovf(n)%num_kmt > max_kmt ) then
           ovf_error = 2
           num_req   = ovf(n)%num_kmt
           goto 10
        endif
        do m=1,ovf(n)%num_kmt
           read(nu,*) ovf(n)%loc_kmt(m)%i,    &
                      ovf(n)%loc_kmt(m)%j,    &
                      ovf(n)%loc_kmt(m)%korg, &
                      ovf(n)%loc_kmt(m)%knew
           write(stdout,*) ovf(n)%loc_kmt(m)%i,    &
                           ovf(n)%loc_kmt(m)%j,    &
                           ovf(n)%loc_kmt(m)%korg, &
                           ovf(n)%loc_kmt(m)%knew
        end do
        call shr_sys_flush(stdout)

        read(nu,*)

! inf,src and ent region limits
        read(nu,*) ovf(n)%reg_inf%imin, &
                   ovf(n)%reg_inf%imax, &
                   ovf(n)%reg_inf%jmin, &
                   ovf(n)%reg_inf%jmax, &
                   ovf(n)%reg_inf%kmin, &
                   ovf(n)%reg_inf%kmax
        read(nu,*) ovf(n)%reg_src%imin, &
                   ovf(n)%reg_src%imax, &
                   ovf(n)%reg_src%jmin, &
                   ovf(n)%reg_src%jmax, &
                   ovf(n)%reg_src%kmin, &
                   ovf(n)%reg_src%kmax
        read(nu,*) ovf(n)%reg_ent%imin, &
                   ovf(n)%reg_ent%imax, &
                   ovf(n)%reg_ent%jmin, &
                   ovf(n)%reg_ent%jmax, &
                   ovf(n)%reg_ent%kmin, &
                   ovf(n)%reg_ent%kmax

        write(stdout,*) ovf(n)%reg_inf%imin, &
                        ovf(n)%reg_inf%imax, &
                        ovf(n)%reg_inf%jmin, &
                        ovf(n)%reg_inf%jmax, &
                        ovf(n)%reg_inf%kmin, &
                        ovf(n)%reg_inf%kmax
        write(stdout,*) ovf(n)%reg_src%imin, &
                        ovf(n)%reg_src%imax, &
                        ovf(n)%reg_src%jmin, &
                        ovf(n)%reg_src%jmax, &
                        ovf(n)%reg_src%kmin, &
                        ovf(n)%reg_src%kmax
        write(stdout,*) ovf(n)%reg_ent%imin, &
                        ovf(n)%reg_ent%imax, &
                        ovf(n)%reg_ent%jmin, &
                        ovf(n)%reg_ent%jmax, &
                        ovf(n)%reg_ent%kmin, &
                        ovf(n)%reg_ent%kmax
        call shr_sys_flush(stdout)

! src points
        read(nu,*) ovf(n)%num_src
        write(stdout,*) ovf(n)%num_src
        if( ovf(n)%num_src <= 1 .or. ovf(n)%num_src > max_src ) then
           ovf_error = 3
           num_req   = ovf(n)%num_src
           goto 10
        endif
        do m=1,ovf(n)%num_src
          read(nu,*) ovf(n)%loc_src(m)%i, &
                     ovf(n)%loc_src(m)%j, &
                     ovf(n)%loc_src(m)%k, &
                     ovf(n)%loc_src(m)%orient
          if ( ovf(n)%loc_src(m)%orient .eq. 1 ) then
             ovf(n)%loc_src(m)%i_adv = ovf(n)%loc_src(m)%i + 1
             ovf(n)%loc_src(m)%j_adv = ovf(n)%loc_src(m)%j
             ovf(n)%loc_src(m)%i_u   = ovf(n)%loc_src(m)%i
             ovf(n)%loc_src(m)%j_u   = ovf(n)%loc_src(m)%j
! some u corner points i_u,j_u zeroed because they are inactive
             if( m == ovf(n)%num_src ) then
                ovf(n)%loc_src(m)%i_u   = 0
                ovf(n)%loc_src(m)%j_u   = 0
             endif
          else if( ovf(n)%loc_src(m)%orient .eq. 2 ) then
             ovf(n)%loc_src(m)%i_adv = ovf(n)%loc_src(m)%i 
             ovf(n)%loc_src(m)%j_adv = ovf(n)%loc_src(m)%j + 1
             ovf(n)%loc_src(m)%i_u   = ovf(n)%loc_src(m)%i - 1
             if(ovf(n)%loc_src(m)%i_u==0) ovf(n)%loc_src(m)%i_u = nx_global
             ovf(n)%loc_src(m)%j_u   = ovf(n)%loc_src(m)%j
             if( m == 1 ) then
                ovf(n)%loc_src(m)%i_u   = 0
                ovf(n)%loc_src(m)%j_u   = 0
             endif
          else if( ovf(n)%loc_src(m)%orient .eq. 3 ) then
             ovf(n)%loc_src(m)%i_adv = ovf(n)%loc_src(m)%i - 1
             if(ovf(n)%loc_src(m)%i_adv==0) ovf(n)%loc_src(m)%i_adv = nx_global
             ovf(n)%loc_src(m)%j_adv = ovf(n)%loc_src(m)%j
             ovf(n)%loc_src(m)%i_u   = ovf(n)%loc_src(m)%i - 1
             if(ovf(n)%loc_src(m)%i_u==0) ovf(n)%loc_src(m)%i_u = nx_global
             ovf(n)%loc_src(m)%j_u   = ovf(n)%loc_src(m)%j - 1
             if( m == 1 ) then
                ovf(n)%loc_src(m)%i_u   = 0
                ovf(n)%loc_src(m)%j_u   = 0
             endif
          else if( ovf(n)%loc_src(m)%orient .eq. 4 ) then
             ovf(n)%loc_src(m)%i_adv = ovf(n)%loc_src(m)%i 
             ovf(n)%loc_src(m)%j_adv = ovf(n)%loc_src(m)%j - 1
             ovf(n)%loc_src(m)%i_u   = ovf(n)%loc_src(m)%i
             ovf(n)%loc_src(m)%j_u   = ovf(n)%loc_src(m)%j - 1
             if( m == ovf(n)%num_src ) then
                ovf(n)%loc_src(m)%i_u   = 0
                ovf(n)%loc_src(m)%j_u   = 0
             endif
          endif
          write(stdout,*) ovf(n)%loc_src(m)%i, &
                          ovf(n)%loc_src(m)%j, &
                          ovf(n)%loc_src(m)%k, &
                          ovf(n)%loc_src(m)%orient
          call shr_sys_flush(stdout)
! check order of ij, constancy of k and range of orient
          if( m==1 ) then
             imin = ovf(n)%loc_src(m)%i 
             jmin = ovf(n)%loc_src(m)%j
             kmin = ovf(n)%loc_src(m)%k
             if( ovf(n)%loc_src(m)%orient < 1 .or.  &
                 ovf(n)%loc_src(m)%orient > 4 ) then 
                ovf_error = 11
                goto 10
             endif
             ornt = ovf(n)%loc_src(m)%orient
          else
             if( ovf(n)%loc_src(m)%i < imin ) then
                ovf_error = 7
                goto 10
             endif
             if( ovf(n)%loc_src(m)%j < jmin ) then
                ovf_error = 7
                goto 10
             endif
             if( ovf(n)%loc_src(m)%i == imin .and. &
                 ovf(n)%loc_src(m)%j == jmin) then
                ovf_error = 8
                goto 10
             endif
             if( ovf(n)%loc_src(m)%i > imin .and. &
                 ovf(n)%loc_src(m)%j > jmin) then
                ovf_error = 9
                goto 10
             endif
             if( ovf(n)%loc_src(m)%k /= kmin ) then
                ovf_error = 10
                goto 10
             endif
             if( ovf(n)%loc_src(m)%orient < 1 .or.  &
                 ovf(n)%loc_src(m)%orient > 4 ) then 
                ovf_error = 11
                goto 10
             endif
             if( ovf(n)%loc_src(m)%orient /= ornt ) then
                ovf_error = 12
                goto 10
             endif
             imin = ovf(n)%loc_src(m)%i 
             jmin = ovf(n)%loc_src(m)%j
             kmin = ovf(n)%loc_src(m)%k
             ornt = ovf(n)%loc_src(m)%orient
          endif
        end do

! ent points
        read(nu,*) ovf(n)%num_ent
        write(stdout,*) ovf(n)%num_ent
        if( ovf(n)%num_ent <= 1 .or. ovf(n)%num_ent > max_ent ) then
           ovf_error = 4
           num_req   = ovf(n)%num_ent
           goto 10
        endif
        do m=1,ovf(n)%num_ent
          read(nu,*) ovf(n)%loc_ent(m)%i, &
                     ovf(n)%loc_ent(m)%j, &
                     ovf(n)%loc_ent(m)%k, &
                     ovf(n)%loc_ent(m)%orient
          if ( ovf(n)%loc_ent(m)%orient .eq. 1 ) then
             ovf(n)%loc_ent(m)%i_adv = ovf(n)%loc_ent(m)%i + 1
             ovf(n)%loc_ent(m)%j_adv = ovf(n)%loc_ent(m)%j
             ovf(n)%loc_ent(m)%i_u   = ovf(n)%loc_ent(m)%i
             ovf(n)%loc_ent(m)%j_u   = ovf(n)%loc_ent(m)%j
! some u corner points i_u,j_u zeroed because they are inactive
             if( m == ovf(n)%num_ent ) then
                ovf(n)%loc_ent(m)%i_u   = 0
                ovf(n)%loc_ent(m)%j_u   = 0
             endif
          else if( ovf(n)%loc_ent(m)%orient .eq. 2 ) then
             ovf(n)%loc_ent(m)%i_adv = ovf(n)%loc_ent(m)%i 
             ovf(n)%loc_ent(m)%j_adv = ovf(n)%loc_ent(m)%j + 1
             ovf(n)%loc_ent(m)%i_u   = ovf(n)%loc_ent(m)%i - 1
             if(ovf(n)%loc_ent(m)%i_u==0) ovf(n)%loc_ent(m)%i_u = nx_global
             ovf(n)%loc_ent(m)%j_u   = ovf(n)%loc_ent(m)%j
             if( m == 1 ) then
                ovf(n)%loc_ent(m)%i_u   = 0
                ovf(n)%loc_ent(m)%j_u   = 0
             endif
          else if( ovf(n)%loc_ent(m)%orient .eq. 3 ) then
             ovf(n)%loc_ent(m)%i_adv = ovf(n)%loc_ent(m)%i - 1
             if(ovf(n)%loc_ent(m)%i_adv==0) ovf(n)%loc_ent(m)%i_adv = nx_global
             ovf(n)%loc_ent(m)%j_adv = ovf(n)%loc_ent(m)%j
             ovf(n)%loc_ent(m)%i_u   = ovf(n)%loc_ent(m)%i - 1
             if(ovf(n)%loc_ent(m)%i_u==0) ovf(n)%loc_ent(m)%i_u = nx_global
             ovf(n)%loc_ent(m)%j_u   = ovf(n)%loc_ent(m)%j - 1
             if( m == 1 ) then
                ovf(n)%loc_ent(m)%i_u   = 0
                ovf(n)%loc_ent(m)%j_u   = 0
             endif
          else if( ovf(n)%loc_ent(m)%orient .eq. 4 ) then
             ovf(n)%loc_ent(m)%i_adv = ovf(n)%loc_ent(m)%i 
             ovf(n)%loc_ent(m)%j_adv = ovf(n)%loc_ent(m)%j - 1
             ovf(n)%loc_ent(m)%i_u   = ovf(n)%loc_ent(m)%i
             ovf(n)%loc_ent(m)%j_u   = ovf(n)%loc_ent(m)%j - 1
             if( m == ovf(n)%num_ent ) then
                ovf(n)%loc_ent(m)%i_u   = 0
                ovf(n)%loc_ent(m)%j_u   = 0
             endif
          endif
          write(stdout,*) ovf(n)%loc_ent(m)%i, &
                          ovf(n)%loc_ent(m)%j, &
                          ovf(n)%loc_ent(m)%k, &
                          ovf(n)%loc_ent(m)%orient
          call shr_sys_flush(stdout)
! check order of ij, constancy of k and range of orient
          if( m==1 ) then
             imin = ovf(n)%loc_ent(m)%i
             jmin = ovf(n)%loc_ent(m)%j
             kmin = ovf(n)%loc_ent(m)%k
             if( ovf(n)%loc_ent(m)%orient < 1 .or.  &
                 ovf(n)%loc_ent(m)%orient > 4 ) then
                ovf_error = 11
                goto 10
             endif
             ornt = ovf(n)%loc_ent(m)%orient
          else
             if( ovf(n)%loc_ent(m)%i < imin ) then
                ovf_error = 7
                goto 10
             endif
             if( ovf(n)%loc_ent(m)%j < jmin ) then
                ovf_error = 7
                goto 10
             endif
             if( ovf(n)%loc_ent(m)%i == imin .and. &
                 ovf(n)%loc_ent(m)%j == jmin) then
                ovf_error = 8
                goto 10
             endif
             if( ovf(n)%loc_ent(m)%i > imin .and. &
                 ovf(n)%loc_ent(m)%j > jmin) then
                ovf_error = 9
                goto 10
             endif
             if( ovf(n)%loc_ent(m)%k /= kmin ) then
                ovf_error = 10
                goto 10
             endif
             if( ovf(n)%loc_ent(m)%orient < 1 .or.  &
                 ovf(n)%loc_ent(m)%orient > 4 ) then
                ovf_error = 11
                goto 10
             endif
             if( ovf(n)%loc_ent(m)%orient /= ornt ) then
                ovf_error = 12
                goto 10
             endif
             imin = ovf(n)%loc_ent(m)%i
             jmin = ovf(n)%loc_ent(m)%j
             kmin = ovf(n)%loc_ent(m)%k
             ornt = ovf(n)%loc_ent(m)%orient
          endif
        end do
        call shr_sys_flush(stdout)

! prd points
        read(nu,*) ovf(n)%num_prd_sets
        write(stdout,*) ovf(n)%num_prd_sets
        if(ovf(n)%num_prd_sets<=0.or.ovf(n)%num_prd_sets>max_prd_sets) then
           ovf_error = 5
           num_req   = ovf(n)%num_prd_sets
           goto 10
        endif
        do m=1,ovf(n)%num_prd_sets
          read(nu,*) ovf(n)%num_prd(m)
          write(stdout,*) ovf(n)%num_prd(m)
          if( ovf(n)%num_prd(m)<=1.or.ovf(n)%num_prd(m)>max_prd) then
             ovf_error = 6
             num_req   = ovf(n)%num_prd(m)
             goto 10
          endif
          do mp=1,ovf(n)%num_prd(m) 
            read(nu,*) ovf(n)%loc_prd(m,mp)%i, &
                       ovf(n)%loc_prd(m,mp)%j, &
                       ovf(n)%loc_prd(m,mp)%k, &
                       ovf(n)%loc_prd(m,mp)%orient
            if ( ovf(n)%loc_prd(m,mp)%orient .eq. 1 ) then
               ovf(n)%loc_prd(m,mp)%i_adv = ovf(n)%loc_prd(m,mp)%i + 1
               ovf(n)%loc_prd(m,mp)%j_adv = ovf(n)%loc_prd(m,mp)%j
               ovf(n)%loc_prd(m,mp)%i_u   = ovf(n)%loc_prd(m,mp)%i
               ovf(n)%loc_prd(m,mp)%j_u   = ovf(n)%loc_prd(m,mp)%j
! some u corner points i_u,j_u zeroed because they are inactive
               if( mp == ovf(n)%num_prd(m) ) then
                  ovf(n)%loc_prd(m,mp)%i_u   = 0
                  ovf(n)%loc_prd(m,mp)%j_u   = 0
               endif
            else if( ovf(n)%loc_prd(m,mp)%orient .eq. 2 ) then
               ovf(n)%loc_prd(m,mp)%i_adv = ovf(n)%loc_prd(m,mp)%i
               ovf(n)%loc_prd(m,mp)%j_adv = ovf(n)%loc_prd(m,mp)%j + 1
               ovf(n)%loc_prd(m,mp)%i_u   = ovf(n)%loc_prd(m,mp)%i - 1
               if(ovf(n)%loc_prd(m,mp)%i_u==0) ovf(n)%loc_prd(m,mp)%i_u = nx_global
               ovf(n)%loc_prd(m,mp)%j_u   = ovf(n)%loc_prd(m,mp)%j
               if( mp == 1 ) then
                  ovf(n)%loc_prd(m,mp)%i_u   = 0
                  ovf(n)%loc_prd(m,mp)%j_u   = 0
               endif
            else if( ovf(n)%loc_prd(m,mp)%orient .eq. 3 ) then
               ovf(n)%loc_prd(m,mp)%i_adv = ovf(n)%loc_prd(m,mp)%i - 1
               if(ovf(n)%loc_prd(m,mp)%i_adv==0) ovf(n)%loc_prd(m,mp)%i_adv = nx_global
               ovf(n)%loc_prd(m,mp)%j_adv = ovf(n)%loc_prd(m,mp)%j
               ovf(n)%loc_prd(m,mp)%i_u   = ovf(n)%loc_prd(m,mp)%i - 1
               if(ovf(n)%loc_prd(m,mp)%i_u==0) ovf(n)%loc_prd(m,mp)%i_u = nx_global
               ovf(n)%loc_prd(m,mp)%j_u   = ovf(n)%loc_prd(m,mp)%j - 1
               if( mp == 1 ) then
                  ovf(n)%loc_prd(m,mp)%i_u   = 0
                  ovf(n)%loc_prd(m,mp)%j_u   = 0
               endif
            else if( ovf(n)%loc_prd(m,mp)%orient .eq. 4 ) then
               ovf(n)%loc_prd(m,mp)%i_adv = ovf(n)%loc_prd(m,mp)%i
               ovf(n)%loc_prd(m,mp)%j_adv = ovf(n)%loc_prd(m,mp)%j - 1
               ovf(n)%loc_prd(m,mp)%i_u   = ovf(n)%loc_prd(m,mp)%i
               ovf(n)%loc_prd(m,mp)%j_u   = ovf(n)%loc_prd(m,mp)%j - 1
               if( mp == ovf(n)%num_prd(m) ) then
                  ovf(n)%loc_prd(m,mp)%i_u   = 0
                  ovf(n)%loc_prd(m,mp)%j_u   = 0
               endif
            endif
            write(stdout,*) ovf(n)%loc_prd(m,mp)%i, &
                            ovf(n)%loc_prd(m,mp)%j, &
                            ovf(n)%loc_prd(m,mp)%k, &
                            ovf(n)%loc_prd(m,mp)%orient
            call shr_sys_flush(stdout)
! check order of ij, constancy of k and range of orient
            if( mp==1 ) then
               imin = ovf(n)%loc_prd(m,mp)%i
               jmin = ovf(n)%loc_prd(m,mp)%j
               kmin = ovf(n)%loc_prd(m,mp)%k
               if( ovf(n)%loc_prd(m,mp)%orient < 1 .or.  &
                   ovf(n)%loc_prd(m,mp)%orient > 4 ) then
                  ovf_error = 11
                  goto 10
               endif
               ornt = ovf(n)%loc_prd(m,mp)%orient
            else
               if( ovf(n)%loc_prd(m,mp)%i < imin ) then
                  ovf_error = 7
                  goto 10
               endif
               if( ovf(n)%loc_prd(m,mp)%j < jmin ) then
                  ovf_error = 7
                  goto 10
               endif
               if( ovf(n)%loc_prd(m,mp)%i == imin .and. &
                   ovf(n)%loc_prd(m,mp)%j == jmin) then
                  ovf_error = 8
                  goto 10
               endif
               if( ovf(n)%loc_prd(m,mp)%i > imin .and. &
                   ovf(n)%loc_prd(m,mp)%j > jmin) then
                  ovf_error = 9
                  goto 10
               endif
               if( ovf(n)%loc_prd(m,mp)%k /= kmin ) then
                  ovf_error = 10
                  goto 10
               endif
               if( ovf(n)%loc_prd(m,mp)%orient < 1 .or.  &
                   ovf(n)%loc_prd(m,mp)%orient > 4 ) then
                  ovf_error = 11
                  goto 10
               endif
               if( ovf(n)%loc_prd(m,mp)%orient /= ornt ) then
                  ovf_error = 12
                  goto 10
               endif
               imin = ovf(n)%loc_prd(m,mp)%i
               jmin = ovf(n)%loc_prd(m,mp)%j
               kmin = ovf(n)%loc_prd(m,mp)%k
               ornt = ovf(n)%loc_prd(m,mp)%orient
            endif
          end do
          call shr_sys_flush(stdout)
        end do

! find src adj limits
        di = 0
        dj = 0
        if( ovf(n)%loc_src(1)%orient .eq. 1 ) di = +1
        if( ovf(n)%loc_src(1)%orient .eq. 2 ) dj = +1
        if( ovf(n)%loc_src(1)%orient .eq. 3 ) di = -1
        if( ovf(n)%loc_src(1)%orient .eq. 4 ) dj = -1
        ovf(n)%adj_src%imin = ovf(n)%loc_src(1)%i+di
        ovf(n)%adj_src%jmin = ovf(n)%loc_src(1)%j+dj
        ovf(n)%adj_src%kmin = ovf(n)%loc_src(1)%k
        ovf(n)%adj_src%imax = ovf(n)%loc_src(1)%i+di
        ovf(n)%adj_src%jmax = ovf(n)%loc_src(1)%j+dj
        ovf(n)%adj_src%kmax = ovf(n)%loc_src(1)%k
        do m=2,ovf(n)%num_src 
          di = 0
          dj = 0
          if( ovf(n)%loc_src(m)%orient .eq. 1 ) di = +1
          if( ovf(n)%loc_src(m)%orient .eq. 2 ) dj = +1
          if( ovf(n)%loc_src(m)%orient .eq. 3 ) di = -1
          if( ovf(n)%loc_src(m)%orient .eq. 4 ) dj = -1
          ovf(n)%adj_src%imin = min &
               (ovf(n)%adj_src%imin,ovf(n)%loc_src(m)%i+di)
          ovf(n)%adj_src%jmin = min &
               (ovf(n)%adj_src%jmin,ovf(n)%loc_src(m)%j+dj)
          ovf(n)%adj_src%kmin = min &
               (ovf(n)%adj_src%kmin,ovf(n)%loc_src(m)%k)
          ovf(n)%adj_src%imax = max &
               (ovf(n)%adj_src%imax,ovf(n)%loc_src(m)%i+di)
          ovf(n)%adj_src%jmax = max &
               (ovf(n)%adj_src%jmax,ovf(n)%loc_src(m)%j+dj)
          ovf(n)%adj_src%kmax = max &
               (ovf(n)%adj_src%kmax,ovf(n)%loc_src(m)%k)
        end do
! print src adj limits
        write(stdout,13)       &
          ovf(n)%adj_src%imin, &
          ovf(n)%adj_src%imax, &
          ovf(n)%adj_src%jmin, &
          ovf(n)%adj_src%jmax, &
          ovf(n)%adj_src%kmin, &
          ovf(n)%adj_src%kmax 
13        format(' Computed source adjacent ijk min/max =',6(i4,2x))

! find ent adj limits
        di = 0
        dj = 0
        if( ovf(n)%loc_ent(1)%orient .eq. 1 ) di = +1
        if( ovf(n)%loc_ent(1)%orient .eq. 2 ) dj = +1
        if( ovf(n)%loc_ent(1)%orient .eq. 3 ) di = -1
        if( ovf(n)%loc_ent(1)%orient .eq. 4 ) dj = -1
        ovf(n)%adj_ent%imin = ovf(n)%loc_ent(1)%i+di
        ovf(n)%adj_ent%jmin = ovf(n)%loc_ent(1)%j+dj
        ovf(n)%adj_ent%kmin = ovf(n)%loc_ent(1)%k
        ovf(n)%adj_ent%imax = ovf(n)%loc_ent(1)%i+di
        ovf(n)%adj_ent%jmax = ovf(n)%loc_ent(1)%j+dj
        ovf(n)%adj_ent%kmax = ovf(n)%loc_ent(1)%k
        do m=2,ovf(n)%num_ent
          di = 0
          dj = 0
          if( ovf(n)%loc_ent(m)%orient .eq. 1 ) di = +1
          if( ovf(n)%loc_ent(m)%orient .eq. 2 ) dj = +1
          if( ovf(n)%loc_ent(m)%orient .eq. 3 ) di = -1
          if( ovf(n)%loc_ent(m)%orient .eq. 4 ) dj = -1
          ovf(n)%adj_ent%imin = min &
               (ovf(n)%adj_ent%imin,ovf(n)%loc_ent(m)%i+di)
          ovf(n)%adj_ent%jmin = min &
               (ovf(n)%adj_ent%jmin,ovf(n)%loc_ent(m)%j+dj)
          ovf(n)%adj_ent%kmin = min &
               (ovf(n)%adj_ent%kmin,ovf(n)%loc_ent(m)%k)
          ovf(n)%adj_ent%imax = max &
               (ovf(n)%adj_ent%imax,ovf(n)%loc_ent(m)%i+di)
          ovf(n)%adj_ent%jmax = max &
               (ovf(n)%adj_ent%jmax,ovf(n)%loc_ent(m)%j+dj)
          ovf(n)%adj_ent%kmax = max &
               (ovf(n)%adj_ent%kmax,ovf(n)%loc_ent(m)%k)
        end do
! print ent adj limits
        write(stdout,14)       &
          ovf(n)%adj_ent%imin, &
          ovf(n)%adj_ent%imax, &
          ovf(n)%adj_ent%jmin, &
          ovf(n)%adj_ent%jmax, &
          ovf(n)%adj_ent%kmin, &
          ovf(n)%adj_ent%kmax 
14        format(' Computed entrainment adjacent ijk min/max =',6(i4,2x))

! find prd adj limits
        do m=1,ovf(n)%num_prd_sets
          di = 0
          dj = 0
          if( ovf(n)%loc_prd(m,1)%orient .eq. 1 ) di = +1
          if( ovf(n)%loc_prd(m,1)%orient .eq. 2 ) dj = +1
          if( ovf(n)%loc_prd(m,1)%orient .eq. 3 ) di = -1
          if( ovf(n)%loc_prd(m,1)%orient .eq. 4 ) dj = -1
          ovf(n)%adj_prd(m)%imin = ovf(n)%loc_prd(m,1)%i+di
          ovf(n)%adj_prd(m)%jmin = ovf(n)%loc_prd(m,1)%j+dj
          ovf(n)%adj_prd(m)%kmin = ovf(n)%loc_prd(m,1)%k
          ovf(n)%adj_prd(m)%imax = ovf(n)%loc_prd(m,1)%i+di
          ovf(n)%adj_prd(m)%jmax = ovf(n)%loc_prd(m,1)%j+dj
          ovf(n)%adj_prd(m)%kmax = ovf(n)%loc_prd(m,1)%k
          do mp=2,ovf(n)%num_prd(m) 
            di = 0
            dj = 0
            if( ovf(n)%loc_prd(m,mp)%orient .eq. 1 ) di = +1
            if( ovf(n)%loc_prd(m,mp)%orient .eq. 2 ) dj = +1
            if( ovf(n)%loc_prd(m,mp)%orient .eq. 3 ) di = -1
            if( ovf(n)%loc_prd(m,mp)%orient .eq. 4 ) dj = -1
            ovf(n)%adj_prd(m)%imin = min &
                 (ovf(n)%adj_prd(m)%imin,ovf(n)%loc_prd(m,mp)%i+di)
            ovf(n)%adj_prd(m)%jmin = min &
                 (ovf(n)%adj_prd(m)%jmin,ovf(n)%loc_prd(m,mp)%j+dj)
            ovf(n)%adj_prd(m)%kmin = min &
                 (ovf(n)%adj_prd(m)%kmin,ovf(n)%loc_prd(m,mp)%k)
            ovf(n)%adj_prd(m)%imax = max &
                 (ovf(n)%adj_prd(m)%imax,ovf(n)%loc_prd(m,mp)%i+di)
            ovf(n)%adj_prd(m)%jmax = max &
                 (ovf(n)%adj_prd(m)%jmax,ovf(n)%loc_prd(m,mp)%j+dj)
            ovf(n)%adj_prd(m)%kmax = max &
                 (ovf(n)%adj_prd(m)%kmax,ovf(n)%loc_prd(m,mp)%k)
          end do
        end do
! print prd adj limits
        do m=1,ovf(n)%num_prd_sets
          write(stdout,15) m,       &
            ovf(n)%adj_prd(m)%imin, &
            ovf(n)%adj_prd(m)%imax, &
            ovf(n)%adj_prd(m)%jmin, &
            ovf(n)%adj_prd(m)%jmax, &
            ovf(n)%adj_prd(m)%kmin, &
            ovf(n)%adj_prd(m)%kmax 
15        format(' Computed product adjacent, set=',i3, &
                 ' ijk min/max =',6(i4,2x))
        end do
      end do  ! ovf loop
      call shr_sys_flush(stdout)

!-----------------------------------------------------------------------
!  end master task section
!-----------------------------------------------------------------------

      close (nu)
   endif  ! master_task
   call release_unit(nu)

! error from goto 10
10 continue
   
   call broadcast_scalar(ovf_error, master_task)
   if (ovf_error /= 0) then
      call broadcast_scalar(num_req, master_task)
      write(stdout,*) 'ERROR on overflow input'
      if( ovf_error == 1 ) then
         write(stdout,*) 'Overflows on but number requested out of range'
         write(stdout,*) 'Number requested = ',num_req
         write(stdout,*) 'Must be > 0 and not greater than ',max_ovf
      else if ( ovf_error == 2 ) then
         write(stdout,*) 'Overflows on with kmt topography changes out of range'
         write(stdout,*) 'Number requested = ',num_req
         write(stdout,*) 'Must be >= 0 and not greater than ',max_kmt
      else if ( ovf_error == 3 ) then
         write(stdout,*) 'Overflows on with number source points out of range'
         write(stdout,*) 'Number requested = ',num_req
         write(stdout,*) 'Must be > 1 and not greater than ',max_src
      else if ( ovf_error == 4 ) then
         write(stdout,*) 'Overflows on with number entrainment points out of range'
         write(stdout,*) 'Number requested = ',num_req
         write(stdout,*) 'Must be > 1 and not greater than ',max_ent
      else if ( ovf_error == 5 ) then
         write(stdout,*) 'Overflows on with number of product sets out of range'
         write(stdout,*) 'Number requested = ',num_req
         write(stdout,*) 'Must be > 0 and not greater than ',max_prd_sets
      else if ( ovf_error == 6 ) then
         write(stdout,*) 'Overflows on with number of product points out of range'
         write(stdout,*) 'Number requested = ',num_req
         write(stdout,*) 'Must be > 1 and not greater than ',max_prd
      else if ( ovf_error == 7 ) then
         write(stdout,*) 'Overflows on with non-monotonic increasing i or j'
      else if ( ovf_error == 8 ) then
         write(stdout,*) 'Overflows on with no change in i and j'
      else if ( ovf_error == 9 ) then
         write(stdout,*) 'Overflows on with both i and j increasing'
      else if ( ovf_error == 10 ) then
         write(stdout,*) 'Overflows on with non-constant level k'
      else if ( ovf_error == 11 ) then
         write(stdout,*) 'Overflows on with orientation either < 0 or > 4'
      else if ( ovf_error == 12 ) then
         write(stdout,*) 'Overflows on with non-constant orientation'
      endif
      call shr_sys_flush(stdout)
      call exit_POP(sigAbort,'ERROR reading overflows_infile')
   endif  ! ovf error

!-----------------------------------------------------------------------
!  broadcast overflows info to all processors
!-----------------------------------------------------------------------

   call broadcast_scalar(num_ovf, master_task)
   do n=1,num_ovf
      call broadcast_scalar(ovf(n)%interactive, master_task)
      call broadcast_scalar(ovf(n)%name, master_task)
! ovf data
      call broadcast_scalar(ovf(n)%ovf_params%lat, master_task)
      call broadcast_scalar(ovf(n)%ovf_params%width, master_task)
      call broadcast_scalar(ovf(n)%ovf_params%source_thick, master_task)
      call broadcast_scalar(ovf(n)%ovf_params%distnc_str_ssb, master_task)
      call broadcast_scalar(ovf(n)%ovf_params%bottom_slope, master_task)
      call broadcast_scalar(ovf(n)%ovf_params%bottom_drag, master_task)
! kmt locations
      call broadcast_scalar(ovf(n)%num_kmt, master_task)
      do m=1,ovf(n)%num_kmt
         call broadcast_scalar(ovf(n)%loc_kmt(m)%i, master_task)
         call broadcast_scalar(ovf(n)%loc_kmt(m)%j, master_task)
         call broadcast_scalar(ovf(n)%loc_kmt(m)%korg, master_task)
         call broadcast_scalar(ovf(n)%loc_kmt(m)%knew, master_task)
      end do
! regional boundaries
!   inflow
      call broadcast_scalar(ovf(n)%reg_inf%imin, master_task)
      call broadcast_scalar(ovf(n)%reg_inf%imax, master_task)
      call broadcast_scalar(ovf(n)%reg_inf%jmin, master_task)
      call broadcast_scalar(ovf(n)%reg_inf%jmax, master_task)
      call broadcast_scalar(ovf(n)%reg_inf%kmin, master_task)
      call broadcast_scalar(ovf(n)%reg_inf%kmax, master_task)
!   source
      call broadcast_scalar(ovf(n)%reg_src%imin, master_task)
      call broadcast_scalar(ovf(n)%reg_src%imax, master_task)
      call broadcast_scalar(ovf(n)%reg_src%jmin, master_task)
      call broadcast_scalar(ovf(n)%reg_src%jmax, master_task)
      call broadcast_scalar(ovf(n)%reg_src%kmin, master_task)
      call broadcast_scalar(ovf(n)%reg_src%kmax, master_task)
!   entrainment
      call broadcast_scalar(ovf(n)%reg_ent%imin, master_task)
      call broadcast_scalar(ovf(n)%reg_ent%imax, master_task)
      call broadcast_scalar(ovf(n)%reg_ent%jmin, master_task)
      call broadcast_scalar(ovf(n)%reg_ent%jmax, master_task)
      call broadcast_scalar(ovf(n)%reg_ent%kmin, master_task)
      call broadcast_scalar(ovf(n)%reg_ent%kmax, master_task)
! src locs and orientation
      call broadcast_scalar(ovf(n)%num_src, master_task)
      do m=1,ovf(n)%num_src
         call broadcast_scalar(ovf(n)%loc_src(m)%i, master_task)
         call broadcast_scalar(ovf(n)%loc_src(m)%j, master_task)
         call broadcast_scalar(ovf(n)%loc_src(m)%i_adv, master_task)
         call broadcast_scalar(ovf(n)%loc_src(m)%j_adv, master_task)
         call broadcast_scalar(ovf(n)%loc_src(m)%i_u, master_task)
         call broadcast_scalar(ovf(n)%loc_src(m)%j_u, master_task)
         call broadcast_scalar(ovf(n)%loc_src(m)%k, master_task)
         call broadcast_scalar(ovf(n)%loc_src(m)%orient, master_task)
      end do
! ent locs and orientation
      call broadcast_scalar(ovf(n)%num_ent, master_task)
      do m=1,ovf(n)%num_ent
         call broadcast_scalar(ovf(n)%loc_ent(m)%i, master_task)
         call broadcast_scalar(ovf(n)%loc_ent(m)%j, master_task)
         call broadcast_scalar(ovf(n)%loc_ent(m)%i_adv, master_task)
         call broadcast_scalar(ovf(n)%loc_ent(m)%j_adv, master_task)
         call broadcast_scalar(ovf(n)%loc_ent(m)%i_u, master_task)
         call broadcast_scalar(ovf(n)%loc_ent(m)%j_u, master_task)
         call broadcast_scalar(ovf(n)%loc_ent(m)%k, master_task)
         call broadcast_scalar(ovf(n)%loc_ent(m)%orient, master_task)
      end do
! prd locs and orientation
      call broadcast_scalar(ovf(n)%num_prd_sets, master_task)
      do m=1,ovf(n)%num_prd_sets
         call broadcast_scalar(ovf(n)%num_prd(m), master_task)
         do mp=1,ovf(n)%num_prd(m)
            call broadcast_scalar(ovf(n)%loc_prd(m,mp)%i, master_task)
            call broadcast_scalar(ovf(n)%loc_prd(m,mp)%j, master_task)
            call broadcast_scalar(ovf(n)%loc_prd(m,mp)%i_adv, master_task)
            call broadcast_scalar(ovf(n)%loc_prd(m,mp)%j_adv, master_task)
            call broadcast_scalar(ovf(n)%loc_prd(m,mp)%i_u, master_task)
            call broadcast_scalar(ovf(n)%loc_prd(m,mp)%j_u, master_task)
            call broadcast_scalar(ovf(n)%loc_prd(m,mp)%k, master_task)
            call broadcast_scalar(ovf(n)%loc_prd(m,mp)%orient, master_task)
         end do
      end do
! adjacent boundaries
      call broadcast_scalar(ovf(n)%adj_src%imin, master_task) 
      call broadcast_scalar(ovf(n)%adj_src%imax, master_task) 
      call broadcast_scalar(ovf(n)%adj_src%jmin, master_task) 
      call broadcast_scalar(ovf(n)%adj_src%jmax, master_task) 
      call broadcast_scalar(ovf(n)%adj_src%kmin, master_task) 
      call broadcast_scalar(ovf(n)%adj_src%kmax, master_task) 
      call broadcast_scalar(ovf(n)%adj_ent%imin, master_task) 
      call broadcast_scalar(ovf(n)%adj_ent%imax, master_task) 
      call broadcast_scalar(ovf(n)%adj_ent%jmin, master_task) 
      call broadcast_scalar(ovf(n)%adj_ent%jmax, master_task) 
      call broadcast_scalar(ovf(n)%adj_ent%kmin, master_task) 
      call broadcast_scalar(ovf(n)%adj_ent%kmax, master_task) 
      do m=1,ovf(n)%num_prd_sets
         call broadcast_scalar(ovf(n)%adj_prd(m)%imin, master_task) 
         call broadcast_scalar(ovf(n)%adj_prd(m)%imax, master_task) 
         call broadcast_scalar(ovf(n)%adj_prd(m)%jmin, master_task) 
         call broadcast_scalar(ovf(n)%adj_prd(m)%jmax, master_task) 
         call broadcast_scalar(ovf(n)%adj_prd(m)%kmin, master_task) 
         call broadcast_scalar(ovf(n)%adj_prd(m)%kmax, master_task) 
      end do
   end do  ! ovf broadcast loop

!-----------------------------------------------------------------------
!  initialize overflow data for all processors, so no need to broadcast
!-----------------------------------------------------------------------

   do n=1,num_ovf
     ovf(n)%Ms            = c0
     ovf(n)%Ms_n          = c0
     ovf(n)%Ms_nm1        = c0
     ovf(n)%Me            = c0
     ovf(n)%Me_n          = c0
     ovf(n)%Me_nm1        = c0
     ovf(n)%phi           = c0
     ovf(n)%Mp            = c0
     ovf(n)%Mp_n          = c0
     ovf(n)%Mp_nm1        = c0
     ovf(n)%wght_reg%inf  = c0
     ovf(n)%wght_reg%src  = c0
     ovf(n)%wght_reg%ent  = c0
     do m=1,ovf(n)%num_prd_sets
       ovf(n)%wght_adj%prd(m)  = c0
     end do
     ovf(n)%rho_reg%inf  = c0
     ovf(n)%rho_reg%src  = c0
     ovf(n)%rho_reg%ent  = c0
     do m=1,ovf(n)%num_prd_sets
       ovf(n)%rho_adj%prd(m)  = c0
     end do
     ovf(n)%prd_set_n = 1
     ovf(n)%prd_set   = 1
     do nn=1,nt
       ovf(n)%trcr_reg%inf(nn) = c0
       ovf(n)%trcr_reg%src(nn) = c0
       ovf(n)%trcr_reg%ent(nn) = c0
       ovf(n)%trcr_adj%src(nn) = c0
       ovf(n)%trcr_adj%ent(nn) = c0
       ovf(n)%trcr_adj%prd(nn) = c0
     end do
     do m=1,ovf(n)%num_src
       do k=1,km
         ovf(n)%loc_src(m)%Utlda(k) = c0
         ovf(n)%loc_src(m)%Vtlda(k) = c0
       end do
       ovf(n)%loc_src(m)%Uovf_nm1   = c0
       ovf(n)%loc_src(m)%Uovf_n     = c0
       ovf(n)%loc_src(m)%Uovf       = c0
       ovf(n)%loc_src(m)%Wovf       = c0
     end do
     do m=1,ovf(n)%num_ent
       do k=1,km
         ovf(n)%loc_ent(m)%Utlda(k) = c0
         ovf(n)%loc_ent(m)%Vtlda(k) = c0
       end do
       ovf(n)%loc_ent(m)%Uovf_nm1   = c0
       ovf(n)%loc_ent(m)%Uovf_n     = c0
       ovf(n)%loc_ent(m)%Uovf       = c0
       ovf(n)%loc_ent(m)%Wovf       = c0
     end do
     do m=1,ovf(n)%num_prd_sets
       do mp=1,ovf(n)%num_prd(m)
         do k=1,km
           ovf(n)%loc_prd(m,mp)%Utlda(k) = c0
           ovf(n)%loc_prd(m,mp)%Vtlda(k) = c0
         end do
         ovf(n)%loc_prd(m,mp)%Uovf_nm1   = c0
         ovf(n)%loc_prd(m,mp)%Uovf_n     = c0
         ovf(n)%loc_prd(m,mp)%Uovf       = c0
         ovf(n)%loc_prd(m,mp)%Wovf       = c0
       end do
     end do
   end do  ! ovf initialization loop for all processors

   else if( overflows_restart_type /= 'ccsm_startup' ) then

     call ovf_read_restart
     call ovf_read_broadcast

   endif






!-----------------------------------------------------------------------
!EOC

 end subroutine init_overflows1

!***********************************************************************
!EOP
! !IROUTINE: init_overflows2
! !INTERFACE:

 subroutine init_overflows2

! !DESCRIPTION:
!  This routine continues the initialization of the overflows by 
!  scattering KMT_G to KMT, then modifying if desired, and finally 
!  computing overflow masks.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  scatter KMT_G to KMT if topography_opt = file 
!
!-----------------------------------------------------------------------

   if (registry_match('topography_opt_file')) then
      if (my_task == master_task) write(stdout,'(a30,a)') &
         ' Reading topography from file:', trim(topography_filename)
      call read_topography(topography_filename,.false.)
   endif
   
   if (.not. overflows_on ) return


! don't initialize the groups yet
  ovf_groups%init = .false.


!-----------------------------------------------------------------------
!
!  modify KMT for overflows if desired and ccsm_startup run
!  make kmt changes regardless of overflows_interactive
!
!-----------------------------------------------------------------------

   call init_overflows_kmt

!-----------------------------------------------------------------------
!
!  set overflow masks for regional averaging
!
!-----------------------------------------------------------------------

   call init_overflows_mask

!-----------------------------------------------------------------------
!EOC

 end subroutine init_overflows2

!***********************************************************************
!EOP
! !IROUTINE: init_overflows_kmt
! !INTERFACE:

 subroutine init_overflows_kmt

! !DESCRIPTION:
!  This routine modifies kmt as required by overflows, if on
!  and interactive.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC

!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (int_kind) :: &
      iblock,i,j,k,m,n,  &  ! dummy loop indices
      ib,ie,jb,je,       &  ! local domain index boundaries
      kmterr                ! error index for kmt changes

   type (block) :: &
      this_block         ! block information for current block

!----------------------------------------------------------------------
!
!  search through kmt and modify for overflows
!
!----------------------------------------------------------------------

      kmterr = 0
      do iblock = 1,nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)
         ib = this_block%ib
         ie = this_block%ie
         jb = this_block%jb
         je = this_block%je
         do j=jb,je
            do i=ib,ie
                do n=1,num_ovf
                   do m=1,ovf(n)%num_kmt
                      if( ovf(n)%loc_kmt(m)%i.eq.this_block%i_glob(i).and.&
                          ovf(n)%loc_kmt(m)%j.eq.this_block%j_glob(j) ) then
                         write(stdout,100) KMT(i,j,iblock),ovf(n)%loc_kmt(m)%i, &
                                ovf(n)%loc_kmt(m)%j,ovf(n)%loc_kmt(m)%knew
                         100 format(' init_overflows_kmt: KMT = ',i5,&
                                    ' at global (i,j) = ',2(i5,1x),&
                                    ' changed to ',i5)
                         if( KMT(i,j,iblock) .ne. ovf(n)%loc_kmt(m)%korg ) then
                            kmterr = kmterr + 1
                         endif
                         KMT(i,j,iblock) = ovf(n)%loc_kmt(m)%knew
                      endif
                   end do
                end do
            enddo
         enddo
      enddo
      if (kmterr > 0) then
         write(stdout,200) kmterr
     200 format(' init_overflows_kmt: kmt inconsistencies for ',i3,' points',/ &
                ' original kmt not equal to actual kmt')
         call shr_sys_flush(stdout)
         call exit_POP(sigAbort,'ERROR kmt inconsistency for overflows')
      endif
      call shr_sys_flush(stdout)

   call POP_HaloUpdate(KMT, POP_haloClinic, POP_gridHorzLocCenter,  &
                       POP_fieldKindScalar, errorCode,              &
                       fillValue = 0_POP_i4)

!----------------------------------------------------------------------
!EOC

 end subroutine init_overflows_kmt

!***********************************************************************
!EOP
! !IROUTINE: init_overflows_mask
! !INTERFACE:

 subroutine init_overflows_mask

! !DESCRIPTION:
!  This routine sets overflow masks for regional and adjacent averaging
!   Need to do this on all procs
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC

!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (int_kind) :: &
      iblock,i,j,n,m,    &  ! dummy loop indices
      ib,ie,jb,je           ! local domain index boundaries

   type (block) :: &
      this_block         ! block information for current block

!----------------------------------------------------------------------
!
!  set masks for regional averaging
!
!----------------------------------------------------------------------

!must go through ALL the overflows
   do n=1,num_ovf
      do iblock = 1,nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)
         ib = this_block%ib
         ie = this_block%ie
         jb = this_block%jb
         je = this_block%je
         do j=jb,je
! inflow region
            if( ovf(n)%reg_inf%jmin  .le. this_block%j_glob(j) .and. &
                this_block%j_glob(j) .le. ovf(n)%reg_inf%jmax ) then
               do i=ib,ie
                  ovf(n)%mask_reg%inf(i,j,iblock) = c0
                  if( ovf(n)%reg_inf%imin  .le. this_block%i_glob(i) .and. &
                      this_block%i_glob(i) .le. ovf(n)%reg_inf%imax ) then
                         ovf(n)%mask_reg%inf(i,j,iblock) = c1
                         write(stdout,30) ovf(n)%name,this_block%i_glob(i), &
                                           this_block%j_glob(j)
                         30 format(' Overflow: ',a24, &
                            ' Inflow region mask at global (ij)=',2(i3,2x))
                  endif
               end do
            endif     ! inflow region
! source region
            if( ovf(n)%reg_src%jmin  .le. this_block%j_glob(j) .and. &
                this_block%j_glob(j) .le. ovf(n)%reg_src%jmax ) then
               do i=ib,ie
                  ovf(n)%mask_reg%src(i,j,iblock) = c0
                  if( ovf(n)%reg_src%imin  .le. this_block%i_glob(i) .and. &
                      this_block%i_glob(i) .le. ovf(n)%reg_src%imax ) then
                         ovf(n)%mask_reg%src(i,j,iblock) = c1
                         write(stdout,31) ovf(n)%name,this_block%i_glob(i), &
                                           this_block%j_glob(j)
                         31 format(' Overflow: ',a24, &
                            ' Source region mask at global (ij)=',2(i3,2x))
                  endif
               end do
            endif     ! source region
! source adjacent
            if( ovf(n)%adj_src%jmin  .le. this_block%j_glob(j) .and. &
                this_block%j_glob(j) .le. ovf(n)%adj_src%jmax ) then
               do i=ib,ie
                  ovf(n)%mask_adj%src(i,j,iblock) = c0
                  if( ovf(n)%adj_src%imin  .le. this_block%i_glob(i) .and. &
                      this_block%i_glob(i) .le. ovf(n)%adj_src%imax ) then
                         ovf(n)%mask_adj%src(i,j,iblock) = c1
                         write(stdout,32) ovf(n)%name,this_block%i_glob(i), &
                                           this_block%j_glob(j)
                         32 format(' Overflow: ',a24, &
                            ' Source adjacent mask at global (ij)=',2(i3,2x))
                  endif
               end do
            endif     ! source adjacent
! entrainment region
            if( ovf(n)%reg_ent%jmin  .le. this_block%j_glob(j) .and. &
                this_block%j_glob(j) .le. ovf(n)%reg_ent%jmax ) then
               do i=ib,ie
                  ovf(n)%mask_reg%ent(i,j,iblock) = c0
                  if( ovf(n)%reg_ent%imin  .le. this_block%i_glob(i) .and. &
                      this_block%i_glob(i) .le. ovf(n)%reg_ent%imax ) then
                         ovf(n)%mask_reg%ent(i,j,iblock) = c1
                         write(stdout,33) ovf(n)%name,this_block%i_glob(i), &
                                           this_block%j_glob(j)
                         33 format(' Overflow: ',a24, &
                            ' Entrainment region mask at global (ij)=',2(i3,2x))
                  endif
               end do
            endif     ! entrainment region
! entrainment adjacent
            if( ovf(n)%adj_ent%jmin  .le. this_block%j_glob(j) .and. &
                this_block%j_glob(j) .le. ovf(n)%adj_ent%jmax ) then
               do i=ib,ie
                  ovf(n)%mask_adj%ent(i,j,iblock) = c0
                  if( ovf(n)%adj_ent%imin  .le. this_block%i_glob(i) .and. &
                      this_block%i_glob(i) .le. ovf(n)%adj_ent%imax ) then
                         ovf(n)%mask_adj%ent(i,j,iblock) = c1
                         write(stdout,34) ovf(n)%name,this_block%i_glob(i), &
                                           this_block%j_glob(j)
                         34 format(' Overflow: ',a24, &
                            ' Entrainment adjacent mask at global (ij)=',2(i3,2x))
                  endif
               end do
            endif     ! entrainment adjacent
         end do  ! j loop
! product adjacent
         do m=1,ovf(n)%num_prd_sets
            do j=jb,je
               if( ovf(n)%adj_prd(m)%jmin  .le. this_block%j_glob(j) .and. &
                   this_block%j_glob(j) .le. ovf(n)%adj_prd(m)%jmax ) then
                  do i=ib,ie
                     ovf(n)%mask_adj%prd(i,j,iblock,m) = c0
                     if( ovf(n)%adj_prd(m)%imin  .le. this_block%i_glob(i) .and. &
                         this_block%i_glob(i) .le. ovf(n)%adj_prd(m)%imax ) then
                            ovf(n)%mask_adj%prd(i,j,iblock,m) = c1
                            write(stdout,35) ovf(n)%name,this_block%i_glob(i), &
                                              this_block%j_glob(j)
                            35 format(' Overflow: ',a24, &
                               ' Product adjacent mask at global (ij)=',2(i3,2x))
                     endif
                  end do
               endif     ! product adjacent
            end do
         end do
      end do
   end do
   call shr_sys_flush(stdout)

!----------------------------------------------------------------------
!EOC

 end subroutine init_overflows_mask

!***********************************************************************
!EOP
! !IROUTINE: init_overflows3
! !INTERFACE:

 subroutine init_overflows3

! !DESCRIPTION:
!  This routine completes the initialization of the overflows by 
!  modifying the 9pt coefficients for the barotropic solution
!  as required for each overflow grid box
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  modify 9pt coefficients for barotropic solver
!
!-----------------------------------------------------------------------

   if( overflows_on .and. overflows_interactive ) then
      call ovf_solvers_9pt
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine init_overflows3

!***********************************************************************
!EOP
! !IROUTINE: init_overflows4
! !INTERFACE:

 subroutine init_overflows4

! !DESCRIPTION:
!  This routine creates the overflow output diagnostics filename, now
!  that the initial model run time is known.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------
   character (char_len) ::  &
      string

   save

!-----------------------------------------------------------------------
!  if overflows off, exit
!-----------------------------------------------------------------------

   if (.not. overflows_on ) return


! don't initialize the groups yet
!  ovf_groups%init = .false.


!-----------------------------------------------------------------------
!  set up output file and unit for overflow diagnostics
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  define ccsm overflow diagnostics output filename
!-----------------------------------------------------------------------
   if (lccsm) then
     call ccsm_date_stamp (ccsm_diag_date, 'ymds')
     string = overflows_diag_outfile
     overflows_diag_outfile = trim(string)/&
                                           &/'.'/&
                                           &/trim(ccsm_diag_date)
   else
!-----------------------------------------------------------------------
!  append runid, initial date to output file names
!  concatenation operator must be split across lines to avoid problems
!    with preprocessors
!-----------------------------------------------------------------------
     if (date_separator == ' ') then
        cdate(1:4) = cyear
        cdate(5:6) = cmonth
        cdate(7:8) = cday
        cdate(9:10)= '  '
     else
        cdate(1:4) = cyear
        cdate(5:5) = date_separator
        cdate(6:7) = cmonth
        cdate(8:8) = date_separator
        cdate(9:10) = cday
     endif
     outfile_tmp = char_blank
     outfile_tmp = trim(overflows_diag_outfile)/&
                                          &/'.'/&
                                          &/trim(runid)/&
                                          &/'.'/&
                                          &/trim(cdate)
     overflows_diag_outfile = trim(outfile_tmp)
   endif ! lccsm


   call get_unit(ovf_diag_unit)
   if (my_task == master_task) then
       open(ovf_diag_unit, file=overflows_diag_outfile, status='unknown')
       write(ovf_diag_unit,*)' '
       close(ovf_diag_unit)

       write(stdout,'(a,a)') &
          'Overflow diagnostics written to file: ', trim(overflows_diag_outfile)
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine init_overflows4

!***********************************************************************
!EOP
! !IROUTINE: init_overflows5
! !INTERFACE:

 subroutine init_overflows5

! !DESCRIPTION:
!  This routine computes regional aveages required at restart
!  for overflow regions, using all available tracers, and also
!  computes regional product values based on source and entrainment.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------
   save

   integer   (int_kind)  :: &
      n         ,& ! index of overflow
     nn            ! ovf tracer index
   real (r8) ::  &
     phi           ! entrainment parameter from actual ratio Me/Mp

   integer (int_kind) ovf_id



!----------------------
!initialize ovf timer (not calling in block loops)
!------------------------

   call get_timer(timer_ovf,'OVF TIMER', 1, &
        distrb_clinic%nprocs)


!-----------------------------------------------------------------------
!
!  compute regional averages.
!
!-----------------------------------------------------------------------

   if( overflows_on .and. overflows_interactive ) then
      call ovf_reg_avgs(oldtime)
! evaluate regional product values based on src,ent averages just computed
! src

  do n=1,ovf_groups%numMyGroups
     ovf_id = ovf_groups%groupIds(n)
        phi = ovf(ovf_id)%phi
        do nn=1,nt
           ovf(ovf_id)%trcr_reg%prd(nn) = ovf(ovf_id)%trcr_reg%src(nn) * (c1 - phi) &
                                   + ovf(ovf_id)%trcr_reg%ent(nn) * phi
        end do
      enddo
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine init_overflows5

!***********************************************************************
!EOP
! !IROUTINE: ovf_write_restart
! !INTERFACE:

 subroutine ovf_write_restart


! !DESCRIPTION:
!  This routine writes the overflow restart file using 
!  selected data from overflow array.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------
   save

   integer   (int_kind)  :: &
      mu,                   &! unit for ovf restart file
      ovf_error,            &! error flag
      n,                    &! ovf loop index
      m,                    &! sub-ovf loop index
      nn,                   &! tracer loop index 
      mp,                   &! sub-ovf sub-loop index
      num_posts,            &
      ierr, count, loc, ovf_id, i

   logical (log_kind), dimension(:), allocatable :: post_array
   integer(i4), dimension(:), allocatable :: myRequests(:)
   integer (i4), dimension(:,:), allocatable:: myStati(:,:)

   integer (int_kind), parameter :: len = 14 !length of buffer for a post

   real (r8), dimension(:), allocatable:: myRecvBuff(:), mySendBuff(:)


   character (char_len)  ::     &
      write_restart_filename,   &! modified file name for restart file
      ovf_restart_pointer_file, &! overflows rpointer filename
      file_suffix,              &! suffix to append to root filename
      char_temp                  ! temporary character string

   character (10)        :: &! for input year,month,day
      cdate

   logical (log_kind) :: post

   logical (log_kind), parameter :: prnt = .false.

   if( prnt .and. my_task == master_task ) then
      write(stdout,*) 'ovf_write_restart called '
      call shr_sys_flush(stdout)
   endif

!-----------------------------------------------------------------------
!  if overflows off, exit
!-----------------------------------------------------------------------

   if( .not. overflows_on ) return

   cdate(1:4) = cyear
   cdate(5:6) = cmonth
   cdate(7:8) = cday
   cdate(9:10)= '  '

   ovf_error = 0
   call get_unit(mu)

   write_restart_filename = char_blank
   file_suffix = char_blank

   if (registry_match('lccsm')) then
     call ccsm_date_stamp(char_temp, 'ymds')
     file_suffix = trim(char_temp)
     !*** must split concatenation operator to avoid preprocessor mangling
     write_restart_filename = trim(overflows_restfile)/&
                                                    &/'.'/&
                                                    &/trim(file_suffix)
   else
     write_restart_filename = trim(overflows_restfile)
   endif


!-----------------------------------------------------------------------
!  master task section
!-----------------------------------------------------------------------

   if (my_task == master_task) then


      !for any ovf groups that the master_task is not in, 
      !must post recv for the info
      allocate(post_array(num_ovf))
      num_posts = num_ovf - ovf_groups%numMyGroups
      if (num_posts > 0) then
         allocate(myRecvBuff(len*num_posts))
      end if

      call ovf_print_init(len, num_posts, myRecvBuff, post_array)

      open(mu, file=write_restart_filename, status='unknown',iostat=ovf_error)

      write(stdout,987) mu,write_restart_filename
987   format(' ovf_write_restart  unit (mu) = ',i5,' file name = ',a64)
      
      write(stdout,99) cdate
99    format(' ovf write restart   cdate yyyymmdd = ',a10)
      call shr_sys_flush(stdout)
      
      write(mu,100) cdate,num_ovf
100   format(30x,'  ! Overflow Restart File for yyyymmdd =',a10/ &
           2x,i10,20x,'! number of overflows')  


      !MAIN ovf loop
      count = 0
      do n=1,num_ovf

         write(mu,101) ovf(n)%name
101      format(2x,a26,'    ! name of overflow')
         
         ! ovf parameters
         write(mu,102) ovf(n)%ovf_params%lat
102      format(2x,1PE27.18,'   ! latitude in degrees')
         write(mu,103) ovf(n)%ovf_params%width
103      format(2x,1PE27.18,'   ! channel width in meters')
         write(mu,105) ovf(n)%ovf_params%source_thick
105      format(2x,1PE27.18,'   ! source thickness in meters')
         write(mu,106) ovf(n)%ovf_params%distnc_str_ssb
106      format(2x,1PE27.18,'   ! strait to shelf-slope break in meters')
         write(mu,107) ovf(n)%ovf_params%bottom_slope
107      format(2x,1PE27.18,'   ! bottom slope dy/dx ')
         write(mu,108) ovf(n)%ovf_params%bottom_drag
108      format(2x,1PE27.18,'   ! bottom drag coefficient')
         
         ! kmt changes, if any
         write(mu,1090) ovf(n)%num_kmt
1090     format(2x,i10,20x,'! number of kmt changes')
         do m=1,ovf(n)%num_kmt
            write(mu,1091) ovf(n)%loc_kmt(m)%i
1091        format(2x,i10,20x,'! i grid box index for kmt change')
            write(mu,1092) ovf(n)%loc_kmt(m)%j
1092        format(2x,i10,20x,'! j grid box index for kmt change')
            write(mu,1093) ovf(n)%loc_kmt(m)%korg
1093        format(2x,i10,20x,'! korg  original grid box k index')
            write(mu,1094) ovf(n)%loc_kmt(m)%knew
1094        format(2x,i10,20x,'! knew  new      grid box k index')
         end do
         
         ! regional boundaries
         !   inflow
         write(mu,110) ovf(n)%reg_inf%imin
110      format(2x,i10,20x,'! inflow region imin')
         write(mu,111) ovf(n)%reg_inf%imax
111      format(2x,i10,20x,'! inflow region imax')
         write(mu,112) ovf(n)%reg_inf%jmin
112      format(2x,i10,20x,'! inflow region jmin')
         write(mu,113) ovf(n)%reg_inf%jmax
113      format(2x,i10,20x,'! inflow region jmax')
         write(mu,114) ovf(n)%reg_inf%kmin
114      format(2x,i10,20x,'! inflow region kmin')
         write(mu,115) ovf(n)%reg_inf%kmax
115      format(2x,i10,20x,'! inflow region kmax')
         !   source
         write(mu,116) ovf(n)%reg_src%imin
116      format(2x,i10,20x,'! source region imin')
         write(mu,117) ovf(n)%reg_src%imax
117      format(2x,i10,20x,'! source region imax')
         write(mu,118) ovf(n)%reg_src%jmin
118      format(2x,i10,20x,'! source region jmin')
         write(mu,119) ovf(n)%reg_src%jmax
119      format(2x,i10,20x,'! source region jmax')
         write(mu,120) ovf(n)%reg_src%kmin
120      format(2x,i10,20x,'! source region kmin')
         write(mu,121) ovf(n)%reg_src%kmax
121      format(2x,i10,20x,'! source region kmax')
         !   entrainment
         write(mu,122) ovf(n)%reg_ent%imin
122      format(2x,i10,20x,'! entrainment region imin')
         write(mu,123) ovf(n)%reg_ent%imax
123      format(2x,i10,20x,'! entrainment region imax')
         write(mu,124) ovf(n)%reg_ent%jmin
124      format(2x,i10,20x,'! entrainment region jmin')
         write(mu,125) ovf(n)%reg_ent%jmax
125      format(2x,i10,20x,'! entrainment region jmax')
         write(mu,126) ovf(n)%reg_ent%kmin
126      format(2x,i10,20x,'! entrainment region kmin')
         write(mu,127) ovf(n)%reg_ent%kmax
127      format(2x,i10,20x,'! entrainment region kmax')
         ! src locs and orientation
         write(mu,128) ovf(n)%num_src
128      format(2x,i10,20x,'! number of source grid boxes')
         do m=1,ovf(n)%num_src
            write(mu,129) ovf(n)%loc_src(m)%i
129         format(2x,i10,20x,'! source box i')
            write(mu,130) ovf(n)%loc_src(m)%j
130         format(2x,i10,20x,'! source box j')
            write(mu,131) ovf(n)%loc_src(m)%i_adv
131         format(2x,i10,20x,'! source box i_adv')
            write(mu,132) ovf(n)%loc_src(m)%j_adv
132         format(2x,i10,20x,'! source box j_adv')
            write(mu,133) ovf(n)%loc_src(m)%i_u
133         format(2x,i10,20x,'! source box i_u')
            write(mu,134) ovf(n)%loc_src(m)%j_u
134         format(2x,i10,20x,'! source box j_u')
            write(mu,135) ovf(n)%loc_src(m)%k
135         format(2x,i10,20x,'! source box k')
            write(mu,136) ovf(n)%loc_src(m)%orient
136         format(2x,i10,20x,'! source box orient')
         end do
         ! ent locs and orientation
         write(mu,137) ovf(n)%num_ent
137      format(2x,i10,20x,'! number of entrainment grid boxes')
         do m=1,ovf(n)%num_ent
            write(mu,138) ovf(n)%loc_ent(m)%i
138         format(2x,i10,20x,'! entrainment box i')
            write(mu,139) ovf(n)%loc_ent(m)%j
139         format(2x,i10,20x,'! entrainment box j')
            write(mu,140) ovf(n)%loc_ent(m)%i_adv
140         format(2x,i10,20x,'! entrainment box i_adv')
            write(mu,141) ovf(n)%loc_ent(m)%j_adv
141         format(2x,i10,20x,'! entrainment box j_adv')
            write(mu,142) ovf(n)%loc_ent(m)%i_u
142         format(2x,i10,20x,'! entrainment box i_u')
            write(mu,143) ovf(n)%loc_ent(m)%j_u
143         format(2x,i10,20x,'! entrainment box j_u')
            write(mu,144) ovf(n)%loc_ent(m)%k
144         format(2x,i10,20x,'! entrainment box k')
            write(mu,145) ovf(n)%loc_ent(m)%orient
145         format(2x,i10,20x,'! entrainment box orient')
         end do
         ! prd locs and orientation
         write(mu,146) ovf(n)%num_prd_sets
146      format(2x,i10,20x,'! number of product sets')
         do m=1,ovf(n)%num_prd_sets
            write(mu,147) ovf(n)%num_prd(m)
147         format(2x,i10,20x, &
                 '! number of product grid boxes for this set')
            do mp=1,ovf(n)%num_prd(m)
               write(mu,148) ovf(n)%loc_prd(m,mp)%i
148            format(2x,i10,20x,'! product box i')
               write(mu,149) ovf(n)%loc_prd(m,mp)%j
149            format(2x,i10,20x,'! product box j')
               write(mu,150) ovf(n)%loc_prd(m,mp)%i_adv
150            format(2x,i10,20x,'! product box i_adv')
               write(mu,151) ovf(n)%loc_prd(m,mp)%j_adv
151            format(2x,i10,20x,'! product box j_adv')
               write(mu,152) ovf(n)%loc_prd(m,mp)%i_u
152            format(2x,i10,20x,'! product box i_u')
               write(mu,153) ovf(n)%loc_prd(m,mp)%j_u
153            format(2x,i10,20x,'! product box j_u')
               write(mu,154) ovf(n)%loc_prd(m,mp)%k
154            format(2x,i10,20x,'! product box k')
               write(mu,155) ovf(n)%loc_prd(m,mp)%orient
155            format(2x,i10,20x,'! product box orient')
          end do
       end do
! adjacent boundaries
       ! src
       write(mu,156) ovf(n)%adj_src%imin 
156    format(2x,i10,20x,'! source adjacent imin')
       write(mu,157) ovf(n)%adj_src%imax 
157    format(2x,i10,20x,'! source adjacent imax')
       write(mu,158) ovf(n)%adj_src%jmin 
158    format(2x,i10,20x,'! source adjacent jmin')
       write(mu,159) ovf(n)%adj_src%jmax 
159    format(2x,i10,20x,'! source adjacent jmax')
       write(mu,160) ovf(n)%adj_src%kmin 
160    format(2x,i10,20x,'! source adjacent kmin')
       write(mu,161) ovf(n)%adj_src%kmax 
161    format(2x,i10,20x,'! source adjacent kmax')
       !ent
       write(mu,162) ovf(n)%adj_ent%imin 
162    format(2x,i10,20x,'! entrainment adjacent imin')
       write(mu,163) ovf(n)%adj_ent%imax 
163    format(2x,i10,20x,'! entrainment adjacent imax')
       write(mu,164) ovf(n)%adj_ent%jmin 
164    format(2x,i10,20x,'! entrainment adjacent jmin')
       write(mu,165) ovf(n)%adj_ent%jmax 
165    format(2x,i10,20x,'! entrainment adjacent jmax')
       write(mu,166) ovf(n)%adj_ent%kmin 
166    format(2x,i10,20x,'! entrainment adjacent kmin')
       write(mu,167) ovf(n)%adj_ent%kmax 
167    format(2x,i10,20x,'! entrainment adjacent kmax')
       !prd
       do m=1,ovf(n)%num_prd_sets
          write(mu,168) ovf(n)%adj_prd(m)%imin 
168       format(2x,i10,20x,'! product adjacent imin')
          write(mu,169) ovf(n)%adj_prd(m)%imax 
169       format(2x,i10,20x,'! product adjacent imax')
          write(mu,170) ovf(n)%adj_prd(m)%jmin 
170       format(2x,i10,20x,'! product adjacent jmin')
          write(mu,171) ovf(n)%adj_prd(m)%jmax 
171       format(2x,i10,20x,'! product adjacent jmax')
          write(mu,172) ovf(n)%adj_prd(m)%kmin 
172       format(2x,i10,20x,'! product adjacent kmin')
          write(mu,173) ovf(n)%adj_prd(m)%kmax 
173       format(2x,i10,20x,'! product adjacent kmax')
       end do
       
       
       ! transports
       !! this information *could* be on other procs that are not the master_task

       if (post_array(n)) then !we posted a recv for this
          count = count+1
          call ovf_print_get(count)
          loc = (count-1)*len+1
          !now unpack
          ovf(n)%Ms = myRecvBuff(loc)
          ovf(n)%Ms_n = myRecvBuff(loc+1)
          ovf(n)%Ms_nm1 = myRecvBuff(loc+2)
          ovf(n)%Me = myRecvBuff(loc+3)
          ovf(n)%Me_n = myRecvBuff(loc+4)
          ovf(n)%Me_nm1 = myRecvBuff(loc+5)
          ovf(n)%phi = myRecvBuff(loc+6)
          ovf(n)%Mp = myRecvBuff(loc+7)
          ovf(n)%Mp_n = myRecvBuff(loc+8)
          ovf(n)%Mp_nm1 = myRecvBuff(loc+9)
          ovf(n)%Tp = myRecvBuff(loc+10)
          ovf(n)%Sp = myRecvBuff(loc+11)
          ovf(n)%prd_set_n = int(myRecvBuff(loc+12), int_kind)
          ovf(n)%prd_set = int(myRecvBuff(loc+13), int_kind)

!          print *, 'CHECK prd_set IAM = ', my_task, 'ovf_id=', n, &
!               'prd_set_n, prd_set', &
!               ovf(n)%prd_set_n, " , ", ovf(n)%prd_set


       end if

       write(mu,174) ovf(n)%Ms
174    format(2x,1PE27.18,'   ! source volume n+1 transport cm3/sec')
       write(mu,175) ovf(n)%Ms_n
175    format(2x,1PE27.18,'   ! source volume n   transport cm3/sec')
       write(mu,176) ovf(n)%Ms_nm1
176    format(2x,1PE27.18,'   ! source volume n-1 transport cm3/sec')
       write(mu,177) ovf(n)%Me
177    format(2x,1PE27.18,'   ! entrainment volume n+1 transport cm3/sec')
       write(mu,178) ovf(n)%Me_n
178    format(2x,1PE27.18,'   ! entrainment volume n   transport cm3/sec')
       write(mu,179) ovf(n)%Me_nm1
179    format(2x,1PE27.18,'   ! entrainment volume n-1 transport cm3/sec')
       write(mu,180) ovf(n)%phi
180    format(2x,1PE27.18,'   ! phi parameter')
       write(mu,181) ovf(n)%Mp
181    format(2x,1PE27.18,'   ! product volume n+1 transport cm3/sec')
       write(mu,182) ovf(n)%Mp_n
182    format(2x,1PE27.18,'   ! product volume n   transport cm3/sec')
       write(mu,183) ovf(n)%Mp_nm1
183    format(2x,1PE27.18,'   ! product volume n-1 transport cm3/sec')
       write(mu,184) ovf(n)%Tp
184    format(2x,1PE27.18,'   ! product temperature C')
       write(mu,185) ovf(n)%Sp
185    format(2x,1PE27.18,'   ! product salinity')
       write(mu,186) ovf(n)%prd_set_n
       write(mu,186) ovf(n)%prd_set
186    format(2x,i10,20x,'! product set index (first is previous time step)')
       
    end do  ! ovf loop 
    
    close(mu)

    ! clean up
    deallocate(post_array)
    if (num_posts > 0) then
       deallocate(myRecvBuff)
    endif
    call ovf_print_finalize(num_posts)

 else    ! my_task \= master_task
    count = 0
    do n=1, ovf_groups%numMyGroups
       if (ovf_groups%amMaster(n)) then
          count = count + 1
       end if
    end do

    call ovf_print_init(len, count)

    if (count > 0) then
       allocate(mySendBuff(count*len))
    end if
    
    num_posts = 0
    do n=1, ovf_groups%numMyGroups
       if (ovf_groups%amMaster(n)) then
          ovf_id = ovf_groups%groupIds(n)
          !send transport info to the master task with isend
          num_posts = num_posts + 1
          !pack send buffer
          loc = (num_posts-1)*len + 1
          mySendBuff(loc) = ovf(ovf_id)%Ms
          mySendBuff(loc+1) = ovf(ovf_id)%Ms_n
          mySendBuff(loc+2) = ovf(ovf_id)%Ms_nm1
          mySendBuff(loc+3) = ovf(ovf_id)%Me
          mySendBuff(loc+4) = ovf(ovf_id)%Me_n
          mySendBuff(loc+5) = ovf(ovf_id)%Me_nm1
          mySendBuff(loc+6) = ovf(ovf_id)%phi
          mySendBuff(loc+7) = ovf(ovf_id)%Mp
          mySendBuff(loc+8) = ovf(ovf_id)%Mp_n
          mySendBuff(loc+9) = ovf(ovf_id)%Mp_nm1
          mySendBuff(loc+10) = ovf(ovf_id)%Tp
          mySendBuff(loc+11) = ovf(ovf_id)%Sp
          mySendBuff(loc+12) = real(ovf(ovf_id)%prd_set_n, r8)
          mySendBuff(loc+13) = real(ovf(ovf_id)%prd_set, r8)

!          print *, 'CHECK prd_set IAM = ', my_task, 'ovf_id=', ovf_id, &
!               'prd_set_n, prd_set', &
!               ovf(ovf_id)%prd_set_n, " , ", ovf(ovf_id)%prd_set

          call ovf_print_send(len, mySendBuff(loc:loc+len-1), num_posts, ovf_id)
       end if
    end do ! n loop
    
    !now clean up
    call ovf_print_finalize(num_posts)
    if (num_posts > 0) then
       deallocate(mySendBuff)
    endif

 endif

call release_unit(mu)
 
 !-----------------------------------------------------------------------
 !
!  if pointer files are used, write filename to pointer file
!
!-----------------------------------------------------------------------

   if (luse_pointer_files) then
      call get_unit(mu)
      if (my_task == master_task) then
       ovf_restart_pointer_file = trim(pointer_filename)/&
                                                       &/'.ovf'
       open(mu, file=ovf_restart_pointer_file, form='formatted', status='unknown')
       write(mu,'(a)') trim(write_restart_filename)
       close(mu)
       write(stdout,blank_fmt)
       write(stdout,*) ' overflow restart pointer file written: ',trim(ovf_restart_pointer_file)
       call POP_IOUnitsFlush(POP_stdout)
     endif
     call release_unit(mu)
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine ovf_write_restart

!***********************************************************************
!EOP
! !IROUTINE: ovf_read_restart
! !INTERFACE:

 subroutine ovf_read_restart

! !DESCRIPTION:
!  This routine reads the overflow restart file for
!  selected data from overflow array.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------
   save

   integer (POP_i4)      :: &
      mu,                   &! unit for ovf restart file
      ovf_error,            &! error flag
      n,                    &! ovf loop index
      m,                    &! sub-ovf loop index
      nn,                   &! tracer loop index 
      mp,                   &! sub-ovf sub-loop index
      ntrcr,                &! number of tracers on read
      cindx,cindx2           ! indices into restart pointer character string

   character (POP_charLength) ::  &
      restart_pointer_file,    &! file name for restart pointer file
      read_overflows_restfile, &! local restart filename
      cdate_label               ! for input year,month,day

   logical (POP_logical), parameter :: prnt = .false.


   if( prnt .and. my_task == master_task ) then
      write(stdout,*) 'ovf_read_restart called '
      call shr_sys_flush(stdout)
   endif

!-----------------------------------------------------------------------
!  if overflows off, exit
!-----------------------------------------------------------------------

   if( .not. overflows_on ) return

   ovf_error = 0

!-----------------------------------------------------------------------
!
!  if pointer files are used, overflows pointer file must be read to get
!  actual filenames - skip this for ccsm_branch initialization
!
!  otherwise use input filename
!-----------------------------------------------------------------------

   errorCode = POP_Success

   read_overflows_restfile = char_blank
   restart_pointer_file = char_blank

   if (luse_pointer_files) then
      call get_unit(mu)
      if (my_task == master_task) then
        restart_pointer_file = pointer_filename
        cindx = len_trim(pointer_filename) + 1
        cindx2= cindx + 3
        restart_pointer_file(cindx:cindx2) = '.ovf'
        write(stdout,*) 'Reading overflow pointer file: ', trim(restart_pointer_file)
        call POP_IOUnitsFlush(POP_stdout)
        open(mu, file=trim(restart_pointer_file), form='formatted', status='old')
        read(mu,'(a)') read_overflows_restfile
        close(mu)
      endif
      call release_unit(mu)
      call broadcast_scalar(read_overflows_restfile, master_task)
   else
      ! use overflows_restfile from namelist 
      read_overflows_restfile = trim(overflows_restfile)
   endif

!-----------------------------------------------------------------------
!  read overflows restart file
!-----------------------------------------------------------------------

   call get_unit(mu)
   if (my_task == master_task) then

    open(mu, file=read_overflows_restfile, status='unknown',iostat=ovf_error)

    write(stdout,987) mu,read_overflows_restfile
    987 format(' ovf_read_restart  unit (mu) = ',i5,' file name = ',a)

    read(mu,99) cdate_label,num_ovf
    99 format(a80/2x,i10)
    write(stdout,100) cdate_label,num_ovf
    100 format(' ovf read restart  label =',/a80/ &
               ' number of overflows = ',i5)
    call shr_sys_flush(stdout)

    do n=1,num_ovf
      read(mu,101) ovf(n)%name
      101 format(2x,a26)

! ovf parameters
      read(mu,102) ovf(n)%ovf_params%lat
      102 format(2x,1PE27.18)
      read(mu,103) ovf(n)%ovf_params%width
      103 format(2x,1PE27.18)
      read(mu,105) ovf(n)%ovf_params%source_thick
      105 format(2x,1PE27.18)
      read(mu,106) ovf(n)%ovf_params%distnc_str_ssb
      106 format(2x,1PE27.18)
      read(mu,107) ovf(n)%ovf_params%bottom_slope
      107 format(2x,1PE27.18)
      read(mu,108) ovf(n)%ovf_params%bottom_drag
      108 format(2x,1PE27.18)
! kmt changes, if any
! GFORTRAN Compiler complains about constants in read format
      read(mu,1090) ovf(n)%num_kmt
1090  format(2x,i10)
      do m=1,ovf(n)%num_kmt
         read(mu,1090) ovf(n)%loc_kmt(m)%i
         read(mu,1090) ovf(n)%loc_kmt(m)%j
         read(mu,1090) ovf(n)%loc_kmt(m)%korg
         read(mu,1090) ovf(n)%loc_kmt(m)%knew
      end do

! regional boundaries
!   inflow
      read(mu,110) ovf(n)%reg_inf%imin
      110 format(2x,i10,20x)
      read(mu,111) ovf(n)%reg_inf%imax
      111 format(2x,i10,20x)
      read(mu,112) ovf(n)%reg_inf%jmin
      112 format(2x,i10,20x)
      read(mu,113) ovf(n)%reg_inf%jmax
      113 format(2x,i10,20x)
      read(mu,114) ovf(n)%reg_inf%kmin
      114 format(2x,i10,20x)
      read(mu,115) ovf(n)%reg_inf%kmax
      115 format(2x,i10,20x)
!   source
      read(mu,116) ovf(n)%reg_src%imin
      116 format(2x,i10,20x)
      read(mu,117) ovf(n)%reg_src%imax
      117 format(2x,i10,20x)
      read(mu,118) ovf(n)%reg_src%jmin
      118 format(2x,i10,20x)
      read(mu,119) ovf(n)%reg_src%jmax
      119 format(2x,i10,20x)
      read(mu,120) ovf(n)%reg_src%kmin
      120 format(2x,i10,20x)
      read(mu,121) ovf(n)%reg_src%kmax
      121 format(2x,i10,20x)
!   entrainment
      read(mu,122) ovf(n)%reg_ent%imin
      122 format(2x,i10,20x)
      read(mu,123) ovf(n)%reg_ent%imax
      123 format(2x,i10,20x)
      read(mu,124) ovf(n)%reg_ent%jmin
      124 format(2x,i10,20x)
      read(mu,125) ovf(n)%reg_ent%jmax
      125 format(2x,i10,20x)
      read(mu,126) ovf(n)%reg_ent%kmin
      126 format(2x,i10,20x)
      read(mu,127) ovf(n)%reg_ent%kmax
      127 format(2x,i10,20x)
! src locs and orientation
      read(mu,128) ovf(n)%num_src
      128 format(2x,i10,20x)
      do m=1,ovf(n)%num_src
         read(mu,129) ovf(n)%loc_src(m)%i
         129 format(2x,i10,20x)
         read(mu,130) ovf(n)%loc_src(m)%j
         130 format(2x,i10,20x)
         read(mu,131) ovf(n)%loc_src(m)%i_adv
         131 format(2x,i10,20x)
         read(mu,132) ovf(n)%loc_src(m)%j_adv
         132 format(2x,i10,20x)
         read(mu,133) ovf(n)%loc_src(m)%i_u
         133 format(2x,i10,20x)
         read(mu,134) ovf(n)%loc_src(m)%j_u
         134 format(2x,i10,20x)
         read(mu,135) ovf(n)%loc_src(m)%k
         135 format(2x,i10,20x)
         read(mu,136) ovf(n)%loc_src(m)%orient
         136 format(2x,i10,20x)
      end do
! ent locs and orientation
      read(mu,137) ovf(n)%num_ent
      137 format(2x,i10,20x)
      do m=1,ovf(n)%num_ent
         read(mu,138) ovf(n)%loc_ent(m)%i
         138 format(2x,i10,20x)
         read(mu,139) ovf(n)%loc_ent(m)%j
         139 format(2x,i10,20x)
         read(mu,140) ovf(n)%loc_ent(m)%i_adv
         140 format(2x,i10,20x)
         read(mu,141) ovf(n)%loc_ent(m)%j_adv
         141 format(2x,i10,20x)
         read(mu,142) ovf(n)%loc_ent(m)%i_u
         142 format(2x,i10,20x)
         read(mu,143) ovf(n)%loc_ent(m)%j_u
         143 format(2x,i10,20x)
         read(mu,144) ovf(n)%loc_ent(m)%k
         144 format(2x,i10,20x)
         read(mu,145) ovf(n)%loc_ent(m)%orient
         145 format(2x,i10,20x)
      end do
! prd locs and orientation
      read(mu,146) ovf(n)%num_prd_sets
      146 format(2x,i10,20x)
      do m=1,ovf(n)%num_prd_sets
         read(mu,147) ovf(n)%num_prd(m)
         147 format(2x,i10,20x)
         do mp=1,ovf(n)%num_prd(m)
            read(mu,148) ovf(n)%loc_prd(m,mp)%i
            148 format(2x,i10,20x)
            read(mu,149) ovf(n)%loc_prd(m,mp)%j
            149 format(2x,i10,20x)
            read(mu,150) ovf(n)%loc_prd(m,mp)%i_adv
            150 format(2x,i10,20x)
            read(mu,151) ovf(n)%loc_prd(m,mp)%j_adv
            151 format(2x,i10,20x)
            read(mu,152) ovf(n)%loc_prd(m,mp)%i_u
            152 format(2x,i10,20x)
            read(mu,153) ovf(n)%loc_prd(m,mp)%j_u
            153 format(2x,i10,20x)
            read(mu,154) ovf(n)%loc_prd(m,mp)%k
            154 format(2x,i10,20x)
            read(mu,155) ovf(n)%loc_prd(m,mp)%orient
            155 format(2x,i10,20x)
         end do
      end do
! adjacent boundaries
! src
      read(mu,156) ovf(n)%adj_src%imin 
      156 format(2x,i10,20x)
      read(mu,157) ovf(n)%adj_src%imax 
      157 format(2x,i10,20x)
      read(mu,158) ovf(n)%adj_src%jmin 
      158 format(2x,i10,20x)
      read(mu,159) ovf(n)%adj_src%jmax 
      159 format(2x,i10,20x)
      read(mu,160) ovf(n)%adj_src%kmin 
      160 format(2x,i10,20x)
      read(mu,161) ovf(n)%adj_src%kmax 
      161 format(2x,i10,20x)
!ent
      read(mu,162) ovf(n)%adj_ent%imin 
      162 format(2x,i10,20x)
      read(mu,163) ovf(n)%adj_ent%imax 
      163 format(2x,i10,20x)
      read(mu,164) ovf(n)%adj_ent%jmin 
      164 format(2x,i10,20x)
      read(mu,165) ovf(n)%adj_ent%jmax 
      165 format(2x,i10,20x)
      read(mu,166) ovf(n)%adj_ent%kmin 
      166 format(2x,i10,20x)
      read(mu,167) ovf(n)%adj_ent%kmax 
      167 format(2x,i10,20x)
!prd
      do m=1,ovf(n)%num_prd_sets
         read(mu,168) ovf(n)%adj_prd(m)%imin 
         168 format(2x,i10,20x)
         read(mu,169) ovf(n)%adj_prd(m)%imax 
         169 format(2x,i10,20x)
         read(mu,170) ovf(n)%adj_prd(m)%jmin 
         170 format(2x,i10,20x)
         read(mu,171) ovf(n)%adj_prd(m)%jmax 
         171 format(2x,i10,20x)
         read(mu,172) ovf(n)%adj_prd(m)%kmin 
         172 format(2x,i10,20x)
         read(mu,173) ovf(n)%adj_prd(m)%kmax 
         173 format(2x,i10,20x)
      end do
! transports
      read(mu,174) ovf(n)%Ms
      174 format(2x,1PE27.18)
      read(mu,175) ovf(n)%Ms_n
      175 format(2x,1PE27.18)
      read(mu,176) ovf(n)%Ms_nm1
      176 format(2x,1PE27.18)
      read(mu,177) ovf(n)%Me
      177 format(2x,1PE27.18)
      read(mu,178) ovf(n)%Me_n
      178 format(2x,1PE27.18)
      read(mu,179) ovf(n)%Me_nm1
      179 format(2x,1PE27.18)
      read(mu,180) ovf(n)%phi
      180 format(2x,1PE27.18)
      read(mu,181) ovf(n)%Mp
      181 format(2x,1PE27.18)
      read(mu,182) ovf(n)%Mp_n
      182 format(2x,1PE27.18)
      read(mu,183) ovf(n)%Mp_nm1
      183 format(2x,1PE27.18)
      read(mu,184) ovf(n)%Tp
      184 format(2x,1PE27.18)
      read(mu,185) ovf(n)%Sp
      185 format(2x,1PE27.18)
      read(mu,186) ovf(n)%prd_set_n
      read(mu,186) ovf(n)%prd_set
      186 format(2x,i10,20x)

   end do  ! ovf loop 

   close(mu)
   endif     ! my_task == master_task

   call release_unit(mu)

!-----------------------------------------------------------------------
!EOC

 end subroutine ovf_read_restart

!***********************************************************************
!EOP
! !IROUTINE: ovf_read_broadcast
! !INTERFACE:

 subroutine ovf_read_broadcast

! !DESCRIPTION:
!  This routine broadcasts selected data in ovf array from the
!  master_task to all processors.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------
   save

   integer   (int_kind)  :: &
      n,                    &! ovf loop index
      m,                    &! sub-ovf loop index
      nn,                   &! tracer loop index 
      mp,                   &! sub-ovf sub-loop index
      k                      ! vertical loop index

!-----------------------------------------------------------------------
!  if overflows off, exit
!-----------------------------------------------------------------------

   if( .not. overflows_on ) return

!-----------------------------------------------------------------------
!  broadcast overflows info to all processors
!-----------------------------------------------------------------------

   call broadcast_scalar(num_ovf, master_task)
   do n=1,num_ovf
      call broadcast_scalar(ovf(n)%name, master_task)
! ovf data
      call broadcast_scalar(ovf(n)%ovf_params%lat, master_task)
      call broadcast_scalar(ovf(n)%ovf_params%width, master_task)
      call broadcast_scalar(ovf(n)%ovf_params%source_thick, master_task)
      call broadcast_scalar(ovf(n)%ovf_params%distnc_str_ssb, master_task)
      call broadcast_scalar(ovf(n)%ovf_params%bottom_slope, master_task)
      call broadcast_scalar(ovf(n)%ovf_params%bottom_drag, master_task)
! kmt changes, if any
      call broadcast_scalar(ovf(n)%num_kmt, master_task)
      do m=1,ovf(n)%num_kmt
         call broadcast_scalar(ovf(n)%loc_kmt(m)%i, master_task)
         call broadcast_scalar(ovf(n)%loc_kmt(m)%j, master_task)
         call broadcast_scalar(ovf(n)%loc_kmt(m)%korg, master_task)
         call broadcast_scalar(ovf(n)%loc_kmt(m)%knew, master_task)
      end do
! regional boundaries
!   inflow
      call broadcast_scalar(ovf(n)%reg_inf%imin, master_task)
      call broadcast_scalar(ovf(n)%reg_inf%imax, master_task)
      call broadcast_scalar(ovf(n)%reg_inf%jmin, master_task)
      call broadcast_scalar(ovf(n)%reg_inf%jmax, master_task)
      call broadcast_scalar(ovf(n)%reg_inf%kmin, master_task)
      call broadcast_scalar(ovf(n)%reg_inf%kmax, master_task)
!   source
      call broadcast_scalar(ovf(n)%reg_src%imin, master_task)
      call broadcast_scalar(ovf(n)%reg_src%imax, master_task)
      call broadcast_scalar(ovf(n)%reg_src%jmin, master_task)
      call broadcast_scalar(ovf(n)%reg_src%jmax, master_task)
      call broadcast_scalar(ovf(n)%reg_src%kmin, master_task)
      call broadcast_scalar(ovf(n)%reg_src%kmax, master_task)
!   entrainment
      call broadcast_scalar(ovf(n)%reg_ent%imin, master_task)
      call broadcast_scalar(ovf(n)%reg_ent%imax, master_task)
      call broadcast_scalar(ovf(n)%reg_ent%jmin, master_task)
      call broadcast_scalar(ovf(n)%reg_ent%jmax, master_task)
      call broadcast_scalar(ovf(n)%reg_ent%kmin, master_task)
      call broadcast_scalar(ovf(n)%reg_ent%kmax, master_task)
! src locs and orientation
      call broadcast_scalar(ovf(n)%num_src, master_task)
      do m=1,ovf(n)%num_src
         call broadcast_scalar(ovf(n)%loc_src(m)%i, master_task)
         call broadcast_scalar(ovf(n)%loc_src(m)%j, master_task)
         call broadcast_scalar(ovf(n)%loc_src(m)%i_adv, master_task)
         call broadcast_scalar(ovf(n)%loc_src(m)%j_adv, master_task)
         call broadcast_scalar(ovf(n)%loc_src(m)%i_u, master_task)
         call broadcast_scalar(ovf(n)%loc_src(m)%j_u, master_task)
         call broadcast_scalar(ovf(n)%loc_src(m)%k, master_task)
         call broadcast_scalar(ovf(n)%loc_src(m)%orient, master_task)
      end do
! ent locs and orientation
      call broadcast_scalar(ovf(n)%num_ent, master_task)
      do m=1,ovf(n)%num_ent
         call broadcast_scalar(ovf(n)%loc_ent(m)%i, master_task)
         call broadcast_scalar(ovf(n)%loc_ent(m)%j, master_task)
         call broadcast_scalar(ovf(n)%loc_ent(m)%i_adv, master_task)
         call broadcast_scalar(ovf(n)%loc_ent(m)%j_adv, master_task)
         call broadcast_scalar(ovf(n)%loc_ent(m)%i_u, master_task)
         call broadcast_scalar(ovf(n)%loc_ent(m)%j_u, master_task)
         call broadcast_scalar(ovf(n)%loc_ent(m)%k, master_task)
         call broadcast_scalar(ovf(n)%loc_ent(m)%orient, master_task)
      end do
! prd locs and orientation
      call broadcast_scalar(ovf(n)%num_prd_sets, master_task)
      do m=1,ovf(n)%num_prd_sets
         call broadcast_scalar(ovf(n)%num_prd(m), master_task)
         do mp=1,ovf(n)%num_prd(m)
            call broadcast_scalar(ovf(n)%loc_prd(m,mp)%i, master_task)
            call broadcast_scalar(ovf(n)%loc_prd(m,mp)%j, master_task)
            call broadcast_scalar(ovf(n)%loc_prd(m,mp)%i_adv, master_task)
            call broadcast_scalar(ovf(n)%loc_prd(m,mp)%j_adv, master_task)
            call broadcast_scalar(ovf(n)%loc_prd(m,mp)%i_u, master_task)
            call broadcast_scalar(ovf(n)%loc_prd(m,mp)%j_u, master_task)
            call broadcast_scalar(ovf(n)%loc_prd(m,mp)%k, master_task)
            call broadcast_scalar(ovf(n)%loc_prd(m,mp)%orient, master_task)
         end do
      end do
! adjacent boundaries
      call broadcast_scalar(ovf(n)%adj_src%imin, master_task) 
      call broadcast_scalar(ovf(n)%adj_src%imax, master_task) 
      call broadcast_scalar(ovf(n)%adj_src%jmin, master_task) 
      call broadcast_scalar(ovf(n)%adj_src%jmax, master_task) 
      call broadcast_scalar(ovf(n)%adj_src%kmin, master_task) 
      call broadcast_scalar(ovf(n)%adj_src%kmax, master_task) 
      call broadcast_scalar(ovf(n)%adj_ent%imin, master_task) 
      call broadcast_scalar(ovf(n)%adj_ent%imax, master_task) 
      call broadcast_scalar(ovf(n)%adj_ent%jmin, master_task) 
      call broadcast_scalar(ovf(n)%adj_ent%jmax, master_task) 
      call broadcast_scalar(ovf(n)%adj_ent%kmin, master_task) 
      call broadcast_scalar(ovf(n)%adj_ent%kmax, master_task) 
      do m=1,ovf(n)%num_prd_sets
         call broadcast_scalar(ovf(n)%adj_prd(m)%imin, master_task) 
         call broadcast_scalar(ovf(n)%adj_prd(m)%imax, master_task) 
         call broadcast_scalar(ovf(n)%adj_prd(m)%jmin, master_task) 
         call broadcast_scalar(ovf(n)%adj_prd(m)%jmax, master_task) 
         call broadcast_scalar(ovf(n)%adj_prd(m)%kmin, master_task) 
         call broadcast_scalar(ovf(n)%adj_prd(m)%kmax, master_task) 
      end do
! transports
      call broadcast_scalar(ovf(n)%Ms, master_task)
      call broadcast_scalar(ovf(n)%Ms_n, master_task)
      call broadcast_scalar(ovf(n)%Ms_nm1, master_task)
      call broadcast_scalar(ovf(n)%Me, master_task)
      call broadcast_scalar(ovf(n)%Me_n, master_task)
      call broadcast_scalar(ovf(n)%Me_nm1, master_task)
      call broadcast_scalar(ovf(n)%phi, master_task)
      call broadcast_scalar(ovf(n)%Mp, master_task)
      call broadcast_scalar(ovf(n)%Mp_n, master_task)
      call broadcast_scalar(ovf(n)%Mp_nm1, master_task)
      call broadcast_scalar(ovf(n)%Tp, master_task)
      call broadcast_scalar(ovf(n)%Sp, master_task)
      call broadcast_scalar(ovf(n)%prd_set_n, master_task)
      call broadcast_scalar(ovf(n)%prd_set, master_task)
      do nn=1,nt
         call broadcast_scalar(ovf(n)%trcr_reg%inf(nn), master_task)
         call broadcast_scalar(ovf(n)%trcr_reg%src(nn), master_task)
         call broadcast_scalar(ovf(n)%trcr_reg%ent(nn), master_task)
         call broadcast_scalar(ovf(n)%trcr_reg%prd(nn), master_task)
         call broadcast_scalar(ovf(n)%trcr_adj%src(nn), master_task)
         call broadcast_scalar(ovf(n)%trcr_adj%ent(nn), master_task)
         call broadcast_scalar(ovf(n)%trcr_adj%prd(nn), master_task)
      end do

   end do  ! ovf broadcast loop for all processors

!-----------------------------------------------------------------------
!EOC

 end subroutine ovf_read_broadcast

!***********************************************************************
!EOP
! !IROUTINE: ovf_advt
! !INTERFACE:

 subroutine ovf_advt(k,TRACER_E,TRACER_N,ntr,this_block, &
                     CE,CW,CN,CS)

! !DESCRIPTION:
!  Modify tracer grid interface value for advection for
!  overflow points; orientation determines if E or N is modified
! !REVISION HISTORY:
!  same as module
!
! ovf t-grid box ij showing advection grid boxes 
! (i_adv,j_adv) set by orientation
!                                        ij+1 
!                                         
!                                     ____2_____ 
!              y ^                   |          | 
!                |                   |          |
!                |            i-1j  3|    ij    |1  i+1j
!                +----->             |          |
!                      x             |__________|
!                                         4
!       
!                                         ij-1   
!
! Note!  Orientations are relative to overflow ij, while
! the advection boxes are offset as in the diagram above.
! Thus, the indices for TRACER_E and TRACER_N are reversed.
! For instance, orient=1 src, the advection box is i+1,j
! above, but when ij is that box (see below), then it is the
! western TRACER_E, or i-1j, that is overwritten. This is
! reversed from the center ij, because of the offset in the
! advection boxes relative to box ij.
!
! Note!  ij loops include ghost points incase advection
! scheme requires them.

!EOP
!BOC
!-----------------------------------------------------------------------
!  input variables
!-----------------------------------------------------------------------

   integer (int_kind), intent(in)                            :: &
      k                     ! vertical index
   real (r8), dimension(nx_block,ny_block), intent(inout)    :: &
      TRACER_E,           & ! east gridbox interface tracer at level k
      TRACER_N              ! north gridbox interface tracer at level k
   integer (int_kind), intent(in)                            :: &
      ntr                   ! tracer index
   type (block), intent(in)                                  :: &
      this_block            ! block information for this block

   real (r8), dimension(nx_block,ny_block), intent(in) :: &
      CN,CS,CE,CW           ! stencil weights based on flux velocities

!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind)  :: &
      n,m,mp,i,j,         & ! dummy loop indices
      ksrc,kent,kprd        ! overflow level indices

   integer (int_kind)  :: &
      iblock                ! local block address for this block

   logical (log_kind), parameter :: prnt = .false.

   integer (int_kind)  :: ovf_id


! turn off print   3 Nov 2008
!   if( prnt .and. my_task == master_task ) then
!      write(stdout,*) 'ovf_advt called '
!      call shr_sys_flush(stdout)
!   endif

   iblock = this_block%local_id

!-----------------------------------------------------------------------
!  overflow loop
!-----------------------------------------------------------------------

   !must have the groups initilaized
   if (.not. ovf_groups%init) then
      call ovf_init_groups()
   end if


! for each of MY overflows
   do n=1,ovf_groups%numMyGroups
      ovf_id = ovf_groups%groupIds(n)
! src
      do m=1,ovf(ovf_id)%num_src  ! source
         ksrc = ovf(ovf_id)%loc_src(m)%k
         if( k == ksrc ) then
            do j=1,ny_block
               if( ovf(ovf_id)%loc_src(m)%j_adv .eq. this_block%j_glob(j) ) then
                  do i=1,nx_block
                     if( ovf(ovf_id)%loc_src(m)%i_adv .eq. this_block%i_glob(i) ) then
                        if( prnt ) then
                           write(stdout,5) nsteps_total,n,ovf(ovf_id)%loc_src(m)%i_adv, &
                           ovf(ovf_id)%loc_src(m)%j_adv,ovf(ovf_id)%loc_src(m)%k, &
                           ovf(ovf_id)%loc_src(m)%orient,ntr, &
                           TRACER_E(i,j),TRACER_E(i-1,j),TRACER_N(i,j),TRACER_N(i,j-1)
                           5  format(' In    ovf_advt src   ',i5,1x,6(i3,1x),4(1pe15.8,1x))
                        endif  ! print
                        if( i > 1 ) then
                           if( ovf(ovf_id)%loc_src(m)%orient .eq. 1 ) then
                              TRACER_E(i-1,j) = ovf(ovf_id)%trcr_reg%src(ntr)
                           endif
                        endif
                        if( j > 1 ) then
                           if( ovf(ovf_id)%loc_src(m)%orient .eq. 2 ) then
                              TRACER_N(i,j-1) = ovf(ovf_id)%trcr_reg%src(ntr)
                           endif
                        endif
                        if( ovf(ovf_id)%loc_src(m)%orient .eq. 3 ) then
                           TRACER_E(i,j)   = ovf(ovf_id)%trcr_reg%src(ntr)
                        endif
                        if( ovf(ovf_id)%loc_src(m)%orient .eq. 4 ) then
                           TRACER_N(i,j)   = ovf(ovf_id)%trcr_reg%src(ntr)
                        endif
                        if( prnt ) then
                           write(stdout,10) nsteps_total,n,ovf(ovf_id)%loc_src(m)%i_adv, &
                           ovf(ovf_id)%loc_src(m)%j_adv,ovf(ovf_id)%loc_src(m)%k, &
                           ovf(ovf_id)%loc_src(m)%orient,ntr, &
                           TRACER_E(i,j),TRACER_E(i-1,j),TRACER_N(i,j),TRACER_N(i,j-1), &
                                            nsteps_total,n,ovf(ovf_id)%loc_src(m)%i_adv, &
                           ovf(ovf_id)%loc_src(m)%j_adv,ovf(ovf_id)%loc_src(m)%k, &
                           ovf(ovf_id)%loc_src(m)%orient,ntr, &
           CE(i,j)*dz(ksrc)*TAREA(i,j,iblock),CW(i,j)*dz(ksrc)*TAREA(i,j,iblock), &
           CN(i,j)*dz(ksrc)*TAREA(i,j,iblock),CS(i,j)*dz(ksrc)*TAREA(i,j,iblock), &
                                            nsteps_total,n,ovf(ovf_id)%loc_src(m)%i_adv, &
                           ovf(ovf_id)%loc_src(m)%j_adv,ovf(ovf_id)%loc_src(m)%k, &
                           ovf(ovf_id)%loc_src(m)%orient,ntr, &
           CE(i,j)*dz(ksrc)*TAREA(i,j,iblock)*TRACER_E(i,j), &
           CW(i,j)*dz(ksrc)*TAREA(i,j,iblock)*TRACER_E(i-1,j)
                           10 format(' Out   ovf_advt src   ',i5,1x,6(i3,1x),4(1pe15.8,1x)/ &
                                     ' Out   ovf_advt src M ',i5,1x,6(i3,1x),4(1pe15.8,1x)/ &
                                     ' Out   ovf_advt src CT',i5,1x,6(i3,1x),2(1pe15.8,1x))
                        endif  ! print
                     endif
                  end do  ! i
               endif
            end do  ! j
         endif  ! k
      end do  ! source
! ent
      do m=1,ovf(ovf_id)%num_ent  ! entrainment
         kent = ovf(ovf_id)%loc_ent(m)%k
         if( k == kent ) then
            do j=1,ny_block
               if( ovf(ovf_id)%loc_ent(m)%j_adv .eq. this_block%j_glob(j) ) then
                  do i=1,nx_block
                     if( ovf(ovf_id)%loc_ent(m)%i_adv .eq. this_block%i_glob(i) ) then
                        if( prnt ) then
                           write(stdout,15) nsteps_total,n,ovf(ovf_id)%loc_ent(m)%i_adv, &
                           ovf(ovf_id)%loc_ent(m)%j_adv,ovf(ovf_id)%loc_ent(m)%k, &
                           ovf(ovf_id)%loc_ent(m)%orient,ntr, &
                           TRACER_E(i,j),TRACER_E(i-1,j),TRACER_N(i,j),TRACER_N(i,j-1)
                           15 format(' In    ovf_advt ent   ',i5,1x,6(i3,1x),4(1pe15.8,1x))
                        endif  ! print
                        if( i > 1 ) then
                           if( ovf(ovf_id)%loc_ent(m)%orient .eq. 1 ) then
                              TRACER_E(i-1,j) = ovf(ovf_id)%trcr_reg%ent(ntr)
                           endif
                        endif
                        if( j > 1 ) then
                           if( ovf(ovf_id)%loc_ent(m)%orient .eq. 2 ) then
                              TRACER_N(i,j-1) = ovf(ovf_id)%trcr_reg%ent(ntr)
                           endif
                        endif
                        if( ovf(ovf_id)%loc_ent(m)%orient .eq. 3 ) then
                           TRACER_E(i,j)   = ovf(ovf_id)%trcr_reg%ent(ntr)
                        endif
                        if( ovf(ovf_id)%loc_ent(m)%orient .eq. 4 ) then
                           TRACER_N(i,j)   = ovf(ovf_id)%trcr_reg%ent(ntr)
                        endif
                        if( prnt ) then
                           write(stdout,20) nsteps_total,n,ovf(ovf_id)%loc_ent(m)%i_adv, &
                           ovf(ovf_id)%loc_ent(m)%j_adv,ovf(ovf_id)%loc_ent(m)%k, &
                           ovf(ovf_id)%loc_ent(m)%orient,ntr, &
                           TRACER_E(i,j),TRACER_E(i-1,j),TRACER_N(i,j),TRACER_N(i,j-1), &
                                            nsteps_total,n,ovf(ovf_id)%loc_ent(m)%i_adv, &
                           ovf(ovf_id)%loc_ent(m)%j_adv,ovf(ovf_id)%loc_ent(m)%k, &
                           ovf(ovf_id)%loc_ent(m)%orient,ntr, &
           CE(i,j)*dz(kent)*TAREA(i,j,iblock),CW(i,j)*dz(kent)*TAREA(i,j,iblock), &
           CN(i,j)*dz(kent)*TAREA(i,j,iblock),CS(i,j)*dz(kent)*TAREA(i,j,iblock), &
                                            nsteps_total,n,ovf(ovf_id)%loc_ent(m)%i_adv, &
                           ovf(ovf_id)%loc_ent(m)%j_adv,ovf(ovf_id)%loc_ent(m)%k, &
                           ovf(ovf_id)%loc_ent(m)%orient,ntr, &
           CE(i,j)*dz(kent)*TAREA(i,j,iblock)*TRACER_E(i,j), &
           CW(i,j)*dz(kent)*TAREA(i,j,iblock)*TRACER_E(i-1,j)
                           20 format(' Out   ovf_advt ent   ',i5,1x,6(i3,1x),4(1pe15.8,1x)/ &
                                     ' Out   ovf_advt ent M ',i5,1x,6(i3,1x),4(1pe15.8,1x)/ &
                                     ' Out   ovf_advt ent CT',i5,1x,6(i3,1x),2(1pe15.8,1x))
                        endif  ! print
                     endif
                  end do  ! i
               endif
            end do  ! j
         endif  ! k
      end do  ! entrainment
! prd
      m = ovf(ovf_id)%prd_set  ! product set for insertion
      do mp=1,ovf(ovf_id)%num_prd(m)  ! product points for insertion
         kprd = ovf(ovf_id)%loc_prd(m,mp)%k
         if( k == kprd ) then
            do j=1,ny_block
               if( ovf(ovf_id)%loc_prd(m,mp)%j_adv .eq. this_block%j_glob(j) ) then
                  do i=1,nx_block
                     if( ovf(ovf_id)%loc_prd(m,mp)%i_adv .eq. this_block%i_glob(i) ) then
                        if( prnt ) then
                           write(stdout,25) nsteps_total,n,ovf(ovf_id)%loc_prd(m,mp)%i_adv, &
                           ovf(ovf_id)%loc_prd(m,mp)%j_adv,ovf(ovf_id)%loc_prd(m,mp)%k, &
                           ovf(ovf_id)%loc_prd(m,mp)%orient,ntr, &
                           TRACER_E(i,j),TRACER_E(i-1,j),TRACER_N(i,j),TRACER_N(i,j-1)
                           25 format(' In    ovf_advt prd   ',i5,1x,6(i3,1x),4(1pe15.8,1x))
                        endif  ! print
                        if( i > 1 ) then
                           if( ovf(ovf_id)%loc_prd(m,mp)%orient .eq. 1 ) then
                              TRACER_E(i-1,j)   = ovf(ovf_id)%trcr_reg%prd(ntr)
                           endif
                        endif
                        if( j > 1 ) then
                           if( ovf(ovf_id)%loc_prd(m,mp)%orient .eq. 2 ) then
                              TRACER_N(i,j-1)   = ovf(ovf_id)%trcr_reg%prd(ntr)
                           endif
                        endif
                        if( ovf(ovf_id)%loc_prd(m,mp)%orient .eq. 3 ) then
                           TRACER_E(i,j)     = ovf(ovf_id)%trcr_reg%prd(ntr)
                        endif
                        if( ovf(ovf_id)%loc_prd(m,mp)%orient .eq. 4 ) then
                           TRACER_N(i,j)     = ovf(ovf_id)%trcr_reg%prd(ntr)
                        endif
                        if( prnt ) then
                           write(stdout,35) nsteps_total,n,ovf(ovf_id)%loc_prd(m,mp)%i_adv, &
                           ovf(ovf_id)%loc_prd(m,mp)%j_adv,ovf(ovf_id)%loc_prd(m,mp)%k, &
                           ovf(ovf_id)%loc_prd(m,mp)%orient,ntr, &
                           TRACER_E(i,j),TRACER_E(i-1,j),TRACER_N(i,j),TRACER_N(i,j-1), &
                                            nsteps_total,n,ovf(ovf_id)%loc_prd(m,mp)%i_adv, &
                           ovf(ovf_id)%loc_prd(m,mp)%j_adv,ovf(ovf_id)%loc_prd(m,mp)%k, &
                           ovf(ovf_id)%loc_prd(m,mp)%orient,ntr, &
           CE(i,j)*dz(kprd)*TAREA(i,j,iblock),CW(i,j)*dz(kprd)*TAREA(i,j,iblock), &
           CN(i,j)*dz(kprd)*TAREA(i,j,iblock),CS(i,j)*dz(kprd)*TAREA(i,j,iblock), &
                                            nsteps_total,n,ovf(ovf_id)%loc_prd(m,mp)%i_adv, &
                           ovf(ovf_id)%loc_prd(m,mp)%j_adv,ovf(ovf_id)%loc_prd(m,mp)%k, &
                           ovf(ovf_id)%loc_prd(m,mp)%orient,ntr, &
           CE(i,j)*dz(kprd)*TAREA(i,j,iblock)*TRACER_E(i,j), &
           CW(i,j)*dz(kprd)*TAREA(i,j,iblock)*TRACER_E(i-1,j)
                           35 format(' Out   ovf_advt prd   ',i5,1x,6(i3,1x),4(1pe15.8,1x)/ &
                                     ' Out   ovf_advt prd M ',i5,1x,6(i3,1x),4(1pe15.8,1x)/ &
                                     ' Out   ovf_advt prd CT',i5,1x,6(i3,1x),2(1pe15.8,1x))
                        endif  ! print
                     endif
                  end do  ! i
               endif
            end do  ! j
         endif  ! k
      end do  ! product points for insertion set
! If prd set just moved and time averaging done previous time step
      if( ovf(ovf_id)%prd_set .ne. ovf(ovf_id)%prd_set_n ) then
       m = ovf(ovf_id)%prd_set_n  ! product set for insertion
       do mp=1,ovf(ovf_id)%num_prd(m)  ! product points for insertion
         kprd = ovf(ovf_id)%loc_prd(m,mp)%k
         if( k == kprd ) then
            do j=1,ny_block
               if( ovf(ovf_id)%loc_prd(m,mp)%j_adv .eq. this_block%j_glob(j) ) then
                  do i=1,nx_block
                     if( ovf(ovf_id)%loc_prd(m,mp)%i_adv .eq. this_block%i_glob(i) ) then
                        if( prnt ) then
                           write(stdout,26) nsteps_total,n,ovf(ovf_id)%loc_prd(m,mp)%i_adv, &
                           ovf(ovf_id)%loc_prd(m,mp)%j_adv,ovf(ovf_id)%loc_prd(m,mp)%k, &
                           ovf(ovf_id)%loc_prd(m,mp)%orient,ntr, &
                           TRACER_E(i,j),TRACER_E(i-1,j),TRACER_N(i,j),TRACER_N(i,j-1)
                           26 format(' In_n  ovf_advt prd   ',i5,1x,6(i3,1x),4(1pe15.8,1x))
                        endif  ! print
                        if( avg_ts_last ) then
                           if( i > 1 ) then
                              if( ovf(ovf_id)%loc_prd(m,mp)%orient .eq. 1 ) then
                                 TRACER_E(i-1,j)   = ovf(ovf_id)%trcr_reg%prd(ntr)
                              endif
                           endif
                           if( j > 1 ) then
                              if( ovf(ovf_id)%loc_prd(m,mp)%orient .eq. 2 ) then
                                 TRACER_N(i,j-1)   = ovf(ovf_id)%trcr_reg%prd(ntr)
                              endif
                           endif
                           if( ovf(ovf_id)%loc_prd(m,mp)%orient .eq. 3 ) then
                              TRACER_E(i,j)     = ovf(ovf_id)%trcr_reg%prd(ntr)
                           endif
                           if( ovf(ovf_id)%loc_prd(m,mp)%orient .eq. 4 ) then
                              TRACER_N(i,j)     = ovf(ovf_id)%trcr_reg%prd(ntr)
                           endif
                        endif
                        if( prnt ) then
                           write(stdout,36) nsteps_total,n,ovf(ovf_id)%loc_prd(m,mp)%i_adv, &
                           ovf(ovf_id)%loc_prd(m,mp)%j_adv,ovf(ovf_id)%loc_prd(m,mp)%k, &
                           ovf(ovf_id)%loc_prd(m,mp)%orient,ntr, &
                           TRACER_E(i,j),TRACER_E(i-1,j),TRACER_N(i,j),TRACER_N(i,j-1), &
                                            nsteps_total,n,ovf(ovf_id)%loc_prd(m,mp)%i_adv, &
                           ovf(ovf_id)%loc_prd(m,mp)%j_adv,ovf(ovf_id)%loc_prd(m,mp)%k, &
                           ovf(ovf_id)%loc_prd(m,mp)%orient,ntr, &
           CE(i,j)*dz(kprd)*TAREA(i,j,iblock),CW(i,j)*dz(kprd)*TAREA(i,j,iblock), &
           CN(i,j)*dz(kprd)*TAREA(i,j,iblock),CS(i,j)*dz(kprd)*TAREA(i,j,iblock), &
                                            nsteps_total,n,ovf(ovf_id)%loc_prd(m,mp)%i_adv, &
                           ovf(ovf_id)%loc_prd(m,mp)%j_adv,ovf(ovf_id)%loc_prd(m,mp)%k, &
                           ovf(ovf_id)%loc_prd(m,mp)%orient,ntr, &
           CE(i,j)*dz(kprd)*TAREA(i,j,iblock)*TRACER_E(i,j), &
           CW(i,j)*dz(kprd)*TAREA(i,j,iblock)*TRACER_E(i-1,j)
                           36 format(' Out_n ovf_advt prd   ',i5,1x,6(i3,1x),4(1pe15.8,1x)/ &
                                     ' Out_n ovf_advt prd M ',i5,1x,6(i3,1x),4(1pe15.8,1x)/ &
                                     ' Out_n ovf_advt prd CT',i5,1x,6(i3,1x),2(1pe15.8,1x))
                        endif  ! print
                     endif
                  end do  ! i
               endif
            end do  ! j
         endif  ! k
       end do  ! product points for insertion set
      endif
   end do  ! each overflow
! special diagnostic  11 nov 2008
   call ovf_UV_check

!----------------------------------------------------------------------
!EOC

 end subroutine ovf_advt

!***********************************************************************
!EOP
! !IROUTINE: ovf_wtkb_check
! !INTERFACE:

 subroutine ovf_wtkb_check(k,WTKB,this_block)

! !DESCRIPTION:
!  Print out wtkb for overflow gridboxes
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!  input variables
!-----------------------------------------------------------------------

   integer (int_kind), intent(in)                            :: &
      k                     ! vertical index
   real (r8), dimension(nx_block,ny_block,nblocks_clinic), intent(in) :: &
      WTKB                  ! WTKB = W at bottom of t-grid box
   type (block), intent(in)                                  :: &
      this_block            ! block information for this block

!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind)  :: &
      n,m,mp,i,j,         & ! dummy loop indices
      ib,ie,jb,je,        & ! local domain index boundaries
      iblock                ! local block address for this block
   logical (log_kind), parameter :: prnt = .false.

   integer (int_kind)  ::ovf_id


   if( prnt .and. my_task == master_task ) then
      write(stdout,*) 'ovf_wtkb_check called '
      call shr_sys_flush(stdout)
   endif

   iblock = this_block%local_id

!-----------------------------------------------------------------------
!  overflow loop
!-----------------------------------------------------------------------

 !must have the groups initilaized
   if (.not. ovf_groups%init) then
      call ovf_init_groups()
   end if  


! for each of MY overflows
   do n=1,ovf_groups%numMyGroups
      ovf_id = ovf_groups%groupIds(n)

! ovf ij
! src
      do m=1,ovf(ovf_id)%num_src  ! source
         ib = this_block%ib
         ie = this_block%ie
         jb = this_block%jb
         je = this_block%je
         do j=jb,je
            if( ovf(ovf_id)%loc_src(m)%j .eq. this_block%j_glob(j) ) then
               do i=ib,ie
                  if( ovf(ovf_id)%loc_src(m)%i .eq. this_block%i_glob(i) ) then
                     if( k == KMT(i,j,iblock) ) then
                        if( prnt ) then
                           write(stdout,10) n,nsteps_total,ovf(ovf_id)%loc_src(m)%i, &
  ovf(ovf_id)%loc_src(m)%j,k,WTKB(i,j,iblock),TAREA(i,j,iblock)*WTKB(i,j,iblock)
  10 format(' ovf_wtkb_ch n=',i3, &
  ' src t,i,j,k         wtkb wtkb*tarea=',4(i4,1x),2(1pe12.5,2x))
                        endif  ! print
                     endif
                  endif
               end do  ! i
            endif
         end do  ! j
      end do  ! source
! ent
      do m=1,ovf(ovf_id)%num_ent  ! entrainment
         ib = this_block%ib
         ie = this_block%ie
         jb = this_block%jb
         je = this_block%je
         do j=jb,je
            if( ovf(ovf_id)%loc_ent(m)%j .eq. this_block%j_glob(j) ) then
               do i=ib,ie
                  if( ovf(ovf_id)%loc_ent(m)%i .eq. this_block%i_glob(i) ) then
                     if( k == KMT(i,j,iblock) ) then
                        if( prnt ) then
                           write(stdout,20) n,nsteps_total,ovf(ovf_id)%loc_ent(m)%i, &
  ovf(ovf_id)%loc_ent(m)%j,k,WTKB(i,j,iblock),TAREA(i,j,iblock)*WTKB(i,j,iblock)
  20 format(' ovf_wtkb_ch n=',i3, &
  ' ent t,i,j,k         wtkb wtkb*tarea=',4(i4,1x),2(1pe12.5,2x))
                        endif  ! print
                     endif
                  endif
               end do  ! i
            endif
         end do  ! j
      end do  ! entrainment
! prd
      m = ovf(ovf_id)%prd_set  ! product set for insertion
      do mp=1,ovf(ovf_id)%num_prd(m)  ! product points for insertion
         ib = this_block%ib
         ie = this_block%ie
         jb = this_block%jb
         je = this_block%je
         do j=jb,je
            if( ovf(ovf_id)%loc_prd(m,mp)%j .eq. this_block%j_glob(j) ) then
               do i=ib,ie
                  if( ovf(ovf_id)%loc_prd(m,mp)%i .eq. this_block%i_glob(i) ) then
                     if( k == KMT(i,j,iblock) ) then
                        if( prnt ) then
                           write(stdout,30) n,nsteps_total,ovf(ovf_id)%loc_prd(m,mp)%i, & 
  ovf(ovf_id)%loc_prd(m,mp)%j,k,WTKB(i,j,iblock),TAREA(i,j,iblock)*WTKB(i,j,iblock)
  30 format(' ovf_wtkb_ch n=',i3, & 
  ' prd t,i,j,k         wtkb wtkb*tarea=',4(i4,1x),2(1pe12.5,2x))
                        endif  ! print
                     endif
                  endif
               end do  ! i
            endif
         end do  ! j
      end do  ! product points for insertion set
! prd
      do m=1,ovf(ovf_id)%num_prd_sets
       do mp=1,ovf(ovf_id)%num_prd(m)  ! product points for insertion
         ib = this_block%ib
         ie = this_block%ie
         jb = this_block%jb
         je = this_block%je
         do j=jb,je
            if( ovf(ovf_id)%loc_prd(m,mp)%j .eq. this_block%j_glob(j) ) then
               do i=ib,ie
                  if( ovf(ovf_id)%loc_prd(m,mp)%i .eq. this_block%i_glob(i) ) then
                     if( k == KMT(i,j,iblock) ) then
                        if( prnt ) then
                           write(stdout,31) n,nsteps_total,ovf(ovf_id)%loc_prd(m,mp)%i, & 
  ovf(ovf_id)%loc_prd(m,mp)%j,k,WTKB(i,j,iblock),TAREA(i,j,iblock)*WTKB(i,j,iblock)
  31 format(' ovf_wtkb_ch n=',i3, & 
  ' all prd t,i,j,k         wtkb wtkb*tarea=',4(i4,1x),2(1pe12.5,2x))
                        endif  ! print
                     endif
                  endif
               end do  ! i
            endif
         end do  ! j
       end do  ! product points for insertion set
      end do  ! product sets
! ovf i_adv j_adv
! src
      do m=1,ovf(ovf_id)%num_src  ! source
         ib = this_block%ib
         ie = this_block%ie
         jb = this_block%jb
         je = this_block%je
         do j=jb,je
            if( ovf(ovf_id)%loc_src(m)%j_adv .eq. this_block%j_glob(j) ) then
               do i=ib,ie
                  if( ovf(ovf_id)%loc_src(m)%i_adv .eq. this_block%i_glob(i) ) then
                     if( k == KMT(i,j,iblock) ) then
                        if( prnt ) then
                           write(stdout,40) n,nsteps_total,ovf(ovf_id)%loc_src(m)%i_adv, &
  ovf(ovf_id)%loc_src(m)%j_adv,k,WTKB(i,j,iblock),TAREA(i,j,iblock)*WTKB(i,j,iblock)
  40 format(' ovf_wtkb_ch n=',i3, &
  ' src t,i_adv,j_adv,k wtkb wtkb*tarea=',4(i4,1x),2(1pe12.5,2x))
                        endif  ! print
                     endif  ! k
                  endif
               end do  ! i
            endif
         end do  ! j
      end do  ! source
! ent
      do m=1,ovf(ovf_id)%num_ent  ! entrainment
         ib = this_block%ib
         ie = this_block%ie
         jb = this_block%jb
         je = this_block%je
         do j=jb,je
            if( ovf(ovf_id)%loc_ent(m)%j_adv .eq. this_block%j_glob(j) ) then
               do i=ib,ie
                  if( ovf(ovf_id)%loc_ent(m)%i_adv .eq. this_block%i_glob(i) ) then
                     if( k == KMT(i,j,iblock) ) then
                        if( prnt ) then
                           write(stdout,50) n,nsteps_total,ovf(ovf_id)%loc_ent(m)%i_adv, &
  ovf(ovf_id)%loc_ent(m)%j_adv,k,WTKB(i,j,iblock),TAREA(i,j,iblock)*WTKB(i,j,iblock)
  50 format(' ovf_wtkb_ch n=',i3, &
  ' ent t,i_adv,j_adv,k wtkb wtkb*tarea=',4(i4,1x),2(1pe12.5,2x))
                        endif  ! print
                     endif  ! k
                  endif
               end do  ! i
            endif
         end do  ! j
      end do  ! entrainment
! prd
      m = ovf(ovf_id)%prd_set  ! product set for insertion
      do mp=1,ovf(ovf_id)%num_prd(m)  ! product points for insertion
         ib = this_block%ib
         ie = this_block%ie
         jb = this_block%jb
         je = this_block%je
         do j=jb,je
            if( ovf(ovf_id)%loc_prd(m,mp)%j_adv .eq. this_block%j_glob(j) ) then
               do i=ib,ie
                  if( ovf(ovf_id)%loc_prd(m,mp)%i_adv .eq. this_block%i_glob(i) ) then
                     if( k == KMT(i,j,iblock) ) then
                        if( prnt ) then
                           write(stdout,60) n,nsteps_total,ovf(ovf_id)%loc_prd(m,mp)%i_adv, & 
  ovf(ovf_id)%loc_prd(m,mp)%j_adv,k,WTKB(i,j,iblock),TAREA(i,j,iblock)*WTKB(i,j,iblock)
  60 format(' ovf_wtkb_ch n=',i3, & 
  ' prd t,i_adv,j_adv,k wtkb wtkb*tarea=',4(i4,1x),2(1pe12.5,2x))
                        endif  ! print
                     endif  ! k
                  endif
               end do  ! i
            endif
         end do  ! j
      end do  ! product points for insertion set
! prd
      do m=1,ovf(ovf_id)%num_prd_sets
       do mp=1,ovf(ovf_id)%num_prd(m)  ! product points for insertion if moved
         ib = this_block%ib
         ie = this_block%ie
         jb = this_block%jb
         je = this_block%je
         do j=jb,je
            if( ovf(ovf_id)%loc_prd(m,mp)%j_adv .eq. this_block%j_glob(j) ) then
               do i=ib,ie
                  if( ovf(ovf_id)%loc_prd(m,mp)%i_adv .eq. this_block%i_glob(i) ) then
                     if( k == KMT(i,j,iblock) ) then
                        if( prnt ) then
                           write(stdout,61) n,nsteps_total,ovf(ovf_id)%loc_prd(m,mp)%i_adv, & 
  ovf(ovf_id)%loc_prd(m,mp)%j_adv,k,WTKB(i,j,iblock),TAREA(i,j,iblock)*WTKB(i,j,iblock)
  61 format(' ovf_wtkb_ch n=',i3, & 
  ' all prd t,i_adv,j_adv,k wtkb wtkb*tarea=',4(i4,1x),2(1pe12.5,2x))
                        endif  ! print
                     endif  ! k
                  endif
               end do  ! i
            endif
         end do  ! j
       end do  ! product points for insertion set
      end do  ! original product set if moved
   end do  ! each overflow

!----------------------------------------------------------------------
!EOC

 end subroutine ovf_wtkb_check

!***********************************************************************
!EOP
! !IROUTINE: ovf_UV_check
! !INTERFACE:

 subroutine ovf_UV_check

! !DESCRIPTION:
!  Print out column UVEL, VVEL for overflow gridboxes
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind)  :: &
      n,m,mp,i,j,k,       & ! dummy loop indices
      ib,ie,jb,je,        & ! local domain index boundaries
      iblock,             & ! local block address for this block
      ksrc,kent,kprd        ! overflow level indices
   type (block)        :: &
      this_block            ! block information for this block
   logical (log_kind), parameter :: prnt = .false.

   integer (int_kind) :: ovf_id

   if( prnt .and. my_task == master_task ) then
      write(stdout,*) 'ovf_UV_check called '
      call shr_sys_flush(stdout)
   endif

   if( prnt ) then
     write(stdout,5) nsteps_total
     5 format(' ovf_UV_check called at nsteps_total=',i6)
!-----------------------------------------------------------------------
!  overflow loop
!-----------------------------------------------------------------------

 !must have the groups initilaized
   if (.not. ovf_groups%init) then
      call ovf_init_groups()
   end if  


! for each of MY overflows
   do n=1,ovf_groups%numMyGroups
      ovf_id = ovf_groups%groupIds(n)
! src
       do m=1,ovf(ovf_id)%num_src  ! source
         ksrc = ovf(ovf_id)%loc_src(m)%k
         do iblock = 1,nblocks_clinic
           this_block = get_block(blocks_clinic(iblock),iblock)
           ib = this_block%ib
           ie = this_block%ie
           jb = this_block%jb
           je = this_block%je
           do j=jb,je
             if( ovf(ovf_id)%loc_src(m)%j_u .eq. this_block%j_glob(j) ) then
               do i=ib,ie
                 if( ovf(ovf_id)%loc_src(m)%i_u .eq. this_block%i_glob(i) ) then
                   write(stdout,15) n,ovf(ovf_id)%loc_src(m)%i_u, &
                   ovf(ovf_id)%loc_src(m)%j_u
                   15 format(' ovf_UV_check n=',i2,' src i_u j_u = ',2(i3,1x))
!                   do k=1,ksrc
                     k=ksrc
!                     write(stdout,10) k,UVEL(i,j,k,oldtime,iblock), &
!                        UVEL(i,j,k,curtime,iblock),UVEL(i,j,k,newtime,iblock), &
!                        VVEL(i,j,k,oldtime,iblock), &
!                        VVEL(i,j,k,curtime,iblock),VVEL(i,j,k,newtime,iblock)
                     10 format('   k old cur new UVEL= ',i2,1x,3(f9.5,1x), &
                     ' VVEL=',3(f9.5,1x))
!                   end do  ! k
                 endif
               end do  ! i
             endif
           end do  ! j
         enddo   ! block
       end do  ! source
! ent
       do m=1,ovf(ovf_id)%num_ent  ! entrainment
         kent = ovf(ovf_id)%loc_ent(m)%k
         do iblock = 1,nblocks_clinic
           this_block = get_block(blocks_clinic(iblock),iblock)
           ib = this_block%ib
           ie = this_block%ie
           jb = this_block%jb
           je = this_block%je
           do j=jb,je
             if( ovf(ovf_id)%loc_ent(m)%j_u .eq. this_block%j_glob(j) ) then
               do i=ib,ie
                 if( ovf(ovf_id)%loc_ent(m)%i_u .eq. this_block%i_glob(i) ) then
                   write(stdout,25) n,ovf(ovf_id)%loc_ent(m)%i_u, &
                   ovf(ovf_id)%loc_ent(m)%j_u
                   25 format(' ovf_UV_check n=',i2,' ent i_u j_u = ',2(i3,1x))
!                   do k=1,kent
                     k=kent
!                     write(stdout,20) k,UVEL(i,j,k,oldtime,iblock), &
!                        UVEL(i,j,k,curtime,iblock),UVEL(i,j,k,newtime,iblock), &
!                        VVEL(i,j,k,oldtime,iblock), &
!                        VVEL(i,j,k,curtime,iblock),VVEL(i,j,k,newtime,iblock)
                     20 format('   k old cur new UVEL= ',i2,1x,3(f9.5,1x), &
                     ' VVEL=',3(f9.5,1x))
!                   end do  ! k
                 endif
               end do  ! i
             endif
           end do  ! j
         enddo   ! block
       end do  ! entrainment
! prd
      do m=1,ovf(ovf_id)%num_prd_sets
       do mp=1,ovf(ovf_id)%num_prd(m)
         kprd = ovf(ovf_id)%loc_prd(m,mp)%k 
         do iblock = 1,nblocks_clinic
           this_block = get_block(blocks_clinic(iblock),iblock)
           ib = this_block%ib
           ie = this_block%ie
           jb = this_block%jb
           je = this_block%je
           do j=jb,je
             if( ovf(ovf_id)%loc_prd(m,mp)%j_u .eq. this_block%j_glob(j) ) then
               do i=ib,ie
                 if( ovf(ovf_id)%loc_prd(m,mp)%i_u .eq. this_block%i_glob(i) ) then
                   write(stdout,35) n,ovf(ovf_id)%loc_prd(m,mp)%i_u, &
                   ovf(ovf_id)%loc_prd(m,mp)%j_u
                   35 format(' ovf_UV_check n=',i2,' prd i_u j_u = ',2(i3,1x))
!                   do k=1,kprd
                     k=kprd
                     write(stdout,30) nsteps_total,n, &
                        ovf(ovf_id)%loc_prd(m,mp)%i_u,ovf(ovf_id)%loc_prd(m,mp)%j_u, &
                        k,UVEL(i,j,k,oldtime,iblock), &
                        UVEL(i,j,k,curtime,iblock),UVEL(i,j,k,newtime,iblock)
                     30 format(' prd t,n,i,j,k old cur new UVEL= ',5(i4,1x),1x,3(f9.5,1x))
!                   end do  ! k
                 endif
               end do  ! i
             endif
           end do  ! j
         enddo   ! block
       end do  ! product
      end do
     end do  ! each overflow
  endif  ! print
!----------------------------------------------------------------------
!EOC

 end subroutine ovf_UV_check

!***********************************************************************
!EOP
! !IROUTINE: ovf_Utlda
! !INTERFACE:

 subroutine ovf_Utlda(iblock)

! !DESCRIPTION:
!  Save ovf sidewall unnormalized baroclinic velocities Utlda. Must be
!  called AFTER the baroclinic solution Utlda is found but BEFORE the
!  baroclinic velocities are normalized (i.e. vertical integral of 
!  baroclinic velocity from surface to bottom topography is zero). 
!
! ij t-grid   i_u,j_u u-grid
!
! assignment of U on u-grid
!   orientation=1  i_u = i    j_u = j
!              =2  i_u = i-1  j_u = j
!              =3  i_u = i-1  j_u = j-1
!              =4  i_u = i    j_u = j-1
!
! ovf t-grid box ij with u-grid
! corners and orientations
!                                         2       (i_u,j_u)
!                                i-1j __________ij
!              y ^                   |          |
!                |                   |          |
!                |                 3 |    ij    | 1
!                +----->             |          |
!                      x             |__________|
!                              i-1j-1            ij-1
!                                         4
! for example, for ovf grid box ij,
! with product orientation 4, the Utlda
! in the above diagram would be ij-1
! lower right corner
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!  input variables
!-----------------------------------------------------------------------

   integer (int_kind),     &
      intent(in)        :: &
      iblock                 ! block index

!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind)    :: &
      i,j,k,n,m,mp,         &  ! dummy loop indices
      ib,ie,jb,je,          &  ! local domain index boundaries
      ksrc,kent,kprd           ! overflow level indices
   type (block)          :: &
      this_block               ! block information for current block
   logical (log_kind), parameter :: prnt = .false.

   integer (int_kind)  :: ovf_id

   if( prnt .and. my_task == master_task ) then
      write(stdout,*) 'ovf_Utlda called '
      call shr_sys_flush(stdout)
   endif

!-----------------------------------------------------------------------
!  overflow loop
!-----------------------------------------------------------------------
!must have the groups initilaized
   if (.not. ovf_groups%init) then
      call ovf_init_groups()
   end if  


! for each of MY overflows
   do n=1,ovf_groups%numMyGroups
      ovf_id = ovf_groups%groupIds(n)
! src
      do m=1,ovf(ovf_id)%num_src  ! source
         ksrc = ovf(ovf_id)%loc_src(m)%k
         this_block = get_block(blocks_clinic(iblock),iblock)
         ib = this_block%ib
         ie = this_block%ie
         jb = this_block%jb
         je = this_block%je
         do j=jb,je
            if( ovf(ovf_id)%loc_src(m)%j_u .eq. this_block%j_glob(j) ) then
               do i=ib,ie
                  if( ovf(ovf_id)%loc_src(m)%i_u .eq. this_block%i_glob(i) ) then
                     do k=1,ksrc-1
                        ovf(ovf_id)%loc_src(m)%Utlda(k) = UVEL(i,j,k,newtime,iblock)
                        ovf(ovf_id)%loc_src(m)%Vtlda(k) = VVEL(i,j,k,newtime,iblock)
                     enddo
                     if(prnt) then
                        write(stdout,10) n,ovf(ovf_id)%loc_src(m)%i_u, &
                                         ovf(ovf_id)%loc_src(m)%j_u, &
                                         ovf(ovf_id)%loc_src(m)%orient,ksrc
                        10 format(' ovf_Utlda n=',i3, &
                        ' src i_u j_u orient k=',4(i4,1x))
                        do k=1,ksrc-1
                           write(stdout,15) k,ovf(ovf_id)%loc_src(m)%Utlda(k), &
                                              ovf(ovf_id)%loc_src(m)%Vtlda(k)
                           15 format('   k=',i3,1x,'Utlda Vtlda= ',2(f9.5,2x))
                        enddo
                     endif
                  endif
               end do  ! i
            endif
         end do  ! j
      end do  ! source
! ent
      do m=1,ovf(ovf_id)%num_ent  ! entrainment
         kent = ovf(ovf_id)%loc_ent(m)%k
         this_block = get_block(blocks_clinic(iblock),iblock)
         ib = this_block%ib
         ie = this_block%ie
         jb = this_block%jb
         je = this_block%je
         do j=jb,je
            if( ovf(ovf_id)%loc_ent(m)%j_u .eq. this_block%j_glob(j) ) then
               do i=ib,ie
                  if( ovf(ovf_id)%loc_ent(m)%i_u .eq. this_block%i_glob(i) ) then
                     do k=1,kent-1
                        ovf(ovf_id)%loc_ent(m)%Utlda(k) = UVEL(i,j,k,newtime,iblock)
                        ovf(ovf_id)%loc_ent(m)%Vtlda(k) = VVEL(i,j,k,newtime,iblock)
                     enddo
                     if(prnt) then
                        write(stdout,20) n,ovf(ovf_id)%loc_ent(m)%i_u, &
                                         ovf(ovf_id)%loc_ent(m)%j_u, &
                                         ovf(ovf_id)%loc_ent(m)%orient,kent
                        20 format(' ovf_Utlda n=',i3, &
                        ' ent i_u j_u orient k=',4(i4,1x))
                        do k=1,kent-1
                           write(stdout,25) k,ovf(ovf_id)%loc_ent(m)%Utlda(k), &
                                              ovf(ovf_id)%loc_ent(m)%Vtlda(k)
                           25 format('   k=',i3,1x,'Utlda Vtlda= ',2(f9.5,2x))
                        enddo
                     endif
                  endif
               end do  ! i
            endif
         end do  ! j
      end do  ! entrainment
! prd
      do m=1,ovf(ovf_id)%num_prd_sets
         do mp=1,ovf(ovf_id)%num_prd(m)  ! product points for each set
            kprd = ovf(ovf_id)%loc_prd(m,mp)%k
            this_block = get_block(blocks_clinic(iblock),iblock)
            ib = this_block%ib
            ie = this_block%ie
            jb = this_block%jb
            je = this_block%je
            do j=jb,je
               if( ovf(ovf_id)%loc_prd(m,mp)%j_u .eq. this_block%j_glob(j) ) then
                  do i=ib,ie
                     if( ovf(ovf_id)%loc_prd(m,mp)%i_u .eq. this_block%i_glob(i) ) then
                        do k=1,kprd-1
                           ovf(ovf_id)%loc_prd(m,mp)%Utlda(k) = UVEL(i,j,k,newtime,iblock)
                           ovf(ovf_id)%loc_prd(m,mp)%Vtlda(k) = VVEL(i,j,k,newtime,iblock)
                        enddo
                        if(prnt) then
                           write(stdout,30) n,ovf(ovf_id)%loc_prd(m,mp)%i_u, &
                                            ovf(ovf_id)%loc_prd(m,mp)%j_u, &
                                            ovf(ovf_id)%loc_prd(m,mp)%orient,kprd
                           30 format(' ovf_Utlda n=',i3, &
                           ' prd i_u j_u orient k=',4(i4,1x))
                           do k=1,kprd-1
                              write(stdout,35) k,ovf(ovf_id)%loc_prd(m,mp)%Utlda(k), &
                                                 ovf(ovf_id)%loc_prd(m,mp)%Vtlda(k)
                              35 format('   k=',i3,1x,'Utlda Vtlda= ',2(f9.5,2x))
                           enddo
                        endif
                     endif
                  end do  ! i
               endif
            end do  ! j
         end do  ! product points for each set
      end do  ! product sets
   end do  ! each overflow

!----------------------------------------------------------------------
!EOC

 end subroutine ovf_Utlda

!***********************************************************************
!EOP
! !IROUTINE: ovf_driver
! !INTERFACE:

 subroutine ovf_driver

! !DESCRIPTION:
!  This routine is the main overflow (ovf) driver, called
!  in step_mod.F90 between baroclinic and barotropic drivers.
!  It calls routines to compute ovf regional means, transports,
!  product locations and sidewall velocity evaluation.
!
! !REVISION HISTORY:
!  same as module

   logical (log_kind), parameter :: prnt = .false.

!EOP
!BOC

   if(prnt) then
      write(stdout,*) ' ovf_driver entered '
      call shr_sys_flush(stdout)
   endif


   call timer_start(timer_ovf)

!----------------------------------------------------------------------
!
!  ovf regional averages
!
!----------------------------------------------------------------------

      call ovf_reg_avgs(curtime)

!----------------------------------------------------------------------
!
!  ovf transports
!
!----------------------------------------------------------------------

      call ovf_transports

!----------------------------------------------------------------------
!
!  ovf location of product
!
!----------------------------------------------------------------------

      call ovf_loc_prd

!----------------------------------------------------------------------
!
!  ovf top W evaluation
!
!----------------------------------------------------------------------

      call ovf_W

!----------------------------------------------------------------------
!
!  ovf sidewall UV evaluation
!
!----------------------------------------------------------------------

      call ovf_UV

!----------------------------------------------------------------------
!EOC


      call timer_stop(timer_ovf)


 end subroutine ovf_driver

!***********************************************************************
!EOP
! !IROUTINE: ovf_reg_avgs
! !INTERFACE:



subroutine ovf_reg_avgs(time_level)

! !DESCRIPTION:
!  Evaluate the ovf regional averages

! note: this is called initially outside of the ovf_driver

!
! !REVISION HISTORY:
!  same as module

!-----------------------------------------------------------------------
!  input variables
!-----------------------------------------------------------------------

   integer (int_kind) :: &! time indices for prognostic arrays
      time_level          ! current time level  (n)

!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (int_kind)      :: &
      iblock,k,n,nn,m             ! dummy loop indices
     

   type (block)            :: &
      this_block                ! block information for current block

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic) :: &
         WRK                    ! temp work array
   real (r8)  vsum_reg_wght,  & ! vertical sum regional weight
              vsum_adj_wght     ! vertical sum adjacent weight

   logical (log_kind), parameter :: prnt = .false.
   
   logical (log_kind) :: post
   logical (log_kind), dimension(:), allocatable :: post_array
   integer (int_kind) ::  num_posts, ierr, count, loc, ovf_id, i
   integer (int_kind), parameter :: len = 15 !length of buffer for a post
   real (r8), dimension(:), allocatable:: myRecvBuff(:), mySendBuff(:)


!EOP
!BOC

   if(prnt) then
      write(stdout,*) ' ovf_reg_avgs called '
      call shr_sys_flush(stdout)
   endif

   !must have the groups initialized
   if (.not. ovf_groups%init) then
      call ovf_init_groups()
   end if  


   !compute the averages
   call ovf_utils_avgs(time_level)

   if (prnt) then 

      if( my_task == master_task ) then
         !for any ovf groups that the master_task is not in, 
         !must post recv for the info
         allocate(post_array(num_ovf))
         num_posts = num_ovf - ovf_groups%numMyGroups
         if (num_posts > 0) then
            allocate(myRecvBuff(len*num_posts))
         endif

         call ovf_print_init(len, num_posts, myRecvBuff, post_array)
         
         count = 0
         do n=1,num_ovf
            if (post_array(n)) then !we posted a recv for this
               count = count+1
               call ovf_print_get(count)
               loc = (count-1)*len+1
               !unpack
               ovf(n)%trcr_reg%inf(1) = myRecvBuff(loc)
               ovf(n)%trcr_reg%inf(2) = myRecvBuff(loc+1)
               ovf(n)%rho_reg%inf = myRecvBuff(loc+2)
               ovf(n)%trcr_reg%src(1) = myRecvBuff(loc+3)
               ovf(n)%trcr_reg%src(2) = myRecvBuff(loc+4)
               ovf(n)%rho_reg%src = myRecvBuff(loc+5)
               ovf(n)%trcr_reg%ent(1) = myRecvBuff(loc+6)
               ovf(n)%trcr_reg%ent(2) = myRecvBuff(loc+7)
               ovf(n)%rho_reg%ent = myRecvBuff(loc+8)

               ovf(n)%trcr_adj%src(1) = myRecvBuff(loc+9)
               ovf(n)%trcr_adj%src(2) = myRecvBuff(loc+10)
               ovf(n)%trcr_adj%ent(1) = myRecvBuff(loc+11)
               ovf(n)%trcr_adj%ent(2) = myRecvBuff(loc+12)
               ovf(n)%rho_adj%prd(1) = myRecvBuff(loc+13)
               ovf(n)%rho_adj%prd(2) = myRecvBuff(loc+14)
            endif

            write(stdout,10) n,ovf(n)%trcr_reg%inf(1),                   &
                 (ovf(n)%trcr_reg%inf(2))*c1000,(ovf(n)%rho_reg%inf-c1)*c1000, & 
                 ovf(n)%trcr_reg%src(1),                                  &
                 (ovf(n)%trcr_reg%src(2))*c1000,(ovf(n)%rho_reg%src-c1)*c1000, & 
                 ovf(n)%trcr_reg%ent(1),                                  &
                 (ovf(n)%trcr_reg%ent(2))*c1000,(ovf(n)%rho_reg%ent-c1)*c1000
10          format(1x,'ovf reg',i3,1x,3(f6.3,1x),3(f6.3,1x),3(f6.3,1x))
            if( n.eq.1 ) then
               write(stdout,11) n,ovf(n)%trcr_adj%src(1),                   &
                    (ovf(n)%trcr_adj%src(2))*c1000,                           & 
                    ovf(n)%trcr_adj%ent(1),                                  &
                    (ovf(n)%trcr_adj%ent(2))*c1000,                           &
                    (ovf(n)%rho_adj%prd(1)-c1)*c1000
11             format(1x,'ovf adj',i3,1x,2(f6.3,1x),1x,2(f6.3,1x),f6.3)
            else
               write(stdout,12) n,ovf(n)%trcr_adj%src(1),                   &
                    (ovf(n)%trcr_adj%src(2))*c1000,                           & 
                    ovf(n)%trcr_adj%ent(1),                                  &
                    (ovf(n)%trcr_adj%ent(2))*c1000,                           &
                    (ovf(n)%rho_adj%prd(1)-c1)*c1000,                         &
                    (ovf(n)%rho_adj%prd(2)-c1)*c1000                         
12             format(1x,'ovf adj',i3,1x,2(f6.3,1x),1x,2(f6.3,1x),     &
                    f6.3,1x,f6.3)
            endif
         end do
         
         !clean up master
         deallocate(post_array)
         if (num_posts > 0) then
            deallocate(myRecvBuff)
         endif
         call ovf_print_finalize(num_posts)
 
      else !not master - send info to master if I am group master
         num_posts = 0
         do n=1, ovf_groups%numMyGroups
            if (ovf_groups%amMaster(n)) then
               num_posts = num_posts + 1
            end if
         end do
         
         call ovf_print_init(len, num_posts)

         if (num_posts > 0) then
            allocate(mySendBuff(num_posts * len))
         end if

         num_posts = 0
         do n=1, ovf_groups%numMyGroups
            if (ovf_groups%amMaster(n)) then
               ovf_id = ovf_groups%groupIds(n)
               !send transport info to the master task with isend
               num_posts = num_posts + 1
               loc = (num_posts-1)*len+1
               !pack send buffer
                mySendBuff(loc) = ovf(ovf_id)%trcr_reg%inf(1)
                mySendBuff(loc+1) = ovf(ovf_id)%trcr_reg%inf(2)
                mySendBuff(loc+2) = ovf(ovf_id)%rho_reg%inf
                mySendBuff(loc+3) = ovf(ovf_id)%trcr_reg%src(1)
                mySendBuff(loc+4) = ovf(ovf_id)%trcr_reg%src(2)
                mySendBuff(loc+5) = ovf(ovf_id)%rho_reg%src
                mySendBuff(loc+6) = ovf(ovf_id)%trcr_reg%ent(1)
                mySendBuff(loc+7) = ovf(ovf_id)%trcr_reg%ent(2)
                mySendBuff(loc+8) = ovf(ovf_id)%rho_reg%ent

                mySendBuff(loc+9) = ovf(ovf_id)%trcr_adj%src(1)
                mySendBuff(loc+10) = ovf(ovf_id)%trcr_adj%src(2)
                mySendBuff(loc+11) = ovf(ovf_id)%trcr_adj%ent(1)
                mySendBuff(loc+12) = ovf(ovf_id)%trcr_adj%ent(2)
                mySendBuff(loc+13) = ovf(ovf_id)%rho_adj%prd(1)
                mySendBuff(loc+14) = ovf(ovf_id)%rho_adj%prd(2)
         
                call ovf_print_send(len, mySendBuff(loc:loc+len-1), num_posts, ovf_id)

             endif
          end do  !n loop

          !now clean up
          call ovf_print_finalize(num_posts)
          if (num_posts > 0) then
             deallocate(mySendBuff)
          endif

       endif !end not master
   endif ! prnt loop

!----------------------------------------------------------------------
!EOC

 end subroutine ovf_reg_avgs

!***********************************************************************
!EOP
! !IROUTINE: ovf_transports
! !INTERFACE:

 subroutine ovf_transports

! !DESCRIPTION:
!  Evaluate the ovf transports. For each overflow, set overflow parameters
!  and evaluate transports.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------
   save

   integer   (int_kind)  :: &
     n          ,& ! ovf loop index
     nn         ,& ! ovf tracer index
     m          ,& ! product level
     k_p           ! product k level
!
   real (r8) ::  &
     lat        ,& ! inflow/source latitude for coriolis parameter (degrees)
     fs            ! coriolis parameter (/s)
!
   real (r8) ::  &
     hu         ,& ! upstream source thickness (cm)
     hs         ,& ! source water vertical thickness (cm)
     Ws         ,& ! source water width (cm)
     xse        ,& ! distance from source to entrainment (cm)
     di         ,& ! depth of inflow (cm)
     ds         ,& ! depth of source (cm)
     de         ,& ! depth of entrainment (cm)
     dp         ,& ! depth of product (cm)
     alpha      ,& ! continental slope between source to entrainment
     cd            ! bottom drag coefficient for spreading, entrained flow
!
   real (r8) ::  &
     T_i        ,& ! inflow mean temperature (C)
     S_i        ,& ! inflow mean salinity
     T_s        ,& ! source mean temperature (C)
     S_s        ,& ! source mean salinity
     T_e        ,& ! entrainment mean temperature (C)
     S_e        ,& ! entrainment mean salinity
     T_p        ,& ! product temperature (C)
     S_p           ! product salinity
!
   real (r8) ::  &
     rho_i      ,& ! inflow mass density (g/cm3)
     rho_s      ,& ! source mass density (g/cm3)
     rho_e      ,& ! entrainment mass density (g/cm3)
     rho_sed    ,& ! source at entrainment depth mass density (g/cm3)
     rho_p         ! product mass density (g/cm3)
!
   real (r8) ::  &
     gp_s       ,& ! source reduced gravity (cm/s2)
     Ms         ,& ! source mass flux (Sv)
     As         ,& ! source cross sectional area (cm2)
     Us            ! source speed (cm/s)
!
   real (r8) ::  &
     gp_e       ,& ! entrainment reduced gravity (cm/s2)
     Me         ,& ! entrainment mass flux (Sv)
     Ue         ,& ! entrainment speed (cm/s)
     Ugeo       ,& ! geostrophic entrainment speed (m/s)
     Uavg       ,& ! average source and geostrophic speed (cm/s)
     a,b,c      ,& ! parameters for quadratic solution
     Wgeo       ,& ! width of geostrophically spread source (cm)
     Kgeo       ,& ! geostrophic Ekman number
     hgeo       ,& ! depth of geostrophically spread source (cm)
     Fgeo       ,& ! Froude number of entrained flow
     phi        ,& ! entrainment parameter from actual ratio Me/Mp
     Mp            ! product mass flux (Sv)

   character (POP_charLength) ::  &
     string

   integer (POP_i4) ::  &
     ier

   integer (int_kind)  :: ovf_id


   integer (int_kind) ::  num_posts, ierr, count, loc, r_n, r_m, indx
   integer (int_kind), parameter :: len = 16 !length of buffer for a post
   real (r8), dimension(:), allocatable:: myRecvBuff(:), mySendBuff(:)


   logical (log_kind) :: post, found
   logical (log_kind), dimension(:), allocatable :: post_array

!EOP
!BOC
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  open overflows_diag_outfile file 
!  append overflows diagnostics to end of overflows diagnostics output file
!
!-----------------------------------------------------------------------
  
  if (my_task == master_task .and. eod) then
     open(ovf_diag_unit, file=overflows_diag_outfile, status='old', position='append')
  endif
 
  !must have the groups initialized
  if (.not. ovf_groups%init) then
     call ovf_init_groups()
  end if

  ! set up
  if( my_task == master_task ) then
     allocate(post_array(num_ovf))
  end if

  if (eod) then !end of day :)
     if (my_task == master_task) then
        num_posts = num_ovf - ovf_groups%numMyGroups
        if (num_posts > 0) then
           allocate(myRecvBuff(len*num_posts))
        end if
        call ovf_print_init(len, num_posts, myRecvBuff, post_array)
     else
        num_posts = 0
        do n=1, ovf_groups%numMyGroups
           if (ovf_groups%amMaster(n)) then
              num_posts = num_posts + 1
           end if
        end do
        call ovf_print_init(len, num_posts)
        if (num_posts > 0) then
           allocate(mySendBuff(num_posts * len))
        end if
     end if
  else ! no printing - just see which are my groups
     if( my_task == master_task ) then
        post_array = .true.
        do n=1, ovf_groups%numMyGroups
           ovf_id = ovf_groups%groupIds(n)
           post_array(ovf_id) = .false.
        end do
     end if
  end if
 
  !do master seperate for easier printing
  if( my_task == master_task ) then
     count = 0
      do n=1,num_ovf
         ovf_id =  n;
         post = post_array(n) ! post is true if master is not in group 
         if ( .not. post) then! master is part of this ovf

            ! set parameters
            lat     = ovf(ovf_id)%ovf_params%lat
            fs      = c2*omega*sin(lat*pi/180.0_r8)
            hu      = ovf(ovf_id)%ovf_params%source_thick
            hs      = hu*(c2/c3)
            xse     = ovf(ovf_id)%ovf_params%distnc_str_ssb
            alpha   = ovf(ovf_id)%ovf_params%bottom_slope
            cd      = ovf(ovf_id)%ovf_params%bottom_drag
            di      = p5*(zt(ovf(ovf_id)%reg_inf%kmin)+zt(ovf(ovf_id)%reg_inf%kmax))
            ds      = zt(ovf(ovf_id)%loc_src(1)%k)
            de      = zt(ovf(ovf_id)%loc_ent(1)%k)
            Ws      = ovf(ovf_id)%ovf_params%width
            ! set region T,S and compute densities
            T_i     = ovf(ovf_id)%trcr_reg%inf(1)
            S_i     = ovf(ovf_id)%trcr_reg%inf(2)
            call state_singlept(T_i,S_i,ds,rho_i)
            T_s     = ovf(ovf_id)%trcr_reg%src(1)
            S_s     = ovf(ovf_id)%trcr_reg%src(2)
            call state_singlept(T_s,S_s,ds,rho_s)
            call state_singlept(T_s,S_s,de,rho_sed)
            T_e     = ovf(ovf_id)%trcr_reg%ent(1)
            S_e     = ovf(ovf_id)%trcr_reg%ent(2)
            call state_singlept(T_e,S_e,de,rho_e)
            ! compute inflow/source reduced gravity and source transport
            gp_s    = grav*(rho_s-rho_i)/rho_sw
            ! if no source overflow, zero out transports
            if( gp_s > c0 ) then
               Ms      = gp_s*hu*hu/(c2*fs)
               As      = hs*Ws
               Us      = Ms/As
               ! compute overflow spreading and entrainment transport
               gp_e    = grav*(rho_sed-rho_e)/rho_sw
               ! zero entrainment transport if gp_e < 0
               if( gp_e > c0 ) then
                  Ugeo    = gp_e*alpha/fs
                  Uavg    = p5*(Us+Ugeo)
                  a       = fs*Ws/c2
                  b       = fs*Ws*hs/c2 + c2*cd*Uavg*xse - Ms*fs/(c2*Ugeo)
                  c       = -fs*Ms*hs/(c2*Ugeo)
                  hgeo    = (-b + sqrt(b*b-c4*a*c))/(c2*a)
                  Fgeo    = Ugeo/sqrt(gp_e*hgeo)
                  phi     = c1-Fgeo**(-c2/c3)
                  Me      = Ms*phi/(c1-phi)
                  ! zero entrainment transport if phi < c0
                  if( phi > c0 ) then
                     Mp  = Ms + Me
                  else
                     Me     = c0
                     Mp     = Ms
                  endif
               else
                  Me     = c0
                  Mp     = Ms
               endif
            else
               Ms     = c0
               Me     = c0
               Mp     = c0
            endif
            ! time shift transports and set output in ovf array
            ovf(ovf_id)%Ms_nm1 = ovf(ovf_id)%Ms_n
            ovf(ovf_id)%Ms_n   = ovf(ovf_id)%Ms
            ovf(ovf_id)%Me_nm1 = ovf(ovf_id)%Me_n
            ovf(ovf_id)%Me_n   = ovf(ovf_id)%Me
            ovf(ovf_id)%Mp_nm1 = ovf(ovf_id)%Mp_n
            ovf(ovf_id)%Mp_n   = ovf(ovf_id)%Mp
            ovf(ovf_id)%Ms     = Ms
            ovf(ovf_id)%Me     = Me
            ovf(ovf_id)%Mp     = Mp
            ! recompute phi based on actual transports
            phi = ovf(ovf_id)%Me / (ovf(ovf_id)%Mp + c1)
            ! if time averaging time step, include last time step 
            if( avg_ts ) then
               phi = (ovf(ovf_id)%Me_n + ovf(ovf_id)%Me) / (ovf(ovf_id)%Mp_n + ovf(ovf_id)%Mp + c1)
            endif
            ovf(ovf_id)%phi = phi
            ! compute product T,S 
            T_p        = T_s*(c1-phi) + T_e*phi
            S_p        = S_s*(c1-phi) + S_e*phi
            ovf(ovf_id)%Tp  = T_p
            ovf(ovf_id)%Sp  = S_p
            do nn=1,nt
               ovf(ovf_id)%trcr_adj%prd(nn) = ovf(ovf_id)%trcr_adj%src(nn) * (c1 - phi) &
                    + ovf(ovf_id)%trcr_adj%ent(nn) * phi
               ovf(ovf_id)%trcr_reg%prd(nn) = ovf(ovf_id)%trcr_reg%src(nn) * (c1 - phi) &
                    + ovf(ovf_id)%trcr_reg%ent(nn) * phi
            end do
            ! product set for insertion
            m = ovf(ovf_id)%prd_set
            
            if (eod) then
               k_p = (ovf(ovf_id)%adj_prd(m)%kmin+ovf(ovf_id)%adj_prd(m)%kmax)/2
               write(ovf_diag_unit,1234) tday,n,phi,1.e-12*Ms,1.e-12*Me,1.e-12*Mp,m,zt(k_p)/100.
1234           format(' ovf_tr: ',f7.1,1x,i2,25x,f7.4,2x,3(f7.4,1x),1x,i2,1x,f8.1)
               write(ovf_diag_unit,1235) tday, n,T_i,S_i*c1000,T_s,S_s*c1000,T_e,S_e*c1000,T_p,S_p*c1000
1235           format(' ovf_TS: ',f7.1,1x,i2,1x,8(f7.4,1x))         
               call shr_sys_flush(ovf_diag_unit)
            endif ! eod - so print

         else !not in this ovf group - but if check if we can print 
            if (eod) then
               count = count+1
               call ovf_print_get(count)
               loc = (n-1)*len+1
               !unpack
               r_n = int(myRecvBuff(loc+1))
               r_m = int(myRecvBuff(loc+6))
               
               write(ovf_diag_unit,1234) myRecvBuff(loc),r_n, myRecvBuff(loc+2), &
                    myRecvBuff(loc+3),myRecvBuff(loc+4),myRecvBuff(loc+5), &
                    r_m, myRecvBuff(loc+7)
1236           format(' ovf_tr: ',f7.1,1x,i2,25x,f7.4,2x,3(f7.4,1x),1x,i2,1x,f8.1)
               write(ovf_diag_unit,1235) myRecvBuff(loc), r_n, myRecvBuff(loc+8), &
                    myRecvBuff(loc+9),myRecvBuff(loc+10), myRecvBuff(loc+11), &
                    myRecvBuff(loc+12), myRecvBuff(loc+13), myRecvBuff(loc+14), &
                    myRecvBuff(loc+15)
1237           format(' ovf_TS: ',f7.1,1x,i2,1x,8(f7.4,1x))         
               call shr_sys_flush(ovf_diag_unit)
            end if ! eod - so print
         end if
      end do ! ovf loop
      deallocate(post_array)
      if (eod) then
         if (num_posts > 0) then
            deallocate(myRecvBuff)
         endif
         call ovf_print_finalize(num_posts)
      end if
  else ! not master
     count = 0
     do n=1,ovf_groups%numMyGroups
        ovf_id = ovf_groups%groupIds(n)
        
        ! set parameters
        lat     = ovf(ovf_id)%ovf_params%lat
        fs      = c2*omega*sin(lat*pi/180.0_r8)
        hu      = ovf(ovf_id)%ovf_params%source_thick
        hs      = hu*(c2/c3)
        xse     = ovf(ovf_id)%ovf_params%distnc_str_ssb
        alpha   = ovf(ovf_id)%ovf_params%bottom_slope
        cd      = ovf(ovf_id)%ovf_params%bottom_drag
        di      = p5*(zt(ovf(ovf_id)%reg_inf%kmin)+zt(ovf(ovf_id)%reg_inf%kmax))
        ds      = zt(ovf(ovf_id)%loc_src(1)%k)
        de      = zt(ovf(ovf_id)%loc_ent(1)%k)
        Ws      = ovf(ovf_id)%ovf_params%width
        ! set region T,S and compute densities
        T_i     = ovf(ovf_id)%trcr_reg%inf(1)
        S_i     = ovf(ovf_id)%trcr_reg%inf(2)
        call state_singlept(T_i,S_i,ds,rho_i)
        T_s     = ovf(ovf_id)%trcr_reg%src(1)
        S_s     = ovf(ovf_id)%trcr_reg%src(2)
        call state_singlept(T_s,S_s,ds,rho_s)
        call state_singlept(T_s,S_s,de,rho_sed)
        T_e     = ovf(ovf_id)%trcr_reg%ent(1)
        S_e     = ovf(ovf_id)%trcr_reg%ent(2)
        call state_singlept(T_e,S_e,de,rho_e)
        ! compute inflow/source reduced gravity and source transport
        gp_s    = grav*(rho_s-rho_i)/rho_sw
        ! if no source overflow, zero out transports
        if( gp_s > c0 ) then
           Ms      = gp_s*hu*hu/(c2*fs)
           As      = hs*Ws
           Us      = Ms/As
           ! compute overflow spreading and entrainment transport
           gp_e    = grav*(rho_sed-rho_e)/rho_sw
           ! zero entrainment transport if gp_e < 0
           if( gp_e > c0 ) then
              Ugeo    = gp_e*alpha/fs
              Uavg    = p5*(Us+Ugeo)
              a       = fs*Ws/c2
              b       = fs*Ws*hs/c2 + c2*cd*Uavg*xse - Ms*fs/(c2*Ugeo)
              c       = -fs*Ms*hs/(c2*Ugeo)
              hgeo    = (-b + sqrt(b*b-c4*a*c))/(c2*a)
              Fgeo    = Ugeo/sqrt(gp_e*hgeo)
              phi     = c1-Fgeo**(-c2/c3)
              Me      = Ms*phi/(c1-phi)
              ! zero entrainment transport if phi < c0
              if( phi > c0 ) then
                 Mp  = Ms + Me
              else
                 Me     = c0
                 Mp     = Ms
              endif
           else
              Me     = c0
              Mp     = Ms
           endif
        else
           Ms     = c0
           Me     = c0
           Mp     = c0
        endif
        ! time shift transports and set output in ovf array
        ovf(ovf_id)%Ms_nm1 = ovf(ovf_id)%Ms_n
        ovf(ovf_id)%Ms_n   = ovf(ovf_id)%Ms
        ovf(ovf_id)%Me_nm1 = ovf(ovf_id)%Me_n
        ovf(ovf_id)%Me_n   = ovf(ovf_id)%Me
        ovf(ovf_id)%Mp_nm1 = ovf(ovf_id)%Mp_n
        ovf(ovf_id)%Mp_n   = ovf(ovf_id)%Mp
        ovf(ovf_id)%Ms     = Ms
        ovf(ovf_id)%Me     = Me
        ovf(ovf_id)%Mp     = Mp
        ! recompute phi based on actual transports
        phi = ovf(ovf_id)%Me / (ovf(ovf_id)%Mp + c1)
        ! if time averaging time step, include last time step 
        if( avg_ts ) then
           phi = (ovf(ovf_id)%Me_n + ovf(ovf_id)%Me) / (ovf(ovf_id)%Mp_n + ovf(ovf_id)%Mp + c1)
        endif
        ovf(ovf_id)%phi = phi
        ! compute product T,S 
        T_p        = T_s*(c1-phi) + T_e*phi
        S_p        = S_s*(c1-phi) + S_e*phi
        ovf(ovf_id)%Tp  = T_p
        ovf(ovf_id)%Sp  = S_p
        do nn=1,nt
           ovf(ovf_id)%trcr_adj%prd(nn) = ovf(ovf_id)%trcr_adj%src(nn) * (c1 - phi) &
                + ovf(ovf_id)%trcr_adj%ent(nn) * phi
           ovf(ovf_id)%trcr_reg%prd(nn) = ovf(ovf_id)%trcr_reg%src(nn) * (c1 - phi) &
                + ovf(ovf_id)%trcr_reg%ent(nn) * phi
        end do
        ! product set for insertion
        m = ovf(ovf_id)%prd_set
        !if iam the group master, need to send the info to be printed
        if (eod .and. ovf_groups%amMaster(n))then

           count = count +1
           loc = (count-1)*len+1
           
           k_p = (ovf(ovf_id)%adj_prd(m)%kmin+ovf(ovf_id)%adj_prd(m)%kmax)/2
           mySendBuff(loc) = tday
           mySendBuff(loc+1) = real(n,r8)
           mySendBuff(loc+2)  = phi
           mySendBuff(loc+3)  = 1.e-12*Ms
           mySendBuff(loc+4)  = 1.e-12*Me
           mySendBuff(loc+5)  = 1.e-12*Mp
           mySendBuff(loc+6)  = real(m,r8)
           mySendBuff(loc+7)  = zt(k_p)/100.
           
           mySendBuff(loc+8)  = T_i
           mySendBuff(loc+9)  = S_i*c1000
           mySendBuff(loc+10)  = T_s
           mySendBuff(loc+11)  = S_s*c1000
           mySendBuff(loc+12)  = T_e
           mySendBuff(loc+13)  = S_e*c1000
           mySendBuff(loc+14)  = T_p
           mySendBuff(loc+15)  = S_p*c1000
           

           call ovf_print_send(len, mySendBuff(loc:loc+len-1), &
                count, ovf_id)
        endif ! printing
     end do   ! n loop over overflows
     if (eod) then
        call ovf_print_finalize(num_posts)
        if (num_posts > 0) then
           deallocate(mySendBuff)
        end if
     end if
  end if !not master_task



!-----------------------------------------------------------------------
!
! close overflows_diag_outfile file 
!
!-----------------------------------------------------------------------

   if (eod .and. my_task == master_task) then
      close(ovf_diag_unit)
   endif ! print_overflows_diag

!----------------------------------------------------------------------
!EOC

 end subroutine ovf_transports

!***********************************************************************
!EOP
! !IROUTINE: ovf_loc_prd
! !INTERFACE:

 subroutine ovf_loc_prd

! !DESCRIPTION:
!  Evaluate the ovf location of product. If product location has moved,
!  set original sidewall velocities on the ugrid to zero and compute 
!  Uovf_n, Uovf_nm1 sidewall velocities on the u-grid at new product 
!  location using Mp_n, Mp_nm1 transport respectively.
!
! !REVISION HISTORY:
!  same as module
!
! ovf t-grid box ij with u-grid
! corners and orientations          
! product moves out of ij ovf box  
! 
!                                   ^
!                                +V | 2       
!                                   | __________ ---> +U
!              y ^                   |          | 1
!                |                   |          |
!                |                   |    ij    | 
!                +----->             |          |
!                      x           3 |__________|
!                             -U <---            |
!                                             4  | -V

!EOP
!BOC
!----------------------------------------------------------------------
!  local variables
!----------------------------------------------------------------------

   integer (int_kind)  :: &
      m_neut_org,         & ! original neutral product density set index
      m_neut                ! neutral product density set index
   real (r8)           :: &
      T_p                ,& ! product temperature (C)
      S_p                ,& ! product salinity
      rho_p              ,& ! product density at each product
      ufrc               ,& ! fraction of ovf velocity for each box
      Uovf_n             ,&   ! U at n
      Uovf_nm1                ! U at n-1
   integer (int_kind)  :: &
      iblock,i,j,n,m,mp,  & ! dummy loop indices
      ib,ie,jb,je,        & ! local domain index boundaries
      k_p,kprd              ! overflow loop and level indices
   type (block)        :: &
      this_block            ! block information for current block
   logical (log_kind), parameter :: prnt = .false.

   integer (int_kind)  :: ovf_id


   integer (int_kind) ::  num_posts, ierr, count, loc, &
        loc_start, num_len, tot_len
   integer (int_kind), parameter :: len = 8 !length of buffer for a post
   real (r8), dimension(:), allocatable:: myRecvBuff(:), mySendBuff(:)
   real (r8) :: rb2, rb3, rb4, rb5, rb6, rb7, rb8
   integer (i4):: ib1, ib4, ib5, indx


   logical (log_kind) :: post, found
   logical (log_kind), dimension(:), allocatable :: post_array

   if(prnt .and. my_task == master_task) then
      write(stdout,*) 'ovf_loc_prd called '
      call shr_sys_flush(stdout)
   endif

!-----------------------------------------------------------------------
!  overflow loop - loop through MY overflows
!-----------------------------------------------------------------------
   


 !must have the groups initialized
   if (.not. ovf_groups%init) then
      call ovf_init_groups()
   end if

   if( my_task == master_task ) then
      allocate(post_array(num_ovf))
   end if

   !set up for printing
   if (prnt) then 
      count = 0
      do n=1, num_ovf
         count = MAX(count, ovf(n)%num_prd_sets)
      end do
      tot_len = len*count + 1 !one spot to tell how long it is
      if( my_task == master_task ) then
         !for any ovf groups that the master_task is not in, 
         !must post recv for the info
         num_posts = num_ovf - ovf_groups%numMyGroups
         if (num_posts > 0) then
            allocate(myRecvBuff(tot_len*num_posts))
         endif
         call ovf_print_init(tot_len, num_posts, myRecvBuff, post_array)
      else !not master
         num_posts = 0
         do n=1, ovf_groups%numMyGroups
            if (ovf_groups%amMaster(n)) then
               num_posts = num_posts + 1
            end if
         end do
         call ovf_print_init(tot_len, num_posts)
         if (num_posts > 0) then
            allocate(mySendBuff(num_posts * tot_len))
         end if
      end if
   else ! not printing
      if( my_task == master_task ) then
         post_array = .true.
         do n=1, ovf_groups%numMyGroups
            ovf_id = ovf_groups%groupIds(n)
            post_array(ovf_id) = .false.
         end do
      end if
   end if


   !seperate the master and slave loops to make printing easier - 
   ! this ends up in some code duplication in the two sections...
   if( my_task == master_task ) then
      count = 0
      do n=1,num_ovf !master goes through all because of possible printing
         ovf_id =  n;
         post = post_array(n) ! post is true if master is not in group 
         if ( .not. post) then! master is part of this ovf
            ! find new product location
            T_p = ovf(ovf_id)%Tp
            S_p = ovf(ovf_id)%Sp
            m_neut_org = ovf(ovf_id)%prd_set
            m_neut     = 0
            if(ovf(ovf_id)%num_prd_sets .eq. 1) then
               m_neut = 1
               k_p = (ovf(ovf_id)%adj_prd(1)%kmin + ovf(ovf_id)%adj_prd(1)%kmax)/2
               call state_singlept(T_p,S_p,zt(k_p),rho_p)
            else
               ! search from deepest to shallowest to allow product water
               ! to go to the deepest possible level
               do m=ovf(ovf_id)%num_prd_sets-1,1,-1
                  k_p = (ovf(ovf_id)%adj_prd(m)%kmin + ovf(ovf_id)%adj_prd(m)%kmax)/2
                  ! get product level for this set
                  
                  !this calculates rho_p based on temp, salinity, and depth
                  call state_singlept(T_p,S_p,zt(k_p),rho_p)

                  if (prnt) then
                     if (m == 1) then
                        write(stdout,5) m, 0.0, &
                             (ovf(ovf_id)%rho_adj%prd(m)-c1)*c1000, &
                             k_p,T_p,S_p,zt(k_p),(rho_p-c1)*c1000
                        else
                           write(stdout,5) m,(ovf(ovf_id)%rho_adj%prd(m-1)-c1)*c1000, &
                                (ovf(ovf_id)%rho_adj%prd(m)-c1)*c1000, &
                                k_p,T_p,S_p,zt(k_p),(rho_p-c1)*c1000
                        end if
5                    format(' neutral lev search- m rho_adj_m-1 rho_adj_m ', &
                          'k_p T_p S_p zt(k_p) rho_p =',/ &
                          2x,i2,2x,2(f12.8,2x),4x,i2,4(f12.8,2x))
                  end if
                  if(rho_p .gt. ovf(ovf_id)%rho_adj%prd(m)) then
                     m_neut = m+1
                     goto 999 !exit loop
                  else
                     m_neut = m
                  endif
               enddo
999            continue
            endif
            ! error check
            if( m_neut .eq. 0 ) then
               write(stdout,10) T_p,S_p,rho_p
10             format(' ovf_loc_prd: no prd lev found for T,S,rho=', &
                    3(f10.5,2x))
               call shr_sys_flush(stdout)
               call exit_POP(sigAbort,'ERROR no product level found')
            endif
            ovf(ovf_id)%prd_set_n = m_neut_org
            ovf(ovf_id)%prd_set   = m_neut
            if (prnt) then
               write(stdout,20) n,T_p,S_p*c1000,(rho_p-c1)*c1000,m_neut
20             format(' For ovf = ',i3,' prd T,S,rho = ',3(f12.8,2x),' prd set =',i5)
            end if
            if( m_neut_org .ne. 0 .and. m_neut_org .ne. m_neut ) then
               ! product point has moved
               if ( overflows_on) then
                  write(stdout,*) 'ovf_loc_prd: nsteps_total=',nsteps_total, &
                       ' ovf=',ovf_id,' swap ovf UV old/new ', &
                       'prd set old/new=',m_neut_org,m_neut
                  call shr_sys_flush(stdout)
               endif
               ! compute Uovf_n, Uovf_nm1 velocities for product sidewall
               m = ovf(ovf_id)%prd_set  ! product set for insertion
               do mp=1,ovf(ovf_id)%num_prd(m)  ! product points for each set
                  kprd = ovf(ovf_id)%loc_prd(m,mp)%k
                  ufrc = c1/real(ovf(ovf_id)%num_prd(m)-1)
                  do iblock = 1,nblocks_clinic
                     this_block = get_block(blocks_clinic(iblock),iblock)
                     ib = this_block%ib
                     ie = this_block%ie
                     jb = this_block%jb
                     je = this_block%je
                     do j=jb,je
                        if( ovf(ovf_id)%loc_prd(m,mp)%j_u .eq. this_block%j_glob(j) ) then
                           do i=ib,ie
                              if( ovf(ovf_id)%loc_prd(m,mp)%i_u .eq. this_block%i_glob(i) ) then
                                 if( ovf(ovf_id)%loc_prd(m,mp)%orient .eq. 1 ) then   
                                    Uovf_nm1 = ovf(ovf_id)%Mp_nm1*ufrc/(dz(kprd)*DYU(i,j,iblock))
                                    Uovf_n   = ovf(ovf_id)%Mp_n  *ufrc/(dz(kprd)*DYU(i,j,iblock))
                                 endif
                                 if( ovf(ovf_id)%loc_prd(m,mp)%orient .eq. 2 ) then  
                                    Uovf_nm1 = ovf(ovf_id)%Mp_nm1*ufrc/(dz(kprd)*DXU(i,j,iblock))
                                    Uovf_n   = ovf(ovf_id)%Mp_n  *ufrc/(dz(kprd)*DXU(i,j,iblock))
                                 endif
                                 if( ovf(ovf_id)%loc_prd(m,mp)%orient .eq. 3 ) then  
                                    Uovf_nm1 = ovf(ovf_id)%Mp_nm1*ufrc/(dz(kprd)*DYU(i,j,iblock))
                                    Uovf_n   = ovf(ovf_id)%Mp_n  *ufrc/(dz(kprd)*DYU(i,j,iblock))
                                 endif
                                 if( ovf(ovf_id)%loc_prd(m,mp)%orient .eq. 4 ) then  
                                    Uovf_nm1 = ovf(ovf_id)%Mp_nm1*ufrc/(dz(kprd)*DXU(i,j,iblock))
                                    Uovf_n   = ovf(ovf_id)%Mp_n  *ufrc/(dz(kprd)*DXU(i,j,iblock))
                                 endif
                                 ovf(ovf_id)%loc_prd(m,mp)%Uovf_nm1 = Uovf_nm1
                                 ovf(ovf_id)%loc_prd(m,mp)%Uovf_n   = Uovf_n
                                 if(prnt) then
                                    write(stdout,30) ovf(ovf_id)%loc_prd(m,mp)%i,ovf(ovf_id)%loc_prd(m,mp)%j, &
                                         ovf(ovf_id)%loc_prd(m,mp)%k,ovf(ovf_id)%Mp_nm1,ufrc,dz(kprd),Uovf_nm1
30                                  format(' loc_prd ijk=',3(i4,1x),'Mp_nm1 uf dz=',3(1pe10.3,1x), &
                                         'Uovf_nm1=',1pe10.3)
                                 endif
                              endif
                           end do  ! i
                        endif
                     end do  ! j
                  end do  ! iblock
               end do  ! product points for each set
               if( overflows_interactive ) then
                  ! zero out original product sidewall U
                  m = m_neut_org
                  do mp=1,ovf(ovf_id)%num_prd(m)  ! product points for each set
                     ! prd  set original Uold sidewalls to zero
                     kprd = ovf(ovf_id)%loc_prd(m,mp)%k
                     do iblock = 1,nblocks_clinic
                        this_block = get_block(blocks_clinic(iblock),iblock)
                        ib = this_block%ib
                        ie = this_block%ie
                        jb = this_block%jb
                        je = this_block%je
                        do j=jb,je
                           if( ovf(ovf_id)%loc_prd(m,mp)%j_u .eq. this_block%j_glob(j) ) then
                              do i=ib,ie
                                 if( ovf(ovf_id)%loc_prd(m,mp)%i_u .eq. this_block%i_glob(i) ) then
                                    UVEL(i,j,kprd,newtime,iblock) = c0
                                    VVEL(i,j,kprd,newtime,iblock) = c0
                                 endif
                              end do  ! i
                           endif
                        end do  ! j
                     end do  ! iblock
                  end do  ! product points for each set
               endif   ! interactive overflows
            endif   ! product point has moved

         else !master not in this ovf - but may need to print (?)
            if (prnt) then
               count = count+1
               call ovf_print_get(count)
               loc = (count-1)*tot_len+1
               num_len = myRecvBuff(count*tot_len)
               !unpack and print
               do i = 1, num_len - 1
                  ib1 = int(myRecvBuff(loc)) 
                  rb2 = myRecvBuff(loc+1)
                  rb3 = myRecvBuff(loc+2)
                  ib4 = int(myRecvBuff(loc+3))
                  rb5 = myRecvBuff(loc+4)
                  rb6 = myRecvBuff(loc+5)
                  rb7 = myRecvBuff(loc+6)
                  rb8 = myRecvBuff(loc+7)
                  write(stdout,50) ib1, rb2, rb3, ib4, rb5, rb6, rb7, rb8
50                format(' neutral lev search- m rho_adj_m-1 rho_adj_m ', &
                       'k_p T_p S_p zt(k_p) rho_p =',/ &
                       2x,i2,2x,2(f12.8,2x),4x,i2,4(f12.8,2x))
                  loc = loc + len !advance for next one
               end do ! i
               ib1 = int(myRecvBuff(loc)) 
               rb2 = myRecvBuff(loc+1)
               rb3 = myRecvBuff(loc+2)
               rb4 = myRecvBuff(loc+3)
               ib5 = int( myRecvBuff(loc+4))
               write(stdout,200) ib1, rb2, rb3, rb4, ib5
200            format(' For ovf = ',i3,' prd T,S,rho = ',3(f12.8,2x),' prd set =',i5)
            end if !print
         end if !end of master not in ovf 
      end do !ovf loop
      !clean up
      deallocate(post_array)
      if (prnt) then
         if (num_posts > 0) then
            deallocate(myRecvBuff)
         endif
         call ovf_print_finalize(num_posts)
      end if

   else ! not master loop
      num_posts = 0
      do n=1, ovf_groups%numMyGroups
         loc_start = -1
         num_len= 0
         ovf_id = ovf_groups%groupIds(n)
         ! find new product location
         T_p = ovf(ovf_id)%Tp
         S_p = ovf(ovf_id)%Sp
         m_neut_org = ovf(ovf_id)%prd_set
         m_neut     = 0
         if(ovf(ovf_id)%num_prd_sets .eq. 1) then
            m_neut = 1
            k_p = (ovf(ovf_id)%adj_prd(1)%kmin + ovf(ovf_id)%adj_prd(1)%kmax)/2
            call state_singlept(T_p,S_p,zt(k_p),rho_p)
         else
            ! search from deepest to shallowest to allow product water
            ! to go to the deepest possible level
            do m=ovf(ovf_id)%num_prd_sets-1,1,-1
               k_p = (ovf(ovf_id)%adj_prd(m)%kmin + ovf(ovf_id)%adj_prd(m)%kmax)/2
               ! get product level for this set
               
               !this calculates rho_p based on temp, salinity, and depth
               call state_singlept(T_p,S_p,zt(k_p),rho_p)
               
               if (prnt .and. ovf_groups%amMaster(n)) then
                  
                  num_len = num_len + 1
                  loc = (num_posts)*tot_len + (num_len - 1)*len + 1
                  
                  if ( loc_start < 0) then !first time for this ovf
                     loc_start = loc
                  endif
                  mySendBuff(loc) = real(m, r8)  !int
                  if (m==1) then
                     mySendBuff(loc+1) =0.0 !can't do m-1
                  else
                     mySendBuff(loc+1) = (ovf(ovf_id)%rho_adj%prd(m-1)-c1)*c1000
                  end if
                  mySendBuff(loc+2) = (ovf(ovf_id)%rho_adj%prd(m)-c1)*c1000
                  mySendBuff(loc+3) = real(k_p,r8) !int
                  mySendBuff(loc+4) = T_p
                  mySendBuff(loc+5) = S_p
                  mySendBuff(loc+6) = zt(k_p)
                  mySendBuff(loc+7) = (rho_p-c1)*c1000
                  
               endif ! print master
  
               if(rho_p .gt. ovf(ovf_id)%rho_adj%prd(m)) then
                  m_neut = m+1
                  goto 998 !exit loop
               else
                  m_neut = m
               endif
            enddo
998         continue
         endif !num prod sets
         !error check
         if( m_neut .eq. 0 ) then
            write(stdout,101) T_p,S_p,rho_p
101         format(' ovf_loc_prd: no prd lev found for T,S,rho=', &
                 3(f10.5,2x))
            call shr_sys_flush(stdout)
            call exit_POP(sigAbort,'ERROR no product level found')
         endif
         ovf(ovf_id)%prd_set_n = m_neut_org
         ovf(ovf_id)%prd_set   = m_neut
         if (prnt .and. ovf_groups%amMaster(n)) then
            num_len = num_len+1
            loc = (num_posts)*tot_len + (num_len - 1)*len + 1
            
            if (loc_start  == -1) then !first time for this ovf
               loc_start = loc;
            endif
            !pack pack buffer and send here
            mySendBuff(loc) = real(ovf_id, r8) !int
            mySendBuff(loc+1) =T_p
            mySendBuff(loc+2) =S_p*c1000
            mySendBuff(loc+3) = (rho_p-c1)*c1000
            mySendBuff(loc+4) = real(m_neut, r8) !int
            
            num_posts = num_posts + 1 !increment here because we do this one once each ovf
            
            !indicate num items sent to master in last location
            mySendBuff(tot_len*num_posts) = real(num_len, r8)
            
            call ovf_print_send(tot_len, mySendBuff(loc_start:loc_start+tot_len-1), &
                 num_posts, ovf_id)
            
         end if ! print master
         
         if( m_neut_org .ne. 0 .and. m_neut_org .ne. m_neut ) then
            ! product point has moved
            if ( overflows_on .and. ovf_groups%amMaster(n)) then
               write(stdout,*) 'ovf_loc_prd: nsteps_total=',nsteps_total, &
                    ' ovf=',ovf_id,' swap ovf UV old/new ', &
                    'prd set old/new=',m_neut_org,m_neut
               call shr_sys_flush(stdout)
            endif
            ! compute Uovf_n, Uovf_nm1 velocities for product sidewall
            m = ovf(ovf_id)%prd_set  ! product set for insertion
            do mp=1,ovf(ovf_id)%num_prd(m)  ! product points for each set
               kprd = ovf(ovf_id)%loc_prd(m,mp)%k
               ufrc = c1/real(ovf(ovf_id)%num_prd(m)-1)
               do iblock = 1,nblocks_clinic
                  this_block = get_block(blocks_clinic(iblock),iblock)
                  ib = this_block%ib
                  ie = this_block%ie
                  jb = this_block%jb
                  je = this_block%je
                  do j=jb,je
                     if( ovf(ovf_id)%loc_prd(m,mp)%j_u .eq. this_block%j_glob(j) ) then
                        do i=ib,ie
                           if( ovf(ovf_id)%loc_prd(m,mp)%i_u .eq. this_block%i_glob(i) ) then
                              if( ovf(ovf_id)%loc_prd(m,mp)%orient .eq. 1 ) then   
                                 Uovf_nm1 = ovf(ovf_id)%Mp_nm1*ufrc/(dz(kprd)*DYU(i,j,iblock))
                                 Uovf_n   = ovf(ovf_id)%Mp_n  *ufrc/(dz(kprd)*DYU(i,j,iblock))
                              endif
                              if( ovf(ovf_id)%loc_prd(m,mp)%orient .eq. 2 ) then  
                                 Uovf_nm1 = ovf(ovf_id)%Mp_nm1*ufrc/(dz(kprd)*DXU(i,j,iblock))
                                 Uovf_n   = ovf(ovf_id)%Mp_n  *ufrc/(dz(kprd)*DXU(i,j,iblock))
                              endif
                              if( ovf(ovf_id)%loc_prd(m,mp)%orient .eq. 3 ) then  
                                 Uovf_nm1 = ovf(ovf_id)%Mp_nm1*ufrc/(dz(kprd)*DYU(i,j,iblock))
                                 Uovf_n   = ovf(ovf_id)%Mp_n  *ufrc/(dz(kprd)*DYU(i,j,iblock))
                              endif
                              if( ovf(ovf_id)%loc_prd(m,mp)%orient .eq. 4 ) then  
                                 Uovf_nm1 = ovf(ovf_id)%Mp_nm1*ufrc/(dz(kprd)*DXU(i,j,iblock))
                                 Uovf_n   = ovf(ovf_id)%Mp_n  *ufrc/(dz(kprd)*DXU(i,j,iblock))
                              endif
                              ovf(ovf_id)%loc_prd(m,mp)%Uovf_nm1 = Uovf_nm1
                              ovf(ovf_id)%loc_prd(m,mp)%Uovf_n   = Uovf_n
                              if(prnt) then
                                 write(stdout,31) ovf(ovf_id)%loc_prd(m,mp)%i,ovf(ovf_id)%loc_prd(m,mp)%j, &
                                      ovf(ovf_id)%loc_prd(m,mp)%k,ovf(ovf_id)%Mp_nm1,ufrc,dz(kprd),Uovf_nm1
31                               format(' loc_prd ijk=',3(i4,1x),'Mp_nm1 uf dz=',3(1pe10.3,1x), &
                                      'Uovf_nm1=',1pe10.3)
                              endif
                           endif
                        end do  ! i
                     endif
                  end do  ! j
               end do  ! iblock
            end do  ! product points for each set
            if( overflows_interactive ) then
               ! zero out original product sidewall U
               m = m_neut_org
               do mp=1,ovf(ovf_id)%num_prd(m)  ! product points for each set
                  ! prd  set original Uold sidewalls to zero
                  kprd = ovf(ovf_id)%loc_prd(m,mp)%k
                  do iblock = 1,nblocks_clinic
                     this_block = get_block(blocks_clinic(iblock),iblock)
                     ib = this_block%ib
                     ie = this_block%ie
                     jb = this_block%jb
                     je = this_block%je
                     do j=jb,je
                        if( ovf(ovf_id)%loc_prd(m,mp)%j_u .eq. this_block%j_glob(j) ) then
                           do i=ib,ie
                              if( ovf(ovf_id)%loc_prd(m,mp)%i_u .eq. this_block%i_glob(i) ) then
                                 UVEL(i,j,kprd,newtime,iblock) = c0
                                 VVEL(i,j,kprd,newtime,iblock) = c0
                              endif
                           end do  ! i
                        endif
                     end do  ! j
                  end do  ! iblock
               end do  ! product points for each set
            endif   ! interactive overflows
         endif   ! product point has moved
      end do  ! overflows

      if (prnt) then
         call ovf_print_finalize(num_posts)
         if (num_posts > 0) then
            deallocate(mySendBuff)
         end if
      end if

   end if ! end of not master part


!----------------------------------------------------------------------
!EOC

 end subroutine ovf_loc_prd

!***********************************************************************
!EOP
! !IROUTINE: ovf_W
! !INTERFACE:

 subroutine ovf_W

! !DESCRIPTION:
!  Evaluate ovf vertical velocity W on the t-grid
!
! !REVISION HISTORY:
!  same as module
!
! ovf t-grid box ij from top
! sidewall transports and top TAREA
! used to compute Wovf at top
! signs of Wovf important!
!                                  
!                     __________
!                    |          |
!          Me  ----> |    ij    | <---- Ms
!          Mp  <---- |  TAREAij |
!                    |__________|
!                                               

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind)  :: &
      iblock,i,j,n,m,mp,  & ! dummy loop indices
      ib,ie,jb,je           ! local domain index boundaries
   type (block)        :: &
      this_block            ! block information for current block
   real (r8)           :: & 
      ufrc                  ! fraction of ovf velocity for each box
   logical (log_kind), parameter :: prnt = .false.

   integer (int_kind)  :: ovf_id 

   if( prnt .and. my_task == master_task ) then
      write(stdout,*) 'ovf_W called '
      call shr_sys_flush(stdout)
   endif

!-----------------------------------------------------------------------
!  overflow loop
!-----------------------------------------------------------------------
 !must have the groups initilaized
   if (.not. ovf_groups%init) then
      call ovf_init_groups()
   end if  

   do n=1, ovf_groups%numMyGroups  ! each overflow
      ovf_id = ovf_groups%groupIds(n)
! src
      do m=1,ovf(ovf_id)%num_src  ! source
         ufrc = c1/real(ovf(ovf_id)%num_src-1)
         if(m==1 .or. m==ovf(ovf_id)%num_src) ufrc = ufrc/c2
         do iblock = 1,nblocks_clinic
            this_block = get_block(blocks_clinic(iblock),iblock)
            ib = this_block%ib
            ie = this_block%ie
            jb = this_block%jb
            je = this_block%je
            do j=jb,je
               if( ovf(ovf_id)%loc_src(m)%j .eq. this_block%j_glob(j) ) then
                  do i=ib,ie
                     if( ovf(ovf_id)%loc_src(m)%i .eq. this_block%i_glob(i) ) then
                        ovf(ovf_id)%loc_src(m)%Wovf = -abs(ovf(ovf_id)%Ms*ufrc &
                                                      /TAREA(i,j,iblock))
                        if(prnt) then
                           write(stdout,10) n,ovf(ovf_id)%loc_src(m)%i, &
                           ovf(ovf_id)%loc_src(m)%j,ovf(ovf_id)%loc_src(m)%k, &
                           ovf(ovf_id)%Ms,ufrc,TAREA(i,j,iblock), &
                           ovf(ovf_id)%loc_src(m)%Wovf
                           10 format(' ovf_W n=',i3,' src ijk=',3(i4,1x), &
                           'Ms uf Ta Wovf=',4(1pe10.3,1x)) 
                        endif
                     endif
                  end do  ! i
               endif
            end do  ! j
         end do  ! iblock
      end do  ! source
! ent
      do m=1,ovf(ovf_id)%num_ent  ! entrainment
         ufrc = c1/real(ovf(ovf_id)%num_ent-1)
         if(m==1 .or. m==ovf(ovf_id)%num_ent) ufrc = ufrc/c2
         do iblock = 1,nblocks_clinic
            this_block = get_block(blocks_clinic(iblock),iblock)
            ib = this_block%ib
            ie = this_block%ie
            jb = this_block%jb
            je = this_block%je
            do j=jb,je
               if( ovf(ovf_id)%loc_ent(m)%j .eq. this_block%j_glob(j) ) then
                  do i=ib,ie
                     if( ovf(ovf_id)%loc_ent(m)%i .eq. this_block%i_glob(i) ) then
                        ovf(ovf_id)%loc_ent(m)%Wovf = -abs(ovf(ovf_id)%Me*ufrc &
                                                      /TAREA(i,j,iblock))
                        if(prnt) then
                           write(stdout,20) n,ovf(ovf_id)%loc_ent(m)%i, &
                           ovf(ovf_id)%loc_ent(m)%j,ovf(ovf_id)%loc_ent(m)%k, &
                           ovf(ovf_id)%Me,ufrc,TAREA(i,j,iblock), &
                           ovf(ovf_id)%loc_ent(m)%Wovf
                           20 format(' ovf_W n=',i3,' ent ijk=',3(i4,1x), &
                           'Me uf Ta Wovf=',4(1pe10.3,1x)) 
                        endif
                     endif
                  end do  ! i
               endif
            end do  ! j
         end do  ! iblock
      end do  ! entrainment
! prd
! set Wovf terms to zero at product points, incase product has moved
      do m=1,ovf(ovf_id)%num_prd_sets
        do mp=1,ovf(ovf_id)%num_prd(m)
          ovf(ovf_id)%loc_prd(m,mp)%Wovf = c0
        end do
      end do
      m = ovf(ovf_id)%prd_set  ! product set for insertion
      do mp=1,ovf(ovf_id)%num_prd(m)  ! product points for each set
         ufrc = c1/real(ovf(ovf_id)%num_prd(m)-1)
         if(mp==1 .or. mp==ovf(ovf_id)%num_prd(m)) ufrc = ufrc/c2
         do iblock = 1,nblocks_clinic
            this_block = get_block(blocks_clinic(iblock),iblock)
            ib = this_block%ib
            ie = this_block%ie
            jb = this_block%jb
            je = this_block%je
            do j=jb,je
               if( ovf(ovf_id)%loc_prd(m,mp)%j .eq. this_block%j_glob(j) ) then
                  do i=ib,ie
                     if( ovf(ovf_id)%loc_prd(m,mp)%i .eq. this_block%i_glob(i) ) then
                        ovf(ovf_id)%loc_prd(m,mp)%Wovf = abs(ovf(ovf_id)%Mp*ufrc &
                                                        /TAREA(i,j,iblock))
                        if(prnt) then
                           write(stdout,30) n,ovf(ovf_id)%loc_prd(m,mp)%i, &
                           ovf(ovf_id)%loc_prd(m,mp)%j,ovf(ovf_id)%loc_prd(m,mp)%k, &
                           ovf(ovf_id)%Mp,ufrc,TAREA(i,j,iblock), &
                           ovf(ovf_id)%loc_prd(m,mp)%Wovf
                           30 format(' ovf_W n=',i3,' prd ijk=',3(i4,1x), &
                           'Mp uf Ta Wovf=',4(1pe10.3,1x)) 
                        endif
                     endif
                  end do  ! i
               endif
            end do  ! j
         end do  ! iblock
      end do  ! product points for insertion set
   end do  ! each overflow
!----------------------------------------------------------------------
!EOC

 end subroutine ovf_W

!***********************************************************************
!EOP
! !IROUTINE: ovf_UV
! !INTERFACE:

 subroutine ovf_UV

! !DESCRIPTION:
!  Evaluate the ovf sidewall velocities UV on the u-grid
! !REVISION HISTORY:
!  same as module
!
! ovf t-grid box ij with U on u-grid
! at corner set by orientation
!
!                                     2
!                                   U __________ U
!              y ^                   |          | 1
!                |                   |          |
!                |                   |    ij    |
!                +----->             |          |
!                      x           3 |__________|
!                                   U            U
!                                               4

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      iblock,i,j,n,m,mp, & ! dummy loop indices
      ib,ie,jb,je,       & ! local domain index boundaries
      ksrc,kent,kprd       ! overflow level indices
   type (block)       :: &
      this_block           ! block information for current block
   real (r8)          :: & 
      ufrc,              & ! fraction of ovf velocity for each box
      Uovf                 ! Uovf at one corner
   logical (log_kind), parameter :: prnt = .false.

   integer (int_kind)    :: ovf_id

   if( prnt .and. my_task == master_task ) then
      write(stdout,*) 'ovf_UV called '
      call shr_sys_flush(stdout)
   endif

!-----------------------------------------------------------------------
!  overflow loop
!-----------------------------------------------------------------------
 !must have the groups initilaized
   if (.not. ovf_groups%init) then
      call ovf_init_groups()
   end if  

   do n=1, ovf_groups%numMyGroups  ! each overflow
      ovf_id = ovf_groups%groupIds(n)
! src
      do m=1,ovf(ovf_id)%num_src  ! source
         ksrc = ovf(ovf_id)%loc_src(m)%k
         ufrc = c1/real(ovf(ovf_id)%num_src-1)
         do iblock = 1,nblocks_clinic
            this_block = get_block(blocks_clinic(iblock),iblock)
            ib = this_block%ib
            ie = this_block%ie
            jb = this_block%jb
            je = this_block%je
            do j=jb,je
               if( ovf(ovf_id)%loc_src(m)%j_u .eq. this_block%j_glob(j) ) then
                  do i=ib,ie
                     if( ovf(ovf_id)%loc_src(m)%i_u .eq. this_block%i_glob(i) ) then
                        if( ovf(ovf_id)%loc_src(m)%orient .eq. 1 ) then
                           Uovf = ovf(ovf_id)%Ms*ufrc/(dz(ksrc)*DYU(i,j,iblock))
                           if( overflows_interactive ) then
                              UVEL(i,j,ksrc,newtime,iblock) = -Uovf
                           endif
                        endif
                        if( ovf(ovf_id)%loc_src(m)%orient .eq. 2 ) then
                           Uovf = ovf(ovf_id)%Ms*ufrc/(dz(ksrc)*DXU(i,j,iblock))
                           if( overflows_interactive ) then
                              VVEL(i,j,ksrc,newtime,iblock) = -Uovf
                           endif
                        endif
                        if( ovf(ovf_id)%loc_src(m)%orient .eq. 3 ) then
                           Uovf = ovf(ovf_id)%Ms*ufrc/(dz(ksrc)*DYU(i,j,iblock))
                           if( overflows_interactive ) then
                              UVEL(i,j,ksrc,newtime,iblock) = +Uovf
                           endif
                        endif
                        if( ovf(ovf_id)%loc_src(m)%orient .eq. 4 ) then
                           Uovf = ovf(ovf_id)%Ms*ufrc/(dz(ksrc)*DXU(i,j,iblock))
                           if( overflows_interactive ) then
                              VVEL(i,j,ksrc,newtime,iblock) = +Uovf
                           endif
                        endif
                        ovf(ovf_id)%loc_src(m)%Uovf = Uovf
                        if( prnt ) then
                           write(stdout,10) n,ovf(ovf_id)%loc_src(m)%i_u, &
                           ovf(ovf_id)%loc_src(m)%j_u,ovf(ovf_id)%loc_src(m)%k,Uovf
                           10 format(' ovf_UV n=',i3,' src i_u j_u k Uovf=', &
                           3(i3,1x),f9.5,2x)
                        endif  ! print
                     endif
                  end do  ! i
               endif
            end do  ! j
         end do  ! iblock
      end do  ! source
! ent
      do m=1,ovf(ovf_id)%num_ent  ! entrainment
         kent = ovf(ovf_id)%loc_ent(m)%k
         ufrc = c1/real(ovf(ovf_id)%num_ent-1)
         do iblock = 1,nblocks_clinic
            this_block = get_block(blocks_clinic(iblock),iblock)
            ib = this_block%ib
            ie = this_block%ie
            jb = this_block%jb
            je = this_block%je
            do j=jb,je
               if( ovf(ovf_id)%loc_ent(m)%j_u .eq. this_block%j_glob(j) ) then
                  do i=ib,ie
                     if( ovf(ovf_id)%loc_ent(m)%i_u .eq. this_block%i_glob(i) ) then
                        if( ovf(ovf_id)%loc_ent(m)%orient .eq. 1 ) then
                           Uovf = ovf(ovf_id)%Me*ufrc/(dz(kent)*DYU(i,j,iblock))
                           if( overflows_interactive ) then
                              UVEL(i,j,kent,newtime,iblock) = -Uovf
                           endif
                        endif
                        if( ovf(ovf_id)%loc_ent(m)%orient .eq. 2 ) then
                           Uovf = ovf(ovf_id)%Me*ufrc/(dz(kent)*DXU(i,j,iblock))
                           if( overflows_interactive ) then
                              VVEL(i,j,kent,newtime,iblock) = -Uovf
                           endif
                        endif
                        if( ovf(ovf_id)%loc_ent(m)%orient .eq. 3 ) then
                           Uovf = ovf(ovf_id)%Me*ufrc/(dz(kent)*DYU(i,j,iblock))
                           if( overflows_interactive ) then
                              UVEL(i,j,kent,newtime,iblock) = +Uovf
                           endif
                        endif
                        if( ovf(ovf_id)%loc_ent(m)%orient .eq. 4 ) then
                           Uovf = ovf(ovf_id)%Me*ufrc/(dz(kent)*DXU(i,j,iblock))
                           if( overflows_interactive ) then
                              VVEL(i,j,kent,newtime,iblock) = +Uovf
                           endif
                        endif
                        ovf(ovf_id)%loc_ent(m)%Uovf = Uovf
                        if( prnt ) then
                           write(stdout,20) n,ovf(ovf_id)%loc_ent(m)%i_u, &
                           ovf(ovf_id)%loc_ent(m)%j_u,ovf(ovf_id)%loc_ent(m)%k,Uovf
                           20 format(' ovf_UV n=',i3,' ent i_u j_u k Uovf=', &
                           3(i3,1x),f9.5,2x)
                        endif  ! print
                     endif
                  end do  ! i
               endif
            end do  ! j
         end do  ! iblock
      end do  ! entrainment
! prd
      m = ovf(ovf_id)%prd_set  ! product set for insertion
      do mp=1,ovf(ovf_id)%num_prd(m)  ! product points for each set
         kprd = ovf(ovf_id)%loc_prd(m,mp)%k
         ufrc = c1/real(ovf(ovf_id)%num_prd(m)-1)
         do iblock = 1,nblocks_clinic
            this_block = get_block(blocks_clinic(iblock),iblock)
            ib = this_block%ib
            ie = this_block%ie
            jb = this_block%jb
            je = this_block%je
            do j=jb,je
               if( ovf(ovf_id)%loc_prd(m,mp)%j_u .eq. this_block%j_glob(j) ) then
                  do i=ib,ie
                     if( ovf(ovf_id)%loc_prd(m,mp)%i_u .eq. this_block%i_glob(i) ) then
                        if( ovf(ovf_id)%loc_prd(m,mp)%orient .eq. 1 ) then
                           Uovf = ovf(ovf_id)%Mp*ufrc/(dz(kprd)*DYU(i,j,iblock))
                           if( overflows_interactive ) then
                              UVEL(i,j,kprd,newtime,iblock) = +Uovf
                           endif
                        endif
                        if( ovf(ovf_id)%loc_prd(m,mp)%orient .eq. 2 ) then
                           Uovf = ovf(ovf_id)%Mp*ufrc/(dz(kprd)*DXU(i,j,iblock))
                           if( overflows_interactive ) then
                              VVEL(i,j,kprd,newtime,iblock) = +Uovf
                           endif
                        endif
                        if( ovf(ovf_id)%loc_prd(m,mp)%orient .eq. 3 ) then
                           Uovf = ovf(ovf_id)%Mp*ufrc/(dz(kprd)*DYU(i,j,iblock))
                           if( overflows_interactive ) then
                              UVEL(i,j,kprd,newtime,iblock) = -Uovf
                           endif
                        endif
                        if( ovf(ovf_id)%loc_prd(m,mp)%orient .eq. 4 ) then
                           Uovf = ovf(ovf_id)%Mp*ufrc/(dz(kprd)*DXU(i,j,iblock))
                           if( overflows_interactive ) then
                              VVEL(i,j,kprd,newtime,iblock) = -Uovf
                           endif
                        endif
                        ovf(ovf_id)%loc_prd(m,mp)%Uovf = Uovf
                        if( prnt ) then
                           write(stdout,30) n,ovf(ovf_id)%loc_prd(m,mp)%i_u, & 
                           ovf(ovf_id)%loc_prd(m,mp)%j_u,ovf(ovf_id)%loc_prd(m,mp)%k,Uovf
                           30 format(' ovf_UV n=',i3,' prd i_u j_u k Uovf=', &
                           3(i3,1x),f9.5,2x)
                        endif  ! print
                     endif
                  end do  ! i
               endif
            end do  ! j
         end do  ! iblock
      end do  ! product points for insertion set
   end do  ! each overflow

!----------------------------------------------------------------------
!EOC

 end subroutine ovf_UV

!***********************************************************************
!EOP
! !IROUTINE: ovf_rhs_brtrpc_momentum
! !INTERFACE:

 subroutine ovf_rhs_brtrpc_momentum(ZX,ZY)

! !DESCRIPTION:
!  Renormalize overflow ZX and ZY vertical integrals of forcing
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!  input variables
!-----------------------------------------------------------------------

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic), &
      intent(inout)      :: &
      ZX, ZY                  ! vertical integrals of forcing

!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind)    :: &
      i,j,k,n,m,mp,         & ! dummy loop indices
      ib,ie,jb,je,          & ! local domain index boundaries
      ksrc,kent,kprd,       & ! level indices
      iblock                  ! block index
   type (block)          :: &
      this_block              ! block information for current block
   real (r8)             :: & 
      dz_sidewall             ! sidewall U-grid depth from top to ovf level
   logical (log_kind), parameter :: prnt = .false.

   integer (int_kind)   :: ovf_id


   if( prnt .and. my_task == master_task ) then
      write(stdout,*) 'ovf_rhs_brtrpc_momentum called '
      call shr_sys_flush(stdout)
   endif

!-----------------------------------------------------------------------
!  overflow loop
!-----------------------------------------------------------------------

 !must have the groups initilaized
   if (.not. ovf_groups%init) then
      call ovf_init_groups()
   end if  

   do n=1, ovf_groups%numMyGroups ! each of my overflows
      ovf_id = ovf_groups%groupIds(n)
! src
      do m=1,ovf(ovf_id)%num_src  ! source
         ksrc = ovf(ovf_id)%loc_src(m)%k
         do iblock = 1,nblocks_clinic
            this_block = get_block(blocks_clinic(iblock),iblock)
            ib = this_block%ib
            ie = this_block%ie
            jb = this_block%jb
            je = this_block%je
            do j=jb,je
               if( ovf(ovf_id)%loc_src(m)%j_u .eq. this_block%j_glob(j) ) then
                  do i=ib,ie
                     if( ovf(ovf_id)%loc_src(m)%i_u .eq. this_block%i_glob(i) ) then
                        dz_sidewall = c0
                        do k=KMU(i,j,iblock)+1,ksrc
                           dz_sidewall = dz_sidewall + dz(k)
                        enddo
                        ZX(i,j,iblock) = (ZX(i,j,iblock)* HU(i,j,iblock)) & 
                           / (HU(i,j,iblock)+dz_sidewall)
                        ZY(i,j,iblock) = (ZY(i,j,iblock)* HU(i,j,iblock)) & 
                           / (HU(i,j,iblock)+dz_sidewall)
                        if(prnt) then
                           write(stdout,10) n,ovf(ovf_id)%loc_src(m)%i_u,ovf(ovf_id)%loc_src(m)%j_u
                           10 format(' ovf_rhs_brtrpc_momentum n=',i3, &
                           ' src ZX,ZY adj at i_u j_u=',2(i3,1x))
                        endif
                     endif
                  end do  ! i
               endif
            end do  ! j
         end do  ! iblock
      end do  ! source
! ent
      do m=1,ovf(ovf_id)%num_ent  ! entrainment
         kent = ovf(ovf_id)%loc_ent(m)%k
         do iblock = 1,nblocks_clinic
            this_block = get_block(blocks_clinic(iblock),iblock)
            ib = this_block%ib
            ie = this_block%ie
            jb = this_block%jb
            je = this_block%je
            do j=jb,je
               if( ovf(ovf_id)%loc_ent(m)%j_u .eq. this_block%j_glob(j) ) then
                  do i=ib,ie
                     if( ovf(ovf_id)%loc_ent(m)%i_u .eq. this_block%i_glob(i) ) then
                        dz_sidewall = c0
                        do k=KMU(i,j,iblock)+1,kent
                           dz_sidewall = dz_sidewall + dz(k)
                        enddo
                        ZX(i,j,iblock) = (ZX(i,j,iblock)* HU(i,j,iblock)) & 
                           / (HU(i,j,iblock)+dz_sidewall)
                        ZY(i,j,iblock) = (ZY(i,j,iblock)* HU(i,j,iblock)) & 
                           / (HU(i,j,iblock)+dz_sidewall)
                        if(prnt) then
                           write(stdout,20) n,ovf(ovf_id)%loc_ent(m)%i_u,ovf(ovf_id)%loc_ent(m)%j_u
                           20 format(' ovf_rhs_brtrpc_momentum n=',i3, &
                           ' ent ZX,ZY adj at i_u j_u=',2(i3,1x))
                        endif
                     endif
                  end do  ! i
               endif
            end do  ! j
         end do  ! iblock
      end do  ! entrainment
! prd
      do m=1,ovf(ovf_id)%num_prd_sets
       do mp=1,ovf(ovf_id)%num_prd(m)  ! product points for each set
         kprd = ovf(ovf_id)%loc_prd(m,mp)%k
         do iblock = 1,nblocks_clinic
            this_block = get_block(blocks_clinic(iblock),iblock)
            ib = this_block%ib
            ie = this_block%ie
            jb = this_block%jb
            je = this_block%je
            do j=jb,je
               if( ovf(ovf_id)%loc_prd(m,mp)%j_u .eq. this_block%j_glob(j) ) then
                  do i=ib,ie
                     if( ovf(ovf_id)%loc_prd(m,mp)%i_u .eq. this_block%i_glob(i) ) then
                        dz_sidewall = c0
                        do k=KMU(i,j,iblock)+1,kprd
                           dz_sidewall = dz_sidewall + dz(k)
                        enddo
                        ZX(i,j,iblock) = (ZX(i,j,iblock)* HU(i,j,iblock)) & 
                           / (HU(i,j,iblock)+dz_sidewall)
                        ZY(i,j,iblock) = (ZY(i,j,iblock)* HU(i,j,iblock)) &
                           / (HU(i,j,iblock)+dz_sidewall)
                        if(prnt) then
                           write(stdout,30) n,ovf(ovf_id)%loc_prd(m,mp)%i_u, &
                                              ovf(ovf_id)%loc_prd(m,mp)%j_u
                           30 format(' ovf_rhs_brtrpc_momentum n=',i3, &
                           ' prd ZX,ZY adj at i_u j_u=',2(i3,1x))
                        endif
                     endif
                  end do  ! i
               endif
            end do  ! j
         end do  ! iblock
       end do  ! product points for each set
      end do  ! all product sets
   end do  ! each overflow

!-----------------------------------------------------------------------
!EOC

 end subroutine ovf_rhs_brtrpc_momentum

!***********************************************************************
!EOP
! !IROUTINE: ovf_brtrpc_renorm
! !INTERFACE:

 subroutine ovf_brtrpc_renorm(WORK3,WORK4,iblock)

! !DESCRIPTION:
!  Renormalize overflow HU for WORK3 and WORK4 in barotropic solution
!  Note- ij limits are 1,nx_block and 1,ny_block to compute ghost values
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!  input variables
!-----------------------------------------------------------------------

   real (r8), dimension(nx_block,ny_block), &
      intent(inout)      :: &
      WORK3,WORK4             ! grid x,y work arrays respectively
   integer (int_kind),      &
      intent(in)         :: &
      iblock                  ! block index   

!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind)    :: &
      i,j,k,n,m,mp,         & ! dummy loop indices
      ksrc,kent,kprd          ! level indices
   type (block)          :: &
      this_block              ! block information for current block
   real (r8)             :: & 
      dz_sidewall             ! sidewall U-grid depth from top to ovf level

   logical (log_kind), parameter :: prnt = .false.

   integer (int_kind) :: ovf_id

   if( prnt .and. my_task == master_task ) then
      write(stdout,*) 'ovf_brtrpc_renorm called '
      call shr_sys_flush(stdout)
   endif

!-----------------------------------------------------------------------
!  overflow loop
!-----------------------------------------------------------------------


 !must have the groups initilaized
   if (.not. ovf_groups%init) then
      call ovf_init_groups()
   end if  

   do n=1, ovf_groups%numMyGroups  ! each overflow
      ovf_id = ovf_groups%groupIds(n)
! src
      do m=1,ovf(ovf_id)%num_src  ! source
         ksrc = ovf(ovf_id)%loc_src(m)%k
         this_block = get_block(blocks_clinic(iblock),iblock)
         do j=1,ny_block
            if( ovf(ovf_id)%loc_src(m)%j_u .eq. this_block%j_glob(j) ) then
               do i=1,nx_block
                  if( ovf(ovf_id)%loc_src(m)%i_u .eq. this_block%i_glob(i) ) then
                     dz_sidewall = c0
                     do k=KMU(i,j,iblock)+1,ksrc
                        dz_sidewall = dz_sidewall + dz(k)
                     enddo
                     WORK3(i,j) = WORK3(i,j)*HUR(i,j,iblock) & 
                                   *(HU(i,j,iblock)+dz_sidewall)
                     WORK4(i,j) = WORK4(i,j)*HUR(i,j,iblock) & 
                                   *(HU(i,j,iblock)+dz_sidewall)
                     if(prnt) then
                        write(stdout,10) n,ovf(ovf_id)%loc_src(m)%i_u, &
                                         ovf(ovf_id)%loc_src(m)%j_u
                        10 format(' ovf_brtrpc_renorm n=',i3, &
                        ' src WORK3/WORK4 adj at i_u j_u=',2(i3,1x))
                     endif
                  endif
               end do  ! i
            endif
         end do  ! j
      end do  ! source
! ent
      do m=1,ovf(ovf_id)%num_ent  ! entrainment
         kent = ovf(ovf_id)%loc_ent(m)%k
         this_block = get_block(blocks_clinic(iblock),iblock)
         do j=1,ny_block
            if( ovf(ovf_id)%loc_ent(m)%j_u .eq. this_block%j_glob(j) ) then
               do i=1,nx_block
                  if( ovf(ovf_id)%loc_ent(m)%i_u .eq. this_block%i_glob(i) ) then
                     dz_sidewall = c0
                     do k=KMU(i,j,iblock)+1,kent
                        dz_sidewall = dz_sidewall + dz(k)
                     enddo
                     WORK3(i,j) = WORK3(i,j)*HUR(i,j,iblock) & 
                                   *(HU(i,j,iblock)+dz_sidewall)
                     WORK4(i,j) = WORK4(i,j)*HUR(i,j,iblock) & 
                                   *(HU(i,j,iblock)+dz_sidewall)
                     if(prnt) then
                        write(stdout,20) n,ovf(ovf_id)%loc_ent(m)%i_u, &
                                         ovf(ovf_id)%loc_ent(m)%j_u
                        20 format(' ovf_brtrpc_renorm n=',i3, &
                        ' ent WORK3/WORK4 adj at i_u j_u=',2(i3,1x))
                     endif
                  endif
               end do  ! i
            endif
         end do  ! j
      end do  ! entrainment
! prd
      do m=1,ovf(ovf_id)%num_prd_sets
       do mp=1,ovf(ovf_id)%num_prd(m)  ! product points for each set
         kprd = ovf(ovf_id)%loc_prd(m,mp)%k
         this_block = get_block(blocks_clinic(iblock),iblock)
         do j=1,ny_block
            if( ovf(ovf_id)%loc_prd(m,mp)%j_u .eq. this_block%j_glob(j) ) then
               do i=1,nx_block
                  if( ovf(ovf_id)%loc_prd(m,mp)%i_u .eq. this_block%i_glob(i) ) then
                     dz_sidewall = c0
                     do k=KMU(i,j,iblock)+1,kprd
                        dz_sidewall = dz_sidewall + dz(k)
                     enddo
                     WORK3(i,j) = WORK3(i,j)*HUR(i,j,iblock) & 
                                   *(HU(i,j,iblock)+dz_sidewall)
                     WORK4(i,j) = WORK4(i,j)*HUR(i,j,iblock) & 
                                   *(HU(i,j,iblock)+dz_sidewall)
                     if(prnt) then
                        write(stdout,30) n,ovf(ovf_id)%loc_prd(m,mp)%i_u, &
                                         ovf(ovf_id)%loc_prd(m,mp)%j_u
                        30 format(' ovf_brtrpc_renorm n=',i3, &
                        ' prd WORK3/WORK4 adj at i_u j_u=',2(i3,1x))
                     endif
                  endif
               end do  ! i
            endif
         end do  ! j
       end do  ! product points for each set
      end do  ! all product sets
   end do  ! each overflow

!-----------------------------------------------------------------------
!EOC

 end subroutine ovf_brtrpc_renorm

!***********************************************************************
!EOP
! !IROUTINE: ovf_rhs_brtrpc_continuity
! !INTERFACE:

 subroutine ovf_rhs_brtrpc_continuity(RHS,iblock)

! !DESCRIPTION:
!  Add overflow vertical velocity to RHS barotropic continuity equation
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!  input variables
!-----------------------------------------------------------------------

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic), &
      intent(inout)     :: &
      RHS                    ! RHS barotropic continuity equation
   integer (int_kind),     &
      intent(in)        :: &
      iblock                 ! block index   

!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind)   :: &
      i,j,n,m,mp,          & ! dummy loop indices
      ib,ie,jb,je            ! local domain index boundaries
   type (block)         :: &
      this_block             ! block information for current block
   logical (log_kind), parameter :: prnt = .false.

   integer (int_kind) :: ovf_id

   if( prnt .and. my_task == master_task ) then
      write(stdout,*) 'ovf_rhs_brtrpc_continuity called '
      call shr_sys_flush(stdout)
   endif

!-----------------------------------------------------------------------
!  overflow loop
!-----------------------------------------------------------------------
 !must have the groups initilaized
   if (.not. ovf_groups%init) then
      call ovf_init_groups()
   end if  
! for each of MY overflows
   do n=1,ovf_groups%numMyGroups
      ovf_id = ovf_groups%groupIds(n)
! src
      do m=1,ovf(ovf_id)%num_src  ! source
         this_block = get_block(blocks_clinic(iblock),iblock)
         ib = this_block%ib
         ie = this_block%ie
         jb = this_block%jb
         je = this_block%je
         do j=jb,je
            if( ovf(ovf_id)%loc_src(m)%j .eq. this_block%j_glob(j) ) then
               do i=ib,ie
                  if( ovf(ovf_id)%loc_src(m)%i .eq. this_block%i_glob(i) ) then
                     RHS(i,j,iblock) = RHS(i,j,iblock) - (ovf(ovf_id)%loc_src(m)%Wovf & 
                                       * TAREA(i,j,iblock)/(beta*c2dtp))
                     if(prnt) then
                        write(stdout,10) n,ovf(ovf_id)%loc_src(m)%i,ovf(ovf_id)%loc_src(m)%j, &
                                         ovf(ovf_id)%loc_src(m)%Wovf,TAREA(i,j,iblock)
                        10 format(' n=',i3,' src RHS adjusted at ij=',2(i3,1x), &
                                  ' Wovf=',f9.6,' TAREA=',1pe11.4)
                     endif
                  endif
               end do  ! i
            endif
         end do  ! j
      end do  ! source
! ent
      do m=1,ovf(ovf_id)%num_ent  ! entrainment
         this_block = get_block(blocks_clinic(iblock),iblock)
         ib = this_block%ib
         ie = this_block%ie
         jb = this_block%jb
         je = this_block%je
         do j=jb,je
            if( ovf(ovf_id)%loc_ent(m)%j .eq. this_block%j_glob(j) ) then
               do i=ib,ie
                  if( ovf(ovf_id)%loc_ent(m)%i .eq. this_block%i_glob(i) ) then
                     RHS(i,j,iblock) = RHS(i,j,iblock) - (ovf(ovf_id)%loc_ent(m)%Wovf & 
                                       * TAREA(i,j,iblock)/(beta*c2dtp))
                     if(prnt) then
                        write(stdout,20) n,ovf(ovf_id)%loc_ent(m)%i,ovf(ovf_id)%loc_ent(m)%j, &
                                         ovf(ovf_id)%loc_ent(m)%Wovf,TAREA(i,j,iblock)
                        20 format(' n=',i3,' ent RHS adjusted at ij=',2(i3,1x), &
                                  ' Wovf=',f9.6,' TAREA=',1pe11.4)
                     endif
                  endif
               end do  ! i
            endif
         end do  ! j
      end do  ! entrainment
! prd
      m = ovf(ovf_id)%prd_set
      do mp=1,ovf(ovf_id)%num_prd(m)  ! product points for each set
         this_block = get_block(blocks_clinic(iblock),iblock)
         ib = this_block%ib
         ie = this_block%ie
         jb = this_block%jb
         je = this_block%je
         do j=jb,je
            if( ovf(ovf_id)%loc_prd(m,mp)%j .eq. this_block%j_glob(j) ) then
               do i=ib,ie
                  if( ovf(ovf_id)%loc_prd(m,mp)%i .eq. this_block%i_glob(i) ) then
                     RHS(i,j,iblock) = RHS(i,j,iblock) - (ovf(ovf_id)%loc_prd(m,mp)%Wovf & 
                                       * TAREA(i,j,iblock)/(beta*c2dtp))
                     if(prnt) then
                        write(stdout,30) n,ovf(ovf_id)%loc_prd(m,mp)%i,ovf(ovf_id)%loc_prd(m,mp)%j, &
                                         ovf(ovf_id)%loc_prd(m,mp)%Wovf,TAREA(i,j,iblock)
                        30 format(' n=',i3,' prd RHS adjusted at ij=',2(i3,1x), &
                                  ' Wovf=',f9.6,' TAREA=',1pe11.4)
                     endif
                  endif
               end do  ! i
            endif
         end do  ! j
      end do  ! product points for each set
   end do  ! each overflow

!-----------------------------------------------------------------------
!EOC

 end subroutine ovf_rhs_brtrpc_continuity

!***********************************************************************
!BOP
! !IROUTINE: ovf_solvers_9pt
! !INTERFACE:

 subroutine ovf_solvers_9pt

! !DESCRIPTION:
!  This routine updates the coefficients of the 9-point stencils for 
!  the barotropic operator for the overflow points
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables:
!
!       {X,Y}{NE,SE,NW,SW} = contribution to {ne,se,nw,sw} coefficients 
!         from {x,y} components of divergence
!       HU = depth at U points
!
!-----------------------------------------------------------------------
   integer (POP_i4) ::  &
      errorCode,        &! error return code
      numBlocksTropic,  &!num local blocks in barotropic distribution
      numBlocksClinic    !num local blocks in baroclinic distribution

   real (POP_r8) ::    &
      xne,xse,xnw,xsw, &! contribution to coefficients from x,y
      yne,yse,ynw,ysw, &!   components of divergence
      ase,anw,asw

   integer (int_kind)   :: &
      i,j,n,               &! dummy counter
      iblock,              &! block counter
      istat

   real (POP_r8), dimension(:,:,:), allocatable :: &
      workNorth,       &!
      workEast,        &!
      workNE,          &!
      HUM               ! HU if no overflows; modified if overflows

!-----------------------------------------------------------------------
!
!  compute nine point operator coefficients: compute on baroclinic
!  decomposition first where grid info defined and redistribute
!  to barotropic distribution
!  leave A0,AC in baroclinic distribution to facilitate easy
!  time-dependent changes in barotropic routine
!
!-----------------------------------------------------------------------

   call POP_DistributionGet(POP_distrbClinic, errorCode, &
                            numLocalBlocks = numBlocksClinic)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ovf_solvers_9pt: error retrieving clinic local block count')
!     activate later, when errorCode is fully supported
!     return
   endif

   allocate(workNorth      (POP_nxBlock,POP_nyBlock,numBlocksClinic), &
            workEast       (POP_nxBlock,POP_nyBlock,numBlocksClinic), &
            workNE         (POP_nxBlock,POP_nyBlock,numBlocksClinic), &
            HUM            (POP_nxBlock,POP_nyBlock,numBlocksClinic), &
            stat=istat)

   if (istat > 0) then
      call POP_ErrorSet(errorCode, &
         'ovf_solvers_9pt: error allocating temporary arrays')
      call exit_POP (sigAbort,'ERROR ovf_solvers_9pt: error allocating temporary arrays')
!     activate later, when errorCode is fully supported
!     return
   endif

   HUM(:,:,:) = HU(:,:,:)
   call ovf_HU(HU,HUM)
   call POP_HaloUpdate(HUM, POP_haloClinic, POP_gridHorzLocNECorner,&
                       POP_fieldKindScalar, errorCode,              &
                       fillValue = 0.0_POP_r8)

   !$OMP PARALLEL DO PRIVATE(iblock,i,j,xne,xse,xnw,xsw,yne,yse,ynw,ysw,ase,anw,asw)
   do iblock = 1,numBlocksClinic

      workNorth            (:,:,iblock) = 0.0_POP_r8
      workEast             (:,:,iblock) = 0.0_POP_r8
      workNE               (:,:,iblock) = 0.0_POP_r8
      centerWgtClinicIndep (:,:,iblock) = 0.0_POP_r8

      do j=2,POP_nyBlock
      do i=2,POP_nxBlock

         xne = 0.25_POP_r8*HUM(i  ,j  ,iblock)*DXUR(i  ,j  ,iblock)* &
                                               DYU (i  ,j  ,iblock)
         xse = 0.25_POP_r8*HUM(i  ,j-1,iblock)*DXUR(i  ,j-1,iblock)* &
                                               DYU (i  ,j-1,iblock)
         xnw = 0.25_POP_r8*HUM(i-1,j  ,iblock)*DXUR(i-1,j  ,iblock)* &
                                               DYU (i-1,j  ,iblock)
         xsw = 0.25_POP_r8*HUM(i-1,j-1,iblock)*DXUR(i-1,j-1,iblock)* &
                                               DYU (i-1,j-1,iblock)

         yne = 0.25_POP_r8*HUM(i  ,j  ,iblock)*DYUR(i  ,j  ,iblock)* &
                                               DXU (i  ,j  ,iblock)
         yse = 0.25_POP_r8*HUM(i  ,j-1,iblock)*DYUR(i  ,j-1,iblock)* &
                                               DXU (i  ,j-1,iblock)
         ynw = 0.25_POP_r8*HUM(i-1,j  ,iblock)*DYUR(i-1,j  ,iblock)* &
                                               DXU (i-1,j  ,iblock)
         ysw = 0.25_POP_r8*HUM(i-1,j-1,iblock)*DYUR(i-1,j-1,iblock)* &
                                               DXU (i-1,j-1,iblock)

         workNE(i,j,iblock) = xne + yne
         ase                = xse + yse
         anw                = xnw + ynw
         asw                = xsw + ysw

         workEast (i,j,iblock)  = xne + xse - yne - yse
         workNorth(i,j,iblock)  = yne + ynw - xne - xnw

         centerWgtClinicIndep(i,j,iblock) = &
                        -(workNE(i,j,iblock) + ase + anw + asw)

      end do
      end do
   end do
 !$OMP END PARALLEL DO


!-----------------------------------------------------------------------
!
!  redistribute operator weights and mask to barotropic distribution
!
!-----------------------------------------------------------------------

   call POP_DistributionGet(POP_distrbTropic, errorCode, &
                            numLocalBlocks = numBlocksTropic)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ovf_solvers_9pt: error retrieving tropic local block count')
      call exit_POP (sigAbort,'ERROR ovf_solvers_9pt: error retrieving tropic local block count')
!     activate later, when errorCode is fully supported
!     return
   endif


   call POP_RedistributeBlocks(btropWgtNorth, POP_distrbTropic, &
                               workNorth,     POP_distrbClinic, errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ovf_solvers_9pt: error redistributing north operator weight')
      call exit_POP (sigAbort,'ERROR ovf_solvers_9pt: error redistributing north operator weight')
!     activate later, when errorCode is fully supported
!     return
   endif


   call POP_RedistributeBlocks(btropWgtEast, POP_distrbTropic, &
                               workEast,     POP_distrbClinic, errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ovf_solvers_9pt: error redistributing east operator weight')
      call exit_POP (sigAbort,'ERROR ovf_solvers_9pt: error redistributing east operator weight')
!     activate later, when errorCode is fully supported
!     return
   endif

   call POP_RedistributeBlocks(btropWgtNE, POP_distrbTropic, &
                               workNE,     POP_distrbClinic, errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ovf_solvers_9pt: error redistributing NE operator weight')
      call exit_POP (sigAbort,'ERROR ovf_solvers_9pt: error redistributing NE operator weight')
!     activate later, when errorCode is fully supported
!     return
   endif

!-----------------------------------------------------------------------
!
!  clean up temporary arrays
!
!-----------------------------------------------------------------------

   deallocate(workNorth, workEast, workNE, HUM, stat=istat)

   if (istat > 0) then
      call POP_ErrorSet(errorCode, &
         'ovf_solvers_9pt: error deallocating temp mask')
      call exit_POP (sigAbort,'ERROR ovf_solvers_9pt: error deallocating temp mask')
!     activate later, when errorCode is fully supported
!     return
   endif


!-----------------------------------------------------------------------
!EOC

 end subroutine ovf_solvers_9pt

!***********************************************************************
!EOP
! !IROUTINE: ovf_HU
! !INTERFACE:

 subroutine ovf_HU(HU,HUM)

! !DESCRIPTION:
!  Modify HU for overflows sidewalls
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!  input variables
!-----------------------------------------------------------------------

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic), &
      intent(in)        :: HU  ! HU 

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic), &
      intent(inout)     :: HUM ! HUM (modified HU)

!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind)   :: &
      i,j,k,n,m,mp,        &   ! dummy loop indices
      ib,ie,jb,je,         &   ! local domain index boundaries
      ksrc,kent,kprd,      &   ! level indices
      iblock                   ! block index
   real (r8)            :: & 
      dz_sidewall              ! sidewall U-grid depth from top to ovf level
   type (block)         :: &
      this_block               ! block information for current block
   logical (log_kind), parameter :: prnt = .false.

   integer (int_kind)  :: ovf_id


   if( prnt .and. my_task == master_task ) then
      write(stdout,*) 'ovf_HU called '
      call shr_sys_flush(stdout)
   endif

!-----------------------------------------------------------------------
!  overflow loop
!-----------------------------------------------------------------------
 !must have the groups initilaized
   if (.not. ovf_groups%init) then
      call ovf_init_groups()
   end if  

   do n=1, ovf_groups%numMyGroups  ! each overflow
      ovf_id = ovf_groups%groupIds(n)
! src
      do m=1,ovf(ovf_id)%num_src  ! source
         ksrc = ovf(ovf_id)%loc_src(m)%k
         do iblock=1,nblocks_clinic
            this_block = get_block(blocks_clinic(iblock),iblock)
            ib = this_block%ib
            ie = this_block%ie
            jb = this_block%jb
            je = this_block%je
            do j=jb,je
               if( ovf(ovf_id)%loc_src(m)%j_u .eq. this_block%j_glob(j) ) then
                  do i=ib,ie
                     if( ovf(ovf_id)%loc_src(m)%i_u .eq. this_block%i_glob(i) ) then
                        dz_sidewall = c0
                        do k=KMU(i,j,iblock)+1,ksrc
                           dz_sidewall = dz_sidewall + dz(k)
                        enddo
                        HUM(i,j,iblock) = HU(i,j,iblock) + dz_sidewall
                        if(prnt) then
                           write(stdout,10) n, &
                           ovf(ovf_id)%loc_src(m)%i_u,ovf(ovf_id)%loc_src(m)%j_u
                           10 format(' n=',i3,' src HU adjusted at i_u j_u =',2(i3,1x))
                        endif
                     endif
                  end do  ! i
               endif
            end do  ! j
         end do  ! iblocks
      end do  ! source
! ent
      do m=1,ovf(ovf_id)%num_ent  ! entrainment
         kent = ovf(ovf_id)%loc_ent(m)%k
         do iblock=1,nblocks_clinic
            this_block = get_block(blocks_clinic(iblock),iblock)
            ib = this_block%ib
            ie = this_block%ie
            jb = this_block%jb
            je = this_block%je
            do j=jb,je
               if( ovf(ovf_id)%loc_ent(m)%j_u .eq. this_block%j_glob(j) ) then
                  do i=ib,ie
                     if( ovf(ovf_id)%loc_ent(m)%i_u .eq. this_block%i_glob(i) ) then
                        dz_sidewall = c0
                        do k=KMU(i,j,iblock)+1,kent
                           dz_sidewall = dz_sidewall + dz(k)
                        enddo
                        HUM(i,j,iblock) = HU(i,j,iblock) + dz_sidewall
                        if(prnt) then
                           write(stdout,20) n, &
                           ovf(ovf_id)%loc_ent(m)%i_u,ovf(ovf_id)%loc_ent(m)%j_u
                           20 format(' n=',i3,' ent HU adjusted at i_u j_u =',2(i3,1x))
                        endif
                     endif
                  end do  ! i
               endif
            end do  ! j
         end do  ! iblocks
      end do  ! entrainment
! prd
      do m=1,ovf(ovf_id)%num_prd_sets
       do mp=1,ovf(ovf_id)%num_prd(m)  ! product points for each set
         kprd = ovf(ovf_id)%loc_prd(m,mp)%k
         do iblock = 1,nblocks_clinic
            this_block = get_block(blocks_clinic(iblock),iblock)
            ib = this_block%ib
            ie = this_block%ie
            jb = this_block%jb
            je = this_block%je
            do j=jb,je
               if( ovf(ovf_id)%loc_prd(m,mp)%j_u .eq. this_block%j_glob(j) ) then
                 do i=ib,ie
                     if( ovf(ovf_id)%loc_prd(m,mp)%i_u .eq. this_block%i_glob(i) ) then
                        dz_sidewall = c0
                        do k=KMU(i,j,iblock)+1,kprd
                           dz_sidewall = dz_sidewall + dz(k)
                        enddo
                        HUM(i,j,iblock) = HU(i,j,iblock) + dz_sidewall
                        if(prnt) then
                           write(stdout,30) n,ovf(ovf_id)%loc_prd(m,mp)%i_u, &
                                            ovf(ovf_id)%loc_prd(m,mp)%j_u
                           30 format(' n=',i3,' prd HU adjusted at i_u j_u =',2(i3,1x))
                        endif
                     endif
                  end do  ! i
               endif
            end do  ! j
         end do  ! iblock
       end do  ! product points for each set
      end do  ! all product sets
   end do  ! each overflow

!-----------------------------------------------------------------------
!EOC

 end subroutine ovf_HU

!***********************************************************************
!EOP
! !IROUTINE: ovf_UV_solution
! !INTERFACE:

 subroutine ovf_UV_solution

! !DESCRIPTION:
!  Evaluate ovf column solution for baroclinic U and V. Should be called
!  BEFORE the final addition of baroclinic and barotropic velocities.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind)    :: &
      iblock,i,j,k,n,m,mp,  & ! dummy loop indices
      ib,ie,jb,je,          & ! local domain index boundaries
      ksrc,kent,kprd          ! overflow level indices

   real (r8)             :: &
      Uovf,                 & ! overflow U
      Uovf_nm1,             & ! overflow U at n-1
      ubar,                 & ! barotropic velocity
      utlda(km)               ! unnormalized baroclinic velocity

   type (block)          :: &
      this_block              ! block information for current block
   logical (log_kind), parameter :: prnt = .false.
   
   integer (int_kind)     :: ovf_id

   if( prnt .and. my_task == master_task ) then
      write(stdout,*) 'ovf_UV_solution called '
      call shr_sys_flush(stdout)
   endif

!-----------------------------------------------------------------------
!  overflow loop
!-----------------------------------------------------------------------
 !must have the groups initilaized
   if (.not. ovf_groups%init) then
      call ovf_init_groups()
   end if  

   do n=1, ovf_groups%numMyGroups  ! each overflow
      ovf_id = ovf_groups%groupIds(n)
! src
      do m=1,ovf(ovf_id)%num_src  ! source
         ksrc = ovf(ovf_id)%loc_src(m)%k
         do iblock = 1,nblocks_clinic
            this_block = get_block(blocks_clinic(iblock),iblock)
            ib = this_block%ib
            ie = this_block%ie
            jb = this_block%jb
            je = this_block%je
            do j=jb,je
               if( ovf(ovf_id)%loc_src(m)%j_u .eq. this_block%j_glob(j) ) then
                  do i=ib,ie
                     if( ovf(ovf_id)%loc_src(m)%i_u .eq. this_block%i_glob(i) ) then
                        if(prnt) then
                           write(stdout,10) n,ovf(ovf_id)%loc_src(m)%i_u, &
                                            ovf(ovf_id)%loc_src(m)%j_u 
                           10 format(' n=',i3,' src iu ju column evaluated=',2(i3,1x))
                        endif
! U
                        do k=1,km
                          utlda(k)  = ovf(ovf_id)%loc_src(m)%Utlda(k)
                        enddo
                        Uovf     = UVEL(i,j,ksrc,newtime,iblock)
                        ubar     = UBTROP(i,j,newtime,iblock)

                        Uovf_nm1 = UVEL(i,j,ksrc,oldtime,iblock)
                        call ovf_U_column(i,j,ksrc,iblock, &
                                 Uovf,ubar,utlda,Uovf_nm1)
! V 
                        do k=1,km
                          utlda(k)  = ovf(ovf_id)%loc_src(m)%Vtlda(k)
                        enddo
                        Uovf     = VVEL(i,j,ksrc,newtime,iblock)
                        ubar     = VBTROP(i,j,newtime,iblock)
                        Uovf_nm1 = VVEL(i,j,ksrc,oldtime,iblock)
                        call ovf_V_column(i,j,ksrc,iblock, &
                                 Uovf,ubar,utlda,Uovf_nm1)
                     endif
                  end do  ! i
               endif
            end do  ! j
         end do  ! iblock
      end do  ! source
! ent
      do m=1,ovf(ovf_id)%num_ent  ! entrainment
         kent = ovf(ovf_id)%loc_ent(m)%k
         do iblock = 1,nblocks_clinic
            this_block = get_block(blocks_clinic(iblock),iblock)
            ib = this_block%ib
            ie = this_block%ie
            jb = this_block%jb
            je = this_block%je
            do j=jb,je
               if( ovf(ovf_id)%loc_ent(m)%j_u .eq. this_block%j_glob(j) ) then
                  do i=ib,ie
                     if( ovf(ovf_id)%loc_ent(m)%i_u .eq. this_block%i_glob(i) ) then
                        if(prnt) then
                           write(stdout,20) n,ovf(ovf_id)%loc_ent(m)%i_u, &
                                            ovf(ovf_id)%loc_ent(m)%j_u 
                           20 format(' n=',i3,' ent iu ju column evaluated=',2(i3,1x))
                        endif
! U
                        do k=1,km
                          utlda(k)  = ovf(ovf_id)%loc_ent(m)%Utlda(k)
                        enddo
                        Uovf     = UVEL(i,j,kent,newtime,iblock)
                        ubar     = UBTROP(i,j,newtime,iblock)
                        Uovf_nm1 = UVEL(i,j,kent,oldtime,iblock)
                        call ovf_U_column(i,j,kent,iblock, &
                                 Uovf,ubar,utlda,Uovf_nm1)
! V 
                        do k=1,km
                          utlda(k)  = ovf(ovf_id)%loc_ent(m)%Vtlda(k)
                        enddo
                        Uovf     = VVEL(i,j,kent,newtime,iblock)
                        ubar     = VBTROP(i,j,newtime,iblock)
                        Uovf_nm1 = VVEL(i,j,kent,oldtime,iblock)
                        call ovf_V_column(i,j,kent,iblock, &
                                 Uovf,ubar,utlda,Uovf_nm1)
                     endif
                  end do  ! i
               endif
            end do  ! j
         end do  ! iblock
      end do  ! entrainment
! prd
      do m=1,ovf(ovf_id)%num_prd_sets
       do mp=1,ovf(ovf_id)%num_prd(m)  ! product points for each set
         kprd = ovf(ovf_id)%loc_prd(m,mp)%k
         do iblock = 1,nblocks_clinic
            this_block = get_block(blocks_clinic(iblock),iblock)
            ib = this_block%ib
            ie = this_block%ie
            jb = this_block%jb
            je = this_block%je
            do j=jb,je
               if( ovf(ovf_id)%loc_prd(m,mp)%j_u .eq. this_block%j_glob(j) ) then
                  do i=ib,ie
                     if( ovf(ovf_id)%loc_prd(m,mp)%i_u .eq. this_block%i_glob(i) ) then
                        if(prnt) then
                           write(stdout,30) n,ovf(ovf_id)%loc_prd(m,mp)%i_u, &
                                            ovf(ovf_id)%loc_prd(m,mp)%j_u 
                           30 format(' n=',i3,' prd iu ju column evaluated=',2(i3,1x))
                        endif
! U
                        do k=1,km
                          utlda(k)  = ovf(ovf_id)%loc_prd(m,mp)%Utlda(k)
                        enddo
                        Uovf     = UVEL(i,j,kprd,newtime,iblock)
                        ubar     = UBTROP(i,j,newtime,iblock)
                        Uovf_nm1 = UVEL(i,j,kprd,oldtime,iblock)
                        call ovf_U_column(i,j,kprd,iblock, &
                                 Uovf,ubar,utlda,Uovf_nm1)
! V 
                        do k=1,km
                          utlda(k)  = ovf(ovf_id)%loc_prd(m,mp)%Vtlda(k)
                        enddo
                        Uovf     = VVEL(i,j,kprd,newtime,iblock)
                        ubar     = VBTROP(i,j,newtime,iblock)
                        Uovf_nm1 = VVEL(i,j,kprd,oldtime,iblock)
                        call ovf_V_column(i,j,kprd,iblock, &
                                 Uovf,ubar,utlda,Uovf_nm1)
                     endif
                  end do  ! i
               endif
            end do  ! j
         end do  ! iblock
       end do  ! product points for each set
      end do  ! all product sets
   end do  ! each overflow

!----------------------------------------------------------------------
!EOC

 end subroutine ovf_UV_solution

!***********************************************************************
!EOP
! !IROUTINE: ovf_U_column
! !INTERFACE:

 subroutine ovf_U_column(i,j,kovf,iblock,Uovf,ubar,utlda,Uovf_nm1)

! !DESCRIPTION:
!  Evaluate ovf column solution for U baroclinic
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!  input variables
!-----------------------------------------------------------------------

   integer (int_kind),       &
      intent(in)          :: &
      i,                     & ! local block i index on u-grid
      j,                     & ! local block j index on u-grid
      kovf,                  & ! k index of overflow
      iblock                   ! block index

   real (r8), intent(in)  :: &
      Uovf,                  & ! overflow U
      ubar,                  & ! barotropic velocity
      utlda(km),             & ! unnormalized baroclinic velocity
      Uovf_nm1                 ! overflow U at n-1

!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind)     :: &
      k                        ! vertical loop index

   real (r8)              :: & 
      uprime,                & ! overflow sidewall baroclinic velocity
      hu,                    & ! HU after accumulation of column dz
      vert_sum,              & ! vertical sum accumulation of utlda*dz
      utlda_bar                ! vertical mean of utlda

   logical (log_kind), parameter :: prnt = .false.

   if( prnt .and. my_task == master_task ) then
     write(stdout,*) 'ovf_U_column called '
     call shr_sys_flush(stdout)
   endif

!  evaluate baroclinic normalization for the overflow column by including 
!  the overflow contributions along the sidewall above the overflow

!  above the topography
   hu       = c0
   vert_sum = c0
   do k=1,KMU(i,j,iblock)
     hu       = hu       +          dz(k)
     vert_sum = vert_sum + utlda(k)*dz(k)
   enddo

!  below the topography but above the overflow
   uprime   = -ubar
   do k=KMU(i,j,iblock)+1,kovf-1
     vert_sum = vert_sum + uprime*dz(k)
   enddo

!  the overflow contribution
   uprime   = Uovf - ubar
   vert_sum = vert_sum + uprime*dz(kovf)

!  adjusted utlda_bar
   utlda_bar = vert_sum/hu

!  evaluate overflow modified baroclinic velocity for the column
   do k=1,KMU(i,j,iblock)
     UVEL(i,j,k,newtime,iblock) = utlda(k) - utlda_bar
   enddo

!  check of zero vertical sum 
   hu       = c0
   vert_sum = c0
   do k=1,KMU(i,j,iblock)
     hu       = hu       +                            dz(k)
     vert_sum = vert_sum + UVEL(i,j,k,newtime,iblock)*dz(k)
   end do
   uprime   = -ubar
   do k=KMU(i,j,iblock)+1,kovf-1
     vert_sum = vert_sum + uprime*dz(k)
   enddo
   uprime   = Uovf - ubar
   vert_sum = vert_sum + uprime*dz(kovf)
   vert_sum = vert_sum / hu

   if( prnt ) then
      write(stdout,*) 'ovf_U_column '
      write(stdout,5) KMU(i,j,iblock),kovf,Uovf,ubar,utlda_bar,Uovf_nm1, &
                      vert_sum
      5 format(' kmu,kovf,Uovf ubar utlda_bar Uovf_nm1 vert_sum='/ &
                 1x,2(i3,1x),2x,4(f10.5,1x),1pe11.4)
      do k=1,kovf
        if( k <= KMU(i,j,iblock) ) then
          write(stdout,10) k,dz(k),utlda(k)-utlda_bar
          10 format('   k dz U_baroclinic =',i3,1x, &
                      f8.2,1x,f9.5)
        else if( k > KMU(i,j,iblock) .and. k < kovf ) then
          write(stdout,15) k,dz(k),-ubar
          15 format('   k dz U_baroclinic =',i3,1x, &
                      f8.2,1x,f9.5)
        else
          write(stdout,20) k,dz(k),Uovf-ubar
          20 format('   k dz U_baroclinic =',i3,1x, &
                      f8.2,1x,f9.5)
        endif
      end do
   endif

!----------------------------------------------------------------------
!EOC

 end subroutine ovf_U_column

!***********************************************************************
!EOP
! !IROUTINE: ovf_V_column
! !INTERFACE:

 subroutine ovf_V_column(i,j,kovf,iblock,Uovf,ubar,utlda,Uovf_nm1)

! !DESCRIPTION:
!  Evaluate ovf column solution for V baroclinic
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!  input variables
!-----------------------------------------------------------------------

   integer (int_kind),       &
      intent(in)          :: &
      i,                     & ! local block i index on u-grid
      j,                     & ! local block j index on u-grid
      kovf,                  & ! k index of overflow
      iblock                   ! block index

   real (r8), intent(in)  :: &
      Uovf,                  & ! overflow U
      ubar,                  & ! barotropic velocity
      utlda(km),             & ! unnormalized baroclinic velocity
      Uovf_nm1                 ! overflow U at n-1

!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind)     :: &
      k                        ! vertical loop index

   real (r8)              :: & 
      uprime,                & ! overflow sidewall baroclinic velocity
      hu,                    & ! HU after accumulation of column dz
      vert_sum,              & ! vertical sum accumulation of utlda*dz
      utlda_bar                ! vertical mean of utlda

   logical (log_kind), parameter :: prnt = .false.

   if( prnt .and. my_task == master_task ) then
     write(stdout,*) 'ovf_V_column called '
     call shr_sys_flush(stdout)
   endif

!  evaluate baroclinic normalization for the overflow column by including
!  the overflow contributions along the sidewall above the overflow

!  above the topography
   hu       = c0
   vert_sum = c0
   do k=1,KMU(i,j,iblock)
     hu       = hu       +          dz(k)
     vert_sum = vert_sum + utlda(k)*dz(k)
   enddo

!  below the topography but above the overflow
   uprime   = -ubar
   do k=KMU(i,j,iblock)+1,kovf-1
     vert_sum = vert_sum + uprime*dz(k)
   enddo

!  the overflow contribution
   uprime   = Uovf - ubar
   vert_sum = vert_sum + uprime*dz(kovf)

!  adjusted utlda_bar
   utlda_bar = vert_sum/hu

!  evaluate overflow modified baroclinic velocity for the column
   do k=1,KMU(i,j,iblock)
     VVEL(i,j,k,newtime,iblock) = utlda(k) - utlda_bar
   enddo

!  check of zero vertical sum 
   hu       = c0
   vert_sum = c0
   do k=1,KMU(i,j,iblock)
     hu       = hu       +                            dz(k)
     vert_sum = vert_sum + VVEL(i,j,k,newtime,iblock)*dz(k)
   end do
   uprime   = -ubar
   do k=KMU(i,j,iblock)+1,kovf-1
     vert_sum = vert_sum + uprime*dz(k)
   enddo
   uprime   = Uovf - ubar
   vert_sum = vert_sum + uprime*dz(kovf)
   vert_sum = vert_sum / hu

   if( prnt ) then
      write(stdout,*) 'ovf_V_column '
      write(stdout,5) KMU(i,j,iblock),kovf,Uovf,ubar,utlda_bar,Uovf_nm1, &
                      vert_sum
      5 format(' kmu,kovf,Uovf ubar utlda_bar Uovf_nm1 vert_sum='/ &
                 1x,2(i3,1x),2x,4(f10.5,1x),1pe11.4)
      do k=1,kovf
        if( k <= KMU(i,j,iblock) ) then
          write(stdout,10) k,dz(k),utlda(k)-utlda_bar
          10 format('   k dz U_baroclinic =',i3,1x, &
                      f8.2,1x,f9.5)
        else if( k > KMU(i,j,iblock) .and. k < kovf ) then
          write(stdout,15) k,dz(k),-ubar
          15 format('   k dz U_baroclinic =',i3,1x, &
                      f8.2,1x,f9.5)
        else
          write(stdout,20) k,dz(k),Uovf-ubar
          20 format('   k dz U_baroclinic =',i3,1x, &
                      f8.2,1x,f9.5)
        endif
      end do
   endif

!----------------------------------------------------------------------
!EOC

 end subroutine ovf_V_column

!***********************************************************************

 end module overflows

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
