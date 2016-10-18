module subcol_tstcp
   !---------------------------------------------------------------------------
   ! Purpose:
   !
   ! Implement the various TestCopy schemes
   !  sub-column schemes
   !
   !---------------------------------------------------------------------------

   use shr_kind_mod,    only: r8=>shr_kind_r8
   use physics_types,   only: physics_state, physics_tend, physics_ptend
   use ppgrid,          only: pcols, psubcols, pver, pverp
   use constituents,    only: pcnst
   use cam_abortutils,  only: endrun
   use spmd_utils,      only: masterproc
   use cam_logfile,     only: iulog

   implicit none

   private
   save

   public :: subcol_gen_tstcp 
   public :: subcol_register_tstcp
   public :: subcol_readnl_tstcp
   public :: subcol_field_avg_tstcp
   public :: subcol_ptend_avg_tstcp

   interface subcol_field_avg_tstcp
      module procedure subcol_field_avg_tstcp_1dr 
      module procedure subcol_field_avg_tstcp_1di 
      module procedure subcol_field_avg_tstcp_2dr 
   end interface

   logical :: subcol_tstcp_noAvg   ! if set, bypasses averaging and assigns back the first subcolumn to grid

   logical :: subcol_tstcp_filter  ! if set, sets up a filter which yields BFB results 
                                   ! (doesn't really excercise the filter arithmetic)

   logical :: subcol_tstcp_weight  ! if set, sets up a weight which yields BFB results 
                                   ! (doesn't really excercise the weight arithmetic)

   logical :: subcol_tstcp_perturb ! if set, turns on the perturbation test which changes the state temperatures
                                   ! to make sure subcolumns differ

   logical :: subcol_tstcp_restart ! if set, sets up weights so that they are more adequately tested in restart,
                                   ! but will not be BFB with non-subcolumnized run

   integer :: tstcpy_scol_idx      ! pbuf index for subcolumn-only test field

contains

   subroutine subcol_register_tstcp()
      use physics_buffer,  only: pbuf_add_field, dtype_i4, col_type_subcol
      use phys_control,    only: phys_getopts

      ! A subcolumn-only test field
      ! pbuf is global so it will show up in restart file
      call pbuf_add_field('TSTCPY_SCOL','global', dtype_i4,                 &
           (/pcols,pver/), tstcpy_scol_idx, col_type_subcol)

   end subroutine subcol_register_tstcp

   subroutine subcol_readnl_tstcp(nlfile)
      use namelist_utils,  only: find_group_name
      use units,           only: getunit, freeunit
      use spmd_utils,      only: masterproc, mpi_logical, masterprocid, mpicom

      character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

      ! Local variables
      integer :: unitn, ierr

      namelist /subcol_tstcp_nl/ subcol_tstcp_noAvg, subcol_tstcp_filter, subcol_tstcp_weight, subcol_tstcp_perturb, &
                                 subcol_tstcp_restart

      !-----------------------------------------------------------------------------

      if (masterproc) then
         unitn = getunit()
         open( unitn, file=trim(nlfile), status='old' )
         call find_group_name(unitn, 'subcol_tstcp_nl', status=ierr)
         if (ierr == 0) then
            read(unitn, subcol_tstcp_nl, iostat=ierr)
            if (ierr /= 0) then
               call endrun('subcol_readnl_tstcp: ERROR reading namelist')
            end if
         end if
         close(unitn)
         call freeunit(unitn)
      end if

#ifdef SPMD
      ! Broadcast namelist variables
      call mpi_bcast(subcol_tstcp_noAvg,   1, mpi_logical, masterprocid, mpicom, ierr)
      call mpi_bcast(subcol_tstcp_filter,  1, mpi_logical, masterprocid, mpicom, ierr)
      call mpi_bcast(subcol_tstcp_weight,  1, mpi_logical, masterprocid, mpicom, ierr)
      call mpi_bcast(subcol_tstcp_perturb, 1, mpi_logical, masterprocid, mpicom, ierr)
      call mpi_bcast(subcol_tstcp_restart, 1, mpi_logical, masterprocid, mpicom, ierr)
#endif
   end subroutine subcol_readnl_tstcp

   subroutine subcol_gen_tstcp(state, tend, state_sc, tend_sc, pbuf)

      use subcol_utils,   only: subcol_set_subcols, subcol_get_nsubcol, subcol_set_weight, subcol_set_filter
      use physics_buffer, only: physics_buffer_desc, pbuf_get_field, col_type_subcol
      use phys_grid,      only: get_gcol_p
      use time_manager,   only: is_first_step, is_first_restart_step


      !-----------------------------------
      ! sub-column generator
      !-----------------------------------
      type(physics_state), intent(inout) :: state
      type(physics_tend),  intent(inout) :: tend
      type(physics_state), intent(inout) :: state_sc        ! sub-column state
      type(physics_tend),  intent(inout) :: tend_sc         ! sub-column tend
      type(physics_buffer_desc), pointer :: pbuf(:)


      !
      ! Local variables
      !
      integer            :: i, k, ngrdcol, indx, indx1, indx2
      integer            :: nsubcol(pcols)
      real(r8)           :: weight(state_sc%psetcols)
      integer            :: filter(state_sc%psetcols)
      integer, pointer   :: test_field(:,:)
      character(len=128) :: errmsg

      ngrdcol    = state%ngrdcol

      !----------------------
      ! Set the number of subcolumns on the 0th time step -- current implementation does not allow
      ! number of subcolumns to vary within a run.  Cannot be done in init as ngrdcol is not known
      ! at init
      !----------------------
      ! Set two subcolumns to 1 to test random number of subcols
      ! Index logic is arithmetic as some compilers balk when hardcoded larger than pcols
      if (is_first_step()) then
         nsubcol(:) = psubcols  ! For test, set all to max number of sub-columns, then reset a few values
         indx1      = anint(pcols/2._r8)
         indx2      = anint(pcols/5._r8)
         if (pcols >= 6) nsubcol(indx1) = 1
         if (pcols >= 11) nsubcol(indx2) = 1
         if (ngrdcol+1 <= pcols) nsubcol(ngrdcol+1:) = 0

         ! Set up the weights once and do not modify - this will test the restart ability to correctly retrieve them
         if(subcol_tstcp_restart) then
           weight=0._r8
           indx=1
           do i=1,ngrdcol
              weight(indx:indx+nsubcol(i)-1)=1._r8/nsubcol(i)
              if (state%lon(i) < -0.5236_r8) then
                if (nsubcol(i) >= 3) then
                  weight(indx) = 2*1._r8/nsubcol(i)
                  weight(indx+1) = 1._r8 - weight(indx)
                  weight(indx+2:indx+nsubcol(i)-1)=0._r8
                end if
              end if
              indx = indx+nsubcol(i)
           end do
           call subcol_set_weight(state%lchnk, weight)
         end if
      else
         call subcol_get_nsubcol(state%lchnk, nsubcol)
         ! Since this is a test generator, check for nsubcol correctness.
         indx1      = anint(pcols/2._r8)
         indx2      = anint(pcols/5._r8)
10       format(a,i3,a,i5)
         do i = 1, pcols
            if ((pcols >= 11) .and. (i == indx2)) then
               if (nsubcol(i) /= 1) then
                  write(errmsg, 10) 'subcol_gen_tstcp: Bad value for nsubcol(',&
                       i,') = ',nsubcol(i)
                  call endrun(errmsg)
               end if
            else if ((pcols >= 6) .and. (i == indx1)) then
               if (nsubcol(i) /= 1) then
                  write(errmsg, 10) 'subcol_gen_tstcp: Bad value for nsubcol(',&
                       i,') = ',nsubcol(i)
                  call endrun(errmsg)
               end if
            else if (i > ngrdcol) then
               if (nsubcol(i) /= 0) then
                  write(errmsg, 10) 'subcol_gen_tstcp: Bad value for nsubcol(',&
                       i,') = ',nsubcol(i)
                  call endrun(errmsg)
               end if
            else
               if (nsubcol(i) /= psubcols) then
                  write(errmsg, 10) 'subcol_gen_tstcp: Bad value for nsubcol(',&
                       i,') = ',nsubcol(i)
                  call endrun(errmsg)
               end if
            end if
         end do
      end if

      call subcol_set_subcols(state, tend, nsubcol, state_sc, tend_sc)

      ! For perturb case, adjust Temperature up and down one degree
      if (subcol_tstcp_perturb) then
        indx=1
        do i=1,ngrdcol
          if (nsubcol(i) >= 2) then
            state_sc%t(indx,:) = state_sc%t(indx,:)+1
            state_sc%t(indx+1,:) = state_sc%t(indx+1,:)-1
          end if
          indx=indx+nsubcol(i)
        end do
      end if
          
      ! Set weight to 1 for first column, 0 for all others -- will be BFB with noUniAv case
      if(subcol_tstcp_filter .and. subcol_tstcp_weight) then
        weight=1._r8
        ! Initialize to 1 - will match doAv_noUni, init to 0 - will match noUniAv
        filter=1
        indx=1
        do i=1,ngrdcol
           weight(indx) = 1.0_r8
           filter(indx) = 1
           indx = indx+nsubcol(i)
        end do
        call subcol_set_weight(state%lchnk, weight)
        call subcol_set_filter(state%lchnk, filter)
      ! Set weight to 1 for first column, 0 for all others -- will be BFB with noUniAv case
      else if(subcol_tstcp_weight) then
        weight=0._r8
        indx=1
        do i=1,ngrdcol
           weight(indx) = 1.0_r8
           indx = indx+nsubcol(i)
        end do
        call subcol_set_weight(state%lchnk, weight)

      ! Set filter to 1 for first column, 0 for all others -- will be BFB with noUniAv case
      else if(subcol_tstcp_filter) then
        filter=0
        indx=1
        do i=1,ngrdcol
           filter(indx) = 1
           indx = indx+nsubcol(i)
        end do
        call subcol_set_filter(state%lchnk, filter)
      end if


      if (is_first_restart_step()) then
         ! Test values for the test pbuf
         call pbuf_get_field(pbuf, tstcpy_scol_idx, test_field,               &
              col_type=col_type_subcol, copy_if_needed=.false.)
         indx = 1
         do i=1,ngrdcol
            do indx1 = 1, nsubcol(i)
               do k = 1, pver
                  indx2 = (get_gcol_p(state%lchnk, i) * 10000)
                  indx2 = k + (100 * (indx1 + indx2))
                  if(test_field(indx, k) /= indx2) then
                     write(iulog, *) 'TSTCPY_SCOL check(',indx,',',k,         &
                          '): expected',indx2,', found',test_field(indx, k)
                     call endrun("Restart check for TSTCPY_SCOL failed")
                  end if
               end do
               indx = indx + 1
            end do
         end do
         ! Unused subcolumn space is not initialized so no check
      else if (is_first_step()) then
         ! Set values for the test pbuf
         call pbuf_get_field(pbuf, tstcpy_scol_idx, test_field,               &
              col_type=col_type_subcol, copy_if_needed=.false.)
         test_field = -1
         indx = 1
         do i=1,ngrdcol
            do indx1 = 1, nsubcol(i)
               do k = 1, pver
                  indx2 = (get_gcol_p(state%lchnk, i) * 10000)
                  indx2 = k + (100 * (indx1 + indx2))
                  test_field(indx, k) = indx2
               end do
               indx = indx + 1
            end do
         end do
      end if

end subroutine subcol_gen_tstcp

subroutine subcol_field_avg_tstcp_1dr (field_sc, ngrdcol, lchnk, field)
      use physics_buffer,   only: physics_buffer_desc
      use subcol_utils,     only: subcol_field_get_firstsubcol, subcol_field_avg_shr, is_filter_set, is_weight_set

      !-----------------------------------
      ! Average the subcolumns dimension (pcols*psubcols) to the grid dimension (pcols)
      !-----------------------------------

      real(r8), intent(in)                        :: field_sc(:)   ! intent in
      integer,  intent(in)                        :: ngrdcol       ! # grid cols
      integer,  intent(in)                        :: lchnk         ! chunk index
      real(r8), intent(out)                       :: field(:)

      ! Unless specialized averaging is needed, most subcolumn schemes will be handled here
      if (subcol_tstcp_noAvg) then
         call subcol_field_get_firstsubcol(field_sc, .true., ngrdcol, lchnk, field)
      else 
         call subcol_field_avg_shr(field_sc, ngrdcol, lchnk, field, is_filter_set(), is_weight_set())
      end if

end subroutine subcol_field_avg_tstcp_1dr

subroutine subcol_field_avg_tstcp_1di (field_sc, ngrdcol, lchnk, field)
      use physics_buffer,   only: physics_buffer_desc
      use subcol_utils,     only: subcol_field_get_firstsubcol, subcol_field_avg_shr, is_filter_set, is_weight_set

      !-----------------------------------
      ! Average the subcolumns dimension (pcols*psubcols) to the grid dimension (pcols)
      !-----------------------------------

      integer, intent(in)                         :: field_sc(:)   ! intent in
      integer, intent(in)                         :: ngrdcol       ! # grid cols
      integer, intent(in)                         :: lchnk         ! chunk index
      integer, intent(out)                        :: field(:)

      ! Unless specialized averaging is needed, most subcolumn schemes will be handled here
      if (subcol_tstcp_noAvg) then
         call subcol_field_get_firstsubcol(field_sc, .true., ngrdcol, lchnk, field)
      else
         call subcol_field_avg_shr(field_sc, ngrdcol, lchnk, field, is_filter_set(), is_weight_set())
      end if

end subroutine subcol_field_avg_tstcp_1di

subroutine subcol_field_avg_tstcp_2dr (field_sc, ngrdcol, lchnk, field)
      use physics_buffer,   only: physics_buffer_desc
      use subcol_utils,     only: subcol_field_get_firstsubcol, subcol_field_avg_shr, is_filter_set, is_weight_set

      !-----------------------------------
      ! Average the subcolumns dimension (pcols*psubcols) to the grid dimension (pcols)
      !-----------------------------------

      real(r8), intent(in)                        :: field_sc(:,:)   ! intent in
      integer,  intent(in)                        :: ngrdcol       ! # grid cols
      integer,  intent(in)                        :: lchnk         ! chunk index
      real(r8), intent(out)                       :: field(:,:)

      ! Unless specialized averaging is needed, most subcolumn schemes will be handled here
      if (subcol_tstcp_noAvg) then
         call subcol_field_get_firstsubcol(field_sc, .true., ngrdcol, lchnk, field)
      else
         call subcol_field_avg_shr(field_sc, ngrdcol, lchnk, field, is_filter_set(), is_weight_set())
      end if

end subroutine subcol_field_avg_tstcp_2dr

subroutine subcol_ptend_avg_tstcp (ptend_sc, ngrdcol, lchnk, ptend)
      use physics_buffer,   only: physics_buffer_desc
      use subcol_utils,     only: subcol_ptend_get_firstsubcol, subcol_ptend_avg_shr, subcol_get_weight, subcol_get_filter, &
                                  is_filter_set, is_weight_set

      !-----------------------------------
      ! Average the subcolumns dimension (pcols*psubcols) to the grid dimension (pcols)
      !-----------------------------------

      type(physics_ptend), intent(in)             :: ptend_sc        ! intent in
      integer,  intent(in)                        :: ngrdcol       ! # grid cols
      integer,  intent(in)                        :: lchnk         ! chunk index
      type(physics_ptend), intent(inout)          :: ptend     

      if (subcol_tstcp_noAvg) then
         call subcol_ptend_get_firstsubcol(ptend_sc, .true., ngrdcol, lchnk, ptend)
      else
         call subcol_ptend_avg_shr(ptend_sc, ngrdcol, lchnk, ptend, is_filter_set(), is_weight_set())
      end if

end subroutine subcol_ptend_avg_tstcp
end module subcol_tstcp
