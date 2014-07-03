module subcol
   !---------------------------------------------------------------------------
   ! Purpose:
   !
   ! Provides the infrastructure for handling the relationship between
   ! grid-box-averaged (GBA) fields and Sub-column (subcol) fields
   !
   ! Different methods for generating and averaging sub-column fields are
   ! called sub-column schemes and are designated with the 'subcol_scheme'
   ! namelist variable.
   !
   ! This module provides several public interfaces (see below) which operate
   ! based on the scheme. In order to implement a new scheme, you need to
   ! follow the following steps:
   ! I)  Implement a sub-column scheme in a separate module (subcol_<scheme>)
   ! I.1) Implement a subcol_register_<scheme> function to register any scheme-
   !      specific fields. Any pbuf_add_field and/or pbuf_register_subcol calls 
   !      need to go here.
   !      This step is optional
   ! I.2) Implement a subcol_init_<scheme> function to initialize any scheme-
   !      specific variables or fields.
   !      This step is optional
   ! I.3) Implement a subcol_gen_<scheme> function to generate the appropriate
   !      subcol fields based on the existing GBA fields and other scheme data.
   ! I.4) Implement a subcol_field_avg_<scheme> function to average subcol fields
   !      back into the appropriate GBA fields.
   !      This step is optional
   ! I.5) Implement a subcol_ptend_avg_<scheme> function to average the 
   !      sub-column ptend to a grid ptend
   !      This step is optional
   ! II) Add necessary cases in the master subcol module (this file)
   ! II.1) Add a case for your scheme name in subcol_register if you are calling
   !       your own subcol_<scheme> registration function.
   ! II.2) Add a case for your scheme name in subcol_init if you are calling your
   !       own subcol_<scheme> initialization function.
   ! II.3) Add a case for your scheme name in subcol_gen and call your
   !       subcol_<scheme> subcol generator.
   ! II.4) Add a case for your scheme name in subcol_field_avg if you are calling your
   !       own subcol_<scheme> field-averaging function.
   ! II.5) Add a case for your scheme name in subcol_ptend_avg if you are calling your
   !       own subcol_<scheme> ptend-averaging function.
   !
   ! New schemes should be implemented in a separate file which is used by
   ! this module (and thus may not 'use' any subcol module variable, function,
   ! or subroutine).
   !
   !---------------------------------------------------------------------------

   use shr_kind_mod,    only: r8=>shr_kind_r8, r4=>shr_kind_r4, i4=>shr_kind_i4
   use physics_types,   only: physics_state, physics_tend, physics_ptend
   use ppgrid,          only: pcols, psubcols, pver, pverp
   use abortutils,      only: endrun
   use subcol_utils,    only: subcol_field_avg_shr, subcol_ptend_avg_shr, &
                              subcol_field_get_firstsubcol, subcol_ptend_get_firstsubcol, &
                              is_filter_set, is_weight_set
   use subcol_tstcp   , only: subcol_gen_tstcp, subcol_register_tstcp, subcol_field_avg_tstcp, subcol_ptend_avg_tstcp

!   CloudObj, SILHS and vamp are currently being developed
!   References are being left in for convenience and demonstration purposes
!   use subcol_CloudObj, only: cloudobj_scheme_name, subcol_register_CloudObj, &
!                              subcol_init_CloudObj, subcol_gen_CloudObj
!   use subcol_CloudObj, only: subcol_ptend_avg_CloudObj
!   use subcol_SILHS,    only: subcol_register_SILHS, subcol_init_SILHS,       &
!                              subcol_gen_SILHS
!   use subcol_vamp,     only: subcol_gen_vamp, subcol_register_vamp, subcol_init_vamp

   implicit none

   private
   save

   !! Public interface functions which implement a sub-column scheme
   public :: subcol_register     ! Read scheme from namelist and initialize any scheme-global variables
   public :: subcol_init         ! Initialize any variables or fields specific to the active scheme
   public :: subcol_gen          ! Generate subcol fields from GBA fields
   public :: subcol_field_avg    ! Average subcol fields back into GBA fields
   public :: subcol_ptend_avg    ! Average sub-column ptend to grid ptend
   public :: subcol_readnl       ! Namelist reader for subcolumns
   public :: subcol_init_restart ! Initialize restart with subcolumn specific fields
   public :: subcol_read_restart ! Read subcolumn specific fields from restart
   public :: subcol_write_restart ! Write subcolumn specific fields for restart


   interface subcol_field_avg
      module procedure subcol_field_avg_1dr
      module procedure subcol_field_avg_1di
      module procedure subcol_field_avg_2dr
   end interface

contains
   subroutine subcol_readnl(nlfile)
      use subcol_utils,    only: subcol_get_scheme, subcol_utils_readnl
      use subcol_tstcp,    only: subcol_readnl_tstcp
!      use subcol_SILHS,    only: subcol_readnl_SILHS
!      use subcol_vamp,     only: subcol_readnl_vamp

      character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

      !
      ! Local variables
      !
      character(len=16) :: subcol_scheme_init          ! Name of subcolumn schem
      !-----------------------------------------------------------------------------

      call subcol_utils_readnl(nlfile)
      subcol_scheme_init = subcol_get_scheme()

      select case(trim(subcol_scheme_init))
      case('tstcp')
         call subcol_readnl_tstcp(nlfile)
!      case ('SILHS')
!         call subcol_readnl_SILHS(nlfile)
!      case (cloudobj_scheme_name)
!         call subcol_readnl_CloudObj(nlfile)
!      case ('vamp')
!         call subcol_readnl_vamp(nlfile)
      case ('off')
         ! No namelist for off 
      case default
         call endrun('subcol_register error: unsupported subcol_scheme specified')
      end select

   end subroutine subcol_readnl

   subroutine subcol_register()
      use phys_control,    only: phys_getopts
      use physics_buffer,  only: pbuf_add_field, dtype_i4
      use subcol_utils,    only: subcol_get_scheme

      select case(subcol_get_scheme())
         case('tstcp')
            call subcol_register_tstcp()
!         case ('SILHS')
!            call subcol_register_SILHS()
!         case (cloudobj_scheme_name)
!            call subcol_register_CloudObj()
!         case ('vamp')
!            call subcol_register_vamp()
         case ('off')
            ! No registration called
         case default
            call endrun('subcol_register error: unsupported subcol_scheme specified')
       end select

   end subroutine subcol_register

   subroutine subcol_init_restart(File, hdimids)
      use subcol_utils, only: subcol_utils_init_restart
      use pio,          only: file_desc_t

      type(file_desc_t),intent(in) :: File
      integer ,intent(in)          :: hdimids(:)

      call subcol_utils_init_restart(File, hdimids)
      ! Put scheme-specific calls here (in select statement)

   end subroutine subcol_init_restart

   subroutine subcol_write_restart(File)
      use subcol_utils, only: subcol_utils_write_restart
      use pio,          only: file_desc_t

      type(file_desc_t), intent(inout) :: File

      call subcol_utils_write_restart(File)
      ! Put scheme-specific calls here (in select statement)

   end subroutine subcol_write_restart

   subroutine subcol_read_restart(File)
      use subcol_utils, only: subcol_utils_read_restart
      use pio,          only: file_desc_t

      type(file_desc_t), intent(inout) :: File

      call subcol_utils_read_restart(File)
      ! Put scheme-specific calls here (in select statement)

   end subroutine subcol_read_restart

   subroutine subcol_init(pbuf2d, subcol_scheme_in)
      use physics_buffer,      only: physics_buffer_desc
      use cam_history_support, only: add_hist_coord
      use subcol_utils,        only: subcol_utils_init, subcol_get_scheme
      use time_manager,        only: is_first_step, is_first_restart_step

      type(physics_buffer_desc), pointer     :: pbuf2d(:,:)
      character(len=*), optional, intent(in) :: subcol_scheme_in ! Name of subcolumn generator

      !
      ! Local variables
      !
      character(len=16) :: subcol_scheme_init          ! Name of subcolumn scheme
  
      ! Set the subcol_scheme_gen to the one passed in , otherwise use the module scheme read from the namelist
      if (present(subcol_scheme_in)) then
         subcol_scheme_init = trim(subcol_scheme_in)
      else
         ! By default, use the module scheme read from the namelist
         subcol_scheme_init = subcol_get_scheme()
      end if
 
      if (is_first_step() .and. .not. is_first_restart_step()) then
         ! Initialize the subcol utility data only at the beginning of the run and not at restart
         call subcol_utils_init(subcol_scheme_init)
      end if

      ! Set the psubcols history coordinate for output
      if (trim(subcol_scheme_init) /= 'off') then
        call add_hist_coord('psubcols', psubcols, 'Subcolumn Index')  
      end if

      ! Call the appropriate subcol init method
      select case(trim(subcol_scheme_init))
        case ('tstcp')
           ! none needed for this scheme
!        case ('SILHS')
!           call subcol_init_SILHS(pbuf2d)
!        case (cloudobj_scheme_name)
!           call subcol_init_CloudObj(pbuf2d)
!        case ('vamp')
!           call subcol_init_vamp()
        case ('off')
           ! No initialization needed for off
        case default
           call endrun('subcol_init error: unsupported subcol_scheme specified')
      end select
   end subroutine subcol_init

   subroutine subcol_gen(state, tend, state_sc, tend_sc, pbuf, subcol_scheme_in)
     use physics_buffer,          only: physics_buffer_desc
     use subcol_utils,            only: subcol_get_scheme


      !-----------------------------------
      ! sub-column generator
      !-----------------------------------
      type(physics_state), intent(inout)     :: state
      type(physics_tend),  intent(inout)     :: tend
      type(physics_state), intent(inout)     :: state_sc         ! sub-column state
      type(physics_tend),  intent(inout)     :: tend_sc          ! sub-column tend
      type(physics_buffer_desc), pointer     :: pbuf(:)
      character(len=*), optional, intent(in) :: subcol_scheme_in ! Name of subcolumn generator

      !
      ! Local variables
      !
      character(len=16) :: subcol_scheme_gen  ! Name of subcolumn scheme

      ! Set the subcol_scheme_gen to the one passed in , otherwise use the module scheme read from the namelist
      if (present(subcol_scheme_in)) then
         subcol_scheme_gen = trim(subcol_scheme_in)
      else
         subcol_scheme_gen = subcol_get_scheme()
      end if

      if (subcol_scheme_gen /= 'off') then
         if (.not. allocated(state_sc%lat)) then
            call endrun('subcol_gen error: state_sc must be allocated before calling subcol_gen')
         end if
      end if

      if (state_sc%psetcols /= (pcols * psubcols)) then
        call endrun('subcol_gen error: state_sc%psetcols must be pcols * psubcols')
      end if

      select case(trim(subcol_scheme_gen))
         case('tstcp')
            call subcol_gen_tstcp(state, tend, state_sc, tend_sc, pbuf)
!         case ('SILHS')
!            call subcol_gen_SILHS(state, tend, state_sc, tend_sc, pbuf)
!         case (cloudobj_scheme_name)
!            call subcol_gen_CloudObj(state, tend, state_sc, tend_sc, pbuf)
!         case ('vamp')
!            call subcol_gen_vamp(state, tend, state_sc, tend_sc, pbuf)
         case default
            call endrun('subcol_gen error: unsupported subcol_scheme specified')
      end select

   end subroutine subcol_gen

   subroutine subcol_field_avg_1dr (field_sc, ngrdcol, lchnk, field, subcol_scheme_in)
      use physics_buffer,   only: physics_buffer_desc
      use subcol_utils,     only: subcol_get_scheme

      !-----------------------------------
      ! Average the subcolumns dimension (pcols*psubcols) to the grid dimension (pcols) 
      !-----------------------------------

      real(r8), intent(in)                        :: field_sc(:)   ! intent in
      integer,  intent(in)                        :: ngrdcol       ! # grid cols
      integer,  intent(in)                        :: lchnk         ! chunk index
      real(r8), intent(out)                       :: field(:)
      character(len=*), optional, intent(in)      :: subcol_scheme_in ! Name of subcolumn generator

      !
      ! Local variables
      !
      character(len=16) :: subcol_scheme_avg           ! Name of subcolumn scheme

      if (present(subcol_scheme_in)) then
         subcol_scheme_avg = trim(subcol_scheme_in)
      else
         ! By default, use the module scheme read from the namelist
         subcol_scheme_avg = subcol_get_scheme()
      end if

      if (size(field_sc,dim=1) .ne. pcols*psubcols) then
         call endrun('subcol_field_avg error: only fields with first dimension pcols*psubcols may use this routine')
      end if

      select case(trim(subcol_scheme_avg))
         ! Example of specialized averaging for specific subcolumn scheme
         case ('tstcp')
            call subcol_field_avg_tstcp(field_sc, ngrdcol, lchnk, field)
         ! Unless specialized averaging is needed, most subcolumn schemes will be handled by the default
         ! If filters and/or weights have been set, they are automatically used by this averager
         case default
            call subcol_field_avg_shr(field_sc, ngrdcol, lchnk, field, is_filter_set(), is_weight_set())
      end select

   end subroutine subcol_field_avg_1dr

   subroutine subcol_field_avg_1di (field_sc, ngrdcol, lchnk, field, subcol_scheme_in)
      use physics_buffer,   only: physics_buffer_desc
      use subcol_utils,     only: subcol_get_scheme

      !-----------------------------------
      ! Average the subcolumns dimension (pcols*psubcols) to the grid dimension (pcols) 
      !-----------------------------------

      integer,  intent(in)                        :: field_sc(:)   ! intent in
      integer,  intent(in)                        :: ngrdcol       ! # grid cols
      integer,  intent(in)                        :: lchnk         ! chunk index
      integer, intent(out)                        :: field(:)
      character(len=*), optional, intent(in)      :: subcol_scheme_in ! Name of subcolumn generator

      !
      ! Local variables
      !
      character(len=16) :: subcol_scheme_avg           ! Name of subcolumn scheme

      if (present(subcol_scheme_in)) then
         subcol_scheme_avg = trim(subcol_scheme_in)
      else
         ! By default, use the module scheme read from the namelist
         subcol_scheme_avg = subcol_get_scheme()
      end if

      if (size(field_sc,dim=1) .ne. pcols*psubcols) then
         call endrun('subcol_field_avg error: only fields with first dimension pcols*psubcols may use this routine')
      end if

      select case(trim(subcol_scheme_avg))
         ! Example of specialized averaging for specific subcolumn scheme
         case ('tstcp')
            call subcol_field_avg_tstcp(field_sc, ngrdcol, lchnk, field)
         ! Unless specialized averaging is needed, most subcolumn schemes will be handled by the default
         ! If filters and/or weights have been set, they are automatically used by this averager
         case default
            call subcol_field_avg_shr(field_sc, ngrdcol, lchnk, field, is_filter_set(), is_weight_set())
      end select

   end subroutine subcol_field_avg_1di

   subroutine subcol_field_avg_2dr (field_sc, ngrdcol, lchnk, field, subcol_scheme_in)
      use physics_buffer,   only: physics_buffer_desc
      use subcol_utils,     only: subcol_get_scheme

      !-----------------------------------
      ! Average the subcolumns dimension (pcols*psubcols) to the grid dimension (pcols) 
      !-----------------------------------

      real(r8), intent(in)                        :: field_sc(:,:) ! intent in
      integer,  intent(in)                        :: ngrdcol       ! # grid cols
      integer,  intent(in)                        :: lchnk         ! chunk index
      real(r8), intent(out)                       :: field(:,:)
      character(len=*), optional, intent(in)      :: subcol_scheme_in ! Name of subcolumn generator

      !
      ! Local variables
      !
      character(len=16) :: subcol_scheme_avg           ! Name of subcolumn scheme

      if (present(subcol_scheme_in)) then
         subcol_scheme_avg = trim(subcol_scheme_in)
      else
         ! By default, use the module scheme read from the namelist
         subcol_scheme_avg = subcol_get_scheme()
      end if

      if (size(field_sc,dim=1) .ne. pcols*psubcols) then
         call endrun('subcol_field_avg error: only fields with first dimension pcols*psubcols may use this routine')
      end if

      select case(trim(subcol_scheme_avg))
         ! Example of specialized averaging for specific subcolumn scheme
         case ('tstcp')
            call subcol_field_avg_tstcp(field_sc, ngrdcol, lchnk, field)
         ! Unless specialized averaging is needed, most subcolumn schemes will be handled with the default
         ! If filters and/or weights have been set, they are automatically used by this averager
         case default
            call subcol_field_avg_shr(field_sc, ngrdcol, lchnk, field, is_filter_set(), is_weight_set())
      end select

   end subroutine subcol_field_avg_2dr

   subroutine subcol_ptend_avg(ptend_sc, ngrdcol, lchnk, ptend, subcol_scheme_in)
      use physics_buffer, only: physics_buffer_desc
      use physics_types,  only: physics_ptend_init
      use subcol_utils,   only: subcol_get_scheme

      !-----------------------------------------------------------------------
      ! Average a sub-column ptend to a grid ptend
      !-----------------------------------------------------------------------

      type(physics_ptend),       intent(in)    :: ptend_sc    ! sub-column ptend 
      integer,                   intent(in)    :: ngrdcol     ! # grid cols
      integer,                   intent(in)    :: lchnk       ! chunk index
      type(physics_ptend),       intent(inout) :: ptend       ! grid ptend
      character(len=*), optional, intent(in)   :: subcol_scheme_in ! Name of subcolumn generator

      !
      ! Local variables
      !
      character(len=16) :: subcol_scheme_avg            ! Name of subcolumn scheme
      integer           :: indx, i, j

      if (present(subcol_scheme_in)) then
         subcol_scheme_avg = trim(subcol_scheme_in)
      else
         ! By default, use the module scheme read from the namelist
         subcol_scheme_avg = subcol_get_scheme()
      end if

      !-----------------------------------------------------------------------

      call physics_ptend_init(ptend, pcols, name=ptend_sc%name, ls=ptend_sc%ls, lu=ptend_sc%lu,  &
           lv=ptend_sc%lv, lq=ptend_sc%lq)

      select case(trim(subcol_scheme_avg))
         case ('tstcp')
            call subcol_ptend_avg_tstcp(ptend_sc, ngrdcol, lchnk, ptend)
         case default
            ! If filters and/or weights have been set, they are automatically used by this averager
            call subcol_ptend_avg_shr(ptend_sc, ngrdcol, lchnk, ptend, is_filter_set(), is_weight_set())
      end select

   end subroutine subcol_ptend_avg

end module subcol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! subcol_field_avg_handler is an external routine used by outfld
!! It is outside the module to prevent a dependency loop
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! MAKE SURE TO UPDATE THE INTERFACE IN OUTFLD IF THIS INTERFACE IS CHANGED
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine subcol_field_avg_handler(idim, field_in, c, field_out)
   use shr_kind_mod,    only: r8 => shr_kind_r8
   use ppgrid,          only: pcols, psubcols
   use phys_grid,       only: get_ncols_p
   use subcol,          only: subcol_field_avg

   implicit none

   !! Dummy arguments
   integer,  intent(in)    :: idim
   real(r8), intent(in)    :: field_in(idim, *)
   integer,  intent(in)    :: c
   real(r8), intent(inout) :: field_out(:,:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !! Local variables
   real(r8), allocatable   :: field_sc(:,:)
   integer                 :: i, j, ngrdcol, dim2

   dim2 = size(field_out, 2)
   allocate(field_sc(pcols*psubcols, dim2))

   do j = 1, dim2
      do i = 1, idim
         field_sc(i, j) = field_in(i, j)
      end do
   end do
   if (idim < (pcols * psubcols)) then
      field_sc(idim+1:pcols*psubcols,:) = 0.0_r8
   end if

   ngrdcol = get_ncols_p(c)
   call subcol_field_avg(field_sc, ngrdcol, c, field_out)

   deallocate(field_sc)
end subroutine subcol_field_avg_handler
