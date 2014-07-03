module mo_extfrc
  !---------------------------------------------------------------
  ! 	... insitu forcing module
  !---------------------------------------------------------------

  use shr_kind_mod, only : r8 => shr_kind_r8
  use ppgrid,       only : pcols, begchunk, endchunk, pver, pverp
  use chem_mods,    only : gas_pcnst, extcnt
  use spmd_utils,   only : masterproc,iam
  use abortutils,   only : endrun
  use cam_history,  only : addfld, outfld, phys_decomp, add_default
  use cam_logfile,  only : iulog
  use tracer_data,  only : trfld,trfile

  implicit none

  type :: forcing
     integer           :: frc_ndx
     real(r8)              :: mw
     character(len=265) :: filename
     real(r8), pointer     :: times(:)
     real(r8), pointer     :: levi(:)
     character(len=8)  :: species
     character(len=8)  :: units
     integer                   :: nsectors
     character(len=32),pointer :: sectors(:)
     type(trfld), pointer      :: fields(:)
     type(trfile)              :: file
  end type forcing

  private
  public  :: extfrc_inti
  public  :: extfrc_set
  public  :: extfrc_timestep_init

  save

  integer, parameter :: time_span = 1

  character(len=256) ::   filename

  logical :: has_extfrc(gas_pcnst)
  type(forcing), allocatable  :: forcings(:)
  integer :: extfrc_cnt = 0

contains

  subroutine extfrc_inti( extfrc_specifier, extfrc_type, extfrc_cycle_yr, extfrc_fixed_ymd, extfrc_fixed_tod)

    !-----------------------------------------------------------------------
    ! 	... initialize the surface forcings
    !-----------------------------------------------------------------------
    use cam_pio_utils, only : cam_pio_openfile
    use pio, only : pio_inq_dimid, pio_inquire, pio_inq_varndims, pio_closefile, &
         pio_inq_varname, pio_nowrite, file_desc_t
    use mo_tracname,   only : solsym
    use mo_chem_utls,  only : get_extfrc_ndx, get_spc_ndx
    use chem_mods,     only : frc_from_dataset
    use tracer_data,   only : trcdata_init
    use phys_control,  only : phys_getopts
    use physics_buffer, only : physics_buffer_desc

    implicit none

    !-----------------------------------------------------------------------
    ! 	... dummy arguments
    !-----------------------------------------------------------------------
    character(len=*), dimension(:), intent(in) :: extfrc_specifier
    character(len=*), intent(in) :: extfrc_type
    integer  , intent(in)        :: extfrc_cycle_yr
    integer  , intent(in)        :: extfrc_fixed_ymd
    integer  , intent(in)        :: extfrc_fixed_tod

    !-----------------------------------------------------------------------
    ! 	... local variables
    !-----------------------------------------------------------------------
    integer :: astat
    integer :: j, l, m, n, i,mm                          ! Indices
    character(len=16)  :: species
    character(len=16)  :: spc_name
    character(len=256) :: locfn
    character(len=256) :: spc_fnames(gas_pcnst)

    integer ::  vid, ndims, nvars, isec, ierr
    type(file_desc_t) :: ncid
    character(len=32)  :: varname

    character(len=1), parameter :: filelist = ''
    character(len=1), parameter :: datapath = ''
    logical         , parameter :: rmv_file = .false.
    logical  :: history_aerosol      ! Output the MAM aerosol tendencies

    !-----------------------------------------------------------------------
 
    call phys_getopts( history_aerosol_out        = history_aerosol   )

    do i = 1, gas_pcnst
       has_extfrc(i) = .false.
       spc_fnames(i) = ''
    enddo

    !-----------------------------------------------------------------------
    ! 	... species has insitu forcing ?
    !-----------------------------------------------------------------------

    !write(iulog,*) 'Species with insitu forcings'

    count_emis: do n=1,gas_pcnst

       if ( len_trim(extfrc_specifier(n) ) == 0 ) then
          exit count_emis
       endif

       i = scan(extfrc_specifier(n),'->')
       spc_name = trim(adjustl(extfrc_specifier(n)(:i-1)))
       filename = trim(adjustl(extfrc_specifier(n)(i+2:)))

       m = get_extfrc_ndx( spc_name )

       if ( m < 1 ) then
          call endrun('extfrc_inti: '//trim(spc_name)// ' does not have an external source')
       endif

       if ( .not. frc_from_dataset(m) ) then
          call endrun('extfrc_inti: '//trim(spc_name)//' cannot have external forcing from additional dataset')
       endif

       mm = get_spc_ndx(spc_name)
       spc_fnames(mm) = filename

       has_extfrc(mm) = .true.
       !write(iulog,*) '   ',  spc_name ,' : filename = ',trim(spc_fnames(mm)),' spc ndx = ',mm

    enddo count_emis

    extfrc_cnt = count( has_extfrc(:) )

    if( extfrc_cnt < 1 ) then
       if (masterproc) write(iulog,*) 'There are no species with insitu forcings'
       return
    end if

    if (masterproc) write(iulog,*) ' '

    !-----------------------------------------------------------------------
    ! 	... allocate forcings type array
    !-----------------------------------------------------------------------
    allocate( forcings(extfrc_cnt), stat=astat )
    if( astat/= 0 ) then
       write(iulog,*) 'extfrc_inti: failed to allocate forcings array; error = ',astat
       call endrun
    end if

    !-----------------------------------------------------------------------
    ! 	... setup the forcing type array
    !-----------------------------------------------------------------------
    n = 0
    species_loop : do m = 1,gas_pcnst
       has_forcing : if( has_extfrc(m) ) then
          spc_name = trim( solsym(m) )
          n        = n + 1
          !-----------------------------------------------------------------------
          ! 	... default settings
          !-----------------------------------------------------------------------
          forcings(n)%frc_ndx          = get_extfrc_ndx( spc_name )
          forcings(n)%species          = spc_name
          forcings(n)%filename         = spc_fnames(m)
          call addfld( trim(spc_name)//'_XFRC',  'molec/cm3/s', pver, 'A', &
                       'external forcing for '//trim(spc_name),   phys_decomp )
          call addfld( trim(spc_name)//'_CLXF',  'molec/cm2/s', 1, 'A', &
                       'vertically intergrated external forcing for '//trim(spc_name),   phys_decomp )
          if ( history_aerosol ) then 
             call add_default( trim(spc_name)//'_XFRC', 1, ' ' )
             call add_default( trim(spc_name)//'_CLXF', 1, ' ' )
          endif
       end if has_forcing
    end do species_loop

    if (masterproc) then
       !-----------------------------------------------------------------------
       ! 	... diagnostics
       !-----------------------------------------------------------------------
       write(iulog,*) ' '
       write(iulog,*) 'extfrc_inti: diagnostics'
       write(iulog,*) ' '
       write(iulog,*) 'extfrc timing specs'
       write(iulog,*) 'type = ',extfrc_type
       if( extfrc_type == 'FIXED' ) then
          write(iulog,*) ' fixed date = ', extfrc_fixed_ymd
          write(iulog,*) ' fixed time = ', extfrc_fixed_tod
       else if( extfrc_type == 'CYCLICAL' ) then
          write(iulog,*) ' cycle year = ',extfrc_cycle_yr
       end if
       write(iulog,*) ' '
       write(iulog,*) 'there are ',extfrc_cnt,' species with external forcing files'
       do m = 1,extfrc_cnt
          write(iulog,*) ' '
          write(iulog,*) 'forcing type ',m
          write(iulog,*) 'species = ',trim(forcings(m)%species)
          write(iulog,*) 'frc ndx = ',forcings(m)%frc_ndx
          write(iulog,*) 'filename= ',trim(forcings(m)%filename)
       end do
       write(iulog,*) ' '
    endif

    !-----------------------------------------------------------------------
    ! read emis files to determine number of sectors
    !-----------------------------------------------------------------------
    frcing_loop: do m = 1, extfrc_cnt

       forcings(m)%nsectors = 0

       call cam_pio_openfile ( ncid, trim(forcings(m)%filename), PIO_NOWRITE)
       ierr = pio_inquire (ncid, nVariables=nvars)

       do vid = 1,nvars

          ierr = pio_inq_varndims (ncid, vid, ndims)

          if( ndims < 4 ) then
             cycle
          elseif( ndims > 4 ) then
             ierr = pio_inq_varname (ncid, vid, varname)
             write(iulog,*) 'extfrc_inti: Skipping variable ', trim(varname),', ndims = ',ndims, &
                  ' , species=',trim(forcings(m)%species)
             cycle
          end if

          forcings(m)%nsectors = forcings(m)%nsectors+1

       enddo

       allocate( forcings(m)%sectors(forcings(m)%nsectors), stat=astat )
       if( astat/= 0 ) then
         write(iulog,*) 'extfrc_inti: failed to allocate forcings(m)%sectors array; error = ',astat
         call endrun
       end if

       isec = 1
       do vid = 1,nvars

          ierr = pio_inq_varndims (ncid, vid, ndims)
          if( ndims == 4 ) then
             ierr = pio_inq_varname(ncid, vid, forcings(m)%sectors(isec))
             isec = isec+1
          endif

       enddo

       call pio_closefile (ncid)

       allocate(forcings(m)%file%in_pbuf(size(forcings(m)%sectors)))
       forcings(m)%file%in_pbuf(:) = .false.
       call trcdata_init( forcings(m)%sectors, &
                          forcings(m)%filename, filelist, datapath, &
                          forcings(m)%fields,  &
                          forcings(m)%file, &
                          rmv_file, extfrc_cycle_yr, extfrc_fixed_ymd, extfrc_fixed_tod, extfrc_type)

    enddo frcing_loop


  end subroutine extfrc_inti

  subroutine extfrc_timestep_init( pbuf2d, state )
    !-----------------------------------------------------------------------
    !       ... check serial case for time span
    !-----------------------------------------------------------------------

    use physics_types,only : physics_state
    use ppgrid,       only : begchunk, endchunk
    use tracer_data,  only : advance_trcdata
    use physics_buffer, only : physics_buffer_desc

    implicit none

    type(physics_state), intent(in):: state(begchunk:endchunk)                 
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    !-----------------------------------------------------------------------
    !       ... local variables
    !-----------------------------------------------------------------------
    integer :: m

    do m = 1,extfrc_cnt
       call advance_trcdata( forcings(m)%fields, forcings(m)%file, state, pbuf2d  )
    end do

  end subroutine extfrc_timestep_init

  subroutine extfrc_set( lchnk, zint, frcing, ncol )
    !--------------------------------------------------------
    !	... form the external forcing
    !--------------------------------------------------------

    implicit none

    !--------------------------------------------------------
    !	... dummy arguments
    !--------------------------------------------------------
    integer,  intent(in)    :: ncol                  ! columns in chunk
    integer,  intent(in)    :: lchnk                 ! chunk index
    real(r8), intent(in)    :: zint(ncol, pverp)                  ! interface geopot above surface (km)
    real(r8), intent(inout) :: frcing(ncol,pver,extcnt)   ! insitu forcings (molec/cm^3/s)

    !--------------------------------------------------------
    !	... local variables
    !--------------------------------------------------------
    integer  ::  i, m, n
    character(len=16) :: xfcname
    real(r8) :: frcing_col(1:ncol)
    integer  :: k, isec
    real(r8),parameter :: km_to_cm = 1.e5_r8

#if (defined MODAL_AERO)
    real(r8) :: fr_to_ait, fr_to_acc, fr_to_pom, fr_to_cor
    real(r8) :: demis_ait, demis_acc, demis_pom, demis_cor
           ! emitted mass-mean diameter (m) for aitken, accum, coarse modes
           ! for log-normal, demis = dgn_emis * exp(1.5*(lnsigmag*2))
    real(r8) :: x_mton_ait, x_mton_acc, x_mton_pom, x_mton_cor
           ! [(number emissions)/(mass emissions)] for aitken, accum, coarse modes
           ! x_mton_xxx = 1/(emitted particle mean mass)  (#/kg)
    real(r8) :: om_to_oc
#endif

    if( extfrc_cnt < 1 .or. extcnt < 1 ) then
       return
    end if

    !--------------------------------------------------------
    !	... set non-zero forcings
    !--------------------------------------------------------
    src_loop : do m = 1,extfrc_cnt

!!$#ifdef DIAGS
!!$          write(iulog,*) ' '
!!$          write(iulog,*) 'src frcing'
!!$          write(iulog,'(1p5g15.7)') frc_data
!!$          write(iulog,*) 'model frcing'
!!$          write(iulog,'(1p5g15.7)') frc_model(:,1)
!!$          data_dz(:)  = forcings(m)%levi(2:nlev+1) - forcings(m)%levi(1:nlev)
!!$          model_dz(:) = model_z(2:pverp) - model_z(1:pver)
!!$          data_sum  = dot_product( frc_data,data_dz )
!!$          model_sum = dot_product( frc_model(:,1),model_dz )
!!$          write(iulog,*) 'data_sum, model_sum = ',data_sum, model_sum
!!$          write(iulog,*) '======================================================'
!!$          if ( abs( (data_sum-model_sum)/max(1.e-20,(data_sum+model_sum)) ) > 0.3 ) then
!!$             call endrun('EXTFRC_SET: data_sum and model_sum disagree')
!!$          endif
!!$#endif

       n = forcings(m)%frc_ndx

       frcing(:ncol,:,n) = 0._r8
       do isec = 1,forcings(m)%nsectors
          if (forcings(m)%file%alt_data) then
             frcing(:ncol,:,n) = frcing(:ncol,:,n) + forcings(m)%fields(isec)%data(:ncol,pver:1:-1,lchnk)
          else
             frcing(:ncol,:,n) = frcing(:ncol,:,n) + forcings(m)%fields(isec)%data(:ncol,:,lchnk)
          endif
       enddo

#if (defined MODAL_AERO)
   !! this is done here because so4_a1 and so4_a2 have the same emis files
   !   fr_to_ait = 0.01             ! fraction of SO4 mass emitted to aitken mode
   !   fr_to_acc = 1.0 - fr_to_ait  ! fraction of SO4 mass emitted to accu. mode
   !   demis_ait = 0.018e-6    ! m
   !   demis_acc = 0.11e-6     ! m
       om_to_oc  = 1.4_r8

   !   x_mton_ait = 6.0 * specmw_so4_amode /      &     ! kmol -> #
   !               (pi*specdens_so4_amode*(demis_ait**3))
   !   x_mton_acc = 6.0 * specmw_so4_amode /      &     ! kmol -> #
   !               (pi*specdens_so4_amode*(demis_acc**3))

       select case( forcings(m)%species )
   !   case( 'so4_a1' )
   !      frcing(:ncol,:,n) = frcing(:ncol,:,n) * fr_to_acc * 1._r8/2.5_r8 ! 0._r8
   !   case( 'so4_a2' )
   !      frcing(:ncol,:,n) = frcing(:ncol,:,n) * fr_to_ait * 1._r8/2.5_r8 ! 0._r8
   !   case( 'num_a1' )
   !      frcing(:ncol,:,n) = frcing(:ncol,:,n) * fr_to_acc * x_mton_acc * 1._r8/2.5_r8 ! 0._r8
   !   case( 'num_a2' )
   !      frcing(:ncol,:,n) = frcing(:ncol,:,n) * fr_to_ait * x_mton_ait * 1._r8/2.5_r8 ! 0._r8
#if ( defined MODAL_AERO_7MODE )
       case( 'pom_a3' )
          frcing(:ncol,:,n) = frcing(:ncol,:,n) * om_to_oc
#elif ( defined MODAL_AERO_3MODE )
       case( 'pom_a1' )
          frcing(:ncol,:,n) = frcing(:ncol,:,n) * om_to_oc
#endif
       end select
#endif

       xfcname = trim(forcings(m)%species)//'_XFRC'
       call outfld( xfcname, frcing(:ncol,:,n), ncol, lchnk )

       frcing_col(:ncol) = 0._r8
       do k = 1,pver
          frcing_col(:ncol) = frcing_col(:ncol) + frcing(:ncol,k,n)*(zint(:ncol,k)-zint(:ncol,k+1))*km_to_cm
       enddo
       xfcname = trim(forcings(m)%species)//'_CLXF'
       call outfld( xfcname, frcing_col(:ncol), ncol, lchnk )

    end do src_loop

  end subroutine extfrc_set


end module mo_extfrc
