!-------------------------------------------------------------------
! Manages reading and interpolation of prescribed aerosol concentrations.
! This places the concentration fields in the physics buffer so that
! radiation package can access them.
!
! This has been generalized so that the field names in the data files
! and the field names in the physics buffer can be arbitrary.
!
! The prescribed_aero_specifier namelist variable specifies a list of 
! variable names of the concentration fields in the netCDF dataset (ncdf_fld_name)
! and the corresponding names used in the physics buffer:
!
! prescribed_aero_specifier = 'pbuf_name1:ncdf_fld_name1','pbuf_name2:ncdf_fld_name2', ...
!
! If there is no ":" then the specified name is used as both the 
! pbuf_name and ncdf_fld_name
!
! Created by: Francis Vitt
!-------------------------------------------------------------------
module prescribed_aero

  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils,   only : endrun
  use spmd_utils,   only : masterproc
  use tracer_data,  only : trfld, trfile
  use cam_logfile,  only : iulog
  use pio,          only : var_desc_t

  implicit none
  private
  save 

  type(trfld), pointer :: fields(:)
  type(trfile)         :: file

  public :: prescribed_aero_init
  public :: prescribed_aero_adv
  public :: write_prescribed_aero_restart
  public :: read_prescribed_aero_restart
  public :: has_prescribed_aero
  public :: prescribed_aero_register
  public :: init_prescribed_aero_restart
  public :: prescribed_aero_readnl

  logical :: has_prescribed_aero = .false.

  ! Decides if its a modal aerosol simulation or not
  logical :: clim_modal_aero     = .false. 

  integer, parameter, public :: N_AERO = 50

  integer :: number_flds

  character(len=256) :: filename = ' '
  character(len=256) :: filelist = ' '
  character(len=256) :: datapath = ' '
  character(len=32)  :: datatype = 'SERIAL'
  logical            :: rmv_file = .false.
  integer            :: cycle_yr = 0
  integer            :: fixed_ymd = 0
  integer            :: fixed_tod = 0
  character(len=32)  :: specifier(N_AERO) = ''

  ! prescribed aerosol names 
  character(len=16), allocatable :: pbuf_names(:)

  integer :: aero_cnt
  integer :: aero_cnt_c = 0 ! clound borne species count for modal aerosols (zero otherwise)

  ! Normal random number which persists from one tiemstep to the next
  real(r8) :: randn_persists

  ! Following definitions are added to
  ! allow randn_persists to persist during restart runs
  type(var_desc_t) :: randn_persists_desc
  character(len=16), parameter :: randn_persists_name = 'prescraero_randn'

contains

!-------------------------------------------------------------------
! registers aerosol fields to the phys buffer
!-------------------------------------------------------------------
  subroutine prescribed_aero_register()

    use ppgrid,         only: pver,pcols
    
    use physics_buffer, only : pbuf_add_field, dtype_r8
    integer :: i,idx

    if (has_prescribed_aero) then
       do i = 1,aero_cnt
          call pbuf_add_field(pbuf_names(i),'physpkg',dtype_r8,(/pcols,pver/),idx)
       enddo
    endif

  endsubroutine prescribed_aero_register

!-------------------------------------------------------------------
! parses prescribed_aero_specifier namelist option
!-------------------------------------------------------------------
  subroutine prescribed_aero_init()

    use tracer_data, only : trcdata_init
    use cam_history, only : addfld, phys_decomp
    use ppgrid,      only : pver
    use error_messages, only: handle_err
    use ppgrid,         only: pcols, pver, begchunk, endchunk
    use physics_buffer, only : physics_buffer_desc

    implicit none

    ! local vars
    character(len=32)  :: spec_a
    integer :: ndx, istat, i, i_c, j
    
    if ( has_prescribed_aero ) then
       if ( masterproc ) then
          write(iulog,*) 'aero is prescribed in :'//trim(filename)
       endif
    else
       return
    endif


    allocate (file%in_pbuf(size(specifier)))
    if (clim_modal_aero) then
      file%in_pbuf(:) = .false.
      do i = 1,N_AERO
          do j=1,aero_cnt
              if(specifier(i) .eq. pbuf_names(j)) then
                  file%in_pbuf(i) = .true.
                  exit
              endif
          enddo
      enddo
    else
      file%in_pbuf(:) = .true.
    endif
    call trcdata_init( specifier, filename, filelist, datapath, fields, file, &
                       rmv_file, cycle_yr, fixed_ymd, fixed_tod, datatype)
        
    number_flds = 0
    if (associated(fields)) number_flds = size( fields )

    if( number_flds < 1 ) then
       if ( masterproc ) then
          write(iulog,*) 'There are no prescribed aerosols'
          write(iulog,*) ' '
       endif
       return
    end if
    
    ! Following loop add fields for output. For modal aerosols, first the cld borne aersols
    ! are added and then their interstitial counterparts are added. The loop exits once all the cld borne
    ! aerosols (and their interstitial counterparts) are added.  Note that the units(fields(i)%units), 
    ! will be same for both interstitial and cloud borne species.
    ! All other aerosol treatments(bulk) are left unchanged
    i_c = 0
    fldloop:do i = 1,number_flds
       if(clim_modal_aero .and. index(trim(fields(i)%fldnam),'_c') > 1) then ! Only cloud borne species
          call addfld(trim(fields(i)%fldnam)//'_D',trim(fields(i)%units), pver, 'A', 'prescribed aero', phys_decomp )
          call spec_c_to_a(trim(fields(i)%fldnam),spec_a)
          call addfld(trim(spec_a)//'_D',trim(fields(i)%units), pver, 'A', 'prescribed aero', phys_decomp )
          i_c = i_c + 1
          if(i_c >= aero_cnt_c) exit fldloop
       else
          call addfld(trim(fields(i)%fldnam)//'_D',trim(fields(i)%units), pver, 'A', 'prescribed aero', phys_decomp )
       endif
    enddo fldloop

  end subroutine prescribed_aero_init

!-------------------------------------------------------------------
! reads namelist options
!-------------------------------------------------------------------
subroutine prescribed_aero_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand
   use rad_constituents, only: rad_cnst_get_info ! Added to query if it is a modal aero sim or not

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr, nmodes, aero_loop_end
   logical :: skip_spec
   character(len=*), parameter :: subname = 'prescribed_aero_readnl'

   character(len=32)  :: prescribed_aero_specifier(N_AERO)
   character(len=256) :: prescribed_aero_file
   character(len=256) :: prescribed_aero_filelist
   character(len=256) :: prescribed_aero_datapath
   character(len=32)  :: prescribed_aero_type
   logical            :: prescribed_aero_rmfile
   integer            :: prescribed_aero_cycle_yr
   integer            :: prescribed_aero_fixed_ymd
   integer            :: prescribed_aero_fixed_tod
   integer :: i,k

   namelist /prescribed_aero_nl/ &
      prescribed_aero_specifier, &
      prescribed_aero_file,      &
      prescribed_aero_filelist,  &
      prescribed_aero_datapath,  &
      prescribed_aero_type,      &
      prescribed_aero_rmfile,    &
      prescribed_aero_cycle_yr,  &
      prescribed_aero_fixed_ymd, &
      prescribed_aero_fixed_tod      
   !-----------------------------------------------------------------------------

   ! Initialize namelist variables from local module variables.
   prescribed_aero_specifier= specifier
   prescribed_aero_file     = filename
   prescribed_aero_filelist = filelist
   prescribed_aero_datapath = datapath
   prescribed_aero_type     = datatype
   prescribed_aero_rmfile   = rmv_file
   prescribed_aero_cycle_yr = cycle_yr
   prescribed_aero_fixed_ymd= fixed_ymd
   prescribed_aero_fixed_tod= fixed_tod

   ! Read namelist
   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'prescribed_aero_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, prescribed_aero_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(prescribed_aero_specifier,len(prescribed_aero_specifier(1))*N_AERO,     mpichar, 0, mpicom)
   call mpibcast(prescribed_aero_file,     len(prescribed_aero_file),     mpichar, 0, mpicom)
   call mpibcast(prescribed_aero_filelist, len(prescribed_aero_filelist), mpichar, 0, mpicom)
   call mpibcast(prescribed_aero_datapath, len(prescribed_aero_datapath), mpichar, 0, mpicom)
   call mpibcast(prescribed_aero_type,     len(prescribed_aero_type),     mpichar, 0, mpicom)
   call mpibcast(prescribed_aero_rmfile,   1, mpilog,  0, mpicom)
   call mpibcast(prescribed_aero_cycle_yr, 1, mpiint,  0, mpicom)
   call mpibcast(prescribed_aero_fixed_ymd,1, mpiint,  0, mpicom)
   call mpibcast(prescribed_aero_fixed_tod,1, mpiint,  0, mpicom)
#endif

   ! Update module variables with user settings.
   specifier  = prescribed_aero_specifier
   filename   = prescribed_aero_file
   filelist   = prescribed_aero_filelist
   datapath   = prescribed_aero_datapath
   datatype   = prescribed_aero_type
   rmv_file   = prescribed_aero_rmfile
   cycle_yr   = prescribed_aero_cycle_yr
   fixed_ymd  = prescribed_aero_fixed_ymd
   fixed_tod  = prescribed_aero_fixed_tod

   ! Turn on prescribed aerosols if user has specified an input dataset.
   has_prescribed_aero = len_trim(filename) > 0 

   if ( .not. has_prescribed_aero) return

   ! Determine whether its a 'modal' aerosol simulation  or not
   call rad_cnst_get_info(0, nmodes=nmodes)
   clim_modal_aero = (nmodes > 0)

   ! For modal aerosols, interstitial species(*_a) are diagnosed from
   ! their *_logm and *_logv counterparts (e.g. soa_a1 is diagnosed from
   ! soa_a1_logm and soa_a1_logv). Therefore, only *_logm and *_logv and cloud
   ! borne (*_c) species are specified in the build-namelist. In the following 
   ! cnt_loop, we will count the cloud borne species and *_logm species(in lieu 
   ! of *_a species). We will skip *_logv species. This will ensure that aero_cnt
   ! variable is the sum of cloud borne and interstitial species (which we will 
   ! manually add in pbuf_names array later). We are also counting cloud borne 
   ! (*_c) species which will help adding the same number of interstitial species
   ! in pbuf_names array
  
   ! determine which prescibed aerosols are specified
   aero_cnt   = 0
   aero_cnt_c = 0 ! cloud borne species count
   cnt_loop: do i = 1,N_AERO
      if ( len_trim(specifier(i))==0 ) exit cnt_loop
      skip_spec = .FALSE.
      if(clim_modal_aero) then
         ! For modal aerosols, skip counting species ending with *_logv 
         if(index(specifier(i),'_c') >= 1) aero_cnt_c = aero_cnt_c + 1
         if(index(specifier(i),'_logv') >= 1) skip_spec = .TRUE.
      endif
      if(.NOT.skip_spec)aero_cnt = aero_cnt+1
   end do cnt_loop

   has_prescribed_aero = aero_cnt>0
   if ( .not. has_prescribed_aero) return

   allocate(pbuf_names(aero_cnt))

   !  For modal aerosols, the following loop will add the cloud borne 
   ! species directly and interstitial species through the "add_interstitial_spec" 
   ! call. Interstitial species are added at the end of the cloud borne species in
   ! pbuf_names array.
   ! Note that aero_cnt_c should be zero for all other aerosol trearments 
   ! except the modal aerosols (e.g. bulk)
   
   if(.NOT. clim_modal_aero) aero_cnt_c = 0
   aero_loop_end = aero_cnt - aero_cnt_c
   
   do i = 1,aero_loop_end
      k = scan( specifier(i),':')
      if (k>1) then
         pbuf_names(i) = trim(specifier(i)(:k-1)) 
         if(clim_modal_aero)call add_interstitial_spec(aero_loop_end,i)
      else
         pbuf_names(i) = trim(specifier(i))
         if(clim_modal_aero)call add_interstitial_spec(aero_loop_end,i)
      endif
   enddo

end subroutine prescribed_aero_readnl

!-------------------------------------------------------------------
! Add interstitial aerosols in pbuf_names array for modal aerosols
!-------------------------------------------------------------------
subroutine add_interstitial_spec(aero_loop_end,i_in)
  implicit none
  
  !Arguments
  integer, intent(in) :: i_in, aero_loop_end
  
  !Local
  character(len=32) :: spec_a

  !  Replace 'c' with 'a' in species name
  call spec_c_to_a(pbuf_names(i_in), spec_a)
  pbuf_names(aero_loop_end+i_in) = spec_a
end subroutine add_interstitial_spec

!-------------------------------------------------------------------
! A generic subroutine which replaces 'c' in the cloud borne 
! species name with 'a' to make it interstital species
!-------------------------------------------------------------------
subroutine spec_c_to_a (spec_c_in,spec_a_out)
  implicit none

  !Arguments
  character(len=*), intent(in)  :: spec_c_in 
  character(len=*), intent(out) :: spec_a_out

  !Local
  character(len=32)   :: name
  character(len=1000) :: errMsg
  integer             :: k_c, k_cp1

  k_c = index(trim(adjustl(spec_c_in)),'_c')
  if(k_c >= 1) then
     k_cp1             = k_c + 1
     name              = trim(adjustl(spec_c_in))
     name(k_cp1:k_cp1) = 'a'
     spec_a_out        = trim(adjustl(name))
  else
     write(errMsg,*) "Input string (",trim(spec_c_in)," is not a cld borne aerosol,", &
          "cannot form interstitial species name"
     call endrun(trim(errMsg))
  endif  
end subroutine spec_c_to_a

!-------------------------------------------------------------------
! advances the aerosol fields to the current time step
!-------------------------------------------------------------------
  subroutine prescribed_aero_adv( state, pbuf2d )

    use tracer_data,  only : advance_trcdata
    use physics_types,only : physics_state
    use ppgrid,       only : begchunk, endchunk
    use ppgrid,       only : pcols, pver
    use string_utils, only : to_lower, GLC
    use cam_history,  only : outfld
    
    use physics_buffer, only : physics_buffer_desc, pbuf_get_field, pbuf_get_chunk, &
         pbuf_get_index

    implicit none

    type(physics_state), intent(in)    :: state(begchunk:endchunk)                 
    
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    type(physics_buffer_desc), pointer :: pbuf_chnk(:)

    character(len=32) :: spec_a
    integer :: c,ncol, i, i_c, spec_a_ndx, errcode
    real(r8),pointer :: outdata(:,:)
    logical :: cld_borne_aero = .FALSE.

    if( .not. has_prescribed_aero ) return

    call advance_trcdata( fields, file, state, pbuf2d )

    ! Diagnose interstital species using random sampling
    if ( clim_modal_aero ) then
      call rand_sample_prescribed_aero(state, pbuf2d )
    endif
    
    ! set the tracer fields with the correct units
    i_c = 0
    fldloop : do i = 1,number_flds
       
!$OMP PARALLEL DO PRIVATE (C, NCOL, OUTDATA,PBUF_CHNK)
       do c = begchunk,endchunk
          ncol = state(c)%ncol
          pbuf_chnk => pbuf_get_chunk(pbuf2d, c)
          if(clim_modal_aero .and. index(trim(fields(i)%fldnam),'_c') > 1) then ! Only cloud borne species
             call pbuf_get_field(pbuf_chnk, fields(i)%pbuf_ndx, outdata)
             call outfld( trim(fields(i)%fldnam)//'_D', outdata(:ncol,:), ncol, state(c)%lchnk )
             
             call spec_c_to_a(trim(fields(i)%fldnam),spec_a)
             spec_a_ndx = pbuf_get_index(trim(fields(i)%fldnam),errcode)
             call pbuf_get_field(pbuf_chnk, spec_a_ndx, outdata)
             call outfld( trim(spec_a)//'_D', outdata(:ncol,:), ncol, state(c)%lchnk )
             cld_borne_aero = .TRUE.
          else
             call pbuf_get_field(pbuf_chnk, fields(i)%pbuf_ndx, outdata)
             call outfld( trim(fields(i)%fldnam)//'_D', outdata(:ncol,:), ncol, state(c)%lchnk )
          endif
       enddo
       if(cld_borne_aero)then
          i_c = i_c + 1
          cld_borne_aero = .FALSE.
          if(i_c >= aero_cnt_c) exit fldloop
       endif
    enddo fldloop

  end subroutine prescribed_aero_adv

!-------------------------------------------------------------------
  subroutine rand_sample_prescribed_aero(state, pbuf2d)

    !Purpose: Generates log normal distribution for the interstitial species.
    !Note that only the interstitial aerosols are diagnosed here
    !
    !Written by Balwinder Singh (12/14/2012)
    !Adapted from Jin-Ho Yoon
    !
    !Update log:

    use physics_types,  only : physics_state
    use ppgrid,         only : begchunk, endchunk, pver
    use physics_buffer, only : physics_buffer_desc, pbuf_get_field, pbuf_get_chunk, &
         pbuf_get_index
    use ppgrid,       only : pcols, pver

    !Arguments
    type(physics_state), intent(in)    :: state(begchunk:endchunk)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
    
    !Local
    real(r8), parameter :: mean_max_val = 5.0_r8
    real(r8), parameter :: std_max_val  = 3.0_r8

    integer  :: c, ncol, i, kc, klog, logv_ndx, logm_ndx, a_ndx, icol, kpver
    real(r8) :: logm2, variance, std, mean_max, std_max, randn
    type(physics_buffer_desc), pointer :: pbuf_chnk(:)
    real(r8), pointer :: data(:,:)
    real(r8) :: logm_data(pcols,pver),logv_data(pcols,pver)

    randn = randn_prescribed_aero()
    do i = 1, aero_cnt
       !Species with '_a' are updated using random sampling.
       !Cloud borne species (ends with _c1,_c2, etc.) have to be skipped
       kc   = index(pbuf_names(i),'_c')
       ! if kc>1, species is cloud borne species
       if( kc > 1) cycle
       
       !If species is ending with _a1, a2 etc., then find indicies of their 
       !logv and lom conterparts
       logv_ndx = logvm_get_index(pbuf_names(i),'v')
       logm_ndx = logvm_get_index(pbuf_names(i),'m')
       a_ndx    = pbuf_get_index(trim(adjustl(pbuf_names(i))))

       do c = begchunk, endchunk
          ncol = state(c)%ncol
          pbuf_chnk => pbuf_get_chunk(pbuf2d, c)

          !Get the species mixing ratio nad its logv and logm values
          call pbuf_get_field(pbuf_chnk, a_ndx, data)
          logv_data = fields(logv_ndx)%data(:,:,c)
          logm_data = fields(logm_ndx)%data(:,:,c)

          do icol = 1, ncol
             do kpver = 1, pver
                logm2    = logm_data(icol,kpver) * logm_data(icol,kpver)

                !Compute variance
                variance = max( 0.0_r8, (logv_data(icol,kpver) - logm2) )

                !Standard deviation
                std      = sqrt(variance)

                !Bounds to keep mixing ratios from going unphysical
                mean_max = exp(logm_data(icol,kpver)) * mean_max_val
                std_max  = exp(logm_data(icol,kpver)  + std_max_val * std )

                data(icol,kpver) = min(exp(logm_data(icol,kpver)+randn*std),mean_max,std_max)

             enddo !pver
          enddo !col
       enddo !chunk
    enddo !flds
   
  end subroutine rand_sample_prescribed_aero
!-------------------------------------------------------------------
  function randn_prescribed_aero()
    !Pupose: This function generates a new normally distributed random
    ! number at end end of each day. This random number then stays the same
    ! for the whole day
    !
    !Written by Balwinder Singh (12/14/2012)
    !Adapted from Jin-Ho Yoon
    !
    !Update log:
    
    use time_manager,   only : is_end_curr_day, is_first_step, get_nstep
    
    integer, parameter :: rconst1_1 = 5000000
    integer, parameter :: rconst1_2 = 50
    integer, parameter :: rconst2_1 = 10000000
    integer, parameter :: rconst2_2 = 10
    
    integer :: i, seed_size, nstep
    integer, allocatable :: seed(:)

    real(r8) :: randn_prescribed_aero
    real(r8) :: randu1, randu2
    
    !Use same random number for the entire day and generate a new normally 
    !distributed random number at the start of the new day
    if(is_first_step() .or. is_end_curr_day()) then
       !Generate two uniformly distributed random numbers (between 0 and 1)
       call random_seed(size=seed_size)
       allocate(seed(seed_size))

       ! Using nstep as a seed to generate same sequence
       nstep = get_nstep()
       do i = 1, seed_size
          seed(i) = rconst1_1*nstep + rconst1_2*(i-1)
       end do
       call random_seed(put=seed)
       call random_number (randu1)

       do i = 1, seed_size
          seed(i) = rconst2_1*nstep + rconst2_2*(i-1)
       end do
       call random_seed(put=seed)
       call random_number (randu2)
       deallocate(seed)
      
       !Normal distribution (Mean = 0, standard dev = 1)
       randn_prescribed_aero = boxMuller(randu1,randu2)
       randn_persists = randn_prescribed_aero
    else
       !Use the previously generated random number
       randn_prescribed_aero = randn_persists
    endif
  end function randn_prescribed_aero
!-------------------------------------------------------------------
  function logvm_get_index(name,type) result(index)
    implicit none
        
    !Args
    character(len=*), intent(in) :: name, type    

    !Loc
    character(len=64)   :: tmp_name
    character(len=1000) :: msgStr
    integer :: index, i
    
    
    index = -1
    tmp_name = trim(adjustl(name))//'_log'//trim(adjustl(type))
           
    fldloop: do i = 1, number_flds
       if(fields(i)%fldnam == tmp_name) then
          index = i       
          exit fldloop
       endif
    enddo fldloop
    
    if(index == -1) then
       write(msgStr,*) "Prescribed_aero.F90: ",tmp_name," doesn't exist in the fields%fldnam"
       call endrun(msgStr)
    endif
    
  end function logvm_get_index
!-------------------------------------------------------------------
  function boxMuller(ru1,ru2) result(rn)
    use physconst,      only : pi 
    implicit none

    !Args
    real(r8), intent(in)  :: ru1, ru2

    !Loc
    real(r8), parameter :: pi2 = 2._r8 * pi    
    real(r8) :: ur, theta, rn

    !Based on Box Muller Method, generate normally distributed random numbers 
    ur    = sqrt(-2._r8 * log(ru1))
    theta = pi2 * ru2
    
    !Normal distribution (Mean = 0, standard dev = 1)
    rn = ur * cos(theta)
    
  end function boxMuller

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine init_prescribed_aero_restart( piofile )
    use pio, only : file_desc_t, pio_def_var, pio_double 
    use tracer_data, only : init_trc_restart
    implicit none
    type(file_desc_t),intent(inout) :: pioFile     ! pio File pointer
    integer :: ierr

    ! For allowing randn_persists to persist during reststarts
    ierr = pio_def_var(piofile, randn_persists_name, pio_double, randn_persists_desc)

    call init_trc_restart( 'prescribed_aero', piofile, file )

  end subroutine init_prescribed_aero_restart
!-------------------------------------------------------------------
  subroutine write_prescribed_aero_restart( piofile )
    use tracer_data, only : write_trc_restart
    use pio, only : file_desc_t, pio_put_var 
    implicit none

    type(file_desc_t) :: piofile
    integer :: ierr

    !  For allowing randn_persists to persist during reststarts
    ierr = pio_put_var(piofile, randn_persists_desc, (/randn_persists/)) 

    call write_trc_restart( piofile, file )

  end subroutine write_prescribed_aero_restart

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine read_prescribed_aero_restart( pioFile )
    use tracer_data, only : read_trc_restart
    use pio, only : file_desc_t, pio_inq_varid, pio_get_var 
    implicit none

    type(file_desc_t) :: piofile    
    integer :: ierr

    ! For allowing randn_persists to persist during reststarts
    ierr = pio_inq_varid(pioFile, randn_persists_name, randn_persists_desc)
    ierr = pio_get_var(pioFile, randn_persists_desc, randn_persists)

    call read_trc_restart( 'prescribed_aero', piofile, file )

  end subroutine read_prescribed_aero_restart

end module prescribed_aero
