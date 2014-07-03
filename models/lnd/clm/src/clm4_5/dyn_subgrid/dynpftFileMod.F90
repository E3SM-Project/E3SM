module dynpftFileMod

#include "shr_assert.h"

  !---------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Handle reading of the pftdyn dataset, which specifies transient areas of natural PFTs
  !
  ! !USES:
  use clmtype
  use shr_kind_mod          , only : r8 => shr_kind_r8
  use decompMod             , only : bounds_type, BOUNDS_LEVEL_PROC
  use dynFileMod            , only : dyn_file_type
  use dynVarTimeInterpMod   , only : dyn_var_time_interp_type
  use clm_varctl            , only : iulog
  use abortutils            , only : endrun
  use spmdMod               , only : masterproc, mpicom
  use shr_log_mod           , only : errMsg => shr_log_errMsg
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  private
  save
  public :: dynpft_init     ! initialize information read from pftdyn dataset
  public :: dynpft_interp   ! interpolate pftdyn information to current time step
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: dynpft_check_consistency ! check consistency with surface dataset
  !
  ! ! PRIVATE TYPES
  type(dyn_file_type), target    :: dynpft_file   ! information for the pftdyn file
  type(dyn_var_time_interp_type) :: wtpft         ! weight of each pft relative to the natural veg landunit

  character(len=*), parameter    :: varname = 'PCT_NAT_PFT'  ! name of variable on file
  !---------------------------------------------------------------------------

contains


  !-----------------------------------------------------------------------
  subroutine dynpft_init(bounds)
    !
    ! !DESCRIPTION:
    ! Initialize dynamic pft dataset (position it to the right time samples
    ! that bound the initial model date)
    !
    ! This also calls dynpft_interp for the initial time
    !
    ! !USES:
    use clm_varctl  , only : fpftdyn
    use clm_varpar  , only : numpft, maxpatch_pft, natpft_size
    use ncdio_pio
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds  ! proc-level bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: wtpft_shape(2)                  ! shape of the wtpft data

    character(len= 32)     :: subname='dynpft_init'! subroutine name
    !-----------------------------------------------------------------------

    SHR_ASSERT(bounds%level == BOUNDS_LEVEL_PROC, subname // ': argument must be PROC-level bounds')

    ! Error check

    if ( maxpatch_pft /= numpft+1 )then
       call endrun(msg=' maxpatch_pft does NOT equal numpft+1 -- this is invalid for dynamic PFT case'//&
            errMsg(__FILE__, __LINE__) )
    end if

    if (masterproc) then
       write(iulog,*) 'Attempting to read pft dynamic landuse data .....'
    end if

    dynpft_file = dyn_file_type(fpftdyn)

    ! Consistency checks
    call check_dim(dynpft_file, 'natpft', natpft_size)
    call dynpft_check_consistency(bounds)

    ! read data PCT_NAT_PFT corresponding to correct year
    !
    ! Note: if you want to change PCT_NAT_PFT so that it is NOT interpolated, but instead
    ! jumps to each year's value on Jan 1 of that year, simply change wtpft to be of type
    ! dyn_var_time_uninterp_type (rather than dyn_var_time_interp_type), and change the
    ! following constructor to construct a variable of dyn_var_time_uninterp_type. That's
    ! all you need to do.

    wtpft_shape = [(bounds%endg-bounds%begg+1), natpft_size]
    wtpft = dyn_var_time_interp_type( &
         dyn_file=dynpft_file, varname=varname, &
         dim1name=grlnd, conversion_factor=100._r8, &
         do_check_sums_equal_1=.true., data_shape=wtpft_shape)

    call dynpft_interp(bounds)

  end subroutine dynpft_init

  !-----------------------------------------------------------------------
  subroutine dynpft_check_consistency(bounds)
    !
    ! !DESCRIPTION:
    ! Check consistency between dynpft file and surface dataset.
    !
    ! This is done by assuming that PCT_NAT_PFT at time 1 in the pftdyn file agrees with
    ! PCT_NAT_PFT on the surface dataset.
    !
    ! !USES:
    use clm_varsur, only : wt_nat_pft
    use clm_varpar, only : natpft_size
    use ncdio_pio
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds  ! proc-level bounds
    !
    ! !LOCAL VARIABLES:
    logical             :: check_dynpft_consistency ! whether to do the consistency check in this routine
    integer             :: g                        ! index
    real(r8), pointer   :: wtpft_time1(:,:)         ! weight of each pft in each grid cell at first time
    logical             :: readvar                  ! whether variable was read
    real(r8), parameter :: tol = 1.e-13_r8          ! tolerance for checking equality

    character(len=*), parameter :: subname = 'dynpft_check_consistency'
    !-----------------------------------------------------------------------
    
    call read_namelist

    if (check_dynpft_consistency) then

       ! Read first time slice of PCT_NAT_PFT

       allocate(wtpft_time1(bounds%begg:bounds%endg, natpft_size))
       call ncd_io(ncid=dynpft_file, varname=varname, flag='read', data=wtpft_time1, &
            dim1name=grlnd, nt=1, readvar=readvar)
       if (.not. readvar) then
          call endrun(msg=' ERROR: ' // trim(varname) // ' NOT on pftdyn file'//&
               errMsg(__FILE__, __LINE__))
       end if

       ! Convert from PCT to weight on grid cell
       wtpft_time1(bounds%begg:bounds%endg,:) = wtpft_time1(bounds%begg:bounds%endg,:) / 100._r8
    
       ! Compare with values read from surface dataset
       do g = bounds%begg, bounds%endg
          if (any(abs(wtpft_time1(g,:) - wt_nat_pft(g,:)) > tol)) then
             write(iulog,*) subname//' mismatch between PCT_NAT_PFT at initial time and that obtained from surface dataset'
             write(iulog,*) 'On pftdyn file: ', wtpft_time1(g,:)
             write(iulog,*) 'On surface dataset: ', wt_nat_pft(g,:)
             write(iulog,*) ' '
             write(iulog,*) 'Confirm that the year of your surface dataset'
             write(iulog,*) 'corresponds to the first year of your pftdyn file'
             write(iulog,*) '(e.g., for a pftdyn file starting at year 1850, which is typical,'
             write(iulog,*) 'you should be using an 1850 surface dataset),'
             write(iulog,*) 'and that your pftdyn file is compatible with the surface dataset.'
             write(iulog,*) ' '
             write(iulog,*) 'If you are confident that you are using the correct pftdyn file'
             write(iulog,*) 'and the correct surface dataset, then you can bypass this check by setting:'
             write(iulog,*) '  check_dynpft_consistency = .false.'
             write(iulog,*) 'in user_nl_clm'
             write(iulog,*) ' '
             call endrun(decomp_index=g, clmlevel=nameg, msg=errMsg(__FILE__, __LINE__))
          end if
       end do

       deallocate(wtpft_time1)

    end if

  contains
    !-----------------------------------------------------------------------
    subroutine read_namelist
      !
      ! !DESCRIPTION:
      ! Read namelist settings related to pftdyn consistency checks
      !
      ! !USES:
      use fileutils      , only : getavu, relavu
      use clm_nlUtilsMod , only : find_nlgroup_name
      use controlMod     , only : NLFilename
      use shr_mpi_mod    , only : shr_mpi_bcast
      !
      ! !ARGUMENTS:
      !
      ! !LOCAL VARIABLES:
      integer :: nu_nml    ! unit for namelist file
      integer :: nml_error ! namelist i/o error flag
      
      character(len=*), parameter :: subname = 'read_namelist'
      !-----------------------------------------------------------------------
      
      namelist /dynpft_consistency_checks/ &
           check_dynpft_consistency

      ! Set default namelist values
      check_dynpft_consistency = .true.

      ! Read namelist
      if (masterproc) then
         nu_nml = getavu()
         open( nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
         call find_nlgroup_name(nu_nml, 'dynpft_consistency_checks', status=nml_error)
         if (nml_error == 0) then
            read(nu_nml, nml=dynpft_consistency_checks,iostat=nml_error)
            if (nml_error /= 0) then
               call endrun(msg='ERROR reading dynpft_consistency_checks namelist'//errMsg(__FILE__, __LINE__))
            end if
         end if
         close(nu_nml)
         call relavu( nu_nml )
      endif

      call shr_mpi_bcast (check_dynpft_consistency, mpicom)

      if (masterproc) then
         write(iulog,*) ' '
         write(iulog,*) 'dynpft_consistency_checks settings:'
         write(iulog,nml=dynpft_consistency_checks)
         write(iulog,*) ' '
      end if

    end subroutine read_namelist


  end subroutine dynpft_check_consistency


  !-----------------------------------------------------------------------
  subroutine dynpft_interp(bounds)
    !
    ! !DESCRIPTION:
    ! Time interpolate dynamic pft data to get pft weights for model time
    !
    ! !USES:
    use clm_varcon      , only : istsoil
    use clm_varpar      , only : natpft_lb, natpft_ub
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds  ! proc-level bounds
    !
    ! !LOCAL VARIABLES:
    integer               :: m,p,l,g          ! indices
    real(r8), allocatable :: wtpft_cur(:,:)   ! current pft weights
    character(len=32) :: subname='dynpft_interp' ! subroutine name
    !-----------------------------------------------------------------------

    ! Interpolate pctpft to current time step - output in pctpft_intp
    ! Map interpolated pctpft to subgrid weights
    ! assumes that maxpatch_pft = numpft + 1, that each landunit has only 1 column, 
    ! SCAM and CNDV have not been defined, and create_croplandunit = .false.

    SHR_ASSERT(bounds%level == BOUNDS_LEVEL_PROC, subname // ': argument must be PROC-level bounds')

    ! Get pft weights for this time step

    call dynpft_file%update_time_info()

    allocate(wtpft_cur(bounds%begg:bounds%endg, natpft_lb:natpft_ub))
    call wtpft%get_current_data(wtpft_cur)

    do p = bounds%begp,bounds%endp
       g = pft%gridcell(p)
       l = pft%landunit(p)

       ! Note that we only deal with the istsoil landunit here, NOT the istcrop landunit
       ! (if there is one)
       ! (However, currently [as of 5-9-13] the code won't let you run with transient
       ! PFTs combined with create_crop_landunit anyway, so it's a moot point.)
       if (lun%itype(l) == istsoil) then
          m = pft%itype(p)

          ! Note that the following assignment assumes that all PFTs share a single column
          pft%wtcol(p) = wtpft_cur(g,m)
       end if

    end do

    deallocate(wtpft_cur)

  end subroutine dynpft_interp

end module dynpftFileMod
