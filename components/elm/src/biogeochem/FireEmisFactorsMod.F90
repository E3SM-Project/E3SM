module FireEmisFactorsMod
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: FireEmisFactorsMod
!
! !DESCRIPTION:
! Manages input of fire emissions factors from netCDF file
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils,   only : endrun
  use elm_varctl,   only : iulog
!
  implicit none
  private
  save
!
! !PUBLIC MEMBERS:
  public :: fire_emis_factors_init
  public :: fire_emis_factors_get

! !PRIVATE MEMBERS:
  integer :: npfts ! number of plant function types
!
  type emis_eff_t
     real(r8), pointer :: eff(:) ! emissions efficiency factor
     real(r8) :: wght            ! molecular weight
  endtype emis_eff_t
!
  type(emis_eff_t), pointer :: comp_factors_table(:)  ! hash table of FireEmis factors (points to an array of pointers)
  integer, pointer :: hash_table_indices(:)           ! pointer to hash table indices
  integer, parameter :: tbl_hash_sz = 2**16           ! hash table size
!
!
! !REVISION HISTORY:
!  28 Oct 2011: Created by Francis Vitt
!
!EOP
!-----------------------------------------------------------------------
contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: fire_emis_factors_get
!
! !INTERFACE:
  subroutine fire_emis_factors_get( comp_name, factors, molecwght )
!
! !DESCRIPTION:
! Method for getting FireEmis information for a named compound 
!
! !USES:
!    use pftconMod , only : nc3crop
    use pftvarcon , only : nc3crop
! !ARGUMENTS:
    character(len=*),intent(in)  :: comp_name      ! FireEmis compound name
    real(r8),        intent(out) :: factors(:)     ! vegetation type factors for the compound of interest
    real(r8),        intent(out) :: molecwght      ! molecular weight of the compound of intrest
!
!EOP
!-----------------------------------------------------------------------
! local vars:
    integer :: hashkey, ndx
    character(len=120) :: errmes

    hashkey = gen_hashkey(comp_name)
    ndx = hash_table_indices(hashkey)

    if (ndx<1) then 
       errmes = 'fire_emis_factors_get: '//trim(comp_name)//' compound not found in FireEmis table'
        write(iulog,*) trim(errmes)
       call endrun(errmes)
    endif

    factors(:npfts) = comp_factors_table( ndx )%eff(:npfts)
    if ( size(factors) > npfts )then
       factors(npfts+1:) = comp_factors_table( ndx )%eff(nc3crop)
    end if
    molecwght  = comp_factors_table( ndx )%wght

  end subroutine fire_emis_factors_get
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: fire_emis_factors_init
!
! !INTERFACE:
  subroutine fire_emis_factors_init( filename )
!
! !DESCRIPTION:
! Initializes the FireEmis factors using data from input file 
!
! !USES:
    use ncdio_pio, only : ncd_pio_openfile,ncd_inqdlen
    use pio, only : pio_inq_varid,pio_get_var,file_desc_t,pio_closefile
    use fileutils   , only : getfil
    use elm_varpar  , only : mxpft
!
! !ARGUMENTS:
    character(len=*),intent(in) :: filename ! FireEmis factors input file

!EOP
!-----------------------------------------------------------------------
!
    character(len=256) :: locfn           ! local file name
    type(file_desc_t) :: ncid             ! netcdf id

    integer :: start(2), count(2)

    integer :: ierr, i, vid
    integer :: dimid, n_comps, n_pfts
    integer :: comp_ef_vid,comp_name_vid,comp_mw_vid

    real(r8),          allocatable :: comp_factors(:)
    character(len=64), allocatable :: comp_names(:)     ! FireEmis compound names
    real(r8),          allocatable :: comp_molecwghts(:)! FireEmis compound molecular weights

    allocate(comp_factors_table(150))
    allocate(hash_table_indices(tbl_hash_sz))

    call getfil(filename, locfn, 0)
    call ncd_pio_openfile (ncid, trim(locfn), 0)

    call ncd_inqdlen( ncid, dimid, n_comps, name='Comp_Num')
    call ncd_inqdlen( ncid, dimid, n_pfts, name='PFT_Num')

    npfts = n_pfts
    if ( npfts /= mxpft .and. npfts /= 16 )then
       call endrun('Number of PFTs on fire emissions file is NOT correct. Its neither the total number of PFTS nor 16')
    end if

    ierr = pio_inq_varid(ncid,'Comp_EF',  comp_ef_vid)
    ierr = pio_inq_varid(ncid,'Comp_Name',comp_name_vid)
    ierr = pio_inq_varid(ncid,'Comp_MW',  comp_mw_vid)

    allocate( comp_factors(n_pfts) )
    allocate( comp_names(n_comps) )
    allocate( comp_molecwghts(n_comps) )

    ierr =  pio_get_var( ncid, comp_name_vid, comp_names ) 
    ierr =  pio_get_var( ncid, comp_mw_vid, comp_molecwghts ) 
 
    ! set up hash table where data is stored
    call  bld_hash_table_indices( comp_names )
    do i=1,n_comps
       start=(/i,1/)
       count=(/1,npfts/)
       ierr = pio_get_var( ncid, comp_ef_vid,  start, count, comp_factors )

       call enter_hash_data( trim(comp_names(i)), comp_factors, comp_molecwghts(i)  )
    enddo

    call pio_closefile(ncid)

    deallocate( comp_factors, comp_names, comp_molecwghts )

  endsubroutine fire_emis_factors_init
!-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  ! Private methods...

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  subroutine bld_hash_table_indices( names )
    character(len=*),intent(in) :: names(:)

    integer :: n, i, hashkey

    hash_table_indices(:) = 0

    n = size(names)
    do i=1,n
       hashkey = gen_hashkey(names(i))
       hash_table_indices(hashkey) = i
    enddo

  endsubroutine bld_hash_table_indices

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  subroutine enter_hash_data( name, data, molec_wght )
    character(len=*), intent(in) :: name
    real(r8), intent(in) :: data(:)
    real(r8), intent(in) :: molec_wght

    integer :: hashkey, ndx
    integer :: nfactors

    hashkey = gen_hashkey(name)
    nfactors = size(data)
    ndx = hash_table_indices(hashkey)

    if(ndx < 1) then
       call endrun('ndx out of bounds '//name)
    endif

    allocate (comp_factors_table(ndx)%eff(nfactors))

    comp_factors_table(ndx)%eff(:) = data(:)
    comp_factors_table(ndx)%wght = molec_wght

  end subroutine enter_hash_data

  !-----------------------------------------------------------------------
  !from cam_history
  !
  ! Purpose: Generate a hash key on the interval [0 .. tbl_hash_sz-1]
  !          given a character string.
  !
  ! Algorithm is a variant of perl's internal hashing function.
  !
  !-----------------------------------------------------------------------
  integer function gen_hashkey(string)

    implicit none
    !
    !  Arguments:
    !
    character(len=*), intent(in) :: string
    !
    !  Local vars
    !
    integer :: hash
    integer :: i
    integer :: strlen
    integer, parameter :: tbl_max_idx = 15  ! 2**N - 1
    integer, parameter :: gen_hash_key_offset = z'000053db'
    integer, dimension(0:tbl_max_idx) :: tbl_gen_hash_key =  (/61,59,53,47,43,41,37,31,29,23,17,13,11,7,3,1/)

    hash = gen_hash_key_offset
    strlen = len_trim(string)
    if ( strlen /= 19 ) then
       !
       ! Process arbitrary string length.
       !
       do i = 1, strlen
          hash = ieor(hash , (ichar(string(i:i)) * tbl_gen_hash_key(iand(i-1,tbl_max_idx))))
       end do
    else
       !
       ! Special case string length = 19
       !
       do i = 1, tbl_max_idx+1
          hash = ieor(hash , ichar(string(i:i))   * tbl_gen_hash_key(i-1)) 
       end do
       do i = tbl_max_idx+2, strlen
          hash = ieor(hash , ichar(string(i:i))   * tbl_gen_hash_key(i-tbl_max_idx-2)) 
       end do
    end if

    gen_hashkey = iand(hash, tbl_hash_sz-1)

    return

  end function gen_hashkey

end module FireEmisFactorsMod


