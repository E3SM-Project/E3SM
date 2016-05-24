module MEGANFactorsMod
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: MEGANFactorsMod
!
! !DESCRIPTION:
! Manages input of MEGAN emissions factors from netCDF file
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils,   only : endrun
  use clm_varctl,   only : iulog
!
  implicit none
  private
  save
!
! !PUBLIC MEMBERS:
  public :: megan_factors_init
  public :: megan_factors_get
  public :: comp_names
!
! !PUBLIC DATA:
  real(r8), public, allocatable :: LDF(:)  ! light dependent fraction
  real(r8), public, allocatable :: Agro(:) ! growing leaf age factor
  real(r8), public, allocatable :: Amat(:) ! mature leaf age factor
  real(r8), public, allocatable :: Anew(:) ! new leaf age factor
  real(r8), public, allocatable :: Aold(:) ! old leaf age factor
  real(r8), public, allocatable :: betaT(:)! temperature factor
  real(r8), public, allocatable :: ct1(:)  ! temperature coefficient 1
  real(r8), public, allocatable :: ct2(:)  ! temperature coefficient 2
  real(r8), public, allocatable :: Ceo(:)  ! Eopt coefficient
!
! !PRIVATE MEMBERS:
  integer :: npfts ! number of plant function types
!
  type emis_eff_t
     real(r8), pointer :: eff(:) ! emissions efficiency factor
     real(r8) :: wght            ! molecular weight
     integer :: class_num        ! MEGAN class number
  endtype emis_eff_t
!
  type(emis_eff_t), pointer :: comp_factors_table(:)  ! hash table of MEGAN factors (points to an array of pointers)
  integer, pointer :: hash_table_indices(:)           ! pointer to hash table indices
  integer, parameter :: tbl_hash_sz = 2**16           ! hash table size
!
  character(len=32), allocatable :: comp_names(:)     ! MEGAN compound names
  real(r8),          allocatable :: comp_molecwghts(:)! MEGAN compound molecular weights
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
! !IROUTINE: megan_factors_get
!
! !INTERFACE:
  subroutine megan_factors_get( comp_name, factors, class_n, molecwght )
!
! !DESCRIPTION:
! Method for getting MEGAN information for a named compound 
!
! !ARGUMENTS:
    character(len=*),intent(in)  :: comp_name      ! MEGAN compound name
    real(r8),        intent(out) :: factors(npfts) ! vegitation type factors for the compound of intrest
    integer,         intent(out) :: class_n        ! MEGAN class number for the compound of intrest
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
       errmes = 'megan_factors_get: '//trim(comp_name)//' compound not found in MEGAN table'
        write(iulog,*) trim(errmes)
       call endrun(errmes)
    endif

    factors(:) = comp_factors_table( ndx )%eff(:)
    class_n    = comp_factors_table( ndx )%class_num
    molecwght  = comp_factors_table( ndx )%wght

  end subroutine megan_factors_get
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: megan_factors_init
!
! !INTERFACE:
  subroutine megan_factors_init( filename )
!
! !DESCRIPTION:
! Initializes the MEGAN factors using data from input file 
!
! !USES:
    use ncdio_pio, only : ncd_pio_openfile,ncd_inqdlen
    use pio, only : pio_inq_varid,pio_get_var,file_desc_t,pio_closefile
    use fileutils   , only : getfil
!
! !ARGUMENTS:
    character(len=*),intent(in) :: filename ! MEGAN factors input file

!EOP
!-----------------------------------------------------------------------
!
    character(len=256) :: locfn           ! local file name
    type(file_desc_t) :: ncid             ! netcdf id

    integer :: start(2), count(2)

    integer :: ierr, i, vid
    integer :: dimid, n_comps, n_classes, n_pfts
    integer :: class_ef_vid,comp_ef_vid,comp_name_vid,class_num_vid
    integer :: comp_mw_vid
    integer, allocatable :: class_nums(:)

    real(r8),allocatable :: factors(:)
    real(r8),allocatable :: comp_factors(:)
    real(r8),allocatable :: class_factors(:)

    allocate(comp_factors_table(150))
    allocate(hash_table_indices(tbl_hash_sz))


    call getfil(filename, locfn, 0)
    call ncd_pio_openfile (ncid, trim(locfn), 0)

    call ncd_inqdlen( ncid, dimid, n_comps, name='Comp_Num')
    call ncd_inqdlen( ncid, dimid, n_classes, name='Class_Num')
    call ncd_inqdlen( ncid, dimid, n_pfts, name='PFT_Num')

    npfts = n_pfts

    ierr = pio_inq_varid(ncid,'Class_EF', class_ef_vid)
    ierr = pio_inq_varid(ncid,'Comp_EF',  comp_ef_vid)
    ierr = pio_inq_varid(ncid,'Comp_Name',comp_name_vid)
    ierr = pio_inq_varid(ncid,'Class_Num',class_num_vid)
    ierr = pio_inq_varid(ncid,'Comp_MW',  comp_mw_vid)

    allocate( factors(n_pfts) )
    allocate( comp_factors(n_pfts) )
    allocate( class_factors(n_pfts) )

    allocate( comp_names(n_comps) )
    allocate( comp_molecwghts(n_comps) )
    allocate( class_nums(n_comps) )

    ierr =  pio_get_var( ncid, comp_name_vid, comp_names ) 
    ierr =  pio_get_var( ncid, comp_mw_vid, comp_molecwghts ) 
    ierr =  pio_get_var( ncid, class_num_vid, class_nums )
 
    ! set up hash table where data is stored
    call  bld_hash_table_indices( comp_names )
    do i=1,n_comps
       start=(/i,1/)
       count=(/1,16/)
       ierr = pio_get_var( ncid, comp_ef_vid,  start, count, comp_factors )
       start=(/class_nums(i),1/)
       ierr = pio_get_var( ncid, class_ef_vid, start, count, class_factors  )
       factors(:) = comp_factors(:)*class_factors(:)
       call enter_hash_data( trim(comp_names(i)), factors, class_nums(i), comp_molecwghts(i) )
    enddo

    allocate( LDF(n_classes) )
    allocate( Agro(n_classes) )
    allocate( Amat(n_classes) )
    allocate( Anew(n_classes) )
    allocate( Aold(n_classes) )
    allocate( betaT(n_classes) )
    allocate( ct1(n_classes) )
    allocate( ct2(n_classes) )
    allocate( Ceo(n_classes) )

    ierr = pio_inq_varid(ncid,'LDF', vid)
    ierr = pio_get_var( ncid, vid, LDF ) 

    ierr = pio_inq_varid(ncid,'Agro', vid)
    ierr = pio_get_var( ncid, vid, Agro ) 

    ierr = pio_inq_varid(ncid,'Amat', vid)
    ierr = pio_get_var( ncid, vid, Amat ) 

    ierr = pio_inq_varid(ncid,'Anew', vid)
    ierr = pio_get_var( ncid, vid, Anew ) 

    ierr = pio_inq_varid(ncid,'Aold', vid)
    ierr = pio_get_var( ncid, vid, Aold ) 

    ierr = pio_inq_varid(ncid,'betaT', vid)
    ierr = pio_get_var( ncid, vid, betaT ) 

    ierr = pio_inq_varid(ncid,'ct1', vid)
    ierr = pio_get_var( ncid, vid, ct1 ) 

    ierr = pio_inq_varid(ncid,'ct2', vid)
    ierr = pio_get_var( ncid, vid, ct2 ) 

    ierr = pio_inq_varid(ncid,'Ceo', vid)
    ierr = pio_get_var( ncid, vid, Ceo ) 

    call pio_closefile(ncid)

    deallocate( class_nums, comp_factors,class_factors,factors )

  endsubroutine megan_factors_init
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
  subroutine enter_hash_data( name, data, class_n, molec_wght )
    character(len=*), intent(in) :: name
    real(r8), intent(in) :: data(:)
    integer,  intent(in) :: class_n
    real(r8), intent(in) :: molec_wght

    integer :: hashkey, ndx
    integer :: nfactors

    hashkey = gen_hashkey(name)
    nfactors = size(data)

    ndx = hash_table_indices(hashkey)

    allocate (comp_factors_table(ndx)%eff(nfactors))

    comp_factors_table(ndx)%eff(:) = data(:)
    comp_factors_table(ndx)%class_num = class_n
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

    integer, parameter :: tbl_max_idx = 15  ! 2**N - 1
    integer, parameter :: gen_hash_key_offset = z'000053db'
    integer, dimension(0:tbl_max_idx) :: tbl_gen_hash_key =  (/61,59,53,47,43,41,37,31,29,23,17,13,11,7,3,1/)

    hash = gen_hash_key_offset

    if ( len(string) /= 19 ) then
       !
       ! Process arbitrary string length.
       !
       do i = 1, len(string)
          hash = ieor(hash , (ichar(string(i:i)) * tbl_gen_hash_key(iand(i-1,tbl_max_idx))))
       end do
    else
       !
       ! Special case string length = 19
       !
       do i = 1, tbl_max_idx+1
          hash = ieor(hash , ichar(string(i:i))   * tbl_gen_hash_key(i-1)) 
       end do
       do i = tbl_max_idx+2, len(string)
          hash = ieor(hash , ichar(string(i:i))   * tbl_gen_hash_key(i-tbl_max_idx-2)) 
       end do
    end if

    gen_hashkey = iand(hash, tbl_hash_sz-1)

    return

  end function gen_hashkey

end module MEGANFactorsMod


