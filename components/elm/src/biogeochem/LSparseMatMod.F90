module LSparseMatMod


  ! !USES:
  ! sparse matrix capability
  use bshr_kind_mod , only : r8 => shr_kind_r8
  use bshr_log_mod  , only : errMsg => shr_log_errMsg
  use elm_varctl    , only : iulog
implicit none
  private
  type, public :: spm_list_type
    real(r8) :: val
    integer  :: icol
    integer  :: irow
    type(spm_list_type), pointer :: next => null()
  end type spm_list_type

  type :: lom_type
  contains
    procedure :: calc_state_pscal
    procedure :: calc_reaction_rscal
    procedure :: apply_reaction_rscal
  end type lom_type
  type, public :: sparseMat_type
    real(r8), pointer :: val(:)      !nonzero val
    integer , pointer :: icol(:)     !col id for each val
    integer , pointer :: pB(:)       !first nonzero column id (as in val) in each row
    integer , pointer :: ncol(:)     !number of nonzero columns in each row
    integer :: szcol                 !number of columns of the matrix
    integer :: szrow
  contains
    procedure, public :: init
  end type sparseMat_type
  real(r8), parameter :: tiny_val=1.e-14_r8
  public :: spm_list_init, spm_list_insert, spm_list_to_mat
  public :: flux_correction
contains

!---------------------------------------------------------------
  function create_spm_type()
    !
    ! ! DESCRIPTION:
    !
    ! create an object of type sparseMat_type.
    ! Right now it is purposely empty
   class(sparseMat_type), pointer :: spm
   class(sparseMat_type), pointer ::  create_spm_type
   allocate(spm)
   create_spm_type => spm
  end function create_spm_type

!---------------------------------------------------------------
  subroutine init(this,nelms,nrows,szcol)
  !
  !DESCRIPTION
  ! initialize the sparse matrix
  implicit none
  class(sparseMat_type), intent(inout) :: this
  integer, intent(in) :: nelms
  integer, intent(in) :: nrows
  integer, intent(in) :: szcol

  allocate(this%val(nelms)) ;  this%val(:)  = 0._r8
  allocate(this%icol(nelms));  this%icol(:) = 0
  allocate(this%pB(nrows));    this%pB(:) = 0
  allocate(this%ncol(nrows));  this%ncol(:) = 0
  this%szcol = szcol
  this%szrow= nrows
  end subroutine init

!---------------------------------------------------------------
  subroutine spm_pack(Amat, nelms)
  !
  ! DESCRIPTION
  ! Pack a full matrix into a sparse matrix

  implicit none
  real(r8), dimension(:,:),intent(in) :: Amat
  integer, intent(in) :: nelms

  class(sparseMat_type), pointer :: spm
  integer :: nrows, szcol
  integer :: ii,jj,kk
  logical :: first
  nrows=size(Amat,1)
  szcol=size(Amat,2)
  allocate(spm, source=create_spm_type())

  call spm%init(nelms, nrows, szcol)

  !pack the nonzero element (defined with absolute values greater than tiny_val) into the sparse matrix
  kk = 0
  do ii = 1, nrows
    first =.true.
    do jj = 1, szcol
      if(abs(amat(ii,jj))>tiny_val)then
        kk = kk + 1
        spm%val(kk) = amat(ii,jj)
        spm%icol(kk)= jj
        if(first)then
          spm%pB(ii) = jj
        endif
        spm%ncol(ii)=spm%ncol(ii)+1
        first=.false.
      endif
    enddo
  enddo

  end subroutine spm_pack

!---------------------------------------------------------------

  subroutine spm_unpack(spm, bmat)

  ! DESCRIPTION
  ! unpack a sparse matrix spm into a full matrix bmat
  implicit none
  class(sparseMat_type), intent(in) :: spm
  real(r8), pointer :: bmat(:,:)

  integer :: nrs, ncs
  integer :: i, j, id

  nrs = size(spm%ncol)
  ncs = spm%szcol
  allocate(bmat(nrs,ncs))
  !unpack the matrix
  bmat(:,:)=0._r8
  do i = 1, nrs
    do j = 1, spm%ncol(i)
      id=j+spm%pB(i)
      bmat(i,spm%icol(id))=spm%val(id)
    enddo
  enddo

  end subroutine spm_unpack

!---------------------------------------------------------------

  subroutine spm_axpy(nx, ny, a, x, spm, y, errinfo)
  !
  ! DESCRIPTION
  ! y=y+a*spm*x
  ! with x, y being vectors, and a a scalar, and spm a sparse matrix. 
  implicit none
  real(r8), intent(in) :: a
  integer , intent(in) :: nx
  integer , intent(in) :: ny
  real(r8), intent(in) :: x(ny)
  type(sparseMat_type), intent(in) :: spm
  real(r8), intent(inout) :: y(nx)
  integer, intent(out) :: errinfo

  integer :: ii, jj, id
  real(r8) :: dy
  errinfo=0
  if(nx /=spm%szrow .or. ny /= spm%szcol)then
    errinfo=-1
    return
  endif

  do ii = 1, spm%szrow
    dy = 0._r8
    do jj = 1, spm%ncol(ii)
      id = spm%pB(ii)+jj-1
      dy=dy+spm%val(id)*x(spm%icol(id))
    enddo
    y(ii) = y(ii) + a*dy
  enddo
  end subroutine spm_axpy
!---------------------------------------------------------------


  subroutine spm_list_init(self, val, icol, irow, nelms)
  !
  ! DESCRIPTION
  ! initialize a sparse matrix
  implicit none
  type(spm_list_type), pointer :: self
  real(r8), intent(in) :: val
  integer, intent(in) :: icol
  integer, intent(in) :: irow
  integer, intent(out):: nelms
  allocate(self)
  nullify(self%next)

  self%val=val
  self%icol=icol
  self%irow=irow
  nelms=1
  end subroutine spm_list_init
!---------------------------------------------------------------

  subroutine spm_list_insert(self, val, icol, irow, nelms)
  ! 
  ! DESCRIPTION
  ! add a specified element to sparse matrix
  implicit none
  type(spm_list_type), pointer :: self
  real(r8), intent(in) :: val
  integer, intent(in) :: icol
  integer, intent(in) :: irow
  integer, intent(inout):: nelms
  type(spm_list_type), pointer :: next

  allocate(next)

  next%val=val
  next%icol=icol
  next%irow=irow
  next%next=> self
  self => next
  nelms=nelms+1
  end subroutine spm_list_insert
!---------------------------------------------------------------

  function spm_list_next(self)result(next)

  implicit none
  type(spm_list_type), pointer :: self
  type(spm_list_type), pointer :: next

  next => self%next
  end function spm_list_next
!---------------------------------------------------------------

  subroutine spm_list_free(self)
  !
  ! DESCRIPTION
  ! Free the allocated memory associated with the sparse matrix
  implicit none
  type(spm_list_type), pointer :: self
  type(spm_list_type), pointer :: current
  type(spm_list_type), pointer :: elem

  elem => self
  do while(associated(elem))
    current => elem
    elem => current%next
    deallocate(current)
  enddo
  end subroutine spm_list_free
!---------------------------------------------------------------

  subroutine spm_list_to_mat(spm_list, spm, nelms, szcol)
  !
  !copy the list to sparse matrix, and free the list
  implicit none
  type(spm_list_type), pointer :: spm_list
  class(sparseMat_type), pointer :: spm
  integer, intent(in) :: nelms    !total number of nonzero elements
  integer, intent(in) :: szcol    !total columns of the matrix

  type(spm_list_type), pointer :: next
  integer :: nrows
  integer :: irow, jj

  nrows=spm_list%irow

  allocate(spm, source=create_spm_type())

  call spm%init(nelms, nrows, szcol)

  jj = nelms
  irow=nrows
  next => spm_list
  do while(associated(next))
    if (irow /= next%irow)then
       spm%pB(irow)=jj+1
       irow=next%irow
    endif
    spm%val(jj) = next%val
    spm%icol(jj)= next%icol
    spm%ncol(irow)=spm%ncol(irow)+1
    next=>spm_list_next(next)
    jj=jj-1
  enddo
  !fill the last row
  spm%pB(irow)=jj+1
  call spm_list_free(spm_list)
  end subroutine spm_list_to_mat



  !-------------------------------------------------------------------------------
  subroutine calc_state_pscal(this, nprimvars, dtime, ystate, p_dt,  d_dt, pscal, lneg, errinfo)
    !
    ! !DESCRIPTION:
    ! calculate limiting factor from each primary state variable
    !
    use BetrstatusType     , only : betr_status_type
    use betr_constants     , only : betr_errmsg_len
    implicit none
    ! !ARGUMENTS:
    class(lom_type), intent(in) :: this
    integer,  intent(in)  :: nprimvars
    real(r8), intent(in)  :: dtime
    real(r8), intent(in)  :: ystate(1:nprimvars)
    real(r8), intent(in)  :: p_dt(1:nprimvars)
    real(r8), intent(in)  :: d_dt(1:nprimvars)
    real(r8), intent(out) :: pscal(1:nprimvars)
    logical,  intent(out) :: lneg
    integer,  intent(out) :: errinfo
    character(len=betr_errmsg_len) :: msg

    ! !LOCAL VARIABLES:
    real(r8) :: yt
    integer  :: j
    real(r8),parameter :: p_par=0.999_r8
    real(r8) :: tmp

    lneg =.false.
    errinfo=0
    do j = 1, nprimvars
       yt = ystate(j) + (p_dt(j)+d_dt(j))*dtime
       if(yt<tiny_val .and. d_dt(j)<0._r8)then
          tmp = dtime*d_dt(j)
          pscal(j) = -(max(p_dt(j),0._r8)*dtime+max(ystate(j),0._r8))/tmp*p_par
          if(abs(pscal(j))<1.e-14)pscal(j)=0._r8
          lneg=.true.
          if(pscal(j)<0._r8)then
             write(iulog,*)'pscal',pscal(j),p_dt(j),ystate(j),d_dt(j)
             errinfo=-1
             return
          endif
       else
          pscal(j) = 1._r8
       endif
    enddo
  end subroutine calc_state_pscal

  !-------------------------------------------------------------------------------
  subroutine calc_reaction_rscal(this, nprimvars, nr, pscal, spm_d, rscal)
    !
    ! !DESCRIPTION:
    ! calculate limiting factor for each reaction
    ! !USES:
    use BetrstatusType     , only : betr_status_type
    implicit none
    ! !ARGUMENTS:
    class(lom_type), intent(in) :: this
    integer , intent(in) :: nprimvars
    integer , intent(in) :: nr
    real(r8), intent(in) :: pscal(1:nprimvars)
    type(sparseMat_type), intent(in) :: spm_d
    real(r8), intent(out):: rscal(1:nr)
    ! !LOCAL VARIABLES:
    integer :: jj, ii, id

    rscal(:)=1._r8
    do ii=1, spm_d%szrow
      do jj = 1, spm_d%ncol(ii)
        id = spm_d%pB(ii)+jj-1
        rscal(spm_d%icol(id)) = min(pscal(ii),rscal(spm_d%icol(id)))
      enddo
    enddo

  end subroutine calc_reaction_rscal

  !-------------------------------------------------------------------------------
  subroutine apply_reaction_rscal(this, nr, rscal, reaction_rates)
    !
    ! !DESCRIPTION:
    ! reduce reaction rates using input scalar
    !
    implicit none
    ! !ARGUMENTS:
    class(lom_type), intent(in) :: this
    integer , intent(in)    :: nr
    real(r8), intent(in)    :: rscal(1:nr)
    real(r8), intent(inout) :: reaction_rates(1:nr)
    ! !LOCAL VARIABLES:
    integer :: j

    do j = 1, nr
       reaction_rates(j) = reaction_rates(j)*rscal(j)
    enddo
  end subroutine  apply_reaction_rscal

  !-------------------------------------------------------------------------------
  subroutine flux_correction(nvars, nreactions, spm_p, spm_d, dtime, ystates, rfluxes)
  !
  ! DESCRIPTION
  ! correcting the fluxes to avoid negative state variables
  ! Reference: Tang and Riley, Biogeosciences, 2016, doi: 10.5194/bg-13-723-2016
  use abortutils             , only : endrun
  implicit none
  integer, intent(in) :: nvars
  integer, intent(in) :: nreactions
  class(sparseMat_type), intent(in) :: spm_p
  class(sparseMat_type), intent(in) :: spm_d
  real(r8), intent(in) :: dtime
  real(r8), intent(in) :: ystates(nvars)
  real(r8), intent(inout):: rfluxes(nreactions)

  integer :: it
  real(r8) :: rscal(nreactions)
  real(r8) :: pscal(nvars)
  real(r8) :: d_dt(nvars)
  real(r8) :: p_dt(nvars)
  integer :: errinfo
  type(lom_type) :: lom
  logical :: lneg
  integer, parameter :: itmax=40

  !initialize the iterator and adjustment scalar vector for the destruction fluxes 
  it=0
  rscal=0._r8
  do
    !obtain the destruction fluxes
    d_dt(:)=0._r8
    call spm_axpy(spm_d%szrow, spm_d%szcol, 1._r8, rfluxes, spm_d, d_dt, errinfo)

    if(errinfo<0)then
      call endrun(msg='ERROR:: in spm_axpy'//errMsg(__FILE__, __LINE__))
    endif
    !obtain the production fluxes
    p_dt(:)=0._r8
    call spm_axpy(spm_p%szrow, spm_d%szcol, 1._r8, rfluxes, spm_p, p_dt, errinfo)

    if(errinfo<0)then
      call endrun(msg='ERROR:: in spm_axpy '//errMsg(__FILE__, __LINE__))
    endif
    !obtain flux reducing scalar for all primary variables
    call lom%calc_state_pscal(spm_d%szrow, dtime, ystates, p_dt,  d_dt, &
      pscal, lneg, errinfo)
    if(errinfo<0)then
      call endrun(msg='ERROR:: in calc_state_pscal '//errMsg(__FILE__, __LINE__))
    endif

    !obtain the adjustment scalar for each destruction flux based on the stoichiometry matrix
    if(lneg .and. it<=itmax)then
      call lom%calc_reaction_rscal(spm_d%szrow, spm_d%szcol,  pscal, &
        spm_d,rscal)

      call lom%apply_reaction_rscal(spm_d%szcol, rscal, rfluxes)
    else
      !terminate the loop and write a warning if maximum iteration is reached.
      if(it==itmax)write(iulog,*)'maximum iterations reached in flux_correction in LSparseMatMod'
      exit
    endif
    it = it + 1
  enddo
  end subroutine flux_correction


end module LSparseMatMod
