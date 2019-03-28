#define DEBUG
Module glm2iac_mod
  
!---------------------------------------------------------------------------
!BOP
!
! !MODULE: glm2iac_mod
!
!  Interface of the integrated assessment component in CCSM
!
! !DESCRIPTION:
!
! !USES:

  use iac_fields_mod
  use shr_cal_mod
  use netcdf

  implicit none
  SAVE
  private                              ! By default make data private

! !PUBLIC MEMBER FUNCTIONS:

  public :: glm2iac_init_mod               ! clm initialization
  public :: glm2iac_run_mod                ! clm run phase
  public :: glm2iac_final_mod              ! clm finalization/cleanup

! !PUBLIC DATA MEMBERS: None


! !REVISION HISTORY:
! Author: T Craig

  real*8, pointer, save :: plodata(:,:)
  real, pointer, save :: plodataf(:,:)
  real, pointer, save :: glmof(:,:)

! !PRIVATE DATA MEMBERS:

!EOP
!===============================================================
contains
!===============================================================

!---------------------------------------------------------------------------
!BOP

! !IROUTINE: glm2iac_init_mod

! !INTERFACE:
  subroutine glm2iac_init_mod( EClock, cdata, glmo, iaci)

! !DESCRIPTION:
! Initialize interface for glm

! !USES:
    implicit none

! !ARGUMENTS:
    integer, pointer :: EClock(:)
    type(iac_cdata_type) :: cdata
    real*8, pointer :: glmo(:,:)
    real*8, pointer :: iaci(:,:)

! !LOCAL VARIABLES:

    integer :: iu
    integer :: nflds, nsize
    character(len=*),parameter :: subname='(glm2iac_init_mod)'

! !REVISION HISTORY:
! Author: T Craig

!EOP
!-----------------------------------------------------------------------

    iu  = cdata%i(iac_cdatai_logunit)

    nflds = iac_iac_npfts + 7
    nsize = cdata%i(iac_cdatai_glm_size)
    ! npfts + extra pft + vh1,vh2,sh1,sh2,sh3,grazing
    allocate(plodata(nflds,nsize))
    allocate(plodataf(nflds,nsize))
    allocate(glmof(size(glmo,dim=1),size(glmo,dim=2)))
    plodata = 0.0

#ifdef DEBUG
    write(iu,*) trim(subname),' allocate plodata ',nflds,nsize
#endif
  end subroutine glm2iac_init_mod

!---------------------------------------------------------------------------
!BOP

! !IROUTINE: glm2iac_run_mod

! !INTERFACE:
  subroutine glm2iac_run_mod( EClock, cdata, glmo, iaci)

! !DESCRIPTION:
! Run interface for glm

! !USES:
    use mksurfdat, only: mksurfdat_run
    implicit none

! !ARGUMENTS:
    integer, pointer :: EClock(:)
    type(iac_cdata_type) :: cdata
    real*8, pointer :: glmo(:,:)
    real*8, pointer :: iaci(:,:)

! !LOCAL VARIABLES:
    logical :: restart_now
    integer :: ymd, tod, dt
    integer :: iu
    integer :: i,j,n,ij,ierr,nmode
    integer :: dimid(3),varid,ncid
    real*8, pointer :: array3d(:,:,:)
    character(len=128) :: fname,casename,hfile
    integer :: myear, mon, day
    character(len=*),parameter :: subname='(glm2iac_run_mod)'

! !REVISION HISTORY:
! Author: T Craig

!EOP
!-----------------------------------------------------------------------

    iu  = cdata%i(iac_cdatai_logunit)

    ymd = EClock(iac_EClock_ymd)
    tod = EClock(iac_EClock_tod)
    dt  = EClock(iac_EClock_dt)

#ifdef DEBUG
    do j = 1,iac_glmo_nflds
       write(iu,*) trim(subname),' glmo minmax ',j,minval(glmo(j,:)),maxval(glmo(j,:))
    enddo
#endif

!--- temporary implementation

!    fname = 'mksurf_landuse_iESM_720x360.nc'
!    myear = ymd
!    write(iu,*) trim(subname),' fname = ',trim(fname)
!    write(iu,*) trim(subname),' myear = ',myear

!--- proper implementation

    call shr_cal_date2ymd(ymd,myear,mon,day)

    casename = trim(cdata%c(iac_cdatac_casename))
    write(hfile,'(a,i4.4,a,i2.2,a,i2.2,a)') trim(casename)//'.iac.hglmo.',myear,'-',mon,'-',day,'.nc'
#ifdef DEBUG
    write(iu,*) trim(subname),' writing history file ',trim(hfile)
#endif
    nmode = ior(NF90_CLOBBER,NF90_64BIT_OFFSET)
    ierr = nf90_create(trim(hfile),nmode,ncid)
    ierr = nf90_def_dim(ncid,'glmo_nx' ,iac_glm_nx,dimid(1))
    ierr = nf90_def_dim(ncid,'glmo_ny' ,iac_glm_ny,dimid(2))
    ierr = nf90_def_dim(ncid,'glmo_nf' ,size(glmo,dim=1),dimid(3))
    ierr = nf90_def_var(ncid,'glmodata',NF90_DOUBLE,dimid,varid)
    ierr = nf90_enddef(ncid)
    allocate(array3d(iac_glm_nx,iac_glm_ny,size(glmo,dim=1)))
    do n = 1,size(glmo,dim=1)
    ij = 0
    do j = 1,iac_glm_ny
    do i = 1,iac_glm_nx
       ij = ij + 1
       array3d(i,j,n) = glmo(n,ij)
    enddo
    enddo
    enddo
    ierr = nf90_put_var(ncid,varid,array3d)
    deallocate(array3d)
    ierr = nf90_close(ncid)

#ifdef DEBUG
    write(iu,*) trim(subname),' date= ',ymd,tod
    write(iu,*) trim(subname),' myear = ',myear
#endif

    plodataf=plodata
    glmof=glmo
    call updateannuallanduse(glmof,plodataf,myear)
!jt    call updateannuallanduse(glmo,plodata,myear)
    plodata=plodataf
    glmo=glmof

#ifdef DEBUG
    do j = 1,size(plodata,dim=1)
       write(iu,*) trim(subname),' plodata minmax ',j,minval(plodata(j,:)),maxval(plodata(j,:))
    enddo
#endif

    casename = trim(cdata%c(iac_cdatac_casename))
    write(hfile,'(a,i4.4,a,i2.2,a,i2.2,a)') trim(casename)//'.iac.hplo.',myear,'-',mon,'-',day,'.nc'
#ifdef DEBUG
    write(iu,*) trim(subname),' writing history file ',trim(hfile)
#endif
    nmode = ior(NF90_CLOBBER,NF90_64BIT_OFFSET)
    ierr = nf90_create(trim(hfile),nmode,ncid)
    ierr = nf90_def_dim(ncid,'plo_nx' ,iac_glm_nx,dimid(1))
    ierr = nf90_def_dim(ncid,'plo_ny' ,iac_glm_ny,dimid(2))
    ierr = nf90_def_dim(ncid,'plo_nf' ,size(plodata,dim=1),dimid(3))
    ierr = nf90_def_var(ncid,'plodata',NF90_DOUBLE,dimid,varid)
    ierr = nf90_enddef(ncid)
    allocate(array3d(iac_glm_nx,iac_glm_ny,size(plodata,dim=1)))
    do n = 1,size(plodata,dim=1)
    ij = 0
    do j = 1,iac_glm_ny
    do i = 1,iac_glm_nx
       ij = ij + 1
       array3d(i,j,n) = plodata(n,ij)
    enddo
    enddo
    enddo
    ierr = nf90_put_var(ncid,varid,array3d)
    deallocate(array3d)
    ierr = nf90_close(ncid)

    call mksurfdat_run(myear,plodata)
!    call mksurfdata(plodata)

  end subroutine glm2iac_run_mod


!---------------------------------------------------------------------------
!BOP

! !IROUTINE: glm2iac_final_mod

! !INTERFACE:
  subroutine glm2iac_final_mod( )

! !DESCRIPTION:
! Finalize glm model
! !USES:
    implicit none

! !ARGUMENTS:

! !LOCAL VARIABLES:
    integer :: iu
    character(len=*),parameter :: subname='(glm2iac_final_mod)'

! !REVISION HISTORY:
! Author: T Craig

!EOP

!---------------------------------------------------------------------------

!    iu  = cdata%i(iac_cdatai_logunit)
!    write(iu,*) trim(subname)
    deallocate(plodata)
    deallocate(plodataf)
    deallocate(glmof)

  end subroutine glm2iac_final_mod

!====================================================================================

end module glm2iac_mod

