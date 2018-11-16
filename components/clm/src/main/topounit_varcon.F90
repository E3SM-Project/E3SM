module topounit_varcon

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing topounit indices and associated variables and routines.
  !
  ! !USES:
!#include "shr_assert.h"
  use ncdio_pio       , only : file_desc_t, var_desc_t, ncd_pio_openfile, ncd_pio_closefile
  use ncdio_pio       , only : ncd_io, check_var, ncd_inqfdims, check_dim, ncd_inqdid, ncd_inqdlen
  use clm_varctl      , only: fsurdat, iulog
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use clm_varcon      , only : grlnd
  !use abortutils      , only : endrun
  use pio
  use spmdMod
  !
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  
  !------------------------------------------------------------------
  ! Initialize topounit type constants
  !------------------------------------------------------------------
  
  !integer, parameter, public :: max_topounits  = 1 ! maximum number of topounits per gridcell
  integer, public :: max_topounits               ! maximum number of topounits per gridcell
  logical, public :: has_topounit                ! true => topounit dimension is on dataset
  !integer, pointer, public :: ntpu_per_grd(:)             ! Number of topounits per grid globally  num_tunits_per_grd
 ! integer, pointer, public :: tpu_lnd(:)             ! Number of topounit per grid for the land domain
  !integer, pointer, public :: tpu_glo_ind(:)         ! global index for number of topounits per grid
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: topounit_varcon_init  ! initialize constants in this module
  !-----------------------------------------------------------------------
  
  contains
  
  !-----------------------------------------------------------------------
  subroutine topounit_varcon_init(begg, endg,lfsurdat, ldomain)
    !
    ! !DESCRIPTION:
    ! Initialize topounit parameters
    !
    ! !USES:
    use fileutils   , only : getfil    
    use domainMod , only : domain_type
    !
    ! !ARGUMENTS:
    integer          ,intent(in)    :: begg, endg
    character(len=*), intent(in) :: lfsurdat    ! surface dataset filename
    type(domain_type),intent(in) :: ldomain     ! land domain
    
    !integer , intent(in) :: amask(:)  
   ! integer :: ni
   ! integer :: nj
    !
    ! !LOCAL VARIABLES:
    type(file_desc_t)     :: ncid         ! netcdf id
    integer :: dimid,varid
    integer :: topounits_size                      ! Size of topounit dimension in the surface data 
    character(len=256):: locfn                ! local file name
    logical :: isgrid2d    
    integer :: n,i,j               ! index 
    integer :: ns,ln,lns
    integer :: ni,nj
    integer :: ier                 ! error status
    integer :: numg
    integer :: t
    type(var_desc_t)   :: vardesc  ! variable descriptor
    character(len=256) :: varname  ! variable name
    logical :: readvar             ! read variable in or not
    integer , allocatable :: idata2d(:,:)
    character(len=*), parameter :: subname = 'topounit_varcon_init'
    !-----------------------------------------------------------------------
     
    if (masterproc) then
       write(iulog,*) 'Attempting to read topounit information from surface data .....'
       if (lfsurdat == ' ') then
          write(iulog,*)'lfsurdat must be specified'
          !call endrun(msg=errMsg(__FILE__, __LINE__))
       endif
    endif
    
    !! count total land gridcells
    !lns = ni*nj
    !numg = 0
    !!ns = len(amask(:))
    !do ln = 1, lns
    !   if (amask(ln) == 1) then
    !      numg = numg + 1
    !   endif
    !enddo

    ! Read surface data
    call getfil( lfsurdat, locfn, 0 )
    call ncd_pio_openfile (ncid, trim(locfn), 0)
    
     !Obtain the maximum number of topounits
    call ncd_inqdid(ncid, 'topounit', dimid, dimexist=has_topounit)
    if (.not. has_topounit) then
       max_topounits = 1
       if (masterproc) then          
          write(iulog,*)'Surface dataset has no topounit dimention; max_topounits is set to 1'
       end if
    else
       call ncd_inqdlen(ncid, dimid, topounits_size)  ! Get the dimension size of topounit from file
       max_topounits = topounits_size
    endif
    write(iulog,*)' TKT max_topounits ', max_topounits
    if(has_topounit .and. max_topounits > 1) then
       !write(iulog,*)' TKT Surface dataset checking max_topounits '
       call check_var(ncid=ncid, varname='topoPerGrid', vardesc=vardesc, readvar=readvar)
       if (readvar) then
          call ncd_io(ncid=ncid, varname= 'topoPerGrid', flag='read', data=ldomain%num_tunits_per_grd, &
             dim1name=grlnd, readvar=readvar)
        endif    
       !write(iulog,*)'TKT Surface dataset chcking max_topounit '
      ! ldomain%num_tunits_per_grd(:) = 0
       !call ncd_io(ncid=ncid, varname='topoPerGrid', flag='read', data=ldomain%num_tunits_per_grd, &
       !     dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          write(iulog,*)' TKT ERROR: Reading number of topounits per grid from lfsurfdat file '
        !   call endrun(msg=errMsg(__FILE__, __LINE__))
        end if
       
       !call ncd_io(ncid=ncid, varname='topoPerGrid', data=idata2d, flag='read', readvar=readvar)
       
       ! Make sure the glc mask is a subset of the land mask
       do n = begg,endg
          !write(iulog,*)'TKT num_tunits_per_grd(n) ', ldomain%num_tunits_per_grd(n)
          if (ldomain%num_tunits_per_grd(n)>1 .and. ldomain%mask(n)==0) then
             write(iulog,*)trim(subname),&
                  'initialize1: landmask/Number of topounits mismatch'
             write(iulog,*)trim(subname),&
                  'More than 1 topounits where landmask = 0, gridcell index', n
             !call endrun(msg=errMsg(__FILE__, __LINE__))
          endif
       enddo
        
       !! Determine dimensions and if grid file is 2d or 1d
       ! call ncd_inqfdims(ncid, isgrid2d, ni, nj, ns)        
       ! allocate(ntpu_per_grd(ns))
       ! 
       ! if (isgrid2d) then
       !    allocate(idata2d(ni,nj))
       !    idata2d(:,:) = 1	
       !    call ncd_io(ncid=ncid, varname='topoPerGrid', data=idata2d, flag='read', readvar=readvar)           
       !    if (readvar) then
       !       do j = 1,nj
       !       do i = 1,ni
       !          n = (j-1)*ni + i	
       !          ntpu_per_grd(n) = idata2d(i,j)
       !       enddo
       !       enddo
       !    end if
       !    deallocate(idata2d)
       ! else
       !    call ncd_io(ncid=ncid, varname='topoPerGrid', data=ntpu_per_grd, flag='read', readvar=readvar)           
       ! end if
       ! if (.not. readvar) then
       !    write(iulog,*)' ERROR: Reading number of topounits per grid from lfsurfdat file '
       ! !   call endrun(msg=errMsg(__FILE__, __LINE__))
       ! end if       
       
        ! count total land gridcells
       !!lns = ni*nj
       !numg = 0
       !!ns = len(amask(:))
       !do ln = 1, ns
       !   if (amask(ln) == 1) then
       !      numg = numg + 1
       !   endif
       !enddo
       !! Extract the number of topounits per grid for the land domain
       !if(has_topounit .and. max_topounits > 1) then
       !   allocate(tpu_glo_ind(numg))
       !   allocate(tpu_lnd(numg)) 
       !   t = 0
       !   do ln = 1,ns
       !      if (amask(ln) == 1) then
       !         t = t + 1
       !         tpu_lnd(t) = ntpu_per_grd(ln)
       !         tpu_glo_ind(t) = ln          
       !      endif
       !   enddo
       !endif    
    end if
    
    call ncd_pio_closefile(ncid)
  end subroutine topounit_varcon_init  
end module topounit_varcon
