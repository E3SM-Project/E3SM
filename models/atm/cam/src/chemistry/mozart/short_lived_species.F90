!---------------------------------------------------------------------
! Manages the storage of non-transported short-lived chemical species
! in the physics buffer.
!
! Created by: Francis Vitt -- 20 Aug 2008
!---------------------------------------------------------------------
module short_lived_species

  use shr_kind_mod, only : r8 => shr_kind_r8
  use chem_mods,    only : slvd_lst, nslvd, gas_pcnst
  use cam_logfile,  only : iulog
  use ppgrid,       only : pcols, pver, begchunk, endchunk
  use spmd_utils,   only : masterproc
  

  implicit none

  save
  private
  public :: map
  public :: register_short_lived_species
  public :: initialize_short_lived_species
  public :: set_short_lived_species
  public :: get_short_lived_species
  public :: slvd_index
  public :: pbf_idx

  integer :: pbf_idx
  integer :: map(nslvd)

  character(len=16), parameter :: pbufname = 'ShortLivedSpecies'

contains

!---------------------------------------------------------------------
!---------------------------------------------------------------------
  subroutine register_short_lived_species
    use physics_buffer, only : pbuf_add_field, dtype_r8

    implicit none

    integer :: m

    if ( nslvd < 1 ) return

    call pbuf_add_field(pbufname,'global',dtype_r8,(/pcols,pver,nslvd/),pbf_idx)

  end subroutine register_short_lived_species

!---------------------------------------------------------------------
!---------------------------------------------------------------------
  subroutine initialize_short_lived_species(ncid_ini, pbuf2d)
    use ioFileMod,      only : getfil
    use error_messages, only : handle_ncerr
    use dycore,         only : dycore_is
    use mo_tracname,    only : solsym
    use ncdio_atm,      only : infld
    use pio,            only : file_desc_t
    use physics_buffer, only : physics_buffer_desc, pbuf_set_field, pbuf_get_chunk, pbuf_get_field

    implicit none

    type(file_desc_t), intent(inout) :: ncid_ini
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)


    integer          :: m,n,lchnk
    character(len=8) :: fieldname
    character(len=4) :: dim1name
    logical          :: found
    real(r8),pointer :: tmpptr(:,:,:)   ! temporary pointer
    real(r8),pointer :: tmpptrout(:,:)  ! temporary pointer

    if ( nslvd < 1 ) return

    found = .false.

    if(dycore_is('se')) then  
       dim1name='ncol'
    else
       dim1name='lon'
    end if

    call pbuf_set_field(pbuf2d, pbf_idx, 0._r8)

    allocate(tmpptr(pcols,pver,begchunk:endchunk))

    do m=1,nslvd
       n = map(m)
       fieldname = solsym(n)
       call infld( fieldname,ncid_ini,dim1name, 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
                   tmpptr, found, grid_map='PHYS')
       call pbuf_set_field(pbuf2d, pbf_idx, tmpptr, start=(/1,1,m/),kount=(/pcols,pver,1/))
    enddo

    deallocate(tmpptr)

  end subroutine initialize_short_lived_species

!---------------------------------------------------------------------
!---------------------------------------------------------------------
  subroutine set_short_lived_species( q, lchnk, ncol, pbuf )

    use physics_buffer, only : physics_buffer_desc, pbuf_set_field

    implicit none 

    real(r8), intent(in)               :: q(pcols,pver,gas_pcnst)
    integer,  intent(in)               :: lchnk, ncol
    type(physics_buffer_desc), pointer :: pbuf(:)

    integer :: m,n

    if ( nslvd < 1 ) return

    do m=1,nslvd
       n = map(m)
       call pbuf_set_field(pbuf, pbf_idx, q(:,:,n), start=(/1,1,m/),kount=(/pcols,pver,1/))
    enddo

  endsubroutine set_short_lived_species

!---------------------------------------------------------------------
!---------------------------------------------------------------------
  subroutine get_short_lived_species( q, lchnk, ncol, pbuf )
    use physics_buffer, only : physics_buffer_desc, pbuf_get_field

    implicit none 

    real(r8), intent(inout)            :: q(pcols,pver,gas_pcnst)
    integer,  intent(in)               :: lchnk, ncol
    type(physics_buffer_desc), pointer :: pbuf(:)
    real(r8),pointer                   :: tmpptr(:,:)


    integer :: m,n 

    if ( nslvd < 1 ) return

    do m=1,nslvd
       n = map(m)
       call pbuf_get_field(pbuf, pbf_idx, tmpptr, start=(/1,1,m/), kount=(/ pcols,pver,1 /))
       q(:ncol,:,n) = tmpptr(:ncol,:)
    enddo

  endsubroutine get_short_lived_species

!---------------------------------------------------------------------
!---------------------------------------------------------------------
  function slvd_index( name )
    implicit none

    character(len=*) :: name
    integer :: slvd_index

    integer :: m

    slvd_index = -1

    if ( nslvd < 1 ) return

    do m=1,nslvd
       if ( name == slvd_lst(m) ) then
          slvd_index = m
          return 
       endif
    enddo

  endfunction slvd_index

end module short_lived_species
