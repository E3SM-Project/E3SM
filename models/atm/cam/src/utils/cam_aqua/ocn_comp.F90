module ocn_comp

! !USES:

  use shr_kind_mod,    only: r8 => shr_kind_r8
  use physconst,       only: tmelt, pi
  use ppgrid,          only: pcols, begchunk, endchunk
  use phys_grid,       only: get_ncols_p, get_rlat_all_p, get_rlon_all_p
  use abortutils,      only: endrun
  use ocn_types
!
! !PUBLIC TYPES:
  implicit none
  private                  ! By default everything private to this module

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public ocn_init          ! Initialization method
  public ocn_run           ! Run method
  save

  real(r8), allocatable :: sst(:,:)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: ocn_init
!
! !DESCRIPTION:
!     initialize docn model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

  subroutine ocn_init( ocn_out )

! !INPUT/OUTPUT PARAMETERS:

    type(ocn_out_t), pointer :: ocn_out(:)

    !--- local variables ---

    integer :: ncol,c,i

!-------------------------------------------------------------------------------

    allocate (ocn_out(begchunk:endchunk))
    do c = begchunk,endchunk
       ocn_out(c)%ts(:)  = 0.0_r8
    end do
    if (.not.allocated(sst)) then
       allocate(sst(pcols,begchunk:endchunk))
    endif

    call prescribed_sst()

    do c = begchunk,endchunk
       ncol = get_ncols_p(c)
       do i = 1, ncol
          ocn_out(c)%ts(i) = sst(i,c) + tmelt
       end do
    end do
  end subroutine ocn_init

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: ocn_run
!
! !DESCRIPTION:
!     run method for ocn model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

  subroutine ocn_run( ocn_out )

! !INPUT/OUTPUT PARAMETERS:
    type(ocn_out_t), intent(inout) :: ocn_out(begchunk:endchunk)
!EOP
   !--- local ---
    integer :: ncol,i,c

!-------------------------------------------------------------------------------

    do c=begchunk,endchunk
       ncol = get_ncols_p(c)
       do i = 1,ncol
          ocn_out(c)%ts(i)   = sst(i,c) + tmelt
       end do
    end do
  end subroutine ocn_run

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: prescribed_sst
!
! !DESCRIPTION:
!     run method for ocn model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

  subroutine prescribed_sst
    implicit none
!EOP
   !--- local ---
    integer,parameter :: sst_option = 1

    real(r8), parameter :: pio180     = pi/180._r8
    !
    ! Parameters for zonally symmetric experiments
    !
    real(r8), parameter ::   t0_max     = 27._r8
    real(r8), parameter ::   t0_min     = 0._r8
    real(r8), parameter ::   maxlat     = 60._r8*pio180
    real(r8), parameter ::   shift      = 5._r8*pio180
    real(r8), parameter ::   shift9     = 10._r8*pio180
    real(r8), parameter ::   shift10    = 15._r8*pio180
    !
    ! Parameters for zonally asymmetric experiments
    !
    real(r8), parameter ::   t0_max6    = 1._r8
    real(r8), parameter ::   t0_max7    = 3._r8
    real(r8), parameter ::   latcen     = 0._r8*pio180
    real(r8), parameter ::   loncen     = 0._r8*pio180
    real(r8), parameter ::   latrad6    = 15._r8*pio180
    real(r8), parameter ::   latrad8    = 30._r8*pio180
    real(r8), parameter ::   lonrad     = 30._r8*pio180

    integer :: lchnk, i, ncols
    real(r8) :: tmp, tmp1, rlat(pcols), rlon(pcols)
!-------------------------------------------------------------------------------
    !
    ! Control
    !
    if (sst_option .lt. 1 .or. sst_option .gt. 10) then
       call endrun ('prescribed_sst: sst_option must be between 1 and 10')
    endif
    if(sst_option == 1 .or. sst_option == 6 .or. &
         sst_option == 7 .or. sst_option == 8     ) then

       do lchnk=begchunk,endchunk
          call get_rlat_all_p(lchnk,pcols,rlat)
          ncols = get_ncols_p(lchnk)
          do i=1,ncols
             if(abs(rlat(i)) > maxlat) then
                sst(i,lchnk)=t0_min
             else
                tmp = sin(rlat(i)*pi*0.5_r8/maxlat)
                tmp = 1._r8 - tmp*tmp
                sst(i,lchnk)=tmp*(t0_max - t0_min) + t0_min
             end if
          end do
       end do
    end if
    !
    ! Flat
    !
    if(sst_option == 2) then
       do lchnk=begchunk,endchunk
          call get_rlat_all_p(lchnk,pcols,rlat)
          ncols = get_ncols_p(lchnk)
          do i=1,ncols
             if(abs(rlat(i)) > maxlat) then
                sst(i,lchnk)=t0_min
             else
                tmp = sin(rlat(i)*pi*0.5_r8/maxlat)
                tmp = 1._r8 - tmp*tmp*tmp*tmp
                sst(i,lchnk)=tmp*(t0_max - t0_min) + t0_min
             end if
          end do
       end do
    end if
    !
    ! Qobs
    !
    if(sst_option == 3) then
       do lchnk=begchunk,endchunk
          call get_rlat_all_p(lchnk,pcols,rlat)
          ncols = get_ncols_p(lchnk)
          do i=1,ncols
             if(abs(rlat(i)) > maxlat) then
                sst(i,lchnk)=t0_min
             else
                tmp = sin(rlat(i)*pi*0.5_r8/maxlat)
                tmp = (2._r8 - tmp*tmp*tmp*tmp - tmp*tmp)*0.5_r8
                sst(i,lchnk)=tmp*(t0_max - t0_min) + t0_min
             end if
          end do
       end do
    end if
    !
    ! Peaked
    !
    if(sst_option == 4) then
       do lchnk=begchunk,endchunk
          call get_rlat_all_p(lchnk,pcols,rlat)
          ncols = get_ncols_p(lchnk)
          do i=1,ncols
             if(abs(rlat(i)) > maxlat) then
                sst(i,lchnk)=t0_min
             else
                tmp = (maxlat - abs(rlat(i)))/maxlat
                tmp1 = 1._r8 - tmp
                sst(i,lchnk)= t0_max*tmp + t0_min*tmp1
             end if
          end do
       end do
    end if
    !
    ! Control-5N
    !
    if(sst_option == 5) then
       do lchnk=begchunk,endchunk
          call get_rlat_all_p(lchnk,pcols,rlat)
          ncols = get_ncols_p(lchnk)
          do i=1,ncols
             if(abs(rlat(i)) > maxlat) then
                sst(i,lchnk)=t0_min
             else if(rlat(i) > shift) then
                tmp = sin((rlat(i)-shift)*pi*0.5_r8/(maxlat-shift))
                tmp = 1._r8 - tmp*tmp
                sst(i,lchnk) = tmp*(t0_max - t0_min) + t0_min
             else
                tmp = sin((rlat(i)-shift)*pi*0.5_r8/(maxlat+shift))	
                tmp = 1._r8 - tmp*tmp
                sst(i,lchnk) = tmp*(t0_max - t0_min) + t0_min
             end if
          end do
       end do
    end if
    !
    ! 1KEQ
    !
    if(sst_option == 6) then
       do lchnk=begchunk,endchunk
          call get_rlat_all_p(lchnk,pcols,rlat)
          call get_rlon_all_p(lchnk,pcols,rlon)
          ncols = get_ncols_p(lchnk)
          do i=1,ncols
             if(abs(rlat(i)-latcen) <= latrad6) then
                tmp1 = cos((rlat(i)-latcen)*pi*0.5_r8/latrad6)
                tmp1 = tmp1*tmp1
                tmp = abs(rlon(i)-loncen)
                tmp = min(tmp , 2._r8*pi-tmp)
                if(tmp <= lonrad) then
                   tmp = cos(tmp*pi*0.5_r8/lonrad)
                   tmp = tmp*tmp
                   sst(i,lchnk) = sst(i,lchnk) + t0_max6*tmp*tmp1
                end if
             end if
          end do
       end do
    end if
    !
    ! 3KEQ
    !
    if(sst_option == 7) then
       do lchnk=begchunk,endchunk
          call get_rlat_all_p(lchnk,pcols,rlat)
          call get_rlon_all_p(lchnk,pcols,rlon)
          ncols = get_ncols_p(lchnk)
          do i=1,ncols
             if(abs(rlat(i)-latcen) <= latrad6) then
                tmp1 = cos((rlat(i)-latcen)*pi*0.5_r8/latrad6)
                tmp1 = tmp1*tmp1
                tmp = abs(rlon(i)-loncen)
                tmp = min(tmp , 2._r8*pi-tmp)
                if(tmp <= lonrad) then
                   tmp = cos(tmp*pi*0.5_r8/lonrad)
                   tmp = tmp*tmp
                   sst(i,lchnk) = sst(i,lchnk) + t0_max7*tmp*tmp1
                end if
             end if
          end do
       end do
    end if
    !
    ! 3KW1
    !
    if(sst_option == 8) then
       do lchnk=begchunk,endchunk
          call get_rlat_all_p(lchnk,pcols,rlat)
          call get_rlon_all_p(lchnk,pcols,rlon)
          ncols = get_ncols_p(lchnk)
          do i=1,ncols
             if(abs(rlat(i)-latcen) <= latrad8) then
                tmp1 = cos((rlat(i)-latcen)*pi*0.5_r8/latrad8)
                tmp1 = tmp1*tmp1
                tmp = cos(rlon(i)-loncen)
                sst(i,lchnk) = sst(i,lchnk) + t0_max7*tmp*tmp1
             end if
          end do
       end do
    end if
    !
    ! Control-10N
    !
    if(sst_option == 9) then
       do lchnk=begchunk,endchunk
          call get_rlat_all_p(lchnk,pcols,rlat)
          ncols = get_ncols_p(lchnk)
          do i=1,ncols
             if(abs(rlat(i)) > maxlat) then
                sst(i,lchnk) = t0_min
             elseif(rlat(i) > shift9) then
                tmp = sin((rlat(i)-shift9)*pi*0.5_r8/(maxlat-shift9))
                tmp = 1._r8 - tmp*tmp
                sst(i,lchnk) = tmp*(t0_max - t0_min) + t0_min
             else
                tmp = sin((rlat(i)-shift9)*pi*0.5_r8/(maxlat+shift9))
                tmp = 1._r8 - tmp*tmp
                sst(i,lchnk) = tmp*(t0_max - t0_min) + t0_min
             end if
          end do
       end do
    end if
    !
    ! Control-15N
    !
    if(sst_option == 10) then
       do lchnk=begchunk,endchunk
          call get_rlat_all_p(lchnk,pcols,rlat)
          ncols = get_ncols_p(lchnk)
          do i=1,ncols
             if(abs(rlat(i)) > maxlat) then
                sst(i,lchnk) = t0_min
             elseif(rlat(i) > shift10) then
                tmp = sin((rlat(i)-shift10)*pi*0.5_r8/(maxlat-shift10))
                tmp = 1._r8 - tmp*tmp
                sst(i,lchnk) = tmp*(t0_max - t0_min) + t0_min
             else
                tmp = sin((rlat(i)-shift10)*pi*0.5_r8/(maxlat+shift10))
                tmp = 1._r8 - tmp*tmp
                sst(i,lchnk) = tmp*(t0_max - t0_min) + t0_min
             end if
          end do
       end do
    endif
  end subroutine prescribed_sst

end module ocn_comp
