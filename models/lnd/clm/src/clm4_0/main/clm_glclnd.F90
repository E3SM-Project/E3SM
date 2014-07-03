module clm_glclnd

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: clm_glclnd
!
! !DESCRIPTION:
! Handle arrays used for exchanging data between glc and land model.
! Based on clm_atmlnd (but without mapping routines because glc data
!  is send and received on the lnd decomposition, at least for now).
!
! The fields sent from the lnd component to the glc component via
!  the coupler are labeled 's2x', or sno to coupler.
! The fields received by the lnd component from the glc component
!  via the coupler are labeled 'x2s', or coupler to sno.
! 'Sno' is a misnomer in that the exchanged data are related to
!  the ice beneath the snow, not the snow itself.  But by CESM convention,
! 'ice' refers to sea ice, not land ice.
!
! !USES:
  use decompMod   , only : get_proc_bounds
  use shr_kind_mod, only : r8 => shr_kind_r8
  use nanMod      , only : nan
  use spmdMod     , only : masterproc
  use clm_varpar  , only : maxpatch_glcmec
  use clm_varctl  , only : iulog, glc_smb
  use abortutils  , only : endrun
!
! !REVISION HISTORY:
! Created by William Lipscomb, Dec. 2007, based on clm_atmlnd.F90.
!
! !PUBLIC TYPES:
  implicit none

!----------------------------------------------------
! glc -> land variables structure
!----------------------------------------------------
  type glc2lnd_type
     real(r8), pointer :: frac(:,:) 
     real(r8), pointer :: topo(:,:) 
     real(r8), pointer :: rofi(:,:) 
     real(r8), pointer :: rofl(:,:) 
     real(r8), pointer :: hflx(:,:) 
  end type glc2lnd_type

!----------------------------------------------------
! land -> glc variables structure
!----------------------------------------------------
  type lnd2glc_type
     real(r8), pointer :: tsrf(:,:) 
     real(r8), pointer :: topo(:,:)
     real(r8), pointer :: qice(:,:)
  end type lnd2glc_type

  type (lnd2glc_type), public, target :: clm_s2x  ! s2x fields on clm grid
  type (glc2lnd_type), public, target :: clm_x2s  ! x2s fields on clm grid

! !PUBLIC MEMBER FUNCTIONS:
  public :: init_glc2lnd_type
  public :: init_lnd2glc_type
!
! !PRIVATE MEMBER FUNCTIONS:

!EOP
!----------------------------------------------------

contains


!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_glc2lnd_type
!
! !INTERFACE:
  subroutine init_glc2lnd_type(beg, end, x2s)
!
! !DESCRIPTION:
! Initialize glc variables required by the land
!
! !ARGUMENTS:
  implicit none
  integer, intent(in) :: beg, end
  type (glc2lnd_type), intent(inout):: x2s
!
! !REVISION HISTORY:
! Created by William Lipscomb, based on init_atm2lnd_type
!EOP
!
! !LOCAL VARIABLES:
  real(r8) :: ival   ! initial value
!------------------------------------------------------------------------

  allocate(x2s%frac(beg:end,maxpatch_glcmec))
  allocate(x2s%topo(beg:end,maxpatch_glcmec))
  allocate(x2s%rofi(beg:end,maxpatch_glcmec))
  allocate(x2s%rofl(beg:end,maxpatch_glcmec))
  allocate(x2s%hflx(beg:end,maxpatch_glcmec))

! ival = nan      ! causes core dump in map_maparray, tcx fix
  ival = 0.0_r8

  x2s%frac(beg:end,:) = ival
  x2s%topo(beg:end,:) = ival
  x2s%rofi(beg:end,:) = ival
  x2s%rofl(beg:end,:) = ival
  x2s%hflx(beg:end,:) = ival

end subroutine init_glc2lnd_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_lnd2glc_type
!
! !INTERFACE:
  subroutine init_lnd2glc_type(beg, end, s2x)
!
! !DESCRIPTION:
! Initialize land variables required by glc
!
! !ARGUMENTS:
  implicit none
  integer, intent(in) :: beg, end
  type (lnd2glc_type), intent(inout):: s2x
!
! !REVISION HISTORY:
! Created by William Lipscomb, based on init_lnd2atm_type
!
!EOP
!
! !LOCAL VARIABLES:
  real(r8) :: ival   ! initial value
!------------------------------------------------------------------------

  allocate(s2x%tsrf(beg:end,maxpatch_glcmec))
  allocate(s2x%topo(beg:end,maxpatch_glcmec))
  allocate(s2x%qice(beg:end,maxpatch_glcmec))

! ival = nan      ! causes core dump in map_maparray, tcx fix
  ival = 0.0_r8

  s2x%tsrf(beg:end,:) = ival
  s2x%topo(beg:end,:) = ival
  s2x%qice(beg:end,:) = ival

end subroutine init_lnd2glc_type

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: create_clm_s2x
!
! !INTERFACE:
  subroutine create_clm_s2x(init)
!
! !DESCRIPTION:
! Assign values to clm_s2x based on the appropriate derived types
!
! !USES:
  use clmtype     
  use domainMod   , only : ldomain
  use clm_varcon  , only : istice_mec
  use clm_atmlnd  , only : clm_l2a, clm_a2l
  use clm_varcon  , only : spval
!
! !ARGUMENTS:
  implicit none

  logical, intent(in)   :: init    ! if true=>only set a subset of fields
!
! !REVISION HISTORY:
! Written by William Lipscomb, Feb. 2009 
!

    integer :: begg, endg              ! per-proc beginning and ending gridcell indices
    integer :: begc, endc              ! per-proc beginning and ending column indices
    integer :: c, l, g, n              ! indices
    integer , pointer :: ityplun(:)    ! landunit type
    integer , pointer :: clandunit(:)  ! column's landunit index
    integer , pointer :: cgridcell(:)  ! column's gridcell index

    ! Assign local pointers to derived type members

    clandunit => col%landunit
    cgridcell => col%gridcell
    ityplun   => lun%itype

    ! Get processor bounds
    call get_proc_bounds(begg, endg, begc=begc, endc=endc)

    ! Initialize qice because otherwise it will remain unset if init=true and
    ! glc_smb=true; note that the value here is the value qice will keep if these
    ! conditions hold

    clm_s2x%qice(:,:) = 0._r8

    ! and initialize the other variables just to be safe

    clm_s2x%tsrf(:,:) = 0._r8
    clm_s2x%topo(:,:) = 0._r8

    ! Fill the clm_s2x vector on the clm grid

    if (glc_smb) then   ! send surface mass balance info
       do c = begc, endc
          l = clandunit(c)
          g = cgridcell(c)

          ! Following assumes all elevation classes are populated
          if (ityplun(l) == istice_mec) then
             n = c - lun%coli(l) + 1    ! elevation class index

             ! t_soisno and glc_topo are valid even in initialization, so tsrf and topo
             ! are set here regardless of the value of init. But qflx_glcice is not valid
             ! until the run loop; thus, in initialization, we will use the default value
             ! for qice, as set above.
             clm_s2x%tsrf(g,n) = ces%t_soisno(c,1)
             clm_s2x%topo(g,n) = cps%glc_topo(c)
             if (.not. init) then
                clm_s2x%qice(g,n) = cwf%qflx_glcice(c)

                ! Check for bad values of qice
                if ( abs(clm_s2x%qice(g,n)) > 1.0_r8 .and. clm_s2x%qice(g,n) /= spval) then
                   write(iulog,*) 'WARNING: qice out of bounds: g, n, qice =', g, n, clm_s2x%qice(g,n)
                endif
             end if

           endif    ! istice_mec
       enddo        ! c
    else  ! Pass PDD info (same info in each elevation class)
       ! Require maxpatch_glcmec = 1 for this case 
       if (maxpatch_glcmec .ne. 1) then
          call endrun('create_clm_s2x error: maxpatch_glcmec must be 1 if glc_smb is false') 
       end if
       n = 1
       do g = begg, endg
          clm_s2x%tsrf(g,n) = clm_l2a%t_ref2m(g)
          clm_s2x%qice(g,n) = clm_a2l%forc_snow(g)   ! Assume rain runs off
          clm_s2x%topo(g,n) = ldomain%topo(g)
          ! Check for bad values of qice
          if (clm_s2x%qice(g,n) > -1.0_r8 .and. clm_s2x%qice(g,n) < 1.0_r8) then
             continue
          else
             write(iulog,*) 'WARNING: qice out of bounds: g, n, qice =', g, n, clm_s2x%qice(g,n)
             write(iulog,*) 'forc_rain =', clm_a2l%forc_rain(g)
             write(iulog,*) 'forc_snow =', clm_a2l%forc_snow(g)
          endif
       enddo
    endif   ! glc_smb

end subroutine create_clm_s2x

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: unpack_clm_x2s
!
! !INTERFACE:
  subroutine unpack_clm_x2s(clm_x2s)
!
! !DESCRIPTION:
! Unpack clm_x2s and update the appropriate derived types
!
! !USES:
  use clm_varcon  , only : istice_mec
  use clmtype     

!
! !ARGUMENTS:
  implicit none

  type(glc2lnd_type), intent(in) :: clm_x2s
!
! !REVISION HISTORY:
! Written by William Lipscomb, Feb. 2009 
!
    integer :: begc, endc              ! per-proc beginning and ending column indices
    integer :: c, l, g, n              ! indices
    integer , pointer :: ityplun(:)    ! landunit type
    integer , pointer :: clandunit(:)  ! column's landunit index
    integer , pointer :: cgridcell(:)  ! column's gridcell index

    logical :: update_glc2sno_fields   ! if true, update glacier_mec fields

    ! Assign local pointers to derived type members

    clandunit     => col%landunit
    cgridcell     => col%gridcell
    ityplun       => lun%itype

    update_glc2sno_fields = .false.
                                           
    if (update_glc2sno_fields) then 
   
       do c = begc, endc
          l = clandunit(c)
          g = cgridcell(c)

          if (ityplun(l) == istice_mec) then
             n = c - lun%coli(l) + 1    ! elevation class index
             cps%glc_frac(c) = clm_x2s%frac(g,n)
             cps%glc_topo(c) = clm_x2s%topo(g,n)
             cwf%glc_rofi(c) = clm_x2s%rofi(g,n)
             cwf%glc_rofl(c) = clm_x2s%rofl(g,n)
             cef%eflx_bot(c) = clm_x2s%hflx(g,n)

          endif
       enddo

    endif   ! update fields

    end subroutine unpack_clm_x2s

!------------------------------------------------------------------------

end module clm_glclnd

