module dynCNDVMod

  !---------------------------------------------------------------------------
  !
  ! !DESCRIPTION:
  ! Handle weight updates associated with prognostic dynamic vegetation (CNDV)
  !
  ! !USES:
  use clmtype
  use decompMod, only : bounds_type
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  private
  save
  public :: dynCNDV_init      ! initialize CNDV weight updates
  public :: dynCNDV_interp    ! interpolate CNDV weight updates to the time step

contains

   !-----------------------------------------------------------------------
   subroutine dynCNDV_init(bounds)
     !
     ! !DESCRIPTION:
     ! Initialize time interpolation of cndv pft weights from annual to time step
     !
     ! Should be called once, in model initialization
     !
     ! !USES:
     use clm_varctl, only : nsrest, nsrStartup
     !
     ! !ARGUMENTS:
     type(bounds_type), intent(in) :: bounds  ! bounds
     !
     ! !LOCAL VARIABLES:
     integer  :: ier, p                        ! error status, do-loop index
     character(len=32) :: subname='dynCNDV_init' ! subroutine name
     !-----------------------------------------------------------------------

     if (nsrest == nsrStartup) then
        do p = bounds%begp,bounds%endp
           pdgvs%fpcgrid(p) = pft%wtcol(p)
           pdgvs%fpcgridold(p) = pft%wtcol(p)
        end do
     end if

  end subroutine dynCNDV_init

  !-----------------------------------------------------------------------
  subroutine dynCNDV_interp( bounds )
    !
    ! !DESCRIPTION:
    ! Time interpolate cndv pft weights from annual to time step
    !
    ! !USES:
    use clm_time_manager, only : get_curr_date, get_step_size, get_nstep, get_curr_yearfrac
    use clm_varcon      , only : istsoil ! CNDV incompatible with dynLU
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds  ! bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: c,g,l,p            ! indices
    real(r8) :: cday               ! current calendar day (1.0 = 0Z on Jan 1)
    real(r8) :: wt1                ! time interpolation weights (weight of time 1)
    real(r8) :: dtime              ! model time step
    real(r8) :: days_per_year      ! days per year
    integer  :: nstep              ! time step number
    integer  :: year               ! year (0, ...) at nstep + 1
    integer  :: mon                ! month (1, ..., 12) at nstep + 1
    integer  :: day                ! day of month (1, ..., 31) at nstep + 1
    integer  :: sec                ! seconds into current date at nstep + 1
    character(len=32) :: subname='dynCNDV_interp' ! subroutine name
    !-----------------------------------------------------------------------

    ! Interpolate pft weight to current time step
    ! Map interpolated pctpft to subgrid weights
    ! assumes maxpatch_pft = numpft + 1, each landunit has 1 column, 
    ! SCAM not defined and create_croplandunit = .false.

    nstep         = get_nstep()
    dtime         = get_step_size()

    wt1 = 1.0_r8 - get_curr_yearfrac(offset = -int(dtime))

    call get_curr_date(year, mon, day, sec, offset=int(dtime))

    do p = bounds%begp,bounds%endp
       g = pft%gridcell(p)
       l = pft%landunit(p)

       if (lun%itype(l) == istsoil .and. lun%wtgcell(l) > 0._r8) then ! CNDV incompatible with dynLU
          pft%wtcol(p)   = pdgvs%fpcgrid(p) + &
                     wt1 * (pdgvs%fpcgridold(p) - pdgvs%fpcgrid(p))

          if (mon==1 .and. day==1 .and. sec==dtime .and. nstep>0) then
             pdgvs%fpcgridold(p) = pdgvs%fpcgrid(p)
          end if
       end if
    end do

  end subroutine dynCNDV_interp

end module dynCNDVMod
