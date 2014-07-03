! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!!  This routine drives the fast microphysics calculations.
!!
!! @author Eric Jensen, Bill McKie
!! @version Sep-1997
subroutine microfast(carma, cstate, iz, scale_threshold, rc)

  ! types
  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carmastate_mod
  use carma_mod
  
  implicit none

  type(carma_type), intent(in)         :: carma       !! the carma object
  type(carmastate_type), intent(inout) :: cstate      !! the carma state object
  integer, intent(in)                  :: iz          !! z index
  real(kind=f)                         :: scale_threshold !! Scaling factor for convergence thresholds
  integer, intent(inout)               :: rc          !! return code, negative indicates failure

  ! Local Variables
  integer                              :: ielem   ! element index
  integer                              :: ibin    ! bin index
  integer                              :: igas    ! gas index
  real(kind=f)                         :: previous_ice(NGAS)      ! total ice at the start of substep
  real(kind=f)                         :: previous_liquid(NGAS)   ! total liquid at the start of substep
  real(kind=f)                         :: previous_supsatl(NGAS)  ! supersaturation wrt ice at the start of substep
  real(kind=f)                         :: previous_supsati(NGAS)  ! supersaturation wrt liquid at the start of substep
  real(kind=f)                         :: supsatold
  real(kind=f)                         :: supsatnew
  real(kind=f)                         :: srat
  real(kind=f)                         :: srat1
  real(kind=f)                         :: srat2
  real(kind=f)                         :: s_threshold

  1 format(/,'microfast::ERROR - excessive change in supersaturation for ',a,' : iz=',i4,',lat=', &
              f7.2,',lon=',f7.2,',srat=',e10.3,',supsatiold=',e10.3,',supsatlold=',e10.3,',supsati=',e10.3, &
              ',supsatl=',e10.3,',t=',f6.2)
  2 format('microfast::ERROR - conditions at beginning of the step : gc=',e10.3,',supsati=',e17.10, &
              ',supsatl=',e17.10,',t=',f6.2,',d_gc=',e10.3,',d_t=',f6.2)
  3 format(/,'microfast::ERROR - excessive change in supersaturation for ',a,' : iz=',i4,',lat=', &
              f7.2,',lon=',f7.2,',supsatiold=',e10.3,',supsatlold=',e10.3,',supsati=',e10.3, &
              ',supsatl=',e10.3,',t=',f6.2)
  
   ! Set production and loss rates to zero.
  call zeromicro(carma, cstate, iz, rc)
  if (rc < RC_OK) return
  

  ! Calculate (implicit) particle loss rates for nucleation, growth,
  ! evaporation, melting, etc.
  if (do_grow) then

    ! Save off the current condensate totals so the gas and latent heating can be
    ! figured out in a way that conserves mass and energy.
    call totalcondensate(carma, cstate, iz, previous_ice, previous_liquid, rc)
    if (rc < RC_OK) return
    
    do igas = 1, NGAS
      call supersat(carma, cstate, iz, igas, rc)
      if (rc < RC_OK) return
    
      previous_supsati(igas) = supsati(iz, igas)
      previous_supsatl(igas) = supsatl(iz, igas)     
    end do
    
    ! Have water vapor and sulfuric acid been defined?
    if ((igash2o /= 0) .and. (igash2so4 /= 0)) then
      
      ! Are both gases avaialble?
      if ((gc(iz, igash2o) > 0._f) .and. (gc(iz,igash2so4) > 0._f)) then

        ! See if any sulfates will form.
        call sulfnuc(carma, cstate, iz, rc) 
      endif 
    end if
    
    call growevapl(carma, cstate, iz, rc)
    if (rc < RC_OK) return

    call actdropl(carma, cstate, iz, rc)
    if (rc < RC_OK) return

    ! The Koop, Tabazadeh and Mohler routines provide different schemes for aerosol freezing.
    ! Only one of these parameterizatons should be active at one time.  However, any
    ! of these routines can be used in conjunction with heterogenous nucleation of glassy
    ! aerosols.
    call freezaerl_tabazadeh2000(carma, cstate, iz, rc)
    if (rc < RC_OK) return

    call freezaerl_koop2000(carma, cstate, iz, rc)
    if (rc < RC_OK) return

    call freezaerl_mohler2010(carma, cstate, iz, rc)
    if (rc < RC_OK) return

    call freezglaerl_murray2010(carma, cstate, iz, rc)
    if (rc < RC_OK) return

    call hetnucl(carma, cstate, iz, rc)
    if (rc < RC_OK) return

    call freezdropl(carma, cstate, iz, rc)
    if (rc < RC_OK) return

    call melticel(carma, cstate, iz, rc)
    if (rc < RC_OK) return    
  endif

  ! Calculate particle production terms and solve for particle 
  ! concentrations at end of time step.
  do ielem = 1,NELEM
    do ibin = 1,NBIN

      if( do_grow )then
        call growp(carma, cstate, iz, ibin, ielem, rc)
        if (rc < RC_OK) return

        call upgxfer(carma, cstate, iz, ibin, ielem, rc)
        if (rc < RC_OK) return
      endif

      call psolve(carma, cstate, iz, ibin, ielem, rc)
     if (rc < RC_OK) return
    enddo
  enddo

  ! Calculate particle production terms for evaporation;
  ! gas loss rates and production terms due to particle nucleation;
  ! growth, and evaporation;
  ! apply evaporation production terms to particle concentrations;
  ! and solve for gas concentrations at end of time step.
  if (do_grow) then
    call evapp(carma, cstate, iz, rc)
    if (rc < RC_OK) return

    call downgxfer(carma, cstate, iz, rc)
    if (rc < RC_OK) return

! NOTE: Not needed because changes in gas concentrations and latent
! heats are now calculated later in gsolve using total condensate.  
!    call gasexchange(carma, cstate, iz, rc)
!    if (rc < RC_OK) return

    call downgevapply(carma, cstate, iz, rc)
    if (rc < RC_OK) return

    call gsolve(carma, cstate, iz, previous_ice, previous_liquid, scale_threshold, rc)
    if (rc /=RC_OK) return
  endif

  ! Update temperature if thermal processes requested
  if (do_thermo) then
    call tsolve(carma, cstate, iz, scale_threshold, rc)
    if (rc /= RC_OK) return
  endif

  !  Update saturation ratios
  if (do_grow .or. do_thermo) then
    do igas = 1, NGAS
      call supersat(carma, cstate, iz, igas, rc)
      if (rc < RC_OK) return
    
      ! Check to see how much the supersaturation changed during this step. If it
      ! has changed to much, then cause a retry.
      if (t(iz) >= T0) then
        supsatold = previous_supsatl(igas)
        supsatnew = supsatl(iz,igas)
      else
        supsatold = previous_supsati(igas)
        supsatnew = supsati(iz,igas)
      end if

      ! If ds_threshold is positive, then it indicates that the criteria should
      ! be based on the percentage change in saturation.
      if (ds_threshold(igas) > 0._f) then
      
        if (supsatold >= 1.e-4_f) then
          srat1 = abs(supsatnew / supsatold - 1._f)
        else
          srat1 = 0._f
        end if
    
        if (supsatnew >= 1.e-4_f) then
          srat2 = abs(supsatold / supsatnew - 1._f)
        else
          srat2 = 0._f
        end if
    
        srat = max(srat1, srat2)
  
        ! Don't let one substep change the supersaturation by too much.
        if (ds_threshold(igas) > 0._f) then
!          if (srat >= ds_threshold(igas)) then
          if ((srat >= ds_threshold(igas)) .and. (abs(supsatold - supsatnew) > 0.1_f)) then
            if (do_substep) then
              if (nretries == maxretries) then 
                if (do_print) write(LUNOPRT,1) trim(gasname(igas)), iz, lat, lon, srat, previous_supsati(igas), previous_supsatl(igas), &
                supsati(iz, igas), supsatl(iz,igas), t(iz)       
                if (do_print) write(LUNOPRT,2) gcl(iz,igas), supsatiold(iz, igas), supsatlold(iz,igas), told(iz), d_gc(iz, igas), d_t(iz)
              end if
              
              rc = RC_WARNING_RETRY
            else
              if (do_print) write(LUNOPRT,1) trim(gasname(igas)), iz, lat, lon, srat, previous_supsati(igas), previous_supsatl(igas), &
              supsati(iz, igas), supsatl(iz,igas), t(iz)       
            end if
          end if
        end if
        
      
      ! If ds_threshold is negative, then it indicates that the criteria is based
      ! upon the supersaturation crossing 0, Indicating a shift from growth to
      ! evaporation and a potential overshoot in the result.
      else if (ds_threshold(igas) < 0._f) then

        ! Adjust the saturation threshold to allow a worse solution if getting a better
        ! solution is taking too much time. The particular solution at any individual
        ! point is probably not going to affect the overall result by too much.
        s_threshold = abs(ds_threshold(igas))
        
        if (nretries >= (0.8_f * maxretries)) then
          s_threshold = 4._f  * s_threshold
        else if (nretries >= (0.7_f * maxretries)) then
          s_threshold = 3.5_f * s_threshold
        else if (nretries >= (0.6_f * maxretries)) then
          s_threshold = 3._f  * s_threshold
        else if (nretries >= (0.5_f * maxretries)) then
          s_threshold = 2.5_f * s_threshold
        else if (nretries >= (0.4_f * maxretries)) then
          s_threshold = 2._f  * s_threshold
        end if
        
        ! If the supersaturation changed signs, then we went from growth to evaporation
        ! or vice versa. Don't let the new supersaturation go too far past 0 in one substep.
        ! This is to prevent overshooting as growth/evaporation should normally stop when
        ! the supersaturation is 0.
        if (((supsatnew * supsatold) < 0._f) .and. (abs(supsatnew) > s_threshold)) then

          if (do_substep) then
            if (nretries == maxretries) then 
              if (do_print) write(LUNOPRT,3) trim(gasname(igas)), iz, lat, lon, previous_supsati(igas), previous_supsatl(igas), &
              supsati(iz, igas), supsatl(iz,igas), t(iz)
              if (do_print) write(LUNOPRT,2) gcl(iz,igas), supsatiold(iz, igas), supsatlold(iz,igas), told(iz), d_gc(iz, igas), d_t(iz)
            end if
          else
            if (do_print) write(LUNOPRT,3) trim(gasname(igas)), iz, lat, lon, previous_supsati(igas), previous_supsatl(igas), &
            supsati(iz, igas), supsatl(iz,igas), t(iz)
          end if
  
          rc = RC_WARNING_RETRY
        end if
      end if
    end do
  endif


  ! Update particle densities
!  if (do_grow) then
!    call rhopart(carma, cstate, iz, rc)
!  end if
  
  ! Return to caller with new particle and gas concentrations.
  return
end
