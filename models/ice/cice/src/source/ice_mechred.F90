!=======================================================================
!BOP
!
! !MODULE: ice_mechred - driver for mechanical redestribution
!
! !DESCRIPTION:
!
! Driver for ice mechanical redistribution (ridging)
!
! See these references:
!
! Flato, G. M., and W. D. Hibler III, 1995: Ridging and strength
!  in modeling the thickness distribution of Arctic sea ice,
!  J. Geophys. Res., 100, 18,611-18,626.
!
! Hibler, W. D. III, 1980: Modeling a variable thickness sea ice
!  cover, Mon. Wea. Rev., 108, 1943-1973, 1980.
!
! Lipscomb, W. H., E. C. Hunke, W. Maslowski, and J. Jakacki, 2006: 
!  Ridging, strength, and stability in sea ice models, submitted 
!  to J. Geophys. Res. 
! 
! Rothrock, D. A., 1975: The energetics of the plastic deformation of
!  pack ice by ridging, J. Geophys. Res., 80, 4514-4519.
!
! Thorndike, A. S., D. A. Rothrock, G. A. Maykut, and R. Colony, 
!  1975: The thickness distribution of sea ice, J. Geophys. Res., 
!  80, 4501-4513. 
!
! !REVISION HISTORY:
!  SVN:$Id: ice_mechred.F90 53 2007-02-08 00:02:16Z dbailey $
!
! authors: William H. Lipscomb, LANL
!          Elizabeth C. Hunke, LANL
!
! 2003: Vectorized by Clifford Chen (Fujitsu) and William Lipscomb
! 2004: Block structure added by William Lipscomb
! 2006: New options for participation and redistribution (WHL)
! 2006: Streamlined for efficiency by Elizabeth Hunke
!       Converted to free source form (F90)
!
! !INTERFACE:
!
      module ice_mechred
!
! !USES:
!
      use ice_kinds_mod
      use ice_domain_size
      use ice_constants
      use ice_fileunits
      use ice_itd, only: hin_max, ilyr1, slyr1, column_sum, &
                         column_conservation_check, compute_tracers
!
!EOP
!
      implicit none
      save

!-----------------------------------------------------------------------
! Ridging parameters
!-----------------------------------------------------------------------

      integer (kind=int_kind) :: & ! defined in namelist 
         kstrength    , & ! 0 for simple Hibler (1979) formulation 
                          ! 1 for Rothrock (1975) pressure formulation 
         krdg_partic  , & ! 0 for Thorndike et al. (1975) formulation 
                          ! 1 for exponential participation function 
         krdg_redist      ! 0 for Hibler (1980) formulation 
                          ! 1 for exponential redistribution function 
 
      real (kind=dbl_kind), parameter :: & 
         Cf = 17._dbl_kind   , & ! ratio of ridging work to PE change in ridging 
         Cs = p25            , & ! fraction of shear energy contrbtng to ridging 
         Cp = p5*gravit*(rhow-rhoi)*rhoi/rhow, & ! proport const for PE 
         fsnowrdg = p5       , & ! snow fraction that survives in ridging 
         Gstar  = p15        , & ! max value of G(h) that participates 
                                 ! (krdg_partic = 0) 
         astar  = p05        , & ! e-folding scale for G(h) participation 
                                 ! (krdg_partic = 1) 
         maxraft= c1         , & ! max value of hrmin - hi = max thickness 
                                 ! of ice that rafts (m) 
         Hstar  = c25        , & ! determines mean thickness of ridged ice (m) 
                                 ! (krdg_redist = 0) 
                                 ! Flato & Hibler (1995) have Hstar = 100 
         mu_rdg = c4         , & ! gives e-folding scale of ridged ice (m^.5) 
                                 ! (krdg_redist = 1) 
         Pstar = 2.75e4_dbl_kind, & ! constant in Hibler strength formula 
                                 ! (kstrength = 0) 
         Cstar = c20             ! constant in Hibler strength formula 
                                 ! (kstrength = 0) 

      logical (kind=log_kind), parameter :: &
         l_conservation_check = .true.  ! if true, check conservation
!         l_conservation_check = .false.  ! if true, check conservation
                                        ! (useful for debugging)

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !ROUTINE: ridge_ice - driver for mechanical redistribution
!
! !DESCRIPTION:
!
! Compute changes in the ice thickness distribution due to divergence
! and shear.
! NOTE: This subroutine operates over a single block.
!
! !REVISION HISTORY:
!
! author: William H. Lipscomb, LANL
!
! !INTERFACE:
!
      subroutine ridge_ice (nx_block,    ny_block,   &
                            dt_dyn,      dt_thm,     &
                            ntrcr,       icells,     &
                            indxi,       indxj,      &
                            rdg_conv,    rdg_shear,  &
                            aicen,       trcrn,      &
                            vicen,       vsnon,      &
                            eicen,       esnon,      &
                            aice0,                   &
                            trcr_depend, l_stop,     &
                            istop,       jstop,      &
                            dardg1dt,    dardg2dt,   &
                            dvirdgdt,    opening,    &
                            fresh,       fhocn,      &
                            fsoot)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells            , & ! number of cells with ice present
         ntrcr                 ! number of tracers in use

      integer (kind=int_kind), dimension (nx_block*ny_block), &
         intent(in) :: &
         indxi, indxj     ! compressed indices for cells with ice

      real (kind=dbl_kind), intent(in) :: &
         dt_dyn  , & ! dynamic time step
         dt_thm      ! thermodynamic time step for diagnostics

      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(in) :: &
         rdg_conv, & ! normalized energy dissipation due to convergence (1/s)
         rdg_shear   ! normalized energy dissipation due to shear (1/s)
 
      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat), &
         intent(inout) :: &
         aicen , & ! concentration of ice
         vicen , & ! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)
 
      real (kind=dbl_kind), dimension (nx_block,ny_block,max_ntrcr,ncat), &
         intent(inout) :: & 
         trcrn     ! ice tracers 
 
      real (kind=dbl_kind), dimension (nx_block,ny_block,ntilyr), &
         intent(inout) :: & 
         eicen     ! energy of melting for each ice layer (J/m^2)
 
      real (kind=dbl_kind), dimension (nx_block,ny_block,ntslyr), &
         intent(inout) :: & 
         esnon     ! energy of melting for each snow layer (J/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(inout) :: & 
         aice0     ! concentration of open water

      integer (kind=int_kind), dimension(max_ntrcr), intent(in) :: &
         trcr_depend

      logical (kind=log_kind), intent(out) :: &
         l_stop   ! if true, abort on return

      integer (kind=int_kind), intent(out) :: &
         istop, jstop ! indices of grid cell where model aborts

      ! optional history fields
      real (kind=dbl_kind), dimension(nx_block,ny_block), &
         intent(inout), optional :: &
         dardg1dt  , & ! rate of fractional area loss by ridging ice (1/s)
         dardg2dt  , & ! rate of fractional area gain by new ridges (1/s)
         dvirdgdt  , & ! rate of ice volume ridged (m/s)
         opening   , & ! rate of opening due to divergence/shear (1/s)
         fresh     , & ! fresh water flux to ocean (kg/m^2/s)
         fhocn         ! net heat flux to ocean (W/m^2)

      real (kind=dbl_kind), dimension(nx_block,ny_block,n_aeromx), &
         intent(inout), optional :: &
         fsoot      ! 
!
!EOP
!
      real (kind=dbl_kind), dimension (icells) :: &
         asum       , & ! sum of ice and open water area
         aksum      , & ! ratio of area removed to area ridged
         msnow_mlt  , & ! mass of snow added to ocean (kg m-2)
         esnow_mlt  , & ! energy needed to melt snow in ocean (J m-2)
         closing_net, & ! net rate at which area is removed    (1/s)
                        ! (ridging ice area - area of new ridges) / dt_dyn
         divu_adv   , & ! divu as implied by transport scheme  (1/s)
         opning     , & ! rate of opening due to divergence/shear
                        ! opning is a local variable;
                        ! opening is the history diagnostic variable
         ardg1      , & ! fractional area loss by ridging ice
         ardg2      , & ! fractional area gain by new ridges
         virdg      , & ! ice volume ridged
         aopen          ! area opening due to divergence/shear

      real (kind=dbl_kind), dimension (icells,n_aeromx) :: &
         msoot       ! mass of soot added to ocean (kg m-2)

      real (kind=dbl_kind), dimension (icells,0:ncat) :: &
         apartic          ! participation function; fraction of ridging
                          ! and closing associated w/ category n

      real (kind=dbl_kind), dimension (icells,ncat) :: &
         hrmin        , & ! minimum ridge thickness
         hrmax        , & ! maximum ridge thickness (krdg_redist = 0)
         hrexp        , & ! ridge e-folding thickness (krdg_redist = 1) 
         krdg             ! mean ridge thickness/thickness of ridging ice

      real (kind=dbl_kind), dimension (icells) :: &
         vice_init, vice_final, & ! ice volume summed over categories
         vsno_init, vsno_final, & ! snow volume summed over categories
         eice_init, eice_final, & ! ice energy summed over layers
         esno_init, esno_final    ! snow energy summed over layers

      integer (kind=int_kind), parameter :: &
         nitermax = 20    ! max number of ridging iterations

      integer (kind=int_kind) :: &
         i,j          , & ! horizontal indices
         n            , & ! thickness category index
         niter        , & ! iteration counter
         ij               ! horizontal index, combines i and j loops

      real (kind=dbl_kind) :: &
         dti              ! 1 / dt_dyn or 1 / dt_thm

      logical (kind=log_kind) :: &
         iterate_ridging, & ! if true, repeat the ridging
         asum_error         ! flag for asum .ne. 1

      character (len=char_len) :: &
         fieldid            ! field identifier

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      l_stop = .false.
      istop = 0
      jstop = 0

      do ij = 1, icells
         msnow_mlt(ij) = c0
         esnow_mlt(ij) = c0
         msoot    (ij,:) = c0
         ardg1    (ij) = c0
         ardg2    (ij) = c0
         virdg    (ij) = c0
!         aopen    (ij) = c0
      enddo

      !-----------------------------------------------------------------
      ! Compute area of ice plus open water before ridging.
      !-----------------------------------------------------------------
      call asum_ridging (nx_block, ny_block,        &
                         icells,   indxi,    indxj, &
                         aicen,    aice0,           &
                         asum)

      !-----------------------------------------------------------------
      ! Compute the area opening and closing.
      !-----------------------------------------------------------------
      call ridge_prep (nx_block, ny_block,      &
                       icells,   indxi,  indxj, &
                       dt_dyn,                  &
                       rdg_conv,  rdg_shear,    &
                       asum,      closing_net,  &
                       divu_adv,  opning)

      !-----------------------------------------------------------------
      ! Compute initial values of conserved quantities. 
      !-----------------------------------------------------------------

      if (l_conservation_check) then

         call column_sum (nx_block,   ny_block,     &
                          icells,   indxi,   indxj, &
                          ncat,                     &
                          vicen,      vice_init)

         call column_sum (nx_block,   ny_block,     &
                          icells,   indxi,   indxj, &
                          ncat,                     &
                          vsnon,      vsno_init)

         call column_sum (nx_block,   ny_block,     &
                          icells,   indxi,   indxj, &
                          ntilyr,                   &
                          eicen,      eice_init)

         call column_sum (nx_block,   ny_block,     &
                          icells,   indxi,   indxj, &
                          ntslyr,                   &
                          esnon,      esno_init)

      endif            

      do niter = 1, nitermax

      !-----------------------------------------------------------------
      ! Compute the thickness distribution of ridging ice
      ! and various quantities associated with the new ridged ice.
      !-----------------------------------------------------------------

         call ridge_itd (nx_block,  ny_block,        &
                         icells,    indxi,    indxj, &
                         aicen,     vicen,           &
                         aice0,                      &
                         aksum,     apartic,         &
                         hrmin,     hrmax,           &
                         hrexp,     krdg)

      !-----------------------------------------------------------------
      ! Redistribute area, volume, and energy.
      !-----------------------------------------------------------------

         call ridge_shift (nx_block,  ny_block,        &
                           icells,    indxi,    indxj, &
                           ntrcr,     dt_dyn,          &
                           aicen,     trcrn,           &
                           vicen,     vsnon,           &
                           eicen,     esnon,           &
                           aice0,     trcr_depend,     &
                           aksum,     apartic,         &
                           hrmin,     hrmax,           &
                           hrexp,     krdg,            &
                           closing_net, opning,        &
                           ardg1,     ardg2,           &
                           virdg,     aopen,           &
                           msnow_mlt, esnow_mlt,       &
                           msoot,                      &
                           l_stop,                     &
                           istop,     jstop)

         if (l_stop) return

      !-----------------------------------------------------------------
      ! Make sure the new area = 1.  If not (because the closing
      ! and opening rates were reduced above), prepare to ridge again
      ! with new rates.
      !-----------------------------------------------------------------

         call asum_ridging (nx_block,  ny_block,      &
                            icells,    indxi,  indxj, &
                            aicen,     aice0,         &
                            asum)

         call ridge_check (nx_block,   ny_block,        &
                           icells,     indxi,    indxj, &
                           dt_dyn,                      &
                           asum,       closing_net,     &
                           divu_adv,   opning,          &
                           iterate_ridging)

      !-----------------------------------------------------------------
      ! If done, exit.  If not, prepare to ridge again.
      !-----------------------------------------------------------------

         if (iterate_ridging) then
            write(nu_diag,*) 'Repeat ridging, niter =', niter
         else
            exit
         endif

         if (niter == nitermax) then
            write(nu_diag,*) ' '
            write(nu_diag,*) 'Exceeded max number of ridging iterations'
            write(nu_diag,*) 'max =',nitermax
            l_stop = .true.
            return
         endif
            
      enddo                     ! niter

      !-----------------------------------------------------------------
      ! Compute final values of conserved quantities. 
      ! Check for conservation (allowing for snow thrown into ocean).
      !-----------------------------------------------------------------

      if (l_conservation_check) then

         call column_sum (nx_block,   ny_block,     &
                          icells,   indxi,   indxj, &
                          ncat,                     &
                          vicen,      vice_final)

         call column_sum (nx_block,   ny_block,     &
                          icells,   indxi,   indxj, &
                          ncat,                     &
                          vsnon,      vsno_final)

         call column_sum (nx_block,   ny_block,     &
                          icells,   indxi,   indxj, &
                          ntilyr,                   &
                          eicen,      eice_final)

         call column_sum (nx_block,   ny_block,     &
                          icells,   indxi,   indxj, &
                          ntslyr,                   &
                          esnon,      esno_final)

         do ij = 1, icells
            vsno_final(ij) = vsno_final(ij) + msnow_mlt(ij)/rhos
            esno_final(ij) = esno_final(ij) + esnow_mlt(ij)
         enddo

         fieldid = 'vice, ridging'
         call column_conservation_check (nx_block,  ny_block,       &
                                         icells,    indxi,   indxj, &
                                         fieldid,                   &
                                         vice_init, vice_final,     &
                                         puny,      l_stop,         &
                                         istop,     jstop)
         if (l_stop) return

         fieldid = 'vsno, ridging'
         call column_conservation_check (nx_block,  ny_block,       &
                                         icells,    indxi,   indxj, &
                                         fieldid,                   &
                                         vsno_init, vsno_final,     &
                                         puny,      l_stop,         &
                                         istop,     jstop)
         if (l_stop) return

         fieldid = 'eice, ridging'
         call column_conservation_check (nx_block,  ny_block,       &
                                         icells,    indxi,   indxj, &
                                         fieldid,                   &
                                         eice_init, eice_final,     &
                                         puny*Lfresh*rhoi,          &
                                         l_stop,                    &
                                         istop,     jstop)
         if (l_stop) return

         fieldid = 'esno, ridging'
         call column_conservation_check (nx_block,  ny_block,       &
                                         icells,    indxi,   indxj, &
                                         fieldid,                   &
                                         esno_init, esno_final,     &
                                         puny*Lfresh*rhos,          &
                                         l_stop,                    &
                                         istop,     jstop)
         if (l_stop) return
         
      endif                     ! l_conservation_check            

      !-----------------------------------------------------------------
      ! Compute ridging diagnostics.
      !-----------------------------------------------------------------

      dti = c1/dt_dyn

      if (present(dardg1dt)) then
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            dardg1dt(i,j) = ardg1(ij)*dti
         enddo
      endif
      if (present(dardg2dt)) then
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            dardg2dt(i,j) = ardg2(ij)*dti
         enddo
      endif
      if (present(dvirdgdt)) then
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            dvirdgdt(i,j) = virdg(ij)*dti
         enddo
      endif
      if (present(opening)) then
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            opening(i,j) = aopen(ij)*dti
         enddo
      endif

      !-----------------------------------------------------------------
      ! Update fresh water and heat fluxes due to snow melt.
      !-----------------------------------------------------------------

      dti = c1/dt_thm

      if (present(fsoot)) then
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            fsoot(i,j,:) = fsoot(i,j,:) + msoot(ij,:)*dti
         enddo
      endif
      if (present(fresh)) then
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            fresh(i,j) = fresh(i,j) + msnow_mlt(ij)*dti
         enddo
      endif
      if (present(fhocn)) then
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            fhocn(i,j) = fhocn(i,j) + esnow_mlt(ij)*dti
         enddo
      endif

      !-----------------------------------------------------------------
      ! Check for fractional ice area > 1.
      !-----------------------------------------------------------------

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
         if (abs(asum(ij) - c1) > puny) then
            l_stop = .true.
            istop = i
            jstop = j

            write(nu_diag,*) ' '
            write(nu_diag,*) 'Ridging error: total area > 1'
            write(nu_diag,*) 'i, j, area:', i, j, asum(ij)
            write(nu_diag,*) 'n, aicen:'
            write(nu_diag,*)  0, aice0(i,j)
            do n = 1, ncat
               write(nu_diag,*) n, aicen(i,j,n)
            enddo
            return
         endif
      enddo

      end subroutine ridge_ice

!=======================================================================
!BOP
!
! !ROUTINE: asum_ridging - find total fractional area
!
! !DESCRIPTION:
!
! Find the total area of ice plus open water in each grid cell.
!
! This is similar to the aggregate_area subroutine except that the
! total area can be greater than 1, so the open water area is
! included in the sum instead of being computed as a residual.
!
! !REVISION HISTORY:
!
! author: William H. Lipscomb, LANL
!
! !INTERFACE:
!
      subroutine asum_ridging (nx_block, ny_block,        &
                               icells,   indxi,    indxj, &
                               aicen,    aice0,           &
                               asum)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells                ! number of cells with ice present

      integer (kind=int_kind), dimension (nx_block*ny_block), &
         intent(in) :: &
         indxi, indxj     ! compressed indices for cells with ice


      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat), &
         intent(in) :: &
         aicen          ! concentration of ice in each category

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         aice0          ! concentration of open water

      real (kind=dbl_kind), dimension (icells), intent(out):: &
         asum           ! sum of ice and open water area
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j, n, &
         ij               ! horizontal index, combines i and j loops

      !-----------------------------------------------------------------
      ! open water
      !-----------------------------------------------------------------
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
         asum(ij) = aice0(i,j)
      enddo

      !-----------------------------------------------------------------
      ! ice categories
      !-----------------------------------------------------------------

      do n = 1, ncat
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            asum(ij) = asum(ij) + aicen(i,j,n)
         enddo
      enddo

      end subroutine asum_ridging

!=======================================================================
!BOP
!
! !ROUTINE: ridge_prep - preparation for ridging
!
! !DESCRIPTION: Initialize arrays, compute area of closing and opening
!
!
! !REVISION HISTORY:
!
! author: William H. Lipscomb, LANL
!
! !INTERFACE:
!
      subroutine ridge_prep (nx_block,   ny_block,        &
                             icells,     indxi,    indxj, &
                             dt,                          &
                             rdg_conv,   rdg_shear,       &
                             asum,       closing_net,     &
                             divu_adv,   opning)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells                ! number of cells with ice present

      integer (kind=int_kind), dimension (nx_block*ny_block), &
         intent(in) :: &
         indxi, indxj     ! compressed indices for cells with ice


      real (kind=dbl_kind), intent(in) :: &
         dt                  ! time step (s)

      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(in) :: &
         rdg_conv, & ! normalized energy dissipation due to convergence (1/s)
         rdg_shear   ! normalized energy dissipation due to shear (1/s)

      real (kind=dbl_kind), dimension(icells), &
         intent(inout):: &
         asum      ! sum of ice and open water area

      real (kind=dbl_kind), dimension(icells), &
         intent(out):: &
         closing_net, & ! net rate at which area is removed    (1/s)
         divu_adv   , & ! divu as implied by transport scheme  (1/s)
         opning         ! rate of opening due to divergence/shear
!
!EOP
!
      real (kind=dbl_kind), parameter :: &
         big = 1.0e+8_dbl_kind

      integer (kind=int_kind) :: &
         i,j, &         ! horizontal indices
         ij             ! horizontal index, combines i and j loops

      ! Set hin_max(ncat) to a big value to ensure that all ridged ice
      ! is thinner than hin_max(ncat).
      hin_max(ncat) = big

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

      !-----------------------------------------------------------------
      ! Compute the net rate of closing due to convergence
      ! and shear, based on Flato and Hibler (1995).
      !
      ! For the elliptical yield curve:
      !    rdg_conv  = -min (divu, 0)
      !    rdg_shear = (1/2) * (Delta - abs(divu))
      ! Note that the shear term also accounts for divergence.
      !
      ! The energy dissipation rate is equal to the net closing rate
      ! times the ice strength.
      !
      ! NOTE: The NET closing rate is equal to the rate that open water
      !  area is removed, plus the rate at which ice area is removed by
      !  ridging, minus the rate at which area is added in new ridges.
      !  The GROSS closing rate is equal to the first two terms (open
      !  water closing and thin ice ridging) without the third term
      !  (thick, newly ridged ice).
      !-----------------------------------------------------------------

         closing_net(ij) = Cs*rdg_shear(i,j) + rdg_conv(i,j)

      !-----------------------------------------------------------------
      ! Compute divu_adv, the divergence rate given by the transport/
      ! advection scheme, which may not be equal to divu as computed
      ! from the velocity field.
      !
      ! If divu_adv < 0, make sure the closing rate is large enough
      ! to give asum = 1.0 after ridging.
      !-----------------------------------------------------------------

         divu_adv(ij) = (c1-asum(ij)) / dt

         if (divu_adv(ij) < c0) &
              closing_net(ij) = max(closing_net(ij), -divu_adv(ij))

      !-----------------------------------------------------------------
      ! Compute the (non-negative) opening rate that will give
      ! asum = 1.0 after ridging.
      !-----------------------------------------------------------------
         opning(ij) = closing_net(ij) + divu_adv(ij)

      enddo

      end subroutine ridge_prep

!=======================================================================
!BOP
!
! !ROUTINE: ridge_itd - thickness distribution of ridging and ridged ice
!
! !DESCRIPTION:
!
! Compute the thickness distribution of the ice and open water
! participating in ridging and of the resulting ridges.
!
! This version includes new options for ridging participation and
!  redistribution.
! The new participation scheme (krdg_partic = 1) improves stability
!  by increasing the time scale for large changes in ice strength.
! The new exponential redistribution function (krdg_redist = 1) improves 
!  agreement between ITDs of modeled and observed ridges.   
!
! !REVISION HISTORY:
!
! author: William H. Lipscomb, LANL
!
! 2006: Changed subroutine name to ridge_itd
!       Added new options for ridging participation and redistribution.  
!
! !INTERFACE:
!
      subroutine ridge_itd (nx_block,    ny_block,        &
                            icells,      indxi,    indxj, &
                            aicen,       vicen,           &
                            aice0,                        &
                            aksum,       apartic,         &
                            hrmin,       hrmax,           &
                            hrexp,       krdg)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells                ! number of cells with ice present

      integer (kind=int_kind), dimension (nx_block*ny_block), &
         intent(in) :: &
         indxi, indxj     ! compressed indices for cells with ice


      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat), &
         intent(in) :: &
         aicen , & ! concentration of ice
         vicen     ! volume per unit area of ice (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         aice0     ! concentration of open water

      real (kind=dbl_kind), dimension (icells), intent(out):: &
         aksum            ! ratio of area removed to area ridged

      real (kind=dbl_kind), dimension (icells,0:ncat), &
         intent(out) :: &
         apartic          ! participation function; fraction of ridging
                          ! and closing associated w/ category n

      real (kind=dbl_kind), dimension (icells,ncat), &
         intent(out) :: &
         hrmin        , & ! minimum ridge thickness
         hrmax        , & ! maximum ridge thickness (krdg_redist = 0)
         hrexp        , & ! ridge e-folding thickness (krdg_redist = 1) 
         krdg             ! mean ridge thickness/thickness of ridging ice
!
!EOP
!
      integer (kind=int_kind) :: &
         i,j          , & ! horizontal indices
         n            , & ! thickness category index
         ij               ! horizontal index, combines i and j loops

      real (kind=dbl_kind), parameter :: &
         Gstari   = c1/Gstar, &
         astari   = c1/astar

      real (kind=dbl_kind), dimension(icells,-1:ncat) :: &
         Gsum             ! Gsum(n) = sum of areas in categories 0 to n

      real (kind=dbl_kind), dimension(icells) :: &
         work             ! temporary work array

      real (kind=dbl_kind) :: &
         hi           , & ! ice thickness for each cat (m)
         hieff        , & ! effective ice thickness (m) (krdg_redist = 2)
         hrmean       , & ! mean ridge thickness (m)
         xtmp             ! temporary variable

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      do ij = 1, icells
         Gsum   (ij,-1) = c0        ! by definition
         Gsum   (ij,0)  = c1         ! to avoid divzero below
         apartic(ij,0)  = c0
      enddo

      do n = 1, ncat
         do ij = 1, icells
            Gsum   (ij,n) = c1    ! to avoid divzero below
            apartic(ij,n) = c0
            hrmin  (ij,n) = c0
            hrmax  (ij,n) = c0
            hrexp  (ij,n) = c0
            krdg   (ij,n) = c1
         enddo
      enddo

      !-----------------------------------------------------------------
      ! Compute the thickness distribution of ice participating in ridging.
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! First compute the cumulative thickness distribution function Gsum,
      !  where Gsum(n) is the fractional area in categories 0 to n.
      ! Ignore categories with very small areas.
      !-----------------------------------------------------------------

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
         if (aice0(i,j) > puny) then
            Gsum(ij,0) = aice0(i,j)
         else
            Gsum(ij,0) = Gsum(ij,-1)
         endif
      enddo

      do n = 1, ncat
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            if (aicen(i,j,n) > puny) then
               Gsum(ij,n) = Gsum(ij,n-1) + aicen(i,j,n)
            else
               Gsum(ij,n) = Gsum(ij,n-1)
            endif
         enddo
      enddo

      ! normalize

      do ij = 1, icells
         work(ij) = c1 / Gsum(ij,ncat)
      enddo
      do n = 0, ncat
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            Gsum(ij,n) = Gsum(ij,n) * work(ij)
         enddo
      enddo

      !-----------------------------------------------------------------
      ! Compute the participation function apartic; this is analogous to
      ! a(h) = b(h)g(h) as defined in Thorndike et al. (1975).
      !
      !                area lost from category n due to ridging/closing
      !  apartic(n) = --------------------------------------------------
      !                  total area lost due to ridging/closing
      !
      !-----------------------------------------------------------------

      if (krdg_partic == 0) then  ! Thornike et al. 1975 formulation

      !-----------------------------------------------------------------
      ! Assume b(h) = (2/Gstar) * (1 - G(h)/Gstar).
      ! The expressions for apartic are found by integrating b(h)g(h) between
      ! the category boundaries.
      !-----------------------------------------------------------------

         do n = 0, ncat
            do ij = 1, icells
               if (Gsum(ij,n) < Gstar) then
                  apartic(ij,n) = Gstari*(Gsum(ij,n)-Gsum(ij,n-1)) * &
                       (c2 - (Gsum(ij,n-1)+Gsum(ij,n))*Gstari)
               elseif (Gsum(ij,n-1) < Gstar) then
                  apartic(ij,n) = Gstari * (Gstar-Gsum(ij,n-1)) * &
                       (c2 - (Gsum(ij,n-1)+Gstar)*Gstari)
               endif
            enddo               ! ij
         enddo                  ! n

      elseif (krdg_partic==1) then   ! exponential dependence on G(h)

      !-----------------------------------------------------------------
      ! b(h) = exp(-G(h)/astar)
      ! apartic(n) = [exp(-G(n-1)/astar - exp(-G(n)/astar] / [1-exp(-1/astar)]. 
      ! The expression for apartic is found by integrating b(h)g(h)
      ! between the category boundaries.
      !-----------------------------------------------------------------

         ! precompute exponential terms using Gsum as work array

         xtmp = c1 / (c1 - exp(-astari))

         do n = -1, ncat
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, icells
               Gsum(ij,n) = exp(-Gsum(ij,n)*astari) * xtmp
            enddo               ! ij
         enddo                  ! n

         do n = 0, ncat
            do ij = 1, icells
               apartic(ij,n) = Gsum(ij,n-1) - Gsum(ij,n)
            enddo               ! ij
         enddo                  ! n

      endif                     ! krdg_partic

      !-----------------------------------------------------------------
      ! Compute variables related to ITD of ridged ice:
      ! 
      ! krdg = mean ridge thickness/ thickness of ridging ice
      ! hrmin  = min ridge thickness
      ! hrmax  = max ridge thickness (krdg_redist = 0)
      ! hrexp  = ridge e-folding scale (krdg_redist = 1)
      !----------------------------------------------------------------

      if (krdg_redist == 0) then  ! Hibler 1980 formulation

      !-----------------------------------------------------------------
      ! Assume ridged ice is uniformly distributed between hrmin and hrmax. 
      ! 
      ! This parameterization is a modified version of Hibler (1980). 
      ! In the original paper the min ridging thickness is hrmin = 2*hi, 
      !  and the max thickness is hrmax = 2*sqrt(hi*Hstar). 
      ! 
      ! Here the min thickness is hrmin = min(2*hi, hi+maxraft), 
      !  so thick ridging ice is not required to raft. 
      !
      !-----------------------------------------------------------------

         do n = 1, ncat
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)

               if (aicen(i,j,n) > puny) then 
                  hi = vicen(i,j,n) / aicen(i,j,n) 
                  hrmin(ij,n) = min(c2*hi, hi + maxraft) 
                  hrmax(ij,n) = c2*sqrt(Hstar*hi) 
                  hrmax(ij,n) = max(hrmax(ij,n), hrmin(ij,n)+puny) 
                  hrmean = p5 * (hrmin(ij,n) + hrmax(ij,n)) 
                  krdg(ij,n) = hrmean / hi 
               endif 

            enddo               ! ij
         enddo                  ! n

      else               ! krdg_redist = 1; exponential redistribution
 
      !----------------------------------------------------------------- 
      ! The ridge ITD is a negative exponential: 
      ! 
      !  g(h) ~ exp[-(h-hrmin)/hrexp], h >= hrmin 
      ! 
      ! where hrmin is the minimum thickness of ridging ice and 
      ! hrexp is the e-folding thickness.
      ! 
      ! Here, assume as above that hrmin = min(2*hi, hi+maxraft).
      ! That is, the minimum ridge thickness results from rafting,
      !  unless the ice is thicker than maxraft.
      !
      ! Also, assume that hrexp = mu_rdg*sqrt(hi).
      ! The parameter mu_rdg is tuned to give e-folding scales mostly
      !  in the range 2-4 m as observed by upward-looking sonar.
      !
      ! Values of mu_rdg in the right column give ice strengths
      !  roughly equal to values of Hstar in the left column
      !  (within ~10 kN/m for typical ITDs):
      !
      !   Hstar     mu_rdg
      !
      !     25        3.0
      !     50        4.0
      !     75        5.0
      !    100        6.0
      !----------------------------------------------------------------- 

         do n = 1, ncat
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)
               if (aicen(i,j,n) > puny) then
                  hi = vicen(i,j,n) / aicen(i,j,n)
                  hi = max(hi,puny)
                  hrmin(ij,n) = min(c2*hi, hi + maxraft)
                  hrexp(ij,n) = mu_rdg * sqrt(hi)
                  krdg(ij,n) = (hrmin(ij,n) + hrexp(ij,n)) / hi
               endif
            enddo
         enddo

      endif                     ! krdg_redist

      !----------------------------------------------------------------
      ! Compute aksum = net ice area removed / total area participating.
      ! For instance, if a unit area of ice with h = 1 participates in
      !  ridging to form a ridge with a = 1/3 and h = 3, then
      !  aksum = 1 - 1/3 = 2/3.
      !---------------------------------------------------------------- 

      do ij = 1, icells
         aksum(ij) = apartic(ij,0) ! area participating = area removed
      enddo

      do n = 1, ncat
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            ! area participating > area removed
            aksum(ij) = aksum(ij)  &
                       + apartic(ij,n) * (c1 - c1/krdg(ij,n)) 
         enddo
      enddo

      end subroutine ridge_itd

!=======================================================================
!BOP
!
! !ROUTINE: ridge_shift - shift ridging ice among thickness categories
!
! !DESCRIPTION:
!
! Remove area, volume, and energy from each ridging category
! and add to thicker ice categories.
!
! !REVISION HISTORY:
!
! author: William H. Lipscomb, LANL
!
! !INTERFACE:
!
      subroutine ridge_shift (nx_block,    ny_block,        &
                              icells,      indxi,    indxj, &
                              ntrcr,       dt,              &
                              aicen,       trcrn,           &
                              vicen,       vsnon,           &
                              eicen,       esnon,           &
                              aice0,       trcr_depend,     &   
                              aksum,       apartic,         &
                              hrmin,       hrmax,           &
                              hrexp,       krdg,            &
                              closing_net, opning,          &
                              ardg1,       ardg2,           &
                              virdg,       aopen,           &
                              msnow_mlt,   esnow_mlt,       &
                              msoot,                        &
                              l_stop,                       &
                              istop,       jstop)
!
! !USES:
!
      use ice_state, only: nt_aero, &
                           nt_alvl, nt_vlvl, tr_lvl
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells            , & ! number of cells with ice present
         ntrcr                 ! number of tracers in use

      integer (kind=int_kind), dimension (nx_block*ny_block), &
         intent(in) :: &
         indxi, indxj     ! compressed indices for cells with ice

      real (kind=dbl_kind), intent(in) :: &
         dt                  ! time step (s)

      integer (kind=int_kind), dimension (max_ntrcr), intent(in) :: &
         trcr_depend ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(inout) :: &
         aice0     ! concentration of open water

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat), &
         intent(inout) :: &
         aicen , & ! concentration of ice
         vicen , & ! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_ntrcr,ncat), &
         intent(inout) :: &
         trcrn     ! ice tracers

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntilyr), &
         intent(inout) :: &
         eicen     ! energy of melting for each ice layer (J/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntslyr), &
         intent(inout) :: &
         esnon     ! energy of melting for each snow layer (J/m^2)

      real (kind=dbl_kind), dimension (icells), intent(in) :: &
         aksum             ! ratio of area removed to area ridged

      real (kind=dbl_kind), dimension (icells,0:ncat), intent(in) :: &
         apartic          ! participation function; fraction of ridging
                          ! and closing associated w/ category n

      real (kind=dbl_kind), dimension (icells,ncat), intent(in) :: &
         hrmin        , & ! minimum ridge thickness
         hrmax        , & ! maximum ridge thickness (krdg_redist = 0)
         hrexp        , & ! ridge e-folding thickness (krdg_redist = 1) 
         krdg             ! mean ridge thickness/thickness of ridging ice

      real (kind=dbl_kind), dimension(icells), intent(inout) :: &
         closing_net, & ! net rate at which area is removed    (1/s)
         opning     , & ! rate of opening due to divergence/shear (1/s)
         ardg1      , & ! fractional area loss by ridging ice
         ardg2      , & ! fractional area gain by new ridges
         virdg      , & ! ice volume ridged (m)
         aopen          ! area opened due to divergence/shear

      real (kind=dbl_kind), dimension(icells), intent(inout) :: &
         msnow_mlt, & ! mass of snow added to ocean (kg m-2)
         esnow_mlt    ! energy needed to melt snow in ocean (J m-2)

      real (kind=dbl_kind), dimension(icells,n_aeromx), intent(inout) :: &
         msoot      ! mass of soot added to ocean (kg m-2)

      logical (kind=log_kind), intent(inout) :: &
         l_stop   ! if true, abort on return

      integer (kind=int_kind), intent(inout) :: &
         istop, jstop ! indices of grid cell where model aborts
!
!EOP
!
      integer (kind=int_kind) :: &
         i,j           , & ! horizontal indices
         n, nr         , & ! thickness category indices
         k             , & ! ice layer index
         it            , & ! tracer index
         ij, m         , & ! horizontal indices, combine i and j loops
         iridge            ! number of cells with nonzero ridging
      integer (kind=int_kind) :: &
         iaero             ! index for number of aerosol tracers

      integer (kind=int_kind), dimension (icells) :: &
         indxii, indxjj  , & ! compressed indices
         indxij              ! compressed indices

      real (kind=dbl_kind), dimension (icells,ncat) :: &
         aicen_init    , & ! ice area before ridging
         vicen_init    , & ! ice volume before ridging
         vsnon_init        ! snow volume before ridging

      real (kind=dbl_kind), dimension (icells,ntilyr) :: &
         eicen_init        ! ice energy before ridging

      real (kind=dbl_kind), dimension (icells,ntslyr) :: &
         esnon_init        ! snow energy before ridging

      real (kind=dbl_kind), dimension(icells,max_ntrcr,ncat) :: &
         atrcrn            ! aicen*trcrn

      real (kind=dbl_kind), dimension (icells) :: &
         closing_gross     ! rate at which area removed, not counting
                           ! area of new ridges

! ECH note:  the following arrays only need be defined on iridge cells
      real (kind=dbl_kind), dimension (icells) :: &
         afrac         , & ! fraction of category area ridged
         ardg1n        , & ! area of ice ridged
         ardg2n        , & ! area of new ridges
         virdgn        , & ! ridging ice volume
         vsrdgn        , & ! ridging snow volume
         dhr           , & ! hrmax - hrmin
         dhr2          , & ! hrmax^2 - hrmin^2
         farea         , & ! fraction of new ridge area going to nr
         fvol              ! fraction of new ridge volume going to nr

      real (kind=dbl_kind), dimension (icells,nilyr) :: &
         eirdgn            ! ridging ice energy

      real (kind=dbl_kind), dimension (icells,nslyr) :: &
         esrdgn            ! ridging snow energy

      real (kind=dbl_kind) :: &
         hi1           , & ! thickness of ridging ice
         hexp          , & ! ridge e-folding thickness
         hL, hR        , & ! left and right limits of integration
         expL, expR    , & ! exponentials involving hL, hR
         tmpfac        , & ! factor by which opening/closing rates are cut
         wk1               ! work variable

      !-----------------------------------------------------------------
      ! Define variables equal to aicen*trcrn, vicen*trcrn, vsnon*trcrn
      !-----------------------------------------------------------------

      do n = 1, ncat
         do it = 1, ntrcr
            if (trcr_depend(it) == 0) then ! ice area tracer
               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)
                  atrcrn(ij,it,n) = aicen(i,j,n)*trcrn(i,j,it,n)
               enddo
            elseif (trcr_depend(it) == 1) then  ! ice volume tracer
               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)
                  atrcrn(ij,it,n) = vicen(i,j,n)*trcrn(i,j,it,n)
               enddo
            elseif (trcr_depend(it) == 2) then  ! snow volume tracer
               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)
                  atrcrn(ij,it,n) = vsnon(i,j,n)*trcrn(i,j,it,n)
               enddo 
            endif
         enddo
      enddo

      !-----------------------------------------------------------------
      ! Based on the ITD of ridging and ridged ice, convert the net
      !  closing rate to a gross closing rate.
      ! NOTE: 0 < aksum <= 1
      !-----------------------------------------------------------------

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         closing_gross(ij) = closing_net(ij) / aksum(ij)

      !-----------------------------------------------------------------
      ! Reduce the closing rate if more than 100% of the open water
      ! would be removed.  Reduce the opening rate proportionately.
      !-----------------------------------------------------------------

         if (apartic(ij,0) > c0) then
            wk1 = apartic(ij,0) * closing_gross(ij) * dt
            if (wk1 > aice0(i,j)) then
               tmpfac = aice0(i,j) / wk1
               closing_gross(ij) = closing_gross(ij) * tmpfac
               opning(ij) = opning(ij) * tmpfac
            endif
         endif

      enddo                     ! ij

      !-----------------------------------------------------------------
      ! Reduce the closing rate if more than 100% of any ice category
      ! would be removed.  Reduce the opening rate proportionately.
      !-----------------------------------------------------------------
      do n = 1, ncat
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            if (aicen(i,j,n) > puny .and. apartic(ij,n) > c0) then
               wk1 = apartic(ij,n) * closing_gross(ij) * dt
               if (wk1 > aicen(i,j,n)) then
                  tmpfac = aicen(i,j,n) / wk1
                  closing_gross(ij) = closing_gross(ij) * tmpfac
                  opning(ij) = opning(ij) * tmpfac
               endif
            endif

         enddo                  ! ij
      enddo                     ! n

      !-----------------------------------------------------------------
      ! Compute change in open water area due to closing and opening.
      !-----------------------------------------------------------------

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
         aice0(i,j) = aice0(i,j) &
                    - apartic(ij,0)*closing_gross(ij)*dt &
                    + opning(ij)*dt
         if (aice0(i,j) < -puny) then
            l_stop = .true.
            istop = i
            jstop = j

            write (nu_diag,*) ' '
            write (nu_diag,*) 'Ridging error: aice0 < 0'
            write (nu_diag,*) 'i, j, aice0:', i, j, aice0(i,j)
            return

         elseif (aice0(i,j) < c0) then    ! roundoff error
            aice0(i,j) = c0
         endif

         aopen(ij) = opning(ij)*dt  ! optional diagnostic

      enddo

      !-----------------------------------------------------------------
      ! Save initial state variables
      !-----------------------------------------------------------------

      do n = 1, ncat
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            aicen_init(ij,n) = aicen(i,j,n)
            vicen_init(ij,n) = vicen(i,j,n)
            vsnon_init(ij,n) = vsnon(i,j,n)
         enddo
      enddo

      do n = 1, ntilyr
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            eicen_init(ij,n) = eicen(i,j,n)
         enddo
      enddo

      do n = 1, ntslyr
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            esnon_init(ij,n) = esnon(i,j,n)
         enddo
      enddo

      !-----------------------------------------------------------------
      ! Compute the area, volume, and energy of ice ridging in each
      !  category, along with the area of the resulting ridge.
      !-----------------------------------------------------------------

      do n = 1, ncat

      !-----------------------------------------------------------------
      ! Identify grid cells with nonzero ridging
      !-----------------------------------------------------------------

         iridge = 0
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            if (aicen_init(ij,n) > puny .and. apartic(ij,n) > c0 &
                 .and. closing_gross(ij) > c0) then
               iridge = iridge + 1
               indxii(iridge) = i
               indxjj(iridge) = j
               indxij(iridge) = ij
            endif
         enddo                  ! ij

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, iridge
            i = indxii(ij)
            j = indxjj(ij)
            m = indxij(ij)

      !-----------------------------------------------------------------
      ! Compute area of ridging ice (ardg1n) and of new ridge (ardg2n).
      ! Make sure ridging fraction <=1.  (Roundoff errors can give
      !  ardg1 slightly greater than aicen.) 
      !-----------------------------------------------------------------

            ardg1n(ij) = apartic(m,n)*closing_gross(m)*dt

            if (ardg1n(ij) > aicen_init(m,n) + puny) then
               l_stop = .true.
               istop = i
               jstop = j

               write (nu_diag,*) ' '
               write (nu_diag,*) 'Ridging error: ardg > aicen'
               write (nu_diag,*) 'i, j, n:', i, j, n
               write (nu_diag,*) 'ardg, aicen:', &
                    ardg1n(ij), aicen_init(m,n)
               return
            else
               ardg1n(ij) = min(aicen_init(m,n), ardg1n(ij))
            endif

            ardg2n(ij) = ardg1n(ij) / krdg(m,n)
            afrac(ij) = ardg1n(ij) / aicen_init(m,n)

      !-----------------------------------------------------------------
      ! Subtract area, volume, and energy from ridging category n.
      ! Note: Tracer values are unchanged.
      !-----------------------------------------------------------------

            virdgn(ij) = vicen_init(m,n) * afrac(ij)
            vsrdgn(ij) = vsnon_init(m,n) * afrac(ij)

            aicen(i,j,n) = aicen(i,j,n) - ardg1n(ij)
            vicen(i,j,n) = vicen(i,j,n) - virdgn(ij)
            vsnon(i,j,n) = vsnon(i,j,n) - vsrdgn(ij)

      !-----------------------------------------------------------------
      ! Increment ridging diagnostics
      !-----------------------------------------------------------------

            ardg1(m) = ardg1(m) + ardg1n(ij)
            ardg2(m) = ardg2(m) + ardg2n(ij)
            virdg(m) = virdg(m) + virdgn(ij)

      !-----------------------------------------------------------------
      ! Decrement level ice area and volume tracers
      !-----------------------------------------------------------------

            if (tr_lvl) then

               ! Assume level and ridged ice both ridge, proportionally.
               ! Subtract the level ice portion of the ridging ice from
               ! the level ice tracers.
               atrcrn(m,nt_alvl,n) = atrcrn(m,nt_alvl,n) * (c1 - afrac(ij))
               atrcrn(m,nt_vlvl,n) = atrcrn(m,nt_vlvl,n) * (c1 - afrac(ij))

            endif

      !-----------------------------------------------------------------
      !  Place part of the snow lost by ridging into the ocean.
      !-----------------------------------------------------------------

            msnow_mlt(m) = msnow_mlt(m) + rhos*vsrdgn(ij)*(c1-fsnowrdg)

      !-----------------------------------------------------------------
      !  Place part of the soot lost by ridging into the ocean.
      !-----------------------------------------------------------------

            if (n_aero >= 1) then
               do iaero=1,n_aero
                msoot(m,iaero) = msoot(m,iaero) &
                        + vsrdgn(ij)*(c1-fsnowrdg) &
                        *(trcrn(i,j,nt_aero  +4*(iaero-1),n)   &
                        + trcrn(i,j,nt_aero+1+4*(iaero-1),n))
               enddo
            endif

      !-----------------------------------------------------------------
      ! Compute quantities used to apportion ice among categories
      ! in the nr loop below
      !-----------------------------------------------------------------

            dhr(ij)  = hrmax(ij,n) - hrmin(m,n)
            dhr2(ij) = hrmax(ij,n) * hrmax(ij,n) - hrmin(m,n) * hrmin(m,n)

         enddo                  ! ij

      !-----------------------------------------------------------------
      ! Subtract ice energy from ridging category n.
      !-----------------------------------------------------------------

         do k = 1, nilyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, iridge
               i = indxii(ij)
               j = indxjj(ij)
               m = indxij(ij)

               eirdgn(ij,k) = eicen_init(m,ilyr1(n)+k-1) * afrac(ij)
               eicen(i,j,ilyr1(n)+k-1) = eicen (i,j,ilyr1(n)+k-1) &
                                       - eirdgn(ij,k)
            enddo
         enddo

      !-----------------------------------------------------------------
      ! Subtract snow energy from ridging category n.
      ! Increment energy needed to melt snow in ocean.
      ! Note that esnow_mlt < 0; the ocean must cool to melt snow.
      !-----------------------------------------------------------------

         do k = 1, nslyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, iridge
               i = indxii(ij)
               j = indxjj(ij)
               m = indxij(ij)

               esrdgn(ij,k) = esnon_init(m,slyr1(n)+k-1) * afrac(ij)
               esnon(i,j,slyr1(n)+k-1) = esnon (i,j,slyr1(n)+k-1) &
                                       - esrdgn(ij,k)
               esnow_mlt(m) = esnow_mlt(m) &
                              + esrdgn(ij,k)*(c1-fsnowrdg)
           enddo
         enddo

      !-----------------------------------------------------------------
      ! Subtract area- and volume-weighted tracers from category n.
      !-----------------------------------------------------------------

         do it = 1, ntrcr
            if (trcr_depend(it) == 0) then ! ice area tracer
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
               do ij = 1, iridge
                  i = indxii(ij)
                  j = indxjj(ij)
                  m = indxij(ij)
                  atrcrn(m,it,n) = atrcrn(m,it,n) &
                                   - ardg1n(ij)*trcrn(i,j,it,n)
               enddo

            elseif (trcr_depend(it) == 1) then ! ice volume tracer
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
               do ij = 1, iridge
                  i = indxii(ij)
                  j = indxjj(ij)
                  m = indxij(ij)
                  atrcrn(m,it,n) = atrcrn(m,it,n) &
                                   - virdgn(ij)*trcrn(i,j,it,n)
               enddo

            elseif (trcr_depend(it) == 2) then ! snow volume tracer
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
               do ij = 1, iridge
                  i = indxii(ij)
                  j = indxjj(ij)
                  m = indxij(ij)
                  atrcrn(m,it,n) = atrcrn(m,it,n) &
                                   - vsrdgn(ij)*trcrn(i,j,it,n)
               enddo
            endif               ! trcr_depend
         enddo                  ! ntrcr


      !-----------------------------------------------------------------
      ! Add area, volume, and energy of new ridge to each category nr.
      !-----------------------------------------------------------------

         do nr = 1, ncat

            if (krdg_redist == 0) then ! Hibler 1980 formulation

               do ij = 1, iridge
                  m = indxij(ij)

      !-----------------------------------------------------------------
      ! Compute the fraction of ridged ice area and volume going to
      !  thickness category nr.
      !-----------------------------------------------------------------

                  if (hrmin(m,n) >= hin_max(nr) .or. &
                      hrmax(ij,n) <= hin_max(nr-1)) then
                     hL = c0
                     hR = c0
                  else
                     hL = max (hrmin(m,n), hin_max(nr-1))
                     hR = min (hrmax(ij,n), hin_max(nr))
                  endif

                  farea(ij) = (hR-hL) / dhr(ij)
                  fvol (ij) = (hR*hR - hL*hL) / dhr2(ij)

               enddo            ! ij

            else         ! krdg_redist = 1; 2005 exponential formulation

      !-----------------------------------------------------------------
      ! Compute the fraction of ridged ice area and volume going to
      !  thickness category nr.
      !-----------------------------------------------------------------

               if (nr < ncat) then

                  do ij = 1, iridge
                     m = indxij(ij)

                     hi1  = hrmin(m,n)
                     hexp = hrexp(m,n)

                     if (hi1 >= hin_max(nr)) then
                        farea(ij) = c0
                        fvol (ij) = c0
                     else
                        hL = max (hi1, hin_max(nr-1))
                        hR = hin_max(nr)
                        expL = exp(-(hL-hi1)/hexp)
                        expR = exp(-(hR-hi1)/hexp)
                        farea(ij) = expL - expR
                        fvol (ij) = ((hL + hexp)*expL  &
                                    - (hR + hexp)*expR) / (hi1 + hexp)
                     endif
                  enddo         ! ij

               else             ! nr = ncat

                  do ij = 1, iridge
                     m = indxij(ij)

                     hi1  = hrmin(m,n)
                     hexp = hrexp(m,n)

                     hL = max (hi1, hin_max(nr-1))
                     expL = exp(-(hL-hi1)/hexp)
                     farea(ij) = expL
                     fvol (ij) = (hL + hexp)*expL / (hi1 + hexp)

                  enddo

               endif            ! nr < ncat

            endif               ! krdg_redist

      !-----------------------------------------------------------------
      ! Transfer ice area, ice volume, and snow volume to category nr.
      !-----------------------------------------------------------------

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, iridge
               i = indxii(ij)
               j = indxjj(ij)
               aicen(i,j,nr) = aicen(i,j,nr) + farea(ij)*ardg2n(ij)
               vicen(i,j,nr) = vicen(i,j,nr) + fvol(ij) *virdgn(ij)
               vsnon(i,j,nr) = vsnon(i,j,nr) &
                             + fvol(ij)*vsrdgn(ij)*fsnowrdg
            enddo

      !-----------------------------------------------------------------
      ! Transfer ice energy to category nr
      !-----------------------------------------------------------------
            do k = 1, nilyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
               do ij = 1, iridge
                  i = indxii(ij)
                  j = indxjj(ij)
                  eicen(i,j,ilyr1(nr)+k-1) = eicen(i,j,ilyr1(nr)+k-1) &
                                           + fvol(ij)*eirdgn(ij,k)
               enddo            ! ij
            enddo               ! k

      !-----------------------------------------------------------------
      ! Transfer snow energy to category nr
      !-----------------------------------------------------------------
            do k = 1, nslyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
               do ij = 1, iridge
                  i = indxii(ij)
                  j = indxjj(ij)
                  esnon(i,j,slyr1(nr)+k-1) = esnon(i,j,slyr1(nr)+k-1) &
                                    + fvol(ij)*esrdgn(ij,k)*fsnowrdg
               enddo            ! ij
            enddo               ! k

      !-----------------------------------------------------------------
      ! Transfer area-weighted and volume-weighted tracers to category nr.
      ! Note: The global sum aicen*trcrn of ice area tracers 
      !       (trcr_depend = 0) is not conserved by ridging.
      !       However, ridging conserves the global sum of volume
      !       tracers (trcr_depend = 1 or 2).
      !-----------------------------------------------------------------

            do it = 1, ntrcr
               if (trcr_depend(it) == 0) then  ! ice area tracer
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
                  do ij = 1, iridge
                     i = indxii(ij)
                     j = indxjj(ij)
                     m = indxij(ij)
                     atrcrn(m,it,nr) = atrcrn(m,it,nr) &
                                + farea(ij)*ardg2n(ij)*trcrn(i,j,it,n)

                  enddo
               elseif (trcr_depend(it) == 1) then ! ice volume tracer
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
                  do ij = 1, iridge
                     i = indxii(ij)
                     j = indxjj(ij)
                     m = indxij(ij)
                     atrcrn(m,it,nr) = atrcrn(m,it,nr) &
                                 + fvol(ij)*virdgn(ij)*trcrn(i,j,it,n)

                  enddo
               elseif (trcr_depend(it) == 2) then ! snow volume tracer
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
                  do ij = 1, iridge
                     i = indxii(ij)
                     j = indxjj(ij)
                     m = indxij(ij)
                     atrcrn(m,it,nr) = atrcrn(m,it,nr) &
                        + fvol(ij)*vsrdgn(ij)*fsnowrdg*trcrn(i,j,it,n)

                  enddo
               endif            ! trcr_depend
            enddo               ! ntrcr

         enddo                  ! nr (new ridges)
      enddo                     ! n (ridging categories)


      !-----------------------------------------------------------------
      ! Compute new tracers
      !-----------------------------------------------------------------

      do n = 1, ncat
         call compute_tracers (nx_block,        ny_block,       &
                               icells,          indxi,   indxj, &
                               ntrcr,           trcr_depend,    &
                               atrcrn(:,:,n),   aicen(:,:,  n), &
                               vicen (:,:,  n), vsnon(:,:,  n), &
                               trcrn(:,:,:,n))
      enddo

      end subroutine ridge_shift

!=======================================================================
!BOP
!
! !ROUTINE: ridge_check - check for ice area > 1
!
! !DESCRIPTION: Make sure ice area <=1.  If not, prepare to repeat ridging.
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!
! !INTERFACE:
!
      subroutine ridge_check (nx_block,  ny_block,        &
                              icells,    indxi,    indxj, &
                              dt,                         &
                              asum,      closing_net,     &
                              divu_adv,  opning,          &
                              iterate_ridging)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells                ! number of cells with ice present

      integer (kind=int_kind), dimension (nx_block*ny_block), &
         intent(in) :: &
         indxi, indxj     ! compressed indices for cells with ice


      real (kind=dbl_kind), intent(in) :: &
         dt               ! time step (s)

      real (kind=dbl_kind), dimension (icells), intent(in) :: &
         asum             ! sum of ice and open water area

      real (kind=dbl_kind), dimension (icells), &
         intent(inout) :: &
         closing_net, & ! net rate at which area is removed    (1/s)
         divu_adv   , & ! divu as implied by transport scheme  (1/s)
         opning         ! rate of opening due to divergence/shear

      logical (kind=log_kind), intent(out) :: &
         iterate_ridging      ! if true, repeat the ridging
!
!EOP
!
      integer (kind=int_kind) :: &
         ij               ! horizontal index, combines i and j loops

      iterate_ridging = .false.

      do ij = 1, icells
         if (abs(asum(ij) - c1) < puny) then
            closing_net(ij) = c0
            opning     (ij) = c0
         else
            iterate_ridging = .true.
            divu_adv(ij) = (c1 - asum(ij)) / dt
            closing_net(ij) = max(c0, -divu_adv(ij))
            opning(ij) = max(c0, divu_adv(ij))
         endif
      enddo

      end subroutine ridge_check

!=======================================================================
!BOP
!
! !ROUTINE: ice_strength - compute ice strength
!
! !DESCRIPTION:
!
! Compute the strength of the ice pack, defined as the energy (J m-2)
! dissipated per unit area removed from the ice pack under compression,
! and assumed proportional to the change in potential energy caused
! by ridging.
!
! See Rothrock (1975) and Hibler (1980).
!
! For simpler strength parameterization, see this reference:
! Hibler, W. D. III, 1979: A dynamic-thermodynamic sea ice model,
!  J. Phys. Oceanog., 9, 817-846.
!
! !REVISION HISTORY:
!
! authors: William H. Lipscomb, LANL
!          Elizabeth C. Hunke, LANL
!
! !INTERFACE:
!
      subroutine ice_strength (nx_block, ny_block, &
                               ilo, ihi, jlo, jhi, &
                               icells,             &
                               indxi,    indxj,    &
                               aice,     vice,     &
                               aice0,    aicen,    &
                               vicen,    strength)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         ilo,ihi,jlo,jhi       ! beg and end of physical domain

      integer (kind=int_kind), intent(in) :: &
         icells       ! no. of cells where icetmask = 1

      integer (kind=int_kind), dimension (nx_block*ny_block), &
         intent(in) :: &
         indxi   , & ! compressed index in i-direction
         indxj       ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         aice   , & ! concentration of ice
         vice   , & ! volume per unit area of ice  (m)
         aice0      ! concentration of open water

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat), &
         intent(in) :: &
         aicen  , & ! concentration of ice
         vicen      ! volume per unit area of ice  (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(out) :: &
         strength   ! ice strength (N/m)
!
!EOP
!
! LOCAL VARIABLES
!
      real (kind=dbl_kind), dimension (icells) :: &
         asum         , & ! sum of ice and open water area
         aksum            ! ratio of area removed to area ridged

      real (kind=dbl_kind), dimension (icells,0:ncat) :: &
         apartic          ! participation function; fraction of ridging
                          ! and closing associated w/ category n

      real (kind=dbl_kind), dimension (icells,ncat) :: &
         hrmin        , & ! minimum ridge thickness
         hrmax        , & ! maximum ridge thickness (krdg_redist = 0)
         hrexp        , & ! ridge e-folding thickness (krdg_redist = 1) 
         krdg             ! mean ridge thickness/thickness of ridging ice

      integer (kind=int_kind) :: &
         i,j             , & ! horizontal indices
         n               , & ! thickness category index
         ij                  ! horizontal index, combines i and j loops

      real (kind=dbl_kind) :: &
         hi              , & ! ice thickness (m)
         h2rdg           , & ! mean value of h^2 for new ridge
         dh2rdg              ! change in mean value of h^2 per unit area
                             ! consumed by ridging 

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      strength(:,:) = c0

      if (kstrength == 1) then  ! Rothrock '75 formulation

      !-----------------------------------------------------------------
      ! Compute thickness distribution of ridging and ridged ice.
      !-----------------------------------------------------------------

         call asum_ridging (nx_block, ny_block,      &
                            icells,   indxi,  indxj, &
                            aicen,    aice0,         &
                            asum)

         call ridge_itd (nx_block,    ny_block,      &
                         icells,      indxi,  indxj, &
                         aicen,       vicen,         &
                         aice0,                      &
                         aksum,       apartic,       &
                         hrmin,       hrmax,         &
                         hrexp,       krdg)

      !-----------------------------------------------------------------
      ! Compute ice strength based on change in potential energy,
      ! as in Rothrock (1975)
      !-----------------------------------------------------------------

         if (krdg_redist==0) then ! Hibler 1980 formulation

            do n = 1, ncat
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)
                  if (aicen(i,j,n) > puny .and. apartic(ij,n) > c0)then
                     hi = vicen(i,j,n) / aicen(i,j,n)
                     h2rdg = p333 * (hrmax(ij,n)**3 - hrmin(ij,n)**3)  &
                                  / (hrmax(ij,n) - hrmin(ij,n)) 
                     dh2rdg = -hi*hi + h2rdg/krdg(ij,n)
                     strength(i,j) = strength(i,j) &
                                   + apartic(ij,n) * dh2rdg
                  endif         ! aicen > puny
               enddo            ! ij
            enddo               ! n

         elseif (krdg_redist==1) then ! exponential formulation

            do n = 1, ncat
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)
                  if (aicen(i,j,n) > puny .and. apartic(ij,n) > c0)then
                     hi = vicen(i,j,n) / aicen(i,j,n)
                     h2rdg =    hrmin(ij,n)*hrmin(ij,n) &
                           + c2*hrmin(ij,n)*hrexp(ij,n) &
                           + c2*hrexp(ij,n)*hrexp(ij,n)
                     dh2rdg = -hi*hi + h2rdg/krdg(ij,n)
                     strength(i,j) = strength(i,j) &
                                   + apartic(ij,n) * dh2rdg
                  endif
               enddo            ! ij
            enddo               ! n

         endif                  ! krdg_redist

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            strength(i,j) = Cf * Cp * strength(i,j) / aksum(ij)
                          ! Cp = (g/2)*(rhow-rhoi)*(rhoi/rhow)
                          ! Cf accounts for frictional dissipation
         enddo            ! ij

      else                      ! kstrength /= 1:  Hibler (1979) form

      !-----------------------------------------------------------------
      ! Compute ice strength as in Hibler (1979)
      !-----------------------------------------------------------------

         do j = jlo, jhi
         do i = ilo, ihi
            strength(i,j) = Pstar*vice(i,j)*exp(-Cstar*(c1-aice(i,j)))
         enddo                  ! j
         enddo                  ! i

      endif                     ! kstrength

      end subroutine ice_strength

!=======================================================================

      end module ice_mechred

!=======================================================================
