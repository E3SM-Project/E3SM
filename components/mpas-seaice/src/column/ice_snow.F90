!  SVN:$Id: ice_snow.F90 972 2015-04-15 19:44:20Z njeffery $
!=======================================================================
!
! authors Elizabeth Hunke, LANL
!         Nicole Jeffery, LANL

      module ice_snow

      use ice_kinds_mod
      use ice_constants_colpkg, only: puny, c0, c1, c10, rhos, Lfresh, &
                                      rhow, rhoi, rhofresh, snwlvlfac, &
                                      rhosmin
      use ice_warnings, only: add_warning

      implicit none
      save

      private
      public :: snow_effective_density, update_snow_radius, snow_redist,&
                drain_snow

      real (kind=dbl_kind), parameter, public :: &
         S_r  = 0.033_dbl_kind, & ! irreducible saturation (Anderson 1976)
         S_wet= 0.422_dbl_kind  ! (um^3/s) wet metamorphism parameters

!=======================================================================

      contains

!=======================================================================

! Compute effective density of snow layers from ice, liquid water mass

      subroutine snow_effective_density(nslyr,     ncat,     &
                                        vsnon,     vsno,     &
                                        smice,     smliq,    &
                                        rhosnew,             &
                                        rhos_effn, rhos_eff, &
                                        rhos_cmpn, rhos_cmp)

      integer (kind=int_kind), intent(in) :: &
         nslyr, & ! number of snow layers
         ncat     ! number of thickness categories

      real (kind=dbl_kind), dimension(:), intent(in) :: &
         vsnon    ! snow volume (m)

      real (kind=dbl_kind), intent(in) :: &
         vsno , & ! total snow volume (m)
         rhosnew  ! new snow density (kg/m^3)

      real (kind=dbl_kind), dimension(:,:), &
         intent(inout) :: &
         smice    , & ! mass of ice in snow (kg/m^3)
         smliq    , & ! mass of liquid in snow (kg/m^3)
         rhos_effn, & ! effective snow density: content (kg/m^3)
         rhos_cmpn    ! effective snow density: compaction (kg/m^3)

      real (kind=dbl_kind), intent(inout) :: &
         rhos_eff , & ! mean effective snow density: content (kg/m^3)
         rhos_cmp     ! mean effective snow density: compaction (kg/m^3)

      integer (kind=int_kind) :: &
         k    , & ! snow layer index
         n    , & ! ice thickness category index
         cnt      ! counter for snow presence

      rhos_eff = c0
      rhos_cmp = c0

      !-----------------------------------------------------------------
      ! Initialize effective snow density (compaction) for new snow
      !-----------------------------------------------------------------

      do n = 1, ncat
            do k = 1, nslyr
               if (rhos_cmpn(k,n) < rhosmin) rhos_cmpn(k,n) = rhosnew
            enddo
      enddo

      !-----------------------------------------------------------------
      ! Compute average effective density of snow
      !-----------------------------------------------------------------

      if (vsno > puny) then

         do n = 1, ncat
            if (vsnon(n) > c0) then
               do k = 1, nslyr
                  rhos_effn(k,n) = rhos_effn(k,n) + smice(k,n) + smliq(k,n)
                  rhos_eff       = rhos_eff + vsnon(n)*rhos_effn(k,n)
                  rhos_cmp       = rhos_cmp + vsnon(n)*rhos_cmpn(k,n)
               enddo
            endif
         enddo
         rhos_eff = rhos_eff/(vsno*real(nslyr,kind=dbl_kind))
         rhos_cmp = rhos_cmp/(vsno*real(nslyr,kind=dbl_kind))

      endif ! vsno

      end subroutine snow_effective_density

!=======================================================================

! Snow redistribution by wind, based on O. Lecomte Ph.D. (2014).
! Namelist option snwredist = 'ITDsd':
! Snow in suspension depends on wind speed, density and the standard
! deviation of the ice thickness distribution. Snow is redistributed
! among ice categories proportionally to the category areas.
! Namelist option snwredist = 'ITDrdg':
! As above, but use the standard deviation of the level and ridged
! ice thickness distribution for snow in suspension, and redistribute
! based on ridged ice area.

! convention:
! volume, mass and energy include factor of ain
! thickness does not

      subroutine snow_redist(dt, nslyr, ncat, wind, ain, vin, vsn, zqsn, &
         snwredist, alvl, vlvl, fresh, fhocn, fsloss, rhos_cmpn, &
         fsnow, rhosmax, windmin, drhosdwind, l_stop, stop_label)

      use ice_therm_vertical, only: adjust_enthalpy

      integer (kind=int_kind), intent(in) :: &
         nslyr     , & ! number of snow layers
         ncat          ! number of thickness categories

      real (kind=dbl_kind), intent(in) :: &
         dt        , & ! time step (s)
         wind      , & ! wind speed (m/s)
         fsnow     , & ! snowfall rate (kg m-2 s-1)
         rhosmax   , & ! maximum snow density (kg/m^3)
         windmin   , & ! minimum wind speed to compact snow (m/s)
         drhosdwind    ! wind compaction factor (kg s/m^4)

      real (kind=dbl_kind), dimension(:), intent(in) :: &
         ain       , & ! ice area fraction
         vin       , & ! ice volume (m)
         alvl      , & ! level ice area tracer
         vlvl          ! level ice volume tracer

      real (kind=dbl_kind), intent(inout) :: &
         fresh     , & ! fresh water flux to ocean (kg/m^2/s)
         fhocn     , & ! net heat flux to ocean (W/m^2)
         fsloss        ! snow loss to leads (kg/m^2/s)

      real (kind=dbl_kind), dimension(:), intent(inout) :: &
         vsn           ! snow volume (m)

      real (kind=dbl_kind), dimension(:,:), intent(inout) :: &
         zqsn      , & ! snow enthalpy (J/m^3)
         rhos_cmpn     ! effective snow density: compaction (kg/m^3)

      character(len=char_len), intent(in) :: &
         snwredist                ! type of snow redistribution

      logical (kind=log_kind), intent(out) :: &
         l_stop          ! if true, print diagnostics and abort on return

      character (len=*), intent(out) :: stop_label

      ! local variables

      integer (kind=int_kind) :: &
         n         , & ! category index
         k             ! layer index

      integer (kind=int_kind), dimension(ncat) :: &
         klyr          ! layer index

      real (kind=dbl_kind), parameter :: &
         refsd   = c1            , & ! standard deviation reference
         gamma   = 1.e-5_dbl_kind    ! tuning coefficient

      real (kind=dbl_kind) :: &
         Vseas     , & ! critical seasonal wind speed (m/s)
         ITDsd     , & ! standard deviation of ITD
         flost     , & ! fraction of snow lost in leads
         alost     , & ! effective lead area for snow lost in leads
         suma      , & ! sum of ice area over categories
         sumv      , & ! sum of ice volume over categories (m)
         summ      , & ! sum of snow mass over categories (kg/m^2)
         sumq      , & ! sum of snow enthalpy over categories (kg/m^2)
         msusp     , & ! potential mass of snow in suspension (kg/m^2)
         msnw_susp , & ! mass of snow in suspension (kg/m^2)
         esnw_susp , & ! energy of snow in suspension (J/m^2)
         asnw_lvl  , & ! mass of snow redeposited on level ice (kg/m^2)
         e_redeptmp, & ! redeposited energy (J/m^2)
         dhsn      , & ! change in snow depth (m)
         dmp       , & ! mass difference in previous layer (kg/m^2)
         hslyr     , & ! snow layer thickness (m)
         hslab     , & ! new snow thickness (m)
         drhos     , & ! change in snow density due to compaction (kg/m^3)
         mlost     , & ! mass of suspended snow lost in leads (kg/m^2)
         elost     , & ! energy of suspended snow lost in leads (J/m^2)
         de        , & ! change in energy (J/m^2)
         al, ar    , & ! areas of level and ridged ice
         hlvl, hrdg, & ! thicknesses of level and ridged ice
         tmp1, tmp2, &          ! temporary values
         tmp3, tmp4, &          ! temporary values
         tmp5      , &          ! temporary values
         work          ! temporary value

      real (kind=dbl_kind), dimension(ncat) :: &
         sfac      , & ! temporary for snwlvlfac
         ardg      , & ! ridged ice area tracer
         m_erosion , & ! eroded mass (kg/m^2)
         e_erosion , & ! eroded energy (J/m^2)
         m_redep   , & ! redeposited mass (kg/m^2)
         e_redep   , & ! redeposited energy (J/m^2)
         vsn_init  , & ! initial volume (m)
         esn_init  , & ! initial energy (J/m^2)
         esn_final , & ! final energy (J/m^2)
         atmp      , & ! temporary variable for ain, for debugging convenience
         hin       , & ! ice thickness (m)
         hsn       , & ! snow depth (m)
         hsn_new       ! new snow depth (m)

      real (kind=dbl_kind), dimension (nslyr) :: &
         dzs             ! snow layer thickness after redistribution (m)

      real (kind=dbl_kind), dimension (nslyr+1) :: &
         zs1         , & ! depth of snow layer boundaries (m)
         zs2             ! adjusted depths, with equal hslyr (m)

      character(len=char_len_long) :: &
         warning  

      !-----------------------------------------------------------------
      ! Conservation checks
      !-----------------------------------------------------------------

      l_stop = .false.
      stop_label = ''
      tmp1 = c0
      tmp3 = c0
      do n = 1, ncat
         ! mass conservation check
         tmp1 = tmp1 + vsn(n)
         vsn_init(n) = vsn(n)
         esn_init(n) = c0
         ! energy conservation check
         do k = 1, nslyr
            tmp3 = tmp3 + vsn(n)*zqsn(k,n)/nslyr
            esn_init(n) = esn_init(n) + vsn(n)*zqsn(k,n)/nslyr
         enddo
      enddo

      !-----------------------------------------------------------------
      ! category thickness and sums
      !-----------------------------------------------------------------

      hin(:) = c0
      hsn(:) = c0
      suma = c0
      sumv = c0
      do n = 1, ncat
         atmp(n) = ain(n)
         if (atmp(n) > puny) then
            hin(n) = vin(n)/atmp(n)
            hsn(n) = vsn(n)/atmp(n)
         endif
         hsn_new(n) = hsn(n)
         suma = suma + atmp(n)
         sumv = sumv + vin(n)
         ! maintain positive definite enthalpy
         do k = 1, nslyr
            zqsn(k,n) = min(zqsn(k,n) + Lfresh*rhos, c0)
         enddo
      enddo ! ncat

      !-----------------------------------------------------------------
      ! standard deviation of ice thickness distribution
      !-----------------------------------------------------------------

      work = c0
      asnw_lvl = c0
      if (trim(snwredist) == 'ITDrdg') then  ! use level and ridged ice
         do n = 1, ncat
            ardg(n) = c1 - alvl(n) ! ridged ice tracer
            al = alvl(n) * atmp(n) ! level
            ar = ardg(n) * atmp(n) ! ridged
            hlvl = c0
            hrdg = c0
            if (al > puny) hlvl = vin(n)*vlvl(n)/al
            if (ar > puny) hrdg = vin(n)*(c1-vlvl(n))/ar
            work = work + al*(hlvl - sumv)**2 + ar*(hrdg - sumv)**2

            ! for redeposition of snow on level ice
            sfac(n) = snwlvlfac
            if (ardg(n) > c0) sfac(n) = min(snwlvlfac, alvl(n)/ardg(n))
            asnw_lvl = asnw_lvl + al - sfac(n)*ar
         enddo
         asnw_lvl = asnw_lvl/suma
      else ! snwredist = 'ITDsd'             ! use standard ITD
         do n = 1, ncat
            work = work + atmp(n)*(hin(n) - sumv)**2
         enddo
      endif
      ITDsd = sqrt(work)

      !-----------------------------------------------------------------
      ! fraction of suspended snow lost in leads
      !-----------------------------------------------------------------

      flost = (c1 - suma) * exp(-ITDsd/refsd)
!echmod      flost = c0
      alost =  c1 - suma  * (c1-flost)

      !-----------------------------------------------------------------
      ! suspended snow
      !-----------------------------------------------------------------

      msusp = c0
      do n = 1, ncat
         ! critical seasonal wind speed needed to compact snow to density rhos
         Vseas = (rhos_cmpn(1,n) - 44.6_dbl_kind)/174.0_dbl_kind ! use top layer
         Vseas = max(Vseas, c0)
         ! maximum mass per unit area of snow in suspension (kg/m^2)
         if (ITDsd > puny) &
            msusp = msusp + atmp(n)*gamma*dt*max(wind-Vseas,c0) &
                  * (rhosmax-rhos_cmpn(1,n))/(rhosmax*ITDsd)
      enddo

      !-----------------------------------------------------------------
      ! erosion
      !-----------------------------------------------------------------

      msnw_susp = c0
      esnw_susp = c0
      klyr(:) = 1
      do n = 1, ncat
         m_erosion(n) = c0                             ! mass
         e_erosion(n) = c0                             ! energy
         if (atmp(n) > puny) then
            m_erosion(n) = min(msusp, rhos*vsn(n))
            if (m_erosion(n) > puny) then
               summ = c0
               dmp = m_erosion(n)
               do k = 1, nslyr
                  if (dmp > c0) then
                     dhsn = min(hsn(n)/nslyr, dmp/(rhos*atmp(n)))
                     msnw_susp  = msnw_susp  + dhsn*rhos*atmp(n) ! total mass in suspension
                     hsn_new(n) = hsn_new(n) - dhsn
                     e_erosion(n) = e_erosion(n) + dhsn*zqsn(k,n)*atmp(n)
                     klyr(n) = k                        ! number of affected layers
                     summ = summ + rhos*vsn(n)/nslyr    ! mass, partial sum
                     dmp = max(m_erosion(n) - summ, c0)
                  endif ! dmp
               enddo
               esnw_susp = esnw_susp + e_erosion(n)     ! total energy in suspension
            endif
         endif
      enddo

      !-----------------------------------------------------------------
      ! redeposition
      !-----------------------------------------------------------------

      do n = 1, ncat
         if (trim(snwredist) == 'ITDrdg') then  ! use level and ridged ice
            work = atmp(n)*(c1-flost)*(ardg(n)*(c1+sfac(n)) + asnw_lvl)
         else                                   ! use standard ITD
            work = atmp(n)*(c1-flost)
         endif
         m_redep(n) = msnw_susp*work    ! mass
         e_redep(n) = c0
         e_redeptmp = esnw_susp*work    ! energy

         ! change in snow depth
         dhsn = c0
         if (atmp(n) > puny) then
            dhsn = m_redep(n) / (rhos*atmp(n))

            if (abs(dhsn) > c0) then

               e_redep(n) = e_redeptmp
               vsn(n) = (hsn_new(n)+dhsn)*atmp(n)

               ! change in snow energy
               de = e_redeptmp / klyr(n)
               ! spread among affected layers
               sumq = c0
               do k = 1, klyr(n)
                  zqsn(k,n) = (atmp(n)*hsn_new(n)*zqsn(k,n) + de) &
                            / (vsn(n))  ! factor of nslyr cancels out

                  if (zqsn(k,n) > c0) then
                     sumq = sumq + zqsn(k,n)
                     zqsn(k,n) = c0
                  endif

               enddo ! klyr
               zqsn(klyr(n),n) = min(zqsn(klyr(n),n) + sumq, c0) ! may lose energy here

      !-----------------------------------------------------------------
      ! Conserving energy, compute the enthalpy of the new equal layers
      !-----------------------------------------------------------------

               if (nslyr > 1) then

                  dzs(:) = hsn(n) / real(nslyr,kind=dbl_kind) ! old layer thickness
                  do k = 1, klyr(n)
                     dzs(k) = dzs(k) + dhsn / klyr(n)         ! old layer thickness (updated)
                  enddo
                  hsn_new(n) = hsn_new(n) + dhsn
                  hslyr  = hsn_new(n) / real(nslyr,kind=dbl_kind) ! new layer thickness

                  zs1(1) = c0
                  zs1(1+nslyr) = hsn_new(n)

                  zs2(1) = c0
                  zs2(1+nslyr) = hsn_new(n)

                  do k = 1, nslyr-1
                     zs1(k+1) = zs1(k) + dzs(k) ! old layer depths (unequal thickness)
                     zs2(k+1) = zs2(k) + hslyr  ! new layer depths (equal thickness)
                  enddo

                  call adjust_enthalpy (nslyr,                &
                                        zs1(:),   zs2(:),     &
                                        hslyr,    hsn_new(n), &
                                        zqsn(:,n))
               endif   ! nslyr > 1
            endif      ! |dhsn| > puny
         endif         ! ain > puny

         ! maintain positive definite enthalpy
         do k = 1, nslyr
            zqsn(k,n) = zqsn(k,n) - Lfresh*rhos
         enddo
      enddo ! ncat

      !-----------------------------------------------------------------
      ! mass of suspended snow lost in leads
      !-----------------------------------------------------------------
      mlost = msnw_susp*alost
      fsloss = fsloss + mlost / dt

      !-----------------------------------------------------------------
      ! mass conservation check
      !-----------------------------------------------------------------

      tmp2 = c0
      do n = 1, ncat
         tmp2 = tmp2 + vsn(n)
      enddo

      if (tmp2 > tmp1) then  ! correct roundoff error
         vsn(:) = vsn(:) * tmp1/tmp2
         tmp2 = c0
         do n = 1, ncat
            tmp2 = tmp2 + vsn(n)
         enddo
      endif

      if (tmp2 < tmp1) fresh = fresh + rhos*(tmp1-tmp2)/dt

      tmp2 = tmp2 + (mlost/rhos)

      if (abs(tmp1-tmp2) > puny) then
         write(warning,*)'mass conservation error in snow_redist', tmp1, tmp2
         call add_warning(warning)
         write(warning,*)'klyr',klyr
         call add_warning(warning)
         write(warning,*)'ain',atmp(:)
         call add_warning(warning)
         write(warning,*)'vsn final',vsn(:)
         call add_warning(warning)
         write(warning,*)'vsn init',vsn_init(:)
         call add_warning(warning)
         write(warning,*)'rhos*vsn init',rhos*vsn_init(:)
         call add_warning(warning)
         write(warning,*)'m_erosion',m_erosion(:)
         call add_warning(warning)
         write(warning,*)'m_redep',m_redep(:)
         call add_warning(warning)
         write(warning,*)'mlost',mlost
         call add_warning(warning)
         write(warning,*)'v_erosion',m_erosion(:)/rhos
         call add_warning(warning)
         write(warning,*)'v_redep',m_redep(:)/rhos
         call add_warning(warning)
         write(warning,*)'v lost',mlost/rhos
         call add_warning(warning)
         write(warning,*)'hsn',hsn(:)
         call add_warning(warning)
         write(warning,*)'hsn_new',hsn_new(:)
         call add_warning(warning)
         write(warning,*)'vsn_new',hsn_new(:)*atmp(:)
         call add_warning(warning)
         write(warning,*)'lost',suma,flost,alost,msnw_susp
         call add_warning(warning)
         stop_label = 'snow redistribution mass conservation error'
         l_stop = .true.
      endif

      !-----------------------------------------------------------------
      ! energy conservation check
      !-----------------------------------------------------------------

      tmp4 = c0
      tmp5 = c0
      esn_final(:) = c0
      do n = 1, ncat
         do k = 1, nslyr
            tmp4 = tmp4 + vsn(n)*zqsn(k,n)/nslyr
            esn_final(n) = esn_final(n) + vsn(n)*zqsn(k,n)/nslyr
         enddo
         tmp5 = tmp5 - e_erosion(n) + e_redep(n)
      enddo
      tmp5 = tmp5 + esnw_susp*alost

      !-----------------------------------------------------------------
      ! energy of suspended snow lost in leads
      !-----------------------------------------------------------------
      elost = tmp3 - tmp4
      fhocn = fhocn + elost / dt

      if (abs(tmp5) > nslyr*Lfresh*puny) then
         write(warning,*)'energy conservation error in snow_redist', tmp3, tmp4, tmp5
         call add_warning(warning)
         write(warning,*)'klyr',klyr
         call add_warning(warning)
         write(warning,*)'ain',atmp(:)
         call add_warning(warning)
         write(warning,*)'vsn final',vsn(:)
         call add_warning(warning)
         write(warning,*)'vsn init',vsn_init(:)
         call add_warning(warning)
         write(warning,*)'rhos*vsn init',rhos*vsn_init(:)
         call add_warning(warning)
         write(warning,*)'m_erosion',m_erosion(:)
         call add_warning(warning)
         write(warning,*)'m_redep',m_redep(:)
         call add_warning(warning)
         write(warning,*)'mlost',mlost
         call add_warning(warning)
         write(warning,*)'v_erosion',m_erosion(:)/rhos
         call add_warning(warning)
         write(warning,*)'v_redep',m_redep(:)/rhos
         call add_warning(warning)
         write(warning,*)'v lost',mlost/rhos
         call add_warning(warning)
         write(warning,*)'hsn',hsn(:)
         call add_warning(warning)
         write(warning,*)'hsn_new',hsn_new(:)
         call add_warning(warning)
         write(warning,*)'vsn_new',hsn_new(:)*atmp(:)
         call add_warning(warning)
         write(warning,*)'lost',suma,flost,alost,msnw_susp
         call add_warning(warning)
         write(warning,*)'tmp3(1)', (vsn(1)*zqsn(k,1)/nslyr,k=1,nslyr)
         call add_warning(warning)
         write(warning,*)'esn init',esn_init(:)
         call add_warning(warning)
         write(warning,*)'esn final',esn_final(:)
         call add_warning(warning)
         write(warning,*)'e_erosion',e_erosion(:)
         call add_warning(warning)
         write(warning,*)'e_redep',e_redep(:)
         call add_warning(warning)
         write(warning,*)'elost',elost,esnw_susp*alost,Lfresh*mlost
         call add_warning(warning)
         write(warning,*)'esnw_susp',esnw_susp
         call add_warning(warning)
         stop_label = 'snow redistribution energy conservation error'
         l_stop = .true.
      endif

      !-----------------------------------------------------------------
      ! wind compaction
      !-----------------------------------------------------------------

      do n = 1, ncat
         if (vsn(n) > puny) then
            ! compact freshly fallen or redistributed snow
            drhos = drhosdwind * max(wind - windmin, c0)
            hslab = c0
            if (fsnow > c0) &
               hslab = max(min(fsnow*dt/(rhos+drhos), hsn_new(n)-hsn(n)), c0)
            hslyr = hsn_new(n) / real(nslyr,kind=dbl_kind)
            do k = 1, nslyr
               work = hslab - hslyr * real(k-1,kind=dbl_kind)
               work = max(c0, min(hslyr, work))
               rhos_cmpn(k,n) = rhos_cmpn(k,n) + drhos*work/hslyr
               rhos_cmpn(k,n) = min(rhos_cmpn(k,n), rhosmax)
            enddo
         endif
      enddo

      end subroutine snow_redist

!=======================================================================

!  Snow grain metamorphism driver

      subroutine update_snow_radius (dt, ncat, nslyr, nilyr, rsnw, hin, &
                                     Tsfc, zTin,  &
                                     hsn, zqsn, smice, smliq, &
                                     rsnw_fall, rsnw_tmax, &
                                     snowage_tau, &
                                     snowage_kappa, &
                                     snowage_drdt0, &
                                     idx_T_max, &
                                     idx_Tgrd_max, &
                                     idx_rhos_max)

      integer (kind=int_kind), intent(in) :: &
         ncat,   & ! number of categories
         nslyr,  & ! number of snow layers
         nilyr, &  ! number of ice layers
         idx_T_max, & ! dimensions of snow parameter matrix
         idx_Tgrd_max, &
         idx_rhos_max

      real (kind=dbl_kind), intent(in) :: &
         dt          ! time step

      real (kind=dbl_kind), dimension(ncat), intent(in) :: &
         zTin        , & ! surface ice temperature (oC)
         Tsfc        , & ! surface temperature (oC)
         hin         , & ! ice thickness (m)
         hsn             ! snow thickness (m)

      real (kind=dbl_kind), dimension(nslyr,ncat), intent(in) :: &
         zqsn            ! enthalpy of snow (J m-3)

      real (kind=dbl_kind), dimension(nslyr,ncat), intent(inout) :: &
         rsnw

      real (kind=dbl_kind), dimension(nslyr,ncat), &
         intent(inout) :: &
         smice, & ! mass of ice in snow (kg/m^2)
         smliq    ! mass of liquid in snow (kg/m^2)

      real (kind=dbl_kind), intent(in) :: &
         rsnw_fall, & ! radius of newly fallen snow (10^-6 m)
         rsnw_tmax    ! maximum grain radius from dry metamorphism (10^-6 m)

      ! dry snow aging parameters
      real (kind=dbl_kind), dimension(idx_rhos_max,idx_Tgrd_max,idx_T_max), intent(in) :: &
         snowage_tau,   & ! (10^-6 m)
         snowage_kappa, & !
         snowage_drdt0    ! (10^-6 m/hr)

      ! local temporary variables

      integer (kind=int_kind) :: k, n

      real (kind=dbl_kind), dimension(nslyr) :: &
         drsnw_wet,    & ! wet metamorphism (10^-6 m)
         drsnw_dry       ! dry (temperature gradient) metamorphism (10^-6 m)

      !-----------------------------------------------------------------
      ! dry metamorphism
      !-----------------------------------------------------------------
       do n = 1, ncat

          if (hsn(n) > puny .and. hin(n) > puny) then

              drsnw_dry(:) = c0
              drsnw_wet(:) = c0

              call snow_dry_metamorph (nslyr, nilyr, dt, rsnw(:,n), drsnw_dry, zqsn(:,n), Tsfc(n), &
                                       zTin(n), hsn(n), hin(n), smice(:,n),smliq(:,n), rsnw_fall, &
                                       snowage_tau, snowage_kappa, snowage_drdt0, &
                                       idx_T_max, idx_Tgrd_max, idx_rhos_max)

      !-----------------------------------------------------------------
      ! wet metamorphism
      !-----------------------------------------------------------------


              do k = 1,nslyr
                    call snow_wet_metamorph  (dt, drsnw_wet(k), rsnw(k,n), smice(k,n),smliq(k,n))
                    rsnw(k,n) = min(rsnw_tmax, rsnw(k,n) + drsnw_dry(k) + drsnw_wet(k))
              enddo
           else
              do k = 1,nslyr
                rsnw(k,n) = max(rsnw_fall,min(rsnw_tmax, rsnw(k,n)))
                smice(k,n) = rhos
                smliq(k,n) = c0
              enddo

           endif
        enddo

      end subroutine update_snow_radius

!=======================================================================

!  Snow grain metamorphism

      subroutine snow_dry_metamorph (nslyr,nilyr, dt, rsnw, drsnw_dry, zqsn, &
                                     Tsfc, zTin1, hsn, hin, smice, smliq, rsnw_fall, &
                                     snowage_tau, snowage_kappa, snowage_drdt0, &
                                     idx_T_max, idx_Tgrd_max, idx_rhos_max)

      use ice_constants_colpkg, only: c0, rhos, Tffresh, Lfresh, cp_ice, p5, puny, c10
      use ice_colpkg_shared, only: idx_T_min, idx_Tgrd_min, idx_rhos_min

    ! Vapor redistribution: Method is to retrieve 3 best-bit parameters that
    ! depend on snow temperature, temperature gradient, and density,
    ! that are derived from the microphysical model described in:
    ! Flanner and Zender (2006), Linking snowpack microphysics and albedo
    ! evolution, J. Geophys. Res., 111, D12208, doi:10.1029/2005JD006834.
    ! The parametric equation has the form:
    ! dr/dt = drdt_0*(tau/(dr_fresh+tau))^(1/kappa), where:
    !   r is the effective radius,
    !   tau and kappa are best-fit parameters,
    !   drdt_0 is the initial rate of change of effective radius, and
    !   dr_fresh is the difference between the current and fresh snow states
    !  (r_current - r_fresh).

      integer (kind=int_kind), intent(in) :: &
         nslyr,  & ! number of snow layers
         nilyr, &  ! number of ice layers
         idx_T_max, & ! dimensions of snow parameter matrix
         idx_Tgrd_max, &
         idx_rhos_max

      real (kind=dbl_kind), intent(in) :: &
         dt                    ! time step (s)

      real (kind=dbl_kind), dimension(nslyr), &
         intent(in) :: &
         smice , & ! mass of ice in snow (kg/m^3)
         smliq , & ! mass of liquid in snow (kg/m^3)
         rsnw,   & ! snow grain radius (10^-6 m)
         zqsn      ! snow enthalpy  (J m-3)

      real (kind=dbl_kind), dimension(nslyr), &
         intent(inout) :: &
         drsnw_dry ! change due to snow aging (10^-6 m)

      real (kind=dbl_kind), intent(in) :: &
         Tsfc,   & ! surface temperature (oC)
         zTin1,  & ! top ice layer temperature (oC)
         hsn,    & ! snow thickness (m)
         hin,    & ! ice thickness (m)
         rsnw_fall

      ! dry snow aging parameters
      real (kind=dbl_kind), dimension(idx_rhos_max,idx_Tgrd_max,idx_T_max), intent(in) :: &
         snowage_tau,   & ! (10^-6 m)
         snowage_kappa, & !
         snowage_drdt0    ! (10^-6 m/hr)

      ! local temporary variables

      integer (kind=int_kind) :: k

      integer (kind=int_kind) :: &
          T_idx,    & ! temperature index
          Tgrd_idx, & ! temperature gradient index
          rhos_idx    ! density index

      real (kind=dbl_kind), dimension(nslyr):: &
         zrhos,   & ! snow density (kg/m^3)  ! for variable snow density
         zdTdz,   & ! temperature gradient (K/s)
         zTsn       ! snow temperature (oC)

      real (kind=dbl_kind) :: &
         bst_tau,   & ! snow aging parameter retrieved from lookup table [hour]
         bst_kappa, & ! snow aging parameter retrieved from lookup table [unitless]
         bst_drdt0, & ! snow aging parameter retrieved from lookup table [um hr-1]
         dr_fresh,  & ! change in snow radius from fresh (10^-6 m)
         dzs,       & ! snow layer thickness (m)
         dzi          ! ice layer thickness (m)

      character(len=char_len_long) :: &
         warning ! warning message

! Needed for variable snow density not currently modeled
! calculate density based on liquid and ice content of snow

         drsnw_dry(:) = c0
         zTsn(:) = c0
         zdTdz(:) = c0
         zrhos(:) = rhos

         dzs = hsn/real(nslyr,kind=dbl_kind)
         dzi = hin/real(nilyr,kind=dbl_kind)

         if (nslyr == 1) then
              zTsn(1)  =(Lfresh + zqsn(1)/rhos)/cp_ice
              zdTdz(1) =  min(c10*idx_Tgrd_max,abs((zTsn(1)*dzi + zTin1*dzs)/(dzs + dzi+puny)- Tsfc)/(hsn+puny))
         else
              zTsn(1)  =(Lfresh + zqsn(1)/rhos)/cp_ice
              do k = 2, nslyr
                 zTsn(k) = (Lfresh + zqsn(k)/rhos)/cp_ice
                 if (k == 2) then
                    zdTdz(k-1) = abs((zTsn(k-1)+zTsn(k))*p5 - Tsfc)/(dzs+puny)
                    zdTdz(k-1) = min(c10*idx_Tgrd_max,zdTdz(k-1))
                 else
                    zdTdz(k-1) = abs(zTsn(k-2)-zTsn(k))*p5/(dzs+puny)
                    zdTdz(k-1) = min(c10*idx_Tgrd_max,zdTdz(k-1))
                 endif
              enddo

              zdTdz(nslyr) = abs((zTsn(nslyr)*dzi + zTin1*dzs)/(dzs + dzi+puny)- &
                            (zTsn(nslyr) + zTsn(nslyr-1))*p5)/(dzs+puny)
              zdTdz(nslyr) = min(c10*idx_Tgrd_max,zdTdz(nslyr))
         endif

         ! best-fit parameters are read from a table
         !  11 temperatures from 225 to 273 K
         !  31 temperature gradients from 0 to 300 K/m
         !   8 snow densities from 0 to 350 kg/m3
         ! pointer snowage_tau, snowage_kappa, snowage_drdt0

         do k = 1, nslyr
          zrhos(k) = smice(k) + smliq(k)

          ! best-fit table indecies:
          T_idx    = nint(abs(zTsn(k)+ Tffresh - 223.15_dbl_kind) / 5.0_dbl_kind, kind=int_kind)
          Tgrd_idx = nint(zdTdz(k) / 10.0_dbl_kind, kind=int_kind)
          !rhos_idx = nint(zrhos(k)-50.0_dbl_kind) / 50.0_dbl_kind, kind=int_kind)   ! variable density
          rhos_idx = nint((rhos-50.0_dbl_kind) / 50.0_dbl_kind, kind=int_kind)        ! fixed density

          ! boundary check:
          T_idx = min(idx_T_max, max(1,T_idx+1))!min(idx_T_max, max(idx_T_min,T_idx))
          Tgrd_idx = min(idx_Tgrd_max, max(1,Tgrd_idx+1))!min(idx_Tgrd_max, max(idx_Tgrd_min,Tgrd_idx))
          rhos_idx = min(idx_rhos_max, max(1,rhos_idx+1)) !min(idx_rhos_max, max(idx_rhos_min,rhos_idx))

          bst_tau   = snowage_tau(rhos_idx,Tgrd_idx,T_idx)
          bst_kappa = snowage_kappa(rhos_idx,Tgrd_idx,T_idx)
          bst_drdt0 = snowage_drdt0(rhos_idx,Tgrd_idx,T_idx)

          ! change in snow effective radius, using best-fit parameters
          dr_fresh = max(c0,rsnw(k)-rsnw_fall)
          drsnw_dry(k) = (bst_drdt0*(bst_tau/(dr_fresh+bst_tau))**(1/bst_kappa))&
                     * (dt/3600.0_dbl_kind)
         enddo

      end subroutine snow_dry_metamorph

!=======================================================================

!  Snow grain metamorphism

      subroutine snow_wet_metamorph (dt, dr_wet, rsnw, smice, smliq)

    use ice_constants_colpkg, only: c0, c1, c4, pi, p1, c100
    !
    ! Liquid water redistribution: Apply the grain growth function from:
    !   Brun, E. (1989), Investigation of wet-snow metamorphism in respect of
    !   liquid-water content, Annals of Glaciology, 13, 22-26.
    !   There are two parameters that describe the grain growth rate as
    !   a function of snow liquid water content (LWC). The "LWC=0" parameter
    !   is zeroed here because we are accounting for dry snowing with a
    !   different representation
    !
      real (kind=dbl_kind), intent(in) :: &
         dt                    ! time step

      real (kind=dbl_kind), &
         intent(in) :: &
         rsnw , & ! snow grain radius (10^-6 m)
         smice, & ! snow ice density (kg/m^3)
         smliq    ! snow liquid density (kg/m^3)

      real (kind=dbl_kind), &
         intent(inout) :: &
         dr_wet

      real (kind=dbl_kind) :: &
         fliq  ! liquid mass fraction

       dr_wet = c0
       fliq = c1
       if (smice + smliq > c0 .and. rsnw > c0) then
         fliq = min(smliq/(smice + smliq),p1)*c100
         dr_wet = S_wet * fliq**3*dt/(c4*pi*rsnw**2)
       endif

      end subroutine snow_wet_metamorph

!=======================================================================

!  Conversions between ice mass, liquid water mass in snow

      subroutine drain_snow (dt, nslyr, vsnon,  aicen, &
                             smice, smliq, meltsliq, use_smliq_pnd)

      integer (kind=int_kind), intent(in) :: &
         nslyr    ! number of snow layers

      real (kind=dbl_kind), intent(in) :: &
         dt,     & ! time step
         vsnon,  & ! snow volume (m)
         aicen     ! aice area

      real (kind=dbl_kind), intent(inout) :: &
         meltsliq  ! total liquid content

      real (kind=dbl_kind), dimension(nslyr), &
         intent(in) :: &
         smice    ! mass of ice in snow (kg/m^2)

      real (kind=dbl_kind), dimension(nslyr), &
         intent(inout) :: &
         smliq    ! mass of liquid in snow (kg/m^2)

      logical (kind=log_kind), intent(in) :: &
         use_smliq_pnd   ! if true, use snow liquid tracer for ponds

      ! local temporary variables

      integer (kind=int_kind) ::  k

      real (kind=dbl_kind) :: &
        hslyr,  & ! snow layer thickness (m)
        hsn,    & ! snow thickness (m)
        meltsliq_tmp  ! temperary snow liquid content

      real (kind=dbl_kind), dimension(nslyr) :: &
         dlin    , & ! liquid into the layer from above (kg/m^2)
         dlout   , & ! liquid out of the layer (kg/m^2)
         phi_liq , & ! volumetric liquid fraction
         phi_ice , & ! volumetric ice fraction
         w_drain    ! flow between layers

      hsn = c0
      meltsliq_tmp = c0
      if (aicen > c0) hsn = vsnon/aicen
      if (hsn > puny) then
        dlin(:) = c0
        dlout(:) = c0
        hslyr    = hsn / real(nslyr,kind=dbl_kind)
        do k = 1,nslyr
            smliq(k)   = smliq(k)  + dlin(k) / hslyr   ! liquid in from above layer
            phi_ice(k) = min(c1, smice(k) / rhoi)
            phi_liq(k) = smliq(k)/rhofresh
            w_drain(k) = max(c0, (phi_liq(k) - S_r*(c1-phi_ice(k))) / dt * rhofresh * hslyr)
            dlout(k)   = w_drain(k) * dt
            smliq(k)   = smliq(k) - dlout(k)/ hslyr
            if (k < nslyr) then
                dlin(k+1) = dlout(k)
            else
                meltsliq_tmp = dlout(nslyr)
            endif
        enddo
      else
        meltsliq_tmp = meltsliq  ! computed in thickness_changes
      endif

      meltsliq = meltsliq
      if (use_smliq_pnd) meltsliq = meltsliq_tmp

      end subroutine drain_snow

!=======================================================================

      end module ice_snow

!=======================================================================
