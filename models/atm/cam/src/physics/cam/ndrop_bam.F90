module ndrop_bam

!---------------------------------------------------------------------------------
!
!   CAM Interface for droplet activation by bulk aerosols.
!
!---------------------------------------------------------------------------------

use shr_kind_mod,     only: r8=>shr_kind_r8
use spmd_utils,       only: masterproc
use ppgrid,           only: pcols, pver, pverp
use physconst,        only: gravit, rair, tmelt, cpair, rh2o, &
     r_universal, mwh2o, rhoh2o, latvap

use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_aer_props

use shr_spfn_mod,     only: erf => shr_spfn_erf, &
                            erfc => shr_spfn_erfc
use wv_saturation,    only: qsat
use cam_history,      only: addfld, add_default, phys_decomp, outfld 
use cam_logfile,      only: iulog
use abortutils,       only: endrun
use ref_pres,         only: top_lev=>trop_cloud_top_lev

implicit none
private
save

public :: ndrop_bam_init, ndrop_bam_run, ndrop_bam_ccn

! activate parameters

real(r8) :: pi     ! pi
integer,  parameter :: psat = 6           ! number of supersaturations to calc ccn concentration
real(r8), parameter :: supersat(psat) = & ! supersaturation (%) to determine ccn concentration
                       (/0.02_r8,0.05_r8,0.1_r8,0.2_r8,0.5_r8,1.0_r8/)
real(r8)              :: super(psat)
real(r8), allocatable :: ccnfact(:,:)

real(r8), allocatable :: alogsig(:)       ! natl log of geometric standard dev of aerosol
real(r8), allocatable :: exp45logsig(:)
real(r8), allocatable :: argfactor(:)
real(r8), allocatable :: amcube(:)        ! cube of dry mode radius (m)
real(r8), allocatable :: smcrit(:)        ! critical supersatuation for activation
real(r8), allocatable :: lnsm(:)          ! ln(smcrit)
real(r8), allocatable :: amcubefactor(:)  ! factors for calculating mode radius
real(r8), allocatable :: smcritfactor(:)  ! factors for calculating critical supersaturation
real(r8), allocatable :: f1(:), f2(:)     ! abdul-razzak functions of width

real(r8) :: aten
real(r8) :: third, sixth
real(r8) :: sq2, sqpi
real(r8) :: alogten, alog2, alog3, alogaten

real(r8) :: pref = 1013.25e2_r8 ! reference pressure (Pa)

! aerosol properties
character(len=20), allocatable :: aername(:)
real(r8), allocatable :: dryrad_aer(:)
real(r8), allocatable :: density_aer(:)
real(r8), allocatable :: hygro_aer(:)
real(r8), allocatable :: dispersion_aer(:)
real(r8), allocatable :: num_to_mass_aer(:)

integer :: naer_all      ! number of aerosols affecting climate
integer :: idxsul   = -1 ! index in aerosol list for sulfate

!===============================================================================
contains
!===============================================================================

subroutine ndrop_bam_init

   !----------------------------------------------------------------------- 
   ! 
   ! Initialize constants for droplet activation by bulk aerosols
   ! 
   !-----------------------------------------------------------------------

   integer  :: l, m, iaer
   real(r8) :: surften       ! surface tension of water w/respect to air (N/m)
   real(r8) :: arg
   !-------------------------------------------------------------------------------

   ! Access the physical properties of the bulk aerosols that are affecting the climate
   ! by using routines from the rad_constituents module.

   call rad_cnst_get_info(0, naero=naer_all)
   allocate( &
      aername(naer_all),        &
      dryrad_aer(naer_all),     &
      density_aer(naer_all),    &
      hygro_aer(naer_all),      &
      dispersion_aer(naer_all), &
      num_to_mass_aer(naer_all) )

   do iaer = 1, naer_all
      call rad_cnst_get_aer_props(0, iaer, &
         aername         = aername(iaer), &
         dryrad_aer      = dryrad_aer(iaer), &
         density_aer     = density_aer(iaer), &
         hygro_aer       = hygro_aer(iaer), &
         dispersion_aer  = dispersion_aer(iaer), &
         num_to_mass_aer = num_to_mass_aer(iaer) )

      ! Look for sulfate aerosol in this list (Bulk aerosol only)
      if (trim(aername(iaer)) == 'SULFATE') idxsul   = iaer

      ! aerosol number concentration
      call addfld(trim(aername(iaer))//'_m3', 'm-3', pver, 'A', 'aerosol number concentration', phys_decomp)

   end do

   if (masterproc) then
      write(iulog,*) 'ndrop_bam_init: iaer, name, dryrad, density, hygro, dispersion, num_to_mass'
      do iaer = 1, naer_all
         write(iulog,*) iaer, aername(iaer), dryrad_aer(iaer), density_aer(iaer), hygro_aer(iaer), &
            dispersion_aer(iaer), num_to_mass_aer(iaer)
      end do
      if (idxsul < 1) then
         write(iulog,*) 'ndrop_bam_init: SULFATE aerosol properties NOT FOUND'
      else
         write(iulog,*) 'ndrop_bam_init: SULFATE aerosol properties FOUND at index ',idxsul
      end if
   end if

   call addfld ('CCN1    ','#/cm3   ',pver, 'A','CCN concentration at S=0.02%',phys_decomp)
   call addfld ('CCN2    ','#/cm3   ',pver, 'A','CCN concentration at S=0.05%',phys_decomp)
   call addfld ('CCN3    ','#/cm3   ',pver, 'A','CCN concentration at S=0.1%',phys_decomp)
   call addfld ('CCN4    ','#/cm3   ',pver, 'A','CCN concentration at S=0.2%',phys_decomp)
   call addfld ('CCN5    ','#/cm3   ',pver, 'A','CCN concentration at S=0.5%',phys_decomp)
   call addfld ('CCN6    ','#/cm3   ',pver, 'A','CCN concentration at S=1.0%',phys_decomp)

   call add_default('CCN3', 1, ' ')

   ! set parameters for droplet activation, following abdul-razzak and ghan 2000, JGR
   
   third   = 1._r8/3._r8
   sixth   = 1._r8/6._r8
   sq2     = sqrt(2._r8)
   pi      = 4._r8*atan(1.0_r8)
   sqpi    = sqrt(pi)
   surften = 0.076_r8
   aten    = 2._r8*mwh2o*surften/(r_universal*tmelt*rhoh2o)
   alogaten= log(aten)
   alog2   = log(2._r8)
   alog3   = log(3._r8)
   super(:)= 0.01_r8*supersat(:)

   allocate(                  &
      alogsig(naer_all),      &
      exp45logsig(naer_all),  &
      argfactor(naer_all),    &
      f1(naer_all),           &
      f2(naer_all),           &
      amcubefactor(naer_all), &
      smcritfactor(naer_all), &
      amcube(naer_all),       &
      smcrit(naer_all),       &
      lnsm(naer_all),         &
      ccnfact(psat,naer_all)  )


   do m = 1, naer_all

      ! Skip aerosols that don't have a dispersion defined.
      if (dispersion_aer(m) == 0._r8) cycle
      
      alogsig(m)     = log(dispersion_aer(m))
      exp45logsig(m) = exp(4.5_r8*alogsig(m)*alogsig(m))
      argfactor(m)   = 2._r8/(3._r8*sqrt(2._r8)*alogsig(m))
      f1(m)          = 0.5_r8*exp(2.5_r8*alogsig(m)*alogsig(m))
      f2(m)          = 1._r8 + 0.25_r8*alogsig(m)
      amcubefactor(m)= 3._r8/(4._r8*pi*exp45logsig(m)*density_aer(m))
      smcritfactor(m)= 2._r8*aten*sqrt(aten/(27._r8*max(1.e-10_r8,hygro_aer(m))))
      amcube(m)      = amcubefactor(m)/num_to_mass_aer(m)

      if (hygro_aer(m) .gt. 1.e-10_r8) then
         smcrit(m) = smcritfactor(m)/sqrt(amcube(m))
      else
         smcrit(m) = 100._r8
      endif
      lnsm(m) = log(smcrit(m))

      do l = 1, psat
         arg = argfactor(m)*log(smcrit(m)/super(l))
         if (arg < 2) then
            if (arg < -2) then
               ccnfact(l,m) = 1.e-6_r8
            else
               ccnfact(l,m) = 1.e-6_r8*0.5_r8*erfc(arg)
            endif
         else
            ccnfact(l,m) = 0._r8
         endif
      enddo

   end do

end subroutine ndrop_bam_init

!===============================================================================

subroutine ndrop_bam_run( &
   wbar, tair, rhoair, na, pmode, &
   nmode, ma, nact)

   ! calculates number fraction of aerosols activated as CCN
   ! assumes an internal mixture within each of up to pmode multiple aerosol modes
   ! a gaussian spectrum of updrafts can be treated.

   !      mks units

   !      Abdul-Razzak and Ghan, A parameterization of aerosol activation.
   !      2. Multiple aerosol types. J. Geophys. Res., 105, 6837-6844.

   ! input
   integer,  intent(in) :: pmode         ! dimension of modes
   integer,  intent(in) :: nmode         ! number of aerosol modes
   real(r8), intent(in) :: wbar          ! grid cell mean vertical velocity (m/s)
   real(r8), intent(in) :: tair          ! air temperature (K)
   real(r8), intent(in) :: rhoair        ! air density (kg/m3)
   real(r8), intent(in) :: na(pmode)     ! aerosol number concentration (1/m3)
   real(r8), intent(in) :: ma(pmode)     ! aerosol mass concentration (kg/m3)

   ! output
   real(r8), intent(out) :: nact         ! number fraction of aerosols activated

   ! local variables
   integer :: maxmodes

   real(r8), allocatable :: volc(:) ! total aerosol volume  concentration (m3/m3)
   real(r8), allocatable :: eta(:)
   real(r8), allocatable :: smc(:)
   real(r8), allocatable :: etafactor2(:)
   real(r8), allocatable :: amcubeloc(:)
   real(r8), allocatable :: lnsmloc(:)

   real(r8) :: pres     ! pressure (Pa)
   real(r8) :: diff0
   real(r8) :: conduct0 ! thermal conductivity (Joule/m/sec/deg)
   real(r8) :: qs       ! water vapor saturation mixing ratio
   real(r8) :: dqsdt    ! change in qs with temperature
   real(r8) :: gloc     ! thermodynamic function (m2/s)
   real(r8) :: zeta
   real(r8) :: lnsmax   ! ln(smax)
   real(r8) :: alpha
   real(r8) :: gammaloc
   real(r8) :: beta
   real(r8) :: sqrtg
   real(r8) :: wnuc
   real(r8) :: alw
   real(r8) :: sqrtalw
   real(r8) :: smax
   real(r8) :: x
   real(r8) :: etafactor1
   real(r8) :: etafactor2max
   real(r8) :: es
   integer  :: m
   !-------------------------------------------------------------------------------

   maxmodes = naer_all
   allocate( &
      volc(maxmodes),       &
      eta(maxmodes),        &
      smc(maxmodes),        &
      etafactor2(maxmodes), &
      amcubeloc(maxmodes),  &
      lnsmloc(maxmodes)     )

   if (maxmodes < pmode) then
      write(iulog,*)'ndrop_bam_run: maxmodes,pmode=',maxmodes,pmode
      call endrun('ndrop_bam_run')
   endif

   nact = 0._r8

   if (nmode .eq. 1  .and.  na(1) .lt. 1.e-20_r8) return

   if (wbar .le. 0._r8) return

   pres     = rair*rhoair*tair
   diff0    = 0.211e-4_r8*(pref/pres)*(tair/tmelt)**1.94_r8
   conduct0 = (5.69_r8 + 0.017_r8*(tair-tmelt))*4.186e2_r8*1.e-5_r8 ! convert to J/m/s/deg
   call qsat(tair, pres, es, qs)
   dqsdt    = latvap/(rh2o*tair*tair)*qs
   alpha    = gravit*(latvap/(cpair*rh2o*tair*tair) - 1._r8/(rair*tair))
   gammaloc = (1 + latvap/cpair*dqsdt)/(rhoair*qs)
   ! growth coefficent Abdul-Razzak & Ghan 1998 eqn 16
   ! should depend on mean radius of mode to account for gas kinetic effects
   gloc     = 1._r8/(rhoh2o/(diff0*rhoair*qs)                               &
              + latvap*rhoh2o/(conduct0*tair)*(latvap/(rh2o*tair) - 1._r8))
   sqrtg    = sqrt(gloc)
   beta     = 4._r8*pi*rhoh2o*gloc*gammaloc
   etafactor2max = 1.e10_r8/(alpha*wbar)**1.5_r8 ! this should make eta big if na is very small.

   do m = 1, nmode
      ! skip aerosols with no dispersion, since they aren't meant to be CCN
      if (dispersion_aer(m) == 0._r8) then
         smc(m)=100._r8
         cycle
      end if
      ! internal mixture of aerosols
      volc(m) = ma(m)/(density_aer(m)) ! only if variable size dist
      if (volc(m) > 1.e-39_r8 .and. na(m) > 1.e-39_r8) then
         etafactor2(m) = 1._r8/(na(m)*beta*sqrtg)  !fixed or variable size dist
         ! number mode radius (m)
         amcubeloc(m) = (3._r8*volc(m)/(4._r8*pi*exp45logsig(m)*na(m)))  ! only if variable size dist
         smc(m) = smcrit(m) ! only for prescribed size dist

         if (hygro_aer(m) > 1.e-10_r8) then   ! loop only if variable size dist
            smc(m) = 2._r8*aten*sqrt(aten/(27._r8*hygro_aer(m)*amcubeloc(m))) 
         else
            smc(m) = 100._r8
         endif
      else
         smc(m) = 1._r8
         etafactor2(m) = etafactor2max ! this should make eta big if na is very small.
      end if
      lnsmloc(m) = log(smc(m)) ! only if variable size dist
   end do

   ! single  updraft
   wnuc    = wbar
   alw     = alpha*wnuc
   sqrtalw = sqrt(alw)
   zeta    = 2._r8*sqrtalw*aten/(3._r8*sqrtg)
   etafactor1 = 2._r8*alw*sqrtalw

   do m = 1, nmode
      ! skip aerosols with no dispersion, since they aren't meant to be CCN
      if (dispersion_aer(m) /= 0._r8) eta(m) = etafactor1*etafactor2(m)
   end do

   call maxsat(zeta, eta, nmode, smc, smax)
   lnsmax = log(smax)

   nact = 0._r8
   do m = 1, nmode
      ! skip aerosols with no dispersion, since they aren't meant to be CCN
      if (dispersion_aer(m) == 0._r8) cycle
      x = 2*(lnsmloc(m) - lnsmax)/(3*sq2*alogsig(m))
      nact = nact + 0.5_r8*(1._r8 - erf(x))*na(m)
   end do
   nact = nact/rhoair ! convert from #/m3 to #/kg

   deallocate( &
      volc,       &
      eta,        &
      smc,        &
      etafactor2, &
      amcubeloc,  &
      lnsmloc     )

end subroutine ndrop_bam_run

!===============================================================================

subroutine ndrop_bam_ccn(lchnk, ncol, maerosol, naer2)

   !-------------------------------------------------------------------------------
   !
   ! Write diagnostic bulk aerosol ccn concentration
   !
   !-------------------------------------------------------------------------------

   ! Input arguments
   integer,  intent(in) :: lchnk           ! chunk identifier
   integer,  intent(in) :: ncol            ! number of columns
   real(r8), intent(in) :: naer2(:,:,:)    ! aerosol number concentration (1/m3)
   real(r8), intent(in) :: maerosol(:,:,:) ! aerosol mass conc (kg/m3)

   ! Local variables
   integer :: i, k, l, m
   real(r8) :: arg

   character*8, parameter :: ccn_name(psat)=(/'CCN1','CCN2','CCN3','CCN4','CCN5','CCN6'/)
   real(r8) :: amcubesulfate(pcols)  ! cube of dry mode radius (m) of sulfate
   real(r8) :: smcritsulfate(pcols)  ! critical supersatuation for activation of sulfate
   real(r8) :: ccnfactsulfate
   real(r8) :: ccn(pcols,pver,psat)  ! number conc of aerosols activated at supersat
   !-------------------------------------------------------------------------------

   ccn(:ncol,:,:) = 0._r8

   do k = top_lev, pver

      do m = 1, naer_all

         if (m == idxsul) then
            ! Lohmann treatment for sulfate has variable size distribution
            do i = 1, ncol
               if (naer2(i,k,m) > 0._r8) then 
                  amcubesulfate(i) = amcubefactor(m)*maerosol(i,k,m)/(naer2(i,k,m))
                  smcritsulfate(i) = smcritfactor(m)/sqrt(amcubesulfate(i))
               else
                  smcritsulfate(i) = 1._r8
               end if
            end do
         end if

         do l = 1, psat

            if (m == idxsul) then
               ! This code is modifying ccnfact for sulfate only.
               do i = 1, ncol
                  arg = argfactor(m)*log(smcritsulfate(i)/super(l))
                  if (arg < 2) then
                     if (arg < -2) then
                        ccnfactsulfate = 1.0e-6_r8
                     else
                        ccnfactsulfate = 0.5e-6_r8*erfc(arg)
                     end if
                  else
                     ccnfactsulfate = 0.0_r8
                  end if
                  ccn(i,k,l) = ccn(i,k,l) + naer2(i,k,m)*ccnfactsulfate
               end do
            else
               ! Non-sulfate species use ccnfact computed by the init routine
               ccn(:ncol,k,l) = ccn(:ncol,k,l) + naer2(:ncol,k,m)*ccnfact(l,m)
            end if

         end do   ! supersaturation
      end do      ! bulk aerosol
   end do         ! level

   do l = 1, psat
      call outfld(ccn_name(l), ccn(1,1,l), pcols, lchnk)
   end do

   do l = 1, naer_all
      call outfld(trim(aername(l))//'_m3', naer2(:,:,l), pcols, lchnk)
   end do

end subroutine ndrop_bam_ccn

!===============================================================================

subroutine maxsat(zeta, eta, nmode, smc, smax)

   ! calculates maximum supersaturation for multiple
   ! competing aerosol modes.

   ! Abdul-Razzak and Ghan, A parameterization of aerosol activation.
   ! 2. Multiple aerosol types. J. Geophys. Res., 105, 6837-6844.

   real(r8), intent(in) :: zeta
   integer,  intent(in) :: nmode ! number of modes
   real(r8), intent(in) :: smc(:) ! critical supersaturation for number mode radius
   real(r8), intent(in) :: eta(:)

   real(r8), intent(out) :: smax ! maximum supersaturation

   integer :: m  ! mode index
   real(r8) :: sum, g1, g2

   do m=1,nmode
      if(zeta.gt.1.e5_r8*eta(m).or.smc(m)*smc(m).gt.1.e5_r8*eta(m))then
         ! weak forcing. essentially none activated
         smax=1.e-20_r8
      else
         ! significant activation of this mode. calc activation all modes.
         go to 1
      endif
   enddo

   return

1  continue

   sum=0
   do m=1,nmode
      if(eta(m).gt.1.e-20_r8)then
         g1=sqrt(zeta/eta(m))
         g1=g1*g1*g1
         g2=smc(m)/sqrt(eta(m)+3*zeta)
         g2=sqrt(g2)
         g2=g2*g2*g2
         sum=sum+(f1(m)*g1+f2(m)*g2)/(smc(m)*smc(m))
      else
         sum=1.e20_r8
      endif
   enddo
   
   smax=1._r8/sqrt(sum)
   
end subroutine maxsat

!===============================================================================

end module ndrop_bam
