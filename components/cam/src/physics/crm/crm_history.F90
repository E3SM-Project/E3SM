module crm_history
!---------------------------------------------------------------------------------------------------
! Purpose: Encapsulate setting up and writing history fields specific to the MMF
!---------------------------------------------------------------------------------------------------

use shr_kind_mod, only: r8 => shr_kind_r8

implicit none 

public :: crm_history_register
public :: crm_history_init
public :: crm_history_out

contains

!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
subroutine crm_history_register()
   use phys_control,        only: phys_getopts
   use cam_history_support, only: add_hist_coord
   use crmdims,             only: crm_nx, crm_ny, crm_nz, crm_nx_rad, crm_ny_rad
#ifdef ECPP
   use ecppvars,            only: NCLASS_CL,ncls_ecpp_in,NCLASS_PR
#endif
   !----------------------------------------------------------------------------
   ! local variables

   logical :: use_ECPP

   !----------------------------------------------------------------------------

   ! Add crm dimensions to cam history 
   call add_hist_coord('crm_nx',     crm_nx,     'CRM NX')
   call add_hist_coord('crm_ny',     crm_ny,     'CRM NY')
   call add_hist_coord('crm_nz',     crm_nz,     'CRM NZ')
   call add_hist_coord('crm_nx_rad', crm_nx_rad, 'Number of x columns for radiation')
   call add_hist_coord('crm_ny_rad', crm_ny_rad, 'Number of y columns for radiation')

#ifdef ECPP
   call phys_getopts(use_ECPP_out = use_ECPP)
   if (use_ECPP) then
      ! Add history coordinate variables for ECPP class dimensions
      call add_hist_coord('NCLASS_CL',    NCLASS_CL,    'NCLASS_CL')
      call add_hist_coord('ncls_ecpp_in', ncls_ecpp_in, 'ncls_ecpp_in')
      call add_hist_coord('NCLASS_PR',    NCLASS_PR,    'NCLASS_PR')
   end if
#endif

end subroutine crm_history_register
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
subroutine crm_history_init(species_class)
   use phys_control,    only: phys_getopts
   use physconst,       only: mwdry, cpair, spec_class_gas
   use ppgrid,          only: pcols, pver, pverp
   use constituents,    only: pcnst, cnst_name
   use cam_history,     only: addfld, add_default, horiz_only
   use crmdims,         only: crm_nx, crm_ny, crm_nz, crm_nx_rad, crm_ny_rad
#ifdef MODAL_AERO
   use cam_history,     only: fieldname_len
   use modal_aero_data, only: cnst_name_cw, lmassptr_amode, &
                              lmassptrcw_amode, nspec_amode, &
                              ntot_amode, numptr_amode, &
                              numptrcw_amode
#endif
   !----------------------------------------------------------------------------
   ! Interface variables

   ! species_class is defined as input so it needs to be outside of 
   ! MODAL_AERO condition for 1-moment micro to work
   integer, intent(in) :: species_class(:)

   !----------------------------------------------------------------------------
   ! local variables
   integer :: m
   logical :: use_ECPP
   character(len=16) :: MMF_microphysics_scheme

#ifdef MODAL_AERO
   integer :: l, lphase, lspec
   character(len=fieldname_len)   :: tmpname
   character(len=fieldname_len+3) :: fieldname
   character(128)                 :: long_name
   character(8)                   :: unit
#endif

   character(len=6), dimension(3) :: dims_crm_3D = (/'crm_nx','crm_ny','crm_nz'/)
   character(len=6), dimension(2) :: dims_crm_2D = (/'crm_nx','crm_ny'/)

   character(len=12), dimension(4) :: ecpp_coord_i = (/'ilev        ','NCLASS_CL   ','ncls_ecpp_in','NCLASS_PR   '/)
   character(len=12), dimension(4) :: ecpp_coord_m = (/'lev         ','NCLASS_CL   ','ncls_ecpp_in','NCLASS_PR   '/)
    
   !----------------------------------------------------------------------------

   call phys_getopts(use_ECPP_out = use_ECPP)
   call phys_getopts(MMF_microphysics_scheme_out = MMF_microphysics_scheme)

   call addfld('CRM_DEI  ', dims_crm_3D, 'A', 'micrometers', 'cloud scale Mitchell ice effective diameter')
   call addfld('CRM_REL  ', dims_crm_3D, 'A', 'micrometers', 'cloud scale droplet effective radius')
   call addfld('CRM_REI  ', dims_crm_3D, 'A', 'micrometers', 'cloud scale ice crystal effective radius')
 
   call addfld('CRM_FSNT',  dims_crm_2D, 'A',  'unitless', 'net TOA shortwave fluxes at CRM grids')
   call addfld('CRM_FSNTC', dims_crm_2D, 'A',  'unitless', 'net TOA clear-sky shortwave fluxes at CRM grids')
   call addfld('CRM_FSNS',  dims_crm_2D, 'A',  'unitless', 'net surface shortwave fluxes at CRM grids')
   call addfld('CRM_FSNSC', dims_crm_2D, 'A',  'unitless',  'net surface clear-sky shortwave fluxes at CRM grids')
   call addfld('CRM_FLNT',  dims_crm_2D, 'A',  'unitless', 'net TOA longwave fluxes at CRM grids')
   call addfld('CRM_FLNTC', dims_crm_2D, 'A',  'unitless', 'net TOA clear-sky longwave fluxes at CRM grids')
   call addfld('CRM_FLNS',  dims_crm_2D, 'A',  'unitless', 'net surface longwave fluxes at CRM grids')
   call addfld('CRM_FLNSC', dims_crm_2D, 'A',  'unitless', 'net surface clear-sky longwave fluxes at CRM grids')
 
   call addfld('CRM_AODVIS', dims_crm_2D,'A', 'unitless', 'Aerosol optical depth at 550nm in CRM grids', flag_xyfill=.true.)
   call addfld('CRM_AOD400', dims_crm_2D,'A', 'unitless', 'Aerosol optical depth at 400nm in CRM grids', flag_xyfill=.true.)
   call addfld('CRM_AOD700', dims_crm_2D,'A', 'unitless', 'Aerosol optical depth at 700nm in CRM grids', flag_xyfill=.true.)
   call addfld('CRM_AODVISZ',dims_crm_3D,'A', 'unitless', 'Aerosol optical depth at 500nm in CRM grids', flag_xyfill=.true.)
   call addfld('AOD400', horiz_only,'A', 'unitless', 'Aerosol optical depth at 400nm', flag_xyfill=.true.)
   call addfld('AOD700', horiz_only,'A', 'unitless', 'Aerosol optical depth at 700nm', flag_xyfill=.true.)
   
   call addfld('CRM_TK',  dims_crm_3D, 'A','m^2/s',   'Eddy viscosity from CRM')
   call addfld('CRM_TKH', dims_crm_3D, 'A','m^2/s',   'Eddy viscosity from CRM')

   if (MMF_microphysics_scheme .eq. 'm2005') then
      call addfld('SPNC    ',(/'lev'/), 'A', '/kg   ','Cloud water dropet number from CRM')
      call addfld('SPNI    ',(/'lev'/), 'A', '/kg   ','Cloud ice crystal number from CRM')
      call addfld('SPNS    ',(/'lev'/), 'A', '/kg   ','Snow particle number from CRM')
      call addfld('SPNG    ',(/'lev'/), 'A', '/kg   ','Graupel particle number from CRM')
      call addfld('SPNR    ',(/'lev'/), 'A', '/kg   ','Rain particle number from CRM')
  
      call addfld('CRM_FLIQ ',dims_crm_3D, 'A', '1      ','Frequency of Occurrence of Liquid' )
      call addfld('CRM_FICE ',dims_crm_3D, 'A', '1      ','Frequency of Occurrence of Ice' )
      call addfld('CRM_FRAIN',dims_crm_3D, 'A', '1      ','Frequency of Occurrence of Rain' )
      call addfld('CRM_FSNOW',dims_crm_3D, 'A', '1      ','Frequency of Occurrence of Snow' )
      call addfld('CRM_FGRAP',dims_crm_3D, 'A', '1      ','Frequency of Occurrence of Graupel' )
  
      call addfld('CRM_QS  ',dims_crm_3D, 'A', 'kg/kg   ','Snow mixing ratio from CRM' )
      call addfld('CRM_QG  ',dims_crm_3D, 'A', 'kg/kg   ','Graupel mixing ratio from CRM' )
      call addfld('CRM_QR  ',dims_crm_3D, 'A', 'kg/kg   ','Rain mixing ratio from CRM' )
  
      call addfld('CRM_NC  ',dims_crm_3D, 'A', '/kg     ','Cloud water dropet number from CRM' )
      call addfld('CRM_NI  ',dims_crm_3D, 'A', '/kg     ','Cloud ice crystal number from CRM' )
      call addfld('CRM_NS  ',dims_crm_3D, 'A', '/kg     ','Snow particle number from CRM' )
      call addfld('CRM_NG  ',dims_crm_3D, 'A', '/kg     ','Graupel particle number from CRM' )
      call addfld('CRM_NR  ',dims_crm_3D, 'A', '/kg     ','Rain particle number from CRM' )
  
      ! below is for *instantaneous* crm output
      call addfld('CRM_AUT  ',dims_crm_3D, 'A', '/s     ','Autoconversion cloud waterfrom CRM' )
      call addfld('CRM_ACC  ',dims_crm_3D, 'A', '/s     ','Accretion cloud water from CRM' )
      call addfld('CRM_EVPC ',dims_crm_3D, 'A', '/s     ','Evaporation cloud water from CRM' )
      call addfld('CRM_EVPR ',dims_crm_3D, 'A', '/s     ','Evaporation rain from CRM' )
      call addfld('CRM_MLT  ',dims_crm_3D, 'A', '/s     ','Melting ice snow graupel from CRM' )
      call addfld('CRM_SUB  ',dims_crm_3D, 'A', '/s     ','Sublimation ice snow graupel from CRM' )
      call addfld('CRM_DEP  ',dims_crm_3D, 'A', '/s     ','Deposition ice snow graupel from CRM' )
      call addfld('CRM_CON  ',dims_crm_3D, 'A', '/s     ','Condensation cloud water from CRM' )
  
      ! *gcm-grid and time-step-avg* process output
      call addfld('A_AUT  ',(/'lev'/), 'A', '/s   ','Avg autoconversion cloud water from CRM' )
      call addfld('A_ACC  ',(/'lev'/), 'A', '/s   ','Avg accretion cloud water from CRM' )
      call addfld('A_EVPC ',(/'lev'/), 'A', '/s   ','Avg evaporation cloud water from CRM' )
      call addfld('A_EVPR ',(/'lev'/), 'A', '/s   ','Avg evaporation rain from CRM' )
      call addfld('A_MLT  ',(/'lev'/), 'A', '/s   ','Avg melting ice snow graupel from CRM' )
      call addfld('A_SUB  ',(/'lev'/), 'A', '/s   ','Avg sublimation ice snow graupel from CRM' )
      call addfld('A_DEP  ',(/'lev'/), 'A', '/s   ','Avg deposition ice snow graupel from CRM' )
      call addfld('A_CON  ',(/'lev'/), 'A', '/s   ','Avg condensation cloud water from CRM' )
  
      call addfld('CRM_DES  ', dims_crm_3D, 'A', 'micrometers', 'cloud scale snow effective diameter')
      call addfld('CRM_MU   ', dims_crm_3D, 'A', 'micrometers', 'cloud scale droplet size distribution shape parameter for radiation')
      call addfld('CRM_LAMBDA',dims_crm_3D, 'A', 'micrometers', 'cloud scale slope of droplet distribution for radiation')
      call addfld('CRM_WVAR' , dims_crm_3D, 'A', 'm/s',         'vertical velocity variance from CRM')
   end if  ! MMF_microphysics_scheme .eq. 'm2005'

#ifdef ECPP
   if (use_ECPP) then
      ! Add ECPP output variables
      call addfld('ABND    ', ecpp_coord_i,'A','fraction','cloud fraction for each sub-sub class for full time period at layer boundary')
      call addfld('ABND_TF ', ecpp_coord_i,'A','fraction','cloud fraction for each sub-sub class for end-portion of time period at layer boundary')
      call addfld('MASFBND ', ecpp_coord_i,'A','kg/m2/s', 'sub-class vertical mass flux (kg/m2/s) at layer boundary')
      call addfld('ACEN    ', ecpp_coord_m,'A','fraction','cloud fraction for each sub-sub class for full time period at layer center')
      call addfld('ACEN_TF ', ecpp_coord_m,'A','fraction','cloud fraction for each sub-sub class for end-portion of time period at layer center')
      call addfld('RHCEN   ', ecpp_coord_m,'A','fraction','relative humidity for each sub-sub calss at layer center')
      call addfld('QCCEN   ', ecpp_coord_m,'A','kg/kg',   'cloud water for each sub-sub class at layer center')
      call addfld('QICEN   ', ecpp_coord_m,'A','kg/kg',   'cloud ice for each sub-sub class at layer center')
      call addfld('QSINK_AFCEN', ecpp_coord_m,'A','/s',   'cloud water loss rate from precip. using cloud water after precip. for each sub-sub class at layer center')
      call addfld('QSINK_BFCEN', ecpp_coord_m,'A','/s',   'cloud water loss rate from precip. using cloud water before precip. for each sub-sub class at layer center')
      call addfld('QSINK_AVGCEN',ecpp_coord_m,'A','/s',   'cloud water loss rate from precip. using avg cloud water and precip. rate for each sub-sub class at layer center')
      call addfld('PRAINCEN', ecpp_coord_m,'A','kg/kg/s', ' cloud water loss rate from precipitation (kg/kg/s) for each sub-sub class at layer center')
      call addfld('PRECRCEN', ecpp_coord_m,'A','kg/m2/s', 'liquid (rain) precipitation rate for each sub-sub class at layer center')
      call addfld('PRECSCEN', ecpp_coord_m,'A','kg/m2/s', 'solid (snow, graupel,...) precipitation rate for each sub-sub class at layer center')
      call addfld('WUPTHRES',     (/'ilev'/),'A','m/s',  'vertical velocity threshold for updraft')
      call addfld('WDNTHRES',     (/'ilev'/),'A','m/s',  'vertical velocity threshold for dndraft')
      call addfld('WWQUI_CEN',    (/ 'lev'/),'A','m2/s2','vertical velocity variance in the quiescent class, layer center')
      call addfld('WWQUI_CLD_CEN',(/ 'lev'/),'A','m2/s2','vertical velocity variance in the cloudy quiescent class, layer center')
      call addfld('WWQUI_BND',    (/'ilev'/),'A','m2/s2','vertical velocity variance in the quiescent class, layer boundary')
      call addfld('WWQUI_CLD_BND',(/'ilev'/),'A','m2/s2','vertical velocity variance in the cloudy quiescent class, layer boundary')
   endif
#endif /* ECPP */

   call addfld('MU_CRM   ',(/'lev'/),'A','Pa/s','mass flux up from CRM')
   call addfld('MD_CRM   ',(/'lev'/),'A','Pa/s','mass flux down from CRM')
   call addfld('DU_CRM   ',(/'lev'/),'A','/s',  'detrainment from updraft from CRM')
   call addfld('EU_CRM   ',(/'lev'/),'A','/s',  'entraiment rate from updraft')
   call addfld('ED_CRM   ',(/'lev'/),'A','/s',  'entraiment rate from downdraft')
 
   do m = 1, pcnst 
     if(cnst_name(m) == 'DMS') then 
        call addfld('DMSCONV',   (/ 'lev' /), 'A', 'kg/kg/s',  'DMS tendency from ZM convection')
     end if
     if(cnst_name(m) == 'SO2') then 
        call addfld('SO2CONV',   (/ 'lev' /), 'A', 'kg/kg/s',  'SO2 tendency from ZM convection')
      end if
   end do
 
   call addfld('SPQRL    ', (/'lev'/), 'A','K/s', 'long-wave heating rate')
   call addfld('SPQRS    ', (/'lev'/), 'A','K/s', 'short-wave heating rate')
   call addfld('LENGC    ', (/'ilev'/),'A','m  ', 'Mixing length scale for the calcuation of vertical difusivity')
 
   call addfld('SPKVH     ',(/'ilev'/),'A','m2/s    ','Vertical diffusivity used in dropmixnuc in the MMF call')
   call addfld('SPWTKE   ', (/'lev'/), 'A','m/s',     'Standard deviation of updraft velocity')
   call addfld('SPLCLOUD  ',(/'lev'/), 'A','        ','Liquid cloud fraction')
#ifdef MAML
   call addfld('CRM_SHF ',dims_crm_2D, 'I', 'W/m2    ', 'CRM Sfc sensible heat flux')
   call addfld('CRM_LHF ',dims_crm_2D, 'I', 'W/m2    ', 'CRM Sfc latent heat flux'  )
   call addfld('CRM_SNOW',dims_crm_2D, 'I', 'm/s     ', 'CRM Snow Rate'             )
   call addfld('CRM_PCP ',dims_crm_2D, 'I', 'm/s     ', 'CRM Precipitation Rate'    )
#endif

   call addfld('SPNDROPMIX',(/ 'lev' /),'A','#/kg/s  ','Droplet number mixing')
   call addfld('SPNDROPSRC',(/ 'lev' /),'A','#/kg/s  ','Droplet number source')
   call addfld('SPNDROPCOL',horiz_only,'A', '#/m2    ','Column droplet number')

   call addfld('PRES    ',(/'lev'/),  'A','Pa      ','Pressure' )
   call addfld('DPRES   ',(/'lev'/),  'A','Pa      ','Pressure thickness of layer' )
   call addfld('SPDT    ',(/'lev'/),  'A','K/s     ','T tendency due to CRM' )
   call addfld('SPDQ    ',(/'lev'/),  'A','kg/kg/s ','Q tendency due to CRM' )
   call addfld('SPDQC   ',(/'lev'/),  'A','kg/kg/s ','QC tendency due to CRM' )
   call addfld('SPDQI   ',(/'lev'/),  'A','kg/kg/s ','QI tendency due to CRM' )
   call addfld('SPMC    ',(/'lev'/),  'A','kg/m2/s ','Total mass flux from CRM' )
   call addfld('SPMCUP  ',(/'lev'/),  'A','kg/m2/s ','Updraft mass flux from CRM' )
   call addfld('SPMCDN  ',(/'lev'/),  'A','kg/m2/s ','Downdraft mass flux from CRM' )
   call addfld('SPMCUUP ',(/'lev'/),  'A','kg/m2/s ','Unsaturated updraft mass flux from CRM' )
   call addfld('SPMCUDN ',(/'lev'/),  'A','kg/m2/s ','Unsaturated downdraft mass flux from CRM' )
   call addfld('SPQC    ',(/'lev'/),  'A','kg/kg   ','Cloud water from CRM' )
   call addfld('SPQI    ',(/'lev'/),  'A','kg/kg   ','Cloud ice from CRM' )
   call addfld('SPQS    ',(/'lev'/),  'A','kg/kg   ','Snow from CRM' )
   call addfld('SPQG    ',(/'lev'/),  'A','kg/kg   ','Graupel from CRM' )
   call addfld('SPQR    ',(/'lev'/),  'A','kg/kg   ','Rain from CRM' )
   call addfld('SPQTFLX ',(/'lev'/),  'A','kg/m2/s ','Nonprecip. water flux from CRM' )
   call addfld('SPUFLX  ',(/'lev'/),  'A','m2/s2   ','x-momentum flux from CRM' )
   call addfld('SPVFLX  ',(/'lev'/),  'A','m2/s2   ','y-momentum flux from CRM' )
   call addfld('SPQTFLXS',(/'lev'/),  'A','kg/m2/s ','SGS Nonprecip. water flux from CRM' )
   call addfld('SPTKE   ',(/'lev'/),  'A','kg/m/s2 ','Total TKE in CRM' )
   call addfld('SPTKES  ',(/'lev'/),  'A','kg/m/s2 ','SGS TKE in CRM' )
   call addfld('SPTK    ',(/'lev'/),  'A','m2/s    ','SGS TK in CRM' )
   call addfld('SPQPFLX ',(/'lev'/),  'A','kg/m2/s ','Precip. water flux from CRM' )
   call addfld('SPPFLX  ',(/'lev'/),  'A','m/s     ','Precipitation flux from CRM' )
   call addfld('SPQTLS  ',(/'lev'/),  'A','kg/kg/s ','L.S. Vapor Tendency from CRM' )
   call addfld('SPQTTR  ',(/'lev'/),  'A','kg/kg/s ','Nonprec. water transport from CRM' )
   call addfld('SPQPTR  ',(/'lev'/),  'A','kg/kg/s ','Prec. water transport from CRM' )
   call addfld('SPQPEVP ',(/'lev'/),  'A','kg/kg/s ','Prec. water evaporation from CRM' )
   call addfld('SPQPFALL',(/'lev'/),  'A','kg/kg/s ','Prec. water fall-out from CRM' )
   call addfld('SPQPSRC ',(/'lev'/),  'A','kg/kg/s ','Prec. water source from CRM' )
   call addfld('SPTLS   ',(/'lev'/),  'A','kg/kg/s ','L.S. LIWSE Tendency from CRM' )
   call addfld('TIMINGF ', horiz_only,'A','        ','CRM CPU usage efficiency: 1 - ideal' )
   call addfld('CLOUDTOP',(/'lev'/),  'A','        ','Cloud Top PDF' )
#if defined(MMF_MOMENTUM_FEEDBACK) || defined(MMF_ESMT)
   call addfld('UCONVMOM',(/'lev'/), 'A', 'm/s2 ','U tendency due to CRM' )
   call addfld('VCONVMOM',(/'lev'/), 'A', 'm/s2 ','V tendency due to CRM' )
#endif
#if defined(MMF_ESMT)
   call addfld('U_ESMT',(/'lev'/), 'A', 'm/s2 ','U tendency due to CRM (ESMT)' )
   call addfld('V_ESMT',(/'lev'/), 'A', 'm/s2 ','V tendency due to CRM (ESMT)' )
#endif

   call addfld('CRM_U   ',dims_crm_3D, 'I', 'm/s     ', 'CRM x-wind' )
   call addfld('CRM_V   ',dims_crm_3D, 'I', 'm/s     ', 'CRM y-wind' )
   call addfld('CRM_W   ',dims_crm_3D, 'I', 'm/s     ', 'CRM z-wind' )
   call addfld('CRM_T   ',dims_crm_3D, 'I', 'K       ', 'CRM Temperature' )
   call addfld('CRM_QV  ',dims_crm_3D, 'I', 'kg/kg   ', 'CRM Water Vapor' )
   call addfld('CRM_QC  ',dims_crm_3D, 'I', 'kg/kg   ', 'CRM Cloud Water' )
   call addfld('CRM_QI  ',dims_crm_3D, 'I', 'kg/kg   ', 'CRM Cloud Ice' )
   call addfld('CRM_QPC ',dims_crm_3D, 'I', 'kg/kg   ', 'CRM Precipitating Water' )
   call addfld('CRM_QPI ',dims_crm_3D, 'I', 'kg/kg   ', 'CRM Precipitating Ice' )
   call addfld('CRM_PREC',dims_crm_2D, 'I', 'm/s     ', 'CRM Precipitation Rate' )

   call addfld('SPTVFLUX ',(/ 'lev' /), 'A', 'W/m2  ','Buoyancy Flux from CRM' )
   call addfld('SPBUOY   ',(/ 'lev' /), 'A', 'W/m3  ','Buoyancy Term from CRM' )
   call addfld('SPBUOYSD ',(/ 'lev' /), 'A', 'W/m3  ','Std Dev of Buoyancy Term from CRM' )
   call addfld('SPMSEF   ',(/ 'lev' /), 'A', 'W/m2  ','Moist Static Energy Flux from CRM' )
   call addfld('SPQVFLUX ',(/ 'lev' /), 'A', 'W/m2  ','Water Wapor Flux from CRM' )

#ifdef MODAL_AERO
   ! add dropmixnuc tendencies for all modal aerosol species
   do m = 1, ntot_amode
      do lphase = 1, 2
         do lspec = 0, nspec_amode(m)+1   ! loop over number + chem constituents + water
            unit = 'kg/m2/s'
            if (lspec == 0) then   ! number
               unit = '#/m2/s'
               if (lphase == 1) then
                  l = numptr_amode(m)
               else
                  l = numptrcw_amode(m)
               endif
            else if (lspec <= nspec_amode(m)) then   ! non-water mass
               if (lphase == 1) then
                  l = lmassptr_amode(lspec,m)
               else
                  l = lmassptrcw_amode(lspec,m)
               endif
            else   ! water mass
               cycle
            end if
            if (lphase == 1) then
               tmpname = cnst_name(l)
            else
               tmpname = cnst_name_cw(l)
            end if

            fieldname = trim(tmpname) // '_mixnuc1sp'
            long_name = trim(tmpname) // ' dropmixnuc mixnuc column tendency in the mmf one '
            call addfld( fieldname,  horiz_only, 'A', unit, long_name)
            call add_default( fieldname, 1, ' ' )

         end do ! lspec
      end do ! lphase
   end do ! m  

   do m = 1, pcnst
      if(species_class(m).eq.spec_class_gas) then
         fieldname = trim(cnst_name(m)) // '_mixnuc1sp'
         long_name = trim(cnst_name(m)) // ' dropmixnuc mixnuc column tendency in the mmf one '
         call addfld( fieldname,  horiz_only, 'A', unit, long_name)
         call add_default( fieldname, 1, ' ' )
      end if
   end do
#endif


   call add_default ('SPDT    ', 1, ' ')
   call add_default ('SPDQ    ', 1, ' ')
   call add_default ('SPDQC   ', 1, ' ')
   call add_default ('SPDQI   ', 1, ' ')
   call add_default ('SPMC    ', 1, ' ')
   call add_default ('SPMCUP  ', 1, ' ')
   call add_default ('SPMCDN  ', 1, ' ')
   call add_default ('SPMCUUP ', 1, ' ')
   call add_default ('SPMCUDN ', 1, ' ')
   call add_default ('SPQC    ', 1, ' ')
   call add_default ('SPQI    ', 1, ' ')
   call add_default ('SPQS    ', 1, ' ')
   call add_default ('SPQG    ', 1, ' ')
   call add_default ('SPQR    ', 1, ' ')
   call add_default ('SPQTFLX ', 1, ' ')
   call add_default ('SPQTFLXS', 1, ' ')
   call add_default ('SPTKE   ', 1, ' ')
   call add_default ('SPTKES  ', 1, ' ')
   call add_default ('SPTK    ', 1, ' ')
   call add_default ('SPQPFLX ', 1, ' ')
   call add_default ('SPPFLX  ', 1, ' ')
   call add_default ('SPQTLS  ', 1, ' ')
   call add_default ('SPQTTR  ', 1, ' ')
   call add_default ('SPQPTR  ', 1, ' ')
   call add_default ('SPQPEVP ', 1, ' ')
   call add_default ('SPQPFALL', 1, ' ')
   call add_default ('SPQPSRC ', 1, ' ')
   call add_default ('SPTLS   ', 1, ' ')
   call add_default ('CLOUDTOP', 1, ' ')
   call add_default ('TIMINGF ', 1, ' ')

   call add_default ('AOD400',  1, ' ')
   call add_default ('AOD700',  1, ' ')

   call add_default ('SPTVFLUX  ', 1, ' ')
   call add_default ('SPBUOY    ', 1, ' ')
   call add_default ('SPBUOYSD  ', 1, ' ')
   call add_default ('SPMSEF    ', 1, ' ')
   call add_default ('SPQVFLUX  ', 1, ' ')

   call add_default ('SPKVH     ', 1, ' ')
   call add_default ('SPWTKE    ', 1, ' ')
   call add_default ('SPLCLOUD  ', 1, ' ')
   call add_default ('SPNDROPSRC', 1, ' ')
   call add_default ('SPNDROPMIX', 1, ' ')
   call add_default ('SPNDROPCOL', 1, ' ')

   if (MMF_microphysics_scheme .eq. 'm2005') then
      call add_default ('SPNC    ', 1, ' ')
      call add_default ('SPNI    ', 1, ' ')
      call add_default ('SPNS    ', 1, ' ')
      call add_default ('SPNG    ', 1, ' ')
      call add_default ('SPNR    ', 1, ' ')
   end if

#if defined(MMF_MOMENTUM_FEEDBACK) || defined(MMF_ESMT)
   call add_default ('UCONVMOM', 1, ' ')
   call add_default ('VCONVMOM', 1, ' ')
#endif

#if defined(MMF_ESMT)
   call add_default ('U_ESMT', 1, ' ')
   call add_default ('V_ESMT', 1, ' ')
#endif

end subroutine crm_history_init
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
subroutine crm_history_out(state, ptend, crm_state, crm_rad, crm_output, crm_ecpp_output, qrs, qrl)
   use physics_types,          only: physics_state, physics_tend, physics_ptend
   use phys_control,           only: phys_getopts
   use crm_state_module,       only: crm_state_type
   use crm_rad_module,         only: crm_rad_type
   use crm_output_module,      only: crm_output_type
   use crm_ecpp_output_module, only: crm_ecpp_output_type
   use constituents,           only: cnst_get_ind
   use ppgrid,                 only: pcols, pver, pverp
   use physconst,              only: cpair
   use cam_history,            only: outfld
   
   !----------------------------------------------------------------------------
   ! interface variables
   type(physics_state),              intent(in) :: state             ! Global model state 
   type(physics_ptend),              intent(in) :: ptend             ! CRM output tendencies
   type(crm_state_type),             intent(in) :: crm_state         ! CRM state 
   type(crm_rad_type),        target,intent(in) :: crm_rad           ! CRM rad variables
   type(crm_output_type),     target,intent(in) :: crm_output        ! CRM output
   type(crm_ecpp_output_type),       intent(in) :: crm_ecpp_output   ! CRM output for ECPP

   real(r8), dimension(:,:), intent(in) :: qrs        ! shortwave radiative heating rate
   real(r8), dimension(:,:), intent(in) :: qrl        ! longwave radiative heating rate

   !----------------------------------------------------------------------------
   ! local variables
   real(r8) :: cwp      (pcols,pver)   ! in-cloud cloud (total) water path (kg/m2)
   real(r8) :: gwp      (pcols,pver)   ! grid-box cloud (total) water path (kg/m2)
   real(r8) :: cicewp   (pcols,pver)   ! in-cloud cloud ice water path (kg/m2)
   real(r8) :: cliqwp   (pcols,pver)   ! in-cloud cloud liquid water path (kg/m2)
   real(r8) :: tgicewp  (pcols)        ! Vertically integrated ice water path (kg/m2
   real(r8) :: tgliqwp  (pcols)        ! Vertically integrated liquid water path (kg/m2)
   real(r8) :: tgwp     (pcols)        ! Vertically integrated (total) cloud water path  (kg/m2)
   real(r8) :: SPDT_out (pcols,pver)   ! CRM heating tendency
   integer :: lchnk                    ! chunk identifier
   integer :: ncol                     ! number of atmospheric columns
   integer :: ixcldliq, ixcldice       ! liquid and ice constituent indices
   integer :: i, k                     ! loop iterators
   logical :: use_ECPP
   character(len=16) :: MMF_microphysics_scheme

   !----------------------------------------------------------------------------

   call cnst_get_ind('CLDLIQ', ixcldliq)
   call cnst_get_ind('CLDICE', ixcldice)

   call phys_getopts(use_ECPP_out = use_ECPP)
   call phys_getopts(MMF_microphysics_scheme_out = MMF_microphysics_scheme)

   lchnk = state%lchnk
   ncol  = state%ncol

   ! Subtract radiative heating for SPDT output
   SPDT_out(:ncol,:pver) = ( ptend%s(:ncol,:pver) - qrs(:ncol,:pver) - qrl(:ncol,:pver) )/cpair

   !----------------------------------------------------------------------------
   ! CRM domain average Tendencies
   call outfld('SPDT    ',SPDT_out,              pcols, lchnk ) 
   call outfld('SPDQ    ',ptend%q(1,1,1),        pcols, lchnk )
   call outfld('SPDQC   ',ptend%q(1,1,ixcldliq), pcols, lchnk )
   call outfld('SPDQI   ',ptend%q(1,1,ixcldice), pcols, lchnk )

   ! CRM radiative heating rate
   ! NOTE: We output the heating rate here because this is the heating that
   ! is applied to the CRM at this GCM timestep, but note that this heating
   ! rate is offset in time from when the radiative heating rate is
   ! calculated in the radiative transfer interface, which comes AFTER this
   ! routine call (because we want it to use the clouds simulated from this
   ! CRM call). Thus, comparing this heating rate with CRM_QRS + CRM_QRL
   ! output in radiation_tend will show a time lag.
   call outfld('CRM_QRAD',crm_rad%qrad, pcols, lchnk )

   ! GCM level rad heating tendencies
   call outfld('SPQRL   ',qrl/cpair, pcols, lchnk )
   call outfld('SPQRS   ',qrs/cpair, pcols, lchnk )

   ! Why do we output this here?
   call outfld('PRES    ',state%pmid, pcols, lchnk )
   call outfld('DPRES   ',state%pdel, pcols, lchnk )

   ! CRM state variables on CRM grid
   call outfld('CRM_U   ',crm_state%u_wind,      pcols, lchnk )
   call outfld('CRM_V   ',crm_state%v_wind,      pcols, lchnk )
   call outfld('CRM_W   ',crm_state%w_wind,      pcols, lchnk )
   call outfld('CRM_T   ',crm_state%temperature, pcols, lchnk )

   if (MMF_microphysics_scheme .eq. 'sam1mom') then
      call outfld('CRM_QV  ',(crm_state%qt-crm_output%qcl-crm_output%qci), pcols, lchnk )
   else if (MMF_microphysics_scheme .eq. 'm2005') then 
      call outfld('CRM_QV  ',crm_state%qt-crm_output%qcl, pcols, lchnk )
   endif

   ! CRM condensate and precipitation on CRM grid
   call outfld('CRM_QC  ',crm_output%qcl,         pcols, lchnk )
   call outfld('CRM_QI  ',crm_output%qci,         pcols, lchnk )
   call outfld('CRM_QPC ',crm_output%qpl,         pcols, lchnk )
   call outfld('CRM_QPI ',crm_output%qpi,         pcols, lchnk )
   call outfld('CRM_PREC',crm_output%prec_crm,    pcols, lchnk )
   call outfld('CRM_TK ', crm_output%tk(:,:,:,:), pcols, lchnk )  
   call outfld('CRM_TKH', crm_output%tkh(:,:,:,:),pcols, lchnk ) 

   ! CRM domain average condensate and precipitation
   call outfld('SPQC    ',crm_output%qc_mean, pcols ,lchnk )
   call outfld('SPQI    ',crm_output%qi_mean, pcols ,lchnk )
   call outfld('SPQS    ',crm_output%qs_mean, pcols ,lchnk )
   call outfld('SPQG    ',crm_output%qg_mean, pcols ,lchnk )
   call outfld('SPQR    ',crm_output%qr_mean, pcols ,lchnk )

   ! CRM domain average fluxes
   call outfld('SPQTFLX ',crm_output%flux_qt,        pcols, lchnk )
   call outfld('SPUFLX  ',crm_output%flux_u,         pcols, lchnk )
   call outfld('SPVFLX  ',crm_output%flux_v,         pcols, lchnk )
   call outfld('SPTKE   ',crm_output%tkez,           pcols, lchnk )
   call outfld('SPTKES  ',crm_output%tkesgsz,        pcols, lchnk )
   call outfld('SPTK    ',crm_output%tkz,            pcols, lchnk )
   call outfld('SPQTFLXS',crm_output%fluxsgs_qt,     pcols, lchnk )
   call outfld('SPQPFLX ',crm_output%flux_qp,        pcols, lchnk )
   call outfld('SPPFLX  ',crm_output%precflux,       pcols, lchnk )
   call outfld('SPQTLS  ',crm_output%qt_ls,          pcols, lchnk )
   call outfld('SPQTTR  ',crm_output%qt_trans,       pcols, lchnk )
   call outfld('SPQPTR  ',crm_output%qp_trans,       pcols, lchnk )
   call outfld('SPQPEVP ',crm_output%qp_evp,         pcols, lchnk )
   call outfld('SPQPFALL',crm_output%qp_fall,        pcols, lchnk )
   call outfld('SPQPSRC ',crm_output%qp_src,         pcols, lchnk )
   call outfld('SPTLS   ',crm_output%t_ls,           pcols, lchnk )

   ! NOTE: these should overwrite cloud outputs from non-MMF routines
   call outfld('CLOUD   ',crm_output%cld,    pcols, lchnk )
   call outfld('CLDTOT  ',crm_output%cltot,  pcols, lchnk )
   call outfld('CLDHGH  ',crm_output%clhgh,  pcols, lchnk )
   call outfld('CLDMED  ',crm_output%clmed,  pcols, lchnk )
   call outfld('CLDLOW  ',crm_output%cllow,  pcols, lchnk )
   call outfld('CLOUDTOP',crm_output%cldtop, pcols, lchnk )

   call outfld('TIMINGF ',crm_output%timing_factor  ,pcols,lchnk)

   ! CRM mass flux
   call outfld('SPMC    ', crm_output%mctot,  pcols, lchnk )
   call outfld('SPMCUP  ', crm_output%mcup,   pcols, lchnk )
   call outfld('SPMCDN  ', crm_output%mcdn,   pcols, lchnk )
   call outfld('SPMCUUP ', crm_output%mcuup,  pcols, lchnk )
   call outfld('SPMCUDN ', crm_output%mcudn,  pcols, lchnk )
   call outfld('MU_CRM  ', crm_output%mu_crm, pcols, lchnk )
   call outfld('MD_CRM  ', crm_output%md_crm, pcols, lchnk )
   call outfld('EU_CRM  ', crm_output%eu_crm, pcols, lchnk )
   call outfld('DU_CRM  ', crm_output%du_crm, pcols, lchnk )
   call outfld('ED_CRM  ', crm_output%ed_crm, pcols, lchnk )

#ifdef MAML
   ! Lower boundary fluxes
   call outfld('CRM_SHF ', cam_in%shf,         pcols, lchnk )
   call outfld('CRM_LHF ', cam_in%lhf,         pcols, lchnk )
   call outfld('CRM_SNOW', crm_output%crm_pcp, pcols, lchnk )
   call outfld('CRM_PCP',  crm_output%crm_snw, pcols, lchnk )
#endif

#ifdef m2005
   if (MMF_microphysics_scheme .eq. 'm2005') then
      ! index is defined in MICRO_M2005/microphysics.F90
      ! Be cautious to use them here. They are defined in crm codes, and these codes are called only 
      ! after the subroutine of crm is called. So they can only be used after the 'crm' subroutine. 
      ! incl, inci, ... can not be used here, for they are defined before we call them???
      call outfld('CRM_NC ',crm_state%nc(:,:,:,:), pcols, lchnk )
      call outfld('CRM_NI ',crm_state%ni(:,:,:,:), pcols, lchnk )
      call outfld('CRM_NR ',crm_state%nr(:,:,:,:), pcols, lchnk )
      call outfld('CRM_NS ',crm_state%ns(:,:,:,:), pcols, lchnk )
      call outfld('CRM_NG ',crm_state%ng(:,:,:,:), pcols, lchnk )
      call outfld('CRM_QR ',crm_state%qr(:,:,:,:), pcols, lchnk )
      call outfld('CRM_QS ',crm_state%qs(:,:,:,:), pcols, lchnk )
      call outfld('CRM_QG ',crm_state%qg(:,:,:,:), pcols, lchnk )
      
      call outfld('CRM_WVAR',crm_output%wvar, pcols, lchnk)

      call outfld('CRM_AUT', crm_output%aut, pcols, lchnk)
      call outfld('CRM_ACC', crm_output%acc, pcols, lchnk)
      call outfld('CRM_MLT', crm_output%mlt, pcols, lchnk)
      call outfld('CRM_SUB', crm_output%sub, pcols, lchnk)
      call outfld('CRM_DEP', crm_output%dep, pcols, lchnk)
      call outfld('CRM_CON', crm_output%con, pcols, lchnk)
      call outfld('CRM_EVPC',crm_output%evpc,pcols, lchnk)
      call outfld('CRM_EVPR',crm_output%evpr,pcols, lchnk)
      
      call outfld('A_AUT', crm_output%aut_a, pcols, lchnk)
      call outfld('A_ACC', crm_output%acc_a, pcols, lchnk)
      call outfld('A_MLT', crm_output%mlt_a, pcols, lchnk)
      call outfld('A_SUB', crm_output%sub_a, pcols, lchnk)
      call outfld('A_DEP', crm_output%dep_a, pcols, lchnk)
      call outfld('A_CON', crm_output%con_a, pcols, lchnk)
      call outfld('A_EVPC',crm_output%evpc_a,pcols, lchnk)
      call outfld('A_EVPR',crm_output%evpr_a,pcols, lchnk)

      call outfld('SPNC    ',crm_output%nc_mean, pcols, lchnk )
      call outfld('SPNI    ',crm_output%ni_mean, pcols, lchnk )
      call outfld('SPNS    ',crm_output%ns_mean, pcols, lchnk )
      call outfld('SPNG    ',crm_output%ng_mean, pcols, lchnk )
      call outfld('SPNR    ',crm_output%nr_mean, pcols, lchnk )
   endif ! m2005
#endif /* m2005 */

   !----------------------------------------------------------------------------
   ! Compute liquid water paths (for diagnostics only)
   !----------------------------------------------------------------------------
   tgicewp(:ncol) = 0.
   tgliqwp(:ncol) = 0.
   do k = 1,pver
      do i = 1,ncol
         cicewp(i,k) = crm_output%gicewp(i,k) * 1.0e-3 / max(0.01_r8,crm_output%cld(i,k)) ! In-cloud ice water path.  g/m2 --> kg/m2
         cliqwp(i,k) = crm_output%gliqwp(i,k) * 1.0e-3 / max(0.01_r8,crm_output%cld(i,k)) ! In-cloud liquid water path. g/m2 --> kg/m2
         tgicewp(i)  = tgicewp(i) + crm_output%gicewp(i,k) *1.0e-3 ! grid cell mean ice water path.  g/m2 --> kg/m2
         tgliqwp(i)  = tgliqwp(i) + crm_output%gliqwp(i,k) *1.0e-3 ! grid cell mean ice water path.  g/m2 --> kg/m2
      end do
   end do
   tgwp(:ncol) = tgicewp(:ncol) + tgliqwp(:ncol)
   gwp(:ncol,:pver) = crm_output%gicewp(:ncol,:pver) + crm_output%gliqwp(:ncol,:pver)
   cwp(:ncol,:pver) = cicewp(:ncol,:pver) + cliqwp(:ncol,:pver)

   call outfld('GCLDLWP' ,gwp,     pcols, lchnk)
   call outfld('TGCLDCWP',tgwp,    pcols, lchnk)
   call outfld('TGCLDLWP',tgliqwp, pcols, lchnk)
   call outfld('TGCLDIWP',tgicewp, pcols, lchnk)
   call outfld('ICLDTWP' ,cwp,     pcols, lchnk)
   call outfld('ICLDIWP' ,cicewp,  pcols, lchnk)

   !----------------------------------------------------------------------------
   !----------------------------------------------------------------------------
#ifdef ECPP
   if (use_ECPP) then
      call outfld('ACEN    ',      crm_ecpp_output%acen,             pcols, lchnk )
      call outfld('ABND    ',      crm_ecpp_output%abnd,             pcols, lchnk )
      call outfld('ACEN_TF ',      crm_ecpp_output%acen_tf,          pcols, lchnk )
      call outfld('ABND_TF ',      crm_ecpp_output%abnd_tf,          pcols, lchnk )
      call outfld('MASFBND ',      crm_ecpp_output%massflxbnd,       pcols, lchnk )
      call outfld('RHCEN   ',      crm_ecpp_output%rhcen,            pcols, lchnk )
      call outfld('QCCEN   ',      crm_ecpp_output%qcloudcen,        pcols, lchnk )
      call outfld('QICEN   ',      crm_ecpp_output%qicecen,          pcols, lchnk )
      call outfld('QSINK_AFCEN',   crm_ecpp_output%qlsink_afcen,     pcols, lchnk )
      call outfld('PRECRCEN',      crm_ecpp_output%precrcen,         pcols, lchnk )
      call outfld('PRECSCEN',      crm_ecpp_output%precsolidcen,     pcols, lchnk )
      call outfld('WUPTHRES',      crm_ecpp_output%wupthresh_bnd,    pcols, lchnk )
      call outfld('WDNTHRES',      crm_ecpp_output%wdownthresh_bnd,  pcols, lchnk )
      call outfld('WWQUI_CEN',     crm_ecpp_output%wwqui_cen,        pcols, lchnk )
      call outfld('WWQUI_CLD_CEN', crm_ecpp_output%wwqui_cloudy_cen, pcols, lchnk )
      call outfld('WWQUI_BND',     crm_ecpp_output%wwqui_cen,        pcols, lchnk )
      call outfld('WWQUI_CLD_BND', crm_ecpp_output%wwqui_cloudy_cen, pcols, lchnk )
      call outfld('QSINK_BFCEN',   crm_ecpp_output%qlsink_bfcen,     pcols, lchnk )
      call outfld('QSINK_AVGCEN',  crm_ecpp_output%qlsink_avgcen,    pcols, lchnk )
      call outfld('PRAINCEN',      crm_ecpp_output%praincen,         pcols, lchnk )
   end if ! use_ECPP
#endif /* ECPP */
   !----------------------------------------------------------------------------
   ! CRM momentum tendencies
   !----------------------------------------------------------------------------

#if defined( MMF_ESMT )
   call outfld('U_ESMT',crm_output%u_tend_esmt, pcols, lchnk )
   call outfld('V_ESMT',crm_output%v_tend_esmt, pcols, lchnk )
#endif /* MMF_ESMT */

#if defined(MMF_MOMENTUM_FEEDBACK)
   ! Make sure these tendencies are set (and rotated) in crm_physics_tend()
   call outfld('UCONVMOM',ptend%u, pcols, lchnk )
   call outfld('VCONVMOM',ptend%v, pcols, lchnk )
#endif /* MMF_MOMENTUM_FEEDBACK */

end subroutine crm_history_out
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
end module crm_history