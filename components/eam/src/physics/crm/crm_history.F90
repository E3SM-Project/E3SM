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

   logical :: use_ECPP                 ! explicit cloud parameterized pollutants

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
   use phys_control,     only: phys_getopts
   use physconst,        only: mwdry, cpair, spec_class_gas
   use ppgrid,           only: pcols, pver, pverp
   use constituents,     only: pcnst, cnst_name
   use cam_history,      only: addfld, add_default, horiz_only
   use crmdims,          only: crm_nx, crm_ny, crm_nz, crm_nx_rad, crm_ny_rad
   use cloud_cover_diags,only: cloud_cover_diags_init
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
   logical :: use_ECPP                 ! explicit cloud parameterized pollutants
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
   ! Register available MMF output variables
   call phys_getopts(use_ECPP_out = use_ECPP)
   call phys_getopts(MMF_microphysics_scheme_out = MMF_microphysics_scheme)

   !----------------------------------------------------------------------------
   ! Instantaneous CRM grid variables
   call addfld('CRM_U    ',dims_crm_3D, 'I','m/s',     'CRM x-wind' )
   call addfld('CRM_V    ',dims_crm_3D, 'I','m/s',     'CRM y-wind' )
   call addfld('CRM_W    ',dims_crm_3D, 'I','m/s',     'CRM z-wind' )
   call addfld('CRM_T    ',dims_crm_3D, 'I','K',       'CRM Temperature' )
   call addfld('CRM_QV   ',dims_crm_3D, 'I','kg/kg',   'CRM Water Vapor' )
   call addfld('CRM_QC   ',dims_crm_3D, 'I','kg/kg',   'CRM Cloud Water' )
   call addfld('CRM_QI   ',dims_crm_3D, 'I','kg/kg',   'CRM Cloud Ice' )
   call addfld('CRM_QPC  ',dims_crm_3D, 'I','kg/kg',   'CRM Precipitating Water' )
   call addfld('CRM_QPI  ',dims_crm_3D, 'I','kg/kg',   'CRM Precipitating Ice' )
   call addfld('CRM_PREC ',dims_crm_2D, 'I','m/s',     'CRM Precipitation Rate' )
   call addfld('CRM_DEI  ',dims_crm_3D, 'A','m^-6',    'cloud scale Mitchell ice effective diameter')
   call addfld('CRM_REL  ',dims_crm_3D, 'A','m^-6',    'cloud scale droplet effective radius')
   call addfld('CRM_REI  ',dims_crm_3D, 'A','m^-6',    'cloud scale ice crystal effective radius') 
   call addfld('CRM_FSNT ',dims_crm_2D, 'A','unitless','net TOA shortwave fluxes at CRM grids')
   call addfld('CRM_FSNTC',dims_crm_2D, 'A','unitless','net TOA clear-sky shortwave fluxes at CRM grids')
   call addfld('CRM_FSNS ',dims_crm_2D, 'A','unitless','net surface shortwave fluxes at CRM grids')
   call addfld('CRM_FSNSC',dims_crm_2D, 'A','unitless','net surface clear-sky shortwave fluxes at CRM grids')
   call addfld('CRM_FLNT ',dims_crm_2D, 'A','unitless','net TOA longwave fluxes at CRM grids')
   call addfld('CRM_FLNTC',dims_crm_2D, 'A','unitless','net TOA clear-sky longwave fluxes at CRM grids')
   call addfld('CRM_FLNS ',dims_crm_2D, 'A','unitless','net surface longwave fluxes at CRM grids')
   call addfld('CRM_FLNSC',dims_crm_2D, 'A','unitless','net surface clear-sky longwave fluxes at CRM grids')
   call addfld('CRM_TK   ',dims_crm_3D, 'A','m^2/s',   'Eddy viscosity from CRM')
   call addfld('CRM_TKH  ',dims_crm_3D, 'A','m^2/s',   'Eddy conductivity/thermal diffusion from CRM')

   !----------------------------------------------------------------------------
   ! Aerosol optical depth
   call addfld('CRM_AODVIS', dims_crm_2D,'A', 'unitless', 'Aerosol optical depth at 550nm in CRM grids', flag_xyfill=.true.)
   call addfld('CRM_AOD400', dims_crm_2D,'A', 'unitless', 'Aerosol optical depth at 400nm in CRM grids', flag_xyfill=.true.)
   call addfld('CRM_AOD700', dims_crm_2D,'A', 'unitless', 'Aerosol optical depth at 700nm in CRM grids', flag_xyfill=.true.)
   call addfld('CRM_AODVISZ',dims_crm_3D,'A', 'unitless', 'Aerosol optical depth at 500nm in CRM grids', flag_xyfill=.true.)
   call addfld('AOD400', horiz_only,'A', 'unitless', 'Aerosol optical depth at 400nm', flag_xyfill=.true.)
   call addfld('AOD700', horiz_only,'A', 'unitless', 'Aerosol optical depth at 700nm', flag_xyfill=.true.)

   !----------------------------------------------------------------------------
   ! 2-moment microphysics variables
   if (MMF_microphysics_scheme .eq. 'm2005') then
      call addfld('MMF_NC    ',(/'lev'/), 'A', '/kg',    'Cloud water dropet number from CRM')
      call addfld('MMF_NI    ',(/'lev'/), 'A', '/kg',    'Cloud ice crystal number from CRM')
      call addfld('MMF_NS    ',(/'lev'/), 'A', '/kg',    'Snow particle number from CRM')
      call addfld('MMF_NG    ',(/'lev'/), 'A', '/kg',    'Graupel particle number from CRM')
      call addfld('MMF_NR    ',(/'lev'/), 'A', '/kg',    'Rain particle number from CRM')
  
      call addfld('CRM_FLIQ ',dims_crm_3D, 'A', '1',   'Frequency of Occurrence of Liquid' )
      call addfld('CRM_FICE ',dims_crm_3D, 'A', '1',   'Frequency of Occurrence of Ice' )
      call addfld('CRM_FRAIN',dims_crm_3D, 'A', '1',   'Frequency of Occurrence of Rain' )
      call addfld('CRM_FSNOW',dims_crm_3D, 'A', '1',   'Frequency of Occurrence of Snow' )
      call addfld('CRM_FGRAP',dims_crm_3D, 'A', '1',   'Frequency of Occurrence of Graupel' )
  
      call addfld('CRM_QS  ',dims_crm_3D, 'A', 'kg/kg','Snow mixing ratio from CRM' )
      call addfld('CRM_QG  ',dims_crm_3D, 'A', 'kg/kg','Graupel mixing ratio from CRM' )
      call addfld('CRM_QR  ',dims_crm_3D, 'A', 'kg/kg','Rain mixing ratio from CRM' )
  
      call addfld('CRM_NC  ',dims_crm_3D, 'A', '/kg',  'Cloud water dropet number from CRM' )
      call addfld('CRM_NI  ',dims_crm_3D, 'A', '/kg',  'Cloud ice crystal number from CRM' )
      call addfld('CRM_NS  ',dims_crm_3D, 'A', '/kg',  'Snow particle number from CRM' )
      call addfld('CRM_NG  ',dims_crm_3D, 'A', '/kg',  'Graupel particle number from CRM' )
      call addfld('CRM_NR  ',dims_crm_3D, 'A', '/kg',  'Rain particle number from CRM' )
  
      ! below is for *instantaneous* crm output
      call addfld('CRM_AUT  ',dims_crm_3D, 'A', '/s',  'Autoconversion cloud waterfrom CRM' )
      call addfld('CRM_ACC  ',dims_crm_3D, 'A', '/s',  'Accretion cloud water from CRM' )
      call addfld('CRM_EVPC ',dims_crm_3D, 'A', '/s',  'Evaporation cloud water from CRM' )
      call addfld('CRM_EVPR ',dims_crm_3D, 'A', '/s',  'Evaporation rain from CRM' )
      call addfld('CRM_MLT  ',dims_crm_3D, 'A', '/s',  'Melting ice snow graupel from CRM' )
      call addfld('CRM_SUB  ',dims_crm_3D, 'A', '/s',  'Sublimation ice snow graupel from CRM' )
      call addfld('CRM_DEP  ',dims_crm_3D, 'A', '/s',  'Deposition ice snow graupel from CRM' )
      call addfld('CRM_CON  ',dims_crm_3D, 'A', '/s',  'Condensation cloud water from CRM' )
  
      ! *gcm-grid and time-step-avg* process output
      call addfld('A_AUT  ',(/'lev'/), 'A', '/s',      'Avg autoconversion cloud water from CRM' )
      call addfld('A_ACC  ',(/'lev'/), 'A', '/s',      'Avg accretion cloud water from CRM' )
      call addfld('A_EVPC ',(/'lev'/), 'A', '/s',      'Avg evaporation cloud water from CRM' )
      call addfld('A_EVPR ',(/'lev'/), 'A', '/s',      'Avg evaporation rain from CRM' )
      call addfld('A_MLT  ',(/'lev'/), 'A', '/s',      'Avg melting ice snow graupel from CRM' )
      call addfld('A_SUB  ',(/'lev'/), 'A', '/s',      'Avg sublimation ice snow graupel from CRM' )
      call addfld('A_DEP  ',(/'lev'/), 'A', '/s',      'Avg deposition ice snow graupel from CRM' )
      call addfld('A_CON  ',(/'lev'/), 'A', '/s',      'Avg condensation cloud water from CRM' )
  
      call addfld('CRM_DES  ', dims_crm_3D,'A','m^-6', 'cloud scale snow effective diameter')
      call addfld('CRM_MU   ', dims_crm_3D,'A','m^-6', 'cloud scale droplet size distribution shape parameter for radiation')
      call addfld('CRM_LAMBDA',dims_crm_3D,'A','m^-6', 'cloud scale slope of droplet distribution for radiation')
      call addfld('CRM_WVAR',  dims_crm_3D,'A','m/s',  'vertical velocity variance from CRM')
   end if  ! MMF_microphysics_scheme .eq. 'm2005'

   !----------------------------------------------------------------------------
   ! ECPP output variables
#ifdef ECPP
   if (use_ECPP) then
      call addfld('ABND    ',    ecpp_coord_i,'A','fraction','cloud fraction for each class for full time period at layer boundary')
      call addfld('ABND_TF ',    ecpp_coord_i,'A','fraction','cloud fraction for each class for end-portion of time period at layer boundary')
      call addfld('MASFBND ',    ecpp_coord_i,'A','kg/m2/s', 'sub-class vertical mass flux (kg/m2/s) at layer boundary')
      call addfld('ACEN    ',    ecpp_coord_m,'A','fraction','cloud fraction for each class for full time period at layer center')
      call addfld('ACEN_TF ',    ecpp_coord_m,'A','fraction','cloud fraction for each class for end-portion of time period at layer center')
      call addfld('RHCEN   ',    ecpp_coord_m,'A','fraction','relative humidity for each calss at layer center')
      call addfld('QCCEN   ',    ecpp_coord_m,'A','kg/kg',   'cloud water for each class at layer center')
      call addfld('QICEN   ',    ecpp_coord_m,'A','kg/kg',   'cloud ice for each class at layer center')
      call addfld('QSINK_AFCEN', ecpp_coord_m,'A','/s',      'cld water loss from precip after precip. for each class at layer center')
      call addfld('QSINK_BFCEN', ecpp_coord_m,'A','/s',      'cld water loss from precip before precip. for each class at layer center')
      call addfld('QSINK_AVGCEN',ecpp_coord_m,'A','/s',      'cld water loss from precip using avg cld water and precip. rate for each class at layer center')
      call addfld('PRAINCEN',    ecpp_coord_m,'A','kg/kg/s', 'cld water loss rate from precip (kg/kg/s) for each class at layer center')
      call addfld('PRECRCEN',    ecpp_coord_m,'A','kg/m2/s', 'liquid (rain) precipitation rate for each class at layer center')
      call addfld('PRECSCEN',    ecpp_coord_m,'A','kg/m2/s', 'solid ice precipitation rate for each class at layer center')
      call addfld('WUPTHRES',     (/'ilev'/), 'A','m/s',     'vertical velocity threshold for updraft')
      call addfld('WDNTHRES',     (/'ilev'/), 'A','m/s',     'vertical velocity threshold for dndraft')
      call addfld('WWQUI_CEN',    (/ 'lev'/), 'A','m2/s2',   'vertical velocity variance, quiescent class, layer center')
      call addfld('WWQUI_CLD_CEN',(/ 'lev'/), 'A','m2/s2',   'vertical velocity variance, cloudy quiescent class, layer center')
      call addfld('WWQUI_BND',    (/'ilev'/), 'A','m2/s2',   'vertical velocity variance, quiescent class, layer boundary')
      call addfld('WWQUI_CLD_BND',(/'ilev'/), 'A','m2/s2',   'vertical velocity variance, cloudy quiescent class, layer boundary')
   endif
#endif /* ECPP */

   !----------------------------------------------------------------------------
   ! MMF variables on the GCM grid
   call addfld('PRES',       (/'lev'/), 'A','Pa',     'Pressure' )
   call addfld('DPRES',      (/'lev'/), 'A','Pa',     'Pressure thickness of layer' )
   call addfld('MMF_DT',     (/'lev'/), 'A','K/s',    'T tendency due to CRM' )
   call addfld('MMF_DQ',     (/'lev'/), 'A','kg/kg/s','Q tendency due to CRM' )
   call addfld('MMF_DQC',    (/'lev'/), 'A','kg/kg/s','QC tendency due to CRM' )
   call addfld('MMF_DQI',    (/'lev'/), 'A','kg/kg/s','QI tendency due to CRM' )
   call addfld('MMF_MC',     (/'lev'/), 'A','kg/m2/s','Total mass flux from CRM' )
   call addfld('MMF_MCUP',   (/'lev'/), 'A','kg/m2/s','Updraft mass flux from CRM' )
   call addfld('MMF_MCDN',   (/'lev'/), 'A','kg/m2/s','Downdraft mass flux from CRM' )
   call addfld('MMF_MCUUP',  (/'lev'/), 'A','kg/m2/s','Unsaturated updraft mass flux from CRM' )
   call addfld('MMF_MCUDN',  (/'lev'/), 'A','kg/m2/s','Unsaturated downdraft mass flux from CRM' )
   call addfld('MU_CRM',     (/'lev'/), 'A','Pa/s',   'mass flux up from CRM')
   call addfld('MD_CRM',     (/'lev'/), 'A','Pa/s',   'mass flux down from CRM')
   call addfld('DU_CRM',     (/'lev'/), 'A','/s',     'detrainment from updraft from CRM')
   call addfld('EU_CRM',     (/'lev'/), 'A','/s',     'entraiment rate from updraft')
   call addfld('ED_CRM',     (/'lev'/), 'A','/s',     'entraiment rate from downdraft')
   call addfld('MMF_QC',     (/'lev'/), 'A','kg/kg',  'Cloud water from CRM' )
   call addfld('MMF_QI',     (/'lev'/), 'A','kg/kg',  'Cloud ice from CRM' )
   call addfld('MMF_QS',     (/'lev'/), 'A','kg/kg',  'Snow from CRM' )
   call addfld('MMF_QG',     (/'lev'/), 'A','kg/kg',  'Graupel from CRM' )
   call addfld('MMF_QR',     (/'lev'/), 'A','kg/kg',  'Rain from CRM' )
   call addfld('MMF_QTFLX',  (/'lev'/), 'A','kg/m2/s','Nonprecip. water flux from CRM' )
   call addfld('MMF_UFLX',   (/'lev'/), 'A','m2/s2',  'x-momentum flux from CRM' )
   call addfld('MMF_VFLX',   (/'lev'/), 'A','m2/s2',  'y-momentum flux from CRM' )
   call addfld('MMF_QTFLXS', (/'lev'/), 'A','kg/m2/s','SGS Nonprecip. water flux from CRM' )
   call addfld('MMF_TKEW',   (/'lev'/), 'A','kg/m/s2','vertical velocity variance in CRM')
   call addfld('MMF_TKE',    (/'lev'/), 'A','kg/m/s2','Total TKE in CRM' )
   call addfld('MMF_TKES',   (/'lev'/), 'A','kg/m/s2','SGS TKE in CRM' )
   call addfld('MMF_TK',     (/'lev'/), 'A','m2/s',   'SGS TK in CRM' )
   call addfld('MMF_QPFLX',  (/'lev'/), 'A','kg/m2/s','Precip. water flux from CRM' )
   call addfld('MMF_PFLX',   (/'lev'/), 'A','m/s',    'Precipitation flux from CRM' )
   call addfld('MMF_QTTR',   (/'lev'/), 'A','kg/kg/s','Nonprec. water transport from CRM' )
   call addfld('MMF_QPTR',   (/'lev'/), 'A','kg/kg/s','Prec. water transport from CRM' )
   call addfld('MMF_QPEVP',  (/'lev'/), 'A','kg/kg/s','Prec. water evaporation from CRM' )
   call addfld('MMF_QPFALL', (/'lev'/), 'A','kg/kg/s','Prec. water fall-out from CRM' )
   call addfld('MMF_QPSRC',  (/'lev'/), 'A','kg/kg/s','Prec. water source from CRM' )
   call addfld('MMF_QTLS',   (/'lev'/), 'A','kg/kg/s','L.S. Total Water Forcing in CRM' )
   call addfld('MMF_TLS',    (/'lev'/), 'A','kg/kg/s','L.S. LIWSE Forcing in CRM' )
   call addfld('MMF_TVFLUX', (/'lev'/), 'A', 'W/m2',  'Buoyancy Flux from CRM' )
   call addfld('MMF_BUOY',   (/'lev'/), 'A', 'W/m3',  'Buoyancy Term from CRM' )
   call addfld('MMF_BUOYSD', (/'lev'/), 'A', 'W/m3',  'Std Dev of Buoyancy Term from CRM' )
   call addfld('MMF_MSEF',   (/'lev'/), 'A', 'W/m2',  'Moist Static Energy Flux from CRM' )
   call addfld('MMF_QVFLUX', (/'lev'/), 'A', 'W/m2',  'Water Wapor Flux from CRM' )
   call addfld('MMF_QRL',    (/'lev'/), 'A','K/s',    'long-wave heating rate')
   call addfld('MMF_QRS',    (/'lev'/), 'A','K/s',    'short-wave heating rate')
   call addfld('MMF_CLDTOP', (/'lev'/), 'A',' ',      'Cloud Top PDF' )
   call addfld('MMF_DU',     (/'lev'/), 'A', 'm/s2 ','U tendency due to CRM' )
   call addfld('MMF_DV',     (/'lev'/), 'A', 'm/s2 ','V tendency due to CRM' )

   call addfld('MMF_SUBCYCLE_FAC', horiz_only,'A',' ', 'CRM subcycle ratio: 1.0 = no subcycling' )

   !----------------------------------------------------------------------------
   ! MMF CRM variance transport
   call addfld('MMF_VT_T'     ,(/'lev'/), 'A',' ','CRM T Variance')
   call addfld('MMF_VT_Q'     ,(/'lev'/), 'A',' ','CRM Q Variance')
   call addfld('MMF_VT_U'     ,(/'lev'/), 'A',' ','CRM U Variance')
   call addfld('MMF_VT_TEND_T',(/'lev'/), 'A',' ','CRM T Variance Tendency')
   call addfld('MMF_VT_TEND_Q',(/'lev'/), 'A',' ','CRM Q Variance Tendency')
   call addfld('MMF_VT_TEND_U',(/'lev'/), 'A',' ','CRM U Variance Tendency')
   call addfld('MMF_VT_TLS',   (/'lev'/), 'A','kg/kg/s','L.S. VT Forcing for LSE' )
   call addfld('MMF_VT_QLS',   (/'lev'/), 'A','kg/kg/s','L.S. VT Forcing for QT' )
   call addfld('MMF_VT_ULS',   (/'lev'/), 'A','kg/kg/s','L.S. VT Forcing for U' )

   !----------------------------------------------------------------------------
   ! mixing diagnostics for dropmixnuc in the GCM grid
   call addfld('MMF_NDROPMIX   ',(/'lev'/), 'A','#/kg/s','Droplet number mixing')
   call addfld('MMF_NDROPSRC   ',(/'lev'/), 'A','#/kg/s','Droplet number source')
   call addfld('MMF_NDROPCOL   ',horiz_only,'A','#/m2',  'Column droplet number')
   call addfld('MMF_NDROPLCLOUD',(/'lev'/), 'A',' ',     'Liquid cloud fraction')
   call addfld('MMF_NDROPKVH   ',(/'ilev'/),'A','m2/s',  'Vertical diffusivity in dropmixnuc for MMF')
   call addfld('MMF_NDROPWTKE  ',(/'lev'/), 'A',' ',     'Vertical velocity variance')

   !----------------------------------------------------------------------------
   ! add dropmixnuc tendencies for all modal aerosol species
#ifdef MODAL_AERO
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

            fieldname = 'MMF_' // trim(tmpname) // '_mixnuc1'
            long_name = trim(tmpname) // ' dropmixnuc mixnuc column tendency in the mmf one '
            call addfld( fieldname,  horiz_only, 'A', unit, long_name)

         end do ! lspec
      end do ! lphase
   end do ! m  

   do m = 1, pcnst
      if(species_class(m).eq.spec_class_gas) then
         fieldname = 'MMF_' // trim(cnst_name(m)) // '_mixnuc1'
         long_name = trim(cnst_name(m)) // ' dropmixnuc mixnuc column tendency in the mmf one '
         call addfld( fieldname,  horiz_only, 'A', unit, long_name)
      end if
   end do
#endif

   !----------------------------------------------------------------------------
   ! Set default MMF fields in monthly h0 files
   call add_default('MMF_DT    ', 1, ' ')
   call add_default('MMF_DQ    ', 1, ' ')
   call add_default('MMF_DQC   ', 1, ' ')
   call add_default('MMF_DQI   ', 1, ' ')
   call add_default('MMF_MC    ', 1, ' ')
   call add_default('MMF_MCUP  ', 1, ' ')
   call add_default('MMF_MCDN  ', 1, ' ')
   call add_default('MMF_MCUUP ', 1, ' ')
   call add_default('MMF_MCUDN ', 1, ' ')
   call add_default('MMF_TKE   ', 1, ' ')
   call add_default('MMF_TKES  ', 1, ' ')
   call add_default('MMF_TK    ', 1, ' ')
   call add_default('MMF_QTLS  ', 1, ' ')
   call add_default('MMF_TLS   ', 1, ' ')
   call add_default('MMF_SUBCYCLE_FAC', 1, ' ')

   if (MMF_microphysics_scheme .eq. 'm2005') then
      call add_default('MMF_NC    ', 1, ' ')
      call add_default('MMF_NI    ', 1, ' ')
      call add_default('MMF_NS    ', 1, ' ')
      call add_default('MMF_NG    ', 1, ' ')
      call add_default('MMF_NR    ', 1, ' ')
   end if

   !----------------------------------------------------------------------------

end subroutine crm_history_init
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
subroutine crm_history_out(state, ptend, crm_state, crm_rad, crm_output, & 
                           crm_ecpp_output, qrs, qrl, icol_beg, icol_end)
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
   use cam_abortutils,         only: endrun
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
   
   integer, intent(in) :: icol_beg ! CRM dimension index range
   integer, intent(in) :: icol_end ! CRM dimension index range

   !----------------------------------------------------------------------------
   ! local variables
   character*15 :: subname='crm_history_out'
   real(r8) :: cwp       (pcols,pver)  ! in-cloud cloud (total) water path (kg/m2)
   real(r8) :: gwp       (pcols,pver)  ! grid-box cloud (total) water path (kg/m2)
   real(r8) :: cicewp    (pcols,pver)  ! in-cloud cloud ice water path (kg/m2)
   real(r8) :: cliqwp    (pcols,pver)  ! in-cloud cloud liquid water path (kg/m2)
   real(r8) :: tgicewp   (pcols)       ! Vertically integrated ice water path (kg/m2
   real(r8) :: tgliqwp   (pcols)       ! Vertically integrated liquid water path (kg/m2)
   real(r8) :: tgwp      (pcols)       ! Vertically integrated (total) cloud water path  (kg/m2)
   real(r8) :: MMF_DT_out(pcols,pver)  ! CRM heating tendency
   integer :: lchnk                    ! chunk identifier
   integer :: ncol                     ! number of atmospheric columns
   integer :: ixcldliq, ixcldice       ! liquid and ice constituent indices
   integer :: i, k, icol               ! loop iterators
   logical :: use_ECPP                 ! explicit cloud parameterized pollutants
   logical :: use_MMF_VT               ! CRM variance transport
   character(len=16) :: MMF_microphysics_scheme
   integer :: idx_vt_t, idx_vt_q, idx_vt_u

   !----------------------------------------------------------------------------

   lchnk = state%lchnk
   ncol  = state%ncol

   ! check that icol bounds are correct
   if ( (icol_end-icol_beg+1) /= ncol ) then
      call endrun(trim(subname)//': icol_beg and icol_end bounds do not match ncol')
   end if

   call cnst_get_ind('CLDLIQ', ixcldliq)
   call cnst_get_ind('CLDICE', ixcldice)

   call phys_getopts(use_ECPP_out = use_ECPP)
   call phys_getopts(use_MMF_VT_out = use_MMF_VT)
   call phys_getopts(MMF_microphysics_scheme_out = MMF_microphysics_scheme)

   ! Subtract radiative heating for MMF_DT output
   MMF_DT_out(1:ncol,:pver) = ( ptend%s(1:ncol,:pver) - qrs(1:ncol,:pver) - qrl(1:ncol,:pver) )/cpair

   !----------------------------------------------------------------------------
   ! CRM domain average Tendencies
   call outfld('MMF_DT    ',MMF_DT_out(1:ncol,:),       ncol, lchnk ) 
   call outfld('MMF_DQ    ',ptend%q(1:ncol,1,1),        ncol, lchnk )
   call outfld('MMF_DQC   ',ptend%q(1:ncol,1,ixcldliq), ncol, lchnk )
   call outfld('MMF_DQI   ',ptend%q(1:ncol,1,ixcldice), ncol, lchnk )

   ! CRM radiative heating rate
   ! NOTE: We output the radiative heating rates (MMF_QRS and MMF_QRL) here 
   ! because this is the heating thatis applied to the CRM at this GCM timestep, 
   ! but note that this heating rate is offset in time from when the radiative 
   ! heating rate is calculated in the radiative transfer interface, which comes 
   ! AFTER this routine call (because we want it to use the clouds simulated 
   ! from this CRM call). Thus, comparing this heating rate with 
   ! CRM_QRS + CRM_QRL output in radiation_tend will show a time lag.

   ! GCM level rad heating tendencies
   call outfld('MMF_QRL   ',qrl(1:ncol,1:pver)/cpair, ncol, lchnk )
   call outfld('MMF_QRS   ',qrs(1:ncol,1:pver)/cpair, ncol, lchnk )

   ! Why do we output this here?
   call outfld('PRES    ',state%pmid(1:ncol,1:pver), ncol, lchnk )
   call outfld('DPRES   ',state%pdel(1:ncol,1:pver), ncol, lchnk )

   ! CRM state variables on CRM grid
   call outfld('CRM_U   ',crm_state%u_wind     (icol_beg:icol_end,:,:,:), ncol, lchnk )
   call outfld('CRM_V   ',crm_state%v_wind     (icol_beg:icol_end,:,:,:), ncol, lchnk )
   call outfld('CRM_W   ',crm_state%w_wind     (icol_beg:icol_end,:,:,:), ncol, lchnk )
   call outfld('CRM_T   ',crm_state%temperature(icol_beg:icol_end,:,:,:), ncol, lchnk )

   if (MMF_microphysics_scheme .eq. 'sam1mom') then
      call outfld('CRM_QV  ',(crm_state%qt(icol_beg:icol_end,:,:,:)    &
                              -crm_output%qcl(icol_beg:icol_end,:,:,:) &
                              -crm_output%qci(icol_beg:icol_end,:,:,:)), ncol, lchnk )
   else if (MMF_microphysics_scheme .eq. 'm2005') then 
      call outfld('CRM_QV  ', crm_state%qt(icol_beg:icol_end,:,:,:)    &
                              -crm_output%qcl(icol_beg:icol_end,:,:,:), ncol, lchnk )
   endif

   !----------------------------------------------------------------------------
   ! CRM condensate and precipitation on CRM grid
   call outfld('CRM_QC  ',crm_output%qcl     (icol_beg:icol_end,:,:,:),ncol, lchnk )
   call outfld('CRM_QI  ',crm_output%qci     (icol_beg:icol_end,:,:,:),ncol, lchnk )
   call outfld('CRM_QPC ',crm_output%qpl     (icol_beg:icol_end,:,:,:),ncol, lchnk )
   call outfld('CRM_QPI ',crm_output%qpi     (icol_beg:icol_end,:,:,:),ncol, lchnk )
   call outfld('CRM_PREC',crm_output%prec_crm(icol_beg:icol_end,:,:),  ncol, lchnk )
   call outfld('CRM_TK ', crm_output%tk      (icol_beg:icol_end,:,:,:),ncol, lchnk )  
   call outfld('CRM_TKH', crm_output%tkh     (icol_beg:icol_end,:,:,:),ncol, lchnk ) 

   !----------------------------------------------------------------------------
   ! CRM domain average condensate and precipitation
   call outfld('MMF_QC    ',crm_output%qc_mean(icol_beg:icol_end,:), ncol ,lchnk )
   call outfld('MMF_QI    ',crm_output%qi_mean(icol_beg:icol_end,:), ncol ,lchnk )
   call outfld('MMF_QS    ',crm_output%qs_mean(icol_beg:icol_end,:), ncol ,lchnk )
   call outfld('MMF_QG    ',crm_output%qg_mean(icol_beg:icol_end,:), ncol ,lchnk )
   call outfld('MMF_QR    ',crm_output%qr_mean(icol_beg:icol_end,:), ncol ,lchnk )

   !----------------------------------------------------------------------------
   ! CRM domain average fluxes
   call outfld('MMF_QTFLX ',crm_output%flux_qt   (icol_beg:icol_end,:), ncol, lchnk )
   call outfld('MMF_UFLX  ',crm_output%flux_u    (icol_beg:icol_end,:), ncol, lchnk )
   call outfld('MMF_VFLX  ',crm_output%flux_v    (icol_beg:icol_end,:), ncol, lchnk )
   call outfld('MMF_TKE   ',crm_output%tkez      (icol_beg:icol_end,:), ncol, lchnk )
   call outfld('MMF_TKEW  ',crm_output%tkew      (icol_beg:icol_end,:), ncol, lchnk )
   call outfld('MMF_TKES  ',crm_output%tkesgsz   (icol_beg:icol_end,:), ncol, lchnk )
   call outfld('MMF_TK    ',crm_output%tkz       (icol_beg:icol_end,:), ncol, lchnk )
   call outfld('MMF_QTFLXS',crm_output%fluxsgs_qt(icol_beg:icol_end,:), ncol, lchnk )
   call outfld('MMF_QPFLX ',crm_output%flux_qp   (icol_beg:icol_end,:), ncol, lchnk )
   call outfld('MMF_PFLX  ',crm_output%precflux  (icol_beg:icol_end,:), ncol, lchnk )
   call outfld('MMF_QTLS  ',crm_output%qt_ls     (icol_beg:icol_end,:), ncol, lchnk )
   call outfld('MMF_QTTR  ',crm_output%qt_trans  (icol_beg:icol_end,:), ncol, lchnk )
   call outfld('MMF_QPTR  ',crm_output%qp_trans  (icol_beg:icol_end,:), ncol, lchnk )
   call outfld('MMF_QPEVP ',crm_output%qp_evp    (icol_beg:icol_end,:), ncol, lchnk )
   call outfld('MMF_QPFALL',crm_output%qp_fall   (icol_beg:icol_end,:), ncol, lchnk )
   call outfld('MMF_QPSRC ',crm_output%qp_src    (icol_beg:icol_end,:), ncol, lchnk )
   call outfld('MMF_TLS   ',crm_output%t_ls      (icol_beg:icol_end,:), ncol, lchnk )

   ! NOTE: these should overwrite cloud outputs from non-MMF routines
   call outfld('CLOUD   ',crm_output%cld     (icol_beg:icol_end,:), ncol, lchnk )
   call outfld('CLDTOT  ',crm_output%cltot   (icol_beg:icol_end),   ncol, lchnk )
   call outfld('CLDHGH  ',crm_output%clhgh   (icol_beg:icol_end),   ncol, lchnk )
   call outfld('CLDMED  ',crm_output%clmed   (icol_beg:icol_end),   ncol, lchnk )
   call outfld('CLDLOW  ',crm_output%cllow   (icol_beg:icol_end),   ncol, lchnk )
   call outfld('MMF_CLDTOP',crm_output%cldtop(icol_beg:icol_end,:), ncol, lchnk )

   call outfld('MMF_SUBCYCLE_FAC',crm_output%subcycle_factor(icol_beg:icol_end), ncol,lchnk)

   !----------------------------------------------------------------------------
   ! CRM mass flux
   call outfld('MMF_MC    ', crm_output%mctot (icol_beg:icol_end,:), ncol, lchnk )
   call outfld('MMF_MCUP  ', crm_output%mcup  (icol_beg:icol_end,:), ncol, lchnk )
   call outfld('MMF_MCDN  ', crm_output%mcdn  (icol_beg:icol_end,:), ncol, lchnk )
   call outfld('MMF_MCUUP ', crm_output%mcuup (icol_beg:icol_end,:), ncol, lchnk )
   call outfld('MMF_MCUDN ', crm_output%mcudn (icol_beg:icol_end,:), ncol, lchnk )
   call outfld('MU_CRM    ', crm_output%mu_crm(icol_beg:icol_end,:), ncol, lchnk )
   call outfld('MD_CRM    ', crm_output%md_crm(icol_beg:icol_end,:), ncol, lchnk )
   call outfld('EU_CRM    ', crm_output%eu_crm(icol_beg:icol_end,:), ncol, lchnk )
   call outfld('DU_CRM    ', crm_output%du_crm(icol_beg:icol_end,:), ncol, lchnk )
   call outfld('ED_CRM    ', crm_output%ed_crm(icol_beg:icol_end,:), ncol, lchnk )

#ifdef m2005
   if (MMF_microphysics_scheme .eq. 'm2005') then
      ! index is defined in MICRO_M2005/microphysics.F90
      ! Be cautious to use them here. They are defined in crm codes, and these codes are called only 
      ! after the subroutine of crm is called. So they can only be used after the 'crm' subroutine. 
      ! incl, inci, ... can not be used here, for they are defined before we call them???
      call outfld('CRM_NC ',crm_state%nc(icol_beg:icol_end,:,:,:), ncol, lchnk )
      call outfld('CRM_NI ',crm_state%ni(icol_beg:icol_end,:,:,:), ncol, lchnk )
      call outfld('CRM_NR ',crm_state%nr(icol_beg:icol_end,:,:,:), ncol, lchnk )
      call outfld('CRM_NS ',crm_state%ns(icol_beg:icol_end,:,:,:), ncol, lchnk )
      call outfld('CRM_NG ',crm_state%ng(icol_beg:icol_end,:,:,:), ncol, lchnk )
      call outfld('CRM_QR ',crm_state%qr(icol_beg:icol_end,:,:,:), ncol, lchnk )
      call outfld('CRM_QS ',crm_state%qs(icol_beg:icol_end,:,:,:), ncol, lchnk )
      call outfld('CRM_QG ',crm_state%qg(icol_beg:icol_end,:,:,:), ncol, lchnk )
      
      call outfld('CRM_WVAR',crm_output%wvar(icol_beg:icol_end,:,:,:), ncol, lchnk)

      call outfld('CRM_AUT', crm_output%aut (icol_beg:icol_end,:,:,:), ncol, lchnk)
      call outfld('CRM_ACC', crm_output%acc (icol_beg:icol_end,:,:,:), ncol, lchnk)
      call outfld('CRM_MLT', crm_output%mlt (icol_beg:icol_end,:,:,:), ncol, lchnk)
      call outfld('CRM_SUB', crm_output%sub (icol_beg:icol_end,:,:,:), ncol, lchnk)
      call outfld('CRM_DEP', crm_output%dep (icol_beg:icol_end,:,:,:), ncol, lchnk)
      call outfld('CRM_CON', crm_output%con (icol_beg:icol_end,:,:,:), ncol, lchnk)
      call outfld('CRM_EVPC',crm_output%evpc(icol_beg:icol_end,:,:,:), ncol, lchnk)
      call outfld('CRM_EVPR',crm_output%evpr(icol_beg:icol_end,:,:,:), ncol, lchnk)
      
      call outfld('A_AUT', crm_output%aut_a (icol_beg:icol_end,:), ncol, lchnk)
      call outfld('A_ACC', crm_output%acc_a (icol_beg:icol_end,:), ncol, lchnk)
      call outfld('A_MLT', crm_output%mlt_a (icol_beg:icol_end,:), ncol, lchnk)
      call outfld('A_SUB', crm_output%sub_a (icol_beg:icol_end,:), ncol, lchnk)
      call outfld('A_DEP', crm_output%dep_a (icol_beg:icol_end,:), ncol, lchnk)
      call outfld('A_CON', crm_output%con_a (icol_beg:icol_end,:), ncol, lchnk)
      call outfld('A_EVPC',crm_output%evpc_a(icol_beg:icol_end,:), ncol, lchnk)
      call outfld('A_EVPR',crm_output%evpr_a(icol_beg:icol_end,:), ncol, lchnk)

      call outfld('MMF_NC    ',crm_output%nc_mean(icol_beg:icol_end,:), ncol, lchnk )
      call outfld('MMF_NI    ',crm_output%ni_mean(icol_beg:icol_end,:), ncol, lchnk )
      call outfld('MMF_NS    ',crm_output%ns_mean(icol_beg:icol_end,:), ncol, lchnk )
      call outfld('MMF_NG    ',crm_output%ng_mean(icol_beg:icol_end,:), ncol, lchnk )
      call outfld('MMF_NR    ',crm_output%nr_mean(icol_beg:icol_end,:), ncol, lchnk )
   endif ! m2005
#endif /* m2005 */

   !----------------------------------------------------------------------------
   ! Compute liquid water paths (for diagnostics only)
   tgicewp(1:ncol) = 0.
   tgliqwp(1:ncol) = 0.
   do k = 1,pver
      do i = 1,ncol
         icol = icol_beg - 1 + i
         cicewp(i,k) = crm_output%gicewp(icol,k) * 1.0e-3 / max(0.01_r8,crm_output%cld(icol,k)) ! In-cloud ice water path.  g/m2 --> kg/m2
         cliqwp(i,k) = crm_output%gliqwp(icol,k) * 1.0e-3 / max(0.01_r8,crm_output%cld(icol,k)) ! In-cloud liquid water path. g/m2 --> kg/m2
         tgicewp(i)  = tgicewp(i) + crm_output%gicewp(icol,k) *1.0e-3 ! grid cell mean ice water path.  g/m2 --> kg/m2
         tgliqwp(i)  = tgliqwp(i) + crm_output%gliqwp(icol,k) *1.0e-3 ! grid cell mean ice water path.  g/m2 --> kg/m2
      end do
   end do
   tgwp(1:ncol) = tgicewp(1:ncol) + tgliqwp(1:ncol)
   gwp(1:ncol,:pver) = crm_output%gicewp(icol_beg:icol_end,:pver) + crm_output%gliqwp(icol_beg:icol_end,:pver)
   cwp(1:ncol,:pver) = cicewp(1:ncol,:pver) + cliqwp(1:ncol,:pver)

   call outfld('GCLDLWP' ,gwp(1:ncol,1:pver),    ncol, lchnk)
   call outfld('TGCLDCWP',tgwp(1:ncol),          ncol, lchnk)
   call outfld('TGCLDLWP',tgliqwp(1:ncol),       ncol, lchnk)
   call outfld('TGCLDIWP',tgicewp(1:ncol),       ncol, lchnk)
   call outfld('ICLDTWP' ,cwp(1:ncol,1:pver),    ncol, lchnk)
   call outfld('ICLDIWP' ,cicewp(1:ncol,1:pver), ncol, lchnk)

   !----------------------------------------------------------------------------
   ! ECPP
#ifdef ECPP
   if (use_ECPP) then
      call outfld('ACEN    ',      crm_ecpp_output%acen            (icol_beg:icol_end,:,:,:,:), ncol, lchnk )
      call outfld('ABND    ',      crm_ecpp_output%abnd            (icol_beg:icol_end,:,:,:,:), ncol, lchnk )
      call outfld('ACEN_TF ',      crm_ecpp_output%acen_tf         (icol_beg:icol_end,:,:,:,:), ncol, lchnk )
      call outfld('ABND_TF ',      crm_ecpp_output%abnd_tf         (icol_beg:icol_end,:,:,:,:), ncol, lchnk )
      call outfld('MASFBND ',      crm_ecpp_output%massflxbnd      (icol_beg:icol_end,:,:,:,:), ncol, lchnk )
      call outfld('RHCEN   ',      crm_ecpp_output%rhcen           (icol_beg:icol_end,:,:,:,:), ncol, lchnk )
      call outfld('QCCEN   ',      crm_ecpp_output%qcloudcen       (icol_beg:icol_end,:,:,:,:), ncol, lchnk )
      call outfld('QICEN   ',      crm_ecpp_output%qicecen         (icol_beg:icol_end,:,:,:,:), ncol, lchnk )
      call outfld('QSINK_AFCEN',   crm_ecpp_output%qlsink_afcen    (icol_beg:icol_end,:,:,:,:), ncol, lchnk )
      call outfld('PRECRCEN',      crm_ecpp_output%precrcen        (icol_beg:icol_end,:,:,:,:), ncol, lchnk )
      call outfld('PRECSCEN',      crm_ecpp_output%precsolidcen    (icol_beg:icol_end,:,:,:,:), ncol, lchnk )
      call outfld('WUPTHRES',      crm_ecpp_output%wupthresh_bnd   (icol_beg:icol_end,:),       ncol, lchnk )
      call outfld('WDNTHRES',      crm_ecpp_output%wdownthresh_bnd (icol_beg:icol_end,:),       ncol, lchnk )
      call outfld('WWQUI_CEN',     crm_ecpp_output%wwqui_cen       (icol_beg:icol_end,:),       ncol, lchnk )
      call outfld('WWQUI_CLD_CEN', crm_ecpp_output%wwqui_cloudy_cen(icol_beg:icol_end,:),       ncol, lchnk )
      call outfld('WWQUI_BND',     crm_ecpp_output%wwqui_cen       (icol_beg:icol_end,:),       ncol, lchnk )
      call outfld('WWQUI_CLD_BND', crm_ecpp_output%wwqui_cloudy_cen(icol_beg:icol_end,:),       ncol, lchnk )
      call outfld('QSINK_BFCEN',   crm_ecpp_output%qlsink_bfcen    (icol_beg:icol_end,:,:,:,:), ncol, lchnk )
      call outfld('QSINK_AVGCEN',  crm_ecpp_output%qlsink_avgcen   (icol_beg:icol_end,:,:,:,:), ncol, lchnk )
      call outfld('PRAINCEN',      crm_ecpp_output%praincen        (icol_beg:icol_end,:,:,:,:), ncol, lchnk )
   end if ! use_ECPP
#endif /* ECPP */

   !----------------------------------------------------------------------------
   ! CRM momentum tendencies
   call outfld('MMF_DU',ptend%u(1:ncol,:), ncol, lchnk )
   call outfld('MMF_DV',ptend%v(1:ncol,:), ncol, lchnk )

   !----------------------------------------------------------------------------
   ! CRM variance transport tracer values and tendencies
   if (use_MMF_VT) then
      call cnst_get_ind( 'VT_T', idx_vt_t )
      call cnst_get_ind( 'VT_Q', idx_vt_q )
      call cnst_get_ind( 'VT_U', idx_vt_u )
      call outfld('MMF_VT_T',      state%q(:,:,idx_vt_t),                   ncol, lchnk )
      call outfld('MMF_VT_Q',      state%q(:,:,idx_vt_q),                   ncol, lchnk )
      call outfld('MMF_VT_U',      state%q(:,:,idx_vt_u),                   ncol, lchnk )
      call outfld('MMF_VT_TEND_T', ptend%q(:,:,idx_vt_t),                   ncol, lchnk )
      call outfld('MMF_VT_TEND_Q', ptend%q(:,:,idx_vt_q),                   ncol, lchnk )
      call outfld('MMF_VT_TEND_U', ptend%q(:,:,idx_vt_u),                   ncol, lchnk )
      call outfld('MMF_VT_TLS',    crm_output%t_vt_ls(icol_beg:icol_end,:), ncol, lchnk )
      call outfld('MMF_VT_QLS',    crm_output%q_vt_ls(icol_beg:icol_end,:), ncol, lchnk )
      call outfld('MMF_VT_ULS',    crm_output%u_vt_ls(icol_beg:icol_end,:), ncol, lchnk )
   end if

   !----------------------------------------------------------------------------

end subroutine crm_history_out
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
end module crm_history
