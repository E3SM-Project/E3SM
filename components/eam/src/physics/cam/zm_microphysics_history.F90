module  zm_microphysics_history
   !----------------------------------------------------------------------------
   ! Purpose: microphysics state structure definition and methods for ZM
   ! Original Author: Xialiang Song and Guang Zhang, June 2010
   !----------------------------------------------------------------------------
   use shr_kind_mod,          only: r8=>shr_kind_r8
   use ppgrid,                only: pcols, pver, pverp
   use zm_microphysics_state, only: zm_microp_st

   public :: zm_microphysics_history_init ! add fields for history output
   public :: zm_microphysics_history_out  ! write history output related to ZM microphysics
  
!===================================================================================================
contains
!===================================================================================================

subroutine zm_microphysics_history_init()
   !----------------------------------------------------------------------------
   ! Purpose: add output history variables for convective microphysics
   !----------------------------------------------------------------------------
   use cam_history, only: addfld, horiz_only, add_default
   !----------------------------------------------------------------------------
   call addfld( 'CLDLIQZM',(/ 'lev' /), 'A', 'g/m3',     'ZM cloud liq water')
   call addfld( 'CLDICEZM',(/ 'lev' /), 'A', 'g/m3',     'ZM cloud ice water')
   call addfld( 'CLIQSNUM',(/ 'lev' /), 'A', '1',        'ZM cloud liq water sample number')
   call addfld( 'CICESNUM',(/ 'lev' /), 'A', '1',        'ZM cloud ice water sample number')
   call addfld( 'QRAINZM' ,(/ 'lev' /), 'A', 'g/m3',     'ZM rain water')
   call addfld( 'QSNOWZM' ,(/ 'lev' /), 'A', 'g/m3',     'ZM snow')
   call addfld( 'QGRAPZM' ,(/ 'lev' /), 'A', 'g/m3',     'ZM graupel')
   call addfld( 'CRAINNUM',(/ 'lev' /), 'A', '1',        'ZM cloud rain water sample number')
   call addfld( 'CSNOWNUM',(/ 'lev' /), 'A', '1',        'ZM cloud snow sample number')
   call addfld( 'CGRAPNUM',(/ 'lev' /), 'A', '1',        'ZM cloud graupel sample number')
   call addfld( 'DIFZM',   (/ 'lev' /), 'A', 'kg/kg/s ', 'ZM detrained ice water')
   call addfld( 'DLFZM',   (/ 'lev' /), 'A', 'kg/kg/s ', 'ZM detrained liq water')
   call addfld( 'DNIFZM',  (/ 'lev' /), 'A', '1/kg/s ',  'ZM detrained ice water num concen')
   call addfld( 'DNLFZM',  (/ 'lev' /), 'A', '1/kg/s ',  'ZM detrained liquid water num concen')
   call addfld( 'WUZM',    (/ 'lev' /), 'A', 'm/s',      'ZM vertical velocity')
   call addfld( 'WUZMSNUM',(/ 'lev' /), 'A', '1',        'ZM vertical velocity sample number')
   call addfld( 'QNLZM',   (/ 'lev' /), 'A', '1/m3',     'ZM cloud liq water number concen')
   call addfld( 'QNIZM',   (/ 'lev' /), 'A', '1/m3',     'ZM cloud ice number concen')
   call addfld( 'QNRZM',   (/ 'lev' /), 'A', '1/m3',     'ZM cloud rain water number concen')
   call addfld( 'QNSZM',   (/ 'lev' /), 'A', '1/m3',     'ZM cloud snow number concen')
   call addfld( 'QNGZM',   (/ 'lev' /), 'A', '1/m3',     'ZM cloud graupel number concen')
   call addfld( 'FRZZM',   (/ 'lev' /), 'A', 'K/s',      'ZM heating tendency due to freezing')
   call addfld( 'AUTOL_M', (/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency due to autoconversion of droplets to rain')
   call addfld( 'ACCRL_M', (/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency due to accretion of droplets by rain')
   call addfld( 'BERGN_M', (/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency due to Bergeron process')
   call addfld( 'FHTIM_M', (/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency due to immersion freezing')
   call addfld( 'FHTCT_M', (/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency due to contact freezing')
   call addfld( 'FHML_M',  (/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency due to homogeneous freezing of droplet')
   call addfld( 'HMPI_M',  (/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency due to HM process')
   call addfld( 'ACCSL_M', (/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency due to accretion of droplet by snow')
   call addfld( 'DLF_M',   (/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency due to detrainment of droplet')
   call addfld( 'COND_M',  (/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency due to condensation')
   call addfld( 'AUTOL_N', (/ 'lev' /), 'A', '1/kg/m',   'ZM num tendency due to autoconversion of droplets to rain')
   call addfld( 'ACCRL_N', (/ 'lev' /), 'A', '1/kg/m',   'ZM num tendency due to accretion of droplets by rain')
   call addfld( 'BERGN_N', (/ 'lev' /), 'A', '1/kg/m',   'ZM num tendency due to Bergeron process')
   call addfld( 'FHTIM_N', (/ 'lev' /), 'A', '1/kg/m',   'ZM num tendency due to immersion freezing')
   call addfld( 'FHTCT_N', (/ 'lev' /), 'A', '1/kg/m',   'ZM num tendency due to contact freezing')
   call addfld( 'FHML_N',  (/ 'lev' /), 'A', '1/kg/m',   'ZM num tendency due to homogeneous freezing of droplet')
   call addfld( 'ACCSL_N', (/ 'lev' /), 'A', '1/kg/m',   'ZM num tendency due to accretion of droplet by snow')
   call addfld( 'ACTIV_N', (/ 'lev' /), 'A', '1/kg/m',   'ZM num tendency due to droplets activation')
   call addfld( 'DLF_N',   (/ 'lev' /), 'A', '1/kg/m',   'ZM num tendency due to detrainment of droplet')
   call addfld( 'AUTOI_M', (/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency due to autoconversion of ice to snow')
   call addfld( 'ACCSI_M', (/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency due to accretion of ice by snow')
   call addfld( 'DIF_M',   (/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency due to detrainment of cloud ice')
   call addfld( 'DEPOS_M', (/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency due to deposition')
   call addfld( 'NUCLI_N', (/ 'lev' /), 'A', '1/kg/m' ,  'ZM num tendency due to ice nucleation')
   call addfld( 'AUTOI_N', (/ 'lev' /), 'A', '1/kg/m' ,  'ZM num tendency due to autoconversion of ice to snow')
   call addfld( 'ACCSI_N', (/ 'lev' /), 'A', '1/kg/m' ,  'ZM num tendency due to accretion of ice by snow')
   call addfld( 'HMPI_N',  (/ 'lev' /), 'A', '1/kg/s' ,  'ZM num tendency due to HM process')
   call addfld( 'DIF_N',   (/ 'lev' /), 'A', '1/kg/m' ,  'ZM num tendency due to detrainment of cloud ice')
   call addfld( 'TRSPC_M', (/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency of droplets due to convective transport')
   call addfld( 'TRSPC_N', (/ 'lev' /), 'A', '1/kg/m' ,  'ZM num tendency of droplets due to convective transport')
   call addfld( 'TRSPI_M', (/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency of ice crystal due to convective transport')
   call addfld( 'TRSPI_N', (/ 'lev' /), 'A', '1/kg/m' ,  'ZM num tendency of ice crystal due to convective transport')
   call addfld( 'ACCGR_M', (/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency due to collection of rain by graupel')
   call addfld( 'ACCGL_M', (/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency due to collection of droplets by graupel')
   call addfld( 'ACCGSL_M',(/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency of graupel due to collection of droplets by snow')
   call addfld( 'ACCGSR_M',(/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency of graupel due to collection of rain by snow')
   call addfld( 'ACCGIR_M',(/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency of graupel due to collection of rain by ice')
   call addfld( 'ACCGRI_M',(/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency of graupel due to collection of ice by rain')
   call addfld( 'ACCGRS_M',(/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency of graupel due to collection of snow by rain')
   call addfld( 'ACCGSL_N',(/ 'lev' /), 'A', '1/kg/m' ,  'ZM num tendency of graupel due to collection of droplets by snow')
   call addfld( 'ACCGSR_N',(/ 'lev' /), 'A', '1/kg/m' ,  'ZM num tendency of graupel due to collection of rain by snow')
   call addfld( 'ACCGIR_N',(/ 'lev' /), 'A', '1/kg/m' ,  'ZM num tendency of graupel due to collection of rain by ice')
   call addfld( 'ACCSRI_M',(/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency of snow due to collection of ice by rain')
   call addfld( 'ACCIGL_M',(/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency of ice mult(splintering) due to acc droplets by graupel')
   call addfld( 'ACCIGR_M',(/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency of ice mult(splintering) due to acc rain by graupel')
   call addfld( 'ACCSIR_M',(/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency of snow due to collection of rain by ice')
   call addfld( 'ACCIGL_N',(/ 'lev' /), 'A', '1/kg/m' ,  'ZM num tendency of ice mult(splintering) due to acc droplets by graupel')
   call addfld( 'ACCIGR_N',(/ 'lev' /), 'A', '1/kg/m' ,  'ZM num tendency of ice mult(splintering) due to acc rain by graupel')
   call addfld( 'ACCSIR_N',(/ 'lev' /), 'A', '1/kg/m' ,  'ZM num tendency of snow due to collection of rain by ice')
   call addfld( 'ACCGL_N', (/ 'lev' /), 'A', '1/kg/m' ,  'ZM num tendency due to collection of droplets by graupel')
   call addfld( 'ACCGR_N', (/ 'lev' /), 'A', '1/kg/m' ,  'ZM num tendency due to collection of rain by graupel')
   call addfld( 'ACCIL_M', (/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency of cloud ice due to collection of droplet by cloud ice')
   call addfld( 'ACCIL_N', (/ 'lev' /), 'A', '1/kg/m' ,  'ZM num tendency of cloud ice due to collection of droplet by cloud ice')
   call addfld( 'FALLR_M', (/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency of rain fallout')
   call addfld( 'FALLS_M', (/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency of snow fallout')
   call addfld( 'FALLG_M', (/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency of graupel fallout')
   call addfld( 'FALLR_N', (/ 'lev' /), 'A', '1/kg/m' ,  'ZM num tendency of rain fallout')
   call addfld( 'FALLS_N', (/ 'lev' /), 'A', '1/kg/m' ,  'ZM num tendency of snow fallout')
   call addfld( 'FALLG_N', (/ 'lev' /), 'A', '1/kg/m' ,  'ZM num tendency of graupel fallout')
   call addfld( 'FHMR_M',  (/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency due to homogeneous freezing of rain')
   call addfld( 'PRECZ_SN',horiz_only , 'A', '#',        'ZM sample num of convective precipitation rate')

   !----------------------------------------------------------------------------

   call add_default( 'CLDLIQZM', 1, ' ')
   call add_default( 'CLDICEZM', 1, ' ')
   call add_default( 'CLIQSNUM', 1, ' ')
   call add_default( 'CICESNUM', 1, ' ')
   call add_default( 'DIFZM',    1, ' ')
   call add_default( 'DLFZM',    1, ' ')
   call add_default( 'DNIFZM',   1, ' ')
   call add_default( 'DNLFZM',   1, ' ')
   call add_default( 'WUZM',     1, ' ')
   call add_default( 'QRAINZM',  1, ' ')
   call add_default( 'QSNOWZM',  1, ' ')
   call add_default( 'QGRAPZM',  1, ' ')
   call add_default( 'CRAINNUM', 1, ' ')
   call add_default( 'CSNOWNUM', 1, ' ')
   call add_default( 'CGRAPNUM', 1, ' ')
   call add_default( 'QNLZM',    1, ' ')
   call add_default( 'QNIZM',    1, ' ')
   call add_default( 'QNRZM',    1, ' ')
   call add_default( 'QNSZM',    1, ' ')
   call add_default( 'QNGZM',    1, ' ')
   call add_default( 'FRZZM',    1, ' ')

end subroutine zm_microphysics_history_init

!===================================================================================================

subroutine zm_microphysics_history_out( lchnk, ncol, microp_st, prec, dlf, dif, dnlf, dnif, frz )
   !----------------------------------------------------------------------------
   ! Purpose: write out history variables for convective microphysics
   !----------------------------------------------------------------------------
   use cam_history, only: outfld
   !----------------------------------------------------------------------------
   ! Arguments
   integer,                         intent(in) :: lchnk     ! chunk identifier
   integer,                         intent(in) :: ncol      ! number of columns in chunk
   type(zm_microp_st),              intent(in) :: microp_st ! ZM microphysics data structure
   real(r8), dimension(pcols),      intent(in) :: prec      ! convective precip rate
   real(r8), dimension(pcols,pver), intent(in) :: dlf       ! detrainment of conv cld liq water mixing ratio
   real(r8), dimension(pcols,pver), intent(in) :: dif       ! detrainment of conv cld ice mixing ratio
   real(r8), dimension(pcols,pver), intent(in) :: dnlf      ! detrainment of conv cld liq water num concen
   real(r8), dimension(pcols,pver), intent(in) :: dnif      ! detrainment of conv cld ice num concen
   real(r8), dimension(pcols,pver), intent(in) :: frz       ! heating rate due to freezing
   !----------------------------------------------------------------------------
   ! Local variables
   integer  :: i,k
   real(r8), dimension(pcols)      :: precz_snum      ! sample num of conv precip rate
   real(r8), dimension(pcols,pver) :: cice_snum       ! convective cloud ice sample number
   real(r8), dimension(pcols,pver) :: cliq_snum       ! convective cloud liquid sample number
   real(r8), dimension(pcols,pver) :: crain_snum      ! convective rain water sample number
   real(r8), dimension(pcols,pver) :: csnow_snum      ! convective snow sample number
   real(r8), dimension(pcols,pver) :: cgraupel_snum   ! convective graupel sample number
   real(r8), dimension(pcols,pver) :: wu_snum         ! vertical velocity sample number
   !----------------------------------------------------------------------------
   do k = 1,pver
      do i = 1,ncol
         if (microp_st%qice(i,k)     >  0) cice_snum(i,k)     = 1
         if (microp_st%qice(i,k)     <= 0) cice_snum(i,k)     = 0
         if (microp_st%qliq(i,k)     >  0) cliq_snum(i,k)     = 1
         if (microp_st%qliq(i,k)     <= 0) cliq_snum(i,k)     = 0
         if (microp_st%qsnow(i,k)    >  0) csnow_snum(i,k)    = 1
         if (microp_st%qsnow(i,k)    <= 0) csnow_snum(i,k)    = 0
         if (microp_st%qrain(i,k)    >  0) crain_snum(i,k)    = 1
         if (microp_st%qrain(i,k)    <= 0) crain_snum(i,k)    = 0
         if (microp_st%qgraupel(i,k) >  0) cgraupel_snum(i,k) = 1
         if (microp_st%qgraupel(i,k) <= 0) cgraupel_snum(i,k) = 0
         if (microp_st%wu(i,k)       >  0) wu_snum(i,k)       = 1
         if (microp_st%wu(i,k)       <= 0) wu_snum(i,k)       = 0
      end do
   end do

   call outfld('CLIQSNUM',cliq_snum          , pcols, lchnk )
   call outfld('CICESNUM',cice_snum          , pcols, lchnk )
   call outfld('CRAINNUM',crain_snum         , pcols, lchnk )
   call outfld('CSNOWNUM',csnow_snum         , pcols, lchnk )
   call outfld('CGRAPNUM',cgraupel_snum      , pcols, lchnk )
   call outfld('WUZMSNUM',wu_snum            , pcols, lchnk )

   call outfld('DIFZM'   ,dif                , pcols, lchnk )
   call outfld('DLFZM'   ,dlf                , pcols, lchnk )
   call outfld('DNIFZM'  ,dnif               , pcols, lchnk )
   call outfld('DNLFZM'  ,dnlf               , pcols, lchnk )
   call outfld('FRZZM'   ,frz                , pcols, lchnk )

   call outfld('WUZM'    ,microp_st%wu       , pcols, lchnk )

   call outfld('CLDLIQZM',microp_st%qliq     , pcols, lchnk )
   call outfld('CLDICEZM',microp_st%qice     , pcols, lchnk )
   call outfld('QRAINZM' ,microp_st%qrain    , pcols, lchnk )
   call outfld('QSNOWZM' ,microp_st%qsnow    , pcols, lchnk )
   call outfld('QGRAPZM' ,microp_st%qgraupel , pcols, lchnk )

   call outfld('QNLZM'   ,microp_st%qnl      , pcols, lchnk )
   call outfld('QNIZM'   ,microp_st%qni      , pcols, lchnk )
   call outfld('QNRZM'   ,microp_st%qnr      , pcols, lchnk )
   call outfld('QNSZM'   ,microp_st%qns      , pcols, lchnk )
   call outfld('QNGZM'   ,microp_st%qng      , pcols, lchnk )

   call outfld('AUTOL_M' ,microp_st%autolm   , pcols, lchnk )
   call outfld('ACCRL_M' ,microp_st%accrlm   , pcols, lchnk )
   call outfld('BERGN_M' ,microp_st%bergnm   , pcols, lchnk )
   call outfld('FHTIM_M' ,microp_st%fhtimm   , pcols, lchnk )
   call outfld('FHTCT_M' ,microp_st%fhtctm   , pcols, lchnk )
   call outfld('FHML_M'  ,microp_st%fhmlm    , pcols, lchnk )
   call outfld('HMPI_M'  ,microp_st%hmpim    , pcols, lchnk )
   call outfld('ACCSL_M' ,microp_st%accslm   , pcols, lchnk )
   call outfld('DLF_M'   ,microp_st%dlfm     , pcols, lchnk )

   call outfld('AUTOL_N' ,microp_st%autoln   , pcols, lchnk )
   call outfld('ACCRL_N' ,microp_st%accrln   , pcols, lchnk )
   call outfld('BERGN_N' ,microp_st%bergnn   , pcols, lchnk )
   call outfld('FHTIM_N' ,microp_st%fhtimn   , pcols, lchnk )
   call outfld('FHTCT_N' ,microp_st%fhtctn   , pcols, lchnk )
   call outfld('FHML_N'  ,microp_st%fhmln    , pcols, lchnk )
   call outfld('ACCSL_N' ,microp_st%accsln   , pcols, lchnk )
   call outfld('ACTIV_N' ,microp_st%activn   , pcols, lchnk )
   call outfld('DLF_N'   ,microp_st%dlfn     , pcols, lchnk )
   call outfld('AUTOI_M' ,microp_st%autoim   , pcols, lchnk )
   call outfld('ACCSI_M' ,microp_st%accsim   , pcols, lchnk )
   call outfld('DIF_M'   ,microp_st%difm     , pcols, lchnk )
   call outfld('NUCLI_N' ,microp_st%nuclin   , pcols, lchnk )
   call outfld('AUTOI_N' ,microp_st%autoin   , pcols, lchnk )
   call outfld('ACCSI_N' ,microp_st%accsin   , pcols, lchnk )
   call outfld('HMPI_N'  ,microp_st%hmpin    , pcols, lchnk )
   call outfld('DIF_N'   ,microp_st%difn     , pcols, lchnk )
   call outfld('COND_M'  ,microp_st%cmel     , pcols, lchnk )
   call outfld('DEPOS_M' ,microp_st%cmei     , pcols, lchnk )

   call outfld('TRSPC_M' ,microp_st%trspcm   , pcols, lchnk )
   call outfld('TRSPC_N' ,microp_st%trspcn   , pcols, lchnk )
   call outfld('TRSPI_M' ,microp_st%trspim   , pcols, lchnk )
   call outfld('TRSPI_N' ,microp_st%trspin   , pcols, lchnk )

   call outfld('ACCGR_M' ,microp_st%accgrm   , pcols, lchnk )
   call outfld('ACCGL_M' ,microp_st%accglm   , pcols, lchnk )
   call outfld('ACCGSL_M',microp_st%accgslm  , pcols, lchnk )
   call outfld('ACCGSR_M',microp_st%accgsrm  , pcols, lchnk )
   call outfld('ACCGIR_M',microp_st%accgirm  , pcols, lchnk )
   call outfld('ACCGRI_M',microp_st%accgrim  , pcols, lchnk )
   call outfld('ACCGRS_M',microp_st%accgrsm  , pcols, lchnk )

   call outfld('ACCGSL_N',microp_st%accgsln  , pcols, lchnk )
   call outfld('ACCGSR_N',microp_st%accgsrn  , pcols, lchnk )
   call outfld('ACCGIR_N',microp_st%accgirn  , pcols, lchnk )

   call outfld('ACCSRI_M',microp_st%accsrim  , pcols, lchnk )
   call outfld('ACCIGL_M',microp_st%acciglm  , pcols, lchnk )
   call outfld('ACCIGR_M',microp_st%accigrm  , pcols, lchnk )
   call outfld('ACCSIR_M',microp_st%accsirm  , pcols, lchnk )

   call outfld('ACCIGL_N',microp_st%accigln  , pcols, lchnk )
   call outfld('ACCIGR_N',microp_st%accigrn  , pcols, lchnk )
   call outfld('ACCSIR_N',microp_st%accsirn  , pcols, lchnk )
   call outfld('ACCGL_N' ,microp_st%accgln   , pcols, lchnk )
   call outfld('ACCGR_N' ,microp_st%accgrn   , pcols, lchnk )

   call outfld('ACCIL_M' ,microp_st%accilm   , pcols, lchnk )
   call outfld('ACCIL_N' ,microp_st%acciln   , pcols, lchnk )

   call outfld('FALLR_M' ,microp_st%fallrm   , pcols, lchnk )
   call outfld('FALLS_M' ,microp_st%fallsm   , pcols, lchnk )
   call outfld('FALLG_M' ,microp_st%fallgm   , pcols, lchnk )
   call outfld('FALLR_N' ,microp_st%fallrn   , pcols, lchnk )
   call outfld('FALLS_N' ,microp_st%fallsn   , pcols, lchnk )
   call outfld('FALLG_N' ,microp_st%fallgn   , pcols, lchnk )

   call outfld('FHMR_M'  ,microp_st%fhmrm    , pcols, lchnk )

   do i = 1,ncol
      if (prec(i) .gt. 0) then
         precz_snum(i) = 1
      else
         precz_snum(i) = 0
      end if
   end do
   call outfld('PRECZ_SN', precz_snum, pcols, lchnk )

end subroutine zm_microphysics_history_out

!===================================================================================================

end module  zm_microphysics_history
