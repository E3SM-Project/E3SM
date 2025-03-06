module  zm_microphysics_history
  !-----------------------------------------------------------------------------
   ! Purpose: microphysics state structure definition and methods for ZM
   ! Original Author: Xialiang Song and Guang Zhang, June 2010
   !-----------------------------------------------------------------------------
   use shr_kind_mod,          only: r8=>shr_kind_r8
   use ppgrid,                only: pcols, pver, pverp
   use zm_microphysics_state, only: zm_microp_st
   use cam_history,           only: outfld, addfld, horiz_only

   public :: zm_microphysics_history_out ! write history output related to ZM microphysics
  
!===================================================================================================
contains
!===================================================================================================

! subroutine zm_microphysics_history_init()
!    !----------------------------------------------------------------------------
! end subroutine zm_microphysics_history_init

!===================================================================================================

subroutine zm_microphysics_history_out(microp_st, dlf, dif, dnlf, dnif, frz, lchnk, ncol)
   !----------------------------------------------------------------------------
   ! Arguments
   type(zm_microp_st),intent(in) :: microp_st   ! ZM microphysics data structure
   real(r8),          intent(in) :: dlf(:,:)    ! detrainment of conv cld liq water mixing ratio
   real(r8),          intent(in) :: dif(:,:)    ! detrainment of conv cld ice mixing ratio
   real(r8),          intent(in) :: dnlf(:,:)   ! detrainment of conv cld liq water num concen
   real(r8),          intent(in) :: dnif(:,:)   ! detrainment of conv cld ice num concen
   real(r8),          intent(in) :: frz(:,:)    ! heating rate due to freezing 
   integer,           intent(in) :: lchnk       ! chunk identifier
   integer,           intent(in) :: ncol        ! number of columns in chunk
   !----------------------------------------------------------------------------
   ! Local variables
   integer  :: i,k
   real(r8) :: cice_snum(pcols,pver)            ! convective cloud ice sample number
   real(r8) :: cliq_snum(pcols,pver)            ! convective cloud liquid sample number
   real(r8) :: crain_snum(pcols,pver)           ! convective rain water sample number
   real(r8) :: csnow_snum(pcols,pver)           ! convective snow sample number
   real(r8) :: cgraupel_snum(pcols,pver)        ! convective graupel sample number
   real(r8) :: wu_snum(pcols,pver)              ! vertical velocity sample number
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

end subroutine zm_microphysics_history_out

!===================================================================================================

end module  zm_microphysics_history