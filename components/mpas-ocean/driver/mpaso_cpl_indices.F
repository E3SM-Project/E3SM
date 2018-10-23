module mpaso_cpl_indices
  
  use seq_flds_mod
  use mct_mod

  implicit none

  SAVE
  public                               ! By default make data private

  ! ocn -> drv

  integer :: index_o2x_So_t      
  integer :: index_o2x_So_u
  integer :: index_o2x_So_v
  integer :: index_o2x_So_s
  integer :: index_o2x_So_dhdx
  integer :: index_o2x_So_dhdy
  integer :: index_o2x_Fioo_meltp
  integer :: index_o2x_Fioo_frazil
  integer :: index_o2x_Faoo_fco2_ocn
  integer :: index_o2x_Faoo_fdms_ocn

  ! ocn -> drv  (BGC)

  integer :: index_o2x_So_algae1
  integer :: index_o2x_So_algae2
  integer :: index_o2x_So_algae3
  integer :: index_o2x_So_doc1
  integer :: index_o2x_So_doc2
  integer :: index_o2x_So_doc3
  integer :: index_o2x_So_dic1
  integer :: index_o2x_So_don1
  integer :: index_o2x_So_no3
  integer :: index_o2x_So_sio3
  integer :: index_o2x_So_nh4
  integer :: index_o2x_So_dms
  integer :: index_o2x_So_dmsp
  integer :: index_o2x_So_docr
  integer :: index_o2x_So_fep1
  integer :: index_o2x_So_fep2
  integer :: index_o2x_So_fed1
  integer :: index_o2x_So_fed2
  integer :: index_o2x_So_zaer1
  integer :: index_o2x_So_zaer2
  integer :: index_o2x_So_zaer3
  integer :: index_o2x_So_zaer4
  integer :: index_o2x_So_zaer5
  integer :: index_o2x_So_zaer6

  ! drv -> ocn

  integer :: index_x2o_Si_ifrac        ! fractional ice wrt ocean
  integer :: index_x2o_Si_bpress       ! ice basal pressure
  integer :: index_x2o_So_duu10n       ! 10m wind speed squared           (m^2/s^2)
  integer :: index_x2o_Sa_pbot         ! atm bottom pressure              (Pa)
  integer :: index_x2o_Sa_co2prog      ! bottom atm level prognostic CO2
  integer :: index_x2o_Sa_co2diag      ! bottom atm level diagnostic CO2
  integer :: index_x2o_Foxx_taux       ! zonal wind stress (taux)         (W/m2   )
  integer :: index_x2o_Foxx_tauy       ! meridonal wind stress (tauy)     (W/m2   )
  integer :: index_x2o_Foxx_swnet      ! net short-wave heat flux         (W/m2   )
  integer :: index_x2o_Foxx_sen        ! sensible heat flux               (W/m2   )
  integer :: index_x2o_Foxx_lat        
  integer :: index_x2o_Foxx_lwup       ! longwave radiation (up)          (W/m2   )
  integer :: index_x2o_Faxa_lwdn       ! longwave radiation (down)        (W/m2   )
  integer :: index_x2o_Fioi_melth      ! heat flux from snow & ice melt   (W/m2   )
  integer :: index_x2o_Fioi_meltw      ! snow melt flux                   (kg/m2/s)
  integer :: index_x2o_Fioi_bergh      ! heat flux from iceberg melt      (W/m2   )
  integer :: index_x2o_Fioi_bergw      ! iceberg melt flux                (kg/m2/s)
  integer :: index_x2o_Fioi_salt       ! salt                             (kg(salt)/m2/s)
  integer :: index_x2o_Foxx_evap       ! evaporation flux                 (kg/m2/s)
  integer :: index_x2o_Faxa_prec         
  integer :: index_x2o_Faxa_snow       ! water flux due to snow           (kg/m2/s)
  integer :: index_x2o_Faxa_rain       ! water flux due to rain           (kg/m2/s)
  integer :: index_x2o_Faxa_bcphidry   ! flux: Black   Carbon hydrophilic dry deposition
  integer :: index_x2o_Faxa_bcphodry   ! flux: Black   Carbon hydrophobic dry deposition
  integer :: index_x2o_Faxa_bcphiwet   ! flux: Black   Carbon hydrophilic wet deposition
  integer :: index_x2o_Faxa_ocphidry   ! flux: Organic Carbon hydrophilic dry deposition
  integer :: index_x2o_Faxa_ocphodry   ! flux: Organic Carbon hydrophobic dry deposition
  integer :: index_x2o_Faxa_ocphiwet   ! flux: Organic Carbon hydrophilic dry deposition
  integer :: index_x2o_Faxa_dstwet1    ! flux: Size 1 dust -- wet deposition
  integer :: index_x2o_Faxa_dstwet2    ! flux: Size 2 dust -- wet deposition
  integer :: index_x2o_Faxa_dstwet3    ! flux: Size 3 dust -- wet deposition
  integer :: index_x2o_Faxa_dstwet4    ! flux: Size 4 dust -- wet deposition
  integer :: index_x2o_Faxa_dstdry1    ! flux: Size 1 dust -- dry deposition
  integer :: index_x2o_Faxa_dstdry2    ! flux: Size 2 dust -- dry deposition
  integer :: index_x2o_Faxa_dstdry3    ! flux: Size 3 dust -- dry deposition
  integer :: index_x2o_Faxa_dstdry4    ! flux: Size 4 dust -- dry deposition
  integer :: index_x2o_Foxx_rofl       ! river runoff flux                (kg/m2/s)
  integer :: index_x2o_Foxx_rofi       ! ice runoff flux                  (kg/m2/s)

  ! drv -> ocn (BGC)

  integer :: index_x2o_Fioi_algae1
  integer :: index_x2o_Fioi_algae2
  integer :: index_x2o_Fioi_algae3
  integer :: index_x2o_Fioi_doc1
  integer :: index_x2o_Fioi_doc2
  integer :: index_x2o_Fioi_doc3
  integer :: index_x2o_Fioi_dic1
  integer :: index_x2o_Fioi_don1
  integer :: index_x2o_Fioi_no3
  integer :: index_x2o_Fioi_sio3
  integer :: index_x2o_Fioi_nh4
  integer :: index_x2o_Fioi_dms
  integer :: index_x2o_Fioi_dmspp
  integer :: index_x2o_Fioi_dmspd
  integer :: index_x2o_Fioi_docr
  integer :: index_x2o_Fioi_fep1
  integer :: index_x2o_Fioi_fep2
  integer :: index_x2o_Fioi_fed1
  integer :: index_x2o_Fioi_fed2
  integer :: index_x2o_Fioi_dust1

contains

  subroutine mpaso_cpl_indices_set( )

    type(mct_aVect) :: o2x      ! temporary
    type(mct_aVect) :: x2o      ! temporary

    ! Determine attribute vector indices

    ! create temporary attribute vectors
    call mct_aVect_init(x2o, rList=seq_flds_x2o_fields, lsize=1)
    call mct_aVect_init(o2x, rList=seq_flds_o2x_fields, lsize=1)

    index_o2x_So_t          = mct_avect_indexra(o2x,'So_t')
    index_o2x_So_u          = mct_avect_indexra(o2x,'So_u')
    index_o2x_So_v          = mct_avect_indexra(o2x,'So_v')
    index_o2x_So_s          = mct_avect_indexra(o2x,'So_s')
    index_o2x_So_dhdx       = mct_avect_indexra(o2x,'So_dhdx')
    index_o2x_So_dhdy       = mct_avect_indexra(o2x,'So_dhdy')
    index_o2x_Fioo_meltp    = mct_avect_indexra(o2x,'Fioo_meltp',perrWith='quiet')
    index_o2x_Fioo_frazil   = mct_avect_indexra(o2x,'Fioo_frazil',perrWith='quiet')
    index_o2x_Faoo_fco2_ocn = mct_avect_indexra(o2x,'Faoo_fco2_ocn',perrWith='quiet')
    index_o2x_Faoo_fdms_ocn = mct_avect_indexra(o2x,'Faoo_fdms_ocn',perrWith='quiet')

    index_o2x_So_algae1     = mct_avect_indexra(o2x,'So_algae1',perrWith='quiet')
    index_o2x_So_algae2     = mct_avect_indexra(o2x,'So_algae2',perrWith='quiet')
    index_o2x_So_algae3     = mct_avect_indexra(o2x,'So_algae3',perrWith='quiet')
    index_o2x_So_doc1       = mct_avect_indexra(o2x,'So_doc1',perrWith='quiet')
    index_o2x_So_doc2       = mct_avect_indexra(o2x,'So_doc2',perrWith='quiet')
    index_o2x_So_doc3       = mct_avect_indexra(o2x,'So_doc3',perrWith='quiet')
    index_o2x_So_dic1       = mct_avect_indexra(o2x,'So_dic1',perrWith='quiet')
    index_o2x_So_don1       = mct_avect_indexra(o2x,'So_don1',perrWith='quiet')
    index_o2x_So_no3        = mct_avect_indexra(o2x,'So_no3',perrWith='quiet')
    index_o2x_So_sio3       = mct_avect_indexra(o2x,'So_sio3',perrWith='quiet')
    index_o2x_So_nh4        = mct_avect_indexra(o2x,'So_nh4',perrWith='quiet')
    index_o2x_So_dms        = mct_avect_indexra(o2x,'So_dms',perrWith='quiet')
    index_o2x_So_dmsp       = mct_avect_indexra(o2x,'So_dmsp',perrWith='quiet')
    index_o2x_So_docr       = mct_avect_indexra(o2x,'So_docr',perrWith='quiet')
    index_o2x_So_fep1       = mct_avect_indexra(o2x,'So_fep1',perrWith='quiet')
    index_o2x_So_fep2       = mct_avect_indexra(o2x,'So_fep2',perrWith='quiet')
    index_o2x_So_fed1       = mct_avect_indexra(o2x,'So_fed1',perrWith='quiet')
    index_o2x_So_fed2       = mct_avect_indexra(o2x,'So_fed2',perrWith='quiet')
    index_o2x_So_zaer1      = mct_avect_indexra(o2x,'So_zaer1',perrWith='quiet')
    index_o2x_So_zaer2      = mct_avect_indexra(o2x,'So_zaer2',perrWith='quiet')
    index_o2x_So_zaer3      = mct_avect_indexra(o2x,'So_zaer3',perrWith='quiet')
    index_o2x_So_zaer4      = mct_avect_indexra(o2x,'So_zaer4',perrWith='quiet')
    index_o2x_So_zaer5      = mct_avect_indexra(o2x,'So_zaer5',perrWith='quiet')
    index_o2x_So_zaer6      = mct_avect_indexra(o2x,'So_zaer6',perrWith='quiet')

    index_x2o_Si_ifrac      = mct_avect_indexra(x2o,'Si_ifrac')
    index_x2o_Si_bpress     = mct_avect_indexra(x2o,'Si_bpress')
    index_x2o_Sa_pbot       = mct_avect_indexra(x2o,'Sa_pbot')
    index_x2o_So_duu10n     = mct_avect_indexra(x2o,'So_duu10n')
    index_x2o_Foxx_tauy     = mct_avect_indexra(x2o,'Foxx_tauy')
    index_x2o_Foxx_taux     = mct_avect_indexra(x2o,'Foxx_taux')
    index_x2o_Foxx_swnet    = mct_avect_indexra(x2o,'Foxx_swnet')
    index_x2o_Foxx_lat      = mct_avect_indexra(x2o,'Foxx_lat')
    index_x2o_Foxx_sen      = mct_avect_indexra(x2o,'Foxx_sen')
    index_x2o_Foxx_lwup     = mct_avect_indexra(x2o,'Foxx_lwup')
    index_x2o_Faxa_lwdn     = mct_avect_indexra(x2o,'Faxa_lwdn')
    index_x2o_Fioi_melth    = mct_avect_indexra(x2o,'Fioi_melth')   
    index_x2o_Fioi_meltw    = mct_avect_indexra(x2o,'Fioi_meltw')
    index_x2o_Fioi_bergh    = mct_avect_indexra(x2o,'PFioi_bergh')
    index_x2o_Fioi_bergw    = mct_avect_indexra(x2o,'PFioi_bergw')
    index_x2o_Fioi_salt     = mct_avect_indexra(x2o,'Fioi_salt')   
    index_x2o_Faxa_prec     = mct_avect_indexra(x2o,'Faxa_prec')   
    index_x2o_Faxa_snow     = mct_avect_indexra(x2o,'Faxa_snow')   
    index_x2o_Faxa_rain     = mct_avect_indexra(x2o,'Faxa_rain')   
    index_x2o_Foxx_evap     = mct_avect_indexra(x2o,'Foxx_evap')
    index_x2o_Foxx_rofl     = mct_avect_indexra(x2o,'Foxx_rofl')
    index_x2o_Foxx_rofi     = mct_avect_indexra(x2o,'Foxx_rofi')
    index_x2o_Faxa_bcphidry = mct_avect_indexra(x2o,'Faxa_bcphidry')
    index_x2o_Faxa_bcphodry = mct_avect_indexra(x2o,'Faxa_bcphodry')
    index_x2o_Faxa_bcphiwet = mct_avect_indexra(x2o,'Faxa_bcphiwet')
    index_x2o_Faxa_ocphidry = mct_avect_indexra(x2o,'Faxa_ocphidry')
    index_x2o_Faxa_ocphodry = mct_avect_indexra(x2o,'Faxa_ocphodry')
    index_x2o_Faxa_ocphiwet = mct_avect_indexra(x2o,'Faxa_ocphiwet')
    index_x2o_Faxa_dstdry1  = mct_avect_indexra(x2o,'Faxa_dstdry1')
    index_x2o_Faxa_dstdry2  = mct_avect_indexra(x2o,'Faxa_dstdry2')
    index_x2o_Faxa_dstdry3  = mct_avect_indexra(x2o,'Faxa_dstdry3')
    index_x2o_Faxa_dstdry4  = mct_avect_indexra(x2o,'Faxa_dstdry4')
    index_x2o_Faxa_dstwet1  = mct_avect_indexra(x2o,'Faxa_dstwet1')
    index_x2o_Faxa_dstwet2  = mct_avect_indexra(x2o,'Faxa_dstwet2')
    index_x2o_Faxa_dstwet3  = mct_avect_indexra(x2o,'Faxa_dstwet3')
    index_x2o_Faxa_dstwet4  = mct_avect_indexra(x2o,'Faxa_dstwet4')
    index_x2o_Sa_co2prog    = mct_avect_indexra(x2o,'Sa_co2prog',perrWith='quiet')
    index_x2o_Sa_co2diag    = mct_avect_indexra(x2o,'Sa_co2diag',perrWith='quiet')

    index_x2o_Fioi_algae1   = mct_avect_indexra(x2o,'Fioi_algae1',perrWith='quiet')
    index_x2o_Fioi_algae2   = mct_avect_indexra(x2o,'Fioi_algae2',perrWith='quiet')
    index_x2o_Fioi_algae3   = mct_avect_indexra(x2o,'Fioi_algae3',perrWith='quiet')
    index_x2o_Fioi_doc1     = mct_avect_indexra(x2o,'Fioi_doc1',perrWith='quiet')
    index_x2o_Fioi_doc2     = mct_avect_indexra(x2o,'Fioi_doc2',perrWith='quiet')
    index_x2o_Fioi_doc3     = mct_avect_indexra(x2o,'Fioi_doc3',perrWith='quiet')
    index_x2o_Fioi_dic1     = mct_avect_indexra(x2o,'Fioi_dic1',perrWith='quiet')
    index_x2o_Fioi_don1     = mct_avect_indexra(x2o,'Fioi_don1',perrWith='quiet')
    index_x2o_Fioi_no3      = mct_avect_indexra(x2o,'Fioi_no3',perrWith='quiet')
    index_x2o_Fioi_sio3     = mct_avect_indexra(x2o,'Fioi_sio3',perrWith='quiet')
    index_x2o_Fioi_nh4      = mct_avect_indexra(x2o,'Fioi_nh4',perrWith='quiet')
    index_x2o_Fioi_dms      = mct_avect_indexra(x2o,'Fioi_dms',perrWith='quiet')
    index_x2o_Fioi_dmspp    = mct_avect_indexra(x2o,'Fioi_dmspp',perrWith='quiet')
    index_x2o_Fioi_dmspd    = mct_avect_indexra(x2o,'Fioi_dmspd',perrWith='quiet')
    index_x2o_Fioi_docr     = mct_avect_indexra(x2o,'Fioi_docr',perrWith='quiet')
    index_x2o_Fioi_fep1     = mct_avect_indexra(x2o,'Fioi_fep1',perrWith='quiet')
    index_x2o_Fioi_fep2     = mct_avect_indexra(x2o,'Fioi_fep2',perrWith='quiet')
    index_x2o_Fioi_fed1     = mct_avect_indexra(x2o,'Fioi_fed1',perrWith='quiet')
    index_x2o_Fioi_fed2     = mct_avect_indexra(x2o,'Fioi_fed2',perrWith='quiet')
    index_x2o_Fioi_dust1    = mct_avect_indexra(x2o,'Fioi_dust1',perrWith='quiet')

    call mct_aVect_clean(x2o)
    call mct_aVect_clean(o2x)

  end subroutine mpaso_cpl_indices_set

end module mpaso_cpl_indices
