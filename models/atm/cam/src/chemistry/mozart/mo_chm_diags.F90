module mo_chm_diags

  use shr_kind_mod, only : r8 => shr_kind_r8
  use chem_mods,    only : gas_pcnst
  use mo_tracname,  only : solsym
  use chem_mods,    only : rxntot, nfs, gas_pcnst, indexm, adv_mass
  use ppgrid,       only : pver
  use mo_constants, only : pi, rgrav, rearth
  use mo_chem_utls, only : get_rxt_ndx, get_spc_ndx
  use cam_history,  only : fieldname_len
  use mo_jeuv,      only : neuv

  private

  public :: chm_diags_inti
  public :: chm_diags
  public :: het_diags

  integer :: id_n,id_no,id_no2,id_no3,id_n2o5,id_hno3,id_ho2no2,id_clono2,id_brono2
  integer :: id_cl,id_clo,id_hocl,id_cl2,id_cl2o2,id_oclo,id_hcl,id_brcl
  integer :: id_ccl4,id_cfc11,id_cfc113,id_ch3ccl3,id_cfc12,id_ch3cl,id_hcfc22,id_cf2clbr
  integer :: id_br,id_bro,id_hbr,id_hobr,id_ch4,id_h2o,id_h2
  integer :: id_o,id_o2,id_h

  integer, parameter :: NJEUV = neuv
  integer :: rid_jeuv(NJEUV), rid_jno_i, rid_jno

  logical :: has_jeuvs, has_jno_i, has_jno

  integer :: nox_species(3),  noy_species(15)
  integer :: clox_species(6), cloy_species(9), tcly_species(17)
  integer :: brox_species(2), broy_species(6)
  integer :: toth_species(3)
  integer :: sox_species(3)
  integer :: nhx_species(3)
  integer :: aer_species(gas_pcnst)

  character(len=fieldname_len) :: dtchem_name(gas_pcnst)
  character(len=fieldname_len) :: depvel_name(gas_pcnst)
  character(len=fieldname_len) :: depflx_name(gas_pcnst)
  character(len=fieldname_len) :: wetdep_name(gas_pcnst)
  character(len=fieldname_len) :: wtrate_name(gas_pcnst)

  real(r8), parameter :: N_molwgt = 14.00674_r8
  real(r8), parameter :: S_molwgt = 32.066_r8

  character(len=32) :: chempkg

contains

  subroutine chm_diags_inti
    !--------------------------------------------------------------------
    !	... initialize utility routine
    !--------------------------------------------------------------------

    use cam_history,  only : addfld, horiz_only, add_default
    use constituents, only : cnst_get_ind, cnst_longname
    use dyn_grid,     only : get_dyn_grid_parm, get_horiz_grid_d
    use phys_control, only: phys_getopts

    implicit none

    integer :: i, j, k, m, n
    character(len=16) :: jname, spc_name, attr
    character(len=2)  :: jchar
    character(len=64) :: lname
    character(len=2)  :: unit_basename  ! Units 'kg' or '1' 

    integer :: id_pan, id_onit, id_mpan, id_isopno3, id_onitr, id_nh4no3
    integer :: id_so2, id_so4, id_h2so4
    integer :: id_nh3, id_nh4
    integer :: id_dst01, id_dst02, id_dst03, id_dst04, id_sslt01, id_sslt02, id_sslt03, id_sslt04
    integer :: id_soa,  id_oc1, id_oc2, id_cb1, id_cb2
    integer :: id_bry, id_cly 
    integer :: id_soam,id_soai,id_soat,id_soab,id_soax

    logical :: history_aerosol      ! Output the MAM aerosol tendencies
    logical :: history_amwg         ! output the variables used by the AMWG diag package
    integer :: bulkaero_species(20)

    !-----------------------------------------------------------------------

    call phys_getopts( history_aerosol_out = history_aerosol, &
                       history_amwg_out    = history_amwg,  &
                       cam_chempkg_out     = chempkg   )

    id_bry     = get_spc_ndx( 'BRY' )
    id_cly     = get_spc_ndx( 'CLY' )

    id_n       = get_spc_ndx( 'N' )
    id_no      = get_spc_ndx( 'NO' )
    id_no2     = get_spc_ndx( 'NO2' )
    id_no3     = get_spc_ndx( 'NO3' )
    id_n2o5    = get_spc_ndx( 'N2O5' )
    id_hno3    = get_spc_ndx( 'HNO3' )
    id_ho2no2  = get_spc_ndx( 'HO2NO2' )
    id_clono2  = get_spc_ndx( 'CLONO2' )
    id_brono2  = get_spc_ndx( 'BRONO2' )
    id_cl      = get_spc_ndx( 'CL' )
    id_clo     = get_spc_ndx( 'CLO' )
    id_hocl    = get_spc_ndx( 'HOCL' )
    id_cl2     = get_spc_ndx( 'CL2' )
    id_cl2o2   = get_spc_ndx( 'CL2O2' )
    id_oclo    = get_spc_ndx( 'OCLO' )
    id_hcl     = get_spc_ndx( 'HCL' )
    id_brcl    = get_spc_ndx( 'BRCL' )
    id_ccl4    = get_spc_ndx( 'CCL4' )
    id_cfc11   = get_spc_ndx( 'CFC11' )
    id_cfc113  = get_spc_ndx( 'CFC113' )
    id_ch3ccl3 = get_spc_ndx( 'CH3CCL3' )
    id_cfc12   = get_spc_ndx( 'CFC12' )
    id_ch3cl   = get_spc_ndx( 'CH3CL' )
    id_hcfc22  = get_spc_ndx( 'HCFC22' )
    id_cf2clbr = get_spc_ndx( 'CF2CLBR' )
    id_br      = get_spc_ndx( 'BR' )
    id_bro     = get_spc_ndx( 'BRO' )
    id_hbr     = get_spc_ndx( 'HBR' )
    id_hobr    = get_spc_ndx( 'HOBR' )
    id_ch4     = get_spc_ndx( 'CH4' )
    id_h2o     = get_spc_ndx( 'H2O' )
    id_h2      = get_spc_ndx( 'H2' )
    id_o       = get_spc_ndx( 'O' )
    id_o2      = get_spc_ndx( 'O2' )
    id_h       = get_spc_ndx( 'H' )

    id_pan     = get_spc_ndx( 'PAN' )
    id_onit    = get_spc_ndx( 'ONIT' )
    id_mpan    = get_spc_ndx( 'MPAN' )
    id_isopno3 = get_spc_ndx( 'ISOPNO3' )
    id_onitr   = get_spc_ndx( 'ONITR' )
    id_nh4no3  = get_spc_ndx( 'NH4NO3' )

    id_so2     = get_spc_ndx( 'SO2' )
    id_so4     = get_spc_ndx( 'SO4' )
    id_h2so4   = get_spc_ndx( 'H2SO4' )

    id_nh3     = get_spc_ndx( 'NH3' )
    id_nh4     = get_spc_ndx( 'NH4' )
    id_nh4no3  = get_spc_ndx( 'NH4NO3' )

    id_dst01   = get_spc_ndx( 'DST01' )
    id_dst02   = get_spc_ndx( 'DST02' )
    id_dst03   = get_spc_ndx( 'DST03' )
    id_dst04   = get_spc_ndx( 'DST04' )
    id_sslt01  = get_spc_ndx( 'SSLT01' )
    id_sslt02  = get_spc_ndx( 'SSLT02' )
    id_sslt03  = get_spc_ndx( 'SSLT03' )
    id_sslt04  = get_spc_ndx( 'SSLT04' )
    id_soa     = get_spc_ndx( 'SOA' )
    id_so4     = get_spc_ndx( 'SO4' )
    id_oc1     = get_spc_ndx( 'OC1' )
    id_oc2     = get_spc_ndx( 'OC2' )
    id_cb1     = get_spc_ndx( 'CB1' )
    id_cb2     = get_spc_ndx( 'CB2' )

    rid_jno   = get_rxt_ndx( 'jno' )
    rid_jno_i = get_rxt_ndx( 'jno_i' )

    id_soam = get_spc_ndx( 'SOAM' )
    id_soai = get_spc_ndx( 'SOAI' )
    id_soat = get_spc_ndx( 'SOAT' )
    id_soab = get_spc_ndx( 'SOAB' )
    id_soax = get_spc_ndx( 'SOAX' )

    nox_species = (/ id_n, id_no, id_no2 /)
    noy_species = (/ id_n, id_no, id_no2, id_no3, id_n2o5, id_hno3, id_ho2no2, id_clono2, &
                     id_brono2, id_pan, id_onit, id_mpan, id_isopno3, id_onitr, id_nh4no3 /)
    clox_species = (/ id_cl, id_clo, id_hocl, id_cl2, id_cl2o2, id_oclo /)
    cloy_species = (/ id_cl, id_clo, id_hocl, id_cl2, id_cl2o2, id_oclo, id_hcl, id_clono2, id_brcl /)
    tcly_species = (/ id_cl, id_clo, id_hocl, id_cl2, id_cl2o2, id_oclo, id_hcl, id_clono2, id_brcl, &
                      id_ccl4, id_cfc11, id_cfc113, id_ch3ccl3, id_cfc12, id_ch3cl, id_hcfc22, id_cf2clbr /)

    brox_species = (/ id_br, id_bro /)
    broy_species = (/ id_br, id_bro, id_hbr, id_brono2, id_brcl, id_hobr /)

    sox_species = (/ id_so2, id_so4, id_h2so4 /)
    nhx_species = (/ id_nh3, id_nh4, id_nh4no3 /)
    bulkaero_species(:) = -1
    bulkaero_species(1:20) = (/ id_dst01, id_dst02, id_dst03, id_dst04, &
                                id_sslt01, id_sslt02, id_sslt03, id_sslt04, &
                                id_soa, id_so4, id_oc1, id_oc2, id_cb1, id_cb2, id_nh4no3, &
                                id_soam,id_soai,id_soat,id_soab,id_soax /)

    aer_species(:) = -1
    n = 1
    do m = 1,gas_pcnst
       k=0
       if ( any(bulkaero_species(:)==m) ) k=1
       if ( k==0 ) k = index(trim(solsym(m)), '_a')
       if ( k==0 ) k = index(trim(solsym(m)), '_c')
       if ( k>0 ) then ! must be aerosol species
          aer_species(n) = m
          n = n+1
       endif
    enddo

    toth_species = (/ id_ch4, id_h2o, id_h2 /)

    call addfld( 'NOX', (/ 'lev' /), 'A', 'mol/mol', 'nox volume mixing ratio' )
    call addfld( 'NOY', (/ 'lev' /), 'A', 'mol/mol', 'noy volume mixing ratio' )
    call addfld( 'BROX', (/ 'lev' /), 'A','mol/mol', 'brox volume mixing ratio' )
    call addfld( 'BROY', (/ 'lev' /), 'A','mol/mol', 'total inorganic bromine (Br+BrO+HOBr+BrONO2+HBr+BrCl)' )
    call addfld( 'CLOX', (/ 'lev' /), 'A','mol/mol', 'clox volume mixing ratio' )
    call addfld( 'CLOY', (/ 'lev' /), 'A','mol/mol', 'total inorganic chlorine (Cl+ClO+2Cl2+2Cl2O2+OClO+HOCl+ClONO2+HCl+BrCl)' )
    call addfld( 'TCLY', (/ 'lev' /), 'A','mol/mol', 'total Cl volume mixing ratio' )
    call addfld( 'TOTH', (/ 'lev' /), 'A','mol/mol', 'total H2 volume mixing ratio' )

    call addfld( 'NOY_mmr', (/ 'lev' /), 'A', 'kg/kg', 'NOy mass mixing ratio' )
    call addfld( 'SOX_mmr', (/ 'lev' /), 'A', 'kg/kg', 'SOx mass mixing ratio' )
    call addfld( 'NHX_mmr', (/ 'lev' /), 'A', 'kg/kg', 'NHx mass mixing ratio' )

    do j = 1,NJEUV
       write( jchar, '(I2)' ) j
       jname = 'jeuv_'//trim(adjustl(jchar))
       rid_jeuv(j) = get_rxt_ndx( trim(jname) )
    enddo

    has_jeuvs = all( rid_jeuv(:) > 0 )
    has_jno_i = rid_jno_i>0
    has_jno   = rid_jno>0

    if ( has_jeuvs ) then
       call addfld( 'PION_EUV', (/ 'lev' /), 'I','/cm^3/s', 'total euv ionization rate' )
       call addfld( 'PEUV1', (/ 'lev' /), 'I',   '/cm^3/s', '(j1+j2+j3)*o' )
       call addfld( 'PEUV1e', (/ 'lev' /), 'I',  '/cm^3/s', '(j14+j15+j16)*o' )
       call addfld( 'PEUV2', (/ 'lev' /), 'I',   '/cm^3/s', 'j4*n' )
       call addfld( 'PEUV3', (/ 'lev' /), 'I',   '/cm^3/s', '(j5+j7+j8+j9)*o2' )
       call addfld( 'PEUV3e', (/ 'lev' /), 'I',  '/cm^3/s', '(j17+j19+j20+j21)*o2' )
       call addfld( 'PEUV4', (/ 'lev' /), 'I',   '/cm^3/s', '(j10+j11)*n2' )
       call addfld( 'PEUV4e', (/ 'lev' /), 'I',  '/cm^3/s', '(j22+j23)*n2' )
       call addfld( 'PEUVN2D', (/ 'lev' /), 'I', '/cm^3/s', '(j11+j13)*n2' )
       call addfld( 'PEUVN2De', (/ 'lev' /), 'I','/cm^3/s', '(j23+j25)*n2' )
    endif
    if ( has_jno ) then
       call addfld( 'PJNO', (/ 'lev' /), 'I', '/cm^3/s', 'jno*no' )
    endif
    if ( has_jno_i ) then
       call addfld( 'PJNO_I', (/ 'lev' /), 'I', '/cm^3/s', 'jno_i*no' )
    endif

    do m = 1,gas_pcnst

       spc_name = trim(solsym(m))

       call cnst_get_ind(spc_name, n, abort=.false. )
       if ( n > 0 ) then
          attr = cnst_longname(n)
       elseif ( trim(spc_name) == 'H2O' ) then
          attr = 'water vapor'
       else
          attr = spc_name
       endif

       depvel_name(m) = 'DV_'//trim(spc_name)
       depflx_name(m) = 'DF_'//trim(spc_name)
       dtchem_name(m) = 'D'//trim(spc_name)//'CHM'

       call addfld( depvel_name(m),   horiz_only,    'A', 'cm/s', 'deposition velocity ' )
       call addfld( depflx_name(m), horiz_only,    'A', 'kg/m2/s', 'dry deposition flux ' )
       call addfld( dtchem_name(m),   (/ 'lev' /), 'A', 'kg/s', 'net tendency from chem' )

       wetdep_name(m) = 'WD_'//trim(spc_name)
       wtrate_name(m) = 'WDR_'//trim(spc_name)

       call addfld( wetdep_name(m),   horiz_only,    'A', 'kg/s', spc_name//' wet deposition' )
       call addfld( wtrate_name(m),   (/ 'lev' /), 'A',   '/s', spc_name//' wet deposition rate' )
       
       if (spc_name(1:3) == 'num') then
          unit_basename = ' 1'
       else
          unit_basename = 'kg'
       endif

       if ( any( aer_species == m ) ) then
          call addfld( spc_name,   (/ 'lev' /), 'A', unit_basename//'/kg ', trim(attr)//' concentration')
          call addfld( trim(spc_name)//'_SRF', horiz_only, 'A', unit_basename//'/kg', trim(attr)//" in bottom layer")    
       else
          call addfld( spc_name, (/ 'lev' /), 'A', 'mol/mol', trim(attr)//' concentration')
          call addfld( trim(spc_name)//'_SRF', horiz_only, 'A', 'mol/mol', trim(attr)//" in bottom layer")
       endif

       if ((m /= id_cly) .and. (m /= id_bry)) then
          if (history_aerosol) then
             call add_default( spc_name, 1, ' ' )
             call add_default( trim(spc_name)//'_SRF', 1, ' ' )
          endif 
          if (history_amwg) then
             call add_default( trim(spc_name)//'_SRF', 1, ' ' )
          endif
       endif

    enddo

    call addfld( 'MASS', (/ 'lev' /), 'A', 'kg', 'mass of grid box' )
    call addfld( 'AREA', horiz_only,    'A', 'm2', 'area of grid box' )

    call addfld( 'WD_NOY', horiz_only, 'A', 'kg/s', 'NOy wet deposition' )
    call addfld( 'DF_NOY', horiz_only, 'I', 'kg/m2/s', 'NOy dry deposition flux ' )

    call addfld( 'WD_SOX', horiz_only, 'A', 'kg/s', 'SOx wet deposition' )
    call addfld( 'DF_SOX', horiz_only, 'I', 'kg/m2/s', 'SOx dry deposition flux ' )

    call addfld( 'WD_NHX', horiz_only, 'A', 'kg/s', 'NHx wet deposition' )
    call addfld( 'DF_NHX', horiz_only, 'I', 'kg/m2/s', 'NHx dry deposition flux ' )

  end subroutine chm_diags_inti

  subroutine chm_diags( lchnk, ncol, vmr, mmr, rxt_rates, invariants, depvel, depflx, mmr_tend, pdel, pbuf )
    !--------------------------------------------------------------------
    !	... utility routine to output chemistry diagnostic variables
    !--------------------------------------------------------------------
    
    use cam_history,  only : outfld
    use constituents, only : pcnst
    use constituents, only : cnst_get_ind
    use phys_grid,    only : get_area_all_p, pcols
    use physics_buffer, only : physics_buffer_desc
    
    implicit none

    !--------------------------------------------------------------------
    !	... dummy arguments
    !--------------------------------------------------------------------
    integer,  intent(in)  :: lchnk
    integer,  intent(in)  :: ncol
    real(r8), intent(in)  :: vmr(ncol,pver,gas_pcnst)
    real(r8), intent(in)  :: mmr(ncol,pver,gas_pcnst)
    real(r8), intent(in)  :: rxt_rates(ncol,pver,rxntot)
    real(r8), intent(in)  :: invariants(ncol,pver,max(1,nfs))
    real(r8), intent(in)  :: depvel(ncol, gas_pcnst)
    real(r8), intent(in)  :: depflx(ncol, gas_pcnst)
    real(r8), intent(in)  :: mmr_tend(ncol,pver,gas_pcnst)
    real(r8), intent(in)  :: pdel(ncol,pver)
    type(physics_buffer_desc), pointer :: pbuf(:)

    !--------------------------------------------------------------------
    !	... local variables
    !--------------------------------------------------------------------
    integer     :: i,j,k, m, n
    integer :: plat
    real(r8)    :: wrk(ncol,pver)
    !      real(r8)    :: tmp(ncol,pver)
    !      real(r8)    :: m(ncol,pver)
    real(r8)    :: un2(ncol)
    
    real(r8), dimension(ncol,pver) :: vmr_nox, vmr_noy, vmr_clox, vmr_cloy, vmr_tcly, vmr_brox, vmr_broy, vmr_toth
    real(r8), dimension(ncol,pver) :: mmr_noy, mmr_sox, mmr_nhx, net_chem
    real(r8), dimension(ncol)      :: df_noy, df_sox, df_nhx

    real(r8) :: area(ncol), mass(ncol,pver)
    real(r8) :: wgt
    character(len=16) :: spc_name

    !--------------------------------------------------------------------
    !	... "diagnostic" groups
    !--------------------------------------------------------------------
    vmr_nox(:ncol,:) = 0._r8
    vmr_noy(:ncol,:) = 0._r8
    vmr_clox(:ncol,:) = 0._r8
    vmr_cloy(:ncol,:) = 0._r8
    vmr_tcly(:ncol,:) = 0._r8
    vmr_brox(:ncol,:) = 0._r8
    vmr_broy(:ncol,:) = 0._r8
    vmr_toth(:ncol,:) = 0._r8
    mmr_noy(:ncol,:) = 0._r8
    mmr_sox(:ncol,:) = 0._r8
    mmr_nhx(:ncol,:) = 0._r8
    df_noy(:ncol) = 0._r8
    df_sox(:ncol) = 0._r8
    df_nhx(:ncol) = 0._r8


    call get_area_all_p(lchnk, ncol, area)
    area = area * rearth**2

    do k = 1,pver
       mass(:ncol,k) = pdel(:ncol,k) * area(:ncol) * rgrav
    enddo

    call outfld( 'AREA', area(:ncol),   ncol, lchnk )
    call outfld( 'MASS', mass(:ncol,:), ncol, lchnk )

    do m = 1,gas_pcnst

       if ( m == id_ch4 .or. m == id_n2o5 .or. m == id_cfc12 .or. m == id_cl2 .or. m == id_cl2o2) then
          wgt = 2._r8
       elseif ( m == id_cfc11 .or. m == id_cfc113 .or. m == id_ch3ccl3 ) then
          wgt = 3._r8
       elseif ( m == id_ccl4 ) then
          wgt = 4._r8
       else
          wgt = 1._r8
       endif

       if ( any( nox_species == m ) ) then
          vmr_nox(:ncol,:) = vmr_nox(:ncol,:) +  wgt * vmr(:ncol,:,m)
       endif
       if ( any( noy_species == m ) ) then
          vmr_noy(:ncol,:) = vmr_noy(:ncol,:) +  wgt * vmr(:ncol,:,m)
       endif

       if ( any( noy_species == m ) ) then
          mmr_noy(:ncol,:) = mmr_noy(:ncol,:) +  wgt * mmr(:ncol,:,m)
       endif
       if ( any( sox_species == m ) ) then
          mmr_sox(:ncol,:) = mmr_sox(:ncol,:) +  wgt * mmr(:ncol,:,m)
       endif
       if ( any( nhx_species == m ) ) then
          mmr_nhx(:ncol,:) = mmr_nhx(:ncol,:) +  wgt * mmr(:ncol,:,m)
       endif

       if ( any( clox_species == m ) ) then
          vmr_clox(:ncol,:) = vmr_clox(:ncol,:) +  wgt * vmr(:ncol,:,m)
       endif
       if ( any( cloy_species == m ) ) then
          vmr_cloy(:ncol,:) = vmr_cloy(:ncol,:) +  wgt * vmr(:ncol,:,m)
       endif
       if ( any( tcly_species == m ) ) then
          vmr_tcly(:ncol,:) = vmr_tcly(:ncol,:) +  wgt * vmr(:ncol,:,m)
       endif

       if ( any( brox_species == m ) ) then
          vmr_brox(:ncol,:) = vmr_brox(:ncol,:) +  wgt * vmr(:ncol,:,m)
       endif
       if ( any( broy_species == m ) ) then
          vmr_broy(:ncol,:) = vmr_broy(:ncol,:) +  wgt * vmr(:ncol,:,m)
       endif

       if ( any ( toth_species == m ) ) then
          vmr_toth(:ncol,:) = vmr_toth(:ncol,:) +  wgt * vmr(:ncol,:,m)
       endif
       
       if ( any( aer_species == m ) ) then
          call outfld( solsym(m), mmr(:ncol,:,m), ncol ,lchnk )
          call outfld( trim(solsym(m))//'_SRF', mmr(:ncol,pver,m), ncol ,lchnk )
       else
          call outfld( solsym(m), vmr(:ncol,:,m), ncol ,lchnk )
          call outfld( trim(solsym(m))//'_SRF', vmr(:ncol,pver,m), ncol ,lchnk )
       endif

       call outfld( depvel_name(m), depvel(:ncol,m), ncol ,lchnk )
       call outfld( depflx_name(m), depflx(:ncol,m), ncol ,lchnk )

       if ( any( noy_species == m ) ) then
          df_noy(:ncol) = df_noy(:ncol) +  wgt * depflx(:ncol,m)*N_molwgt/adv_mass(m)
       endif
       if ( any( sox_species == m ) ) then
          df_sox(:ncol) = df_sox(:ncol) +  wgt * depflx(:ncol,m)*S_molwgt/adv_mass(m)
       endif
       if ( any( nhx_species == m ) ) then
          df_nhx(:ncol) = df_nhx(:ncol) +  wgt * depflx(:ncol,m)*N_molwgt/adv_mass(m)
       endif

       do k=1,pver
          do i=1,ncol
             net_chem(i,k) = mmr_tend(i,k,m) * mass(i,k) 
          end do
       end do
       call outfld( dtchem_name(m), net_chem(:ncol,:), ncol, lchnk )

    enddo


    call outfld( 'NOX',  vmr_nox(:ncol,:),  ncol, lchnk )
    call outfld( 'NOY',  vmr_noy(:ncol,:),  ncol, lchnk )
    call outfld( 'CLOX', vmr_clox(:ncol,:), ncol, lchnk )
    call outfld( 'CLOY', vmr_cloy(:ncol,:), ncol, lchnk )
    call outfld( 'BROX', vmr_brox(:ncol,:), ncol, lchnk )
    call outfld( 'BROY', vmr_broy(:ncol,:), ncol, lchnk )
    call outfld( 'TCLY', vmr_tcly(:ncol,:), ncol, lchnk )
    call outfld( 'NOY_mmr', mmr_noy(:ncol,:), ncol ,lchnk )
    call outfld( 'SOX_mmr', mmr_sox(:ncol,:), ncol ,lchnk )
    call outfld( 'NHX_mmr', mmr_nhx(:ncol,:), ncol ,lchnk )
    call outfld( 'DF_NOY', df_noy(:ncol), ncol ,lchnk )
    call outfld( 'DF_SOX', df_sox(:ncol), ncol ,lchnk )
    call outfld( 'DF_NHX', df_nhx(:ncol), ncol ,lchnk )

    !--------------------------------------------------------------------
    !	... euv ion production
    !--------------------------------------------------------------------

    jeuvs: if ( has_jeuvs ) then
       do k = 1,pver
          un2(:)   = 1._r8 - (vmr(:,k,id_o) + vmr(:,k,id_o2) + vmr(:,k,id_h))
          wrk(:,k) = vmr(:,k,id_o)*(rxt_rates(:,k,rid_jeuv(1)) + rxt_rates(:,k,rid_jeuv(2)) &
               + rxt_rates(:,k,rid_jeuv(3)) + rxt_rates(:,k,rid_jeuv(14)) &
               + rxt_rates(:,k,rid_jeuv(15)) + rxt_rates(:,k,rid_jeuv(16))) &
               + vmr(:,k,id_n)*rxt_rates(:,k,rid_jeuv(4)) &
               + vmr(:,k,id_o2)*(rxt_rates(:,k,rid_jeuv(5)) + rxt_rates(:,k,rid_jeuv(7)) &
               + rxt_rates(:,k,rid_jeuv(8)) + rxt_rates(:,k,rid_jeuv(9)) &
               + rxt_rates(:,k,rid_jeuv(17)) + rxt_rates(:,k,rid_jeuv(19)) &
               + rxt_rates(:,k,rid_jeuv(20)) + rxt_rates(:,k,rid_jeuv(21))) &
               + un2(:)*(rxt_rates(:,k,rid_jeuv(6)) + rxt_rates(:,k,rid_jeuv(10)) &
               + rxt_rates(:,k,rid_jeuv(11)) + rxt_rates(:,k,rid_jeuv(18)) &
               + rxt_rates(:,k,rid_jeuv(22)) + rxt_rates(:,k,rid_jeuv(23)))
          wrk(:,k) = wrk(:,k) * invariants(:,k,indexm)
       end do
       call outfld( 'PION_EUV', wrk, ncol, lchnk )

       do k = 1,pver
          wrk(:,k) = vmr(:,k,id_o)*(rxt_rates(:,k,rid_jeuv(1)) + rxt_rates(:,k,rid_jeuv(2)) &
               + rxt_rates(:,k,rid_jeuv(3)))
          wrk(:,k) = wrk(:,k) * invariants(:,k,indexm)
       end do
       call outfld( 'PEUV1', wrk, ncol, lchnk )
       do k = 1,pver
          wrk(:,k) = vmr(:,k,id_o)*(rxt_rates(:,k,rid_jeuv(14)) + rxt_rates(:,k,rid_jeuv(15)) &
               + rxt_rates(:,k,rid_jeuv(16)))
          wrk(:,k) = wrk(:,k) * invariants(:,k,indexm)
       end do
       call outfld( 'PEUV1e', wrk, ncol, lchnk )
       do k = 1,pver
          wrk(:,k) = vmr(:,k,id_n)*rxt_rates(:,k,rid_jeuv(4))
          wrk(:,k) = wrk(:,k) * invariants(:,k,indexm)
       end do
       call outfld( 'PEUV2', wrk, ncol, lchnk )
       do k = 1,pver
          wrk(:,k) = vmr(:,k,id_o2)*(rxt_rates(:,k,rid_jeuv(5)) + rxt_rates(:,k,rid_jeuv(7)) &
               + rxt_rates(:,k,rid_jeuv(8)) + rxt_rates(:,k,rid_jeuv(9)))
          wrk(:,k) = wrk(:,k) * invariants(:,k,indexm)
       end do
       call outfld( 'PEUV3', wrk, ncol, lchnk )
       do k = 1,pver
          wrk(:,k) = vmr(:,k,id_o2)*(rxt_rates(:,k,rid_jeuv(17)) + rxt_rates(:,k,rid_jeuv(19)) &
               + rxt_rates(:,k,rid_jeuv(20)) + rxt_rates(:,k,rid_jeuv(21)))
          wrk(:,k) = wrk(:,k) * invariants(:,k,indexm)
       end do
       call outfld( 'PEUV3e', wrk, ncol, lchnk )
       do k = 1,pver
          un2(:)   = 1._r8 - (vmr(:,k,id_o) + vmr(:,k,id_o2) + vmr(:,k,id_h))
          wrk(:,k) = un2(:)*(rxt_rates(:,k,rid_jeuv(6)) + rxt_rates(:,k,rid_jeuv(10)) + rxt_rates(:,k,rid_jeuv(11)))
          wrk(:,k) = wrk(:,k) * invariants(:,k,indexm)
       end do
       call outfld( 'PEUV4', wrk, ncol, lchnk )
       do k = 1,pver
          un2(:)   = 1._r8 - (vmr(:,k,id_o) + vmr(:,k,id_o2) + vmr(:,k,id_h))
          wrk(:,k) = un2(:)*(rxt_rates(:,k,rid_jeuv(18)) + rxt_rates(:,k,rid_jeuv(22)) + rxt_rates(:,k,rid_jeuv(23)))
          wrk(:,k) = wrk(:,k) * invariants(:,k,indexm)
       end do
       call outfld( 'PEUV4e', wrk, ncol, lchnk )
       do k = 1,pver
          un2(:)   = 1._r8 - (vmr(:,k,id_o) + vmr(:,k,id_o2) + vmr(:,k,id_h))
          wrk(:,k) = un2(:)*(rxt_rates(:,k,rid_jeuv(11)) + rxt_rates(:,k,rid_jeuv(13)))
          wrk(:,k) = wrk(:,k) * invariants(:,k,indexm)
       end do
       call outfld( 'PEUVN2D', wrk, ncol, lchnk )
       do k = 1,pver
          un2(:)   = 1._r8 - (vmr(:,k,id_o) + vmr(:,k,id_o2) + vmr(:,k,id_h))
          wrk(:,k) = un2(:)*(rxt_rates(:,k,rid_jeuv(23)) + rxt_rates(:,k,rid_jeuv(25)))
          wrk(:,k) = wrk(:,k) * invariants(:,k,indexm)
       end do
       call outfld( 'PEUVN2De', wrk, ncol, lchnk )
    endif jeuvs

    if ( has_jno_i ) then
       do k = 1,pver
          wrk(:,k) = vmr(:,k,id_no)*rxt_rates(:,k,rid_jno_i)
          wrk(:,k) = wrk(:,k) * invariants(:,k,indexm)
       end do
       call outfld( 'PJNO_I', wrk, ncol, lchnk )
    endif
    if ( has_jno ) then
       do k = 1,pver
          wrk(:,k) = vmr(:,k,id_no)*rxt_rates(:,k,rid_jno)
          wrk(:,k) = wrk(:,k) * invariants(:,k,indexm)
       end do
       call outfld( 'PJNO', wrk, ncol, lchnk )
    endif

  end subroutine chm_diags

  subroutine het_diags( het_rates, mmr, pdel, lchnk, ncol )

    use cam_history,  only : outfld
    use phys_grid,    only : get_wght_all_p
    implicit none

    integer,  intent(in)  :: lchnk
    integer,  intent(in)  :: ncol
    real(r8), intent(in)  :: het_rates(ncol,pver,max(1,gas_pcnst))
    real(r8), intent(in)  :: mmr(ncol,pver,gas_pcnst)
    real(r8), intent(in)  :: pdel(ncol,pver)

    real(r8), dimension(ncol) :: noy_wk, sox_wk, nhx_wk, wrk_wd
    integer :: m, k, j
    integer :: plat
    real(r8) :: wght(ncol)
    !
    ! output integrated wet deposition field
    !
    noy_wk(:) = 0._r8
    sox_wk(:) = 0._r8
    nhx_wk(:) = 0._r8

    call get_wght_all_p(lchnk, ncol, wght)

    do m = 1,gas_pcnst
       !
       ! compute vertical integral
       !
       wrk_wd(:ncol) = 0._r8
       do k = 1,pver
          wrk_wd(:ncol) = wrk_wd(:ncol) + het_rates(:ncol,k,m) * mmr(:ncol,k,m) * pdel(:ncol,k) 
       end do
       !
       wrk_wd(:ncol) = wrk_wd(:ncol) * rgrav * wght(:ncol) * rearth**2
       !

       call outfld( wetdep_name(m), wrk_wd(:ncol),               ncol, lchnk )
       call outfld( wtrate_name(m), het_rates(:ncol,:,m), ncol, lchnk )

       if ( any(noy_species == m ) ) then
          noy_wk(:ncol) = noy_wk(:ncol) + wrk_wd(:ncol)*N_molwgt/adv_mass(m)
       endif
       if ( m == id_n2o5 ) then  ! 2 NOy molecules in N2O5
          noy_wk(:ncol) = noy_wk(:ncol) + wrk_wd(:ncol)*N_molwgt/adv_mass(m)
       endif
       if ( any(sox_species == m ) ) then
          sox_wk(:ncol) = sox_wk(:ncol) + wrk_wd(:ncol)*S_molwgt/adv_mass(m)
       endif
       if ( any(nhx_species == m ) ) then
          nhx_wk(:ncol) = nhx_wk(:ncol) + wrk_wd(:ncol)*N_molwgt/adv_mass(m)
       endif

    end do
    
    call outfld( 'WD_NOY', noy_wk(:ncol), ncol, lchnk )
    call outfld( 'WD_SOX', sox_wk(:ncol), ncol, lchnk )
    call outfld( 'WD_NHX', nhx_wk(:ncol), ncol, lchnk )

  end subroutine het_diags

end module mo_chm_diags
