module output_aerocom_aie 
!------------------------------------------------------------------------------------------
! 
!  Purpose:  to output variables requried by AEROCOM AIE intercomparison. 
! 
!  Method:  to sample cloud properties at the cloud top to facilliate the comparions 
!           between the model and satellite observations. 
!           Also to sample cloud properties only at the time of satellite overpass. 
!
!  Author: Minghuai Wang (2008-04)
! 
!  Updated by Minghuai Wang on August, 2013 for the third AeroCOM AIE intercomparison (IND3)
!
!------------------------------------------------------------------------------------------  

   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid,       only: pcols, pver, begchunk, endchunk
   use cam_history,      only: addfld, add_default, horiz_only,  outfld
   use cam_history_support,   only : fillvalue
   use cam_abortutils,    only: endrun

   implicit none

   private 

   public ::  output_aerocom_aie_init 
   public ::  output_aerocom_aie_register
   public ::  aerocom_calc
   public :: cloud_top_aerocom
   public ::  output_aerocom

   logical, public :: do_aerocom_ind3 

   real(r8) :: p0 = 1.0e5    

!  variable index for pbuf
   integer :: angstrm_idx   ! Angstrom coefficients
   integer :: cdr_idx      ! droplet effective radius at top of liquid water clouds (meter)
   integer :: cdnc_idx     ! droplet number concentration in top layer of liquid water clouds (#/m3)
   integer :: cdnum_idx    ! column-integrated droplet number concentrations
   integer :: icnum_idx    ! column-integrated ice crystal number concentrations
   integer :: clt_idx      ! fractional cover by all clouds
   integer :: lcc_idx      ! fractional cover by liquid water clouds
   integer :: lwp_idx      ! in-cloud liquid water path for liquid clouds (kg/m^2)
   integer :: iwp_idx      !  in-cloud ice water path for ice clouds (kg/m^2)
   integer :: icr_idx      ! effective radius of crystals at top of ice clouds (meter)
   integer :: icc_idx      ! fractional cover by ice clouds
   integer :: cod_idx      ! in-cloud optical depth 
   integer :: ccn_idx      ! cloud condensation nuclei number concentration for liquid water
                                       ! clouds where activation corresponding to CDR and CDN
   integer :: ttop_idx     ! Temperature at top of clouds (K)
   integer :: ptop_idx     ! Pressure at top fo clouds (Pa)
   integer :: autoconv_idx      ! Column-integrated autoconversion rate
   integer :: accretn_idx      ! Column-integrated autoconversion rate
   integer :: rh700_idx      ! relative humidity at 700 hPa
   integer :: icnc_idx     ! ice crystal number concentration in top layer of ice clouds (#/m3)

   integer :: intccn_idx   
   integer :: colrv_idx
   integer :: rwp_idx
   integer :: lwp2_idx
   integer :: iwp2_idx

   integer :: ccn3d_idx    ! 
   integer :: autocl_idx   ! 3D autoconversion rate
   integer :: accretl_idx  ! 3D accreation rate

   integer :: cldo_idx     ! pbuf index cloud fraction 
   integer :: qrain_idx    ! pbuf index for grid-mean rain water
   integer :: cld_tau_idx  ! cloud optical depth
   integer :: rel_idx      ! droplet effective radius
   integer :: rei_idx      ! ice crystal effective radius
   
   integer :: cldliqbf_idx  ! liquid water before microphyiscs
   integer :: cldicebf_idx  ! ice water before microphyiscs
   integer :: numliqbf_idx  ! liquid droplet number before microphyiscs
   integer :: numicebf_idx  ! ice crystal number before microphyiscs

   integer :: ixcldice, ixcldliq, ixnumliq, ixnumice

CONTAINS

   subroutine  output_aerocom_aie_register ()
   use physics_buffer,      only : pbuf_add_field, dtype_r8

   call pbuf_add_field('angstrm',  'physpkg', dtype_r8, (/pcols/),      angstrm_idx)
   call pbuf_add_field('cdr',  'physpkg', dtype_r8, (/pcols/),      cdr_idx)
   call pbuf_add_field('cdnc',  'physpkg', dtype_r8, (/pcols/),      cdnc_idx)
   call pbuf_add_field('cdnum',  'physpkg', dtype_r8, (/pcols/),      cdnum_idx)
   call pbuf_add_field('icnum',  'physpkg', dtype_r8, (/pcols/),      icnum_idx)
   call pbuf_add_field('clt',  'physpkg', dtype_r8, (/pcols/),      clt_idx)
   call pbuf_add_field('lcc',  'physpkg', dtype_r8, (/pcols/),      lcc_idx)
   call pbuf_add_field('lwp',  'physpkg', dtype_r8, (/pcols/),      lwp_idx)
   call pbuf_add_field('iwp',  'physpkg', dtype_r8, (/pcols/),      iwp_idx)
   call pbuf_add_field('icr',  'physpkg', dtype_r8, (/pcols/),      icr_idx)
   call pbuf_add_field('icc',  'physpkg', dtype_r8, (/pcols/),      icc_idx)
   call pbuf_add_field('cod',  'physpkg', dtype_r8, (/pcols/),      cod_idx)
   call pbuf_add_field('ccn',  'physpkg', dtype_r8, (/pcols/),      ccn_idx)
   call pbuf_add_field('ttop', 'physpkg', dtype_r8, (/pcols/),      ttop_idx)
   call pbuf_add_field('ptop', 'physpkg', dtype_r8, (/pcols/),      ptop_idx)

   call pbuf_add_field('autoconv', 'physpkg', dtype_r8, (/pcols/),      autoconv_idx)
   call pbuf_add_field('accretn', 'physpkg', dtype_r8, (/pcols/),      accretn_idx)
   call pbuf_add_field('icnc', 'physpkg', dtype_r8, (/pcols/),      icnc_idx)
   call pbuf_add_field('rh700', 'physpkg', dtype_r8, (/pcols/),      rh700_idx)

   call pbuf_add_field('intccn', 'physpkg', dtype_r8, (/pcols/),      intccn_idx)
   call pbuf_add_field('colrv', 'physpkg', dtype_r8, (/pcols/),      colrv_idx)
   call pbuf_add_field('rwp', 'physpkg', dtype_r8, (/pcols/),      rwp_idx)
   call pbuf_add_field('lwp2', 'physpkg', dtype_r8, (/pcols/),      lwp2_idx)
   call pbuf_add_field('iwp2', 'physpkg', dtype_r8, (/pcols/),      iwp2_idx)

   call pbuf_add_field('ccn3d', 'physpkg', dtype_r8, (/pcols, pver/),     ccn3d_idx)
   call pbuf_add_field('autocl', 'physpkg', dtype_r8, (/pcols, pver/),     autocl_idx)
   call pbuf_add_field('accretl', 'physpkg', dtype_r8, (/pcols, pver/),     accretl_idx)

   call pbuf_add_field('QRAIN', 'physpkg', dtype_r8, (/pcols, pver/),     qrain_idx)
   call pbuf_add_field('cld_tau', 'physpkg', dtype_r8, (/pcols, pver/),     cld_tau_idx)

   call pbuf_add_field('cldliqbf', 'physpkg', dtype_r8, (/pcols, pver/),     cldliqbf_idx)
   call pbuf_add_field('cldicebf', 'physpkg', dtype_r8, (/pcols, pver/),     cldicebf_idx)
   call pbuf_add_field('numliqbf', 'physpkg', dtype_r8, (/pcols, pver/),     numliqbf_idx)
   call pbuf_add_field('numicebf', 'physpkg', dtype_r8, (/pcols, pver/),     numicebf_idx)

   return
   end subroutine output_aerocom_aie_register
!----------------------------------------------------------------------------------------------------

!====================================================================================================
   subroutine  output_aerocom_aie_init ()
   use infnan
   use physics_buffer,      only : pbuf_get_index, dtype_r8
   use constituents,  only: cnst_get_ind

   call addfld ('angstrm', horiz_only  ,'A' ,  '#', 'Angstrom coefficient', flag_xyfill=.true.)
   call addfld ('aerindex', horiz_only ,'A' ,  '#', 'Aerosol Index (Angstrom coefficient * AOD)', flag_xyfill=.true.)
   call addfld ('cdr', horiz_only  ,'A'   ,   'meter', &
               'Grid-cell mean droplet effective radius at top of liquid water clouds',  flag_xyfill=.true.)
   call addfld ('cdnc', horiz_only ,'A'    ,  '#/m3', &
               'Grid-cell mean droplet number concentration at top of liquid water clouds',  flag_xyfill=.true.)
   call addfld ('cdnum', horiz_only ,'A'    ,  '#/m2', &
               'Grid-cell mean column-integrated droplet number concentrations',  flag_xyfill=.true.)
   call addfld ('icnum', horiz_only ,'A'    ,  '#/m2',  &
               'Grid-cell mean column-integrated ice crystal number concentrations',  flag_xyfill=.true.)
   call addfld ('clt', horiz_only ,'A'    ,   'fraction','Fractional cover by all clouds',   flag_xyfill=.true.) 
   call addfld ('lcc', horiz_only ,'A'    ,   'fraction','Fractional cover by liquid water clouds',   flag_xyfill=.true.)
   call addfld ('lwp', horiz_only ,'A'    ,   'kg/m2',  &
               'Grid-cell mean liquid water path for liquid water clouds',  flag_xyfill=.true.)
   call addfld ('iwp', horiz_only ,'A'    ,   'kg/m2',  'Grid-cell mean ice water path for ice clouds',  flag_xyfill=.true.)
   call addfld ('icr', horiz_only ,'A'    ,   'meter',  &
               'Grid-cell mean effective radius of crystals at top of ice clouds',  flag_xyfill=.true.)
   call addfld ('icc', horiz_only ,'A'    ,   'fracton','Fractional cover by ice clouds',  flag_xyfill=.true.)
   call addfld ('cod', horiz_only ,'A'    ,   'amount',  'Grid-cell mean cloud optical depth',  flag_xyfill=.true.)
   call addfld ('ccn', horiz_only ,'A'    ,   '#/m3',    &
               'CCN number concentration at 0.3% at the top layer of liquid water clouds',  flag_xyfill=.true.)
   call addfld ('ttop', horiz_only ,'A'    ,  'K',      'Temperature at top of clouds',  flag_xyfill=.true.)
   call addfld ('ptop', horiz_only ,'A'    ,  'Pa',      'Pressure at top of clouds',  flag_xyfill=.true.)
   call addfld ('autoconv', horiz_only ,'A'    ,  'kg/m2/s',   'Grid-mean surface precipitation rate',  flag_xyfill=.true.)
   call addfld ('accretn', horiz_only ,'A'    ,  'kg/m2/s',    'Grid-mean surface precipitation rate',  flag_xyfill=.true.)
   call addfld ('icnc', horiz_only ,'A'    ,  '#/m3',   &
               'Ice crystal number concentration at top of ice clouds',  flag_xyfill=.true.)
   call addfld ('rh700', horiz_only ,'A'    , 'fraction', 'Relative humidity at 700 hPa',  flag_xyfill=.true.)

   call addfld ('rwp', horiz_only ,'A'    , 'kg/m2',  'Rain water path',  flag_xyfill=.true.)
   call addfld ('intccn', horiz_only ,'A'    , '#/m2', 'Column-integrated CCN number concentration',  flag_xyfill=.true.)
   call addfld ('colrv', horiz_only ,'A'    , 'm',  'Column-integrated volume-mean droplet effective radius',  flag_xyfill=.true.)
   call addfld ('lwp2', horiz_only ,'A'    ,   'kg/m2',  &
               'Grid-cell mean liquid water path for liquid water clouds (new)',  flag_xyfill=.true.)
   call addfld ('iwp2', horiz_only ,'A'    ,   'kg/m2',  'Grid-cell mean ice water path for ice clouds (new)',  flag_xyfill=.true.)

   call addfld ('lwpbf', horiz_only ,'A'    ,   'kg/m2',  &
               'Grid-cell mean liquid water path for liquid water clouds before microphysics',  flag_xyfill=.true.)
   call addfld ('iwpbf', horiz_only ,'A'    ,   'kg/m2',  &
               'Grid-cell mean ice water path for ice clouds before microphysics',  flag_xyfill=.true.)
   call addfld ('cdnumbf', horiz_only ,'A'  ,  '#/m2',  &
              'Grid-cell mean column-integrated droplet number concentrations before microphysics',  flag_xyfill=.true.)
   call addfld ('icnumbf', horiz_only ,'A'    ,  '#/m2',  &
              'Grid-cell mean column-integrated ice crystal number concentrations before microphysics',  flag_xyfill=.true.)

   call addfld ('aod400',horiz_only ,'A'    ,  '#',   'Aerosol optical depth at 400 nm',  flag_xyfill=.true.)
   call addfld ('aod700', horiz_only ,'A'    , '#',   'Aerosol optical depth at 700 nm',  flag_xyfill=.true.)

   call addfld('colccn.1',horiz_only ,'A'    ,'#/m2',   'Column-integrated CCN concentration at S=0.1%')
   call addfld('colccn.3',horiz_only ,'A'    ,'#/m2',   'Column-integrated CCN concentration at S=0.3%')
   call addfld('ccn.1bl',horiz_only ,'A'    ,'#/m3',   'CCN concentration at S=0.1% at 1km above surface')
   call addfld('ccn.3bl',horiz_only ,'A'    ,'#/m3',   'CCN concentration at S=0.3% at 1km above surface')

   call addfld('lwc', (/ 'lev' /) ,'A'    ,   'kg/m3  ', 'Cloud liquid water content')
   call addfld('iwc', (/ 'lev' /) ,'A'    ,   'kg/m3  ', 'Cloud ice water content')
   call addfld('nc',  (/ 'lev' /) ,'A'    ,  '#/m3  ', 'Cloud liquid droplet number concentration')
   call addfld('ni', (/ 'lev' /) ,'A'    ,   '#/m3  ', 'Cloud ice crystal number concentration')
   call addfld('airmass',(/ 'lev' /) ,'A'    ,   'kg/m2', 'Atmosphere mass content of air ')
   call addfld('zaltitude', (/ 'lev' /) ,'A'    ,   'm', 'Altitiude from ground')
   call addfld('dz', (/ 'lev' /) ,'A',   'm', 'Layer thickness')

   cldo_idx     = pbuf_get_index('AST')
   rel_idx      = pbuf_get_index('REL') 
   rei_idx      = pbuf_get_index('REI')

   call cnst_get_ind('CLDLIQ',ixcldliq)
   call cnst_get_ind('CLDICE',ixcldice)
   call cnst_get_ind('NUMLIQ',ixnumliq)
   call cnst_get_ind('NUMICE',ixnumice)

   return
   end subroutine output_aerocom_aie_init
!=================================================================================

!---------------------------------------------------------------------------------
   subroutine cloud_top_aerocom(state, pbuf)
!-----------------------------------------------------------------------------------
!  Purpose: to sample the variables at the cloud top for both ice and liquid clouds. 
!
!  Source: the original codes is  provided by AEROCOM (AIE intercomparison)
!
!------------------------------------------------------------------------------------

   use physics_types,   only: physics_state
   use physics_buffer, only : physics_buffer_desc, pbuf_get_field
   use physconst,     only: gravit, rair, cappa
   use cam_logfile,      only: iulog

   type(physics_state), intent(in), target :: state
   type(physics_buffer_desc), pointer      :: pbuf(:)
   

!  Local variables
   integer, parameter  ::  iovl = 3   ! Overlap assumption for clouds: 1, maximum; 2, random; 3, max/random.
   
   real(r8) :: cld     (pcols,pver)          ! cloud fractin, the same with cldo, except that it is 0 
   real(r8) :: cldliq  (pcols, pver) 
   real(r8) :: cldice  (pcols, pver)
   real(r8) :: numliq  (pcols, pver)
   real(r8) :: numice  (pcols, pver)
   real(r8) :: rhoair(pcols, pver)   ! air density (kg/m3)
   real(r8) :: lnd  (pcols, pver)    ! in-cloud droplet number concentrations 
   real(r8) :: ind  (pcols, pver)    ! in-cloud ice crystal number concentrations 
   real(r8) :: lwpbf  (pcols)
   real(r8) :: iwpbf  (pcols)
   real(r8) :: cdnumbf  (pcols)
   real(r8) :: icnumbf  (pcols)

   real(r8) :: airmass (pcols, pver) 
   real(r8) :: zaltitude (pcols, pver)
   real(r8) :: dz(pcols, pver)
   real(r8) :: lwc (pcols, pver)
   real(r8) :: iwc (pcols, pver)
   real(r8) :: nc(pcols, pver)
   real(r8) :: ni(pcols, pver)

   real(r8), pointer, dimension(:, :) :: cldo 
   real(r8), pointer, dimension(:, :) :: rel    ! Liquid cloud particle effective radius
   real(r8), pointer, dimension(:, :) :: rei    ! Ice effective drop size (microns)
   real(r8), pointer, dimension(:, :) :: tau    ! cloud optical depth
   real(r8), pointer, dimension(:, :) :: qrain  ! grid-mean rain water mixing ratio

   real(r8), pointer, dimension(:) ::  clt 
   real(r8), pointer, dimension(:) ::  icc
   real(r8), pointer, dimension(:) ::  lcc
   real(r8), pointer, dimension(:) ::  ttop 
   real(r8), pointer, dimension(:) ::  ptop
   real(r8), pointer, dimension(:) ::  cdr 
   real(r8), pointer, dimension(:) ::  icr
   real(r8), pointer, dimension(:) ::  cdnc 
   real(r8), pointer, dimension(:) ::  icnc 
   real(r8), pointer, dimension(:) ::  ccn

   real(r8), pointer, dimension(:) :: cod
   real(r8), pointer, dimension(:) :: lwp
   real(r8), pointer, dimension(:) :: iwp
   real(r8), pointer, dimension(:) :: cdnum
   real(r8), pointer, dimension(:) :: icnum
   real(r8), pointer, dimension(:) :: autoconv
   real(r8), pointer, dimension(:) :: accretn
   
   real(r8), pointer, dimension(:) :: intccn
   real(r8), pointer, dimension(:) :: colrv
   real(r8), pointer, dimension(:) :: rwp
   real(r8), pointer, dimension(:) :: lwp2
   real(r8), pointer, dimension(:) :: iwp2

   real(r8), pointer, dimension(:, :) :: autocl
   real(r8), pointer, dimension(:, :) :: accretl 
   real(r8), pointer, dimension(:, :) :: ccn3d

   real(r8), pointer, dimension(:, :) :: cldliqbf, cldicebf, numliqbf, numicebf


   real(r8) :: ftmp(pcols)
   real(r8) :: zi(pver+1) 
   real(r8) :: flag_max 
   real(r8) :: thres_cld, thres_cod, thres_cwp, max_cld, fr, cwp
   real(r8) :: rv, t700
   real(r8) :: iiwp, ilwp, iphase, lphase
   
   integer :: itrue
   integer :: lchnk, ncol
   integer :: i, k
   
   thres_cld = 0.001
   thres_cwp = 1.0e-5  
   max_cld = 1.0_r8 - 1.0e-9
!   thres_cod = 0.3
   thres_cod = 1.0e-5 

   lchnk = state%lchnk
   ncol  = state%ncol

   cldliq(:ncol, :) = state%q(:ncol,:,ixcldliq)
   cldice(:ncol, :) = state%q(:ncol,:,ixcldice)
   numliq(:ncol, :) = state%q(:ncol,:,ixnumliq)
   numice(:ncol, :) = state%q(:ncol,:,ixnumice)
   rhoair(:ncol, :) = state%pmid(:ncol,:) / (rair * state%t(:ncol,:))

   call pbuf_get_field(pbuf, cldo_idx, cldo)
   call pbuf_get_field(pbuf, qrain_idx, qrain)
   call pbuf_get_field(pbuf, cld_tau_idx, tau)
   call pbuf_get_field(pbuf, rel_idx, rel)
   call pbuf_get_field(pbuf, rei_idx, rei)
   call pbuf_get_field(pbuf, ccn3d_idx, ccn3d)
   call pbuf_get_field(pbuf, autocl_idx, autocl)
   call pbuf_get_field(pbuf, accretl_idx, accretl)

   call pbuf_get_field(pbuf, cldliqbf_idx, cldliqbf)
   call pbuf_get_field(pbuf, cldicebf_idx, cldicebf)
   call pbuf_get_field(pbuf, numliqbf_idx, numliqbf)
   call pbuf_get_field(pbuf, numicebf_idx, numicebf)

   call pbuf_get_field(pbuf, clt_idx, clt)
   call pbuf_get_field(pbuf, icc_idx, icc)
   call pbuf_get_field(pbuf, lcc_idx, lcc)
   call pbuf_get_field(pbuf, ttop_idx, ttop)
   call pbuf_get_field(pbuf, ptop_idx, ptop)
   call pbuf_get_field(pbuf, cdr_idx, cdr)
   call pbuf_get_field(pbuf, icr_idx, icr)
   call pbuf_get_field(pbuf, cdnc_idx, cdnc)
   call pbuf_get_field(pbuf, icnc_idx, icnc)
   call pbuf_get_field(pbuf, ccn_idx, ccn)
   call pbuf_get_field(pbuf, cod_idx, cod)
   call pbuf_get_field(pbuf, lwp_idx, lwp)
   call pbuf_get_field(pbuf, iwp_idx, iwp)
   call pbuf_get_field(pbuf, cdnum_idx, cdnum)
   call pbuf_get_field(pbuf, icnum_idx, icnum)
   call pbuf_get_field(pbuf, autoconv_idx, autoconv) 
   call pbuf_get_field(pbuf, accretn_idx, accretn)
   call pbuf_get_field(pbuf, intccn_idx, intccn)
   call pbuf_get_field(pbuf, colrv_idx, colrv)
   call pbuf_get_field(pbuf, rwp_idx, rwp)
   call pbuf_get_field(pbuf, lwp2_idx, lwp2)
   call pbuf_get_field(pbuf, iwp2_idx, iwp2)

! calcluate in-cloud droplet number concentrations
! and ice crystal number concentrations
   do i=1, ncol
    do k=1, pver
      if(cldo(i,k).gt.thres_cld) then
         lnd(i,k) = numliq(i,k) * rhoair(i,k)/cldo(i,k)
         ind(i,k) = numice(i,k) * rhoair(i,k)/cldo(i,k)
      else
         lnd(i,k) = 0.0
         ind(i,k) = 0.0
      end if
    end do
   end do

   cld(:, :) = cldo(:, :)
   if ( iovl.eq.2 .or.iovl.eq.3 ) then 
    clt(:) = 1._r8
   else 
    clt(:) = 0.0_r8
   end if 
   icc(:) = 0.0_r8
   lcc(:) = 0.0_r8
   ttop(:) = 0.0_r8
   ptop(:) = 0.0_r8
   cdr(:) = 0.0_r8
   icr(:) = 0.0_r8
   cdnc(:) = 0.0_r8
   icnc(:) = 0.0_r8
   ccn(:) = 0.0_r8
   lwp2(:) = 0.0_r8
   iwp2(:) = 0.0_r8

   do i=1, ncol
    do k=2, pver ! assumption: uppermost layer is cloud-free (k=1)
       cwp =  (cldliq(i,k)+cldice(i,k)) * state%pdel(i,k)/gravit 
       cwp = cwp/max(0.001_r8, cld(i,k))
!       if ( tau(i,k).ge.thres_cod.and.cld(i,k).ge.thres_cld ) then ! visible, not-too-small cloud 
       if ( cwp.ge.thres_cwp.and.cld(i,k).ge.thres_cld ) then ! not-too-small cloud, use cwp intead of tau, as at nightly grids, tau is not defined. 
!       flag_max is needed since the vertical integration for maximum overlap is different from 
!       the two others: for maximum, clt is the actual cloud cover in the level, for the two others, the actual cloud cover is 1 - clt
!       ftmp is total cloud cover seen from above down to the current level
!       clt is ftmp from the level just above
!       ftmp - clt is thus the additional cloud fraction seen from above in this level
        if (iovl.eq.1 ) then 
          flag_max = -1._r8
          ftmp(i) = max(clt(i), cld(i,k))  ! maximum overlap	
        elseif ( iovl.eq.2 ) then 
          flag_max = 1._r8
          ftmp(i) = clt(i) * ( 1._r8 - cld(i,k) ) ! random overlap	
        elseif ( iovl.eq. 3 ) then 
          flag_max = 1._r8
!          ftmp(i) = clt(i) * ( 1._r8 - min( max( cld(i,k), cld(i,k-1) ), max_cld) ) /                &
!            ( 1.0_r8 - min( cld(i,k-1), max_cld ) )  ! maximum-random overlap
          fr= ( 1._r8 - min( max( cld(i,k), cld(i,k-1) ), max_cld) ) /                &
            ( 1.0_r8 - min( cld(i,k-1), max_cld ) )  ! maximum-random overlap
          fr = min(1.0_r8, fr)
          if(fr.gt.1.0_r8) then
            write(6, *) 'cloud overlap',  cld(i,k), cld(i, k-1), ( 1._r8 - min( max( cld(i,k), cld(i,k-1) ), max_cld) ), &
                       ( 1.0_r8 - min( cld(i,k-1), max_cld ) ) 
            call endrun('cloud overlap')
          endif
          fr = min(1.0_r8, fr)
          ftmp(i) = clt(i) * fr
        endif 

        ttop(i) = ttop(i) + state%t(i,k) * ( clt(i) - ftmp(i) )*flag_max 
        ptop(i) = ptop(i) + state%pmid(i,k) *( clt(i) - ftmp(i) )*flag_max

        if(cldice(i,k).gt.1.0e-8) then
          iphase = 1.0
          iiwp = cldice(i,k) * state%pdel(i,k)/gravit 
        else 
          iphase = 0.0
          iiwp = 0.0
        end if

        if(cldliq(i,k).gt.1.0e-8) then
          lphase = 1.0
          ilwp = cldliq(i,k) * state%pdel(i,k)/gravit
        else 
          lphase = 0.0
          ilwp = 0.0
        end if

!       ice clouds
        icr(i) = icr(i) + rei(i,k) * iphase * ( clt(i) - ftmp(i) )*flag_max  * 1.0e-6  ! um --> m
        icc(i) = icc(i) + iphase * ( clt(i) - ftmp(i) )*flag_max 
        icnc(i) = icnc(i) + ind(i,k) * iphase * ( clt(i) - ftmp(i) )*flag_max
        iwp2(i) = iwp2(i) + iiwp 

!       liquid water clouds
        cdr(i) = cdr(i) + rel(i,k) * lphase * ( clt(i) - ftmp(i) )*flag_max  * 1.0e-6  ! um -> m
        cdnc(i) = cdnc(i) + lnd(i,k) * lphase * ( clt(i) - ftmp(i) )*flag_max 
        ccn(i)  = ccn(i) + ccn3d(i,k) * lphase * ( clt(i) - ftmp(i)) * flag_max * 1.0e6  ! #/cm3 --> #/m3
        lcc(i) = lcc(i) + lphase * ( clt(i) - ftmp(i) )*flag_max
        lwp2(i) = lwp2(i) + ilwp
		
        clt(i) = ftmp(i)
       else 
        cld(i,k) = 0.0_r8   ! reset the cloud fraction to be zero. This is necessary for maximum-random overlap  
                            ! for cld(i, k-1) is used there. Otherwise you will find maximum-random overlap gives
                            ! smaller cloud fraction in some cases. 
       end if ! is there a visible, not-too-small cloud?
    end do ! loop over k

    if ( iovl.eq.2 .or. iovl.eq.3 ) then 
      clt(i) = 1._r8 - clt(i)
    end if 
!    if ( clt(i).le.thres_cld) then
!      ttop(i) = fillvalue
!      ptop(i) = fillvalue
!    else 
!      ttop(i) = ttop(i)/clt(i)
!      ptop(i) = ptop(i)/clt(i)
!    end if
!    if (lcc(i).le.thres_cld) then
!      cdr(i)  = fillvalue
!      cdnc(i) = fillvalue
!      ccn(i)  = fillvalue
!      lwp2(i) = fillvalue
!    else
!      cdr(i) = cdr(i)/lcc(i) * 1.0e-6  ! micron meter -> meter
!      cdnc(i,lchnk) = cdnc(i)/lcc(i) * 1.0e6  ! #/cm3 --> #/m3
!      ccn(i) = ccn(i)/lcc(i)
!      lwp2(i) = lwp2(i)/lcc(i) 
!    end if 
!    if (icc(i).le.thres_cld) then
!      icr(i)  = fillvalue
!      icnc(i) = fillvalue
!      iwp2(i) = fillvalue
!    else
!      icr(i) = icr(i)/icc(i) * 1.0e-6    ! micron meter -> meter
!      icnc(i) = icnc(i)/icc(i)  * 1.0e6  ! #/cm3 --> #/m3
!      iwp2(i) = iwp2(i)/icc(i) 
!    end if
   end do ! loop over I

! Diagnose other variables 
   do i=1, ncol 
     lwp(i) = 0.0
     iwp(i) = 0.0
     colrv(i) = 0.0
     cod(i) = 0.0
     intccn(i) = 0.0
     rwp(i) = 0.0
     cdnum(i) = 0.0
     icnum(i) = 0.0
     autoconv(i) = 0.0
     accretn(i) = 0.0
     lwpbf(i)=0.0
     iwpbf(i)=0.0
     cdnumbf(i)=0.0
     icnumbf(i)=0.0
     do k=1, pver
       lwp(i) = lwp(i) + cldliq(i,k) * state%pdel(i,k)/gravit
       iwp(i) = iwp(i) + cldice(i,k) * state%pdel(i,k)/gravit
       rwp(i) = rwp(i) + qrain(i,k) * state%pdel(i,k)/gravit
       cdnum(i) = cdnum(i) + numliq(i, k) * state%pdel(i,k)/gravit             
       icnum(i) = icnum(i) + numice(i, k) * state%pdel(i,k)/gravit  
       rv = min(25._r8, max(1.0_r8, (3 * cldliq(i,k)/max(4*3.14159*numliq(i,k)*1.0e3, 1.0_r8))**(1/3.) * 1.0e6))  ! micrometer
       colrv(i) = colrv(i) + rv * cldliq(i,k) * state%pdel(i,k)/gravit  

       cod(i) = cod(i) + tau(i,k) * cld(i,k)
       intccn(i) = intccn(i)+ccn3d(i,k) * (state%pdel(i,k)/gravit/rhoair(i,k)) * 1.0e6 ! #/cm3 -->#/m2
       autoconv(i) = autoconv(i) + autocl(i,k) * state%pdel(i,k)/gravit
       accretn(i) = accretn(i) + accretl(i,k) * state%pdel(i,k)/gravit

       lwpbf(i) = lwpbf(i) + cldliqbf(i,k) * state%pdel(i,k)/gravit
       iwpbf(i) = iwpbf(i) + cldicebf(i,k) * state%pdel(i,k)/gravit
       cdnumbf(i) = cdnumbf(i) + numliqbf(i, k) * state%pdel(i,k)/gravit
       icnumbf(i) = icnumbf(i) + numicebf(i, k) * state%pdel(i,k)/gravit

     end do
!     if(lcc(i).gt.thres_cld) then

!       colrv(i) = colrv(i) / max(0.00001, lwp(i))
!       cdnum(i) = cdnum(i) / max(0.00001, lwp(i))

!       lwp(i) = lwp(i) / max(clt(i), 0.01)
!       rwp(i) = rwp(i) / max(clt(i), 0.01) 
!     else
!       lwp(i) = fillvalue
!       rwp(i) = fillvalue
!       colrv(i) = fillvalue
!       cdnum(i) = fillvalue
!     end if
!     if(icc(i).gt.thres_cld) then
!       iwp(i) = iwp(i) / max(clt(i), 0.01)
!     else
!       iwp(i) = fillvalue
!     end if
!     if(clt(i).gt.thres_cld) then
!       cod(i) = cod(i) / max(clt(i), 0.01)
!     else
!       cod(i) = fillvalue
!     end if
   end do

   call outfld('cdr', cdr, pcols, lchnk) 
   call outfld('cdnc', cdnc, pcols, lchnk)
   call outfld('cdnum', cdnum, pcols, lchnk)
   call outfld('icnum', icnum, pcols, lchnk)
   call outfld('clt', clt, pcols, lchnk)
   call outfld('lcc', lcc, pcols, lchnk)
   call outfld('lwp', lwp, pcols, lchnk)
   call outfld('iwp', iwp, pcols, lchnk)
   call outfld('icc', icc, pcols, lchnk)
   call outfld('icnc', icnc, pcols, lchnk)
   call outfld('icr', icr, pcols, lchnk)
   call outfld('cod', cod, pcols, lchnk)
   call outfld('ccn', ccn, pcols, lchnk)
   call outfld('ptop', ptop, pcols, lchnk)
   call outfld('ttop', ttop, pcols, lchnk)
   call outfld('intccn', intccn, pcols, lchnk)
   call outfld('colrv', colrv, pcols, lchnk)
   call outfld('rwp', rwp, pcols, lchnk)
   call outfld('lwp2', lwp2, pcols, lchnk)
   call outfld('iwp2', iwp2, pcols, lchnk)

   call outfld('autoconv', autoconv, pcols, lchnk)
   call outfld('accretn', accretn, pcols, lchnk)

   call outfld('lwpbf', lwpbf, pcols, lchnk)
   call outfld('iwpbf', iwpbf, pcols, lchnk)
   call outfld('cdnumbf', cdnumbf, pcols, lchnk)
   call outfld('icnumbf', icnumbf, pcols, lchnk)

   zaltitude = 0.0
   do i=1, ncol
     zi = 0.0
     do k=pver, 1, -1
       airmass(i,k) = state%pdel(i,k)/gravit     ! kg/m2
       dz(i, k) = airmass(i,k)/rhoair(i,k) ! layer thickness in m
       zi(k) = zi(k+1) + dz(i,k) ! layer thickness in m
       zaltitude(i, k) = (zi(k+1)+zi(k))/2._r8
       lwc(i,k) = cldliq(i,k) * rhoair(i,k)  ! kg/kg --> kg/m3
       iwc(i,k) = cldice(i,k) * rhoair(i,k)  ! kg/kg --> kg/m3
       nc(i,k) = numliq(i,k) * rhoair(i,k)  ! #/kg --> #/m3
       ni(i,k) = numice(i,k) * rhoair(i,k)  ! #/kg --> #/m3
     end do
   end do
   call outfld('lwc', lwc, pcols, lchnk)
   call outfld('iwc', iwc, pcols, lchnk)
   call outfld('nc',  nc,  pcols, lchnk)
   call outfld('ni',  ni,  pcols, lchnk)
   call outfld('airmass', airmass, pcols, lchnk)
   call outfld('zaltitude', zaltitude, pcols, lchnk)
   call outfld('dz', dz, pcols, lchnk)

   return
   end subroutine cloud_top_aerocom
!==========================================================================

!--------------------------------------------------------------------------
   subroutine aerocom_calc(state, cld, rel, rei, lnd, ind, tau, cliqwp, cicewp, coszrs)
!-------------------------------------------------------------------
!  calculate required variables.
!------------------------------------------------------------------
   use physics_types,   only: physics_state

   implicit none
!!#include <comctl.h>

   type(physics_state), intent(in) :: state
   real(r8), intent(in) :: cld(pcols,pver)     ! cloud cover
   real(r8), intent(in) :: rel(pcols,pver)     ! Liquid cloud particle effective radius (microns)
   real(r8), intent(in) :: rei(pcols,pver)     ! Ice effective drop size (microns) 
   real(r8), intent(in) :: lnd(pcols, pver)    ! Liquid cloud number concentration (#/cm3)
   real(r8), intent(in) :: ind(pcols, pver)    ! Liquid cloud number concentration (#/cm3)
   real(r8), intent(in) :: tau(pcols, pver)    ! cloud optical depth 
   real(r8), intent(in) :: cliqwp(pcols,pver)  ! in-cloud liquid water path  (g/m2)
   real(r8), intent(in) :: cicewp(pcols,pver)  ! in-cloud ice water path (g/m2)
   real(r8), intent(in) :: coszrs(pcols)       ! cosine solar zenith angle (to tell if day or night)

!  Local variables
   integer :: lchnk, ncol
   integer :: i, k

   lchnk = state%lchnk
   ncol  = state%ncol

!   call cloud_top_aerocom(state, pbuf)

   return
   end subroutine aerocom_calc
!-----------------------------------------------------------------------------

!=============================================================================
   subroutine output_aerocom (state, pbuf)
   use physics_types,   only: physics_state
   use physics_buffer, only : physics_buffer_desc, pbuf_get_field

   type(physics_state), intent(in), target :: state
   type(physics_buffer_desc), pointer      :: pbuf(:)

!  Local variables
   
!   if(dosw) then
!     call outfld('OD550',   od550(:),   pcols, lchnk)
!     call outfld('ALBS',    albs(:),    pcols, lchnk)
!     call outfld('RST',     rst(:),     pcols, lchnk)
!     call outfld('RSTCS',   rstcs(:),   pcols, lchnk)
!     call outfld('RSS',     rss(:),     pcols, lchnk)
!     call outfld('RSSCS',   rsscs(:),     pcols, lchnk)
!     call outfld('RSDS',    rsds(:),    pcols, lchnk)
!   end if
!   if (dolw) then
!     call outfld('RLT',     rlt(:),     pcols, lchnk)
!     call outfld('RLTCS',   rltcs(:),   pcols, lchnk)
!     call outfld('RLS',     rls(:),     pcols, lchnk)
!     call outfld('RLSCS',     rlscs(:),     pcols, lchnk)
!   end if

!   call outfld('HFLS',    hfls(:),    pcols, lchnk)
!   call outfld('HFSS',    hfss(:),    pcols, lchnk)

   return
   end subroutine output_aerocom
!--------------------------------------------------------------------------

end module  output_aerocom_aie 
