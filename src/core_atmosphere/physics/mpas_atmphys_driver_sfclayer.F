! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================
 module mpas_atmphys_driver_sfclayer
 use mpas_kind_types
 use mpas_grid_types

 use mpas_atmphys_constants
 use mpas_atmphys_vars

!wrf physics:
 use module_sf_sfclay

 implicit none
 private
 public:: init_sfclayer,       &
          allocate_sfclayer,   &
          deallocate_sfclayer, &
          driver_sfclayer

 integer,parameter,private:: isfflx   = 1        !=1 for surface heat and moisture fluxes.
 integer,parameter,private:: isftcflx = 0        !=0,(Charnock and Carlson-Boland).
 integer,parameter,private:: iz0tlnd  = 0        !=0,(Carlson-Boland).
 integer,parameter,private:: scm_force_flux = 0  !SCM surface forcing by surface fluxes.
                                                 !0=no 1=yes (WRF single column model option only).

 integer,private:: i,j


!>\brief MPAS driver for parameterization of the surface layer.
!>\author Laura D. Fowler (send comments to laura@ucar.edu).
!>\date 2013-05-01.
!>
!>\details
!>
!> subroutines in mpas_atmphys_driver_sfclayer:
!> --------------------------------------------
!> allocate_sfclayer    : allocate local arrays for parameterization of surface layer.
!> deallocate_sfclayer  : deallocate local arrays for parameterization of surface layer.
!> init_sfclayer        : initialization of individual surface layer schemes.
!> driver_sfclayer      : main driver (called from subroutine physics_driver).
!> sfclayer_from_MPAS   : initialize local arrays.
!> sfclayer_to_MPAS     : copy local arrays to MPAS arrays.
!>
!> WRF physics called from driver_sfclayer:
!> ----------------------------------------
!>    * module_sf_sfclay: Monin-Obukhov surface layer scheme.
!>
!> add-ons and modifications to sourcecode:
!> ----------------------------------------
!>    * removed the pre-processor option "do_hydrostatic_pressure" before call to the subroutine
!>      sfclay.
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2013-05-29.
!>    * updated the definition of the horizontal resolution to the actual mean distance between
!>      cell centers.
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2013-08-23.
!>    * in call to subroutine sfclay, replaced the variable g (that originally pointed to gravity)
!>      with gravity, for simplicity.
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2014-03-21.
!>    * in subroutine sfclayer_from_MPAS, added initialization of ustm, cd, cda, ck, and cka. in
!>      subroutine sfclayer_to_MPAS, filled diag_physics%ustm with ustm_p after call to subroutine
!>      sfclay.
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2014-04-16. 
!>    * throughout the sourcecode, replaced all "var_struct" defined arrays by local pointers.
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2014-04-22.
!>    * Modified sourcecode to use pools.
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2014-05-15.


 contains


!==================================================================================================
 subroutine allocate_sfclayer
!==================================================================================================

 if(.not.allocated(dx_p)     ) allocate(dx_p(ims:ime,jms:jme)     )
 if(.not.allocated(br_p)     ) allocate(br_p(ims:ime,jms:jme)     )
 if(.not.allocated(cd_p)     ) allocate(cd_p(ims:ime,jms:jme)     )
 if(.not.allocated(cda_p)    ) allocate(cda_p(ims:ime,jms:jme)    )
 if(.not.allocated(chs_p)    ) allocate(chs_p(ims:ime,jms:jme)    )
 if(.not.allocated(chs2_p)   ) allocate(chs2_p(ims:ime,jms:jme)   )
 if(.not.allocated(ck_p)     ) allocate(ck_p(ims:ime,jms:jme)     )
 if(.not.allocated(cka_p)    ) allocate(cka_p(ims:ime,jms:jme)    )
 if(.not.allocated(cpm_p)    ) allocate(cpm_p(ims:ime,jms:jme)    )
 if(.not.allocated(cqs2_p)   ) allocate(cqs2_p(ims:ime,jms:jme)   )
 if(.not.allocated(gz1oz0_p) ) allocate(gz1oz0_p(ims:ime,jms:jme) )
 if(.not.allocated(flhc_p)   ) allocate(flhc_p(ims:ime,jms:jme)   )
 if(.not.allocated(flqc_p)   ) allocate(flqc_p(ims:ime,jms:jme)   )
 if(.not.allocated(fh_p)     ) allocate(fh_p(ims:ime,jms:jme)     )
 if(.not.allocated(fm_p)     ) allocate(fm_p(ims:ime,jms:jme)     )
 if(.not.allocated(hfx_p)    ) allocate(hfx_p(ims:ime,jms:jme)    )
 if(.not.allocated(hpbl_p)   ) allocate(hpbl_p(ims:ime,jms:jme)   )
 if(.not.allocated(lh_p)     ) allocate(lh_p(ims:ime,jms:jme)     )
 if(.not.allocated(mavail_p) ) allocate(mavail_p(ims:ime,jms:jme) )
 if(.not.allocated(mol_p)    ) allocate(mol_p(ims:ime,jms:jme)    )
 if(.not.allocated(psih_p)   ) allocate(psih_p(ims:ime,jms:jme)   )
 if(.not.allocated(psim_p)   ) allocate(psim_p(ims:ime,jms:jme)   )
 if(.not.allocated(q2_p)     ) allocate(q2_p(ims:ime,jms:jme)     )
 if(.not.allocated(qfx_p)    ) allocate(qfx_p(ims:ime,jms:jme)    )
 if(.not.allocated(qgh_p)    ) allocate(qgh_p(ims:ime,jms:jme)    )
 if(.not.allocated(qsfc_p)   ) allocate(qsfc_p(ims:ime,jms:jme)   )
 if(.not.allocated(regime_p) ) allocate(regime_p(ims:ime,jms:jme) )
 if(.not.allocated(rmol_p)   ) allocate(rmol_p(ims:ime,jms:jme)   )
 if(.not.allocated(t2m_p)    ) allocate(t2m_p(ims:ime,jms:jme)    )
 if(.not.allocated(tsk_p)    ) allocate(tsk_p(ims:ime,jms:jme)    )
 if(.not.allocated(th2m_p)   ) allocate(th2m_p(ims:ime,jms:jme)   )
 if(.not.allocated(u10_p)    ) allocate(u10_p(ims:ime,jms:jme)    )
 if(.not.allocated(ust_p)    ) allocate(ust_p(ims:ime,jms:jme)    )
 if(.not.allocated(ustm_p)   ) allocate(ustm_p(ims:ime,jms:jme)   )
 if(.not.allocated(v10_p)    ) allocate(v10_p(ims:ime,jms:jme)    )
 if(.not.allocated(wspd_p)   ) allocate(wspd_p(ims:ime,jms:jme)   )
 if(.not.allocated(xland_p)  ) allocate(xland_p(ims:ime,jms:jme)  )
 if(.not.allocated(zol_p)    ) allocate(zol_p(ims:ime,jms:jme)    )
 if(.not.allocated(znt_p)    ) allocate(znt_p(ims:ime,jms:jme)    )

 end subroutine allocate_sfclayer

!==================================================================================================
 subroutine deallocate_sfclayer
!==================================================================================================

 if(allocated(dx_p)     ) deallocate(dx_p     )
 if(allocated(br_p)     ) deallocate(br_p     )
 if(allocated(cd_p)     ) deallocate(cd_p     )
 if(allocated(cda_p)    ) deallocate(cda_p    )
 if(allocated(chs_p)    ) deallocate(chs_p    )
 if(allocated(chs2_p)   ) deallocate(chs2_p   )
 if(allocated(ck_p)     ) deallocate(ck_p     )
 if(allocated(cka_p)    ) deallocate(cka_p    )
 if(allocated(cpm_p)    ) deallocate(cpm_p    )
 if(allocated(cqs2_p)   ) deallocate(cqs2_p   )
 if(allocated(gz1oz0_p) ) deallocate(gz1oz0_p )
 if(allocated(flhc_p)   ) deallocate(flhc_p   )
 if(allocated(flqc_p)   ) deallocate(flqc_p   )
 if(allocated(fh_p)     ) deallocate(fh_p     )
 if(allocated(fm_p)     ) deallocate(fm_p     )
 if(allocated(hfx_p)    ) deallocate(hfx_p    )
 if(allocated(hpbl_p)   ) deallocate(hpbl_p   )
 if(allocated(lh_p)     ) deallocate(lh_p     )
 if(allocated(mavail_p) ) deallocate(mavail_p )
 if(allocated(mol_p)    ) deallocate(mol_p    )
 if(allocated(psih_p)   ) deallocate(psih_p   )
 if(allocated(psim_p)   ) deallocate(psim_p   )
 if(allocated(q2_p)     ) deallocate(q2_p     )
 if(allocated(qfx_p)    ) deallocate(qfx_p    )
 if(allocated(qgh_p)    ) deallocate(qgh_p    )
 if(allocated(qsfc_p)   ) deallocate(qsfc_p   )
 if(allocated(regime_p) ) deallocate(regime_p )
 if(allocated(rmol_p)   ) deallocate(rmol_p   )
 if(allocated(t2m_p)    ) deallocate(t2m_p    )
 if(allocated(tsk_p)    ) deallocate(tsk_p    )
 if(allocated(th2m_p)   ) deallocate(th2m_p   )
 if(allocated(u10_p)    ) deallocate(u10_p    )
 if(allocated(ust_p)    ) deallocate(ust_p    )
 if(allocated(ustm_p)   ) deallocate(ustm_p   )
 if(allocated(v10_p)    ) deallocate(v10_p    )
 if(allocated(wspd_p)   ) deallocate(wspd_p   )
 if(allocated(xland_p)  ) deallocate(xland_p  )
 if(allocated(zol_p)    ) deallocate(zol_p    )
 if(allocated(znt_p)    ) deallocate(znt_p    )

 end subroutine deallocate_sfclayer

!==================================================================================================
 subroutine sfclayer_from_MPAS(mesh,diag_physics,sfc_input)
!==================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: sfc_input
 type(mpas_pool_type),intent(inout):: diag_physics

!local pointers:
 real(kind=RKIND),dimension(:),pointer:: skintemp,xland
 real(kind=RKIND),dimension(:),pointer:: dcEdge_m,hpbl,mavail
 real(kind=RKIND),dimension(:),pointer:: br,cpm,chs,chs2,cqs2,fh,fm,flhc,flqc,gz1oz0,hfx,    &
                                         qfx,qgh,qsfc,lh,mol,psim,psih,regime,rmol,ust,ustm, &
                                         wspd,znt,zol

!--------------------------------------------------------------------------------------------------

 call mpas_pool_get_array(sfc_input   ,'skintemp',skintemp)
 call mpas_pool_get_array(sfc_input   ,'xland'   ,xland   )
 call mpas_pool_get_array(diag_physics,'dcEdge_m',dcEdge_m)
 call mpas_pool_get_array(diag_physics,'hpbl'    ,hpbl    )
 call mpas_pool_get_array(diag_physics,'mavail'  ,mavail  )
 call mpas_pool_get_array(diag_physics,'br'      ,br      )
 call mpas_pool_get_array(diag_physics,'cpm'     ,cpm     )
 call mpas_pool_get_array(diag_physics,'chs'     ,chs     )
 call mpas_pool_get_array(diag_physics,'chs2'    ,chs2    )
 call mpas_pool_get_array(diag_physics,'cqs2'    ,cqs2    )
 call mpas_pool_get_array(diag_physics,'fh'      ,fh      )
 call mpas_pool_get_array(diag_physics,'fm'      ,fm      )
 call mpas_pool_get_array(diag_physics,'flhc'    ,flhc    )
 call mpas_pool_get_array(diag_physics,'flqc'    ,flqc    )
 call mpas_pool_get_array(diag_physics,'gz1oz0'  ,gz1oz0  )
 call mpas_pool_get_array(diag_physics,'hfx'     ,hfx     )
 call mpas_pool_get_array(diag_physics,'qfx'     ,qfx     )
 call mpas_pool_get_array(diag_physics,'qgh'     ,qgh     )
 call mpas_pool_get_array(diag_physics,'qsfc'    ,qsfc    )
 call mpas_pool_get_array(diag_physics,'lh'      ,lh      )
 call mpas_pool_get_array(diag_physics,'mol'     ,mol     )
 call mpas_pool_get_array(diag_physics,'psim'    ,psim    )
 call mpas_pool_get_array(diag_physics,'psih'    ,psih    )
 call mpas_pool_get_array(diag_physics,'regime'  ,regime  )
 call mpas_pool_get_array(diag_physics,'rmol'    ,rmol    )
 call mpas_pool_get_array(diag_physics,'ust'     ,ust     )
 call mpas_pool_get_array(diag_physics,'ustm'    ,ustm    )
 call mpas_pool_get_array(diag_physics,'wspd'    ,wspd    )
 call mpas_pool_get_array(diag_physics,'znt'     ,znt     )
 call mpas_pool_get_array(diag_physics,'zol'     ,zol     )

 do j = jts,jte
 do i = its,ite
    !input variables:
    dx_p(i,j)     = dcEdge_m(i)
    hpbl_p(i,j)   = hpbl(i)
    mavail_p(i,j) = mavail(i)
    tsk_p(i,j)    = skintemp(i)
    xland_p(i,j)  = xland(i)       
    !inout variables:
    br_p(i,j)     = br(i)
    cpm_p(i,j)    = cpm(i)
    chs_p(i,j)    = chs(i)
    chs2_p(i,j)   = chs2(i)
    cqs2_p(i,j)   = cqs2(i)
    fh_p(i,j)     = fh(i)
    fm_p(i,j)     = fm(i)
    flhc_p(i,j)   = flhc(i)
    flqc_p(i,j)   = flqc(i)
    gz1oz0_p(i,j) = gz1oz0(i)
    hfx_p(i,j)    = hfx(i)
    qfx_p(i,j)    = qfx(i)
    qgh_p(i,j)    = qgh(i)
    qsfc_p(i,j)   = qsfc(i) 
    lh_p(i,j)     = lh(i)
    mol_p(i,j)    = mol(i) 
    psim_p(i,j)   = psim(i)
    psih_p(i,j)   = psih(i)
    regime_p(i,j) = regime(i)
    rmol_p(i,j)   = rmol(i)
    ust_p(i,j)    = ust(i)
    ustm_p(i,j)   = ustm(i)
    wspd_p(i,j)   = wspd(i)
    znt_p(i,j)    = znt(i) 
    zol_p(i,j)    = zol(i) 
    !output variables:
    cd_p(i,j)     = 0._RKIND
    cda_p(i,j)    = 0._RKIND
    ck_p(i,j)     = 0._RKIND
    cka_p(i,j)    = 0._RKIND

    q2_p(i,j)     = 0._RKIND
    t2m_p(i,j)    = 0._RKIND
    th2m_p(i,j)   = 0._RKIND
    u10_p(i,j)    = 0._RKIND
    v10_p(i,j)    = 0._RKIND
 enddo
 enddo

 end subroutine sfclayer_from_MPAS

!==================================================================================================
 subroutine sfclayer_to_MPAS(diag_physics)
!==================================================================================================

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics

!local pointers:
 real(kind=RKIND),dimension(:),pointer:: br,cpm,chs,chs2,cqs2,fh,fm,flhc,flqc,gz1oz0,hfx,    &
                                         qfx,qgh,qsfc,lh,mol,psim,psih,regime,rmol,ust,ustm, &
                                         wspd,znt,zol
 real(kind=RKIND),dimension(:),pointer:: q2,t2m,th2m,u10,v10

!--------------------------------------------------------------------------------------------------

 call mpas_pool_get_array(diag_physics,'br'    , br    )
 call mpas_pool_get_array(diag_physics,'cpm'   , cpm   )
 call mpas_pool_get_array(diag_physics,'chs'   , chs   )
 call mpas_pool_get_array(diag_physics,'chs2'  , chs2  )
 call mpas_pool_get_array(diag_physics,'cqs2'  , cqs2  )
 call mpas_pool_get_array(diag_physics,'fh'    , fh    )
 call mpas_pool_get_array(diag_physics,'fm'    , fm    )
 call mpas_pool_get_array(diag_physics,'flhc'  , flhc  )
 call mpas_pool_get_array(diag_physics,'flqc'  , flqc  )
 call mpas_pool_get_array(diag_physics,'gz1oz0', gz1oz0)
 call mpas_pool_get_array(diag_physics,'hfx'   , hfx   )
 call mpas_pool_get_array(diag_physics,'qfx'   , qfx   )
 call mpas_pool_get_array(diag_physics,'qgh'   , qgh   )
 call mpas_pool_get_array(diag_physics,'qsfc'  , qsfc  )
 call mpas_pool_get_array(diag_physics,'lh'    , lh    )
 call mpas_pool_get_array(diag_physics,'mol'   , mol   )
 call mpas_pool_get_array(diag_physics,'psim'  , psim  )
 call mpas_pool_get_array(diag_physics,'psih'  , psih  )
 call mpas_pool_get_array(diag_physics,'regime', regime)
 call mpas_pool_get_array(diag_physics,'rmol'  , rmol  )
 call mpas_pool_get_array(diag_physics,'ust'   , ust   )
 call mpas_pool_get_array(diag_physics,'ustm'  , ustm  )
 call mpas_pool_get_array(diag_physics,'wspd'  , wspd  )
 call mpas_pool_get_array(diag_physics,'znt'   , znt   )
 call mpas_pool_get_array(diag_physics,'zol'   , zol   )

 call mpas_pool_get_array(diag_physics,'q2'    , q2    )
 call mpas_pool_get_array(diag_physics,'t2m'   , t2m   )
 call mpas_pool_get_array(diag_physics,'th2m'  , th2m  )
 call mpas_pool_get_array(diag_physics,'u10'   , u10   )
 call mpas_pool_get_array(diag_physics,'v10'   , v10   )

 do j = jts,jte
 do i = its,ite
    br(i)     = br_p(i,j)
    cpm(i)    = cpm_p(i,j)
    chs(i)    = chs_p(i,j)
    chs2(i)   = chs2_p(i,j)
    cqs2(i)   = cqs2_p(i,j)
    fh(i)     = fh_p(i,j)
    fm(i)     = fm_p(i,j)
    flhc(i)   = flhc_p(i,j)
    flqc(i)   = flqc_p(i,j)
    gz1oz0(i) = gz1oz0_p(i,j)
    hfx(i)    = hfx_p(i,j)
    lh(i)     = lh_p(i,j)
    mol(i)    = mol_p(i,j)
    qfx(i)    = qfx_p(i,j)
    qgh(i)    = qgh_p(i,j)
    qsfc(i)   = qsfc_p(i,j)
    psim(i)   = psim_p(i,j)
    psih(i)   = psih_p(i,j)
    regime(i) = regime_p(i,j)
    rmol(i)   = rmol_p(i,j)
    ust(i)    = ust_p(i,j)
    ustm(i)   = ustm_p(i,j)
    wspd(i)   = wspd_p(i,j)
    zol(i)    = zol_p(i,j)
    znt(i)    = znt_p(i,j)
    !diagnostics:
    q2(i)     = q2_p(i,j)
    t2m(i)    = t2m_p(i,j)
    th2m(i)   = th2m_p(i,j)
    u10(i)    = u10_p(i,j)
    v10(i)    = v10_p(i,j)
 enddo
 enddo

 end subroutine sfclayer_to_MPAS

!==================================================================================================
 subroutine init_sfclayer
!==================================================================================================

!local variables:
 logical:: allowed_to_read

!--------------------------------------------------------------------------------------------------
 write(0,*)
 write(0,*) '--- enter sfclayer_initialization:'
 sfclayer_select: select case (trim(sfclayer_scheme))

    case("monin_obukhov")
       write(0,*) '    enter monin_obukhov initialization:'
       call sfclayinit(allowed_to_read)
       write(0,*) '    end monin_obukhov initialization'

    case default

 end select sfclayer_select
 write(0,*) '--- end sfclayer_initialization'

 end subroutine init_sfclayer

!==================================================================================================
 subroutine driver_sfclayer(mesh,diag_physics,sfc_input)
!==================================================================================================

!input and inout arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: sfc_input
 type(mpas_pool_type),intent(inout):: diag_physics

!local pointers:
 real(kind=RKIND),dimension(:),pointer:: areaCell

!local variables:
 real(kind=RKIND):: dx

!--------------------------------------------------------------------------------------------------
 write(0,*)
 write(0,*) '--- enter subroutine driver_sfclayer:'

 call mpas_pool_get_array(mesh,'areaCell',areaCell)

!copy all MPAS arrays to rectanguler grid:
 call sfclayer_from_MPAS(mesh,diag_physics,sfc_input)

 sfclayer_select: select case (trim(sfclayer_scheme))

    case("monin_obukhov")
       dx = sqrt(maxval(areaCell))
       call sfclay( &
                   p3d      = pres_hyd_p , psfc     = psfc_p     , t3d      = t_p        , &
                   u3d      = u_p        , v3d      = v_p        , qv3d     = qv_p       , &
                   dz8w     = dz_p       , cp       = cp         , g        = gravity    , &
                   rovcp    = rcp        , R        = R_d        , xlv      = xlv        , & 
                   chs      = chs_p      , chs2     = chs2_p     , cqs2     = cqs2_p     , &
                   cpm      = cpm_p      , znt      = znt_p      , ust      = ust_p      , &
                   pblh     = hpbl_p     , mavail   = mavail_p   , zol      = zol_p      , &
                   mol      = mol_p      , regime   = regime_p   , psim     = psim_p     , &
                   psih     = psih_p     , fm       = fm_p       , fh       = fh_p       , &
                   xland    = xland_p    , hfx      = hfx_p      , qfx      = qfx_p      , &
                   lh       = lh_p       , tsk      = tsk_p      , flhc     = flhc_p     , &
                   flqc     = flqc_p     , qgh      = qgh_p      , qsfc     = qsfc_p     , &
                   rmol     = rmol_p     , u10      = u10_p      , v10      = v10_p      , &
                   th2      = th2m_p     , t2       = t2m_p      , q2       = q2_p       , &
                   gz1oz0   = gz1oz0_p   , wspd     = wspd_p     , br       = br_p       , &
                   isfflx   = isfflx     , dx       = dx         , svp1     = svp1       , &
                   svp2     = svp2       , svp3     = svp3       , svpt0    = svpt0      , &
                   ep1      = ep_1       , ep2      = ep_2       , karman   = karman     , &
                   eomeg    = eomeg      , stbolt   = stbolt     , P1000mb  = P0         , &
                   dxCell   = dx_p       , ustm     = ustm_p     , ck       = ck_p       , & 
                   cka      = cka_p      , cd       = cd_p       , cda      = cda_p      , & 
                   isftcflx = isftcflx   , iz0tlnd  = iz0tlnd    ,                         &
                   scm_force_flux = scm_force_flux               ,                         & 
                   ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde , &
                   ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme , &
                   its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte   &
                  )
        
    case default

 end select sfclayer_select

!copy local arrays to MPAS grid:
 call sfclayer_to_MPAS(diag_physics)

 write(0,*) '--- end subroutine driver_sfclayer'

 end subroutine driver_sfclayer

!==================================================================================================
 end module mpas_atmphys_driver_sfclayer
!==================================================================================================
