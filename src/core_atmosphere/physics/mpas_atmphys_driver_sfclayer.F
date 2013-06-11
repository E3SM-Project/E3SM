! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================
 module mpas_atmphys_driver_sfclayer
 use mpas_grid_types

 use mpas_atmphys_constants
 use mpas_atmphys_vars

!from wrf physics:
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

 contains

!==================================================================================================
 subroutine allocate_sfclayer
!==================================================================================================

 if(.not.allocated(area_p)   ) allocate(area_p(ims:ime,jms:jme)   )
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

 if(allocated(area_p)   ) deallocate(area_p   )
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
 type(mesh_type),intent(in):: mesh
 type(sfc_input_type),intent(in):: sfc_input
 type(diag_physics_type),intent(inout):: diag_physics

!--------------------------------------------------------------------------------------------------

 do j = jts,jte
 do i = its,ite
    !input variables:
    area_p(i,j)   = mesh % areaCell % array(i)
    hpbl_p(i,j)   = diag_physics % hpbl   % array(i)
    mavail_p(i,j) = diag_physics % mavail % array(i)
    tsk_p(i,j)    = sfc_input    % skintemp % array(i)
    xland_p(i,j)  = sfc_input    % xland    % array(i)       
    !inout variables:
    br_p(i,j)     = diag_physics % br     % array(i)
    cpm_p(i,j)    = diag_physics % cpm    % array(i)
    chs_p(i,j)    = diag_physics % chs    % array(i)
    chs2_p(i,j)   = diag_physics % chs2   % array(i)
    cqs2_p(i,j)   = diag_physics % cqs2   % array(i)
    fh_p(i,j)     = diag_physics % fh     % array(i)
    fm_p(i,j)     = diag_physics % fm     % array(i)
    flhc_p(i,j)   = diag_physics % flhc   % array(i)
    flqc_p(i,j)   = diag_physics % flqc   % array(i)
    gz1oz0_p(i,j) = diag_physics % gz1oz0 % array(i)
    hfx_p(i,j)    = diag_physics % hfx    % array(i)
    qfx_p(i,j)    = diag_physics % qfx    % array(i)
    qgh_p(i,j)    = diag_physics % qgh    % array(i)
    qsfc_p(i,j)   = diag_physics % qsfc   % array(i) 
    lh_p(i,j)     = diag_physics % lh     % array(i)
    mol_p(i,j)    = diag_physics % mol    % array(i) 
    psim_p(i,j)   = diag_physics % psim   % array(i)
    psih_p(i,j)   = diag_physics % psih   % array(i)
    regime_p(i,j) = diag_physics % regime % array(i)
    rmol_p(i,j)   = diag_physics % rmol   % array(i)
    ust_p(i,j)    = diag_physics % ust    % array(i)
    wspd_p(i,j)   = diag_physics % wspd   % array(i)
    znt_p(i,j)    = diag_physics % znt    % array(i) 
    zol_p(i,j)    = diag_physics % zol    % array(i) 
    !output variables:
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
 type(diag_physics_type),intent(inout):: diag_physics

!--------------------------------------------------------------------------------------------------

 do j = jts,jte
 do i = its,ite
    diag_physics % br     % array(i) = br_p(i,j)
    diag_physics % cpm    % array(i) = cpm_p(i,j)
    diag_physics % chs    % array(i) = chs_p(i,j)
    diag_physics % chs2   % array(i) = chs2_p(i,j)
    diag_physics % cqs2   % array(i) = cqs2_p(i,j)
    diag_physics % fh     % array(i) = fh_p(i,j)
    diag_physics % fm     % array(i) = fm_p(i,j)
    diag_physics % flhc   % array(i) = flhc_p(i,j)
    diag_physics % flqc   % array(i) = flqc_p(i,j)
    diag_physics % gz1oz0 % array(i) = gz1oz0_p(i,j)
    diag_physics % hfx    % array(i) = hfx_p(i,j)
    diag_physics % lh     % array(i) = lh_p(i,j)
    diag_physics % mol    % array(i) = mol_p(i,j)
    diag_physics % qfx    % array(i) = qfx_p(i,j)
    diag_physics % qgh    % array(i) = qgh_p(i,j)
    diag_physics % qsfc   % array(i) = qsfc_p(i,j)
    diag_physics % psim   % array(i) = psim_p(i,j)
    diag_physics % psih   % array(i) = psih_p(i,j)
    diag_physics % regime % array(i) = regime_p(i,j)
    diag_physics % rmol   % array(i) = rmol_p(i,j)
    diag_physics % ust    % array(i) = ust_p(i,j)
    diag_physics % wspd   % array(i) = wspd_p(i,j)
    diag_physics % zol    % array(i) = zol_p(i,j)
    diag_physics % znt    % array(i) = znt_p(i,j)
    !diagnostics:
    diag_physics % q2     % array(i) = q2_p(i,j)
    diag_physics % t2m    % array(i) = t2m_p(i,j)
    diag_physics % th2m   % array(i) = th2m_p(i,j)
    diag_physics % u10    % array(i) = u10_p(i,j)
    diag_physics % v10    % array(i) = v10_p(i,j)
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
!--------------------------
 type(mesh_type)        ,intent(in):: mesh
 type(sfc_input_type)   ,intent(in):: sfc_input
 type(diag_physics_type),intent(inout):: diag_physics

!local variables:
!----------------
 real(kind=RKIND):: dx

!--------------------------------------------------------------------------------------------------
 write(0,*)
 write(0,*) '--- enter subroutine driver_sfclayer:'

!copy all MPAS arrays to rectanguler grid:
 call sfclayer_from_MPAS(mesh,diag_physics,sfc_input)

 sfclayer_select: select case (trim(sfclayer_scheme))

    case("monin_obukhov")
       dx = sqrt(maxval(mesh % areaCell % array))
#if defined(do_hydrostatic_pressure)
!... REARRANGED CALL USING HYDROSTATIC PRESSURE:
       call sfclay( &
                   p3d      = pres_hyd_p , psfc     = psfc_p     , t3d      = t_p        , &
                   u3d      = u_p        , v3d      = v_p        , qv3d     = qv_p       , &
                   dz8w     = dz_p       , cp       = cp         , g        = g          , &
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
                   areaCell = area_p     , ustm     = ustm_p     , ck       = ck_p       , & 
                   cka      = cka_p      , cd       = cd_p       , cda      = cda_p      , & 
                   isftcflx = isftcflx   , iz0tlnd  = iz0tlnd    ,                         &
                   scm_force_flux = scm_force_flux               ,                         & 
                   ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde , &
                   ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme , &
                   its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte   &
                 )
#else
!... REARRANGED CALL:
       call sfclay( &
                   p3d      = pres_p     , psfc     = psfc_p     , t3d      = t_p        , &
                   u3d      = u_p        , v3d      = v_p        , qv3d     = qv_p       , &
                   dz8w     = dz_p       , cp       = cp         , g        = g          , &
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
                   areaCell = area_p     , ustm     = ustm_p     , ck       = ck_p       , & 
                   cka      = cka_p      , cd       = cd_p       , cda      = cda_p      , & 
                   isftcflx = isftcflx   , iz0tlnd  = iz0tlnd    ,                         &
                   scm_force_flux = scm_force_flux               ,                         & 
                   ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde , &
                   ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme , &
                   its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte   &
                 )
#endif
        
    case default

 end select sfclayer_select

!copy all arrays back to the MPAS grid:
 call sfclayer_to_MPAS(diag_physics)

 write(0,*) '--- end subroutine driver_sfclayer'

 end subroutine driver_sfclayer

!==================================================================================================
 end module mpas_atmphys_driver_sfclayer
!==================================================================================================
