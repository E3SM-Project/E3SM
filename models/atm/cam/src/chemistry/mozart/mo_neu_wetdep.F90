!
! code written by J.-F. Lamarque, S. Walters and F. Vitt
! based on the original code from J. Neu developed for UC Irvine
! model
!
module mo_neu_wetdep
!
  use shr_kind_mod, only : r8 => shr_kind_r8
  use cam_logfile,  only : iulog
  use constituents, only : pcnst
  use spmd_utils,   only : masterproc
  use abortutils,   only : endrun
  use seq_drydep_mod, only : n_species_table, species_name_table, dheff
  use gas_wetdep_opts,only : gas_wetdep_method, gas_wetdep_list, gas_wetdep_cnt
!
  implicit none
!
  private
  public :: neu_wetdep_init
  public :: neu_wetdep_tend
  public :: do_neu_wetdep
!
  save
!
  integer, allocatable, dimension(:) :: mapping_to_heff,mapping_to_mmr
  real(r8),allocatable, dimension(:) :: mol_weight
  logical ,allocatable, dimension(:) :: ice_uptake
  integer                     :: index_cldice,index_cldliq,nh3_ndx,co2_ndx
  logical                     :: debug   = .false.
  integer                     :: hno3_ndx = 0
!
! diagnostics
!
  logical                     :: do_diag = .false.
  integer, parameter          :: kdiag = 18
!
  real(r8), parameter :: zero = 0._r8
  real(r8), parameter :: one  = 1._r8
!
  logical :: do_neu_wetdep
!
contains

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine neu_wetdep_init
!
  use constituents, only : cnst_get_ind,cnst_mw
  use cam_history,  only : addfld, add_default, phys_decomp
  use ppgrid,       only : pver
!
  integer :: m,l
  character*20 :: test_name

  do_neu_wetdep = gas_wetdep_method == 'NEU' .and. gas_wetdep_cnt>0

  if (.not.do_neu_wetdep) return

  allocate( mapping_to_heff(gas_wetdep_cnt) )
  allocate( mapping_to_mmr(gas_wetdep_cnt) )
  allocate( ice_uptake(gas_wetdep_cnt) )
  allocate( mol_weight(gas_wetdep_cnt) )

!
! find mapping to heff table
!
  if ( debug ) then
    print '(a,i4)','gas_wetdep_cnt=',gas_wetdep_cnt
    print '(a,i4)','n_species_table=',n_species_table
  end if
  mapping_to_heff = -99
  do m=1,gas_wetdep_cnt
!
    test_name = gas_wetdep_list(m)
    if ( debug ) print '(i4,a)',m,trim(test_name)
!
! mapping based on the MOZART4 wet removal subroutine;
! this might need to be redone (JFL: Sep 2010)
!
    select case( trim(test_name) )
       
      case( 'HYAC', 'CH3COOH' , 'HCOOH', 'EOOH' )
         test_name = 'CH2O'
      case ( 'SO2','SOGB','SOGI','SOGM','SOGT','SOGX' )
         test_name = 'H2O2'
      case ( 'CLONO2','BRONO2','HCL','HOCL','HOBR','HBR', 'Pb', 'MACROOH', 'ISOPOOH', 'XOOH', 'H2SO4' )
         test_name = 'HNO3'
      case ( 'ALKOOH', 'MEKOOH', 'TOLOOH', 'TERPOOH' )
         test_name = 'CH3OOH'        

    end select
!
    do l = 1,n_species_table
!
!      if ( debug ) print '(i4,a)',l,trim(species_name_table(l))
!
       if( trim(test_name) == trim( species_name_table(l) ) ) then
          mapping_to_heff(m)  = l
          if ( debug ) print '(a,a,i4)','mapping to heff of ',trim(species_name_table(l)),l
          exit
       end if
    end do
    if ( mapping_to_heff(m) == -99 ) then
      if (masterproc) print *,'problem with mapping_to_heff of ',trim(test_name)
!      call endrun()
    end if
!
! special cases for NH3 and CO2
!
    if ( trim(test_name) == 'NH3' ) then
      nh3_ndx = m
    end if
    if ( trim(test_name) == 'CO2' ) then
      co2_ndx = m
    end if
    if ( trim(test_name) == 'HNO3' ) then
      hno3_ndx = m
    end if
!
  end do
   
   if (any ( mapping_to_heff(:) == -99 ))  call endrun('mo_neu_wet->depwetdep_init: unmapped species error' )
!
  if ( debug ) then
    print '(a,i4)','co2_ndx',co2_ndx
    print '(a,i4)','nh3_ndx',nh3_ndx
  end if
!
! find mapping to species
!
  mapping_to_mmr = -99
  do m=1,gas_wetdep_cnt
    if ( debug ) print '(i4,a)',m,trim(gas_wetdep_list(m))
    call cnst_get_ind(gas_wetdep_list(m), mapping_to_mmr(m), abort=.false. )
    if ( debug ) print '(a,i4)','mapping_to_mmr ',mapping_to_mmr(m)
    if ( mapping_to_mmr(m) <= 0 ) then
      print *,'problem with mapping_to_mmr of ',gas_wetdep_list(m)
      call endrun('problem with mapping_to_mmr of '//trim(gas_wetdep_list(m)))
    end if
  end do
!
! define specie-dependent arrays
!
  do m=1,gas_wetdep_cnt
!
    mol_weight     (m) = cnst_mw(mapping_to_mmr(m))
    if ( debug ) print '(i4,a,f8.4)',m,' mol_weight ',mol_weight(m)
    ice_uptake(m) = .false.
    if ( trim(gas_wetdep_list(m)) == 'HNO3' ) then
      ice_uptake(m) = .true.
    end if
!
  end do
!
! indices for cloud quantities
!
  call cnst_get_ind( 'CLDICE', index_cldice )
  call cnst_get_ind( 'CLDLIQ', index_cldliq )
!
! define output
!
  do m=1,gas_wetdep_cnt
    call addfld     ('DTWR_'//trim(gas_wetdep_list(m)),'mol/mol/s',pver, 'A','wet removal Neu scheme tendency',phys_decomp)
    call add_default('DTWR_'//trim(gas_wetdep_list(m)), 1, ' ')
  end do
!
  if ( do_diag ) then
    call addfld     ('QT_RAIN_HNO3','mol/mol/s',pver, 'A','wet removal Neu scheme rain tendency',phys_decomp)
    call addfld     ('QT_RIME_HNO3','mol/mol/s',pver, 'A','wet removal Neu scheme rain tendency',phys_decomp)
    call addfld     ('QT_WASH_HNO3','mol/mol/s',pver, 'A','wet removal Neu scheme rain tendency',phys_decomp)
    call addfld     ('QT_EVAP_HNO3','mol/mol/s',pver, 'A','wet removal Neu scheme rain tendency',phys_decomp)
    call add_default('QT_RAIN_HNO3',1,' ')
    call add_default('QT_RIME_HNO3',1,' ')
    call add_default('QT_WASH_HNO3',1,' ')
    call add_default('QT_EVAP_HNO3',1,' ')
  end if
!
  return
!
end subroutine neu_wetdep_init
!
subroutine neu_wetdep_tend(lchnk,ncol,mmr,pmid,pdel,zint,tfld,delt, &
     prain, nevapr, cld, cmfdqr, wd_tend)
!
  use ppgrid,           only : pcols, pver
  
  use phys_grid,        only : get_area_all_p
  use shr_const_mod,    only : SHR_CONST_REARTH,SHR_CONST_G
  use cam_history,      only : outfld
!
  implicit none
!
  integer,        intent(in)    :: lchnk,ncol
  real(r8),       intent(in)    :: mmr(pcols,pver,pcnst)    ! mass mixing ratio (kg/kg)
  real(r8),       intent(in)    :: pmid(pcols,pver)         ! midpoint pressures (Pa)
  real(r8),       intent(in)    :: pdel(pcols,pver)         ! pressure delta about midpoints (Pa)
  real(r8),       intent(in)    :: zint(pcols,pver+1)       ! interface geopotential height above the surface (m)
  real(r8),       intent(in)    :: tfld(pcols,pver)         ! midpoint temperature (K)
  real(r8),       intent(in)    :: delt                     ! timestep (s)
!  

  real(r8),       intent(in)    :: prain(ncol, pver)
  real(r8),       intent(in)    :: nevapr(ncol, pver)
  real(r8),       intent(in)    :: cld(ncol, pver)
  real(r8),       intent(in)    :: cmfdqr(ncol, pver)

  real(r8),       intent(inout) :: wd_tend(pcols,pver,pcnst)



!
! local arrays and variables
!
  integer :: i,k,l,kk,m,id
  real(r8), parameter                       :: rearth = SHR_CONST_REARTH    ! radius earth (m)
  real(r8), parameter                       :: gravit = SHR_CONST_G         ! m/s^2
  real(r8), dimension(ncol)                 :: area
  real(r8), dimension(ncol,pver)            :: cldice,cldliq,cldfrc,totprec,totevap,delz,delp,p
  real(r8), dimension(ncol,pver)            :: rls,evaprate,mass_in_layer,temp
  real(r8), dimension(ncol,pver,gas_wetdep_cnt) :: trc_mass,heff,dtwr
  real(r8), dimension(ncol,pver,gas_wetdep_cnt) :: wd_mmr
  logical , dimension(gas_wetdep_cnt)           :: tckaqb
  integer , dimension(ncol)                 :: test_flag
!
! arrays for HNO3 diagnostics
!
  real(r8), dimension(ncol,pver)            :: qt_rain,qt_rime,qt_wash,qt_evap
!
! for Henry's law calculations
!
  real(r8), parameter       :: t0     = 298._r8
  real(r8), parameter       :: ph     = 1.e-5_r8
  real(r8), parameter       :: ph_inv = 1._r8/ph
  real(r8)                  :: e298, dhr
  real(r8), dimension(ncol) :: dk1s,dk2s,wrk
!
! from cam/src/physics/cam/stratiform.F90
!

  if (.not.do_neu_wetdep) return
!
! don't do anything if there are no species to be removed
!
  if ( gas_wetdep_cnt == 0 ) return
!
! reset output variables
!
!  wd_tend = 0._r8
!
! get area (in radians square)
!
  call get_area_all_p(lchnk, ncol, area)
  area = area * rearth**2                     ! in m^2
!
! reverse order along the vertical before calling
! J. Neu's wet removal subroutine
!
  do k=1,pver
    kk = pver - k + 1
    do i=1,ncol
!
      mass_in_layer(i,k) = area(i) * pdel(i,kk)/gravit          ! kg
!
      cldice (i,k) = mmr(i,kk,index_cldice)                     ! kg/kg
      cldliq (i,k) = mmr(i,kk,index_cldliq)                     ! kg/kg
      cldfrc (i,k) = cld(i,kk)                                  ! unitless
!
      totprec(i,k) = (prain(i,kk)+cmfdqr(i,kk)) &
                                  * mass_in_layer(i,k)          ! kg/s
      totevap(i,k) = nevapr(i,kk) * mass_in_layer(i,k)          ! kg/s
!
      delz(i,k) = zint(i,kk) - zint(i,kk+1)                     ! in m
!
      temp(i,k) = tfld(i,kk)
!
! convert tracer mass to kg
!
      trc_mass(i,k,:) = mmr(i,kk,mapping_to_mmr(:)) * mass_in_layer(i,k)
!
      delp(i,k) = pdel(i,kk) * 0.01_r8          ! in hPa
      p   (i,k) = pmid(i,kk) * 0.01_r8          ! in hPa
!
    end do
  end do
!
! define array for tendency calculation (on model grid)
!
  dtwr(1:ncol,:,:) = mmr(1:ncol,:,mapping_to_mmr(:))
!
! compute 1) integrated precipitation flux across the interfaces (rls)
!         2) evaporation rate
!
  rls      (:,pver) = 0._r8
  evaprate (:,pver) = 0._r8
  do k=pver-1,1,-1
    rls     (:,k) = max(0._r8,totprec(:,k)-totevap(:,k)+rls(:,k+1))
    !evaprate(:,k) = min(1._r8,totevap(:,k)/(rls(:,k+1)+totprec(:,k)+1.e-36_r8))
    evaprate(:,k) = min(1._r8,totevap(:,k)/(rls(:,k+1)+1.e-36_r8)) 
  end do
!
! compute effective Henry's law coefficients
! code taken from models/drv/shr/seq_drydep_mod.F90
!
  heff = 0._r8
  do k=1,pver
!
    kk = pver - k + 1
!
    wrk(:) = (t0-tfld(1:ncol,kk))/(t0*tfld(1:ncol,kk))
!
    do m=1,gas_wetdep_cnt
!
      l    = mapping_to_heff(m)
      id   = 6*(l - 1)
      e298 = dheff(id+1)
      dhr  = dheff(id+2)
      heff(:,k,m) = e298*exp( dhr*wrk(:) )
      test_flag = -99
      if( dheff(id+3) /= 0._r8 .and. dheff(id+5) == 0._r8 ) then
        e298 = dheff(id+3)
        dhr  = dheff(id+4)
        dk1s(:) = e298*exp( dhr*wrk(:) )
        where( heff(:,k,m) /= 0._r8 )
          heff(:,k,m) = heff(:,k,m)*(1._r8 + dk1s(:)*ph_inv)
        elsewhere
          test_flag = 1
          heff(:,k,m) = dk1s(:)*ph_inv
        endwhere
      end if
!
      if (k.eq.1 .and. maxval(test_flag) > 0 .and. debug ) print '(a,i4)','heff for m=',m
!
      if( dheff(id+5) /= 0._r8 ) then
        if( nh3_ndx > 0 .or. co2_ndx > 0 ) then
          e298 = dheff(id+3)
          dhr  = dheff(id+4)
          dk1s(:) = e298*exp( dhr*wrk(:) )
          e298 = dheff(id+5)
          dhr  = dheff(id+6)
          dk2s(:) = e298*exp( dhr*wrk(:) )
          if( m == co2_ndx ) then
             heff(:,k,m) = heff(:,k,m)*(1._r8 + dk1s(:)*ph_inv)*(1._r8 + dk2s(:)*ph_inv)
          else if( m == nh3_ndx ) then
             heff(:,k,m) = heff(:,k,m)*(1._r8 + dk1s(:)*ph/dk2s(:))
          else
             write(iulog,*) 'error in assigning henrys law coefficients'
             write(iulog,*) 'species ',m
          end if
        end if
      end if
!
    end do
  end do
!
! define flag for high effective Henry's law
!
  do m=1,gas_wetdep_cnt
    if ( maxval(heff(:,:,m)) > 1.e4_r8 ) then
      tckaqb(m) = .true.
    else
      tckaqb(m) = .false.
    end if
  end do
!
  if ( debug ) then
    print '(a,50f8.2)','tckaqb     ',tckaqb
    print '(a,50e12.4)','heff      ',heff(1,1,:)
    print '(a,50i4)'  ,'ice_uptake ',ice_uptake
    print '(a,50f8.2)','mol_weight ',mol_weight(:)
    print '(a,50f8.2)','temp       ',temp(1,:)
    print '(a,50f8.2)','p          ',p   (1,:)
  end if
!
! call J. Neu's subroutine
!
  do i=1,ncol
!
    call washo(pver,gas_wetdep_cnt,delt,trc_mass(i,:,:),mass_in_layer(i,:),p(i,:),delz(i,:) &
              ,rls(i,:),cldliq(i,:),cldice(i,:),cldfrc(i,:),temp(i,:),evaprate(i,:) &
              ,area(i),heff(i,:,:),mol_weight(:),tckaqb(:),ice_uptake(:) &
              ,qt_rain(i,:),qt_rime(i,:),qt_wash(i,:),qt_evap(i,:) )
!
  end do
!
! compute tendencies and convert back to mmr
! on original vertical grid
!
  do k=1,pver
    kk = pver - k + 1
    do i=1,ncol
!
! convert tracer mass from kg
!
      wd_mmr(i,kk,:) = trc_mass(i,k,:) / mass_in_layer(i,k)
!
    end do
  end do
!
! tendency calculation (on model grid)
!
  dtwr(1:ncol,:,:) = wd_mmr(1:ncol,:,:) - dtwr(1:ncol,:,:)
  dtwr(1:ncol,:,:) = dtwr(1:ncol,:,:) / delt 
!
! output tendencies
!
  do m=1,gas_wetdep_cnt
    wd_tend(1:ncol,:,mapping_to_mmr(m)) = wd_tend(1:ncol,:,mapping_to_mmr(m)) + dtwr(1:ncol,:,m)
    call outfld( 'DTWR_'//trim(gas_wetdep_list(m)),dtwr(:,:,m),ncol,lchnk )
  end do
!
  if ( do_diag ) then
    call outfld('QT_RAIN_HNO3', qt_rain, ncol, lchnk )
    call outfld('QT_RIME_HNO3', qt_rime, ncol, lchnk )
    call outfld('QT_WASH_HNO3', qt_wash, ncol, lchnk )
    call outfld('QT_EVAP_HNO3', qt_evap, ncol, lchnk )
  end if
!
  return
end subroutine neu_wetdep_tend

!-----------------------------------------------------------------------
!
! Original code from Jessica Neu
! Updated by S. Walters and J.-F. Lamarque (March-April 2011)
!
!-----------------------------------------------------------------------

      subroutine WASHO(LPAR,NTRACE,DTSCAV,QTTJFL,QM,POFL,DELZ,  &
      RLS,CLWC,CIWC,CFR,TEM,EVAPRATE,GAREA,HSTAR,TCMASS,TCKAQB, &
      TCNION, qt_rain, qt_rime, qt_wash, qt_evap)
!
      implicit none

!-----------------------------------------------------------------------
!---p-conde 5.4 (2007)   -----called from main-----
!---called from pmain to calculate rainout and washout of tracers
!---revised by JNEU 8/2007
!---
!-LAER has been removed - no scavenging for aerosols
!-LAER could be used as LWASHTYP
!---WILL THIS WORK FOR T42->T21???????????
!-----------------------------------------------------------------------
      
      integer LPAR, NTRACE
      real(r8),  intent(inout) ::  QTTJFL(LPAR,NTRACE)
      real(r8),  intent(in) :: DTSCAV, QM(LPAR),POFL(LPAR),DELZ(LPAR),GAREA
      real(r8),  intent(in) :: RLS(LPAR),CLWC(LPAR),CIWC(LPAR),CFR(LPAR),TEM(LPAR),      &
                               EVAPRATE(LPAR)
      real(r8),  intent(in) :: HSTAR(LPAR,NTRACE),TCMASS(NTRACE)
      logical ,  intent(in) :: TCKAQB(NTRACE),TCNION(NTRACE) 
!
      real(r8),  intent(inout) :: qt_rain(lpar)
      real(r8),  intent(inout) :: qt_rime(lpar)
      real(r8),  intent(inout) :: qt_wash(lpar)
      real(r8),  intent(inout) :: qt_evap(lpar)
!
      integer I,J,L,N,LE, LM1
      real(r8), dimension(LPAR) :: CFXX
      real(r8), dimension(LPAR) :: QTT, QTTNEW

      real(r8) WRK, RNEW_TST
      real(r8) CLWX
      real(r8) RNEW,RPRECIP,DELTARIMEMASS,DELTARIME,RAMPCT
      real(r8) MASSLOSS
      real(r8) DOR,DNEW,DEMP,COLEFFSNOW,RHOSNOW
      real(r8) WEMP,REMP,RRAIN,RWASH
      real(r8) QTPRECIP,QTRAIN,QTCXA,QTAX,QTOC

      real(r8) FAMA,RAMA,DAMA,FCA,RCA,DCA
      real(r8) FAX,RAX,DAX,FCXA,RCXA,DCXA,FCXB,RCXB,DCXB
      real(r8) RAXADJ,FAXADJ,RAXADJF
      real(r8) QTDISCF,QTDISRIME,QTDISCXA
      real(r8) QTEVAPAXP,QTEVAPAXW,QTEVAPAX
      real(r8) QTWASHAX
      real(r8) QTEVAPCXAP,QTEVAPCXAW,QTEVAPCXA
      real(r8) QTWASHCXA,QTRIMECXA
      real(r8) QTRAINCXA,QTRAINCXB
      real(r8) QTTOPCA,QTTOPAA,QTTOPCAX,QTTOPAAX

      real(r8) AMPCT,AMCLPCT,CLNEWPCT,CLNEWAMPCT,CLOLDPCT,CLOLDAMPCT
      real(r8) RAXLOC,RCXALOC,RCXBLOC,RCALOC,RAMALOC,RCXPCT

      real(r8) QTNETLCXA,QTNETLCXB,QTNETLAX,QTNETL
      real(r8) QTDISSTAR
      

      real(r8), parameter  :: TICE=273._r8
      real(r8), parameter  :: CFMIN=0.1_r8
      real(r8), parameter  :: CWMIN=1.0e-5_r8
      real(r8), parameter  :: DMIN=1.0e-1_r8       !mm
      real(r8), parameter  :: VOLPOW=1._r8/3._r8
      real(r8), parameter  :: RHORAIN=1.0e3_r8     !kg/m3
      real(r8), parameter  :: RHOSNOWFIX=1.0e2_r8     !kg/m3
      real(r8), parameter  :: COLEFFRAIN=0.7_r8
      real(r8), parameter  :: TMIX=258._r8
      real(r8), parameter  :: TFROZ=240._r8
      real(r8), parameter  :: COLEFFAER=0.05_r8
!
! additional work arrays and diagnostics
!
      real(r8) :: rls_wrk(lpar)
      real(r8) :: rnew_wrk(lpar)
      real(r8) :: rca_wrk(lpar)
      real(r8) :: fca_wrk(lpar)
      real(r8) :: rcxa_wrk(lpar)
      real(r8) :: fcxa_wrk(lpar)
      real(r8) :: rcxb_wrk(lpar)
      real(r8) :: fcxb_wrk(lpar)
      real(r8) :: rax_wrk(lpar,2)
      real(r8) :: fax_wrk(lpar,2)
      real(r8) :: rama_wrk(lpar)
      real(r8) :: fama_wrk(lpar)
      real(r8) :: deltarime_wrk(lpar)
      real(r8) :: clwx_wrk(lpar)
      real(r8) :: frc(lpar,3)
      real(r8) :: rlsog(lpar)
!
      logical :: is_hno3
      logical :: rls_flag(lpar)
      logical :: rnew_flag(lpar)
      logical :: cf_trigger(lpar)
      logical :: freezing(lpar)
!
      real(r8), parameter :: four = 4._r8
      real(r8), parameter :: adj_factor = one + 10._r8*epsilon( one )
!
      integer :: LWASHTYP,LICETYP
!
      if ( debug ) then
        print '(a,50f8.2)','tckaqb     ',tckaqb
        print '(a,50e12.4)','hstar     ',hstar(1,:)
        print '(a,50i4)'  ,'ice_uptake ',TCNION
        print '(a,50f8.2)','mol_weight ',TCMASS(:)
        print '(a,50f8.2)','temp       ',tem(:)
        print '(a,50f8.2)','p          ',pofl(:)
      end if

!-----------------------------------------------------------------------
      LE = LPAR-1   
!
      rls_flag(1:le) = rls(1:le) > zero 
      freezing(1:le) = tem(1:le) < tice
      rlsog(1:le) = rls(1:le)/garea
!
species_loop : &
     do N = 1,NTRACE
       QTT(:lpar)    = QTTJFL(:lpar,N)
       QTTNEW(:lpar) = QTTJFL(:lpar,N)
       is_hno3 = n == hno3_ndx
       if( is_hno3 ) then
         qt_rain(:lpar) = zero
         qt_rime(:lpar) = zero
         qt_wash(:lpar) = zero
         qt_evap(:lpar) = zero
         rca_wrk(:lpar) = zero
         fca_wrk(:lpar) = zero
         rcxa_wrk(:lpar) = zero
         fcxa_wrk(:lpar) = zero
         rcxb_wrk(:lpar) = zero
         fcxb_wrk(:lpar) = zero
         rls_wrk(:lpar) = zero
         rnew_wrk(:lpar) = zero
         cf_trigger(:lpar) = .false.
         clwx_wrk(:lpar) = -9999._r8
         deltarime_wrk(:lpar) = -9999._r8
         rax_wrk(:lpar,:) = zero
         fax_wrk(:lpar,:) = zero
       endif
!-----------------------------------------------------------------------
!  calculate scavenging by large-scale stratiform precipitation
!  check whether mass-limited or henry's law
!-----------------------------------------------------------------------
       if( TCKAQB(N) ) then
         LWASHTYP = 1
       else
         LWASHTYP = 2
       end if
!-----------------------------------------------------------------------
!  check whether soluble in ice
!-----------------------------------------------------------------------
       if( TCNION(N) ) then
         LICETYP = 1
       else
         LICETYP = 2
       end if

!-----------------------------------------------------------------------
!  initialization
!-----------------------------------------------------------------------
       QTTOPAA = zero
       QTTOPCA = zero

       RCA  = zero
       FCA  = zero
       DCA  = zero
       RAMA = zero
       FAMA = zero
       DAMA = zero

       AMPCT      = zero
       AMCLPCT    = zero
       CLNEWPCT   = zero
       CLNEWAMPCT = zero
       CLOLDPCT   = zero
       CLOLDAMPCT = zero
!-----------------------------------------------------------------------
!  Check whether precip in top layer - if so, require CF ge 0.2
!-----------------------------------------------------------------------
       if( RLS(LE) > zero ) then
         CFXX(LE) = max( CFMIN,CFR(LE) )
       else
         CFXX(LE) = CFR(LE)
       endif

       rnew_flag(1:le) = .false.

level_loop : &
       do L = LE,1,-1
         LM1  = L - 1
         FAX  = zero
         RAX  = zero
         DAX  = zero
         FCXA = zero
         FCXB = zero
         DCXA = zero
         DCXB = zero
         RCXA = zero
         RCXB = zero

         QTDISCF   = zero
         QTDISRIME = zero
         QTDISCXA  = zero

         QTEVAPAXP = zero
         QTEVAPAXW = zero
         QTEVAPAX  = zero
         QTWASHAX  = zero

         QTEVAPCXAP = zero
         QTEVAPCXAW = zero
         QTEVAPCXA  = zero
         QTRIMECXA  = zero
         QTWASHCXA  = zero
         QTRAINCXA  = zero
         QTRAINCXB  = zero
         
         RAMPCT = zero
         RCXPCT = zero

         RCXALOC = zero
         RCXBLOC = zero
         RAXLOC  = zero
         RAMALOC = zero
         RCALOC  = zero

         RPRECIP       = zero
         DELTARIMEMASS = zero
         DELTARIME     = zero
         DOR           = zero
         DNEW          = zero

         QTTOPAAX = zero
         QTTOPCAX = zero

has_rls : &
         if( rls_flag(l) ) then
!-----------------------------------------------------------------------
!-----Evaporate ambient precip and decrease area-------------------------
!-----If ice, diam=diam falling from above  If rain, diam=4mm (not used)
!-----Evaporate tracer contained in evaporated precip
!-----Can't evaporate more than we start with-----------------------------
!-----Don't do washout until we adjust ambient precip to match Rbot if needed
!------(after RNEW if statements)
!-----------------------------------------------------------------------
           FAX = max( zero,FAMA*(one - evaprate(l)) )
           RAX = RAMA								     !kg/m2/s
           if ( debug ) then
             if( (l == 3 .or. l == 2) ) then
               write(*,*) 'washout: l,rls,fax = ',l,rls(l),fax
             endif
           endif
           if( FAMA > zero ) then
             if( freezing(l) ) then
               DAX = DAMA      !mm
             else
               DAX = four    !mm - not necessary
             endif
           else
             DAX = zero
           endif

           if( RAMA > zero ) then
             QTEVAPAXP = min( QTTOPAA,EVAPRATE(L)*QTTOPAA )
           else
             QTEVAPAXP = zero
           endif
           if( is_hno3 ) then
             rax_wrk(l,1) = rax
             fax_wrk(l,1) = fax
           endif


!-----------------------------------------------------------------------
!  Determine how much the in-cloud precip rate has increased------
!-----------------------------------------------------------------------
           WRK = RAX*FAX + RCA*FCA
           if( WRK > 0._r8 ) then
             RNEW_TST = RLS(L)/(GAREA * WRK)
           else
             RNEW_TST = 10._r8
           endif
           RNEW = RLSOG(L) - (RAX*FAX + RCA*FCA)     !GBA*CF
           rnew_wrk(l) = rnew_tst
           if ( debug ) then
             if( is_hno3 .and. l == kdiag-1 ) then
               write(*,*) ' '
               write(*,*) 'washout: rls,rax,fax,rca,fca'
               write(*,'(1p,5g15.7)') rls(l),rax,fax,rca,fca
               write(*,*) ' '
             endif
           endif
!-----------------------------------------------------------------------
!  if RNEW>0, there is growth and/or new precip formation
!-----------------------------------------------------------------------
has_rnew:  if( rlsog(l) > adj_factor*(rax*fax + rca*fca) ) then
!-----------------------------------------------------------------------
!  Min cloudwater requirement for cloud with new precip
!  Min CF is set at top for LE, at end for other levels
!  CWMIN is only needed for new precip formation - do not need for RNEW<0
!-----------------------------------------------------------------------
             if( cfxx(l) == zero ) then
               if ( do_diag ) then
                 write(*,*) 'cfxx(l) == zero',l
                 write(*,*) qttjfl(:,n)
                 write(*,*) qm(:)
                 write(*,*) pofl(:)
                 write(*,*) delz(:)
                 write(*,*) rls(:)
                 write(*,*) clwc(:)
                 write(*,*) ciwc(:)
                 write(*,*) cfr(:)
                 write(*,*) tem(:)
                 write(*,*) evaprate(:)
                 write(*,*) hstar(:,n)
               end if
!
! if we are here,, that means that there is
! a inconsistency and this will lead to a division
! by 0 later on! This column should then be skipped
!
               QTTJFL(:lpar,n) = QTT(:lpar)
               cycle species_loop
!
!              call endrun()
!
             endif
             rnew_flag(l) = .true.
             CLWX = max( CLWC(L)+CIWC(L),CWMIN*CFXX(L) )
             if( is_hno3 ) then
               clwx_wrk(l) = clwx
             endif
!-----------------------------------------------------------------------
!  Area of old cloud and new cloud
!-----------------------------------------------------------------------
             FCXA = FCA
             FCXB = max( zero,CFXX(L)-FCXA )
!-----------------------------------------------------------------------
!                           ICE
!  For ice and mixed phase, grow precip in old cloud by riming
!  Use only portion of cloudwater in old cloud fraction
!  and rain above old cloud fraction
!  COLEFF from Lohmann and Roeckner (1996), Loss rate from Rotstayn (1997)
!-----------------------------------------------------------------------
is_freezing : &
             if( freezing(l) ) then
               COLEFFSNOW = exp( 2.5e-2_r8*(TEM(L) - TICE) )
               if( TEM(L) <= TFROZ ) then
                 RHOSNOW = RHOSNOWFIX
               else
                 RHOSNOW = 0.303_r8*(TEM(L) - TFROZ)*RHOSNOWFIX
               endif
               if( FCXA > zero ) then
                 if( DCA > zero ) then
                   DELTARIMEMASS = CLWX*QM(L)*(FCXA/CFXX(L))* &
                     (one - exp( (-COLEFFSNOW/(DCA*1.e-3_r8))*((RCA)/(2._r8*RHOSNOW))*DTSCAV ))   !uses GBA R
                 else
                   DELTARIMEMASS = zero
                 endif
               else
                 DELTARIMEMASS = zero
               endif
!-----------------------------------------------------------------------
!  Increase in precip rate due to riming (kg/m2/s):
!  Limit to total increase in R in cloud
!-----------------------------------------------------------------------
               if( FCXA > zero ) then
                 DELTARIME = min( RNEW/FCXA,DELTARIMEMASS/(FCXA*GAREA*DTSCAV) ) !GBA
               else
                 DELTARIME = zero
               endif
               if( is_hno3 ) then
                 deltarime_wrk(l) = deltarime
               endif
!-----------------------------------------------------------------------
!  Find diameter of rimed precip, must be at least .1mm
!-----------------------------------------------------------------------
               if( RCA > zero ) then
                 DOR = max( DMIN,(((RCA+DELTARIME)/RCA)**VOLPOW)*DCA )
               else
                 DOR = zero
               endif
!-----------------------------------------------------------------------
!  If there is some in-cloud precip left, we have new precip formation
!  Will be spread over whole cloud fraction 
!-----------------------------------------------------------------------
!  Calculate precip rate in old and new cloud fractions
!-----------------------------------------------------------------------
               RPRECIP = (RNEW-(DELTARIME*FCXA))/CFXX(L) !kg/m2/s    !GBA
!-----------------------------------------------------------------------
!  Calculate precip rate in old and new cloud fractions
!-----------------------------------------------------------------------
               RCXA = RCA + DELTARIME + RPRECIP          !kg/m2/s GBA
               RCXB = RPRECIP                            !kg/m2/s GBA

!-----------------------------------------------------------------------
!  Find diameter of new precip from empirical relation using Rprecip
!  in given area of box- use density of water, not snow, to convert kg/s
!  to mm/s -> as given in Field and Heymsfield
!  Also calculate diameter of mixed precip,DCXA, from empirical relation
!  using total R in FCXA - this will give larger particles than averaging DOR and
!  DNEW in the next level
!  DNEW and DCXA must be at least .1mm
!-----------------------------------------------------------------------
               if( RPRECIP > zero ) then
                 WEMP = (CLWX*QM(L))/(GAREA*CFXX(L)*DELZ(L)) !kg/m3
                 REMP = RPRECIP/((RHORAIN/1.e3_r8))             !mm/s local
                 DNEW = DEMPIRICAL( WEMP, REMP )
                 if ( debug ) then
                   if( is_hno3 .and. l >= 15 ) then
                     write(*,*) ' '
                     write(*,*) 'washout: wemp,remp.dnew @ l = ',l
                     write(*,'(1p,3g15.7)') wemp,remp,dnew
                     write(*,*) ' '
                   endif
                 endif
                 DNEW = max( DMIN,DNEW )
                 if( FCXB > zero ) then
                   DCXB = DNEW
                 else
                   DCXB = zero
                 endif
               else
                 DCXB = zero
               endif

               if( FCXA > zero ) then
                 WEMP = (CLWX*QM(L)*(FCXA/CFXX(L)))/(GAREA*FCXA*DELZ(L)) !kg/m3
                 REMP = RCXA/((RHORAIN/1.e3_r8))                         !mm/s local
                 DEMP = DEMPIRICAL( WEMP, REMP )
                 DCXA = ((RCA+DELTARIME)/RCXA)*DOR + (RPRECIP/RCXA)*DNEW
                 DCXA = max( DEMP,DCXA )
                 DCXA = max( DMIN,DCXA )
               else
                 DCXA = zero
               endif
               if ( debug ) then
                 if( is_hno3 .and. l >= 15 ) then
                   write(*,*) ' '
                   write(*,*) 'washout: rca,rcxa,deltarime,dor,rprecip,dnew @ l = ',l
                   write(*,'(1p,6g15.7)') rca,rcxa,deltarime,dor,rprecip,dnew 
                   write(*,*) 'washout: dcxa,dcxb,wemp,remp,demp'
                   write(*,'(1p,5g15.7)') dcxa,dcxb,wemp,remp,demp
                   write(*,*) ' '
                 end if
               endif

               if( QTT(L) > zero ) then   
!-----------------------------------------------------------------------
!                       ICE SCAVENGING
!-----------------------------------------------------------------------
!  For ice, rainout only hno3/aerosols using new precip
!  Tracer dissolved given by Kaercher and Voigt (2006) for T<258K
!  For T>258K, use Henry's Law with Retention coefficient
!  Rain out in whole CF
!-----------------------------------------------------------------------
                 if( RPRECIP > zero ) then
                   if( LICETYP == 1 ) then
                     RRAIN = RPRECIP*GAREA                                  !kg/s local
                     call DISGAS( CLWX, CFXX(L), TCMASS(N), HSTAR(L,N), &
                                  TEM(L),POFL(L),QM(L),                 &
                                  QTT(L)*CFXX(L),QTDISCF )
                     call RAINGAS( RRAIN, DTSCAV, CLWX, CFXX(L),        &
                                   QM(L), QTT(L), QTDISCF, QTRAIN )
                     WRK       = QTRAIN/CFXX(L)
                     QTRAINCXA = FCXA*WRK
                     QTRAINCXB = FCXB*WRK
                   elseif( LICETYP == 2 ) then
                     QTRAINCXA = zero
                     QTRAINCXB = zero
                   endif
                   if( debug .and. is_hno3 .and. l == kdiag ) then
                     write(*,*) ' '
                     write(*,*) 'washout: Ice Scavenging'
                     write(*,*) 'washout: qtraincxa, qtraincxb, fcxa, fcxb, qt_rain, cfxx(l), wrk @ level = ',l
                     write(*,'(1p,7g15.7)') qtraincxa, qtraincxb, fcxa, fcxb, qt_rain(l), cfxx(l), wrk
                     write(*,*) ' '
                   endif
                 endif
!-----------------------------------------------------------------------
!  For ice, accretion removal for hno3 and aerosols is propotional to riming, 
!  no accretion removal for gases
!  remove only in mixed portion of cloud
!  Limit DELTARIMEMASS to RNEW*DTSCAV for ice - evaporation of rimed ice to match
!  RNEW precip rate would result in HNO3 escaping from ice (no trapping) 
!-----------------------------------------------------------------------
                 if( DELTARIME > zero ) then
                   if( LICETYP == 1 ) then
                     if( TEM(L) <= TFROZ ) then
                       RHOSNOW = RHOSNOWFIX
                     else
                       RHOSNOW = 0.303_r8*(TEM(L) - TFROZ)*RHOSNOWFIX
                     endif
                     QTCXA = QTT(L)*FCXA
                     call DISGAS( CLWX*(FCXA/CFXX(L)), FCXA, TCMASS(N),   &
                                  HSTAR(L,N), TEM(L), POFL(L),            &
                                  QM(L), QTCXA, QTDISRIME )       
                     QTDISSTAR = (QTDISRIME*QTCXA)/(QTDISRIME + QTCXA)
                     if ( debug ) then
                       if( is_hno3 .and. l >= 15 ) then
                         write(*,*) ' '
                         write(*,*) 'washout: fcxa,dca,rca,qtdisstar @ l = ',l
                         write(*,'(1p,4g15.7)') fcxa,dca,rca,qtdisstar
                         write(*,*) ' '
                       endif
                     endif
                     QTRIMECXA = QTCXA*                             &
                        (one - exp((-COLEFFSNOW/(DCA*1.e-3_r8))*       &
                        (RCA/(2._r8*RHOSNOW))*                         &  !uses GBA R    
                        (QTDISSTAR/QTCXA)*DTSCAV))
                     QTRIMECXA = min( QTRIMECXA, &               
                        ((RNEW*GAREA*DTSCAV)/(CLWX*QM(L)*(FCXA/CFXX(L))))*QTDISSTAR)
                   elseif( LICETYP == 2 ) then
                     QTRIMECXA = zero
                   endif
                 endif
               else
                 QTRAINCXA = zero
                 QTRAINCXB = zero
                 QTRIMECXA = zero
               endif
!-----------------------------------------------------------------------
!  For ice, no washout in interstitial cloud air
!-----------------------------------------------------------------------
               QTWASHCXA = zero
               QTEVAPCXA = zero

!-----------------------------------------------------------------------
!                      RAIN
!  For rain, accretion increases rain rate but diameter remains constant
!  Diameter is 4mm (not used)
!-----------------------------------------------------------------------
             else is_freezing
               if( FCXA > zero ) then
                 DELTARIMEMASS = (CLWX*QM(L))*(FCXA/CFXX(L))*           &
                   (one - exp( -0.24_r8*COLEFFRAIN*((RCA)**0.75_r8)*DTSCAV ))  !local
               else
                 DELTARIMEMASS = zero
               endif
!-----------------------------------------------------------------------
!  Increase in precip rate due to riming (kg/m2/s):
!  Limit to total increase in R in cloud
!-----------------------------------------------------------------------
               if( FCXA > zero ) then
                 DELTARIME = min( RNEW/FCXA,DELTARIMEMASS/(FCXA*GAREA*DTSCAV) ) !GBA
               else
                 DELTARIME = zero
               endif
!-----------------------------------------------------------------------
!  If there is some in-cloud precip left, we have new precip formation 
!-----------------------------------------------------------------------
               RPRECIP = (RNEW-(DELTARIME*FCXA))/CFXX(L)       !GBA

               RCXA = RCA + DELTARIME + RPRECIP            !kg/m2/s GBA
               RCXB = RPRECIP                              !kg/m2/s GBA
               DCXA = FOUR  
               if( FCXB > zero ) then
                 DCXB = FOUR
               else
                 DCXB = zero
               endif
!-----------------------------------------------------------------------
!                         RAIN SCAVENGING
!  For rain, rainout both hno3/aerosols and gases using new precip
!-----------------------------------------------------------------------
               if( QTT(L) > zero ) then
                 if( RPRECIP > zero ) then
                   RRAIN = (RPRECIP*GAREA) !kg/s local
                   call DISGAS( CLWX, CFXX(L), TCMASS(N), HSTAR(L,N), &
                                TEM(L), POFL(L), QM(L),               &
                                QTT(L)*CFXX(L), QTDISCF )
                   call RAINGAS( RRAIN, DTSCAV, CLWX, CFXX(L),        &
                                 QM(L), QTT(L), QTDISCF, QTRAIN )
                   WRK       = QTRAIN/CFXX(L)
                   QTRAINCXA = FCXA*WRK
                   QTRAINCXB = FCXB*WRK
                   if( debug .and. is_hno3 .and. l == kdiag ) then
                     write(*,*) ' '
                     write(*,*) 'washout: Rain Scavenging'
                     write(*,*) 'washout: qtraincxa, qtraincxb, fcxa, fcxb, qt_rain, cfxx(l), wrk @ level = ',l
                     write(*,'(1p,7g15.7)') qtraincxa, qtraincxb, fcxa, fcxb, qt_rain(l), cfxx(l), wrk
                     write(*,*) ' '
                   endif
                 endif
!-----------------------------------------------------------------------
!  For rain, accretion removal is propotional to riming
!  caclulate for hno3/aerosols and gases
!  Remove only in mixed portion of cloud
!  Limit DELTARIMEMASS to RNEW*DTSCAV
!-----------------------------------------------------------------------
                 if( DELTARIME > zero ) then
                   QTCXA = QTT(L)*FCXA
                   call DISGAS( CLWX*(FCXA/CFXX(L)), FCXA, TCMASS(N),    &
                                HSTAR(L,N), TEM(L), POFL(L),             &
                                QM(L), QTCXA, QTDISRIME )
                   QTDISSTAR = (QTDISRIME*QTCXA)/(QTDISRIME + QTCXA)
                   QTRIMECXA = QTCXA*                              &
                      (one - exp(-0.24_r8*COLEFFRAIN*                 &
                      ((RCA)**0.75_r8)*                               & !local 
                      (QTDISSTAR/QTCXA)*DTSCAV))               
                   QTRIMECXA = min( QTRIMECXA, &
                      ((RNEW*GAREA*DTSCAV)/(CLWX*QM(L)*(FCXA/CFXX(L))))*QTDISSTAR)
                 else
                   QTRIMECXA = zero
                 endif
               else
                 QTRAINCXA = zero
                 QTRAINCXB = zero
                 QTRIMECXA = zero
               endif
!-----------------------------------------------------------------------
!  For rain, washout gases and HNO3/aerosols using rain from above old cloud
!  Washout for HNO3/aerosols is only on non-dissolved portion, impaction-style
!  Washout for gases is on non-dissolved portion, limited by QTTOP+QTRIME
!-----------------------------------------------------------------------
               if( RCA > zero ) then
                 QTPRECIP = FCXA*QTT(L) - QTDISRIME
                 if( LWASHTYP == 1 ) then
                   if( QTPRECIP > zero ) then
                     QTWASHCXA = QTPRECIP*(one - exp( -0.24_r8*COLEFFAER*((RCA)**0.75_r8)*DTSCAV ))   !local
                   else
                     QTWASHCXA = zero
                   endif
                   QTEVAPCXA = zero
                 elseif( LWASHTYP == 2 ) then
                   RWASH = RCA*GAREA                                !kg/s local
                   if( QTPRECIP > zero ) then
                     call WASHGAS( RWASH, FCA, DTSCAV, QTTOPCA+QTRIMECXA, &
                                   HSTAR(L,N), TEM(L), POFL(L),           &
                                   QM(L), QTPRECIP, QTWASHCXA, QTEVAPCXA )
                   else
                     QTWASHCXA = zero
                     QTEVAPCXA = zero
                   endif
                 endif
               endif
             endif is_freezing
!-----------------------------------------------------------------------
!  If RNEW<O, confine precip to area of cloud above
!  FCXA does not require a minimum (could be zero if R(L).le.what
!  evaporated in ambient)
!-----------------------------------------------------------------------
           else has_rnew
             CLWX = CLWC(L) + CIWC(L)
             if( is_hno3 ) then
               clwx_wrk(l) = clwx
             endif
             FCXA = FCA
             FCXB = max( zero,CFXX(L)-FCXA )
             RCXB = zero
             DCXB = zero
             QTRAINCXA = zero
             QTRAINCXB = zero
             QTRIMECXA = zero

!-----------------------------------------------------------------------
!  Put rain into cloud up to RCA so that we evaporate
!  from ambient first
!  Adjust ambient to try to match RLS(L)
!  If no cloud, RAX=R(L)
!-----------------------------------------------------------------------
             if( FCXA > zero ) then
               RCXA = min( RCA,RLS(L)/(GAREA*FCXA) )     !kg/m2/s  GBA
               if( FAX > zero .and. ((RCXA+1.e-12_r8) < RLS(L)/(GAREA*FCXA)) ) then
                 RAXADJF = RLS(L)/GAREA - RCXA*FCXA
                 RAMPCT = RAXADJF/(RAX*FAX)
                 FAXADJ = RAMPCT*FAX
                 if( FAXADJ > zero ) then
                   RAXADJ = RAXADJF/FAXADJ
                 else
                   RAXADJ = zero
                 endif
               else
                 RAXADJ = zero
                 RAMPCT = zero
                 FAXADJ = zero
               endif
             else
               RCXA = zero
               if( FAX > zero ) then
                 RAXADJF = RLS(L)/GAREA
                 RAMPCT = RAXADJF/(RAX*FAX)
                 FAXADJ = RAMPCT*FAX
                 if( FAXADJ > zero ) then
                   RAXADJ = RAXADJF/FAXADJ
                 else
                   RAXADJ = zero
                 endif              
               else
                 RAXADJ = zero
                 RAMPCT = zero
                 FAXADJ = zero
               endif
             endif
  
             QTEVAPAXP = min( QTTOPAA,QTTOPAA - (RAMPCT*(QTTOPAA-QTEVAPAXP)) )
             FAX = FAXADJ
             RAX = RAXADJ
             if ( debug ) then
               if( (l == 3 .or. l == 2) ) then
                 write(*,*) 'washout: l,fcxa,fax = ',l,fcxa,fax
               endif
             endif

!-----------------------------------------------------------------------
!                IN-CLOUD EVAPORATION/WASHOUT
!  If precip out the bottom of the cloud is 0, evaporate everything
!  If there is no cloud, QTTOPCA=0, so nothing happens
!-----------------------------------------------------------------------
             if( RCXA <= zero ) then
               QTEVAPCXA = QTTOPCA
               RCXA = zero
               DCXA = zero
             else
!-----------------------------------------------------------------------
!  If rain out the bottom of the cloud is >0 (but .le. RCA):
!  For ice, decrease particle size,
!  no washout
!  no evap for non-ice gases (b/c there is nothing in ice)
!  T<Tmix,release hno3& aerosols
!  release is amount dissolved in ice mass released
!  T>Tmix, hno3&aerosols are incorporated into ice structure:
!  do not release
!  For rain, assume full evaporation of some raindrops
!  proportional evaporation for all species 
!  washout for gases using Rbot 
!  impact washout for hno3/aerosol portion in gas phase              
!-----------------------------------------------------------------------
!              if (TEM(L) < TICE ) then
is_freezing_a : &
               if( freezing(l) ) then
                 QTWASHCXA = zero
                 DCXA = ((RCXA/RCA)**VOLPOW)*DCA
                 if( LICETYP == 1 ) then
                   if( TEM(L) <= TMIX ) then
                     MASSLOSS = (RCA-RCXA)*FCXA*GAREA*DTSCAV
!-----------------------------------------------------------------------
!  note-QTT doesn't matter b/c T<258K
!-----------------------------------------------------------------------
                     call DISGAS( (MASSLOSS/QM(L)), FCXA, TCMASS(N),   &
                                   HSTAR(L,N), TEM(L), POFL(L),        &
                                   QM(L), QTT(L), QTEVAPCXA )
                     QTEVAPCXA = min( QTTOPCA,QTEVAPCXA )
                   else
                     QTEVAPCXA = zero
                   endif
                 elseif( LICETYP == 2 ) then   
                   QTEVAPCXA = zero
                 endif
               else is_freezing_a
                 QTEVAPCXAP = (RCA - RCXA)/RCA*QTTOPCA
                 DCXA = FOUR
                 QTCXA = FCXA*QTT(L)
                 if( LWASHTYP == 1 ) then
                   if( QTT(L) > zero ) then
                     call DISGAS( CLWX*(FCXA/CFXX(L)), FCXA, TCMASS(N),   &
                                  HSTAR(L,N), TEM(L), POFL(L),            &
                                  QM(L), QTCXA, QTDISCXA )
                     if( QTCXA > QTDISCXA ) then
                       QTWASHCXA = (QTCXA - QTDISCXA)*(one - exp( -0.24_r8*COLEFFAER*((RCXA)**0.75_r8)*DTSCAV )) !local
                     else
                       QTWASHCXA = zero
                     endif
                     QTEVAPCXAW = zero
                   else
                     QTWASHCXA  = zero
                     QTEVAPCXAW = zero
                   endif
                 elseif (LWASHTYP == 2 ) then
                   RWASH = RCXA*GAREA                         !kg/s local
                   call WASHGAS( RWASH, FCXA, DTSCAV, QTTOPCA, HSTAR(L,N), &
                                 TEM(L), POFL(L), QM(L),                   &
                                 QTCXA-QTDISCXA, QTWASHCXA, QTEVAPCXAW )
                 endif
                 QTEVAPCXA = QTEVAPCXAP + QTEVAPCXAW
               endif is_freezing_a
             endif
           endif has_rnew

!-----------------------------------------------------------------------
!                 AMBIENT WASHOUT
!  Ambient precip is finalized - if it is rain, washout
!  no ambient washout for ice, since gases are in vapor phase
!-----------------------------------------------------------------------
           if( RAX > zero ) then
             if( .not. freezing(l) ) then
               QTAX = FAX*QTT(L)
               if( LWASHTYP == 1 ) then
                 QTWASHAX = QTAX*                        &
                    (one - exp(-0.24_r8*COLEFFAER*       &
                   ((RAX)**0.75_r8)*DTSCAV))  !local
                 QTEVAPAXW = zero
               elseif( LWASHTYP == 2 ) then
                 RWASH = RAX*GAREA   !kg/s local
                 call WASHGAS( RWASH, FAX, DTSCAV, QTTOPAA, HSTAR(L,N), &
                               TEM(L), POFL(L), QM(L), QTAX,            &
                               QTWASHAX, QTEVAPAXW )
               endif
             else
               QTEVAPAXW = zero
               QTWASHAX  = zero
             endif
           else
             QTEVAPAXW = zero
             QTWASHAX  = zero
           endif
           QTEVAPAX = QTEVAPAXP + QTEVAPAXW

!-----------------------------------------------------------------------
!                  END SCAVENGING
!  Require CF if our ambient evaporation rate would give less 
!  precip than R from model.
!-----------------------------------------------------------------------
           if( do_diag .and. is_hno3 ) then
             rls_wrk(l) = rls(l)/garea
             rca_wrk(l) = rca
             fca_wrk(l) = fca
             rcxa_wrk(l) = rcxa
             fcxa_wrk(l) = fcxa
             rcxb_wrk(l) = rcxb
             fcxb_wrk(l) = fcxb
             rax_wrk(l,2) = rax
             fax_wrk(l,2) = fax
           endif
upper_level : &
           if( L > 1 ) then
             FAMA = max( FCXA + FCXB + FAX - CFR(LM1),zero )
             if( FAX > zero ) then
               RAXLOC = RAX/FAX
             else
               RAXLOC = zero
             endif
             if( FCXA > zero ) then
               RCXALOC = RCXA/FCXA
             else
               RCXALOC = zero
             endif
             if( FCXB > zero ) then
               RCXBLOC = RCXB/FCXB
             else
               RCXBLOC = zero
             endif

             if( CFR(LM1) >= CFMIN ) then
               CFXX(LM1) = CFR(LM1)
             else
               if( adj_factor*RLSOG(LM1) >= (RCXA*FCXA + RCXB*FCXB + RAX*FAX)*(one - EVAPRATE(LM1)) ) then
                 CFXX(LM1) = CFMIN
                 cf_trigger(lm1) = .true.
               else
                 CFXX(LM1) = CFR(LM1)
               endif
               if( is_hno3 .and. lm1 == kdiag .and. debug ) then
                 write(*,*) ' '
                 write(*,*) 'washout: rls,garea,rcxa,fcxa,rcxb,fcxb,rax,fax'
                 write(*,'(1p,8g15.7)') rls(lm1),garea,rcxa,fcxa,rcxb,fcxb,rax,fax
                 write(*,*) ' '
               endif
             endif
!-----------------------------------------------------------------------
!  Figure out what will go into ambient and cloud below
!  Don't do for lowest level
!-----------------------------------------------------------------------
             if( FAX > zero ) then
               RAXLOC = RAX/FAX
               AMPCT = max( zero,min( one,(CFXX(L) + FAX - CFXX(LM1))/FAX ) )
               AMCLPCT = one - AMPCT
             else
               RAXLOC  = zero
               AMPCT   = zero
               AMCLPCT = zero
             endif
             if( FCXB > zero ) then
               RCXBLOC = RCXB/FCXB
               CLNEWPCT = max( zero,min( (CFXX(LM1) - FCXA)/FCXB,one ) )
               CLNEWAMPCT = one - CLNEWPCT
             else
               RCXBLOC    = zero
               CLNEWPCT   = zero
               CLNEWAMPCT = zero
             endif
             if( FCXA > zero ) then
               RCXALOC = RCXA/FCXA
               CLOLDPCT = max( zero,min( CFXX(LM1)/FCXA,one ) )
               CLOLDAMPCT = one - CLOLDPCT
             else
               RCXALOC    = zero
               CLOLDPCT   = zero
               CLOLDAMPCT = zero
             endif
!-----------------------------------------------------------------------
!  Remix everything for the next level
!-----------------------------------------------------------------------
             FCA = min( CFXX(LM1),FCXA*CLOLDPCT + CLNEWPCT*FCXB + AMCLPCT*FAX )
             if( FCA > zero ) then
!-----------------------------------------------------------------------
!  Maintain cloud core by reducing NC and AM area going into cloud below
!-----------------------------------------------------------------------
               RCA = (RCXA*FCXA*CLOLDPCT + RCXB*FCXB*CLNEWPCT + RAX*FAX*AMCLPCT)/FCA
               if ( debug ) then
                 if( is_hno3 ) then
                   write(*,*) ' '
                   write(*,*) 'washout: rcxa,fcxa,cloldpctrca,rca,fca,dcxa @ l = ',l
                   write(*,'(1p,6g15.7)') rcxa,fcxa,cloldpct,rca,fca,dcxa
                   write(*,*) 'washout: rcxb,fcxb,clnewpct,dcxb'
                   write(*,'(1p,4g15.7)') rcxb,fcxb,clnewpct,dcxb
                   write(*,*) 'washout: rax,fax,amclpct,dax'
                   write(*,'(1p,4g15.7)') rax,fax,amclpct,dax
                   write(*,*) ' '
                 endif
               endif
               DCA = (RCXA*FCXA*CLOLDPCT)/(RCA*FCA)*DCXA + & 
                     (RCXB*FCXB*CLNEWPCT)/(RCA*FCA)*DCXB + &
                     (RAX*FAX*AMCLPCT)/(RCA*FCA)*DAX
             else
               FCA = zero
               DCA = zero
               RCA = zero
             endif

             FAMA = FCXA + FCXB + FAX - CFXX(LM1)
             if( FAMA > zero ) then
               RAMA = (RCXA*FCXA*CLOLDAMPCT + RCXB*FCXB*CLNEWAMPCT + RAX*FAX*AMPCT)/FAMA
	       if( RAMA > zero ) then
                 DAMA = (RCXA*FCXA*CLOLDAMPCT)/(RAMA*FAMA)*DCXA +  &
                        (RCXB*FCXB*CLNEWAMPCT)/(RAMA*FAMA)*DCXB +  &
                        (RAX*FAX*AMPCT)/(RAMA*FAMA)*DAX
	       else
		  FAMA = zero
                  DAMA = zero
	       endif
             else
               FAMA = zero
               DAMA = zero
               RAMA = zero
             endif
           else upper_level
             AMPCT      = zero
             AMCLPCT    = zero
             CLNEWPCT   = zero
             CLNEWAMPCT = zero
             CLOLDPCT   = zero
             CLOLDAMPCT = zero
           endif upper_level
         else has_rls
	   RNEW = zero
           QTEVAPCXA = QTTOPCA
           QTEVAPAX = QTTOPAA
           if( L > 1 ) then
             if( RLS(LM1) > zero ) then
               CFXX(LM1) = max( CFMIN,CFR(LM1) )
!              if( CFR(LM1) >= CFMIN ) then
!                CFXX(LM1) = CFR(LM1)
!              else
!                CFXX(LM1) = CFMIN
!              endif
             else
               CFXX(LM1) = CFR(LM1)
             endif
           endif
           AMPCT      = zero
           AMCLPCT    = zero
           CLNEWPCT   = zero
           CLNEWAMPCT = zero
           CLOLDPCT   = zero
           CLOLDAMPCT = zero
           RCA        = zero
           RAMA       = zero
           FCA        = zero
           FAMA       = zero
           DCA        = zero
           DAMA       = zero
         endif has_rls

         if( do_diag .and. is_hno3 ) then
           fama_wrk(l) = fama
           rama_wrk(l) = rama
         endif
!-----------------------------------------------------------------------
!  Net loss can not exceed QTT in each region
!-----------------------------------------------------------------------
         QTNETLCXA = QTRAINCXA + QTRIMECXA + QTWASHCXA - QTEVAPCXA
         QTNETLCXA = min( QTT(L)*FCXA,QTNETLCXA )
   
         QTNETLCXB =QTRAINCXB
         QTNETLCXB = min( QTT(L)*FCXB,QTNETLCXB )

         QTNETLAX = QTWASHAX - QTEVAPAX
         QTNETLAX = min( QTT(L)*FAX,QTNETLAX )
              
         QTTNEW(L) = QTT(L) - (QTNETLCXA + QTNETLCXB + QTNETLAX)

         if( do_diag .and. is_hno3 ) then
           qt_rain(l) = qtraincxa + qtraincxb
           qt_rime(l) = qtrimecxa
           qt_wash(l) = qtwashcxa + qtwashax
           qt_evap(l) = qtevapcxa + qtevapax
           frc(l,1) = qtnetlcxa
           frc(l,2) = qtnetlcxb
           frc(l,3) = qtnetlax
         endif
         if( debug .and. is_hno3 .and. l == kdiag ) then
           write(*,*) ' '
           write(*,*) 'washout: qtraincxa, qtraincxb, qtrimecxa @ level = ',l
           write(*,'(1p,3g15.7)') qtraincxa, qtraincxb, qtrimecxa
           write(*,*) ' '
         endif
         if ( debug ) then
           if( (l == 3 .or. l == 2) ) then
             write(*,*) 'washout: hno3, hno3, qtnetlca,b, qtnetlax @ level = ',l
             write(*,'(1p,5g15.7)') qttnew(l), qtt(l), qtnetlcxa, qtnetlcxb, qtnetlax
             write(*,*) 'washout: qtwashax, qtevapax,fax,fama'
             write(*,'(1p,5g15.7)') qtwashax, qtevapax, fax, fama
           endif
         endif

         QTTOPCAX = (QTTOPCA + QTNETLCXA)*CLOLDPCT + QTNETLCXB*CLNEWPCT + (QTTOPAA + QTNETLAX)*AMCLPCT
         QTTOPAAX = (QTTOPCA + QTNETLCXA)*CLOLDAMPCT + QTNETLCXB*CLNEWAMPCT + (QTTOPAA + QTNETLAX)*AMPCT
         QTTOPCA = QTTOPCAX
         QTTOPAA = QTTOPAAX
       end do level_loop

       if ( debug ) then
         if( is_hno3 ) then
           write(*,*) ' '
           write(*,*) 'washout: clwx_wrk'
           write(*,'(1p,5g15.7)') clwx_wrk(1:le)
           write(*,*) 'washout: cfr'
           write(*,'(1p,5g15.7)') cfr(1:le)
           write(*,*) 'washout: cfxx'
           write(*,'(1p,5g15.7)') cfxx(1:le)
           write(*,*) 'washout: cf trigger'
           write(*,'(10l4)') cf_trigger(1:le)
           write(*,*) 'washout: evaprate'
           write(*,'(1p,5g15.7)') evaprate(1:le)
           write(*,*) 'washout: rls'
           write(*,'(1p,5g15.7)') rls(1:le)
           write(*,*) 'washout: rls/garea'
           write(*,'(1p,5g15.7)') rls_wrk(1:le)
           write(*,*) 'washout: rnew_wrk'
           write(*,'(1p,5g15.7)') rnew_wrk(1:le)
           write(*,*) 'washout: rnew_flag'
           write(*,'(10l4)') rnew_flag(1:le)
           write(*,*) 'washout: deltarime_wrk'
           write(*,'(1p,5g15.7)') deltarime_wrk(1:le)
           write(*,*) 'washout: rama_wrk'
           write(*,'(1p,5g15.7)') rama_wrk(1:le)
           write(*,*) 'washout: fama_wrk'
           write(*,'(1p,5g15.7)') fama_wrk(1:le)
           write(*,*) 'washout: rca_wrk'
           write(*,'(1p,5g15.7)') rca_wrk(1:le)
           write(*,*) 'washout: fca_wrk'
           write(*,'(1p,5g15.7)') fca_wrk(1:le)
           write(*,*) 'washout: rcxa_wrk'
           write(*,'(1p,5g15.7)') rcxa_wrk(1:le)
           write(*,*) 'washout: fcxa_wrk'
           write(*,'(1p,5g15.7)') fcxa_wrk(1:le)
           write(*,*) 'washout: rcxb_wrk'
           write(*,'(1p,5g15.7)') rcxb_wrk(1:le)
           write(*,*) 'washout: fcxb_wrk'
           write(*,'(1p,5g15.7)') fcxb_wrk(1:le)
           write(*,*) 'washout: rax1_wrk'
           write(*,'(1p,5g15.7)') rax_wrk(1:le,1)
           write(*,*) 'washout: fax1_wrk'
           write(*,'(1p,5g15.7)') fax_wrk(1:le,1)
           write(*,*) 'washout: rax2_wrk'
           write(*,'(1p,5g15.7)') rax_wrk(1:le,2)
           write(*,*) 'washout: fax2_wrk'
           write(*,'(1p,5g15.7)') fax_wrk(1:le,2)
           write(*,*) 'washout: rls_flag'
           write(*,'(1p,10l4)') rls_flag(1:le)
           write(*,*) 'washout: freezing'
           write(*,'(1p,10l4)') freezing(1:le)
           write(*,*) 'washout: qtnetlcxa'
           write(*,'(1p,5g15.7)') frc(1:le,1)
           write(*,*) 'washout: qtnetlcxb'
           write(*,'(1p,5g15.7)') frc(1:le,2)
           write(*,*) 'washout: qtnetlax'
           write(*,'(1p,5g15.7)') frc(1:le,3)
           write(*,*) ' '
         endif
       endif
!-----------------------------------------------------------------------
!  reload new tracer mass and rescale moments: check upper limits (LE) 
!-----------------------------------------------------------------------
       QTTJFL(:le,N) = QTTNEW(:le)

     end do species_loop
!
     return
   end subroutine washo
!---------------------------------------------------------------------
      subroutine DISGAS (CLWX,CFX,MOLMASS,HSTAR,TM,PR,QM,QT,QTDIS)
!---------------------------------------------------------------------
      implicit none
      real(r8), intent(in) :: CLWX,CFX    !cloud water,cloud fraction 
      real(r8), intent(in) :: MOLMASS     !molecular mass of tracer
      real(r8), intent(in) :: HSTAR       !Henry's Law coeffs A*exp(-B/T)
      real(r8), intent(in) :: TM          !temperature of box (K)
      real(r8), intent(in) :: PR          !pressure of box (hPa)
      real(r8), intent(in) :: QM          !air mass in box (kg)
      real(r8), intent(in) :: QT          !tracer in box (kg)
      real(r8), intent(out) :: QTDIS      !tracer dissolved in aqueous phase 
 
      real(r8)  MUEMP
      real(r8), parameter :: INV298 = 1._r8/298._r8
      real(r8), parameter  :: TICE=273._r8
      real(r8), parameter  :: TMIX=258._r8
      real(r8), parameter  :: RETEFF=0.5_r8
!---Next calculate rate of uptake of tracer

!---effective Henry's Law constant: H* = moles-T / liter-precip / press(atm-T)
!---p(atm of tracer-T) = (QT/QM) * (.029/MolWt-T) * pressr(hPa)/1000
!---limit temperature effects to T above freezing
!----MU from fit to Kaercher and Voigt (2006)

      if(TM .ge. TICE) then
         QTDIS=(HSTAR*(QT/(QM*CFX))*0.029_r8*(PR/1.0e3_r8))*(CLWX*QM)
      elseif (TM .le. TMIX) then
         MUEMP=exp(-14.2252_r8+(1.55704e-1_r8*TM)-(7.1929e-4_r8*(TM**2.0_r8)))
         QTDIS=MUEMP*(MOLMASS/18._r8)*(CLWX*QM)
      else
       QTDIS=RETEFF*((HSTAR*(QT/(QM*CFX))*0.029_r8*(PR/1.0e3_r8))*(CLWX*QM))
      endif

      return
      end subroutine DISGAS

!-----------------------------------------------------------------------
      subroutine RAINGAS (RRAIN,DTSCAV,CLWX,CFX,QM,QT,QTDIS,QTRAIN)
!-----------------------------------------------------------------------
!---New trace-gas rainout from large-scale precip with two time scales,
!---one based on precip formation from cloud water and one based on 
!---Henry's Law solubility: correct limit for delta-t
!---    
!---NB this code does not consider the aqueous dissociation (eg, C-q) 
!---   that makes uptake of HNO3 and H2SO4 so complete.  To do so would
!---   require that we keep track of the pH of the falling rain.
!---THUS the Henry's Law coefficient KHA needs to be enhanced to incldue this!
!---ALSO the possible formation of other soluble species from, eg, CH2O, H2O2
!---   can be considered with enhanced values of KHA.
!---
!---Does NOT now use RMC (moist conv rain) but could, assuming 30% coverage
!-----------------------------------------------------------------------
      implicit none
      real(r8), intent(in) :: RRAIN       !new rain formation in box (kg/s)
      real(r8), intent(in) :: DTSCAV      !time step (s)
      real(r8), intent(in) :: CLWX,CFX !cloud water and cloud fraction
      real(r8), intent(in) :: QM          !air mass in box (kg)
      real(r8), intent(in) :: QT          !tracer in box (kg) 
      real(r8), intent(in) :: QTDIS          !tracer in aqueous phase (kg) 
      real(r8), intent(out) :: QTRAIN      !tracer picked up by new rain  

      real(r8)   QTLF,QTDISSTAR





      QTDISSTAR=(QTDIS*(QT*CFX))/(QTDIS+(QT*CFX))
 
!---Tracer Loss frequency (1/s) within cloud fraction:
      QTLF = (RRAIN*QTDISSTAR)/(CLWX*QM*QT*CFX)

!---in time = DTSCAV, the amount of QTT scavenged is calculated 
!---from CF*AMOUNT OF UPTAKE 
      QTRAIN = QT*CFX*(1._r8 - exp(-DTSCAV*QTLF))

      return
      end subroutine RAINGAS


!-----------------------------------------------------------------------
      subroutine WASHGAS (RWASH,BOXF,DTSCAV,QTRTOP,HSTAR,TM,PR,QM, &
                            QT,QTWASH,QTEVAP)
!-----------------------------------------------------------------------
!---for most gases below-cloud washout assume Henry-Law equilib with precip
!---assumes that precip is liquid, if frozen, do not call this sub
!---since solubility is moderate, fraction of box with rain does not matter
!---NB this code does not consider the aqueous dissociation (eg, C-q) 
!---   that makes uptake of HNO3 and H2SO4 so complete.  To do so would
!---   require that we keep track of the pH of the falling rain.
!---THUS the Henry's Law coefficient KHA needs to be enhanced to incldue this!
!---ALSO the possible formation of other soluble species from, eg, CH2O, H2O2
!---   can be considered with enhanced values of KHA.
!-----------------------------------------------------------------------
      implicit none
      real(r8), intent(in)  :: RWASH   ! precip leaving bottom of box (kg/s)
      real(r8), intent(in)  :: BOXF   ! fraction of box with washout
      real(r8), intent(in)  :: DTSCAV  ! time step (s)
      real(r8), intent(in)  :: QTRTOP  ! tracer-T in rain entering top of box 
!                                              over time step (kg)
      real(r8), intent(in)  :: HSTAR ! Henry's Law coeffs A*exp(-B/T)
      real(r8), intent(in)  :: TM      ! temperature of box (K)
      real(r8), intent(in)  :: PR      ! pressure of box (hPa)
      real(r8), intent(in)  :: QT      ! tracer in box (kg)
      real(r8), intent(in)  :: QM      ! air mass in box (kg)
      real(r8), intent(out) :: QTWASH  ! tracer picked up by precip (kg)
      real(r8), intent(out) :: QTEVAP  ! tracer evaporated from precip (kg)
      
      real(r8), parameter :: INV298 = 1._r8/298._r8
      real(r8)            :: FWASH, QTMAX, QTDIF

!---effective Henry's Law constant: H* = moles-T / liter-precip / press(atm-T)
!---p(atm of tracer-T) = (QT/QM) * (.029/MolWt-T) * pressr(hPa)/1000
!---limit temperature effects to T above freezing

!
! jfl
!
! added test for BOXF = 0.
!
      if ( BOXF == 0._r8 ) then
        QTWASH = 0._r8
        QTEVAP = 0._r8
        return
      end if

!---effective washout frequency (1/s):
        FWASH = (RWASH*HSTAR*29.e-6_r8*PR)/(QM*BOXF)
!---equilib amount of T (kg) in rain thru bottom of box over time step
        QTMAX = QT*FWASH*DTSCAV
      if (QTMAX .gt. QTRTOP) then
!---more of tracer T can go into rain
         QTDIF = min (QT, QTMAX-QTRTOP)
         QTWASH = QTDIF * (1._r8 - exp(-DTSCAV*FWASH))
         QTEVAP=0._r8
      else
!--too much of T in rain, must degas/evap T
         QTWASH = 0._r8
         QTEVAP = QTRTOP - QTMAX
      endif
     
      return
      end subroutine WASHGAS

!-----------------------------------------------------------------------
      function DEMPIRICAL (CWATER,RRATE)
!-----------------------------------------------------------------------

      implicit none
      real(r8), intent(in)  :: CWATER   
      real(r8), intent(in)  :: RRATE

      real(r8) :: DEMPIRICAL
 
      real(r8) RRATEX,WX,THETA,PHI,ETA,BETA,ALPHA,BEE
      real(r8) GAMTHETA,GAMBETA
  


      RRATEX=RRATE*3600._r8       !mm/hr
      WX=CWATER*1.0e3_r8  !g/m3

      if(RRATEX .gt. 0.04_r8) then
         THETA=exp(-1.43_r8*dlog10(7._r8*RRATEX))+2.8_r8
      else
         THETA=5._r8
      endif
      PHI=RRATEX/(3600._r8*10._r8) !cgs units
      ETA=exp((3.01_r8*THETA)-10.5_r8)
      BETA=THETA/(1._r8+0.638_r8)
      ALPHA=exp(4._r8*(BETA-3.5_r8))
      BEE=(.638_r8*THETA/(1._r8+.638_r8))-1.0_r8
      GAMTHETA = GAMMA(THETA)
      GAMBETA  = GAMMA(BETA+1._r8)
      DEMPIRICAL=(((WX*ETA*GAMTHETA)/(1.0e6_r8*ALPHA*PHI*GAMBETA))** &
                 (-1._r8/BEE))*10._r8      ! in mm (wx/1e6 for cgs)
      

      return
      end function DEMPIRICAL
function GAMMA( X )
!-----------------------------------------------------------------------
!       Purpose: Compute the gamma function â(x)
!       Input :  x  --- Argument of â(x)
!                       ( x is not equal to 0,-1,-2,úúú )
!       Output:  GA --- â(x)
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  dummy arguments
!-----------------------------------------------------------------------
        real(r8), intent(in)  :: X

!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------
        real(r8), parameter :: PI = 3.141592653589793e0_r8

        integer :: k, M, M1
        real(r8)    :: GR, R, Z
        real(r8)    :: G(26)

!-----------------------------------------------------------------------
!  function definition
!-----------------------------------------------------------------------
        real(r8) :: GAMMA

        DATA G/1.0e0_r8,0.5772156649015329_r8,                     &
            -0.6558780715202538e0_r8, -0.420026350340952e-1_r8,    &
            0.1665386113822915e0_r8,-.421977345555443e-1_r8,       &
            -.96219715278770e-2_r8, .72189432466630e-2_r8,         &
            -.11651675918591e-2_r8, -.2152416741149e-3_r8,         &
            .1280502823882e-3_r8, -.201348547807e-4_r8,            &
            -.12504934821e-5_r8, .11330272320e-5_r8,               &
            -.2056338417e-6_r8, .61160950e-8_r8,                   &
            .50020075e-8_r8, -.11812746e-8_r8,                     &
            .1043427e-9_r8, .77823e-11_r8,                         &
            -.36968e-11_r8, .51e-12_r8,                            &
            -.206e-13_r8, -.54e-14_r8, .14e-14_r8, .1e-15_r8/

is_integer : &
        IF( x == real( int(x) ) ) then
          IF( X > zero ) THEN
            GAMMA = ONE
            M1 = INT(X) - 1
            DO K = 2,M1
              GAMMA = GAMMA*real(K)
            END DO
          ELSE
            GAMMA = 1.0e36_r8
          ENDIF
        ELSE is_integer
          IF( ABS(X) > ONE ) THEN
            Z = ABS(X)
            M = INT(Z)
            R = ONE
            DO K = 1,M
              R = R*(Z - real(k))
            END DO
            Z = Z - real(M)
          ELSE
            Z = X
          ENDIF
          GR = G(26)
          DO K = 25,1,-1
            GR = GR*Z + G(K)
          end DO
          GAMMA = ONE/(GR*Z)
          IF( ABS(X) > ONE ) THEN
            GAMMA = GAMMA*R
            IF( X < zero ) then
              GAMMA = -PI/(X*GAMMA*SIN( PI*X ))
            ENDIF
          ENDIF
        ENDIF is_integer

END function GAMMA

!
end module mo_neu_wetdep
