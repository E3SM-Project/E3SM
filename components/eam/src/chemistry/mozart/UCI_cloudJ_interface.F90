module UCI_cloudJ_interface
  !>>>>>>>>  Cloud-J version 7.7   (Adapted to E3SM by Philip Cameron-Smith, 2020)
  !
  ! It just uses read in SZA, and LatxMONTH for a climatology of T & O3
  !       CLDFLAG = 1  :  Clear sky J's
  !       CLDFLAG = 2  :  Averaged cloud cover
  !       CLDFLAG = 3  :  cloud-fract**3/2, then average cloud cover
  !       CLDFLAG = 4  :  ****not used
  !       CLDFLAG = 5  :  Random select NRANDO ICA's from all(Independent Column Atmos.)
  !       CLDFLAG = 6  :  Use all (up to 4) quadrature cloud cover QCAs (mid-pts of bin)
  !       CLDFLAG = 7  :  Use all (up to 4) QCAs (average clouds within each Q-bin)
  !       CLDFLAG = 8  :  Calculate J's for ALL ICAs (up to 20,000 per cell!)

  ! NOTE: Fast-JX works 'upsidedown', ie the Fast-JX code indexes
  ! grid-levels from the ground up, while E3SM indexes grid-levels
  ! from the top of the atmosphere downwards.  Hence, this interface
  ! routine flips the E3SM levels before putting them into the Fast-JX
  ! variables, and then flips the Fast-JX result back before returning
  ! it to E3SM.

  
  use ppgrid,        only : pcols, pver, pverp, begchunk, endchunk
  use chem_mods,     only : nabscol, phtcnt, gas_pcnst, nfs
  use dust_model,    only : ndst => dust_nbin
  use mo_constants,  only : pi,r2d,boltz,d2r
  use shr_kind_mod,  only : r8 => shr_kind_r8
  use mo_chem_utls,  only : get_spc_ndx, get_rxt_ndx, get_inv_ndx
!  use physconst,     only : mwdry, mwch4, mwn2o, mwf11, mwf12, mwh2o, mwo3  ! pjc needed?
  use chem_surfvals, only : chem_surfvals_get
  use mo_tracname,   only : solsym
  
  use cam_history,      only : fieldname_len
  use cam_logfile,      only : iulog
  use cam_abortutils,   only : endrun
  use spmd_utils,       only : iam, masterproc ! PJC for debugging only.
  use shr_sys_mod,      only : shr_sys_flush ! PJC for debugging only.
  use thread_mod,       only : omp_get_thread_num  ! PJC for debugging only.
  
!!! I/O utilities.  Used by Xactive_Photo.   Use for FastJ?
! use pio
! use cam_pio_utils,only : cam_pio_openfile

  USE FJX_CMN_MOD
  
  implicit none
  private


!!!  Variables common to this module
  logical, parameter          :: LPRTJ  = .FALSE. !Turn On/off diagnostic output to atm.log and fort.7
  logical, parameter          :: LPRTJ7 = .FALSE. !Turn On/off diagnostic to fort.7 (which only makes sense if run on a single CPU).
  character*6,  dimension(JVN_)  ::  TITLJXX
  integer                     :: JVNU,ANU,L1U
  integer                     :: NJXX


  integer ::  ox_ndx, o3_ndx, o3_inv_ndx, o3rad_ndx    ! O3 related indicies into concentration arrays
  integer ::  ch4_ndx, inv_ndx_ch4          ! CH4 index into concentration arrays
  integer ::  inv_ndx_M, inv_ndx_h2o        ! index in invariants array of airmass
  real(r8)::  CH4_SURF                      ! Surface concentration of CH4 (vmr)

  ! flag for the source used  (1=advection_array, 2=invariants_array, 3=surface_concentration)
  integer ::  CH4_source_flag      ! CH4 

  integer ::  E3SM_FastJ_aerosol_index(2,gas_pcnst) ! 1st element = E3SM index; 2nd element = Fast-J aerosol index
  integer ::  num_E3SM_FastJ_aerosols               ! # of aerosols in E3SM_FastJ_aerosol_index
  
!---------------------------------------------------------------------------------
! Public interfaces
!---------------------------------------------------------------------------------
  public :: cloudj_init, cloudJ_interface
 
!================================================================================================
contains

!---------------------------------------------------------------------------------
! cloudJ_init routine
!---------------------------------------------------------------------------------

subroutine cloudJ_init()  !pmid, pint, zmid, zint, rlats, rlons, ncol)

  USE FJX_INIT_MOD

  integer :: tracer_num,kk     ! local variables for mapping from E3SM tracers to Fast-J aerosols
  character(len=16) :: aerosol_string
!  character(len=*)  :: error_string
  
!  real(r8), intent(in)    :: pmid(pcols,pver)             ! midpoint pressure (Pa)
!  real(r8), intent(in)    :: pint(pcols,pver+1)           ! interface pressures (Pa)
!  real(r8), intent(in)    :: zmid(ncol,pver)              ! midpoint height (m)
!  real(r8), intent(in)    :: zint(ncol,pver+1)            ! interface geopotential in km
!  real(r8), intent(in)    :: rlats(ncol), rlons(ncol)     ! lats and longs (radians)

  
  ANU = AN_     ! max # FJX aerosols in layer
  JVNU = JVN_   ! max # of J-values.  Currently set in UCI_fjx_cmn_mod.f90.  Can it be set by E3SM?
  L1U = L1_     ! model levels +1

!!! Fast-J sets the number of atmospheric levels via the LPAR
!!! parameter in the FJX_CMN_MOD module. This must be consistent with
!!! E3SM.  The easiest solution is to just abort if the # of levels
!!! are inconsistent.   [PJC]
  IF ( pver .NE. LPAR ) THEN
     write (iulog,*) 'cloudj_interface: #E3SM_atm_levels = pver =',pver,',  LPAR =',LPAR
     write (iulog,*) 'cloudj_interface: LPAR may need to be modified to match model levels.'
     call endrun('cloudj_interface: ERROR. # of model layers in E3SM and Cloud-J are different.')
  ENDIF

  IF ( pver .NE. LWEPAR ) THEN
     write (iulog,*) 'cloudj_interface: #E3SM_atm_levels = pver =',pver,',  LWEPAR =',LWEPAR
     write (iulog,*) 'cloudj_interface: LWEPAR may need to be modified to match model levels.'
     write (iulog,*) 'cloudj_interface: NOTE: The interface routine could be modified to use LWEPAR, but need to be careful about indexing.'
     call endrun('cloudj_interface: ERROR. # of cloud layers in E3SM and Cloud-J are different.')
  ENDIF

!!! Determine tracer indicies within concentration arrays
    o3_ndx     = get_spc_ndx( 'O3' )
    if (LPRTJ) write(iulog,'(a,i3)')'FAST-JX: o3_ndx =',o3_ndx

    ch4_ndx     = get_spc_ndx( 'CH4' )
    CH4_source_flag = 1
    if (LPRTJ) write(iulog,'(a,i3)')'FAST-JX: ch4_ndx =',ch4_ndx
    if (ch4_ndx < 0) then
       CH4_source_flag = 2
       if (LPRTJ) write(iulog,'(a)')'Prognostic tracer for CH4 is not available.  Trying invariant field...'
       inv_ndx_ch4 = get_inv_ndx( 'CH4' )
       inv_ndx_M   = get_inv_ndx( 'M' )
       if (LPRTJ) write(iulog,'(a,2i3)')'FAST-JX: inv_ndx_ch4, inv_ndx_M =',inv_ndx_ch4,inv_ndx_M
       if ( (inv_ndx_ch4 < 0)  .OR. (inv_ndx_M < 0) ) then
       CH4_source_flag = 3
          if (LPRTJ) write(iulog,'(a)')'Invariant fields needed for CH4 are not available.  Trying surface concentration...'
          CH4_surf = chem_surfvals_get('CH4VMR')
          if (LPRTJ) write(iulog,'(a,E12.3)') 'FAST-JX: Set CH4 to surface value.  CH4_surf(vmr) =',CH4_surf
       endif
    endif
    if (LPRTJ) write(iulog,'(a,i3)')'FAST-JX: CH4_source_flag =',CH4_source_flag

    inv_ndx_h2o = get_inv_ndx( 'H2O' )

    !-----------------------------------------------------------------------
    !---read in & store all fast-JX data:   single call at set up
    !-----------------------------------------------------------------------
    
    !!!  JVNU passed in to INIT_FJX
    !!!  NJXX passed back from INIT_FJX. Number of species to calculate J-values for.
    !!!  TITLJXX passed back from INIT_FJX.  

!!!    !$omp master
    call INIT_FJX (TITLJXX,JVNU,NJXX)
!!!    !$omp end master          
!!!    !$omp barrier                                                                                                                               
    !  if (LPRTJ) write(iulog,'(a)')     'FAST-JX: NOTE: number of E3SM photolysis reactions will usually be similar to the number of entries in the j2j data file.' 
    if (LPRTJ) write(iulog,'(a,2i5)') 'FAST-JX: phtcnt, NRATJ =',phtcnt, NRATJ

    IF ( phtcnt .NE. NRATJ ) then
     !NOTE: The UCI-CTM may be able to handle cases in which the j2j file doesn't match the mechanism exactly. However, for E3SM, the current requirement is that the j2j file much match the chemical mechanism. 
       call endrun('The number of photolysis reaction specified in j2j file must be the same as thenumber of photolysis reactions in E3SM chemical mechanism.')  ! The lack of spaces before 'number' is needed to line-break properly.
    endif

    !     open (77,file='tables/atmos_PTClds.dat',status='old',SHARED) !UNIX

    !! Set up mapping from E3SM tracer number to Fast-J aerosol index
    !! The Fast-J aerosol properties are not based on the E3SM
    !! aerosols, so this mapping is just a place holder.
    !! NOTE: this is based on code from chm_diags_inti_ac

    !! This code assumes we are using MAM.
    E3SM_FastJ_aerosol_index(:,:) = 0
    num_E3SM_FastJ_aerosols = 0
    do tracer_num = 1,gas_pcnst
       kk=0
       if ( kk==0 ) kk = index(trim(solsym(tracer_num)), '_a')    ! indicates this is an aerosol
!?!       if ( kk==0 ) kk = index(trim(solsym(tracer_num)), '_c') ! indicates this is an aerosol in cloud phase (I think)
       if ( kk>0 ) then ! must be aerosol species
          num_E3SM_FastJ_aerosols = num_E3SM_FastJ_aerosols + 1
          E3SM_FastJ_aerosol_index(1,num_E3SM_FastJ_aerosols) = tracer_num
          aerosol_string = solsym(tracer_num)
          aerosol_string = aerosol_string(1:kk-1)         ! isolate the name of the aerosol
!          write(iulog,'(A)') 'PJC Fast-J: '//trim(aerosol_string)
          select case (trim(aerosol_string))
          case ('so4','pom','soa','ncl','mom')
             E3SM_FastJ_aerosol_index(2,num_E3SM_FastJ_aerosols) =  5    ! 05 = UT-sulfM  in FJX_scat-aer.dat
          case ('dst')
             E3SM_FastJ_aerosol_index(2,num_E3SM_FastJ_aerosols) = 14    ! 14 = HDust.80  in FJX_scat-aer.dat
          case ('bc')
             E3SM_FastJ_aerosol_index(2,num_E3SM_FastJ_aerosols) =  6    ! 06 = UM-BC1    in FJX_scat-aer.dat
          case ('num')
             ! number density is not needed by Fast-J, so ignore for now.
          case default
!             error_string = 'FAST-J: Unknown E3SM aerosol type = '//trim(aerosol_string)
             call endrun('FAST-J: Unknown E3SM aerosol type = '//trim(aerosol_string))
          end select
       endif
    enddo
    IF (num_E3SM_FastJ_aerosols .GT. ANU) THEN
       call endrun('FAST-J: ERROR, num_E3SM_FastJ_aerosols > ANU')
    ENDIF
    
    
end subroutine cloudJ_init

  
!---------------------------------------------------------------------------------
! cloudJ_interface routine
!---------------------------------------------------------------------------------

subroutine cloudJ_interface(photos, vmr, mmr, invariants, temper, cldfr, cldwat, cldice, &
                 pmid, pint, zmid, zint, rlats, rlons, col_dens, zen_angle, srf_alb, &
                 tdens, ps, ts, esfact, relhum, dust_vmr, &
 !                           dt_diag, fracday, &
                 ncol, lchnk, do_clouds, do_aerosols)
     
  USE FJX_SUB_MOD
  USE CLD_SUB_MOD, ONLY : CLOUD_JX
  USE OSA_SUB_MOD
  
     !----------- CloudJ_interface arguments -------------------

  real(r8), intent(out)   :: photos(ncol,pver,phtcnt)     ! photodissociation rates (1/s)
  integer,  intent(in)    :: ncol, lchnk
  real(r8), intent(in)    :: esfact                       ! earth sun distance factor (currently not used by Fast-J)
  real(r8), intent(in)    :: ps(pcols)                    ! surface pressure (Pa)
  real(r8), intent(in)    :: ts(ncol)                     ! surface temperature (K)
  real(r8), intent(in)    :: col_dens(ncol,pver,nabscol) ! column densities (molecules/cm^2)
  real(r8), intent(in)    :: zen_angle(ncol)              ! solar zenith angle (radians)
  real(r8), intent(in)    :: srf_alb(pcols)               ! surface albedo
  real(r8), intent(in)    :: tdens(ncol,pver)             ! total atms density (molecules/cm^3)
  real(r8), intent(in)    :: vmr(ncol,pver,gas_pcnst)     ! species concentration (mol/mol)
  real(r8), intent(in)    :: mmr(pcols,pver,gas_pcnst)    ! species concentration (kg/kg)
  real(r8), intent(in)    :: invariants(ncol,pver,nfs)    ! Fixed species, including N2 and species from files (cgs density)
  real(r8), intent(in)    :: pmid(pcols,pver)             ! midpoint pressure (Pa)
  real(r8), intent(in)    :: pint(pcols,pver+1)           ! interface pressures (Pa)
  real(r8), intent(in)    :: zmid(ncol,pver)              ! midpoint height (m)
  real(r8), intent(in)    :: zint(ncol,pver+1)            ! interface geopotential in km
  real(r8), intent(in)    :: rlats(ncol), rlons(ncol)     ! lats and longs (radians)
  real(r8), intent(in)    :: temper(pcols,pver)           ! midpoint temperature (K)
  real(r8), intent(in)    :: relhum(ncol,pver)            ! relative humidity
! real(r8), intent(in)    :: cwat(ncol,pver)              ! cloud water (kg/kg)
  real(r8), intent(in)    :: cldfr(ncol,pver)             ! cloud fraction
  real(r8), intent(in)    :: cldwat(ncol,pver)            ! cloud water (kg/kg)
  real(r8), intent(in)    :: cldice(ncol,pver)            ! cloud ice (kg/kg)
  real(r8), intent(in)    :: dust_vmr(ncol,pver,ndst)     ! dust concentration (mol/mol)
  logical,  intent(in)    :: do_clouds                    ! Should Cloud-J use cloud optical depths?
  logical,  intent(in)    :: do_aerosols                  ! Should Cloud-J use aerosol optical depths?
!     real(r8), intent(out)   :: dt_diag(pcols,8)             ! od diagnostics
!     real(r8), intent(out)   :: fracday(pcols)               ! fraction of day
     
  !---------------key params in/out of CLOUD_J-------------------------
  logical                     :: LDARK
  integer                     :: IRAN
  integer                     :: NICA,JCOUNT
  integer                     :: JP
  real*8                      :: U0,SZA,SOLF
  real*8,  dimension(L2_  )   :: PPP,ZZZ
  real*8,  dimension(L1_  )   :: TTT,HHH,DDD,RRR,OOO,CCC
  real*8,  dimension(L1_  )   :: O3,CH4,H2O, OD18
  real*8,  dimension(S_+2,L1_):: SKPERD
  real*8,  dimension(6)       :: SWMSQ
  real*8,  dimension(L2_)     :: CLF,LWP,IWP,REFFL,REFFI
  integer, dimension(L2_)     :: CLDIW
  real*8,  dimension(L2_,AN_) :: AERSP
  integer, dimension(L2_,AN_) :: NDXAER
  real*8,  dimension(L_,JVN_) :: VALJXX
  real*8,  dimension(5,S_)    :: RFL
  real*8,  dimension(NQD_)    :: WTQCA
  
!-------------local use-----------------------
!  integer :: NSZA,MONTH,ILAT
! beware the OSA code uses single R*4 variables
  real*4  :: OWAVEL,OWIND,OCHLR,OSA_dir(5)
  real*8  :: YLAT,PSURF,ALBEDO(5),WIND,CHLR
  real*8, dimension(L2_) :: ETAA,ETAB, ZOFL,RI,TI,AER1,AER2,PPPX
  integer,dimension(L2_) :: NAA1,NAA2
  real*8, dimension(L_)  :: WLC,WIC
  real*8, dimension(LWEPAR) :: CLDFRW,CLDIWCW,CLDLWCW
  real*8  SCALEH,CF,PMID_UCI,PDEL,ZDEL,ICWC,F1
  integer I,J,K,L,N
  integer LTOP,JP04,JP09
  character*11, dimension(4)     ::  TITJX
  real*8 VJOSA(L2_,2),VJSTD(L2_,2)
  
  integer :: column_number
  integer :: level_number
  integer, parameter :: cloud_freezing_temperature = 273 !PJC: temperature below which cloud is assumed to be ice
  integer :: aerosol_num    ! local variable for mapping from E3SM to FAST-J aerosols

  photos(:,:,:) = 0.   ! This makes sure that values on the dark side (which Fast-JX doesn't calculate) are zero. 

!  write(iulog,'(A,2i3)') '*** CloudJ_interface, iam, thread_num = ',iam,OMP_get_thread_num()   !pjc

!if ( .not.masterproc ) return     !pjc: limit to one task, to make diagnostic output manageable.

!  if ( iam .GT. 2) return   ! pjc: limit to a few tasks.
  
!!!$OMP MASTER      ! to activate this line, only use 1 exclamation mark
!   call sleep(1*iam)
  if (LPRTJ) then
     write(iulog,'(A,i3,A,i3,A,i6)') '*** CloudJ_interface, iam =',iam,'   ncol =',ncol,'   lchnk =',lchnk    !pjc
     call shr_sys_flush(iulog)   ! pjc: flush buffer to iulog
  endif
     


!-----------------------------------------------------------------------

  DO column_number=1,ncol    !!!! Loop over columns in chunk
     if (LPRTJ) then
        WRITE (iulog,'(A,I3)') '*** PJC, FAST-JX: column_number =',column_number
     endif
     
!     MONTH = 8       ! month
!     YLAT  = 20.     ! Latitude
!     PSURF = 1013.   ! Surface Pressure

     ALBEDO(1) = srf_alb(column_number)   ! Default albedo for quad angle = 86 deg zenith angle
     ALBEDO(2) = srf_alb(column_number)   ! Default albedo for quad angle = 71 deg zenith angle
     ALBEDO(3) = srf_alb(column_number)   ! Default albedo for quad angle = 48 deg zenith angle
     ALBEDO(4) = srf_alb(column_number)   ! Default albedo for quad angle = 22 deg zenith angle
     ALBEDO(5) = srf_alb(column_number)   ! Default albedo for SZA incident ray

     WIND = 6.0     ! Surface wind in m/s (at what altitude?)
     CHLR = 0.08    ! Chlorophyl mg/m3
     
!     if (LPRTJ) write(iulog,'(a,2i5,5x,a,i5)') 'Atmosphere:',LPAR,LWEPAR, 'LPAR / LWEPAR', L1_
!     if (LPRTJ) write(iulog,'(a,f10.4)') 'P surface', PSURF
!     if (LPRTJ) write(iulog,'(a,3i4)') 'MONTH/ LAT',MONTH,NINT(YLAT)
     if (LPRTJ) then
        write(iulog,'(a,2F10.3)') 'lat, lon (degrees) =',rlats(column_number)/CPI180,rlons(column_number)/CPI180
        write(iulog,'(a,5f8.4)') 'Albedos 1:4 & 5=SZA', ALBEDO
        write(iulog,'(a,2f8.3)') 'OSA: wind & chlor-a',WIND,CHLR
     endif
!     read (77,*)     ! Read 7 header lines that aren't used anymore.
!     read (77,*)
!     read (77,*)
!     read (77,*)
!     read (77,*)
!     read (77,*)
!     read (77,*)
!     do L = 1,LPAR+1      ! Level interfaces
!        read (77,'(i3,1x,2f11.7,2x,f5.1,f5.2,f11.2,2(f7.3,i4))') &
!             J,ETAA(L),ETAB(L),TI(L),RI(L),ZOFL(L) &     
!             ,AER1(L),NAA1(L),AER2(L),NAA2(L)
!        ! eta-A   eta-B   Temp(K)   RH(fraction)  equivalent_altitude(m)   AER-P NDA  AER-P NDA
!     enddo

     TTT(1:pver) = temper(column_number,pver:1:-1)    ! Layer mid-point temperature (K)
     RRR(1:pver) = relhum(column_number,pver:1:-1)    ! Mid-point Relative Humidity (fraction)

!     AER1(:) = 0.     !!! PJC: temporarily set to zero   [Now set to MAM4 aerosols below]
!     NAA1(:) = 0.
!     AER2(:) = 0.
!     NAA2(:) = 0.

!     read (77,*)
!     do L = LWEPAR,1,-1   ! Clouds in levels   
!        read (77,'(i3,1p,e14.5,28x,2e14.5)') &
!             J,CLDFRW(L),CLDLWCW(L),CLDIWCW(L)
!        ! Cloud_fraction?    WLC(g/g)    WIC(g/g)
!     enddo
!     close(77)

     if ( do_clouds ) then
        CLDFRW(:)       = cldfr (column_number,pver:1:-1)   ! Cloud fraction (fraction) [[?]]
        CLDLWCW(1:pver) = cldwat(column_number,pver:1:-1)   ! Cloud liquid water (g/g)
        CLDIWCW(1:pver) = cldice(column_number,pver:1:-1)   ! Cloud ice water (g/g)
     else
        CLDFRW (:) = 0.
        CLDLWCW(:) = 0.
        CLDIWCW(:) = 0.
     endif
        
     PPP(1:pver+1) = pint(column_number,pver+1:1:-1) / 100D0      ! Interface pressures (mbar)
     PPP(pver+2)   = 0.
!     do L = 1,L1_
!        PPP(L) = ETAA(L) + ETAB(L)*PSURF    ! Pressure on level interfaces (mbar)
!     enddo
!!!! just for print out and levels          ! PJC: Not used in this code
!!!     do L = 1,L1_
!!!        PPPX(L) = 0.5d0*(PPP(L)+PPP(L+1))   ! Mid-layer pressure (mbar)
!!!     enddo
!---sets climatologies for O3, T, D & Z
!-----------------------------------------------------------------------
!     call ACLIM_FJX (YLAT,MONTH,PPP, TTT,O3,CH4, L1_)
        ! Returns Temperature(K), O3(ppm), CH4(ppb) (in layers plus above model top?)
!-----------------------------------------------------------------------
!     do L = 1,L_
!!!!       TTT(L) = TI(L)  keep climatology T's and O3's
!        RRR(L) = RI(L)
!     enddo

!     ZZZ(1)  = 16.d5*log10(1013.25d0/PPP(1))        ! zzz in cm
     ZZZ(1:pver+1) = zint(column_number,pver+1:1:-1) *1D5            ! Altitude in cm

     O3(:pver)  = vmr(column_number,pver:1:-1,o3_ndx) * 1D6     ! convert from vmr to ppm
     O3(pver+1) = O3(pver) 

     select case (CH4_source_flag)
     CASE(1)     ! CH4 from advected/prognostic array
        CH4(:pver) = vmr(column_number,pver:1:-1,ch4_ndx) * 1D9  !convert from vmr to ppb
     CASE(2)     ! CH4 from invariant/fixed array
        CH4(:pver) = invariants(column_number,pver:1:-1,inv_ndx_CH4) /   &
                       invariants(column_number,pver:1:-1,inv_ndx_M)  * 1D9  ! convert to ppb
     CASE(3)     ! CH4 from surface concentration
        CH4(:pver) = CH4_surf *1D9                               !convert to ppb        
     CASE default
        call endrun('FAST-JX: CH4_source_flag is not a recognized case.')
     end select
     CH4(pver+1) = CH4(pver)        !set concentration above model top be same as top layer

     do L = 1,L_
        DDD(L)  = (PPP(L)-PPP(L+1))*MASFAC
!!! geopotential since assumes g = constant
!        SCALEH      = 1.3806d-19*MASFAC*TTT(L)
!        ZZZ(L+1) = ZZZ(L) -( log(PPP(L+1)/PPP(L)) * SCALEH )
        OOO(L) = DDD(L)*O3(L)*1.d-6
        CCC(L) = DDD(L)*CH4(L)*1.d-9
     enddo
     L = L_+1
     ZZZ(L+1)= ZZZ(L) + ZZHT
     DDD(L)  = (PPP(L)-PPP(L+1))*MASFAC
     OOO(L) = DDD(L)*O3(L)*1.d-6
     CCC(L) = DDD(L)*CH4(L)*1.d-9
!-----------------------------------------------------------------------
!       call ACLIM_RH (PL, TL, QL, HHH, L1U)
!-----------------------------------------------------------------------
!     HHH(:) = 0.50d0    ! This is actually water vapor column density (molecules/cm2)
!!! set up clouds and aerosols
     HHH(:pver)  = DDD(:pver) * invariants(column_number,pver:1:-1,inv_ndx_H2O) /   &
                                invariants(column_number,pver:1:-1,inv_ndx_M)
     HHH(pver+1) = HHH(pver)
     
     AERSP(:,:)  = 0.d0
     NDXAER(:,:) = 0

     IF ( do_aerosols ) THEN
        do aerosol_num = 1, num_E3SM_FastJ_aerosols 
           do L = 1,L_
              NDXAER(L,aerosol_num)= E3SM_FastJ_aerosol_index(2,aerosol_num)
              AERSP(L,aerosol_num) = mmr(column_number,pver-L+1,E3SM_FastJ_aerosol_index(1,aerosol_num)) * (PPP(L)-PPP(L+1)) * 100d0 * 1000d0 / 9.8d0  !P(mbar) 100(Pa/mbar) *1000(g/kg) / g(m/s2) => g_air/m2
           enddo
        enddo
     ENDIF
  
     LTOP  = LWEPAR
     if (maxval(CLDFRW) .le. 0.005d0) then
        IWP(:) = 0.d0
        REFFI(:) = 0.d0
        LWP(:) = 0.d0
        REFFL(:) = 0.d0
     endif
     do L = 1,LTOP
        CLDIW(L) = 0
        CF  = CLDFRW(L)
        if (CF .gt. 0.005d0) then
           CLF(L) = CF
           WLC(L) = CLDLWCW(L) / CF
           WIC(L) = CLDIWCW(L) / CF
!  CLDIW is an integer flag: 1 = water cld, 2 = ice cloud, 3 = both
           if (WLC(L) .gt. 1.d-11) CLDIW(L) = 1
           if (WIC(L) .gt. 1.d-11) CLDIW(L) = CLDIW(L) + 2
        else
           CLF(L) = 0.d0
           WLC(L) = 0.d0
           WIC(L) = 0.d0
        endif
     enddo
!---derive R-effective for clouds:  the current UCI algorithm - use your own
     do L = 1,LTOP
!---ice clouds
        if (WIC(L) .gt. 1.d-12) then
           PDEL = PPP(L) - PPP(L+1)
           ZDEL = (ZZZ(L+1) - ZZZ(L))*0.01d0  ! m
           IWP(L) = 1000.d0*WIC(L)*PDEL*G100 /CLF(L)   ! g/m2
           ICWC =        IWP(L) / ZDEL          ! g/m3
           REFFI(L) = 164.d0 * (ICWC**0.23d0)
        else
           IWP(L) = 0.d0
           REFFI(L) = 0.d0
        endif
!---water clouds
        if (WLC(L) .gt. 1.d-12) then
           PMID_UCI = 0.5d0*(PPP(L)+PPP(L+1))
           PDEL = PPP(L) - PPP(L+1)
           F1   = 0.005d0 * (PMID_UCI - 610.d0)
           F1   = min(1.d0, max(0.d0, F1))
           LWP(L) = 1000.d0*WLC(L)*PDEL*G100     ! g/m2
           REFFL(L) = 9.6d0*F1 + 12.68d0*(1.d0-F1)
        else
           LWP(L) = 0.d0
           REFFL(L) = 0.d0
        endif
     enddo
     do L = 1,LTOP
        CLDFRW(L) = CLF(L)
     enddo
!!!  end of atmosphere setup


!!! begin call to Cloud_J

        SZA = zen_angle(column_number)/CPI180              ! solar zenith angle (degrees)
        
        if (LPRTJ) then
           write(iulog,'(a,i3)') 'CLDFLAG =',CLDFLAG
        endif
        
        do L = 1,LTOP
           CLF(L) = CLDFRW(L)
        enddo
        IRAN = 1
        SOLF = 1.d0
        U0 = cos(SZA*CPI180)       !  U0 = cos(sza)

! beware the OSA code uses single R*4 variables
        ANGLES(1) = sngl(EMU(1))
        ANGLES(2) = sngl(EMU(2))
        ANGLES(3) = sngl(EMU(3))
        ANGLES(4) = sngl(EMU(4))
        ANGLES(5) = sngl(U0)
        OWIND = sngl(WIND)
        OCHLR = sngl(CHLR)
        do K = 1,NS2
           OWAVEL = sngl(WL(K))
           call OSA(OWAVEL,OWIND,OCHLR, ANGLES,OSA_dir)
           do J = 1,5
              RFL(J,K) = dble(OSA_dir(J))
! this overwrite the OSA with the readin values above
              RFL(J,K) = ALBEDO(J)
           enddo
        enddo

        if (LPRTJ) then
           write(iulog,'(a,f8.3,3f8.5)')'SZA SOLF U0 albedo =' &
                ,SZA,SOLF,U0,RFL(5,18)
           call JP_ATM0(PPP,TTT,DDD,OOO,ZZZ, L_)
           write(iulog,'(A5,2A12)') 'level','O3(ppm)','CH4(ppb)'
           do L=pver+1,1,-1
              write(iulog,'(I5,2G12.3)') L, O3(L), CH4(L)
           enddo
           write(iulog,'(5X,A8,5A8)') 'wvl','u1','u2','u3','u4','u0'
           do K=1,NS2
              write(iulog,'(i5,f8.1,5f8.4)') K,WL(K), (RFL(J,K), J=1,5)
           enddo

           write(iulog,'(a)')''
        endif

        SKPERD(:,:)=0.d0
        SWMSQ(:)= 0.d0
        OD18(:) =0.d0
        WTQCA(:)= 0.d0
        
!=======================================================================
        call CLOUD_JX (U0,SZA,RFL,SOLF,LPRTJ,PPP,ZZZ,TTT,HHH,DDD,       &
               RRR,OOO,CCC,  LWP,IWP,REFFL,REFFI, CLF,CLDCOR,CLDIW,    &
               AERSP,NDXAER,L1U,ANU,JVNU, VALJXX,SKPERD,SWMSQ,OD18,    &
               CLDFLAG,NRANDO,IRAN,LNRG,NICA, JCOUNT,LDARK,WTQCA)
!=======================================================================

! summary results. This outputs to the file 'fort.7'. This may not
! work on all systems, but it works on Cori, and makes it easier to
! compare to the Fast-JX standalone output.
        if (LPRTJ7) then
           N=7
           write(N,'(a)')'===================== FAST-JX(v7.7) column ====================='
           write(N,'(a,i3)') 'column number in chunk =',column_number
           write(N,'(a,2F10.3)') ' lat, lon (degrees) =',rlats(column_number)/CPI180,rlons(column_number)/CPI180
           write(N,'(a,i5,F10.3)') ' CLDFLAG, SZA =', CLDFLAG,SZA
           write(N,*) 'LDARK =',LDARK
           write(N,*) 'WTQCA =',WTQCA
           write(N,'(a)') 'Average J-values. NOTE: levels are indexed with 1 being the ground!!!'
           write(N,'(1x,a,72(a6,3x))') 'Lev ',(TITLEJX(K), K=1,NJX)
           do L = L_,1,-1
              write(N,'(i3,1p, 72e9.2)') L,(VALJXX(L,K),K=1,NJX)
           enddo
           
           SKPERD(:,:)=NINT( SKPERD(:,:)*100 ) / 100. !Rounding to 2 decimal places for consistent output.
           
           write(N,'(a)') 'heating rate profiles in K/day v7.6  180-778nm '
           write(N,'(a4, 32f7.1)')'wvl ',(WL(K),K=NW1,NW2)
           do L = L_,1,-1
              write(N,'(i4,32f7.2)') L,(SKPERD(K,L), K=NW1,NW2)
           enddo
           write(N,'(a)') 'heating rate profiles in K/day v7.6 778-...nm plus 1:18 19:27 1:27'
           write(N,'(a4,32f7.1)')'wvl ',(WL(K),K=NW2+1,NS2)
           do L = L_,1,-1
              write(N,'(i4,35f7.2)') L,(SKPERD(K,L), K=NW2+1,NS2+2),SKPERD(S_+1,L)+SKPERD(S_+2,L)
           enddo
           
           write(N,'(a)') '---PHOTO_JX internal print: Solar fluxes (W/m2)--'
           write(N,'(a11,f12.4)')    ' inc TOTAL ',SWMSQ(1)
           write(N,'(a11,f12.4)')    ' rfl outtop',SWMSQ(2)
           write(N,'(a11,f12.4)')    ' abs in atm',SWMSQ(3)
           write(N,'(a11,f12.4)')    ' abs at srf',SWMSQ(4)
           write(N,'(a11,1p,e12.4)') ' PAR direct',SWMSQ(5)
           write(N,'(a11,1p,e12.4)') ' PAR diffus',SWMSQ(6)
           write(N,'(a)')''
        endif
!---map the J-values from fast-JX onto CTM (ZPJQUAD) using JIND & JFACTA
!--- from the 'FJX_j2j.dat' tables
!      do J = 1,NRATJ
!         JP = JIND(J)
!         if (JP .gt. 0) then
!            do L = 1,L_
!               ZPJQUAD(L,J) = ZPJQUAD(L,J) + VALJXX(L,JP)*JFACTA(J)
!          enddo
!        endif
!      enddo


!  real(r8), intent(inout) :: photos(ncol,pver,phtcnt)     ! photodissociation rates (1/s)

      IF (LPRTJ) WRITE(iulog,'(a,i4)')'NRATJ =',NRATJ  
      IF (LPRTJ) WRITE(iulog,'(a,200i3)')'JIND = ',JIND  
      do J = 1,NRATJ
         JP = JIND(J)
         IF (LPRTJ) write(iulog,'(a,3i4)')'J, JIND(J) =',J, JIND(J)
         if (JP .gt. 0) then
            IF (LPRTJ) write(iulog,'(a,F10.4)')'  JFACTA(J) =',JFACTA(J)
            photos(column_number,:,J) = VALJXX(L_:1:-1,JP)*JFACTA(J)           
         endif
      enddo
      IF (LPRTJ) then
         write(iulog,*)''
         write(iulog,'(a,i4)') 'Average J-values. NOTE: E3SM levels are indexed with the ground =',pver
         write(iulog,'(1x,a,72(a6,3x))') 'Lev ',(TITLEJX(JIND(K)),K=1,NRATJ)
         write(iulog,*)''
         do level_number = 1,pver
            write(iulog,'(i3,1p,72e9.2)')level_number,(photos(column_number,level_number,K),K=1,phtcnt)
         enddo
      endif

     enddo !!!! end of loop over columns in chunk


     
!pjc      goto 92
!pjc   91 stop 'error in opening .dat file'
!pjc   92 stop

!!!$OMP END MASTER    ! to activate this line, only use 1 exclamation mark

   end subroutine cloudJ_interface

 end module UCI_cloudJ_interface

    
