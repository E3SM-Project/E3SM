!===============================================================================
! SVN $Id: seq_diag_mct.F90 59750 2014-05-01 15:17:20Z sacks $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/drv/seq_mct/trunk_tags/drvseq5_0_12/driver/seq_diag_mct.F90 $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: seq_diag_mod -- computes spatial \& time averages of fluxed quatities
!
! !DESCRIPTION:
!    The coupler is required to do certain diagnostics, those calculations are
!    located in this module.
!
! !REMARKS:
!    CESM sign convention for fluxes is positive downward with hierarchy being
!       atm/glc/lnd/rof/ice/ocn
!    Sign convention:
!       positive value <=> the model is gaining water, heat, momentum, etc.
!    Unit convention:
!       heat flux     ~ W/m^2
!       momentum flux ~ N/m^2
!       water flux    ~ (kg/s)/m^2
!       salt  flux    ~ (kg/s)/m^2
!
! !REVISION HISTORY:
!    2012-aug-20 - T. Craig    - add rof component
!    2008-jul-10 - T. Craig    - updated budget implementation
!    2007-may-07 - B. Kauffman - initial port to cpl7. 
!    2002-nov-21 - R. Jacob    - initial port to cpl6. 
!    199x-mmm-dd - B. Kauffman - original version in cpl4.
!
! !INTERFACE: ------------------------------------------------------------------

module seq_diag_mct
  
! !USES:

   use shr_kind_mod, only: r8 => shr_kind_r8, in=>shr_kind_in
   use shr_kind_mod, only: i8 => shr_kind_i8,  cl=>shr_kind_cl
   use shr_sys_mod       ! system calls
   use shr_mpi_mod       ! mpi wrappers
   use shr_const_mod     ! shared constants
   use mct_mod           ! mct wrappers
   use esmf

   use seq_comm_mct  ! mpi comm groups & related
   use seq_timemgr_mod
   use component_type_mod

   implicit none
   save
   private

! !PUBLIC TYPES:

   ! none

!PUBLIC MEMBER FUNCTIONS:

   public seq_diag_zero_mct
   public seq_diag_atm_mct
   public seq_diag_lnd_mct
   public seq_diag_rof_mct
   public seq_diag_glc_mct
   public seq_diag_ocn_mct
   public seq_diag_ice_mct
   public seq_diag_accum_mct
   public seq_diag_sum0_mct
   public seq_diag_print_mct
   public seq_diag_avect_mct
   public seq_diag_avdiff_mct

!EOP

   !----------------------------------------------------------------------------
   ! Local data
   !----------------------------------------------------------------------------

   !----- local constants -----
   real(r8),parameter :: HFLXtoWFLX = & ! water flux implied by latent heat of fusion
   &  - (shr_const_ocn_ref_sal-shr_const_ice_ref_sal) / &
   &    (shr_const_ocn_ref_sal*shr_const_latice)

   !--- C for component ---
   !--- "r" is recieve in the coupler, "s" is send from the coupler

   integer(in),parameter :: c_size = 22

   integer(in),parameter :: c_atm_as   = 1 ! model index: atm
   integer(in),parameter :: c_atm_ar   = 2 ! model index: atm
   integer(in),parameter :: c_inh_is   = 3 ! model index: ice, northern
   integer(in),parameter :: c_inh_ir   = 4 ! model index: ice, northern
   integer(in),parameter :: c_ish_is   = 5 ! model index: ice, southern
   integer(in),parameter :: c_ish_ir   = 6 ! model index: ice, southern
   integer(in),parameter :: c_lnd_ls   = 7 ! model index: lnd
   integer(in),parameter :: c_lnd_lr   = 8 ! model index: lnd
   integer(in),parameter :: c_ocn_os   = 9 ! model index: ocn
   integer(in),parameter :: c_ocn_or   =10 ! model index: ocn
   integer(in),parameter :: c_rof_rs   =11 ! model index: rof
   integer(in),parameter :: c_rof_rr   =12 ! model index: rof
   integer(in),parameter :: c_glc_gs   =13 ! model index: glc
   integer(in),parameter :: c_glc_gr   =14 ! model index: glc
   ! --- on atm grid ---
   integer(in),parameter :: c_inh_as   =15 ! model index: ice, northern
   integer(in),parameter :: c_inh_ar   =16 ! model index: ice, northern
   integer(in),parameter :: c_ish_as   =17 ! model index: ice, southern
   integer(in),parameter :: c_ish_ar   =18 ! model index: ice, southern
   integer(in),parameter :: c_lnd_as   =19 ! model index: lnd
   integer(in),parameter :: c_lnd_ar   =20 ! model index: lnd
   integer(in),parameter :: c_ocn_as   =21 ! model index: ocn
   integer(in),parameter :: c_ocn_ar   =22 ! model index: ocn

   character(len=8),parameter :: cname(c_size) = &
      (/' c2a_atm',' a2c_atm',' c2i_inh',' i2c_inh',' c2i_ish',' i2c_ish', &
        ' c2l_lnd',' l2c_lnd',' c2o_ocn',' o2c_ocn',' c2r_rof',' r2c_rof', &
        ' c2g_glc',' g2c_glc', &
        ' c2a_inh',' a2c_inh',' c2a_ish',' a2c_ish', &
        ' c2a_lnd',' a2c_lnd',' c2a_ocn',' a2c_ocn' /)

   !--- F for field ---

   integer(in),parameter :: f_size = 17
   integer(in),parameter :: f_a    = 1    ! index for area
   integer(in),parameter :: f_h    = 2    ! 1st index for heat
   integer(in),parameter :: f_w    = 11   ! 1st index for water

   integer(in),parameter :: f_area    = 1 ! area (wrt to unit sphere)
   integer(in),parameter :: f_hfrz    = 2 ! heat : latent, freezing
   integer(in),parameter :: f_hmelt   = 3 ! heat : latent, melting
   integer(in),parameter :: f_hswnet  = 4 ! heat : short wave, net
   integer(in),parameter :: f_hlwdn   = 5 ! heat : longwave down
   integer(in),parameter :: f_hlwup   = 6 ! heat : longwave up
   integer(in),parameter :: f_hlatv   = 7 ! heat : latent, vaporization
   integer(in),parameter :: f_hlatf   = 8 ! heat : latent, fusion, snow       
   integer(in),parameter :: f_hioff   = 9 ! heat : latent, fusion, frozen runoff
   integer(in),parameter :: f_hsen    =10 ! heat : sensible
   integer(in),parameter :: f_wfrz    =11 ! water: freezing
   integer(in),parameter :: f_wmelt   =12 ! water: melting
   integer(in),parameter :: f_wrain   =13 ! water: precip, liquid
   integer(in),parameter :: f_wsnow   =14 ! water: precip, frozen
   integer(in),parameter :: f_wevap   =15 ! water: evaporation
   integer(in),parameter :: f_wroff   =16 ! water: runoff/flood
   integer(in),parameter :: f_wioff   =17 ! water: frozen runoff

   character(len=8),parameter :: fname(f_size) = &
      (/'    area',' hfreeze','   hmelt','  hnetsw','   hlwdn', &
        '   hlwup',' hlatvap',' hlatfus','  hiroff','    hsen', &
        ' wfreeze','   wmelt','   wrain','   wsnow', &
        '   wevap',' wrunoff',' wfrzrof' /)

   !--- P for period ---

   integer(in),parameter :: p_size = 5

   integer(in),parameter :: p_inst = 1
   integer(in),parameter :: p_day  = 2
   integer(in),parameter :: p_mon  = 3
   integer(in),parameter :: p_ann  = 4
   integer(in),parameter :: p_inf  = 5

   character(len=8),parameter :: pname(p_size) = &
      (/'    inst','   daily',' monthly','  annual','all_time' /)

! !PUBLIC DATA MEMBERS

   !--- time-averaged (annual?) global budge diagnostics ---
   !--- note: call sum0 then save budg_dataG and budg_ns on restart from/to root pe ---
   real(r8),public :: budg_dataL(f_size,c_size,p_size) ! local sum, valid on all pes
   real(r8),public :: budg_dataG(f_size,c_size,p_size) ! global sum, valid only on root pe
   real(r8),public :: budg_ns   (f_size,c_size,p_size) ! counter, valid only on root pe

   character(len=*),parameter :: afldname  = 'aream'
   character(len=*),parameter :: latname   = 'lat'
   character(len=*),parameter :: afracname = 'afrac'
   character(len=*),parameter :: lfracname = 'lfrac'
   character(len=*),parameter :: ofracname = 'ofrac'
   character(len=*),parameter :: ifracname = 'ifrac'

   character(*),parameter :: modName = "(seq_diag_mct) "

   integer(in),parameter :: debug = 0 ! internal debug level

! !PRIVATE DATA MEMBERS

   integer :: index_a2x_Faxa_swnet
   integer :: index_a2x_Faxa_lwdn
   integer :: index_a2x_Faxa_rainc
   integer :: index_a2x_Faxa_rainl
   integer :: index_a2x_Faxa_snowc
   integer :: index_a2x_Faxa_snowl

   integer :: index_x2a_Faxx_lwup
   integer :: index_x2a_Faxx_lat
   integer :: index_x2a_Faxx_sen
   integer :: index_x2a_Faxx_evap

   integer :: index_l2x_Fall_swnet
   integer :: index_l2x_Fall_lwup
   integer :: index_l2x_Fall_lat
   integer :: index_l2x_Fall_sen
   integer :: index_l2x_Fall_evap
   integer :: index_l2x_Flrl_rofl
   integer :: index_l2x_Flrl_rofi

   integer :: index_x2l_Faxa_lwdn
   integer :: index_x2l_Faxa_rainc
   integer :: index_x2l_Faxa_rainl
   integer :: index_x2l_Faxa_snowc
   integer :: index_x2l_Faxa_snowl
   integer :: index_x2l_Flrr_flood

   integer :: index_r2x_Forr_rofl
   integer :: index_r2x_Forr_rofi
   integer :: index_r2x_Firr_rofi
   integer :: index_r2x_Flrr_flood

   integer :: index_x2r_Flrl_rofl
   integer :: index_x2r_Flrl_rofi

   integer :: index_o2x_Fioo_q

   integer :: index_xao_Faox_lwup
   integer :: index_xao_Faox_lat
   integer :: index_xao_Faox_sen
   integer :: index_xao_Faox_evap

   integer :: index_x2o_Foxx_lwup
   integer :: index_x2o_Foxx_lat
   integer :: index_x2o_Foxx_sen
   integer :: index_x2o_Foxx_evap
   integer :: index_x2o_Foxx_swnet
   integer :: index_x2o_Foxx_rofl
   integer :: index_x2o_Foxx_rofi
   integer :: index_x2o_Faxa_lwdn
   integer :: index_x2o_Faxa_rain
   integer :: index_x2o_Faxa_snow
   integer :: index_x2o_Fioi_melth
   integer :: index_x2o_Fioi_meltw

   integer :: index_i2x_Fioi_melth
   integer :: index_i2x_Fioi_meltw
   integer :: index_i2x_Faii_swnet
   integer :: index_i2x_Fioi_swpen
   integer :: index_i2x_Faii_lwup
   integer :: index_i2x_Faii_lat
   integer :: index_i2x_Faii_sen
   integer :: index_i2x_Faii_evap

   integer :: index_x2i_Faxa_lwdn
   integer :: index_x2i_Faxa_rain
   integer :: index_x2i_Faxa_snow
   integer :: index_x2i_Fioo_q
   integer :: index_x2i_Fixx_rofi

   integer :: index_g2x_Fogg_rofl
   integer :: index_g2x_Fogg_rofi
   integer :: index_g2x_Figg_rofi

!===============================================================================
contains
!===============================================================================

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: seq_diag_zero_mct - zero out global budget diagnostic data.
!
! !DESCRIPTION:
!    Zero out global budget diagnostic data.
!
! !REVISION HISTORY:
!    2008-jul-11 - T. Craig - update
!
! !INTERFACE: ------------------------------------------------------------------

subroutine seq_diag_zero_mct(EClock,mode)

! !INPUT/OUTPUT PARAMETERS:

   type(ESMF_Clock), intent(in),optional :: EClock
   character(len=*), intent(in),optional :: mode

!EOP

   integer(IN) :: ip,yr,mon,day,sec
   !----- formats -----
   character(*),parameter :: subName = '(seq_diag_zero_mct) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   if (.not. present(EClock) .and. .not. present(mode)) then
      call shr_sys_abort(subName//' ERROR EClock or mode should be present')
   endif

   if (present(EClock)) then
      call seq_timemgr_EClockGetData(EClock,curr_yr=yr, &
         curr_mon=mon,curr_day=day,curr_tod=sec)

      do ip = 1,p_size
         if (ip == p_inst) then
            budg_dataL(:,:,ip) = 0.0_r8
            budg_dataG(:,:,ip) = 0.0_r8
            budg_ns(:,:,ip) = 0.0_r8
         endif
         if (ip==p_day .and. sec==0) then
            budg_dataL(:,:,ip) = 0.0_r8
            budg_dataG(:,:,ip) = 0.0_r8
            budg_ns(:,:,ip) = 0.0_r8
         endif
         if (ip==p_mon .and. day==1 .and. sec==0) then
            budg_dataL(:,:,ip) = 0.0_r8
            budg_dataG(:,:,ip) = 0.0_r8
            budg_ns(:,:,ip) = 0.0_r8
         endif
         if (ip==p_ann .and. mon==1 .and. day==1 .and. sec==0) then
            budg_dataL(:,:,ip) = 0.0_r8
            budg_dataG(:,:,ip) = 0.0_r8
            budg_ns(:,:,ip) = 0.0_r8
         endif
      enddo
   endif

   if (present(mode)) then
      if (trim(mode) == 'inst') then
         budg_dataL(:,:,p_inst) = 0.0_r8
         budg_dataG(:,:,p_inst) = 0.0_r8
         budg_ns(:,:,p_inst) = 0.0_r8
      elseif (trim(mode) == 'day') then
         budg_dataL(:,:,p_day) = 0.0_r8
         budg_dataG(:,:,p_day) = 0.0_r8
         budg_ns(:,:,p_day) = 0.0_r8
      elseif (trim(mode) == 'mon') then
         budg_dataL(:,:,p_mon) = 0.0_r8
         budg_dataG(:,:,p_mon) = 0.0_r8
         budg_ns(:,:,p_mon) = 0.0_r8
      elseif (trim(mode) == 'ann') then
         budg_dataL(:,:,p_ann) = 0.0_r8
         budg_dataG(:,:,p_ann) = 0.0_r8
         budg_ns(:,:,p_ann) = 0.0_r8
      elseif (trim(mode) == 'inf') then
         budg_dataL(:,:,p_inf) = 0.0_r8
         budg_dataG(:,:,p_inf) = 0.0_r8
         budg_ns(:,:,p_inf) = 0.0_r8
      elseif (trim(mode) == 'all') then
         budg_dataL(:,:,:) = 0.0_r8
         budg_dataG(:,:,:) = 0.0_r8
         budg_ns(:,:,:) = 0.0_r8
      else
         call shr_sys_abort(subname//' ERROR in mode '//trim(mode))
      endif
   endif

end subroutine seq_diag_zero_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: seq_diag_accum_mct - accum out global budget diagnostic data.
!
! !DESCRIPTION:
!    Accum out global budget diagnostic data.
!
! !REVISION HISTORY:
!    2008-jul-11 - T. Craig - update
!
! !INTERFACE: ------------------------------------------------------------------

subroutine seq_diag_accum_mct()

! !INPUT/OUTPUT PARAMETERS:

!EOP

   integer(in) :: ip

   !----- formats -----
   character(*),parameter :: subName = '(seq_diag_accum_mct) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   do ip = p_inst+1,p_size
      budg_dataL(:,:,ip) = budg_dataL(:,:,ip) + budg_dataL(:,:,p_inst)
   enddo
   budg_ns(:,:,:) = budg_ns(:,:,:) + 1.0_r8

end subroutine seq_diag_accum_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: seq_diag_sum0_mct - sum local to global on root
!
! !DESCRIPTION:
!    Sum local values to global on root
!
! !REVISION HISTORY:
!    2008-jul-19 - T. Craig - update
!
! !INTERFACE: ------------------------------------------------------------------

subroutine seq_diag_sum0_mct()

! !INPUT/OUTPUT PARAMETERS:

!EOP

   real(r8) :: budg_dataGtmp(f_size,c_size,p_size) ! temporary sum
   integer(in)      :: mpicom      ! mpi comm
   !----- formats -----
   character(*),parameter :: subName = '(seq_diag_sum0_mct) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   call seq_comm_setptrs(CPLID,mpicom=mpicom)
   budg_dataGtmp = 0.0_r8
   call shr_mpi_sum(budg_dataL,budg_dataGtmp,mpicom,subName)
   budg_dataG = budg_dataG + budg_dataGtmp
   budg_dataL = 0.0_r8

end subroutine seq_diag_sum0_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: seq_diag_atm_mct - compute global atm input/output flux diagnostics
!
! !DESCRIPTION:
!     Compute global atm input/output flux diagnostics
!
! !REVISION HISTORY:
!    2008-jul-10 - T. Craig - update
!
! !INTERFACE: ------------------------------------------------------------------

subroutine seq_diag_atm_mct( atm, frac_a, do_a2x, do_x2a )

! !INPUT/OUTPUT PARAMETERS:

   type(component_type), intent(in) :: atm    ! component type for instance1
   type(mct_aVect)     , intent(in) :: frac_a ! frac bundle
   logical, optional   , intent(in) :: do_a2x             
   logical, optional   , intent(in) :: do_x2a             

!EOP

   !----- local -----
   type(mct_aVect), pointer :: a2x_a        ! model to drv bundle
   type(mct_aVect), pointer :: x2a_a        ! drv to model bundle
   type(mct_ggrid), pointer :: dom_a
   integer(in)              :: k,n,ic,if,ip      ! generic index
   integer(in)              :: kArea             ! index of area field in aVect
   integer(in)              :: kLat              ! index of lat field in aVect
   integer(in)              :: kl,ka,ko,ki       ! fraction indices
   integer(in)              :: lSize             ! size of aVect
   real(r8)                 :: da,di,do,dl       ! area of a grid cell
   logical,save             :: first_time = .true.

   !----- formats -----
   character(*),parameter :: subName = '(seq_diag_atm_mct) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   dom_a => component_get_dom_cx(atm)
   a2x_a => component_get_c2x_cx(atm)  
   x2a_a => component_get_x2c_cx(atm)  

   kArea = mct_aVect_indexRA(dom_a%data,afldname)
   kLat  = mct_aVect_indexRA(dom_a%data,latname)
   ka    = mct_aVect_indexRA(frac_a,afracname)
   kl    = mct_aVect_indexRA(frac_a,lfracname)
   ko    = mct_aVect_indexRA(frac_a,ofracname)
   ki    = mct_aVect_indexRA(frac_a,ifracname)

   !---------------------------------------------------------------------------
   ! add values found in this bundle to the budget table
   !---------------------------------------------------------------------------

   ip = p_inst

   if (present(do_a2x)) then
      if (first_time) then
         index_a2x_Faxa_swnet  = mct_aVect_indexRA(a2x_a,'Faxa_swnet')
         index_a2x_Faxa_lwdn   = mct_aVect_indexRA(a2x_a,'Faxa_lwdn')
         index_a2x_Faxa_rainc  = mct_aVect_indexRA(a2x_a,'Faxa_rainc')
         index_a2x_Faxa_rainl  = mct_aVect_indexRA(a2x_a,'Faxa_rainl')
         index_a2x_Faxa_snowc  = mct_aVect_indexRA(a2x_a,'Faxa_snowc')
         index_a2x_Faxa_snowl  = mct_aVect_indexRA(a2x_a,'Faxa_snowl')
      end if

      lSize = mct_avect_lSize(a2x_a)
      do n=1,lSize
      do k=1,4

         if (k == 1) then
            ic = c_atm_ar
            da = -dom_a%data%rAttr(kArea,n) * frac_a%rAttr(ka,n)
         elseif (k == 2) then
            ic = c_lnd_ar
            da =  dom_a%data%rAttr(kArea,n) * frac_a%rAttr(kl,n)
         elseif (k == 3) then
            ic = c_ocn_ar
            da =  dom_a%data%rAttr(kArea,n) * frac_a%rAttr(ko,n)
         elseif (k == 4) then
            if (dom_a%data%rAttr(kLat,n) > 0.0_r8) then
               ic = c_inh_ar
            else
               ic = c_ish_ar
            endif
            da = dom_a%data%rAttr(kArea,n) * frac_a%rAttr(ki,n)
         endif

         if = f_area  ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + da
         if = f_hswnet; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + da*a2x_a%rAttr(index_a2x_Faxa_swnet,n)
         if = f_hlwdn ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + da*a2x_a%rAttr(index_a2x_Faxa_lwdn,n)
         if = f_wrain ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + da*a2x_a%rAttr(index_a2x_Faxa_rainc,n) &
                                                                    + da*a2x_a%rAttr(index_a2x_Faxa_rainl,n)
         if = f_wsnow ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + da*a2x_a%rAttr(index_a2x_Faxa_snowc,n) &
                                                                    + da*a2x_a%rAttr(index_a2x_Faxa_snowl,n)
      enddo
      enddo
      ! --- heat implied by snow flux ---
      ic = c_atm_ar;  budg_dataL(f_hlatf,ic,ip) = -budg_dataL(f_wsnow,ic,ip)*shr_const_latice
      ic = c_lnd_ar;  budg_dataL(f_hlatf,ic,ip) = -budg_dataL(f_wsnow,ic,ip)*shr_const_latice
      ic = c_ocn_ar;  budg_dataL(f_hlatf,ic,ip) = -budg_dataL(f_wsnow,ic,ip)*shr_const_latice
      ic = c_inh_ar;  budg_dataL(f_hlatf,ic,ip) = -budg_dataL(f_wsnow,ic,ip)*shr_const_latice
      ic = c_ish_ar;  budg_dataL(f_hlatf,ic,ip) = -budg_dataL(f_wsnow,ic,ip)*shr_const_latice
   end if

   if (present(do_x2a)) then
      if (first_time) then
         index_x2a_Faxx_lwup   = mct_aVect_indexRA(x2a_a,'Faxx_lwup')
         index_x2a_Faxx_lat    = mct_aVect_indexRA(x2a_a,'Faxx_lat')
         index_x2a_Faxx_sen    = mct_aVect_indexRA(x2a_a,'Faxx_sen')
         index_x2a_Faxx_evap   = mct_aVect_indexRA(x2a_a,'Faxx_evap')
      end if

      lSize = mct_avect_lSize(x2a_a)
      do n=1,lSize
      do k=1,4

         if (k == 1) then
            ic = c_atm_as
            da = -dom_a%data%rAttr(kArea,n) * frac_a%rAttr(ka,n)
         elseif (k == 2) then
            ic = c_lnd_as
            da =  dom_a%data%rAttr(kArea,n) * frac_a%rAttr(kl,n)
         elseif (k == 3) then
            ic = c_ocn_as
            da =  dom_a%data%rAttr(kArea,n) * frac_a%rAttr(ko,n)
         elseif (k == 4) then
            if (dom_a%data%rAttr(kLat,n) > 0.0_r8) then
               ic = c_inh_as
            else
               ic = c_ish_as
            endif
            da = dom_a%data%rAttr(kArea,n) * frac_a%rAttr(ki,n)
         endif

         if = f_area ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + da
         if = f_hlwup; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + da*x2a_a%rAttr(index_x2a_Faxx_lwup,n)
         if = f_hlatv; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + da*x2a_a%rAttr(index_x2a_Faxx_lat,n)
         if = f_hsen ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + da*x2a_a%rAttr(index_x2a_Faxx_sen,n)
         if = f_wevap; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + da*x2a_a%rAttr(index_x2a_Faxx_evap,n)

      enddo
      enddo
   end if

   first_time = .false.

end subroutine seq_diag_atm_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: seq_diag_lnd_mct - compute global lnd input/output flux diagnostics
!
! !DESCRIPTION:
!     Compute global lnd input/output flux diagnostics
!
! !REVISION HISTORY:
!    2008-jul-10 - T. Craig - update
!
! !INTERFACE: ------------------------------------------------------------------

subroutine seq_diag_lnd_mct( lnd, frac_l, do_l2x, do_x2l)

   type(component_type), intent(in) :: lnd    ! component type for instance1
   type(mct_aVect)     , intent(in) :: frac_l ! frac bundle
   logical, optional   , intent(in) :: do_l2x             
   logical, optional   , intent(in) :: do_x2l             

!EOP

   !----- local -----
   type(mct_aVect), pointer :: l2x_l        ! model to drv bundle
   type(mct_aVect), pointer :: x2l_l        ! drv to model bundle
   type(mct_ggrid), pointer :: dom_l
   integer(in)              :: k,n,ic,if,ip ! generic index
   integer(in)              :: kArea        ! index of area field in aVect
   integer(in)              :: kLat         ! index of lat field in aVect
   integer(in)              :: kl,ka,ko,ki  ! fraction indices
   integer(in)              :: lSize        ! size of aVect
   real(r8)                 :: da,di,do,dl  ! area of a grid cell
   logical,save             :: first_time = .true.

   !----- formats -----
   character(*),parameter :: subName = '(seq_diag_lnd_mct) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   !---------------------------------------------------------------------------
   ! add values found in this bundle to the budget table
   !---------------------------------------------------------------------------

   dom_l => component_get_dom_cx(lnd)
   l2x_l => component_get_c2x_cx(lnd)  
   x2l_l => component_get_x2c_cx(lnd)  

   ip = p_inst

   kArea = mct_aVect_indexRA(dom_l%data,afldname)
   kl    = mct_aVect_indexRA(frac_l,lfracname)

   if (present(do_l2x)) then
      if (first_time) then
         index_l2x_Fall_swnet  = mct_aVect_indexRA(l2x_l,'Fall_swnet')
         index_l2x_Fall_lwup   = mct_aVect_indexRA(l2x_l,'Fall_lwup')
         index_l2x_Fall_lat    = mct_aVect_indexRA(l2x_l,'Fall_lat')
         index_l2x_Fall_sen    = mct_aVect_indexRA(l2x_l,'Fall_sen')
         index_l2x_Fall_evap   = mct_aVect_indexRA(l2x_l,'Fall_evap')
         index_l2x_Flrl_rofl   = mct_aVect_indexRA(l2x_l,'Flrl_rofl')
         index_l2x_Flrl_rofi   = mct_aVect_indexRA(l2x_l,'Flrl_rofi')
      end if

      lSize = mct_avect_lSize(l2x_l)
      ic = c_lnd_lr
      do n=1,lSize
         dl =  dom_l%data%rAttr(kArea,n) * frac_l%rAttr(kl,n)
         if = f_area  ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + dl
         if = f_hswnet; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + dl*l2x_l%rAttr(index_l2x_Fall_swnet,n)
         if = f_hlwup ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + dl*l2x_l%rAttr(index_l2x_Fall_lwup,n)
         if = f_hlatv ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + dl*l2x_l%rAttr(index_l2x_Fall_lat,n)
         if = f_hsen  ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + dl*l2x_l%rAttr(index_l2x_Fall_sen,n)
         if = f_wevap ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + dl*l2x_l%rAttr(index_l2x_Fall_evap,n)
         if = f_wroff ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) - dl*l2x_l%rAttr(index_l2x_Flrl_rofl,n)
         if = f_wioff ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) - dl*l2x_l%rAttr(index_l2x_Flrl_rofi,n)
      end do
      budg_dataL(f_hioff,ic,ip) = -budg_dataL(f_wioff,ic,ip)*shr_const_latice
   end if

   if (present(do_x2l)) then
      if (first_time) then
         index_x2l_Faxa_lwdn   = mct_aVect_indexRA(x2l_l,'Faxa_lwdn')
         index_x2l_Faxa_rainc  = mct_aVect_indexRA(x2l_l,'Faxa_rainc')
         index_x2l_Faxa_rainl  = mct_aVect_indexRA(x2l_l,'Faxa_rainl')
         index_x2l_Faxa_snowc  = mct_aVect_indexRA(x2l_l,'Faxa_snowc')
         index_x2l_Faxa_snowl  = mct_aVect_indexRA(x2l_l,'Faxa_snowl')
         index_x2l_Flrr_flood  = mct_aVect_indexRA(x2l_l,'Flrr_flood')
      end if

      lSize = mct_avect_lSize(x2l_l)
      ic = c_lnd_ls
      do n=1,lSize
         dl =  dom_l%data%rAttr(kArea,n) * frac_l%rAttr(kl,n)
         if = f_area ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + dl
         if = f_hlwdn; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + dl*x2l_l%rAttr(index_x2l_Faxa_lwdn,n)
         if = f_wrain; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + dl*x2l_l%rAttr(index_x2l_Faxa_rainc,n) &
                                                                   + dl*x2l_l%rAttr(index_x2l_Faxa_rainl,n)
         if = f_wsnow; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + dl*x2l_l%rAttr(index_x2l_Faxa_snowc,n) &
                                                                   + dl*x2l_l%rAttr(index_x2l_Faxa_snowl,n)
         if = f_wroff; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) - dl*x2l_l%rAttr(index_x2l_Flrr_flood,n)
      end do
      budg_dataL(f_hlatf,ic,ip) = -budg_dataL(f_wsnow,ic,ip)*shr_const_latice
   end if

   first_time = .false.

end subroutine seq_diag_lnd_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: seq_diag_rof_mct - compute global rof input/output flux diagnostics
!
! !DESCRIPTION:
!     Compute global rof input/output flux diagnostics
!
! !REVISION HISTORY:
!    2008-jul-10 - T. Craig - update
!
! !INTERFACE: ------------------------------------------------------------------

subroutine seq_diag_rof_mct( rof, frac_r)

   type(component_type), intent(in) :: rof    ! component type for instance1
   type(mct_aVect)     , intent(in) :: frac_r ! frac bundle

!EOP

   !----- local -----
   type(mct_aVect), pointer :: r2x_r
   type(mct_aVect), pointer :: x2r_r
   type(mct_ggrid), pointer :: dom_r
   integer(in)              :: k,n,ic,if,ip      ! generic index
   integer(in)              :: kArea             ! index of area field in aVect
   integer(in)              :: kLat              ! index of lat field in aVect
   integer(in)              :: kl,ka,ko,ki,kr    ! fraction indices
   integer(in)              :: lSize             ! size of aVect
   real(r8)                 :: da,di,do,dl,dr    ! area of a grid cell
   logical,save             :: first_time = .true.

   !----- formats -----
   character(*),parameter :: subName = '(seq_diag_rof_mct) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   !---------------------------------------------------------------------------
   ! add values found in this bundle to the budget table
   !---------------------------------------------------------------------------

   dom_r => component_get_dom_cx(rof)
   r2x_r => component_get_c2x_cx(rof)  
   x2r_r => component_get_x2c_cx(rof)  

   if (first_time) then
      index_x2r_Flrl_rofl  = mct_aVect_indexRA(x2r_r,'Flrl_rofl')
      index_x2r_Flrl_rofi  = mct_aVect_indexRA(x2r_r,'Flrl_rofi')
   end if

   ip = p_inst
   ic = c_rof_rr
   kArea = mct_aVect_indexRA(dom_r%data,afldname)
   lSize = mct_avect_lSize(x2r_r)
   do n=1,lSize
      dr =  dom_r%data%rAttr(kArea,n)
      if = f_wroff; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + dr*x2r_r%rAttr(index_x2r_Flrl_rofl,n)
      if = f_wioff; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + dr*x2r_r%rAttr(index_x2r_Flrl_rofi,n)
   end do
   budg_dataL(f_hioff,ic,ip) = -budg_dataL(f_wioff,ic,ip)*shr_const_latice

   if (first_time) then
      index_r2x_Forr_rofl   = mct_aVect_indexRA(r2x_r,'Forr_rofl')
      index_r2x_Forr_rofi   = mct_aVect_indexRA(r2x_r,'Forr_rofi')
      index_r2x_Firr_rofi   = mct_aVect_indexRA(r2x_r,'Firr_rofi')
      index_r2x_Flrr_flood  = mct_aVect_indexRA(r2x_r,'Flrr_flood')
   end if

   ip = p_inst
   ic = c_rof_rs
   kArea = mct_aVect_indexRA(dom_r%data,afldname)
   lSize = mct_avect_lSize(r2x_r)
   do n=1,lSize
      dr =  dom_r%data%rAttr(kArea,n)
      if = f_wroff; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) - dr*r2x_r%rAttr(index_r2x_Forr_rofl,n) &
                                                                + dr*r2x_r%rAttr(index_r2x_Flrr_flood,n)
      if = f_wioff; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) - dr*r2x_r%rAttr(index_r2x_Forr_rofi,n) &
                                                                - dr*r2x_r%rAttr(index_r2x_Firr_rofi,n)
   end do
   budg_dataL(f_hioff,ic,ip) = -budg_dataL(f_wioff,ic,ip)*shr_const_latice

   first_time = .false.

end subroutine seq_diag_rof_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: seq_diag_glc_mct - compute global glc input/output flux diagnostics
!
! !DESCRIPTION:
!     Compute global glc input/output flux diagnostics
!
! !REVISION HISTORY:
!    2008-jul-10 - T. Craig - update
!
! !INTERFACE: ------------------------------------------------------------------

subroutine seq_diag_glc_mct( glc, frac_g)

   type(component_type), intent(in) :: glc    ! component type for instance1
   type(mct_aVect)     , intent(in) :: frac_g ! frac bundle

!EOP

   !----- local -----
   type(mct_aVect), pointer :: g2x_g
   type(mct_aVect), pointer :: x2g_g
   type(mct_ggrid), pointer :: dom_g
   integer(in)              :: k,n,ic,if,ip      ! generic index
   integer(in)              :: kArea             ! index of area field in aVect
   integer(in)              :: kLat              ! index of lat field in aVect
   integer(in)              :: kl,ka,ko,ki,kr,kg ! fraction indices
   integer(in)              :: lSize             ! size of aVect
   real(r8)                 :: da,di,do,dl,dr,dg ! area of a grid cell
   logical,save             :: first_time = .true.

   !----- formats -----
   character(*),parameter :: subName = '(seq_diag_glc_mct) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   !---------------------------------------------------------------------------
   ! add values found in this bundle to the budget table
   !---------------------------------------------------------------------------

   dom_g => component_get_dom_cx(glc)
   g2x_g => component_get_c2x_cx(glc)  
   x2g_g => component_get_x2c_cx(glc)  

   if (first_time) then
      index_g2x_Fogg_rofl   = mct_aVect_indexRA(g2x_g,'Fogg_rofl')
      index_g2x_Fogg_rofi   = mct_aVect_indexRA(g2x_g,'Fogg_rofi')
      index_g2x_Figg_rofi   = mct_aVect_indexRA(g2x_g,'Figg_rofi')
   end if

   ip = p_inst
   ic = c_glc_gs
   kArea = mct_aVect_indexRA(dom_g%data,afldname)
   lSize = mct_avect_lSize(g2x_g)
   do n=1,lSize
      dg =  dom_g%data%rAttr(kArea,n)
      if = f_wroff; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) - dg*g2x_g%rAttr(index_g2x_Fogg_rofl,n)
      if = f_wioff; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) - dg*g2x_g%rAttr(index_g2x_Fogg_rofi,n) &
                                                                - dg*g2x_g%rAttr(index_g2x_Figg_rofi,n)
   end do
   budg_dataL(f_hioff,ic,ip) = -budg_dataL(f_wioff,ic,ip)*shr_const_latice

   first_time = .false.

end subroutine seq_diag_glc_mct

!BOP ===========================================================================
!
! !IROUTINE: seq_diag_ocn_mct - compute global ocn input/output flux diagnostics
!
! !DESCRIPTION:
!     Compute global ocn input/output flux diagnostics
!
! !REVISION HISTORY:
!    2008-jul-10 - T. Craig - update
!
! !INTERFACE: ------------------------------------------------------------------

subroutine seq_diag_ocn_mct( ocn, xao_o, frac_o, do_o2x, do_x2o, do_xao)

   type(component_type) , intent(in)          :: ocn    ! component type for instance1
   type(mct_aVect)      , intent(in)          :: frac_o ! frac bundle
   type(mct_aVect)      , intent(in)          :: xao_o  
   logical              , intent(in),optional :: do_o2x
   logical              , intent(in),optional :: do_x2o
   logical              , intent(in),optional :: do_xao

!EOP

   !----- local -----
   type(mct_aVect), pointer :: o2x_o        ! model to drv bundle
   type(mct_aVect), pointer :: x2o_o        ! drv to model bundle
   type(mct_ggrid), pointer :: dom_o
   integer(in)              :: k,n,if,ic,ip ! generic index
   integer(in)              :: kArea        ! index of area field in aVect
   integer(in)              :: kLat         ! index of lat field in aVect
   integer(in)              :: kl,ka,ko,ki  ! fraction indices
   integer(in)              :: lSize        ! size of aVect
   real(r8)                 :: da,di,do,dl  ! area of a grid cell
   logical,save             :: first_time = .true.

   !----- formats -----
   character(*),parameter :: subName = '(seq_diag_ocn_mct) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   if (.not. present(do_o2x) .and. &
       .not. present(do_x2o) .and. &
       .not. present(do_xao)) then
      call shr_sys_abort(subName//"ERROR: must input a bundle")
   end if

   !---------------------------------------------------------------------------
   ! add values found in this bundle to the budget table
   !---------------------------------------------------------------------------

   dom_o => component_get_dom_cx(ocn)
   o2x_o => component_get_c2x_cx(ocn)  
   x2o_o => component_get_x2c_cx(ocn)  

   ip = p_inst

   kArea = mct_aVect_indexRA(dom_o%data,afldname)
   ko    = mct_aVect_indexRA(frac_o,ofracname)
   ki    = mct_aVect_indexRA(frac_o,ifracname)

   if (present(do_o2x)) then
      if (first_time) then
         index_o2x_Fioo_q      = mct_aVect_indexRA(o2x_o,'Fioo_q')
      end if

      lSize = mct_avect_lSize(o2x_o)
      ic = c_ocn_or
      do n=1,lSize
         do =  dom_o%data%rAttr(kArea,n) * frac_o%rAttr(ko,n)
         di =  dom_o%data%rAttr(kArea,n) * frac_o%rAttr(ki,n)
         if = f_area; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + do
         if = f_hfrz; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + (do+di)*max(0.0_r8,o2x_o%rAttr(index_o2x_Fioo_q,n))
      end do
      budg_dataL(f_wfrz,ic,ip) = budg_dataL(f_hfrz,ic,ip) * HFLXtoWFLX
   end if

   if (present(do_xao)) then
      if (first_time) then
         index_xao_Faox_lwup   = mct_aVect_indexRA(xao_o,'Faox_lwup') 
         index_xao_Faox_lat    = mct_aVect_indexRA(xao_o,'Faox_lat')  
         index_xao_Faox_sen    = mct_aVect_indexRA(xao_o,'Faox_sen') 
         index_xao_Faox_evap   = mct_aVect_indexRA(xao_o,'Faox_evap')  
      end if

      lSize = mct_avect_lSize(xao_o)
      ic = c_ocn_or
      do n=1,lSize
         do =  dom_o%data%rAttr(kArea,n) * frac_o%rAttr(ko,n)
         if = f_hlwup; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + do*xao_o%rAttr(index_xao_Faox_lwup,n)
         if = f_hlatv; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + do*xao_o%rAttr(index_xao_Faox_lat,n)
         if = f_hsen ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + do*xao_o%rAttr(index_xao_Faox_sen,n)
         if = f_wevap; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + do*xao_o%rAttr(index_xao_Faox_evap,n)
      end do
   end if

   if (present(do_x2o)) then
      if (first_time) then
         index_x2o_Fioi_melth  = mct_aVect_indexRA(x2o_o,'Fioi_melth')  
         index_x2o_Fioi_meltw  = mct_aVect_indexRA(x2o_o,'Fioi_meltw') 
         index_x2o_Foxx_swnet  = mct_aVect_indexRA(x2o_o,'Foxx_swnet')
         index_x2o_Faxa_lwdn   = mct_aVect_indexRA(x2o_o,'Faxa_lwdn')
         index_x2o_Faxa_rain   = mct_aVect_indexRA(x2o_o,'Faxa_rain') 
         index_x2o_Faxa_snow   = mct_aVect_indexRA(x2o_o,'Faxa_snow')  
         index_x2o_Foxx_lwup   = mct_aVect_indexRA(x2o_o,'Foxx_lwup') 
         index_x2o_Foxx_lat    = mct_aVect_indexRA(x2o_o,'Foxx_lat')  
         index_x2o_Foxx_sen    = mct_aVect_indexRA(x2o_o,'Foxx_sen') 
         index_x2o_Foxx_evap   = mct_aVect_indexRA(x2o_o,'Foxx_evap')  
         index_x2o_Foxx_rofl   = mct_aVect_indexRA(x2o_o,'Foxx_rofl')
         index_x2o_Foxx_rofi   = mct_aVect_indexRA(x2o_o,'Foxx_rofi')
      end if

      if (.not. present(do_xao)) then
         ! these are in x2o but they really are the atm/ocean flux 
         ! computed in the coupler and are "like" an o2x
         lSize = mct_avect_lSize(x2o_o)
         ic = c_ocn_or
         do n=1,lSize
            do =  dom_o%data%rAttr(kArea,n) * frac_o%rAttr(ko,n)
            di =  dom_o%data%rAttr(kArea,n) * frac_o%rAttr(ki,n)
            if = f_hlwup; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + (do+di)*x2o_o%rAttr(index_x2o_Foxx_lwup,n)
            if = f_hlatv; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + (do+di)*x2o_o%rAttr(index_x2o_Foxx_lat,n)
            if = f_hsen ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + (do+di)*x2o_o%rAttr(index_x2o_Foxx_sen,n)
            if = f_wevap; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + (do+di)*x2o_o%rAttr(index_x2o_Foxx_evap,n)
         end do
      endif

      lSize = mct_avect_lSize(x2o_o)
      ic = c_ocn_os
      do n=1,lSize
         do =  dom_o%data%rAttr(kArea,n) * frac_o%rAttr(ko,n)
         di =  dom_o%data%rAttr(kArea,n) * frac_o%rAttr(ki,n)
         if = f_area  ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + do
         if = f_wmelt ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + (do+di)*x2o_o%rAttr(index_x2o_Fioi_meltw,n)
         if = f_hmelt ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + (do+di)*x2o_o%rAttr(index_x2o_Fioi_melth,n)
         if = f_hswnet; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + (do+di)*x2o_o%rAttr(index_x2o_Foxx_swnet,n)
         if = f_hlwdn ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + (do+di)*x2o_o%rAttr(index_x2o_Faxa_lwdn,n)
         if = f_wrain ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + (do+di)*x2o_o%rAttr(index_x2o_Faxa_rain,n)
         if = f_wsnow ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + (do+di)*x2o_o%rAttr(index_x2o_Faxa_snow,n)
         if = f_wroff ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + (do+di)*x2o_o%rAttr(index_x2o_Foxx_rofl,n)
         if = f_wioff ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + (do+di)*x2o_o%rAttr(index_x2o_Foxx_rofi,n)
      end do
      budg_dataL(f_hlatf,ic,ip) = -budg_dataL(f_wsnow,ic,ip)*shr_const_latice
      budg_dataL(f_hioff,ic,ip) = -budg_dataL(f_wioff,ic,ip)*shr_const_latice
   end if

   first_time = .false.

end subroutine seq_diag_ocn_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: seq_diag_ice_mct - compute global ice input/output flux diagnostics
!
! !DESCRIPTION:
!     Compute global ice input/output flux diagnostics
!
! !REVISION HISTORY:
!    2008-jul-10 - T. Craig - update
!
! !INTERFACE: ------------------------------------------------------------------

subroutine seq_diag_ice_mct( ice, frac_i, do_i2x, do_x2i)

   type(component_type), intent(in)           :: ice    ! component type for instance1
   type(mct_aVect)     , intent(in)           :: frac_i ! frac bundle
   logical             , intent(in), optional :: do_i2x
   logical             , intent(in), optional :: do_x2i

!EOP

   !----- local -----
   type(mct_aVect), pointer :: i2x_i        ! model to drv bundle
   type(mct_aVect), pointer :: x2i_i        ! drv to model bundle
   type(mct_ggrid), pointer :: dom_i
   integer(in)              :: k,n,ic,if,ip ! generic index
   integer(in)              :: kArea        ! index of area field in aVect
   integer(in)              :: kLat         ! index of lat field in aVect
   integer(in)              :: kl,ka,ko,ki  ! fraction indices
   integer(in)              :: lSize        ! size of aVect
   real(r8)                 :: da,di,do,dl  ! area of a grid cell
   logical,save             :: first_time = .true.

   !----- formats -----
   character(*),parameter :: subName = '(seq_diag_ice_mct) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   !---------------------------------------------------------------------------
   ! add values found in this bundle to the budget table
   !---------------------------------------------------------------------------

   dom_i => component_get_dom_cx(ice)
   i2x_i => component_get_c2x_cx(ice)  
   x2i_i => component_get_x2c_cx(ice)  

   ip = p_inst

   kArea = mct_aVect_indexRA(dom_i%data,afldname)
   kLat  = mct_aVect_indexRA(dom_i%data,latname)
   ki    = mct_aVect_indexRA(frac_i,ifracname)
   ko    = mct_aVect_indexRA(frac_i,ofracname)

   if (present(do_i2x)) then
         index_i2x_Fioi_melth  = mct_aVect_indexRA(i2x_i,'Fioi_melth')
         index_i2x_Fioi_meltw  = mct_aVect_indexRA(i2x_i,'Fioi_meltw')
         index_i2x_Fioi_swpen  = mct_aVect_indexRA(i2x_i,'Fioi_swpen')
         index_i2x_Faii_swnet  = mct_aVect_indexRA(i2x_i,'Faii_swnet')
         index_i2x_Faii_lwup   = mct_aVect_indexRA(i2x_i,'Faii_lwup')
         index_i2x_Faii_lat    = mct_aVect_indexRA(i2x_i,'Faii_lat')
         index_i2x_Faii_sen    = mct_aVect_indexRA(i2x_i,'Faii_sen')
         index_i2x_Faii_evap   = mct_aVect_indexRA(i2x_i,'Faii_evap')

      lSize = mct_avect_lSize(i2x_i)
      do n=1,lSize
         if (dom_i%data%rAttr(kLat,n) > 0.0_r8) then
            ic = c_inh_ir
         else
            ic = c_ish_ir
         endif
         do =  dom_i%data%rAttr(kArea,n) * frac_i%rAttr(ko,n)
         di =  dom_i%data%rAttr(kArea,n) * frac_i%rAttr(ki,n)
         if = f_area  ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + di
         if = f_hmelt ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) - di*i2x_i%rAttr(index_i2x_Fioi_melth,n)
         if = f_wmelt ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) - di*i2x_i%rAttr(index_i2x_Fioi_meltw,n)
         if = f_hswnet; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + di*i2x_i%rAttr(index_i2x_Faii_swnet,n) &
                                                                    - di*i2x_i%rAttr(index_i2x_Fioi_swpen,n)
         if = f_hlwup ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + di*i2x_i%rAttr(index_i2x_Faii_lwup,n)
         if = f_hlatv ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + di*i2x_i%rAttr(index_i2x_Faii_lat,n)
         if = f_hsen  ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + di*i2x_i%rAttr(index_i2x_Faii_sen,n)
         if = f_wevap ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + di*i2x_i%rAttr(index_i2x_Faii_evap,n)
      end do
   end if

   if (present(do_x2i)) then
      if (first_time) then
         index_x2i_Faxa_lwdn   = mct_aVect_indexRA(x2i_i,'Faxa_lwdn') 
         index_x2i_Faxa_rain   = mct_aVect_indexRA(x2i_i,'Faxa_rain')  
         index_x2i_Faxa_snow   = mct_aVect_indexRA(x2i_i,'Faxa_snow')  
         index_x2i_Fioo_q      = mct_aVect_indexRA(x2i_i,'Fioo_q')  
         index_x2i_Fixx_rofi   = mct_aVect_indexRA(x2i_i,'Fixx_rofi')
      end if

      lSize = mct_avect_lSize(x2i_i)
      do n=1,lSize
         if (dom_i%data%rAttr(kLat,n) > 0.0_r8) then
            ic = c_inh_is
         else
            ic = c_ish_is
         endif
         do =  dom_i%data%rAttr(kArea,n) * frac_i%rAttr(ko,n)
         di =  dom_i%data%rAttr(kArea,n) * frac_i%rAttr(ki,n)
         if = f_area ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + di
         if = f_hlwdn; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + di*x2i_i%rAttr(index_x2i_Faxa_lwdn,n)
         if = f_wrain; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + di*x2i_i%rAttr(index_x2i_Faxa_rain,n)
         if = f_wsnow; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + di*x2i_i%rAttr(index_x2i_Faxa_snow,n)
         if = f_wioff; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) + di*x2i_i%rAttr(index_x2i_Fixx_rofi,n)
         if = f_hfrz ; budg_dataL(if,ic,ip) = budg_dataL(if,ic,ip) - (do+di)*max(0.0_r8,x2i_i%rAttr(index_x2i_Fioo_q,n))
      end do
      ic = c_inh_is  
      budg_dataL(f_hlatf,ic,ip) = -budg_dataL(f_wsnow,ic,ip)*shr_const_latice
      budg_dataL(f_hioff,ic,ip) = -budg_dataL(f_wioff,ic,ip)*shr_const_latice
      budg_dataL(f_wfrz ,ic,ip) =  budg_dataL(f_hfrz ,ic,ip)*HFLXtoWFLX
      ic = c_ish_is
      budg_dataL(f_hlatf,ic,ip) = -budg_dataL(f_wsnow,ic,ip)*shr_const_latice
      budg_dataL(f_hioff,ic,ip) = -budg_dataL(f_wioff,ic,ip)*shr_const_latice
      budg_dataL(f_wfrz ,ic,ip) =  budg_dataL(f_hfrz ,ic,ip)*HFLXtoWFLX
   end if

   first_time = .false.

end subroutine seq_diag_ice_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: seq_diag_print_mct - print global budget diagnostics
!
! !DESCRIPTION:
!   Print global budget diagnostics.
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE seq_diag_print_mct(EClock, stop_alarm, &
     budg_print_inst,  budg_print_daily,  budg_print_month,  &
     budg_print_ann,  budg_print_ltann,  budg_print_ltend)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(ESMF_Clock) , intent(in) :: EClock
   logical          , intent(in) :: stop_alarm
   integer          , intent(in) :: budg_print_inst
   integer          , intent(in) :: budg_print_daily
   integer          , intent(in) :: budg_print_month
   integer          , intent(in) :: budg_print_ann
   integer          , intent(in) :: budg_print_ltann
   integer          , intent(in) :: budg_print_ltend

!EOP

   !--- local ---
   integer(in)      :: ic,if,ip    ! data array indicies
   integer(in)      :: ica,icl,icn,ics,ico
   integer(in)      :: icar,icxs,icxr,icas
   integer(in)      :: n           ! loop counter
   integer(in)      :: nday        ! number of days in time avg
   integer(in)      :: cdate,sec   ! coded date, seconds
   integer(in)      :: yr,mon,day  ! date
   integer(in)      :: iam         ! pe number
   integer(in)      :: plev        ! print level
   logical          :: sumdone     ! has a sum been computed yet
   character(len=40):: str         ! string
   real(r8) :: dataGpr (f_size,c_size,p_size) ! values to print, scaled and such

   !----- formats -----
   character(*),parameter :: subName = '(seq_diag_print_mct) '
   character(*),parameter :: F00   = "('(seq_diag_print_mct) ',4a)"

   !----- formats -----
   character(*),parameter :: FAH="(4a,i9,i6)"
   character(*),parameter :: FA0= "('    ',8x,6(6x,a8,1x))"
   character(*),parameter :: FA1= "('    ',a8,6f15.8)"
   character(*),parameter :: FA0r="('    ',8x,8(6x,a8,1x))"
   character(*),parameter :: FA1r="('    ',a8,8f15.8)"

!-------------------------------------------------------------------------------
! print instantaneous budget data
!-------------------------------------------------------------------------------

   sumdone = .false.
   call seq_comm_setptrs(CPLID,iam=iam)
   call seq_timemgr_EClockGetData(EClock,curr_yr=yr, &
      curr_mon=mon,curr_day=day,curr_tod=sec)
   cdate = yr*10000+mon*100+day

   do ip = 1,p_size
      plev = 0
      if (ip == p_inst) then
         plev = max(plev,budg_print_inst)
      endif
      if (ip==p_day .and. sec==0) then
         plev = max(plev,budg_print_daily)
      endif
      if (ip==p_mon .and. day==1 .and. sec==0) then
         plev = max(plev,budg_print_month)
      endif
      if (ip==p_ann .and. mon==1 .and. day==1 .and. sec==0) then
         plev = max(plev,budg_print_ann)
      endif
      if (ip==p_inf .and. mon==1 .and. day==1 .and. sec==0) then
         plev = max(plev,budg_print_ltann)
      endif
      if (ip==p_inf .and. stop_alarm) then
         plev = max(plev,budg_print_ltend)
      endif

   if (plev > 0) then
! ---- doprint ---- doprint ---- doprint ----

   if (.not.sumdone) then
      call seq_diag_sum0_mct()
      dataGpr = budg_dataG
      sumdone = .true.

   !  old budget normalizations (global area and 1e6 for water)
      dataGpr = dataGpr/(4.0_r8*shr_const_pi)
      dataGpr(f_w:f_size,:,:) = dataGpr(f_w:f_size,:,:) * 1.0e6_r8
      dataGpr = dataGpr/budg_ns

      if (iam /= 0) return
   endif

   ! ---------------------------------------------------------
   ! ---- detail atm budgets and breakdown into components ---
   ! ---------------------------------------------------------

   if (plev >= 3) then
   do ic = 1,2
      if (ic == 1) then
         ica = c_atm_ar
         icl = c_lnd_ar
         icn = c_inh_ar
         ics = c_ish_ar
         ico = c_ocn_ar
         str = "ATM_to_CPL"
      elseif (ic == 2) then
         ica = c_atm_as
         icl = c_lnd_as
         icn = c_inh_as
         ics = c_ish_as
         ico = c_ocn_as
         str = "CPL_TO_ATM"
      else
         call shr_sys_abort(subname//' ERROR in ic index code 411')
      endif

      write(logunit,*) ' '
      write(logunit,FAH) subname,trim(str)//' AREA BUDGET (m2/m2): period = ',trim(pname(ip)),': date = ',cdate,sec
      write(logunit,FA0) cname(ica),cname(icl),cname(icn),cname(ics),cname(ico),' *SUM*  '
      do if = f_a, f_h-1
         write(logunit,FA1)    fname(if),dataGpr(if,ica,ip),dataGpr(if,icl,ip), &
                   dataGpr(if,icn,ip),dataGpr(if,ics,ip),dataGpr(if,ico,ip), &
                                         dataGpr(if,ica,ip)+dataGpr(if,icl,ip)+ &
                   dataGpr(if,icn,ip)+dataGpr(if,ics,ip)+dataGpr(if,ico,ip) 
      enddo

      write(logunit,*) ' '
      write(logunit,FAH) subname,trim(str)//' HEAT BUDGET (W/m2): period = ',trim(pname(ip)),': date = ',cdate,sec
      write(logunit,FA0) cname(ica),cname(icl),cname(icn),cname(ics),cname(ico),' *SUM*  '
      do if = f_h, f_w-1
         write(logunit,FA1)    fname(if),dataGpr(if,ica,ip),dataGpr(if,icl,ip), &
                   dataGpr(if,icn,ip),dataGpr(if,ics,ip),dataGpr(if,ico,ip), &
                                         dataGpr(if,ica,ip)+dataGpr(if,icl,ip)+ &
                   dataGpr(if,icn,ip)+dataGpr(if,ics,ip)+dataGpr(if,ico,ip) 
      enddo
      write(logunit,FA1)    '   *SUM*',sum(dataGpr(f_h:f_w-1,ica,ip)),sum(dataGpr(f_h:f_w-1,icl,ip)), &
         sum(dataGpr(f_h:f_w-1,icn,ip)),sum(dataGpr(f_h:f_w-1,ics,ip)),sum(dataGpr(f_h:f_w-1,ico,ip)), &
                                       sum(dataGpr(f_h:f_w-1,ica,ip))+sum(dataGpr(f_h:f_w-1,icl,ip))+ &
         sum(dataGpr(f_h:f_w-1,icn,ip))+sum(dataGpr(f_h:f_w-1,ics,ip))+sum(dataGpr(f_h:f_w-1,ico,ip)) 

      write(logunit,*) ' '
      write(logunit,FAH) subname,trim(str)//' WATER BUDGET (kg/m2s*1e6): period = ',trim(pname(ip)),': date = ',cdate,sec
      write(logunit,FA0) cname(ica),cname(icl),cname(icn),cname(ics),cname(ico),' *SUM*  '
      do if = f_w, f_size
         write(logunit,FA1)    fname(if),dataGpr(if,ica,ip),dataGpr(if,icl,ip), &
                   dataGpr(if,icn,ip),dataGpr(if,ics,ip),dataGpr(if,ico,ip), &
                                         dataGpr(if,ica,ip)+dataGpr(if,icl,ip)+ &
                   dataGpr(if,icn,ip)+dataGpr(if,ics,ip)+dataGpr(if,ico,ip) 
      enddo
      write(logunit,FA1)    '   *SUM*',sum(dataGpr(f_w:f_size,ica,ip)),sum(dataGpr(f_w:f_size,icl,ip)), &
         sum(dataGpr(f_w:f_size,icn,ip)),sum(dataGpr(f_w:f_size,ics,ip)),sum(dataGpr(f_w:f_size,ico,ip)), &
                                       sum(dataGpr(f_w:f_size,ica,ip))+sum(dataGpr(f_w:f_size,icl,ip))+ &
         sum(dataGpr(f_w:f_size,icn,ip))+sum(dataGpr(f_w:f_size,ics,ip))+sum(dataGpr(f_w:f_size,ico,ip)) 
   enddo
   endif   ! plev

   ! ---------------------------------------------------------
   ! ---- detail lnd/ocn/ice component budgets ----
   ! ---------------------------------------------------------

   if (plev >= 2) then
   do ic = 1,4
      if (ic == 1) then
         icar = c_lnd_ar
         icxs = c_lnd_ls
         icxr = c_lnd_lr
         icas = c_lnd_as
         str = "LND"
      elseif (ic == 2) then
         icar = c_ocn_ar
         icxs = c_ocn_os
         icxr = c_ocn_or
         icas = c_ocn_as
         str = "OCN"
      elseif (ic == 3) then
         icar = c_inh_ar
         icxs = c_inh_is
         icxr = c_inh_ir
         icas = c_inh_as
         str = "ICE_NH"
      elseif (ic == 4) then
         icar = c_ish_ar
         icxs = c_ish_is
         icxr = c_ish_ir
         icas = c_ish_as
         str = "ICE_SH"
      else
         call shr_sys_abort(subname//' ERROR in ic index code 412')
      endif

      write(logunit,*) ' '
      write(logunit,FAH) subname,trim(str)//' HEAT BUDGET (W/m2): period = ',trim(pname(ip)),': date = ',cdate,sec
      write(logunit,FA0) cname(icar),cname(icxs),cname(icxr),cname(icas),' *SUM*  '
      do if = f_h, f_w-1
         write(logunit,FA1)    fname(if),-dataGpr(if,icar,ip),dataGpr(if,icxs,ip), &
                                          dataGpr(if,icxr,ip),-dataGpr(if,icas,ip), &
                                         -dataGpr(if,icar,ip)+dataGpr(if,icxs,ip)+ &
                                          dataGpr(if,icxr,ip)-dataGpr(if,icas,ip)
      enddo
      write(logunit,FA1)    '   *SUM*',-sum(dataGpr(f_h:f_w-1,icar,ip)),sum(dataGpr(f_h:f_w-1,icxs,ip)), &
                                       sum(dataGpr(f_h:f_w-1,icxr,ip)),-sum(dataGpr(f_h:f_w-1,icas,ip)), &
                                       -sum(dataGpr(f_h:f_w-1,icar,ip))+sum(dataGpr(f_h:f_w-1,icxs,ip))+ &
                                       sum(dataGpr(f_h:f_w-1,icxr,ip))-sum(dataGpr(f_h:f_w-1,icas,ip))

      write(logunit,*) ' '
      write(logunit,FAH) subname,trim(str)//' WATER BUDGET (kg/m2s*1e6): period = ',trim(pname(ip)),': date = ',cdate,sec
      write(logunit,FA0) cname(icar),cname(icxs),cname(icxr),cname(icas),' *SUM*  '
      do if = f_w, f_size
         write(logunit,FA1)    fname(if),-dataGpr(if,icar,ip),dataGpr(if,icxs,ip), &
                                         dataGpr(if,icxr,ip),-dataGpr(if,icas,ip), &
                                         -dataGpr(if,icar,ip)+dataGpr(if,icxs,ip)+ &
                                         dataGpr(if,icxr,ip)-dataGpr(if,icas,ip)
      enddo
      write(logunit,FA1)    '   *SUM*',-sum(dataGpr(f_w:f_size,icar,ip)),sum(dataGpr(f_w:f_size,icxs,ip)), &
                                       sum(dataGpr(f_w:f_size,icxr,ip)),-sum(dataGpr(f_w:f_size,icas,ip)), &
                                       -sum(dataGpr(f_w:f_size,icar,ip))+sum(dataGpr(f_w:f_size,icxs,ip))+ &
                                       sum(dataGpr(f_w:f_size,icxr,ip))-sum(dataGpr(f_w:f_size,icas,ip))

   enddo
   endif   ! plev

   ! ---------------------------------------------------------
   ! ---- net summary budgets ----
   ! ---------------------------------------------------------

   if (plev >= 1) then

      write(logunit,*) ' '
      write(logunit,FAH) subname,'NET AREA BUDGET (m2/m2): period = ',trim(pname(ip)),': date = ',cdate,sec
      write(logunit,FA0) '     atm','     lnd','     ocn','  ice nh','  ice sh',' *SUM*  '
      do if = 1,f_h-1
         write(logunit,FA1)    fname(if),dataGpr(if,c_atm_ar,ip), &
                                         dataGpr(if,c_lnd_lr,ip), &
                                         dataGpr(if,c_ocn_or,ip), &
                                         dataGpr(if,c_inh_ir,ip), &
                                         dataGpr(if,c_ish_ir,ip), &
                                         dataGpr(if,c_atm_ar,ip)+ &
                                         dataGpr(if,c_lnd_lr,ip)+ &
                                         dataGpr(if,c_ocn_or,ip)+ &
                                         dataGpr(if,c_inh_ir,ip)+ &
                                         dataGpr(if,c_ish_ir,ip)
      enddo

      write(logunit,*) ' '
      write(logunit,FAH) subname,'NET HEAT BUDGET (W/m2): period = ',trim(pname(ip)),': date = ',cdate,sec
      write(logunit,FA0r) '     atm','     lnd','     rof','     ocn','  ice nh','  ice sh','     glc',' *SUM*  '
      do if = f_h, f_w-1
         write(logunit,FA1r)   fname(if),dataGpr(if,c_atm_ar,ip)+dataGpr(if,c_atm_as,ip), &
                                         dataGpr(if,c_lnd_lr,ip)+dataGpr(if,c_lnd_ls,ip), &
                                         dataGpr(if,c_rof_rr,ip)+dataGpr(if,c_rof_rs,ip), &
                                         dataGpr(if,c_ocn_or,ip)+dataGpr(if,c_ocn_os,ip), &
                                         dataGpr(if,c_inh_ir,ip)+dataGpr(if,c_inh_is,ip), &
                                         dataGpr(if,c_ish_ir,ip)+dataGpr(if,c_ish_is,ip), &
                                         dataGpr(if,c_glc_gr,ip)+dataGpr(if,c_glc_gs,ip), &
                                         dataGpr(if,c_atm_ar,ip)+dataGpr(if,c_atm_as,ip)+ &
                                         dataGpr(if,c_lnd_lr,ip)+dataGpr(if,c_lnd_ls,ip)+ &
                                         dataGpr(if,c_rof_rr,ip)+dataGpr(if,c_rof_rs,ip)+ &
                                         dataGpr(if,c_ocn_or,ip)+dataGpr(if,c_ocn_os,ip)+ &
                                         dataGpr(if,c_inh_ir,ip)+dataGpr(if,c_inh_is,ip)+ &
                                         dataGpr(if,c_ish_ir,ip)+dataGpr(if,c_ish_is,ip)+ &
                                         dataGpr(if,c_glc_gr,ip)+dataGpr(if,c_glc_gs,ip)
      enddo
      write(logunit,FA1r)'   *SUM*',sum(dataGpr(f_h:f_w-1,c_atm_ar,ip))+sum(dataGpr(f_h:f_w-1,c_atm_as,ip)), &
                                    sum(dataGpr(f_h:f_w-1,c_lnd_lr,ip))+sum(dataGpr(f_h:f_w-1,c_lnd_ls,ip)), &
                                    sum(dataGpr(f_h:f_w-1,c_rof_rr,ip))+sum(dataGpr(f_h:f_w-1,c_rof_rs,ip)), &
                                    sum(dataGpr(f_h:f_w-1,c_ocn_or,ip))+sum(dataGpr(f_h:f_w-1,c_ocn_os,ip)), &
                                    sum(dataGpr(f_h:f_w-1,c_inh_ir,ip))+sum(dataGpr(f_h:f_w-1,c_inh_is,ip)), &
                                    sum(dataGpr(f_h:f_w-1,c_ish_ir,ip))+sum(dataGpr(f_h:f_w-1,c_ish_is,ip)), &
                                    sum(dataGpr(f_h:f_w-1,c_glc_gr,ip))+sum(dataGpr(f_h:f_w-1,c_glc_gs,ip)), &
                                    sum(dataGpr(f_h:f_w-1,c_atm_ar,ip))+sum(dataGpr(f_h:f_w-1,c_atm_as,ip))+ &
                                    sum(dataGpr(f_h:f_w-1,c_lnd_lr,ip))+sum(dataGpr(f_h:f_w-1,c_lnd_ls,ip))+ &
                                    sum(dataGpr(f_h:f_w-1,c_rof_rr,ip))+sum(dataGpr(f_h:f_w-1,c_rof_rs,ip))+ &
                                    sum(dataGpr(f_h:f_w-1,c_ocn_or,ip))+sum(dataGpr(f_h:f_w-1,c_ocn_os,ip))+ &
                                    sum(dataGpr(f_h:f_w-1,c_inh_ir,ip))+sum(dataGpr(f_h:f_w-1,c_inh_is,ip))+ &
                                    sum(dataGpr(f_h:f_w-1,c_ish_ir,ip))+sum(dataGpr(f_h:f_w-1,c_ish_is,ip))+ &
                                    sum(dataGpr(f_h:f_w-1,c_glc_gr,ip))+sum(dataGpr(f_h:f_w-1,c_glc_gs,ip))

      write(logunit,*) ' '
      write(logunit,FAH) subname,'NET WATER BUDGET (kg/m2s*1e6): period = ',trim(pname(ip)),': date = ',cdate,sec
      write(logunit,FA0r) '     atm','     lnd','     rof','     ocn','  ice nh','  ice sh','     glc',' *SUM*  '
      do if = f_w, f_size
         write(logunit,FA1r)   fname(if),dataGpr(if,c_atm_ar,ip)+dataGpr(if,c_atm_as,ip), &
                                         dataGpr(if,c_lnd_lr,ip)+dataGpr(if,c_lnd_ls,ip), &
                                         dataGpr(if,c_rof_rr,ip)+dataGpr(if,c_rof_rs,ip), &
                                         dataGpr(if,c_ocn_or,ip)+dataGpr(if,c_ocn_os,ip), &
                                         dataGpr(if,c_inh_ir,ip)+dataGpr(if,c_inh_is,ip), &
                                         dataGpr(if,c_ish_ir,ip)+dataGpr(if,c_ish_is,ip), &
                                         dataGpr(if,c_glc_gr,ip)+dataGpr(if,c_glc_gs,ip), &
                                         dataGpr(if,c_atm_ar,ip)+dataGpr(if,c_atm_as,ip)+ &
                                         dataGpr(if,c_lnd_lr,ip)+dataGpr(if,c_lnd_ls,ip)+ &
                                         dataGpr(if,c_rof_rr,ip)+dataGpr(if,c_rof_rs,ip)+ &
                                         dataGpr(if,c_ocn_or,ip)+dataGpr(if,c_ocn_os,ip)+ &
                                         dataGpr(if,c_inh_ir,ip)+dataGpr(if,c_inh_is,ip)+ &
                                         dataGpr(if,c_ish_ir,ip)+dataGpr(if,c_ish_is,ip)+ &
                                         dataGpr(if,c_glc_gr,ip)+dataGpr(if,c_glc_gs,ip)
      enddo
      write(logunit,FA1r)'   *SUM*',sum(dataGpr(f_w:f_size,c_atm_ar,ip))+sum(dataGpr(f_w:f_size,c_atm_as,ip)), &
                                    sum(dataGpr(f_w:f_size,c_lnd_lr,ip))+sum(dataGpr(f_w:f_size,c_lnd_ls,ip)), &
                                    sum(dataGpr(f_w:f_size,c_rof_rr,ip))+sum(dataGpr(f_w:f_size,c_rof_rs,ip)), &
                                    sum(dataGpr(f_w:f_size,c_ocn_or,ip))+sum(dataGpr(f_w:f_size,c_ocn_os,ip)), &
                                    sum(dataGpr(f_w:f_size,c_inh_ir,ip))+sum(dataGpr(f_w:f_size,c_inh_is,ip)), &
                                    sum(dataGpr(f_w:f_size,c_ish_ir,ip))+sum(dataGpr(f_w:f_size,c_ish_is,ip)), &
                                    sum(dataGpr(f_w:f_size,c_glc_gr,ip))+sum(dataGpr(f_w:f_size,c_glc_gs,ip)), &
                                    sum(dataGpr(f_w:f_size,c_atm_ar,ip))+sum(dataGpr(f_w:f_size,c_atm_as,ip))+ &
                                    sum(dataGpr(f_w:f_size,c_lnd_lr,ip))+sum(dataGpr(f_w:f_size,c_lnd_ls,ip))+ &
                                    sum(dataGpr(f_w:f_size,c_rof_rr,ip))+sum(dataGpr(f_w:f_size,c_rof_rs,ip))+ &
                                    sum(dataGpr(f_w:f_size,c_ocn_or,ip))+sum(dataGpr(f_w:f_size,c_ocn_os,ip))+ &
                                    sum(dataGpr(f_w:f_size,c_inh_ir,ip))+sum(dataGpr(f_w:f_size,c_inh_is,ip))+ &
                                    sum(dataGpr(f_w:f_size,c_ish_ir,ip))+sum(dataGpr(f_w:f_size,c_ish_is,ip))+ &
                                    sum(dataGpr(f_w:f_size,c_glc_gr,ip))+sum(dataGpr(f_w:f_size,c_glc_gs,ip))

   endif

   write(logunit,*) ' '
! ---- doprint ---- doprint ---- doprint ----
   endif  ! plev > 0
   enddo  ! ip = 1,p_size

end subroutine seq_diag_print_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: seq_diag_avect_mct - print global budget diagnostics
!
! !DESCRIPTION:
!   Print global diagnostics for AV/ID.
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE seq_diag_avect_mct(infodata, id, av, dom, gsmap, comment)

   use seq_infodata_mod

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(seq_infodata_type) , intent(in)           :: infodata
   integer(in)             , intent(in)           :: ID
   type(mct_aVect)         , intent(in)           :: av
   type(mct_gGrid)         , pointer              :: dom
   type(mct_gsMap)         , pointer              :: gsmap
   character(len=*)        , intent(in), optional :: comment

!EOP

   !--- local ---
   logical                          :: bfbflag
   integer(in)                      :: n,k         ! counters
   integer(in)                      :: npts,nptsg  ! number of local/global pts in AV
   integer(in)                      :: kflds       ! number of fields in AV
   real(r8),                pointer :: sumbuf (:)  ! sum buffer
   real(r8),                pointer :: minbuf (:)  ! min buffer
   real(r8),                pointer :: maxbuf (:)  ! max buffer
   real(r8),                pointer :: sumbufg(:)  ! sum buffer reduced
   real(r8),                pointer :: minbufg(:)  ! min buffer reduced
   real(r8),                pointer :: maxbufg(:)  ! max buffer reduced
   integer(i8),             pointer :: isumbuf (:) ! integer local sum
   integer(i8),             pointer :: isumbufg(:) ! integer global sum
   integer(i8)                      :: ihuge       ! huge
   integer(in)                      :: mpicom      ! mpi comm
   integer(in)                      :: iam         ! pe number
   integer(in)                      :: km,ka       ! field indices
   integer(in)                      :: ns          ! size of local AV
   integer(in)                      :: rcode       ! error code
   real(r8),                pointer :: weight(:)   ! weight
   type(mct_string)                 :: mstring     ! mct char type
   character(CL)                    :: lcomment    ! should be long enough
   character(CL)                    :: itemc       ! string converted to char

   type(mct_avect)                  :: AV1         ! local avect with one field
   type(mct_avect)                  :: AVr1        ! avect on root with one field
   type(mct_avect)                  :: AVr2        ! avect on root with one field

   !----- formats -----
   character(*),parameter :: subName = '(seq_diag_avect_mct) '
   character(*),parameter :: F00   = "('(seq_diag_avect_mct) ',4a)"

!-------------------------------------------------------------------------------
! print instantaneous budget data
!-------------------------------------------------------------------------------

   call seq_comm_setptrs(ID,&
        mpicom=mpicom, iam=iam)

   call seq_infodata_GetData(infodata,&
        bfbflag=bfbflag)

   lcomment = ''
   if (present(comment)) then
      lcomment=trim(comment)
   endif

   ns = mct_aVect_lsize(AV)
   npts = mct_aVect_lsize(dom%data)
   if (ns /= npts) call shr_sys_abort(trim(subname)//' ERROR: size of AV,dom')
   km = mct_aVect_indexRA(dom%data,'mask')
   ka = mct_aVect_indexRA(dom%data,afldname)
   kflds = mct_aVect_nRattr(AV)
   allocate(sumbuf(kflds),sumbufg(kflds))

   sumbuf =       0.0_r8

   if (bfbflag) then

      npts = mct_aVect_lsize(AV)
      allocate(weight(npts))
      weight(:) = 1.0_r8 
      do n = 1,npts
         if (dom%data%rAttr(km,n) <= 1.0e-06_R8) then
            weight(n) = 0.0_r8
         else
            weight(n) = dom%data%rAttr(ka,n)*shr_const_rearth*shr_const_rearth
         endif
      enddo

      allocate(maxbuf(kflds),maxbufg(kflds))
      maxbuf = 0.0_r8

      do n = 1,npts
      do k = 1,kflds
         if (AV%rAttr(k,n) > 1.01_r8*shr_const_spval .or. &
             AV%rAttr(k,n) < 0.99_r8*shr_const_spval) then
             maxbuf(k) = max(maxbuf(k),abs(AV%rAttr(k,n)*weight(n)))
         endif
      enddo
      enddo

      call shr_mpi_max(maxbuf,maxbufg,mpicom,subname,all=.true.)
      call shr_mpi_sum(npts,nptsg,mpicom,subname,all=.true.)

      do k = 1,kflds
         if (maxbufg(k) < 1000.0*TINY(maxbufg(k)) .or. &
             maxbufg(k) > HUGE(maxbufg(k))/(2.0_r8*nptsg)) then
            maxbufg(k) = 0.0_r8
         else
            maxbufg(k) = (1.1_r8) * maxbufg(k) * nptsg
         endif
      enddo

      allocate(isumbuf(kflds),isumbufg(kflds))
      isumbuf = 0
      ihuge = HUGE(isumbuf)

      do n = 1,npts
      do k = 1,kflds
         if (AV%rAttr(k,n) > 1.01_r8*shr_const_spval .or. &
             AV%rAttr(k,n) < 0.99_r8*shr_const_spval) then
             if (abs(maxbufg(k)) > 1000.0_r8 * TINY(maxbufg)) then
                isumbuf(k) = isumbuf(k) + int((AV%rAttr(k,n)*weight(n)/maxbufg(k))*ihuge,i8)
             endif
         endif
      enddo
      enddo

      call shr_mpi_sum(isumbuf,isumbufg,mpicom,subname)

      do k = 1,kflds
         sumbufg(k) = isumbufg(k)*maxbufg(k)/ihuge
      enddo

      deallocate(weight)
      deallocate(maxbuf,maxbufg)
      deallocate(isumbuf,isumbufg)

   else

      npts = mct_aVect_lsize(AV)
      allocate(weight(npts))
      weight(:) = 1.0_r8 
      do n = 1,npts
         if (dom%data%rAttr(km,n) <= 1.0e-06_R8) then
            weight(n) = 0.0_r8
         else
            weight(n) = dom%data%rAttr(ka,n)*shr_const_rearth*shr_const_rearth
         endif
      enddo

      do n = 1,npts
      do k = 1,kflds
         if (AV%rAttr(k,n) > 1.01_r8*shr_const_spval .or. &
             AV%rAttr(k,n) < 0.99_r8*shr_const_spval) then
             sumbuf(k) = sumbuf(k) + AV%rAttr(k,n)*weight(n)
         endif
      enddo
      enddo

      !--- global reduction ---
      call shr_mpi_sum(sumbuf,sumbufg,mpicom,subname)

      deallocate(weight)

   endif

   if (iam == 0) then
      !      write(logunit,*) 'sdAV: *** writing ',trim(lcomment),': k fld min/max/sum ***'
      do k = 1,kflds
         call mct_aVect_getRList(mstring,k,AV)
         itemc = mct_string_toChar(mstring)
         call mct_string_clean(mstring)
         if (len_trim(lcomment) > 0) then
            write(logunit,100) 'xxx','sorr',k,sumbufg(k),trim(lcomment),trim(itemc)
         else
            write(logunit,101) 'xxx','sorr',k,sumbufg(k),trim(itemc)
         endif
      enddo
      call shr_sys_flush(logunit)
   endif

   deallocate(sumbuf,sumbufg)

100  format('comm_diag ',a3,1x,a4,1x,i3,es26.19,1x,a,1x,a)
101  format('comm_diag ',a3,1x,a4,1x,i3,es26.19,1x,a)

end subroutine seq_diag_avect_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: seq_diag_avdiff_mct - print global budget diagnostics
!
! !DESCRIPTION:
!   Print global diagnostics for AV/ID.
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE seq_diag_avdiff_mct(AV1,AV2,ID,comment)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(mct_aVect) , intent(in) :: AV1
   type(mct_aVect) , intent(in) :: AV2
   integer         , intent(in) :: ID
   character(len=*), intent(in), optional :: comment

!EOP

   !--- local ---
   integer(in)      :: n,k,n1,k1,n2,k2         ! counters
   integer(in)      :: iam         ! pe number
   integer(in)      :: cnt         ! counter
   real(r8)         :: adiff,rdiff ! diff values
   type(mct_string) :: mstring     ! mct char type
   character(len=64):: lcomment    ! should be long enough

   !----- formats -----
   character(*),parameter :: subName = '(seq_diag_avdiff_mct) '
   character(*),parameter :: F00   = "('(seq_diag_avdiff_mct) ',4a)"

!-------------------------------------------------------------------------------
! print instantaneous budget data
!-------------------------------------------------------------------------------

   call seq_comm_setptrs(ID,iam=iam)

   lcomment = ''
   if (present(comment)) then
      lcomment=trim(comment)
   endif

   n1 = mct_aVect_lsize(AV1)
   k1 = mct_aVect_nRattr(AV1)
   n2 = mct_aVect_lsize(AV2)
   k2 = mct_aVect_nRattr(AV2)

   if (n1 /= n2 .or. k1 /= k2) then
      write(s_logunit,*) subname,trim(lcomment),' AV sizes different ',n1,n2,k1,k2
      return
   endif

   do k = 1,k1
      cnt = 0
      adiff = 0.
      rdiff = 0.
      do n = 1,n1
         if (AV1%rAttr(k,n) /= AV2%rAttr(k,n)) then
            cnt = cnt + 1
            adiff = max(adiff, abs(AV1%rAttr(k,n)-AV2%rAttr(k,n)))
            rdiff = max(rdiff, abs(AV1%rAttr(k,n)-AV2%rAttr(k,n))/(abs(AV1%rAttr(k,n))+abs(AV2%rAttr(k,n))))
         endif
      enddo
      if (cnt > 0) then
         call mct_aVect_getRList(mstring,k,AV1)
         write(s_logunit,*) subname,trim(lcomment),' AVs fld k diff ', &
            iam,mct_string_toChar(mstring),cnt,adiff,rdiff, &
            minval(AV1%rAttr(k,:)),minval(AV1%rAttr(k,:)), &
            maxval(AV1%rAttr(k,:)),maxval(AV2%rAttr(k,:))
         call mct_string_clean(mstring)       
      endif
   enddo

end subroutine seq_diag_avdiff_mct

!===============================================================================
end module seq_diag_mct
