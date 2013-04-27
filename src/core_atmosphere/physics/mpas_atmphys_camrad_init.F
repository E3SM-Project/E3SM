! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
#define DM_BCAST_MACRO(A) call mpas_dmpar_bcast_reals(dminfo,size(A),A)

!==================================================================================================
 module mpas_atmphys_camrad_init
 use mpas_dmpar
 use mpas_kind_types
 use mpas_grid_types
 
 use mpas_atmphys_constants,only: cp,degrad,ep_2,g,R_d,R_v,stbolt
 use mpas_atmphys_utilities

!wrf physics:
 use module_ra_cam_support

 implicit none
 private
 public:: camradinit


!>\brief Initialization of CAM radiation codes using MPAS MPI decomposition.
!>\author Laura D. Fowler (send comments to laura@ucar.edu).
!>\date 2013-05-01.
!>
!>\details
!> subroutine camradinit calls the main subroutines needed to initialize the long- and short-wave
!> CAM radiation codes, and read input data from auxillary files.
!>
!> subroutines called in mpas_atmphys_camrad_init:
!> -----------------------------------------------
!> radini      :initialization of radiation constants.
!> esinti      :initialization of saturation vapor pressures.
!> oznini      :initialization of climatological monthly-mean ozone profiles.
!> aerosol_init:initialization of aerosol optical properties.
!>
!> add-ons and modifications to sourcecode:
!> ----------------------------------------
!>    * added initialization of variable mxaerl which is the number of layers below 900 hPa in
!>      which background aerosols are present. mxaerl is computed using the pressure-base array.
!>      -> added diag in the argument list of subroutines camradinit and aerosol_init.
!>      -> in subroutine aerosol_init, added initialization of variable mxaerl.
!>      Laura D. Fowler (birch.ucar.edu) / 2013-07-01.
!>    * moved the arrays pin and ozmixm from the mesh structure to the atm_input structure in
!>      subroutine oznini.
!>      Laura D. Fowler (birch.ucar.edu) / 2013-07-08.

!local parameters:
 integer,parameter:: latsiz = 64
 integer,parameter:: lonsiz = 1

 contains

!==================================================================================================
 subroutine camradinit(dminfo,mesh,atm_input,diag,state_1,state_2)
!==================================================================================================

!input arguments:
 type(dm_info),intent(in):: dminfo
 type(mesh_type),intent(in):: mesh
 type(diag_type),intent(in):: diag

!inout arguments:
 type(atm_input_type),intent(inout):: atm_input
 type(state_type),intent(inout):: state_1,state_2 

!local variables:
 real(r8):: pstd
 real(r8):: rh2o, cpair

!--------------------------------------------------------------------------------------------------

!...these were made allocatable 20090612 to save static memory allocation. JM:
 if ( .not. allocated( ksul   ) ) allocate( ksul(nrh,nspint)   )
 if ( .not. allocated( wsul   ) ) allocate( wsul(nrh,nspint)   )
 if ( .not. allocated( gsul   ) ) allocate( gsul(nrh,nspint)   )
 if ( .not. allocated( ksslt  ) ) allocate( ksslt(nrh,nspint)  )
 if ( .not. allocated( wsslt  ) ) allocate( wsslt(nrh,nspint)  )
 if ( .not. allocated( gsslt  ) ) allocate( gsslt(nrh,nspint)  )
 if ( .not. allocated( kcphil ) ) allocate( kcphil(nrh,nspint) )
 if ( .not. allocated( wcphil ) ) allocate( wcphil(nrh,nspint) )
 if ( .not. allocated( gcphil ) ) allocate( gcphil(nrh,nspint) )

 if (.not. allocated(ah2onw  ) ) allocate( ah2onw(n_p,n_tp,n_u,n_te,n_rh)   )
 if (.not. allocated(eh2onw  ) ) allocate( eh2onw(n_p,n_tp,n_u,n_te,n_rh)   )
 if (.not. allocated(ah2ow   ) ) allocate( ah2ow (n_p,n_tp,n_u,n_te,n_rh)   )
 if (.not. allocated(cn_ah2ow) ) allocate( cn_ah2ow(n_p,n_tp,n_u,n_te,n_rh) )
 if (.not. allocated(cn_eh2ow) ) allocate( cn_eh2ow(n_p,n_tp,n_u,n_te,n_rh) )
 if (.not. allocated(ln_ah2ow) ) allocate( ln_ah2ow(n_p,n_tp,n_u,n_te,n_rh) )
 if (.not. allocated(ln_eh2ow) ) allocate( ln_eh2ow(n_p,n_tp,n_u,n_te,n_rh) )

 ozncyc   = .true.
 indirect = .true.
 ixcldliq = 2
 ixcldice = 3

 pstd = 101325.0

!...from physconst module:
 mwdry = 28.966            ! molecular weight dry air ~ kg/kmole (shr_const_mwdair)
 mwco2 =  44.              ! molecular weight co2
 mwh2o = 18.016            ! molecular weight water vapor (shr_const_mwwv)
 mwch4 =  16.              ! molecular weight ch4
 mwn2o =  44.              ! molecular weight n2o
 mwf11 = 136.              ! molecular weight cfc11
 mwf12 = 120.              ! molecular weight cfc12
 cappa = R_D/CP
 rair = R_D
 tmelt = 273.16            ! freezing T of fresh water ~ K 
 r_universal = 6.02214e26 * stbolt   ! Universal gas constant ~ J/K/kmole
 latvap = 2.501e6          ! latent heat of evaporation ~ J/kg
 latice = 3.336e5          ! latent heat of fusion ~ J/kg
 zvir = R_V/R_D - 1.
 rh2o = R_V
 cpair = CP

 epsqs = EP_2

!initialization of some constants:
 call radini(dminfo,g,cp,ep_2,stbolt,pstd*10.0)
 write(0,*) '    end subroutine radini'

!initialization of saturation vapor pressures:
 call esinti(epsqs,latvap,latice,rh2o,cpair,tmelt)
 write(0,*) '    end subroutine esinti'
 
!initialization of ozone mixing ratios:
 call oznini(mesh,atm_input)
 write(0,*) '    end subroutine oznini'

!initialization of aerosol concentrations:
 call aerosol_init(dminfo,mesh,diag,state_1,state_2)
 write(0,*) '    end subroutine aerosol_init'

 end subroutine camradinit

!==================================================================================================
 subroutine radini(dminfo,gravx,cpairx,epsilox,stebolx,pstdx)
!--------------------------------------------------------------------------------------------------
! 
! Purpose: 
! Initialize various constants for radiation scheme; note that
! the radiation scheme uses cgs units.
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: W. Collins (H2O parameterization) and J. Kiehl
! 
!-----------------------------------------------------------------------
!  use shr_kind_mod, only: r8 => shr_kind_r8
!  use ppgrid,       only: pver, pverp
!  use comozp,       only: cplos, cplol
!  use pmgrid,       only: masterproc, plev, plevp
!  use radae,        only: radaeini
!  use physconst,    only: mwdry, mwco2
#if ( defined SPMD )
!   use mpishorthand
#endif
   implicit none

!------------------------------Arguments--------------------------------
!
! Input arguments
!
   type(dm_info),intent(in):: dminfo

   real, intent(in) :: gravx      ! Acceleration of gravity (MKS)
   real, intent(in) :: cpairx     ! Specific heat of dry air (MKS)
   real, intent(in) :: epsilox    ! Ratio of mol. wght of H2O to dry air
   real, intent(in) :: stebolx    ! Stefan-Boltzmann's constant (MKS)
   real(r8), intent(in) :: pstdx      ! Standard pressure (Pascals)
!
!---------------------------Local variables-----------------------------
!
   integer k       ! Loop variable

   real(r8) v0         ! Volume of a gas at stp (m**3/kmol)
   real(r8) p0         ! Standard pressure (pascals)
   real(r8) amd        ! Effective molecular weight of dry air (kg/kmol)
   real(r8) goz        ! Acceleration of gravity (m/s**2)
!
!-----------------------------------------------------------------------
!
! Set general radiation consts; convert to cgs units where appropriate:
!
   gravit  =  100.*gravx
   rga     =  1./gravit
   gravmks =  gravx
   cpair   =  1.e4*cpairx
   epsilo  =  epsilox
   sslp    =  1.013250e6
   stebol  =  1.e3*stebolx
   rgsslp  =  0.5/(gravit*sslp)
   dpfo3   =  2.5e-3
   dpfco2  =  5.0e-3
   dayspy  =  365.
   pie     =  4.*atan(1.)
!
! Initialize ozone data.
!
   v0  = 22.4136         ! Volume of a gas at stp (m**3/kmol)
   p0  = 0.1*sslp        ! Standard pressure (pascals)
   amd = 28.9644         ! Molecular weight of dry air (kg/kmol)
   goz = gravx           ! Acceleration of gravity (m/s**2)
!
! Constants for ozone path integrals (multiplication by 100 for unit
! conversion to cgs from mks):
!
   cplos = v0/(amd*goz)       *100.0
   cplol = v0/(amd*goz*p0)*0.5*100.0
!
! Derived constants
! If the top model level is above ~90 km (0.1 Pa), set the top level to compute
! longwave cooling to about 80 km (1 Pa)
! WRF: assume top level > 0.1 mb
!  if (hypm(1) .lt. 0.1) then
!     do k = 1, pver
!        if (hypm(k) .lt. 1.) ntoplw  = k
!     end do
!  else
      ntoplw = 1
!  end if
!   if (masterproc) then
!     write (6,*) 'RADINI: ntoplw =',ntoplw, ' pressure:',hypm(ntoplw)
!   endif

   call radaeini(dminfo,pstdx,mwdry,mwco2)

 end subroutine radini

!==================================================================================================
 subroutine radaeini(dminfo,pstdx,mwdryx,mwco2x)
!==================================================================================================

!input arguments:
 type(dm_info),intent(in):: dminfo

 real(r8), intent(in) :: pstdx   ! Standard pressure (dynes/cm^2)
 real(r8), intent(in) :: mwdryx  ! Molecular weight of dry air 
 real(r8), intent(in) :: mwco2x  ! Molecular weight of carbon dioxide

!local variables:

!variables for loading absorptivity/emissivity:
 integer:: ncid_ae               ! NetCDF file id for abs/ems file
 integer:: pdimid                ! pressure dimension id
 integer:: psize                 ! pressure dimension size
 integer:: tpdimid               ! path temperature dimension id
 integer:: tpsize                ! path temperature size
 integer:: tedimid               ! emission temperature dimension id
 integer:: tesize                ! emission temperature size

 integer:: udimid                ! u (H2O path) dimension id
 integer:: usize                 ! u (H2O path) dimension size

 integer:: rhdimid               ! relative humidity dimension id
 integer:: rhsize                ! relative humidity dimension size

 integer::    ah2onwid           ! var. id for non-wndw abs.
 integer::    eh2onwid           ! var. id for non-wndw ems.
 integer::    ah2owid            ! var. id for wndw abs. (adjacent layers)
 integer:: cn_ah2owid            ! var. id for continuum trans. for wndw abs.
 integer:: cn_eh2owid            ! var. id for continuum trans. for wndw ems.
 integer:: ln_ah2owid            ! var. id for line trans. for wndw abs.
 integer:: ln_eh2owid            ! var. id for line trans. for wndw ems.
   
!character*(NF_MAX_NAME) tmpname! dummy variable for var/dim names
 character(len=StrKIND):: locfn      ! local filename
 integer:: tmptype                 ! dummy variable for variable type
 integer:: ndims                   ! number of dimensions
!integer dims(NF_MAX_VAR_DIMS)   ! vector of dimension ids
 integer:: natt                    ! number of attributes

!Variables for setting up H2O table:
 integer:: t                     ! path temperature
 integer:: tmin                  ! mininum path temperature
 integer:: tmax                  ! maximum path temperature
 integer:: itype                 ! type of sat. pressure (=0 -> H2O only)
 real(r8):: tdbl

 integer:: i,istat,cam_abs_unit
 logical:: opened
 character(len=StrKIND):: errmess

 integer:: i_te,i_rh

!--------------------------------------------------------------------------------------------------

!... constants to set:
 p0     = pstdx
 amd    = mwdryx
 amco2  = mwco2x

!... coefficients for h2o emissivity and absorptivity for overlap of H2O and trace gases:
 c16  = coefj(3,1)/coefj(2,1)
 c17  = coefk(3,1)/coefk(2,1)
 c26  = coefj(3,2)/coefj(2,2)
 c27  = coefk(3,2)/coefk(2,2)
 c28  = .5
 c29  = .002053
 c30  = .1
 c31  = 3.0e-5

!... initialize further longwave constants referring to far wing correction for overlap of H2O
!    and trace gases; R&D refers to:
!    Ramanathan,V., and P.Downey, 1986:A Nonisothermal Emissivity and Absorptivity Formulation
!            for Water Vapor Journal of Geophysical Research, vol. 91., D8, pp 8649-8666.

 fwcoef = .1           ! See eq(33) R&D
 fwc1   = .30          ! See eq(33) R&D
 fwc2   = 4.5          ! See eq(33) and eq(34) in R&D
 fc1    = 2.6          ! See eq(34) R&D

 istat = -999
 if(dminfo % my_proc_id == IO_NODE) then
    do i = 10,99
       inquire(i,opened = opened,iostat=istat)
       if(.not. opened) then
          cam_abs_unit = i
          exit
       endif
    enddo
    if(istat /= 0) &
       call physics_error_fatal('module_ra_cam: radaeinit: Cannot find unused '//&
                                'fortran unit to read in lookup table.')
 endif

!distribute unit to other processors:
 call mpas_dmpar_bcast_int(dminfo,cam_abs_unit)

!open init file:
 if(dminfo % my_proc_id == IO_NODE) then

#ifdef SINGLE_PRECISION
    open(cam_abs_unit,file='CAM_ABS_DATA.TBL',form='UNFORMATTED',status='OLD',iostat=istat)
#else
    open(cam_abs_unit,file='CAM_ABS_DATA.DBL',form='UNFORMATTED',status='OLD',iostat=istat)
#endif
    if(istat /= 0) then
       write(errmess,'(A,I4)') 'module_ra_cam: error reading CAM_ABS_DATA on unit', &
             cam_abs_unit
       call physics_error_fatal(errmess)
    endif

    read(cam_abs_unit,iostat=istat) ah2onw
    read(cam_abs_unit,iostat=istat) eh2onw 
    read(cam_abs_unit,iostat=istat) ah2ow 
    read(cam_abs_unit,iostat=istat) cn_ah2ow 
    read(cam_abs_unit,iostat=istat) cn_eh2ow 
    read(cam_abs_unit,iostat=istat) ln_ah2ow 
    read(cam_abs_unit,iostat=istat) ln_eh2ow

 endif
 
 DM_BCAST_MACRO(ah2onw)
 DM_BCAST_MACRO(eh2onw)
 DM_BCAST_MACRO(ah2ow)
 DM_BCAST_MACRO(cn_ah2ow)
 DM_BCAST_MACRO(cn_eh2ow)
 DM_BCAST_MACRO(ln_ah2ow)
 DM_BCAST_MACRO(ln_eh2ow)

 if(dminfo % my_proc_id == IO_NODE) close(cam_abs_unit)

! Set up table of H2O saturation vapor pressures for use in calculation effective path RH.  
! Need separate table from table in wv_saturation because:
! (1. Path temperatures can fall below minimum of that table; and
! (2. Abs/Emissivity tables are derived with RH for water only.

 tmin = nint(min_tp_h2o)
 tmax = nint(max_tp_h2o)+1
 itype = 0
 do t = tmin, tmax
!   call gffgch(dble(t),estblh2o(t-tmin),itype)
    tdbl = t
    call gffgch(tdbl,estblh2o(t-tmin),itype)
 enddo

 end subroutine radaeini

!==================================================================================================
 subroutine aerosol_init(dminfo,mesh,diag,state_1,state_2)
!==================================================================================================

!This subroutine assumes a uniform aerosol distribution in both time and space. It should be
!modified if aerosol data are available from WRF-CHEM or other sources.

!input arguments:!
 type(dm_info),intent(in)  :: dminfo
 type(mesh_type),intent(in):: mesh
 type(diag_type),intent(in):: diag

!inout arguments:
 type(state_type),intent(inout):: state_1,state_2

!local variables:
 integer:: max_mxaerl
 integer,dimension(:),allocatable:: mxaerl_local
 integer:: iCell,k,kk,nAerLevels,nCells,nCellsSolve,nVertLevels

 real(kind=RKIND):: psurf
 real(kind=RKIND),dimension(:),pointer:: m_psp,m_psn
 real(kind=RKIND),dimension(:,:),pointer:: m_hybi
 real(kind=RKIND),dimension(:,:),pointer:: pressure_b
 real(kind=RKIND),dimension(:,:,:),pointer:: aerosolcn,aerosolcp

 real(kind=RKIND),dimension(29) :: hybi  
 data hybi/0, 0.0065700002014637, 0.0138600002974272, 0.023089999333024 , &
              0.0346900001168251, 0.0491999983787537, 0.0672300010919571, &
              0.0894500017166138, 0.116539999842644 , 0.149159997701645 , &
              0.187830001115799 , 0.232859998941422 , 0.284209996461868 , &
              0.341369986534119 , 0.403340011835098 , 0.468600004911423 , &
              0.535290002822876 , 0.601350009441376 , 0.66482001543045  , &
              0.724009990692139 , 0.777729988098145 , 0.825269997119904 , & 
              0.866419970989227 , 0.901350021362305 , 0.930540025234222 , & 
              0.954590022563934 , 0.974179983139038 , 0.990000009536743 , 1/

!--------------------------------------------------------------------------------------------------

!initialization:
 nCells      = mesh % nCells
 nCellsSolve = mesh % nCellsSolve
 nAerLevels  = mesh % nAerLevels
 nVertLevels = mesh % nVertLevels
 m_hybi => mesh % m_hybi % array

 pressure_b => diag % pressure_base % array

 m_psp => state_1 % m_ps % array
 m_psn => state_2 % m_ps % array

 aerosolcp => state_1 % aerosols % array
 aerosolcn => state_2 % aerosols % array 

!initialization of aerosol levels:
 do k = 1, nAerLevels
 do iCell = 1, nCells
    m_hybi(k,iCell) = hybi(k)
 enddo
 enddo

 psurf = 1.e05
 do iCell = 1, nCells
    m_psp(iCell) = psurf
    m_psn(iCell) = psurf
 enddo

!mxaerl = max number of levels (from bottom) for background aerosol. Limit background aerosol
!height below 900 mb:
 if(.not.allocated(mxaerl_local)) allocate(mxaerl_local(1:nCellsSolve))
 mxaerl_local(1:nCellsSolve) = 0
 do k = 1, nVertLevels
 do iCell = 1, nCellsSolve
    if(pressure_b(k,iCell) .ge. 9.e04) mxaerl_local(iCell) = mxaerl_local(iCell) + 1
 enddo
 enddo 

 max_mxaerl = maxval(mxaerl_local(1:nCellsSolve))
 call mpas_dmpar_max_int(dminfo,max_mxaerl,mxaerl)
 write(0,101) max_mxaerl,mxaerl
 101 format(/,1x,'--- AEROSOLS:mxaerl on node = ',i5, &
            /,1x,'    AEROSOLS:Background aerosol limited to bottom ',i5,' model interfaces')
 if(allocated(mxaerl_local)) deallocate(mxaerl_local)

!initialize indices for water species:
 ozncyc = .true.
 indirect = .true.
 ixcldliq = 2
 ixcldice = 3

!initialize indices for aerosol species:
 idxSUL         = state_1 % index_sul
 idxSSLT        = state_1 % index_sslt
 idxDUSTfirst   = state_1 % index_dust1
 idxOCPHO       = state_1 % index_ocpho
 idxCARBONFIRST = state_1 % index_ocpho
 idxBCPHO       = state_1 % index_bcpho
 idxOCPHI       = state_1 % index_ocphi
 idxBCPHI       = state_1 % index_bcphi
 idxBG          = state_1 % index_bg
 idxVOLC        = state_1 % index_volc
 
 write(0,*) '    idxSUL         =',idxSUL
 write(0,*) '    idxSSLT        =',idxSSLT
 write(0,*) '    idxDUSTfirst   =',idxDUSTfirst
 write(0,*) '    idxOCPHO       =',idxOCPHO
 write(0,*) '    idxCARBONfirst =',idxCARBONfirst
 write(0,*) '    idxBCPHO       =',idxBCPHO
 write(0,*) '    idxOCPHI       =',idxOCPHI
 write(0,*) '    idxBCPHI       =',idxBCPHI
 write(0,*) '    idxBG          =',idxBG
 write(0,*) '    idxVOLC        =',idxVOLC

 do iCell = 1, nCells
 do k = 1, nAerLevels
!aerosolc arrays are upward cumulative (kg/m2) at each level. Here we assume uniform vertical
!distribution (aerosolc linear with hybi)
    aerosolcp(idxSUL,k,iCell)=1.e-7*(1.-hybi(k))
    aerosolcn(idxSUL,k,iCell)=1.e-7*(1.-hybi(k))
    aerosolcp(idxSSLT,k,iCell)=1.e-22*(1.-hybi(k))
    aerosolcn(idxSSLT,k,iCell)=1.e-22*(1.-hybi(k))
    aerosolcp(idxDUSTfirst,k,iCell)=1.e-7*(1.-hybi(k))
    aerosolcn(idxDUSTfirst,k,iCell)=1.e-7*(1.-hybi(k))
    aerosolcp(idxDUSTfirst+1,k,iCell)=1.e-7*(1.-hybi(k))
    aerosolcn(idxDUSTfirst+1,k,iCell)=1.e-7*(1.-hybi(k))
    aerosolcp(idxDUSTfirst+2,k,iCell)=1.e-7*(1.-hybi(k))
    aerosolcn(idxDUSTfirst+2,k,iCell)=1.e-7*(1.-hybi(k))
    aerosolcp(idxDUSTfirst+3,k,iCell)=1.e-7*(1.-hybi(k))
    aerosolcn(idxDUSTfirst+3,k,iCell)=1.e-7*(1.-hybi(k))
    aerosolcp(idxOCPHO,k,iCell)=1.e-7*(1.-hybi(k))
    aerosolcn(idxOCPHO,k,iCell)=1.e-7*(1.-hybi(k))
    aerosolcp(idxBCPHO,k,iCell)=1.e-9*(1.-hybi(k))
    aerosolcn(idxBCPHO,k,iCell)=1.e-9*(1.-hybi(k))
    aerosolcp(idxOCPHI,k,iCell)=1.e-7*(1.-hybi(k))
    aerosolcn(idxOCPHI,k,iCell)=1.e-7*(1.-hybi(k))
    aerosolcp(idxBCPHI,k,iCell)=1.e-8*(1.-hybi(k))
    aerosolcn(idxBCPHI,k,iCell)=1.e-8*(1.-hybi(k))
 enddo
 enddo

 call aer_optics_initialize(dminfo)

 end subroutine aerosol_init

!==================================================================================================
 subroutine aer_optics_initialize(dminfo)
!==================================================================================================

!input arguments:
 type(dm_info):: dminfo

!local variables:
 integer:: nrh_opac  ! number of relative humidity values for OPAC data
 integer:: nbnd      ! number of spectral bands, should be identical to nspint
 integer:: krh_opac  ! rh index for OPAC rh grid
 integer:: krh       ! another rh index
 integer:: ksz       ! dust size bin index
 integer:: kbnd      ! band index

 real(r8), parameter :: wgt_sscm = 6.0 / 7.0
 real(r8):: rh       ! local relative humidity variable

 integer, parameter :: irh=8
 real(r8):: rh_opac(irh)             ! OPAC relative humidity grid
 real(r8):: ksul_opac(irh,nspint)    ! sulfate  extinction
 real(r8):: wsul_opac(irh,nspint)    !          single scattering albedo
 real(r8):: gsul_opac(irh,nspint)    !          asymmetry parameter
 real(r8):: ksslt_opac(irh,nspint)   ! sea-salt
 real(r8):: wsslt_opac(irh,nspint)
 real(r8):: gsslt_opac(irh,nspint)
 real(r8):: kssam_opac(irh,nspint)   ! sea-salt accumulation mode
 real(r8):: wssam_opac(irh,nspint)
 real(r8):: gssam_opac(irh,nspint)
 real(r8):: ksscm_opac(irh,nspint)   ! sea-salt coarse mode
 real(r8):: wsscm_opac(irh,nspint)
 real(r8):: gsscm_opac(irh,nspint)
 real(r8):: kcphil_opac(irh,nspint)  ! hydrophilic organic carbon
 real(r8):: wcphil_opac(irh,nspint)
 real(r8):: gcphil_opac(irh,nspint)
 real(r8):: dummy(nspint)

 integer:: i,istat,cam_aer_unit
 logical:: opened
 character(len=StrKIND):: errmess
 
!--------------------------------------------------------------------------------------------------

!write(0,*) '--- enter subroutine aer_optics_initialize:'

!READ AEROSOL OPTICS DATA:
 istat = -999
 if(dminfo % my_proc_id == IO_NODE) then
    do i = 10,99
       inquire(i,opened = opened,iostat=istat)
       if(.not. opened) then
          cam_aer_unit = i
          exit
       endif
    enddo
    if(istat /= 0) &
       call physics_error_fatal('module_ra_cam: aer_optics_initialize: Cannot find unused '//&
                                'fortran unit to read in lookup table.')
 endif

!distribute unit to other processors:
 call mpas_dmpar_bcast_int(dminfo,cam_aer_unit)

!open init file:
 if(dminfo % my_proc_id == IO_NODE) then

#ifdef SINGLE_PRECISION
    open(cam_aer_unit,file='CAM_AEROPT_DATA.TBL',form='UNFORMATTED',status='OLD',iostat=istat)
#else
    open(cam_aer_unit,file='CAM_AEROPT_DATA.DBL',form='UNFORMATTED',status='OLD',iostat=istat)
#endif
    if(istat /= 0) then
       write(errmess,'(A,I4)') 'module_ra_cam: error reading CAM_AEROPT_DATA on unit', &
             cam_aer_unit
       call physics_error_fatal(errmess)
    endif

    read(cam_aer_unit,iostat=istat) dummy
    read(cam_aer_unit,iostat=istat) rh_opac 
    read(cam_aer_unit,iostat=istat) ksul_opac 
    read(cam_aer_unit,iostat=istat) wsul_opac 
    read(cam_aer_unit,iostat=istat) gsul_opac 
    read(cam_aer_unit,iostat=istat) kssam_opac 
    read(cam_aer_unit,iostat=istat) wssam_opac 
    read(cam_aer_unit,iostat=istat) gssam_opac 
    read(cam_aer_unit,iostat=istat) ksscm_opac 
    read(cam_aer_unit,iostat=istat) wsscm_opac 
    read(cam_aer_unit,iostat=istat) gsscm_opac
    read(cam_aer_unit,iostat=istat) kcphil_opac 
    read(cam_aer_unit,iostat=istat) wcphil_opac 
    read(cam_aer_unit,iostat=istat) gcphil_opac 
    read(cam_aer_unit,iostat=istat) kcb 
    read(cam_aer_unit,iostat=istat) wcb 
    read(cam_aer_unit,iostat=istat) gcb 
    read(cam_aer_unit,iostat=istat) kdst 
    read(cam_aer_unit,iostat=istat) wdst 
    read(cam_aer_unit,iostat=istat) gdst 
    read(cam_aer_unit,iostat=istat) kbg 
    read(cam_aer_unit,iostat=istat) wbg 
    read(cam_aer_unit,iostat=istat) gbg
    read(cam_aer_unit,iostat=istat) kvolc 
    read(cam_aer_unit,iostat=istat) wvolc 
    read(cam_aer_unit,iostat=istat) gvolc

 endif

 DM_BCAST_MACRO(rh_opac)
 DM_BCAST_MACRO(ksul_opac)
 DM_BCAST_MACRO(wsul_opac)
 DM_BCAST_MACRO(gsul_opac)
 DM_BCAST_MACRO(kssam_opac)
 DM_BCAST_MACRO(wssam_opac)
 DM_BCAST_MACRO(gssam_opac)
 DM_BCAST_MACRO(ksscm_opac)
 DM_BCAST_MACRO(wsscm_opac)
 DM_BCAST_MACRO(gsscm_opac)
 DM_BCAST_MACRO(kcphil_opac)
 DM_BCAST_MACRO(wcphil_opac)
 DM_BCAST_MACRO(gcphil_opac)
 DM_BCAST_MACRO(kcb)
 DM_BCAST_MACRO(wcb)
 DM_BCAST_MACRO(gcb)
 DM_BCAST_MACRO(kdst)
 DM_BCAST_MACRO(wdst)
 DM_BCAST_MACRO(gdst)
 DM_BCAST_MACRO(kbg)
 DM_BCAST_MACRO(wbg)
 DM_BCAST_MACRO(gbg)
 DM_BCAST_MACRO(kvolc)
 DM_BCAST_MACRO(wvolc)
 DM_BCAST_MACRO(gvolc)

 if(dminfo % my_proc_id == IO_NODE) close(cam_aer_unit)

! map OPAC aerosol species onto CAM aerosol species
! CAM name             OPAC name
! sul   or SO4         = suso                  sulfate soluble
! sslt  or SSLT        = 1/7 ssam + 6/7 sscm   sea-salt accumulation/coagulation mode
! cphil or CPHI        = waso                  water soluble (carbon)
! cphob or CPHO        = waso @ rh = 0
! cb    or BCPHI/BCPHO = soot

 ksslt_opac(:,:) = (1.0 - wgt_sscm) * kssam_opac(:,:) + wgt_sscm * ksscm_opac(:,:)

 wsslt_opac(:,:) = ( (1.0 - wgt_sscm) * kssam_opac(:,:) * wssam_opac(:,:) &
                 + wgt_sscm * ksscm_opac(:,:) * wsscm_opac(:,:) ) &
                 / ksslt_opac(:,:)

 gsslt_opac(:,:) = ( (1.0 - wgt_sscm) * kssam_opac(:,:) * wssam_opac(:,:) * gssam_opac(:,:) &
                 + wgt_sscm * ksscm_opac(:,:) * wsscm_opac(:,:) * gsscm_opac(:,:) ) &
                 / ( ksslt_opac(:,:) * wsslt_opac(:,:) )

 do i = 1, nspint
    kcphob(i) = kcphil_opac(1,i)
    wcphob(i) = wcphil_opac(1,i)
    gcphob(i) = gcphil_opac(1,i)
 end do

!interpolate optical properties of hygrospopic aerosol species onto a uniform
!relative humidity grid:

 nbnd = nspint

 do krh = 1, nrh
    rh = 1.0_r8 / nrh * (krh - 1)
    do kbnd = 1, nbnd
       ksul(krh,kbnd) = exp_interpol(rh_opac, &
                 ksul_opac(:,kbnd) / ksul_opac(1,kbnd),rh) * ksul_opac(1,kbnd)
       wsul(krh,kbnd) = lin_interpol(rh_opac, &
                 wsul_opac(:,kbnd) / wsul_opac(1,kbnd),rh) * wsul_opac(1,kbnd)
       gsul(krh,kbnd) = lin_interpol(rh_opac, &
                 gsul_opac(:,kbnd) / gsul_opac(1,kbnd),rh) * gsul_opac(1,kbnd)
       ksslt(krh,kbnd) = exp_interpol(rh_opac, &
                 ksslt_opac(:,kbnd) / ksslt_opac(1,kbnd),rh) * ksslt_opac(1,kbnd)
       wsslt(krh,kbnd) = lin_interpol(rh_opac, &
                 wsslt_opac(:,kbnd) / wsslt_opac(1,kbnd),rh) * wsslt_opac(1,kbnd)
       gsslt(krh,kbnd) = lin_interpol(rh_opac, &
          gsslt_opac(:,kbnd) / gsslt_opac(1,kbnd),rh) * gsslt_opac(1,kbnd)
       kcphil(krh,kbnd) = exp_interpol(rh_opac, &
          kcphil_opac(:,kbnd) / kcphil_opac(1,kbnd),rh) * kcphil_opac(1,kbnd)
       wcphil(krh,kbnd) = lin_interpol(rh_opac, &
          wcphil_opac(:,kbnd) / wcphil_opac(1,kbnd),rh) * wcphil_opac(1,kbnd)
       gcphil(krh,kbnd) = lin_interpol(rh_opac, &
          gcphil_opac(:,kbnd) / gcphil_opac(1,kbnd),rh) * gcphil_opac(1,kbnd)
    enddo
 enddo

 write(0,*) '    end subroutine aer_optics_initialize:'

 end subroutine aer_optics_initialize

!==================================================================================================
 subroutine oznini(mesh,atm_input)
!==================================================================================================

!This subroutine assumes a uniform distribution of ozone concentration. It should be replaced
!with monthly climatology varying ozone distribution.

!input arguments:
 type(mesh_type),intent(in):: mesh

!inout arguments:
 type(atm_input_type),intent(inout):: atm_input

!local variables:
 integer,parameter:: pin_unit = 27
 integer,parameter:: lat_unit = 28
 integer,parameter:: oz_unit  = 29
 integer,parameter:: open_ok  = 0

 integer:: i,i1,i2,istat,k,j,m
 integer:: iCell,nCells,num_months,levsiz
 
 real(kind=RKIND),dimension(:),pointer:: latCell,lonCell
 real(kind=RKIND),dimension(:),pointer:: pin
 real(kind=RKIND),dimension(:,:,:),pointer:: ozmixm

 real(kind=RKIND):: lat,lon,dlat,dlatCell
 real(kind=RKIND),dimension(latsiz):: lat_ozone
!real(Kind=RKIND),dimension(lonsiz,levsiz,latsiz,num_months):: ozmixin
 real(Kind=RKIND),dimension(:,:,:,:),allocatable:: ozmixin

!--------------------------------------------------------------------------------------------------

 nCells     = mesh % nCells
 num_months = mesh % nMonths
 levsiz     = mesh % nOznLevels
 latCell => mesh % latCell % array
 lonCell => mesh % lonCell % array
 pin     => atm_input % pin % array
 ozmixm  => atm_input % ozmixm % array

!-- read in ozone pressure data:
 open(pin_unit,file='OZONE_PLEV.TBL',action='READ',status='OLD',iostat=istat)
 if(istat /= open_ok) &
    call physics_error_fatal('subroutine oznini: ' // &
                             'failure opening OZONE_PLEV.TBL')
 do k = 1,levsiz
    read(pin_unit,*) pin(k)
    pin(k) = pin(k)*100.
!   write(0,*) pin(k)
 enddo
 close(pin_unit)

!-- read in ozone lat data:
 open(lat_unit, file='OZONE_LAT.TBL',action='READ',status='OLD',iostat=istat) 
 if(istat /= open_ok) &
    call physics_error_fatal('subroutine oznini: ' // &
                             'failure opening OZONE_LAT.TBL')
 do j = 1, latsiz
    read(lat_unit,*) lat_ozone(j)
!   write(0,101) j,lat_ozone(j)
 enddo
 close(lat_unit)

!-- read in ozone data:
 open(oz_unit,file='OZONE_DAT.TBL',action='READ',status='OLD',iostat=istat)
 if(istat /= open_ok) &
    call physics_error_fatal('subroutine oznini: ' // &
                                'failure opening OZONE_DAT.TBL')

 allocate(ozmixin(lonsiz,levsiz,latsiz,num_months))
 do m=1,num_months
 do j=1,latsiz ! latsiz=64
 do k=1,levsiz ! levsiz=59
 do i=1,lonsiz ! lonsiz=1
    read(oz_unit,*) ozmixin(i,k,j,m)
 enddo
 enddo
 enddo
 enddo
 close(oz_unit)

!INTERPOLATION OF INPUT OZONE DATA TO MPAS GRID:
!write(0,*) 'max latCell=', maxval(latCell)/degrad
!write(0,*) 'min latCell=', minval(latCell)/degrad
!write(0,*) 'max lonCell=', maxval(lonCell)/degrad
!write(0,*) 'min lonCell=', minval(lonCell)/degrad
!write(0,*)
!write(0,*) 'max lat_ozone=',maxval(lat_ozone)
!write(0,*) 'min lat_ozone=',minval(lat_ozone)
 do iCell = 1,nCells
    lat = latCell(iCell)/degrad
    lon = lonCell(iCell)/degrad
    if(lat .gt. lat_ozone(latsiz)) then
     i1 = latsiz
     i2 = latsiz
    elseif(lat .lt. lat_ozone(1)) then
       i1 = 1
       i2 = 1
   else
       do i = 1, latsiz
          if(lat.ge.lat_ozone(i) .and. lat.lt.lat_ozone(i+1)) exit
       enddo
       i1 = i
       i2 = i+1
    endif

    do m = 1,num_months
    do k = 1,levsiz
    do j = 1,lonsiz
       dlat     = lat_ozone(i2)-lat_ozone(i1)
       dlatCell = lat-lat_ozone(i1)
       if(dlat == 0.) then
          ozmixm(m,k,iCell) = ozmixin(j,k,i1,m)
       else
          ozmixm(m,k,iCell) = ozmixin(j,k,i1,m) &
                     + (ozmixin(j,k,i2,m)-ozmixin(j,k,i1,m))*dlatCell/dlat
       endif
    enddo 
    enddo       
    enddo
!   do k = 1, levsiz
!      write(0,102) iCell,i1,i2,lat_ozone(i1),lat,lat_ozone(i2),ozmixin(1,k,i1,1), &
!                   ozmixm(1,k,iCell),ozmixin(1,k,i2,1)
!   enddo
 enddo
 deallocate(ozmixin)

!formats:
 101 format(i3,12(1x,e15.8))
 102 format(i6,i6,i6,6(1x,e15.8))

 end subroutine oznini

!==================================================================================================
 end module mpas_atmphys_camrad_init
!==================================================================================================
