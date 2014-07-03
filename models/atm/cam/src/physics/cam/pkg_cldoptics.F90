module pkg_cldoptics

!---------------------------------------------------------------------------------
! Purpose:
!
! Compute cloud optical properties: liquid and ice partical size; emissivity
!
! Author: Byron Boville  Sept 06, 2002, assembled from existing subroutines
!
!---------------------------------------------------------------------------------

  use shr_kind_mod,  only: r8=>shr_kind_r8
  use ppgrid,        only: pcols, pver, pverp

  implicit none
  private
  save

  public :: cldefr, cldems, cldovrlap, cldclw, reitab, reltab

contains

!===============================================================================
  subroutine cldefr(lchnk   ,ncol    , &
       landfrac,t       ,rel     ,rei     ,ps      ,pmid    , landm, icefrac, snowh)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute cloud water and ice particle size 
! 
! Method: 
! use empirical formulas to construct effective radii
! 
! Author: J.T. Kiehl, B. A. Boville, P. Rasch
! 
!-----------------------------------------------------------------------

!------------------------------Arguments--------------------------------
!
! Input arguments
!
    integer, intent(in) :: lchnk                 ! chunk identifier
    integer, intent(in) :: ncol                  ! number of atmospheric columns

    real(r8), intent(in) :: landfrac(pcols)      ! Land fraction
    real(r8), intent(in) :: icefrac(pcols)       ! Ice fraction
    real(r8), intent(in) :: t(pcols,pver)        ! Temperature
    real(r8), intent(in) :: ps(pcols)            ! Surface pressure
    real(r8), intent(in) :: pmid(pcols,pver)     ! Midpoint pressures
    real(r8), intent(in) :: landm(pcols)
    real(r8), intent(in) :: snowh(pcols)         ! Snow depth over land, water equivalent (m)
!
! Output arguments
!
    real(r8), intent(out) :: rel(pcols,pver)      ! Liquid effective drop size (microns)
    real(r8), intent(out) :: rei(pcols,pver)      ! Ice effective drop size (microns)
!
! following Kiehl
         call reltab(ncol, t, landfrac, landm, icefrac, rel, snowh)

! following Kristjansson and Mitchell
         call reitab(ncol, t, rei)

    return
  end subroutine cldefr

!===============================================================================
  subroutine cldems(lchnk   ,ncol    ,clwp    ,fice    ,rei     ,emis    ,cldtau)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute cloud emissivity using cloud liquid water path (g/m**2)
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: J.T. Kiehl
! 
!-----------------------------------------------------------------------

    use phys_control,    only: phys_getopts

!------------------------------Parameters-------------------------------
!
    real(r8) kabsl                  ! longwave liquid absorption coeff (m**2/g)
    parameter (kabsl = 0.090361_r8)
!
!------------------------------Arguments--------------------------------
!
! Input arguments
!
    integer, intent(in) :: lchnk                   ! chunk identifier
    integer, intent(in) :: ncol                    ! number of atmospheric columns

    real(r8), intent(in) :: clwp(pcols,pver)       ! cloud liquid water path (g/m**2)
    real(r8), intent(in) :: rei(pcols,pver)        ! ice effective drop size (microns)
    real(r8), intent(in) :: fice(pcols,pver)       ! fractional ice content within cloud
!
! Output arguments
!
    real(r8), intent(out) :: emis(pcols,pver)       ! cloud emissivity (fraction)
    real(r8), intent(out) :: cldtau(pcols,pver)     ! cloud optical depth
!
!---------------------------Local workspace-----------------------------
!
    integer i,k                 ! longitude, level indices
    real(r8) kabs                   ! longwave absorption coeff (m**2/g)
    real(r8) kabsi                  ! ice absorption coefficient

    character(len=16) :: microp_scheme  ! microphysics scheme
!-----------------------------------------------------------------------
!
    call phys_getopts(microp_scheme_out=microp_scheme)

    do k=1,pver
       do i=1,ncol

          !note that optical properties for ice valid only
          !in range of 13 > rei > 130 micron (Ebert and Curry 92)
    	  if ( microp_scheme == 'RK' ) then
             kabsi = 0.005_r8 + 1._r8/rei(i,k)
          else
             kabsi = 0.005_r8 + 1._r8/min(max(13._r8,rei(i,k)),130._r8)
          end if
          kabs = kabsl*(1._r8-fice(i,k)) + kabsi*fice(i,k)
          emis(i,k) = 1._r8 - exp(-1.66_r8*kabs*clwp(i,k))
          cldtau(i,k) = kabs*clwp(i,k)
       end do
    end do
!
    return
  end subroutine cldems

!===============================================================================
  subroutine cldovrlap(lchnk   ,ncol    ,pint    ,cld     ,nmxrgn  ,pmxrgn  )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Partitions each column into regions with clouds in neighboring layers.
! This information is used to implement maximum overlap in these regions
! with random overlap between them.
! On output,
!    nmxrgn contains the number of regions in each column
!    pmxrgn contains the interface pressures for the lower boundaries of
!           each region! 
! Method: 

! 
! Author: W. Collins
! 
!-----------------------------------------------------------------------

!
! Input arguments
!
    integer, intent(in) :: lchnk                ! chunk identifier
    integer, intent(in) :: ncol                 ! number of atmospheric columns

    real(r8), intent(in) :: pint(pcols,pverp)   ! Interface pressure
    real(r8), intent(in) :: cld(pcols,pver)     ! Fractional cloud cover
!
! Output arguments
!
    integer,  intent(out) :: nmxrgn(pcols)      ! Number of maximally overlapped regions
    real(r8), intent(out) :: pmxrgn(pcols,pverp)! Maximum values of pressure for each
!    maximally overlapped region.
!    0->pmxrgn(i,1) is range of pressure for
!    1st region,pmxrgn(i,1)->pmxrgn(i,2) for
!    2nd region, etc
!
!---------------------------Local variables-----------------------------
!
    integer i                    ! Longitude index
    integer k                    ! Level index
    integer n                    ! Max-overlap region counter

    real(r8) pnm(pcols,pverp)    ! Interface pressure

    logical cld_found            ! Flag for detection of cloud
    logical cld_layer(pver)      ! Flag for cloud in layer
!
!------------------------------------------------------------------------
!

    do i = 1, ncol
       cld_found = .false.
       cld_layer(:) = cld(i,:) > 0.0_r8
       pmxrgn(i,:) = 0.0_r8
       pnm(i,:)=pint(i,:)*10._r8
       n = 1
       do k = 1, pver
          if (cld_layer(k) .and.  .not. cld_found) then
             cld_found = .true.
          else if ( .not. cld_layer(k) .and. cld_found) then
             cld_found = .false.
             if (count(cld_layer(k:pver)) == 0) then
                exit
             endif
             pmxrgn(i,n) = pnm(i,k)
             n = n + 1
          endif
       end do
       pmxrgn(i,n) = pnm(i,pverp)
       nmxrgn(i) = n
    end do

    return
  end subroutine cldovrlap

!===============================================================================
  subroutine cldclw(lchnk   ,ncol    ,zi      ,clwp    ,tpw     ,hl      )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Evaluate cloud liquid water path clwp (g/m**2)
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: J.T. Kiehl
! 
!-----------------------------------------------------------------------


!
! Input arguments
!
    integer, intent(in) :: lchnk                 ! chunk identifier
    integer, intent(in) :: ncol                  ! number of atmospheric columns

    real(r8), intent(in) :: zi(pcols,pverp)      ! height at layer interfaces(m)
    real(r8), intent(in) :: tpw(pcols)           ! total precipitable water (mm)
!
! Output arguments
!
    real(r8), intent(out) :: clwp(pcols,pver)     ! cloud liquid water path (g/m**2)
    real(r8), intent(out) :: hl(pcols)            ! liquid water scale height

!
!---------------------------Local workspace-----------------------------
!
    integer  :: i,k                  ! longitude, level indices
    real(r8) :: clwc0                ! reference liquid water concentration (g/m**3)
    real(r8) :: emziohl(pcols,pverp) ! exp(-zi/hl)
    real(r8) :: rhl(pcols)           ! 1/hl
!
!-----------------------------------------------------------------------
!
! Set reference liquid water concentration
!
    clwc0 = 0.21_r8
!
! Diagnose liquid water scale height from precipitable water
!
    do i=1,ncol
       hl(i)  = 700.0_r8*log(max(tpw(i)+1.0_r8,1.0_r8))
       rhl(i) = 1.0_r8/hl(i)
    end do
!
! Evaluate cloud liquid water path (vertical integral of exponential fn)
!
    do k=1,pverp
       do i=1,ncol
          emziohl(i,k) = exp(-zi(i,k)*rhl(i))
       end do
    end do
    do k=1,pver
       do i=1,ncol
          clwp(i,k) = clwc0*hl(i)*(emziohl(i,k+1) - emziohl(i,k))
       end do
    end do
!
    return
  end subroutine cldclw


!===============================================================================
  subroutine reltab(ncol, t, landfrac, landm, icefrac, rel, snowh)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute cloud water size
! 
! Method: 
! analytic formula following the formulation originally developed by J. T. Kiehl
! 
! Author: Phil Rasch
! 
!-----------------------------------------------------------------------
    use physconst,          only: tmelt
!------------------------------Arguments--------------------------------
!
! Input arguments
!
    integer, intent(in) :: ncol
    real(r8), intent(in) :: landfrac(pcols)      ! Land fraction
    real(r8), intent(in) :: icefrac(pcols)       ! Ice fraction
    real(r8), intent(in) :: snowh(pcols)         ! Snow depth over land, water equivalent (m)
    real(r8), intent(in) :: landm(pcols)         ! Land fraction ramping to zero over ocean
    real(r8), intent(in) :: t(pcols,pver)        ! Temperature

!
! Output arguments
!
    real(r8), intent(out) :: rel(pcols,pver)      ! Liquid effective drop size (microns)
!
!---------------------------Local workspace-----------------------------
!
    integer i,k               ! Lon, lev indices
    real(r8) rliqland         ! liquid drop size if over land
    real(r8) rliqocean        ! liquid drop size if over ocean
    real(r8) rliqice          ! liquid drop size if over sea ice
!
!-----------------------------------------------------------------------
!
    rliqocean = 14.0_r8
    rliqice   = 14.0_r8
    rliqland  = 8.0_r8
    do k=1,pver
       do i=1,ncol
! jrm Reworked effective radius algorithm
          ! Start with temperature-dependent value appropriate for continental air
          ! Note: findmcnew has a pressure dependence here
          rel(i,k) = rliqland + (rliqocean-rliqland) * min(1.0_r8,max(0.0_r8,(tmelt-t(i,k))*0.05_r8))
          ! Modify for snow depth over land
          rel(i,k) = rel(i,k) + (rliqocean-rel(i,k)) * min(1.0_r8,max(0.0_r8,snowh(i)*10._r8))
          ! Ramp between polluted value over land to clean value over ocean.
          rel(i,k) = rel(i,k) + (rliqocean-rel(i,k)) * min(1.0_r8,max(0.0_r8,1.0_r8-landm(i)))
          ! Ramp between the resultant value and a sea ice value in the presence of ice.
          rel(i,k) = rel(i,k) + (rliqice-rel(i,k)) * min(1.0_r8,max(0.0_r8,icefrac(i)))
! end jrm
       end do
    end do
  end subroutine reltab

!===============================================================================
  subroutine reitab(ncol, t, re)

    use phys_control,    only: do_waccm_phys

    integer, intent(in) :: ncol
    real(r8), intent(out) :: re(pcols,pver)
    real(r8), intent(in) :: t(pcols,pver)
    integer , parameter :: len_retab = 138
    real(r8), parameter :: min_retab = 136._r8
    real(r8) retab(len_retab)
    real(r8) corr
    integer i
    integer k
    integer index
    integer min_index
    !
    !       Tabulated values of re(T) in the temperature interval
    !       180 K -- 274 K; hexagonal columns assumed:
    !
    !       Modified for pmc formation: 136K -- 274K    
    !
    data retab /                                                &
         0.05_r8,    0.05_r8,    0.05_r8,    0.05_r8,    0.05_r8,   0.05_r8,      &
         0.055_r8,   0.06_r8,    0.07_r8,    0.08_r8,    0.09_r8,    0.1_r8,      &
         0.2_r8,     0.3_r8,     0.40_r8,    0.50_r8,    0.60_r8,    0.70_r8,     &
         0.8_r8 ,    0.9_r8,     1.0_r8,     1.1_r8,     1.2_r8,     1.3_r8,      &
         1.4_r8,     1.5_r8,     1.6_r8,     1.8_r8,     2.0_r8,     2.2_r8,      &
         2.4_r8,     2.6_r8,     2.8_r8,     3.0_r8,     3.2_r8,     3.5_r8,      &
         3.8_r8,     4.1_r8,     4.4_r8,     4.7_r8,     5.0_r8,     5.3_r8,      &
         5.6_r8, & ! If WACCM physics is not on, only the values below this line are used.
         5.92779_r8, 6.26422_r8, 6.61973_r8, 6.99539_r8, 7.39234_r8,           &
         7.81177_r8, 8.25496_r8, 8.72323_r8, 9.21800_r8, 9.74075_r8, 10.2930_r8,	&
         10.8765_r8, 11.4929_r8, 12.1440_r8, 12.8317_r8, 13.5581_r8, 14.2319_r8, 	&
         15.0351_r8, 15.8799_r8, 16.7674_r8, 17.6986_r8, 18.6744_r8, 19.6955_r8,	&
         20.7623_r8, 21.8757_r8, 23.0364_r8, 24.2452_r8, 25.5034_r8, 26.8125_r8,	&
         27.7895_r8, 28.6450_r8, 29.4167_r8, 30.1088_r8, 30.7306_r8, 31.2943_r8, 	&
         31.8151_r8, 32.3077_r8, 32.7870_r8, 33.2657_r8, 33.7540_r8, 34.2601_r8, 	&
         34.7892_r8, 35.3442_r8, 35.9255_r8, 36.5316_r8, 37.1602_r8, 37.8078_r8,	&
         38.4720_r8, 39.1508_r8, 39.8442_r8, 40.5552_r8, 41.2912_r8, 42.0635_r8,	&
         42.8876_r8, 43.7863_r8, 44.7853_r8, 45.9170_r8, 47.2165_r8, 48.7221_r8,	&
         50.4710_r8, 52.4980_r8, 54.8315_r8, 57.4898_r8, 60.4785_r8, 63.7898_r8,	&
         65.5604_r8, 71.2885_r8, 75.4113_r8, 79.7368_r8, 84.2351_r8, 88.8833_r8,	&
         93.6658_r8, 98.5739_r8, 103.603_r8, 108.752_r8, 114.025_r8, 119.424_r8, 	&
         124.954_r8, 130.630_r8, 136.457_r8, 142.446_r8, 148.608_r8, 154.956_r8,	&
         161.503_r8, 168.262_r8, 175.248_r8, 182.473_r8, 189.952_r8, 197.699_r8,	&
         205.728_r8, 214.055_r8, 222.694_r8, 231.661_r8, 240.971_r8, 250.639_r8/	
    !
    save retab
    !

    ! This is done purely to preserve answers that were set when WACCM and
    ! CAM had separate tables.
    if (do_waccm_phys()) then
       min_index = 1
    else
       min_index = 44
    end if

    do k=1,pver
       do i=1,ncol
          index = int(t(i,k)-min_retab)
          index = min(max(index,min_index),len_retab-1)
          corr = t(i,k) - int(t(i,k))
          re(i,k) = retab(index)*(1._r8-corr)		&
               +retab(index+1)*corr
          !           re(i,k) = amax1(amin1(re(i,k),30.),10.)
       end do
    end do
    !
    return
  end subroutine reitab

end module pkg_cldoptics
