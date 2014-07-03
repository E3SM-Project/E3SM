module spehox

!----------------------------------------------------------------------- 
!
! BOP
!
! !MODULE: spedata
!
! !DESCRIPTION
! Determines the HOx production factor assoctioned with
! solar proton ionization.
!
! !USES
  use shr_kind_mod, only: r8 => shr_kind_r8,r4 => shr_kind_r4
  use abortutils,   only: endrun

  implicit none

  private  ! all unless made public
  save

  public hox_prod_factor

! !REVISION HISTORY:
!   17 Nov 2005  Francis Vitt     Creation
!
! EOP
!----------------------------------------------------------------------- 
! $Id: spehox.F90,v 1.1.2.1 2006/05/03 20:53:09 stacy Exp $
! $Author: stacy $
!----------------------------------------------------------------------- 
! HOx production per ion pair (cm^-3 s-1) from Figure 2 
! of Solomon et al. (1981)
! 
! Data source: 
!           Dr. Charles H. Jackman, Code 916                       
!           NASA/Goddard Space Flight Center                       
!                Greenbelt, MD  20771                              
!  ph:301-614-6053  fx:301-614-5903  Charles.H.Jackman@nasa.gov    
! 
!   Alt(km)     10        100      1000      10000    100000
!  
!    40.0      2.00      2.00      2.00      1.99      1.99
!    42.5      2.00      2.00      1.99      1.99      1.99
!    45.0      2.00      2.00      1.99      1.99      1.99
!    47.5      2.00      2.00      1.99      1.99      1.98
!    50.0      2.00      1.99      1.99      1.98      1.98
!    52.5      2.00      1.99      1.99      1.98      1.95
!    55.0      2.00      1.99      1.98      1.97      1.93
!    57.5      2.00      1.99      1.98      1.95      1.89
!    60.0      1.99      1.98      1.97      1.94      1.85
!    62.5      1.99      1.98      1.96      1.90      1.81
!    65.0      1.99      1.98      1.94      1.87      1.77
!    67.5      1.98      1.96      1.91      1.82      1.72
!    70.0      1.98      1.94      1.87      1.77      1.64
!    72.5      1.96      1.90      1.80      1.70      1.50
!    75.0      1.93      1.84      1.73      1.60      1.30
!    77.5      1.84      1.72      1.60      1.40      0.93
!    80.0      1.60      1.40      1.20      0.95      0.40
!    82.5      0.80      0.60      0.40      0.15      0.00
!    85.0      0.30      0.15      0.10      0.00      0.00
!    87.5      0.00      0.00      0.00      0.00      0.00
!    90.0      0.00      0.00      0.00      0.00      0.00
!----------------------------------------------------------------------- 

  integer, parameter :: nalts = 21
  integer, parameter :: nprods = 5
  real(r8) :: alts(nalts)
  real(r8) :: log_ion_prod(nprods)
  real(r8) :: factor_tbl(nalts,nprods)

  data alts(1:21) / 40.0_r8, 42.5_r8, 45.0_r8, 47.5_r8, 50.0_r8, 52.5_r8, 55.0_r8, 57.5_r8, &
                    60.0_r8, 62.5_r8, 65.0_r8, 67.5_r8, 70.0_r8, 72.5_r8, 75.0_r8, 77.5_r8, &
                    80.0_r8, 82.5_r8, 85.0_r8, 87.5_r8, 90.0_r8 / 

  data log_ion_prod(1:5)  /  1._r8,    2._r8,    3._r8,    4._r8,    5._r8 /

  data factor_tbl( 1,1:5) /  2.00_r8,  2.00_r8,  2.00_r8,  1.99_r8,  1.99_r8 /
  data factor_tbl( 2,1:5) /  2.00_r8,  2.00_r8,  2.00_r8,  1.99_r8,  1.99_r8 /
  data factor_tbl( 3,1:5) /  2.00_r8,  2.00_r8,  2.00_r8,  1.99_r8,  1.99_r8 /
  data factor_tbl( 4,1:5) /  2.00_r8,  2.00_r8,  2.00_r8,  1.99_r8,  1.99_r8 /
  data factor_tbl( 5,1:5) /  2.00_r8,  2.00_r8,  2.00_r8,  1.99_r8,  1.99_r8 /
  data factor_tbl( 6,1:5) /  2.00_r8,  1.99_r8,  1.99_r8,  1.98_r8,  1.95_r8 /
  data factor_tbl( 7,1:5) /  2.00_r8,  1.99_r8,  1.98_r8,  1.97_r8,  1.93_r8 /
  data factor_tbl( 8,1:5) /  2.00_r8,  1.99_r8,  1.98_r8,  1.95_r8,  1.89_r8 /
  data factor_tbl( 9,1:5) /  1.99_r8,  1.98_r8,  1.97_r8,  1.94_r8,  1.85_r8 /
  data factor_tbl(10,1:5) /  1.99_r8,  1.98_r8,  1.96_r8,  1.90_r8,  1.81_r8 /
  data factor_tbl(11,1:5) /  1.99_r8,  1.98_r8,  1.94_r8,  1.87_r8,  1.77_r8 /
  data factor_tbl(12,1:5) /  1.98_r8,  1.96_r8,  1.91_r8,  1.82_r8,  1.72_r8 /
  data factor_tbl(13,1:5) /  1.98_r8,  1.94_r8,  1.87_r8,  1.77_r8,  1.64_r8 /
  data factor_tbl(14,1:5) /  1.96_r8,  1.90_r8,  1.80_r8,  1.70_r8,  1.50_r8 /
  data factor_tbl(15,1:5) /  1.93_r8,  1.84_r8,  1.73_r8,  1.60_r8,  1.30_r8 /
  data factor_tbl(16,1:5) /  1.84_r8,  1.72_r8,  1.60_r8,  1.40_r8,  0.93_r8 /
  data factor_tbl(17,1:5) /  1.60_r8,  1.40_r8,  1.20_r8,  0.95_r8,  0.40_r8 /
  data factor_tbl(18,1:5) /  0.80_r8,  0.60_r8,  0.40_r8,  0.15_r8,  0.00_r8 /
  data factor_tbl(19,1:5) /  0.30_r8,  0.15_r8,  0.10_r8,  0.00_r8,  0.00_r8 /
  data factor_tbl(20,1:5) /  0.00_r8,  0.00_r8,  0.00_r8,  0.00_r8,  0.00_r8 /
  data factor_tbl(21,1:5) /  0.00_r8,  0.00_r8,  0.00_r8,  0.00_r8,  0.00_r8 /

contains

!----------------------------------------------------------------------- 
! Returns the HOx production factor for the ionization profile
!----------------------------------------------------------------------- 
  function hox_prod_factor( ion_pairs, zmid )

    use ppgrid,  only : pcols, pver
    ! for each level interpolate the factor table

    implicit none

    real(r8),intent(in) :: ion_pairs(pver)
    real(r8),intent(in) :: zmid(pver)
    real(r8) :: hox_prod_factor(pver)

    integer :: k
    integer :: lastk


    lastk = 1

    ! start at the bottom since the table goes from bottom to top
    do k = pver,1,-1
       hox_prod_factor(k) = interp_factor_tbl( ion_pairs(k), zmid(k), lastk )
    enddo

  end function hox_prod_factor

!----------------------------------------------------------------------- 
! bilinear interpolates the above table of factors
!----------------------------------------------------------------------- 
  function interp_factor_tbl( ionp, z, lastk )

    implicit none

    real(r8),intent(in) :: ionp
    real(r8),intent(in) :: z
    integer, intent(inout) :: lastk
    real(r8) :: interp_factor_tbl

    integer :: i,j
    real(r8) :: logp
    real(r8) :: atlwgt1,atlwgt2
    real(r8) :: prodwgt1,prodwgt2
    real(r8) :: fact(nprods)

    if (ionp <= 0._r8 ) then 
       interp_factor_tbl = 0._r8
       return
    endif

    ! interpolate log10 of the ionization rate since the table is 
    ! on a log scale.
    logp = log10(ionp)

    if ( z <= alts(1) ) then 
       fact(:) = factor_tbl(1,:)
    else if ( z >= alts(nalts) ) then
       fact(:) = factor_tbl(nalts,:)
    else
       do i = lastk,nalts
          if ( z > alts(i) .and. z <= alts(i+1) ) then 
             atlwgt1 = (alts(i+1) - z)/(alts(i+1) - alts(i))
             atlwgt2 =   (z - alts(i))/(alts(i+1) - alts(i))
             fact(:) = atlwgt1*factor_tbl(i,:) + atlwgt2*factor_tbl(i+1,:)
             lastk = i
             exit
          endif
       enddo
    endif

    if ( logp <= log_ion_prod(1) ) then 
       interp_factor_tbl = fact(1)
    else if ( logp >= log_ion_prod(nprods) ) then
       interp_factor_tbl = fact(nprods)
    else 
       do i = 1,nprods
          if ( logp > log_ion_prod(i) .and. logp <= log_ion_prod(i+1) ) then 
             prodwgt1 = (log_ion_prod(i+1) - logp)/(log_ion_prod(i+1) - log_ion_prod(i))
             prodwgt2 =   (logp - log_ion_prod(i))/(log_ion_prod(i+1) - log_ion_prod(i))
             interp_factor_tbl = prodwgt1*fact(i) + prodwgt2*fact(i+1)
             exit
          endif
       enddo
    endif

  end function interp_factor_tbl

end module spehox
