! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!!  This routine evaluates the coagulation kernels, ckernel(k,j1,j2,i1,i2)
!!  [cm^3 s^-1] and pkernel. Indices correspond to aritrary array of columns <ic, iy>
!!  vertical level <k>, aerosol groups <j1,j2> and bins <i1,i2> of colliding particles.
!!
!!  ckernel is calculated as a static array for use each timestep
!!  ckern0 is also created for a basis to calculate new ckernels each timestep, if desired. (coagwet.f)
!!
!!  This routine requires that vertical profiles of temperature <T>,
!!  air density <rhoa>, and viscosity <rmu> are defined.
!!
!!  @version Oct-1995  
!!  @author  Andy Ackerman 
subroutine setupckern(carma, cstate, rc)

  ! types
  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carmastate_mod
  use carma_mod
  
  implicit none

  type(carma_type), intent(in)         :: carma   !! the carma object
  type(carmastate_type), intent(inout) :: cstate  !! the carma state object
  integer, intent(inout)               :: rc      !! return code, negative indicates failure
  
  ! Local declarations
  ! 2-D collision efficiency for current group pair under
  ! consideration (for extrapolation of input data)
  real(kind=f)            :: e_coll2(NBIN,NBIN)
  integer, parameter      :: NP_DATA = 21       ! number of collector/collected pairs in input data
  integer, parameter      :: NR_DATA = 12       ! number of radius bins in input data
  real(kind=f), parameter :: e_small = 0.0001_f   ! smallest collision efficiency
  logical, save           :: init_data = .FALSE. ! did data_p and data_r get initialized?
  real(kind=f), save      :: data_p(NP_DATA)          ! radius ratios (collected/collector)
  real(kind=f), save      :: data_r(NR_DATA)          ! collector drop radii (um)
  real(kind=f), save      :: data_e(NP_DATA, NR_DATA) ! geometric collection efficiencies
      
  integer :: ip
  integer :: ig, jg
  real(kind=f) :: cstick_calc    ! the probability that two particles that collide through thermal coagulation will stick to each other.
  integer :: i1, i2, j1, j2, k
  integer :: i, j
  integer :: igrp
  integer :: ibin

  real(kind=f) :: rhoa_cgs
  real(kind=f) :: temp1, temp2

  real(kind=f) :: r1
  real(kind=f) :: di
  real(kind=f) :: gi
  real(kind=f) :: rlbi
  real(kind=f) :: dti1
  real(kind=f) :: dti2
  real(kind=f) :: dti

  real(kind=f) :: r2
  real(kind=f) :: dj
  real(kind=f) :: gj
  real(kind=f) :: rlbj 
  real(kind=f) :: dtj1
  real(kind=f) :: dtj2 
  real(kind=f) :: dtj 

  real(kind=f) :: rp
  real(kind=f) :: dp
  real(kind=f) :: gg
  real(kind=f) :: delt
  real(kind=f) :: term1
  real(kind=f) :: term2
  real(kind=f) :: cbr

  real(kind=f) :: r_larg
  real(kind=f) :: r_smal
  integer :: i_larg
  integer :: i_smal
  integer :: ig_larg
  integer :: ig_smal
  real(kind=f) :: d_larg

  real(kind=f) :: re_larg
  real(kind=f) :: pe 
  real(kind=f) :: pe3 
  real(kind=f) :: ccd 

  real(kind=f) :: e_coll
  real(kind=f) :: vfc_smal
  real(kind=f) :: vfc_larg 
  real(kind=f) :: sk
  real(kind=f) :: e1
  real(kind=f) :: e3
  real(kind=f) :: e_langmuir
  real(kind=f) :: re60

  real(kind=f) :: pr 
  real(kind=f) :: e_fuchs

  integer :: jp, jj, jr

  real(kind=f) :: pblni
  real(kind=f) :: rblni 

  real(kind=f) :: term3
  real(kind=f) :: term4

  real(kind=f) :: beta
  real(kind=f) :: b_coal
  real(kind=f) :: a_coal 
  real(kind=f) :: x_coal 
  real(kind=f) :: e_coal
  real(kind=f) :: vfc_1
  real(kind=f) :: vfc_2
  real(kind=f) :: cgr 
  

!  Add constants for calculating effect of Van Der Waal's forces on coagulation
!  See Chan and Mozurkewich, J. Atmos. Sci., June 2001  
  real(kind=f), parameter :: vwa1 = 0.0757_f
  real(kind=f), parameter :: vwa3 = 0.0015_f
  real(kind=f), parameter :: vwb0 = 0.0151_f
  real(kind=f), parameter :: vwb1 = -0.186_f
  real(kind=f), parameter :: vwb3 = -0.0163_f
  real(kind=f), parameter :: ham  = 6.4e-13_f   ! erg, Hamaker constant
  real(kind=f) :: hp, hpln, Enot, Einf
  logical      :: use_vw(NGROUP, NGROUP)
  integer      :: ielem
  

!  Initialization of input data for gravitational collection.
!  The data were compiled by Hall (J. Atmos. Sci. 37, 2486-2507, 1980).

  data data_p/0.00_f,0.05_f,0.10_f,0.15_f,0.20_f,0.25_f,0.30_f,0.35_f,0.40_f,0.45_f, &
    0.50_f,0.55_f,0.60_f,0.65_f,0.70_f,0.75_f,0.80_f,0.85_f,0.90_f,0.95_f,1.00_f/

  data data_r( 1), (data_e(ip, 1),ip=1,NP_DATA) /   10.0, &
    0.0001, 0.0001, 0.0001, 0.0001, 0.0140, 0.0170, 0.0190, 0.0220, &
    0.0270, 0.0300, 0.0330, 0.0350, 0.0370, 0.0380, 0.0380, 0.0370, &
    0.0360, 0.0350, 0.0320, 0.0290, 0.0270 /
  data data_r( 2), (data_e(ip, 2),ip=1,NP_DATA) /   20.0, &
    0.0001, 0.0001, 0.0001, 0.0050, 0.0160, 0.0220, 0.0300, 0.0430, &
    0.0520, 0.0640, 0.0720, 0.0790, 0.0820, 0.0800, 0.0760, 0.0670, &
    0.0570, 0.0480, 0.0400, 0.0330, 0.0270 /
  data data_r( 3), (data_e(ip, 3),ip=1,NP_DATA) /   30.0, &
    0.0001, 0.0001, 0.0020, 0.0200, 0.0400, 0.0850, 0.1700, 0.2700, &
    0.4000, 0.5000, 0.5500, 0.5800, 0.5900, 0.5800, 0.5400, 0.5100, &
    0.4900, 0.4700, 0.4500, 0.4700, 0.5200 /
  data data_r( 4), (data_e(ip, 4),ip=1,NP_DATA) /   40.0, &
    0.0001, 0.0010, 0.0700, 0.2800, 0.5000, 0.6200, 0.6800, 0.7400, &
    0.7800, 0.8000, 0.8000, 0.8000, 0.7800, 0.7700, 0.7600, 0.7700, &
    0.7700, 0.7800, 0.7900, 0.9500, 1.4000 /
  data data_r( 5), (data_e(ip, 5),ip=1,NP_DATA) /   50.0, &
    0.0001, 0.0050, 0.4000, 0.6000, 0.7000, 0.7800, 0.8300, 0.8600, &
    0.8800, 0.9000, 0.9000, 0.9000, 0.9000, 0.8900, 0.8800, 0.8800, &
    0.8900, 0.9200, 1.0100, 1.3000, 2.3000 /
  data data_r( 6), (data_e(ip, 6),ip=1,NP_DATA) /   60.0, &
    0.0001, 0.0500, 0.4300, 0.6400, 0.7700, 0.8400, 0.8700, 0.8900, &
    0.9000, 0.9100, 0.9100, 0.9100, 0.9100, 0.9100, 0.9200, 0.9300, &
    0.9500, 1.0000, 1.0300, 1.7000, 3.0000 /
  data data_r( 7), (data_e(ip, 7),ip=1,NP_DATA) /   70.0, &
    0.0001, 0.2000, 0.5800, 0.7500, 0.8400, 0.8800, 0.9000, 0.9200, &
    0.9400, 0.9500, 0.9500, 0.9500, 0.9500, 0.9500, 0.9500, 0.9700, &
    1.0000, 1.0200, 1.0400, 2.3000, 4.0000 /
  data data_r( 8), (data_e(ip, 8),ip=1,NP_DATA) /  100.0, &
    0.0001, 0.5000, 0.7900, 0.9100, 0.9500, 0.9500, 1.0000, 1.0000, &
    1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, &
    1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
  data data_r( 9), (data_e(ip, 9),ip=1,NP_DATA) /  150.0, &
    0.0001, 0.7700, 0.9300, 0.9700, 0.9700, 1.0000, 1.0000, 1.0000, &
    1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, &
    1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
  data data_r(10), (data_e(ip,10),ip=1,NP_DATA) /  200.0, &
    0.0001, 0.8700, 0.9600, 0.9800, 1.0000, 1.0000, 1.0000, 1.0000, &
    1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, &
    1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
  data data_r(11), (data_e(ip,11),ip=1,NP_DATA) /  300.0, &
    0.0001, 0.9700, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, &
    1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, &
    1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
  data data_r(12), (data_e(ip,12),ip=1,NP_DATA) / 1000.0, &
    0.0001, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, &
    1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, &
    1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /


  ! Use constant kernel if <icoagop> = I_COAGOP_CONST
  if( icoagop .eq. I_COAGOP_CONST )then
    ckernel(:,:,:,:,:) = ck0
  else
 
    if( icollec .eq. I_COLLEC_DATA )then
    
      ! Convert <data_r> from um to cm and take logarithm of <data_e>;
      ! however, we only want to do this once.
      !
      ! If we are using Open/MP, we only want one thread to do this
      ! operation once. This is a kludge, and this table should probably
      ! get set up a different way.
      !$OMP CRITICAL(CARMA_HALL)
      if (.not. init_data) then
        init_data = .TRUE.
        
        do i = 1, NR_DATA
          data_r(i) = data_r(i)/1.e4_f
          do ip = 1, NP_DATA
            data_e(ip,i) = log(data_e(ip,i))
          enddo
        enddo
      endif
      !$OMP END CRITICAL(CARMA_HALL)
    endif
   
    ! Loop over the grid
    do k = 1, NZ
    
      ! This is <rhoa> in cartesian coordinates.
      rhoa_cgs = rhoa(k) / (xmet(k)*ymet(k)*zmet(k))
  
      temp1 = BK*t(k)
      temp2 = 6._f*PI*rmu(k)
      
      do j1 = 1, NGROUP
        do j2 = j1, NGROUP
          use_vw(j1, j2) = is_grp_sulfate(j1) .and. is_grp_sulfate(j2)
        end do
      end do
  
      ! Loop over groups!
      do j1 = 1, NGROUP
        do j2 = 1, NGROUP
    
          if( icoag(j1,j2) .ne. 0 )then
  
            ! First particle
            do i1 = 1, NBIN
    
              r1 = r_wet(k,i1,j1)
              di = temp1*bpm(k,i1,j1)/(temp2*r1)
              gi  = sqrt( 8._f*temp1/(PI*rmass(i1,j1)) )
              rlbi = 8._f*di/(PI*gi)
              dti1= (2._f*r1 + rlbi)**3
              dti2= (4._f*r1*r1 + rlbi*rlbi)**1.5_f
              dti = 1._f/(6._f*r1*rlbi)
              dti = dti*(dti1 - dti2) - 2._f*r1
  
              !  Second particle
              do i2 = 1, NBIN
                r2  = r_wet(k,i2,j2)
                dj  = temp1*bpm(k,i2,j2)/(temp2*r2)
                gj  = sqrt( 8._f*temp1/(PI*rmass(i2,j2)) )
                rlbj = 8._f*dj/(PI*gj)
                dtj1= (2._f*r2 + rlbj)**3
                dtj2= (4._f*r2*r2 + rlbj*rlbj)**1.5_f
                dtj = 1._f/(6._f*r2*rlbj)
                dtj = dtj*(dtj1 - dtj2) - 2._f*r2
                
                !  Account for the charging effect of small particles (Van Der Waal's forces).  
                !  Set cstick to E_infinity/Eo, then multiply cbr kernel by Eo
                !  See Chan and Mozurkewich, J. Atmos. Sci., June 2001
                !  Only applicable to groups with sulfate elements    
                if (use_vw(j1,j2)) then
                  hp = ham / temp1 * (4._f * r1 * r2 / (r1 + r2)**2)
                  hpln = log(1._f + hp)
                  Enot = 1._f + vwa1 * hpln + vwa3 * hpln**3
                  Einf = 1._f + sqrt(hp / 3._f) / (1._f + vwb0*sqrt(hp)) + vwb1 * hpln + vwb3 * hpln**3
                  cstick_calc =  Einf / Enot
                else
                  cstick_calc = cstick
                end if

                !  First calculate thermal coagulation kernel
                rp  = r1 + r2
                dp  = di + dj
                gg  = sqrt(gi*gi + gj*gj)*cstick_calc
                delt= sqrt(dti*dti + dtj*dtj)
                term1 = rp/(rp + delt)
                term2 = 4._f*dp/(gg*rp)
  
                ! <cbr> is thermal (brownian) coagulation coefficient
                cbr = 4._f*PI*rp*dp/(term1 + term2)
  
                ! Determine indices of larger and smaller particles (of the pair)
                if (r2 .ge. r1) then
                  r_larg = r2
                  r_smal = r1
                  i_larg = i2
                  i_smal = i1
                  ig_larg = j2
                  ig_smal = j1
                  d_larg  = dj
                else
                  r_larg = r1
                  r_smal = r2
                  i_larg = i1
                  i_smal = i2
                  ig_larg = j1
                  ig_smal = j2
                  d_larg  = di
                endif
                
                ! Calculate enhancement of coagulation due to convective diffusion 
                ! as described in Pruppacher and Klett (Eqs. 17-12 and 17-14).
  
                ! Enhancement applies to larger particle.
                re_larg = re(k,i_larg,ig_larg)
  
                ! <pe> is Peclet number.
                pe  = re_larg*rmu(k) / (rhoa_cgs*d_larg)
                pe3 = pe**(1._f/3._f)
                
                ! <ccd> is convective diffusion coagulation coefficient
                if( re_larg .lt. 1._f )then
                  ccd = 0.45_f*cbr*pe3
                else 
                  ccd = 0.45_f*cbr*pe3*re_larg**(ONE/6._f)
                endif
  
                ! Next calculate gravitational collection kernel.  
                
                ! First evaluate collection efficiency <e>.
                if( icollec .eq. I_COLLEC_CONST )then
                  !   constant value
                  e_coll = grav_e_coll0
                else if( icollec .eq. I_COLLEC_FUCHS )then
                  ! Find maximum of Langmuir's formulation and Fuchs' value.
                  ! First calculate Langmuir's efficiency <e_langmuir>.
  
                  ! <sk> is stokes number.
                  ! <vfc_{larg,smal}> is the fallspeed in cartesian coordinates.!
                  vfc_smal = vf(k,i_smal,ig_smal) * zmet(k)
                  vfc_larg = vf(k,i_larg,ig_larg) * zmet(k)
  
                  sk = vfc_smal * (vfc_larg - vfc_smal) / (r_larg*GRAV)
     
                  if( sk .lt. 0.08333334_f )then
                    e1 = 0._f
                  else 
                    e1 = (sk/(sk + 0.25_f))**2
                  endif
     
                  if( sk .lt. 1.214_f )then
                    e3  = 0._f
                  else
                    e3  = 1._f/(1._f+.75_f*log(2._f*sk)/(sk-1.214_f))**2
                  endif
     
                  if( re_larg .lt. 1._f )then
                    e_langmuir = e3
                  else if( re_larg .gt. 1000._f )then
                    e_langmuir = e1
                  else if( re_larg .le. 1000._f )then
                    re60 = re_larg/60._f
                    e_langmuir = (e3  + re60*e1)/(1._f + re60)
                  endif
  
                  ! Next calculate Fuchs' efficiency (valid for r < 10 um).
                  pr = r_smal/r_larg
                  e_fuchs   = (pr/(1.414_f*(1. + pr)))**2
    
                  e_coll = max( e_fuchs, e_langmuir )
  
                else if( icollec .eq. I_COLLEC_DATA )then
  
                  ! Interpolate input data (from data statment at beginning of subroutine).
                  pr = r_smal/r_larg
  
                  ! First treat cases outside the data range
                  if( pr .lt. data_p(2) )then
  
                    ! Radius ratio is smaller than lowest nonzero ratio in input data --
                    ! use constant values (as in Beard and Ochs, 1984) if available,
                    ! otherwise use very small efficiencty
                    if( i2 .eq. i_larg )then
                      if( i2.eq.1 )then
                        e_coll = e_small
                      else
                        e_coll = e_coll2(i1,i2-1)
                      endif
                    else
                      if( i1.eq.1 )then
                        e_coll = e_small
                      else
                        e_coll = e_coll2(i1-1,i2)
                      endif
                    endif
    
                  elseif( r_larg .lt. data_r(1) )then
                    ! Radius of larger particle is smaller than smallest radius in input data -- 
                    ! assign very small efficiency.
                    e_coll = e_small
                  else
  
                    ! Both droplets are either within grid (interpolate) or larger
                    ! droplet is larger than maximum on grid (extrapolate) -- in both cases 
                    ! will interpolate on ratio of droplet radii.
  
                    ! Find <jp> such that data_p(jp) <= pr <= data_p(jp+1) and calculate
                    ! <pblni> = fractional distance of <pr> between points in <data_p> 
                    jp = NP_DATA
                    do jj = NP_DATA-1, 2, -1
                      if( pr .le. data_p(jj+1) ) jp = jj
                    enddo
                    
                    ! should not need this if-stmt
                    if( jp .lt. NP_DATA )then
                      pblni = (pr - data_p(jp)) / (data_p(jp+1) - data_p(jp))
                    else
                      ! nor this else-stmt
                      if (do_print) write(LUNOPRT, *) 'setupckern::ERROR NP_DATA < jp = ', jp
                      return
                    endif
  
                    if( r_larg .gt. data_r(NR_DATA) )then
  
                      ! Extrapolate on R and interpolate on p 
                      !
                      ! NOTE: This expression has a bugin it, since jr won't
                      ! be defined.
                      e_coll = (1._f-pblni)*data_e(jp  ,jr) + &
                               (   pblni)*data_e(jp+1,jr)
  
                    else
  
                      ! Find <jr> such that data_r(jr) <= r_larg <= data_r(jr+1) and calculate
                      ! <rblni> = fractional distance of <r_larg> between points in <data_r>
                      jr = NR_DATA
                      do jj = NR_DATA-1, 1, -1
                        if( r_larg .le. data_r(jj+1) ) jr = jj
                      enddo
                      rblni = (r_larg - data_r(jr)) / (data_r(jr+1) - data_r(jr))
   
                      ! Bilinear interpolation of logarithm of data.
                      e_coll = (1._f-pblni)*(1._f-rblni)*data_e(jp  ,jr  ) + &
                               (   pblni)*(1._f-rblni)*data_e(jp+1,jr  ) + &
                               (1._f-pblni)*(   rblni)*data_e(jp  ,jr+1) + &
                               (   pblni)*(   rblni)*data_e(jp+1,jr+1)  
  
                      ! (since data_e is logarithm of efficiencies)
                      term1 = (1._f-rblni)*(1._f-pblni)*data_e(jp,jr)
                      
                      if( jp .lt. NP_DATA )then
                        term2 = pblni*(1.-rblni)*data_e(jp+1,jr)
                      else
                        term2 = -100._f
                      endif
                      
                      if( jr .lt. NR_DATA )then
                        term3 = (1._f-pblni)*rblni*data_e(jp,jr+1)
                      else
                        term3 = -100._f
                      endif
                      
                      if( jr .lt. NR_DATA .and. jp .lt. NP_DATA )then
                        term4 = pblni*rblni*data_e(jp+1,jr+1)
                      else
                        term4 = -100._f
                      endif
  
                      e_coll = exp(term1 + term2 + term3 + term4)
                    endif
                  endif
  
                  e_coll2(i1,i2) = e_coll
                endif
  
                ! Now calculate coalescence efficiency from Beard and Ochs 
                ! (J. Geophys. Res. 89, 7165-7169, 1984).
                beta = log(r_smal*1.e4_f) + 0.44_f*log(r_larg*50._f)
                b_coal = 0.0946_f*beta - 0.319_f
                a_coal = sqrt(b_coal**2 + 0.00441)
                x_coal = (a_coal-b_coal)**(ONE/3._f) &
                       - (a_coal+b_coal)**(ONE/3._f)
                x_coal = x_coal + 0.459_f
  
                ! Limit extrapolated values to no less than 50% and no more than 100%
                x_coal = max(x_coal,.5_f)
                e_coal = min(x_coal,1._f)
  
                ! Now use coalescence efficiency and collision efficiency in definition
                ! of (geometric) gravitational collection efficiency <cgr>.
                vfc_1 = vf(k,i1,j1) * zmet(k)
                vfc_2 = vf(k,i2,j2) * zmet(k)
                cgr = e_coal * e_coll *  PI * rp**2 * abs( vfc_1 - vfc_2 )
  
                ! Long's (1974) kernel that only depends on size of larger droplet
  !                 if( r_larg .le. 50.e-4_f )then
  !                   cgr = 1.1e10_f * vol(i_larg,ig_larg)**2
  !                 else
  !                   cgr = 6.33e3_f * vol(i_larg,ig_larg)
  !                 endif
  
                ! Now combine all the coagulation and collection kernels into the
                ! overall kernel.
                ckernel(k,i1,i2,j1,j2) = cbr + ccd + cgr
                
                ! To avoid generation of large, non-physical hydrometeors by
                ! coagulation, cut down ckernel for large radii
  !                 if( ( r1 .gt. 0.18_f .and. r2 .gt. 10.e-4_f ) .or. &
  !                     ( r2 .gt. 0.18_f .and. r1 .gt. 10.e-4_f ) ) then
  !                   ckernel(k,i1,i2,j1,j2) = ckernel(k,i1,i2,j1,j2) / 1.e6_f
  !                 endif
  
              enddo    ! second particle bin
            enddo    ! first particle bin
          endif     ! icoag ne 0 
        enddo     ! second particle group
      enddo     ! first particle group
    enddo     ! vertical level
  endif     ! not constant
  
  ! return to caller with coagulation kernels evaluated.
  return
end
