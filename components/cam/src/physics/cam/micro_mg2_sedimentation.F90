module micro_mg2_sedimentation
  use micro_mg_utils, only: r8
  implicit none
  private
  save

  public :: sed_CalcFallRate, sed_AdvanceOneStep
  integer, parameter, public  :: MG_LIQUID = 1
  integer, parameter, public  :: MG_ICE = 2
  integer, parameter, public  :: MG_RAIN = 3
  integer, parameter, public  :: MG_SNOW = 4
  real(r8), parameter, public :: CFL = 0.9_r8

  real(r8) :: brDeff_dim = -1._r8
  real(r8) :: eff_dimDbr = -1._r8
  real(r8) :: oneDshape_coeff = -1._r8
  real(r8) :: lamPbr_bounds(2) = (/ -1._r8, -1._r8 /)

contains

  subroutine sed_CalcFallRate(q,n,cloud_frac,rho,pdel,nlev,i, &
    mg_type,deltat,g,an,rhof,alphaq,alphan,cfl,ncons,nnst,gamma_b_plus1,gamma_b_plus4)
    use micro_mg_utils, only: size_dist_param_liq, &
                                           size_dist_param_basic, &
                                           mg_ice_props, mg_liq_props, &
                                           mg_rain_props, mg_snow_props, &
                                           qsmall, bi, bc, br, bs, &
                                           limiter_is_on
    implicit none
    real(r8), intent(in)            :: q(:,:), n(:,:)
    real(r8), intent(in)            :: rho(:,:), pdel(:,:), an(:,:), rhof(:,:)
    real(r8), intent(in)            :: cloud_frac(:,:), deltat, g
    integer, intent(in)             :: nlev, i, mg_type
    real(r8), intent(in), optional  :: nnst, gamma_b_plus1, gamma_b_plus4
    logical, intent(in), optional   :: ncons
    real(r8), intent(out)           :: alphaq(:), alphan(:)
    real(r8), intent(out)           :: cfl

    real(r8) :: qic(nlev), nic(nlev), lam(nlev), pgam(nlev), cq, cn, lamPbr(nlev)
    integer :: k, ngptl

    ! compute module lambda quantities for rain if not already computed
    if (brDeff_dim .eq. -1._r8) then
      brDeff_dim = br/mg_rain_props%eff_dim
      eff_dimDbr = mg_rain_props%eff_dim/br
      oneDshape_coeff = 1._r8/mg_rain_props%shape_coef
      lamPbr_bounds(1) = mg_rain_props%lambda_bounds(1)**br
      lamPbr_bounds(2) = mg_rain_props%lambda_bounds(2)**br
    end if

    ! use quantity in cloud
    qic = q(i,:)/cloud_frac(i,:)
    nic = n(i,:)/cloud_frac(i,:)

    ! switch for specification of droplet and crystal number or ice cloud number
    if ((mg_type == MG_LIQUID  .or. mg_type == MG_ICE) .and. ncons) then
      nic = nnst/rho(i,:)
    end if

    ! calculate lambda
    select case (mg_type)

      case (MG_ICE)
        call size_dist_param_basic(mg_ice_props, qic(:), nic(:), lam(:))

      case (MG_LIQUID)
        call size_dist_param_liq(mg_liq_props, qic(:), nic(:), rho(i,:), pgam(:), lam(:))

      case (MG_RAIN)
        do k=1,nlev
          if (qic(k) > qsmall) then
            ! add upper limit to in-cloud number concentration to prevent
            ! numerical error
            if (limiter_is_on(mg_rain_props%min_mean_mass)) then
              nic(k) = min(nic(k), qic(k) / mg_rain_props%min_mean_mass)
            end if
            ! lambda^b = (c nic/qic)^(b/d)
            lamPbr(k) = (mg_rain_props%shape_coef * nic(k)/qic(k))**brDeff_dim
            ! check for slope
            ! adjust vars
            if (lamPbr(k) < lamPbr_bounds(1)) then
              lamPbr(k) = lamPbr_bounds(1)
              nic(k) = lamPbr(k)**eff_dimDbr * qic(k)*oneDshape_coeff
            else if (lamPbr(k) > lamPbr_bounds(2)) then
              lamPbr(k) = lamPbr_bounds(2)
              nic(k) = lamPbr(k)**eff_dimDbr * qic(k)*oneDshape_coeff
            end if
          end if
        end do

      case (MG_SNOW)
        call size_dist_param_basic(mg_snow_props, qic(:), nic(:), lam(:))

      case default
        print *, "Invalid mg_type to mg2_sedimentation"
        stop

    end select

    ! Loop through levels to compute alphaq, alphan, and possibly eigenspaces
    select case (mg_type)

      case (MG_ICE)
        do k=1,nlev
          if (qic(k) .ge. qsmall) then
            alphaq(k) = g*rho(i,k)*min(an(i,k)*gamma_b_plus4/(6._r8*lam(k)**bi), &
                  1.2_r8*rhof(i,k))
            alphan(k) = g*rho(i,k)* &
                  min(an(i,k)*gamma_b_plus1/lam(k)**bi,1.2_r8*rhof(i,k))
          end if
        end do

      case (MG_LIQUID)
        do k=1,nlev
          if (qic(k) .ge. qsmall) then
            alphaq(k) = g*rho(i,k)*an(i,k)*gamma(4._r8+bc+pgam(k))/ &
                        (lam(k)**bc*gamma(pgam(k)+4._r8))

            alphan(k) = g*rho(i,k)* &
                        an(i,k)*gamma(1._r8+bc+pgam(k))/ &
                        (lam(k)**bc*gamma(pgam(k)+1._r8))
          end if
        end do

      case (MG_RAIN)
        do k=1,nlev
          if (qic(k) > qsmall) then
            cq = g*rho(i,k)*an(i,k)*gamma_b_plus4/6._r8
            cn = g*rho(i,k)*an(i,k)*gamma_b_plus1
            alphaq(k) = min(cq/lamPbr(k), 9.1_r8*g*rho(i,k)*rhof(i,k))
            alphan(k) = min(cn/lamPbr(k), 9.1_r8*g*rho(i,k)*rhof(i,k))
          end if
        end do

      case (MG_SNOW)
        do k=1,nlev
          if (lam(k) .ge. qsmall) then
            alphaq(k) = g*rho(i,k)* &
              min(an(i,k)*gamma_b_plus4/(6._r8*lam(k)**bs),1.2_r8*rhof(i,k))
            alphan(k) = g*rho(i,k)* &
              min(an(i,k)*gamma_b_plus1/lam(k)**bs,1.2_r8*rhof(i,k))
          end if
        end do

      case default
        print *, "Invalid mg_type to mg2_sedimentation"
        stop

    end select

    ! compute CFL number
    cfl = max(maxval(alphaq(:)*deltat/pdel(i,:)),maxval(alphan(:)*deltat/pdel(i,:)))

  end subroutine sed_CalcFallRate

  subroutine sed_AdvanceOneStep(q,n,alphaq,alphan,pdel,deltat,deltat_sed,nlev,&
    i,mg_type,g,qtend,ntend,prect,qsedtend,cloud_frac,qvlat,tlat,xxl,preci,qsevap)
    use micro_mg_utils, only: qsmall
    implicit none
    real(r8), intent(in)              :: pdel(:,:), alphaq(:), alphan(:)
    real(r8), intent(in)              :: deltat, deltat_sed, g
    integer, intent(in)               :: nlev, i, mg_type
    real(r8), intent(inout)           :: q(:,:), qtend(:,:), n(:,:), ntend(:,:)
    real(r8), intent(inout)           :: prect(:)
    real(r8), intent(inout)           :: qsedtend(:,:)
    real(r8), intent(in), optional    :: cloud_frac(:,:), xxl
    real(r8), intent(inout), optional :: qvlat(:,:), tlat(:,:), preci(:)
    real(r8), intent(inout), optional :: qsevap(:,:)

    real(r8) :: fq(0:nlev), fn(0:nlev), ratio(nlev), deltaFluxQ, deltaFluxN
    real(r8) :: deltafluxQ_evap
    integer :: k


    ! Compute flux
    fq(0) = 0._r8
    fq(1:nlev) = alphaq(:)*q(i,:)
    fn(0) = 0._r8
    fn(1:nlev) = alphan(:)*n(i,:)

    ! for cloud liquid and ice, if cloud fraction increases with height
    ! then add flux from above to both vapor and cloud water of current level
    ! this means that flux entering clear portion of cell from above evaporates
    ! instantly
    ! note: this is not an issue with precip, since we assume max overlap
    if (mg_type == MG_ICE .or. mg_type == MG_LIQUID) then
      ratio(1) = 1._r8
      ratio(2:nlev) = cloud_frac(i,2:nlev)/cloud_frac(i,1:(nlev-1))
    else
      ratio = 1._r8
    end if

    do k = 1,nlev
      ratio(k)=min(ratio(k),1._r8)
      deltafluxQ = (fq(k)-ratio(k)*fq(k-1))/pdel(i,k)
      deltafluxN = (fn(k)-ratio(k)*fn(k-1))/pdel(i,k)
      ! add fallout terms to eulerian tendencies
      qtend(i,k) = qtend(i,k) - deltafluxQ*deltat_sed/deltat
      ntend(i,k) = ntend(i,k) - deltafluxN*deltat_sed/deltat
      ! sedimentation tendency for output
      qsedtend(i,k) = qsedtend(i,k) - deltafluxQ*deltat_sed/deltat
      if (mg_type == MG_ICE .or. mg_type == MG_LIQUID) then
        ! add terms to to evap/sub of cloud water
        deltafluxQ_evap = (ratio(k)-1._r8)*fq(k-1)/pdel(i,k)
        qvlat(i,k) = qvlat(i,k) - deltafluxQ_evap*deltat_sed/deltat
        qsevap(i,k) = qsevap(i,k) - deltafluxQ_evap*deltat_sed/deltat
        tlat(i,k) = tlat(i,k) + deltafluxQ_evap*xxl*deltat_sed/deltat
      end if

      q(i,k) = q(i,k) - deltat_sed*deltafluxQ
      n(i,k) = n(i,k) - deltat_sed*deltafluxN

      if (q(i,k) < -1.d-10) then
         if (mg_type == MG_ICE) then
            write(iulog,1001) 'qi',i,k,q(i,k)
            stop "qi negative"
         else if (mg_type == MG_SNOW) then
            write(iulog,1001) 'qs',i,k,q(i,k)
            stop "qs negative"
         if (mg_type == MG_LIQUID) then
            write(iulog,1001) 'qc',i,k,q(i,k)
            stop "qc negative" 
        else if (mg_type == MG_RAIN) then
            write(iulog,1001) 'qr',i,k,q(i,k)
            stop "qr negative"
         end if

      else if (n(i,k) < -1.d-10) then
        print *, n(i,k)
        stop "n negative"
      end if
    end do

!if index i or k has more than 6 digits, will fail.
1001 format ( a2,'(', i6,',', i6,') = ',e12.3) 

    ! units below are m/s
    ! sedimentation flux at surface is added to precip flux at surface
    ! to get total precip (cloud + precip water) rate

    prect(i) = prect(i) + deltat_sed/deltat*fq(nlev)/g/1000._r8
    if (mg_type == MG_ICE .or. mg_type == MG_SNOW) then
      preci(i) = preci(i) + deltat_sed/deltat*fq(nlev)/g/1000._r8
    end if

  end subroutine

end module micro_mg2_sedimentation
