module Sedimentation
  use micro_mg_utils, only: r8
  implicit none
  private
  save
  public :: sed_CalcCFL, sed_CalcFlux, sed_AdvanceOneStep
  integer, parameter, public  :: MG_LIQUID = 1
  integer, parameter, public  :: MG_ICE = 2
  integer, parameter, public  :: MG_RAIN = 3
  integer, parameter, public  :: MG_SNOW = 4
  real(r8), parameter, public :: CFL = 0.9_r8

contains

  subroutine sed_CalcCFL(q,n,cloud_frac,rho,pdel,nlev,i, &
    mg_type,deltat,g,an,rhof,alphaq,alphan,s1,s2,cfl,ncons,nnst,gamma_b_plus1,gamma_b_plus4)

    use micro_mg_utils, only: size_dist_param_liq, &
         size_dist_param_basic, &
         mg_ice_props, mg_liq_props, &
         mg_snow_props, mg_rain_props, &
         qsmall, bi, bc, br, bs

    implicit none
    real(r8), intent(in)            :: q(:,:), n(:,:)
    real(r8), intent(in)            :: rho(:,:), pdel(:,:), an(:,:), rhof(:,:)
    real(r8), intent(in)            :: cloud_frac(:,:), deltat, g
    integer,  intent(in)            :: nlev, i, mg_type
    real(r8), intent(in), optional  :: nnst, gamma_b_plus1, gamma_b_plus4
    logical, intent(in), optional   :: ncons
    real(r8), intent(out)           :: cfl, alphaq(:), alphan(:)
    real(r8), intent(out)           :: s1(:), s2(:)

    real(r8) :: qic(nlev), nic(nlev)
    real(r8) :: lambda_bounds(2)
    real(r8) :: lam(nlev), pgam(nlev), cq, cn, clam, tr, sqrtdisc, lambr(nlev)
    integer :: k
    logical :: limited(nlev)

    ! use quantity in cloud
    qic = q(i,:)/cloud_frac(i,:)
    nic = n(i,:)/cloud_frac(i,:)

    ! switch for specification of droplet and crystal number or ice cloud number
    if ((mg_type == MG_LIQUID  .or. mg_type == MG_ICE) .and. ncons) then
      nic = nnst/rho(i,:)
    end if

    if (mg_type == MG_LIQUID) then
      call size_dist_param_liq(mg_liq_props, qic, nic, rho(i,:), &
                               pgam, lam)

    else if (mg_type == MG_ICE) then
      call size_dist_param_basic(mg_ice_props, qic, nic, lam)

    else if (mg_type == MG_RAIN) then
      lambda_bounds(1) = mg_rain_props%lambda_bounds(1)**br
      lambda_bounds(2) = mg_rain_props%lambda_bounds(2)**br
      do k=1,nlev
        if (qic(k) > qsmall) then
          ! add upper limit to in-cloud number concentration to prevent
          ! numerical error
           nic(k) = min(nic(k), qic(k) / mg_rain_props%min_mean_mass)

          ! lambda^b = (c nic/qic)^(b/d)
          lambr(k) = (mg_rain_props%shape_coef * nic(k)/qic(k))**(br/mg_rain_props%eff_dim)
          ! check for slope
          ! adjust vars
          limited(k) = .false.
          if (lambr(k) < lambda_bounds(1)) then
            lambr(k) = lambda_bounds(1)
            nic(k) = lambr(k)**(mg_rain_props%eff_dim/br) * qic(k)/mg_rain_props%shape_coef
            limited(k) = .true.
          else if (lambr(k) > lambda_bounds(2)) then
            lambr(k) = lambda_bounds(2)
            nic(k) = lambr(k)**(mg_rain_props%eff_dim/br) * qic(k)/mg_rain_props%shape_coef
            limited(k) = .true.
          end if
        end if
      end do

    else if (mg_type == MG_SNOW) then
      call size_dist_param_basic(mg_snow_props, qic, nic, lam)
    else
      print *, "Invalid mg_type to mg2_sedimentation"
      stop
    end if

    ! Loop through levels to compute alphaq, alphan, and the CFL number
    ! initialize values to zero
    cfl = 0._r8
    alphaq = 0._r8
    alphan = 0._r8
    s1 = 0._r8
    s2 = 0._r8

    do k=1,nlev

      if(mg_type == MG_LIQUID) then
        if (qic(k) .ge. qsmall) then
          alphaq(k) = g*rho(i,k)*an(i,k)*gamma(4._r8+bc+pgam(k))/ &
               (lam(k)**bc*gamma(pgam(k)+4._r8))

          alphan(k) = g*rho(i,k)* &
               an(i,k)*gamma(1._r8+bc+pgam(k))/ &
               (lam(k)**bc*gamma(pgam(k)+1._r8))
          s1(k) = alphaq(k)
          s2(k) = alphan(k)
        end if

      else if (mg_type == MG_ICE) then
        if (qic(k) .ge. qsmall) then
          alphaq(k) = g*rho(i,k)*min(an(i,k)*gamma_b_plus4/(6._r8*lam(k)**bi), &
                1.2_r8*rhof(i,k))
          alphan(k) = g*rho(i,k)* &
                min(an(i,k)*gamma_b_plus1/lam(k)**bi,1.2_r8*rhof(i,k))
          s1(k) = alphaq(k)
          s2(k) = alphan(k)
        end if

      else if (mg_type == MG_RAIN) then
        if (qic(k) > qsmall) then
          cq = g*rho(i,k)*an(i,k)*gamma_b_plus4/6._r8
          cn = g*rho(i,k)*an(i,k)*gamma_b_plus1
          alphaq(k) = min(cq/lambr(k), 9.1_r8*g*rho(i,k)*rhof(i,k))
          alphan(k) = min(cn/lambr(k), 9.1_r8*g*rho(i,k)*rhof(i,k))

          ! obtain eigenvalues
          s1(k) = alphaq(k)
          s2(k) = alphan(k)
        end if

      else if (mg_type == MG_SNOW) then
        if (lam(k) .ge. qsmall) then
          alphaq(k) = g*rho(i,k)* &
            min(an(i,k)*gamma_b_plus4/(6._r8*lam(k)**bs),1.2_r8*rhof(i,k))
          alphan(k) = g*rho(i,k)* &
            min(an(i,k)*gamma_b_plus1/lam(k)**bs,1.2_r8*rhof(i,k))
          s1(k) = alphaq(k)
          s2(k) = alphan(k)
        end if
      end if

      ! Update CFL number
      cfl = max(cfl,s1(k)*deltat/pdel(i,k),s2(k)*deltat/pdel(i,k))

      ! TODO: detect shock formation

    end do

  end subroutine sed_CalcCFL

  subroutine sed_CalcFlux(q,n,alphaq,alphan,s1,s2,nlev,i,mg_type,fq,fn)
    implicit none
    real(r8), intent(in)            :: q(:,:), n(:,:), alphaq(:), alphan(:)
    real(r8), intent(in)            :: s1(:), s2(:)
    integer,  intent(in)            :: nlev, i, mg_type
    real(r8), intent(out)           :: fq(:,:), fn(:,:)

    real(r8) :: qtemp(0:nlev), ntemp(0:nlev)
    real(r8) :: a, b, c, d, dq, dn, coeff(2)
    integer :: k

    ! initialize values to zero
    fq = 0._r8
    fn = 0._r8



    ! initialize qtemp and ntemp with ghost cell
    qtemp(0) = 0._r8
    qtemp(1:nlev) = q(i,:)
    ntemp(0) = 0._r8
    ntemp(1:nlev) = n(i,:)

    ! Loop through levels to compute flux
    do k=1,nlev

      if(mg_type == MG_LIQUID) then
        fq(k,:) = alphaq(k)*qtemp(k)
        fn(k,:) = alphan(k)*ntemp(k)

      else if (mg_type == MG_ICE) then
        fq(k,:) = alphaq(k)*qtemp(k)
        fn(k,:) = alphan(k)*ntemp(k)

      else if (mg_type == MG_RAIN) then
        fq(k,:) = alphaq(k)*qtemp(k)
        fn(k,:) = alphan(k)*ntemp(k)

      else if (mg_type == MG_SNOW) then
        fq(k,:) = alphaq(k)*qtemp(k)
        fn(k,:) = alphan(k)*ntemp(k)
      end if

    end do

  end subroutine sed_CalcFlux


  subroutine sed_AdvanceOneStep(q,fq,n,fn,pdel,deltat,deltat_sed,nlev,i,mg_type,g, &
    qtend,ntend,prect,qsedtend,cloud_frac,qvlat,tlat,xxl,preci,qsevap)
    implicit none
    real(r8), intent(in)           :: pdel(:,:)
    real(r8), intent(in)              :: deltat, deltat_sed, g
    integer, intent(in)               :: nlev, i, mg_type
    real(r8), intent(inout)           :: q(:,:), qtend(:,:), n(:,:), ntend(:,:)
    real(r8), intent(inout)           :: fq(:,:), fn(:,:), prect(:)
    real(r8), intent(inout)           :: qsedtend(:,:)
    real(r8), intent(in), optional    :: cloud_frac(:,:), xxl
    real(r8), intent(inout), optional :: qvlat(:,:), tlat(:,:), preci(:)
    real(r8), intent(inout), optional :: qsevap(:,:)

    real(r8) :: ratio(nlev), deltaFluxQ, deltaFluxN, deltafluxQ_evap
    integer :: k

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
      deltafluxQ = (fq(k,1)-ratio(k)*fq(k-1,2))/pdel(i,k)
      deltafluxN = (fn(k,1)-ratio(k)*fn(k-1,2))/pdel(i,k)
      ! add fallout terms to eulerian tendencies
      qtend(i,k) = qtend(i,k) - deltafluxQ*deltat_sed/deltat
      ntend(i,k) = ntend(i,k) - deltafluxN*deltat_sed/deltat
      ! sedimentation tendency for output
      qsedtend(i,k) = qsedtend(i,k) - deltafluxQ*deltat_sed/deltat
      if (mg_type == MG_ICE .or. mg_type == MG_LIQUID) then
        ! add terms to to evap/sub of cloud water
        deltafluxQ_evap = (ratio(k)-1._r8)*fq(k-1,2)/pdel(i,k)
        qvlat(i,k) = qvlat(i,k) - deltafluxQ_evap*deltat_sed/deltat
        qsevap(i,k) = qsevap(i,k) - deltafluxQ_evap*deltat_sed/deltat
        tlat(i,k) = tlat(i,k) + deltafluxQ_evap*xxl*deltat_sed/deltat
      end if

      q(i,k) = q(i,k) - deltat_sed*deltafluxQ
      n(i,k) = n(i,k) - deltat_sed*deltafluxN

      if (q(i,k) < -1.d-16) then
        print *, q(i,k)
        stop "q negative"
      else if (n(i,k) < -1.d-16) then
        print *, n(i,k)
        stop "n negative"
      end if
    end do

    ! units below are m/s
    ! sedimentation flux at surface is added to precip flux at surface
    ! to get total precip (cloud + precip water) rate

    prect(i) = prect(i) + deltat_sed/deltat*fq(nlev,2)/g/1000._r8
    if (mg_type == MG_ICE .or. mg_type == MG_SNOW) then
      preci(i) = preci(i) + deltat_sed/deltat*fq(nlev,2)/g/1000._r8
    end if

  end subroutine

end module Sedimentation
