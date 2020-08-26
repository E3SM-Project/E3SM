module precip_init_mod
  use params, only: asyncid
  use task_util_mod
  implicit none

contains

  subroutine precip_init(ncrms)
    ! Initialize precipitation related stuff
    use vars
    use micro_params
    use params
    use sat_mod
    implicit none
    integer, intent(in) :: ncrms
    real(crm_rknd) pratio, coef1, coef2,estw,esti,rrr1,rrr2
    real*4 :: gammafff
    external :: gammafff
    integer :: k,icrm

    gam3 = 3.
    gamr1 = 3.+b_rain
    gamr2 = (5.+b_rain)/2.
    gamr3 = 4.+b_rain
    gams1 = 3.+b_snow
    gams2 = (5.+b_snow)/2.
    gams3 = 4.+b_snow
    gamg1 = 3.+b_grau
    gamg2 = (5.+b_grau)/2.
    gamg3 = 4.+b_grau
    gam3 = gammafff(gam3)
    gamr1 = gammafff(gamr1)
    gamr2 = gammafff(gamr2)
    gamr3 = gammafff(gamr3)
    gams1 = gammafff(gams1)
    gams2 = gammafff(gams2)
    gams3 = gammafff(gams3)
    gamg1 = gammafff(gamg1)
    gamg2 = gammafff(gamg2)
    gamg3 = gammafff(gamg3)

    if(nint(gam3).ne.2) then
      if(masterproc)print*,'cannot compute gamma-function in precip_init. Exiting...'
      call task_abort
    end if


    !$acc parallel loop collapse(2)  async(asyncid)
    do k=1,nzm
      do icrm = 1 , ncrms
        pratio = sqrt(1.29 / rho(icrm,k))
        rrr1=393./(tabs0(icrm,k)+120.)*(tabs0(icrm,k)/273.)**1.5
        rrr2=(tabs0(icrm,k)/273.)**1.94*(1000./pres(icrm,k))
        estw = 100.*esatw_crm(tabs0(icrm,k))
        esti = 100.*esati_crm(tabs0(icrm,k))

        ! accretion by snow:
        coef1 = 0.25 * pi * nzeros * a_snow * gams1 * pratio/(pi * rhos * nzeros/rho(icrm,k) ) ** ((3+b_snow)/4.)
        coef2 = exp(0.025*(tabs0(icrm,k) - 273.15))
        accrsi(icrm,k) =  coef1 * coef2 * esicoef
        accrsc(icrm,k) =  coef1 * esccoef
        coefice(icrm,k) =  coef2

        ! evaporation of snow:
        coef1  =(lsub/(tabs0(icrm,k)*rv)-1.)*lsub/(therco*rrr1*tabs0(icrm,k))
        coef2  = rv*tabs0(icrm,k)/(diffelq*rrr2*esti)
        evaps1(icrm,k)  =  0.65*4.*nzeros/sqrt(pi*rhos*nzeros)/(coef1+coef2)/sqrt(rho(icrm,k))
        evaps2(icrm,k)  =  0.49*4.*nzeros*gams2*sqrt(a_snow/(muelq*rrr1))/(pi*rhos*nzeros)**((5+b_snow)/8.) / (coef1+coef2) * rho(icrm,k)**((1+b_snow)/8.)*sqrt(pratio)

        ! accretion by graupel:
        coef1 = 0.25*pi*nzerog*a_grau*gamg1*pratio/(pi*rhog*nzerog/rho(icrm,k))**((3+b_grau)/4.)
        coef2 = exp(0.025*(tabs0(icrm,k) - 273.15))
        accrgi(icrm,k) =  coef1 * coef2 * egicoef
        accrgc(icrm,k) =  coef1 * egccoef

        ! evaporation of graupel:
        coef1  =(lsub/(tabs0(icrm,k)*rv)-1.)*lsub/(therco*rrr1*tabs0(icrm,k))
        coef2  = rv*tabs0(icrm,k)/(diffelq*rrr2*esti)
        evapg1(icrm,k)  = 0.65*4.*nzerog/sqrt(pi*rhog*nzerog)/(coef1+coef2)/sqrt(rho(icrm,k))
        evapg2(icrm,k)  = 0.49*4.*nzerog*gamg2*sqrt(a_grau/(muelq*rrr1))/(pi * rhog * nzerog)**((5+b_grau)/8.) / (coef1+coef2) * rho(icrm,k)**((1+b_grau)/8.)*sqrt(pratio)

        ! accretion by rain:
        accrrc(icrm,k)=  0.25 * pi * nzeror * a_rain * gamr1 * pratio/(pi * rhor * nzeror / rho(icrm,k)) ** ((3+b_rain)/4.)* erccoef

        ! evaporation of rain:
        coef1  =(lcond/(tabs0(icrm,k)*rv)-1.)*lcond/(therco*rrr1*tabs0(icrm,k))
        coef2  = rv*tabs0(icrm,k)/(diffelq * rrr2 * estw)
        evapr1(icrm,k)  =  0.78 * 2. * pi * nzeror / &
        sqrt(pi * rhor * nzeror) / (coef1+coef2) / sqrt(rho(icrm,k))
        evapr2(icrm,k)  =  0.31 * 2. * pi  * nzeror * gamr2 * 0.89 * sqrt(a_rain/(muelq*rrr1))/(pi * rhor * nzeror)**((5+b_rain)/8.) / (coef1+coef2) * rho(icrm,k)**((1+b_rain)/8.)*sqrt(pratio)
      end do
    end do

  end subroutine precip_init


end module precip_init_mod
