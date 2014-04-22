! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================
 module mpas_atmphys_driver_cloudiness
 use mpas_grid_types

 use mpas_atmphys_constants, only: ep_2
 use mpas_atmphys_vars
 
 implicit none
 private
 public:: allocate_cloudiness,   &
          deallocate_cloudiness, &
          driver_cloudiness

 integer,private:: i,j,k


!>\brief MPAS driver for parameterization of the diagnostic cloud fraction.
!>\author Laura D. Fowler (send comments to laura@ucar.edu).
!>\date 2013-05-01.
!>
!>\details
!>
!> subroutines in mpas_atmphys_driver_cloudiness:
!> ----------------------------------------------
!> allocate_cloudiness  : allocate local arrays for parameterization of diagnostic cloudiness.
!> deallocate_cloudiness: deallocate local arrays for parameterization of diagnostic cloudiness.
!> cloudiness_from_MPAS : initialize local arrays.
!> cloudiness_to_MPAS   : copy local arrays to MPAS arrays.
!> driver_cloudiness    : main driver (called from subroutine physics_driver).
!> calc_cldincidence    : calculates the cloud fraction as 0 or 1, depending on cloud condensates.
!> calc_cldfraction     : calculates the cloud fraction as a function of the relative humidity.
!>
!> add-ons and modifications to sourcecode:
!> ----------------------------------------
!>    * throughout the sourcecode, replaced all "var_struct" defined arrays by local pointers.
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2014-04-22.


 contains

!==================================================================================================
 subroutine allocate_cloudiness
!==================================================================================================

 if(.not.allocated(cldfrac_p) ) allocate(cldfrac_p(ims:ime,kms:kme,jms:jme) )

 end subroutine allocate_cloudiness

!==================================================================================================
 subroutine deallocate_cloudiness
!==================================================================================================

 if(allocated(cldfrac_p) ) deallocate(cldfrac_p )

 end subroutine deallocate_cloudiness

!==================================================================================================
 subroutine cloudiness_from_MPAS(diag_physics)
!==================================================================================================

!inout arguments:
 type(diag_physics_type),intent(inout):: diag_physics

!--------------------------------------------------------------------------------------------------

 do j = jts,jte
 do k = kts,kte
 do i = its,ite
    cldfrac_p(i,k,j) = 0.
 enddo
 enddo
 enddo

 end subroutine cloudiness_from_MPAS

!==================================================================================================
 subroutine cloudiness_to_MPAS(diag_physics)
!==================================================================================================

!inout arguments:
 type(diag_physics_type),intent(inout):: diag_physics

!local variables:
 real(kind=RKIND),dimension(:,:),pointer:: cldfrac

!--------------------------------------------------------------------------------------------------

 cldfrac => diag_physics % cldfrac % array

 do j = jts,jte
 do k = kts,kte
 do i = its,ite
    cldfrac(k,i) = cldfrac_p(i,k,j)
 enddo
 enddo
 enddo
 
 end subroutine cloudiness_to_MPAS

!==================================================================================================
 subroutine driver_cloudiness(diag_physics)
!==================================================================================================

!inout arguments:
 type(diag_physics_type),intent(inout):: diag_physics

!--------------------------------------------------------------------------------------------------
 write(0,*)
 write(0,*) '--- enter subroutine driver_cloudiness:'

!copy MPAS arrays to local arrays:
 call cloudiness_from_MPAS(diag_physics)

 cld_fraction_select: select case (trim(radt_cld_scheme))
    case("cld_incidence")

      !calculate the incidence of clouds:
      write(0,*) '--- enter subroutine calc_cldincidence:'
      call calc_cldincidence(cldfrac_p,qc_p,qi_p,f_qc,f_qi)
      write(0,*) '--- exit subroutine calc_cldincidence'

    case("cld_fraction")

      !calculate the cloud fraction based on the relative humidity:
      write(0,*) '--- enter subroutine calc_cldfraction:'
      call calc_cldfraction(cldfrac_p,t_p,pres_p,qv_p,qc_p,qi_p,qs_p)
      write(0,*) '--- exit subroutine calc_cldfraction'

    case default

 end select cld_fraction_select

!copy local arrays to MPAS grid:
 call cloudiness_to_MPAS(diag_physics)

 write(0,*) '--- exit subroutine driver_cloudiness'

 end subroutine driver_cloudiness

!==================================================================================================
 subroutine calc_cldincidence(cldfrac,qc,qi,f_qc,f_qi)
!==================================================================================================

!input arguments:
 logical,intent(in):: f_qc,f_qi
 real(kind=RKIND),intent(in),dimension(ims:ime,kms:kme,jms:jme):: qc,qi

!output arguments:
 real(kind=RKIND),intent(out),dimension(ims:ime,kms:kme,jms:jme):: cldfrac

!local variables:
 real(kind=RKIND),parameter:: thresh = 1.e-06

!--------------------------------------------------------------------------------------------------

 do j = jts,jte
 do k = kts,kte
 do i = its,ite
    cldfrac(i,k,j) = 0._RKIND
 enddo
 enddo
 enddo

 if(f_qc .and. f_qi) then
    do j = jts,jte
    do k = kts,kte
    do i = its,ite
       if(qc(i,k,j)+qi(i,k,j) .gt. thresh) cldfrac(i,k,j) = 1.0_RKIND
    enddo
    enddo
    enddo
 elseif(f_qc) then
    do j = jts,jte
    do k = kts,kte
    do i = its,ite
       if(qc(i,k,j) .gt. thresh) cldfrac(i,k,j) = 1.0_RKIND
    enddo
    enddo
    enddo
 endif
    
 end subroutine calc_cldincidence

!==================================================================================================
 subroutine calc_cldfraction(cldfrac,t_p,pres_p,qv,qc,qi,qs)
!==================================================================================================

!input arguments:
 real(kind=RKIND),intent(in),dimension(ims:ime,kms:kme,jms:jme):: qv,qc,qi,qs
 real(kind=RKIND),intent(in),dimension(ims:ime,kms:kme,jms:jme):: t_p,pres_p

!output arguments:
 real(kind=RKIND),intent(out),dimension(ims:ime,kms:kme,jms:jme):: cldfrac

!local variables:
 real(kind=RKIND),parameter:: alpha0  = 100.
 real(kind=RKIND),parameter:: gamma   = 0.49
 real(kind=RKIND),parameter:: qcldmin = 1.e-12
 real(kind=RKIND),parameter:: pexp    = 0.25
 real(kind=RKIND),parameter:: rhgrid  = 1.0

 real(kind=RKIND),parameter:: svp1  = 0.61078
 real(kind=RKIND),parameter:: svp2  = 17.2693882
 real(kind=RKIND),parameter:: svpi2 = 21.8745584
 real(kind=RKIND),parameter:: svp3  = 35.86
 real(kind=RKIND),parameter:: svpi3 = 7.66
 real(kind=RKIND),parameter:: svpt0 = 273.15

 real(kind=RKIND):: esi,esw,qvsi,qvsw
 real(kind=RKIND):: arg,denom,qcld,qvs,rhum,subsat,weight

!--------------------------------------------------------------------------------------------------

 do j = jts,jte
 do k = kts,kte
 do i = its,ite
    cldfrac(i,k,j) = 0._RKIND
 enddo
 enddo
 enddo
 
 do j = jts,jte
 do k = kts,kte
 do i = its,ite

!... calculation of the saturation mixing ratios over water and over ice (Murray, 1966):
    esw = 1000. * svp1 * exp(svp2 * (t_p(i,k,j) - svpt0) / (t_p(i,k,j) - svp3))
    esi = 1000. * svp1 * exp(svpi2 * (t_p(i,k,j) - svpt0) / (t_p(i,k,j) - svpi3))

    qvsw = ep_2 * esw / (pres_p(i,k,j) - esw)
    qvsi = ep_2 * esi / (pres_p(i,k,j) - esi)

    qcld = qc(i,k,j) + qi(i,k,j) + qs(i,k,j)
    if(qcld .lt. qcldmin) then
       weight = 0.
    else
       weight = (qi(i,k,j) + qs(i,k,j)) / qcld
    endif

    qvs = (1-weight) * qvsw + weight * qvsi
    rhum = qv(i,k,j) / qvs

    if(qcld .lt. qcldmin) then

       !assume that the cloud fraction is equal to 0. when the cloudy mixing ratio equals 0.
       cldfrac(i,k,j) = 0.

    elseif(rhum .ge. rhgrid) then
       !assume that the cloud fraction is equal to 1. when the relative humidity equal 100%.
       cldfrac(i,k,j) = 1.

    else
       !computation of the cloud fraction:
       subsat = max(1.e-10,rhgrid*qvs-qv(i,k,j))
       denom  = subsat**gamma
       arg    = max(-6.9,-alpha0*qcld/denom) ! exp(-6.9) = 0.001

       rhum = max(1.e-10,rhum)
       cldfrac(i,k,j) = (rhum/rhgrid)**pexp*(1.-exp(arg))
       if(cldfrac(i,k,j) .lt. 0.01) cldfrac(i,k,j) = 0.

    endif

!   if(qcld .ge. qcldmin) write(0,101) i,k,qcld,rhum,cldfrac(i,k,j) 

 enddo
 enddo
 enddo

 101 format(i8,1x,i3,3(1x,e15.8))

 end subroutine calc_cldfraction

!==================================================================================================
 end module mpas_atmphys_driver_cloudiness
!==================================================================================================
