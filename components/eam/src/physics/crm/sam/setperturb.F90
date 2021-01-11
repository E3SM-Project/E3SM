module setperturb_mod
   use random_mod
   implicit none

contains

   subroutine setperturb(ncrms,icrm,iseed)

      ! Add random noise near the surface to help turbulence develop

      ! This surboutine has been updated for SPCAM5 (Minghuai.Wang@pnnl.gov, April, 2012).
      ! Now the random generator is seeded based on the global column id, which gets rid
      ! of the dependence of the SPCAM results on pcols.

      ! This module was updated to use a Mersenne Twister algorithm, because compiler deependent
      ! issues were identified with the intrinisic random number routines (e.g. random_number())
      ! Walter Hannah - LLNL - Mar 2018

      use vars
      use sgs,    only: setperturb_sgs
      use params, only: crm_rknd
      use RNG_MT

      implicit none
      integer, intent(in) :: ncrms,icrm
      integer, intent(in) :: iseed
      
      integer i,j,k
      real(crm_rknd)     :: rand_perturb              ! variable to hold random number generator output
      real(crm_rknd)     :: t02                       ! new average liquid static energy (LSE) for energy conservation scaling
      real(crm_rknd)     :: factor_xy                 ! 1/(nx*ny)
      real(crm_rknd)     :: perturb_k_scaling         ! scaling factor so perturbation magnitudes decrease with altitude
      integer, parameter :: perturb_num_layers  = 5   ! Number of levels to perturb
      integer, parameter :: perturb_t_magnitude = 1.0 ! perturbation LSE amplitube [K]

      factor_xy = 1.D0/real((nx*ny),crm_rknd)

      ! set the sub-grid scale (SGS) turbulence fields
      call setperturb_sgs(ncrms,icrm,0)  

      ! set the seed
      call RNG_MT_set_seed(iseed)

      !--------------------------------------------------------
      ! Apply random liquid static energy (LSE) perturbations
      !--------------------------------------------------------
      do k = 1,perturb_num_layers

         ! set perturb_k_scaling so that perturbation magnitude decreases with altitude
         perturb_k_scaling = real( perturb_num_layers+1-k ,crm_rknd) &
                            /real( perturb_num_layers ,crm_rknd)

         t02 = 0.0
         do j = 1,ny
            do i = 1,nx

               ! Generate a uniform random number in interval (0,1)
               call RNG_MT_gen_rand(rand_perturb)

               ! convert perturbation range from (0,1) to (-1,1)
               rand_perturb = 1.-2.*rand_perturb

               ! apply perturbation 
               t(icrm,i,j,k) = t(icrm,i,j,k) + perturb_t_magnitude &
                                             * rand_perturb * perturb_k_scaling
               
               ! Calculate new average LSE for energy conservation scaling below
               t02 = t02 + t(icrm,i,j,k)*factor_xy

            end do ! i
         end do ! j

         ! enforce energy conservation
         do j = 1,ny
            do i = 1,nx
               t(icrm,i,j,k) = t(icrm,i,j,k) *  t0(icrm,k)/t02
            end do ! i
         end do ! j

      end do ! k
      !--------------------------------------------------------
      !--------------------------------------------------------

   end subroutine setperturb

end module setperturb_mod
