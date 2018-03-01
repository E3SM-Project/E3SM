!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glide_stress.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2013
!   Glimmer-CISM contributors - see AUTHORS file for list of contributors
!
!   This file is part of Glimmer-CISM.
!
!   Glimmer-CISM is free software: you can redistribute it and/or modify it
!   under the terms of the Lesser GNU General Public License as published
!   by the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   Glimmer-CISM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   Lesser GNU General Public License for more details.
!
!   You should have received a copy of the Lesser GNU General Public License
!   along with Glimmer-CISM. If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! *sfp* module to hold subroutines for calculation of stress components from converged, higher-order
! stress and effective viscosity fields. To be called at the end of HO vel calculation.

!TODO - Combine with glam_strs2?  Not sure it needs to be its own module.

module glide_stress

    use glimmer_paramets, only : dp
    use glide_types
    use parallel

    implicit none

    private
    public :: glide_calcstrsstr  

    contains

    subroutine glide_calcstrsstr( model )

        type(glide_global_type) :: model

        call calcstrsstr(model%general%ewn,  model%general%nsn,  model%general%upn, &
                         model%numerics%dew,       model%numerics%dns,              &
                         model%numerics%sigma,     model%numerics%stagsigma,        & 
                         model%geometry%thck,                                       &
                         model%geomderv%dusrfdew,   model%geomderv%dusrfdns,        &
                         model%geomderv%dthckdew,   model%geomderv%dthckdns,        &
                         model%velocity%uvel,       model%velocity%vvel,            &
                         model%stress%efvs,                                         &
                         model%stress%tau%xx,      model%stress%tau%yy,             &
                         model%stress%tau%xy,      model%stress%tau%scalar,         &
                         model%stress%tau%xz,      model%stress%tau%yz )

    end subroutine glide_calcstrsstr

    subroutine calcstrsstr( ewn,  nsn,  upn,  &
                            dew,        dns,       &
                            sigma,      stagsigma, &  
                            thck,                  &
                            dusrfdew,   dusrfdns,  &
                            dthckdew,   dthckdns,  &
                            uvel,       vvel,      &
                            efvs,                  &
                            tauxx,      tauyy,     &
                            tauxy,      tau,       &
                            tauxz,      tauyz )

!TODO - Remove scaling.
        use glimmer_paramets, only : len0, thk0
!!        use glimmer_horiz_bcs, only: horiz_bcs_unstag_scalar

        implicit none

        integer, intent(in) :: ewn, nsn, upn

        real(dp), intent(in) :: dew, dns 
        real(dp), intent(in), dimension(:)     :: sigma, stagsigma
        real(dp), intent(in), dimension(:,:,:) :: efvs, uvel, vvel
        real(dp), intent(in), dimension(:,:) :: thck, dusrfdew, &
                                                        dusrfdns, dthckdew, dthckdns

        real(dp), intent(out), dimension(:,:,:) :: tauxx, tauyy, tauxy, &
                                                           tauxz, tauyz, tau
        !*sfp* local vars
        integer :: ew, ns, up
        real(dp), parameter :: f1 = len0 / thk0
        real(dp) :: dew2, dew4, dns2, dns4
        real(dp), dimension(upn-1) :: dup, dupm        

        !*sfp* note that these are already defined and used in glam_strs2. If needed by PB&J 

        ! stress calc routine as well, consider moving these up-scope 
        dup(1:upn-1) = sigma(2:upn) - sigma(1:upn-1)
        dupm(:) = - 0.25_dp / dup(:)
        dew2 = 2.0_dp * dew; dns2 = 2.0_dp * dns        ! *sp* 2x the standard grid spacing
        dew4 = 4.0_dp * dew; dns4 = 4.0_dp * dns        ! *sp* 4x the standard grid spacing

!TODO - I think this loop should be over locally owned cells only.
        do ns = 2,nsn-1
            do ew = 2,ewn-1;

            if (thck(ew,ns) > 0.0_dp) then

                tauxz(:,ew,ns) = vertideriv(upn, hsum(uvel(:,ew-1:ew,ns-1:ns)), &
                                                  thck(ew,ns), dupm(1:upn-1))
                tauyz(:,ew,ns) = vertideriv(upn, hsum(vvel(:,ew-1:ew,ns-1:ns)), &
                                                  thck(ew,ns), dupm(1:upn-1))
                tauxx(:,ew,ns) = horizderiv(upn,  stagsigma,       &
                              sum(uvel(:,ew-1:ew,ns-1:ns),3), &
                              dew4, tauxz(:,ew,ns),           &
                              sum(dusrfdew(ew-1:ew,ns-1:ns)), &
                              sum(dthckdew(ew-1:ew,ns-1:ns)))
                tauyy(:,ew,ns) = horizderiv(upn,  stagsigma,       &
                              sum(vvel(:,ew-1:ew,ns-1:ns),2), &
                              dns4, tauyz(:,ew,ns),           &
                              sum(dusrfdns(ew-1:ew,ns-1:ns)), &
                              sum(dthckdns(ew-1:ew,ns-1:ns)))
                tauxy(:,ew,ns) = horizderiv(upn,  stagsigma,       &
                              sum(uvel(:,ew-1:ew,ns-1:ns),2), &
                              dns4, tauxz(:,ew,ns),           &
                              sum(dusrfdns(ew-1:ew,ns-1:ns)), &
                              sum(dthckdns(ew-1:ew,ns-1:ns))) + &
                              horizderiv(upn,  stagsigma,                &
                              sum(vvel(:,ew-1:ew,ns-1:ns),3), &
                              dew4, tauyz(:,ew,ns),           &
                              sum(dusrfdew(ew-1:ew,ns-1:ns)), &
                              sum(dthckdew(ew-1:ew,ns-1:ns)))
            else
                tauxz(:,ew,ns) = 0.0_dp
                tauyz(:,ew,ns) = 0.0_dp
                tauxx(:,ew,ns) = 0.0_dp
                tauyy(:,ew,ns) = 0.0_dp
                tauxy(:,ew,ns) = 0.0_dp
            end if

            end do
        end do

        !TODO - This should be over locally owned cells only.  Move into loop above?
        tauxz = f1 * efvs * tauxz     
        tauyz = f1 * efvs * tauyz     
        tauxx = 2.0_dp * efvs * tauxx 
        tauyy = 2.0_dp * efvs * tauyy 
        tauxy = efvs * tauxy          

        !*sfp* expanding this in terms of viscosity and velocity gradients, I've confirmed that 
        ! one gets the same thing as if one took Tau_eff = N_eff * Eps_eff, where Eps_eff is the 
        ! 1st order approx. to the 2nd strain-rate invariant (outlined in model description document).
        tau = sqrt(tauxz**2 + tauyz**2 + tauxx**2 + tauyy**2 + tauxx*tauyy + tauxy**2)

!TODO - I don't think these halo updates are needed.  
!       (If they are, they should be moved up to the glissade driver level.)

        call parallel_halo(tauxx)
!!        call horiz_bcs_unstag_scalar(tauxx)
        call parallel_halo(tauyy)
!!        call horiz_bcs_unstag_scalar(tauyy)
        call parallel_halo(tauxy)
!!        call horiz_bcs_unstag_scalar(tauxy)
        call parallel_halo(tauxz)
!!        call horiz_bcs_unstag_scalar(tauxz)
        call parallel_halo(tauyz)
!!        call horiz_bcs_unstag_scalar(tauyz)
        call parallel_halo(tau)
!!        call horiz_bcs_unstag_scalar(tau)
        return

    end subroutine calcstrsstr


    function vertideriv( upn, varb, thck, dupm )

        implicit none

        integer, intent(in) :: upn
        real(dp), intent(in), dimension(:) :: varb
        real(dp), intent(in) :: thck
        real(dp), intent(in), dimension(:) :: dupm            

        real(dp), dimension(size(varb)-1) :: vertideriv

        !*sfp* 'dupm' defined as -1/(2*del_sigma), in which case it seems like 
        !there should be a '-' in front of this expression ... or, negative sign
        !may be implicit in the vert indices ( "arb(2:upn) - varb(1:upn-1)" ) and
        !the fact that up=1 at the sfc and up=upn at the bed ??? 
        vertideriv(1:upn-1) = dupm(1:upn-1) * (varb(2:upn) - varb(1:upn-1)) / thck

        return

   end function vertideriv

   function horizderiv( upn,     stagsigma,   &
                         varb,    grid,        &
                         dvarbdz, dusrfdx, dthckdx)

        implicit none

        integer, intent(in) :: upn
        real(dp), dimension(:), intent(in) :: stagsigma
        real(dp), dimension(:,:), intent(in) :: varb
        real(dp), dimension(:), intent(in) :: dvarbdz
        real(dp), intent(in) :: dusrfdx, dthckdx, grid

        real(dp) :: horizderiv(size(varb,1)-1)

        ! *sfp* where does this factor of 1/4 come from ... averaging? 
        horizderiv = (varb(1:upn-1,2) + varb(2:upn,2) - varb(1:upn-1,1) - varb(2:upn,1)) / grid - &
                   dvarbdz * (dusrfdx - stagsigma * dthckdx) / 4.0_dp

        return

   end function horizderiv

   function hsum(inp)

      implicit none

      real(dp), dimension(:,:,:), intent(in) :: inp
      real(dp), dimension(size(inp,dim=1)) :: hsum

      hsum = sum(sum(inp(:,:,:),dim=3),dim=2)

      return

   end function hsum

end module glide_stress
