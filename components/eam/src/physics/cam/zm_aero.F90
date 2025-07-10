module zm_aero
   !----------------------------------------------------------------------------
   ! Purpose: microphysics state structure definition and methods for ZM
   ! Original Author: Xialiang Song and Guang Zhang, June 2010
   !----------------------------------------------------------------------------
   use shr_kind_mod,     only: r8=>shr_kind_r8
   use ppgrid,           only: pcols, pver, pverp
   use cam_abortutils,   only: endrun
   use cam_logfile,      only: iulog
   use zm_aero_type,     only: zm_aero_t

   public :: zm_aero_init  ! aerosol stype initialization

!===================================================================================================
contains
!===================================================================================================

subroutine zm_aero_init(nmodes, nbulk, aero)
   !----------------------------------------------------------------------------
   ! Purpose: Initialize the zm_aero_t object for modal aerosols
   !----------------------------------------------------------------------------
   use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_mode_props, rad_cnst_get_aer_props
   use physconst,        only: pi
   !----------------------------------------------------------------------------
   ! Arguments
   integer,         intent(in   ) :: nmodes
   integer,         intent(in   ) :: nbulk
   type(zm_aero_t), intent(inout) :: aero
   !----------------------------------------------------------------------------
   ! Local variables
   character(len=*), parameter :: routine = 'zm_aero_init'
   character(len=20), allocatable :: aername(:)
   character(len=32) :: str32
   integer  :: iaer, l, m
   integer  :: nspecmx ! max number of species in a mode
   real(r8) :: sigmag
   real(r8) :: dgnumlo
   real(r8) :: dgnumhi
   real(r8) :: alnsg
   !----------------------------------------------------------------------------
   aero%nmodes = nmodes
   aero%nbulk  = nbulk

   if (nmodes > 0) then

      ! Initialize the modal aerosol information
      aero%scheme = 'modal'

      ! Get number of species in each mode, and find max.
      allocate(aero%nspec(aero%nmodes))
      nspecmx = 0
      do m = 1, aero%nmodes
         call rad_cnst_get_info(0, m, nspec=aero%nspec(m), mode_type=str32)
         nspecmx = max(nspecmx, aero%nspec(m))
         ! save mode index for specified mode types
         select case (trim(str32))
         case ('accum')
            aero%mode_accum_idx = m
         case ('aitken')
            aero%mode_aitken_idx = m
         case ('coarse')
            aero%mode_coarse_idx = m
         end select
      end do

      ! Check that required mode types were found
      if (aero%mode_accum_idx == -1 .or. &
          aero%mode_aitken_idx == -1 .or. &
          aero%mode_coarse_idx == -1) then
         write(iulog,*) routine//': ERROR required mode type not found - mode idx:', &
            aero%mode_accum_idx, aero%mode_aitken_idx, aero%mode_coarse_idx
         call endrun(routine//': ERROR required mode type not found')
      end if

      ! find indices for the dust and seasalt species in the coarse mode
      do l = 1, aero%nspec(aero%mode_coarse_idx)
         call rad_cnst_get_info(0, aero%mode_coarse_idx, l, spec_type=str32)
         select case (trim(str32))
         case ('dust')
            aero%coarse_dust_idx = l
         case ('seasalt')
            aero%coarse_nacl_idx = l
         case ('sulfate')
            aero%coarse_so4_idx = l
#if ( defined MODAL_AERO_4MODE_MOM  || defined MODAL_AERO_5MODE )
         case ('m-organic')
            aero%coarse_mom_idx  = l
#endif
#if ( defined RAIN_EVAP_TO_COARSE_AERO )
         case ('black-c')
            aero%coarse_bc_idx   = l
         case ('p-organic')
            aero%coarse_pom_idx  = l
         case ('s-organic')
            aero%coarse_soa_idx  = l
#endif
         end select
      end do

      ! Check that required modal species types were found
      if (aero%coarse_dust_idx == -1 .or. &
          aero%coarse_nacl_idx == -1 .or. &
          aero%coarse_so4_idx == -1) then
         write(iulog,*) routine//': ERROR required mode-species type not found - indicies:', &
            aero%coarse_dust_idx, aero%coarse_nacl_idx, aero%coarse_so4_idx
         call endrun(routine//': ERROR required mode-species type not found')
      end if

#if ( defined MODAL_AERO_4MODE_MOM  || defined MODAL_AERO_5MODE )
      if (aero%coarse_mom_idx == -1) then
         write(iulog,*) routine//': ERROR required mode-species type not found - indicies:', &
            aero%coarse_mom_idx
         call endrun(routine//': ERROR required mode-species type not found')
      end if
#endif

#if ( defined RAIN_EVAP_TO_COARSE_AERO )
      if (aero%coarse_bc_idx == -1 .or. &
          aero%coarse_pom_idx == -1 .or. &
          aero%coarse_soa_idx == -1) then
         write(iulog,*) routine//': ERROR required mode-species type not found - indicies:', &
            aero%coarse_bc_idx, aero%coarse_pom_idx, aero%coarse_soa_idx
         call endrun(routine//': ERROR required mode-species type not found')
      end if
#endif

      allocate( &
         aero%num_a(nmodes), &
         aero%mmr_a(nspecmx,nmodes), &
         aero%numg_a(pcols,pver,nmodes), &
         aero%mmrg_a(pcols,pver,nspecmx,nmodes), &
         aero%voltonumblo(nmodes), &
         aero%voltonumbhi(nmodes), &
         aero%specdens(nspecmx,nmodes), &
         aero%spechygro(nspecmx,nmodes), &
         aero%dgnum(nmodes), &
         aero%dgnumg(pcols,pver,nmodes) )

      do m = 1, nmodes
         ! Properties of modes
         call rad_cnst_get_mode_props( 0, m, sigmag=sigmag, dgnumlo=dgnumlo, dgnumhi=dgnumhi )
         alnsg               = log(sigmag)
         aero%voltonumblo(m) = 1 / ( (pi/6.0_r8)*(dgnumlo**3)*exp(4.5_r8*alnsg**2) )
         aero%voltonumbhi(m) = 1 / ( (pi/6.0_r8)*(dgnumhi**3)*exp(4.5_r8*alnsg**2) )
         ! save sigmag of aitken mode
         if (m == aero%mode_aitken_idx) aero%sigmag_aitken = sigmag
         ! Properties of modal species
         do l = 1, aero%nspec(m)
            call rad_cnst_get_aer_props(0, m, l, &
               density_aer = aero%specdens(l,m), &
               hygro_aer   = aero%spechygro(l,m))
         end do
      end do

   else if (nbulk > 0) then

      aero%scheme = 'bulk'

      ! Properties needed for BAM number concentration calcs.
      allocate( &
         aername(nbulk),                   &
         aero%num_to_mass_aer(nbulk),      &
         aero%mmr_bulk(nbulk),             &
         aero%mmrg_bulk(pcols,pver,nbulk)  )

      do iaer = 1, aero%nbulk
         call rad_cnst_get_aer_props(0, iaer, &
            aername         = aername(iaer), &
            num_to_mass_aer = aero%num_to_mass_aer(iaer) )
         ! Look for sulfate aerosol in this list (Bulk aerosol only)
         if (trim(aername(iaer)) == 'SULFATE') aero%idxsul   = iaer
         if (trim(aername(iaer)) == 'DUST1')   aero%idxdst1  = iaer
         if (trim(aername(iaer)) == 'DUST2')   aero%idxdst2  = iaer
         if (trim(aername(iaer)) == 'DUST3')   aero%idxdst3  = iaer
         if (trim(aername(iaer)) == 'DUST4')   aero%idxdst4  = iaer
         if (trim(aername(iaer)) == 'BCPHI')   aero%idxbcphi = iaer
      end do

   end if

end subroutine zm_aero_init

!===================================================================================================

end module zm_aero
