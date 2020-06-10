!**********************************************************************************  
! This computer software was prepared by Battelle Memorial Institute, hereinafter
! the Contractor, under Contract No. DE-AC05-76RL0 1830 with the Department of 
! Energy (DOE). NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY,
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
! MOSAIC module: see module_mosaic_driver.F for information and terms of use
!**********************************************************************************  
	module module_data_mosaic_asect


	implicit none


!-----------------------------------------------------------------------
!
!   The variables in this module provide a means of organizing and accessing
!   aerosol species in the "chem" array by their chemical component, 
!   size bin (or mode), "type", and "phase"
!
!   Their purpose is to allow flexible coding of process modules, 
!   compared to "hard-coding" using the chem array p_xxx indices
!   (e.g., p_so4_a01, p_so4_a02, ...; p_num_a01, ...)
!
!-----------------------------------------------------------------------
!
!   rce & sg 2004-dec-03 - added phase and type capability,
!	which changed this module almost completely
!
!-----------------------------------------------------------------------
!
!   maxd_atype = maximum allowable number of aerosol types
!   maxd_asize = maximum allowable number of aerosol size bins
!   maxd_acomp = maximum allowable number of chemical components
!	in each aerosol size bin
!   maxd_aphase = maximum allowable number of aerosol phases 
!	(gas, cloud, ice, rain, ...)
!
!   ntype_aer = number of aerosol types
!	The aerosol type will allow treatment of an externally mixed 
!	aerosol.  The current MOSAIC code has only 1 type, with the implicit
!	assumption of internal mixing.  Eventually, multiple types 
!	could treat fresh primary BC/OC, fresh SO4 from nucleation, 
!	aged BC/OC/SO4/... mixture, soil dust, sea salt, ... 
!
!   nphase_aer = number of aerosol phases
!
!   ai_phase = phase (p) index for interstitial (unactivated) aerosol particles
!   cw_phase = phase (p) index for aerosol particles in cloud water
!   ci_phase = phase (p) index for aerosol particles in cloud ice
!   rn_phase = phase (p) index for aerosol particles in rain
!   sn_phase = phase (p) index for aerosol particles in snow
!   gr_phase = phase (p) index for aerosol particles in graupel
!   [Note:  the value of "xx_phase" will be between 1 and nphase_aer 
!	for phases that are active in a simulation.  The others
!	will have non-positive values.]
!
!   nsize_aer(t) = number of aerosol size bins for aerosol type t
!
!   ncomp_aer(t) = number of "regular" chemical components for aerosol type t
!
!   massptr_aer(c,s,t,p) = the position/index in the chem array for mixing- 
!	ratio for chemical component c, size bin s, type t, and phase p.
!
!   numptr_aer(s,t,p) = the position/index in the chem array for mixing- 
!	ratio of particle number for size bin s, type t, and phase p.
!
!-----------------------------------------------------------------------
!
!   dens_aer(c,t) = dry density (g/cm^3) of aerosol chemical component 
!	c of type t
!   [Note:  dens_aer(c,t) == dens_mastercomp_aer(mastercompptr_aer(c,t))
!	The dens_mastercomp_aer is used in some initialization routines.
!	The dens_aer is used in most other places because of convenience.]
!
!-----------------------------------------------------------------------
!
!   volumlo_sect(s,t) = 1-particle volume (cm^3) at lower boundary of section m
!   volumhi_sect(s,t) = 1-particle volume (cm^3) at upper boundary of section m
!   volumcen_sect(s,t)= 1-particle volume (cm^3) at "center" of section m
!
!   [Note:  the "center" values are defined as follows:
!       volumcen_sect == 0.5*(volumlo_sect + volumhi_sect)
!                     == (pi/6) * (dcen_sect**3) ]
!
!
!-----------------------------------------------------------------------

	integer, save :: maxd_atype = 0
	integer, save :: maxd_asize = 0
	integer, save :: maxd_acomp = 0
	integer, save :: maxd_aphase = 0 

	integer, save :: ai_phase = -999888777
	integer, save :: cw_phase = -999888777
!	integer, save :: ci_phase = -999888777
!	integer, save :: rn_phase = -999888777
!	integer, save :: sn_phase = -999888777
!	integer, save :: gr_phase = -999888777

	integer, save :: ntype_aer = 0 ! number of types
	integer, save :: nphase_aer = 0 ! number of phases

        integer, allocatable ::   &
           nsize_aer (:),     & ! number of size bins
           ncomp_aer (:),     & ! number of chemical components
           massptr_aer( :, :, :, :), &
                ! index for mixing ratio
           numptr_aer( :, :, :) ! index for the number mixing ratio

        real, allocatable  ::  dens_aer( :, :)   ! aerosol density
        real, allocatable  ::  hygro_aer( :, :)  ! hygroscopicity 
        real, allocatable  ::  sigmag_aer(:, :)  ! geometric standard deviation for aerosol

!  added by Yang Zhang
        real, allocatable ::   &
          volumhi_sect( :, :),   &
          volumlo_sect( :, :),   &
          dcen_sect( :, : ),   &
          dlo_sect( :, : ),   &
          dhi_sect( :, : ) 

! flag for aerosols +++mhwang
        logical, allocatable :: is_aerosol(:) ! true if field is aerosol (any phase)

        integer, allocatable ::   &
                iphase_of_aerosol(:), isize_of_aerosol(:), itype_of_aerosol(:),   &
                inmw_of_aerosol(:), laicwpair_of_aerosol(:)

	end module module_data_mosaic_asect
