!     path:      $Source: /storm/rc1/cvsroot/rc/rrtmg_lw/src/mcica_subcol_gen_lw.f90,v $
!     author:    $Author: mike $
!     revision:  $Revision: 1.3 $
!     created:   $Date: 2007/08/28 22:38:11 $
!

      module mcica_subcol_gen_lw

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2006-2007, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------

! Purpose: Create McICA stochastic arrays for cloud physical or optical properties.
! Two options are possible:
! 1) Input cloud physical properties: cloud fraction, ice and liquid water
!    paths, ice fraction, and particle sizes.  Output will be stochastic
!    arrays of these variables.  (inflag = 1)
! 2) Input cloud optical properties directly: cloud optical depth, single
!    scattering albedo and asymmetry parameter.  Output will be stochastic
!    arrays of these variables.  (inflag = 0; longwave scattering is not
!    yet available, ssac and asmc are for future expansion)

! --------- Modules ----------

      use shr_kind_mod, only: r8 => shr_kind_r8
      use cam_abortutils,   only: endrun

!      use parkind, only : jpim, jprb
      use parrrtm, only : nbndlw, ngptlw
      use rrlw_con, only: grav
      use rrlw_wvn, only: ngb
      use rrlw_vsn

      implicit none
      private

! public interfaces/functions/subroutines
      public :: mcica_subcol_lw, generate_stochastic_clouds 

      contains

!------------------------------------------------------------------
! Public subroutines
!------------------------------------------------------------------

      subroutine mcica_subcol_lw(lchnk, ncol, nlay, icld, permuteseed, play, &
                       cldfrac, ciwp, clwp, rei, rel, tauc, cldfmcl, &
                       ciwpmcl, clwpmcl, reicmcl, relqmcl, taucmcl, rnglw, pergro)

! ----- Input -----
! Control
      integer, intent(in) :: lchnk           ! chunk identifier
      integer, intent(in) :: ncol            ! number of columns
      integer, intent(in) :: nlay            ! number of model layers
      integer, intent(in) :: icld            ! clear/cloud, cloud overlap flag
      integer, intent(in) :: permuteseed     ! if the cloud generator is called multiple times, 
                                                        ! permute the seed between each call.
                                                        ! between calls for LW and SW, recommended
                                                        ! permuteseed differes by 'ngpt'
      logical, intent(in) :: pergro

! Atmosphere
      real(kind=r8), intent(in) :: play(:,:)          ! layer pressures (mb) 
                                                        !    Dimensions: (ncol,nlay)

! Atmosphere/clouds - cldprop
      real(kind=r8), intent(in) :: cldfrac(:,:)       ! layer cloud fraction
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: tauc(:,:,:)        ! cloud optical depth
                                                        !    Dimensions: (nbndlw,ncol,nlay)
!      real(kind=r8), intent(in) :: ssac(:,:,:)       ! cloud single scattering albedo
                                                        !    Dimensions: (nbndlw,ncol,nlay)
!      real(kind=r8), intent(in) :: asmc(:,:,:)       ! cloud asymmetry parameter
                                                        !    Dimensions: (nbndlw,ncol,nlay)
      real(kind=r8), intent(in) :: ciwp(:,:)          ! cloud ice water path
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: clwp(:,:)          ! cloud liquid water path
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: rei(:,:)           ! cloud ice particle size
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: rel(:,:)           ! cloud liquid particle size
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: rnglw(:,:,:)           ! rand lw

! ----- Output -----
! Atmosphere/clouds - cldprmc [mcica]
      real(kind=r8), intent(out) :: cldfmcl(:,:,:)    ! cloud fraction [mcica]
                                                        !    Dimensions: (ngptlw,ncol,nlay)
      real(kind=r8), intent(out) :: ciwpmcl(:,:,:)    ! cloud ice water path [mcica]
                                                        !    Dimensions: (ngptlw,ncol,nlay)
      real(kind=r8), intent(out) :: clwpmcl(:,:,:)    ! cloud liquid water path [mcica]
                                                        !    Dimensions: (ngptlw,ncol,nlay)
      real(kind=r8), intent(out) :: relqmcl(:,:)      ! liquid particle size (microns)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(out) :: reicmcl(:,:)      ! ice partcle size (microns)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(out) :: taucmcl(:,:,:)    ! cloud optical depth [mcica]
                                                        !    Dimensions: (ngptlw,ncol,nlay)
!      real(kind=r8), intent(out) :: ssacmcl(:,:,:)   ! cloud single scattering albedo [mcica]
                                                        !    Dimensions: (ngptlw,ncol,nlay)
!      real(kind=r8), intent(out) :: asmcmcl(:,:,:)   ! cloud asymmetry parameter [mcica]
                                                        !    Dimensions: (ngptlw,ncol,nlay)

! ----- Local -----

! Stochastic cloud generator variables [mcica]
      integer, parameter :: nsubclw = ngptlw ! number of sub-columns (g-point intervals)
      integer :: km, im, nm                  ! loop indices

      real(kind=r8) :: pmid(ncol, nlay)               ! layer pressures (Pa) 
!      real(kind=r8) :: pdel(ncol, nlay)              ! layer pressure thickness (Pa) 
!      real(kind=r8) :: qi(ncol, nlay)                ! ice water (specific humidity)
!      real(kind=r8) :: ql(ncol, nlay)                ! liq water (specific humidity)


! Return if clear sky; or stop if icld out of range
      if (icld.eq.0) return
      if (icld.lt.0.or.icld.gt.3) then 
         call endrun('MCICA_SUBCOL: INVALID ICLD')
      endif 

! NOTE: For GCM mode, permuteseed must be offset between LW and SW by at least the number of subcolumns


! Pass particle sizes to new arrays, no subcolumns for these properties yet
! Convert pressures from mb to Pa

      reicmcl(:ncol,:nlay) = rei(:ncol,:nlay)
      relqmcl(:ncol,:nlay) = rel(:ncol,:nlay)
      pmid(:ncol,:nlay)    = play(:ncol,:nlay)*1.e2_r8

! Convert input ice and liquid cloud water paths to specific humidity ice and liquid components 

!      cwp =  (q * pdel * 1000.) / gravit)
!           = (kg/kg * kg m-1 s-2 *1000.) / m s-2
!           = (g m-2)
!
!      q  = (cwp * gravit) / (pdel *1000.)
!         = (g m-2 * m s-2) / (kg m-1 s-2 * 1000.)
!         =  kg/kg

!      do km = 1, nlay
!         qi(km) = (ciwp(km) * grav) / (pdel(km) * 1000._r8)
!         ql(km) = (clwp(km) * grav) / (pdel(km) * 1000._r8)
!      enddo

!  Generate the stochastic subcolumns of cloud optical properties for the longwave;
      call generate_stochastic_clouds (ncol, nlay, nsubclw, icld, pmid, cldfrac, clwp, ciwp, tauc, &
                               cldfmcl, clwpmcl, ciwpmcl, taucmcl, permuteseed, rnglw, pergro)!BSINGH

      end subroutine mcica_subcol_lw


!-------------------------------------------------------------------------------------------------
      subroutine generate_stochastic_clouds(ncol, nlay, nsubcol, icld, pmid, cld, clwp, ciwp, tauc, &
                                   cld_stoch, clwp_stoch, ciwp_stoch, tauc_stoch, changeSeed,rnglw, pergro)!BSINGH  
!-------------------------------------------------------------------------------------------------

  !----------------------------------------------------------------------------------------------------------------
  ! ---------------------
  ! Contact: Cecile Hannay (hannay@ucar.edu)
  ! 
  ! Original code: Based on Raisanen et al., QJRMS, 2004.
  ! 
  ! Modifications: Generalized for use with RRTMG and added Mersenne Twister as the default
  !   random number generator, which can be changed to the optional kissvec random number generator
  !   with flag 'irnd' below. Some extra functionality has been commented or removed.  
  !   Michael J. Iacono, AER, Inc., February 2007
  !
  ! Given a profile of cloud fraction, cloud water and cloud ice, we produce a set of subcolumns.
  ! Each layer within each subcolumn is homogeneous, with cloud fraction equal to zero or one 
  ! and uniform cloud liquid and cloud ice concentration.
  ! The ensemble as a whole reproduces the probability function of cloud liquid and ice within each layer 
  ! and obeys an overlap assumption in the vertical.   
  ! 
  ! Overlap assumption:
  !  The cloud are consistent with 4 overlap assumptions: random, maximum, maximum-random and exponential. 
  !  The default option is maximum-random (option 3)
  !  The options are: 1=random overlap, 2=max/random, 3=maximum overlap, 4=exponential overlap
  !  This is set with the variable "overlap" 
  !mji - Exponential overlap option (overlap=4) has been deactivated in this version
  !  The exponential overlap uses also a length scale, Zo. (real,    parameter  :: Zo = 2500. ) 
  ! 
  ! Seed:
  !  If the stochastic cloud generator is called several times during the same timestep, 
  !  one should change the seed between the call to insure that the subcolumns are different.
  !  This is done by changing the argument 'changeSeed'
  !  For example, if one wants to create a set of columns for the shortwave and another set for the longwave ,
  !  use 'changeSeed = 1' for the first call and'changeSeed = 2' for the second call 
  !
  ! PDF assumption:
  !  We can use arbitrary complicated PDFS. 
  !  In the present version, we produce homogeneuous clouds (the simplest case).  
  !  Future developments include using the PDF scheme of Ben Johnson. 
  !
  ! History file:
  !  Option to add diagnostics variables in the history file. (using FINCL in the namelist)
  !  nsubcol = number of subcolumns
  !  overlap = overlap type (1-3)
  !  Zo = length scale 
  !  CLOUD_S = mean of the subcolumn cloud fraction ('_S" means Stochastic)
  !  CLDLIQ_S = mean of the subcolumn cloud water
  !  CLDICE_S = mean of the subcolumn cloud ice 
  !
  ! Note:
  !   Here: we force that the cloud condensate to be consistent with the cloud fraction 
  !   i.e we only have cloud condensate when the cell is cloudy. 
  !   In CAM: The cloud condensate and the cloud fraction are obtained from 2 different equations 
  !   and the 2 quantities can be inconsistent (i.e. CAM can produce cloud fraction 
  !   without cloud condensate or the opposite).
  !---------------------------------------------------------------------------------------------------------------

      use mcica_random_numbers
! The Mersenne Twister random number engine
      use MersenneTwister, only: randomNumberSequence, &   
                                 new_RandomNumberSequence, getRandomReal
      type(randomNumberSequence) :: randomNumbers

! -- Arguments

      integer, intent(in) :: ncol            ! number of columns
      integer, intent(in) :: nlay            ! number of layers
      integer, intent(in) :: icld            ! clear/cloud, cloud overlap flag
      integer, intent(in) :: nsubcol         ! number of sub-columns (g-point intervals)      
      integer, optional, intent(in) :: changeSeed     ! allows permuting seed
      logical, intent(in) :: pergro

! Column state (cloud fraction, cloud water, cloud ice) + variables needed to read physics state 
      real(kind=r8), intent(in) :: pmid(:,:)          ! layer pressure (Pa)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: cld(:,:)           ! cloud fraction 
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: clwp(:,:)          ! cloud liquid water path
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: ciwp(:,:)          ! cloud ice water path
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: tauc(:,:,:)        ! cloud optical depth
                                                        !    Dimensions: (nbndlw,ncol,nlay)
!      real(kind=r8), intent(in) :: ssac(:,:,:)       ! cloud single scattering albedo
                                                        !    Dimensions: (nbndlw,ncol,nlay)
                                                        !   inactive - for future expansion
       real(kind=r8), intent(in) :: rnglw(:,:,:)           ! rand #lw
!      real(kind=r8), intent(in) :: asmc(:,:,:)       ! cloud asymmetry parameter
                                                        !    Dimensions: (nbndlw,ncol,nlay)
                                                        !   inactive - for future expansion

      real(kind=r8), intent(out) :: cld_stoch(:,:,:)  ! subcolumn cloud fraction 
                                                        !    Dimensions: (ngptlw,ncol,nlay)
      real(kind=r8), intent(out) :: clwp_stoch(:,:,:) ! subcolumn cloud liquid water path
                                                        !    Dimensions: (ngptlw,ncol,nlay)
      real(kind=r8), intent(out) :: ciwp_stoch(:,:,:) ! subcolumn cloud ice water path
                                                        !    Dimensions: (ngptlw,ncol,nlay)
      real(kind=r8), intent(out) :: tauc_stoch(:,:,:) ! subcolumn cloud optical depth
                                                        !    Dimensions: (ngptlw,ncol,nlay)
!      real(kind=r8), intent(out) :: ssac_stoch(:,:,:)! subcolumn cloud single scattering albedo
                                                        !    Dimensions: (ngptlw,ncol,nlay)
                                                        !   inactive - for future expansion
!      real(kind=r8), intent(out) :: asmc_stoch(:,:,:)! subcolumn cloud asymmetry parameter
                                                        !    Dimensions: (ngptlw,ncol,nlay)
                                                        !   inactive - for future expansion

! -- Local variables
      real(kind=r8) :: cldf(ncol,nlay)                ! cloud fraction 
    
! Mean over the subcolumns (cloud fraction, cloud water , cloud ice) - inactive
!      real(kind=r8) :: mean_cld_stoch(ncol, nlay)    ! cloud fraction 
!      real(kind=r8) :: mean_clwp_stoch(ncol, nlay)   ! cloud water
!      real(kind=r8) :: mean_ciwp_stoch(ncol, nlay)   ! cloud ice
!      real(kind=r8) :: mean_tauc_stoch(ncol, nlay)   ! cloud optical depth
!      real(kind=r8) :: mean_ssac_stoch(ncol, nlay)   ! cloud single scattering albedo
!      real(kind=r8) :: mean_asmc_stoch(ncol, nlay)   ! cloud asymmetry parameter

! Set overlap
      integer :: overlap                     ! 1 = random overlap, 2 = maximum/random,
                                                        ! 3 = maximum overlap, 
!      real(kind=r8), parameter  :: Zo = 2500._r8   ! length scale (m) 
!      real(kind=r8) :: zm(ncol,nlay)                 ! Height of midpoints (above surface)
!      real(kind=r8), dimension(nlay) :: alpha=0.0_r8    ! overlap parameter  

! Constants (min value for cloud fraction and cloud water and ice)
      real(kind=r8), parameter :: cldmin = 1.0e-80_r8     ! min cloud fraction
!      real(kind=r8), parameter :: qmin   = 1.0e-10_r8   ! min cloud water and cloud ice (not used)

! Variables related to random number and seed 
      integer :: irnd                        ! flag for random number generator
                                                        !  0 = kissvec
                                                        !  1 = Mersenne Twister

      real(kind=r8), dimension(nsubcol, ncol, nlay) :: CDF, CDF2      ! random numbers
      integer, dimension(ncol) :: seed1, seed2, seed3, seed4 ! seed to create random number (kissvec)
      real(kind=r8), dimension(ncol) :: rand_num      ! random number (kissvec)
      integer :: iseed                       ! seed to create random number (Mersenne Teister)
      real(kind=r8) :: rand_num_mt                    ! random number (Mersenne Twister)

! Flag to identify cloud fraction in subcolumns
      logical,  dimension(nsubcol, ncol, nlay) :: iscloudy   ! flag that says whether a gridbox is cloudy

! Indices
      integer :: ilev, isubcol, i, n         ! indices

!------------------------------------------------------------------------------------------ 

! Set randum number generator to use (0 = kissvec; 1 = mersennetwister)
      irnd = 0
!      irnd = 1

! Pass input cloud overlap setting to local variable
      overlap = icld

! ensure that cloud fractions are in bounds 
      cldf(:,:) = cld(:ncol,:nlay)
      where (cldf(:,:) < cldmin)
          cldf(:,:) = 0._r8
      end where

! ----- Create seed  --------
   
! Advance randum number generator by changeseed values
      if (irnd.eq.0) then   
! For kissvec, create a seed that depends on the state of the columns. Maybe not the best way, but it works.  
! Must use pmid from bottom four layers. 
         do i=1,ncol
            if (pmid(i,nlay).lt.pmid(i,nlay-1)) then 
               call endrun('MCICA_SUBCOL: KISSVEC SEED GENERATOR REQUIRES PMID FROM BOTTOM FOUR LAYERS.')
            endif 
            seed1(i) = (pmid(i,nlay) - int(pmid(i,nlay)))  * 1000000000
            seed2(i) = (pmid(i,nlay-1) - int(pmid(i,nlay-1)))  * 1000000000
            seed3(i) = (pmid(i,nlay-2) - int(pmid(i,nlay-2)))  * 1000000000
            seed4(i) = (pmid(i,nlay-3) - int(pmid(i,nlay-3)))  * 1000000000
          enddo
         do i=1,changeSeed
            call kissvec(seed1, seed2, seed3, seed4, rand_num)
         enddo
      elseif (irnd.eq.1) then
         randomNumbers = new_RandomNumberSequence(seed = changeSeed)
      endif 


! ------ Apply overlap assumption --------

! generate the random numbers  

      select case (overlap)

      case(1) 
! Random overlap
! i) pick a random value at every level
  
         if (irnd.eq.0) then 
            do isubcol = 1,nsubcol
               do ilev = 1,nlay
                  call kissvec(seed1, seed2, seed3, seed4, rand_num)  ! we get different random number for each level
                  CDF(isubcol,:,ilev) = rand_num
               enddo
            enddo
         elseif (irnd.eq.1) then
            do isubcol = 1, nsubcol
               do i = 1, ncol
                  do ilev = 1, nlay
                     rand_num_mt = getRandomReal(randomNumbers)
                     CDF(isubcol,i,ilev) = rand_num_mt
                  enddo
               enddo
             enddo
         endif

      case(2) 
! Maximum-Random overlap
! i) pick a random number for top layer.
! ii) walk down the column: 
!    - if the layer above is cloudy, we use the same random number than in the layer above
!    - if the layer above is clear, we use a new random number 

         if (irnd.eq.0) then 
            do isubcol = 1,nsubcol
               do ilev = 1,nlay
                  call kissvec(seed1, seed2, seed3, seed4, rand_num) 
                  CDF(isubcol,:,ilev) = rand_num !BSINGH -commented this line
                  if(pergro)CDF(isubcol,:,ilev) = rnglw(isubcol,1:ncol,ilev) !BSINGH - added this line
               enddo
            enddo
         elseif (irnd.eq.1) then
            do isubcol = 1, nsubcol
               do i = 1, ncol
                  do ilev = 1, nlay
                     rand_num_mt = getRandomReal(randomNumbers)
                     CDF(isubcol,i,ilev) = rand_num_mt
                  enddo
               enddo
             enddo
         endif

         do ilev = 2,nlay
            do i = 1, ncol
               do isubcol = 1, nsubcol
                  if (CDF(isubcol, i, ilev-1) > 1._r8 - cldf(i,ilev-1) ) then
                     CDF(isubcol,i,ilev) = CDF(isubcol,i,ilev-1) 
                  else
                     CDF(isubcol,i,ilev) = CDF(isubcol,i,ilev) * (1._r8 - cldf(i,ilev-1))
                  end if
               end do
            end do
         enddo
       
      case(3) 
! Maximum overlap
! i) pick the same random numebr at every level  

         if (irnd.eq.0) then 
            do isubcol = 1,nsubcol
               call kissvec(seed1, seed2, seed3, seed4, rand_num)
               do ilev = 1,nlay
                  CDF(isubcol,:,ilev) = rand_num
               enddo
            enddo
         elseif (irnd.eq.1) then
            do isubcol = 1, nsubcol
               do i = 1, ncol
                  rand_num_mt = getRandomReal(randomNumbers)
                  do ilev = 1, nlay
                     CDF(isubcol,i,ilev) = rand_num_mt
                  enddo
               enddo
             enddo
         endif

!    case(4) - inactive
!       ! Exponential overlap: weighting between maximum and random overlap increases with the distance. 
!       ! The random numbers for exponential overlap verify:
!       ! j=1   RAN(j)=RND1
!       ! j>1   if RND1 < alpha(j,j-1) => RAN(j) = RAN(j-1)
!       !                                 RAN(j) = RND2
!       ! alpha is obtained from the equation
!       ! alpha = exp(- (Zi-Zj-1)/Zo) where Zo is a characteristic length scale    


!       ! compute alpha
!       zm    = state%zm     
!       alpha(:, 1) = 0.
!       do ilev = 2,nlay
!          alpha(:, ilev) = exp( -( zm (:, ilev-1) -  zm (:, ilev)) / Zo)
!       end do
       
!       ! generate 2 streams of random numbers
!       do isubcol = 1,nsubcol
!          do ilev = 1,nlay
!             call kissvec(seed1, seed2, seed3, seed4, rand_num)
!             CDF(isubcol, :, ilev) = rand_num
!             call kissvec(seed1, seed2, seed3, seed4, rand_num)
!             CDF2(isubcol, :, ilev) = rand_num
!          end do
!       end do

!       ! generate random numbers
!       do ilev = 2,nlay
!          where (CDF2(:, :, ilev) < spread(alpha (:,ilev), dim=1, nCopies=nsubcol) )
!             CDF(:,:,ilev) = CDF(:,:,ilev-1) 
!          end where
!       end do

      end select

 
! -- generate subcolumns for homogeneous clouds -----
      do ilev = 1,nlay
         iscloudy(:,:,ilev) = (CDF(:,:,ilev) >= 1._r8 - spread(cldf(:,ilev), dim=1, nCopies=nsubcol) )
      enddo

! where the subcolumn is cloudy, the subcolumn cloud fraction is 1;
! where the subcolumn is not cloudy, the subcolumn cloud fraction is 0

      do ilev = 1,nlay
         do i = 1, ncol
            do isubcol = 1, ngptlw
               if (iscloudy(isubcol,i,ilev) ) then
                  cld_stoch(isubcol,i,ilev) = 1._r8
               else
                  cld_stoch(isubcol,i,ilev) = 0._r8
               endif
            end do
         end do
      enddo

! where there is a cloud, set the subcolumn cloud properties;
! incoming clwp, ciwp and tauc should be in-cloud quantites and not grid-averaged quantities

      do ilev = 1,nlay
         do i = 1, ncol
            do isubcol = 1, nsubcol
               if ( iscloudy(isubcol,i,ilev) .and. (cldf(i,ilev) > 0._r8) ) then
                  clwp_stoch(isubcol,i,ilev) = clwp(i,ilev)
                  ciwp_stoch(isubcol,i,ilev) = ciwp(i,ilev)
               else
                  clwp_stoch(isubcol,i,ilev) = 0._r8
                  ciwp_stoch(isubcol,i,ilev) = 0._r8
               end if
            end do
         end do
      enddo
      do ilev = 1,nlay
         do i = 1,ncol
            do isubcol = 1,nsubcol
               if ( iscloudy(isubcol,i,ilev) .and. (cldf(i,ilev) > 0._r8) ) then
                  n = ngb(isubcol)
                  tauc_stoch(isubcol,i,ilev) = tauc(n,i,ilev)
!                  ssac_stoch(isubcol,i,ilev) = ssac(n,i,ilev)
!                  asmc_stoch(isubcol,i,ilev) = asmc(n,i,ilev)
               else
                  tauc_stoch(isubcol,i,ilev) = 0._r8
!                  ssac_stoch(isubcol,i,ilev) = 1._r8
!                  asmc_stoch(isubcol,i,ilev) = 0._r8
               endif
            enddo
         enddo
      enddo

! -- compute the means of the subcolumns ---
!      mean_cld_stoch(:,:) = 0._r8
!      mean_clwp_stoch(:,:) = 0._r8
!      mean_ciwp_stoch(:,:) = 0._r8
!      mean_tauc_stoch(:,:) = 0._r8
!      mean_ssac_stoch(:,:) = 0._r8
!      mean_asmc_stoch(:,:) = 0._r8
!      do i = 1, nsubcol
!         mean_cld_stoch(:,:) =  cld_stoch(i,:,:) + mean_cld_stoch(:,:) 
!         mean_clwp_stoch(:,:) =  clwp_stoch( i,:,:) + mean_clwp_stoch(:,:) 
!         mean_ciwp_stoch(:,:) =  ciwp_stoch( i,:,:) + mean_ciwp_stoch(:,:) 
!         mean_tauc_stoch(:,:) =  tauc_stoch( i,:,:) + mean_tauc_stoch(:,:) 
!         mean_ssac_stoch(:,:) =  ssac_stoch( i,:,:) + mean_ssac_stoch(:,:) 
!         mean_asmc_stoch(:,:) =  asmc_stoch( i,:,:) + mean_asmc_stoch(:,:) 
!      end do
!      mean_cld_stoch(:,:) = mean_cld_stoch(:,:) / nsubcol
!      mean_clwp_stoch(:,:) = mean_clwp_stoch(:,:) / nsubcol
!      mean_ciwp_stoch(:,:) = mean_ciwp_stoch(:,:) / nsubcol
!      mean_tauc_stoch(:,:) = mean_tauc_stoch(:,:) / nsubcol
!      mean_ssac_stoch(:,:) = mean_ssac_stoch(:,:) / nsubcol
!      mean_asmc_stoch(:,:) = mean_asmc_stoch(:,:) / nsubcol

      end subroutine generate_stochastic_clouds


!------------------------------------------------------------------
! Private subroutines
!------------------------------------------------------------------

!-------------------------------------------------------------------------------------------------- 
      subroutine kissvec(seed1,seed2,seed3,seed4,ran_arr)
!-------------------------------------------------------------------------------------------------- 

! public domain code
! made available from http://www.fortran.com/
! downloaded by pjr on 03/16/04 for NCAR CAM
! converted to vector form, functions inlined by pjr,mvr on 05/10/2004

! safeguard against integer overflow, statement function changed to
! internal function by santos, Nov. 2012

! The  KISS (Keep It Simple Stupid) random number generator. Combines:
! (1) The congruential generator x(n)=69069*x(n-1)+1327217885, period 2^32.
! (2) A 3-shift shift-register generator, period 2^32-1,
! (3) Two 16-bit multiply-with-carry generators, period 597273182964842497>2^59
!  Overall period>2^123; 
!
      use shr_kind_mod, only: i8 => shr_kind_i8

      real(kind=r8), dimension(:), intent(inout)  :: ran_arr
      integer, dimension(:), intent(inout) :: seed1,seed2,seed3,seed4
      integer(i8) :: kiss
      integer :: i

      logical :: big_endian

      big_endian = (transfer(1_i8, 1) == 0)

      do i = 1, size(ran_arr)
         kiss = 69069_i8 * seed1(i) + 1327217885
         seed1(i) = low_byte(kiss)
         seed2(i) = m (m (m (seed2(i), 13), - 17), 5)
         seed3(i) = 18000 * iand (seed3(i), 65535) + ishft (seed3(i), - 16)
         seed4(i) = 30903 * iand (seed4(i), 65535) + ishft (seed4(i), - 16)
         kiss = int(seed1(i), i8) + seed2(i) + ishft (seed3(i), 16) + seed4(i)
         ran_arr(i) = low_byte(kiss)*2.328306e-10_r8 + 0.5_r8
      end do
    
      contains

        pure integer function m(k, n)
          integer, intent(in) :: k
          integer, intent(in) :: n

          m = ieor (k, ishft (k, n) )

        end function m

        pure integer function low_byte(i)
          integer(i8), intent(in) :: i

          if (big_endian) then
             low_byte = transfer(ishft(i,bit_size(1)),1)
          else
             low_byte = transfer(i,1)
          end if

        end function low_byte

      end subroutine kissvec

      end module mcica_subcol_gen_lw

