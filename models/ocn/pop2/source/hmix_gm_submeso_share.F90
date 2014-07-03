!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module hmix_gm_submeso_share 

!BOP
! !MODULE: hmix_gm_submeso_share 

! !DESCRIPTION:
!  This module contains routines for computing tracer and density differences for
!  use in the hmix_gm and mix_submeso routines. In addition, isopycnal slopes are
!  computed if necessary.

! !REVISION HISTORY:
!  SVN:$Id: hmix_gm_submeso_share.F90

! !USES:

   use registry
   use blocks
   use kinds_mod
   use grid
   use constants
   use state_mod
   use time_management
   use domain_size, only: km, nt
   use domain, only: nblocks_clinic

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: init_meso_mixing,   &
             tracer_diffs_and_isopyc_slopes 

!-----------------------------------------------------------------------
!
!  variables to save from one call to next
!
!-----------------------------------------------------------------------

   real (r8), dimension(:,:,:,:,:), allocatable, public :: &
      RX,RY,            &     ! Dx(rho), Dy(rho)
      TX,TY,TZ               ! tracer differences in each direction
   real (r8), dimension(:,:,:,:,:,:), allocatable, public :: &
      SLX, SLY                ! slope of isopycnal sfcs in x,y-direction
   real (r8), dimension(:,:,:,:), allocatable, public :: &
      RZ_SAVE                 ! Dz(rho)
   real (r8), dimension(:,:,:), allocatable, public :: &
      HXY,              &     ! dx/dy for y-z plane
      HYX                     ! dy/dx for x-z plane



!***********************************************************************

 contains

!***********************************************************************

! !IROUTINE: init_meso_mixing
! !INTERFACE:

   subroutine init_meso_mixing(hmix_tracer_itype,hmix_tracer_type_gm)

! !DESCRIPTION:
!  Initializes various submesoscale and GM mixing options and allocates necessary
!  space. Also computes some time-independent arrays.
!

   integer (int_kind) :: &
      iblock,            &  ! block index
      hmix_tracer_itype,hmix_tracer_type_gm
 
!-----------------------------------------------------------------------
!
!  allocate GM and submeso_flux arrays
!
!-----------------------------------------------------------------------    
    
   allocate (HXY (nx_block,ny_block,nblocks_clinic),    &
             HYX (nx_block,ny_block,nblocks_clinic))
   if(hmix_tracer_itype == hmix_tracer_type_gm) then
    allocate (SLX   (nx_block,ny_block,2,2,km,nblocks_clinic),  &
              SLY   (nx_block,ny_block,2,2,km,nblocks_clinic))
   endif
   allocate (TX(nx_block,ny_block,km,nt,nblocks_clinic),  &
             TY(nx_block,ny_block,km,nt,nblocks_clinic),  &
             TZ(nx_block,ny_block,km,nt,nblocks_clinic))

   allocate (RX(nx_block,ny_block,2,km,nblocks_clinic),  &
             RY(nx_block,ny_block,2,km,nblocks_clinic))

   allocate (RZ_SAVE(nx_block,ny_block,km,nblocks_clinic))

 
   HXY      = c0
   HYX      = c0
   SLX      = c0
   SLY      = c0   
   TX       = c0
   TY       = c0
   TZ       = c0
   RX       = c0
   RY       = c0
   RZ_SAVE  = c0
   

!-----------------------------------------------------------------------
!
!  register init_meso_mixing
!
!-----------------------------------------------------------------------

   call register_string ('init_meso_mixing')
   

!-----------------------------------------------------------------------
!
!  initialize various time-independent arrays
!
!-----------------------------------------------------------------------

  do iblock = 1,nblocks_clinic
     
     !*** Hyx = dy/dx for x-z plane

     HYX(:,:,iblock) = HTE(:,:,iblock) / HUS(:,:,iblock)

     !*** Hxy = dx/dy for y-z plane

     HXY(:,:,iblock) = HTN(:,:,iblock) / HUW(:,:,iblock)

  enddo
  
   
!-----------------------------------------------------------------------

   end subroutine init_meso_mixing

!-----------------------------------------------------------------------


! !IROUTINE: tracer_diffs_and_isopyc_slopes
! !INTERFACE:

   subroutine tracer_diffs_and_isopyc_slopes (TMIX, this_block)

! !DESCRIPTION:
!  Calculates common variables used in hmix_gm and mix_submeso.
!
! !INPUT PARAMETERS:


      real (r8), dimension(nx_block,ny_block,km,nt), intent(in) :: &
         TMIX                  ! tracers at all vertical levels
                               !   at mixing time level
      type (block), intent(in) :: &
         this_block            ! block info for this sub block


!-----------------------------------------------------------------------
!
!      local variables
!
!-----------------------------------------------------------------------

      integer (int_kind), parameter :: &
         ieast  = 1, iwest  = 2,       &
	 jnorth = 1, jsouth = 2
      integer (int_kind) :: &
         i,j,n,kk,k,        &! dummy loop counters
	 ktmp,              &! array indices
         kn, ks,            &! cyclic pointers for 2-level local arrays
         bid                 ! local block address for this sub block
      real (r8), dimension(nx_block,ny_block) :: &
         KMASK, KMASKE, KMASKN,   &! ocean mask
	 DRDT, DRDS                ! expansion coefficients d(rho)/dT,S
      real (r8), dimension(nx_block,ny_block,2) :: &
         TXP, TYP, TZP, TEMP
      real (r8), dimension(nx_block,ny_block) :: & 
         RZ                  ! Dz(rho)
      integer (int_kind), parameter :: &
         ktp = 1, kbt = 2     ! refer to the top and bottom halves of a 
                              ! grid cell, respectively

!-----------------------------------------------------------------------
!
!  register tracer_diffs_and_isopyc_slopes
!
!-----------------------------------------------------------------------

   call register_string ('tracer_diffs_and_isopyc_slopes')

!-----------------------------------------------------------------------
!
!     initialize various quantities
!
!-----------------------------------------------------------------------

      bid = this_block%local_id

      DRDT   = c0
      DRDS   = c0
      TXP    = c0
      TYP    = c0
      TZP    = c0
      TEMP   = c0


        kn = 1
        ks = 2

        do kk=1,km

          KMASK = merge(c1, c0, kk < KMT(:,:,bid))

!-----------------------------------------------------------------------
!
!     compute RX=Dx(rho) and RY=Dy(rho) for all vertical levels. 
!
!-----------------------------------------------------------------------

          if ( kk == 1 ) then

            do j=1,ny_block
              do i=1,nx_block
                if ( kk <= KMT(i,j,bid) .and. kk <= KMTE(i,j,bid) ) then
                  KMASKE(i,j) = c1
                else
                  KMASKE(i,j) = c0
                endif
                if ( kk <= KMT(i,j,bid) .and. kk <= KMTN(i,j,bid) ) then
                  KMASKN(i,j) = c1
                else
                  KMASKN(i,j) = c0
                endif
                TEMP(i,j,kn) = max(-c2, TMIX(i,j,kk,1))
              enddo
            enddo

            do j=1,ny_block
              do i=1,nx_block-1
                TXP(i,j,kn) = KMASKE(i,j) * (TEMP(i+1,j,kn)  &
                                            -TEMP(i,  j,kn))
              enddo
            enddo

            do j=1,ny_block-1
              do i=1,nx_block
                TYP(i,j,kn) = KMASKN(i,j) * (TEMP(i,j+1,kn)  &
                                            -TEMP(i,j,  kn))
              enddo
            enddo

            do n=1,nt
              do j=1,ny_block
                do i=1,nx_block-1
                  TX(i,j,kk,n,bid) = KMASKE(i,j)  &
                              * (TMIX(i+1,j,kk,n) - TMIX(i,j,kk,n))
                enddo
              enddo

              do j=1,ny_block-1
                do i=1,nx_block
                  TY(i,j,kk,n,bid) = KMASKN(i,j)  &
                              * (TMIX(i,j+1,kk,n) - TMIX(i,j,kk,n))
                enddo
              enddo
            enddo

!     D_T(rho) & D_S(rho) at level 1

            call state (kk, kk, TMIX(:,:,kk,1), TMIX(:,:,kk,2),  &
                        this_block, DRHODT=DRDT, DRHODS=DRDS) 

!     RX = Dx(rho) = DRDT*Dx(T) + DRDS*Dx(S)
!     RY = Dy(rho) = DRDT*Dy(T) + DRDS*Dy(S)

            RX(:,:,ieast ,kk,bid) = DRDT * TXP(:,:,kn)  &
                                  + DRDS * TX(:,:,kk,2,bid) 
            RY(:,:,jnorth,kk,bid) = DRDT * TYP(:,:,kn)  &
                                  + DRDS * TY(:,:,kk,2,bid) 

            do j=1,ny_block
              do i=2,nx_block
                RX(i,j,iwest,kk,bid) = DRDT(i,j) * TXP(i-1,j,kn)  &
                                     + DRDS(i,j) * TX (i-1,j,kk,2,bid)
              enddo
            enddo

            do j=2,ny_block
              do i=1,nx_block
                RY(i,j,jsouth,kk,bid) = DRDT(i,j) * TYP(i,j-1,kn)  &
                                      + DRDS(i,j) * TY (i,j-1,kk,2,bid)
              enddo
            enddo

          endif  ! end of kk == 1 if statement

!-----------------------------------------------------------------------
!
!     compute RZ=Dz(rho) and
!     SLX = RX / RZ = slope of isopycnal surfaces in x-direction
!     SLY = RY / RZ = slope of isopycnal surfaces in y-direction
!
!-----------------------------------------------------------------------

          if ( kk < km ) then

            TEMP(:,:,ks) = max(-c2, TMIX(:,:,kk+1,1))

            TZ(:,:,kk+1,1,bid) = TMIX(:,:,kk  ,1) - TMIX(:,:,kk+1,1)
            TZ(:,:,kk+1,2,bid) = TMIX(:,:,kk  ,2) - TMIX(:,:,kk+1,2) 
            TZP(:,:,ks) = TEMP(:,:,kn) - TEMP(:,:,ks)

!     RZ = Dz(rho) = DRDT*Dz(T) + DRDS*Dz(S)

            RZ = DRDT * TZP(:,:,ks) + DRDS * TZ (:,:,kk+1,2,bid) 
            RZ = min(RZ,-eps2)

            if (registry_match('init_gm')) then
              SLX(:,:,ieast ,kbt,kk,bid) = KMASK * RX(:,:,ieast ,kk,bid) / RZ
              SLX(:,:,iwest ,kbt,kk,bid) = KMASK * RX(:,:,iwest ,kk,bid) / RZ
              SLY(:,:,jnorth,kbt,kk,bid) = KMASK * RY(:,:,jnorth,kk,bid) / RZ
              SLY(:,:,jsouth,kbt,kk,bid) = KMASK * RY(:,:,jsouth,kk,bid) / RZ
            endif

!-----------------------------------------------------------------------
!
!     compute Dx(rho), Dy(rho) at level kk+1
!
!-----------------------------------------------------------------------

            KMASKE = merge(c1, c0, kk+1 <= KMT(:,:,bid) .and.  &
                           kk+1 <= KMTE(:,:,bid))
            KMASKN = merge(c1, c0, kk+1 <= KMT(:,:,bid) .and.  &
                           kk+1 <= KMTN(:,:,bid))

            do j=1,ny_block
              do i=1,nx_block-1
                TXP(i,j,ks) = KMASKE(i,j)*(TEMP(i+1,j,ks)  &
                                         - TEMP(i,j,ks)) 
              enddo
            enddo

            do j=1,ny_block-1
              do i=1,nx_block
                TYP(i,j,ks) = KMASKN(i,j)*(TEMP(i,j+1,ks)  &
                                         - TEMP(i,j,ks))
              enddo
            enddo

            do n=1,nt
              do j=1,ny_block
                do i=1,nx_block-1
                  TX(i,j,kk+1,n,bid) = KMASKE(i,j)  &
                            * (TMIX(i+1,j,kk+1,n) - TMIX(i,j,kk+1,n))
                enddo
              enddo

              do j=1,ny_block-1
                do i=1,nx_block
                  TY(i,j,kk+1,n,bid) = KMASKN(i,j)  &
                            * (TMIX(i,j+1,kk+1,n) - TMIX(i,j,kk+1,n))
                enddo
              enddo
            enddo

!     D_T(rho) & D_S(rho) at level kk+1

            call state (kk+1, kk+1, TMIX(:,:,kk+1,1),  &
                        TMIX(:,:,kk+1,2), this_block,  &
                        DRHODT=DRDT, DRHODS=DRDS)

            RX(:,:,ieast ,kk+1,bid) = DRDT * TXP(:,:,ks)  &
                                    + DRDS * TX(:,:,kk+1,2,bid) 
            RY(:,:,jnorth,kk+1,bid) = DRDT * TYP(:,:,ks)  &
                                    + DRDS * TY(:,:,kk+1,2,bid) 

            do j=1,ny_block
              do i=2,nx_block
                RX(i,j,iwest,kk+1,bid) = DRDT(i,j) * TXP(i-1,j,ks)  &
                                       + DRDS(i,j) * TX (i-1,j,kk+1,2,bid)
              enddo
            enddo

            do j=2,ny_block
              do i=1,nx_block
                RY(i,j,jsouth,kk+1,bid) = DRDT(i,j) * TYP(i,j-1,ks)  &
                                        + DRDS(i,j) * TY (i,j-1,kk+1,2,bid)
              enddo
            enddo

            RZ = DRDT * TZP(:,:,ks) + DRDS * TZ(:,:,kk+1,2,bid) 
            RZ_SAVE(:,:,kk+1,bid) = min(RZ,c0)
            RZ = min(RZ,-eps2)
	    
            if (registry_match('init_gm')) then

!-----------------------------------------------------------------------
!
!     compute slope of isopycnal surfaces at level kk+1
!
!-----------------------------------------------------------------------

              where ( kk+1 <= KMT(:,:,bid) )
                SLX(:,:,ieast, ktp,kk+1,bid) = RX(:,:,ieast ,kk+1,bid) / RZ
                SLX(:,:,iwest, ktp,kk+1,bid) = RX(:,:,iwest ,kk+1,bid) / RZ
                SLY(:,:,jnorth,ktp,kk+1,bid) = RY(:,:,jnorth,kk+1,bid) / RZ
                SLY(:,:,jsouth,ktp,kk+1,bid) = RY(:,:,jsouth,kk+1,bid) / RZ
              end where
            endif

!-----------------------------------------------------------------------
!
!     end of kk < km if block
!
!-----------------------------------------------------------------------

          endif

          ktmp = kn
          kn   = ks
          ks   = ktmp

        enddo   ! end of kk-loop
	

!-----------------------------------------------------------------------
!
     end subroutine tracer_diffs_and_isopyc_slopes
!
!***********************************************************************


 end module hmix_gm_submeso_share

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
