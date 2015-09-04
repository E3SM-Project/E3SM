! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!!  This routine evaluates the derived mapping arrays and sets up
!!  the particle size bins.
!!
!!  @author Eric Jensen
!!  @ version Oct-1995
subroutine setupbins(carma, rc)

  ! types
  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carma_mod

  implicit none

  type(carma_type), intent(inout) :: carma   !! the carma object
  integer, intent(inout)          :: rc      !! return code, negative indicates failure

  ! Local declarations
  integer :: ielem, ibin, i, j, ix, iy, iz, ie, ig, ip, igrp, jgrp
  real(kind=f)  :: tmp_rhop(NBIN, NGROUP)
  real(kind=f)  :: vrfact
  real(kind=f)  :: cpi
  ! Local declarations needed for creation of fractal bin structure
  real(kind=f)  :: rf, rp
  real(kind=f)  :: vpor, upor, gamma, happel, perm, brinkman, epsil, omega

  !  Define formats
  !
  1 format(a,':  ',12i6)
  2 format(a,':  ',i6)
  3 format(a,':  ',f12.2)
  4 format(a,':  ',12f12.2)
  5 format(/,'Particle grid structure (setupbins):')
  6 format(a,':  ',1p12e12.3)
  7 format(a,':  ',12l6)


  !  Determine which elements are particle number concentrations
  !  <ienconc(igroup)> is the element corresponding to particle number 
  !  concentration in group <igroup>
  !
  igrp = 0
  do ielem = 1, NELEM
    if( itype(ielem) .eq. I_INVOLATILE .or. &
        itype(ielem) .eq. I_VOLATILE )then

      igrp = igrp + 1
      ienconc(igrp) = ielem
    endif
  enddo
  
  if( igrp .gt. NGROUP )then
    if (do_print) write(LUNOPRT,'(/,a)') 'CARMA_setupbin:: ERROR - bad itype array'
    rc = -1
    return
  endif

  !  Determine which group each element belongs to
  !  i.e., <igelem(ielem)> is the group to which element <ielem> belongs!
  igrp = 0
  do ielem = 1, NELEM
    if( itype(ielem) .eq. I_INVOLATILE .or. &
       itype(ielem) .eq. I_VOLATILE )then
      igrp = igrp + 1
    endif
    igelem(ielem) = igrp
  enddo

  !  Determine how many cores are in each group <ncore>.
  !  The core elements in a group are given by <icorelem(1:ncore,igroup)>.
  !
  !  Also evaluate whether or not second moment is used <if_sec_mom> for each group.
  ielem = 0
  
  do igrp = 1, NGROUP
  
    ncore(igrp) = 0
    if_sec_mom(igrp) = .false.
    imomelem(igrp) = 0
  
    do j = 1, nelemg(igrp)
  
      ielem = ielem + 1
  
      if( itype(ielem) .eq. I_COREMASS .or. &
          itype(ielem) .eq. I_VOLCORE )then
  
        ncore(igrp) = ncore(igrp) + 1
        icorelem(ncore(igrp),igrp) = ielem
  
      elseif( itype(ielem) .eq. I_CORE2MOM )then
  
        if_sec_mom(igrp) = .true.
        imomelem(igrp) = ielem
  
      endif
  
    enddo
  enddo

  !  Particle mass densities (NBIN for each group) -- the user might want
  !  to modify this (this code segment does not appear in setupaer subroutine
  !  because <igelem> is not defined until this subroutine).
  do ig = 1,NGROUP
    ie = ienconc(ig)
    do ibin = 1,NBIN
      tmp_rhop(ibin, ig) = rhoelem(ibin, ie)

        !  Set initial density of all hydrometeor groups to 1 such that nucleation
        !  mapping arrays are calculated correctly.
        !  or not
!           if( itype(ie) .ne. I_INVOLATILE ) then
!             rhop3(ixyz,ibin,ig) = 1.
!           endif
    enddo
  enddo
  
  !  Set up the particle bins.
  !  For each particle group, the mass of a particle in
  !  bin j is <rmrat> times that in bin j-1
  !
  !    rmass(NBIN,NGROUP)     =  bin center mass [g]
  !    r(NBIN,NGROUP)         =  bin mean (volume-weighted) radius [cm]
  !    vol(NBIN,NGROUP)       =  bin center volume [cm^3]
  !    dr(NBIN,NGROUP)        =  bin width in radius space [cm]
  !    dv(NBIN,NGROUP)        =  bin width in volume space [cm^3]
  !    dm(NBIN,NGROUP)        =  bin width in mass space [g]
  cpi = 4._f/3._f*PI

  do igrp = 1, NGROUP

    vrfact = ( (3._f/2._f/PI/(rmrat(igrp)+1._f))**(ONE/3._f) )* &
            ( rmrat(igrp)**(ONE/3._f) - 1._f )

    ! If rmassmin wasn't specified, then use rmin to determine the mass
    ! of the first bin.
    if (rmassmin(igrp) == 0._f) then
      rmassmin(igrp) = cpi*tmp_rhop(1,igrp)*rmin(igrp)**3
    else
      
      ! Just for internal consistency, recalculate rmin based on the rmass
      ! that is being used.
      rmin(igrp) = (rmassmin(igrp) / cpi / tmp_rhop(1,igrp)) ** (1._f / 3._f)
    end if
    
    do j = 1, NBIN
      rmass(j,igrp)   = rmassmin(igrp) * rmrat(igrp)**(j-1)
      rmassup(j,igrp) = 2._f*rmrat(igrp)/(rmrat(igrp)+1._f)*rmass(j,igrp)
      dm(j,igrp)      = 2._f*(rmrat(igrp)-1._f)/(rmrat(igrp)+1._f)*rmass(j,igrp)
      vol(j,igrp) = rmass(j,igrp) / tmp_rhop(j,igrp)
      r(j,igrp)   = ( rmass(j,igrp)/tmp_rhop(j,igrp)/cpi )**(ONE/3._f)
      rup(j,igrp) = ( rmassup(j,igrp)/tmp_rhop(j,igrp)/cpi )**(ONE/3._f)
      dr(j,igrp)  = vrfact*(rmass(j,igrp)/tmp_rhop(j,igrp))**(ONE/3._f)
      rlow(j,igrp) = rup(j,igrp) - dr(j,igrp)
 
      if (is_grp_fractal(igrp)) then
      ! fractal flag is true

        if (r(j,igrp) .le. rmon(igrp)) then   ! if the bin radius is less than the monomer size
                                     
          nmon(j,igrp) = 1.0_f
          rrat(j,igrp) = 1.0_f
          arat(j,igrp) = 1.0_f
          rprat(j,igrp) = 1.0_f
          df(j,igrp) = 3.0_f  ! Reset fractal dimension to 3 (this is a formality)

        else   ! if bin radius is greater than the monomer size

          rf = (1.0_f/falpha(igrp))**(1.0_f/df(j,igrp))*r(j,igrp)**(3.0_f/df(j,igrp))*rmon(igrp)**(1.0_f-3.0_f/df(j,igrp))
          nmon(j,igrp) = falpha(igrp)*(rf/rmon(igrp))**df(j,igrp)

          rrat(j,igrp) = rf/r(j,igrp)           
                                                                                         
          ! Calculate mobility radius for permeable aggregates
          ! using Vainshtein (2003) formulation     
          vpor = 1.0_f - (nmon(j,igrp))**(1.0_f-3.0_f/df(j,igrp))           ! Volume average porosity (eq. 3.2)
          upor = 1.0_f-(1.0_f - vpor)*sqrt(df(j,igrp)/3.0_f)                ! Uniform poroisty (eq. 3.10)
          gamma = (1.0_f - upor)**(1.0_f/3.0_f)
          happel = 2.0_f/(9.0_f*(1.0_f-upor))*   &                          ! Happel permeability model
                  (3.0_f-4.5_f*gamma+4.5_f*gamma**5.0_f-3.0_f*gamma**6.0_f)/  &
                  (3.0_f+2.0_f*gamma**5.0_f)    
          perm = happel*rmon(igrp)**2.0_f                                   ! Permeability (eq. 3.3)
          brinkman = nmon(j,igrp)**(1.0_f/df(j,igrp))*1.0_f/sqrt(happel)    ! Brinkman parameter (eq. 3.9) 
          epsil = 1.0_f - brinkman**(-1.)*tanh(brinkman)                    !
          omega = 2.0_f/3.0_f*epsil/(2.0_f/3.0_f+epsil/brinkman**2.0_f)     ! drag coefficient (eq. 2.7)
          rp = rf * omega
          rprat(j,igrp) = rp/r(j,igrp)

          arat(j,igrp) = (rprat(j,igrp) / rrat(j, igrp))**2.0_f
        endif
      else
         ! Not a fractal.
         nmon(j,igrp) = 1.0_f
         rprat(j,igrp) = 1.0_f
         df(j,igrp) = 3.0_f
      endif
   enddo
  enddo
  
  !  Evaluate differences between valuse of <rmass> in different bins.
  do igrp = 1, NGROUP
   do jgrp = 1, NGROUP
    do i = 1, NBIN
     do j = 1, NBIN
       diffmass(i,igrp,j,jgrp) = rmass(i,igrp) - rmass(j,jgrp)
     enddo
    enddo
   enddo
  enddo
  
  !  Report some initialization values
  if (do_print_init) then
    write(LUNOPRT,5)
    write(LUNOPRT,2) 'NGROUP ',NGROUP
    write(LUNOPRT,2) 'NELEM  ',NELEM
    write(LUNOPRT,2) 'NBIN   ',NBIN
    write(LUNOPRT,6) 'Massmin',(rmassmin(i),i=1,NGROUP)
    write(LUNOPRT,4) 'Mrat   ',(rmrat(i),i=1,NGROUP)
    write(LUNOPRT,1) 'nelemg ',(nelemg(i),i=1,NGROUP)
    write(LUNOPRT,1) 'itype  ',(itype(i),i=1,NELEM)
    write(LUNOPRT,1) 'ienconc',(ienconc(i),i=1,NGROUP)
    write(LUNOPRT,1) 'igelem ',(igelem(i),i=1,NELEM)
    write(LUNOPRT,1) 'ncore  ',(ncore(i),i=1,NGROUP)
    write(LUNOPRT,7) 'fractal',(is_grp_fractal(i),i=1,NGROUP)
  end if
 
  !  Return to caller with particle grid initialized
  return
end
