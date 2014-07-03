! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! This routine sets up mapping arrays for coagulation. It only computes varaibles that
!! are independent of the model state. The calculation of factors needed for coagulation
!! that depend on state are calculated in <i>setupckern</i>.
!!
!! @author Eric Jensen
!! @ version Oct-1995
subroutine setupcoag(carma, rc)

  ! types
  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carma_mod

  implicit none

  type(carma_type), intent(inout) :: carma   !! the CARMA object
  integer, intent(inout)         :: rc       !! return code, negative indicates failure

  ! Local declarations
  integer      :: ielem, isolto, icompto, igto, ig, iepart 
  integer      :: icompfrom, ic, iecore  
  integer      :: isolfrom
  integer      :: igrp, jg, i, j , ipair
  real(kind=f) :: rmsum
  integer      :: ibin
  real(kind=f) :: rmkbin
  integer      :: kb, ncg 
  real(kind=f) :: rmk 
  logical      :: fill_bot           ! used for filling <icoag>
  integer      :: irow, icol
  logical      :: isCoag
  integer      :: igtest
  real(kind=f) :: pkernl, pkernu


  ! NOTE: Moved this section from from setupckern.f, since it is not dependent on the
  ! model's state.
  !
  ! Fill <icoag>, maintaining diagonal symmetry
  ! -------------------------------------------
  ! Fill bottom of matrix if non-zero term(s) in upper half;
  ! also check for non-zero, non-matching, non-diagonal terms.
  fill_bot = .true.
  do irow = 2, NGROUP
    do icol = 1, irow-1
      if( icoag(irow,icol) .ne. 0 )then
        fill_bot = .false.
        if( icoag(icol,irow) .ne. 0 .and. &
            icoag(icol,irow) .ne. icoag(irow,icol) )then
          if (do_print) write(LUNOPRT, *) 'setupcoag::ERROR bad icoag array'
          rc = -1
          return
        endif
      endif
    enddo
  enddo

  do ig = 2, NGROUP
    do jg = 1, ig-1
      if( fill_bot )then
        irow = ig
        icol = jg
      else
        irow = jg
        icol = ig
      endif
      icoag(irow,icol) = icoag(icol,irow)
    enddo
  enddo

  ! Initialize <icoagelem> with zeros
  do ielem = 1,NELEM
    do ig = 1,NGROUP
      icoagelem(ielem,ig) = 0
      icoagelem_cm(ielem,ig) = 0
    enddo
  enddo

  ! For each element <ielem> and each group <ig>, determine which element in <ig>
  ! contributes to production  in <ielem>: <icoagelem(ielem,ig)>.
  ! If no elements in <ig> are transfered into element <ielem> during coagulation,
  ! then set <icoagelem(ielem,ig)> to 0.
  do ielem = 1,NELEM
    isolto = isolelem(ielem)           ! target solute type
    icompto = icomp(ielem)             ! target element compound
    igto = igelem(ielem)               ! target group

    do ig = 1, NGROUP                 ! source group
      ! source particle number concentration element
      iepart = ienconc(ig)

      ! source element compound
      icompfrom = icomp(iepart) 
      
      ! Check to see if the target group is produced by coagulation of any
      ! group with the source group.
      isCoag = .FALSE.
      
      do igtest = 1, NGROUP
        if (icoag(ig, igtest) .eq. igto .or. icoag(igtest, ig) .eq. igto) then
          isCoag = .TRUE.
        endif
      end do
      
      ! Only find the source production element if the group igto can
      ! be produced by coagulation from group ig.
      if (isCoag) then
      
        ! If <ig> only has no cores, then the only way to make particles
        ! would be if the one element <iepart> is the same type as the
        ! source.
        if( ncore(ig) .eq. 0 ) then
        
          if( icompfrom .eq. icompto )then
            icoagelem(ielem,ig) = iepart
          endif
        else
        
          ! Search the elements in the group to see if one has the same
          ! type as the source.
          
          ! First check the particle number concentration element of the group.
          !
          ! NOTE: No matter what else happens, you need to adjust the total
          ! particle mass.
          if( icompfrom .eq. icompto )then
            icoagelem(ielem,ig) = iepart
          else
          
            ! Now check the other cores for a match.
            do ic = 1,ncore(ig)
              iecore = icorelem(ic,ig)       ! absolute element number of core
              icompfrom = icomp(iecore)    ! source element compound
              
              if( icompfrom .eq. icompto ) then
                        
                ! For core second moment elements, we need additional pairs of source
                ! elements c  to account for core moment production due to products
                ! of source particle core mass.
                if( itype(ielem) .eq. I_CORE2MOM )then
                  icoagelem_cm(ielem,ig) = iecore
                  icoagelem(ielem,ig) = imomelem(ig)
                else
                  icoagelem(ielem,ig) = iecore
                endif
              endif
            enddo
          endif
        endif

        ! If <ielem> is a core mass type and <ig> is a pure CN group and the
        ! solutes don't match, then set <icoagelem> to zero to make sure no
        ! coag production occurs.
        if( itype(ielem) .eq. I_COREMASS .and. &
            itype(ienconc(ig)).eq. I_INVOLATILE &
            .and. ncore(ig) .eq. 0 ) then
          isolfrom = isolelem(ienconc(ig))
          if( isolfrom .ne. isolto ) then
            icoagelem(ielem,ig) = 0
          endif
        endif
        
        ! If there is a source and this is a multi-component group,
        ! then we need to make sure that the particle concentration
        ! of the group also gets updated, since this keeps track of
        ! the total mass.
        if (icoagelem(ielem,ig) .ne. 0) then
          if (ncore(igto) .ne. 0 .and. ielem .ne. ienconc(igto)) then
            icoagelem(ienconc(igto), ig) = iepart
          endif
        endif
        
      endif
    enddo          ! end of (ig = 1, NGROUP)
  enddo            ! end of (ielem = 1,NELEM)
  

  ! Coagulation won't work properly if any of the elements are produced by
  ! items that come later in the element list than themselves. Report an
  ! error if that is the case.
  do ielem = 1, NELEM
    do ig = 1, NGROUP
      if (icoagelem(ielem, ig) .gt. ielem) then
         if (do_print) write(LUNOPRT, '(a,i3,a,i3,a)') &
           'setupcoag::ERROR For coagulation, element (', &
           icoagelem(ielem,ig), ') must come before (', ielem, &
           ') in the element list.'
         rc = -1
         return
      endif
    enddo
  enddo
  
   
  !  Calculate lower bin <kbin> which coagulated particle goes into
  !  and make sure it is less than <NBIN>+1
  !
  !  Colliding particles come from group <ig>, bin <i> and group <jg>, bin <j>
  !  Resulting particle lands in group <igrp>, between <ibin> and <ibin> + 1
  do igrp = 1, NGROUP
    do ig = 1, NGROUP
      do jg = 1, NGROUP
        do i = 1, NBIN
          do j = 1, NBIN

            rmsum = rmass(i,ig) + rmass(j,jg)

            do ibin = 1, NBIN-1
              if( rmsum .ge. rmass(ibin,igrp) .and. rmsum .lt. rmass(ibin+1,igrp) ) then
                kbin(igrp,ig,jg,i,j) = ibin
              endif
            enddo

            ibin = NBIN
            if( rmsum .ge. rmass(ibin,igrp) ) kbin(igrp,ig,jg,i,j) = NBIN
          enddo
        enddo
      enddo
    enddo
  enddo
  
  ! Calculate partial loss fraction
  !
  ! This fraction is needed because when a particle in bin <i> collides
  ! with a particle in bin <j> resulting in a particle whose mass falls
  ! between <i> and <i>+1, only partial loss occurs from bin <i>.
  !
  ! Since different particle groups have different radius grids, this
  ! fraction is a function of the colliding groups and the resulting group.
  do igrp = 1, NGROUP
    do ig = 1, NGROUP
      do jg = 1, NGROUP

        if( igrp .eq. icoag(ig,jg) ) then

          do i = 1, NBIN
            do j = 1,NBIN
              volx(igrp,ig,jg,i,j) = 1.

              if(kbin(igrp,ig,jg,i,j).eq.i) then

                ibin = kbin(igrp,ig,jg,i,j)
                rmkbin = rmass(ibin,igrp)
                volx(igrp,ig,jg,i,j) = 1. - &
                   (rmrat(igrp)*rmkbin-rmass(i,ig)-rmass(j,jg)) &
                   /(rmrat(igrp)*rmkbin-rmkbin)* &
                   rmass(i,ig)/(rmass(i,ig) + rmass(j,jg))
              endif
            enddo
          enddo
        endif
      enddo
    enddo
  enddo
  
  ! Calculate mapping functions that specify sets of quadruples
  ! (group pairs and bin pairs) that contribute to production
  ! in each bin. Mass transfer from <ig,i> to <igrp,ibin> occurs due to
  ! collisions between particles in <ig,i> and particles in <jg,j>.
  ! 2 sets of quadruples must be generated:
  !    low: k = ibin and (k != i or ig != igrp)  and  icoag(ig,jg) = igrp
  !     up: k+1 = ibin        and  icoag(ig,jg) = igrp
  !
  ! npair#(igrp,ibin) is the number of pairs in each set (# = l,u)
  ! i#, j#, ig#, and jg# are the bin pairs and group pairs in each
  ! set (# = low, up)
  do igrp = 1, NGROUP
    do ibin = 1, NBIN

      npairl(igrp,ibin) = 0
      npairu(igrp,ibin) = 0

      do ig = 1, NGROUP
        do jg = 1, NGROUP
          do i = 1, NBIN
            do j = 1, NBIN
              kb = kbin(igrp,ig,jg,i,j)
              ncg = icoag(ig,jg)
    
              if( kb+1.eq.ibin .and. ncg.eq.igrp ) then
                npairu(igrp,ibin) = npairu(igrp,ibin) + 1
                iup(igrp,ibin,npairu(igrp,ibin)) = i
                jup(igrp,ibin,npairu(igrp,ibin)) = j
                igup(igrp,ibin,npairu(igrp,ibin)) = ig
                jgup(igrp,ibin,npairu(igrp,ibin)) = jg
              endif
    
              if( kb.eq.ibin .and. ncg.eq.igrp .and. (i.ne.ibin .or. ig.ne.igrp) ) then
                npairl(igrp,ibin) = npairl(igrp,ibin) + 1
                ilow(igrp,ibin,npairl(igrp,ibin)) = i
                jlow(igrp,ibin,npairl(igrp,ibin)) = j
                iglow(igrp,ibin,npairl(igrp,ibin)) = ig
                jglow(igrp,ibin,npairl(igrp,ibin)) = jg
              endif
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo


! NOTE: Split ckernel out of pkernel, so that it can be made independent of model state.
! It also reduces the size of the tables and should improve the intialization time.

!  Calculate variables needed in routine coagp.f
  do igrp = 1, NGROUP
    do jg = 1, NGROUP
      do ig = 1, NGROUP

        if( igrp .eq. icoag(ig,jg) ) then
        
          do j = 1, NBIN
            do i = 1, NBIN

              ibin = kbin(igrp,ig,jg,i,j)
              rmk = rmass(ibin,igrp)
              rmsum = rmass(i,ig) + rmass(j,jg)

              pkernl = (rmrat(igrp)*rmk - rmsum) / (rmrat(igrp)*rmk - rmk)
                        
              pkernu = (rmsum - rmk) / (rmrat(igrp)*rmk - rmk)

              if( ibin .eq. NBIN )then
                pkernl = rmsum / rmass(ibin,igrp)
                pkernu = 0._f
              endif
  
              pkernel(i,j,ig,jg,igrp,1) = pkernu * rmass(i,ig)/rmsum
              pkernel(i,j,ig,jg,igrp,2) = pkernl * rmass(i,ig)/rmsum
              pkernel(i,j,ig,jg,igrp,3) = pkernu * rmk*rmrat(igrp)/rmsum
              pkernel(i,j,ig,jg,igrp,4) = pkernl * rmk/rmsum
              pkernel(i,j,ig,jg,igrp,5) = pkernu * ( rmk*rmrat(igrp)/rmsum )**2
              pkernel(i,j,ig,jg,igrp,6) = pkernl * ( rmk/rmsum )**2
            enddo
          enddo
        endif
      enddo
    enddo
  enddo

  ! Do some extra debugging reports  (normally commented)
  if (do_print_init) then
    write(LUNOPRT,*) ' '
    write(LUNOPRT,*) 'Coagulation group mapping:'
    do ig = 1, NGROUP
      do jg = 1, NGROUP
        write(LUNOPRT,*) 'ig jg icoag = ', ig, jg, icoag(ig,jg)
      enddo
    enddo
    write(LUNOPRT,*) ' '
    write(LUNOPRT,*) 'Coagulation element mapping:'
    do ielem = 1, NELEM
      do ig = 1, NGROUP
        write(LUNOPRT,*) 'ielem ig icoagelem icomp(ielem) = ', &
          ielem, ig, icoagelem(ielem,ig), icomp(ielem)
      enddo
    enddo
    write(LUNOPRT,*) ' '
    write(LUNOPRT,*) 'Coagulation bin mapping arrays'
    do igrp = 1, NGROUP
      do ibin = 1,3
        write(LUNOPRT,*) 'igrp, ibin = ',igrp, ibin
        do ipair = 1,npairl(igrp,ibin)
          write(LUNOPRT,*) 'low:np,ig,jg,i,j ', &
              ipair,iglow(igrp,ibin,ipair), &
          jglow(igrp,ibin,ipair), ilow(igrp,ibin,ipair), &
                jlow(igrp,ibin,ipair)
        enddo
        do ipair = 1,npairu(igrp,ibin)
          write(LUNOPRT,*) 'up:np,ig,jg,i,j ', &
             ipair,igup(igrp,ibin,ipair), &
         jgup(igrp,ibin,ipair), iup(igrp,ibin,ipair), &
              jup(igrp,ibin,ipair)
        enddo
      enddo
    enddo
  endif
  
  ! Return to caller with coagulation mapping arrays defined
  return
end
