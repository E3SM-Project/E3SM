module mkglcmecMod
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mkglcmecMod
!
! !DESCRIPTION:
! Make glacier multi-elevation class  data
!
! !REVISION HISTORY:
! Author: Erik Kluzek, Mariana Vertenstein
!
!-----------------------------------------------------------------------
!!USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod , only : shr_sys_flush
  use mkdomainMod , only : domain_checksame
  implicit none

  private           ! By default make data private
!
! !PUBLIC MEMBER FUNCTIONS:
!
  public mkglcmecInit  ! Initialization
  public mkglcmec      ! Set glacier multi-elevation class
  public mkglacier     ! Set percent glacier
!
! !PUBLIC DATA MEMBERS: 
!
  integer, public       :: nglcec         = 10   ! number of elevation classes for glaciers
  real(r8), pointer     :: elevclass(:)          ! elevation classes
!EOP
!===============================================================
contains
!===============================================================

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkglcmecInit
!
! !INTERFACE:
subroutine mkglcmecInit( elevclass_o )
!
! !DESCRIPTION:
! Initialize of Make glacier multi-elevation class data
! !USES:
!
! !ARGUMENTS:
  implicit none
  real(r8), intent(OUT) :: elevclass_o(:)          ! elevation classes
!
! !CALLED FROM:
! subroutine mksrfdat in module mksrfdatMod
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
!
! !LOCAL VARIABLES:
!EOP
  character(len=32) :: subname = 'mkglcmecInit:: '
!-----------------------------------------------------------------------
  allocate( elevclass(nglcec+1) )

  ! -----------------------------------------------------------------
  ! Define elevation classes, represents lower boundary of each class
  ! -----------------------------------------------------------------

  if (      nglcec == 36 )then
     elevclass(:) = (/ 0.,   200.,   400.,   600.,   800.,  &
                    1000.,  1200.,  1400.,  1600.,  1800.,  &
                    2000.,  2200.,  2400.,  2600.,  2800.,  &
                    3000.,  3200.,  3400.,  3600.,  3800.,  &
                    4000.,  4200.,  4400.,  4600.,  4800.,  &
                    5000.,  5200.,  5400.,  5600.,  5800.,  &
                    6000.,  6200.,  6400.,  6600.,  6800.,  &
                    7000., 10000./)
  else if ( nglcec == 10 )then
     elevclass(1)  =     0.
     elevclass(2)  =   200.
     elevclass(3)  =   400.
     elevclass(4)  =   700.
     elevclass(5)  =  1000.
     elevclass(6)  =  1300.
     elevclass(7)  =  1600.
     elevclass(8)  =  2000.
     elevclass(9)  =  2500.
     elevclass(10) =  3000.
     elevclass(11) = 10000.
  else if ( nglcec == 5  )then
     elevclass(1)  =     0.
     elevclass(2)  =   500.
     elevclass(3)  =  1000.
     elevclass(4)  =  1500.
     elevclass(5)  =  2000.
     elevclass(6)  = 10000.
  else if ( nglcec == 3  )then
     elevclass(1)  =     0.
     elevclass(2)  =  1000.
     elevclass(3)  =  2000.
     elevclass(4)  = 10000.
  else if ( nglcec == 1  )then
     elevclass(1)  =     0.
     elevclass(2)  = 10000.
  else if ( nglcec == 0  )then
     elevclass(1)  = 10000.
  else
     write(6,*) subname//"ERROR:: nglcec must be 0, 1, 3, 5, 10 or 36",&
          " to work with CLM: "
     call abort()
  end if

  elevclass_o(:) = elevclass(:)

end subroutine mkglcmecInit

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkglcmec
!
! !INTERFACE:
subroutine mkglcmec(ldomain, mapfname, &
                    datfname_fglacier, ndiag, &
                    pctglac_o, pctglac_o_uncorrected, &
                    pctglcmec_o, topoglcmec_o, &
                    pctglcmec_gic_o, pctglcmec_icesheet_o, &
                    pctglc_gic_o, pctglc_icesheet_o)
!
! !DESCRIPTION: 
! make percent glacier on multiple elevation classes, mean elevation for each
! elevation class, and associated fields
! 
! Note that the raw glacier data are specified by level, and thus implicitly include the
! necessary topo data for breaking pct glacier into elevation classes. Each level in the
! input data is assigned to an elevation (given by BIN_CENTERS in the input data). Thus,
! all of the input glacier in level 1 is treated as being at the same elevation, and
! likewise for each other level. These elevations are then used in assigning pct_glacier
! to the appropriate elevation class in the output data, as well as determining the mean
! topographic height of each elevation class in the output data.
!
! Does nothing if nglcec==0.
!
! !USES:
  use mkdomainMod, only : domain_type, domain_clean, domain_read
  use mkgridmapMod
  use mkvarpar	
  use mkutilsMod, only : slightly_below, slightly_above
  use mkncdio
!
! !ARGUMENTS:
  implicit none
  type(domain_type) , intent(in) :: ldomain
  character(len=*)  , intent(in) :: mapfname           ! input mapping file name
  character(len=*)  , intent(in) :: datfname_fglacier  ! raw glacier data
  integer           , intent(in) :: ndiag              ! unit number for diag out
  real(r8)          , intent(in) :: pctglac_o(:)       ! % glac on output glacier grid
  real(r8)          , intent(in) :: pctglac_o_uncorrected(:) ! % glac on output glacier grid, before any corrections were done
  real(r8)          , intent(out):: pctglcmec_o (:,:)  ! % for each elevation class on output glacier grid
  real(r8)          , intent(out):: topoglcmec_o(:,:)  ! mean elevation for each elevation classs on output glacier grid
  real(r8)          , intent(out):: pctglcmec_gic_o(:,:)      ! % glc gic on output grid, by elevation class
  real(r8)          , intent(out):: pctglcmec_icesheet_o(:,:) ! % glc ice sheet on output grid, by elevation class
  real(r8)          , intent(out):: pctglc_gic_o(:)           ! % glc gic on output grid, summed across elevation classes
  real(r8)          , intent(out):: pctglc_icesheet_o(:)      ! % glc ice sheet on output grid, summed across elevation classes
!
! !CALLED FROM:
! subroutine mksrfdat in module mksrfdatMod
!
! !REVISION HISTORY:
! Author: David Lawrence
! 7/12/11: Bill Sacks: substantial rewrite to use input topo and % glacier at same resolution
! 9/25/12: Bill Sacks: substantial rewrite to use new format of fglacier, which provides
!          percent by elevation bin (thus the separate topo dataset is no longer needed
!          in this routine)
!
!
! !LOCAL VARIABLES:
!EOP
  type(domain_type)     :: tdomain              ! local domain
  type(gridmap_type)    :: tgridmap             ! local gridmap
  real(r8), allocatable :: pctglc_gic_i(:)      ! input GIC percentage for a single level
  real(r8), allocatable :: pctglc_icesheet_i(:) ! input icesheet percentage for a single level
  real(r8), allocatable :: topoglcmec_unnorm_o(:,:) ! same as topoglcmec_o, but unnormalized
  real(r8) :: topoice_i                         ! topographic height of this level
  real(r8) :: pctglc_i                          ! input total pct glacier for a single level & single point
  real(r8) :: wt, frac                          ! weighting factors for remapping
  integer  :: ndims                             ! number of dimensions in input variables
  integer  :: dim_lengths(nf_max_var_dims)      ! lengths of dimensions in input variables
  integer, allocatable :: starts(:), counts(:)  ! start indices & counts for reading variable slices
  integer  :: ni,no,ns_o,nst,lev                ! indices
  integer  :: n,m                               ! indices
  integer  :: ncid,dimid,varid                  ! input netCDF id's
  integer  :: nlev                              ! number of levels in input file
  real(r8) :: glc_sum                           ! temporary
  integer  :: ier                               ! error status
  logical  :: errors                            ! error status
  real(r8), parameter :: minglac = 1.e-6_r8     ! Minimum glacier amount
  character(len=32) :: subname = 'mkglcmec'
!-----------------------------------------------------------------------

  ! Initialize all output fields to zero

  pctglcmec_o(:,:)          = 0.
  topoglcmec_o(:,:)         = 0.
  pctglcmec_gic_o(:,:)      = 0.
  pctglcmec_icesheet_o(:,:) = 0.
  pctglc_gic_o(:)           = 0.
  pctglc_icesheet_o(:)      = 0.

  ! Set number of output points

  ns_o = ldomain%ns

  ! -----------------------------------------------------------------
  ! Exit early, if no glaciers exist
  ! -----------------------------------------------------------------
  if ( all(pctglac_o < minglac ) )then
     write (6,*) 'No glaciers exist, set glcmec to zero as well'
     return
  end if

  if ( nglcec == 0 )then
     write (6,*) 'Number of glacier elevation classes is zero ',&
          '-- set glcmec to zero as well'
     call shr_sys_flush(6)
     return
  end if

  write (6,*) 'Attempting to make percent elevation class ',&
       'and mean elevation for glaciers .....'
  call shr_sys_flush(6)

  ! -----------------------------------------------------------------
  ! Read domain and dimension information from glacier raw data file
  ! -----------------------------------------------------------------

  call domain_read(tdomain,datfname_fglacier)
  nst = tdomain%ns

  ! Read z dimension size
  write (6,*) 'Open glacier file: ', trim(datfname_fglacier)
  call check_ret(nf_open(datfname_fglacier, 0, ncid), subname)
  ier = nf_inq_dimid (ncid, 'z', dimid)
  if (ier /= NF_NOERR) then
     write (6,*) trim(subname), ' ERROR: z dimension not found on glacier file:'
     write (6,*) trim(datfname_fglacier)
     write (6,*) 'Perhaps you are trying to use an old-format glacier file?'
     write (6,*) '(prior to Sept., 2012)'
     call abort()
  end if
  call check_ret(nf_inq_dimlen (ncid, dimid, nlev), subname)

  ! -----------------------------------------------------------------
  ! Read mapping data, check for consistency with domains
  ! -----------------------------------------------------------------

  ! Mapping for raw glacier -> model output grid
  call gridmap_mapread(tgridmap, mapfname )

  ! Error checks for domain and map consistencies
  call domain_checksame( tdomain, ldomain, tgridmap )

  ! -----------------------------------------------------------------
  ! Determine dimension lengths and create start & count arrays
  ! for later reading one level at a time
  ! -----------------------------------------------------------------

  call get_dim_lengths(ncid, 'PCT_GLC_GIC', ndims, dim_lengths)

  allocate(starts(ndims), counts(ndims), stat=ier)
  if (ier/=0) call abort()

  starts(1:ndims) = 1

  ! We assume that the last dimension is the level dimension
  counts(1:ndims-1) = dim_lengths(1:ndims-1)
  counts(ndims) = 1

  ! -------------------------------------------------------------------- 
  ! Compute fields on the output grid
  ! -------------------------------------------------------------------- 

  allocate(pctglc_gic_i(nst), pctglc_icesheet_i(nst), stat=ier)
  if (ier/=0) call abort()

  allocate(topoglcmec_unnorm_o(ns_o,nglcec), stat=ier)
  if (ier/=0) call abort()

  topoglcmec_unnorm_o(:,:) = 0.

  write(6,'(a)',advance='no') 'Level: '
  do lev = 1, nlev
     write(6,'(i4)',advance='no') lev

     ! Read this level's data
     ! We assume that the last dimension is the level dimension
     starts(ndims) = lev
     call check_ret(nf_inq_varid (ncid, 'BIN_CENTERS', varid), subname)
     call check_ret(nf_get_vara_double (ncid, varid, (/lev/), (/1/), topoice_i), subname)
     call check_ret(nf_inq_varid (ncid, 'PCT_GLC_GIC', varid), subname)
     call check_ret(nf_get_vara_double (ncid, varid, starts, counts, pctglc_gic_i), subname)
     call check_ret(nf_inq_varid (ncid, 'PCT_GLC_ICESHEET', varid), subname)
     call check_ret(nf_get_vara_double (ncid, varid, starts, counts, pctglc_icesheet_i), subname)
     
     ! Determine elevation class
     m = get_elevclass(topoice_i)
     if (m < 1 .or. m > nglcec) then 
        call abort()
     end if

     do n = 1,tgridmap%ns
        ni = tgridmap%src_indx(n)
        no = tgridmap%dst_indx(n)
        wt = tgridmap%wovr(n)

        ! fraction of this destination cell that is covered by source cells that are within the source landmask
        frac = tgridmap%frac_dst(no)

        ! We don't bother with the following if pctglac_o(no) <= minglac, because the
        ! pctglcmec_o values will end up being <= minglac.
        ! Also, if frac == 0, then we can't do this, to avoid divide by 0. In this case, the
        ! outputs remain equal to 0 (their initialized value).
        if (pctglac_o(no) > minglac .and. frac > 0) then
           pctglc_i = pctglc_gic_i(ni) + pctglc_icesheet_i(ni)
           pctglcmec_o(no,m)          = pctglcmec_o(no,m)          + wt*pctglc_i / frac
           pctglcmec_gic_o(no,m)      = pctglcmec_gic_o(no,m)      + wt*pctglc_gic_i(ni) / frac
           pctglcmec_icesheet_o(no,m) = pctglcmec_icesheet_o(no,m) + wt*pctglc_icesheet_i(ni) / frac

           ! note that, by weighting the following by pctglc_i, we are getting something
           ! like the average topographic height over glaciated areas - NOT the average
           ! topographic height of the entire grid cell
           topoglcmec_unnorm_o(no,m) = topoglcmec_unnorm_o(no,m) + wt*pctglc_i*topoice_i / frac
        end if
     end do
  end do

  ! advance to next line (needed because of 'advance=no' writes above)
  write(6,*) ' ' 

  ! Close glacier input file
  call check_ret(nf_close(ncid), subname)

  ! Now pctglcmec_o, pctglcmec_gic_o and pctglcmec_icesheet_o are already normalized
  ! appropriately, because for destarea normalization, the sum of the wovr values for a
  ! given destination grid cell should equal frac_dst, so sum(wt/frac) = 1 for each
  ! destination grid cell.
  !
  ! However, we need to normalize topoglcmec_o. To do this, note that pctglcmec_o(n,m) is
  ! equal to the sum of the weights used in doing the weighted average of topoice_i
  ! (weight = wt*pctglc_i/frac); hence pctglcmec_o(n,m) is the correct normalization
  ! factor
  do no = 1,ns_o
     do m = 1,nglcec
        if (pctglcmec_o(no,m) > 0) then
           topoglcmec_o(no,m) = topoglcmec_unnorm_o(no,m) / pctglcmec_o(no,m)
        end if

        ! Correct for rounding errors that put topoglcmec_o(no,m) slightly outside the
        ! allowed bounds for this elevation class
        if (slightly_below(topoglcmec_o(no,m), elevclass(m))) then
           write(6,*) 'Warning: topoglcmec_o was slightly lower than lower bound; setting equal&
                & to lower bound; for: ', no, m, topoglcmec_o(no,m), elevclass(m)
           write(6,*) '(this is informational only, and probably just indicates rounding error)'
           topoglcmec_o(no,m) = elevclass(m)
        else if (slightly_above(topoglcmec_o(no,m), elevclass(m+1))) then
           write(6,*) 'Warning: topoglcmec_o was slightly higher than upper bound; setting equal&
                & to upper bound; for: ', no, m, topoglcmec_o(no,m), elevclass(m+1)
           write(6,*) '(this is informational only, and probably just indicates rounding error)'
           topoglcmec_o(no,m) = elevclass(m+1)
        end if
     end do
  end do

  ! We want the sum of pctglcmec_o across elevation classes to equal pctglac_o for each
  ! grid cell. This should be approximately true at this point, but there have been a
  ! number of adjustments to pctglac_o since its original remapping (e.g., setting values
  ! < 1 to 0, and rescaling all percentages so that the total landcover is equal to
  ! 100%). We need to essentially apply these same adjustments to pctglcmec_o.
  ! 
  ! To do this, we use the stored, uncorrected values of pctglac_o. In theory, we could do
  ! this renormalization without reference to these uncorrected values, since
  ! pctglac_o_uncorrected should equal the sum of pctglcmec_o across elevation classes for
  ! each grid cell. However, pctglac_o_uncorrected is helpful for error checking: it
  ! allows us to verify that the pctglcmec sums are correct (this check is done below).
  ! 
  ! Here we also rescale pctglcmec_gic_o and pctglcmec_icesheet_o similarly, to
  ! achieve the same thing (i.e., we want sums to add up to pct_glacier)
  do no = 1,ns_o
     if (pctglac_o_uncorrected(no) > 0) then
        pctglcmec_o(no,:)          = pctglcmec_o(no,:)          * (pctglac_o(no) / pctglac_o_uncorrected(no))
        pctglcmec_gic_o(no,:)      = pctglcmec_gic_o(no,:)      * (pctglac_o(no) / pctglac_o_uncorrected(no))
        pctglcmec_icesheet_o(no,:) = pctglcmec_icesheet_o(no,:) * (pctglac_o(no) / pctglac_o_uncorrected(no))

     else if (pctglac_o_uncorrected(no) == 0 .and. pctglac_o(no) > 0) then
        ! There isn't a clear way to handle the case when pctglac_o_uncorrected==0 (and
        ! hence pctglcmec_o==0 for all elevation classes) but pctglac_o==0. Fortunately,
        ! this should never happen: the only place where there is a potential for this is
        ! over the south pole, where pctglac_o is set to 100%.

        write(6,*) 'ERROR in ', subname
        write(6,*) 'pctglac_o_uncorrected==0 but pctglac_o > 0'
        write(6,*) 'no = ', no
        write(6,*) 'lon = ', tgridmap%xc_dst(no)
        write(6,*) 'lat = ', tgridmap%yc_dst(no)
        call abort()
     end if

     ! The only other possibility is that (pctglac_o_uncorrected(no) == 0 .and.
     ! pctglac_o(no) == 0); but in this case, all pctglcmec_o(no,:) values should be 0
     ! and can remain 0 (and similarly for pctglc_gic, pctglc_icesheet, etc.)
  end do

  ! Set pctglc_gic_o to sum of pctglcmec_gic_o across elevation classes, and similarly for pctglc_icesheet_o
  pctglc_gic_o      = sum(pctglcmec_gic_o, dim=2)
  pctglc_icesheet_o = sum(pctglcmec_icesheet_o, dim=2)

  ! -------------------------------------------------------------------- 
  ! Perform various sanity checks
  ! -------------------------------------------------------------------- 

  errors = .false.

  ! Check that sum over pctglcmec_o (from 1 to nglcec) is equal to pctglac_o(no)
  !
  ! Note: We use a threshold of 2e-5 rather than 1e-6 or smaller because there is single-
  ! precision roundoff error in the storage of the flat (i.e., non-multi-level)
  ! PCT_GLACIER variable on the input dataset (which is used to compute pctglac_o) (since
  ! the input dataset uses single-precision floats)
  do no = 1,ns_o
     glc_sum = 0.
     do m = 1,nglcec
        glc_sum = glc_sum + pctglcmec_o(no,m)
     end do
     if (abs(glc_sum - pctglac_o(no)) > 2.e-5) then
        write(6,*)'no,pctglc,pctglac= ',no,glc_sum,pctglac_o(no)
        errors = .true.
     end if
  end do

  ! Check that GIC + ICESHEET = total glacier, for variables summed over elevation classes
  !
  ! Note: We use a threshold of 2e-5 rather than 1e-6 or smaller because there is single-
  ! precision roundoff error in the storage of the flat (i.e., non-multi-level)
  ! PCT_GLACIER variable on the input dataset (which is used to compute pctglac_o) (since
  ! the input dataset uses single-precision floats)
  do no = 1,ns_o
     if (abs((pctglc_gic_o(no) + pctglc_icesheet_o(no)) - pctglac_o(no)) > 2.e-5) then
        write(6,*)'no,pctglc_gic,pctglc_icesheet,pctglac,pctglac_uncorrected,lon,lat=',no,pctglc_gic_o(no),&
             pctglc_icesheet_o(no),pctglac_o(no),pctglac_o_uncorrected(no),tgridmap%&
             xc_dst(no),tgridmap%yc_dst(no)
        errors = .true.
     end if
  end do

  ! Check that GIC + ICESHEET = total glacier, for each elevation class
  do m = 1,nglcec
     do no = 1,ns_o
        if (abs((pctglcmec_gic_o(no,m) + pctglcmec_icesheet_o(no,m)) - pctglcmec_o(no,m)) > 1.e-6) then
           write(6,*)'no,m,pctglcmec_gic,pctglcmec_icesheet,pctglcmec,lat,lon=',no,m,&
                pctglcmec_gic_o(no,m),pctglcmec_icesheet_o(no,m),pctglcmec_o(no,m),&
                tgridmap%xc_dst(no),tgridmap%yc_dst(no)
           errors = .true.
        end if
     end do
  end do

  ! Error check: are all elevations within elevation class range
  do no = 1,ns_o
     if (pctglac_o(no) .gt. minglac) then 
        do m = 1,nglcec
           ! WJS (7-12-11): the original check was ((topo outside range) and (topo .ne.
           ! 0)). I think that the second condition is meant to essentially check whether
           ! pctglcmec_o(no,m) > 0. So I think we could just check whether ((topo outside
           ! range) and (pctglcmec_o(no,m) > 0). However, I am keeping a warning message
           ! if ((topo outside range) and (topo .ne. 0)) as well, because I don't want to
           ! get rid of error checks.
           if ( (topoglcmec_o(no,m) .lt. elevclass(m) .or. topoglcmec_o(no,m) .gt. elevclass(m+1)) &
                .and. (pctglcmec_o(no,m) .gt. 0 .or. topoglcmec_o(no,m) .ne. 0)) then
              write(6,*) 'Error: mean elevation does not fall within elevation class '
              write(6,*) elevclass(m),elevclass(m+1),topoglcmec_o(no,m),pctglcmec_o(no,m),m,no
              errors = .true.
           endif
        end do
     end if
  end do

  if (errors) then
     call abort()
  end if

  ! Deallocate dynamic memory

  call domain_clean(tdomain)
  call gridmap_clean(tgridmap)
  deallocate(pctglc_gic_i, pctglc_icesheet_i)
  deallocate(topoglcmec_unnorm_o)
  deallocate(starts, counts)

  write (6,*) 'Successfully made percent elevation class and mean elevation for glaciers'
  write (6,*)
  call shr_sys_flush(6)

end subroutine mkglcmec

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkglacier
!
! !INTERFACE:
subroutine mkglacier(ldomain, mapfname, datfname, ndiag, zero_out, glac_o, glac_uncorrected)
!
! !DESCRIPTION:
! make percent glacier
!
! In contrast to mkglcmec, this uses a "flat" PCT_GLACIER field (not separated by
! elevation class, and not separated into icesheet vs GIC). 
!
! This simpler routine is sufficient for cases when we run without multiple elevation
! classes. This routine is also used when running with multiple elevation classes: we
! first regrid the flat PCT_GLACIER field, then later create the multiple elevation class
! data. This multi-step process makes it easier to do corrections on the total
! PCT_GLACIER, and make sure these corrections apply appropriately to the multi-level
! output. The assumption is that PCT_GLACIER is the sum of both PCT_GLC_GIC and
! PCT_GLC_ICESHEET across all elevation bins.
!
! !USES:
  use mkdomainMod , only : domain_type, domain_clean, domain_read
  use mkgridmapMod
  use mkvarpar	
  use mkvarctl    
  use mkncdio
!
! !ARGUMENTS:
  implicit none
  type(domain_type), intent(in) :: ldomain
  character(len=*)  , intent(in) :: mapfname  ! input mapping file name
  character(len=*)  , intent(in) :: datfname  ! input data file name
  integer           , intent(in) :: ndiag     ! unit number for diag out
  logical           , intent(in) :: zero_out  ! if should zero glacier out
  real(r8)          , intent(out):: glac_o(:) ! output grid: %glacier
  real(r8)          , intent(out):: glac_uncorrected(:) ! output grid: %glacier before any
                                                        ! corrections are done
!
! !CALLED FROM:
! subroutine mksrfdat in module mksrfdatMod
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
  type(gridmap_type)   :: tgridmap
  type(domain_type)    :: tdomain            ! local domain
  real(r8), allocatable :: glac_i(:)          ! input grid: percent glac
  real(r8) :: sum_fldi                        ! global sum of dummy input fld
  real(r8) :: sum_fldo                        ! global sum of dummy output fld
  real(r8) :: gglac_i                         ! input  grid: global glac
  real(r8) :: garea_i                         ! input  grid: global area
  real(r8) :: gglac_o                         ! output grid: global glac
  real(r8) :: garea_o                         ! output grid: global area
  integer  :: ni,no,k,n,m,ns                  ! indices
  integer  :: ncid,dimid,varid                ! input netCDF id's
  integer  :: ier                             ! error status
  real(r8) :: relerr = 0.00001                ! max error: sum overlap wts ne 1
  character(len=32) :: subname = 'mkglacier'
!-----------------------------------------------------------------------

  write (6,*) 'Attempting to make %glacier .....'
  call shr_sys_flush(6)

  ! -----------------------------------------------------------------
  ! Read input file
  ! -----------------------------------------------------------------

  ! Obtain input grid info, read local fields

  call domain_read(tdomain,datfname)
  ns = tdomain%ns
  allocate(glac_i(ns), stat=ier)
  if (ier/=0) call abort()

  write (6,*) 'Open glacier file: ', trim(datfname)
  call check_ret(nf_open(datfname, 0, ncid), subname)
  call check_ret(nf_inq_varid (ncid, 'PCT_GLACIER', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, glac_i), subname)
  call check_ret(nf_close(ncid), subname)

  ! Area-average percent cover on input grid to output grid 
  ! and correct according to land landmask
  ! Note that percent cover is in terms of total grid area.

  if ( zero_out )then

     do no = 1, ldomain%ns
        glac_o(no) = 0.
     enddo

  else

     call gridmap_mapread(tgridmap, mapfname )

     ! Error checks for domain and map consistencies
     call domain_checksame( tdomain, ldomain, tgridmap )
     
     ! Determine glac_o on output grid

     call gridmap_areaave(tgridmap, glac_i, glac_o)
     
     ! Save a copy of glac_o before any corrections are done. This is needed for
     ! normalization in mkglcmec
     glac_uncorrected(:) = glac_o(:)

     do no = 1, ldomain%ns
        if (glac_o(no) < 1.) glac_o(no) = 0.
     enddo
  end if

  ! Check for conservation

  do no = 1, ldomain%ns
     if ((glac_o(no)) > 100.000001_r8) then
        write (6,*) 'MKGLACIER error: glacier = ',glac_o(no), &
                ' greater than 100.000001 for column, row = ',no
        call shr_sys_flush(6)
        stop
     end if
  enddo

  ! Some error checking and writing of global values before and after the regrid

  if ( .not. zero_out )then

     ! Global sum of output field -- must multiply by fraction of
     ! output grid that is land as determined by input grid

     sum_fldi = 0.0_r8
     do ni = 1, tdomain%ns
        sum_fldi = sum_fldi + tgridmap%area_src(ni) * tgridmap%frac_src(ni)
     enddo

     sum_fldo = 0.
     do no = 1, ldomain%ns
        sum_fldo = sum_fldo + tgridmap%area_dst(no) * tgridmap%frac_dst(no)
     end do

     ! -----------------------------------------------------------------
     ! Error check1
     ! Compare global sum fld_o to global sum fld_i.
     ! -----------------------------------------------------------------

     if ( trim(mksrf_gridtype) == 'global') then
        if ( abs(sum_fldo/sum_fldi-1.) > relerr ) then
           write (6,*) 'MKGLACIER error: input field not conserved'
           write (6,'(a30,e20.10)') 'global sum output field = ',sum_fldo
           write (6,'(a30,e20.10)') 'global sum input  field = ',sum_fldi
           stop
        end if
     end if

     ! -----------------------------------------------------------------
     ! Error check2
     ! Compare global areas on input and output grids
     ! -----------------------------------------------------------------

     ! Input grid

     gglac_i = 0.
     garea_i = 0.
     do ni = 1, tdomain%ns
        garea_i = garea_i + tgridmap%area_src(ni)*re**2
        gglac_i = gglac_i + glac_i(ni)*(tgridmap%area_src(ni)/100.)*&
                                        tgridmap%frac_src(ni)*re**2
     end do

     ! Output grid

     gglac_o = 0.
     garea_o = 0.
     do no = 1, ldomain%ns
        garea_o = garea_o + tgridmap%area_dst(no)*re**2
        gglac_o = gglac_o + glac_o(no)*(tgridmap%area_dst(no)/100.)*&
                                        tgridmap%frac_dst(no)*re**2
     end do

     ! Diagnostic output

     write (ndiag,*)
     write (ndiag,'(1x,70a1)') ('=',k=1,70)
     write (ndiag,*) 'Glacier Output'
     write (ndiag,'(1x,70a1)') ('=',k=1,70)

     write (ndiag,*)
     write (ndiag,'(1x,70a1)') ('.',k=1,70)
     write (ndiag,2001)
2001 format (1x,'surface type   input grid area  output grid area'/ &
             1x,'                 10**6 km**2      10**6 km**2   ')
     write (ndiag,'(1x,70a1)') ('.',k=1,70)
     write (ndiag,*)
     write (ndiag,2002) gglac_i*1.e-06,gglac_o*1.e-06
     write (ndiag,2004) garea_i*1.e-06,garea_o*1.e-06
2002 format (1x,'glaciers    ',f14.3,f17.3)
2004 format (1x,'all surface ',f14.3,f17.3)

  end if

  ! Deallocate dynamic memory

  call domain_clean(tdomain) 
  if ( .not. zero_out )then
     call gridmap_clean(tgridmap)
     deallocate (glac_i)
  end if

  write (6,*) 'Successfully made %glacier'
  write (6,*)
  call shr_sys_flush(6)

end subroutine mkglacier

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_elevclass
!
! !INTERFACE:
integer function get_elevclass(topo, writewarn)
!
! !DESCRIPTION:
! Returns elevation class index (1..nglcec) given the topographic height.
! If topo is lower than the lowest elevation class, returns 0.
! If topo is higher than the highest elevation class, returns (nglcec+1).
! In either of the two latter cases, the function also writes a warning message, unless
! writewarn is present and false.
!
! !ARGUMENTS:
   implicit none
   real(r8), intent(in) :: topo  ! topographic height (m)
   logical, intent(in), optional :: writewarn  ! should warning messages be written? (default: true)
!
! !REVISION HISTORY:
! Author: Bill Sacks
!
! !LOCAL VARIABLES:
!EOP
   integer :: m
   logical :: my_writewarn
   character(len=32) :: subname = 'mkglcmec'
!-----------------------------------------------------------------------
   
   if (present(writewarn)) then
      my_writewarn = writewarn
   else
      my_writewarn = .true.
   end if

   if (topo < elevclass(1)) then
      if (my_writewarn) then
         write(6,*) 'WARNING in ', trim(subname)
         write(6,*) 'topo out of bounds'
         write(6,*) 'topo = ', topo
         write(6,*) 'elevclass(1) = ', elevclass(1)
      end if
      get_elevclass = 0
      return
   end if
   
   do m = 1, nglcec
      if (topo < elevclass(m+1)) then
         ! note that we already know that topo >= elevclass(m), otherwise we would have
         ! returned earlier
         get_elevclass = m
         return
      end if
   end do

   if (my_writewarn) then
      write(6,*) 'WARNING in ', trim(subname)
      write(6,*) 'topo out of bounds'
      write(6,*) 'topo = ', topo
      write(6,*) 'elevclass(nglcec+1) = ', elevclass(nglcec+1)
   end if
   get_elevclass = nglcec+1

end function get_elevclass
         
!-----------------------------------------------------------------------

end module mkglcmecMod
