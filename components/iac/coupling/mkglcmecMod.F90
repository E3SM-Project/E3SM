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

  if (      nglcec == 10 )then
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
  else
     write(6,*) subname//"ERROR:: nglcec must be 1, 3, 5, or 10 to work with CLM: "
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
subroutine mkglcmec(ldomain, mapfname_t2g, mapfname_g2g, datfname_fglctopo, datfname_fglacier, &
                    ndiag, pctglac_o, pctglcmec_o, topoglcmec_o, thckglcmec_o )
!
! !DESCRIPTION:
! make percent glacier on multiple elevation classes and mean elevation for each elevation class
!
! !USES:
  use mkfileutils, only : getfil
  use mkdomainMod, only : domain1_type, domain1_clean, domain1_read
  use mkgridmapMod
  use mkvarpar	
  use mkncdio
!
! !ARGUMENTS:
  implicit none
  type(domain1_type), intent(in) :: ldomain
  character(len=*)  , intent(in) :: mapfname_t2g       ! map raw topo -> raw glacier
  character(len=*)  , intent(in) :: mapfname_g2g       ! map raw glacier -> output glacier
  character(len=*)  , intent(in) :: datfname_fglctopo  ! raw glc topo data
  character(len=*)  , intent(in) :: datfname_fglacier  ! raw glacier data
  integer           , intent(in) :: ndiag              ! unit number for diag out
  real(r8)          , intent(in) :: pctglac_o(:)       ! % glac on output glacier grid
  real(r8)          , intent(out):: pctglcmec_o (:,:)  ! % for each elevation class on output glacier
  real(r8)          , intent(out):: topoglcmec_o(:,:)  ! mean elevation for each elevation classs on output glacier 
  real(r8)          , intent(out):: thckglcmec_o(:,:)  ! mean ice thickness for each elevation classs on input glacier 
!
! !CALLED FROM:
! subroutine mksrfdat in module mksrfdatMod
!
! !REVISION HISTORY:
! Author: David Lawrence
!
!
! !LOCAL VARIABLES:
!EOP
  type(domain1_type)    :: tdomain_topo        ! local domain: topo_ice , topo_bedrock
  type(domain1_type)    :: tdomain_glac        ! local domain: fracdata
  type(gridmap_type)    :: tgridmap            ! local gridmap
  real(r8), allocatable :: pctglcmec_i(:,:)    ! % for each elevation class on input glacier grid
  real(r8), allocatable :: topoglcmec_i(:,:)   ! mean elevation for each elevation classs on input glacier grid
  real(r8), allocatable :: thckglcmec_i(:,:)   ! mean ice thickness for each elevation classs on input glacier grid
  real(r8), allocatable :: topoice_i(:)        ! topo of ice surface
  real(r8), allocatable :: topobedrock_i(:)    ! topo of bedrock surface
  real(r8), allocatable :: glac_i(:)           ! input glacier fraction
  integer , allocatable :: nptsec(:,:)         ! number of points in an elevation class 
  real(r8), allocatable :: topoec(:,:)         ! sum of all elevations, weighted by area in elevation class
  real(r8), allocatable :: thckec(:,:)         ! sum of all elevations, weighted by area in elevation class
  real(r8), allocatable :: areaec(:,:)         ! total area for all points within elevation class
  real(r8), allocatable :: pcttotec(:,:)       ! accumulator for total percentage area in grid 
  real(r8), allocatable :: areatot(:)          ! total area for glacier grid cell
  real(r8), allocatable :: pctareaec(:)        ! accumulator for pct area for elevation class
  real(r8), allocatable :: pctectot(:)         ! accumulator for total percentage area of elevation classes
  real(r8), allocatable :: pcttot(:)           ! accumulator for total percentage area in grid 
  integer , allocatable :: nptsectot(:)        ! total number of points within glacier grid cell that are not ocean
  integer , allocatable :: po(:)               ! temporary flag
  real(r8) :: wt                               ! temporary weight
  integer  :: ni,no,ns_o,nst_i,nsg_i           ! indices
  integer  :: k,l,n,m                          ! indices
  integer  :: ncid,dimid,varid                 ! input netCDF id's
  integer  :: ier                              ! error status
  real(r8), parameter :: minglac = 1.e-6_r8    ! Minimum glacier amount
  character(len=256) locfn                     ! local dataset file name
  character(len=32) :: subname = 'mkglcmec'
!-----------------------------------------------------------------------

  ! Initialize output to zero

  pctglcmec_o = 0.
  topoglcmec_o = 0.
  thckglcmec_o = 0.

  ns_o = ldomain%ns

  ! -----------------------------------------------------------------
  ! Exit early, if no glaciers exist
  ! -----------------------------------------------------------------
  if ( all(pctglac_o < minglac ) )then
     write (6,*) 'No glaciers exist, set glcmec to zero as well'
     return
  end if

  write (6,*) 'Attempting to make percent elevation class and mean elevation for glaciers .....'
  call shr_sys_flush(6)

  ! -----------------------------------------------------------------
  ! Read input file
  ! -----------------------------------------------------------------

  ! Get raw topo data 

  call getfil (datfname_fglctopo, locfn, 0)
  call domain1_read(tdomain_topo,locfn)
  nst_i = tdomain_topo%ns

  allocate(topoice_i(nst_i), topobedrock_i(nst_i), stat=ier)
  if (ier/=0) call abort()

  call check_ret(nf_open(locfn, 0, ncid), subname)
  call check_ret(nf_inq_varid (ncid, 'TOPO_ICE', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, topoice_i), subname)
  call check_ret(nf_inq_varid (ncid, 'TOPO_BEDROCK', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, topobedrock_i), subname)
  call check_ret(nf_close(ncid), subname)

  ! Get raw glacier data 

  call getfil (datfname_fglacier, locfn, 0)
  call domain1_read(tdomain_glac,locfn)
  nsg_i = tdomain_glac%ns
  allocate(glac_i(nsg_i), stat=ier)
  if (ier/=0) call abort()

  call check_ret(nf_open(locfn, 0, ncid), subname)
  call check_ret(nf_inq_varid (ncid, 'PCT_GLACIER', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, glac_i), subname)
  call check_ret(nf_close(ncid), subname)

  ! Mapping for raw topo -> raw glacier data (e.g. 5 min -> 1/2 degree)

  call gridmap_mapread(tgridmap, mapfname_t2g )

  ! Error checks for domain and map consistencies
  ! Note that the topo dataset has no landmask - so a unit landmask is assumed

  if (tdomain_topo%ns /= tgridmap%na) then
     write(6,*)'input domain size and gridmap source size are not the same size'
     write(6,*)' domain size = ',tdomain_topo%ns
     write(6,*)' map src size= ',tgridmap%na
     stop
  end if
  do n = 1,tgridmap%ns
     ni = tgridmap%src_indx(n)
     if (tgridmap%mask_src(ni) /= 1) then
        write(6,*)'input domain mask and gridmap mask are not the same at ni = ',ni
        write(6,*)' domain  mask= ',tdomain_topo%mask(ni)
        write(6,*)' gridmap mask= 1'
        stop
     end if
     if (tdomain_topo%lonc(ni) /= tgridmap%xc_src(ni)) then
        write(6,*)'input domain lon and gridmap lon not the same at ni = ',ni
        write(6,*)' domain  lon= ',tdomain_topo%lonc(ni)
        write(6,*)' gridmap lon= ',tgridmap%xc_src(ni)
        stop
     end if
     if (tdomain_topo%latc(ni) /= tgridmap%yc_src(ni)) then
        write(6,*)'input domain lat and gridmap lat not the same at ni = ',ni
        write(6,*)' domain  lat= ',tdomain_topo%latc(ni)
        write(6,*)' gridmap lat= ',tgridmap%yc_src(ni)
        stop
     end if
  end do

  ! First calculate elevation classes on input glacier grid
  ! Output is topoec,  
  ! Note that glac_i is on raw glacier grid, so will have index no

  allocate(pctglcmec_i(nsg_i,nglcec),  &
           topoglcmec_i(nsg_i,nglcec), &
           thckglcmec_i(nsg_i,nglcec), &
           nptsec(nsg_i,nglcec), &
           topoec(nsg_i,nglcec), &    
           thckec(nsg_i,nglcec), &
           areaec(nsg_i,nglcec), &
           areatot(nsg_i), &   
           pctareaec(nsg_i), & 
           pctectot(nsg_i), &  
           nptsectot(nsg_i), &
           po(nsg_i), stat=ier) 
  if (ier/=0) call abort()

  pctglcmec_i(:,:)  = 0.
  topoglcmec_i(:,:) = 0.
  thckglcmec_i(:,:) = 0.

  nptsec(:,:)  = 0
  topoec(:,:)  = 0.
  thckec(:,:)  = 0.
  areaec(:,:)  = 0.
  areatot(:)   = 0.
  pctareaec(:) = 0.
  pctectot(:)  = 0.
  po(:)        = 0
  nptsectot(:) = 0
     
  do n = 1,tgridmap%ns
     ni = tgridmap%src_indx(n)
     no = tgridmap%dst_indx(n)
     wt = tgridmap%wovr(n)

     if (glac_i(no) > minglac) then 
        areatot(no) = areatot(no) + tgridmap%area_src(ni)
        if (tgridmap%frac_src(ni) > 0) then
           do m = 1,nglcec
              if (topoice_i(ni) .ge. elevclass(m)  .and. &
                  topoice_i(ni) .lt. elevclass(m+1)) then
                 nptsec(no,m) = nptsec(no,m)  + 1
                 nptsectot(no)= nptsectot(no) + nptsec(no,m)
                 topoec(no,m) = topoec(no,m)  + wt*tgridmap%area_src(ni)*topoice_i(ni)
                 
                 ! bedrock cannot be below mean sea level; required to avoid overly thick ice sheets 
                 ! in ice shelf terrain (?)
                 if ( (topoice_i(ni) - topobedrock_i(ni)) .gt. elevclass(m+1) ) then
                    topobedrock_i(ni) = 0
                 endif
                 thckec(no,m) = thckec(no,m) + wt*tgridmap%area_src(ni)*(topoice_i(ni)-topobedrock_i(ni))
                 areaec(no,m) = areaec(no,m) + wt*tgridmap%area_src(ni)
              endif
           end do
        end if
     end if
  end do

  do no = 1,nsg_i
     if (glac_i(no) > minglac .and. areatot(no) > 0.) then 
        do m = nglcec,1,-1
           pctareaec(no) = pctareaec(no) + 100.*areaec(no,m)/areatot(no)
           if (pctareaec(no) .le. glac_i(no) .and. areaec(no,m) .gt. 0) then 
              pctglcmec_i(no,m)  = areaec(no,m)/areatot(no)*100.
              topoglcmec_i(no,m) = topoec(no,m)/areaec(no,m)
              thckglcmec_i(no,m) = thckec(no,m)/areaec(no,m)
           else if (pctareaec(no) .gt. glac_i(no) .and. po(no) .eq. 0) then 
              pctglcmec_i(no,m)  = areaec(no,m)/areatot(no)*100.
              topoglcmec_i(no,m) = topoec(no,m)/areaec(no,m)
              thckglcmec_i(no,m) = thckec(no,m)/areaec(no,m)
              po(no) = 1
           endif
           
           ! if all topo points within a glacier grid point are zero, then glacier is ice 
           ! shelf with an assumed elevation of 5m
           if (nptsectot(no) .eq. 0) then
              pctglcmec_i(no,1) = 100.
              topoglcmec_i(no,1) = 5
              thckglcmec_i(no,1) = 5
           endif
           pctectot(no) = pctectot(no) + pctglcmec_i(no,m)
        enddo

        ! error check: are all elevations within elevation class range
        do m = 1,nglcec
           if ((topoglcmec_i(no,m) .lt. elevclass(m)  .or. &
                topoglcmec_i(no,m) .gt. elevclass(m+1)) .and. &
                topoglcmec_i(no,m) .ne. 0) then
              write(6,*) 'Warning: mean elevation does not fall within elevation class '
              write(6,*) elevclass(m),elevclass(m+1),topoglcmec_i(no,m),m,no
           endif
        enddo
        
        ! normalize so that sum of pctglcmec_i adds up to one
        if ( pctectot(no) /= 0.0_r8 ) then
           do m = 1,nglcec
              pctglcmec_i(no,m) = (pctglcmec_i(no,m)/pctectot(no)) * 100.
           end do
        end if

     end if
  end do

  call gridmap_clean(tgridmap)

  deallocate(nptsec,    &
             topoec,    &
             thckec,    &
             areaec,    &
             areatot,   &
             pctareaec, &
             pctectot,  &
             nptsectot, &
             po)

  ! Average from input pct_glacier to output grid
  ! Note that glac_i(no) above is now glaci_i(ni) here

  call gridmap_mapread(tgridmap, mapfname_g2g )

  ! Error checks for domain and map consistencies

  if (tdomain_glac%ns /= tgridmap%na) then
     write(6,*)'input domain size and gridmap source size are not the same size'
     write(6,*)' domain size = ',tdomain_glac%ns
     write(6,*)' map src size= ',tgridmap%na
     stop
  end if
  do n = 1,tgridmap%ns
     ni = tgridmap%src_indx(n)
     if (tdomain_glac%mask(ni) /= tgridmap%mask_src(ni)) then
        write(6,*)'input domain mask and gridmap mask are not the same at ni = ',ni
        write(6,*)' domain  mask= ',tdomain_glac%mask(ni)
        write(6,*)' gridmap mask= ',tgridmap%mask_src(ni)
        stop
     end if
     if (tdomain_glac%lonc(ni) /= tgridmap%xc_src(ni)) then
        write(6,*)'input domain lon and gridmap lon not the same at ni = ',ni
        write(6,*)' domain  lon= ',tdomain_glac%lonc(ni)
        write(6,*)' gridmap lon= ',tgridmap%xc_src(ni)
        stop
     end if
     if (tdomain_glac%latc(ni) /= tgridmap%yc_src(ni)) then
        write(6,*)'input domain lat and gridmap lat not the same at ni = ',ni
        write(6,*)' domain  lat= ',tdomain_glac%latc(ni)
        write(6,*)' gridmap lat= ',tgridmap%yc_src(ni)
        stop
     end if
  end do

  allocate(pcttot(ns_o), pcttotec(ns_o,nglcec))

  pcttot(:)   = 0.
  pcttotec(:,:) = 0.

  do n = 1,tgridmap%ns
     ni = tgridmap%src_indx(n)
     no = tgridmap%dst_indx(n)
     wt = tgridmap%wovr(n)

     if (pctglac_o(no) .gt. minglac) then 
        pcttot(no) = pcttot(no) + wt*glac_i(ni)
        do m = 1,nglcec
           pcttotec(no,m) = pcttotec(no,m) + wt*glac_i(ni)*pctglcmec_i(ni,m)
        enddo
     end if
  end do

  do n = 1,tgridmap%ns
     ni = tgridmap%src_indx(n)
     no = tgridmap%dst_indx(n)
     wt = tgridmap%wovr(n)

     if (pctglac_o(no) .gt. minglac) then 
        do m = 1,nglcec
           pctglcmec_o(no,m) = pctglcmec_o(no,m) + wt*glac_i(ni)*pctglcmec_i(ni,m)/pcttot(no)
           if (pcttotec(no,m) .ne. 0) then
              topoglcmec_o(no,m) = topoglcmec_o(no,m) + wt*pctglcmec_i(ni,m)*glac_i(ni)*&
                                                        topoglcmec_i(ni,m)/pcttotec(no,m)
              thckglcmec_o(no,m) = thckglcmec_o(no,m) + wt*pctglcmec_i(ni,m)*glac_i(ni)*&
                                                        thckglcmec_i(ni,m)/pcttotec(no,m)
           endif
        end do
     end if
  end do
     
  do no = 1,ns_o
     if (pctglac_o(no) .gt. minglac) then 
        ! Scale according to grid cell pct_glacier
        do m = 1,nglcec
           pctglcmec_o(no,m) = pctglac_o(no)*pctglcmec_o(no,m)/100._r8
        end do
        
        ! Error check: are all elevations within elevation class range
        do m = 1,nglcec
           if ( (topoglcmec_o(no,m) .lt. elevclass(m) .or. topoglcmec_o(no,m) .gt. elevclass(m+1)) &
                .and. topoglcmec_o(no,m) .ne. 0) then
              write(6,*) 'Warning: mean elevation does not fall within elevation class '
              write(6,*) elevclass(m),elevclass(m+1),topoglcmec_o(no,m),m,no
           endif
        end do
     end if
  end do
     
  ! Deallocate dynamic memory

  call domain1_clean(tdomain_topo)
  call domain1_clean(tdomain_glac)
  call gridmap_clean(tgridmap)

  deallocate (topoice_i, topobedrock_i)
  deallocate (glac_i)
  deallocate (pctglcmec_i, topoglcmec_i, thckglcmec_i)

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
subroutine mkglacier(ldomain, mapfname, datfname, ndiag, zero_out, glac_o)
!
! !DESCRIPTION:
! make percent glacier
!
! !USES:
  use mkfileutils, only : getfil
  use mkdomainMod, only : domain1_type, domain1_clean, domain1_read
  use mkgridmapMod
  use mkvarpar	
  use mkvarctl    
  use mkncdio
!
! !ARGUMENTS:
  implicit none
  type(domain1_type), intent(in) :: ldomain
  character(len=*)  , intent(in) :: mapfname  ! input mapping file name
  character(len=*)  , intent(in) :: datfname  ! input data file name
  integer           , intent(in) :: ndiag     ! unit number for diag out
  logical           , intent(in) :: zero_out  ! if should zero glacier out
  real(r8)          , intent(out):: glac_o(:) ! output grid: %glacier
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
  type(domain1_type)    :: tdomain            ! local domain
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
  character(len=256) locfn                    ! local dataset file name
  character(len=32) :: subname = 'mkglacier'
!-----------------------------------------------------------------------

  write (6,*) 'Attempting to make %glacier .....'
  call shr_sys_flush(6)

  ! -----------------------------------------------------------------
  ! Read input file
  ! -----------------------------------------------------------------

  ! Obtain input grid info, read local fields

  call getfil (datfname, locfn, 0)

  call domain1_read(tdomain,locfn)
  ns = tdomain%ns
  allocate(glac_i(ns), stat=ier)
  if (ier/=0) call abort()

  call check_ret(nf_open(locfn, 0, ncid), subname)
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
     
     if (tdomain%ns /= tgridmap%na) then
        write(6,*)'input domain size and gridmap source size are not the same size'
        write(6,*)' domain size = ',tdomain%ns
        write(6,*)' map src size= ',tgridmap%na
        stop
     end if
     do n = 1,tgridmap%ns
        ni = tgridmap%src_indx(n)
        if (tdomain%mask(ni) /= tgridmap%mask_src(ni)) then
           write(6,*)'input domain mask and gridmap mask are not the same at ni = ',ni
           write(6,*)' domain  mask= ',tdomain%mask(ni)
           write(6,*)' gridmap mask= ',tgridmap%mask_src(ni)
           stop
        end if
        if (tdomain%lonc(ni) /= tgridmap%xc_src(ni)) then
           write(6,*)'input domain lon and gridmap lon not the same at ni = ',ni
           write(6,*)' domain  lon= ',tdomain%lonc(ni)
           write(6,*)' gridmap lon= ',tgridmap%xc_src(ni)
           stop
        end if
        if (tdomain%latc(ni) /= tgridmap%yc_src(ni)) then
           write(6,*)'input domain lat and gridmap lat not the same at ni = ',ni
           write(6,*)' domain  lat= ',tdomain%latc(ni)
           write(6,*)' gridmap lat= ',tgridmap%yc_src(ni)
           stop
        end if
     end do

     ! Determine glac_o on output grid

     call gridmap_areaave(tgridmap, glac_i, glac_o)

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

  call domain1_clean(tdomain) 
  call gridmap_clean(tgridmap)
  deallocate (glac_i)

  write (6,*) 'Successfully made %glacier'
  write (6,*)
  call shr_sys_flush(6)

end subroutine mkglacier

!-----------------------------------------------------------------------

end module mkglcmecMod
