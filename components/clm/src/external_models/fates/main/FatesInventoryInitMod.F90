module FatesInventoryInitMod

   !-----------------------------------------------------------------------------------------------
   ! This module contains the majority of the code used to read forest inventory data and
   ! initialize a simulation from that data.  Some important points:
   ! - This procedure is called through the host model's "cold-start" initialization and not a 
   !   restart type of simulation.
   ! - This procedure, if called from a massive grid, is both probably inappropriate, and probably
   !   time consuming.
   ! - This procedure is assumed to be called over a small subset of sites, for instance a single
   !   site, or a small collection of sparse/irregularly spaced group of sites
   !
   ! Created: Ryan Knox June 2017
   ! This code borrows heavily in concept from what is done in ED2.  We will also do our best to
   ! maintain compatibility with the PSS/CSS file formats that were used in ED2.
   ! See: https://github.com/EDmodel/ED2/blob/master/ED/src/io/ed_read_ed10_20_history.f90
   ! At the time of writing this ED2 is unlicensed, and only concepts were borrowed with no direct
   ! code copied.
   !-----------------------------------------------------------------------------------------------

   ! CIME GLOBALS

   use shr_log_mod      , only : errMsg => shr_log_errMsg

   ! FATES GLOBALS
   use FatesConstantsMod, only : r8 => fates_r8
   use FatesConstantsMod, only : pi_const
   use FatesGlobals     , only : endrun => fates_endrun
   use FatesGlobals     , only : fates_log
   use FatesInterfaceMod, only : bc_in_type
   use FatesInterfaceMod, only : hlm_inventory_ctrl_file
   use EDTypesMod       , only : ed_site_type
   use EDTypesMod       , only : ed_patch_type
   use EDTypesMod       , only : ed_cohort_type 
   use EDTypesMod       , only : area
   use EDPftvarcon      , only : EDPftvarcon_inst

   implicit none
   private

   ! This derived type is to allow an array of pointers to the LL patch structure
   ! This is different than allocating a vector of patches. This is needed for
   ! quickly matching cohort string identifiers, the indices that match thos identifiers
   ! with a patch.  BY having a vector of patch pointers that lines up with the string
   ! identifier array, this can be done quickly.
   type pp_array
      type(ed_patch_type), pointer :: cpatch
   end type pp_array

   character(len=*), parameter, private :: sourcefile =  __FILE__

   logical, parameter :: debug_inv = .false.       ! Debug flag for devs

   ! String length specifiers
   integer, parameter :: patchname_strlen = 64   
   integer, parameter :: cohortname_strlen = 64
   integer, parameter :: line_strlen = 512
   integer, parameter :: path_strlen = 256


   real(r8), parameter :: max_site_adjacency_deg = 0.05_r8 ! The maximum distance in degrees
                                                           ! allowed between a site's coordinate
                                                           ! defined in model memory and a physical
                                                           ! site listed in the file

   public :: initialize_sites_by_inventory

contains

   ! ==============================================================================================

   subroutine initialize_sites_by_inventory(nsites,sites,bc_in)

      ! !USES:
      use shr_file_mod, only        : shr_file_getUnit
      use shr_file_mod, only        : shr_file_freeUnit
      use EDTypesMod, only          : nclmax
      use EDTypesMod, only          : maxpft
      use EDTypesMod, only          : ncwd
      use EDPatchDynamicsMod, only  : create_patch
      use EDPatchDynamicsMod, only  : fuse_patches
      use EDCohortDynamicsMod, only : fuse_cohorts
      use EDCohortDynamicsMod, only : sort_cohorts
      use EDcohortDynamicsMod, only : count_cohorts

      ! Arguments
      integer,            intent(in)               :: nsites
      type(ed_site_type), intent(inout), target    :: sites(nsites)
      type(bc_in_type),   intent(in)               :: bc_in(nsites)

      ! Locals
      type(ed_patch_type), pointer                 :: currentpatch
      type(ed_cohort_type), pointer                :: currentcohort
      type(ed_patch_type), pointer                 :: newpatch
      type(ed_patch_type), pointer                 :: olderpatch
      integer                                      :: sitelist_file_unit   ! fortran file unit for site list
      integer                                      :: pss_file_unit        ! fortran file unit for the pss file
      integer                                      :: css_file_unit        ! fortran file unit for the css file
      integer                                      :: nfilesites           ! number of sites in file list
      logical                                      :: lopen                ! logical, file is open
      logical                                      :: lexist               ! logical, file exists
      integer                                      :: ios                  ! integer, "IO" status
      character(len=line_strlen)                   :: header_str           ! large string for whole lines
      real(r8)                                     :: age_init             ! dummy value for creating a patch
      real(r8)                                     :: area_init            ! dummy value for creating a patch
      real(r8)                                     :: cwd_ag_init(ncwd)    ! dummy value for creating a patch
      real(r8)                                     :: cwd_bg_init(ncwd)    ! dummy value for creating a patch
      real(r8)                                     :: leaf_litter_init(maxpft) ! dummy value for creating a patch
      real(r8)                                     :: root_litter_init(maxpft) ! dummy value for creating a patch
      integer                                      :: s                    ! site index
      integer                                      :: ipa                  ! patch index
      integer                                      :: total_cohorts        ! cohort counter for error checking
      integer,                         allocatable :: inv_format_list(:)   ! list of format specs
      character(len=path_strlen),      allocatable :: inv_css_list(:)      ! list of css file names
      character(len=path_strlen),      allocatable :: inv_pss_list(:)      ! list of pss file names
      real(r8),                        allocatable :: inv_lat_list(:)      ! list of lat coords
      real(r8),                        allocatable :: inv_lon_list(:)      ! list of lon coords
      integer                                      :: invsite              ! index of inventory site 
                                                                           ! closest to actual site
      character(len=patchname_strlen)              :: patch_name           ! patch ID string in the file
      integer                                      :: npatches             ! number of patches found in PSS
      type(pp_array),                  allocatable :: patch_pointer_vec(:) ! vector of pointers to patch LL
      character(len=patchname_strlen), allocatable :: patch_name_vec(:)    ! vector of patch ID strings
      real(r8)                                     :: basal_area_postf     ! basal area before fusion (m2/ha)
      real(r8)                                     :: basal_area_pref      ! basal area after fusion (m2/ha)

      real(r8), parameter                          :: max_ba_diff = 1.0e-2 ! 1% is the maximum allowable
                                                                           ! change in BA due to fusion

      ! I. Load the inventory list file, do some file handle checks
      ! ------------------------------------------------------------------------------------------

      sitelist_file_unit = shr_file_getUnit()
      inquire(file=trim(hlm_inventory_ctrl_file),exist=lexist,opened=lopen)
      if( .not.lexist ) then   ! The inventory file list DNE
         write(fates_log(), *) 'An inventory Initialization was requested.'
         write(fates_log(), *) 'However the inventory file: ',trim(hlm_inventory_ctrl_file),' DNE'
         write(fates_log(), *) 'Aborting'
         call endrun(msg=errMsg(sourcefile, __LINE__))
      end if
      if( lopen ) then        ! The inventory file should not be open
         write(fates_log(), *) 'The inventory list file is open but should not be.'
         write(fates_log(), *) 'Aborting.'
         call endrun(msg=errMsg(sourcefile, __LINE__))
      end if

      open(unit=sitelist_file_unit,file=trim(hlm_inventory_ctrl_file),status='OLD',action='READ',form='FORMATTED')
      rewind(sitelist_file_unit)

      ! There should be at least 1 line
      read(sitelist_file_unit,fmt='(A)',iostat=ios) header_str
      read(sitelist_file_unit,fmt='(A)',iostat=ios) header_str
      if( ios /= 0 ) then
         write(fates_log(), *) 'The inventory file does not contain at least two lines'
         write(fates_log(), *) 'of data, ie a header and 1 site.  Aborting.'
         call endrun(msg=errMsg(sourcefile, __LINE__))
      end if
      rewind(unit=sitelist_file_unit)


      ! Count the number of sites that are listed in this file, and allocate storage arrays
      ! ------------------------------------------------------------------------------------------

      nfilesites = count_inventory_sites(sitelist_file_unit)

      allocate(inv_format_list(nfilesites))
      allocate(inv_pss_list(nfilesites))
      allocate(inv_css_list(nfilesites))
      allocate(inv_lat_list(nfilesites))
      allocate(inv_lon_list(nfilesites))


      ! Check through the sites that are listed and do some sanity checks
      ! ------------------------------------------------------------------------------------------
      call assess_inventory_sites(sitelist_file_unit, nfilesites, inv_format_list, &
            inv_pss_list, inv_css_list, &
            inv_lat_list, inv_lon_list)

      ! We can close the list file now
      close(sitelist_file_unit, iostat = ios)
      if( ios /= 0 ) then
         write(fates_log(), *) 'The inventory file needed to be closed, but was still open'
         write(fates_log(), *) 'aborting'
         call endrun(msg=errMsg(sourcefile, __LINE__))
      end if
      call shr_file_freeUnit(sitelist_file_unit)


      ! For each site, identify the most proximal PSS/CSS couplet, read-in the data
      ! allocate linked lists and assign to memory
      do s = 1, nsites
         invsite = &
               minloc( (sites(s)%lat-inv_lat_list(:))**2.0_r8 + &
               (sites(s)%lon-inv_lon_list(:))**2.0_r8 , dim=1)

         ! Do a sanity check on the distance separation between physical site and model site
         if ( sqrt( (sites(s)%lat-inv_lat_list(invsite))**2.0_r8 + &
               (sites(s)%lon-inv_lon_list(invsite))**2.0_r8 ) > max_site_adjacency_deg ) then
            write(fates_log(), *) 'Model site at lat:',sites(s)%lat,' lon:',sites(s)%lon
            write(fates_log(), *) 'has no reasonably proximal site in the inventory site list.'
            write(fates_log(), *) 'Closest is at lat:',inv_lat_list(invsite),' lon:',inv_lon_list(invsite)
            write(fates_log(), *) 'Separation must be less than ',max_site_adjacency_deg,' degrees'
            write(fates_log(), *) 'Exiting'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if


         ! Open the PSS/CSS couplet and initialize the ED data structures.
         ! Lets start withe the PSS
         ! ---------------------------------------------------------------------------------------

         pss_file_unit = shr_file_getUnit()
         open(unit=pss_file_unit,file=trim(inv_pss_list(invsite)), &
               status='OLD',action='READ',form='FORMATTED')
         rewind(pss_file_unit)
         read(pss_file_unit,fmt=*) header_str

         ! Do one quick pass through just to count lines
         ipa = 0
         countpatchloop: do
            read(pss_file_unit,fmt=*,iostat=ios) header_str
            if(ios/=0) exit
            ipa = ipa + 1
         end do countpatchloop
         rewind(pss_file_unit)
         read(pss_file_unit,fmt=*) header_str

         npatches = ipa
         allocate(patch_name_vec(npatches))
         allocate(patch_pointer_vec(npatches))


         do ipa=1,npatches

            allocate(newpatch)

            newpatch%patchno = ipa
            newpatch%younger => null()
            newpatch%older   => null()

            ! This call doesn't do much asside from initializing the patch with
            ! nominal values, NaNs, zero's and allocating some vectors. We should
            ! be able to get the following values from the patch files. But on 
            ! the patch creation step, we don't have that information.

            age_init            = 0.0_r8
            area_init           = 0.0_r8 
            cwd_ag_init(:)      = 0.0_r8
            cwd_bg_init(:)      = 0.0_r8
            leaf_litter_init(:) = 0.0_r8
            root_litter_init(:) = 0.0_r8

            call create_patch(sites(s), newpatch, age_init, area_init, &
                  cwd_ag_init, cwd_bg_init, &
                  leaf_litter_init, root_litter_init )

            if( inv_format_list(invsite) == 1 ) then
               call set_inventory_edpatch_type1(newpatch,pss_file_unit,ipa,ios,patch_name)
            end if

            ! Add it to the site's patch list
            ! ------------------------------------------------------------------------------------

            patch_name_vec(ipa)           = trim(patch_name)
            patch_pointer_vec(ipa)%cpatch => newpatch

            if(ipa == 1) then
               ! This is the first patch to be added
               ! It starts off as the oldest and youngest patch in the list
               sites(s)%youngest_patch => newpatch
               sites(s)%oldest_patch   => newpatch
            else

               ! Insert this patch into the patch LL
               ! First check the two end-cases

               ! Youngest Patch
               if(newpatch%age<=sites(s)%youngest_patch%age)then
                  newpatch%older                  => sites(s)%youngest_patch
                  newpatch%younger                => null()
                  sites(s)%youngest_patch%younger => newpatch
                  sites(s)%youngest_patch         => newpatch

                  ! Oldest Patch
               else if(newpatch%age>sites(s)%oldest_patch%age)then
                  newpatch%older              => null()
                  newpatch%younger            => sites(s)%oldest_patch
                  sites(s)%oldest_patch%older => newpatch
                  sites(s)%oldest_patch       => newpatch

                  ! Somewhere in the middle
               else
                  currentpatch => sites(s)%youngest_patch
                  do while(associated(currentpatch))
                     olderpatch => currentpatch%older
                     if(associated(currentpatch%older)) then
                        if(newpatch%age >= currentpatch%age .and. &
                              newpatch%age < olderpatch%age) then
                           ! Set the new patches pointers
                           newpatch%older   => currentpatch%older
                           newpatch%younger => currentpatch
                           ! Fix the patch's older pointer
                           currentpatch%older => newpatch
                           ! Fix the older patch's younger pointer
                           olderpatch%younger => newpatch
                        end if
                     end if
                     currentPatch => olderpatch
                  enddo
               end if
            end if
         end do

         close(pss_file_unit,iostat=ios)
         if( ios /= 0 ) then
            write(fates_log(), *) 'The pss file: ',inv_pss_list(invsite),' could not be closed'
            write(fates_log(), *) 'aborting'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if
         call shr_file_freeUnit(pss_file_unit)

         if(debug_inv) then
            write(fates_log(),*) 'Raw List of Inventory Patches, Age Sorted:'
            currentpatch => sites(s)%youngest_patch
            do while(associated(currentpatch))
               write(fates_log(),*) ' AGE: ',currentpatch%age,' AREA: ',currentpatch%area
               currentPatch => currentpatch%older
            enddo
         end if
         
         ! OPEN THE CSS FILE
         ! ---------------------------------------------------------------------------------------
         css_file_unit = shr_file_getUnit()
         open(unit=css_file_unit,file=trim(inv_css_list(invsite)), &
               status='OLD',action='READ',form='FORMATTED')
         rewind(css_file_unit)
         read(css_file_unit,fmt=*) header_str

         ! Read in each cohort line. Each line is associated with a patch from the PSS
         ! file via a patch name identification string.  We pass the whole site pointer
         ! to this routine, because inside the routine we identify the patch by making
         ! comparisons with patch_name_vec and identifying the patch pointer 
         ! from patch_pointer_vec

         invcohortloop: do
            if ( inv_format_list(invsite) == 1 ) then
               call set_inventory_edcohort_type1(sites(s),bc_in(s),css_file_unit, &
                     npatches, patch_pointer_vec,patch_name_vec, ios)
            end if
            if ( ios/=0 ) exit
         end do invcohortloop

         close(css_file_unit,iostat=ios)
         if( ios/=0 ) then
            write(fates_log(), *) 'The css file: ',inv_css_list(invsite),' could not be closed'
            write(fates_log(), *) 'aborting'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if
         call shr_file_freeUnit(css_file_unit)

         deallocate(patch_pointer_vec,patch_name_vec)

         ! Report Basal Area (as a check on if things were read in)
         ! ------------------------------------------------------------------------------
         basal_area_pref = 0.0_r8
         currentpatch => sites(s)%youngest_patch
         do while(associated(currentpatch))
            currentcohort => currentpatch%tallest
            do while(associated(currentcohort))
               basal_area_pref = basal_area_pref + &
                     currentcohort%n*0.25*((currentcohort%dbh/100.0_r8)**2.0_r8)*pi_const
               currentcohort => currentcohort%shorter
            end do
            currentPatch => currentpatch%older
         enddo

         write(fates_log(),*) '-------------------------------------------------------'
         write(fates_log(),*) 'Basal Area from inventory, BEFORE fusion'
         write(fates_log(),*) 'Lat: ',sites(s)%lat,' Lon: ',sites(s)%lon
         write(fates_log(),*) basal_area_pref,' [m2/ha]'
         write(fates_log(),*) '-------------------------------------------------------'

         ! Update the patch index numbers and fuse the cohorts in the patches
         ! ----------------------------------------------------------------------------------------
         ipa=1
         total_cohorts = 0
         currentpatch => sites(s)%youngest_patch
         do while(associated(currentpatch))
            currentpatch%patchno = ipa
            ipa=ipa+1
            
            ! Perform Cohort Fusion
            call fuse_cohorts(currentpatch,bc_in(s))
            call sort_cohorts(currentpatch)
            total_cohorts = total_cohorts + count_cohorts(currentpatch)

            currentPatch => currentpatch%older
         enddo

         if(total_cohorts .eq. 0)then
            write(fates_log(), *) 'The inventory initialization produced no cohorts.'
            write(fates_log(), *) 'aborting'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         ! Fuse patches
         ! ----------------------------------------------------------------------------------------
         call fuse_patches(sites(s), bc_in(s) ) 

         ! Report Basal Area (as a check on if things were read in)
         ! ----------------------------------------------------------------------------------------
         basal_area_postf = 0.0_r8
         currentpatch => sites(s)%youngest_patch
         do while(associated(currentpatch))
            currentcohort => currentpatch%tallest
            do while(associated(currentcohort))
               basal_area_postf = basal_area_postf + &
                     currentcohort%n*0.25*((currentcohort%dbh/100.0_r8)**2.0_r8)*pi_const
               currentcohort => currentcohort%shorter
            end do
            currentPatch => currentpatch%older
         enddo

         write(fates_log(),*) '-------------------------------------------------------'
         write(fates_log(),*) 'Basal Area from inventory, AFTER fusion'
         write(fates_log(),*) 'Lat: ',sites(s)%lat,' Lon: ',sites(s)%lon
         write(fates_log(),*) basal_area_postf,' [m2/ha]'
         write(fates_log(),*) '-------------------------------------------------------'

         ! Check to see if the fusion process has changed too much
         ! We are sensitive to fusion in inventories because we may be asking for a massive amount
         ! of fusion. For instance some init files are directly from inventory, where a cohort
         ! is synomomous with a single plant.

         if( abs(basal_area_postf-basal_area_pref)/basal_area_pref > max_ba_diff ) then
            write(fates_log(),*) 'Inventory Fusion Changed total biomass beyond reasonable limit'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if
         
      end do

      deallocate(inv_format_list, inv_pss_list, inv_css_list, inv_lat_list, inv_lon_list)

      return
   end subroutine initialize_sites_by_inventory

   ! ==============================================================================================

   function count_inventory_sites(sitelist_file_unit) result(nsites)

      ! Simple function that counts the number of lines in the inventory descriptor file

      ! Arguments
      integer, intent(in)        :: sitelist_file_unit

      ! Locals
      character(len=line_strlen) :: header_str
      character(len=line_strlen) :: site_str
      integer                    :: ios
      integer                    :: nsites


      ! Set the file position to the top of the file
      ! Read in the header line
      ! Read through sites and check coordinates and file existence
      rewind(unit=sitelist_file_unit)
      read(sitelist_file_unit,fmt='(A)') header_str
      nsites = 0
      do
         read(sitelist_file_unit,fmt='(A)',iostat=ios) site_str
         if (ios/=0) exit
         nsites = nsites + 1
      end do

      return
   end function count_inventory_sites

   ! ==============================================================================================

   subroutine assess_inventory_sites(sitelist_file_unit,nsites, inv_format_list, &
         inv_pss_list,inv_css_list, &
         inv_lat_list,inv_lon_list)

      ! -------------------------------------------------------------------------------------------
      ! This subroutine looks through the inventory descriptor file
      ! and line by line reads information about the available inventory
      ! sites, and saves their information (such as location and file path)
      ! to arrays.  This routine also does some simple checks to make
      ! sure it is not reading nonsense
      !
      ! File Format for the inventory site file:
      ! 1 line header
      ! 1 line listing each available inventory site with the following fields:
      ! type     latitude    longitude   pss-name   css-name
      !
      ! The fields for each site are described as follows:
      !
      ! <short-name>    <value-kind>     <description>
      !
      ! type            integer          We will accomodate different file format with different
      !                                  field values as the need arises. format 1 will read in 
      !                                  datasets via "set_inventory_edpatch_type1()", 
      !                                  "set_inventory_edcohort_type1()"
      ! 
      ! latitude        float            The geographic latitude coordinate of the site
      ! longitude       float            The geogarphic longitude coordinate of the site
      ! pss-name        string           The full path to the patch descriptor file (PSS)
      ! css-name        string           The full path to the cohort descriptor file (CSS)
      ! -------------------------------------------------------------------------------------------


      ! Arguments
      integer, intent(in)                      :: sitelist_file_unit       ! file unit for sitelist
      integer, intent(in)                      :: nsites                   ! number of inventory sites
      integer, intent(inout)                   :: inv_format_list(nsites)  ! array of formats for each inventory site
      character(len=path_strlen),intent(inout) :: inv_pss_list(nsites)     ! array of pss file paths for each site
      character(len=path_strlen),intent(inout) :: inv_css_list(nsites)     ! array of css file paths for each site
      real(r8),intent(inout)                   :: inv_lat_list(nsites)     ! array of latitudes for each site
      real(r8),intent(inout)                   :: inv_lon_list(nsites)     ! array of longitudes for each site

      ! Locals
      character(len=line_strlen)               :: header_str  ! a string to hold the header information
      character(len=line_strlen)               :: site_str    ! a string to hold each site-line in the file
      integer                                  :: isite       ! site index
      integer                                  :: ios         ! fortran read status flag
      character(len=path_strlen)               :: pss_file
      character(len=path_strlen)               :: css_file
      real(r8)                                 :: site_lat    ! inventory site latitude
      real(r8)                                 :: site_lon    ! site longitude
      integer                                  :: iblnk       ! Index used for string parsing
      integer                                  :: file_format ! format type (1=legacy ED pss/css)
      logical                                  :: lexist      ! file existence flag

      rewind(unit=sitelist_file_unit)
      read(sitelist_file_unit,fmt='(4A)') header_str

      do isite=1,nsites

         ! Read in the whole line
         read(sitelist_file_unit,fmt='(a)',iostat=ios) site_str

         ! Parse the format identifier
         read(site_str,*) file_format

         ! Parse the latitude
         site_str = adjustl(site_str)
         iblnk = index(site_str,' ')
         site_str = adjustl(site_str(iblnk:))
         read(site_str,*) site_lat

         ! Parse the longitude
         site_str = adjustl(site_str)
         iblnk = index(site_str,' ')
         site_str = adjustl(site_str(iblnk:))
         read(site_str,*) site_lon

         ! Parse the pss file name
         site_str = adjustl(site_str)
         iblnk    = index(site_str,' ')
         site_str = adjustl(site_str(iblnk:))
         iblnk    = index(site_str,' ')
         read(site_str(:iblnk),fmt='(1A)') pss_file

         ! Parse the css file name
         site_str = adjustl(site_str)
         iblnk    = index(site_str,' ')
         site_str = adjustl(site_str(iblnk:))
         read(site_str,fmt='(1A)') css_file


         if ( site_lat < -90.0_r8 .or. site_lat > 90.0_r8 ) then
            write(fates_log(), *) 'read invalid latitude: ',site_lat,' from inventory site list'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         ! Longitude should be converted to positive coordinate

         if( site_lon<0.0_r8 ) site_lon = 360.0_r8 + site_lon

         if ( site_lon < 0.0_r8 .or. site_lon > 360.0_r8 ) then
            write(fates_log(), *) 'read invalid longitude: ',site_lon,' from inventory site list'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         inquire(file=trim(pss_file),exist=lexist)
         if( .not.lexist ) then
            write(fates_log(), *) 'the following pss file could not be found:'
            write(fates_log(), *) trim(pss_file)
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         inquire(file=trim(css_file),exist=lexist)
         if( .not.lexist ) then
            write(fates_log(), *) 'the following css file could not be found:'
            write(fates_log(), *) trim(css_file)
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         ! If we have made it to this point, then in all likelihood, the PSS/CSS
         ! File has probably been correctly identified

         inv_format_list(isite) = file_format
         inv_pss_list(isite) = pss_file
         inv_css_list(isite) = css_file
         inv_lat_list(isite) = site_lat
         inv_lon_list(isite) = site_lon

      end do

      return
   end subroutine assess_inventory_sites

   ! ==============================================================================================

   subroutine set_inventory_edpatch_type1(newpatch,pss_file_unit,ipa,ios,patch_name)

      ! --------------------------------------------------------------------------------------------
      ! This subroutine reads in a line of an inventory patch file (pss)
      ! And populates a new patch with that information.
      ! This routine specifically reads PSS files that are "Type 1" formatted
      ! 
      ! The file is formatted text, which contains 1 header line to label columns
      ! and then 1 line for each patch containing the following fields:
      !
      ! time	(year)     year of measurement
      ! patch	(string)   patch id string
      ! trk	(integer)  LU type index (0 non-forest, 1 secondary, 2 primary
      ! age	(years)    Time since this patch was disturbed (created)
      ! area	(fraction) Fraction of the site occupied by this patch
      ! water	(NA)       Water content of soil (NOT USED)
      ! fsc	(kg/m2)    Fast Soil Carbon
      ! stsc	(kg/m2)    Structural Soil Carbon
      ! stsl	(kg/m2)    Structural Soil Lignan
      ! ssc	(kg/m2)    Slow Soil Carbon
      ! psc	(NA)       Passive Soil Carbon (NOT USED)
      ! msn	(kg/m2)    Mineralized Soil Nitrogen
      ! fsn     (kg/m2)    Fast Soil Nitrogen
      ! --------------------------------------------------------------------------------------------

      use FatesSizeAgeTypeIndicesMod, only: get_age_class_index
      use EDtypesMod, only: AREA
      use EDTypesMod, only: ncwd
      use SFParamsMod , only : SF_val_CWD_frac

      ! Arguments
      type(ed_patch_type),intent(inout), target   :: newpatch      ! Patch structure
      integer,intent(in)                          :: pss_file_unit ! Self explanatory
      integer,intent(in)                          :: ipa           ! Patch index (line number)
      integer,intent(out)                         :: ios           ! Return flag
      character(len=patchname_strlen),intent(out) :: patch_name    ! unique string identifier of patch

      ! Locals
      real(r8)                                    :: p_time     ! Time patch was recorded
      real(r8)                                    :: p_trk      ! Land Use index (see above descriptions)
      character(len=patchname_strlen)             :: p_name     ! unique string identifier of patch
      real(r8)                                    :: p_age      ! Patch age [years]
      real(r8)                                    :: p_area     ! Patch area [fraction]
      real(r8)                                    :: p_water    ! Patch water (unused)
      real(r8)                                    :: p_fsc      ! Patch fast soil carbon
      real(r8)                                    :: p_stsc     ! Patch structural soil carbon
      real(r8)                                    :: p_stsl     ! Patch structural soil lignans
      real(r8)                                    :: p_ssc      ! Patch slow soil carbon
      real(r8)                                    :: p_psc      ! Patch P soil carbon
      real(r8)                                    :: p_msn      ! Patch mean soil nitrogen
      real(r8)                                    :: p_fsn      ! Patch fast soil nitrogen
      integer                                     :: icwd       ! index for counting CWD pools
      integer                                     :: ipft       ! index for counting PFTs
      real(r8)                                    :: pftfrac    ! the inverse of the total number of PFTs

      character(len=128),parameter    :: wr_fmt = &
            '(F5.2,2X,A4,2X,F5.2,2X,F5.2,2X,F5.2,2X,F5.2,2X,F5.2,2X,F5.2,2X,F5.2,2X,F5.2,2X,F5.2,2X,F5.2,2X,F5.2)'


      read(pss_file_unit,fmt=*,iostat=ios) p_time, p_name, p_trk, p_age, p_area, &
            p_water,p_fsc, p_stsc, p_stsl, p_ssc, &
            p_psc, p_msn, p_fsn

      if (ios/=0) return

      patch_name = trim(p_name)

      if( debug_inv) then
         write(*,fmt=wr_fmt) &
               p_time, p_name, p_trk, p_age, p_area, &
               p_water,p_fsc, p_stsc, p_stsl, p_ssc, &
               p_psc, p_msn, p_fsn
      end if

      ! Fill in the patch's memory structures

      newpatch%age       = p_age
      newpatch%age_class = get_age_class_index(newpatch%age)
      newpatch%area      = p_area*AREA

      ! The non-litter patch soil variables need to be sent to the HLM
      ! This is not being sent as of this message (RGK 06-2017)

      ! ---------------------------------------------------------------------
      ! The litter and CWD could be estimated from at least two methods
      ! 1) after reading in the cohort data, assuming a SS turnover rate
      ! 2) again assuming SS, back out the CWD and Litter flux rates into
      !    the SSC, STSC and FSC pools that balance with their decomp rates
      !    and then use those flux rates, to calculate the CWD and litter
      !    pool sizes based on another SS model where flux out balances with
      !    mortality and litter fluxes into the non-decomposing pool
      !
      ! This is significant science modeling and does not have a simple
      ! first hack solution. (RGK 06-2017)
      ! ----------------------------------------------------------------------

      do icwd = 1, ncwd
         newpatch%cwd_ag(icwd) = 0.0_r8
         newpatch%cwd_bg(icwd) = 0.0_r8
      end do


      newpatch%leaf_litter(:) = 0.0_r8
      newpatch%root_litter(:) = 0.0_r8


      return
   end subroutine set_inventory_edpatch_type1


   ! ==============================================================================================

   subroutine set_inventory_edcohort_type1(csite,bc_in,css_file_unit,npatches, &
                                           patch_pointer_vec,patch_name_vec,ios)

      ! --------------------------------------------------------------------------------------------
      ! This subroutine reads in a line of an inventory cohort file (css)
      ! And populates a new cohort with that information.
      ! This routine specifically reads CSS files that are "Type 1" formatted
      ! 
      ! The file formatted text, which contains 1 header line to label columns
      ! and then 1 line for each cohort containing the following fields:
      ! 
      ! FILE FORMAT:
      ! time	(year)     year of measurement
      ! patch	(string)   patch id string associated with this cohort
      ! index	(integer)  cohort index
      ! dbh	(cm)       diameter at breast height
      ! height   (m)        height of the tree
      ! pft      (integer)  the plant functional type index (must be consistent with param file)
      ! n 	(/m2)      The plant number density
      ! bdead    (kgC/plant)The dead biomass per indiv of this cohort (NOT USED)
      ! balive   (kgC/plant)The live biomass per indiv of this cohort (NOT USED)
      ! avgRG    (cm/yr?)   Average Radial Growth (NOT USED)
      ! --------------------------------------------------------------------------------------------

      use EDGrowthFunctionsMod, only : hite
      use EDGrowthFunctionsMod, only : bleaf
      use EDGrowthFunctionsMod, only : bdead
      use EDCohortDynamicsMod , only : create_cohort
      use FatesInterfaceMod   , only : numpft

      ! Arguments
      type(ed_site_type),intent(inout), target    :: csite         ! current site
      type(bc_in_type),intent(in)                 :: bc_in         ! boundary conditions
      integer, intent(in)                         :: css_file_unit     ! Self explanatory
      integer, intent(in)                         :: npatches      ! number of patches
      type(pp_array), intent(in)                  :: patch_pointer_vec(npatches)
      character(len=patchname_strlen), intent(in) :: patch_name_vec(npatches)
      integer,intent(out)                         :: ios           ! Return flag
      
      ! Locals
      real(r8)                                    :: c_time        ! Time patch was recorded
      character(len=patchname_strlen)             :: p_name        ! The patch associated with this cohort
      character(len=cohortname_strlen)            :: c_name        ! cohort index
      real(r8)                                    :: c_dbh         ! diameter at breast height (cm)
      real(r8)                                    :: c_height      ! tree height (m)
      integer                                     :: c_pft         ! plant functional type index
      real(r8)                                    :: c_nplant      ! plant density (/m2)
      real(r8)                                    :: c_bdead       ! dead biomass (kg)
      real(r8)                                    :: c_balive      ! live biomass (kg)
      real(r8)                                    :: c_avgRG       ! avg radial growth (NOT USED)
      integer                                     :: cstatus       ! 
      type(ed_patch_type), pointer                :: cpatch        ! current patch pointer
      type(ed_cohort_type), pointer               :: temp_cohort   ! temporary patch (needed for allom funcs)
      integer                                     :: ipa           ! patch idex
      logical                                     :: matched_patch ! check if cohort was matched w/ patch

      character(len=128),parameter    :: wr_fmt = &
           '(F7.1,2X,A20,2X,A20,2X,F5.2,2X,F5.2,2X,I4,2X,F5.2,2X,F5.2,2X,F5.2,2X,F5.2)'

      real(r8), parameter :: abnormal_large_nplant = 1000.0_r8  ! Used to catch bad values
      real(r8), parameter :: abnormal_large_dbh    = 500.0_r8   ! I've never heard of a tree > 3m

      read(css_file_unit,fmt=*,iostat=ios) c_time, p_name, c_name, c_dbh, c_height, &
            c_pft, c_nplant, c_bdead, c_balive, c_avgRG
      
      if( debug_inv) then
         write(*,fmt=wr_fmt) &
              c_time, p_name, c_name, c_dbh, c_height, &
              c_pft, c_nplant, c_bdead, c_balive, c_avgRG
      end if

      if (ios/=0) return

      ! Identify the patch based on the patch_name
      matched_patch = .false.
      do ipa=1,npatches
         if( trim(p_name) == trim(patch_name_vec(ipa))) then
            cpatch => patch_pointer_vec(ipa)%cpatch
            matched_patch = .true.
         end if
      end do

      if(.not.matched_patch)then
         write(fates_log(), *) 'could not match a cohort with a patch'
         write(fates_log(),fmt=wr_fmt) &
               c_time, p_name, c_name, c_dbh, c_height, &
               c_pft, c_nplant, c_bdead, c_balive, c_avgRG
         call endrun(msg=errMsg(sourcefile, __LINE__))
      end if

      ! Run some sanity checks on the input data
      ! pft, nplant and dbh are the critical ones in this format specification
      ! -------------------------------------------------------------------------------------------

      if (c_pft > numpft ) then
         write(fates_log(), *) 'inventory pft: ',c_pft
         write(fates_log(), *) 'An inventory cohort file specified a pft index'
         write(fates_log(), *) 'greater than the maximum specified pfts numpft'
         call endrun(msg=errMsg(sourcefile, __LINE__))
      end if

      if (c_pft <= 0 ) then
         write(fates_log(), *) 'inventory pft: ',c_pft
         write(fates_log(), *) 'The inventory produced a cohort with <=0 pft index'
         call endrun(msg=errMsg(sourcefile, __LINE__))
      end if

      if (c_dbh <=0 ) then
         write(fates_log(), *) 'inventory dbh: ', c_dbh
         write(fates_log(), *) 'The inventory produced a cohort with <= 0 dbh'
         call endrun(msg=errMsg(sourcefile, __LINE__))
      end if
      
      if (c_dbh > abnormal_large_dbh ) then
         write(fates_log(), *) 'inventory dbh: ', c_nplant
         write(fates_log(), *) 'The inventory produced a cohort with very large diameter [cm]'
         call endrun(msg=errMsg(sourcefile, __LINE__))
      end if

      if (c_nplant <=0 ) then
         write(fates_log(), *) 'inventory nplant: ', c_nplant
         write(fates_log(), *) 'The inventory produced a cohort with <= 0 density /m2'
         call endrun(msg=errMsg(sourcefile, __LINE__))
      end if

      if (c_nplant > abnormal_large_nplant ) then
         write(fates_log(), *) 'inventory nplant: ', c_nplant
         write(fates_log(), *) 'The inventory produced a cohort with very large density /m2'
         call endrun(msg=errMsg(sourcefile, __LINE__))
      end if
      
      allocate(temp_cohort)   ! A temporary cohort is needed because we want to make
                              ! use of the allometry functions
      
     

      temp_cohort%pft         = c_pft
      temp_cohort%n           = c_nplant * cpatch%area
      temp_cohort%hite        = Hite(temp_cohort)
      temp_cohort%dbh         = c_dbh
      temp_cohort%canopy_trim = 1.0_r8
      temp_cohort%bdead       = Bdead(temp_cohort)
      temp_cohort%balive      = Bleaf(temp_cohort)*(1.0_r8 + EDPftvarcon_inst%allom_l2fr(c_pft) &
            + EDpftvarcon_inst%allom_latosa_int(c_pft)*temp_cohort%hite)
      temp_cohort%b           = temp_cohort%balive + temp_cohort%bdead
      
      if( EDPftvarcon_inst%evergreen(c_pft) == 1) then
         temp_cohort%bstore = Bleaf(temp_cohort) * EDPftvarcon_inst%cushion(c_pft)
         temp_cohort%laimemory = 0._r8
         cstatus = 2
      endif
      
      if( EDPftvarcon_inst%season_decid(c_pft) == 1 ) then !for dorment places
         temp_cohort%bstore = Bleaf(temp_cohort) * EDPftvarcon_inst%cushion(c_pft) !stored carbon in new seedlings.
         if(csite%status == 2)then 
            temp_cohort%laimemory = 0.0_r8
         else
            temp_cohort%laimemory = Bleaf(temp_cohort)
         endif
         ! reduce biomass according to size of store, this will be recovered when elaves com on.
         temp_cohort%balive = temp_cohort%balive - temp_cohort%laimemory
         cstatus = csite%status
      endif
      
      if ( EDPftvarcon_inst%stress_decid(c_pft) == 1 ) then
         temp_cohort%bstore = Bleaf(temp_cohort) * EDPftvarcon_inst%cushion(c_pft)
         temp_cohort%laimemory = Bleaf(temp_cohort)
         temp_cohort%balive = temp_cohort%balive - temp_cohort%laimemory
         cstatus = csite%dstatus
      endif
      
      call create_cohort(cpatch, c_pft, temp_cohort%n, temp_cohort%hite, temp_cohort%dbh, &
            temp_cohort%balive, temp_cohort%bdead, temp_cohort%bstore, &
            temp_cohort%laimemory,  cstatus, temp_cohort%canopy_trim, 1, bc_in)
      
      deallocate(temp_cohort) ! get rid of temporary cohort

      return
   end subroutine set_inventory_edcohort_type1

end module FatesInventoryInitMod
