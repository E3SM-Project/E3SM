PROGRAM main

   use shr_timer_mod
   use kind_mod

   use map_mod
   use mapsort_mod
   use fixroff_mod
   use smooth_mod

   implicit none

   !--- mapping matricies ---
   type(sMatrix),target :: map_orig    ! orig matrix -- produced by scrip, may have runoff over land
   type(sMatrix),target :: map_corr    ! orig matrix, corrected -- no runoff over land
   type(sMatrix),target :: map_new     ! orig matrix, corrected, smoothed -- unsorted & sorted
   type(sMatrix),target :: map_smooth  ! smoothing matrix

   !--- for correction (runoff relocations) and smoothing ---
   integer(IN)       :: n               !
   real(R8)          :: eFold           ! smoothing e-fold curve
   real(R8)          :: rMax            ! max smoothing radius

   integer(IN)       :: t01,t02,t03,t04,t05,t06,t07,t08,t09,t10 ! timer numbers
   integer(IN)       :: rCode           ! return code
   character( 8)     :: dstr            ! wall clock date
   character(10)     :: tstr            ! wall clock time

   !--- namelist vars ---
   character(180) :: gridtype                  ! type of run-off grid
   character(180) :: file_roff                 ! file name: rtm rdirc file
   character(180) :: file_roff_out             ! file name: rtm rdirc file
   character(180) :: file_ocn_global_mask      ! file name: ocn scrip grid file (all ocean points masked as such)
   character(180) :: file_ocn_coastal_mask     ! file name: ocn scrip grid file (only coastal points masked as ocean)
   character(180) :: file_nn                   ! file name: orig matrix, corrected
   character(180) :: file_new                  ! file name: orig matrix, corrected, smoothed, sorted -- done
   character(180) :: file_smooth               ! file name: smoothing matrix
   character(180) :: title                     ! netCDF title attribute
   logical        :: nn_dest_is_smooth_src     ! use step1 dest cells as
                                               ! step2 source cells
   logical        :: step1                     ! gen nn
   logical        :: step2                     ! gen smooth
   logical        :: step3                     ! mat mult
   logical        :: lmake_rSCRIP              ! .true. => convert runoff grid to SCRIP

   namelist / input_nml /   &
      gridtype              &
   ,  file_roff             &
   ,  file_roff_out         &
   ,  file_ocn_global_mask  &
   ,  file_ocn_coastal_mask &
   ,  file_nn               &
   ,  file_new              &
   ,  file_smooth           &
   ,  title                 &
   ,  eFold                 &
   ,  rMax                  &
   ,  lmake_rSCRIP          &
   ,  nn_dest_is_smooth_src &
   ,  step1, step2, step3

   !--- formats ---
   character(*),parameter :: F00 = "('(main) ',6a)"
   character(*),parameter :: F01 = "('(main) ',a,i7)"
   character(*),parameter :: F02 = "('(main) ',a,es13.6)"
   character(*),parameter :: F03 = "('(main) ',a,l1)"
   character(*),parameter :: F10 = "('(main) ',73('-')/,'(main) ',a/,'(main) ',73('-'))"
   character(*),parameter :: F12 = "('(main) date & time:',1x,a4,2('-',a2),2x,a2,2(':',a2))"

!-------------------------------------------------------------------------------
! PURPOSE:
!   Take a scrip conservative map matrix, runoff->ocean (unmasked ocean)
!   and produce a runoff->ocean map suitable for use in CCSM
!
! Step 1: Read raw grids and genearate nearest neighbor mapping file
!         (or read nearest neighbor map)
! Step 2: Create the smoothing map (and sort it)
!         smoothing map determines how runoff gets spread to more ocean cells
! Step 3: Multiply corrected original with smoothing (and sort new map)
!-------------------------------------------------------------------------------

   write(6,F10) "correct/smooth/sort runoff -> ocean map"

   call shr_timer_init()

   !----------------------------------------------------------------------------
   write(6,F10) "Step 0:  read input namelist data"
   !----------------------------------------------------------------------------

   gridtype              = 'unset'
   file_roff             = 'unset'
   file_ocn_global_mask  = 'unset'
   file_ocn_coastal_mask = 'unset'
   file_nn               = 'unset'
   file_smooth           = 'unset'
   file_new              = 'unset'
   title                 = 'unset'
   eFold                 = 1000000.00000 ! smoothing eFold distance in meters
   rMax                  =  500000.00000 ! max smoothing radius in meters
   nn_dest_is_smooth_src = .true.
   step1                 = .true.
   step2                 = .true.
   step3                 = .true.

   ! These two variables typically don't appear in namelist
   lmake_rSCRIP  = .false.
   file_roff_out = "runoff.nc"

   read (*,nml=input_nml,iostat=rCode)
   if (trim(file_ocn_coastal_mask) == 'unset') then
     file_ocn_coastal_mask = file_ocn_global_mask
   end if

   write(6,F00) "Namelist values..."
   write(6,F00) "   gridtype              = ",trim(gridtype      )
   write(6,F00) "   file_roff             = ",trim(file_roff     )
   write(6,F00) "   file_ocn (global)     = ",trim(file_ocn_global_mask)
   write(6,F00) "   file_ocn (coast)      = ",trim(file_ocn_coastal_mask)
   write(6,F00) "   file_nn               = ",trim(file_nn       )
   write(6,F00) "   file_smooth           = ",trim(file_smooth   )
   write(6,F00) "   file_new              = ",trim(file_new      )
   write(6,F00) "   title                 = ",trim(title         )
   write(6,F02) "   eFold distance        = ",eFold
   write(6,F02) "   rMax  distance        = ",rMax
   write(6,F03) "   nn_dest_is_smooth_src = ",nn_dest_is_smooth_src
   write(6,F03) "   step1                 = ",step1
   write(6,F03) "   step2                 = ",step2
   write(6,F03) "   step3                 = ",step3
!  if (rCode > 0) then
!     write(6,F01) 'ERROR: reading input namelist, iostat=',rCode
!     stop
!  end if

   call date_and_time(dstr,tstr)
   write(6,F12) dstr(1:4),dstr(5:6),dstr(7:8) ,tstr(1:2),tstr(3:4),tstr(5:6)

   if (lmake_rSCRIP) then
     write(6,F10) "Generating SCRIP file from runoff input"
     ! note that this overloads ofilename in map_gridRead... use as
     ! ocean scrip grid file in normal use, but output scrip grid file
     ! when generating SCRIP file from runoff data.
     call map_gridRead(map_orig, trim(file_roff), trim(file_roff_out),        &
                       gridtype, .true.)
     write(6,*) "Successfully generated ", trim(file_roff_out)
     stop
   end if

   call shr_timer_get  (t01,"Step 1: Grid Read, Gen NN")
   call shr_timer_start(t01)
if (step1) then
   !----------------------------------------------------------------------------
   write(6,F10) "Step 1: read grid info & create nearest neighbor map"
   !----------------------------------------------------------------------------
   call map_gridRead(map_orig , trim(file_roff), trim(file_ocn_coastal_mask), gridtype)

   call date_and_time(dstr,tstr)
   call map_print(map_orig)
!  call map_gennn(map_orig)     ! optimized nn map generation -- has bugs?
   call map_gennn0(map_orig)    ! non-optimized nn map generation
   call mapsort_sort(map_orig)  ! sort map
   call map_check(map_orig)
   call map_write(map_orig, trim(file_nn))
else if (step2.or.step3) then
   !----------------------------------------------------------------------------
   write(6,F10) "Step 1: read nearest neighbor map"
   !----------------------------------------------------------------------------
   call map_read (map_orig, trim(file_nn))  ! read corrected r2o map
   call mapsort_sort(map_orig)  ! sort map
   call map_check(map_orig)
endif
   call shr_timer_stop (t01)
   call shr_timer_print(t01)


   call shr_timer_get  (t06,"Step 2: Create the smoothing map")
   call shr_timer_start(t06)
if (step2) then
   !----------------------------------------------------------------------------
   write(6,F10) "Step 2:  create smoothing map"
   !----------------------------------------------------------------------------
   write(6,F02) "   eFold distance = ",eFold
   write(6,F02) "   rMax  distance = ",rMax
   call smooth_init(file_ocn_global_mask, nn_dest_is_smooth_src, map_orig,map_smooth)
!  call restrictSources(map_smooth,trim(file_sources)) ! restrict cells subject to smoothing
   call smooth(map_smooth,eFold,rMax)
   call mapsort_sort(map_smooth)
   call map_check(map_smooth)
   call map_write(map_smooth, trim(file_smooth))
else if (step3) then
   !----------------------------------------------------------------------------
   write(6,F10) "Step 2:  read smoothing map"
   !----------------------------------------------------------------------------
   call map_read (map_smooth, trim(file_smooth))  ! read corrected r2o map
   call mapsort_sort(map_smooth)  ! sort map
   call map_check(map_smooth)
endif
   call shr_timer_stop (t06)
   call shr_timer_print(t06)

   call shr_timer_get  (t07,"Step 3:  Matrix Multiply")
   call shr_timer_start(t07)


if (step3) then
   !----------------------------------------------------------------------------
   write(6,F10) "Step 3:  create smoothed version of nearest neighbor  map"
   !----------------------------------------------------------------------------
   !--- create new map datatype to hold result of matrix-matrix multiply ---
   call map_dup(map_orig,map_new)
   map_new%title  = trim(title)
   call map_matMatMult(map_orig,map_new,map_smooth) ! mult(A,B,S): B=S*A
   call mapsort_sort(map_new)
   call map_check(map_new)
   call map_write(map_new, trim(file_new))
else if (.not.(step1.or.step2)) then
   ! If we want to generate the rof->ocn nearest-neighbor or ocn->ocn smooth
   ! map, we won't necessarily have the rof->ocn nnsm map to sort and check.
   ! It only makes sense to call this step if we are skipping all three
   ! steps.
   call map_read(map_new, trim(file_new))
   call mapsort_sort(map_new)
   call map_check(map_new)
endif
   call shr_timer_stop (t07)
   call shr_timer_print(t07)

   call shr_timer_print_all()
   stop

END PROGRAM main

!===============================================================================
