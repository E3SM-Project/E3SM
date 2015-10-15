#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module repro_sum_mod
#define SPMD
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute reproducible global sums of a set of arrays across an MPI 
! subcommunicator
!
! Methods: 
! Compute using both a scalable, reproducible algorithm and a
! scalable, nonreproducible algorithm:
! * Reproducible (scalable): 
!    Convert to fixed point (integer vector representation) to enable
!    reproducibility when using MPI_Allreduce
! * Alternative reproducible (scalable):
!    Use parallel double-double algorithm due to Helen He and 
!    Chris Ding, based on David Bailey's/Don Knuth's DDPDD algorithm
! * Nonreproducible (scalable): 
!    Floating point and MPI_Allreduce based.
! Compare these and report relative difference (if absolute difference
! less than sum) or absolute difference back to calling routine.
!
! Author: P. Worley (based on suggestions from J. White for fixed 
!                    point algorithm and on He/Ding paper for ddpdd 
!                    algorithm)
! 
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!- use statements ------------------------------------------------------
!-----------------------------------------------------------------------
#if ( defined noI8 )
   ! Workaround for when shr_kind_i8 is not supported.
   use shr_kind_mod,  only: r8 => shr_kind_r8, i8 => shr_kind_i4 
#else
   use shr_kind_mod,  only: r8 => shr_kind_r8, i8 => shr_kind_i8 
#endif
   use shr_sys_mod,   only: shr_sys_abort  
   use perf_mod  ! _EXTERNAL
   use parallel_mod, only : abortmp
!-----------------------------------------------------------------------
!- module boilerplate --------------------------------------------------
!-----------------------------------------------------------------------
   implicit none
   private
!-----------------------------------------------------------------------
!- include statements --------------------------------------------------
!-----------------------------------------------------------------------
#if ( defined SPMD )
#include <mpif.h>  
#endif

   save

!-----------------------------------------------------------------------
! Public interfaces ----------------------------------------------------
!-----------------------------------------------------------------------
   public :: &
      repro_sum_defaultopts, &! get default runtime options
      repro_sum_setopts,     &! set runtime options
      repro_sum,             &! calculate distributed sum
      repro_sum_tol_exceeded  ! utility function to check relative 
                              !  differences against the tolerance

!-----------------------------------------------------------------------
! Public data ----------------------------------------------------------
!-----------------------------------------------------------------------
   logical, parameter  :: def_repro_sum_recompute = .false.
   logical, public     :: repro_sum_recompute = def_repro_sum_recompute

   real(r8), parameter :: def_rel_diff_max = -1.0_r8
   real(r8), public    :: repro_sum_rel_diff_max = def_rel_diff_max

!-----------------------------------------------------------------------
! Private interfaces ----------------------------------------------------
!-----------------------------------------------------------------------
   private :: &
      ddpdd             ! double-double sum routine

!-----------------------------------------------------------------------
! Private data ----------------------------------------------------------
!-----------------------------------------------------------------------

   integer, parameter :: def_repro_logunit = 6               ! default

   !----------------------------------------------------------------------------
   ! repro_sum_mod options
   !----------------------------------------------------------------------------
   logical, parameter :: def_use_ddpdd = .false.             ! default
   logical            :: repro_sum_use_ddpdd = def_use_ddpdd

   CONTAINS

!
!========================================================================
!
   subroutine repro_sum_defaultopts(repro_sum_use_ddpdd_out,    &
                                    repro_sum_rel_diff_max_out, &
                                    repro_sum_recompute_out     )
!----------------------------------------------------------------------- 
! Purpose: Return default runtime options
! Author: P. Worley
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
      ! Use DDPDD algorithm instead of fixed precision algorithm
      logical, intent(out), optional :: repro_sum_use_ddpdd_out
      ! maximum permissible difference between reproducible and 
      ! nonreproducible sums
      real(r8), intent(out), optional :: repro_sum_rel_diff_max_out
      ! recompute using different algorithm when difference between 
      ! reproducible and nonreproducible sums is too great 
      logical, intent(out), optional :: repro_sum_recompute_out
!-----------------------------------------------------------------------
      if ( present(repro_sum_use_ddpdd_out) ) then
         repro_sum_use_ddpdd_out = def_use_ddpdd
      endif
      if ( present(repro_sum_rel_diff_max_out) ) then
         repro_sum_rel_diff_max_out = def_rel_diff_max
      endif
      if ( present(repro_sum_recompute_out) ) then
         repro_sum_recompute_out = def_repro_sum_recompute
      endif
   end subroutine repro_sum_defaultopts
!
!========================================================================
!
   subroutine repro_sum_setopts(repro_sum_use_ddpdd_in,    &
                                repro_sum_rel_diff_max_in, &
                                repro_sum_recompute_in,    &
                                repro_sum_master,          &
                                repro_sum_logunit          )
                                
!----------------------------------------------------------------------- 
! Purpose: Set runtime options
! Author: P. Worley
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
      ! Use DDPDD algorithm instead of fixed precision algorithm
      logical, intent(in), optional :: repro_sum_use_ddpdd_in
      ! maximum permissible difference between reproducible and 
      ! nonreproducible sums
      real(r8), intent(in), optional :: repro_sum_rel_diff_max_in
      ! recompute using different algorithm when difference between 
      ! reproducible and nonreproducible sums is too great 
      logical, intent(in), optional  :: repro_sum_recompute_in
      ! flag indicating whether this process should output
      ! log messages
      logical, intent(in), optional  :: repro_sum_master
      ! unit number for log messages
      integer, intent(in), optional  :: repro_sum_logunit
!---------------------------Local Workspace-----------------------------
      integer logunit            ! unit number for log messages
      logical master             ! local master?
!-----------------------------------------------------------------------
      if ( present(repro_sum_master) ) then
         master = repro_sum_master
      else
         master = .false.
      endif

      if ( present(repro_sum_logunit) ) then
         logunit = repro_sum_logunit
      else
         logunit = def_repro_logunit
      endif

      if ( present(repro_sum_use_ddpdd_in) ) then
         repro_sum_use_ddpdd = repro_sum_use_ddpdd_in
      endif
      if ( present(repro_sum_rel_diff_max_in) ) then
         repro_sum_rel_diff_max = repro_sum_rel_diff_max_in
      endif
      if ( present(repro_sum_recompute_in) ) then
         repro_sum_recompute = repro_sum_recompute_in
      endif
      if (master) then
         if ( repro_sum_use_ddpdd ) then
            write(logunit,*) 'REPRO_SUM_SETOPTS: ',&
              'Using double-double-based (scalable) usually reproducible ', &
              'distributed sum algorithm'
         else
            write(logunit,*) 'REPRO_SUM_SETOPTS: ',&
              'Using fixed-point-based (scalable) reproducible ', &
              'distributed sum algorithm'
         endif

         if (repro_sum_rel_diff_max >= 0._r8) then
            write(logunit,*) '                    ',&
              'with a maximum relative error tolerance of ', &
              repro_sum_rel_diff_max
            if (repro_sum_recompute) then
               write(logunit,*) '                    ',&
                 'If tolerance exceeded, sum is recomputed using ', &
                 'a serial algorithm.'
            else
               write(logunit,*) '                    ',&
                 'If tolerance exceeded, fixed-precision is sum used ', &
                 'but a warning is output.'
            endif
         else
            write(logunit,*) '                    ',&
              'and not comparing with floating point algorithms.'
         endif

      endif
   end subroutine repro_sum_setopts

!
!========================================================================
!

   subroutine repro_sum (arr, arr_gsum, nsummands, dsummands, nflds, &
                         ddpdd_sum,                                  &
                         arr_max_levels, gbl_count, gbl_count_out,   &
                         arr_gbl_max, repro_sum_validate,            &
                         repro_sum_stats, rel_diff, commid           )
!----------------------------------------------------------------------
! 
! Purpose: 
! Compute the global sum of each field in "arr" using the indicated 
! communicator with a reproducible yet scalable implementation based 
! on a fixed point algorithm. An alternative is to use an "almost
! always" reproducible floating point algorithm, as described below.
!
! The accuracy of the fixed point algorithm is controlled by the
! number of "levels" of integer expansion. The algorithm will calculate
! the number of levels that is required for the sum to be essentially
! exact. The optional parameter arr_max_levels can be used to override
! the calculated value. 
!
! The algorithm also requires an upper bound on
! the maximum summand (in absolute value) for each field, and will 
! calculate this internally. However, if the optional parameters
! arr_max_levels and arr_gbl_max are both set, then the algorithm will
! use the values in arr_gbl_max for the upper bounds instead. If these
! are not upper bounds, or if the upper bounds are not tight enough
! to achieve the requisite accuracy, the algorithm will repeat the
! computation with appropriate upper bounds, returning the calculated
! upper bounds in arr_gbl_max. If only arr_gbl_max is present, then 
! just the calculated upper bounds are returned.
!
! Finally, the algorithm requires the total number of summands in the
! distributed sum. This will be calculated internally, using an MPI
! collective, but the value in the optional argument "gbl_count" will 
! be used instead if it is present and if >= nsummands. The accuracy 
! of this value is not checked. However, if set to < nsummands the 
! value will instead be calculated. If the optional parameter 
! gbl_count_out is present, then the value used 
! (gbl_count if >= nsummands; calculated otherwise) will be 
! returned.
!
! If requested (by setting repro_sum_rel_diff_max >= 0.0 and passing in
! the optional rel_diff parameter), results are compared with a 
! nonreproducible floating point algorithm. 
!
! The alternative algorithm is a minor modification of a parallel
! implementation of David Bailey's routine DDPDD by Helen He
! and Chris Ding. Bailey uses the Knuth trick to implement quadruple
! precision summation of double precision values with 10 double 
! precision operations. The advantage of this algorithm is that 
! it requires a single MPI_Allreduce. The disadvantage is that it 
! is not guaranteed to be reproducible (though it is in all 
! practical applications). This alternative is used when the optional 
! parameter ddpdd_sum is set to .true. It is also used if the 
! fixed precision algorithm radix assumption does not hold. 
! 
!----------------------------------------------------------------------
!
! Arguments
!
      integer,  intent(in) :: nsummands  ! number of local summands
      integer,  intent(in) :: dsummands  ! declared first dimension
      integer,  intent(in) :: nflds      ! number of fields
      real(r8), intent(in) :: arr(dsummands,nflds) 
                                         ! input array

      real(r8), intent(out):: arr_gsum(nflds)      
                                         ! global means

      logical,  intent(in),    optional :: ddpdd_sum
                                         ! use ddpdd algorithm instead
                                         ! of fixed precision algorithm

      integer,  intent(in),    optional :: arr_max_levels(nflds)
                                         ! maximum number of levels of
                                         ! integer expansion to use
      integer,  intent(in),    optional :: commid     
                                         ! MPI communicator

      integer,  intent(in),    optional :: gbl_count  
                                         ! total number of summands

      integer,  intent(out),   optional :: gbl_count_out
                                         ! calculated total number 
                                         ! of summands

      integer,  intent(inout), optional :: repro_sum_stats(5)
                                   ! increment running totals for
                                   !  (1) one-reduction repro_sum
                                   !  (2) two-reduction repro_sum
                                   !  (3) both types in one call
                                   !  (4) nonrepro_sum
                                   !  (5) global count reduction

      real(r8), intent(inout), optional :: arr_gbl_max(nflds)
                                         ! upper bound on max(abs(arr))

      real(r8), intent(out),   optional :: rel_diff(2,nflds)
                                         ! relative and absolute
                                         !  differences between fixed 
                                         !  and floating point sums

      logical,  intent(in),    optional :: repro_sum_validate
         ! flag enabling/disabling testing that gmax and  max_levels are 
         ! accurate/sufficient. Default is enabled.

!
! Local workspace
!
      logical :: use_ddpdd_sum           ! flag indicating whether to
                                         !  use repro_sum_ddpdd or not
      logical :: recompute               ! flag indicating need to 
                                         !  determine gmax/gmin before 
                                         !  computing sum
      logical :: validate                ! flag indicating need to 
                                         !  verify gmax and max_levels 
                                         !  are accurate/sufficient
      integer :: mpi_comm                ! MPI subcommunicator
      integer :: tot_summands            ! total number of summands
      integer :: ierr                    ! MPI error return
      integer :: ifld, isum              ! loop variables

      integer :: arr_lmax_exp(nflds)     ! local exponent maxima
      integer :: arr_lmin_exp(nflds)     ! local exponent minima
      integer :: arr_lextremes(nflds,2)  ! local exponent extrema
      integer :: arr_gextremes(nflds,2)  ! global exponent extrema

      integer :: arr_gmax_exp(nflds)     ! global exponents maxima 
      integer :: arr_gmin_exp(nflds)     ! global exponents minima
      integer :: arr_max_shift           ! maximum safe exponent for 
                                         !  value < 1 (so that sum does
                                         !  not overflow)
      integer :: max_levels(nflds)       ! maximum number of levels of 
                                         !  integer expansion to use
      integer :: max_level               ! maximum value in max_levels
      integer :: gbl_cnt_red             ! global count reduction? (0/1)
      integer :: repro_sum_fast          ! 1 reduction repro_sum? (0/1)
      integer :: repro_sum_slow          ! 2 reduction repro_sum? (0/1)
      integer :: repro_sum_both          ! both fast and slow? (0/1)
      integer :: nonrepro_sum            ! nonrepro_sum? (0/1)

      real(r8) :: xtot_summands          ! r8 version of tot_summands

      real(r8) :: arr_lsum(nflds)        ! local sums
      real(r8) :: arr_gsum_fast(nflds)   ! global sum calculated using 
                                         !  fast, nonreproducible, 
                                         !  floating point alg.
      real(r8) :: abs_diff               ! absolute difference between 
                                         !  fixed and floating point 
                                         !  sums
!
!-----------------------------------------------------------------------
!
! check whether should use repro_sum_ddpdd algorithm
      use_ddpdd_sum = repro_sum_use_ddpdd
      if ( present(ddpdd_sum) ) then
         use_ddpdd_sum = ddpdd_sum
      endif

! check whether intrinsic-based algorithm will work on this system
      if (.not. use_ddpdd_sum) then
         if (radix(0_r8) .ne. radix(0_i8)) then
            use_ddpdd_sum = .true.
!           call shr_sys_abort('repro_sum failed: mathematical models for '// &
!                              'integer and floating point use different bases.' )
         call abortmp("use_ddpdd_sum cannot be TRUE in homme-standalone")
         endif
      endif

! initialize local statistics variables
      gbl_cnt_red = 0
      repro_sum_fast = 0
      repro_sum_slow = 0
      repro_sum_both = 0
      nonrepro_sum = 0

#if ( defined SPMD )
! set MPI communicator
      if ( present(commid) ) then
         mpi_comm = commid
      else
         mpi_comm = MPI_COMM_WORLD
      endif
      call t_barrierf('sync_repro_sum',mpi_comm)
#else
      mpi_comm = 0
#endif

      if ( use_ddpdd_sum ) then

         call t_startf('repro_sum_ddpdd')

         call repro_sum_ddpdd(arr, arr_gsum, nsummands, dsummands, &
                              nflds, mpi_comm)
         repro_sum_fast = 1

         call t_stopf('repro_sum_ddpdd')

      else

         call t_startf('repro_sum_int')

! determine global number of summands
#if ( defined SPMD )
         if ( present(gbl_count) ) then
            if (gbl_count < nsummands) then
               call mpi_allreduce (nsummands, tot_summands, 1, &
                                   MPI_INTEGER, MPI_SUM, mpi_comm, ierr)
               gbl_cnt_red = 1
            else
               tot_summands = gbl_count
            endif
         else
            call mpi_allreduce (nsummands, tot_summands, 1, &
                                MPI_INTEGER, MPI_SUM, mpi_comm, ierr)
            gbl_cnt_red = 1
         endif
#else
         tot_summands = nsummands
#endif
         if ( present( gbl_count_out) ) then
            gbl_count_out = tot_summands
         endif

! determine maximum shift
!  (shift needs to be small enough that summation does not exceed 
!   maximum number of digits in i8)
         xtot_summands = tot_summands
         arr_max_shift = digits(0_i8) - (exponent(xtot_summands) + 1)
         if (arr_max_shift < 2) then
            call shr_sys_abort('repro_sum failed: number of summands too '// &
                               'large for fixed precision algorithm' )
         endif

! see if have sufficient information to not require max/min allreduce
         recompute = .true.
         validate = .false.
         if ( present(arr_gbl_max) .and. present(arr_max_levels) ) then
            recompute = .false.

            max_level = 0
            do ifld=1,nflds
               if ((arr_gbl_max(ifld) .ge. 0.0_r8) .and. &
                   (arr_max_levels(ifld) > 0)) then

                  arr_gmax_exp(ifld)  = exponent(arr_gbl_max(ifld))
                  if (max_level < arr_max_levels(ifld)) &
                     max_level = arr_max_levels(ifld)

               else
                  recompute = .true.
               endif
            enddo

            if (.not. recompute) then
! calculate sum
               if (present(repro_sum_validate)) then
                  validate = repro_sum_validate
               else
                  validate = .true.
               endif
               call repro_sum_int(arr, arr_gsum, nsummands, dsummands, &
                                  nflds, arr_max_shift, arr_gmax_exp, &
                                  arr_max_levels, max_level, validate, &
                                  recompute, mpi_comm)
! record statistics
               repro_sum_fast = 1
               if (recompute) then
                  repro_sum_both = 1
               endif
            endif
         endif

! do not have sufficient information; calculate global max/min and 
! use to compute required number of levels
         if (recompute) then

! record statistic
            repro_sum_slow = 1

! determine maximum and minimum (non-zero) summand values
            arr_lmax_exp(:) = MINEXPONENT(1._r8)
            arr_lmin_exp(:) = MAXEXPONENT(1._r8)
!$omp parallel do      &
!$omp default(shared)  &
!$omp private(ifld, isum)
            do ifld=1,nflds
               do isum=1,nsummands
                  if (arr(isum,ifld) .ne. 0.0) then
                     arr_lmax_exp(ifld) = &
                        max(exponent(arr(isum,ifld)),arr_lmax_exp(ifld))
                     arr_lmin_exp(ifld) = &
                        min(exponent(arr(isum,ifld)),arr_lmin_exp(ifld))
                  endif
               end do
            end do

#if ( defined SPMD )
            arr_lextremes(:,1) = -arr_lmax_exp(:)
            arr_lextremes(:,2) = arr_lmin_exp(:)
            call mpi_allreduce (arr_lextremes, arr_gextremes, 2*nflds, &
                                MPI_INTEGER, MPI_MIN, mpi_comm, ierr)
            arr_gmax_exp(:) = -arr_gextremes(:,1)
            arr_gmin_exp(:) =  arr_gextremes(:,2)
#else 
            arr_gmax_exp(:) = arr_lmax_exp(:)
            arr_gmin_exp(:) = arr_lmin_exp(:)
#endif
! if a field is identically zero, arr_gmin_exp still equals MAXEXPONENT
!   and arr_gmax_exp still equals MINEXPONENT. In this case, set 
!   arr_gmin_exp = arr_gmax_exp = MINEXPONENT
            do ifld=1,nflds
               arr_gmin_exp(ifld) = min(arr_gmax_exp(ifld),arr_gmin_exp(ifld))
            enddo

! determine maximum number of levels required for each field
! ((digits(0_i8) + (arr_gmax_exp(ifld)-arr_gmin_exp(ifld))) / arr_max_shift)
! + 1 because first truncation probably does not involve a maximal shift
! + 1 to guarantee that the integer division rounds up (not down)
            max_level = 0
            do ifld=1,nflds
               max_levels(ifld) = 2 + &
                ((digits(0_i8) + (arr_gmax_exp(ifld)-arr_gmin_exp(ifld))) &
                / arr_max_shift)
               if ( present(arr_max_levels) .and. (.not. validate) ) then
! if validate true, then computation with arr_max_levels failed
!  previously
                  if ( arr_max_levels(ifld) > 0 ) then
                     max_levels(ifld) = &
                        min(arr_max_levels(ifld),max_levels(ifld))
                  endif
               endif
               if (max_level < max_levels(ifld)) &
                  max_level = max_levels(ifld)
            enddo

! calculate sum
            validate = .false.
            call repro_sum_int(arr, arr_gsum, nsummands, dsummands, nflds, &
                               arr_max_shift, arr_gmax_exp, max_levels, &
                               max_level, validate, recompute, mpi_comm)

! return upper bounds on observed maxima 
            if ( present(arr_gbl_max) ) then
               do ifld=1,nflds
                  arr_gbl_max(ifld) = scale(1.0_r8,arr_gmax_exp(ifld))
               enddo
            endif

         endif

         call t_stopf('repro_sum_int')

      endif

! compare fixed and floating point results
      if ( present(rel_diff) ) then
         if (repro_sum_rel_diff_max >= 0.0_r8) then

#if ( defined SPMD )
            call t_barrierf('sync_nonrepro_sum',mpi_comm)
#endif
            call t_startf('nonrepro_sum')
! record statistic
            nonrepro_sum = 1
! compute nonreproducible sum
            arr_lsum(:) = 0._r8
!$omp parallel do      &
!$omp default(shared)  &
!$omp private(ifld, isum)
            do ifld=1,nflds
               do isum=1,nsummands
                  arr_lsum(ifld) = arr(isum,ifld) + arr_lsum(ifld)
               end do
            end do

#if ( defined SPMD )
            call mpi_allreduce (arr_lsum, arr_gsum_fast, nflds, &
                                MPI_REAL8, MPI_SUM, mpi_comm, ierr)
#else 
            arr_gsum_fast(:) = arr_lsum(:)
#endif
            call t_stopf('nonrepro_sum')

! determine differences
!$omp parallel do      &
!$omp default(shared)  &
!$omp private(ifld, abs_diff)
            do ifld=1,nflds
               abs_diff = abs(arr_gsum_fast(ifld)-arr_gsum(ifld))
               if (abs(arr_gsum(ifld)) > abs_diff) then
                  rel_diff(1,ifld) = abs_diff/abs(arr_gsum(ifld))
               else
                  rel_diff(1,ifld) = abs_diff
               endif
               rel_diff(2,ifld) = abs_diff
            enddo
         else
            rel_diff(:,:) = 0.0_r8
         endif
      endif

! return statistics
      if ( present(repro_sum_stats) ) then
         repro_sum_stats(1) = repro_sum_stats(1) + repro_sum_fast
         repro_sum_stats(2) = repro_sum_stats(2) + repro_sum_slow
         repro_sum_stats(3) = repro_sum_stats(3) + repro_sum_both
         repro_sum_stats(4) = repro_sum_stats(4) + nonrepro_sum
         repro_sum_stats(5) = repro_sum_stats(5) + gbl_cnt_red
      endif

      return

   end subroutine repro_sum

!
!========================================================================
!

   subroutine repro_sum_int (arr, arr_gsum, nsummands, dsummands, nflds, &
                             arr_max_shift, arr_gmax_exp, max_levels, &
                             max_level, validate, recompute, mpi_comm)
!----------------------------------------------------------------------
! 
! Purpose: 
! Compute the global sum of each field in "arr" using the indicated 
! communicator with a reproducible yet scalable implementation based 
! on a fixed point algorithm. The accuracy of the fixed point algorithm
! is controlled by the number of "levels" of integer expansion, the
! maximum value of which is specified by max_level. 
! 
!----------------------------------------------------------------------
!
! Arguments
!
      integer,  intent(in) :: nsummands     ! number of local summands
      integer,  intent(in) :: dsummands     ! declared first dimension
      integer,  intent(in) :: nflds         ! number of fields
      integer,  intent(in) :: arr_max_shift ! maximum safe exponent for 
                                            !  value < 1 (so that sum 
                                            !  does not overflow)
      integer,  intent(in) :: arr_gmax_exp(nflds)  
                                            ! exponents of global maxima
      integer,  intent(in) :: max_levels(nflds) 
                                            ! maximum number of levels 
                                            !  of integer expansion
      integer,  intent(in) :: max_level     ! maximum value in 
                                            !  max_levels
      integer,  intent(in) :: mpi_comm      ! MPI subcommunicator

      real(r8), intent(in) :: arr(dsummands,nflds) 
                                             ! input array

      logical,  intent(in):: validate
         ! flag indicating that accuracy of solution generated from
         ! arr_gmax_exp and max_levels should be tested

      logical,  intent(out):: recompute            
         ! flag indicating that either the upper bounds are inaccurate,
         !  or max_levels and arr_gmax_exp do not generate accurate 
         !  enough sums

      real(r8), intent(out):: arr_gsum(nflds)      ! global means
!
! Local workspace
!
      integer, parameter  :: max_jlevel = &
                                1 + (digits(0_i8)/digits(0.0_r8))

      integer(i8) :: i8_arr_lsum_level((max_level+2)*nflds) 
                                   ! integer vector representing local 
                                   !  sum
      integer(i8) :: i8_arr_level  ! integer part of summand for current
                                   !  expansion level
      integer(i8) :: i8_arr_gsum_level((max_level+2)*nflds) 
                                   ! integer vector representing global
                                   !  sum
      integer(i8) :: IX_8          ! integer representation of current 
                                   !  jlevels of X_8 ('part' of 
                                   !  i8_arr_gsum_level)

      integer :: max_error(nflds)  ! accurate of upper bound on data?
      integer :: not_exact(nflds)  ! max_levels sufficient to 
                                   !  capture all digits?
      integer :: ifld, isum        ! loop variables
      integer :: arr_exp           ! exponent of summand
      integer :: arr_shift         ! exponent used to generate integer 
                                   !  for current expansion level
      integer :: ilevel            ! current integer expansion level
      integer :: offset(nflds)     ! beginning location in 
                                   !  i8_arr_{g,l}sum_level for integer 
                                   !  expansion of current ifld
      integer :: voffset           ! modification to offset used to 
                                   !  include validation metrics 
      integer :: ioffset           ! offset(ifld)
      integer :: jlevel            ! number of floating point 'pieces' 
                                   !  extracted from a given i8 integer
      integer :: ierr              ! MPI error return
      integer :: LX(max_jlevel)    ! exponent of X_8 (see below)
      integer :: veclth            ! total length of i8_arr_lsum_level
      integer :: sum_digits        ! lower bound on number of significant
                                   !  in integer expansion of sum
      integer :: curr_exp          ! exponent of partial sum during
                                   !  reconstruction from integer vector
      integer :: corr_exp          ! exponent of current summand in
                                   !  reconstruction from integer vector

      real(r8) :: arr_frac         ! fraction of summand
      real(r8) :: arr_remainder    ! part of summand remaining after 
                                   !  current level of integer expansion
      real(r8) :: X_8(max_jlevel)  ! r8 vector representation of current 
                                   !  i8_arr_gsum_level
      real(r8) :: RX_8             ! r8 representation of difference 
                                   !  between current i8_arr_gsum_level
                                   !  and current jlevels of X_8 
                                   !  (== IX_8). Also used in final
                                   !  scaling step

      logical :: first             ! flag used to indicate that just
                                   !  beginning reconstruction of sum
                                   !  from integer vector

!
!-----------------------------------------------------------------------
! if validating upper bounds, reserve space for validation metrics
      if (validate) then
        voffset = 2
      else
        voffset = 0
      endif

! compute offsets for each field
      offset(1) = voffset
      do ifld=2,nflds
         offset(ifld) = offset(ifld-1) &
                        + (max_levels(ifld-1) + voffset)
      enddo
      veclth = offset(nflds) + max_levels(nflds)

! convert local summands to vector of integers and sum
! (Using scale instead of set_exponent because arr_remainder may not be
!  "normal" after level 1 calculation)
      i8_arr_lsum_level(:) = 0_i8
      max_error(:) = 0
      not_exact(:) = 0

!$omp parallel do      &
!$omp default(shared)  &
!$omp private(ifld, ioffset, isum, arr_frac, arr_exp, arr_shift, &
!$omp         ilevel, i8_arr_level, arr_remainder)
      do ifld=1,nflds
         ioffset = offset(ifld)

         do isum=1,nsummands
            arr_remainder = 0.0_r8

            if (arr(isum,ifld) .ne. 0.0_r8) then
               arr_exp   = exponent(arr(isum,ifld))
               arr_frac  = fraction(arr(isum,ifld))

! test that global maximum upper bound is an upper bound
               if (arr_exp > arr_gmax_exp(ifld)) then
                  max_error(ifld) = 1
                  exit
               endif

! calculate first shift
               arr_shift = arr_max_shift - (arr_gmax_exp(ifld)-arr_exp)

! determine first (probably) nonzero level (assuming initial fraction is 
!  'normal' - algorithm still works if this is not true)
!  NOTE: this is critical; scale will set to zero if min exponent is too small.
               if (arr_shift < 1) then
                  ilevel = (1 + (arr_gmax_exp(ifld)-arr_exp))/arr_max_shift
                  arr_shift = ilevel*arr_max_shift - (arr_gmax_exp(ifld)-arr_exp)

                  do while (arr_shift < 1)
                     arr_shift = arr_shift + arr_max_shift
                     ilevel = ilevel + 1
                  enddo
               else
                  ilevel = 1
               endif

               if (ilevel .le. max_levels(ifld)) then
! apply first shift/truncate, add it to the relevant running 
! sum, and calculate the remainder.
                  arr_remainder = scale(arr_frac,arr_shift)
                  i8_arr_level = int(arr_remainder,i8)
                  i8_arr_lsum_level(ioffset+ilevel) = &
                     i8_arr_lsum_level(ioffset+ilevel) + i8_arr_level
                  arr_remainder = arr_remainder - i8_arr_level

! while the remainder is non-zero, continue to shift, truncate, 
! sum, and calculate new remainder
                  do while ((arr_remainder .ne. 0.0_r8) &
                     .and. (ilevel < max_levels(ifld)))
                     ilevel = ilevel + 1
                     arr_remainder = scale(arr_remainder,arr_max_shift)
                     i8_arr_level = int(arr_remainder,i8)
                     i8_arr_lsum_level(ioffset+ilevel) = &
                        i8_arr_lsum_level(ioffset+ilevel) + i8_arr_level
                     arr_remainder = arr_remainder - i8_arr_level
                  enddo

               endif
            endif

            if (arr_remainder .ne. 0.0_r8) then
               not_exact(ifld) = 1
            endif

         enddo

! record if upper bound was inaccurate or if level expansion stopped
! before full accuracy was achieved
         if (validate) then
            i8_arr_lsum_level(ioffset-voffset+1) = max_error(ifld)
            i8_arr_lsum_level(ioffset-voffset+2) = not_exact(ifld)
         endif

      enddo

#if ( defined SPMD )
! sum integer vector element-wise
#if ( defined noI8 )
     ! Workaround for when shr_kind_i8 is not supported.
      call mpi_allreduce (i8_arr_lsum_level, i8_arr_gsum_level, & 
                          veclth, MPI_INTEGER, MPI_SUM, mpi_comm, ierr)
#else
      call mpi_allreduce (i8_arr_lsum_level, i8_arr_gsum_level, & 
                          veclth, MPI_INTEGER8, MPI_SUM, mpi_comm, ierr)
#endif
#else
      i8_arr_gsum_level(:) = i8_arr_lsum_level(:)
#endif

! Construct global sum from integer vector representation:
!  1) arr_max_shift is the shift applied to fraction(arr_gmax) .
!   When shifting back, need to "add back in" true arr_gmax exponent. This was
!   removed implicitly by working only with the fraction .
!  2) want to add levels into sum in reverse order (smallest to largest)
!  3) assignment to X_8 will usually lose accuracy since maximum integer 
!   size is greater than the max number of 'digits' in r8 value (if tot_summands 
!   correction not very large). Calculate remainder and add in first (since 
!   smaller). One correction is sufficient for r8 (53 digits) and i8 (63 digits). 
!   For r4 (24 digits) may need to correct twice. Code is written in a general
!   fashion, to work no matter how many corrections are necessary (assuming
!   max_jlevel parameter calculation is correct).

      recompute = .false.
      do ifld=1,nflds
         arr_gsum(ifld) = 0.0_r8
         ioffset = offset(ifld)

! if validate is .true., test whether the summand upper bound
!  was exceeded on any of the processes
         if (validate) then
            if (i8_arr_gsum_level(ioffset-voffset+1) .ne. 0_i8) then
               recompute = .true.
            endif
         endif

         if (.not. recompute) then

! start with maximum shift, and work up to larger values
            arr_shift = arr_gmax_exp(ifld) &
                        - max_levels(ifld)*arr_max_shift
            curr_exp = 0
            first = .true.
            do ilevel=max_levels(ifld),1,-1

               if (i8_arr_gsum_level(ioffset+ilevel) .ne. 0_i8) then
                  jlevel = 1

! r8 representation of higher order bits in integer 
                  X_8(jlevel) = i8_arr_gsum_level(ioffset+ilevel)
                  LX(jlevel)  = exponent(X_8(jlevel))

! calculate remainder
                  IX_8 = int(X_8(jlevel),i8)
                  RX_8 = (i8_arr_gsum_level(ioffset+ilevel) - IX_8)

! repeat using remainder
                  do while (RX_8 .ne. 0.0_r8)
                     jlevel = jlevel + 1
                     X_8(jlevel) = RX_8
                     LX(jlevel) = exponent(RX_8)
                     IX_8 = IX_8 + int(RX_8,i8)
                     RX_8 = (i8_arr_gsum_level(ioffset+ilevel) - IX_8)
                  enddo

! add in contributions, smaller to larger, rescaling for each
! addition to guarantee that exponent of working summand is always 
! larger than minexponent
                  do while (jlevel > 0)
                     if (first) then
                        curr_exp = LX(jlevel) + arr_shift
                        arr_gsum(ifld) = fraction(X_8(jlevel))
                        first = .false.
                     else
                        corr_exp = curr_exp - (LX(jlevel) + arr_shift)
                        arr_gsum(ifld) = fraction(X_8(jlevel)) &
                                       + scale(arr_gsum(ifld),corr_exp)
                        curr_exp = LX(jlevel) + arr_shift
                     endif
                     jlevel = jlevel - 1
                  enddo

               endif

               arr_shift = arr_shift + arr_max_shift
            enddo

! apply final exponent correction, scaling first if exponent is too small
! to apply directly
            corr_exp = curr_exp + exponent(arr_gsum(ifld))
            if (corr_exp .ge. MINEXPONENT(1._r8)) then 
               arr_gsum(ifld) = set_exponent(arr_gsum(ifld),corr_exp)
            else
               RX_8 = set_exponent(arr_gsum(ifld), &
                                   corr_exp-MINEXPONENT(1._r8))
               arr_gsum(ifld) = scale(RX_8,MINEXPONENT(1._r8))
            endif

! if validate is .true. and some precision lost, test whether 'too much' 
!  was lost, due to too loose an upper bound, too stringent a limit on number
!  of levels of expansion, cancellation, .... Calculated by comparing lower 
!  bound on number of sigificant digits with number of digits in 1.0_r8 .
            if (validate) then
               if (i8_arr_gsum_level(ioffset-voffset+2) .ne. 0_i8) then

! find first nonzero level and use exponent for this level, then assume all
! subsequent levels contribute arr_max_shift digits (ignoring probable
! increase in number of digits due to summing over tot_summands).
                  sum_digits = 0
                  do ilevel=1,max_levels(ifld)
                     if (i8_arr_gsum_level(ioffset+ilevel) .ne. 0_i8) then
                        X_8(1) = i8_arr_gsum_level(ioffset+ilevel)
                        LX(1)  = exponent(X_8(1))
                        if (sum_digits .eq. 0) then
                           sum_digits = LX(1)
                        else
                           sum_digits = sum_digits + arr_max_shift
                        endif
                     endif
                  enddo

                  if (sum_digits < digits(1.0_r8)) then
                     recompute = .true.
                  endif
               endif
            endif

         endif

      enddo

      return

   end subroutine repro_sum_int

!
!========================================================================
!

   logical function repro_sum_tol_exceeded (name, nflds, master, &
                                            logunit, rel_diff    )
!----------------------------------------------------------------------
! 
! Purpose: 
! Test whether distributed sum exceeds tolerance and print out a
! warning message.
! 
!----------------------------------------------------------------------
!
! Arguments
!
      character(len=*), intent(in) :: name    ! distributed sum identifier
      integer,  intent(in) :: nflds           ! number of fields
      logical,  intent(in) :: master          ! process that will write
                                              !  warning messages?
      integer,  intent(in) :: logunit         ! unit warning messages 
                                              !  written to
      real(r8), intent(in) :: rel_diff(2,nflds)
                                              ! relative and absolute
                                              !  differences between fixed 
                                              !  and floating point sums

!
! Local workspace
!
      integer  :: ifld                        ! field index
      integer  :: exceeds_limit               ! number of fields whose
                                              !  sum exceeds tolerance
      real(r8) :: max_rel_diff                ! maximum relative difference
      integer  :: max_rel_diff_idx            ! field index for max. rel. diff.
      real(r8) :: max_abs_diff                ! maximum absolute difference
      integer  :: max_abs_diff_idx            ! field index for max. abs. diff.
!
!-----------------------------------------------------------------------
!
      repro_sum_tol_exceeded = .false.
      if (repro_sum_rel_diff_max < 0.0_r8) return

      ! check that "fast" reproducible sum is accurate enough.
      exceeds_limit = 0
      max_rel_diff = 0.0_r8
      max_abs_diff = 0.0_r8
      do ifld=1,nflds
         if (rel_diff(1,ifld) > repro_sum_rel_diff_max) then
            exceeds_limit = exceeds_limit + 1
            if (rel_diff(1,ifld) > max_rel_diff) then
               max_rel_diff = rel_diff(1,ifld)
               max_rel_diff_idx = ifld
            endif
            if (rel_diff(2,ifld) > max_abs_diff) then
               max_abs_diff = rel_diff(2,ifld)
               max_abs_diff_idx = ifld
            endif
         endif
      enddo

      if (exceeds_limit > 0) then
         if (master) then
            write(logunit,*) trim(name), &
                            ': difference in fixed and floating point sums ', &
                            ' exceeds tolerance in ', exceeds_limit, &
                            ' fields.'
            write(logunit,*) '  Maximum relative diff: (rel)', &
                            rel_diff(1,max_rel_diff_idx), ' (abs) ', &
                            rel_diff(2,max_rel_diff_idx)
            write(logunit,*) '  Maximum absolute diff: (rel)', &
                            rel_diff(1,max_abs_diff_idx), ' (abs) ', &
                            rel_diff(2,max_abs_diff_idx)
         endif
         repro_sum_tol_exceeded = .true.
       endif

      return

   end function repro_sum_tol_exceeded

!
!========================================================================
!

   subroutine repro_sum_ddpdd (arr, arr_gsum, nsummands, dsummands,   &
                               nflds, mpi_comm                        )
!----------------------------------------------------------------------
! 
! Purpose: 
! Compute the global sum of each field in "arr" using the indicated 
! communicator with a reproducible yet scalable implementation based 
! on He and Ding's implementation of the double-double algorithm.
! 
!----------------------------------------------------------------------
!
! Arguments
!
      integer,  intent(in) :: nsummands  ! number of local summands
      integer,  intent(in) :: dsummands  ! declared first dimension
      integer,  intent(in) :: nflds      ! number of fields
      real(r8), intent(in) :: arr(dsummands,nflds) 
                                         ! input array
      integer,  intent(in) :: mpi_comm   ! MPI subcommunicator

      real(r8), intent(out):: arr_gsum(nflds)      
                                         ! global sums

!
! Local workspace
!
      integer :: old_cw                  ! for x86 processors, save
                                         !  current arithmetic mode
      integer :: ifld, isum              ! loop variables
      integer :: ierr                    ! MPI error return

      real(r8)    :: e, t1, t2           ! temporaries
      complex(r8) :: arr_lsum_dd(nflds)  ! local sums (in double-double
                                         !  format)
      complex(r8) :: arr_gsum_dd(nflds)  ! global sums (in double-double
                                         !  format)

      integer, save :: mpi_sumdd
      logical, save :: first_time = .true.

!
!-----------------------------------------------------------------------
!
!     call x86_fix_start (old_cw)

      if (first_time) then
#if ( defined SPMD )
         call mpi_op_create(ddpdd, .true., mpi_sumdd, ierr)
#endif
         first_time = .false.
      endif

      do ifld=1,nflds
         arr_lsum_dd(ifld) = (0.0_r8,0.0_r8)

         do isum=1,nsummands
                        
            ! Compute arr(isum,ifld) + arr_lsum_dd(ifld) using Knuth''s 
            ! trick.
            t1 = arr(isum,ifld) + real(arr_lsum_dd(ifld))
            e  = t1 - arr(isum,ifld)
            t2 = ((real(arr_lsum_dd(ifld)) - e) &
                  + (arr(isum,ifld) - (t1 - e))) &
                 + aimag(arr_lsum_dd(ifld))

            ! The result is t1 + t2, after normalization.
            arr_lsum_dd(ifld) = cmplx ( t1 + t2, t2 - ((t1 + t2) - t1), r8 )
         enddo

      enddo

#if ( defined SPMD )
      call mpi_allreduce (arr_lsum_dd, arr_gsum_dd, nflds, &
                          MPI_COMPLEX16, mpi_sumdd, mpi_comm, ierr)
      do ifld=1,nflds
         arr_gsum(ifld) = real(arr_gsum_dd(ifld))
      enddo
#else
      do ifld=1,nflds
         arr_gsum(ifld) = real(arr_lsum_dd(ifld))
      enddo
#endif

!     call x86_fix_end (old_cw)

      return

   end subroutine repro_sum_ddpdd
!
!-----------------------------------------------------------------------
!
   subroutine DDPDD (dda, ddb, len, itype)
!----------------------------------------------------------------------
! 
! Purpose: 
! Modification of original codes written by David H. Bailey    
! This subroutine computes ddb(i) = dda(i)+ddb(i)
! 
!----------------------------------------------------------------------
!
! Arguments
!
      integer, intent(in)        :: len       ! array length
      complex(r8), intent(in)    :: dda(len)  ! input
      complex(r8), intent(inout) :: ddb(len)  ! result
      integer, intent(in)        :: itype     ! unused
!
! Local workspace
!
      real(r8) e, t1, t2
      integer i
!
!-----------------------------------------------------------------------
!
      do i = 1, len
!   Compute dda + ddb using Knuth's trick.
         t1 = real(dda(i)) + real(ddb(i))
         e  = t1 - real(dda(i))
         t2 = ((real(ddb(i)) - e) + (real(dda(i)) - (t1 - e))) &
              + aimag(dda(i)) + aimag(ddb(i))

!   The result is t1 + t2, after normalization.
         ddb(i) = cmplx ( t1 + t2, t2 - ((t1 + t2) - t1), r8 )
      enddo

      return

   end subroutine DDPDD
!
!========================================================================
!
end module repro_sum_mod
