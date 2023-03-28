module shr_reprosum_mod
!------------------------------------------------------------------------
!
! Purpose:
! Compute reproducible global sums of a set of arrays across an MPI
! subcommunicator
!
! Methods:
! Compute using either or both a scalable, reproducible algorithm and a
! scalable, nonreproducible algorithm:
! * Reproducible (scalable):
!    Convert each floating point summand to an integer vector
!    representation, to enable reproducibility when using
!    MPI_Allreduce, then convert the resulting global sum back to a
!    floating point representation locally;
! * Alternative usually reproducible (scalable):
!    Use parallel double-double algorithm due to Helen He and
!    Chris Ding, based on David Bailey's/Don Knuth's DDPDD algorithm;
! * Nonreproducible (scalable):
!    Floating point and MPI_Allreduce based.
! If computing both reproducible and nonreproducible sums, compare
! these and report relative difference (if absolute difference
! less than sum) or absolute difference back to calling routine.
!
! Author: P. Worley (based on suggestions from J. White for integer
!                    vector algorithm and on He/Ding paper for DDPDD
!                    algorithm)
!
!------------------------------------------------------------------------

!------------------------------------------------------------------------
!- use statements -------------------------------------------------------
!------------------------------------------------------------------------
#if ( defined noI8 )
   ! Workaround for when shr_kind_i8 is not supported.
   use shr_kind_mod,  only: r8 => shr_kind_r8, i8 => shr_kind_i4
#else
   use shr_kind_mod,  only: r8 => shr_kind_r8, i8 => shr_kind_i8
#endif
   use shr_log_mod,   only: s_loglev  => shr_log_Level
   use shr_log_mod,   only: s_logunit => shr_log_Unit
   use shr_sys_mod,   only: shr_sys_abort
   use shr_infnan_mod,only: shr_infnan_inf_type, assignment(=), &
                            shr_infnan_posinf, shr_infnan_neginf, &
                            shr_infnan_nan, &
                            shr_infnan_isnan, shr_infnan_isinf, &
                            shr_infnan_isposinf, shr_infnan_isneginf
   use perf_mod

!------------------------------------------------------------------------
!- module boilerplate ---------------------------------------------------
!------------------------------------------------------------------------
   implicit none
   private

!------------------------------------------------------------------------
!- include statements ---------------------------------------------------
!------------------------------------------------------------------------
#include <mpif.h>

   save

!------------------------------------------------------------------------
! Public interfaces -----------------------------------------------------
!------------------------------------------------------------------------
   public :: &
      shr_reprosum_setopts,        &! set runtime options
      shr_reprosum_calc,           &! calculate distributed sum
      shr_reprosum_tolExceeded      ! utility function to check relative
                                    !  differences against the tolerance

!------------------------------------------------------------------------
! Public data -----------------------------------------------------------
!------------------------------------------------------------------------
   logical, public     :: shr_reprosum_recompute = .false.

   real(r8), public    :: shr_reprosum_reldiffmax = -1.0_r8

!------------------------------------------------------------------------
! Private interfaces ----------------------------------------------------
!------------------------------------------------------------------------
   private :: &
      ddpdd,           &! double-double sum routine
      split_indices     ! split indices among OMP threads

!------------------------------------------------------------------------
! Private data ----------------------------------------------------------
!------------------------------------------------------------------------

   !---------------------------------------------------------------------
   ! shr_reprosum_mod options
   !---------------------------------------------------------------------
   logical            :: repro_sum_use_ddpdd = .false.

   logical            :: repro_sum_allow_infnan = .false.

   CONTAINS

!
!========================================================================
!
   subroutine shr_reprosum_setopts(repro_sum_use_ddpdd_in,    &
                                   repro_sum_allow_infnan_in, &
                                   repro_sum_rel_diff_max_in, &
                                   repro_sum_recompute_in,    &
                                   repro_sum_master,          &
                                   repro_sum_logunit          )

!------------------------------------------------------------------------
! Purpose: Set runtime options
! Author: P. Worley
!------------------------------------------------------------------------
!------------------------------Arguments---------------------------------
      ! Use DDPDD algorithm instead of integer vector algorithm
      logical, intent(in), optional :: repro_sum_use_ddpdd_in
      ! Allow INF or NaN in summands
      logical, intent(in), optional :: repro_sum_allow_infnan_in
      ! maximum permissible difference between reproducible and
      ! nonreproducible sums
      real(r8), intent(in), optional :: repro_sum_rel_diff_max_in
      ! recompute using different algorithm when difference between
      ! reproducible and nonreproducible sums is too great
      logical, intent(in), optional  :: repro_sum_recompute_in
      ! flag indicating whether this MPI task should output
      ! log messages
      logical, intent(in), optional  :: repro_sum_master
      ! unit number for log messages
      integer, intent(in), optional  :: repro_sum_logunit
!---------------------------Local Workspace------------------------------
      integer logunit            ! unit number for log messages
      logical master             ! local master?
      logical,save :: firstcall = .true.  ! first call
!------------------------------------------------------------------------

      if ( present(repro_sum_master) ) then
         master = repro_sum_master
      else
         master = .false.
      endif

      if ( present(repro_sum_logunit) ) then
         logunit = repro_sum_logunit
      else
         logunit = s_logunit
      endif

      if (.not. firstcall) then
         write(logunit,*) 'shr_reprosum_setopts: ERROR can only be called once'
         call shr_sys_abort('shr_reprosum_setopts ERROR: multiple calls')
      endif
      firstcall = .false.

      if ( present(repro_sum_use_ddpdd_in) ) then
         repro_sum_use_ddpdd = repro_sum_use_ddpdd_in
      endif
      if ( present(repro_sum_allow_infnan_in) ) then
         repro_sum_allow_infnan = repro_sum_allow_infnan_in
      endif
      if ( present(repro_sum_rel_diff_max_in) ) then
         shr_reprosum_reldiffmax = repro_sum_rel_diff_max_in
      endif
      if ( present(repro_sum_recompute_in) ) then
         shr_reprosum_recompute = repro_sum_recompute_in
      endif
      if (master) then
         if ( repro_sum_use_ddpdd ) then
            write(logunit,*) 'SHR_REPROSUM_SETOPTS: ',&
              'Using double-double-based (scalable) usually reproducible ', &
              'distributed sum algorithm'
         else
            write(logunit,*) 'SHR_REPROSUM_SETOPTS: ',&
              'Using integer-vector-based (scalable) reproducible ', &
              'distributed sum algorithm'
         endif

         if ( repro_sum_allow_infnan ) then
            write(logunit,*) 'SHR_REPROSUM_SETOPTS: ',&
              'Will calculate sum when INF or NaN are included in summands'
         else
            write(logunit,*) 'SHR_REPROSUM_SETOPTS: ',&
              'Will abort if INF or NaN are included in summands'
         endif

         if (shr_reprosum_reldiffmax >= 0.0_r8) then
            write(logunit,*) '                    ',&
              'with a maximum relative error tolerance of ', &
              shr_reprosum_reldiffmax
            if (shr_reprosum_recompute) then
               write(logunit,*) '                    ',&
                 'If tolerance exceeded, sum is recomputed using ', &
                 'a serial algorithm.'
            else
               write(logunit,*) '                    ',&
                 'If tolerance exceeded, integer-vector-based sum is used ', &
                 'but a warning is output.'
            endif
         else
            write(logunit,*) '                    ',&
              'and not comparing with floating point algorithms.'
         endif

      endif
   end subroutine shr_reprosum_setopts

!
!========================================================================
!

   subroutine shr_reprosum_calc (arr, arr_gsum, nsummands, dsummands,     &
                                 nflds, allow_infnan, ddpdd_sum,          &
                                 arr_gbl_max, arr_gbl_max_out,            &
                                 arr_max_levels, arr_max_levels_out,      &
                                 gbl_max_nsummands, gbl_max_nsummands_out,&
                                 gbl_count, repro_sum_validate,           &
                                 repro_sum_stats, rel_diff, commid        )
!------------------------------------------------------------------------
!
! Purpose:
! Compute the global sum of each field in 'arr' using the indicated
! communicator with a reproducible yet scalable implementation based
! on first converting each floating point summand into an equivalent
! representation using a vector of integers, summing the integer
! vectors, then converting the resulting sum back to a floating point
! representation. An alternative is to use an 'almost always
! reproducible' floating point algorithm (DDPDD), as described below.
!
! Description of integer vector algorithm:
!-----------------------------------------
! The basic idea is to represent the mantissa of each floating point
! value as an integer, add these integers, and then convert back to a
! floating point value. For a real*8 value, there are enough digits in
! an integer*8 variable to not lose any information (in the
! mantissa). However, each of these integers would have a different
! implicit exponent if done in a naive way, and so the sum would not
! be accurate. Also, even with the same 'normalization', the sum might
! exceed the maximum value representable by an integer*8, causing
! an overflow. Instead, a vector of integers is generated, where a
! given element (or level using the terminology used in the code) of
! the vector is associated with a particular exponent. The mantissa
! for a given floating point value is then converted to some number of
! integer values, depending on the exponent of the floating point
! value, the normalization of its mantissa, the maximum number of
! summands, the number of participating MPI tasks and of OpenMP
! threads, and the exponents associated with the levels of the integer
! vector, and added into the appropriate levels of the integer
! vector. Each MPI task has its own integer vector representing the
! local sum. This is then summed across all participating MPI tasks
! using an MPI_Allreduce, and, lastly, converted back to a floating
! point value. Note that the same approach works for a vector of
! integer*4 variables, simply requiring more levels, both for the full
! summation vector and for each individual real*8 summand. This is a
! compile time option in the code, in support of systems for which the
! compiler or MPI library has issues when using integer*8. As
! implemented, this algorithm should work for any floating point and
! integer type as long as they share the same base. The code is
! written as if for real*8 and integer*8 variables, but the only
! dependence is on the types 'r8' and 'i8', which are defined in the
! code, currently with reference to the corresponding types in
! shr_kind_mod. This is how integer*4 support is implemented, by
! defining i8 to be shr_kind_i4 instead of shr_kind_i8.
!
! For this to work, each MPI task must have the same number of levels
! and same implicit exponent for each level. These levels must be
! sufficient to represent the smallest and largest nonzero individual
! summands (in absolute value) and the largest possible intermediate
! sum, including the final sum. Most of the complexity in the
! algorithm is in identifying the number of levels, the exponent
! associated with each level, and the appropriate levels to target
! when converting a floating point value into its integer vector
! representation. There are also some subtleties in reconstructing the
! final sum from the integer vector, as described below. For each
! floating point value, the exponent and mantissa are extracted using
! the fortran intrinsics 'exponent' and 'fraction'. The mantissa is
! then 'shifted' to match the exponent for a target level in the
! integer vector using the 'scale' intrinsic. 'int(X,i8)' is used
! for the conversion for the given level, and subtraction between
! this integer and the original 'shifted' value identifies the
! remainder that will be converted to an integer for the next level
! in the vector. The logic continues until the remainder is zero. As
! mentioned above, the only requirement, due to the implementation
! using these fortran intrinsics, is that floating point and integer
! models use the same base, e.g.  
!  radix(1.0_r8) == radix(1_i8)
! for real*8 and integer*8. If not, then the alternative algorithm
! DDPDD mentioned above and described below is used instead. The
! integer representation must also have enough digits for the
! potential growth of the sum for each level, so could conceivably be
! too small for a large number of summands. 
!
! Upper bounds on the total number of summands and on all intermediate
! sums are calculated as
!  <number of MPI tasks>*<max number of summands per MPI task>
! and
!  <number of MPI tasks>*<max number of summands per MPI task>
!   *<max absolute value over all nonzero summands>
! respectively. The maximum number of summands per MPI task and the
! maximum absolute value over all nonzero summands are global
! information that need to be determined with additional MPI
! collectives. The minimum nonzero absolute value summand is also
! global information. Fortunately, all of these can be determined with
! a single MPI_Allreduce call, so only one more than that required for
! the sum itself. (Note that, in actuality, the exponents of max and
! min summands are determined, and these are used to calculate bounds
! on the maximum and minimum, allowing the use of an MPI_INTEGER
! vector in the MPI_Allreduce call.)
!
! The actual code is made a little messier by (a) supporting summation
! of multiple fields without increasing the number of MPI_Allreduce
! calls, (b) supporting OpenMP threading of the local arithmetic, (c)
! allowing the user to specify estimates for the global information
! (to avoid the additional MPI_Allreduce), (d) including a check of
! whether user specified bounds were sufficient and, if not,
! determining the actual bounds and recomputing the sum, and (e)
! allowing the user to specify the maximum number of levels to use,
! potentially losing accuracy but still preserving reproducibility and
! being somewhat cheaper to compute.
!
! The conversion of the local summands to vectors of integers, the
! summation of the local vectors of integers, and the summation of the
! distributed vectors of integers will be exact (if optional parameters
! are not used to decrease the accuracy - see below). However, the
! conversion of the vector of integer representation to a floating
! point value may be subject to rounding errors. Before the
! conversion, the vector of integers is adjusted so that all elements
! have the same sign, and so that the value, in absolute value, at a
! given level is strictly less than what can be represented at the
! next lower level (larger exponent) and strictly greater than what
! can represented at the next higher level (smaller exponent). Since
! all elements have the same sign, the sign is set to positive
! temporarily and then restored when the conversion to floating point
! is complete. These are all integer operations, so no accuracy is
! lost. These adjustments eliminate the possibility of catastrophic
! cancellation. Also, when converting the individual elements to
! floating point values and summing them, the summation is now
! equivalent to concatenating the digits in the mantissas for the
! component summands. In consequence, in the final step when each
! element of this modified vector of integers is converted to a
! floating point value and added into the intermediate sum, any
! rounding is limited to the least significant digit representable
! in the final floating point sum.
!
! Any such rounding error will be sensitive to the particular floating
! values generated from the integer vector, and so will be 
! sensitive to the number of levels in the vector and the implicit
! exponent associated with each level, which are themselves functions
! of the numbers of MPI tasks and OpenMP threads and the number of
! digits representable in an integer. To avoid this sensitivity,
! (effectively) generate a new integer vector in which each component
! integer has a fixed number of significant digits (e.g.,
! digits(1.0_r8)) and generate the floating point values from these
! before summing. (See comments in code for more details.) This
! creates a sequence of floating point values to be summed that are
! independent of, for example, the numbers of MPI tasks and OpenMP
! threads or whether using integer*8 or integer*4 internal
! representations in the integer vector, and thus ensure
! reproducibility with respect to these options.
!
! Description of optional parameters for integer vector algorithm:
!-----------------------------------------------------------------
! The accuracy of the integer vector algorithm is controlled by the
! total number of levels of integer expansion. The algorithm
! calculates the number of levels that is required for the sum to be
! essentially exact. (The sum as represented by the integer expansion
! is exact, but roundoff may perturb the least significant digit of
! the returned floating point representation of the sum.) The optional
! parameter arr_max_levels can be used to override the calculated
! value for each field. The optional parameter arr_max_levels_out can
! be used to return the values used.
!
! The algorithm requires an upper bound on the maximum summand
! (in absolute value) for each field, and will calculate this internally
! using an MPI_Allreduce. However, if the optional parameters
! arr_max_levels and arr_gbl_max are both set, then the algorithm will
! use the values in arr_gbl_max for the upper bounds instead. If only
! arr_gbl_max is present, then the maxima are computed internally
! (and the specified values are ignored). The optional parameter
! arr_gbl_max_out can be used to return the values used.
!
! The algorithm also requires an upper bound on the number of
! local summands across all MPI tasks. (By definition, the number of
! local summands is the same for each field on a given MPI task, i.e.,
! the input parameter nsummands.) This will be calculated internally,
! using an MPI_Allreduce, but the value in the optional argument
! gbl_max_nsummands will be used instead if (1) it is present,
! (2) the value is > 0, and (3) the maximum values and required number
! of levels are also specified. (If the maximum values are calculated,
! then the same MPI_Allreduce is used to determine the maximum numbers
! of local summands.) The accuracy of the user-specified value is not
! checked. However, if set to < 1, the value will instead be calculated.
! If the optional parameter gbl_max_nsummands_out is present,
! then the value used (gbl_max_nsummands if >= 1; calculated otherwise)
! will be returned.
!
! If the user-specified upper bounds on maximum summands are
! inaccurate or if the user-specified upper bounds (maximum summands
! and number of local summands) and numbers of levels causes
! any of the global sums to have fewer than the expected
! number of significant digits, and if the optional parameter
! repro_sum_validate is NOT set to .false., then the algorithm will
! repeat the computations with internally calculated values for
! arr_max_levels, arr_gbl_max, and gbl_max_nsummands.
!
! If requested (by setting shr_reprosum_reldiffmax >= 0.0 and passing in
! the optional rel_diff parameter), results are compared with a
! nonreproducible floating point algorithm.
!
! Note that the cost of the algorithm is not strongly correlated with
! the number of levels, which primarily shows up as a (modest) increase
! in the cost of the MPI_Allreduce as a function of vector length.
! Rather the cost is more a function of (a) the number of integers
! required to represent an individual summand and (b) the number of
! MPI_Allreduce calls. The number of integers required to represent an
! individual summand is 1 or 2 when using 8-byte integers for 8-byte
! real summands when the number of local summands and number of MPI
! tasks are not too large. As the magnitude of either of these increase,
! the number of integers required increases. The number of
! MPI_Allreduce calls is either 2 (specifying nothing or just
! arr_max_levels and arr_gbl_max correctly) or 1 (specifying
! gbl_max_nsummands, arr_max_levels, and arr_gbl_max correctly).
! When specifying arr_max_nsummands, arr_max_levels, or arr_gbl_max
! incorrectly, 3 or 4 MPI_Allreduce calls will be required.
!
! Description of alternative (DDPDD) algorithm:
!----------------------------------------------
! The alternative algorithm is a minor modification of a parallel
! implementation of David Bailey's routine DDPDD by Helen He
! and Chris Ding. See, for example,
!  Y. He, and C. Ding, 'Using Accurate Arithmetics to Improve
!  Numerical Reproducibility and Stability in Parallel Applications,'
!  J. Supercomputing, vol. 18, no. 3, 2001, pp. 259–277
! and the citations therein. Bailey uses the Knuth trick to implement
! quadruple precision summation of double precision values with 10
! double precision operations. The advantage of this algorithm is that
! it requires a single MPI_Allreduce and is less expensive per summand
! than is the integer vector algorithm. The disadvantage is that it
! is not guaranteed to be reproducible (though it is reproducible
! much more often than is the standard floating point algorithm).
! This alternative is used when the optional parameter ddpdd_sum is
! set to .true. It is also used if the integer vector algorithm radix
! assumption does not hold.
!
!------------------------------------------------------------------------
!
! Arguments
!
      integer,  intent(in) :: nsummands  ! number of local summands
      integer,  intent(in) :: dsummands  ! declared first dimension
      integer,  intent(in) :: nflds      ! number of fields
      real(r8), intent(in) :: arr(dsummands,nflds)
                                         ! input array

      real(r8), intent(out):: arr_gsum(nflds)
                                         ! global sums

      logical,  intent(in),    optional :: ddpdd_sum
                                         ! use ddpdd algorithm instead
                                         ! of integer vector algorithm

      logical,  intent(in),    optional :: allow_infnan
         ! if .true., allow INF or NaN input values.
         ! if .false. (the default), then abort.

      real(r8), intent(in),    optional :: arr_gbl_max(nflds)
                                         ! upper bound on max(abs(arr))

      real(r8), intent(out),   optional :: arr_gbl_max_out(nflds)
                                         ! calculated upper bound on
                                         ! max(abs(arr))

      integer,  intent(in),    optional :: arr_max_levels(nflds)
                                         ! maximum number of levels of
                                         ! integer expansion to use

      integer,  intent(out),   optional :: arr_max_levels_out(nflds)
                                         ! output of number of levels of
                                         ! integer expansion to used

      integer,  intent(in),    optional :: gbl_max_nsummands
                                         ! maximum of nsummand over all
                                         ! MPI tasks

      integer,  intent(out),   optional :: gbl_max_nsummands_out
                                         ! calculated maximum nsummands
                                         ! over all MPI tasks

      integer,  intent(in),    optional :: gbl_count
                                         ! was total number of summands;
                                         ! now is ignored; use
                                         ! gbl_max_nsummands instead

      logical,  intent(in),    optional :: repro_sum_validate
         ! flag enabling/disabling testing that gmax and  max_levels are
         ! accurate/sufficient. Default is enabled.

      integer,  intent(inout), optional :: repro_sum_stats(6)
                                   ! increment running totals for
                                   !  (1) one-reduction repro_sum
                                   !  (2) two-reduction repro_sum
                                   !  (3) both types in one call
                                   !  (4) nonrepro_sum
                                   !  (5) global max nsummands reduction
                                   !  (6) global lor 3*nflds reduction

      real(r8), intent(out),   optional :: rel_diff(2,nflds)
                                         ! relative and absolute
                                         !  differences between integer
                                         !  vector and floating point sums

      integer,  intent(in),    optional :: commid
                                         ! MPI communicator

!
! Local workspace
!
      logical :: abort_inf_nan           ! flag indicating whether to
                                         !  abort if INF or NaN found in input
      logical :: use_ddpdd_sum           ! flag indicating whether to
                                         !  use shr_reprosum_ddpdd or not
      logical :: recompute               ! flag indicating need to
                                         !  determine gmax/gmin before
                                         !  computing sum
      logical :: validate                ! flag indicating need to
                                         !  verify gmax and max_levels
                                         !  are accurate/sufficient
      logical :: nan_check, inf_check    ! flag on whether there are
                                         !  NaNs and INFs in input array
      logical :: inf_nan_lchecks(3,nflds)! flags on whether there are
                                         !  NaNs, positive INFs, or negative INFs
                                         !  for each input field locally
      logical :: inf_nan_gchecks(3,nflds)! flags on whether there are
                                         !  NaNs, positive INFs, or negative INFs
                                         !  for each input field
      logical :: arr_gsum_infnan(nflds)  ! flag on whether field sum is a
                                         !  NaN or INF

      integer :: gbl_lor_red             ! global lor reduction? (0/1)
      integer :: gbl_max_red             ! global max reduction? (0/1)
      integer :: repro_sum_fast          ! 1 reduction repro_sum? (0/1)
      integer :: repro_sum_slow          ! 2 reduction repro_sum? (0/1)
      integer :: repro_sum_both          ! both fast and slow? (0/1)
      integer :: nonrepro_sum            ! nonrepro_sum? (0/1)

      integer :: nan_count, inf_count    ! local count of NaNs and INFs in
                                         !  input array
      integer :: omp_nthreads            ! number of OpenMP threads
      integer :: mpi_comm                ! MPI subcommunicator
      integer :: mypid                   ! MPI task ID (COMM_WORLD)
      integer :: tasks                   ! number of MPI tasks
      integer :: ierr                    ! MPI error return
      integer :: ifld, isum, ithread     ! loop variables
      integer :: max_nsummands           ! max nsummands over all MPI tasks
                                         !  or threads (used in both ways)

      integer, allocatable :: isum_beg(:), isum_end(:)
                                         ! range of summand indices for each
                                         !  OpenMP thread
      integer, allocatable :: arr_tlmin_exp(:,:)
                                         ! per thread local exponent minima
      integer, allocatable :: arr_tlmax_exp(:,:)
                                         ! per thread local exponent maxima
      integer :: arr_exp, arr_exp_tlmin, arr_exp_tlmax
                                         ! summand exponent and working min/max
      integer :: arr_lmin_exp(nflds)     ! local exponent minima
      integer :: arr_lmax_exp(nflds)     ! local exponent maxima
      integer :: arr_lextremes(0:nflds,2)! local exponent extrema
      integer :: arr_gextremes(0:nflds,2)! global exponent extrema

      integer :: arr_gmax_exp(nflds)     ! global exponents maxima
      integer :: arr_gmin_exp(nflds)     ! global exponents minima
      integer :: arr_max_shift           ! maximum safe exponent for
                                         !  value < 1 (so that sum does
                                         !  not overflow)
      integer :: max_levels(nflds)       ! maximum number of levels of
                                         !  integer expansion to use
      integer :: max_level               ! maximum value in max_levels
      integer :: extra_levels            ! number of extra levels needed
                                         !  to guarantee that sum over threads
                                         !  or tasks does not cause overflow

      real(r8) :: xmax_nsummands         ! real(max_nsummands,r8)
      real(r8) :: arr_lsum(nflds)        ! local sums
      real(r8) :: arr_gsum_fast(nflds)   ! global sum calculated using
                                         !  fast, nonreproducible,
                                         !  floating point alg.
      real(r8) :: abs_diff               ! absolute difference between
                                         !  integer vector and floating point
                                         !  sums
#ifdef _OPENMP
      integer omp_get_max_threads
      external omp_get_max_threads
#endif
!
!------------------------------------------------------------------------
!
! Initialize local statistics variables
      gbl_lor_red = 0
      gbl_max_red = 0
      repro_sum_fast = 0
      repro_sum_slow = 0
      repro_sum_both = 0
      nonrepro_sum = 0

! Set MPI communicator
      if ( present(commid) ) then
         mpi_comm = commid
      else
         mpi_comm = MPI_COMM_WORLD
      endif
      call t_barrierf('sync_repro_sum',mpi_comm)

! Check whether should abort if input contains NaNs or INFs
      abort_inf_nan = .not. repro_sum_allow_infnan
      if ( present(allow_infnan) ) then
         abort_inf_nan = .not. allow_infnan
      endif

! With Fujitsu always abort on NaNs or INFs in input
#ifdef CPRFJ
      abort_inf_nan = .true.
#endif

      call t_startf('shr_reprosum_INF_NaN_Chk')

! Initialize flags to indicate that no NaNs or INFs are present in the input data
      inf_nan_gchecks = .false.
      arr_gsum_infnan = .false.

      if (abort_inf_nan) then

! Check whether input contains NaNs or INFs, and abort if so
         nan_check = any(shr_infnan_isnan(arr))
         inf_check = any(shr_infnan_isinf(arr))

         if (nan_check .or. inf_check) then

            nan_count = count(shr_infnan_isnan(arr))
            inf_count = count(shr_infnan_isinf(arr))

            if ((nan_count > 0) .or. (inf_count > 0)) then
               call mpi_comm_rank(MPI_COMM_WORLD, mypid, ierr)
               write(s_logunit,37) real(nan_count,r8), real(inf_count,r8), mypid
37 format("SHR_REPROSUM_CALC: Input contains ",e12.5, &
          " NaNs and ", e12.5, " INFs on MPI task ", i7)
               call shr_sys_abort("shr_reprosum_calc ERROR: NaNs or INFs in input")
            endif

         endif

#ifndef CPRFJ
      else

! Determine whether any fields contain NaNs or INFs, and avoid processing them
! via integer expansions
         inf_nan_lchecks = .false.

         do ifld=1,nflds
            inf_nan_lchecks(1,ifld) = any(shr_infnan_isnan(arr(:,ifld)))
            inf_nan_lchecks(2,ifld) = any(shr_infnan_isposinf(arr(:,ifld)))
            inf_nan_lchecks(3,ifld) = any(shr_infnan_isneginf(arr(:,ifld)))
         end do

         call t_startf("repro_sum_allr_lor")
         call mpi_allreduce (inf_nan_lchecks, inf_nan_gchecks, 3*nflds, &
                             MPI_LOGICAL, MPI_LOR, mpi_comm, ierr)
         gbl_lor_red = 1
         call t_stopf("repro_sum_allr_lor")

         do ifld=1,nflds
            arr_gsum_infnan(ifld) = any(inf_nan_gchecks(:,ifld))
         enddo
#endif

      endif

      call t_stopf('shr_reprosum_INF_NaN_Chk')

! Check whether should use shr_reprosum_ddpdd algorithm
      use_ddpdd_sum = repro_sum_use_ddpdd
      if ( present(ddpdd_sum) ) then
         use_ddpdd_sum = ddpdd_sum
      endif

! Check whether intrinsic-based algorithm will work on this system
! (requires floating point and integer bases to be the same)
! If not, always use ddpdd.
      use_ddpdd_sum = use_ddpdd_sum .or. (radix(1.0_r8) /= radix(1_i8))

      if ( use_ddpdd_sum ) then

         call t_startf('shr_reprosum_ddpdd')

         call shr_reprosum_ddpdd(arr, arr_gsum, nsummands, dsummands, &
                              nflds, mpi_comm)
         repro_sum_fast = 1

         call t_stopf('shr_reprosum_ddpdd')

      else

         call t_startf('shr_reprosum_int')

! Get number of MPI tasks
         call mpi_comm_size(mpi_comm, tasks, ierr)

! Get number of OpenMP threads
#ifdef _OPENMP
         omp_nthreads = omp_get_max_threads()
#else
         omp_nthreads = 1
#endif

! See if have sufficient information to not require max/min allreduce
         recompute = .true.
         validate = .false.
         if ( present(arr_gbl_max) .and. present(arr_max_levels) ) then
            recompute = .false.

! Setting lower bound on max_level*nflds to be 64 to improve OpenMP
! performance for loopb in shr_reprosum_int
            max_level = (64/nflds) + 1
            do ifld=1,nflds
               if ((arr_gbl_max(ifld) >= 0.0_r8) .and. &
                   (arr_max_levels(ifld) > 0)) then

                  arr_gmax_exp(ifld)  = exponent(arr_gbl_max(ifld))
                  if (max_level < arr_max_levels(ifld)) &
                     max_level = arr_max_levels(ifld)

               else
                  recompute = .true.
               endif
            enddo

            if (.not. recompute) then

! Determine maximum number of summands in local phases of the
! algorithm
               call t_startf("repro_sum_allr_max")
               if ( present(gbl_max_nsummands) ) then
                  if (gbl_max_nsummands < 1) then
                     call mpi_allreduce (nsummands, max_nsummands, 1, &
                                         MPI_INTEGER, MPI_MAX, mpi_comm, ierr)
                     gbl_max_red = 1
                  else
                     max_nsummands = gbl_max_nsummands
                  endif
               else
                  call mpi_allreduce (nsummands, max_nsummands, 1, &
                                      MPI_INTEGER, MPI_MAX, mpi_comm, ierr)
                  gbl_max_red = 1
               endif
               call t_stopf("repro_sum_allr_max")

! Determine maximum shift. Shift needs to be small enough that summation,
! in absolute value, does not exceed maximum value representable by i8.

! If requested, return max_nsummands before it is redefined
               if ( present( gbl_max_nsummands_out) ) then
                  gbl_max_nsummands_out = max_nsummands
               endif

! Summing within each thread first (adding 1 to max_nsummands
! to ensure that integer division rounds up)
               max_nsummands = (max_nsummands/omp_nthreads) + 1
! then over threads and tasks
               max_nsummands = max(max_nsummands, tasks*omp_nthreads)
! A 'max' is used in the above calculation because the partial sum for
! each thread, calculated in shr_reprosum_int, is postprocessed so that
! each integer in the corresponding vector of integers is reduced in
! magnitude to be less than (radix(1_i8)**arr_max_shift). Therefore,
! the maximum shift can be calculated separately for per thread sums
! and sums over threads and tasks, and the smaller value used. This is
! equivalent to using max_nsummands as defined above.

               xmax_nsummands = real(max_nsummands,r8)
               arr_max_shift = digits(1_i8) - (exponent(xmax_nsummands) + 1)
               if (arr_max_shift < 2) then
                  call shr_sys_abort('repro_sum failed: number of summands too '// &
                                     'large for integer vector algorithm' )
               endif
! Note: by construction, each floating point value will be decomposed
! into a vector of integers each component of which will be strictly
! less than radix(1_i8)**arr_max_shift in absolute value, and the
! summation of max_nsummands of these, again in absolute value, will
! then be less than 
!  radix(1_i8)**(arr_max_shift + exponent(xmax_nsummands))
! or radix(1_i8)**(digits(1_i8) - 1). This is more conservative than
! necessary, but it also allows the postprocessing mentioned above
! (and described later) to proceed without danger of introducing
! overflow. 

! Determine additional number of levels needed to support the
! postprocessing that reduces the magnitude of each component
! of the integer vector of the partial sum for each thread
! to be less than (radix(1_i8)**arr_max_shift).
               extra_levels = (digits(1_i8) - 1)/arr_max_shift
! Extra levels are indexed by (-(extra_levels-1):0)
! Derivation of this is described in the comments in
! shr_reprosum_int.

! Calculate sum
               if (present(repro_sum_validate)) then
                  validate = repro_sum_validate
               else
                  validate = .true.
               endif
               call shr_reprosum_int(arr, arr_gsum, nsummands, dsummands, &
                                     nflds, arr_max_shift, arr_gmax_exp, &
                                     arr_max_levels, max_level, extra_levels, &
                                     arr_gsum_infnan, validate, recompute, &
                                     omp_nthreads, mpi_comm)

! Record statistics, etc.
               repro_sum_fast = 1
               if (recompute) then
                  repro_sum_both = 1
               else
! If requested, return specified levels and upper bounds on maxima
                  if ( present(arr_max_levels_out) ) then
                     do ifld=1,nflds
                        arr_max_levels_out(ifld) = arr_max_levels(ifld)
                     enddo
                  endif
                  if ( present(arr_gbl_max_out) ) then
                     do ifld=1,nflds
                       arr_gbl_max_out(ifld) = arr_gbl_max(ifld)
                     enddo
                  endif
               endif
            endif
         endif

! Do not have sufficient information; calculate global max/min and
! use to compute required number of levels
         if (recompute) then

! Record statistic
            repro_sum_slow = 1

! Determine maximum and minimum (non-zero) summand values and
! maximum number of local summands

! Allocate thread-specific work space
            allocate(arr_tlmax_exp(nflds,omp_nthreads))
            allocate(arr_tlmin_exp(nflds,omp_nthreads))
            allocate(isum_beg(omp_nthreads))
            allocate(isum_end(omp_nthreads))

! Split summand index range over OpenMP threads
            call split_indices(nsummands, omp_nthreads, isum_beg, isum_end)

!$omp parallel do      &
!$omp default(shared)  &
!$omp private(ithread, ifld, isum, arr_exp, arr_exp_tlmin, arr_exp_tlmax)
            do ithread=1,omp_nthreads
               call t_startf('repro_sum_loopa')
               do ifld=1,nflds
                  arr_exp_tlmin = MAXEXPONENT(1.0_r8)
                  arr_exp_tlmax = MINEXPONENT(1.0_r8)
                  if (.not. arr_gsum_infnan(ifld)) then
                     do isum=isum_beg(ithread),isum_end(ithread)
                        if (arr(isum,ifld) /= 0.0_r8) then
                           arr_exp = exponent(arr(isum,ifld))
                           arr_exp_tlmin = min(arr_exp,arr_exp_tlmin)
                           arr_exp_tlmax = max(arr_exp,arr_exp_tlmax)
                        endif
                     end do
                  endif
                  arr_tlmin_exp(ifld,ithread) = arr_exp_tlmin
                  arr_tlmax_exp(ifld,ithread) = arr_exp_tlmax
               end do
               call t_stopf('repro_sum_loopa')
            end do

            do ifld=1,nflds
               arr_lmax_exp(ifld) = maxval(arr_tlmax_exp(ifld,:))
               arr_lmin_exp(ifld) = minval(arr_tlmin_exp(ifld,:))
            end do
            deallocate(arr_tlmin_exp,arr_tlmax_exp,isum_beg,isum_end)

            arr_lextremes(0,:) = -nsummands
            arr_lextremes(1:nflds,1) = -arr_lmax_exp(:)
            arr_lextremes(1:nflds,2) = arr_lmin_exp(:)
            call t_startf("repro_sum_allr_minmax")
            call mpi_allreduce (arr_lextremes, arr_gextremes, 2*(nflds+1), &
                                MPI_INTEGER, MPI_MIN, mpi_comm, ierr)
            call t_stopf("repro_sum_allr_minmax")
            max_nsummands   = -arr_gextremes(0,1)
            arr_gmax_exp(:) = -arr_gextremes(1:nflds,1)
            arr_gmin_exp(:) =  arr_gextremes(1:nflds,2)

! If a field is identically zero or contains INFs or NaNs, arr_gmin_exp
! still equals MAXEXPONENT and arr_gmax_exp still equals MINEXPONENT.
! In this case, set arr_gmin_exp = arr_gmax_exp = MINEXPONENT
            do ifld=1,nflds
               arr_gmin_exp(ifld) = min(arr_gmax_exp(ifld),arr_gmin_exp(ifld))
            enddo

! If requested, return upper bounds on observed maxima
            if ( present(arr_gbl_max_out) ) then
               do ifld=1,nflds
                  arr_gbl_max_out(ifld) = scale(1.0_r8,arr_gmax_exp(ifld))
               enddo
            endif

! If requested, return max_nsummands before it is redefined
            if ( present( gbl_max_nsummands_out) ) then
               gbl_max_nsummands_out = max_nsummands
            endif

! Determine maximum shift (same as in previous branch, but with calculated
! max_nsummands). Shift needs to be small enough that summation, in absolute
! value, does not exceed maximum value representable by i8.

! Summing within each thread first (adding 1 to max_nsummands
! to ensure that integer division rounds up)
            max_nsummands = (max_nsummands/omp_nthreads) + 1
! then over threads and tasks
            max_nsummands = max(max_nsummands, tasks*omp_nthreads)
! A 'max' is used in the above calculation because the partial sum for
! each thread, calculated in shr_reprosum_int, is postprocessed so that
! each integer in the corresponding vector of integers is reduced in
! magnitude to be less than (radix(1_i8)**arr_max_shift). Therefore,
! the maximum shift can be calculated separately for per thread sums
! and sums over threads and tasks, and the smaller value used. This is
! equivalent to using max_nsummands as defined above.

            xmax_nsummands = real(max_nsummands,r8)
            arr_max_shift = digits(1_i8) - (exponent(xmax_nsummands) + 1)
            if (arr_max_shift < 2) then
               call shr_sys_abort('repro_sum failed: number of summands too '// &
                                  'large for integer vector algorithm' )
            endif
! Note: by construction, each floating point value will be decomposed
! into a vector of integers each component of which will be strictly
! less than radix(1_i8)**arr_max_shift in absolute value, and the
! summation of max_nsummands of these, again in absolute value, will
! then be less than 
!  radix(1_i8)**(arr_max_shift + exponent(xmax_nsummands))
! or radix(1_i8)**(digits(1_i8) - 1). This is more conservative than
! necessary, but it also allows the postprocessing mentioned above
! (and described later) to proceed without danger of introducing
! overflow. 

! Determine maximum number of levels required for each field.
! Need enough levels to represent both the smallest and largest
! nonzero summands (in absolute value), and any values in between.
! The number of digits from the most significant digit in the
! largest summand to the most significant digit in the smallest
! summand is (arr_gmax_exp(ifld)-arr_gmin_exp(ifld)), and the maximum
! number of digits needed to represent the smallest value is
! digits(1.0_r8). Divide this total number of digits by the number of
! digits per level (arr_max_shift) to get the number of levels
!  ((digits(1.0_r8) + (arr_gmax_exp(ifld)-arr_gmin_exp(ifld))) / arr_max_shift)
! with some tweaks:
!  + 1 because first truncation for any given summand probably does
!  not involve a maximal shift (but this adds only one to the total)
!  + 1 to guarantee that the integer division rounds up (not down)
! (setting lower bound on max_level*nflds to be 64 to improve OpenMP
!  performance for loopb in shr_reprosum_int)
            max_level = (64/nflds) + 1
            do ifld=1,nflds
               max_levels(ifld) = 2 + &
                ((digits(1.0_r8) + (arr_gmax_exp(ifld)-arr_gmin_exp(ifld))) &
                / arr_max_shift)
               if ( present(arr_max_levels) .and. (.not. validate) ) then
! If validate true, then computation with arr_max_levels failed
! previously
                  if ( arr_max_levels(ifld) > 0 ) then
                     max_levels(ifld) = &
                        min(arr_max_levels(ifld),max_levels(ifld))
                  endif
               endif
               if (max_level < max_levels(ifld)) &
                  max_level = max_levels(ifld)
            enddo

! If requested, return calculated levels
            if ( present(arr_max_levels_out) ) then
               do ifld=1,nflds
                  arr_max_levels_out(ifld) = max_levels(ifld)
               enddo
            endif

! Determine additional number of levels needed to support the
! postprocessing that reduces the magnitude of each component
! of the integer vector of the partial sum for each thread
! to be less than (radix(1_i8)**arr_max_shift).
            extra_levels = (digits(1_i8) - 1)/arr_max_shift
! Extra levels are indexed by (-(extra_levels-1):0)
! Derivation of this is described in the comments in
! shr_reprosum_int.

! Calculate sum
            validate = .false.
            call shr_reprosum_int(arr, arr_gsum, nsummands, dsummands, &
                                  nflds, arr_max_shift, arr_gmax_exp, &
                                  max_levels, max_level, extra_levels, &
                                  arr_gsum_infnan, validate, recompute, &
                                  omp_nthreads, mpi_comm)

         endif

         call t_stopf('shr_reprosum_int')

      endif

! Compare integer vector and floating point results
      if ( present(rel_diff) ) then
         if (shr_reprosum_reldiffmax >= 0.0_r8) then

            call t_barrierf('sync_nonrepro_sum',mpi_comm)
            call t_startf('nonrepro_sum')
! Record statistic
            nonrepro_sum = 1
! Compute nonreproducible sum
            arr_lsum(:) = 0.0_r8
!$omp parallel do      &
!$omp default(shared)  &
!$omp private(ifld, isum)
            do ifld=1,nflds
               if (.not. arr_gsum_infnan(ifld)) then
                  do isum=1,nsummands
                     arr_lsum(ifld) = arr(isum,ifld) + arr_lsum(ifld)
                  end do
               endif
            end do

            call t_startf("nonrepro_sum_allr_r8")
            call mpi_allreduce (arr_lsum, arr_gsum_fast, nflds, &
                                MPI_REAL8, MPI_SUM, mpi_comm, ierr)
            call t_stopf("nonrepro_sum_allr_r8")

            call t_stopf('nonrepro_sum')

! Determine differences
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

! Set field sums to NaN and INF, as needed
      do ifld=1,nflds
         if (arr_gsum_infnan(ifld)) then
            if (inf_nan_gchecks(1,ifld)) then
               ! NaN => NaN
               arr_gsum(ifld) = shr_infnan_nan
            else if (inf_nan_gchecks(2,ifld) .and. inf_nan_gchecks(3,ifld)) then
               ! posINF and negINF => NaN
               arr_gsum(ifld) = shr_infnan_nan
            else if (inf_nan_gchecks(2,ifld)) then
               ! posINF only => posINF
               arr_gsum(ifld) = shr_infnan_posinf
            else if (inf_nan_gchecks(3,ifld)) then
               ! negINF only => negINF
               arr_gsum(ifld) = shr_infnan_neginf
            endif
         endif
      end do

! Return statistics
      if ( present(repro_sum_stats) ) then
         repro_sum_stats(1) = repro_sum_stats(1) + repro_sum_fast
         repro_sum_stats(2) = repro_sum_stats(2) + repro_sum_slow
         repro_sum_stats(3) = repro_sum_stats(3) + repro_sum_both
         repro_sum_stats(4) = repro_sum_stats(4) + nonrepro_sum
         repro_sum_stats(5) = repro_sum_stats(5) + gbl_max_red
         repro_sum_stats(6) = repro_sum_stats(6) + gbl_lor_red
      endif

   end subroutine shr_reprosum_calc

!
!========================================================================
!

   subroutine shr_reprosum_int (arr, arr_gsum, nsummands, dsummands, nflds, &
                                arr_max_shift, arr_gmax_exp, max_levels,    &
                                max_level, extra_levels, skip_field,        &
                                validate, recompute, omp_nthreads, mpi_comm )
!------------------------------------------------------------------------
!
! Purpose:
! Compute the global sum of each field in 'arr' using the indicated
! communicator with a reproducible yet scalable implementation based
! on first converting each floating point summand into an equivalent
! representation using a vector of integers, summing the integer
! vectors, then converting the resulting sum back to a floating point
! representation. The accuracy of the integer vector algorithm is
! controlled by the number of 'levels' of integer expansion, the maximum
! value of which is specified by max_level.
!
!------------------------------------------------------------------------
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
      integer,  intent(in) :: extra_levels  ! number of extra levels
                                            ! needed to guarantee that
                                            ! sum over threads or tasks
                                            ! does not cause overflow
      integer,  intent(in) :: omp_nthreads  ! number of OpenMP threads
      integer,  intent(in) :: mpi_comm      ! MPI subcommunicator

      real(r8), intent(in) :: arr(dsummands,nflds)
                                            ! input array

      logical,  intent(in) :: skip_field(nflds)
         ! flag indicating whether the sum for this field should be
         ! computed or not (used to skip over fields containing
         ! NaN or INF summands)

      logical,  intent(in) :: validate
         ! flag indicating that accuracy of solution generated from
         ! arr_gmax_exp and max_levels should be tested

      logical,  intent(out):: recompute
         ! flag indicating that either the upper bounds are inaccurate,
         !  or max_levels and arr_gmax_exp do not generate accurate
         !  enough sums

      real(r8), intent(out):: arr_gsum(nflds)      ! global sums
!
! Local workspace
!
      integer, parameter  :: max_svlevel_factor = &
                                1 + (digits(1_i8)/digits(1.0_r8))

      integer(i8) :: i8_arr_tlsum_level(-(extra_levels-1):max_level,nflds,omp_nthreads)
                                   ! integer vector representing local
                                   !  sum (per thread, per field)
      integer(i8) :: i8_arr_lsum_level((max_level+extra_levels+2)*nflds)
                                   ! integer vector representing local
                                   !  sum
      integer(i8) :: i8_arr_level  ! integer part of summand for current
                                   !  expansion level
      integer(i8) :: i8_arr_gsum_level((max_level+extra_levels+2)*nflds)
                                   ! integer vector representing global
                                   !  sum
      integer(i8) :: i8_gsum_level(-(extra_levels-1):max_level)
                                   ! integer vector representing global
                                   !  sum for one field
      integer(i8) :: IX_8          ! integer representation of r8 value
      integer(i8) :: i8_sign       ! sign global sum
      integer(i8) :: i8_radix      ! radix for i8 variables (and r8
                                   !  variables by earlier if-test)

      integer :: max_error(nflds,omp_nthreads)
                                   ! accurate upper bound on data?
      integer :: not_exact(nflds,omp_nthreads)
                                   ! max_levels sufficient to
                                   !  capture all digits?
      integer :: isum_beg(omp_nthreads), isum_end(omp_nthreads)
                                   ! range of summand indices for each
                                   !  OpenMP thread
      integer :: ifld, isum, ithread, jlevel
                                   ! loop variables
      integer :: arr_exp           ! exponent of summand
      integer :: arr_shift         ! exponent used to generate integer
                                   !  for current expansion level
      integer :: ilevel            ! current integer expansion level
      integer :: offset(nflds)     ! beginning location in
                                   !  i8_arr_{g,l}sum_level for integer
                                   !  expansion of current ifld
      integer :: voffset           ! modification to offset used to
                                   !  include validation metrics
      integer :: min_level         ! index of minimum levels (including
                                   !  extra levels) for i8_arr_tlsum_level
      integer :: ioffset           ! offset(ifld)
      integer :: svlevel           ! number of summands in summand_vector
      integer :: ierr              ! MPI error return
      integer :: LX                ! exponent of X_8 (see below)
      integer :: veclth            ! total length of i8_arr_lsum_level
      integer :: i8_digit_count    ! number of digits in integer
                                   !  expansion of sum
      integer :: i8_begin_level    ! level starting from in
                                   !  creating next 'exactly representable'
                                   !  floating point value from modified
                                   !  integer expansion of the sum
      integer :: i8_trunc_level    ! level at which the number of digits in
                                   !  the modified integer expansion of the
                                   !  sum exceeds the number of representable
                                   !  digits in the floating point sum
      integer :: i8_trunc_loc      ! location of last digit at i8_trunc_level
                                   !  in the modified integer expansion of the
                                   !  sum that is representable in the floating
                                   !  point sum
      integer(i8) :: i8_trunc_level_rem
                                   ! truncated digits at i8_trunc_level
                                   !  in the modified integer expansion
                                   !  of the sum
      integer :: curr_exp          ! exponent of partial sum during
                                   !  reconstruction from integer vector
      integer :: corr_exp          ! exponent of current summand in
                                   !  reconstruction from integer vector

      real(r8) :: arr_frac         ! fraction of summand
      real(r8) :: arr_remainder    ! part of summand remaining after
                                   !  current level of integer expansion
      real(r8) :: X_8              ! r8 representation of current
                                   !  i8_arr_gsum_level
      real(r8) :: RX_8             ! r8 representation of (other)
                                   !  integers used in calculation.
      real(r8) :: summand_vector((max_level+extra_levels)*max_svlevel_factor)
                                   ! vector of r8 values generated from
                                   !  integer vector representation to be
                                   !  summed to generate global sum

      logical :: first_stepd_iteration
                                   ! flag used to indicate whether first
                                   !  time through process of converting
                                   !  vector of integers into a floating
                                   !  point value, as it requires
                                   !  special logic
!
!------------------------------------------------------------------------
! Save radix of i8 variables in an i8 variable
      i8_radix = radix(IX_8)

! If validating upper bounds, reserve space for validation metrics
! In both cases, reserve extra levels for overflows from the top level
      if (validate) then
        voffset = extra_levels + 2
      else
        voffset = extra_levels
      endif

! For convenience, define minimum level index for i8_arr_tlsum_level
      min_level = -(extra_levels-1)

! Compute offsets for each field
      offset(1) = voffset
      do ifld=2,nflds
         offset(ifld) = offset(ifld-1) &
                        + (max_levels(ifld-1) + voffset)
      enddo
      veclth = offset(nflds) + max_levels(nflds)

! Split summand index range over OpenMP threads
      call split_indices(nsummands, omp_nthreads, isum_beg, isum_end)

! Convert local summands to vector of integers and sum
! (Using scale instead of set_exponent because arr_remainder may not be
! 'normal' after level 1 calculation)
      i8_arr_lsum_level(:) = 0_i8

!$omp parallel do      &
!$omp default(shared)  &
!$omp private(ithread, ifld, ioffset, isum, arr_frac, arr_exp, &
!$omp         arr_shift, ilevel, i8_arr_level, arr_remainder, RX_8, IX_8)
      do ithread=1,omp_nthreads
       call t_startf('repro_sum_loopb')
       do ifld=1,nflds
          ioffset = offset(ifld)

          max_error(ifld,ithread) = 0
          not_exact(ifld,ithread) = 0
          i8_arr_tlsum_level(:,ifld,ithread) = 0_i8

          if (skip_field(ifld)) cycle

          do isum=isum_beg(ithread),isum_end(ithread)
            arr_remainder = 0.0_r8

            if (arr(isum,ifld) /= 0.0_r8) then
               arr_exp   = exponent(arr(isum,ifld))
               arr_frac  = fraction(arr(isum,ifld))

! Test that global maximum upper bound is an upper bound
               if (arr_exp > arr_gmax_exp(ifld)) then
                  max_error(ifld,ithread) = 1
                  exit
               endif

! Calculate first shift
               arr_shift = arr_max_shift - (arr_gmax_exp(ifld)-arr_exp)

! Determine first (probably) nonzero level (assuming initial fraction is
! 'normal' - algorithm still works if this is not true)
! NOTE: this is critical; scale will set to zero if min exponent is too small.
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

               if (ilevel <= max_levels(ifld)) then
! Apply first shift/truncate, add it to the relevant running
! sum, and calculate the remainder.
                  arr_remainder = scale(arr_frac,arr_shift)
                  i8_arr_level = int(arr_remainder,i8)
                  i8_arr_tlsum_level(ilevel,ifld,ithread) = &
                     i8_arr_tlsum_level(ilevel,ifld,ithread) + i8_arr_level
                  arr_remainder = arr_remainder - i8_arr_level

! While the remainder is non-zero, continue to shift, truncate,
! sum, and calculate new remainder
                  do while ((arr_remainder /= 0.0_r8) &
                     .and. (ilevel < max_levels(ifld)))
                     ilevel = ilevel + 1
                     arr_remainder = scale(arr_remainder,arr_max_shift)
                     i8_arr_level = int(arr_remainder,i8)
                     i8_arr_tlsum_level(ilevel,ifld,ithread) = &
                        i8_arr_tlsum_level(ilevel,ifld,ithread) + i8_arr_level
                     arr_remainder = arr_remainder - i8_arr_level
                  enddo

               endif
            endif

            if (arr_remainder /= 0.0_r8) then
               not_exact(ifld,ithread) = 1
            endif

          enddo
! Postprocess integer vector to eliminate possibility of overflow
! during subsequent sum over threads and tasks, as per earlier
! comment on logic behind definition of max_nsummands. If value at a
! given level is larger than or equal to
! (radix(1_i8)**arr_max_shift), subtract this 'overlap' from the
! current value and add it (appropriately shifted) to the value at
! the next smaller level in the vector.
! (a) As described earlier, prior to this postprocessing the integer
!     components are each strictly less than
!     radix(1_i8)**(digits(1_i8) - 1) in absolute value. So, after
!     shifting, the absolute value of the amount added to level
!     max_levels(ifld)-1 from level max_levels(ifld) is less than
!     radix(1_i8)**(digits(1_i8) - 1 - arr_max_shift) with the
!     resulting sum, in absolute value, being less than
!      (radix(1_i8)**(digits(1_i8) - 1))*(1 + radix(1_i8)**(-arr_max_shift)).
!     Any overlap from this component is then added to the level
!     max_levels(ifld)-2, etc., with resulting intermediate sums, in
!     absolute value, for levels 1 to max_levels(ifld) being bounded
!     from above by 
!      (radix(1_i8)**(digits(1_i8) - 1))*sum{i=0,inf}(radix(1_i8)**(-i*arr_max_shift)).
!     Since radix(1_i8) >= 2 and arr_max_shift is also required to be
!     >= 2 (otherwise the code exits with an error) this is less than
!     or equal to 
!      (radix(1_i8)**(digits(1_i8) - 1))*sum{i=0,inf}(2**(-2i)),
!     or
!      (radix(1_i8)**(digits(1_i8) - 1))*(4/3).
!     In summary, this shows that no absolute value generated during
!     this process will exceed the maximum value representable in i8,
!     i.e. (radix(1_i8)**(digits(1_i8)) - 1), as long as
!     digits(1_i8) >= 2. 
! (b) 'ilevel==0,...,-(extra_levels-1)' correspond to extra levels
!     used to continue the above process until values at all levels
!     are less than radix(1_i8)**arr_max_shift in absolute value
!     (except level -(extra_levels-1), as described below). The
!     result of shifting the overlap from level 1 to level 0, which
!     is initially zero, is bounded in absolute value by 
!      (radix(1_i8)**(digits(1_i8) - 1 - arr_max_shift))*(4/3).
!     After removing any overlap from level 0, the upper bound for
!     level -1, which is also initially zero, is
!      (radix(1_i8)**(digits(1_i8) - 1 - 2*arr_max_shift))*(4/3).
!     Continuing the process, when get to level -(extra_levels-1),
!     the upper bound is
!      (radix(1_i8)**(digits(1_i8) - 1 - extra_levels*arr_max_shift))*(4/3).
!     If we define
!      extra_levels = ceiling[(digits(1_i8) - 1)/arr_max_shift - 1]
!     then the upper bound is
!      (radix(1_i8)**(arr_max_shift))*(4/3).
!     Setting 
!      extra_levels = (digits(1_i8) - 1)/arr_max_shift
!     is then a slightly conservative estimate that achieves the same
!     upper bound. While the above upper bound at level
!     -(extra_levels-1)is a factor of (4/3) larger than the target
!     radix(1_i8)**arr_max_shift, it is still small enough so that
!     the sum over threads and tasks, bounded from above in absolute
!     value by 
!      (radix(1_i8)**(digits(1_i8) - 1))*(4/3),
!     will not cause an overflow at level -(extra_levels-1) as long as
!     digits(1_i8) >= 2.
          do ilevel=max_levels(ifld),min_level+1,-1
             if (abs(i8_arr_tlsum_level(ilevel,ifld,ithread)) >= &
                    (i8_radix**arr_max_shift)) then
                
                IX_8 = i8_arr_tlsum_level(ilevel,ifld,ithread) &
                       / (i8_radix**arr_max_shift)
                i8_arr_tlsum_level(ilevel-1,ifld,ithread) = &
                   i8_arr_tlsum_level(ilevel-1,ifld,ithread) + IX_8
                
                IX_8 = IX_8*(i8_radix**arr_max_shift)
                i8_arr_tlsum_level(ilevel,ifld,ithread)   = &
                   i8_arr_tlsum_level(ilevel,ifld,ithread) - IX_8
             endif
          enddo
       enddo
       call t_stopf('repro_sum_loopb')
      enddo

! Sum contributions from different threads
      do ifld=1,nflds
         ioffset = offset(ifld)
         do ithread = 1,omp_nthreads
            do ilevel = min_level,max_levels(ifld)
               i8_arr_lsum_level(ioffset+ilevel) = &
                  i8_arr_lsum_level(ioffset+ilevel) &
                  + i8_arr_tlsum_level(ilevel,ifld,ithread)
            enddo
         enddo
      enddo

! Record if upper bound was inaccurate or if level expansion stopped
! before full accuracy was achieved
      if (validate) then
         do ifld=1,nflds
            ioffset = offset(ifld)
            i8_arr_lsum_level(ioffset-voffset+1) = maxval(max_error(ifld,:))
            i8_arr_lsum_level(ioffset-voffset+2) = maxval(not_exact(ifld,:))
         enddo
      endif

! Sum integer vector element-wise
#if ( defined noI8 )
     ! Workaround for when shr_kind_i8 is not supported.
      call t_startf("repro_sum_allr_i4")
      call mpi_allreduce (i8_arr_lsum_level, i8_arr_gsum_level, &
                          veclth, MPI_INTEGER, MPI_SUM, mpi_comm, ierr)
      call t_stopf("repro_sum_allr_i4")
#else
      call t_startf("repro_sum_allr_i8")
      call mpi_allreduce (i8_arr_lsum_level, i8_arr_gsum_level, &
                          veclth, MPI_INTEGER8, MPI_SUM, mpi_comm, ierr)
      call t_stopf("repro_sum_allr_i8")
#endif

      call t_startf('repro_sum_finalsum')
! Construct global sum from integer vector representation:
!  1) arr_max_shift is the shift applied to fraction(arr_gmax) .
!   When shifting back, need to 'add back in' the true arr_gmax exponent.
!   This was removed implicitly by working only with the fraction.
!  2) To avoid the possibility of catastrophic cancellation, and
!   an unacceptable floating point rounding error, can do some arithmetic
!   with the integer vector so that all components have the same sign.
!  3) If convert each integer in the integer vector to a floating
!   point value and then add these together, smallest to largest, to
!   calculate the final sum, there may be roundoff error in the least
!   significant digit. This error will be sensitive to the particular
!   floating values generated from the integer vector, and so will be
!   sensitive to the number of levels in the vector and the implicit
!   exponent associated with each level. So this approach is not
!   guaranteed to be reproducible with respect to a change in the
!   number of MPI tasks and OpenMP threads (as this changes the
!   definition of max_nsummands, and thus also arr_max_shift). It is
!   also not guaranteed to be reproducible with respect to changing
!   the integer size, e.g. from i8 to i4, as this also changes
!   arr_max_shift. However, can eliminate this potential loss of
!   reproducibility by taking the following steps. 
!   a) Manipulate the integer vector so that
!    i)  the component values do not 'overlap', that is, the value
!        represented by a component is strictly less than the value
!        represented by the least significant digit in the previous
!        component, and 
!    ii) all components are positive (saving the sign to be restored
!        to the final result).
!   b) Identify the digit in the resulting integer vector that is the
!      last representable in the floating point representation, then
!      truncate the vector at this point, i.e., all digits of lesser
!      significance in the given component and all components
!      representing digits of lesser significance (call this the
!      remainder). 
!   c) Convert each integer component in the modified integer vector
!      to its corresponding floating point value and sum the
!      sequence. (Order is unimportant, as explained below, but here
!      add largest to smallest.) 
!   d) Repeat (b) and (c) for the remainder (recursively, as
!      necessary).
!   e) Sum all floating point numbers generated by step (c), smallest
!      to largest. 
!   f) Restore the sign.
!   With the manipulations in (a) and (b), the summation in (c) is
!   equivalent to concatenating the digits in the mantissas for the
!   component summands, so rounding is irrelevant (so far). Repeating
!   this with the remainder(s) generates a sequence of 'exact'
!   floating point numbers. Summing these can still generate a
!   rounding error in the least significant digit in the largest
!   floating point value (which is the last representable digit in the
!   final result), but the floating point values being summed and
!   order of summation are independent of the number of levels and
!   implicit exponents, so reproducibility is ensured.
!  
!  Note that assignment of an i8 integer value to an r8 floating point
!  variable in step (c) can lead to a loss of accuracy because the
!  maximum number of digits in the i8 integer can be greater than the
!  maximum number of digits representable in the r8 variable (if the
!  xmax_nsummands correction is not very large). With the same sign
!  and nonoverlapping properties of the integer components, these lost
!  digits will also not be representable in the final sum. The process
!  described above of truncating at this last representable digit, and
!  then separately generating floating point value(s) for the
!  remainder, takes care of this automatically. Similar reasoning
!  applies to r4 floating point values with either i8 or i4 integer
!  components.

      recompute = .false.
      do ifld=1,nflds
         arr_gsum(ifld) = 0.0_r8
         ioffset = offset(ifld)
         svlevel = 0

! If validate is .true., test whether the summand upper bound
! was exceeded on any of the MPI tasks
         if (validate) then
            if (i8_arr_gsum_level(ioffset-voffset+1) /= 0_i8) then
               recompute = .true.
            endif
         endif

         if (.not. recompute) then
! Copy integer vector for current field from i8_arr_gsum_level, so that
! can be modified without changing i8_arr_gsum_level. (Preserving
! i8_arr_gsum_level unchanged is not necessary, but is convenient for debugging
! and makes indexing clearer and less error prone.)
            i8_gsum_level(:) = 0_i8
            do ilevel=min_level,max_levels(ifld)
               i8_gsum_level(ilevel) = i8_arr_gsum_level(ioffset+ilevel)
            enddo

! Preprocess integer vector (as described in 3(a) above):
!  i) If value larger than or equal to (radix(1_i8)**arr_max_shift),
!     add this 'overlap' to the value at the next smaller level
!     in the vector, resulting in nonoverlapping ranges for each
!     component.
!
!     As before, no intermediate sums for levels
!     max_levels(ifld) to -(extra_levels-2), in absolute value,
!     will exceed the the maximum value representable in i8, but the
!     upper bound on the final sum, in absolute value, at
!     level -(extra_levels-1) is now 
!      (radix(1_i8)**(digits(1_i8) - 1))*(4/3) + 
!            + sum{i=1,inf}(radix(1_i8)**(-i*arr_max_shift))
!      = (radix(1_i8)**(digits(1_i8) - 1))*
!            ((4/3) + sum{i=1,inf}(radix(1_i8)**(-i*arr_max_shift)).
!     which is less than or equal to 
!      (radix(1_i8)**(digits(1_i8) - 1))*((4/3) + (1/3))
!     or
!      (radix(1_i8)**(digits(1_i8) - 1))*(5/3)
!     which will not cause an overflow at level -(extra_levels-1)
!     as long as digits(1_i8) >= 3.
!
!     Since the exponents associated with each successive level
!     differ by arr_max_shift, monotonically decreasing with
!     increasing level, the absolute value at each level after this
!     preprocessing is strictly less than what can be represented at
!     the next lower level (larger exponent). If nonzero, it is also
!     strictly greater than what is represented at the next higher
!     level (smaller exponent). Note that the smallest level,
!     -(extra_levels-1), does not have to be less than
!     (radix(1_i8)**arr_max_shift) for this 'nonoverlap' property to
!     hold. 
            do ilevel=max_levels(ifld),min_level+1,-1
               if (abs(i8_gsum_level(ilevel)) >= &
                      (i8_radix**arr_max_shift)) then
                  
                  IX_8 = i8_gsum_level(ilevel) &
                         / (i8_radix**arr_max_shift)
                  i8_gsum_level(ilevel-1) = &
                     i8_gsum_level(ilevel-1) + IX_8
                  
                  IX_8 = IX_8*(i8_radix**arr_max_shift)
                  i8_gsum_level(ilevel)   = &
                     i8_gsum_level(ilevel) - IX_8
               endif
            enddo

!  ii) Working consecutively from the first level with a nonzero value
!      up to level max_levels(ifld), subtract +/- 1 from level with
!      larger exponent (e.g., ilevel) and add  +/-
!      (i8_radix**arr_max_shift) to level with smaller exponent
!      (ilevel+1), when necessary, so that the value at ilevel+1
!      has the same sign as the value at ilevel. Treat a zero value at
!      ilevel+1 as always a different sign from the value at ilevel so
!      that the process always makes this nonzero. (Otherwise, the
!      wrong sign could be reintroduced by subtracting from a zero
!      value at the next step.) When finished with the process values
!      at all levels are either greater than or equal to zero or all
!      are less than or equal to zero. Note that this can decrease
!      (but not increase) the absolute value at level
!      -(extra_levels-1) by 1. All other levels are now less than or
!      equal to (radix(1_i8)**arr_max_shift) in absolute value rather
!      than strictly less than. 
            ilevel = min_level
            do while ((i8_gsum_level(ilevel) == 0_i8) &
                      .and. (ilevel < max_levels(ifld)))
               ilevel = ilevel + 1
            enddo
!
            if (i8_gsum_level(ilevel) < 0_i8) then
               i8_sign = -1_i8
            else
               i8_sign = 1_i8
            endif
!            
            if (ilevel < max_levels(ifld)) then
               do jlevel=ilevel,max_levels(ifld)-1
                  if ((sign(1_i8,i8_gsum_level(jlevel)) &
                       /= sign(1_i8,i8_gsum_level(jlevel+1)))&
                      .or. (i8_gsum_level(jlevel+1) == 0_i8)) then
                     i8_gsum_level(jlevel)   = &
                       i8_gsum_level(jlevel) - i8_sign
                     i8_gsum_level(jlevel+1) = &
                        i8_gsum_level(jlevel+1) &
                        + i8_sign*(i8_radix**arr_max_shift)
                  endif
               enddo
            endif

!  iii) If 'same sign' is negative, then change to positive
!       temporarily. 
            if (i8_sign < 0_i8) then
               do jlevel=ilevel,max_levels(ifld)
                  i8_gsum_level(jlevel) = -i8_gsum_level(jlevel)
               enddo
            endif

!  iv) Nonoverlap property can be lost after imposition of same sign
!      over components. Reintroduce this property (retaining same sign
!      property). Note that carryover is never more than '1' to the
!      next smaller level, so, again, no intermediate or final sums
!      will exceed the maximum value representable in i8, including
!      level -(extra_levels-1) as long as digits(1_i8) >= 4. 
            do ilevel=max_levels(ifld),min_level+1,-1
               if (abs(i8_gsum_level(ilevel)) >= &
                      (i8_radix**arr_max_shift)) then
                  
                  IX_8 = i8_gsum_level(ilevel)/ &
                         (i8_radix**arr_max_shift)
                  i8_gsum_level(ilevel-1) = &
                     i8_gsum_level(ilevel-1) + IX_8
                  
                  IX_8 = IX_8*(i8_radix**arr_max_shift)
                  i8_gsum_level(ilevel)   = &
                     i8_gsum_level(ilevel) - IX_8
               endif
            enddo

! Step 3(d): iterate over steps 3(b) and 3(c), truncating integer
! vector to 'fit' into a floating point value, then repeating with
! remainder
            first_stepd_iteration = .true.
            arr_shift = arr_gmax_exp(ifld) - (min_level)*arr_max_shift
            i8_digit_count = 0
            i8_begin_level = min_level
            do while (i8_begin_level <= max_levels(ifld))
   
! Determine at which level the total number of integer digits equals
! or exceeds the number of digits representable in the floating point
! sum. Then determine which digit at this level is the last
! representable in the floating point sum. Note that this location
! (i8_trunc_loc) is zero-based, i.e. smallest digit is at location
! 0. Note that the exponent is a count of the number of digits for the
! first nonzero level. All subsequent levels contribute arr_max_shift
! digits.
               i8_trunc_loc = 0
               i8_trunc_level = max_levels(ifld)
               do ilevel=i8_begin_level,max_levels(ifld)
                  if (first_stepd_iteration) then
! Special logic for first time through. Subsequent iterations treat
! leading zeroes as significant.                  
                     if (i8_digit_count == 0) then
                        if (i8_gsum_level(ilevel) /= 0_i8) then
                           X_8 = i8_gsum_level(ilevel)
                           LX  = exponent(X_8)
! Note that even if i8_gsum_level(ilevel) is truncated when assigned
! to X_8, the exponent LX will still capture the original number of
! digits.
                        else
                           LX = 0
                        endif
                     else
                        LX = arr_max_shift
                     endif
                  else
! If i8_digit_count /= 0 during the first iteration
! (ilevel == i8_begin_level), then there is a remainder left at the
! previous i8_trunc_level and LX should be set to zero for this
! iteration.
                     if ((ilevel == i8_begin_level) .and. (i8_digit_count /= 0)) then
                        LX = 0
                     else
                        LX = arr_max_shift
                     endif
                  endif
                  if (i8_digit_count + LX >= digits(1.0_r8)) then
                     i8_trunc_level = ilevel
                     i8_trunc_loc = (i8_digit_count + LX) - digits(1.0_r8)
                     exit
                  else
                     i8_digit_count = i8_digit_count + LX
                  endif
               enddo
               first_stepd_iteration = .false.

! Truncate at i8_trunc_loc as needed and determine what the remainder
! is. 
               if (i8_trunc_loc == 0) then
! No truncation is necessary, and remainder is just the components
! for the remaining levels
                  i8_trunc_level_rem = 0
               else
! Shift right to identify the digits to be preserved and truncate
! there
                  IX_8 = i8_gsum_level(i8_trunc_level)/ &
                         (i8_radix**i8_trunc_loc)
! Shift left to put digits in the correct location (right fill with
! zeroes)
                  IX_8 = IX_8*(i8_radix**i8_trunc_loc)
! Calculate local remainder
                  i8_trunc_level_rem = (i8_gsum_level(i8_trunc_level) - IX_8)
! Update level with the truncated value
                  i8_gsum_level(i8_trunc_level) = IX_8
               endif
                  
! Calculate floating point value corresponding to modified integer
! vector. Note that, by construction, i8 integer value will fit into
! r8 floating point value, so do not need to test for this.
               svlevel = svlevel + 1
               summand_vector(svlevel) = 0.0_r8
               do ilevel=i8_begin_level,i8_trunc_level
                  if (i8_gsum_level(ilevel) /= 0_i8) then

! Convert integer to floating point representation
                     X_8 = i8_gsum_level(ilevel)
                     LX  = exponent(X_8)

! Add to vector of floating point summands, scaling first if exponent
! is too small to apply directly
                     curr_exp = LX + arr_shift
                     if (curr_exp >= MINEXPONENT(1.0_r8)) then
                        summand_vector(svlevel) = &
                           summand_vector(svlevel) + set_exponent(X_8,curr_exp)
                     else
                        RX_8 = set_exponent(X_8, &
                                            curr_exp-MINEXPONENT(1.0_r8))
                        summand_vector(svlevel) = &
                           summand_vector(svlevel) + scale(RX_8,MINEXPONENT(1.0_r8))
                     endif
                     
                  endif

! Note that same arr_shift should be used for next 'step 3(d)'
! iteration if i8_trunc_loc > 0.
                  if ((ilevel < i8_trunc_level) .or. (i8_trunc_loc == 0)) then
                     arr_shift = arr_shift - arr_max_shift
                  endif

               enddo

               if (i8_trunc_loc == 0) then
                  i8_digit_count = 0
                  i8_begin_level = i8_trunc_level + 1
               else
                  i8_digit_count = i8_trunc_loc
                  i8_begin_level = i8_trunc_level
                  i8_gsum_level(i8_trunc_level) = i8_trunc_level_rem
               endif
            
            enddo

! Step 3(e): sum vector of floating point values, smallest to largest
            arr_gsum(ifld) = 0.0_r8
            do jlevel=svlevel,1,-1
               arr_gsum(ifld) = arr_gsum(ifld) + summand_vector(jlevel)
            enddo

! Step 3(f): restore the sign
            arr_gsum(ifld) = i8_sign*arr_gsum(ifld)

! If validate is .true. and some precision lost, test whether 'too
! much' was lost, due to too loose an upper bound, too stringent a
! limit on number of levels of expansion, cancellation, ...
! Calculated by comparing lower bound on number of significant digits
! with number of digits in 1.0_r8 .
            if (validate) then
               if (i8_arr_gsum_level(ioffset-voffset+2) /= 0_i8) then

! Find first nonzero level and use exponent for this level, then
! assume all subsequent levels contribute arr_max_shift digits.
                  i8_digit_count = 0
                  do ilevel=min_level,max_levels(ifld)
                     if (i8_digit_count == 0) then
                        if (i8_arr_gsum_level(ioffset+ilevel) /= 0_i8) then
                           X_8 = i8_arr_gsum_level(ioffset+ilevel)
                           LX  = exponent(X_8)
                           i8_digit_count = LX
                        endif
                     else
                        i8_digit_count = i8_digit_count + arr_max_shift
                     endif
                  enddo

                  if (i8_digit_count < digits(1.0_r8)) then
                     recompute = .true.
                  endif
               endif
            endif

         endif

      enddo
      call t_stopf('repro_sum_finalsum')

   end subroutine shr_reprosum_int

!
!========================================================================
!

   logical function shr_reprosum_tolExceeded (name, nflds, master, &
                                              logunit, rel_diff    )
!------------------------------------------------------------------------
!
! Purpose:
! Test whether distributed sum exceeds tolerance and print out a
! warning message.
!
!------------------------------------------------------------------------
!
! Arguments
!
      character(len=*), intent(in) :: name    ! distributed sum identifier
      integer,  intent(in) :: nflds           ! number of fields
      logical,  intent(in) :: master          ! MPI task that will write
                                              !  warning messages?
      integer, optional, intent(in) :: logunit! unit warning messages
                                              !  written to
      real(r8), intent(in) :: rel_diff(2,nflds)
                                              ! relative and absolute
                                              !  differences between integer
                                              !  vector and floating point sums

!
! Local workspace
!
      integer  :: llogunit                    ! local log unit
      integer  :: ifld                        ! field index
      integer  :: exceeds_limit               ! number of fields whose
                                              !  sum exceeds tolerance
      real(r8) :: max_rel_diff                ! maximum relative difference
      integer  :: max_rel_diff_idx            ! field index for max. rel. diff.
      real(r8) :: max_abs_diff                ! maximum absolute difference
      integer  :: max_abs_diff_idx            ! field index for max. abs. diff.
!
!------------------------------------------------------------------------
!
      shr_reprosum_tolExceeded = .false.
      if (shr_reprosum_reldiffmax < 0.0_r8) return

      if ( present(logunit) ) then
         llogunit = logunit
      else
         llogunit = s_logunit
      endif

! Check that 'fast' reproducible sum is accurate enough.
      exceeds_limit = 0
      max_rel_diff = 0.0_r8
      max_abs_diff = 0.0_r8
      max_rel_diff_idx = 0
      do ifld=1,nflds
         if (rel_diff(1,ifld) > shr_reprosum_reldiffmax) then
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
            write(llogunit,*) trim(name), &
                            ': difference between integer vector and floating point sums ', &
                            ' exceeds tolerance in ', exceeds_limit, &
                            ' fields.'
            write(llogunit,*) '  Maximum relative diff: (rel)', &
                            rel_diff(1,max_rel_diff_idx), ' (abs) ', &
                            rel_diff(2,max_rel_diff_idx)
            write(llogunit,*) '  Maximum absolute diff: (rel)', &
                            rel_diff(1,max_abs_diff_idx), ' (abs) ', &
                            rel_diff(2,max_abs_diff_idx)
         endif
         shr_reprosum_tolExceeded = .true.
       endif

   end function shr_reprosum_tolExceeded

!
!========================================================================
!

   subroutine shr_reprosum_ddpdd (arr, arr_gsum, nsummands, dsummands,  &
                                  nflds, mpi_comm                       )
!------------------------------------------------------------------------
!
! Purpose:
! Compute the global sum of each field in 'arr' using the indicated
! communicator with a reproducible yet scalable implementation based
! on He and Ding's implementation of the double-double algorithm.
!
!------------------------------------------------------------------------
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
!------------------------------------------------------------------------
!
      call shr_reprosumx86_fix_start (old_cw)

      if (first_time) then
         call mpi_op_create(ddpdd, .true., mpi_sumdd, ierr)
         first_time = .false.
      endif

      do ifld=1,nflds
         arr_lsum_dd(ifld) = (0.0_r8,0.0_r8)

         do isum=1,nsummands

! Compute arr(isum,ifld) + arr_lsum_dd(ifld) using Knuth''s trick.
            t1 = arr(isum,ifld) + real(arr_lsum_dd(ifld))
            e  = t1 - arr(isum,ifld)
            t2 = ((real(arr_lsum_dd(ifld)) - e) &
                  + (arr(isum,ifld) - (t1 - e))) &
                 + aimag(arr_lsum_dd(ifld))

            ! The result is t1 + t2, after normalization.
            arr_lsum_dd(ifld) = cmplx ( t1 + t2, t2 - ((t1 + t2) - t1), r8 )
         enddo

      enddo

      call t_startf("repro_sum_allr_c16")
      call mpi_allreduce (arr_lsum_dd, arr_gsum_dd, nflds, &
                          MPI_COMPLEX16, mpi_sumdd, mpi_comm, ierr)
      call t_stopf("repro_sum_allr_c16")

      do ifld=1,nflds
         arr_gsum(ifld) = real(arr_gsum_dd(ifld))
      enddo

      call shr_reprosumx86_fix_end (old_cw)

   end subroutine shr_reprosum_ddpdd
!
!------------------------------------------------------------------------
!
   subroutine DDPDD (dda, ddb, len, itype)
!------------------------------------------------------------------------
!
! Purpose:
! Modification of original codes written by David H. Bailey
! This subroutine computes ddb(i) = dda(i)+ddb(i)
!
!------------------------------------------------------------------------
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
!------------------------------------------------------------------------
!
      do i = 1, len

! Compute dda + ddb using Knuth's trick.
         t1 = real(dda(i)) + real(ddb(i))
         e  = t1 - real(dda(i))
         t2 = ((real(ddb(i)) - e) + (real(dda(i)) - (t1 - e))) &
              + aimag(dda(i)) + aimag(ddb(i))

! The result is t1 + t2, after normalization.
         ddb(i) = cmplx ( t1 + t2, t2 - ((t1 + t2) - t1), r8 )
      enddo

   end subroutine DDPDD
!
!------------------------------------------------------------------------
!
  subroutine split_indices(total,num_pieces,ibeg,iend)
!------------------------------------------------------------------------
!
! Purpose:
! Split range into 'num_pieces'
!
!------------------------------------------------------------------------
!
! Arguments
!
    integer, intent(in) :: total
    integer, intent(in) :: num_pieces
    integer, intent(out) :: ibeg(num_pieces), iend(num_pieces)
!
! Local workspace
!
    integer :: itmp1, itmp2, ioffset, i
!
!------------------------------------------------------------------------
!
    itmp1 = total/num_pieces
    itmp2 = mod(total,num_pieces)
    ioffset = 0
    do i=1,itmp2
       ibeg(i) = ioffset + 1
       iend(i) = ioffset + (itmp1+1)
       ioffset = iend(i)
    enddo
    do i=itmp2+1,num_pieces
       ibeg(i) = ioffset + 1
       if (ibeg(i) > total) then
          iend(i) = ibeg(i) - 1
       else
          iend(i) = ioffset + itmp1
          ioffset = iend(i)
       endif
    enddo

  end subroutine split_indices
!
!========================================================================
!
end module shr_reprosum_mod
