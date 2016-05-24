
      module mo_grid
!---------------------------------------------------------------------
! 	... Basic grid point resolution parameters
!---------------------------------------------------------------------
      implicit none

      save

      integer, parameter :: &
                pcnst    = PCNST+1, &     ! number of advected constituents including cloud water
                pcnstm1  = PCNST, &     ! number of advected constituents excluding cloud water
                plev     = PLEV, &         ! number of vertical levels
                plevp    = plev+1, &      ! plev plus 1
                plevm    = plev-1, &      ! plev minus 1
                plon     = PLON, &         ! number of longitudes
                plat     = PLAT            ! number of latitudes

      integer, parameter :: &
                pnats    = GRPCNT    ! number of non-advected trace species

#ifdef STRAT_CHEM
      integer, parameter :: &
                phmu     = PCNST      ! number of long-lived species 
#endif

      integer :: nodes                ! mpi task count
      integer :: plonl                ! longitude tile dimension
      integer :: pplon                ! longitude tile count
      integer :: plnplv               ! plonl * plev

      end module mo_grid
