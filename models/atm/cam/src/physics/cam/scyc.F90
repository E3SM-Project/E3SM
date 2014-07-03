
module scyc

   !----------------------------------------------------------------------- 
   ! Purpose: 
   ! Provide logical control functions for the prognostic sulfur cycle
   ! parameterization.
   !
   ! Author: B. Eaton
   !-----------------------------------------------------------------------

   use cam_logfile, only: iulog
   implicit none
   save
   private
   public :: &
      scycini,       &! initialize control functions for sulfur cycle
      doscyc,        &! returns .true. when prognostic sulfur cycle is enabled
      useGEIA,       &! returns .true. when using 1985 GEIA sulfur emissions data
      add_volc_emis, &! returns .true. when including non-eruptive volcano emissions
      scyc_ncomp,    &! returns number of sulfur cycle components
      scyc_idx1,     &! returns constituent index of 1st sulfur cycle component
      scyc_idx_set    ! set constituent index of 1st sulfur cycle component

   integer, parameter :: ncomp = 4            ! number of sulfur cycle components
   integer      sulfur_idx1                   ! start constituent index for first sulfur cycle component

   logical :: &
      progsul,     &! true => prognostic sulfur (full sulfur cycle)
      geia_emis,   &! true => use the GEIA emission data from M. Barth.
      volc_emis     ! true => include non-eruptive volcano emissions

!##############################################################################
contains
!##############################################################################

   subroutine scycini( xprogsul, emis_type, xvolc_emis )

      !----------------------------------------------------------------------- 
      ! Purpose: Initialize the logical control functions.
      ! 
      ! Author: B. Eaton
      !-----------------------------------------------------------------------

      implicit none

      logical, intent(in) :: &
         xprogsul,     &! true => prognostic sulfur (full sulfur cycle)
         xvolc_emis     ! true => include non-eruptive volcano emissions

      character(len=*) ::&
         emis_type     ! emission dataset type.  Either SMITH or GEIA.  SMITH is
                       ! the default.
      !-----------------------------------------------------------------------

      progsul = xprogsul
      geia_emis = .false.
      if ( trim( emis_type ) == 'GEIA' ) geia_emis = .true.
      volc_emis = xvolc_emis

      if ( progsul ) then
         write(iulog,*)'Sulfur cycle configuration:'
         if ( geia_emis ) then
            write(iulog,*)'  Emissions from 1985 GEIA data prepared by M. Barth.'
         else
            write(iulog,*)'  Emissions from GEIA/SMITH multiyear dataset.'
         end if
         if ( volc_emis ) then
            write(iulog,*)'  Include non-eruptive volcano emissions.'
         else
            write(iulog,*)'  No volcano emissions.'
         end if
      else
         write(iulog,*)'Sulfur cycle disabled.'
      end if

   end subroutine scycini

!#######################################################################

   logical function doscyc ()

      !----------------------------------------------------------------------- 
      ! Purpose: Return .true. when prognostic sulfur cycle is being used.
      ! 
      ! Author: B. Eaton
      !-----------------------------------------------------------------------

      implicit none

      doscyc = progsul

   end function doscyc

!#######################################################################

   logical function useGEIA()

      !----------------------------------------------------------------------- 
      ! Purpose: Return .true. when using the 1985 GEIA sulfur emissions data
      !          prepared by M. Barth.
      ! 
      ! Author: B. Eaton
      !-----------------------------------------------------------------------

      implicit none

      useGEIA = geia_emis

   end function useGEIA

!#######################################################################

   integer function scyc_ncomp()
      implicit none
      scyc_ncomp = ncomp
   end function scyc_ncomp

!#######################################################################

   subroutine scyc_idx_set(m)
     ! store the start index of the sulfur species
     implicit none
     integer m
     write(iulog,*) ' scyc: scyc_idx_set to ', m
     sulfur_idx1 = m
   end subroutine scyc_idx_set

   integer function scyc_idx1()
      implicit none
      scyc_idx1 = sulfur_idx1
   end function scyc_idx1

!#######################################################################

   logical function add_volc_emis()
      implicit none
      add_volc_emis = volc_emis
   end function add_volc_emis

!#######################################################################

end module scyc
