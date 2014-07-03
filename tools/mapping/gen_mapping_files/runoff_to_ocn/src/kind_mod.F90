!===============================================================================
! SVN $Id: kind_mod.F90 35698 2012-03-22 23:59:57Z kauff $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/tools/mapping/trunk_tags/mapping_140422a/gen_mapping_files/runoff_to_ocn/src/kind_mod.F90 $
!===============================================================================

MODULE kind_mod

   use shr_kind_mod
   !----------------------------------------------------------------------------
   ! precision/kind constants add data public
   !----------------------------------------------------------------------------
   public
   integer,parameter :: R8 = SHR_KIND_R8  ! 8 byte real
   integer,parameter :: R4 = SHR_KIND_R4  ! 4 byte real
   integer,parameter :: RN = SHR_KIND_RN  ! native real
   integer,parameter :: I8 = SHR_KIND_I8  ! 8 byte integer
   integer,parameter :: I4 = SHR_KIND_I4  ! 4 byte integer
   integer,parameter :: IN = SHR_KIND_IN  ! native integer
   integer,parameter :: CL = SHR_KIND_CL  ! long char
   integer,parameter :: CS = SHR_KIND_CS  ! short char

END MODULE kind_mod
