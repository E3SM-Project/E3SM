
      module mo_files_hdr

      use io, only : temp_path, procfiles_path

      contains

      subroutine files_hdr

      implicit none

!-----------------------------------------------------------------------
!        ... The local variables
!-----------------------------------------------------------------------
      integer           :: slen
      logical           :: lexist
      
      inquire( file = 'files.h', exist = lexist )
      if( lexist ) then
         call system( 'rm files.h' )
      end if
      open( unit = 30, file = 'files.h' )

      write(30,'(''#define SETRXTFILE '',a)') TRIM( temp_path ) // 'mo_setrxt.F'
      write(30,'(''#define ADJRXTFILE '',a)') TRIM( temp_path ) // 'mo_adjrxt.F'
      write(30,'(''#define PHTADJFILE '',a)') TRIM( temp_path ) // 'mo_phtadj.F'
      write(30,'(''#define RXTMODFILE '',a)') TRIM( temp_path ) // 'mo_rxt_mod.F'
      write(30,'(''#define GRPVMRFILE '',a)') TRIM( temp_path ) // 'mo_make_grp_vmr.F'
      write(30,'(''#define SETDATFILE '',a)') TRIM( temp_path ) // 'mo_sim_dat.F'
      slen = len_trim( procfiles_path )
      write(30,'(''#define EXPSLVPATH '',a)') procfiles_path(:slen-1)
      write(30,'(''#define IMPSLVPATH '',a)') procfiles_path(:slen-1)
      write(30,'(''#define MODSPATH '',a)') procfiles_path(:slen-1)

      close(30)
      
      end subroutine files_hdr

      end module mo_files_hdr
