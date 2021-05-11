!=======================================================================
! Copyright 2008, LANSLLC. All rights reserved.
! Unless otherwise indicated, this information has been authored by an 
! employee or employees of the Los Alamos National Security, LLC (LANS), 
! operator of the Los Alamos National Laboratory under Contract No. 
! DE-AC52-06NA25396 with the U.S. Department of Energy. The U.S. Government 
! has rights to use, reproduce, and distribute this information. The public 
! may copy and use this information without charge, provided that this 
! Notice and any statement of authorship are reproduced on all copies. 
! Neither the Government nor LANS makes any warranty, express or implied, 
! or assumes any liability or responsibility for the use of this 
! information.
!
! CICE is developed and maintained by Elizabeth C. Hunke (eclare@lanl.gov)
! and William H. Lipscomb (lipscomb@lanl.gov) of Group T-3 (Fluid 
! Dynamics), Los Alamos National Laboratory, with support from the 
! Climate Change Prediction Program (CCPP) and the Scientific 
! Discovery through Advanced Computing (SciDAC) program of the U.S. 
! Department of Energy.  We thank John Dukowicz (T-3), Phil Jones (T-3), 
! and Robert Malone (CCS-2) for their support of the sea ice modeling 
! effort at LANL.
!   
! CICE has been developed in close collaboration with the NCAR CCSM
! climate modeling project and includes ideas and efforts from 
! members the CCSM Polar Climate Working Group (PCWG).  We especially 
! thank the following members of the PCWG code development team:
!
! Cecilia Bitz, UW
! Bruce Briegleb, NCAR
! Tony Craig, NCAR
! John Dennis, NCAR
! Marika Holland, NCAR
! Bonnie Light, UW
! Julie Schramm, NCAR
! David Bailey, NCAR
!
! Assistance from researchers from U.K. Met Office Hadley Centre is
! also much appreciated, especially our collaboration with Alison 
! McLaren and Ann Keen. Numerous others have contributed to this 
! effort--thanks to all! 
!=======================================================================
!BOP
!
! !ROUTINE: debug_ice - wrapper for print_state
!
! !DESCRIPTION:
!
! Wrapper for the print_state debugging routine.
! Useful for debugging in the main driver (see ice.F_debug)
! ip, jp, mtask are set in ice_diagnostics.F
!
! !REVISION HISTORY:
!
! author Elizabeth C. Hunke, LANL
!
! !INTERFACE:
!
      subroutine debug_ice(plabeld)
!
! !USES:
!
      use ice_kinds_mod
      use ice_diagnostics
      use ice_domain, only: nblocks
      use ice_blocks, only: nx_block, ny_block
!
! !INPUT/OUTPUT PARAMETERS:
!
      character (char_len), intent(in) :: plabeld
!
!EOP
!
      integer (kind=int_kind) :: i, j, iblk

      do iblk = 1, nblocks
      do j = 1, ny_block
      do i = 1, nx_block
         if (iblk==iblkp .and. i==ip .and. j==jp .and. my_task==mtask) &
              call print_state(plabeld,i,j,iblk)
      enddo
      enddo
      enddo

      end subroutine debug_ice

!=======================================================================
