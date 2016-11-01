! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module atmphys_build_tables_thompson
 use module_mp_thompson

 implicit none
 private
 public:: build_tables_thompson

!builds the files containing the look-up tables for the Thompson cloud microphysics scheme.
!Laura D. Fowler (send comments to laura@ucar.edu).
!2016-11-01.


 contains


!=================================================================================================================
 subroutine build_tables_thompson
!=================================================================================================================

!local variables:
 logical, parameter:: l_mp_tables = .false.

!-----------------------------------------------------------------------------------------------------------------

!--- partial initialization before building the look-up tables:
 call thompson_init(l_mp_tables)

!--- building look-up table for rain collecting graupel:
 write(0,*)
 write(0,*) '--- enter subroutine qr_acr_qg:'
 call qr_acr_qg
 open(unit=11,file='MP_THOMPSON_QRacrQG_DATA.DBL',form='unformatted')
 write(11) tcg_racg
 write(11) tmr_racg
 write(11) tcr_gacr
 write(11) tmg_gacr
 write(11) tnr_racg
 write(11) tnr_gacr
 close(unit=11)
 write(0,*) '--- end subroutine qr_acr_qg.'

!--- building look-up table for rain collecting snow:
 write(0,*)
 write(0,*) '--- enter subroutine qr_acr_qs:'
 call qr_acr_qs
 open(unit=11,file='MP_THOMPSON_QRacrQS_DATA.DBL',form='unformatted',status='new')
 write(11)tcs_racs1
 write(11)tmr_racs1
 write(11)tcs_racs2
 write(11)tmr_racs2
 write(11)tcr_sacr1
 write(11)tms_sacr1
 write(11)tcr_sacr2
 write(11)tms_sacr2
 write(11)tnr_racs1
 write(11)tnr_racs2
 write(11)tnr_sacr1
 write(11)tnr_sacr2
 close(unit=11)
 write(0,*) '--- end subroutine qr_acr_qs'

!--- building look-up table for freezing of cloud droplets:
 write(0,*)
 write(0,*) '--- enter subroutine freezeH2O:'
 call freezeH2O
 open(unit=11,file='MP_THOMPSON_freezeH2O_DATA.DBL',form='unformatted')
 write(11) tpi_qrfz
 write(11) tni_qrfz
 write(11) tpg_qrfz
 write(11) tnr_qrfz
 write(11) tpi_qcfz
 write(11) tni_qcfz
 close(unit=11)
 write(0,*) '--- end subroutine freezeH2O'
 
!--- building look-up table for autoconversion of cloud ice to snow:
 write(0,*)
 write(0,*) '--- enter subroutine qi_aut_qs:'
 call qi_aut_qs
 open(unit=11,file='MP_THOMPSON_QIautQS_DATA.DBL',form='unformatted')
 write(11) tpi_ide
 write(11) tps_iaus
 write(11) tni_iaus
 close(unit=11)
 write(0,*) '--- end subroutine qi_aut_qs'

 end subroutine build_tables_thompson

!=================================================================================================================
 end module atmphys_build_tables_thompson 
!=================================================================================================================
