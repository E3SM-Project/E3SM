!#**********************************************************************************  
! This computer software was prepared by Battelle Memorial Institute, hereinafter
! the Contractor, under Contract No. DE-AC05-76RL0 1830 with the Department of 
! Energy (DOE). NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY,
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
! miscellaneous debuging routines for CBMZ and MOSAIC
!**********************************************************************************  
	module module_ecpp_util
!==Guangxing Lin
        !use abortutils, only: endrun
        use cam_abortutils, only: endrun
!==Guangxing Lin

	contains

!-----------------------------------------------------------------------
	subroutine ecpp_debugmsg( lun, level, str )
!
! when lun >  0, writes "str" to unit "lun"
! when lun <= 0, passes "str" on to wrf_debug
!
	implicit none
! subr arguments
	integer, intent(in) :: lun, level
	character(len=*), intent(in) :: str
! local variables
	integer n

	n = max( 1, len_trim(str) )
	if (lun .ge. 0) then
	    write(lun,'(a)') str(1:n)
	else
	    call endrun( str(1:n) )
	end if
	return
	end subroutine ecpp_debugmsg


!-----------------------------------------------------------------------
	subroutine ecpp_message( lun, str )
!
! when lun >  0, writes "str" to unit "lun"
! when lun <= 0, passes "str" on to wrf_message
!
	implicit none
! subr arguments
	integer, intent(in) :: lun
	character(len=*), intent(in) :: str
! local variables
	integer n

	n = max( 1, len_trim(str) )
	if (lun .ge. 0) then
	    write(lun,'(a)') str(1:n)
	else
	    call endrun( str(1:n) )
	end if
	return
	end subroutine ecpp_message


!-----------------------------------------------------------------------
	subroutine ecpp_error_fatal( lun, str )
!
! when lun >  0, writes "str" to unit "lun"
! then (always)  passes "str" on to wrf_error_fatal
!
	implicit none
! subr arguments
	integer, intent(in) :: lun
	character(len=*), intent(in) :: str
! local variables
	integer n

	n = max( 1, len_trim(str) )
        call  endrun( str(1:n) )
	return
	end subroutine ecpp_error_fatal


!-----------------------------------------------------------------------
	subroutine parampollu_1clm_set_opts(               &
		xppopt_updn_prof_aa,                       &
		xppopt_quiescn_mf, xppopt_quiescn_sosi,    &
		xppopt_chemtend_wq, xppopt_chemtend_dtsub, &
		xppopt_chemtend_updnfreq                   )

	use module_data_ecpp1

	implicit none


!   subr arguments
	integer, intent(in) ::   &
		xppopt_updn_prof_aa,   &
		xppopt_quiescn_mf, xppopt_quiescn_sosi,   &
		xppopt_chemtend_wq, xppopt_chemtend_dtsub,   &
		xppopt_chemtend_updnfreq
     

	ppopt_updn_prof_aa      = xppopt_updn_prof_aa
	ppopt_quiescn_mf        = xppopt_quiescn_mf
	ppopt_quiescn_sosi      = xppopt_quiescn_sosi
	ppopt_chemtend_wq       = xppopt_chemtend_wq
	ppopt_chemtend_dtsub    = xppopt_chemtend_dtsub
	ppopt_chemtend_updnfreq = xppopt_chemtend_updnfreq


	return
	end subroutine parampollu_1clm_set_opts

!-----------------------------------------------------------------------
	end module module_ecpp_util
