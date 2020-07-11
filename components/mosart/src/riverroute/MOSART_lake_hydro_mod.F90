!
MODULE MOSART_lake_hydro_mod
! Description: core code of lake module in MOSART framework
! A simple wier formula (Kindsvater and Carter (1959))is used to estimate outflow from lakes 
! Developed by Wondie Yigzaw, July 2020. 
! REVISION HISTORY:
! 
!-----------------------------------------------------------------------
    
! !USES:
    use shr_kind_mod  , only : r8 => shr_kind_r8
    use RunoffMod , only : Tctl, TUnit, TRunoff, THeat, TPara, rtmCTL
    use rof_cpl_indices, only : nt_nliq, nt_nice
	use RtmVar         , only : iulog, ngeom, nlayers, rstraflag, lakeflag
	use MOSART_heat_mod
	use MOSART_lake_mod
	use RtmTimeManager
    implicit none
	
    public mosart_lake_hydro
	
! !PUBLIC MEMBER FUNCTIONS:
    contains  
		
	subroutine mosart_lake_hydro(iunit,nt,localDeltaT)
	! !DESCRIPTION: calculate the water temperature of reservoir.
	
		use shr_sys_mod , only : shr_sys_flush
	
		implicit none
		integer,  intent(in) :: iunit, nt
		real(r8), intent(in) :: localDeltaT
		character(len=*),parameter :: subname = '(mosart_lake_hydro)'
		
        integer :: j,i                   ! indices
		real(r8) :: rsd_time = 2592000   ! residence time (s)
		real(r8) :: inflow, outflow      ! Lake inflow, outflow (m3/s)
		
		
!**************************************************************************************************************************************************
		
		if (TRunoff%lake_flg(iunit) ==1) then	! Lake module active if there is natural lake
			do i=2,nlayers
               if (TRunoff%v_zt(iunit,i)<TRunoff%v_zt(iunit,i-1) .or. TRunoff%d_v(iunit,j)<0.0_r8)return
			end do
			do j = 2, ngeom+1	 
               if (TRunoff%a_di(iunit,j)<0._r8 .or. TRunoff%v_zti(iunit,j)<0._r8) return
			end do
			
			TRunoff%lake_inflow(iunit) = -TRunoff%erout(iunit,nt)
			
			if (TRunoff%V_str(iunit)/TRunoff%V_max(iunit)>= 0.75_r8 .and. TRunoff%V_str(iunit)/TRunoff%V_max(iunit)<= 1._r8) then
			    TRunoff%lake_outflow(iunit) = - (TRunoff%V_max(iunit) - TRunoff%V_str(iunit))*1e6/rsd_time  ! Storage converted from mcm to m3
			elseif (TRunoff%V_str(iunit)/TRunoff%V_max(iunit)>1._r8) then
			    TRunoff%lake_outflow(iunit) = - (TRunoff%V_str(iunit) - TRunoff%V_max(iunit))*1e6/rsd_time - TRunoff%lake_inflow(iunit) ! Storage converted from mcm to m3
			else
			    TRunoff%lake_outflow(iunit) = 0._r8
			end if
			TRunoff%erout(iunit,nt) = TRunoff%lake_outflow(iunit) ! now the lake outflow is discharged to downstream as the outflow from the current grid
			
			TRunoff%dV_str(iunit) = TRunoff%lake_inflow(iunit) + TRunoff%lake_outflow(iunit)
			
		end if ! lake hydro
		
    end subroutine mosart_lake_hydro
  
end MODULE MOSART_lake_hydro_mod