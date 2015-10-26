!
MODULE WRM_read_print
! Description: module to read and write, update from previous version in WRM_subw_IO.F90
! Developed by Nathalie Voisin 12/17/2013
! REVISION HISTORY:
!-----------------------------------------------------------------------

! !USES:
  use shr_kind_mod  , only : r8 => shr_kind_r8, SHR_KIND_CL
  use shr_const_mod , only : SHR_CONST_REARTH, SHR_CONST_PI
  use RunoffMod, only : Trunoff, Tctl, Tunit
  use rof_cpl_indices, only : nt_rtm
  use WRM_type_mod, only : ctlSubwWRM, WRMUnit, StorWater
  use WRM_start_op_year
  
  implicit none
  private

!-----------------------------------------------------------------------
  contains
!-----------------------------------------------------------------------
!______________________________________________________________________________________________________  
  subroutine read_file_grid(FileName, value)
     implicit none
     character(len=250), intent(in) :: FileName
     integer :: iunit 
     integer :: ios
     real(r8), intent(out) :: value(:)

     open (unit=1, file=FileName, status="old", action="read", iostat=ios)
     if(ios /= 0) then
        print*, "Cannot find file ", FileName
        stop
     end if
     iunit=0
     do iunit=1, Tctl%NUnit
!tcx        read(unit=1, *) value(iunit)
        if ( value(iunit) < 0._r8 ) then
           print*, "negative value, iunit ", iunit, FileName
           value(iunit) = 0._r8
        end if
     end do
     close(unit=1)
  end subroutine read_file_grid 
!_______________________________________________________________________________________________________
  subroutine readPotentialEvap(theTime)
  !! DESCRIPTION: read the simulated potential evaporation to adjust the storage of reservoir
     implicit none
     character(len=*), intent(in) :: theTime
     character(len=250) :: petFileName
     integer :: ios, iunit, ilat, ilon    ! flag of IO status, local indices
     real(r8) :: ftemp1            ! tempory array

     petFileName = adjustl(trim(Tctl%runoffPath))//'pet/'//theTime//'.pet'
     CALL read_file_grid(petFileName,StorWater%pot_evap)
     StorWater%pot_evap = 0.75_r8 * StorWater%pot_evap * 0.001_r8/Tctl%DATAH  !mm/day-->m/s, or mm/hr-->m/s, TO DO
     ! the 0.75 is due to the fact that raw potential evap overestimate the evaporation from large bodies of water
     ! USGS assume between 0.65 to 0.85 woith 0.7 when air temp = water temp
  end subroutine readPotentialEvap
!______________________________________________________________________________________________________
  
  subroutine readDemand(theTime)
  ! !DESCRIPTION: read in the irrigation demand data for each time step
     implicit none
     character(len=*), intent(in) :: theTime
     character(len=250) :: demFileName  ! water demand file names
     integer :: ios, iunit, ilat, ilon    ! flag of IO status, local indices
     real(r8) :: ftemp1            ! tempory array

     if ( ctlSubwWRM%GroundwaterFlag > 0 ) then
        if ( ctlSubwWRM%TotalDemandFlag > 0 ) then
           demFileName = adjustl(trim(ctlSubwWRM%demandPath))//trim(Tctl%baseName)//'.gw_nonirr.txt'
           CALL read_file_grid(demFileName,StorWater%GWShareNonIrrig)
        endif
        demFileName = adjustl(trim(ctlSubwWRM%demandPath))//trim(Tctl%baseName)//'.gw_irr.txt'
        CALL read_file_grid(demFileName,StorWater%GWShareIrrig)
     endif

     demFileName = adjustl(trim(ctlSubwWRM%demandPath))//trim(Tctl%baseName)//'_'//theTime//'.ConIrrig'
     CALL read_file_grid(demFileName,StorWater%ConDemIrrig)
     StorWater%demand = StorWater%ConDemIrrig * (1._r8 - StorWater%GWShareIrrig)

     ! Toal demand means differentiation between irrigation and non irrigation demand
     if ( ctlSubwWRM%TotalDemandFlag > 0 ) then
        demFileName = adjustl(trim(ctlSubwWRM%demandPath))//trim(Tctl%baseName)//'_'//theTime//'.ConNonIrrig'
        CALL read_file_grid(demFileName,StorWater%ConDemNonIrrig) 
        StorWater%demand = StorWater%ConDemIrrig*(1._r8-StorWater%GWShareIrrig) + StorWater%ConDemNonIrrig*(1._r8-StorWater%GWShareNonIrrig)
     endif
             
     ! Return flow option means difference between consumptive use and withdrawal
     if ( ctlSubwWRM%ReturnFlowFlag > 0 ) then
        demFileName = adjustl(trim(ctlSubwWRM%demandPath))//trim(Tctl%baseName)//'_'//theTime//'.WithIrrig'
        CALL read_file_grid(demFileName,StorWater%WithDemIrrig)
        StorWater%demand = StorWater%WithDemIrrig*(1._r8-StorWater%GWShareIrrig)

        if ( ctlSubwWRM%TotalDemandFlag > 0 ) then
           demFileName = adjustl(trim(ctlSubwWRM%demandPath))//trim(Tctl%baseName)//'_'//theTime//'.WithNonIrrig'
           CALL read_file_grid(demFileName,StorWater%WithDemNonIrrig)
           StorWater%demand = StorWater%WithDemIrrig*(1._r8-StorWater%GWShareIrrig) + StorWater%WithDemNonIrrig*(1._r8-StorWater%GWShareNonIrrig)
        endif
     endif

  end subroutine readDemand
!________________________________________________________________________________________________________________
  subroutine print_grid_file(theTime, strLen, nio, value)
     ! !DESCRIPTION: output the simulation results into external files
     implicit none
     character(len=*), intent(in) :: theTime  ! the time step to output
     integer, intent(in) :: strLen     ! length of each line to print
     integer, intent(in) :: nio        ! unit of the file to print
     real(r8), intent(in) :: value(:)
     integer :: ios                    ! flag of io status
     integer :: iunit          ! local index
     character(len=strLen) :: strLine
     character(len=20) :: stemp

     strLine = ''
     do iunit=1, Tctl%NUnit
        stemp = ''
        call num2str(value(iunit), stemp, 'e20.10')
        strLine = trim(strLine)//adjustr(stemp)
     end do
     !print*, "StorWater%demand(137)", StorWater%demand(137)
     write(unit=nio,fmt="((a10), (a))") theTime, strLine

  end subroutine print_grid_file
!________________________________________________________________________________________________________________

  subroutine print_dam_file(theTime, strLen, nio, value)
     ! !DESCRIPTION: output the simulation results into external files
     implicit none
     character(len=*), intent(in) :: theTime  ! the time step to output
     integer, intent(in) :: strLen     ! length of each line to print
     integer, intent(in) :: nio        ! unit of the file to print
     real(r8), intent(in) :: value(:)
     integer :: ios                    ! flag of io status
     integer :: iunit          ! local index
     character(len=strLen) :: strLine
     character(len=20) :: stemp

     strLine = ''
     do iunit=1, ctlSubwWRM%NDam
        stemp = ''
        call num2str(value(iunit), stemp, 'e20.10')
        strLine = trim(strLine)//adjustr(stemp)
     end do
     write(unit=nio,fmt="((a10), (a))") theTime, strLine
     !write(unit=nio,fmt="((a10)") "supply"

  end subroutine print_dam_file
!________________________________________________________________________________________________________________

  subroutine printTest2(theTime, theUnit, nio)
     ! !DESCRIPTION: output the simulation results into external files
     implicit none
     character(len=*), intent(in) :: theTime  ! the time step to output
     integer, intent(in) :: theUnit    ! the index of the unit to print 
     integer, intent(in) :: nio        ! unit of the file to print
        
     integer :: ios,ii,nt                    ! flag of io status

     ii=theUnit !167 !191 !

     do nt = 1,nt_rtm
        !write(unit=nio,fmt="((a10),(e20.11))") theTime, Trunoff%flow(ii)
        write(unit=nio,fmt="((a10),6(e20.11))") theTime, nt, Trunoff%qsur(ii,nt), Trunoff%qsub(ii,nt), Trunoff%etin(ii,nt)/(TUnit%area(ii)*TUnit%frac(ii)), Trunoff%erlateral(ii,nt), Trunoff%erin(ii,nt), Trunoff%flow(ii,nt)
     enddo

  end subroutine printTest2   
!________________________________________________________________________________________________________________

end MODULE WRM_read_print
