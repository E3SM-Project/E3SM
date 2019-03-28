SUBROUTINE AG2INIT( CALP2 )

! Initialization subroutine for the AgLU module

! Created by: Kenny Gillingham 7/10/01

USE Ag2Global8

IMPLICIT NONE

real(8), intent( inout ) :: CALP2(NLP,10)

INTEGER T,i

CLIMATE = 1.0_8 ! initialize the climate yield multipliers

! Populate the Ag2 arrays with data
CALL Ag2control

! Bring in GDP Data from MiniCAM to overwrite read-in data
! Could calculate GDP here using MiniCAM's labor prod numbers
! GDP in trillions of 1990 dollars
DO i=1, NLP
  CALP(i,8) = CALP2(i,8) ! Only brought period 1 in, clear the rest. Could be very wrong.
END DO

! Procedure to calculate GDP per capita using MiniCAM data(starting in 1990)
gdpcap(:,2:9) = 1.0d9 * gdp(:,2:9) / popu(:,2:NMP)

! Extend GDP, GDPpercap and Population data for 3 more time steps for fwd forests
CALL Ag2Extend

! Calibrate the AgLU module
CALL Ag2Calibrate

CALP2 = CALP

! Initialize tree yields
T = 2
TreeYield(:,T+1) = TreeYield(:,T) * KJ(2,:,T+1)
TreeYield(:,T+2) = TreeYield(:,T) * KJ(2,:,T+2)

END SUBROUTINE AG2INIT