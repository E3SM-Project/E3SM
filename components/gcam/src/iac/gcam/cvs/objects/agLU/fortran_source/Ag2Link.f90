subroutine SETGNP( regionNumber, gnps )

    use Ag2Global8

    integer(4), intent( in ) :: regionNumber
    real(8), intent( in ) :: gnps( 1:NMP )

    gdp( regionNumber + 1, 1:NMP ) = gnps( 1:NMP )
    ! write(*, * ) gdp( regionNumber+1, : )

end subroutine SETGNP

function GETGNP( regionNumber, yearNumber ) result( retGNP )

    use Ag2Global8

    integer, intent( in ) :: regionNumber
    integer, intent( in ) :: yearNumber

    real(8) :: retGNP

    retGNP = gdp( regionNumber + 1, yearNumber + 1 )

end function GETGNP

subroutine SETPOP( regionNumber, pops )

    use Ag2Global8

    integer(4), intent( in ) :: regionNumber
    real(8), intent( in ) :: pops( 1:NMP + 1 )

    ! Fill the AgLU's pop array with the MiniCAM's pop data
    popu( regionNumber + 1, 1:NMP ) = pops( 2:NMP + 1 )

end subroutine SETPOP

    
function GETPOP( regionNumber, yearNumber ) result( retPop )

    use Ag2Global8

    integer, intent( in ) :: regionNumber
    integer, intent( in ) :: yearNumber

    real(8) :: retPop

    retPop = popu( regionNumber + 1, yearNumber )

end function GETPOP

subroutine SETBIOMASSPRICE( biomassPrice )
    
    use Ag2Global8
    
    real(8), intent( in ) :: biomassPrice(:)
    write( *,* ) "Setting biomass price."

    CALP( 1, 8 ) = biomassPrice( 1 ) * GJperGcal * CVRT90

end subroutine SETBIOMASSPRICE

! Do we need to initialize other years?
function GETBIOMASSPRICE( ) result( retBio )

    use Ag2Global8

    real(8) :: retBio

    retBio = CALP( 1, 8 )

end function GETBIOMASSPRICE

! SUBROUTINE AG2INITC( P, GNPinv )
SUBROUTINE AG2INITC( P )
    
    
	! written 7/01 ktg
	! rewritten 1/03 jpl

	use Ag2Global8

	implicit none

	! local arrays: prices to be read into AgLU module (10 is max number of market
	! indicies in AgLU model), global supply, global demand,calibration prices,
	! MiniCAM GDP function which returns YLM (needed for GDP), Land use change emissions
	
	! Arguments
	real(8), intent( inout ) :: P( 12, 14 )
   
    ! There is a weird issue with the 0:NMP in YLM. Unneccessary I believe

	! Variables
	integer(4) :: i, j, t
     real(8) :: MC_CALP(NLP,10)
    P = 0.0_8

	! Initialize AgLU market indicies
	JWood = 1
	JFWood = 2
	JBeef = 3
	JFoodGr = 4
	JCoarseGr = 5
	JOilCrops = 6
	JMiscCrops = 7
	JBio = 8
	JPast = 10	! only prices are set to JPast - pasture sup & dem are INPAST=3
	
	! Initialize miniCAM indices
	INOIL = 1			! Oil
	INGAS = 2			! Gas
	INCOAL = 3			! Coal
	INBMASS = 4 		! Biomass
	INCARB = 5			! Carbon
	INFOREST = 6		! Wood
	INFFOREST = 7		! Forward Wood
	INFOODGR = 8		! Food Grains
	INCOARSEGR = 9		! Coarse Grains
	INOILCROPS = 10 	! Oil Crops
	INMISCCROPS = 11	! Misc Crops
	INPAST = 12 		! Pasture

    MC_CALP = 0.0_8

	! Initialize YLM for Mkt_GDP function within Ag2init
	do i=1, NLP
	  ! Initialize the AgLU price of biomass and convert to 1990$/Gcal
	  MC_CALP(i,JBio) = 23.153_8 ! 2.5_8 * CVRT90 * GJperGCAL
	end do

	! Initialize the AgLU model, passing in population and passing out initial prices
	! Demog module doesn't currently produce forward looking pop estimates
	! which AgLU needs, so the Demog cannot be used with AgLU
	! Also pass the size of the population array
	! call Ag2init( ZLM, MC_CALP, NNLPMax )
    call AG2INIT( MC_CALP )

	! Initialize MiniCAM Ag sector calibration prices
	do i=1, NLP
	  P(INFOREST, i ) = MC_CALP(i,JWood)
	  P(INFFOREST, i ) = MC_CALP(i,JFWood)
	  P(INPAST, i ) = MC_CALP(i,JPast)
	  P(INFOODGR, i ) = MC_CALP(i,JFoodGr)
	  P(INCOARSEGR, i ) = MC_CALP(i,JCoarseGr)
	  P(INOILCROPS, i ) = MC_CALP(i,JOilCrops)
	  P(INMISCCROPS, i ) = MC_CALP(i,JMiscCrops)
	enddo

end subroutine AG2INITC

subroutine AG2RUN( P, Region, Period, AGDEM, AGSUP )
	
	use Ag2Global8
	implicit none
	
	! Arguments
	real(8) , intent( inout ) :: P( NINP ) 
	integer(4), intent( inout ) :: Region, Period
	real(8), intent( inout ) :: AGDEM( NINP ), AGSUP( NINP )

	! Variables
	integer(4)	:: i, in
	real(8) 	:: agPriceLocal( 10, NLP ), Sup( 8, NLP, NMP ), Dem( 8, NLP, NMP )
    
    period = period + 1
    region = region + 1

	! Don't run if period 1 (1975).  Return 1s.
	IF ( period == 1) THEN 
	  AGSUP = 1.0_8
      AGDEM = 1.0_8
	ELSE		 
	  AGSUP = 0.0_8
      AGDEM = 0.0_8
	  agPriceLocal = 0.0_8
      Sup = 0.0_8
      Dem = 0.0_8

	  ! Assign prices from MiniCAM solution to ag model arrays		
	  ! Convert price of biomass first to 1990$/GJ and then to 1990$/Gcal
	  ! Region and Period are passed in.
	  agPriceLocal( JWood, Region ) = P( INFOREST )
	  agPriceLocal( JFWood, Region ) = P( INFFOREST )
	  agPriceLocal( JPast, Region ) = P( INPAST )
	  agPriceLocal( JFoodGr, Region ) = P( INFOODGR )
	  agPriceLocal( JCoarseGr, Region ) = P( INCOARSEGR )
	  agPriceLocal( JOilCrops, Region ) = P( INOILCROPS )
	  agPriceLocal( JMiscCrops, Region ) = P( INMISCCROPS )
	  agPriceLocal( JBio, Region ) = P( INBMASS ) * CVRT90 * GJperGCAL !
       
	! Now adjust biomass price for carbon price. 
	! The market "price" for market INLUCEm is the emissions per unit biomass production in Tg/EJ
	! So conversion is the same as for TXUILM0 (except LU model takes prices in $1990)
	 

	  ! Call the AgLU model, passing in prices, passing out supply and demand
	  call Ag2model( Region, Period, agPriceLocal, Dem, Sup )

	  ! Return supply back to the MiniCAM
	  AGSUP( INFOREST ) = Sup( JWood, Region, Period )
	  AGSUP( INFFOREST ) = Sup( JFWood, Region, Period)
	  AGSUP( INPAST ) = Sup( JBeef, Region, Period)
	  AGSUP( INFOODGR ) = Sup( JFoodGr, Region, Period)
	  AGSUP( INCOARSEGR ) = Sup( JCoarseGr, Region, Period)
	  AGSUP( INOILCROPS ) = Sup( JOilCrops, Region, Period)
	  AGSUP( INMISCCROPS ) = Sup( JMiscCrops, Region, Period)
	  AGSUP( INBMASS ) = Sup( JBio, Region, Period) * GJperGCAL / (100000d0) ! Convert from 10^13cals to EJ


	  ! Return demand back to the MiniCAM (not for biomass)
	  AGDEM( INFOREST  ) = Dem( JWood, Region, Period )
	  AGDEM( INFFOREST ) = Dem( JFWood, Region, Period )
	  AGDEM( INPAST ) = Dem( JBeef, Region, Period )
	  AGDEM( INFOODGR ) = Dem( JFoodGr, Region, Period )
	  AGDEM( INCOARSEGR ) = Dem( JCoarseGr, Region, Period )
	  AGDEM( INOILCROPS ) = Dem( JOilCrops, Region, Period )
	  AGDEM( INMISCCROPS ) = Dem( JMiscCrops, Region, Period )
	END IF

end subroutine AG2RUN

subroutine AG2LINKOUT
    call AG2OUTPUT()
end subroutine AG2LINKOUT

function AG2CO2EMISSIONS( period, region ) result( CARBLAND )
	
	use Ag2Global8
	implicit none
	
	! Arguments   
	integer(4), intent( inout ) :: period
    integer(4), intent( inout ) :: region
	real(8) :: CARBLAND
	
	! Variables

	real(8) :: LUCEmiss( NLP,NMP ) 
	! Body of subroutine
    
    period = period + 1
    region = region + 1

	! Retrieve CO2 emissions and other gas activities at end of each time period
	call Ag2CH4N2O( period, AGCH4ACT( 1:3, 1:NLP ), AGN2OACT( 1:3, 1:NLP ) )
    
    CARBLAND = 0.0_8

	if ( period > 1 ) then
		call Ag2Emiss( period, LUCEmiss )
	    CARBLAND = LUCEmiss( region, period )	! Fill MiniCAM array with LUC emissions
    endif

end function AG2CO2EMISSIONS
