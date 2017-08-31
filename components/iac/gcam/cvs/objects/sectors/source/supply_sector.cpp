/*
 * LEGAL NOTICE
 * This computer software was prepared by Battelle Memorial Institute,
 * hereinafter the Contractor, under Contract No. DE-AC05-76RL0 1830
 * with the Department of Energy (DOE). NEITHER THE GOVERNMENT NOR THE
 * CONTRACTOR MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY
 * LIABILITY FOR THE USE OF THIS SOFTWARE. This notice including this
 * sentence must appear on any copies of this computer software.
 * 
 * EXPORT CONTROL
 * User agrees that the Software will not be shipped, transferred or
 * exported into any country or used in any manner prohibited by the
 * United States Export Administration Act or any other applicable
 * export laws, restrictions or regulations (collectively the "Export Laws").
 * Export of the Software may require some form of license or other
 * authority from the U.S. Government, and failure to obtain such
 * export control license may result in criminal liability under
 * U.S. laws. In addition, if the Software is identified as export controlled
 * items under the Export Laws, User represents and warrants that User
 * is not a citizen, or otherwise located within, an embargoed nation
 * (including without limitation Iran, Syria, Sudan, Cuba, and North Korea)
 *     and that User is not otherwise prohibited
 * under the Export Laws from receiving the Software.
 * 
 * Copyright 2011 Battelle Memorial Institute.  All Rights Reserved.
 * Distributed as open-source under the terms of the Educational Community 
 * License version 2.0 (ECL 2.0). http://www.opensource.org/licenses/ecl2.php
 * 
 * For further details, see: http://www.globalchange.umd.edu/models/gcam/
 * 
 */

/*! 
* \file supply_sector.cpp
* \ingroup Objects
* \brief SupplySector class source file.
* \author James Blackwood, Sonny Kim
*/

#include "util/base/include/definitions.h"
#include <string>
#include <cassert>

// xml headers
#include <xercesc/dom/DOMNode.hpp>

#include "util/base/include/xml_helper.h"
#include "sectors/include/supply_sector.h"
#include "containers/include/scenario.h"
#include "marketplace/include/marketplace.h"
#include "sectors/include/subsector.h"
#include "util/logger/include/ilogger.h"
#include "marketplace/include/imarket_type.h"
#include "util/base/include/configuration.h"
#include "containers/include/iinfo.h"
#include "sectors/include/sector_utils.h"
#include "sectors/include/cal_quantity_tabulator.h"
#include "reporting/include/indirect_emissions_calculator.h"
#include "util/base/include/summary.h"

using namespace std;
using namespace xercesc;

extern Scenario* scenario;

// static initialize.
const string SupplySector::XML_NAME = "supplysector";

/* \brief Constructor
* \param aRegionName The name of the region.
*/
SupplySector::SupplySector( const string& aRegionName ):
Sector( aRegionName ),mHasTrialSupplyMarket( false ),
mBiomassAdder( scenario->getModeltime()->getmaxper() ),
// The default price for a trial supply market is 0.001
mPriceTrialSupplyMarket( scenario->getModeltime()->getmaxper(), 0.001 )
{
}

/*! \brief Get the XML node name for output to XML.
*
* This public function accesses the private constant string, XML_NAME.
* This way the tag is always consistent for both read-in and output and can be easily changed.
* This function may be virtual to be overridden by derived class pointers.
* \author Josh Lurz, James Blackwood
* \return The constant XML_NAME.
*/
const string& SupplySector::getXMLName() const {
    return XML_NAME;
}

/*! \brief Get the XML node name in static form for comparison when parsing XML.
*
* This public function accesses the private constant string, XML_NAME.
* This way the tag is always consistent for both read-in and output and can be easily changed.
* The "==" operator that is used when parsing, required this second function to return static.
* \note A function cannot be static and virtual.
* \author Josh Lurz, James Blackwood
* \return The constant XML_NAME as a static.
*/
const std::string& SupplySector::getXMLNameStatic() {
    return XML_NAME;
}

/*! \brief Parses any child nodes specific to derived classes
*
* Method parses any input data from child nodes that are specific to the classes derived from this class. Since Sector is the generic base class, there are no values here.
*
* \author Josh Lurz, Steve Smith, Sonny Kim
* \param nodeName name of current node
* \param curr pointer to the current node in the XML input tree
*/
bool SupplySector::XMLDerivedClassParse( const string& nodeName, const DOMNode* curr ) {
    // Temporary hack for CCTP.
    if( nodeName == "biomass-price-adder" ){
        XMLHelper<double>::insertValueIntoVector( curr, mBiomassAdder, scenario->getModeltime() );
    }
    else if( nodeName == "has-trial-supply-market" ){
        mHasTrialSupplyMarket = XMLHelper<bool>::getValue( curr );
    }
    else if( nodeName == "price-trial-supply" ){
        XMLHelper<double>::insertValueIntoVector( curr, mPriceTrialSupplyMarket, scenario->getModeltime() );
    }        
    else {
        return false;
    }
    return true;
}

/*! \brief XML output stream for derived classes
*
* Function writes output due to any variables specific to derived classes to XML.
* This function is called by toInputXML in the base Sector class.
*
* \author Steve Smith, Josh Lurz, Sonny Kim
* \param out reference to the output stream
* \param tabs A tabs object responsible for printing the correct number of tabs. 
*/
void SupplySector::toInputXMLDerived( ostream& aOut, Tabs* aTabs ) const {  

    const Modeltime* modeltime = scenario->getModeltime();
    // Temporary CCTP hack.
    XMLWriteVector( mBiomassAdder, "biomass-price-adder", aOut, aTabs, modeltime, 0.0 );
    if( mHasTrialSupplyMarket ){
        XMLWriteElement( mHasTrialSupplyMarket, "has-trial-supply-market", aOut, aTabs );
    }
    for( int period = 0; period < modeltime->getmaxper(); ++period ) {
        XMLWriteElementCheckDefault( mPriceTrialSupplyMarket[ period ], "price-trial-supply",
                                     aOut, aTabs, 0.001, modeltime->getper_to_yr( period ) );
    }
}

/*! \brief XML debugging output stream for derived classes
*
* Function writes output due to any variables specific to derived classes to XML.
* This function is called by toInputXML in the base Sector class.
*
* \author Steve Smith, Josh Lurz, Sonny Kim
* \param out reference to the output stream
* \param tabs A tabs object responsible for printing the correct number of tabs. 
*/
void SupplySector::toDebugXMLDerived( const int aPeriod, ostream& aOut, Tabs* aTabs ) const {

    XMLWriteElement( mBiomassAdder[ aPeriod ], "biomass-price-adder", aOut, aTabs );
    if( mHasTrialSupplyMarket ){
        XMLWriteElement( mHasTrialSupplyMarket, "has-trial-supply-market", aOut, aTabs );
        XMLWriteElement( mPriceTrialSupplyMarket[ aPeriod ], "price-trial-supply", aOut, aTabs );
    }
}

/*! \brief Complete the initialization of the supply sector.
* \param aRegionInfo Regional information object.
* \param aDependencyFinder Regional dependency finder.
*/
void SupplySector::completeInit( const IInfo* aRegionInfo,
                                 DependencyFinder* aDependencyFinder,
                                 ILandAllocator* aLandAllocator )
{
	// default unit to EJ
	if ( mOutputUnit.empty() ) {
		mOutputUnit = "EJ"; 
	}
	// default unit to EJ
	if ( mInputUnit.empty() ) {
		mInputUnit = "EJ"; 
	}
	// default unit to $/GJ
	if ( mPriceUnit.empty() ) {
		mPriceUnit = "75$/GJ"; 
	}
	Sector::completeInit( aRegionInfo, aDependencyFinder, aLandAllocator );	
    setMarket();
}

/*! \brief Create new market for this Sector
*
* Sets up the appropriate market within the marketplace for this Sector. Note that the type of market is NORMAL -- 
* signifying that this market is a normal market that is solved (if necessary).
*
* \author Sonny Kim, Josh Lurz, Steve Smith
*/
void SupplySector::setMarket() {    
    Marketplace* marketplace = scenario->getMarketplace();
    // Creates a regional market. MiniCAM supply sectors are not independent and 
    // cannot be members of multi-region markets.
    if( marketplace->createMarket( regionName, regionName, name, IMarketType::NORMAL ) ) {
        // Initialize prices for markets
        marketplace->setPriceVector( name, regionName, mPrice );

        // Set price and output units for period 0 market info
        IInfo* marketInfo = marketplace->getMarketInfo( name, regionName, 0, true );
        marketInfo->setString( "price-unit", mPriceUnit );
        marketInfo->setString( "output-unit", mOutputUnit );
    }
    // Create trial supply market.
    if( mHasTrialSupplyMarket ){
        if( SectorUtils::createTrialSupplyMarket(regionName,name,
            marketplace->getMarketInfo(name,regionName,0,true)) ){
            // Initialize "prices" (quantities) for trial supply market
            marketplace->setPriceVector( SectorUtils::getTrialMarketName(name),
                regionName, mPriceTrialSupplyMarket );
        }
    }
    else {
        // If this sector does not explicitly create a trial supply it may still
        // have had trial demands because the dependency finder choose to make one
        // for this sector.  In that case we will store the initial trial demand
        // in the market info and if the dependency finder decides to make a trial
        // market again for this sector it can utilize this value.
        const Modeltime* modeltime = scenario->getModeltime();
        for( int period = 1; period < modeltime->getmaxper(); ++period ) {
            if( mPriceTrialSupplyMarket[ period ] != 0.001 ) {
                marketplace->getMarketInfo( name, regionName, period, true )
                    ->setDouble( "initial-trial-demand", mPriceTrialSupplyMarket[ period ] );
            }
        }
    }
}

/*! \brief Initialize the SupplySector.
* \details Currently only calls the base class initCalc.
* \param aNationalAccount National accounts container.
* \param aDemographics Regional demographics object.
* \param aPeriod Period for which to initialize the SupplySector.
*/
void SupplySector::initCalc( NationalAccount* aNationalAccount,
                            const Demographic* aDemographics,
                            const int aPeriod )
{
    Sector::initCalc( aNationalAccount, aDemographics, aPeriod );

    // Check if the sector should create a trial supply market or energy final
    // demand supply object. First check if the flag is already set. This is
    // only done in period 1 so that other markets have a chance to set the
    // flag.
    if( aPeriod == 1 ){
        if( SectorUtils::isFinalEnergySector( regionName, name ) ){
            mFinalEnergySupplier.reset( new FinalEnergySupplier( name ) );
        }
    }
}

/*! \brief returns Sector output.
*
* Returns the total amount of the SupplySector. 
*
* \author Sonny Kim
* \param period Model period
* \todo make year 1975 regular model year so that logic below can be removed
* \return total output
*/
double SupplySector::getOutput( const int aPeriod ) const {
    double output = 0;
    for ( unsigned int i = 0; i < subsec.size(); ++i ) {
        double subsecOutput = subsec[ i ]->getOutput( aPeriod );
        // error check.
        if ( !util::isValidNumber( subsecOutput ) ){
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::ERROR );
            mainLog << "Output for subsector " << subsec[ i ]->getName() << " in Sector " << name 
                << " in region " << regionName <<" is not valid." << endl;
            continue;
        }
        output += subsecOutput;
    }

    // In the base period return a read in output if there is none.
    if( aPeriod == 0 && output == 0 ){
        return mBaseOutput;
    }

    return output;
}

/*! \brief Return the price of the SupplySector.
* \details The price of a SupplySector is the weighted average subsector price.
* \param aPeriod Model period.
* \return Price.
* \todo Move entire calculation here once demand sectors are rewritten.
*/
double SupplySector::getPrice( const GDP* aGDP, const int aPeriod ) const {
    return Sector::getPrice( aGDP, aPeriod );
}

/*! \brief Calculate the final supply price.
* \details Calculates shares for the sector and price for the supply sector, and
*          then sets the price of the good into the marketplace.
* \param aGDP The regional GDP container.
* \param aPeriod The period in which to calculate the final supply price.
*/
void SupplySector::calcFinalSupplyPrice( const GDP* aGDP, const int aPeriod ){
    // Instruct all subsectors to calculate their costs. This must be done
    // before prices can be calculated.
    calcCosts( aPeriod );

    // Set the price into the market.
    Marketplace* marketplace = scenario->getMarketplace();

    double avgMarginalPrice = getPrice( aGDP, aPeriod );

    // Temporary hack for CCTP.
    if( name == "regional biomass" ){
        // Adjust for biomass carbon value adder.
        // Prices must be adjusted by the price multiplier if there is a
        // carbon tax in place.
        double carbonPrice = marketplace->getPrice( "CO2", regionName, aPeriod, false );
        // Retrieve proportional tax rate.
        const IInfo* marketInfo = marketplace->getMarketInfo( "CO2", regionName, aPeriod, false );
        // Note: the key includes the region name.
        const double proportionalTaxRate = 
            ( marketInfo && marketInfo->hasValue( "proportional-tax-rate" + regionName ) ) 
            ? marketInfo->getDouble( "proportional-tax-rate" + regionName, true )
            : 1.0;
        if( carbonPrice != Marketplace::NO_MARKET_PRICE && carbonPrice > 0 ){
            // Adjust carbon price with the proportional tax rate.
            carbonPrice *= proportionalTaxRate;

            // The carbon price and all crops except biomass are in 1990
            // dollars. Biomass must be adjusted by a 1975 carbon price.
            const double CONVERT_90_TO_75 = 2.212;
            
            avgMarginalPrice += mBiomassAdder[ aPeriod ] * carbonPrice / CONVERT_90_TO_75;
        }
    }
    
    marketplace->setPrice( name, regionName, avgMarginalPrice, aPeriod, true );
}

/*! \brief Set supply Sector output
* \details This routine takes the market demand and propagates that through the
*          supply subsectors where it is shared out (and subsequently passed to
*          the technology level within each sub-Sector to be shared out).
* \author Sonny Kim
* \param aGDP GDP object uses to calculate various types of GDPs.
* \param aPeriod Model period
*/
void SupplySector::supply( const GDP* aGDP, const int aPeriod ) {
	Marketplace* marketplace = scenario->getMarketplace();
	// demand for the good produced by this Sector
	double marketDemand = marketplace->getDemand( name, regionName, aPeriod );

	// Determine if fixed output must be scaled because fixed supply
	// exceeded demand.
	double fixedOutput = getFixedOutput( aPeriod );
	double scaleFactor = SectorUtils::calcFixedOutputScaleFactor( marketDemand, fixedOutput );

	// Calculate the demand for new investment.
	double newInvestment = max( marketDemand - fixedOutput, 0.0 );
	const vector<double> subsecShares = calcSubsectorShares( aGDP, aPeriod );

	// This is where subsector and technology outputs are set
	for( unsigned int i = 0; i < subsec.size(); ++i ){
		// set subsector output from Sector demand
		subsec[ i ]->setOutput( subsecShares[ i ] * newInvestment, scaleFactor, aGDP, aPeriod );
	}    

	// Set the final energy for the calibration market.
    if( mFinalEnergySupplier.get() ){
        mFinalEnergySupplier->setFinalEnergy( regionName,
                                              getEnergyInput( aPeriod ),
                                              aPeriod );
    }

    // Add demand to trial supply market.
    if( mHasTrialSupplyMarket ){
        SectorUtils::addToTrialDemand( regionName, name, marketDemand, aPeriod );
    }

	const static bool debugChecking = Configuration::getInstance()->getBool( "debugChecking" );
	if ( debugChecking ) {
		// If the model is working correctly this should never give an error
		// An error here means that the supply summed up from the supply sectors 
		// is not equal to the demand that was passed in 
		double mrksupply = getOutput( aPeriod );

		// if demand identically = 1 then must be in initial iteration so is not an error
		if ( aPeriod > 0 && fabs(mrksupply - marketDemand ) > 0.01 && marketDemand != 1 ) {
			ILogger& mainLog = ILogger::getLogger( "main_log" );
			mainLog.setLevel( ILogger::WARNING );
			mainLog << regionName << " Market "<<  name << " demand and derived supply are not equal by: ";
			mainLog << fabs( mrksupply - marketDemand ) << ": ";
			mainLog << "S: " << mrksupply << " D: " << marketDemand << " Fixed-Supply: " << getFixedOutput( aPeriod ) << endl;
		}
	}
}

/*!
 * \brief Get the energy input for the SupplySector.
 * \todo If there is a DemandSupplySector, move this there.
 * \param aPeriod Period.
 * \return Total energy input.
 */
double SupplySector::getEnergyInput( const int aPeriod ) const {
    double totalEnergy = 0;
    for( unsigned int i = 0; i < subsec.size(); ++i ){
        totalEnergy += subsec[ i ]->getEnergyInput( aPeriod );
    }
    return totalEnergy;
}

/*!
 * \brief Constructor.
 * \param aSectorName Name of the parent sector.
 */
SupplySector::FinalEnergySupplier::FinalEnergySupplier( const string& aSectorName ) {
    mTFEMarketName = SectorUtils::createTFEMarketName( aSectorName );
}

/*!
 * \brief Set the quantity of final energy into the final energy market for the
 *        final demand sector.
 * \param aRegionName Region name.
 * \param aFinalEnergy Quantity of final energy.
 * \param aPeriod Period.
 */
void SupplySector::FinalEnergySupplier::setFinalEnergy( const string& aRegionName,
                                                        const double aFinalEnergy,
                                                        const int aPeriod )
{
	if( aPeriod > 1 ){
        Marketplace* marketplace = scenario->getMarketplace();
		marketplace->addToDemand( mTFEMarketName, aRegionName,
                                  aFinalEnergy, aPeriod, false );
	}
}

/*! \brief Function to finalize objects after a period is solved.
* \details This function is used to calculate and store variables which are only needed after the current
* period is complete.
* \param aPeriod The period to finalize.
* \todo Finish this function.
* \author Josh Lurz, Sonny Kim
*/
void SupplySector::postCalc( const int aPeriod ){
    Sector::postCalc( aPeriod );
    const Marketplace* marketplace = scenario->getMarketplace();
    // If trial supply market exists, get solved trial "prices" and set to member
    // price vector.
    if( aPeriod > 0 && mHasTrialSupplyMarket ){
        mPriceTrialSupplyMarket[ aPeriod ] = marketplace->getPrice(
            SectorUtils::getTrialMarketName(name),
            regionName, aPeriod, true );
    }
    else if( aPeriod > 0 && marketplace->getMarketInfo( name, regionName, aPeriod, true )
             ->getBoolean( "has-split-market", false ) )
    {
        // We may still have trial "prices" if the dependency finder generated a
        // trial market for this sector.
        mPriceTrialSupplyMarket[ aPeriod ] = marketplace->getDemand( name, regionName, aPeriod );
    }
}

//! Write MiniCAM style Sector output to database.
void SupplySector::dbOutput( const GDP* aGDP,
                             const IndirectEmissionsCalculator* aIndEmissCalc ) const
{
    const Modeltime* modeltime = scenario->getModeltime();
    // function protocol
    void dboutput4(string var1name,string var2name,string var3name,string var4name,
        string uname,vector<double> dout);

    // total Sector output
    int maxper = modeltime->getmaxper();
    vector<double> temp(maxper);
    for( int per = 0; per < maxper; ++per ){
        temp[ per ] = getOutput( per );
    }
    dboutput4( regionName,"Secondary Energy Prod","by Sector",name,mOutputUnit, temp );
    dboutput4( regionName,"Secondary Energy Prod",name,"zTotal",mOutputUnit, temp );


    string str; // temporary string

    // Sector fuel consumption by fuel type
    typedef map<string,double>:: const_iterator CI;
    map<string,double> tfuelmap = summary[0].getfuelcons();
    for (CI fmap=tfuelmap.begin(); fmap!=tfuelmap.end(); ++fmap) {
        for (int m=0;m<maxper;m++) {
            temp[m] = summary[m].get_fmap_second(fmap->first);
        }
        if( fmap->first == "" ){
            dboutput4( regionName,"Fuel Consumption",name, "No Fuelname", mInputUnit,temp);
        }
        else {
            dboutput4( regionName,"Fuel Consumption",name,fmap->first,mInputUnit,temp);
        }
    }

    // Sector emissions for all greenhouse gases
    map<string,double> temissmap = summary[0].getemission(); // get gases for per 0
    for (CI gmap=temissmap.begin(); gmap!=temissmap.end(); ++gmap) {
        for (int m=0;m<maxper;m++) {
            temp[m] = summary[m].get_emissmap_second(gmap->first);
        }
        dboutput4(regionName,"Emissions","Sec-"+name,gmap->first,"MTC",temp);
    }
    // CO2 emissions by Sector
    for ( int m=0;m<maxper;m++) {
        temp[m] = summary[m].get_emissmap_second("CO2");
    }
    dboutput4( regionName,"CO2 Emiss","by Sector",name,"MTC",temp);
    dboutput4( regionName,"CO2 Emiss",name,"zTotal","MTC",temp);

    // CO2 indirect emissions by Sector
    for ( int m=0;m<maxper;m++) {
        temp[m] = aIndEmissCalc->getIndirectEmissions( name, m );
    }
    dboutput4( regionName,"CO2 Emiss(ind)",name,"zTotal","MTC",temp);

    // Sector price
    for ( int m=0;m<maxper;m++) {
        temp[m] = getPrice( aGDP, m );
    }
    dboutput4( regionName,"Price",name,"zSectorAvg",mPriceUnit, temp );
    // for electricity Sector only
    if (name == "electricity") {
        for ( int m=0;m<maxper;m++) {
            temp[m] = getPrice( aGDP, m ) * 2.212 * 0.36;
        }
        dboutput4( regionName,"Price","electricity C/kWh","zSectorAvg","90C/kWh",temp);
    }

    // Sector price
    for ( int m = 0; m < maxper; m++ ) {
        temp[m] = getPrice( aGDP, m );
    }
    dboutput4( regionName,"Price","by Sector",name,mPriceUnit, temp );

    // do for all sub sectors in the Sector
    for( int m = 0; m < maxper; m++ ) {
        temp[ m ] = getOutput( m );
    }

    for( unsigned int i = 0; i < subsec.size(); ++i ){
        // output or demand for each technology
        subsec[ i ]->MCoutputSupplySector( aGDP );
        subsec[ i ]->MCoutputAllSectors( aGDP, aIndEmissCalc, temp );
    }

    // do for all sub sectors in the Sector
    for( unsigned int i = 0; i < subsec.size(); ++i ){
        // output or demand for each technology
        subsec[ i ]->csvOutputFile( aGDP, aIndEmissCalc );
    }
}
