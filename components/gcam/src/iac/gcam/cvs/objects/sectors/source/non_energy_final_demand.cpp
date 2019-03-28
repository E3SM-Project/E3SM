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
* \file non_energy_final_demand.cpp
* \ingroup Objects
* \brief NonNonEnergyFinalDemand class source file.
* \author Josh Lurz
*/

#include "util/base/include/definitions.h"
#include <string>

#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>

#include "util/base/include/xml_helper.h"
#include "sectors/include/non_energy_final_demand.h"
#include "util/base/include/model_time.h"
#include "marketplace/include/marketplace.h"
#include "containers/include/scenario.h"
#include "containers/include/gdp.h"
#include "util/base/include/configuration.h"
#include "util/base/include/ivisitor.h"
#include "containers/include/iinfo.h"
#include "sectors/include/sector_utils.h"

using namespace std;
using namespace xercesc;

extern Scenario* scenario;

/*! \brief Constructor.
* \author Sonny Kim, Steve Smith, Josh Lurz
*/
NonEnergyFinalDemand::NonEnergyFinalDemand(){
    const Modeltime* modeltime = scenario->getModeltime();
    const int maxper = modeltime->getmaxper();
    
    mIsPerCapitaBased = false;
    // mIncomeElasticity = 0;
    mIncomeElasticity.resize( maxper );

    // mPriceElasticity = 0;
    mPriceElasticity.resize( maxper );

    // resize vectors

    mTechnicalChange.resize( maxper );
    mServiceDemands.resize( maxper );
    mAEEI.resize( maxper );
}

/*! \brief Destructor.
*/
NonEnergyFinalDemand::~NonEnergyFinalDemand(){
}

//! Complete the initialization of the object.
void NonEnergyFinalDemand::completeInit( const string& aRegionName, const IInfo* aRegionInfo ) {
    if( mBaseService <= 0 ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::WARNING );
        mainLog << "Zero base service for demand sector " << mName << " in region " << aRegionName << "." << endl;
        mBaseService = 1;
    }
    // Copy the base service to the base year of the service demand output.
    mServiceDemands[ 0 ] = mBaseService;
}

//! init calc
void NonEnergyFinalDemand::initCalc( const string& aRegionName, const GDP* aGDP, const int aPeriod ){
    calcTechChange( aPeriod );
}

/*! \brief Set data members from XML input
*
* \author Josh Lurz
* \param aNode pointer to the current node in the XML input tree
*/
bool NonEnergyFinalDemand::XMLParse( const DOMNode* aNode ) {
    /*! \pre make sure we were passed a valid node. */
    assert( aNode );

    // get the name attribute.
    mName = XMLHelper<string>::getAttr( aNode, "name" );

    // get all child nodes.
    DOMNodeList* nodeList = aNode->getChildNodes();
    const Modeltime* modeltime = scenario->getModeltime();

    // loop through the child nodes.
    for( unsigned int i = 0; i < nodeList->getLength(); ++i ){
        DOMNode* curr = nodeList->item( i );
        const string nodeName = XMLHelper<string>::safeTranscode( curr->getNodeName() );

        if( nodeName == "#text" ) {
            continue;
        }
        else if( nodeName == "price-elasticity" ) {
            XMLHelper<double>::insertValueIntoVector( curr, mPriceElasticity, modeltime );
            // mPriceElasticity = XMLHelper<double>::getValue( curr );
        }
        else if( nodeName == "service-output" ){
            mBaseService = XMLHelper<double>::getValue( curr );
        }
        else if( nodeName == "income-elasticity" ){
            // mIncomeElasticity = XMLHelper<double>::getValue( curr );
            XMLHelper<double>::insertValueIntoVector( curr, mIncomeElasticity, modeltime );
        }
        else if( nodeName == "aeei" ) {
            XMLHelper<double>::insertValueIntoVector( curr, mAEEI, modeltime );
        }
        else if( nodeName == "perCapitaBased" ) {
            mIsPerCapitaBased = XMLHelper<bool>::getValue( curr );
        }
        else {
            return false;
        }
    }
    return true;
}

/*! \brief Write object to xml output stream for use as a future input file.
* \details Method writes the contents of this object to the XML output stream.
* \author Steve Smith, Josh Lurz
* \param aOut reference to the output stream
* \param aTabs A aTabs object responsible for printing the correct number of aTabs. 
*/
void NonEnergyFinalDemand::toInputXML( ostream& aOut, Tabs* aTabs ) const {
	XMLWriteOpeningTag( getXMLNameStatic(), aOut, aTabs, mName );
    const Modeltime* modeltime = scenario->getModeltime();
   
    // write the xml for the class members.
    XMLWriteElementCheckDefault( mIsPerCapitaBased, "perCapitaBased", aOut, aTabs, false );
    for( unsigned int i = 0; i < mPriceElasticity.size(); ++i ){
        XMLWriteElementCheckDefault( mPriceElasticity[ i ], "priceelasticity", aOut, aTabs, 0.0, modeltime->getper_to_yr( i ) );
    }
    for( unsigned int i = 0; i < mServiceDemands.size(); ++i ){
        XMLWriteElementCheckDefault( mServiceDemands[ i ], "serviceoutput", aOut, aTabs, 0.0, modeltime->getper_to_yr( i ) );
    }
    for( unsigned int i = 0; i < mIncomeElasticity.size(); ++i ){
        XMLWriteElementCheckDefault( mIncomeElasticity[ i ], "incomeelasticity", aOut, aTabs, 0.0, modeltime->getper_to_yr( i ) );
    }
    for( unsigned int i = 0; i < mAEEI.size(); ++i ){
        XMLWriteElementCheckDefault( mAEEI[ i ], "aeei", aOut, aTabs, 0.0, modeltime->getper_to_yr( i ) );
    }

    XMLWriteClosingTag( getXMLNameStatic(), aOut, aTabs );
}	

/*! \brief Write information useful for debugging to XML output stream
*
* Function writes market and other useful info to XML. Useful for debugging.
*
* \author Josh Lurz
* \param aPeriod model period
* \param aOut reference to the output stream
* \param aTabs A tabs object responsible for printing the correct number of tabs.
*/
void NonEnergyFinalDemand::toDebugXML( const int aPeriod, ostream& aOut, Tabs* aTabs ) const {
	XMLWriteOpeningTag ( getXMLNameStatic(), aOut, aTabs, mName );
    // write the xml for the class members.
	XMLWriteElement( mIsPerCapitaBased, "perCapitaBased", aOut, aTabs );
    XMLWriteElement( mTechnicalChange[ aPeriod ], "cumm-tech-change", aOut, aTabs );

    // Now write out own members.
    XMLWriteElement( mServiceDemands[ aPeriod ], "service", aOut, aTabs );
    XMLWriteElement( mIncomeElasticity[ aPeriod ], "iElasticity", aOut, aTabs );
    XMLWriteElement( mPriceElasticity[ aPeriod ], "pElasticity", aOut, aTabs );
    XMLWriteElement( mAEEI[ aPeriod ], "aeei", aOut, aTabs );

    XMLWriteClosingTag( getXMLNameStatic(), aOut, aTabs );
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
const string& NonEnergyFinalDemand::getXMLNameStatic() {
    const static string XML_NAME = "non-energy-final-demand";
	return XML_NAME;
}

//! Get the name of the NonEnergyFinalDemand. This may not be needed.
const string& NonEnergyFinalDemand::getName() const {
    return mName;
}

/* \brief Tabulate the fixed demands of the energy final demand.
*
* For now this sets -1 to flag that final demands are not fixed.
* \author Steve Smith
* \param period model period.
*/
void NonEnergyFinalDemand::tabulateFixedDemands( const string& aRegionName,
                                                 const Demographic* aDemographics,
                                                 const GDP* aGDP,
                                                 const int aPeriod ) const
{
	Marketplace* marketplace = scenario->getMarketplace();
	marketplace->getMarketInfo( mName, aRegionName, aPeriod, true )->setDouble( "calDemand", 1 );
}

void NonEnergyFinalDemand::scaleCalibratedValues( const string& aFuelName,
                                                  const double aScaleValue,
                                                  const int aPeriod )
{
    // Only scale if the input is the input of the final demand, or if all inputs should be scaled.
    if( aFuelName == mName || aFuelName == "allInputs" ){
        // Adjust the base service level.
        mBaseService *= aScaleValue;
    }
}

/*! \brief Calculate the final demand and set it into the marketplace.
* \details Calls the internal method to calculate final demand, stores the value for the period, and sets it into the marketplace.
* \param aRegionName Region name.
* \param aDemographics Regional demographics.
* \param aGDP Regional GDP.
* \param aPeriod model period.
*/
void NonEnergyFinalDemand::setFinalDemand( const string& aRegionName,
                                        const Demographic* aDemographics,
                                        const GDP* aGDP,
                                        const int aPeriod )
{
    mServiceDemands[ aPeriod ] = calcDemand( aRegionName, aDemographics, aGDP, aPeriod );
    
    // Set the service demand into the marketplace.
    Marketplace* marketplace = scenario->getMarketplace();
    marketplace->addToDemand( mName, aRegionName, mServiceDemands[ aPeriod ], aPeriod );
}

/*! \brief Aggrgate sector energy service demand function
*
* Function calculates the aggregate demand for energy services and passes that down to the sub-sectors. 
* Demand is proportional to either GDP (to a power) or GDP per capita (to a power) times population.
*
* \author Sonny Kim, Josh Lurz
* \param gdp aGDP object for calculating various types of gdps.
* \param aPeriod Model aPeriod
*/
double NonEnergyFinalDemand::calcDemand( const string& aRegionName,
                                      const Demographic* aDemographics,
                                      const GDP* aGDP,
                                      const int aPeriod ) const
{
    // Calculate demand
    double adjustedDemand;
    if( aPeriod == 0 ){
        adjustedDemand = mBaseService;
    }
    else {
        const Modeltime* modeltime = scenario->getModeltime();
        const int normPeriod = modeltime->getyr_to_per( 1990 );
        double priceRatio = SectorUtils::calcPriceRatio( aRegionName, mName, normPeriod, aPeriod );

        const int basePer = modeltime->getyr_to_per( modeltime->getStartYear() );
        double gdpRatio = aGDP->getGDP( aPeriod ) / aGDP->getGDP( basePer );
        // If perCapitaBased, service_demand = B * P^r * GDPperCap^r * Population.
        // All values are relative to the base year
        double serviceDemand;
        if ( mIsPerCapitaBased ) { // demand based on per capita GDP
            double scaledGDPperCap = aGDP->getScaledGDPperCap( aPeriod );
            serviceDemand = mServiceDemands[ 0 ] * pow( priceRatio, mPriceElasticity[ aPeriod ] ) 
                * pow( scaledGDPperCap, mIncomeElasticity[ aPeriod ] );
            // need to multiply above by population ratio (aNodeent population/base year
            // population).  This ratio provides the population ratio.
            serviceDemand *= gdpRatio / scaledGDPperCap;

        }
        // If not perCapitaBased, service_demand = B * P^r * GDP^r
        else { // demand based on scale of GDP    
            serviceDemand = mServiceDemands[ 0 ] * pow( priceRatio, mPriceElasticity[ aPeriod ] ) 
                * pow( gdpRatio, mIncomeElasticity[aPeriod] );
        }

        // demand sector output is total end-use sector demand for service
        // adjust demand using cummulative technical change.
        assert( mTechnicalChange[ aPeriod ] != -1 );
        adjustedDemand = serviceDemand / mTechnicalChange[ aPeriod ];
    }

    // Return the final demand adjusted for technical change.
    return adjustedDemand;
}

double NonEnergyFinalDemand::getWeightedEnergyPrice( const string& aRegionName, const int aPeriod ) const {
    return 0;
}

//! Calculate technical change
void NonEnergyFinalDemand::calcTechChange( const int aPeriod ){
    if( aPeriod > 0 ){
        // calculate cummulative technical change using AEEI, autonomous end-use energy intensity
        mTechnicalChange[ aPeriod ] = mTechnicalChange[ aPeriod - 1 ] * pow( 1 + mAEEI[ aPeriod ], 
            scenario->getModeltime()->gettimestep( aPeriod ) );
    }
    else {
        // There is no tech change in the base period.
        mTechnicalChange[ aPeriod ] = 1;
    }
}

// Documentation is inherited.
void NonEnergyFinalDemand::accept( IVisitor* aVisitor, const int aPeriod ) const {
    aVisitor->startVisitFinalDemand( this, aPeriod );
    aVisitor->endVisitFinalDemand( this, aPeriod );
}

//! Write sector output to database.
void NonEnergyFinalDemand::csvOutputFile( const string& aRegionName ) const {
    // function protocol
    void fileoutput3( string var1name,string var2name,string var3name,
        string var4name,string var5name,string uname,vector<double> dout);
    
    // function arguments are variable name, double array, db name, table name
    // the function writes all years
    // total Sector output
    fileoutput3( aRegionName, mName, " ", " ", "demand", "SerUnit", mServiceDemands);
}

//! Write MiniCAM style demand sector output to database.
void NonEnergyFinalDemand::dbOutput( const string& aRegionName ) const {
    const Modeltime* modeltime = scenario->getModeltime();
    const int maxper = modeltime->getmaxper();
    vector<double> temp(maxper);
    
    // function protocol
    void dboutput4(string var1name,string var2name,string var3name,string var4name,
        string uname,vector<double> dout);
    
    // total sector output
    dboutput4( aRegionName,"End-Use Service","by Sector", mName, "Ser Unit", mServiceDemands );

    // End-use service price elasticity
    dboutput4( aRegionName,"End-Use Service","Elasticity", mName + "_price" ," ", mPriceElasticity );
    // End-use service income elasticity
    dboutput4( aRegionName,"End-Use Service","Elasticity", mName + "_income"," ", mIncomeElasticity );
}
