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
* \file tran_subsector.cpp
* \ingroup Objects
* \brief transporation technology class source file.
* \author Sonny Kim, Josh Lurz, Steve Smith, Marshall Wise
*/

#include "util/base/include/definitions.h"
#include <string>
#include <iostream>
#include <cassert>
#include <vector>
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>

#include "sectors/include/tran_subsector.h"
#include "containers/include/scenario.h"
#include "util/base/include/model_time.h"
#include "containers/include/info_factory.h"
#include "containers/include/iinfo.h"
#include "util/base/include/xml_helper.h"
#include "marketplace/include/marketplace.h"
#include "util/base/include/summary.h"
#include "containers/include/gdp.h"
#include "demographics/include/demographic.h"
#include "technologies/include/itechnology_container.h"
#include "technologies/include/itechnology.h"

using namespace std;
using namespace xercesc;

extern Scenario* scenario;
const string TranSubsector::XML_NAME = "tranSubsector";

/*  Begin TranSubsector Method Definitions */

/*! \brief Default constructor for TranSubsector.
*
* Default constructor takes region name, sector name and units as arguments
* and resizes and initializes vectors.
* \param regionName The name of the region.
* \param sectorName The name of the sector.
* \param aUnit The sector output unit.
* \author Josh Lurz, Sonny Kim
*/
TranSubsector::TranSubsector( const string& regionName, const string& sectorName ): Subsector( regionName, sectorName ) {
	// resize vectors
	const Modeltime* modeltime = scenario->getModeltime();
	const int maxper = modeltime->getmaxper();
	speed.resize( maxper ); // average speed of mode
	mPopulation.resize( maxper ); // copy of population since demog object not available
	popDenseElasticity.resize( maxper );
	popDensity = 1; // initialize to 1 for now
	mAddTimeValue = false; // initialize to false
}

/*! \brief Get the XML node name for output to XML.
*
* This public function accesses the private constant string, XML_NAME.
* This way the tag is always consistent for both read-in and output and can be easily changed.
* This function may be virtual to be overridden by derived class pointers.
* \author Josh Lurz, James Blackwood
* \return The constant XML_NAME.
*/
const std::string& TranSubsector::getXMLName() const {
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
const std::string& TranSubsector::getXMLNameStatic() {
    return XML_NAME;
}

/*! \brief Function Parses any input variables specific to derived classes.
* \param nodeName The name of the XML node.
* \param curr A pointer to the XML DOM node.
* \author Josh Lurz, Sonny Kim
* \return Boolean for node match.
*/
bool TranSubsector::XMLDerivedClassParse( const string& nodeName, const DOMNode* curr ) {    
	// additional read in for transportation
	const Modeltime* modeltime = scenario->getModeltime();
	if( nodeName == "addTimeValue" ){
		mAddTimeValue = XMLHelper<bool>::getValue( curr );
	}
	else if( nodeName == "speed" ){
		XMLHelper<double>::insertValueIntoVector( curr, speed, modeltime );
	}
	else if( nodeName == "popDenseElasticity" ){
		XMLHelper<double>::insertValueIntoVector( curr, popDenseElasticity, modeltime );
	}
	else {
		return false;
	}
	return true;
}

/*! \brief XML output stream for derived classes
*
* Function writes output due to any variables specific to derived classes to XML
* \author Josh Lurz, Sonny Kim
* \param out reference to the output stream
* \param tabs A tabs object responsible for printing the correct number of tabs. 
*/
void TranSubsector::toInputXMLDerived( ostream& out, Tabs* tabs ) const {
	XMLWriteElementCheckDefault( mAddTimeValue, "addTimeValue", out, tabs, false );
    const Modeltime* modeltime = scenario->getModeltime();
    XMLWriteVector( speed, "speed", out, tabs, modeltime, 0.0 );
    XMLWriteVector( popDenseElasticity, "popDenseElasticity", out, tabs, modeltime, 0.0 );
    XMLWriteVector( mServiceOutputs, "serviceoutput", out, tabs, modeltime, 0.0 );
}

/*! \brief XML output for debugging.
* Function writes output to debugging XML
* \author Josh Lurz, Sonny Kim
* \param out reference to the output stream
* \param tabs A tabs object responsible for printing the correct number of tabs. 
*/
void TranSubsector::toDebugXMLDerived( const int period, ostream& out, Tabs* tabs ) const {
	XMLWriteElement( mAddTimeValue, "addTimeValue", out, tabs );
	XMLWriteElement( popDenseElasticity[ period ], "popDenseElasticity", out, tabs );
	XMLWriteElement( popDensity, "popDensity", out, tabs );
	XMLWriteElement( speed[ period ], "speed", out, tabs );
    
    // Write out useful debugging info
	XMLWriteElement( mTimeValue, "timeValue", out, tabs );
}

/*! \brief Perform any initializations needed for each period.
* \author Sonny Kim, Steve Smith, Josh Lurz
* \param aNationalAccount National accounts information
* \param aDemographics Regional demographics information
* \param aMoreSectorInfo Additional sector information required for subsector
* \param aPeriod Model period
*/
void TranSubsector::completeInit( const IInfo* aSectorInfo,
                                  DependencyFinder* aDependencyFinder,
                                  ILandAllocator* aLandAllocator )
{
    // Only call base class completeInit.
    Subsector::completeInit( aSectorInfo, aDependencyFinder, aLandAllocator );
}

/*!
* \brief Perform any initializations needed for each period.
* \details Perform any initializations or calculations that only need to be done
*          once per period (instead of every iteration) should be placed in this
*          function.
* \warning The ghg part of this routine assumes the existence of technologies in
*          the previous and future periods
* \author Steve Smith, Sonny Kim
* \param aNationalAccount National accounts container.
* \param aDemographics Regional demographics container.
* \param aMoreSectorInfo sector info object.
* \param aPeriod Model period
*/
void TranSubsector::initCalc( NationalAccount* aNationalAccount,
							 const Demographic* aDemographics,
							 const MoreSectorInfo* aMoreSectorInfo,
							 const int aPeriod )
{
	// Check if illegal values have been read in
	if ( speed[ aPeriod ] <= 0 ) {
		speed[ aPeriod ] = 1;
		ILogger& mainLog = ILogger::getLogger( "main_log" );
		mainLog.setLevel( ILogger::ERROR );
		mainLog << "Speed was zero or negative in subsector: " << name << " in region " 
			<< regionName << ". Reset to 1." << endl;
	}
	// time in transit
	// initialize vector to hold population (thousands)
	// TODO: revise access to population to avoid statement below
	mPopulation[ aPeriod ] = aDemographics->getTotal( aPeriod );

	Subsector::initCalc( aNationalAccount, aDemographics, aMoreSectorInfo, aPeriod );
}

/*! \brief returns the subsector price.
* \details Calculates and returns share-weighted total price (subsectorprice)
*          with or without value of time.
* \author Sonny Kim
* \param aGDP Regional GDP object.
* \param aPeriod Model period
* \return The subsector price with or without value of time. 
*/
double TranSubsector::getPrice( const GDP* aGDP, const int aPeriod ) const {
	// mAddTimeValue is a boolean that determines whether the service price includes
	// the value of time
	if (mAddTimeValue) {
		return getGeneralizedPrice( aGDP, aPeriod );
	}
	// normal share-weighted total technology cost only
	return Subsector::getPrice( aGDP, aPeriod );
}

/*! \brief Get the time value for the period.
* \param aGDP The regional GDP container.
* \param aPeriod The model period.
* \author Sonny Kim
* \return The time value.
*/
double TranSubsector::getTimeValue( const GDP* aGDP, const int aPeriod ) const {
	const double WEEKS_PER_YEAR = 50;
	const double HOURS_PER_WEEK = 40;
	// calculate time value based on hours worked per year Convert GDPperCap
	// into dollars (instead of 1000's of $'s) GDP value at this point in the
	// code does not include energy feedback calculation for this year, so is,
	// therefore, approximate
	return aGDP->getApproxGDPperCap( aPeriod ) * 1000 / ( HOURS_PER_WEEK * WEEKS_PER_YEAR ) / speed[ aPeriod ];
}

/*! \brief Calculate the generalized service price for the mode that includes time value.
* \author Sonny Kim
* \param aGDP The regional GDP container.
* \param aPeriod The model period.
* \return The the generalized price.
*/
double TranSubsector::getGeneralizedPrice( const GDP* aGDP, const int aPeriod ) const {
	// add cost of time spent on travel by converting gdp/cap into an hourly
	// wage and multiplying by average speed.
	// The price unit is $ per service, e.g. $/pass-mi or $/ton-mi
    
    // Save time value so can print out
    // Maybe also write to XML DB?
    mTimeValue =  getTimeValue( aGDP, aPeriod );
	return Subsector::getPrice( aGDP, aPeriod ) + mTimeValue;
}

/*! \brief Get the time in transit per day per person for the period.
*  Currently used for reporting only.
* \author Sonny Kim
* \param aPeriod The model period.
* \return The time in transit.
*/
double TranSubsector::getTimeInTransit( const int aPeriod ) const {
	const double DAYS_PER_YEAR = 365;
	const double POP_MILE_CONV = 1000;
	// calculate time in transit per day for each person using total population
	return getOutput( aPeriod ) / mPopulation[ aPeriod ] * POP_MILE_CONV  
		/ speed[ aPeriod ] / DAYS_PER_YEAR ;
}

/*! \brief Get service per day per capita for the period.
* \author Sonny Kim
* \param aPeriod The model period.
* \return The service per day per capita.
*/
double TranSubsector::getServicePerCapita( const int aPeriod ) const {
	const double DAYS_PER_YEAR = 365;
	const double POP_MILE_CONV = 1000;
	// units: million pass or ton mi / thousand persons
	return getOutput( aPeriod ) / mPopulation[ aPeriod ] * POP_MILE_CONV
		/ DAYS_PER_YEAR ;
}

void TranSubsector::setOutput( const double aVariableSubsectorDemand,
							  const double aFixedOutputScaleFactor,
							  const GDP* aGDP,
							  const int aPeriod )
{
    Subsector::setOutput( aVariableSubsectorDemand,
                          aFixedOutputScaleFactor,
                          aGDP,
                          aPeriod );
}

/*! \brief Write supply sector MiniCAM style Subsector output to database.
*
* Writes outputs with titles and units appropriate to supply sectors.
*
* \author Sonny Kim
*/
void TranSubsector::MCoutputSupplySector( const GDP* aGDP ) const {
    // function protocol
    void dboutput4(string var1name,string var2name,string var3name,string var4name,
        string uname,vector<double> dout);
    
    const Modeltime* modeltime = scenario->getModeltime();
    const int maxper = modeltime->getmaxper();
	const string& outputUnit = mSubsectorInfo->getString( "output-unit", true );
	const string& priceUnit = mSubsectorInfo->getString( "price-unit", true );
    vector<double> temp(maxper);
    
    // total Subsector output
    for( int per = 0; per < maxper; ++per ){
        temp[ per ] = getOutput( per );
    }
    dboutput4(regionName,"Secondary Energy Prod",sectorName,name,outputUnit,temp);
    // Subsector price
	for( int m = 0; m < maxper; m++ ){
        temp[ m ] = getPrice( aGDP, m );
    }
    dboutput4(regionName,"Price",sectorName,name,priceUnit,temp);
    
    // do for all technologies in the Subsector
    for( unsigned int i = 0; i < mTechContainers.size(); ++i ){

        // secondary energy and price output by tech
        // output or demand for each Technology
        for ( int m=0;m<maxper;m++) {
            temp[m] = mTechContainers[i]->getNewVintageTechnology(m)->getOutput( m );
        }

        dboutput4( regionName, "Secondary Energy Prod", sectorName + "_tech-new-investment", 
            mTechContainers[i]->getName(), outputUnit, temp );
		
		// Output for all vintages.
        for ( int m=0; m < maxper;m++) {
			temp[ m ] = 0;
            // Only sum output to the current period.
			for( int j = 0; j <= m; ++j ){
				temp[m] += mTechContainers[i]->getNewVintageTechnology(j)->getOutput( m );
			}
        }
        dboutput4( regionName, "Secondary Energy Prod", sectorName + "_tech-total", mTechContainers[i]->getName(), outputUnit, temp );
        // Transportation technology cost already in 1990 $.
        for ( int m=0;m<maxper;m++) {
            temp[m] = mTechContainers[i]->getNewVintageTechnology(m)->getCost( m );
        }
        dboutput4( regionName, "Price", sectorName + "_tech", mTechContainers[i]->getName(), priceUnit, temp );
    }
}

/*! \brief Writes variables specific to transportation class to database.
*
* \author Sonny Kim, Steve Smith
* \param aGDP The GDP info object.
* \param aSectorOutput The vector of sector outputs.
*/
void TranSubsector::MCoutputAllSectors( const GDP* aGDP,
                                        const IndirectEmissionsCalculator* aIndirectEmissionsCalc,
                                        const vector<double> aSectorOutput ) const
{
	Subsector::MCoutputAllSectors( aGDP, aIndirectEmissionsCalc, aSectorOutput );

	// function protocol
	void dboutput4(string var1name,string var2name,string var3name,string var4name,
		string uname,vector<double> dout);
	const int maxPeriod = scenario->getModeltime()->getmaxper();
	const string& priceUnit = mSubsectorInfo->getString( "price-unit", true );
	vector<double> temp( maxPeriod );
	// Subsector timeValue price
	for( int per = 0; per < maxPeriod; ++per ){
		temp [ per ] = getTimeValue( aGDP, per );
	}
	dboutput4( regionName, "General", "TimeValue", sectorName + name, priceUnit, temp );
	// Subsector speed
	dboutput4( regionName, "General", "Speed", sectorName + name, "Miles/hr", speed );
	// time in transit (hours/day/person)
	for( int per = 0; per < maxPeriod; ++per ){
		temp[ per ] = getTimeInTransit( per );
	}
	dboutput4( regionName, "End-Use Service", sectorName+" TimeInTransit", name, "hrs/day/per", temp);
	// service per day per person (service/day/person)
	for( int per = 0; per < maxPeriod; ++per ){
		temp[ per ] = getServicePerCapita( per );
	}
}
