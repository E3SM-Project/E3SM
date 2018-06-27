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
* \file final_demand_sector.cpp
* \ingroup Objects
* \brief The Final Demand Sector class source file.
*
*  Detailed Description.
*
* \author Pralit Patel
* \author Sonny Kim
*/

#include "util/base/include/definitions.h"
#include <string>
#include <xercesc/dom/DOMNode.hpp>

#include "sectors/include/final_demand_sector.h"
#include "sectors/include/subsector.h"

using namespace std;

//! Default constructor
FinalDemandSector::FinalDemandSector( const string& aRegionName ):Sector( aRegionName ){
}

//! Destructor
FinalDemandSector::~FinalDemandSector() {}

//! Setup the markets for the FinalDemandSector. This currently does nothing.
void FinalDemandSector::setMarket(){
}

//! Parse xml file for data
bool FinalDemandSector::XMLDerivedClassParse( const string& nodeName, const xercesc::DOMNode* curr ) {
    return false;
}

//! For derived classes to output XML data
void FinalDemandSector::toInputXMLDerived( ostream& out, Tabs* tabs ) const {
	//do something maybe
}

//! Output debug info for derived class
void FinalDemandSector::toDebugXMLDerived( const int period, ostream& out, Tabs* tabs ) const {
}

const string& FinalDemandSector::getXMLName() const {
	return getXMLNameStatic();
}

const string& FinalDemandSector::getXMLNameStatic() {
    const static string XML_NAME = "finalDemandSector";
	return XML_NAME;
}

/*!
 * \brief Get the FinalDemandSector price.
 * \param aGDP Regional GDP container(null for SGM).
 * \param aPeriod Model period.
 * \return Final demand sectors do not have a price, so return 0.
 */
double FinalDemandSector::getPrice( const GDP* aGDP,
                                    const int aPeriod ) const {
    return 0;
}

//! Operate the consumers.
void FinalDemandSector::operate( NationalAccount& aNationalAccount, const Demographic* aDemographic,
                                 const int aPeriod )
{
    for( SubsectorIterator currSub = subsec.begin(); currSub != subsec.end(); ++currSub ){
        // flag tells the subsector to operate all capital, old and new.
        (*currSub)->operate( aNationalAccount, aDemographic, moreSectorInfo.get(), true, aPeriod );
    }
}

/*! \brief Perform any initializations needed for each period.
*
* Any initializations or calculations that only need to be done once per period (instead of every iteration) should be placed in this function.
*
* \author Sonny Kim
* \param aNationalAccount NationalAccount object
* \param aDemographics Demographics object
* \param aPeriod Period to initialize
*/
void FinalDemandSector::initCalc( NationalAccount* aNationalAccount,
                                  const Demographic* aDemographics,
                                  const int aPeriod )
{
    // do any sub-Sector initializations
    for( SubsectorIterator currSub = subsec.begin(); currSub != subsec.end(); ++currSub ){
        (*currSub)->initCalc( aNationalAccount, aDemographics,
                              moreSectorInfo.get(), aPeriod );
    }
}

/*! \brief Complete the initialization of the final demand sector.
* \param aRegionInfo The regional information object.
* \param aDependencyFinder Region's dependency finder, should be null for CGE regions.
* \param aLandAllocator Region's land allocator, should be null for CGE regions.
*/
void FinalDemandSector::completeInit( const IInfo* aRegionInfo,
                                      DependencyFinder* aDependencyFinder,
                                      ILandAllocator* aLandAllocator )
{
    Sector::completeInit( aRegionInfo, aDependencyFinder, aLandAllocator );
}
