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
* \file subsector_add_techcosts.cpp
* \ingroup Objects
* \brief SubsectorAddTechCostsAddTechCosts class source file.
* \author Sonny Kim
*/

#include "util/base/include/definitions.h"
#include <string>
#include <cassert>
#include <xercesc/dom/DOMNode.hpp>

#include "sectors/include/subsector_add_techcosts.h"
#include "technologies/include/itechnology_container.h"
#include "technologies/include/itechnology.h"

using namespace std;
using namespace xercesc;

/*! \brief Default constructor.
*
* Constructor initializes member variables with default values, sets vector sizes, etc.
*
* \author Sonny Kim
*/
SubsectorAddTechCosts::SubsectorAddTechCosts( const string& aRegionName, const string& aSectorName )
: Subsector(aRegionName, aSectorName)
{
}

//! Parses any input variables specific to derived classes
bool SubsectorAddTechCosts::XMLDerivedClassParse( const string nodeName, const DOMNode* curr ) {
    return false;
}

/*! \brief Get the XML node name for output to XML.
*
* This public function accesses the private constant string, XML_NAME.
* This way the tag is always consistent for both read-in and output and can be easily changed.
* This function may be virtual to be overriden by derived class pointers.
* \author Josh Lurz, James Blackwood
* \return The constant XML_NAME.
*/
const string& SubsectorAddTechCosts::getXMLName() const {
	return getXMLNameStatic();
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
const string& SubsectorAddTechCosts::getXMLNameStatic() {
    const static string XML_NAME = "subsectorAddTechCosts";
	return XML_NAME;
}

/*! \brief Computes total costs of all technologies in SubsectorAddTechCosts.
*
* Called from calcShare after technology shares are determined. 
* Calculates total price (subsectorprice) and cost of fuel (fuelprice). 
*
* Price function separated to allow different weighting for SubsectorAddTechCosts price
* changed to void return maw
*
* \author Sonny Kim
* \param regionName region name
* \param period Model period
*/
double SubsectorAddTechCosts::getPrice( const GDP* aGDP, const int aPeriod ) const {
    double subsectorPrice = 0;
    for ( unsigned int i = 0; i < mTechContainers.size(); ++i) {
        const ITechnology* newVintageTech = mTechContainers[ i ]->getNewVintageTechnology( aPeriod );
		double currCost = newVintageTech->getCost( aPeriod );
        // calculate total price for SubsectorAddTechCosts
		// subsector price is additive of all technology costs
		// **** Beware of Share Weights
		// **** Share Weights are used here to turn off technology(or not include) and
		// **** to adjust total technology cost and its additive contribution to
		// **** subsector price.
		if( currCost > 0 ){
			subsectorPrice += newVintageTech->getShareWeight() * currCost;
		}
    }
	// Check for the condition where all technologies were fixed.
	return ( subsectorPrice > 0 ) ? subsectorPrice : -1;
}
