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
* \file accumulated_grade.cpp
* \ingroup Objects
* \brief AccumulatedGrade class source file.
* \author Sonny Kim
*/

#include "util/base/include/definitions.h"
#include "resources/include/accumulated_grade.h"
#include "marketplace/include/marketplace.h"
#include "containers/include/scenario.h"
#include "containers/include/iinfo.h"
#include "util/base/include/model_time.h"

using namespace std;
using namespace xercesc;

extern Scenario* scenario;

//! Default constructor
AccumulatedGrade::AccumulatedGrade(){
}

/*! \brief Get the XML node name for output to XML.
* \details This public function accesses the private constant string, XML_NAME.
* This way the tag is always consistent for both read-in and output and can be easily changed.
* This function may be virtual to be overriden by derived class pointers.
* \author Josh Lurz, James Blackwood
* \return The constant XML_NAME.
*/
const string& AccumulatedGrade::getXMLName() const {
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
const string& AccumulatedGrade::getXMLNameStatic() {
    const static string XML_NAME = "accumulatedGrade";
    return XML_NAME;
}

/*! \brief Perform any initializations needed for each period.
* Any initializations or calcuations that only need to be done once per period 
* (instead of every iteration) should be placed in this function.
*
* \author Sonny Kim
* \param period Model period
*/
void AccumulatedGrade::initCalc( const string& aRegionName, const string& aResourceName, const int aPeriod ) {
    // This is called only once in the beginning of each period.
    // Add to available resource the amount of accumulated resource from a by-product   // SHK 10/28/04
    Marketplace* marketplace = scenario->getMarketplace();

    // none accumulated for period = 0
    // for period = 1, take slope of 0 and last period accumulated amount
    double lastPerAmount = 0; // annual amount last period
    double lastLastPerAmount = 0; // annual amount two periods before
    if ( aPeriod >= 1 ) {
        const IInfo* resourceInfo = marketplace->getMarketInfo( aResourceName, aRegionName, aPeriod - 1, true );
        lastPerAmount = resourceInfo->getDouble( "AccumulatedRsc", true );
    }
    // for period > 1, take slope of last and the last two period accumulated amount
    if ( aPeriod > 1 ) {
        const IInfo* prevResourceInfo = marketplace->getMarketInfo( aResourceName, aRegionName, aPeriod - 2, true );
        lastLastPerAmount = prevResourceInfo->getDouble( "AccumulatedRsc", true );
    }
    double diffAmount = lastPerAmount - lastLastPerAmount;
    const Modeltime* modeltime = scenario->getModeltime();
    double accumulatedAmount = ( lastLastPerAmount + 0.5 * diffAmount ) * modeltime->gettimestep( aPeriod );
    // if assuming steady-state for nuclear breeder reactors
    //accumulatedAmount = ( lastPerAmount ) * modeltime->gettimestep( period );
    mAvailable += accumulatedAmount; // add total accumulated from current period
}

/*! \brief Perform any initializations needed after each period.
* \details Any initializations or calcuations that only need to be done once per
*          period (instead of every iteration) should be placed in this
*          function.
* \author Sonny Kim
* \param aRegionName Region name.
* \param aResourceName Resource name.
* \param aPeriod Model period.
*/
void AccumulatedGrade::postCalc( const string& aRegionName,
                                 const string& aResourceName,
                                 const int aPeriod )
{
}
