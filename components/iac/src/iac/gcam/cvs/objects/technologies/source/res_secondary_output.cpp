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
* \file secondary_output.cpp
* \ingroup Objects
* \brief RESSecondaryOutput class source file.
* \author Josh Lurz
*/
#include "util/base/include/definitions.h"
#include <string>
#include "technologies/include/res_secondary_output.h"
#include "containers/include/scenario.h"
#include "util/base/include/model_time.h"
#include "marketplace/include/marketplace.h"



using namespace std;
using namespace xercesc;

extern Scenario* scenario;

// static initialize.
const string RESSecondaryOutput::XML_REPORTING_NAME = "res-output-secondary";

/*! \brief Get the XML name for reporting to XML file.
*
* This public function accesses the private constant string, XML_NAME. This way
* the tag is always consistent for reporting outputs and can be easily
* changed.
* \author Sonny Kim
* \return The constant XML_NAME.
*/
const string& RESSecondaryOutput::getXMLReportingName() const{
    return XML_REPORTING_NAME;
}

const string& RESSecondaryOutput::getXMLNameStatic()
{
    const static string XML_NAME = "res-secondary-output";
    return XML_NAME;
}

RESSecondaryOutput::RESSecondaryOutput()
    : SecondaryOutput()
{
}

RESSecondaryOutput* RESSecondaryOutput::clone() const
{
    return new RESSecondaryOutput( *this );
}

bool RESSecondaryOutput::isSameType( const string& aType ) const
{
    return aType == getXMLNameStatic();
}

void RESSecondaryOutput::setPhysicalOutput( const double aPrimaryOutput,
                                        const string& aRegionName,
                                        ICaptureComponent* aCaptureComponent,
                                        const int aPeriod )
{
    // Secondary output is the primary output multiplied by the output ratio.
    mPhysicalOutputs[ aPeriod ] = calcPhysicalOutputInternal( aPrimaryOutput );
    
    // Remove the secondary output from demand instead of adding to supply
    // because the sector which has this output as a primary will attempt to
    // fill all of demand. If this technology also added to supply, supply would
    // not equal demand.
    Marketplace* marketplace = scenario->getMarketplace();
    marketplace->addToSupply( mName, aRegionName, mPhysicalOutputs[ aPeriod ], aPeriod, true );

}
