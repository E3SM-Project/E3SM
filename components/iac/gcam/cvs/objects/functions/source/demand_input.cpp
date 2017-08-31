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
* \file demand_input.cpp
* \ingroup Objects
* \brief The Demand Input source file.
*
*  Detailed description.
*
* \author Pralit Patel
* \author Sonny Kim
*/

#include "util/base/include/definitions.h"
#include "functions/include/demand_input.h"
#include "util/base/include/xml_helper.h"
#include "util/base/include/ivisitor.h"

using namespace std;

//! Default Constructor
DemandInput::DemandInput() {
}

/*! \brief Creates a clone of this class
 *
 * \author Pralit Patel
 * \warning Does not copy everything, only non calculated values to pass on to the next period
 * \return Pointer to the new class created
 */
DemandInput* DemandInput::clone() const {
    return new DemandInput( *this );
}

void DemandInput::copyParam( const IInput* aInput,
                             const int aPeriod ) {
    SGMInput::copyParam( aInput, aPeriod );
    aInput->copyParamsInto( *this, aPeriod );
}

void DemandInput::copyParamsInto( DemandInput& aDemandInput,
                                  const int aPeriod ) const {
    aDemandInput.mPriceElasticity.init( mPriceElasticity );
    aDemandInput.mIncomeElasticity.init( mIncomeElasticity );
}

void DemandInput::copyParamsInto( ProductionInput& aProductionInput,
                                  const int aPeriod ) const {
	/*! \invariant This should never be called. */
	assert( false );
}

//! Get Price Elasticity
double DemandInput::getPriceElasticity() const {
    return mPriceElasticity;
}

//! Get Income Elasticity
double DemandInput::getIncomeElasticity() const {
    return mIncomeElasticity;
}

//! XML parsing for derived class
bool DemandInput::XMLDerivedClassParse( const string& nodeName, const xercesc::DOMNode* curr ) {
    if ( nodeName == "incomeElasticity" ) {
        mIncomeElasticity = XMLHelper<double>::getValue( curr );
    }
    else if ( nodeName == "priceElasticity" ) {
        mPriceElasticity = XMLHelper<double>::getValue( curr );
    }
    else {
        return false;
    }
    return true;
}

//! Output XML for derived class
void DemandInput::toInputXMLDerived( ostream& out, Tabs* tabs ) const {
    XMLWriteElement( mIncomeElasticity, "incomeElasticity", out, tabs );
    XMLWriteElement( mPriceElasticity, "priceElasticity", out, tabs );
}

//! Output debug info to XML for derived class
void DemandInput::toDebugXMLDerived( const int period, ostream& out, Tabs* tabs ) const {
    XMLWriteElement( mIncomeElasticity, "incomeElasticity", out, tabs );
    XMLWriteElement( mPriceElasticity, "priceElasticity", out, tabs );
}

/*! \brief Get the XML node name for output to XML.
*
* This public function accesses the private constant string, XML_NAME.
* This way the tag is always consistent for both read-in and output and can be easily changed.
* This function may be virtual to be overridden by derived class pointers.
* \author Josh Lurz, James Blackwood
* \return The constant XML_NAME.
*/
const string& DemandInput::getXMLName() const {
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
const string& DemandInput::getXMLNameStatic() {
    const static string XML_NAME = "demandInput";
    return XML_NAME;
}

void DemandInput::accept( IVisitor* aVisitor, const int aPeriod ) const {
    aVisitor->startVisitDemandInput( this, aPeriod );
    SGMInput::accept( aVisitor, aPeriod );
    aVisitor->endVisitDemandInput( this, aPeriod );
}
