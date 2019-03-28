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
* \file production_input.cpp
* \ingroup Objects
* \brief The Production Input class source file.
*
*  Detailed description.
*
* \author Pralit Patel
* \author Sonny Kim
*/

#include "util/base/include/definitions.h"
#include "functions/include/production_input.h"
#include "util/base/include/xml_helper.h"
#include "util/base/include/ivisitor.h"

using namespace std;
using namespace xercesc;

// static initialize.
const string ProductionInput::XML_REPORTING_NAME = "input-production";

//! Default Constructor
ProductionInput::ProductionInput() {
}

/*! \brief Creates a clone of this class
 *
 * \author Pralit Patel
 * \warning Does not copy everything, only non calculated values to pass on to the next period
 * \return Pointer to the new class created
 */
ProductionInput* ProductionInput::clone() const {
    return new ProductionInput( *this );
}

void ProductionInput::copyParam( const IInput* aInput,
                                 const int aPeriod ) {
    SGMInput::copyParam( aInput, aPeriod );
    aInput->copyParamsInto( *this, aPeriod );
}

void ProductionInput::copyParamsInto( ProductionInput& aProductionInput,
                                      const int aPeriod ) const {
    aProductionInput.mSalesTaxRate.init( mSalesTaxRate );
}

void ProductionInput::copyParamsInto( DemandInput& aDemandInput,
                                      const int aPeriod ) const {
	/*! \invariant This should never be called. */
	assert( false );
}

//! XML parsing for derived class
bool ProductionInput::XMLDerivedClassParse( const string& nodeName, const DOMNode* curr ) {
    return false;
}

//! Output XML for derived class
void ProductionInput::toInputXMLDerived( ostream& out, Tabs* tabs ) const {
}

//! Output debug info to XML for derived class
void ProductionInput::toDebugXMLDerived( const int period, ostream& out, Tabs* tabs ) const {
}

//! Get Price Elasticity
double ProductionInput::getPriceElasticity() const {
	return 0;
}

//! Get Income Elasticity
double ProductionInput::getIncomeElasticity() const {
	return 0;
}

/*! \brief Get the XML node name for output to XML.
*
* This public function accesses the private constant string, XML_NAME.
* This way the tag is always consistent for both read-in and output and can be easily changed.
* This function may be virtual to be overridden by derived class pointers.
* \author Josh Lurz, James Blackwood
* \return The constant XML_NAME.
*/
const string& ProductionInput::getXMLName() const {
    return getXMLNameStatic();
}

/*! \brief Get the XML node name in static form for comparison when parsing XML.
* \details This public function accesses the private constant string, XML_NAME.
*          This way
* the tag is always consistent for both read-in and output and can be easily
* changed. The "==" operator that is used when parsing, required this second
* function to return static.
* \note A function cannot be static and virtual.
* \author Josh Lurz, James Blackwood
* \return The constant XML_NAME as a static.
*/
const string& ProductionInput::getXMLNameStatic() {
    const static string XML_NAME = "productionInput";
    return XML_NAME;
}

/*! \brief Get the XML name for reporting to XML file.
*
* This public function accesses the private constant string, XML_NAME. This way
* the tag is always consistent for reporting outputs and can be easily
* changed.
* \author Sonny Kim
* \return The constant XML_NAME.
*/
const string& ProductionInput::getXMLReportingName() const{
    return XML_REPORTING_NAME;
}

void ProductionInput::accept( IVisitor* aVisitor, const int aPeriod ) const {
    aVisitor->startVisitProductionInput( this, aPeriod );
    SGMInput::accept( aVisitor, aPeriod );
    aVisitor->endVisitProductionInput( this, aPeriod );
}
