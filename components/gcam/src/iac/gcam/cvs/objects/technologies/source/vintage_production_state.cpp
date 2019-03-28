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
* \file vintage_production_state.cpp
* \ingroup Objects
* \brief VintageProductionState class source file.
* \author Josh Lurz
*/

#include "util/base/include/definitions.h"
#include <cassert>
#include "technologies/include/vintage_production_state.h"
#include "technologies/include/ishutdown_decider.h"
#include "util/base/include/xml_helper.h"
#include "technologies/include/marginal_profit_calculator.h"

using namespace std;

VintageProductionState::VintageProductionState():
mInitialYear( -1 ){
}

VintageProductionState* VintageProductionState::clone() const {
    return new VintageProductionState( *this );
}

bool VintageProductionState::isSameType( const string& aType ) const {
    return aType == getXMLNameStatic();
}

const string& VintageProductionState::getName() const {
    return getXMLNameStatic();
}

void VintageProductionState::toDebugXML( const int aPeriod,
                                        ostream& aOut,
                                        Tabs* aTabs ) const
{
    XMLWriteOpeningTag( getXMLNameStatic(), aOut, aTabs );
    XMLWriteElement( mBaseOutput, "base-output", aOut, aTabs, mInitialYear );
    XMLWriteClosingTag( getXMLNameStatic(), aOut, aTabs );
}

const string& VintageProductionState::getXMLNameStatic() {
    const static string XML_NAME = "vintage-production-state";
    return XML_NAME;
}

double VintageProductionState::calcProduction( const string& aRegionName,
                                              const string& aSectorName,
                                              const double aVariableOutput,
                                              const MarginalProfitCalculator* aMarginalProfitCalc,
                                              const double aFixedOutputScaleFactor,
                                              const vector<IShutdownDecider*>& aShutdownDeciders,
                                              const int aPeriod ) const
{
    // Check that setBaseOutput has been called to initialize the base output
    // and initial output year. 
    assert( mBaseOutput.isInited() && mInitialYear != 1 );

    double shutdownCoef = calcShutdownCoefficient( aRegionName, aSectorName,
        aShutdownDeciders,
        aMarginalProfitCalc,
        aPeriod );
    // Get the amount of output from the Technology in a period after the initial
    // investment period. This is calculated based on the aggregate shutdown
    // coefficient and the amount of output in the initial investment period.
    return mBaseOutput * shutdownCoef * aFixedOutputScaleFactor;
}

/*!
* \brief Return the aggregate shutdown coefficient which shuts down Technology
*          production.
* \details The aggregate shutdown coefficient is multiplied by the output in
*          any operating year to determine the Technology output. The aggregate
*          coefficient is the product of all the Technology's shutdown
*          decisions.
* \param aRegionName Region name.
* \param aSectorName Sector name.
* \param aShutdownDeciders Set of shutdown decision makers.
* \param aMarginalProfitCalc Calculator of the marginal profit rate for the
*        Technology.
* \param aPeriod Model period.
* \return The aggregate shutdown coefficient.
*/
double VintageProductionState::calcShutdownCoefficient( const string& aRegionName,
                                                       const string& aSectorName,
                                                       const vector<IShutdownDecider*>& aShutdownDeciders,
                                                       const MarginalProfitCalculator* aMarginalProfitCalc,
                                                       const int aPeriod ) const
{
    // Loop through the shutdown decision makers and calculate the product of
    // their coefficients.
    double shutdownCoef = 1;

    // Avoid expensive marginal profit calculation if there are no shutdown deciders.
    if( aShutdownDeciders.empty() ){
        return shutdownCoef;
    }

    double marginalProfit = aMarginalProfitCalc->calcShortTermMarginalProfit( aRegionName,
        aSectorName,
        aPeriod );

    for( unsigned int i = 0; i < aShutdownDeciders.size(); ++i ){
        shutdownCoef *= aShutdownDeciders[ i ]->calcShutdownCoef( 0, marginalProfit, aRegionName,
            aSectorName, mInitialYear, aPeriod );
    }
    return shutdownCoef;
}

void VintageProductionState::setBaseOutput( const double aBaseOutput,
                                           const int aBaseYear )
{
    assert( aBaseOutput >= 0 );

    // Store the base level of output and the period in which it was produced.
    mBaseOutput.set( aBaseOutput );
    mInitialYear = aBaseYear;
}

bool VintageProductionState::isNewInvestment() const {
    return false;
}

bool VintageProductionState::isOperating() const {
    return true;
}
