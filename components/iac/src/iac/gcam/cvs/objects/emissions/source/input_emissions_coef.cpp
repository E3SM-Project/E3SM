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
 * \file input_emissions_coef.h
 * \ingroup Objects
 * \brief InputEmissionsCoef header file.
 * \author Jim Naslund
 */

#include "util/base/include/definitions.h"

#include "emissions/include/input_emissions_coef.h"
#include "util/base/include/util.h"

using namespace std;

//! Constructor that initializes input emissions.
InputEmissionsCoef::InputEmissionsCoef( const double aInputEmissions ):
AEmissionsCoef(),
mInputEmissions( aInputEmissions ){
    // Cannot use member init. list because mOverridesFuturePeriods is a member of AEmissionsCoef
    mOverridesFuturePeriods = true;
}

//! Clone operator.
InputEmissionsCoef* InputEmissionsCoef::clone() const {
    return new InputEmissionsCoef( *this );
}

void InputEmissionsCoef::updateCoef( const double adjEmissDriver ){
    // Updates the emissions coefficient to be proportional to emissions divided by driver.
    // Check for divide by zero.
    mEmissionsCoef = adjEmissDriver > util::getSmallNumber() ? mInputEmissions / adjEmissDriver : 0;
}

void InputEmissionsCoef::initCalc( const IInfo* aSubsectorInfo, const string& aName, const int aPeriod ){
}

double InputEmissionsCoef::getInputEmissions() const {
    return mInputEmissions;
}

/*!
* \brief Calculate the maxCntrl parameter of the control function
* \details maxCntrl is the most a non-CO2 gas can be reduced. The control
*          function approaches this limit as gdp per capita grows. This 
*          method calculates and returns that parameter
* \param aFinalEmissCoef The minimum emissions per unit of input allowed
* \param aB Parameter used to convert aFinalEmissCoef to maxCntrl
* \param aMultiplier The input driver ( e.g., energy consumed ) adjusted
*                    for any MAC or emAdjust applicable
* \return The maximum amount of control allowed for this non-CO2 gas 
*/
double InputEmissionsCoef::calcMaxCntrl( const double aFinalEmissCoef, const double aB,
                                         const double aMultiplier ) const {
    double maxControl = 0.0;
    // Guard against divide by 0
    if ( aMultiplier > util::getSmallNumber() ) {
        double numerator = aFinalEmissCoef - ( mInputEmissions / aMultiplier );
        double denominator = ( aFinalEmissCoef / aB ) - ( mInputEmissions / aMultiplier );
        maxControl = 100.0 * ( numerator / denominator );
    }

    return maxControl;
}

bool InputEmissionsCoef::needsCalcForAdjustment() const {
    return true;
}

const string& InputEmissionsCoef::getXMLName() const{
    static const string XML_NAME = "input-emissions";
    return XML_NAME;
}

double InputEmissionsCoef::getXMLValue() const{
    return mInputEmissions;
}
