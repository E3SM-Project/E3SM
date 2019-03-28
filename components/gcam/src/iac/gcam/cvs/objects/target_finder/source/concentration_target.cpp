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
* \file concentration_target.cpp
* \ingroup Objects
* \brief ConcentrationTarget class source file.
* \author Josh Lurz
*/

#include "util/base/include/definitions.h"
#include <cassert>
#include <string>
#include "containers/include/scenario.h"
#include "util/base/include/model_time.h"
#include "climate/include/iclimate_model.h"
#include "target_finder/include/concentration_target.h"
#include "util/base/include/configuration.h"
#include "util/logger/include/ilogger.h"

using namespace std;

extern Scenario* scenario;

/*!
 * \brief Constructor
 * \param aClimateModel The climate model.
 * \param aTargetValue The target value.
 * \param aFirstTaxYear The first tax year.
 */
ConcentrationTarget::ConcentrationTarget( const IClimateModel* aClimateModel,
                                          const double aTargetValue,
                                          const int aFirstTaxYear ):
mClimateModel( aClimateModel ),
mTargetValue( aTargetValue ),
mFirstTaxYear( aFirstTaxYear )
{
    // Store configuration variables.
    const Configuration* conf = Configuration::getInstance();
    mTargetGas = conf->getString( "concentration-target-gas", "CO2" );
}

/*! \brief Return the static name of the object.
 * \return The static name of the object.
 */
const string& ConcentrationTarget::getXMLNameStatic(){
	static const string XML_NAME = "concentration-target";
	return XML_NAME;
}

/*!
 * \brief Get the status of the last trial with respect to the target.
 * \param aTolerance Solution tolerance.
 * \param aYear Year in which to get the status.
 * \return Status of the last trial.
 */
double ConcentrationTarget::getStatus( const int aYear ) const
{
    // Make sure we are using the correct year.
    const int year = aYear == ITarget::getUseMaxTargetYearFlag() ? getYearOfMaxTargetValue()
        : aYear;
    /*!
     * \pre year must be greater than mFirstTaxYear otherwise we will have no
     *      ability to change the status in that year.
     */
    assert( year >= mFirstTaxYear );
    const double currConcentration = mClimateModel->getConcentration( mTargetGas, year );

    // Determine how how far away from the target the current estimate is.
    double percentOff = ( currConcentration - mTargetValue ) / mTargetValue * 100;
    
    // Print an information message.
    ILogger& mainLog = ILogger::getLogger( "main_log" );
    mainLog.setLevel( ILogger::NOTICE );
    mainLog << "Currently " << percentOff << " percent away from the concentration target." << endl;
    mainLog << "Current: " << currConcentration << " Target: " << mTargetValue
            << " In year " << year << endl;
    // Print information to target log since target info is not available elsewhere
    ILogger& targetLog = ILogger::getLogger( "target_finder_log" );
    targetLog.setLevel( ILogger::NOTICE );
    targetLog << "Target Values - Current: " << currConcentration << " Target: " << mTargetValue
              << " In year: " << year << endl;

    return percentOff;
}

int ConcentrationTarget::getYearOfMaxTargetValue() const {
    const int finalYearToCheck = scenario->getModeltime()->getEndYear();
    double maxConc = 0;
    int maxYear = mFirstTaxYear - 1;
    
    // Loop over possible year and find the max concentration and the year it occurs in.
    for( int year = mFirstTaxYear; year <= finalYearToCheck; ++year ) {
        if( maxConc < mClimateModel->getConcentration( mTargetGas, year ) ) {
            maxConc = mClimateModel->getConcentration( mTargetGas, year );
            maxYear = year;
        }
    }
    return maxYear;
}
