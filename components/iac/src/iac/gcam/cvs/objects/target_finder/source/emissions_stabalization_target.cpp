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
 * \file emissions_stabalization_target.cpp
 * \ingroup Objects
 * \brief EmissionsStabalizationTarget class source file.
 * \author Josh Lurz
 */

#include "util/base/include/definitions.h"
#include <cassert>
#include <string>
#include "containers/include/scenario.h"
#include "util/base/include/model_time.h"
#include "climate/include/iclimate_model.h"
#include "target_finder/include/emissions_stabalization_target.h"
#include "util/logger/include/ilogger.h"
#include "util/base/include/configuration.h"
#include "util/base/include/util.h"

extern Scenario* scenario;

using namespace std;

/*!
 * \brief Constructor
 * \param aClimateModel The climate model.
 * \param aTargetValue The value of the target.
 * \param aFirstTaxYear The first tax year.
 */
EmissionsStabalizationTarget::EmissionsStabalizationTarget(
    const IClimateModel* aClimateModel,
    const double aTargetValue,
    const int aFirstTaxYear ):
mClimateModel( aClimateModel ),
mFirstTaxYear( aFirstTaxYear )
{
    // Store configuration variables.
    const Configuration* conf = Configuration::getInstance();
    mTargetGas = conf->getString( "concentration-target-gas", "CO2" );
}

/*!
 * \brief Return the static name of the object.
 * \return The static name of the object.
 */
const string& EmissionsStabalizationTarget::getXMLNameStatic(){
	static const string XML_NAME = "stabilization-target";
	return XML_NAME;
}

double EmissionsStabalizationTarget::getStatus( const int aYear ) const
{
    // Determine the total system emissions. NOTE: Net terrestrial uptake
    // includes land use emissions as a negative, so they are not added here as
    // that would double account.
    
    // NOTE: Assumes that net terrestrial uptake is not removed from industrial
    // emissions.

    // Make sure we are using the correct year.
    const int year = aYear == ITarget::getUseMaxTargetYearFlag() ? getYearOfMaxTargetValue()
        : aYear;
    /*!
     * \pre year must be greater than mFirstTaxYear otherwise we will have no
     *      ability to change the status in that year.
     */
    assert( year >= mFirstTaxYear );

    double totalEmissions = mClimateModel->getEmissions( "CO2", year );

    // Year can be between two years
    double netOceanUp = mClimateModel->getNetOceanUptake( year );

    double netTerrUp = mClimateModel->getNetTerrestrialUptake( year );

    double totalNetUptake = netOceanUp + netTerrUp;

    // Determine how how far away from the target the current estimate is.
    double percentOff = objects::percentDiff( totalNetUptake, totalEmissions );
    
    // Print an information message.
    ILogger& targetLog = ILogger::getLogger( "target_finder_log" );
    targetLog.setLevel( ILogger::NOTICE );
    targetLog << "Currently " << percentOff
              << " percent away from the stabilization target. "
              << "Emissions: " << totalEmissions << " Total net uptake: "
              << totalNetUptake << endl;

    return percentOff;
}

int EmissionsStabalizationTarget::getYearOfMaxTargetValue() const {
    // TODO: I think realistically this is always the last year.
    
    double minDiff = numeric_limits<double>::max();
    const int finalYearToCheck = scenario->getModeltime()->getEndYear();
    int minYear = mFirstTaxYear - 1;
    
    // Loop over possible year and find the min difference in emisssions and the
    // year it occurs in.
    for( int year = mFirstTaxYear; year <= finalYearToCheck; ++year ) {
        double totalEmissions = mClimateModel->getEmissions( "CO2", year );
        double netOceanUp = mClimateModel->getNetOceanUptake( year );
        double netTerrUp = mClimateModel->getNetTerrestrialUptake( year );
        double totalNetUptake = netOceanUp + netTerrUp;
        double percentOff = objects::percentDiff( totalNetUptake, totalEmissions );
        if( percentOff < minDiff ) {
            minDiff = percentOff;
            minYear = year;
        }
    }
    return minYear;
}
