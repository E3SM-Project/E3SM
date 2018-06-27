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
 * \file indirect_emissions_calculator.cpp
 * \ingroup Objects
 * \brief The IndirectEmissionsCalculator source file.
 * \author Josh Lurz
 */

#include "util/base/include/definitions.h"
#include <cassert>
#include "reporting/include/indirect_emissions_calculator.h"
#include "sectors/include/sector.h"
#include "technologies/include/technology.h"
#include "util/logger/include/ilogger.h"
#include "util/base/include/util.h"

// Remove when modeltime is no longer required.
#include "containers/include/scenario.h"
#include "util/base/include/model_time.h"
extern Scenario* scenario;

using namespace std;

/*! 
 * \brief Constructor
 */
IndirectEmissionsCalculator::IndirectEmissionsCalculator():
mCurrTotalEmissions( 0 ),
mCurrIndirectEmissions( 0 ),
mCurrOutput( 0 ){
}

/*!
 * \brief Get the upstream emissions coefficient for a Sector.
 * \details Returns the upstream emissions coefficient for the Sector. The
 *          upstream emissions coefficient includes both direct emissions and
 *          upstream emissions.
 * \param aSector Name of the sector for which to get the upstream emissions
 *        coefficient.
 * \return The upstream emissions coefficient for a sector.
 */
double IndirectEmissionsCalculator::getUpstreamEmissionsCoefficient( const string& aSector,
                                                                     const int aPeriod ) const
{
    DoubleMap::const_iterator upstreamEmissionsCoef = mUpstreamEmissionsCoefficients.find( aSector );
    if( mUpstreamEmissionsCoefficients.end() == upstreamEmissionsCoef
        || !upstreamEmissionsCoef->second[ aPeriod ].isInited() )
    {
        // Upstream emissions coefficients do not exist for resources, renewable
        // fuels and the 'none' fuel.
        return 0;
    }
    return upstreamEmissionsCoef->second[ aPeriod ];
}

/*!
 * \brief Get the total indirect emissions for a Sector.
 * \details Returns the total indirect emissions for the Sector. The indirect
 *          emissions include all upstream emissions, but does not include
 *          direct emissions for the sector.
 * \param aSector Name of the sector for which to get the total indirect
 *        emissions.
 * \return The total indirect emissions.
 * \todo Remove this function once the Access database is removed.
 */
double IndirectEmissionsCalculator::getIndirectEmissions( const string& aSector,
                                                          const int aPeriod ) const
{
    DoubleMap::const_iterator indirectEmissions = mIndirectEmissions.find( aSector );
    if( mIndirectEmissions.end() == indirectEmissions || !indirectEmissions->second[ aPeriod ].isInited() ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::WARNING );
        mainLog << "No indirect emissions were calculated for sector " << aSector << "." << endl;
        return 0;
    }
    return indirectEmissions->second[ aPeriod ];
}

void IndirectEmissionsCalculator::startVisitSector( const Sector* aSector,
                                                    const int aPeriod )
{
    // Check that aPeriod is not all period mode.
    assert( aPeriod != -1 );

    // Ensure all current variables have been cleared.
    assert( mCurrSectorName.empty() && mCurrTotalEmissions == 0 && mCurrIndirectEmissions == 0
            && mCurrOutput == 0 );

    // Ensure values have not already been calculated.
    assert( mUpstreamEmissionsCoefficients.find( mCurrSectorName ) == mUpstreamEmissionsCoefficients.end() );
    assert( mIndirectEmissions.find( mCurrSectorName ) == mIndirectEmissions.end() );

    // Set the current sector name.
    mCurrSectorName = aSector->getName();
}

void IndirectEmissionsCalculator::endVisitSector( const Sector* aSector, const int aPeriod ){
    // Check that aPeriod is not all period mode.
    assert( aPeriod != -1 );

    // Ensure the sector name has been set.
    assert( !mCurrSectorName.empty() );

    // Calculate the indirect emissions coefficient.
    mUpstreamEmissionsCoefficients[ mCurrSectorName ][ aPeriod ] = mCurrOutput > util::getSmallNumber()
                                                                   ? mCurrTotalEmissions / mCurrOutput : 0;

    // Store the indirect emissions.
    mIndirectEmissions[ mCurrSectorName ][ aPeriod ] = mCurrIndirectEmissions;

    // Clear the state variables.
    mCurrSectorName.clear();
    mCurrTotalEmissions = 0;
    mCurrIndirectEmissions = 0;
    mCurrOutput = 0;
}

void IndirectEmissionsCalculator::startVisitTechnology( const Technology* aTechnology,
                                                        const int aPeriod )
{
    // TODO: Fix this for vintaging.
    if( scenario->getModeltime()->getper_to_yr( aPeriod ) != aTechnology->year ){
        return;
    }

    // Ensure the sector name has been set.
    assert( !mCurrSectorName.empty() );

    // TODO: FIX THIS CLASS FOR MULTI-INPUT.
    double upstreamCoef = 0; // getUpstreamEmissionsCoefficient( aTechnology->getFuelName(), aPeriod );

    // Calculate indirect and total emissions. Total emissions are required to
    // calculate the indirect emissions coefficient for this sector which will
    // include direct emissions.
    // TODO: FIX THIS
    double indirectEmissions = upstreamCoef; // * aTechnology->getInput( aPeriod );

    double totalEmissions = aTechnology->getEmissionsByGas( "CO2", aPeriod )
                            + indirectEmissions;

    // Add to the regional lists.
    mCurrIndirectEmissions += indirectEmissions;
    mCurrTotalEmissions += totalEmissions;
    mCurrOutput += aTechnology->getOutput( aPeriod );
}
