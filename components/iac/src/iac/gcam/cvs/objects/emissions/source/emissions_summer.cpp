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
* \file emissions_summer.cpp
* \ingroup Objects
* \brief The EmissionsSummer class source file.
*
* \author Josh Lurz
*/

#include "util/base/include/definitions.h"
#include <cassert>
#include "emissions/include/emissions_summer.h"
#include "emissions/include/aghg.h"

using namespace std;

/*! \brief Constructor
* \param aGHG GHG that is being summed.
*/
EmissionsSummer::EmissionsSummer( const string& aGHGName ):
mGHGName( aGHGName ){
}

/*! \brief Add emissions from a GHG to the stored emissions.
* \param aGHG GHG to update emissions from.
* \param aPeriod Period in which to update.
*/
void EmissionsSummer::startVisitGHG( const AGHG* aGHG, const int aPeriod ){
    if( aGHG->getName() == mGHGName ){
        mEmissionsByPeriod[ aPeriod ] += aGHG->getEmission( aPeriod );
    }
}

/*! \brief Get the current emissions sum.
* \param aPeriod Model period for which to get emissions.
* \return The emissions sum.
*/
double EmissionsSummer::getEmissions( const int aPeriod ) const {
    // The value may not be initialized if there were no GHGs.  
    // The default zero will be correct though.

    // The emissions sum may be negative if uptake is occurring.
    return mEmissionsByPeriod[ aPeriod ];
}

/*! \brief Return whether any emissions were set for the period.
* \param aPeriod Model period.
* \return Whether any emissions were set.
*/
double EmissionsSummer::areEmissionsSet( const int aPeriod ) const {
    return mEmissionsByPeriod[ aPeriod ].isInited();
}

/*! \brief Get the GHG name.
 * \return The name of the GHG that is summed by this object.
 */
const string& EmissionsSummer::getGHGName() const {
    return mGHGName;
}

/*!
 * \brief Add an EmissionsSummer to the group.
 * \details The given EmissionsSummer will be updated for all model periods.  The
 *          memory for the given EmissionsSummer will not be managed by this object.
 * \param A refernce to an EmissionsSummer to update when this group is updated.
 */
void GroupedEmissionsSummer::addEmissionsSummer( EmissionsSummer* aEmissionsSummer ) {
    mEmissionsSummers[ aEmissionsSummer->getGHGName() ] = aEmissionsSummer;
}

void GroupedEmissionsSummer::startVisitGHG( const AGHG* aGHG, const int aPeriod ) {
    // We are currently assuming all periods should be updated.
    assert( aPeriod == -1 );
    
    CSummerIterator it = mEmissionsSummers.find( aGHG->getName() );
    if( it != mEmissionsSummers.end() ) {
        for( int period = 1; period < scenario->getModeltime()->getmaxper(); ++period ) {
            (*it).second->startVisitGHG( aGHG, period );
        }
    }
}
