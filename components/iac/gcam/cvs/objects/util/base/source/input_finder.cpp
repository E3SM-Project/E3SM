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
* \file input_finder.cpp
* \ingroup Objects
* \brief The InputFinder class source file.
*
* \author Josh Lurz
*/

#include "util/base/include/definitions.h"
#include <cassert>
#include <algorithm>
#include "util/base/include/input_finder.h"
#include "technologies/include/technology.h"
#include "functions/include/minicam_input.h"
#include "sectors/include/afinal_demand.h"

using namespace std;

/*! \brief Constructor
*/
InputFinder::InputFinder(){
}

/*!
* \brief Add an input used by the technology to the list of inputs.
* \details Adds the technology's input to the list of inputs if it is unique.
* \param aTechnology Technology from which to update inputs.
* \param aPeriod Period in which to update.
*/
void InputFinder::startVisitMiniCAMInput( const MiniCAMInput* aInput, const int aPeriod ){
    // TODO: Add other types of inputs.
    if( !aInput->hasTypeFlag( IInput::ENERGY ) ){
        return;
    }

    // Check if the input is already known.
    const string input = aInput->getName();
    if( !input.empty() && find( mInputs.begin(), mInputs.end(), input ) == mInputs.end() ){
        mInputs.push_back( input );
    }
}

/*!
* \brief Adds the single input of a final demand to the list of inputs.
* \details Final demands have a single input which is the good they demand. This is the same as the final demand name.
* \param aFinalDemand Final demand from which to update inputs.
* \param aPeriod Period in which to update.
*/
void InputFinder::startVisitFinalDemand( const AFinalDemand* aFinalDemand, const int aPeriod ){
    // Check if the input is already known.
    const string input = aFinalDemand->getName();
    if( !input.empty() && find( mInputs.begin(), mInputs.end(), input ) == mInputs.end() ){
        mInputs.push_back( input );
    }
}


/*! \brief Get the list of inputs used.
* \pre The InputFinder has already visited technologies.
* \return The list of inputs used.
*/
const list<string>& InputFinder::getInputs() const {
    return mInputs;
}

void InputFinder::finish() const {
}
