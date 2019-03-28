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
* \file ademand_function.cpp
* \ingroup Objects
* \brief The ADemandFunction class source file.
*
*  Detailed description.
*
* \author Pralit Patel
* \author Sonny Kim
* \author Josh Lurz
* \date $Date: 2005/06/01 21:23:59 $
* \version $Revision: 1.2 $
*/

#include "util/base/include/definitions.h"
#include "functions/include/ademand_function.h"
#include "functions/include/function_utils.h"

using namespace std;

/*! \brief Apply technical change to demand functions.
* \details 
* \note This function currently makes a call to
*       FunctionUtils::applyTechChangeInternal so that it can share code with
*       AProductionFunction::applyTechnicalChange. In the future these
*       implementations may vary so that function is not called directly.
* \param input Vector of inputs for the demand function.
* \param aTechChange A structure containing the various possible types of
*        technical change.
* \param regionName Name of the region containing the function.
* \param sectorName Name of the sector containing the function.
* \param alphaZero The up-front scaler.
* \param sigma Sigma coefficient.
* \return The new alpha zero.
*/
double ADemandFunction::applyTechnicalChange( InputSet& input, const TechChange& aTechChange,
                                              const string& regionName, const string& sectorName,
                                              const int aPeriod, double alphaZero, double sigma ) const 
{
    return FunctionUtils::applyTechnicalChangeInternal( input, aTechChange, regionName, sectorName,
                                                        aPeriod, alphaZero, sigma );
}
