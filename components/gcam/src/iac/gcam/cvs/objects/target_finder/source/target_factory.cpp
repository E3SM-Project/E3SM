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
* \file target_factory.cpp
* \ingroup Objects
* \brief TargetFactory source file.
* \author Josh Lurz
*/

#include "util/base/include/definitions.h"
#include <string>
#include <cassert>

#include "target_finder/include/target_factory.h"
#include "util/logger/include/ilogger.h"

// Add new types here.
#include "target_finder/include/concentration_target.h"
#include "target_finder/include/forcing_target.h"
#include "target_finder/include/temperature_target.h"
#include "target_finder/include/emissions_stabalization_target.h"
#include "target_finder/include/kyoto_forcing_target.h"

using namespace std;

/*! \brief Returns whether the requested type is a type the factory knows how to
*          create.
* \param aType Type to determine if the factory can create.
* \return Whether the factory can create the type.
*/
bool TargetFactory::isOfType( const string& aType ) {
	// Search the list of known types.
	return ( ( aType == ConcentrationTarget::getXMLNameStatic() )
		|| ( aType == ForcingTarget::getXMLNameStatic() )
		|| ( aType == TemperatureTarget::getXMLNameStatic() )
        || ( aType == KyotoForcingTarget::getXMLNameStatic() )
        || ( aType == EmissionsStabalizationTarget::getXMLNameStatic() ) );
}

/*!
 * \brief Return a new instance of a component of the requested type.
 * \param aType Type of ITarget to return.
 * \param aClimateModel Scenario's climate model.
 * \param aTargetValue The target value.
 * \param aFirstTaxYear The first year in which the target could be checked.
 * \return A newly created ITarget wrapped in an auto_ptr. The pointer
 *         is null if the type is unknown.
 */
auto_ptr<ITarget> TargetFactory::create( const string& aType,
                                         const IClimateModel* aClimateModel,
                                         double aTargetValue,
                                         int aFirstTaxYear )
{
	// Search the list of known types.
	if( aType == ConcentrationTarget::getXMLNameStatic() ) {
		return auto_ptr<ITarget>( new ConcentrationTarget( aClimateModel,
                                                           aTargetValue,
                                                           aFirstTaxYear ) );
	}
	if( aType == ForcingTarget::getXMLNameStatic() ){
		return auto_ptr<ITarget>( new ForcingTarget( aClimateModel,
                                                     aTargetValue,
                                                     aFirstTaxYear ) );
	}
	if( aType == TemperatureTarget::getXMLNameStatic() ){
		return auto_ptr<ITarget>( new TemperatureTarget( aClimateModel,
                                                         aTargetValue,
                                                         aFirstTaxYear ) );
	}
    
    if( aType == KyotoForcingTarget::getXMLNameStatic() ){
        return auto_ptr<ITarget>( new KyotoForcingTarget( aClimateModel,
                                                          aTargetValue,
                                                          aFirstTaxYear ) );
    }

    if( aType == EmissionsStabalizationTarget::getXMLNameStatic() ){
        return auto_ptr<ITarget>( new EmissionsStabalizationTarget( aClimateModel,
                                                                    aTargetValue,
                                                                    aFirstTaxYear ) );
    }

    // Make sure this create is in sync with isOfType.
    assert( !isOfType( aType ) );

	// Unknown type.
	ILogger& mainLog = ILogger::getLogger( "main_log" );
	mainLog.setLevel( ILogger::ERROR );
	mainLog << "Could not create Target of type " << aType << "." << endl;
	return auto_ptr<ITarget>();
}
