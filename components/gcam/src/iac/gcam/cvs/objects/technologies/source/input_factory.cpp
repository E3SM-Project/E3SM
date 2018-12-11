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
 * \file input_factory.cpp
 * \ingroup Objects
 * \brief InputFactory source file.
 * \author Josh Lurz
 */

#include "util/base/include/definitions.h"
#include <string>
#include "technologies/include/input_factory.h"
#include "functions/include/iinput.h"
#include "util/logger/include/ilogger.h"

// Add new types here.
#include "functions/include/energy_input.h"
#include "functions/include/non_energy_input.h"
#include "functions/include/input_capital.h"
#include "functions/include/input_OM_fixed.h"
#include "functions/include/input_OM_var.h"
#include "functions/include/renewable_input.h"
#include "functions/include/building_demand_input.h"
#include "functions/include/input_subsidy.h"
#include "functions/include/input_tax.h"

using namespace std;

/*!
 * \brief Returns whether the requested type is a type the factory knows how to
 *        create.
 * \param aType Type to determine if the factory can create.
 * \return Whether the factory can create the type.
 */
bool InputFactory::isOfType( const std::string& aType ) {
    if( aType == EnergyInput::getXMLNameStatic() ){
        return true;
    }
    if( aType == NonEnergyInput::getXMLNameStatic() ){
        return true;
    }
    if( aType == RenewableInput::getXMLNameStatic() ){
        return true;
    }
    if( aType == InputCapital::getXMLNameStatic() ){
        return true;
    }
    if( aType == InputOMFixed::getXMLNameStatic() ){
        return true;
    }
    if( aType == InputOMVar::getXMLNameStatic() ){
        return true;
    }
    if( aType == BuildingDemandInput::getXMLNameStatic() ) {
        return true;
    }
    if( aType == InputSubsidy::getXMLNameStatic() ) {
        return true;
    }
    if( aType == InputTax::getXMLNameStatic() ) {
        return true;
    }
    return false;
}

/*!
 * \brief Return a new instance of an input of the requested type.
 * \param aType Type of input to return.
 * \return A newly created input wrapped in an auto_ptr. The pointer is null if
 *         the type is unknown.
 */
auto_ptr<IInput> InputFactory::create( const std::string& aType ) {
    if( aType == EnergyInput::getXMLNameStatic() ){
        return auto_ptr<IInput>( new EnergyInput );
    }
    if( aType == NonEnergyInput::getXMLNameStatic() ){
        return auto_ptr<IInput>( new NonEnergyInput );
    }
    if( aType == InputCapital::getXMLNameStatic() ){
        return auto_ptr<IInput>( new InputCapital );
    }
    if( aType == InputOMFixed::getXMLNameStatic() ){
        return auto_ptr<IInput>( new InputOMFixed );
    }
    if( aType == InputOMVar::getXMLNameStatic() ){
        return auto_ptr<IInput>( new InputOMVar );
    }
    if( aType == RenewableInput::getXMLNameStatic() ){
        return auto_ptr<IInput>( new RenewableInput );
    }
    if( aType == BuildingDemandInput::getXMLNameStatic() ){
        return auto_ptr<IInput>( new BuildingDemandInput );
    }
    if( aType == InputSubsidy::getXMLNameStatic() ){
        return auto_ptr<IInput>( new InputSubsidy );
    }
    if( aType == InputTax::getXMLNameStatic() ){
        return auto_ptr<IInput>( new InputTax );
    }

    // Check for consistency.
    assert( !isOfType( aType ) );

    // Unknown type.
    ILogger& mainLog = ILogger::getLogger( "main_log" );
    mainLog.setLevel( ILogger::ERROR );
    mainLog << "Could not create input of type " << aType << "." << endl;
    return auto_ptr<IInput>();
}
