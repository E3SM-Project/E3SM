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
 * \file interpolation_function_factory.cpp
 * \ingroup Objects
 * \brief InterpolationFunctionFactory class source file.
 * \author Pralit Patel
 */

#include "util/base/include/definitions.h"
#include <string>
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>

#include "util/base/include/interpolation_function_factory.h"
#include "util/logger/include/ilogger.h"

// IInterpolationFunction subclasses
#include "util/base/include/linear_interpolation_function.h"
#include "util/base/include/fixed_interpolation_function.h"
#include "util/base/include/s_curve_interpolation_function.h"

using namespace std;
using namespace xercesc;

/*!
 * \brief Returns whether this factory can create a function with the given xml
 *        name attribute value.
 * \param aXMLAttrNameValue The name attribute value of an xml element to check.
 * \return True if this factory has a function with the given name, false otherwise.
 * \note The list of known function here needs to be kept in sync with
 *       the ones found in createAndParseFunction.
 */
bool InterpolationFunctionFactory::hasInterpolationFunction( const string& aXMLAttrNameValue ) {
    return LinearInterpolationFunction::getXMLAttrNameStatic() == aXMLAttrNameValue
        || FixedInterpolationFunction::getXMLAttrNameStatic() == aXMLAttrNameValue
        || SCurveInterpolationFunction::getXMLAttrNameStatic() == aXMLAttrNameValue;
}

/*!
 * \brief Creates and parses the function with the given xml attribute name value.
 * \details Creates the function and calls XMLParse on it before returning it,
 *          if there are no known functions which match the given name null
 *          is returned.
 * \param aXMLAttrNameValue The name attribute value of an xml element to check.
 * \param aNode The xml which defines the function to be created.
 * \return The newly created and parsed function or null if given an unknown type.
 * \note The list of known functions here must be kept in sync with
 *       the ones found in hasInterpolationFunction.
 */
IInterpolationFunction* InterpolationFunctionFactory::createAndParseFunction( const string& aXMLAttrNameValue,
                                                                              const DOMNode* aNode )
{
    // make sure we know about this function
    if( !hasInterpolationFunction( aXMLAttrNameValue ) ) {
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::WARNING );
        mainLog << "Could not create unknown IInterpolationFunction: " << aXMLAttrNameValue << endl;
        return 0;
    }
    
    // create the requested function
    IInterpolationFunction* retFunction;
    if( LinearInterpolationFunction::getXMLAttrNameStatic() == aXMLAttrNameValue ) {
        retFunction = new LinearInterpolationFunction();
    }
    else if( FixedInterpolationFunction::getXMLAttrNameStatic() == aXMLAttrNameValue ) {
        retFunction = new FixedInterpolationFunction();
    }
    else if( SCurveInterpolationFunction::getXMLAttrNameStatic() == aXMLAttrNameValue ) {
        retFunction = new SCurveInterpolationFunction();
    }
    else {
        // this must mean createAndParseFunction and hasInterpolationFunction
        // are out of sync with known functions
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::WARNING );
        mainLog << "Could not create unknown IInterpolationFunction: " << aXMLAttrNameValue
            << ", createAndParseFunction may be out of sync with hasInterpolationFunction." << endl;
        return 0;
    }
    
    // parse the created function
    if( aNode ) {
        retFunction->XMLParse( aNode );
    }
    return retFunction;
}
