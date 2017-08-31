#ifndef _INTERPOLATION_FUNCTION_FACTORY_H_
#define _INTERPOLATION_FUNCTION_FACTORY_H_
#if defined(_MSC_VER)
#pragma once
#endif

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
 * \file interpolation_function_factory.h
 * \ingroup Objects
 * \brief Header file for the InterpolationFunctionFactory class.
 * \author Pralit Patel
 */
#include <xercesc/dom/DOMNode.hpp>
#include <string>

class IInterpolationFunction;

/*!
 * \ingroup Objects
 * \brief A factory which can be used to create instances of an interpolation function.
 * \details There are two static methods to determine if this factory can create an
 *          interpolation function.  The XML setup for an interpolation function is
 *          a little different so that generating the XML tags would be easier.  The
 *          node name will always be "interpolation-function" and the name attribute
 *          will be used to determine which function should be created by this factory.
 *
 * \author Pralit Patel
 * \author Sonny Kim
 */
class InterpolationFunctionFactory {
public:
    static bool hasInterpolationFunction( const std::string& aXMLAttrNameValue );
    
    static IInterpolationFunction* createAndParseFunction( const std::string& aXMLAttrNameValue,
        const xercesc::DOMNode* aNode );
};

#endif // _INTERPOLATION_FUNCTION_FACTORY_H_
