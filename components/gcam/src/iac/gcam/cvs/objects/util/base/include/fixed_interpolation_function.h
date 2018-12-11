#ifndef _FIXED_INTERPOLATION_FUNCTION_H_
#define _FIXED_INTERPOLATION_FUNCTION_H_
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
 * \file fixed_interpolation_function.h
 * \ingroup Objects
 * \brief Header file for the FixedInterpolationFunction class.
 * \author Pralit Patel
 */
#include <xercesc/dom/DOMNode.hpp>
#include <string>

#include "util/base/include/iinterpolation_function.h"

class DataPoint;

/*!
 * \ingroup Objects
 * \brief A fixed interpolation function.
 * \details This function simply returns the y-value of the left data point
 *          regardless of the values given for the right data point as well as
 *          the x-value to interpolate for.  Note that no error-checking is
 *          done by this interpolation function.
 *          <b>XML specification for FixedInterpolationFunction</b>
 *          - XML name: \c interpolation-function
 *          - Contained by:
 *          - Parsing inherited from class: None.
 *          - Attributes:
 *              - \c name = fixed
 *                      The XML name attribute value which differentiates this
 *                      IInterpolationFunction from the others.
 *          - Elements:
 *
 * \author Pralit Patel
 * \author Sonny Kim
 */
class FixedInterpolationFunction : public IInterpolationFunction {
public:
    FixedInterpolationFunction();
    ~FixedInterpolationFunction();
    
    static const std::string& getXMLAttrNameStatic();
    
    // IInterpolationFunction methods
    virtual double interpolate( const DataPoint* aLeftPoint, const DataPoint* aRightPoint,
        const double aXValue ) const;
    
    // IParsable methods
    virtual bool XMLParse( const xercesc::DOMNode* aNode );

    // IRoundTrippable methods
    virtual void toInputXML( std::ostream& aOut, Tabs* aTabs ) const;
};

#endif // _FIXED_INTERPOLATION_FUNCTION_H_
