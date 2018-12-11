#ifndef _S_CURVE_INTERPOLATION_FUNCTION_H_
#define _S_CURVE_INTERPOLATION_FUNCTION_H_
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
 * \file s_curve_interpolation_function.h
 * \ingroup Objects
 * \brief Header file for the SCurveInterpolationFunction class.
 * \author Pralit Patel
 */
#include <xercesc/dom/DOMNode.hpp>
#include <string>

#include "util/base/include/iinterpolation_function.h"

class DataPoint;

/*!
 * \ingroup Objects
 * \brief A S-curve interpolation function.
 * \details Interpolate a value by constructing an S-curve which connects
 *          the left and right data points and evaluating that curve at the
 *          given x-value.  Tuning parameters to shape the S-curve are parsed.
 *          TODO: what do we do about errors?
 *          <b>XML specification for SCurveInterpolationFunction</b>
 *          - XML name: \c interpolation-function
 *          - Contained by:
 *          - Parsing inherited from class: None.
 *          - Attributes:
 *              - \c name = s-curve
 *                      The XML name attribute value which differentiates this
 *                      IInterpolationFunction from the others.
 *          - Elements:
 *              - \c steepness double SCurveInterpolationFunction::mSteepness
 *                      Adjustment parameter to control steepness of the s-curve,
 *                      the higher the number the steeper the curve.
 *              - \c median-x-value double SCurveInterpolationFunction::mMedianXValue
 *                      Adjustment parameter to control at which x-value will the s-curve
 *                      equal the median of the y-values from the left and right
 *                      data points when interpolating.
 *
 * \author Pralit Patel
 * \author Sonny Kim
 */
class SCurveInterpolationFunction : public IInterpolationFunction {
public:
    SCurveInterpolationFunction();
    SCurveInterpolationFunction( const double aSteepness,
        const double aMedianXValue );
    ~SCurveInterpolationFunction();
    
    static const std::string& getXMLAttrNameStatic();
    
    // IInterpolationFunction methods
    virtual double interpolate( const DataPoint* aLeftPoint, const DataPoint* aRightPoint,
        const double aXValue ) const;
    
    // IParsable methods
    virtual bool XMLParse( const xercesc::DOMNode* aNode );

    // IRoundTrippable methods
    virtual void toInputXML( std::ostream& aOut, Tabs* aTabs ) const;

private:
    //! Parameter for the steepness of the s-curve. Higher number means steeper ascent.
    double mSteepness;

    //! Parameter for at which x-value the curve will have exactly the median of
    //! the y-values from the left and right data points.
    double mMedianXValue;
};

#endif // _S_CURVE_INTERPOLATION_FUNCTION_H_
