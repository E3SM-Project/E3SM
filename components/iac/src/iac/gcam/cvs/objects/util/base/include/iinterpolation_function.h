#ifndef _IINTERPOLATION_FUNCTION_H_
#define _IINTERPOLATION_FUNCTION_H_
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
 * \file iinterpolation_function.h
 * \ingroup Objects
 * \brief Header file for the IInterpolationFunction interface.
 * \author Pralit Patel
 */
#include "util/base/include/iparsable.h"
#include "util/base/include/iround_trippable.h"

class DataPoint;

/*!
 * \ingroup Objects
 * \brief An interface for a function which can be used to interpolate a value
 *        between two data points.
 * \author Pralit Patel
 * \author Sonny Kim
 */
class IInterpolationFunction : public IParsable, public IRoundTrippable {
public:
    //! Virtual destructor so that instances of the interface may be deleted
    //! correctly through a pointer to the interface.
    inline virtual ~IInterpolationFunction();
    
    /*!
     * \brief Interpolate a y-value at the given x-value that is between the
     *        given left and right data points.
     * \param aLeftPoint The left data point to bracket this interpolation.
     * \param aRightPoint The right data point to bracket this interpolation.
     * \param aXValue The x value at which we want to interpolate.
     * \return The interpolated y-value at the given aXValue.
     */
    virtual double interpolate( const DataPoint* aLeftPoint, const DataPoint* aRightPoint,
        const double aXValue ) const = 0;

    /*!
     * \brief The the XML element name which will be shared by interpolation
     *        functions.
     * \details The XML attribute name will be used to differentiate the subclasses
     *          of IInterpolationRule.  This setup was choosen to ease generating
     *          the XML tags for interpolation functions.  See getXMLAttrNameStatic()
     *          from the subclasses.
     */
    static const std::string& getXMLNameStatic() {
        const static std::string XML_NAME = "interpolation-function";
        return XML_NAME;
    }
};

// Inline function definitions.
IInterpolationFunction::~IInterpolationFunction(){
}

#endif // _IINTERPOLATION_FUNCTION_H_
