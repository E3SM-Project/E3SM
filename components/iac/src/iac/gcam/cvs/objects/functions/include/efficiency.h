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
#ifndef _EFFICIENCY_H_
#define _EFFICIENCY_H_
#if defined(_MSC_VER)
#pragma once
#endif

/*! 
 * \file efficiency.h
 * \ingroup Objects
 * \brief Efficiency class header file.
 * \author Josh Lurz
 */

#include <string>
#include "util/base/include/value.h"
#include "functions/include/icoefficient.h"
class Tabs;

/*! 
 * \ingroup Objects
 * \brief Defines an efficiency for an input in a production function.
 * \details The efficiency for an input in a production function is equal to the
 *          amount of output produced per unit of input(output/input). This
 *          class defines an efficiency which is internally converted to an
 *          intensity for use in production functions. Efficiencies are intended
 *          to more easily represent single energy-input production functions.
 *          Coefficients(intensities) should be used for multi-energy-input
 *          production functions.
 * \author Josh Lurz
 */
class Efficiency: public ICoefficient { 
public:
    const static std::string& getXMLNameStatic();

    Efficiency( const double aEfficiency );

    virtual Efficiency* clone() const;

    virtual bool isSameType( const std::string& aType ) const;

    virtual const std::string& getName() const;

    virtual void toInputXML( std::ostream& aOut,
                             Tabs* aTabs ) const;

    virtual void toDebugXML( const int aPeriod,
                             std::ostream& aOut,
                             Tabs* aTabs ) const;

    virtual void completeInit();

    virtual double getCoefficient() const;
private:
    //! The read-in efficiency.
    Value mReadInEfficiency;
};

#endif // _EFFICIENCY_H_
