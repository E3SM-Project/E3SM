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

#ifndef _INTENSITY_H_
#define _INTENSITY_H_
#if defined(_MSC_VER)
#pragma once
#endif

/*! 
 * \file intensity.h
 * \ingroup Objects
 * \brief Intensity class header file.
 * \author Josh Lurz
*/

#include <string>
#include "util/base/include/value.h"
#include "functions/include/icoefficient.h"

class Tabs;

/*! 
 * \ingroup Objects
 * \brief Represents an intensity.
 * \details The intensity for an input is equal to the amount required to
 *          produce one unit of output. This class defines an intensity for use
 *          in a production function.
 * \author Josh Lurz
 */
class Intensity : public ICoefficient { 
public:
    const static std::string& getXMLNameStatic();

    Intensity( const double aIntensity );

    virtual Intensity* clone() const;
    
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
    //! The read-in intensity.
    Value mReadInIntensity;
};

#endif // _INTENSITY_H_
