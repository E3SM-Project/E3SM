#ifndef _FORCING_TARGET_H_
#define _FORCING_TARGET_H_
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
* \file forcing_target.h
* \ingroup Objects
* \brief The ForcingTarget class header file.
* \author Josh Lurz
*/

#include "target_finder/include/itarget.h"
#include <string>

class IClimateModel;

/*! \brief A class representing a forcing target.
* \details A target for a total climate forcing for a set of gases.
*/
class ForcingTarget: public ITarget {
public:
    ForcingTarget( const IClimateModel* aClimateModel,
                   const double aTargetValue,
                   const int aFirstTaxYear );
    
    static const std::string& getXMLNameStatic();

    // ITarget methods
    virtual double getStatus( const int aYear ) const;
    
    virtual int getYearOfMaxTargetValue() const;
private:
    //! The climate model.
    const IClimateModel* mClimateModel;

    //! The target value.
    const double mTargetValue;
    
    //! The first policy year which would be the first valid year to check
    //! getStatus in.
    const int mFirstTaxYear;
};

#endif // _FORCING_TARGET_H_
