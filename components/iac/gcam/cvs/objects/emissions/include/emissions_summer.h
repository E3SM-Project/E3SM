#ifndef _EMISSIONS_SUMMER_H_
#define _EMISSIONS_SUMMER_H_
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
* \file emissions_summer.h
* \ingroup Objects
* \brief EmissionsSummer class header file.
* \author Josh Lurz
*/

#include "util/base/include/time_vector.h"
#include "util/base/include/default_visitor.h"
#include "util/base/include/value.h"

/*! 
* \ingroup Objects
* \brief A class which sums emissions for a particular gas.
* \details 
* \author Josh Lurz
*/

class EmissionsSummer : public DefaultVisitor {
public:
    explicit EmissionsSummer( const std::string& aGHGName );

    virtual void startVisitGHG( const AGHG* aGHG,
                                const int aPeriod );

    // Non-IVisitor interface methods.
    double getEmissions( const int aPeriod ) const;

    double areEmissionsSet( const int aPeriod ) const;
    
    const std::string& getGHGName() const;
private:
    //! The name of the GHG being summed.
    const std::string mGHGName;

    //! The current sum.
    objects::PeriodVector<Value> mEmissionsByPeriod;
};

/*!
 * \brief A container for EmissionsSummer objects so that the model can only be
 *        visited once and have all of the EmissionsSummers updated.
 * \details Visitors tend to be slow performance wise due to the fact that they
 *          have to visit all objects within the model.  This becomes noticeable
 *          when a large number of EmissionsSummers need to be updated for all
 *          model periods.  Using this class to visit just once for all gasses
 *          and all periods drastically reducing runtime.
 * \note Unlike the EmissionsSummer visitor this class assumes it is visited with
 *       the aPeriod -1 flag to indicate that it will be updated for all periods.
 *       Although this behavior could easily be changed.
 * \author Pralit Patel
 */
class GroupedEmissionsSummer : public DefaultVisitor {
public:
    void addEmissionsSummer( EmissionsSummer* aEmissionsSummer );
    
    // DefaultVisitor methods
    virtual void startVisitGHG( const AGHG* aGHG,
                                const int aPeriod );
    
    /*virtual void startVisitCarbonCalc( const ICarbonCalc* aCarbonCalc,
                                       const int aPeriod );*/
    
private:
    //! A map of emissions summer by GHG name.  The memory for the EmissionsSummer
    //! is not managed by this class.
    std::map<std::string, EmissionsSummer*> mEmissionsSummers;
    
    typedef std::map<std::string, EmissionsSummer*>::const_iterator CSummerIterator;
};

#endif // _EMISSIONS_SUMMER_H_
