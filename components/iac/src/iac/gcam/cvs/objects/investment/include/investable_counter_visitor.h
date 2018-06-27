#ifndef _INVESTABLE_COUNTER_VISITOR_H_
#define _INVESTABLE_COUNTER_VISITOR_H_
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
* \file investable_counter_visitor.h
* \ingroup Objects
* \brief InvestableCounterVisitor class header file.
* \author Sonny Kim
*/

#include "util/base/include/default_visitor.h"
#include <vector>
#include <string>

/*! 
* \ingroup Objects
* \brief A visitor which visits all subsectors and technologies for each sector 
*          for determining the number of calibration markets to create.
* \details This visitor when passed to a Region will tabulate individually the
*          fixed and calibrated supplies of all goods in the region. It will
*          also keep track of whether there are exist variable supplies of each
*          good in the region. This should be used instead of
*          Sector::getCalOutput if secondary outputs need to be included in the
*          sums. The visitor should be used by passing this object as an
*          argument to the accept function, and then access the information
*          through the getSupplyInfo function.  This Visitor will be responsible
*          for creating the calibration markets and setting them to solve.  Note 
*          that an investable that is 100% fixed investment will not be set to
*          solve.  InvestableCounterVisitor will use a naming convention 
*          "SubSec-"+Subsector::getName() or "Tech-"+BaseTechnology::getName().
* \warning When subsectors or technologies have the same name accross sectors/subsectors
*          this calibration naming scheme breaks.  Using the full path in the naming scheme
*          is a somewhat better solution.
* \see     GetDistributedInvestmentVisitor
* \see     SetShareWeightVisitor
* \see     MarketBasedInvestor
*
* \author Sonny Kim
*/
class InvestableCounterVisitor : public DefaultVisitor {
public:

    InvestableCounterVisitor( const std::string& aRegionName, const int aInvestableCount );

    virtual void startVisitSubsector( const Subsector* aSubsector,
                                      const int aPeriod );

    virtual void endVisitSubsector( const Subsector* aSubsector,
                                    const int aPeriod );

    virtual void startVisitBaseTechnology( const BaseTechnology* aTechnology,
                                       const int aPeriod );

    virtual void endVisitBaseTechnology( const BaseTechnology* aTechnology,
                                     const int aPeriod );
    
private:
    //! Name of the Region the calibration and fixed outputs are currently being
    //! summed for.
    std::string mCurrentRegionName;

    //! Name of the subsector currently being visited.
    std::string mCurrentSubSectorName;

    //! A counter for the numer of subsectors that have been processed
    //! used to determine which to hold at a constant share weight
    int mSubSectorCount;

    //! A counter for the number of technologies that have been processed
    //! used to determine which to hold at a constant share weight
    int mTechCount;

    //! A flag which will be determined during construction
    //! of this object.  If it is false then no calibration markets
    //! are necessary at the subsector level.
    const bool mMultipleSubsectors;

    //! A flag set at the subsector level to determine if calibration markets
    //! are necessary at the technology level.
    bool mMultipleTechnologies;
};

#endif // _INVESTABLE_COUNTER_VISITOR_H_
