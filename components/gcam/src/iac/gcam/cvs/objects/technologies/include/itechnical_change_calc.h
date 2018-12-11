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

#ifndef _ITECHNICAL_CHANGE_CALC_H_
#define _ITECHNICAL_CHANGE_CALC_H_
#if defined(_MSC_VER)
#pragma once
#endif

/*! 
 * \file itechnical_change_calc.h
 * \ingroup Objects
 * \brief ITechnicalChangeCalc interface header file.
 * \author Josh Lurz
 */
#include <vector>
#include "util/base/include/istandard_component.h"

class IInput;
struct PreviousPeriodInfo;
class IFunction;

/*! 
 * \ingroup Objects
 * \brief This object is responsible for determining the quantity of technical
 *        change to apply to each input and to the production function on
 *        aggregate.
 * \details This object can be added on to Technologies so that they can apply
 *          technical change to inputs and to the production function. The
 *          technical change calculator is responsible for calculating and
 *          adjusting inputs for technical change. It is also responsible for
 *          calculating Hicks Neutral technical change, which applies to the
 *          entire production function.
 * \author Josh Lurz
*/
class ITechnicalChangeCalc : public IParsedComponent { 
public:
    // Clone operator must be declared explicitly even though it is inherited
    // from IStandardComponent so that the return type can be changed. Since
    // this class is a subtype of IStandardComponent, this is legal and referred
    // to as a covariant return type.
    virtual ITechnicalChangeCalc* clone() const = 0;
    
    /*!
     * \brief Complete the initialization of the technical change calculator.
     */
    virtual void completeInit() = 0;

    /*!
     * \brief Calculate the total Hick's neutral technical change and adjust the
     *        technology's inputs for the technical change.
     * \param aInputs List of inputs to adjust.
     * \param aProductionFunc The Technology's production function.
     * \param aPeriod Model period.
     * \return The cumulative Hicks neutral technical change.
     */
    virtual double calcAndAdjustForTechChange( std::vector<IInput*>& aInputs,
                                               PreviousPeriodInfo& aPreviousPeriodInfo,
                                               const IFunction* aProductionFunc,
                                               const std::string& aRegionName,
                                               const std::string& aSectorName,
                                               const int aPeriod ) const = 0;
};

#endif // _ITECHNICAL_CHANGE_CALC_H_
