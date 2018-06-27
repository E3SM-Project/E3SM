#ifndef _VARIABLE_COST_SHUTDOWN_DECIDER_H_
#define _VARIABLE_COST_SHUTDOWN_DECIDER_H_
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
* \file variable_cost_shutdown_decider.h
* \ingroup Objects
* \brief The VariableCostShutdownDecider header file.
* \author Josh Lurz
* \date $Date: 2005/06/03 22:11:39 $
* \version $Revision: 1.2 $
*/
#include "technologies/include/ishutdown_decider.h"

#include <string>
struct ProductionFunctionInfo;
/*! 
* \ingroup Objects
* \brief This object makes the shutdown decision for a vintage based on its
*        variable cost exceeding its price received.
* \details This object makes the shutdown decision for a vintage by determining
*          the point at which its variable costs, defined as the costs per unit
*          of all inputs other than capital, are equal to the price received for
*          the product. This is the classic short-term shutdown decision. This
*          calculation is done as a ratio with the denominator as the capital
*          for the vintage, and when the vintages reaches a minimum level the
*          vintage is smoothly shutdown.
* \author Josh Lurz
*/
class VariableCostShutdownDecider: public IShutdownDecider
{
public:
    VariableCostShutdownDecider();
	VariableCostShutdownDecider* clone() const;
    double calcShutdownCoef( const ProductionFunctionInfo* aFuncInfo,
							 const double aCalculatedProfits,
                             const std::string& aRegionName,
                             const std::string& aSectorName,
                             const int aPeriod ) const;
};


#endif // _VARIABLE_COST_SHUTDOWN_DECIDER_H_
