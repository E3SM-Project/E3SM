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
* \file calc_capital_good_price_visitor.cpp
* \ingroup Objects
* \brief The CalcCapitalGoodPriceVisitor class source file.
* \author Pralit Patel
*/

#include "util/base/include/definitions.h"
#include <cassert>
#include "consumers/include/calc_capital_good_price_visitor.h"
#include "consumers/include/invest_consumer.h"
#include "containers/include/scenario.h"
#include "util/base/include/model_time.h"
#include "functions/include/node_input.h"
#include "functions/include/function_utils.h"

using namespace std;

extern Scenario* scenario;

CalcCapitalGoodPriceVisitor::CalcCapitalGoodPriceVisitor( std::string& aRegionName )
:mRegionName( aRegionName )
{
}

void CalcCapitalGoodPriceVisitor::startVisitInvestConsumer( const InvestConsumer* aInvestConsumer, 
                                                                const int aPeriod )
{
    const Modeltime* modelTime = scenario->getModeltime();
    // only calc for the current year consumer
    if( aInvestConsumer->year == modelTime->getper_to_yr( aPeriod ) ) {
        // calculate the price of the capital good which is the CES aggregate of
        // the investment consumer's inputs
        aInvestConsumer->calcPricePaid( 0, mRegionName, "", aPeriod, 
                    modelTime->gettimestep( aPeriod ) );
        aInvestConsumer->mNestedInputRoot->calcLevelizedCost( mRegionName, "", aPeriod, 
                aInvestConsumer->mNestedInputRoot->getCoefficient( aPeriod ) ); // alpha zero is the root's alpha
        
        // Set the price of the capital good so the production sectors can get to it
        FunctionUtils::setCapitalGoodPrice( mRegionName, aPeriod, 
            aInvestConsumer->mNestedInputRoot->getPricePaid( mRegionName, aPeriod ) );
    }
}
