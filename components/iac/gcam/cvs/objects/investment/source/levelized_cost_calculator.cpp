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
 * \file levelized_cost_calculator.cpp
 * \ingroup Objects
 * \brief LevelizedCostCalculator class source file.
 * \author Josh Lurz
 */

#include "util/base/include/definitions.h"
#include <cassert>
#include <cmath>
#include "investment/include/levelized_cost_calculator.h"
#include "util/base/include/util.h"
#include "investment/include/investment_utils.h"
#include "investment/include/iinvestable.h"
#include "functions/include/ifunction.h"
#include "containers/include/national_account.h"
#include "functions/include/iinput.h"
#include "functions/include/inested_input.h"
#include "functions/include/function_utils.h"
#include "util/base/include/xml_helper.h"

using namespace std;

//! Constructor
LevelizedCostCalculator::LevelizedCostCalculator(){
}

/*! \brief Return the XML name for this object as a string.
* \return The XML name for the object.
*/
const string& LevelizedCostCalculator::getXMLNameStatic(){
    const static string XML_NAME = "levelized-cost-calculator";
    return XML_NAME;
}

/*! \brief Write the object to an XML output stream.
* \param aOut The output stream to write to.
* \param aTabs The object which tracks the number of tabs to write.
*/
void LevelizedCostCalculator::toInputXML( ostream& aOut, Tabs* aTabs ) const {
    XMLWriteOpeningTag( getXMLNameStatic(), aOut, aTabs );
    XMLWriteClosingTag( getXMLNameStatic(), aOut, aTabs );
}

/*! \brief Write the object to an XML output stream for debugging.
* \param aPeriod Period to write debugging information for.
* \param aOut The output stream to write to.
* \param aTabs The object which tracks the number of tabs to write.
*/
void LevelizedCostCalculator::toDebugXML( const int aPeriod, ostream& aOut, Tabs* aTabs ) const {
    XMLWriteOpeningTag( getXMLNameStatic(), aOut, aTabs );
    XMLWriteClosingTag( getXMLNameStatic(), aOut, aTabs );
}

/*! \brief Calculate the parent level levelized cost.
* \details TODO
* \param aInvestables Children to for which to calculate a levelized cost.
* \param aNationalAccount National accounts container.
* \param aRegionName Name of the region in which investment is occurring.
* \param aGoodName Name of the sector in which investment is occurring.
* \param aInvestmentLogitExp The investment logit exponential.
* \param aIsShareCalc Whether this expected profit rate is being used to
*        calculate shares.
* \param aIsDistributing Whether this expected profit rate is being used
*        to distribute investment.
* \param aPeriod The period in which to calculate the expected profit rate.
* \return The parent level levelized cost.
*/
double LevelizedCostCalculator::calcSectorExpectedProfitRate( const vector<IInvestable*>& aInvestables,
                                                              const NationalAccount& aNationalAccount,
                                                              const string& aRegionName,
                                                              const string& aGoodName,
                                                              const double aInvestmentLogitExp,
                                                              const bool aIsShareCalc,
                                                              const bool aIsDistributing,
                                                              const int aPeriod ) const
{
    // Sum levelized cost for the subsector/sector with a logit distribution.
    // TODO: shouldn't there be a way to calculate shares other that rate logit?
    double levelizedCostDenom = 0;
    vector<double> levelizedCosts( aInvestables.size(), 0 );
    int levelizedCostPos = 0;
    
    // beta is equivalent to the share weight of the investable
    for( InvestmentUtils::CInvestableIterator currInv = aInvestables.begin();
        currInv != aInvestables.end(); ++currInv )
    {
        double currLevelizedCost = (*currInv)->getExpectedProfitRate( aNationalAccount,
                                                                      aRegionName,
                                                                      aGoodName,
                                                                      this,
                                                                      // TODO: what is the reason for this hack?
                                                                      -1 * aInvestmentLogitExp, // HACK
                                                                      aIsShareCalc,
                                                                      aIsDistributing,
                                                                      aPeriod );
        if( currLevelizedCost > 0 ){
            levelizedCosts[ levelizedCostPos ] = currLevelizedCost;
            levelizedCostDenom += (*currInv)->getShareWeight( aPeriod ) *
                pow( currLevelizedCost, aInvestmentLogitExp );
        }
        ++levelizedCostPos;
    }

    // If this is the share calc return only the sum of the profit rate to the
    // logit.
    if( aIsShareCalc ){
        return levelizedCostDenom;
    }

    levelizedCostPos = 0;
    double sectorLevelizedCost = 0;
    for( InvestmentUtils::CInvestableIterator currInv = aInvestables.begin();
        currInv != aInvestables.end(); ++currInv )
    {
        double currLevelizedCost = levelizedCosts[ levelizedCostPos ];
        if( currLevelizedCost > 0 ){
            // maybe put this share calc in a utility
            double share = (*currInv)->getShareWeight( aPeriod ) *
                    pow( currLevelizedCost, aInvestmentLogitExp ) / levelizedCostDenom;
            sectorLevelizedCost += share * currLevelizedCost;
        }
        ++levelizedCostPos;
    }

    return sectorLevelizedCost;
}

/*! \brief Calculate the levelized cost for a technology.
* \details This function calculates the levelized cost for a production
*          technology. Currently this calls into the passed in
*          ProductionDemandFunction and uses the its levelized cost function to
*          calculate a base levelized cost. This levelized cost is then adjusted
*          for the investment tax credit. Adjustments should also be made for
*          technologies with longer lifetimes or which have delays before they
*          come on-line, but this has not been implemented.
* \param aTechProdFuncInfo A structure containing the necessary data items to
*        call the production technology's levelized cost function.
* \param aRegionName Name of the region in which investment is occurring.
* \param aGoodName Name of the sector in which investment is occurring.
* \param aDelayedInvestmentTime The lag before this technology will come on-line.
* \param aLifetime Nameplate lifetime of the technology.
* \param aTimestep Length in years of the time step.
* \param aPeriod Period in which to calculate the levelized cost.
* \return Levelized cost per unit of output for the technology.
* \todo Need to handle lifetime and delayed investment time.
*/
double LevelizedCostCalculator::calcTechnologyExpectedProfitRate( const ProductionFunctionInfo& aTechProdFuncInfo,
                                                                  const NationalAccount& aNationalAccount,
                                                                  const string& aRegionName,
                                                                  const string& aSectorName,
                                                                  const double aDelayedInvestmentTime,
                                                                  const int aLifetime,
                                                                  const int aTimeStep,
                                                                  const int aPeriod ) const
{
    // Compute the expected profit rate.
    double levelizedCost = aTechProdFuncInfo.mNestedInputRoot->getLevelizedCost( aRegionName,
                                                                                  aSectorName,
                                                                                  aPeriod );
    // Decrease the raw levelized cost by the investment tax credit rate.
    // Check if this is right.
    levelizedCost /= ( 1 + aNationalAccount.getAccountValue( NationalAccount::INVESTMENT_TAX_CREDIT ) );
    
    // Need to handle lifetime and delayed investment time.

    /*! \post expected profit is greater than or equal to zero.*/
    assert( levelizedCost >= 0 );
    /*! \post expected profit is a valid number. */
    assert( util::isValidNumber( levelizedCost ) );
    return levelizedCost;
}
