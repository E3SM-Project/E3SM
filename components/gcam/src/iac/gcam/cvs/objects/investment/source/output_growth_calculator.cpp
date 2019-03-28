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
 * \file output_growth_calculator.cpp
 * \ingroup Objects
 * \brief OutputGrowthCalculator class source file.
 * \author Josh Lurz
 */

#include "util/base/include/definitions.h"
#include <string>
#include <cassert>
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#include "investment/include/output_growth_calculator.h"
#include "util/base/include/xml_helper.h"
#include "investment/include/investment_utils.h"
#include "containers/include/scenario.h"
#include "util/base/include/model_time.h"
#include "marketplace/include/marketplace.h"
#include "investment/include/iinvestable.h"
#include "investment/include/simple_expected_profit_calculator.h"
#include "investment/include/rate_logit_distributor.h"
#include "util/logger/include/ilogger.h"

using namespace std;

extern Scenario* scenario;

//! Constructor which initializes all member variables to default values.
OutputGrowthCalculator::OutputGrowthCalculator ():
mAggregateInvestmentFraction( 0.01 ),
mOutputGrowthRate( scenario->getModeltime()->getmaxper() ),
mTrialCapital( scenario->getModeltime()->getmaxper() )
{
}

//! Return the XML name of this object statically.
const string& OutputGrowthCalculator::getXMLNameStatic(){
    const static string XML_NAME = "output-growth-calculator";
    return XML_NAME;
}

/*! \brief Parses all data associated with the class.
* \author Josh Lurz
* \param aNode pointer to the current node in the XML input tree
*/
void OutputGrowthCalculator::XMLParse( const xercesc::DOMNode* aNode ) {
    /*! \pre make sure we were passed a valid node. */
    assert( aNode );

    // get all child nodes.
    const xercesc::DOMNodeList* nodeList = aNode->getChildNodes();

    // loop through the child nodes.
    for( unsigned int i = 0; i < nodeList->getLength(); i++ ){
        const xercesc::DOMNode* curr = nodeList->item( i );
        if( curr->getNodeType() == xercesc::DOMNode::TEXT_NODE ){
            continue;
        }

        const string nodeName = XMLHelper<string>::safeTranscode( curr->getNodeName() );
        if( nodeName == "aggregate-investment-fraction" ){
            mAggregateInvestmentFraction = XMLHelper<double>::getValue( curr );
        }
        else if( nodeName == "output-growth-rate" ){
            XMLHelper<double>::insertValueIntoVector( curr, mOutputGrowthRate, scenario->getModeltime() );
        }
        else {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::DEBUG );
            mainLog << "Unknown node " << nodeName << " found while parsing " << getXMLNameStatic() << endl;
        }
    }
}

//! Write out debugging information.
void OutputGrowthCalculator::toDebugXML( const int aPeriod, ostream& aOut, Tabs* aTabs ) const {
    XMLWriteOpeningTag( getXMLNameStatic(), aOut, aTabs );
    XMLWriteElement( mAggregateInvestmentFraction, "aggregate-investment-fraction", aOut, aTabs );
    XMLWriteElement( mOutputGrowthRate[ aPeriod ], "output-growth-rate", aOut, aTabs );
    XMLWriteElement( mTrialCapital[ aPeriod ], "trial-capital", aOut, aTabs );
    XMLWriteClosingTag( getXMLNameStatic(), aOut, aTabs );
}

//! Write out input XML information.
void OutputGrowthCalculator::toInputXML( ostream& aOut, Tabs* aTabs ) const {
    XMLWriteOpeningTag( getXMLNameStatic(), aOut, aTabs );
    XMLWriteElementCheckDefault( mAggregateInvestmentFraction, "aggregate-investment-fraction", aOut, aTabs );
    
    const Modeltime* modeltime = scenario->getModeltime();
    // 0 isn't really the default.
    XMLWriteVector( mOutputGrowthRate, "output-growth-rate", aOut, aTabs, modeltime, 0.0 );
    XMLWriteClosingTag( getXMLNameStatic(), aOut, aTabs );
}

/*! \brief Calculate an overall scalar used to grow investment from the previous
*          period.
* \details Calculates a scalar which accelerates investment from the previous
*          period by growing output at a set rate. To calculate the investment
*          growth the function first calculates a desired level of output and
*          then backs out the trial investment using an average sector capital
*          to output ratio.
* \note This calculation is only performed on the first iteration of a period,
*       so the final equilibrium output will not be exactly the last period's
*       output multiplied by the growth rate as the capital to output ratio will
*       have adjusted during solution.
* \param aInvestables The investable children of this investment object.
* \param aDemographic The Demographics object.
* \param aNationalAccount The national accounts container.
* \param aGoodName The name of the sector in which investment is occurring.
* \param aRegionName The name of the region containing the sector being invested
*        in.
* \param aPrevInvestment Total sector investment for the previous period.
* \param aInvestmentLogitExp The investment logit exponential.
* \param aPeriod The period in which to calculate the scalar.
* \return An overall scalar used to grow investment from the previous period.
* \author Josh Lurz
*/
double OutputGrowthCalculator::calcInvestmentDependencyScalar( const vector<IInvestable*>& aInvestables,
                                                               const Demographic* aDemographic,
                                                               const NationalAccount& aNationalAccount,
                                                               const string& aGoodName,
                                                               const string& aRegionName,
                                                               const double aPrevInvestment,
                                                               const double aInvestmentLogitExp,
                                                               const int aPeriod ) 
{
    // TODO: We should store the entire scalar for the period, not just TESTK
    // like Fortran does.
    
    // Calculate the starting level of capital.
    const double baseCapital = InvestmentUtils::calcBaseCapital( aRegionName, aPrevInvestment,
                                                                 mAggregateInvestmentFraction, 
                                                                 aPeriod );
    // Calculate trial capital from the subsectors.
    const double trialCapital = calcTrialCapital( aInvestables, aNationalAccount, aGoodName,
                                                  aRegionName, aInvestmentLogitExp, aPeriod );
    
    // Calculate the trial investment. Performs a capital interpolation.
    const double trialInvestment = ( trialCapital - ( 2 * baseCapital ) ) / 3;

    // Calculate the scalar.
    const double invDepScalar = max( trialInvestment, trialCapital / 10 );

    /*! \post The investment scalar is positive */
    assert( invDepScalar > 0 );
    return invDepScalar;
}

/*! \brief Calculate the trial amount of capital required to meet the desired
*          output level.
* \details This function will calculate the amount of capital required to meet
*          the desired level of output, but only in the first iteration of each
*          period. In latter iterations, it will return the value stored in the
*          first iteration. To calculate this capital quantity, the function
*          first calculates the sector average capital to output ratio and the
*          output gap. The function uses these to calculate the amount of
*          capital required to meet the output gap, and adds the amounts
*          required by fixed investment and interpolated annual investment flows
*          to determine the level of trial capital.
* \note The Fortran equivalent of this return value of this function is TESTK.
* \param aInvestables The investable objects
* \param aNationalAccount The national accounts container.
* \param aGoodName The sector name.
* \param aRegionName The region name.
* \param aInvestmentLogitExp The investment logit function exponential.
* \param aPeriod The period in which to calculate trial capital.
* \return The amount of trial capital.
*/
double OutputGrowthCalculator::calcTrialCapital( const vector<IInvestable*>& aInvestables, 
                                                 const NationalAccount& aNationalAccount,
                                                 const string& aGoodName,
                                                 const string& aRegionName,
                                                 const double aInvestmentLogitExp,
                                                 const int aPeriod )
{
    // Fortran code special cases period 0, but you shouldn't calculate investment in period 0.
    assert( aPeriod > 0 );
    
    // If the trial capital for this period has already been calculated, return
    // the stored value. This ensures this calculation is only performed once
    // per period.
    if( mTrialCapital[ aPeriod ] > 0 ){
        return mTrialCapital[ aPeriod ];
    }
    
    // Calculate the additional output required to meet the new level of output.
    const double additionalOutput = calcOutputGap( aInvestables, aNationalAccount, aGoodName,
                                                   aRegionName, aInvestmentLogitExp, aPeriod );
    
    // Create a simple expected profit rate calculator and distributor which are
    // needed to calculate the average sector capital to output ratio.
    SimpleExpectedProfitCalculator expProfitRateCalc;
    RateLogitDistributor distributor;

    // Calculate the amount of capital required to produce one unit of output.
    const double capitalOutputRatio = distributor.calcCapitalOutputRatio( aInvestables,
                                                                          &expProfitRateCalc,
                                                                          aNationalAccount,
                                                                          aRegionName,
                                                                          aGoodName,
                                                                          aPeriod );
    
    // Calculate the trial level of capital needed to reach the desired output.
    mTrialCapital[ aPeriod ] = capitalOutputRatio * additionalOutput +
                               InvestmentUtils::sumInvestment( aInvestables, aPeriod - 1 ) * 2 +
                               InvestmentUtils::sumFixedInvestment( aInvestables, aPeriod ) * 3;
    
    assert( mTrialCapital[ aPeriod ] > 0 );
    return mTrialCapital[ aPeriod ];
}

/*! \brief Calculate the output gap between what would be produced with no new
*          investment and what is desired based on an output growth rate and the
*          previous period's output.
* \details To calculate the output gap, the function first scales previous
*          period's supply by the read-in growth rate for the period. This is
*          the projected sales. The function then deducts all output produced by
*          existing vintages to determine output that must be produced by new
*          vintages. Next the function iterates through the children of this
*          investor and calculates a child level capital to output ratio. This
*          is used to calculate output expected to be produced by fixed
*          investment in the current period, and produced by last period's
*          annual investment scaling down to zero. These outputs are removed
*          from the desired new vintage output to determine the output gap,
*          which the function ensures is zero or greater.
* \param aInvestables The investable objects
* \param aNationalAccount The national accounts container.
* \param aGoodName The sector name.
* \param aRegionName The region name.
* \param aInvestmentLogitExp The investment logit function exponential.
* \param aPeriod The period in which to calculate the output gap.
* \return The output gap.
*/
double OutputGrowthCalculator::calcOutputGap( const vector<IInvestable*>& aInvestables,
                                              const NationalAccount& aNationalAccount,
                                              const string& aGoodName,
                                              const string& aRegionName,
                                              const double aInvestmentLogitExp,
                                              const int aPeriod ) const 
{
    const Marketplace* marketplace = scenario->getMarketplace();
    
    // Get the sales of the good and accelerate the growth.
    const double projSales = marketplace->getSupply( aGoodName, aRegionName, aPeriod - 1 ) 
                             * mOutputGrowthRate[ aPeriod ];
    
    // Get the supply of the good. This is called after only operating old vintages(QOLD)
    const double oldVintageSupply = marketplace->getSupply( aGoodName, aRegionName, aPeriod );
    
    // Calculate the quantity the new vintages will need to produce(QNEW)
    const double desiredNewVintageOutput = projSales - oldVintageSupply;

    // Calculate the amount of investment that would be done if the accelerator was one.
    // This uses the capital interpolation formula to determine a total.
    double totalFixedInvestment = 0;
    double totalPrevAnnualInvestment = 0;

    // Create a simple expected profit rate calculator and distributor which are
    // needed to calculate the average sector capital to output ratio.
    SimpleExpectedProfitCalculator expProfitRateCalc;
    RateLogitDistributor distributor;
    
    for( unsigned int i = 0; i < aInvestables.size(); ++i ){
        // Determine the capital to output ratio for the subsector.
        const double capOutputRatio = aInvestables[ i ]->getCapitalOutputRatio( &distributor,
                                                                                &expProfitRateCalc,
                                                                                aNationalAccount,
                                                                                aRegionName,
                                                                                aGoodName,
                                                                                aPeriod );
        // Invert the ratio to find output per unit of capital.
        const double capQuotient = ( capOutputRatio > 0 ) ? 1 / capOutputRatio : 0;

        // Add the previous periods investment.
        // I think these are only right for a 5 year time step.
        totalPrevAnnualInvestment += aInvestables[ i ]->getAnnualInvestment( aPeriod - 1 ) *
                                     capQuotient * 2;
        // Add any fixed investment.
        totalFixedInvestment += aInvestables[ i ]->getFixedInvestment( aPeriod ) *
                                capQuotient * 3;
    }

    // Determine the gap between desired and actual output(RGAP). Don't let the
    // output gap drop below zero.
    return  max( desiredNewVintageOutput - totalFixedInvestment - totalPrevAnnualInvestment, 0.0 );
}
