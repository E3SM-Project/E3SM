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
 * \file accelerator.cpp
 * \ingroup Objects
 * \brief Accelerator class source file.
 * \author Josh Lurz
 */

#include "util/base/include/definitions.h"
#include <string>
#include <cassert>
#include <cmath>
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>

#include "investment/include/accelerator.h"
#include "util/base/include/xml_helper.h"
#include "containers/include/scenario.h"
#include "util/base/include/model_time.h"
#include "investment/include/investment_utils.h"
#include "util/base/include/util.h"
#include "util/logger/include/ilogger.h"
#include "investment/include/simple_expected_profit_calculator.h"
#include "investment/include/levelized_cost_calculator.h"
#include "investment/include/rate_logit_distributor.h"

// Replace these with a factory method.
#include "investment/include/igrowth_calculator.h"
#include "investment/include/investment_growth_calculator.h"
#include "investment/include/output_growth_calculator.h"

extern Scenario* scenario;

using namespace std;

Accelerator::Accelerator ():
mInvestmentLogitExp( 1 ),
mProfitElasExp( 1 ),
mInvestments( scenario->getModeltime()->getmaxper() ),
mFixedInvestments( scenario->getModeltime()->getmaxper(), -1.0 )
{
}

/*! \brief Destructor
* \note Needed so that that destructor is not in-lined in the header before the
*       IGrowthCalculator class is fully defined.
*/
Accelerator::~Accelerator(){
}

/*! \brief Get the XML node name in static form for comparison when parsing XML.
* \details This public function accesses the private constant string, XML_NAME.
*          This way the tag is always consistent for both read-in and output and
*          can be easily changed. The "==" operator that is used when parsing,
*          required this second function to return static.
* \note A function cannot be static and virtual.
* \author Josh Lurz
* \return The constant XML_NAME as a static.
*/
const std::string& Accelerator::getXMLNameStatic() {
    const static string XML_NAME = "Accelerator";
    return XML_NAME;
}

/*! \brief Parses all data for the class
* \author Josh Lurz
* \param aCurr pointer to the current node in the XML input tree
*/
void Accelerator::XMLParse( const xercesc::DOMNode* aCurr ) {
    /*! \pre make sure we were passed a valid node. */
    assert( aCurr );

    // get all child nodes.
    const xercesc::DOMNodeList* nodeList = aCurr->getChildNodes();

    // loop through the child nodes.
    for( unsigned int i = 0; i < nodeList->getLength(); ++i ){
        const xercesc::DOMNode* curr = nodeList->item( i );
        // Skip any text nodes.
        if( curr->getNodeType() == xercesc::DOMNode::TEXT_NODE ){
            continue;
        }
        const string nodeName = XMLHelper<string>::safeTranscode( curr->getNodeName() );
        if( nodeName == "FixedInvestment" ){
            XMLHelper<double>::insertValueIntoVector( curr, mFixedInvestments, scenario->getModeltime() );
        }
        else if( nodeName == "InvestmentLogitExp" ){
            mInvestmentLogitExp = XMLHelper<double>::getValue( curr );
        }
        else if( nodeName == "ProfitElasExp" ){
            mProfitElasExp = XMLHelper<double>::getValue( curr );
        }
        // THIS NEEDS TO BE FIXED
        else if( nodeName == InvestmentGrowthCalculator::getXMLNameStatic() ){
            // Check if the current type of growth calculator is equal to the new
            // type. If it isn't, which will also happen if the current type is
            // blank, create a new one.
            if( mGrowthCalculatorType != InvestmentGrowthCalculator::getXMLNameStatic() ){
                mGrowthCalculator.reset( new InvestmentGrowthCalculator );
                mGrowthCalculatorType = InvestmentGrowthCalculator::getXMLNameStatic();
            }
            assert( mGrowthCalculator.get() );
            mGrowthCalculator->XMLParse( curr );
        }
        else if( nodeName == OutputGrowthCalculator::getXMLNameStatic() ){
            // Check if the current type of growth calculator is equal to the new
            // type. If it isn't, which will also happen if the current type is
            // blank, create a new one.
            if( mGrowthCalculatorType != OutputGrowthCalculator::getXMLNameStatic() ){
                mGrowthCalculator.reset( new OutputGrowthCalculator );
                mGrowthCalculatorType = OutputGrowthCalculator::getXMLNameStatic();
            }
            assert( mGrowthCalculator.get() );
            mGrowthCalculator->XMLParse( curr );
        }

        else if( nodeName == SimpleExpectedProfitCalculator::getXMLNameStatic() ){
            // Check if the current type of growth calculator is equal to the new
            // type. If it isn't, which will also happen if the current type is
            // blank, create a new one.
            if( mProfitRateCalculatorType != SimpleExpectedProfitCalculator::getXMLNameStatic() ){
                mProfitRateCalculator.reset( new SimpleExpectedProfitCalculator );
                mProfitRateCalculatorType = SimpleExpectedProfitCalculator::getXMLNameStatic();
            }
            assert( mProfitRateCalculator.get() );
        }
        else if( nodeName == LevelizedCostCalculator::getXMLNameStatic() ){
            // Check if the current type of growth calculator is equal to the new
            // type. If it isn't, which will also happen if the current type is
            // blank, create a new one.
            if( mProfitRateCalculatorType != LevelizedCostCalculator::getXMLNameStatic() ){
                mProfitRateCalculator.reset( new LevelizedCostCalculator );
                mProfitRateCalculatorType = LevelizedCostCalculator::getXMLNameStatic();
            }
            assert( mProfitRateCalculator.get() );
        }
        // Add other types of growth calculators here for now until we have a
        // factory method.
        else {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::WARNING );
            mainLog << "Unrecognized node " << nodeName << " found while parsing " << getXMLNameStatic() 
                    << "." << endl;
        }
    }
}

/*! \brief Write the object to an XML output stream for debugging.
* \param aPeriod Period to write debugging information for.
* \param aOut The output stream to write to.
* \param aTabs The object which tracks the number of tabs to write.
*/
void Accelerator::toDebugXML( const int aPeriod, ostream& aOut, Tabs* aTabs ) const {
    XMLWriteOpeningTag( getXMLNameStatic(), aOut, aTabs );
    
    XMLWriteElement( mInvestments[ aPeriod ], "investment", aOut, aTabs,
        scenario->getModeltime()->getper_to_yr( aPeriod ) );
    
    XMLWriteElement( mFixedInvestments[ aPeriod ], "FixedInvestment", aOut, aTabs,
        scenario->getModeltime()->getper_to_yr( aPeriod ) );
    
    XMLWriteElement( mInvestmentLogitExp, "InvestmentLogitExp", aOut, aTabs );
    XMLWriteElement( mProfitElasExp, "ProfitElasExp", aOut, aTabs );
    assert( mGrowthCalculator.get() );
    mGrowthCalculator->toDebugXML( aPeriod, aOut, aTabs );
    assert( mProfitRateCalculator.get() );
    mProfitRateCalculator->toDebugXML( aPeriod, aOut, aTabs );
    XMLWriteClosingTag( getXMLNameStatic(), aOut, aTabs );
}

/*! \brief Write the object to an XML output stream.
* \param aOut The output stream to write to.
* \param aTabs The object which tracks the number of tabs to write.
*/
void Accelerator::toInputXML( ostream& aOut, Tabs* aTabs ) const {
    XMLWriteOpeningTag( getXMLNameStatic(), aOut, aTabs );
    XMLWriteVector( mFixedInvestments, "FixedInvestment", aOut, aTabs, scenario->getModeltime(), -1.0 );

    XMLWriteElementCheckDefault( mInvestmentLogitExp, "InvestmentLogitExp", aOut, aTabs, 1.0 );
    XMLWriteElementCheckDefault( mProfitElasExp, "ProfitElasExp", aOut, aTabs, 1.0 );
    assert( mGrowthCalculator.get() );
    mGrowthCalculator->toInputXML( aOut, aTabs );
    assert( mProfitRateCalculator.get() );
    mProfitRateCalculator->toInputXML( aOut, aTabs );
    XMLWriteClosingTag( getXMLNameStatic(), aOut, aTabs );
}

/*! \brief Complete the initialization of the Accelerator before it is used.
* \details This function stores the region and sector name internally, and then
*          creates default growth and profit calculation objects if they were
*          not read in. The default growth calculator is the
*          InvestmentGrowthCalculator and the default profit rate calculator is
*          the SimpleExpectedProfitRateCalculator. These objects are created
*          with their default parameters.
* \param aRegionName Name of the region containing the accelerator.
* \param aSectorName Name of the sector containing the accelerator.
*/
void Accelerator::completeInit( const string& aRegionName, const string& aSectorName ){
    mRegionName = aRegionName;
    mSectorName = aSectorName;
    // Create a default type of growth calculator if one was not created during
    // parsing.
    if( !mGrowthCalculator.get() ){
        mGrowthCalculator.reset( new InvestmentGrowthCalculator );
        mGrowthCalculatorType = InvestmentGrowthCalculator::getXMLNameStatic();
    }
    if( !mProfitRateCalculator.get() ){
        mProfitRateCalculator.reset( new SimpleExpectedProfitCalculator() );
        mProfitRateCalculatorType = SimpleExpectedProfitCalculator::getXMLNameStatic();
    }

    // Warn if fixed investment for the base period was read in, as it will be ignored.
    if( mFixedInvestments[ 0 ] != -1 ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::NOTICE );
        mainLog << "Will ignore fixed investment in the base period for sector "
                << mSectorName << " in region " << mRegionName << endl;
    }
}

/*! \brief Initialization of the accelerator based investor for each period.
* \param aRegionName Name of the region containing this investor.
* \param aSectorName Name of the sector the investor is investing in.
* \author Josh Lurz
*/
void Accelerator::initCalc( vector<IInvestable*>& aInvestables,
                            NationalAccount& aNationalAccount, 
                            const Demographic* aDemographic,
                            const int aPeriod )
{
    /*! \pre Check that the period is not nonsensical */
    assert( aPeriod >= 0 );
    if( aPeriod == 0 ){
        // Don't need to determine investment for the base period.
        mInvestments[ aPeriod ] = InvestmentUtils::sumInvestment( aInvestables, aPeriod );
        assert( mInvestments[ aPeriod ] >= 0 );
    }        
}
/*! \brief Calculates a lump sum investment using an acceleration or scaling and
*          distributes that amount.
* \details This function is the main operation of the Accelerator. In the base
*          period, this function returns the sum of all read in investment for
*          the subsectors. In later periods, it uses the read-in
*          IGrowthCalculator object to accelerate investment, and then
*          determines a lump sum of new investment. It then uses the read-in
*          IExpectedProfitRateCalculator object and a IDistributor object to
*          distribute the investment among the subsectors. Finally, it checks
*          that the investment distributed is equal to the initial investment
*          calculated, stores the value internally, and returns it.
* \param aInvestables The vector of children which will receive investment.
* \param aNationalAccount The national accounts container.
* \param aDemographic A pointer to the Demographics object.
* \param aPeriod The period in which investment is calculated and distributed.
* \return The total investment which occurred.
* \author Josh Lurz
*/

double Accelerator::calcAndDistributeInvestment( vector<IInvestable*>& aInvestables,
                                                 NationalAccount& aNationalAccount, 
                                                 const Demographic* aDemographic,
                                                 const int aPeriod )
{   
    // For the base period use read-in investments. Do not calculate new investments.
    // Revise to calculate investment is base period is to be solved.
    // SHK debugging
    double currInvestment = 0;
    if( aPeriod == 0 ){
        //return mInvestments[ aPeriod ];
        currInvestment = mInvestments[ aPeriod ];
    }
    else{
        currInvestment = mInvestments[ aPeriod - 1 ];
    }

    // Calculate the investment dependency scalar.
    double invDepScalar = mGrowthCalculator->calcInvestmentDependencyScalar( aInvestables,
                          aDemographic, aNationalAccount, mSectorName, mRegionName, 
//                          mInvestments[ aPeriod - 1 ], mInvestmentLogitExp, aPeriod );
                          currInvestment, mInvestmentLogitExp, aPeriod );
    // Calculate total investment
    const double newInvestment = calcNewInvestment( aInvestables, aNationalAccount,
                                                    invDepScalar, aPeriod );
    assert( newInvestment >= 0 );

    // Create a logit based investment distributor.
    RateLogitDistributor invDistributor;

    // Use the investment distributor to distribute the investment.
    mInvestments[ aPeriod ] = invDistributor.distribute( mProfitRateCalculator.get(),
                              aInvestables, aNationalAccount, mInvestmentLogitExp,
                              mRegionName, mSectorName, newInvestment, aPeriod );

    // Check that total investment and distributed investment are equal.
    if( !util::isEqual( newInvestment, mInvestments[ aPeriod ] ) ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::WARNING );
        mainLog << fabs( newInvestment - mInvestments[ aPeriod ] )
                << " difference between desired investment and distributed investment in "
                << mSectorName << " in " << mRegionName << endl;
    }
    // Return the total amount of investment actually distributed
    return mInvestments[ aPeriod ];
}

/*! \brief A function to calculate the amount of new investment given a growth
*          scalar.
* \details This function initially checks for any level of fixed investment,
*          which it uses as the investment quantity if it was read-in. It then
*          calculated a sector level expected profit rate to use to scale the
*          investment. This does not use the read-in IExpectedProfitRate object,
*          it always uses the SimpleExpectedProfitRate calculator. This is
*          because an expected profit rate is always desired to scale
*          investment. It then uses this expected profit rate and the passed in
*          scalar to calculate a new total investment, based on the previous
*          periods. It finally checks to make sure this value exceeds the fixed
*          investment which must occur in this period, and returns the higher of
*          the two values.
* \param aInvestables The vector of children which will receive investment.
* \param aNationalAccount The national accounts container.
* \param aCapDependencyScalar The previous calculator growth scalar.
* \param aPeriod The period in which to calculate new investment.
* \return The total investment which should occur.
* \author Josh Lurz
*/
double Accelerator::calcNewInvestment( vector<IInvestable*>& aInvestables,
                                       NationalAccount& aNationalAccount, 
                                       const double aCapDependencyScalar,
                                       const int aPeriod ) const 
{
    // Check if parent level investment is fixed.
    if( mFixedInvestments[ aPeriod ] != -1 ){
        // Use sector level fixed investment instead of calculating a scaled
        // level. If all children are fixed and that sum does not equal this
        // value, there will be a warning later when the investment is
        // distributed.
        double newInvestment = mFixedInvestments[ aPeriod ];
        // Check if total child investment is greater than this amount, as that
        // will always override this amount.
        double childSumFixed = InvestmentUtils::sumFixedInvestment( aInvestables, aPeriod );
        if( childSumFixed > mFixedInvestments[ aPeriod ] ){
            newInvestment = childSumFixed;
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::NOTICE );
            mainLog << "Overriding parent level investment with child level investment sum. " << endl;
        }
        return newInvestment;
    }

    // Create a standard profit rate calculator.
    SimpleExpectedProfitCalculator expProfitRateCalc;

    // Calculate the sector level expected profit rate.
    const double expProfitRate = expProfitRateCalc.calcSectorExpectedProfitRate( aInvestables,
                                                                                 aNationalAccount,
                                                                                 mRegionName,
                                                                                 mSectorName,
                                                                                 mInvestmentLogitExp,
                                                                                 false,
                                                                                 false,
                                                                                 aPeriod );
    
    // Calculate the total new investment using the expected profit and the
    // profit rate elasticity.
    const double newInvestment = ( expProfitRate > 0 ) ?
                                 aCapDependencyScalar * pow( expProfitRate, mProfitElasExp ) : 0;
    // Determine if all subsectors have fixed investment and the total amount of
    // fixed investment.
    double sumFixed = InvestmentUtils::sumFixedInvestment( aInvestables, aPeriod );
    // Ensure that we always return at least the amount of fixed investment
    // regardless of the profit rate.
    return max( newInvestment, sumFixed );
}   

/*! \brief A function to sets the efficienty conditions for the investment.
* \details Do nothing the investment based on the Accelerator.
* \param aInvestables The vector of children which will receive investment.
* \param aNationalAccount The national accounts container.
* \param aCapDependencyScalar The previous calculator growth scalar.
* \param aPeriod The period in which to calculate new investment.
* \author Sonny Kim
*/
void Accelerator::setEfficiencyConditions( vector<IInvestable*>& aInvestables,
                                           NationalAccount& aNationalAccount, 
                                           const Demographic* aDemographic,
                                           const int aPeriod ) const
{
    // Do nothing.  Only applicable for market based investment calculation.
}

