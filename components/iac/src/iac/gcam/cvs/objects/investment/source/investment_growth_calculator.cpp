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
 * \file investment_growth_calculator.cpp
 * \ingroup Objects
 * \brief InvestmentGrowthCalculator class source file.
 * \author Josh Lurz
 */

#include "util/base/include/definitions.h"
#include <string>
#include <cassert>
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#include "investment/include/investment_growth_calculator.h"
#include "util/base/include/xml_helper.h"
#include "demographics/include/demographic.h"
#include "investment/include/investment_utils.h"

using namespace std;

//! Constructor which initializes all values to their defaults.
InvestmentGrowthCalculator::InvestmentGrowthCalculator ():
mAggregateInvestmentFraction( 0.01 ),
mInvestmentAcceleratorScalar( 1.2 ),
mEconomicGrowthExp( 1 ),
mMarginalValueDollar( 1 )
{
}

//! Return the XML name of this object statically.
const string& InvestmentGrowthCalculator::getXMLNameStatic(){
    const static string XML_NAME = "investment-growth-calculator";
    return XML_NAME;
}

/*! \brief Parses all data associated with the class.
* \author Josh Lurz
* \param aNode pointer to the current node in the XML input tree
*/
void InvestmentGrowthCalculator::XMLParse( const xercesc::DOMNode* aNode ) {
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
        if( nodeName == "AggregateInvestmentFraction" ){
            mAggregateInvestmentFraction = XMLHelper<double>::getValue( curr );
        }
        else if( nodeName == "InvestmentAcceleratorScalar" ){
            mInvestmentAcceleratorScalar = XMLHelper<double>::getValue( curr );
        }
        else if( nodeName == "EconomicGrowthExp" ){
            mEconomicGrowthExp = XMLHelper<double>::getValue( curr );
        }
        else if( nodeName == "MarginalValueDollar" ){
            mMarginalValueDollar = XMLHelper<double>::getValue( curr );
        }
        else {
            cout << "Warning unknown node " << nodeName << " found while parsing " << getXMLNameStatic() << endl;
        }
    }
}

//! Write out debugging information.
void InvestmentGrowthCalculator::toDebugXML( const int aPeriod, ostream& aOut, Tabs* aTabs ) const {
    XMLWriteOpeningTag( getXMLNameStatic(), aOut, aTabs );
    XMLWriteElement( mAggregateInvestmentFraction, "AggregateInvestmentFraction", aOut, aTabs );
    XMLWriteElement( mInvestmentAcceleratorScalar, "InvestmentAcceleratorScalar", aOut, aTabs );
    XMLWriteElement( mEconomicGrowthExp, "EconomicGrowthExp", aOut, aTabs );
    XMLWriteElement( mMarginalValueDollar, "MarginalValueDollar", aOut, aTabs );
    XMLWriteClosingTag( getXMLNameStatic(), aOut, aTabs );
}

//! Write out input XML information.
void InvestmentGrowthCalculator::toInputXML( ostream& aOut, Tabs* aTabs ) const {
    XMLWriteOpeningTag( getXMLNameStatic(), aOut, aTabs );
    XMLWriteElementCheckDefault( mAggregateInvestmentFraction, "AggregateInvestmentFraction", aOut, aTabs );
    XMLWriteElementCheckDefault( mInvestmentAcceleratorScalar, "InvestmentAcceleratorScalar", aOut, aTabs );
    XMLWriteElementCheckDefault( mEconomicGrowthExp, "EconomicGrowthExp", aOut, aTabs, 1.0 );
    XMLWriteElementCheckDefault( mMarginalValueDollar, "MarginalValueDollar", aOut, aTabs, 1.0 );
    XMLWriteClosingTag( getXMLNameStatic(), aOut, aTabs );
}

/*! \brief Calculate an overall scalar used to grow investment from the previous
*          aPeriod.
* \details Calculates a scalar to grow investment which is based on several
*          read-in parameters and the regional demographics. The parameter is
*          calculated as the previous period's capital multiplied by the read-in
*          growth rate, the economic growth rate, and an adjustment for the
*          marginal value of the dollar.
* \param aInvestables The investable children contained by the object currently
*        calculating investment.
* \param aDemographic The Demographics object needed for calculating the
*        increase in working age population.
* \param aNationalAccount The national account container.
* \param aGoodName The name of the sector calculating investment.
* \param aRegionName The name of the region.
* \param aPrevInvestment Total sector investment for the previous period.
* \param aInvestmentLogitExp The investment logit exponential.
* \param aPeriod The aPeriod in which to calculate the scalar.
* \return An overall scalar used to grow investment from the previous aPeriod.
* \author Josh Lurz
*/
double InvestmentGrowthCalculator::calcInvestmentDependencyScalar( const vector<IInvestable*>& aInvestables,
                                                                   const Demographic* aDemographic,
                                                                   const NationalAccount& aNationalAccount,
                                                                   const string& aGoodName,
                                                                   const string& aRegionName,
                                                                   const double aPrevInvestment,
                                                                   const double aInvestmentLogitExp,
                                                                   const int aPeriod ) 
{   
    // Calculate the starting level of investment based on the previous period.
    const double baseCapital = InvestmentUtils::calcBaseCapital( aRegionName, aPrevInvestment,
                                                                 mAggregateInvestmentFraction, aPeriod );

    // Calculate the scalar based on economic growth. 
    const double economicGrowthScalar = calcEconomicGrowthScalar( aDemographic, aPeriod );
    assert( economicGrowthScalar > 0 );

    // Calculate the scalar for the investment.
    double invDepScalar = baseCapital * mInvestmentAcceleratorScalar * economicGrowthScalar *
                          pow( mMarginalValueDollar, -1 * aInvestmentLogitExp );
    assert( invDepScalar > 0 );
    return invDepScalar;
}

/*! \brief Calculate a growth scalar based on the increase in economic activity
*          in the region.
* \brief Calculates a parameter to approximate economic growth in the region.
*        This is equal to the growth rate in working age population for the
*        region raised to a read-in exponent.
* \param aDemographic The Demographic object used to calculate the change in
*        working age population.
* \param The aPeriod in which to calculate the economic growth scalar.
* \return The economic growth scalar.
* \author Josh Lurz
*/
double InvestmentGrowthCalculator::calcEconomicGrowthScalar( const Demographic* aDemographic,
                                                             const int aPeriod ) const
{
    /*! \pre aPeriod is greater than the base aPeriod. */
    assert( aPeriod > 0 );
    
    // Calculate the change in the working age population.
    double workingAgeRateChange = aDemographic->getWorkingAgePopulation( aPeriod ) 
                                  / aDemographic->getWorkingAgePopulation( aPeriod - 1 );
    assert( workingAgeRateChange > 0 );
    
    // Calculate the economic growth scalar. 
    double econGrowthScalar = pow( workingAgeRateChange, mEconomicGrowthExp );
    
    /*! \post Economic Growth Scalar is greater than zero. */
    assert( econGrowthScalar > 0 );
    
    return econGrowthScalar;
}
