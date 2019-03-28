#ifndef _IEXPECTED_PROFIT_CALCULATOR_H_
#define _IEXPECTED_PROFIT_CALCULATOR_H_
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
 * \file iexpected_profit_calculator.h
 * \ingroup Objects
 * \brief The IExpectedProfitCalculator interface header file.
 * \author Josh Lurz
 */

#include <vector>
#include <iosfwd>
class NationalAccount;
class IInvestable;
class IInput;
class IFunction;
class Tabs;
struct ProductionFunctionInfo;

/*! 
 * \ingroup Objects
 * \brief This is the interface to an object responsible for calculating the
 *        expected profit rate of a ProductionSector, Subsector or
 *        ProductionTechnology.
 * \details TODO
 * \author Josh Lurz
 */
class IExpectedProfitRateCalculator
{
public:
    IExpectedProfitRateCalculator();

	virtual ~IExpectedProfitRateCalculator();

    // TODO: Inherit so that documentation is inherited.
    // TODO: Add an XMLParse since there is a toInputXML.
    virtual void toInputXML( std::ostream& aOut, Tabs* aTabs ) const = 0;

    virtual void toDebugXML( const int aPeriod,
                             std::ostream& aOut,
                             Tabs* aTabs ) const = 0;

    /*!
     * \brief Calculate the expected profit rate for a ProductionSector or
     *        Subsector.
     * \details Determines the expected profit rate for a Sector or Subsector
     *          using the given parameters and the passed in investable
     *          children.
     * \param aInvestables Children of the current sector or subsector for which
     *        to calculate the average expected profit rate.
     * \param aNationalAccount Regional national accounts container.
     * \param aRegionName Name of the region in which expected profit rates are
     *        being calculated.
     * \param aSectorName Name of the sector for which expected profit rates are
     *        being calculated.
     * \param aInvestmentLogitExp The investment logit exponential.
     * \param aIsShareCalc Whether this calculation is for the investment share
     *        calculation.
     * \param aIsDistributing Whether this expected profit rate is being used
     *        to distribute investment.
     * \param aPeriod Period in which to calculate expected profits.
     * \return The expected profit rate for the ProductionSector or Subsector.
     */
    virtual double calcSectorExpectedProfitRate( const std::vector<IInvestable*>& aInvestables,
                                                 const NationalAccount& aNationalAccount,
                                                 const std::string& aRegionName,
                                                 const std::string& aSectorName,
                                                 const double aInvestmentLogitExp,
                                                 const bool aIsShareCalc,
                                                 const bool aIsDistributing,
                                                 const int aPeriod ) const = 0;

    /*!
     * \brief Calculate the expected profit rate for a ProductionTechnology
     * \details Determines the expected profit rate for a ProductionTechnology
     *          using the passed in parameters.
     * \param aTechProdFuncInfo A structure containing the necessary information
     *        to call the production function for the given
     *        ProductionTechnology. This includes the production function
     *        itself.
     * \param aNationalAccount Regional national accounts container.
     * \param aRegionName Name of the region in which expected profit rates are
     *        being calculated.
     * \param aSectorName Name of the sector for which expected profit rates are
     *        being calculated.
     * \param aDelayedInvestmentTime Amount of time between the investment
     *        occurring and the ProductionTechnology being brought on-line.
     * \param aLifetime Lifetime of the ProductionTechnology.
     * \param aTimeStep Time step for the period in which investment is
     *        occurring.
     * \param aPeriod Period in which to calculate expected profits.
     * \return The expected profit rate for the ProductionTechnology.
     */
    virtual double calcTechnologyExpectedProfitRate( const ProductionFunctionInfo& aTechProdFuncInfo,
                                                     const NationalAccount& aNationalAccount,
                                                     const std::string& aRegionName,
                                                     const std::string& aSectorName,
                                                     const double aDelayedInvestmentTime,
                                                     const int aLifetime,
                                                     const int aTimeStep,
                                                     const int aPeriod ) const = 0;
};

// Define empty inline methods.
//! Constructor
inline IExpectedProfitRateCalculator::IExpectedProfitRateCalculator(){
}

//! Destructor
inline IExpectedProfitRateCalculator::~IExpectedProfitRateCalculator(){
}

#endif // _IEXPECTED_PROFIT_CALCULATOR_H_
