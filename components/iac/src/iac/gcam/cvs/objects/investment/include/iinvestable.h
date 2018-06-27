#ifndef _IINVESTABLE_H_
#define _IINVESTABLE_H_
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
 * \file iinvestable.h
 * \ingroup Objects
 * \brief The IInvestable interface header file.
 * \author Josh Lurz
 */
#include "util/base/include/ivisitable.h"

class NationalAccount;
class IExpectedProfitRateCalculator;
class IDistributor;

/*! 
 * \ingroup Objects
 * \brief An interface to any object which may be invested in.
 * \details This interface is used to represent an object that can accept
 *          investment. This interface is currently implemented by SGM
 *          subsectors and ProductionTechnologies. This interface allows
 *          investment decision making objects to treat sets of Subsectors and
 *          ProductionTechnologies equivalently.
 * \todo Fix argument ordering of these functions to match.
 * \author Josh Lurz
 */
class IInvestable: public IVisitable
{
public:
    IInvestable();
    virtual ~IInvestable();
    
    /*!
     * \brief Get the expected profit rate for the object.
     * \details Uses the IExpectedProfitRateCalculator object to calculate the
     *          expected profit rate for the object in the given period.
     * \param aNationalAccount Regional national accounts container.
     * \param aRegionName Region name.
     * \param aSectorName Sector name.
     * \param aExpProfitRateCalc An object responsible for calculating expected
     *        profit rates. This object must be used by the implementing class
     *        to determine the expected profit rate.
     * \param aInvestmentLogitExp The investment logit exponential.
     * \param aIsShareCalculation Whether this is expected profit rate is being
     *        used for a share calculation.
     * \param aIsDistributing Whether this expected profit rate is being used
     *        to distribute investment.
     * \param aPeriod Period.
     * \return The expected profit rate for the current level.
     */
    virtual double getExpectedProfitRate( const NationalAccount& aNationalAccount,
                                          const std::string& aRegionName,
                                          const std::string& aSectorName,
                                          const IExpectedProfitRateCalculator* aExpProfitRateCalc,
                                          const double aInvestmentLogitExp,
                                          const bool aIsShareCalc,
                                          const bool aIsDistributing,
                                          const int aPeriod ) const = 0;

    /*!
     * \brief Get the quantity of fixed investment for the period.
     * \details Objects which implement this interface may choose to read in or
     *          otherwise set a quantity of fixed investment by period. This
     *          fixed investment quantity will be used as the investment amount
     *          for the object in that period, as long as the overall investment
     *          at that level is enough to satisfy all fixed investment.
     * \param Period for which to get the fixed investment.
     * \return Fixed investment for the period.
     */
    virtual double getFixedInvestment( const int aPeriod ) const = 0;

    /*!
     * \brief Get the quantity of fixed investment for the period.
     * \details Objects which implement this interface may choose to read in or
     *          otherwise set a quantity of fixed investment by period. This
     *          fixed investment quantity will be used as the investment amount
     *          for the object in that period, as long as the overall investment
     *          at that level is enough to satisfy all fixed investment.
     * \param Period for which to get the fixed investment.
     * \return Fixed investment for the period.
     */
    virtual bool hasCalibrationMarket() const = 0;

    /*!
     * \brief Get the annual investment for a period.
     * \details Returns the annual investment for the period. This not available
     *          until investment has been calculated and distributed for the
     *          period.
     * \param aPeriod Period for which to return annual investment. If this is
     *        -1 then return any annual investment.
     * \todo Don't allow period to be -1.
     * \return Annual investment for the period.
     */
    virtual double getAnnualInvestment( const int aPeriod ) const = 0;

    /*!
     * \brief Return the amount of capital required to produce one unit of
     *        output.
     * \details Determines the amount of capital required to produce one unit of
     *          output, given that the IDistributor is used to distribute the
     *          investment with the given parameters. If there was only one
     *          subsector and one technology, this would be equal to the capital
     *          coefficient multiplied by the overall scalar(alpha zero). This
     *          function must use the same methods for calculating profits and
     *          distributing investment as will later be used once the
     *          investment level is known for this function to be correct.
     * \param aDistributor The object responsible for distributing investment.
     * \param aExpProfitRateCalculator The object responsible for calculating
     *        expected profits.
     * \param aNationalAccount The regional national accounts container.
     * \param aRegionName Region name.
     * \param aSectorName Sector name.
     * \param aPeriod Model period.
     * \return The amount of capital required to produce one unit of output.
     */
    virtual double getCapitalOutputRatio( const IDistributor* aDistributor,
                                          const IExpectedProfitRateCalculator* aExpProfitRateCalc,
                                          const NationalAccount& aNationalAccount,
                                          const std::string& aRegionName,
                                          const std::string& aSectorName, 
                                          const int aPeriod ) const = 0;

    /*! 
     * \brief Get the total output.
     * \details Returns the total quantity of output for the given period by the
     *          object and all its' children. This must be called after output
     *          is calculated for the period.
     * \return Total output for the period.
     */
    virtual double getOutput( const int aPeriod ) const = 0;
    
    /*!
     * \brief Distribute investment to the object.
     * \details Distributes investment to the object. If this object itself has
     *          more investable children, it must use the IDistributor and
     *          aExpProfitRateCalc to further distribute the investment to its'
     *          children.
     * \param aDistributor Object responsible for distributing investment.
     * \param aNationalAccount Regional national accounts container.
     * \param aExpProfitRateCalc Object responsible for calculating expected
     *        profits.
     * \param aRegionName Region name.
     * \param aSectorName Sector name.
     * \param aNewInvestment Amount of investment to distribute.
     * \param aPeriod Model period.
     * \return The total amount of investment actually distributed. This may
     *         differ from aNewInvestment if there were not enough profitable
     *         children to which to distribute investment or all children had
     *         fixed investment.
     */
    virtual double distributeInvestment( const IDistributor* aDistributor,
                                         NationalAccount& aNationalAccount,
                                         const IExpectedProfitRateCalculator* aExpProfitRateCalc,
                                         const std::string& aRegionName,
                                         const std::string& aSectorName,
                                         const double aNewInvestment,
                                         const int aPeriod ) = 0;

    /*!
     * \brief Gets the share weight from the object.
     * \details Gets the value used to bais the investment to an Investable.
     *          This share weight should be calibrated such that the all of 
     *          the investment gets distributed.
     * \param aPeriod The period for which to get the share weight.
     * \return The share weight for the given period.
     */
    virtual double getShareWeight( const int aPeriod ) const = 0;
};

// Define empty inline methods.
//! Constructor
inline IInvestable::IInvestable(){
}

//! Destructor
inline IInvestable::~IInvestable(){
}

#endif // _IINVESTABLE_H_
