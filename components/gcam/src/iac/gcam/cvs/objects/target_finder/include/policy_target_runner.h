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

#ifndef _POLICY_TARGET_RUNNER_H_
#define _POLICY_TARGET_RUNNER_H_
#if defined(_MSC_VER)
#pragma once
#endif

/*!
 * \file policy_target_runner.h
 * \ingroup Objects
 * \brief The PolicyTargetRunner class header file.
 * \author Josh Lurz
 */

#include <memory>
#include <vector>
#include "containers/include/iscenario_runner.h"
#include "util/base/include/value.h"

class Timer;
class TotalPolicyCostCalculator;
class ITarget;
class Modeltime;

/*! 
 * \ingroup Objects
 * \brief A scenario runner that determines an optimal Peck-Wan carbon price
 *        pathway.
 * \details This scenario runner runs the scenario multiple times to determine a
 *          near optimal carbon price pathway. The input parameters for the
 *          pathway are an initial tax year, an interest rate, and a final
 *          climate parameter level, such as temperature or concentration. The
 *          model will construct a Hotelling price path such that the climate
 *          parameter level increases until the target peaks at the specified
 *          level, and then the climate parameter is held constant.
 *
 *          An "overshoot" is also allowed where the hotelling price path is 
 *          increased until the climate parameter hits the target in the desired
 *          year.
 *
 *          This scenario is controlled by the "find-path" boolean or directly
 *          by the BatchRunner. If it is run independently and not from the
 *          BatchRunner, it reads its configuration from the filename specified
 *          by "policy-target-file". If it is run from the BatchRunner the
 *          configuration values are parsed from that file directly.
 *
 * \note The LUE feedback must be disconnected for this to work using the bisection
 *       solver, otherwise at high carbon prices the land use change emissions cause
 *       the sign of the derivative on total emissions with respect to the carbon price
 *       to changes signs, which defeats the bisection solution mechanism.  Using
 *       the Secant method should work.
 *
 *          <b>XML specification for PolicyTargetRunner</b>
 *          - XML name: \c policy-target-runner
 *          - Contained by: None or BatchRunner.
 *          - Parsing inherited from class: None.
 *          - Attributes:
 *              - \c name PolicyTargetRunner::mName
 *          - Elements:
 *              - \c target-value PolicyTargetRunner::mTargetValue
 *              - \c target-type PolicyTargetRunner::mTargetType
 *                   (optional) The default is concentration.
 *              - \c tax-name PolicyTargetRunner::mTaxName
 *                   (optional) The default is CO2.
 *              - \c target-tolerance PolicyTargetRunner::mTolerance
 *                   (optional) The default is 0.01.
 *              - \c path-discount-rate PolicyTargetRunner::mPathDiscountRate
 *                   (optional) The default is 0.05.
 *              - \c first-tax-year PolicyTargetRunner::mFirstTaxYear
 *                   (optional) The default is 2020.
 *              - \c max-iterations PolicyTargetRunner::mMaxIterations
 *                   (optional) The default is 100.
 *              - \c stabilization PolicyTargetRunner::mInitialTargetYear
 *                   (optional) Set the initial target year to the flag
 *                   ITarget::getUseMaxTargetYearFlag(), this is the default.
 *              - \c overshoot PolicyTargetRunner::mInitialTargetYear
 *                   (optional) Set the initial target year to the value of the
 *                   year attribute or the last model year if that attribute is
 *                   not specified.
 *
 * \author Josh Lurz
 * \author Pralit Patel
 */
class PolicyTargetRunner: public IScenarioRunner {
    friend class ScenarioRunnerFactory;
public:
    // IParsable interface
    virtual bool XMLParse( const xercesc::DOMNode* aRoot );

    virtual ~PolicyTargetRunner();

    virtual const std::string& getName() const;

    virtual bool setupScenarios( Timer& timer,
        const std::string aName = "",
        const std::list<std::string> aScenComponents =
            std::list<std::string>() );
    
    virtual bool runScenarios( const int aSingleScenario,
                              const bool aPrintDebugging,
                              Timer& timer );

    virtual void printOutput( Timer& timer,
        const bool aCloseDB = true ) const;

    virtual Scenario* getInternalScenario();
    virtual const Scenario* getInternalScenario() const;
private:
    //! The scenario runner which controls running the initial scenario, and all
    //! fixed taxed scenarios after.
    std::auto_ptr<IScenarioRunner> mSingleScenario;

    //! The delegate object which calculates total costs.
    std::auto_ptr<TotalPolicyCostCalculator> mPolicyCostCalculator;

    //! The name of the policy target runner.
    std::string mName;

    //! The type of policy target.
    std::string mTargetType;

    //! The name of the tax.
    std::string mTaxName;

    //! Path discount rate. This is the interest rate used to determine the
    //! Hotelling tax path.
    Value mPathDiscountRate;

    //! The tolerance as a percent.
    Value mTolerance;

    //! The target for the climate parameter.
    Value mTargetValue;

    //! The first year to tax.
    unsigned int mFirstTaxYear;
    
    //! The maximum number of bisection iterations to perform when determining the
    //! initial tax or the trial tax in a single future period past the
    //! stabilization year.
    unsigned int mMaxIterations;
    
    //! Initial target year if we are doing an overshoot or the flag
    //! ITarget::getUseMaxTargetYearFlag() if we are doing a stabilization.
    int mInitialTargetYear;

    //! Whether the target runner has already parsed its data. The XML parse
    //! can be called directly from the BatchRunner and in that case the object
    //! should not parse data from its separate configuration file.
    bool mHasParsedConfig;
    
    //! The number of periods to forward look when trying to stay on target
    int mNumForwardLooking;

    static std::vector<double>
        calculateHotellingPath( const double aIntialTax,
                                const double aHotellingRate,
                                const Modeltime* aModeltime,
                                const int aInitialYear,
                                const int aFinalYear );

    void setTrialTaxes( const std::vector<double> aTaxes );
    
    bool solveInitialTarget( std::vector<double>& aTaxes,
                             const ITarget* aPolicyTarget,
                             const unsigned int aLimitIterations,
                             const double aTolerance,
                             Timer& aTimer );
    
    bool solveFutureTarget( std::vector<double>& aTaxes,
                            const ITarget* aPolicyTarget,
                            const unsigned int aLimitIterations,
                            const double aTolerance,
                            const int aPeriod,
                            Timer& aTimer );
    
    bool skipFuturePeriod( std::vector<double>& aTaxes,
                           const ITarget* aPolicyTarget,
                           const unsigned int aLimitIterations,
                           const double aTolerance,
                           const int aFirstSkippedPeriod,
                           const int aPeriod,
                           Timer& aTimer );
    PolicyTargetRunner();
    static const std::string& getXMLNameStatic();
};
#endif // _POLICY_TARGET_RUNNER_H_
