#ifndef _SIMPLE_POLICY_TARGET_RUNNER_H_
#define _SIMPLE_POLICY_TARGET_RUNNER_H_
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
 * \file simple_policy_target_runner.h
 * \ingroup Objects
 * \brief The SimplePolicyTargetRunner class header file.
 * \author Jim Naslund
 */

#include <memory>
#include <xercesc/dom/DOMNode.hpp>
#include <boost/noncopyable.hpp>
#include "containers/include/iscenario_runner.h"

class Timer;
class ITarget;
class Curve;
class TotalPolicyCostCalculator;

typedef std::vector<std::pair<double,double> > VectorOfPairs;

/*!
 * \ingroup Objects
 * \brief A scenario runner which allows the user to specify two emissions
 *        pathways to interpolate between to find a climate target.
 * \details This class runs a scenario multiple times. It reads in two emissions
 *          curves and creates a third curve by performing a interpolating
 *          between the two pathways using a constant. The constant is varied
 *          until the target is reached.
 *
 *          This scenario is controlled by the "simple-find-path" boolean or
 *          directly by the BatchRunner. If it is run independently and not from
 *          the BatchRunner, it reads its configuration from the filename
 *          specified by "sPolicyInputFileName". If it is run from the
 *          BatchRunner the configuration values are parsed from that file
 *          directly. It outputs to the filename specified by
 *          "sPolicyOutputFileName". The default output filename is
 *          sPolicyFinalEmissionsCurve.xml.
 *
 *          <b>XML specification for SimplePolicyTargetRunner</b>
 *          - XML name: \c simple-policy-target-runner
 *          - Contained by: SimplePolicyTargetRunner
 *          - Parsing inherited from class: None.
 *          - Elements:
 *              - \c target-year SimplePolicyTargetRunner::mTargetYear
 *              - \c target-value SimplePolicyTargetRunner::mTargetValue
 *              - \c target-type SimplePolicyTargetRunner::mTargetType
 *              - \c target-tolerance SimplePolicyTargetRunner::mTargetTolerance
 *                (optional)
 *                                    The default is 0.005.
 *              - \c tax-name SimplePolicyTargetRunner::mTaxName
 *                (optional)
 *                                    The default is CO2.
 *              - \c Curve SimplePolicyTargetRunner::mLowerBound
 *                  -Attributes:
 *                      - \c type PointSetCurve
 *                      - \c name lower-bound
 *                      - \c wre wre level of curve
 *                  -Elements:
 *                      - \c PointSet
 *                          -Attributes
 *                              - \c type ExplicitPointSet
 *                          -Elements
 *                              - \c DataPoint (there can be up to 8 of these)
 *                                  -Attributes
 *                                      - \c type XYDataPoint
 *                                  -Elements
 *                                      - \c x x point
 *                                      - \c y y point
 *              - \c Curve SimplePolicyTargetRunner::mUpperBound
 *                  -Attributes:
 *                      - \c type PointSetCurve
 *                      - \c name upper-bound
 *                      - \c wre wre level of curve (this has no function, for
 *                        reference only)
 *                  -Elements:
 *                      - \c PointSet
 *                          -Attributes
 *                              - \c type ExplicitPointSet
 *                          -Elements
 *                              - \c DataPoint (there can be a datapoint for
 *                                each period)
 *                                  -Attributes
 *                                      - \c type XYDataPoint
 *                                  -Elements
 *                                      - \c x x point
 *                                      - \c y y point
 *                           
 *
 * \author Jim Naslund
 */
class SimplePolicyTargetRunner: public IScenarioRunner, protected boost::noncopyable {
    friend class ScenarioRunnerFactory;
public:
    virtual ~SimplePolicyTargetRunner();

    virtual const std::string& getName() const;

    virtual bool setupScenarios( Timer& timer,
        const std::string aName = "",
        const std::list<std::string> aScenComponents = std::list<std::string>() );
    
    virtual bool runScenarios( const int aSingleScenario,
                               const bool aPrintDebugging,
                               Timer& timer );

    virtual void printOutput( Timer& timer,
                              const bool aCloseDB ) const;

    virtual Scenario* getInternalScenario();
    virtual const Scenario* getInternalScenario() const;

    // IParsable Interface.
    bool XMLParse( const xercesc::DOMNode* aRoot );

protected:
    //! The scenario runner which controls running the initial scenario, and all
    //! fixed taxed scenarios after.
    std::auto_ptr<IScenarioRunner> mSingleScenario;

    //! The policy target.
    std::auto_ptr<ITarget> mPolicyTarget;

    //! The delegate object which calculates total costs.
    std::auto_ptr<TotalPolicyCostCalculator> mPolicyCostCalculator;

    //! The type of policy target.
    std::string mTargetType;
    
    //! The name of the tax to modify.
    std::string mTaxName;

    //! The year in which to reach the target.
    unsigned int mTargetYear;

    //! The target value
    double mTargetValue;

    //! Lower Bound Curve
    std::auto_ptr<Curve> mLowerBound;

    //! Upper Bound Curve
    std::auto_ptr<Curve> mUpperBound;

    //! Interpolated Curve
    std::auto_ptr<Curve> mInterpolatedCurve;

    //! Tolerance as a percent
    double mTolerance;

    //! Whether the target runner has already parsed it's data. The XML parse
    //! can be called directly from the BatchRunner and in that case the object
    //! should not parse data from its separate configuration file.
    bool mHasParsedConfig;

    SimplePolicyTargetRunner();
    static const std::string& getXMLNameStatic();

private:

    std::vector<double> curveToConstraintVector( const Curve* aCurve ) const;
    void combineCurves( const std::vector<double>& aEmissionValues ) const;
    void setTrialTaxes( const std::vector<double>& aEmissions );
    static std::vector<double> preComputeDifferences( const VectorOfPairs&,
                                                      const VectorOfPairs& );
    static std::vector<double> preComputeEmissions( const VectorOfPairs& aLower,
                                                    const std::vector<double>& aDifferences,
                                                    const double aConstant );

};

#endif // _SIMPLE_POLICY_TARGET_RUNNER_H_
