#ifndef _TOTAL_POLICY_COST_CALCULATOR_H_
#define _TOTAL_POLICY_COST_CALCULATOR_H_
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
* \file total_policy_cost_calculator.h
* \ingroup Objects
* \brief The TotalPolicyCostCalculator class header file.
* \author Josh Lurz
*/

#include <map>
#include <memory>
#include <vector>

class IScenarioRunner;
class Curve;

/*! 
* \ingroup Objects
* \brief A delegate class which calculates the total cost of a policy given an
*        already run scenario.
* \details This class runs a scenario multiple times while varying a fixed
*          carbon price, to determine the MAC curve and total cost for the
*          scenario.
* \author Josh Lurz
*/
class TotalPolicyCostCalculator {
public:
    explicit TotalPolicyCostCalculator( IScenarioRunner* aSingleScenario );
    ~TotalPolicyCostCalculator();
    bool calculateAbatementCostCurve();
    void printOutput() const;
private:
    //! The total global cost of the policy.
    double mGlobalCost;

    //! The total global cost of the policy discounted at the read-in rate.
    double mGlobalDiscountedCost;

    //! Whether costs have been successfully run.
    bool mRanCosts;

    //! The number of points to use to calculate the marginal abatement curve.
    unsigned int mNumPoints;

    //! The name of the GHG for which to calculate the marginal abatement curve.
    std::string mGHGName;

    //! The scenario runner which controls running the initial scenario, and all
    //! fixed taxed scenarios after. This is a weak reference.
    IScenarioRunner* mSingleScenario;

    typedef std::map<const std::string, double> RegionalCosts;
    typedef RegionalCosts::const_iterator CRegionalCostsIterator;

    //! Total costs indexed by region name.
    RegionalCosts mRegionalCosts;

    //! Total discounted costs indexed by region name.
    RegionalCosts mRegionalDiscountedCosts;

    typedef std::vector<std::map<const std::string, const Curve*> > VectorRegionCurves;
    typedef VectorRegionCurves::iterator VectorRegionCurvesIterator;
    typedef VectorRegionCurves::const_iterator CVectorRegionCurvesIterator;
    typedef std::map<const std::string, const Curve* > RegionCurves;
    typedef RegionCurves::const_iterator CRegionCurvesIterator;
    typedef RegionCurves::iterator RegionCurvesIterator;
    VectorRegionCurves mEmissionsQCurves;
    VectorRegionCurves mEmissionsTCurves;
    VectorRegionCurves mPeriodCostCurves;
    RegionCurves mRegionalCostCurves;

    bool runTrials();
    void createCostCurvesByPeriod();
    void createRegionalCostCurves();
    const std::string createXMLOutputString() const;
    void writeToDB() const;
    void writeToCSV() const;
};
#endif // _TOTAL_POLICY_COST_CALCULATOR_H_
