#ifndef _MAC_GENERATOR_SCENARIO_RUNNER_H_
#define _MAC_GENERATOR_SCENARIO_RUNNER_H_
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
* \file mac_generator_scenario_runner.h
* \ingroup Objects
* \brief The MACGeneratorScenarioRunner class header file.
* \author Josh Lurz
*/

#include <memory>
#include "containers/include/iscenario_runner.h"

class Timer;
class TotalPolicyCostCalculator;

/*! 
* \ingroup Objects
* \brief A derived ScenarioRunner class that runs a scenario multiple times in
*        order to generate a marginal abatement cost (MAC) curve for each time
*        period.
* \details This class runs a scenario multiple times while varying a fixed
*          carbon price, to determine the MAC curve and total cost for the
*          scenario.
* \author Josh Lurz
*/
class MACGeneratorScenarioRunner: public IScenarioRunner {
    friend class ScenarioRunnerFactory;
public:
    virtual ~MACGeneratorScenarioRunner();

    virtual const std::string& getName() const;

    // IParsable interface
    virtual bool XMLParse( const xercesc::DOMNode* aRoot );

    virtual bool setupScenarios( Timer& timer,
        const std::string aName = "",
        const std::list<std::string> aScenComponents = std::list<std::string>() );
    
    virtual bool runScenarios( const int aSinglePeriod,
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

    MACGeneratorScenarioRunner();
    static const std::string& getXMLNameStatic();
};
#endif // _MAC_GENERATOR_SCENARIO_RUNNER_H_
