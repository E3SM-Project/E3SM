#ifndef _ISCENARIO_RUNNER_H_
#define _ISCENARIO_RUNNER_H_
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
 * \file iscenario_runner.h
 * \ingroup Objects
 * \brief The IScenarioRunner interface header file.
 * \author Josh Lurz
 */
#include <list>
#include <string>
#include "util/base/include/iparsable.h"

class Timer;
class Scenario;

/*! 
 * \ingroup Objects
 * \brief An interface to a class which is responsible for running a scenario or
 *        set of scenarios.
 * \details An object implementing this interface defines a series of methods
 *          for setting up, running, and printing results for a set of
 *          scenarios. The method for generating the scenarios is implemented by
 *          the object, and not defined by the interface.
 * \author Josh Lurz
 */
class IScenarioRunner: public IParsable {
public:
    /*!
     * \brief Virtual destructor so that derived instances can be deleted
     *        through the base class pointer.
     */
    virtual ~IScenarioRunner();

    /*!
     * \brief Get the name of the scenario runner.
     */
    virtual const std::string& getName() const = 0;

    /*!
     * \brief Setup the ScenarioRunner before running a single or series of
     *        scenarios.
     * \details This method must be called to setup the scenario runner before
     *          the runScenarios method is called. This performs initialization
     *          for the scenario runner.
     * \param aTimer A reference to the global timer.
     * \param aName Name of the scenario or set of scenarios. Defaults to the
     *        empty string.
     * \param aScenComponents A list of locations of scenario components.
     *        Defaults to an empty list.
     * \return Whether the ScenarioRunner was setup successfully.
     */
    virtual bool setupScenarios( Timer& aTimer,
                                 const std::string aName = std::string(),
                                 const std::list<std::string> aScenComponents = std::list<std::string>() ) = 0;
    
    /*!
     * \brief Run the scenario or set of scenarios.
     * \details Run the scenarios in a manner defined by the type of the
     *          scenario runner.
     * \param aSinglePeriod The single period to run or
     *        Scenario::RUN_ALL_PERIODS to run all periods.
     * \param aPrintDebugging Whether to print debugging information during the
     *        periods. This is only respected for the initial run if the same
     *        input set is being run multiple times.
     * \param aTimer Reference to the global timer.
     * \return Whether the scenario or set of scenarios ran successfully.
     */
    virtual bool runScenarios( const int aSinglePeriod,
                               const bool aPrintDebugging,
                               Timer& aTimer ) = 0;
    
    /*!
     * \brief Print the output from the set of scenarios run.
     * \details Print the output of the scenario runs.
     * \param aTimer Reference to the global timer.
     * \param aCloseDB Whether to close the database. Defaults to true.
     * \todo Rename the aCloseDB parameter to something more general.
     */
    virtual void printOutput( Timer& aTimer,
                              const bool aCloseDB = true ) const = 0;

    /*!
     * \brief Get the a mutable reference to the internal scenario object.
     * \details At any given time there is only one Scenario object within the
     *          model that is performing runs. This function returns a pointer
     *          to that internal scenario. Since scenario runners may contain
     *          other scenario runners, this often requires searching down
     *          through several levels.
     * \return The internal scenario.
     */
    virtual Scenario* getInternalScenario() = 0;


    /*!
     * \brief Get a constant reference to the internal scenario object.
     * \details At any given time there is only one Scenario object within the
     *          model that is performing runs. This function returns a pointer
     *          to that internal scenario. Since scenario runners may contain
     *          other scenario runners, this often requires searching down
     *          through several levels.
     * \return Constant pointer to the internal scenario.
     */
    virtual const Scenario* getInternalScenario() const = 0;
};

// Inline destructor.
inline IScenarioRunner::~IScenarioRunner(){}

#endif // _SCENARIO_RUNNER_H_
