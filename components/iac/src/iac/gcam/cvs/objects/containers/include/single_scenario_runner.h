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

#ifndef _SINGLE_SCENARIO_RUNNER_H_
#define _SINGLE_SCENARIO_RUNNER_H_
#if defined(_MSC_VER)
#pragma once
#endif

/*! 
 * \file single_scenario_runner.h
 * \ingroup Objects
 * \brief The SingleScenarioRunner class header file.
 * \author Josh Lurz
 */

#include <memory>
#include <list>
#include "containers/include/iscenario_runner.h"

class Timer;
class Scenario;

/*! 
 * \ingroup Object
 * \brief A class which is responsible for running a single scenario.
 * \details This class is responsible for running a single scenario. The
 *          SingleScenarioRunner will initialize a new Scenario each time
 *          setupScenarios is called. This will reset all input data and create
 *          a new run identifier for use in the output databases. Inputs read in
 *          last will override any earlier inputs. Inputs are read in the
 *          following order:
 *          -# Base input file as specified in the configuration file.
 *          -# Scenario components as specified in the configuration file, in
 *             top to bottom order in the scenario components section of the
 *             file.
 *          -# Scenario components passed into the function, in the order they
 *             are passed in.
 *
 *          setupScenarios must be called before runScenarios. runScenarios may
 *          be called multiple times, as in the case when total policy costs are
 *          calculated. Care should be taken to ensure that operations that
 *          manipulate the structure of the model after setupScenarios is called
 *          correctly reinitialize any existing data correctly, such as
 *          resetting markets.
 *
 *          printOutput is called after runScenarios to print output to any
 *          configured databases.
 *          
 *          The getInternalScenarios functions only return a valid scenario
 *          after setupScenarios is called.
 *
 * \todo What should be documented here vs. the IScenarioRunner interface.
 *
 * \author Josh Lurz
 */
class SingleScenarioRunner: public IScenarioRunner {
    friend class ScenarioRunnerFactory;
public:
    virtual const std::string& getName() const;

    // IParsable interface
    virtual bool XMLParse( const xercesc::DOMNode* aRoot );

    virtual ~SingleScenarioRunner();

    virtual bool setupScenarios( Timer& timer, const std::string aName = "",
                                 const std::list<std::string> aScenComponents =
                                   std::list<std::string>() );
    
    virtual bool runScenarios( const int aSinglePeriod,
                               const bool aPrintDebugging,
                               Timer& aTimer );
    
    virtual void printOutput( Timer& timer, const bool aCloseDB = true ) const;

    virtual Scenario* getInternalScenario();
    virtual const Scenario* getInternalScenario() const;

protected:    
    SingleScenarioRunner();
    static const std::string& getXMLNameStatic();
    //! The scenario which will be run.
    std::auto_ptr<Scenario> mScenario;
};
#endif // _SINGLE_SCENARIO_RUNNER_H_
