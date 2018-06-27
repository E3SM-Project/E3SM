#ifndef _CCSM_SCENARIO_H_
#define _CCSM_SCENARIO_H_
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
* \file ccsm_scenario.h
* \ingroup Objects
* \brief The Ccsm_Scenario class header file.
* \author Sonny Kim
*/
#include <vector>
#include <map>
#include <memory>
#include <string>
#include <boost/shared_ptr.hpp>
#include "containers/include/scenario.h"

// Forward declarations
class Modeltime;
class Scenario;
class Marketplace;
class World;
class Curve;
class Tabs;
class Solver;
class GHGPolicy;
class IClimateModel;
class OutputMetaData;
class SolutionInfoParamParser;

/*!
* \ingroup Objects
* \brief A class which defines a model ccsm_scenario.
* \details The Ccsm_Scenario class object inherits and redefines the run method
*          in the base class.
* \author John Truesdale
*/

class Ccsm_Scenario: public Scenario
{
public:
    Ccsm_Scenario();
    ~Ccsm_Scenario();
    bool run( const int aSinglePeriod, const bool aPrintDebugging, const std::string& aFilenameEnding = "" );
private:
    Tabs *scenarioTabs;
    std::ofstream XMLDebugFile;
    std::ofstream SGMDebugFile;
};

#endif // _CCSM_SCENARIO_H_

