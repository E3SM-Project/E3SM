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
 * \file solver_factory.cpp
 * \ingroup Objects
 * \brief SolverFactory class source file.
 * \author Pralit Patel
 */

#include "util/base/include/definitions.h"
#include <string>
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>

#include "solution/solvers/include/solver_factory.h"
#include "util/logger/include/ilogger.h"

// Solver subclasses
#include "solution/solvers/include/bisection_nr_solver.h"
#include "solution/solvers/include/bisect_policy_nr_solver.h"
#include "solution/solvers/include/user_configurable_solver.h"

using namespace std;
using namespace xercesc;

/*!
 * \brief Returns whether this factory can create a solver with the given xml
 *        name.
 * \param aXMLName The name of an xml element to check.
 * \return True if this factory has a solver with the given xml name, false otherwise.
 * \note The list of known solvers here needs to be kept in sync with
 *       the ones found in createAndParseSolver.
 */
bool SolverFactory::hasSolver( const string& aXMLName ) {
    return BisectPolicyNRSolver::getXMLNameStatic() == aXMLName
        || BisectionNRSolver::getXMLNameStatic() == aXMLName
        || UserConfigurableSolver::getXMLNameStatic() == aXMLName;
}

/*!
 * \brief Creates and parses the solver with the given xml name.
 * \details Creates the solver and calls XMLParse on it before returning it,
 *          if there are no known solvers which match the given xml name null
 *          is returned.
 * \param aXMLName The element name of the given xml node.
 * \param aNode The xml which defines the solver to be created.
 * \return The newly created and parsed solver or null if given an unknown type.
 * \note The list of known solvers here must be kept in sync with
 *       the ones found in hasSolver.
 */
Solver* SolverFactory::createAndParseSolver( const string& aXMLName, Marketplace* aMarketplace,
                                             World* aWorld, const DOMNode* aNode )
{
    // make sure we know about this solver
    if( !hasSolver( aXMLName ) ) {
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::WARNING );
        mainLog << "Could not create unknown Solver: " << aXMLName << endl;
        return 0;
    }
    
    // create the requested solver
    Solver* retSolver;
    if( BisectPolicyNRSolver::getXMLNameStatic() == aXMLName ) {
        retSolver = new BisectPolicyNRSolver( aMarketplace, aWorld );
    }
    else if( BisectionNRSolver::getXMLNameStatic() == aXMLName ) {
        retSolver = new BisectionNRSolver( aMarketplace, aWorld );
    }
    else if( UserConfigurableSolver::getXMLNameStatic() == aXMLName ) {
        retSolver = new UserConfigurableSolver( aMarketplace, aWorld );
    }
    else {
        // this must mean createAndParseSolver and hasSolver are out of
        // sync with known solvers
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::WARNING );
        mainLog << "Could not create unknown Solver: " << aXMLName
            << ", createAndParseSolver may be out of sync with hasSolver." << endl;
        return 0;
    }
    
    // parse the created solver
    retSolver->XMLParse( aNode );
    return retSolver;
}
