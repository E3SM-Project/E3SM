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
 * \file solver_component_factory.cpp
 * \ingroup Objects
 * \brief SolverComponentFactory class source file.
 * \author Pralit Patel
 */

#include "util/base/include/definitions.h"
#include <string>
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>

#include "solution/solvers/include/solver_component_factory.h"
#include "solution/solvers/include/solver_component.h"
#include "util/logger/include/ilogger.h"

// SolverComponent subclasses
#include "solution/solvers/include/log_newton_raphson.h"
#include "solution/solvers/include/log_newton_raphson_sd.h"
#include "solution/solvers/include/bisect_all.h"
#include "solution/solvers/include/bisect_one.h"
#include "solution/solvers/include/bisect_policy.h"

using namespace std;
using namespace xercesc;

/*!
 * \brief Returns whether this factory can create a solver component with the given
 *        xml name.
 * \param aXMLName The name of an xml element to check.
 * \return True if this factory has a solver component with the given xml name,
 *         false otherwise.
 * \note The list of known solvers components here needs to be kept in sync with
 *       the ones found in createAndParseSolverComponent.
 */
bool SolverComponentFactory::hasSolverComponent( const string& aXMLName ) {
    return LogNewtonRaphson::getXMLNameStatic() == aXMLName
        || BisectAll::getXMLNameStatic() == aXMLName
        || LogNewtonRaphsonSaveDeriv::getXMLNameStatic() == aXMLName
        || BisectOne::getXMLNameStatic() == aXMLName
        || BisectPolicy::getXMLNameStatic() == aXMLName;
}

/*!
 * \brief Creates and parses the solver component with the given xml name.
 * \details Creates the solver component and calls XMLParse on it before returning
 *          it,  if there are no known solver components which match the given xml
 *          name null is returned.
 * \param aXMLName The element name of the given xml node.
 * \param aNode The xml which defines the solver component to be created.
 * \return The newly created and parsed solver component or null if given an unknown type.
 * \note The list of known solvers components here must be kept in sync with
 *       the ones found in hasSolverComponent.
 */
SolverComponent* SolverComponentFactory::createAndParseSolverComponent( const string& aXMLName,
                                                                        Marketplace* aMarketplace,
                                                                        World* aWorld,
                                                                        CalcCounter* aCalcCounter,
                                                                        const DOMNode* aNode )
{
    // make sure we know about this solver component
    if( !hasSolverComponent( aXMLName ) ) {
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::WARNING );
        mainLog << "Could not create unknown SolverComponent: " << aXMLName << endl;
        return 0;
    }
    
    // create the requested solver component
    SolverComponent* retSolverComponent;
    if( LogNewtonRaphson::getXMLNameStatic() == aXMLName ) {
        retSolverComponent = new LogNewtonRaphson( aMarketplace, aWorld, aCalcCounter );
    }
    else if( BisectAll::getXMLNameStatic() == aXMLName ) {
        retSolverComponent = new BisectAll( aMarketplace, aWorld, aCalcCounter );
    }
    else if( LogNewtonRaphsonSaveDeriv::getXMLNameStatic() == aXMLName ) {
        retSolverComponent = new LogNewtonRaphsonSaveDeriv( aMarketplace, aWorld, aCalcCounter );
    }
    else if( BisectOne::getXMLNameStatic() == aXMLName ) {
        retSolverComponent = new BisectOne( aMarketplace, aWorld, aCalcCounter );
    }
    else if( BisectPolicy::getXMLNameStatic() == aXMLName ) {
        retSolverComponent = new BisectPolicy( aMarketplace, aWorld, aCalcCounter );
    }
    else {
        // this must mean createAndParseSolverComponent and hasSolverComponent
        // are out of sync with known solver components
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::WARNING );
        mainLog << "Could not create unknown SolverComponent: " << aXMLName
            << ", createAndParseSolverComponent may be out of sync with hasSolverComponent." << endl;
        return 0;
    }
    
    // parse the created solver component if we have something to parse
    if( aNode ) {
        retSolverComponent->XMLParse( aNode );
    }
    return retSolverComponent;
}
