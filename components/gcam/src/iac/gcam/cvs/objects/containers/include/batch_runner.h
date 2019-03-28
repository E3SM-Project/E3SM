#ifndef _BATCH_RUNNER_H_
#define _BATCH_RUNNER_H_
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
 * \file batch_runner.h
 * \ingroup Objects
 * \brief The BatchRunner class header file.
 * \author Josh Lurz
 */

#include <string>
#include <list>
#include <memory>
#include <xercesc/dom/DOMNode.hpp>
#include "containers/include/iscenario_runner.h"
class Timer;

/*! 
 * \ingroup Objects
 * \brief A scenario runner which reads in set of input files and runs multiple
 *        scenarios by determining all possible permutations of these sets of
 *        files.
 * \details This class runs multiple scenarios by creating combinations of the
 *          sets of files within each component. For each scenario, one file set
 *          from each component will be chosen. It will continue running
 *          scenarios until all possible combinations have been run. This may
 *          also be used to run a simple set of scenarios by defining a single
 *          component with one file set for each run.
 *
 *          A special type of component set may also be read which is the
 *          runner-set. Only a single runner set may be read, and each element
 *          in the runner set represents a type of ScenarioRunner to use.
 *          ScenarioRunners in the runner-set will be used similarly to
 *          FileSets, permutations will be created for each member. Only one
 *          runner-set may be read, if none are read a ScenarioRunner will be
 *          created from the values in the configuration file as would occur in
 *          a single scenario run.
 *          
 *          The batch runner is turned on using the boolean configuration value
 *          "BatchMode". The name of the configuration file is determined by the
 *          file configuration value "BatchFileName".
 *
 *          <b>XML specification for BatchRunner</b>
 *          - XML name: \c BatchRunner
 *          - Contained by: None.
 *          - Parsing inherited from class: None.
 *          - Elements:
 *              - \c ComponentSet BatchRunner::ComponentSet
 *                  - Attributes:
 *                      - \c name Name of the ComponentSet.
 *                  - Elements:
 *                      - \c FileSet BatchRunner::FileSet
 *                          - Attributes:
 *                              - \c name Name of the FileSet.
 *                          - Elements:
 *                              - \c %Value Path to a single file.
 *                                  - Attributes:
 *                                      - \c name Name of the file.
 *              - \c runner-set A set of ScenarioRunners to use. The contents of
 *                              all runner-set elements are combined into a single
 *                              set of scenario runners.
 *                  - Elements:
 *                      - \c %Value Path to a scenario runner configuration file.
 *                      - \c merge-runner MergeRunner
 *                      - \c single-scenario-runner SingleScenarioRunner
 *                      - \c mac-generator-scenario-runner MACGeneratorScenarioRunner
 *                      - \c policy-target-runner PolicyTargetRunner
 *                      - \c simple-policy-target-runner SimplePolicyTargetRunner
 *
 * \author Josh Lurz
 */
class BatchRunner: public IScenarioRunner {
	friend class ScenarioRunnerFactory;
public:
    virtual ~BatchRunner();

    virtual const std::string& getName() const;

    virtual bool setupScenarios( Timer& aTimer,
                                 const std::string aName = "",
                                 const std::list<std::string> aScenComponents = std::list<std::string>() );

    virtual bool runScenarios( const int aSinglePeriod,
                               const bool aPrintDebugging,
                               Timer& aTimer );

    virtual void printOutput( Timer& aTimer, const bool aCloseDB ) const;

	virtual Scenario* getInternalScenario();

	virtual const Scenario* getInternalScenario() const;

	// IParsable Interface.
    bool XMLParse( const xercesc::DOMNode* aRoot );
protected:
    //! A structure which defines a single file.
    struct File {
        //! The name of the file.
        std::string mName;

        //! The path to the file.
        std::string mPath;
    };

    //! A structure which defines a single named set of files.
    struct FileSet {
        //! The set of files.
        std::list<File> mFiles;

        //! The name for the set of files.
        std::string mName;
    };

    //! A structure which defines a single named component consisting of
    //! multiple file sets.
    struct Component {
        //! The list of contained file sets.
        std::list<FileSet> mFileSets;

        //! An iterator to the current file set.
        std::list<FileSet>::const_iterator mFileSetIterator;

        //! The name of the component set.
        std::string mName;
    };
    
    //! A vector containing a series of Components.
    typedef std::vector<Component> ComponentSet;

    //! The set of all components.
    ComponentSet mComponentSet;

    //! A set of IScenarioRunners to use. Permutation will be created to use
    //! each IScenarioRunner with each scenario.
    std::list<IScenarioRunner*> mScenarioRunners;
    
    //! The names of all scenarios which did not solve.
    std::list<std::string> mUnsolvedNames;

    //! The current scenario runner.
    IScenarioRunner* mInternalRunner;

	BatchRunner();
	bool runSingleScenario( IScenarioRunner* aScenarioRunner,
                            const Component& aCurrComponent,
                            Timer& aTimer );

    bool XMLParseComponentSet( const xercesc::DOMNode* aNode );

    bool XMLParseRunnerSet( const xercesc::DOMNode* aNode );

    bool XMLParseFileSet( const xercesc::DOMNode* aNode,
                          Component& aCurrComponentSet );

	static const std::string& getXMLNameStatic();

    /*! 
     * \brief Assists with creating and parsing XML data for IScenarioRunners.
     * \details Used to peek into an XML file and create the type of
     *          IScenarioRunner requested. XML parsing for the class is then
     *          forwarded to the newly created IScenarioRunner.
     */
    class ParseHelper: public IParsable {
    public:
        // IParsableInterface
        bool XMLParse( const xercesc::DOMNode* aRoot );

        std::auto_ptr<IScenarioRunner>& getParsedScenarioRunner();
    private:
        //! The IScenarioRunner created by the XMLParse method.
        std::auto_ptr<IScenarioRunner> mScenarioRunner;
    };

};
#endif // _BATCH_RUNNER_H_
