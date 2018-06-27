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
* \file policy_target_runner.cpp
* \ingroup Objects
* \brief SimplePolicyTargetRunner class source file.
* \author Jim Naslund
*/

#include <cassert>
#include <utility>

#include "util/base/include/definitions.h"

#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>

#include "target_finder/include/simple_policy_target_runner.h"
#include "target_finder/include/itarget.h"
#include "util/base/include/configuration.h"
#include "util/logger/include/ilogger.h"
#include "util/base/include/xml_helper.h"
#include "util/curves/include/point_set_curve.h"
#include "containers/include/scenario_runner_factory.h"
#include "policy/include/policy_ghg.h"
#include "target_finder/include/target_factory.h"
#include "target_finder/include/bisecter.h"
#include "util/curves/include/explicit_point_set.h"
#include "util/curves/include/xy_data_point.h"
#include "containers/include/total_policy_cost_calculator.h"
#include "util/base/include/auto_file.h"

using namespace std;
using namespace xercesc;

extern void closeDB();
extern ofstream outFile;
extern void createMCvarid();

class Curve;

//! Constructor that initializes needed variables
SimplePolicyTargetRunner::SimplePolicyTargetRunner()
: mLowerBound( new PointSetCurve ),
  mUpperBound( new PointSetCurve ),
  mInterpolatedCurve( new PointSetCurve ),
  mTargetValue( -1 ),
  mTargetYear( 0 ),
  mTolerance( 0.005 ),
  mHasParsedConfig( false ),
  mTaxName( "CO2" ) // Default to only taxing CO2.
{
    mSingleScenario = ScenarioRunnerFactory::create( "single-scenario-runner" );

    // Check to make sure calibration is off.
    const Configuration* conf = Configuration::getInstance();
    if( conf->getBool( "debugChecking" ) && conf->getBool( "CalibrationActive" ) ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::WARNING );
        mainLog << "Calibration may be incompatible with target finding." << endl;
    }

    // Initialize the total policy cost calculator if the user requested that
    // total costs should be calculated.
    if( conf->getBool( "createCostCurve" ) ){
        mPolicyCostCalculator.reset( new TotalPolicyCostCalculator( mSingleScenario.get() ) );
    }
}
//! Destructor.
SimplePolicyTargetRunner::~SimplePolicyTargetRunner(){
}

const string& SimplePolicyTargetRunner::getName() const {
    return getXMLNameStatic();
}

/*! \brief Setup the Scenario to be run.
* \details This function sets up the contained SingleScenarioRunner.
* \param aTimer The timer used to print out the amount of time spent performing
*        operations.
* \param aName The name to add on to the name read in in the Configuration file.
* \param aScenComponents A list of additional scenario components to read in.
* \return Whether the setup completed successfully.
*/
bool SimplePolicyTargetRunner::setupScenarios( Timer& aTimer, const string aName, const list<string> aScenComponents ){
    bool success = mSingleScenario->setupScenarios( aTimer, aName, aScenComponents );

    // Get the name of the input file from the Configuration.
    const string fileName = Configuration::getInstance()->getFile( "sPolicyInputFileName", "" );

    if( fileName.empty() ){
        return false;
    }

    // Only read from the configuration file if the data has not already been
    // directly parsed from the BatchRunner configuration file.
    if( !mHasParsedConfig ){
        // Add note so that if XML read fails here user knows what happened
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::NOTICE );
        mainLog << "Reading simple target finder configuration file " << fileName << endl;

        // Parse the file.
        success &= XMLHelper<void>::parseXML( fileName, this );
    }


    // Vector of strings to hold input errors.
    vector<string> errorMsgs;

    // Evaluates to true if all the variables that should have been defined were actually defined.
    // False otherwise.
    if( !( mLowerBound.get() != 0 && mUpperBound.get() != 0 && !mTargetType.empty()
        && mTargetValue != -1 && mTargetYear != 0 ) ){
        errorMsgs.push_back("Missing input variables.");
        success = false;
    }
    // Evaluates to true if both curves have the same number of points.
    if( mLowerBound->getSortedPairs().size() 
        != mUpperBound->getSortedPairs().size() ){
        errorMsgs.push_back("Curves do not have the same number of points.");
        success = false;
    }
    // Evaluates to true if the lower bound has less than or equal to the number of periods
    // It is unnecessary to check the upper bound because the two curves are guaranteed to have
    // the same number of points.
    if( mLowerBound->getSortedPairs().size() 
        > (unsigned int)getInternalScenario()->getModeltime()->getmaxper() ){
        errorMsgs.push_back("Curves have more points than the number of periods.");
        success = false;
    }
    // Evaluates to true if both curves have the same starting and ending point.
    if( ( mLowerBound->getMinX() != mUpperBound->getMinX() )
        || ( mLowerBound->getMaxX() != mUpperBound->getMaxX() ) ){
        errorMsgs.push_back("Curves do not have the same starting and ending point.");
        success = false;
    }

    if( !success ){
        // Output all errors in input to the user.
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::ERROR );
        for(unsigned int i = 0; i < errorMsgs.size(); i++){
            mainLog << errorMsgs[ i ] << endl;
        }
    }

    /*! \pre The target year is initialized. */
    assert( mTargetYear != 0 );

    // Create a a policy target.
    mPolicyTarget = TargetFactory::create( mTargetType,
                                           getInternalScenario()->getClimateModel(),
                                           mTargetValue,
                                           mLowerBound->getMinX() );

    // Initialize the interpolated curve to have the same x points as the
    // lower-bound.
    // Create a PointSet from the lower bound.
    VectorOfPairs points = mLowerBound->getSortedPairs();
    PointSet* xPoints = new ExplicitPointSet();
    for(unsigned int i = 0; i < points.size(); i++){
        xPoints->addPoint( new XYDataPoint( points[ i ].first, 0 ) );
    }
    mInterpolatedCurve.reset( new PointSetCurve( xPoints ) );

    return success;
}

/*! \brief Runs the scenario and adjusts the emissions until a target is reached.
* \details Performs a bisection search on the two curves that are read in.  Determines
*          an emissions pathway that will cause the specified target to be reached in
*          the end of the target period.
* \param aSinglePeriod This parameter is currently ignored.
* \param aPrintDebugging Whether to print debugging information.
* \param aTimer The timer used to print out the amount of time spent performing
*        operations.
* \return Whether all model runs solved successfully.
*/
bool SimplePolicyTargetRunner::runScenarios( const int aSingleScenario,
                                             const bool aPrintDebugging,
                                             Timer& timer ){
    // Search until a limit is reach or the solution is found.
    const unsigned int LIMIT_ITERATIONS = 100;

    // Run the model for all periods
    bool success = mSingleScenario->runScenarios( Scenario::RUN_ALL_PERIODS, false, timer );

    // Construct a bisecter which has an initial trial equal to 1/2
    Bisecter bisecter( mPolicyTarget.get(),
                       mTolerance,
                       0,
                       1,
                       .5,
                       1,
                       mTargetYear );

    // Declare some variables we will be using in the loop.
    vector<double> differences = preComputeDifferences( mLowerBound->getSortedPairs(),
                                                        mUpperBound->getSortedPairs());
    double scalingValue = .5; // initial value

    while( bisecter.getIterations() < LIMIT_ITERATIONS ){

        // Compute the emissions based on the constant c and set those emissions as a 
        // constraint in the model.
        vector<double> emissions = preComputeEmissions( mLowerBound->getSortedPairs(), differences, scalingValue );
        combineCurves( emissions );
        setTrialTaxes( curveToConstraintVector( mInterpolatedCurve.get() ) );

        // Run the applicable periods
        const Modeltime* modeltime = getInternalScenario()->getModeltime();
        // The first x point on the either bound curve is the first year
        unsigned int start = modeltime->getyr_to_per( static_cast<unsigned int>(mLowerBound->getMinX()) );
        unsigned int max = getInternalScenario()->getModeltime()->getmaxper();
        for( unsigned int i = start; i < max; i++ ){
            success &= mSingleScenario->runScenarios( i, false, timer );
        }

        pair<double, bool> trial = bisecter.getNextValue();

        // Check for solution.
        if( trial.second ){
            break;
        }

        scalingValue = trial.first;
    }

    if( bisecter.getIterations() >= LIMIT_ITERATIONS ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::ERROR );
        mainLog << "Exiting target finding search as the iterations limit was reached." << endl;
        success = false;
    }
    else {
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::NOTICE );
        mainLog << "Target value was found by search algorithm in "
                << bisecter.getIterations() << " iterations." << endl;
    }

    mSingleScenario->printOutput( timer, false );

    // Now calculate the abatement curves.
    if( mPolicyCostCalculator.get() ){
        success &= mPolicyCostCalculator->calculateAbatementCostCurve();
    }

    return success;
}

/*! \brief Print the output.
* \details Call the various types of printing routines.
* \param aTimer Scenario timer.
* \param aCloseDB Whether to close the Access database when complete.
*/
void SimplePolicyTargetRunner::printOutput( Timer& timer, const bool aCloseDB ) const {
    if( mPolicyCostCalculator.get() ){
        mPolicyCostCalculator->printOutput();
    }

    // Open the XML output file and write the final emissions pathway to it.  It will
    // automatically be closed.
    AutoOutputFile out( "sPolicyOutputFileName", "sPolicyFinalEmissionsCurve.xml" );
    Tabs tabs;
    mInterpolatedCurve->toInputXMLDerived( *out, &tabs );

    static const bool printDB = Configuration::getInstance()->getBool( "write-access-db", true );
    
    // Close the database.
    if( printDB && aCloseDB ){
        createMCvarid();
        closeDB();
        outFile.close();
    }
}

/*! \brief Get the internal scenario.
* \return The internal scenario.
*/
Scenario* SimplePolicyTargetRunner::getInternalScenario(){
    return mSingleScenario->getInternalScenario();
}

/*! \brief Get the internal scenario.
* \return Constant pointer to the internal scenario.
*/
const Scenario* SimplePolicyTargetRunner::getInternalScenario() const {
    return mSingleScenario->getInternalScenario();
}

// IParsable Interface callback function
bool SimplePolicyTargetRunner::XMLParse( const xercesc::DOMNode* aRoot ){
    // Check for double initialization.
    assert( !mHasParsedConfig );

    // Set the configuration has been parsed.
    mHasParsedConfig = true;

    // assume we were passed a valid node.
    assert( aRoot );

    // get the children of the node.
    DOMNodeList* nodeList = aRoot->getChildNodes();
    bool success = true;
    // loop through the children
    for ( unsigned int i = 0; i < nodeList->getLength(); i++ ){
        DOMNode* curr = nodeList->item( i );
        string nodeName = XMLHelper<string>::safeTranscode( curr->getNodeName() );

        if( nodeName == XMLHelper<void>::text() ) {
            continue;
        }

        // Read the target year.
        else if ( nodeName == "target-year" ){
            mTargetYear = XMLHelper<int>::getValue( curr );
        }

        // Read the target value.
        else if ( nodeName == "target-value" ){
            mTargetValue = XMLHelper<double>::getValue( curr );
        }

        // Read the target type.
        else if ( nodeName == "target-type" ){
            mTargetType = XMLHelper<string>::getValue( curr ) + "-target";
        }

        else if ( nodeName == "target-tolerance" ){
            mTolerance = XMLHelper<double>::getValue( curr );
        }
        else if( nodeName == "tax-name" ){
            mTaxName == XMLHelper<string>::getValue( curr );
        }
        // Read the lower-bound and upper-bound curves
        else if ( nodeName == Curve::getXMLNameStatic() ){
            const string nameOfCurve = XMLHelper<string>::getAttr( curr, "name" );
            if( nameOfCurve == "lower-bound"){
                mLowerBound->XMLParse( curr );
            }
            else if( nameOfCurve == "upper-bound"){
                mUpperBound->XMLParse( curr );
            }
            else {
                ILogger& mainLog = ILogger::getLogger( "main_log" );
                mainLog.setLevel( ILogger::WARNING );
                mainLog << "Unrecognized name of curve: " << nameOfCurve
                        << " found while parsing SimplePolicyTargetRunner." << endl
                        << "Looking for \"lower-bound\" or \"upper-bound\"" << endl;
                success = false;
            }
        }

        // Handle unknown nodes.
        else {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::WARNING );
            mainLog << "Unrecognized string: " << nodeName
                    << " found while parsing SimplePolicyTargetRunner." << endl;
            success = false;
        }
    }
    return success;
}

const string& SimplePolicyTargetRunner::getXMLNameStatic(){
    static const string XML_NAME = "simple-policy-target-runner";
    return XML_NAME;
}

/*! \brief Sets the emissions values passed in into mInterpolatedCurve
*   \param emissionValues precomputed emissions values to use
*/
void SimplePolicyTargetRunner::combineCurves( const vector<double>& aEmissionValues ) const {
    VectorOfPairs points = mInterpolatedCurve->getSortedPairs();

    // Loop through this curve and change the Y values
    unsigned int i = 0;
    for( VectorOfPairs::iterator newCurveIt = points.begin();
         newCurveIt != points.end(); newCurveIt++, i++ ){
            mInterpolatedCurve->setY( newCurveIt->first, aEmissionValues[ i ] );
    }
}


/*! 
* \brief Helper function to compute the difference between the upper and lower
*        bounds.
* \param aLowerPts lower bound points
* \param aUpperPts upper bound points
* \return Vector of differences between upper and lower bound.
*/
vector<double> SimplePolicyTargetRunner::preComputeDifferences(const VectorOfPairs& aLowerPts,
                                            const VectorOfPairs& aUpperPts){
    vector<double> differences;
    for(unsigned int i = 0; i < aLowerPts.size(); i++){
        differences.push_back( aUpperPts[ i ].second - aLowerPts[ i ].second );
    }
    return differences;
}

/*! \brief Helper function to compute the emissions based on the upper and lower bound
*   \param aLower lower bound points
*   \param differences vector of differences between upper and lower bound
*   \param c constant by which to combine the curves
*/
vector<double> SimplePolicyTargetRunner::preComputeEmissions( const VectorOfPairs& aLower,
                                                              const vector<double>& aDifferences,
                                                              const double aConstant ){
    vector<double> emissions;
    for( unsigned int i = 0; i < aDifferences.size(); i++ ){
        emissions.push_back( aLower[ i ].second + aConstant * aDifferences[ i ] );
    }
    return emissions;
}

/*! \brief Converts the Y values of a curve into a vector of doubles.  The x index of
*          the curve will correspond to the appropriate period.
* \param aCurve the curve to convert.
* \return vector of doubles
*/
vector<double> SimplePolicyTargetRunner::curveToConstraintVector( const Curve* aCurve ) const {
    // Get a pointer to the modelTime
    const Modeltime* modelTime = mSingleScenario->getInternalScenario()->getModeltime();

    VectorOfPairs points = aCurve->getSortedPairs();
    vector<double> constraint( modelTime->getmaxper() );

    // Set the values of the non-computed periods to -1 because they are not being used.
    // GHGPolicy.completeInit will not do anything with constraints that
    // are -1.
    // The first x point on the either bound curve is the first year
    unsigned int max = modelTime->getyr_to_per( static_cast<unsigned int>(mLowerBound->getMinX()) );
    for(unsigned int i = 0; i < max; i++){
        constraint[ i ] = -1;
    }

    // Loop through the vector of pairs
    for(unsigned int i = 0; i < points.size(); i++){
        // Assign constraint[period of this points year] to this points value.
        constraint[ modelTime->getyr_to_per( ( int)points[i].first ) ] = points[ i ].second;
    }
    return constraint;
}

/*! \brief Set the trial emissions.
* \param aEmmissions emissions vector
*/
void SimplePolicyTargetRunner::setTrialTaxes( const vector<double>& aEmissions ) {
    // Set the fixed constraint into the world. The world will clone this tax object,
    // this object retains ownership of the original.
    GHGPolicy tax( mTaxName, "global" );
    tax.setConstraint( aEmissions );
    mSingleScenario->getInternalScenario()->setTax( &tax );
}
