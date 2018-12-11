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
* \file node_input.cpp
* \ingroup Objects
* \brief The NodeInput class source file.
* \author Pralit Patel
* \author Ron Sands
*/

#include "util/base/include/definitions.h"
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>

#include "functions/include/node_input.h"
#include "functions/include/ifunction.h"
#include "functions/include/production_input.h"
#include "functions/include/demand_input.h"
#include "functions/include/trade_input.h"
#include "containers/include/scenario.h"
#include "util/base/include/xml_helper.h"
#include "functions/include/function_utils.h"
#include "functions/include/function_manager.h"
#include "util/base/include/ivisitor.h"

// I could get rid of these two if I can figure out a better way to
// create the utility demand fn markets
#include "marketplace/include/marketplace.h"
#include "util/base/include/configuration.h"

using namespace std;
using namespace xercesc;

extern Scenario* scenario;

//! Default Constructor
NodeInput::NodeInput(){
}

//! Destructor
NodeInput::~NodeInput() {
    for( CNestedInputIterator it = mNestedInputs.begin(); it != mNestedInputs.end(); ++it ) {
        delete *it;
    }
}

void NodeInput::XMLParse( const xercesc::DOMNode* node ) {
    /*! \pre make sure we were passed a valid node. */
    assert( node );

    // get the name attribute.
    mName = XMLHelper<string>::getAttr( node, "name" );

    // get all child nodes.
    const DOMNodeList* nodeList = node->getChildNodes();

    // loop through the child nodes.
    for( unsigned int i = 0; i < nodeList->getLength(); i++ ){
        const DOMNode* curr = nodeList->item( i );
        if( curr->getNodeType() == DOMNode::TEXT_NODE ){
            continue;
        }
        const string nodeName = XMLHelper<string>::safeTranscode( curr->getNodeName() );

        if ( nodeName == ProductionInput::getXMLNameStatic() ) {
            parseContainerNode( curr, mNestedInputs, new ProductionInput() );
        }
        else if ( nodeName == DemandInput::getXMLNameStatic() ) {
            parseContainerNode( curr, mNestedInputs, new DemandInput() );
        }
        else if ( nodeName == TradeInput::getXMLNameStatic() ) {
            parseContainerNode( curr, mNestedInputs, new TradeInput() );
        }
        else if ( nodeName == NodeInput::getXMLNameStatic() ) {
            parseContainerNode( curr, mNestedInputs, new NodeInput() );
        }
        else if ( nodeName == "prodDmdFnType" ) {
            mProdDmdFnType = XMLHelper<string>::getValue( curr );
        }
        else if ( nodeName == "price-recieved" ) {
            mPricePaid.set( XMLHelper<double>::getValue( curr ) );
        }
        else if ( nodeName == "Sigma1" ) {
            mSigmaNewCapital.set( XMLHelper<double>::getValue( curr ) );
        }
        else if ( nodeName == "Sigma2" ) {
            mSigmaOldCapital.set( XMLHelper<double>::getValue( curr ) );
        }
        // TODO: these shouldn't really be in here, they are specific to the UtilityDemandFunction
        else if ( nodeName == "alpha-utility-param" ) {
            mAlphaUtilityParam.set( XMLHelper<double>::getValue( curr ) );
        }
        else if ( nodeName == "beta-utility-param" ) {
            mBetaUtilityParam.set( XMLHelper<double>::getValue( curr ) );
        }
        else if ( nodeName == "gamma-utility-param" ) {
            // putting gamma in here due to hack in UtilityDemandFunction
            mAlphaCoef.set( XMLHelper<double>::getValue( curr ) );
        }
        else if ( nodeName == "technicalChange" ) {
            mTechChange.set( XMLHelper<double>::getValue( curr ) );
        }
        // TODO: do I need a derived?
        else /*if( !XMLDerivedClassParse( nodeName, curr ) )*/{
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::WARNING );
            mainLog << "Unrecognized text string: " << nodeName << " found while parsing " 
                << getXMLReportingName() << "." << endl;
        }
    }
}

void NodeInput::completeInit( const string& aRegionName,
                             const string& aSectorName,
                             const string& aSubsectorName,
                             const string& aTechName,
                             DependencyFinder* aDependencyFinder , const IInfo* aTechInfo)
{
    // This is a hack to create the trial utility market when that is the function set
    // it would ideally be created somwhere else such as with in the function itself
    // however it does not have a completeInit method.
    if( !mProdDmdFnType.empty() ) {
        mProdDmdFn = FunctionManager::getFunction( mProdDmdFnType );
        // TODO: a better way to create this market that is less hackish
        if( mProdDmdFnType == "UtilityDemandFunction" ) {
            Marketplace* marketplace = scenario->getMarketplace();
            const string utilityMarketName = aSectorName+"-utility";
            bool createdUtilityMarket = marketplace->createMarket( aRegionName, aRegionName, utilityMarketName,
                IMarketType::NORMAL );
            // it may not have created the market if it was created by another input in say an earlier
            // consumer
            if( createdUtilityMarket && !Configuration::getInstance()->getBool( "CalibrationActive" ) ){
                const Modeltime* modeltime = scenario->getModeltime();
                for( int period = 0; period < modeltime->getmaxper(); ++period ){
                    marketplace->setMarketToSolve( utilityMarketName, aRegionName, period );
                }
            }
        }
    }
    // Initially set the current sigma to new capital even though we technically have old vintage
    // technologies in the base year.  Theoritically the base year should be read in balanced
    // and thus the sigmas we use there should not matter.
    mCurrentSigma = mSigmaNewCapital;

    // we sort the child inputs now to make things easier on us when it comes time to merge
    // node inputs in copyParamsInto.
    InputNameComparator comp;
    sort( mNestedInputs.begin(), mNestedInputs.end(), comp );

    // have all contained inputs do completeInit as well
    for( NestedInputIterator it = mNestedInputs.begin(); it != mNestedInputs.end(); ++it ) {
        (*it)->completeInit( aRegionName, aSectorName, aSubsectorName, aTechName, 
            aDependencyFinder, aTechInfo );
    }
}

void NodeInput::initCalc( const string& aRegionName,
                         const string& aSectorName,
                         const bool aIsNewInvestmentPeriod,
                         const bool aIsTrade,
                         const int aPeriod )
{
    /*!
     * \pre We must have a demand fn by now.
     */
    assert( !mProdDmdFnType.empty() && mProdDmdFn );

    // do I need isRoot?
    /*if( mName == "root" ) {
        mPricePaid.init( FunctionUtils::getPriceReceived( aRegionName, aSectorName, aPeriod ) );
    }*/
    // create a cache of the child INestedInput* as IInput* since the production function
    // needs a vector of those and the compiler can not be conviced it is safe to use a
    // vector of the subclass INestedInput* instead of IInput*
    mChildInputsCache.clear();
    mChildInputsCache.reserve( mNestedInputs.size() );
    for( NestedInputIterator nestedInputIter = mNestedInputs.begin();
        nestedInputIter != mNestedInputs.end(); ++nestedInputIter )
    {
        (*nestedInputIter)->initCalc( aRegionName, aSectorName, aIsNewInvestmentPeriod, 
                                      aIsTrade, aPeriod );
        mChildInputsCache.push_back( *nestedInputIter );
    }
    // initialized the hack
    mNodePriceSet = false;
}

void NodeInput::copyParam( const IInput* aInput,
                          const int aPeriod )
{
    aInput->copyParamsInto( *this, aPeriod );
}

void NodeInput::copyParamsInto( NodeInput& aInput, const int aPeriod ) const {
    /*!
     * \pre We should be copying into inputs of the same name.
     */
    assert( mName == aInput.mName );

    aInput.mAlphaCoef.init( mAlphaCoef );
    aInput.mPricePaid.init( mPricePaid );
    aInput.mSigmaNewCapital.init( mSigmaNewCapital );
    aInput.mSigmaOldCapital.init( mSigmaOldCapital );
    aInput.mCurrentSigma.set( mCurrentSigma );
    aInput.mProdDmdFnType = mProdDmdFnType;
    aInput.mProdDmdFn = mProdDmdFn;

    aInput.mBasePricePaid.set( mBasePricePaid );

    aInput.mAlphaUtilityParam = mAlphaUtilityParam;
    aInput.mBetaUtilityParam = mBetaUtilityParam;

    // copy children
    CNestedInputIterator itThis = mNestedInputs.begin();
    NestedInputIterator itArg = aInput.mNestedInputs.begin();
    while( itArg != aInput.mNestedInputs.end() && itThis != mNestedInputs.end() ) {
        if( (*itThis)->getName() == (*itArg)->getName() ) {
            // exists in both so just copy the params
            (*itArg)->copyParam( *itThis, aPeriod );
            ++itArg;
            ++itThis;
        }
        else if ( (*itArg)->getName() < (*itThis)->getName() ) {
            // did not exist previously so erase it
            delete *itArg;
            itArg = aInput.mNestedInputs.erase( itArg );
        }
        else {
            // does not have this input so clone it
            itArg = 1 +
                aInput.mNestedInputs.insert( itArg, static_cast<INestedInput*>( (*itThis)->clone() ) );
            ++itThis;
        }
    }
    // anything left in the argument's children gets erased
    aInput.mNestedInputs.erase( itArg, aInput.mNestedInputs.end() );
    
    // anything left in this input's children gets copied
    for( ; itThis != mNestedInputs.end(); ++itThis ) {
        aInput.mNestedInputs.push_back( static_cast<INestedInput*>( (*itThis)->clone() ) );
    }
}

IInput* NodeInput::clone() const {
    NodeInput* retNodeInput = new NodeInput;
    retNodeInput->copy( *this );
    return retNodeInput;
}

void NodeInput::copy( const NodeInput& aNodeInput ) {
    mName = aNodeInput.mName;
    mSigmaNewCapital.init( aNodeInput.mSigmaNewCapital );
    mSigmaOldCapital.init( aNodeInput.mSigmaOldCapital );
    mCurrentSigma.set( aNodeInput.mCurrentSigma );
    mAlphaCoef.init( aNodeInput.mAlphaCoef );
    mPricePaid.init( aNodeInput.mPricePaid ); 
    mProdDmdFnType = aNodeInput.mProdDmdFnType;
    mProdDmdFn = aNodeInput.mProdDmdFn;
    mBasePricePaid.set( aNodeInput.mBasePricePaid );

    mAlphaUtilityParam = aNodeInput.mAlphaUtilityParam;
    mBetaUtilityParam = aNodeInput.mBetaUtilityParam;

    // copy children
    for( CNestedInputIterator it = aNodeInput.mNestedInputs.begin(); it != aNodeInput.mNestedInputs.end(); ++it ) {
        mNestedInputs.push_back( static_cast<INestedInput*>( (*it)->clone() ) );
    }
}

bool NodeInput::isSameType( const string& aType ) const {
    return aType == getXMLReportingName();
}

bool NodeInput::hasTypeFlag( const int aTypeFlag ) const {
    // TODO: do node inputs have a type? Just doing false for now
    return false;
}

//! Output to XML data
void NodeInput::toInputXML( ostream& aOut, Tabs* aTabs ) const {
    // write the beginning tag.
    XMLWriteOpeningTag ( getXMLReportingName(), aOut, aTabs, mName );
    
    XMLWriteElement( mSigmaNewCapital, "Sigma1", aOut, aTabs );
    XMLWriteElement( mSigmaOldCapital, "Sigma2", aOut, aTabs );

    //toInputXMLDerived( out, tabs );
    XMLWriteElement( mProdDmdFnType, "prodDmdFnType", aOut, aTabs );
    for( CNestedInputIterator it = mNestedInputs.begin(); it != mNestedInputs.end(); ++it ) {
        (*it)->toInputXML( aOut, aTabs );
    }

    // write the closing tag.
    XMLWriteClosingTag( getXMLReportingName(), aOut, aTabs );
}

//! Output debug info to XML
void NodeInput::toDebugXML( const int aPeriod, ostream& aOut, Tabs* aTabs ) const {
    // write the beginning tag.
    XMLWriteOpeningTag ( getXMLReportingName(), aOut, aTabs, mName );

    XMLWriteElement( mSigmaNewCapital, "Sigma1", aOut, aTabs );
    XMLWriteElement( mSigmaOldCapital, "Sigma2", aOut, aTabs );
    XMLWriteElement( mAlphaCoef, "coefficient", aOut, aTabs );
    XMLWriteElement( mNodeCurrencyDemand, "demandCurrency", aOut, aTabs );
    XMLWriteElement( mPricePaid, "pricePaid", aOut, aTabs );

    //toDebugXMLDerived( aPeriod, aOut, aTabs );
    XMLWriteElement( mProdDmdFnType, "prodDmdFnType", aOut, aTabs );
    for( CNestedInputIterator it = mNestedInputs.begin(); it != mNestedInputs.end(); ++it ) {
        (*it)->toDebugXML( aPeriod, aOut, aTabs );
    }

    // write the closing tag.
    XMLWriteClosingTag( getXMLReportingName(), aOut, aTabs );
}

/*! \brief Get the XML node name for output to XML.
*
* This public function accesses the private constant string, XML_NAME.
* This way the tag is always consistent for both read-in and output and can be easily changed.
* This function may be virtual to be overridden by derived class pointers.
* \author Josh Lurz, James Blackwood
* \return The constant XML_NAME.
*/
const string& NodeInput::getXMLReportingName() const {
    return getXMLNameStatic();
}

/*! \brief Get the XML node name in static form for comparison when parsing XML.
* \details This public function accesses the private constant string, XML_NAME.
*          This way
* the tag is always consistent for both read-in and output and can be easily
* changed. The "==" operator that is used when parsing, required this second
* function to return static.
* \note A function cannot be static and virtual.
* \author Josh Lurz, James Blackwood
* \return The constant XML_NAME as a static.
*/
const string& NodeInput::getXMLNameStatic() {
    const static string XML_NAME = "nodeInput";
    return XML_NAME;
}

//! Get the name of the NodeInput
const string& NodeInput::getName() const {
    return mName;
}

void NodeInput::removeEmptyInputs() {
    const Modeltime* modeltime = scenario->getModeltime();
    const int BASE_PERIOD = modeltime->getBasePeriod();
    double currNodeDemand = 0;

    for( NestedInputIterator nestedInputIter = mNestedInputs.begin(); 
        nestedInputIter != mNestedInputs.end(); )
    {
        (*nestedInputIter)->removeEmptyInputs();
        double currLeafDemand = (*nestedInputIter)->getPhysicalDemand( BASE_PERIOD );
        if( currLeafDemand == 0 ){
            // Remove the empty input object.
            delete *nestedInputIter;
            // Erase the empty input from the vector.
            nestedInputIter = mNestedInputs.erase( nestedInputIter );
        }
        // If demand is valid, check next input.
        else {
            currNodeDemand += currLeafDemand;
            ++nestedInputIter;
        }
    }
    // Set the currency demand for the node so that
    // the node is not removed.  If the node demand is 0, the node
    // will be removed as well as the leaves.
    mNodeCurrencyDemand.init( currNodeDemand );
}

void NodeInput::initialize() {
    double valueSum = 0;
    for( NestedInputIterator it = mNestedInputs.begin(); it != mNestedInputs.end(); ++it ) {
        (*it)->initialize();
        // note, this is really currency during initialization
        valueSum += (*it)->getPhysicalDemand( 0 );
    }
    mNodeCurrencyDemand.init( valueSum );

    // should I create a isRoot()?
    if( mName != "root" ) {
        // note the use of init here implies do not overwrite an already read in value
        mPricePaid.init( 1.0 );
        // TODO: hack for consumers
        mAlphaCoef.init( 1.0 );
    }
}

void NodeInput::calcCoefficient( const std::string& aRegionName, const std::string& aSectorName,
        const int aTechPeriod )
{
    // We need to create a vector which includes this input as the first element since the production
    // function will need access to it.  Really it just needs the totalDemand and output price at this
    // node.
    vector<IInput*> temp;
    temp.push_back( this );
    // have child inputs calculate their children's coefficients first
    for( NestedInputIterator it = mNestedInputs.begin(); it != mNestedInputs.end(); ++it ) {
        (*it)->calcCoefficient( aRegionName, aSectorName, aTechPeriod );
        temp.push_back( *it );
    }

    // now have the production function calculate the coefficient for our direct children
    double retAlpha = mProdDmdFn->calcCoefficient( temp, 0, aRegionName, aSectorName, aTechPeriod, mCurrentSigma, 0, 0 );
    // TODO: check this
    if( mName == "root" ) {
        mAlphaCoef = retAlpha;
    }
}

void NodeInput::changeElasticity( const std::string& aRegionName, const int aPeriod, const double aAlphaZero ) {
    // have children adjust their children's coefficients first
    for( NestedInputIterator it = mNestedInputs.begin(); it != mNestedInputs.end(); ++it ) {
        (*it)->changeElasticity( aRegionName, aPeriod, aAlphaZero );
    }
    // set the current sigma to old capital
    mCurrentSigma = mSigmaOldCapital;
    
    // calculate new alpha coefficients for our direct children
    // TODO: could save some time by not calculating these if old and new sigmas
    // are the same
    mProdDmdFn->changeElasticity( mChildInputsCache, aRegionName, mBasePricePaid, 0, 0, aPeriod, aAlphaZero,
        mSigmaNewCapital, mSigmaOldCapital );
    if( mName == "root" ) {
        // set alpha zero to 1 as we have pulled the difference into the other alphas
        mAlphaCoef = 1.0;
    }
}
void NodeInput::changeSigma( const string& aRegionName, const int aPeriod,
        const double aAlphaZero )
{
    // have children adjust their children's coefficients first
    for( NestedInputIterator it = mNestedInputs.begin(); it != mNestedInputs.end(); ++it ) {
        (*it)->changeSigma( aRegionName, aPeriod, aAlphaZero );
    }
    
    // calculate new alpha coefficients for our direct children
    // Notice we are changing sigmas from the current (which would have been
    // copied forward from the previous year's nesting) to new capital which allows the
    // user to properly change sigmas over time.
    // TODO: could save some time by not calculating these if current and new sigmas
    // are the same
    mProdDmdFn->changeElasticity( mChildInputsCache, aRegionName, mBasePricePaid, 1, 0, aPeriod, aAlphaZero,
        mCurrentSigma, mSigmaNewCapital );

    // set the current sigma to old capital
    mCurrentSigma = mSigmaNewCapital;
    if( mName == "root" ) {
        // set alpha zero to 1 as we have pulled the difference into the other alphas
        mAlphaCoef = 1.0;
    }
}

void NodeInput::calcLevelizedCost( const std::string& aRegionName, const std::string& aSectorName,
        const int aPeriod, const double aAlphaZero ) 
{
    /*!
     * \warning Using the mNodePriceSet hack here to avoid excessive calls to calcLevelizedCost
     *          this means we are relying on the technology to call resetCalcLevelizedCostFlag
     *          at the end of it's operate to ensure that prices will be recalculated the next
     *          time the solver changes prices.  Note calling calcVariableCost will also implicitly
     *          reset the flag.
     */
    if( mNodePriceSet ) {
        return;
    }
    // have children calculate their levelized costs first
    // the leaves are assumed to already have calculated their appropriate price paid
    for( NestedInputIterator it = mNestedInputs.begin(); it != mNestedInputs.end(); ++it ) {
        (*it)->calcLevelizedCost( aRegionName, aSectorName, aPeriod, aAlphaZero );
    }

    // use the function to calculate our levelized costs
    double tempPrice = mProdDmdFn->calcLevelizedCost( mChildInputsCache, aRegionName, aSectorName, aPeriod,
        aAlphaZero, mCurrentSigma );

    setPricePaid( tempPrice, aPeriod );

    // we only set the hack for the root so that we don't have to recurse through
    // the nest when it comes time to reset this flag
    if( mName == "root" ) {
        mNodePriceSet = true;
    }

    // We need to store the base year price paids since they are require to adjust our coefficients
    // with new sigmas in the future.  Note this may be inconsistent if the base year was not read in
    // balanced.
    if ( aPeriod == 0 ){
        mBasePricePaid = tempPrice;
    }
}

/*!
 * \brief Reset the hack.
 * \details Calling this will force the nesting structure to recalculate node
 *          prices next time calcLeveizedCost is called.
 */
void NodeInput::resetCalcLevelizedCostFlag() {
    mNodePriceSet = false;
}

/*
void NodeInput::calcCapitalPrice( const std::string& aRegionName, const std::string& aSectorName,
        const int aPeriod, const double aAlphaZero ) {
    // we must have a price by now to calculate the price for the capital path
    assert( mPricePaid.isInited() );
    assert( hasCapitalInput() );

    //ILogger& sgmLog = ILogger::getLogger( "sgm_debug_log" );
    //sgmLog.setLevel( ILogger::DEBUG );
    // TODO: this math should be in a production function
    double r = 1 - mCurrentSigma;
    double siblingSum = 0.0;
    for( NestedInputIterator it = mNestedInputs.begin(); it != mNestedInputs.end(); ++it ) {
        if( it != mCapitalPath ) {
            siblingSum += pow( (*it)->getPricePaid( aRegionName, aPeriod ) / 
                (*it)->getCoefficient( aPeriod ), r );
        }
    }
    double capitalPathPricePaid = pow( getPricePaid( aRegionName, aPeriod ) * aAlphaZero, r );
    capitalPathPricePaid -= siblingSum;
    if( capitalPathPricePaid < 0 ) {
        //cout << "went < 0 here " << endl;
        capitalPathPricePaid = 0;
    }
    capitalPathPricePaid = pow( capitalPathPricePaid, 1 / r ) * (*mCapitalPath)->getCoefficient( aPeriod );
    //sgmLog << "calc cap  p: " << capitalPathPricePaid << endl;
    (*mCapitalPath)->setPricePaid( capitalPathPricePaid, aPeriod );
    (*mCapitalPath)->calcCapitalPrice( aRegionName, aSectorName, aPeriod, aAlphaZero );
}
*/

double NodeInput::calcInputDemand( const std::string& aRegionName, const std::string& aSectorName,
        const int aPeriod, const double aPhysicalOutput, const double aUtilityParameterA,
        const double aAlphaZero )
{
    // first calculate the demands for the direct children
    double retDemand = mProdDmdFn->calcDemand( mChildInputsCache, aPhysicalOutput, aRegionName, aSectorName, 1,
        aPeriod, aUtilityParameterA, aAlphaZero, mCurrentSigma, mPricePaid );

    // have all of the children calculate demands for their children
    for( NestedInputIterator it = mNestedInputs.begin(); it != mNestedInputs.end(); ++it ) {
        // note we are passing the just calculated demand for the child to it's calcInputDemand 
        // so that we can ensure that the correct amout of demand will be distrubuted by that node
        // rather than relying on having the input output ratio from that node to the root
        (*it)->calcInputDemand( aRegionName, aSectorName, aPeriod, (*it)->getPhysicalDemand( aPeriod ),
            aUtilityParameterA, aAlphaZero );
    }
    return retDemand;
}

double NodeInput::calcCapitalOutputRatio( const std::string& aRegionName, const std::string& aSectorName,
        const int aPeriod, const double aAlphaZero ) {
    /*
    assert( hasCapitalInput() );
    // this is a hack because the leaf does not have a prodDmdFn so the node will
    // calculate the io coef for it's child and the root will do itself and the child
    vector<IInput*> temp;
    temp.push_back( *mCapitalPath );
    double childIOCoef = mProdDmdFn->getCapitalOutputRatio( temp, aRegionName, aSectorName, mPricePaid, aPeriod,
        aAlphaZero, mCurrentSigma );
    return childIOCoef * (*mCapitalPath)->calcCapitalOutputRatio( aRegionName, aSectorName, aPeriod, aAlphaZero );
    */
    // TODO: might be able to get this to work
    return 0;
}

void NodeInput::calcVariableLevelizedCost(const std::string& aRegionName, const std::string& aSectorName,
        const int aPeriod, const double aAlphaZero ) 
{
    // have children calculate their levelized costs first
    vector<IInput*> temp;
    for( NestedInputIterator it = mNestedInputs.begin(); it != mNestedInputs.end(); ++it ) {
        // variable means don't include capital
        if( !(*it)->hasTypeFlag( IInput::CAPITAL ) ) {
            (*it)->calcVariableLevelizedCost( aRegionName, aSectorName, aPeriod, aAlphaZero );
            temp.push_back( *it );
        }
    }

    // now calculate the cost for this node, note that any capital nodes will be omitted from
    // this calculation
    double tempPrice = mProdDmdFn->calcLevelizedCost( temp, aRegionName, aSectorName, aPeriod,
        aAlphaZero, mCurrentSigma );
    
    setPricePaid( tempPrice, aPeriod );

    // since we just overwrote all of the node prices with varible costs that implies we
    // had better recalculate the node prices the next time calcLevelizedCost is called
    resetCalcLevelizedCostFlag();
}

void NodeInput::applyTechnicalChange( const std::string& aRegionName, const std::string& aSectorName,
        const int aPeriod, const TechChange& aTechChange )
{
    // have children adjust their children's coefficients for tech chagne first
    for( NestedInputIterator it = mNestedInputs.begin(); it != mNestedInputs.end(); ++it ) {
        (*it)->applyTechnicalChange( aRegionName, aSectorName, aPeriod, aTechChange );
    }

    // use the production function to apply tech change to direct children now
    // TODO: I don't think alpha zero and sigma are necessary but I could pass them in I guess
    mProdDmdFn->applyTechnicalChange( mChildInputsCache, aTechChange, aRegionName, aSectorName, aPeriod, 1, 1 );
}

const IFunction* NodeInput::getFunction() const {
    assert( mProdDmdFn );

    return mProdDmdFn;
}

double NodeInput::getLevelizedCost( const std::string& aRegionName, const std::string& aSectorName,
        const int aPeriod ) const
{
    // this is not the best solution but since it must be const I have to calc seperatly
    // before you can get the levelized cost
    return mPricePaid;
}

double NodeInput::getCurrencyDemand( const int aPeriod ) const {
    assert( mNodeCurrencyDemand.isInited() );

    return mNodeCurrencyDemand;
}

void NodeInput::setCurrencyDemand( const double aCurrencyDemand,
                                    const std::string& aRegionName, 
                                    const int aPeriod ) {
    // setting a currency demand for a node makes no sense
    assert( false );
}

double NodeInput::getPhysicalDemand( const int aPeriod ) const {
    // does not make sense for nodes but for the sake of consistency
    // for removeEmptyInputs and initialize return the currency amout
    return getCurrencyDemand( aPeriod );
}
    
double NodeInput::getCarbonContent( const int aPeriod ) const {
    // does not make sense for nodes
    assert( false );
    return 0;
}
    
void NodeInput::setPhysicalDemand( const double aPhysicalDemand,
                                    const std::string& aRegionName, 
                                    const int aPeriod )
{
    // allowing this now so that calcInputDemand can be
    // a top down operation.
    // TODO: should make this mNodeValue
    mNodeCurrencyDemand.set( aPhysicalDemand );
}

double NodeInput::getPrice( const std::string& aRegionName,
                             const int aPeriod ) const
{
    return mPricePaid;
}

void NodeInput::setPrice( const std::string& aRegionName,
                           const double aPrice,
                           const int aPeriod )
{
    mPricePaid = aPrice;
}

double NodeInput::getPriceAdjustment() const {
    // TODO:
    return 0;
}

double NodeInput::getPricePaid( const std::string& aRegionName,
                                 const int aPeriod ) const
{
    return mPricePaid;
}

void NodeInput::setPricePaid( const double aPricePaid,
                               const int aPeriod )
{
    mPricePaid = aPricePaid;
}

double NodeInput::getCoefficient( const int aPeriod ) const {
    assert( mAlphaCoef.isInited() );

    return mAlphaCoef;
}

void NodeInput::setCoefficient( const double aCoefficient,
                                 const int aPeriod )
{
    mAlphaCoef = aCoefficient;
}

double NodeInput::getConversionFactor( const int aPeriod ) const {
    // TODO:
    return 0;
}

double NodeInput::getCO2EmissionsCoefficient( const std::string& aGHGName,
                                             const int aPeriod ) const
{
    // TODO:
    return 0;
}

void NodeInput::tabulateFixedQuantity( const std::string& aRegionName,
                                        const double aFixedOutput,
                                        const bool aIsInvestmentPeriod,
                                        const int aPeriod )
{
    // TODO:
}

void NodeInput::scaleCalibrationQuantity( const double aScaleFactor ) {
    // TODO:
}

double NodeInput::getCalibrationQuantity( const int aPeriod ) const {
    // TODO:
    return 0;
}

double NodeInput::getPriceElasticity() const {
    // TODO: should not be here
    return mAlphaUtilityParam;
}

double NodeInput::getIncomeElasticity() const {
    // TODO: should not be here
    return mBetaUtilityParam;
}

double NodeInput::getTechChange( const int aPeriod ) const {
    return mTechChange;
}

void NodeInput::accept( IVisitor* aVisitor, const int aPeriod ) const {
    aVisitor->startVisitNodeInput( this, aPeriod );

    // visit the nested inputs
    for( CNestedInputIterator it = mNestedInputs.begin(); it != mNestedInputs.end(); ++it ) {
        (*it)->accept( aVisitor, aPeriod );
    }

    aVisitor->endVisitNodeInput( this, aPeriod );
}
