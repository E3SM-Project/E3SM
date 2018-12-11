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
* \file sgm_input.cpp
* \ingroup Objects
* \brief The SGMInput class source file.
* \author Pralit Patel
* \author Sonny Kim
*/

#include "util/base/include/definitions.h"
#include <iostream>
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#include <cmath>

#include "functions/include/sgm_input.h"
#include "containers/include/scenario.h"
#include "marketplace/include/marketplace.h"
#include "util/base/include/xml_helper.h"
#include "util/base/include/configuration.h"
#include "functions/include/function_utils.h"
#include "containers/include/iinfo.h"
#include "util/base/include/ivisitor.h"
#include "emissions/include/aghg.h"
#include "sectors/include/more_sector_info.h"
#include "containers/include/national_account.h"
#include "technologies/include/expenditure.h"
#include "marketplace/include/cached_market.h"

using namespace std;
using namespace xercesc;

extern Scenario* scenario;

// static initialize.
const string SGMInput::XML_REPORTING_NAME = "input-SGM";

//! Default Constructor
SGMInput::SGMInput(): mTypeFlags( 0 ){
}

/*! Copy Constructor
* We need to define a copy constructor so that values which are stored
* for reporting purposes such as currency demand do not get copied to
* a new vintage since it really didn't have that currency demand.
* \param aInput The input to copy.
* \author Pralit Patel
*/
SGMInput::SGMInput( const SGMInput& aInput ):mName( aInput.mName),
mPriceAdjustFactor( aInput.mPriceAdjustFactor ),
mTechnicalChange( aInput.mTechnicalChange ),
mConversionFactor( aInput.mConversionFactor ),
mCO2Coefficient( aInput.mCO2Coefficient ),
mTypeFlags( aInput.mTypeFlags ),
mSalesTaxRate( aInput.mSalesTaxRate )
{
    // we do need to pass forward the coeff however we don't
    // know which and to where since we don't have period info
    // so we will do the last period that had an inited value
    // and pass it into the next period
    int currPeriod = aInput.mCoefficient.size() - 1;
    for( ; currPeriod >= 0 
            && !aInput.mCoefficient[ currPeriod ].isInited(); --currPeriod ) {
        // nothing to do in the loop
    }
    if( currPeriod >= 0 ) {
        // if there was a valid coef to copy put it in as at currPer + 1
        // which is the next period after the last period that had a valid coef
        // note that trying to copy the last period forward is undefined and
        // will likely crash here
        mCoefficient[ currPeriod + 1] = aInput.mCoefficient[ currPeriod ];
        
        // this is a hack because we may need to copy the coeff back to
        // vintages before period 0 note that this will also leave period
        // 1 vintages with a coef in period 0 which is not entirely correct
        if( currPeriod == 0 ) {
            mCoefficient[ 0 ] = aInput.mCoefficient[ 0 ];
        }
        
        // always copy the base year quantities and prices
        mPhysicalDemand[ 0 ] = aInput.mPhysicalDemand[ 0 ];
        mPricePaid[ 0 ] = aInput.mPricePaid[ 0 ];
    }
}

//! Destructor
SGMInput::~SGMInput() {}

/*! \brief Get the XML name for reporting to XML file.
*
* This public function accesses the private constant string, XML_NAME. This way
* the tag is always consistent for reporting outputs and can be easily
* changed.
* \author Sonny Kim
* \return The constant XML_NAME.
*/
const string& SGMInput::getXMLReportingName() const{
    return XML_REPORTING_NAME;
}

//! Parse the SGMInput's XML.
void SGMInput::XMLParse( const xercesc::DOMNode* node ) {
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

        if ( nodeName == "coefficient" ) {
            // TODO: assuming period zero here
            mCoefficient[ 0 ].init( XMLHelper<double>::getValue( curr ) );
        }
        else if( nodeName == "demandCurrency" ) {
            // TODO: assuming period zero here if no year attr exists,
            //       it would be better to just force the year attr
            /*!
             * \warning Placing initial currency demand in physical demand,
             *          this means coef calculation routines will need to use
             *          getPhysicalDemand.
             */
            if( XMLHelper<int>::getAttr( curr, "year" ) == 0 ) {
                mPhysicalDemand[ 0 ].init( XMLHelper<double>::getValue( curr ) );
            }
            else {
                XMLHelper<Value>::insertValueIntoVector( curr, mPhysicalDemand, scenario->getModeltime() );
            }
        }
        else if( nodeName == "priceAdjustFactor" ) {
            mPriceAdjustFactor.init( XMLHelper<double>::getValue( curr ) );
        }
        else if( nodeName == "technicalChange" ) {
            mTechnicalChange.init( XMLHelper<double>::getValue( curr ) );
        }
        else if( nodeName == "flag" ) {
            setFlagsByName( XMLHelper<string>::getAttr( curr, "type" ) );
        }
        else if( nodeName == "sales-tax-rate" ) {
            mSalesTaxRate.init( XMLHelper<double>::getValue( curr ) );
        }
        else if( !XMLDerivedClassParse( nodeName, curr ) ){
            cout << "Unrecognized text string: " << nodeName << " found while parsing " << getXMLName() << "." << endl;
        }
    }
}

//! Complete the initialization of the SGMInput.
void SGMInput::completeInit( const string& aRegionName,
                             const string& aSectorName,
                             const string& aSubsectorName,
                             const string& aTechName,
                             DependencyFinder* aDependencyFinder , const IInfo* aTechInfo)
{
}

void SGMInput::initCalc( const string& aRegionName,
                         const string& aSectorName,
                         const bool aIsNewInvestmentPeriod,
                         const bool aIsTrade,
                         const int aPeriod )
{
    /*! \pre There must be a valid region name. */
    assert( !aRegionName.empty() );

    mCachedMarket = scenario->getMarketplace()->locateMarket( mName, aRegionName, aPeriod );

    initializeCachedCoefficients( aRegionName );
    initializeTypeFlags( aRegionName );
    
    // copy forward the coeff only if it has not already been inited
    if( aPeriod > 0 && !mCoefficient[ aPeriod ].isInited() ) {
        mCoefficient[ aPeriod ].init( getCoefficient( aPeriod - 1 ) );
    }
    
    // trade does not use coeffs but to avoid assertions init them
    if( aPeriod == 0 && aIsTrade ) {
        mCoefficient[ aPeriod ].init( 0 );
    }

}

void SGMInput::copyParam( const IInput* aInput,
                          const int aPeriod )
{
    // Name must be already defined.
    // ensure that inputs before the base period just copy from the base period
    int prevPeriod = aPeriod > 0 ? aPeriod - 1 : 0;
    mCoefficient[ aPeriod ].set( aInput->getCoefficient( prevPeriod ) );
    
    // Parameters which should not get copied forward if they didn't exist.
    // using the init method of Value will ensure this
    mPriceAdjustFactor.init( aInput->getPriceAdjustment() );
    mTechnicalChange.init( aInput->getTechChange( prevPeriod ) );

    // Always copy the base year price paid and quantity
    mPricePaid[ 0 ].set( aInput->getPricePaid( "", 0 ) );
    mPhysicalDemand[ 0 ].set( aInput->getPhysicalDemand( 0 ) );
    
    // warning: had to copy sales tax rate in production input
}

bool SGMInput::isSameType( const string& aType ) const {
    return aType == getXMLName();
}

//! Output to XML data
void SGMInput::toInputXML( ostream& out, Tabs* tabs ) const {
    // write the beginning tag.
    XMLWriteOpeningTag ( getXMLName(), out, tabs, mName );

    XMLWriteElement( mCoefficient[ 0 ], "coefficient", out, tabs );
    XMLWriteElement( mPhysicalDemand[ 0 ], "demandCurrency", out, tabs );
    XMLWriteElement( mConversionFactor, "conversionFactor", out, tabs );
    XMLWriteElement( mPriceAdjustFactor, "priceAdjustFactor", out, tabs );
    XMLWriteElement( mSalesTaxRate, "sales-tax-rate", out, tabs );

    toInputXMLDerived( out, tabs );

    // write the closing tag.
    XMLWriteClosingTag( getXMLName(), out, tabs );
}

//! Output debug info to XML
void SGMInput::toDebugXML( const int period, ostream& out, Tabs* tabs ) const {
    // write the beginning tag.
    XMLWriteOpeningTag ( getXMLName(), out, tabs, mName );

    XMLWriteElement( mCoefficient[ period ], "coefficient", out, tabs );
    XMLWriteElement( mPhysicalDemand[ period ], "demandPhysical", out, tabs );
    XMLWriteElement( mPricePaid[ period ], "pricePaid", out, tabs );
    XMLWriteElement( mConversionFactor, "conversionFactor", out, tabs );
    XMLWriteElement( mPriceAdjustFactor, "priceAdjustFactor", out, tabs );
    XMLWriteElement( mSalesTaxRate, "sales-tax-rate", out, tabs );

    toDebugXMLDerived( period, out, tabs );

    // write the closing tag.
    XMLWriteClosingTag( getXMLName(), out, tabs );
}

//! Get the name of the SGMInput
const string& SGMInput::getName() const {
    return mName;
}

/*! \brief Get the currency to physical conversion factor for the SGMInput.
* \details This function returns the conversion factor used to convert the SGMInput
*          from a currency to a physical demand.
* \param aPeriod Model period.
* \return The currency to physical conversion factor.
* \warning SGM now relies on prices to convert from currency to physical and this
*          method should likely be removed or always return one.
* \author Josh Lurz
*/
double SGMInput::getConversionFactor( const int aPeriod ) const {
    // TODO: what to do about this as it no longer makes sense
    /*! \pre The conversion factor is initialized. */
    assert( mConversionFactor.isInited() );
    return mConversionFactor;
}

/*! \brief Get the emissions gas coefficient for an SGMInput.
* \details This function returns the emissions gas coefficent for an SGMInput.
* \param aGHGName Name of the GHG to return the coefficient for.
* \param aPeriod Model period.
* \return The emissions coefficient.
* \todo This doesn't work for multiple gases.
* \author Josh Lurz
*/
double SGMInput::getCO2EmissionsCoefficient( const string& aGHGName, const int aPeriod ) const {
    // Currently this assumes the only gas is CO2.
    assert( aGHGName == "CO2" );
    
    /*! \pre The CO2 coefficient is initialized. */
    assert( mCO2Coefficient.isInited() );
    return mCO2Coefficient;
}

/*! \brief Get the currency demand for the SGMInput.
* \details Return the value of demand for this SGMInput in currency units. This
*          is calculated by multiplying the physical demand by the price
*          paid for the SGMInput.
* \param aPeriod Model period.
* \return The demand for this SGMInput in currency units.
* \warning The price paid must be calculated before using this method.
* \author Josh Lurz
*/
double SGMInput::getCurrencyDemand( const int aPeriod ) const {
    // TODO: the region name isn't necessary for getPricePaid
    // however using "" doesn't seem like a good idea
    // TODO: check the price paid for capital it may well be the price of the capital good
    const double pricePaid = !hasTypeFlag( IInput::CAPITAL ) ? getPricePaid( "", aPeriod ) : 1;
    return getPhysicalDemand( aPeriod ) * pricePaid;
}

/*! \brief Get the carbon content for the SGMInput.
* \param aPeriod Model period.
* \return Carbon content of input.
* \author Sonny Kim
*/
double SGMInput::getCarbonContent( const int aPeriod ) const {
    return 0;
}

/*! \brief Get the Physical Demand.
* \param aPeriod Model period.
* \return Physical demand.
*/
double SGMInput::getPhysicalDemand( const int aPeriod ) const {
    // can't check this due to remove empty inputs
    //assert( mPhysicalDemand[ aPeriod ].isInited() );
    
    return mPhysicalDemand[ aPeriod ];
}

/*!
 * \brief Get the input specific technical change.
 * \param aPeriod Model period.
 * \return The input specific technical change.
 * \author Josh Lurz
 */
double SGMInput::getTechChange( const int aPeriod ) const {
    return mTechnicalChange;
}

//! Set Physical Demand.
void SGMInput::setPhysicalDemand( double aPhysicalDemand, const string& aRegionName, const int aPeriod )
{
    mPhysicalDemand[ aPeriod ].set( aPhysicalDemand );
    
    // never add capital demand directly to the marketplace since
    // this physical demand is a quantity of capital and the Capital market
    // works in annual dollar amounts
    if( !hasTypeFlag( IInput::CAPITAL ) && aPhysicalDemand > util::getSmallNumber() ) {
        mCachedMarket->addToDemand( mName, aRegionName, aPhysicalDemand, aPeriod );
    }
}

//! Set currency Demand.
void SGMInput::setCurrencyDemand( double aCurrencyDemand, const string& aRegionName,
                               const int aPeriod )
{
    // SGM inputs cannot directly set currency demands.
    // could possibly get this to work
    assert( false );
}

/*! \brief Get the SGMInput-output coefficient.
* \param aPeriod Model period.
* \return The IO coefficient.
*/
double SGMInput::getCoefficient( const int aPeriod ) const {
    // TODO: it would be nice to assert this however it will cause the 
    // xmldb outputter to crash since it checks all periods
    //assert( mCoefficient[ aPeriod ].isInited() );
    
    return mCoefficient[ aPeriod ];
}

/*! \brief Set the IO coefficient.
 * \author Pralit Patel
 * \param aCoefficient new coefficient value
 * \param aPeriod Model period.
 */
void SGMInput::setCoefficient( const double aCoefficient, const int aPeriod ) {
    assert( aCoefficient != 0 ); // Can't set coefficients to zero.
    mCoefficient[ aPeriod ].set( aCoefficient );
}

/*! \brief Return the market price, or unadjusted price, for the SGMInput.
* \param aRegionName Region containing the SGMInput.
* \param aPeriod Period to find the price in.
* \return The market or unadjusted price.
* \author Josh Lurz
*/
double SGMInput::getPrice( const string& aRegionName, const int aPeriod ) const {
    return mCachedMarket->getPrice( mName, aRegionName, aPeriod );
}

void SGMInput::setPrice( const string& aRegionName,
                         const double aPrice,
                         const int aPeriod )
{
    // SGM markets are solved so the input price cannot be set.
}

/*! \brief Returns the price paid for each SGMInput.
* \param aRegionName Name of the containing region.-
* \param aPeriod Model period.
* \author Sonny Kim
*/
double SGMInput::getPricePaid( const string& aRegionName, const int aPeriod ) const{
    // TODO: it would be nice to assert this however it will cause the 
    // xmldb outputter to crash since it checks all periods
    //assert( mPricePaid[ aPeriod ].isInited() );
    return mPricePaid[ aPeriod ];
}
/*! \brief Set the price paid for each SGMInput.
*
* \param aPricePaid new price paid value
* \param aPeriod Model period.
* \author Sonny Kim
*/
void SGMInput::setPricePaid( double aPricePaid, const int aPeriod ) {
    mPricePaid[ aPeriod ].set( aPricePaid );
}

void SGMInput::calcPricePaid( const string& aRegionName, const string& aSectorName, const MoreSectorInfo* aMoreSectorInfo,
                              const vector<AGHG*>& aGhgs, const ICaptureComponent* aSequestrationDevice,
                              const int aLifetimeYears, const int aPeriod )
{
    // initialize so that aMoreSectorInfo does not have any impact on price 
    // and override only if aMoreSectorInfo is not null
    double transportationAdder = 0;
    double transportationMult = 1;
    double proportionalTax = 1;
    double additiveTax = 0;

    // if pointer is not null
    if( aMoreSectorInfo ) {
        transportationAdder = aMoreSectorInfo->getValue( MoreSectorInfo::TRANSPORTATION_COST );
        transportationMult = aMoreSectorInfo->getValue( MoreSectorInfo::TRAN_COST_MULT );
        proportionalTax = aMoreSectorInfo->getValue( MoreSectorInfo::PROPORTIONAL_TAX_RATE );
        additiveTax = aMoreSectorInfo->getValue( MoreSectorInfo::ADDITIVE_TAX );
    }
    double tempPricePaid = 0;
    if( hasTypeFlag( IInput::CAPITAL ) ){
        double intrestRate = getPrice( aRegionName, aPeriod ) + getPriceAdjustment();
        double capitalGoodPrice = FunctionUtils::getCapitalGoodPrice( aRegionName, aPeriod );
        double depreciatedIntrestRate = intrestRate / 
            ( 1 - pow( 1 / ( 1 + intrestRate), aLifetimeYears ) );
        tempPricePaid = depreciatedIntrestRate * capitalGoodPrice * ( 1 + mSalesTaxRate );
    }
    else {
        // Calculate GHG taxes.
        double carbonTax = 0;
        for( vector<AGHG*>::const_iterator ghg = aGhgs.begin(); ghg != aGhgs.end(); ++ghg ){
            carbonTax += (*ghg)->getGHGValue( this, aRegionName, aSectorName,
                aSequestrationDevice, aPeriod );
        }
        tempPricePaid = ( ( getPrice( aRegionName, aPeriod ) * ( 1 + mSalesTaxRate ) + 
            ( transportationAdder * transportationMult ) ) * proportionalTax + additiveTax + carbonTax ) * getPriceAdjustment();
    }
    setPricePaid( tempPricePaid, aPeriod );
}

double SGMInput::calcTaxes( const string& aRegionName, NationalAccount* aNationalAccount, 
                         Expenditure* aExpenditure, const int aPeriod ) const
{
    // TODO: I am not sure if this is properly taking into accoutn transportation adder/mult or price
    // adjustments
    double salesTax;
    if( !hasTypeFlag( IInput::CAPITAL ) ) {
        salesTax = mSalesTaxRate * getPrice( aRegionName, aPeriod ) * getPhysicalDemand( aPeriod );
    }
    else {
        // since we do not have the depreciation we have to back out the correct price
        salesTax = ( getPricePaid( aRegionName, aPeriod ) / ( 1 + mSalesTaxRate ) ) * mSalesTaxRate * getPhysicalDemand( aPeriod );
    }
    // TODO: factor inputs could possibly go into another account but would have to create
    // that accounting category
    if( aNationalAccount && aExpenditure ) {
        if( !hasTypeFlag( IInput::LABOR ) ) {
            aNationalAccount->addToAccount( NationalAccount::INDIRECT_BUSINESS_TAX, salesTax );
            aExpenditure->addToType( Expenditure::INDIRECT_TAXES, salesTax );
        } else {
            aNationalAccount->addToAccount( NationalAccount::SOCIAL_SECURITY_TAX, salesTax );
            aExpenditure->addToType( Expenditure::SOCIAL_SECURITY_TAX, salesTax );
        }
    }
    return salesTax;
}

/*! \brief Returns the price received for an SGMInput.
* \details Queries the marketplace and get the price received from the market info.
* \param aRegionName Name of the region containing the SGMInput.
* \param aPeriod Period
* \return The price received for the SGMInput.
* \author Josh Lurz
*/
double SGMInput::getPriceReceived( const string& aRegionName, const int aPeriod ) const {
    return FunctionUtils::getPriceReceived( aRegionName, mName, aPeriod );
}

/*! \brief Returns the price adjustment.
* \note For all inputs other than capital this is a multiplier, for capital it is an adder.
* \author Sonny Kim
*/
double SGMInput::getPriceAdjustment() const {
    return mPriceAdjustFactor;
}

bool SGMInput::hasTypeFlag( const int aTypeFlag ) const {

    /*! \pre The type flags must be initialized. */
    //assert( mTypeFlags != INITIALIZED );
    //assert( mTypeFlags != 0 );

    return ( ( aTypeFlag & ~mTypeFlags ) == 0 );
}

double SGMInput::getCalibrationQuantity( const int aPeriod ) const
{
    return -1;
}

/*! \brief Initialize the cached emissions coefficient and physical to energy conversion factor.
* \param aRegionName Region name.
*/
void SGMInput::initializeCachedCoefficients( const string& aRegionName ){
    // If the conversion factor has not been initialized, store the value from
    // the marketplace internally.
    if( !mConversionFactor.isInited() ){
        // Get the conversion factor from the marketplace.
        const IInfo* marketInfo = scenario->getMarketplace()->getMarketInfo( mName, aRegionName, 0, false );
        const double convFactor = marketInfo ? marketInfo->getDouble( "ConversionFactor", false ) : 0;
        if( convFactor == 0 && isEnergyGood( aRegionName ) ){
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::WARNING );
            mainLog << "Conversion factor of zero for energy SGMInput " << mName << "." << endl;
        }
        mConversionFactor.init( convFactor );
    }
    
    // If the coefficient has not been initialized, store the value from the
    // marketplace internally.
    if( !mCO2Coefficient.isInited() ){
        mCO2Coefficient.init( FunctionUtils::getCO2Coef( aRegionName, mName, 0 ) );
    }
}

/*! \brief Initialize the type flags.
* \details Sets any additional type flags based on the input name so
*          that we can categorize inputs.
* \see setFlagsByName
* \see isEnergyGood
* \see isPrimaryEnergyGood
* \see isSecondaryEnergyGood
* \param aRegionName Name of the region.
*/
void SGMInput::initializeTypeFlags( const string& aRegionName ) {

    // Set that the type flags are initialized.
    mTypeFlags |= IInput::INITIALIZED;

    // TODO: it may be a good idea to no longer rely on input names
    // and rather just read in flags by reading a flag XML tag.
    setFlagsByName( mName );

    if( !hasTypeFlag( IInput::FACTOR ) ){
        if( isEnergyGood( aRegionName ) ){
            mTypeFlags |= IInput::ENERGY;
        }
        else {
            mTypeFlags |= IInput::MATERIAL;
        }
    }
    // Initialize primary and secondary energy flags.
    if( isPrimaryEnergyGood( aRegionName ) ){
        mTypeFlags |= IInput::PRIMARY;
    }
    else if( isSecondaryEnergyGood( aRegionName ) ){
        mTypeFlags |= IInput::SECONDARY;
    }
}

/*! \brief Add the correct flags to the mTypeFlags for the given type name.
 * \details This allows us to flexibly initialize the type flags based on
 *          input names or read in flags.
 * \param aTypeName The name to use when determining the type flags.
 * \todo Expand for possible read in flags.
 */
void SGMInput::setFlagsByName( const string& aTypeName ) {
    // Initialize the type.
    // TODO: just read in the flag rather than hard coding the input names here
    if( aTypeName == "Land" || aTypeName == "NatRes" || aTypeName == "ColRes" || aTypeName == "OilRes" || aTypeName == "GasRes" ){
        mTypeFlags |= IInput::LAND;
        mTypeFlags |= IInput::FACTOR;
    }
    else if( aTypeName == "Labor" || aTypeName == "UnSkLab" || aTypeName == "SkLab" ){
        mTypeFlags |=  IInput::LABOR;
        mTypeFlags |=  IInput::FACTOR;
    }
    else if( aTypeName == "Capital" ){
        mTypeFlags |= IInput::CAPITAL;
        mTypeFlags |= IInput::FACTOR;
    }

    // Inititialize the numeraire flag.
    const static string numeraireInputName = Configuration::getInstance()->getString( "numeraire-good", "SVS" );
    if( aTypeName == numeraireInputName ){
        mTypeFlags |= IInput::NUMERAIRE;
    }
}

/*! \brief Return whether a good is an energy good.
* \param aGoodName Good name.
* \return Whether the good is an energy price good.
*/
bool SGMInput::isEnergyGood( const string& aRegionName ) const
{
    const IInfo* marketInfo = scenario->getMarketplace()->getMarketInfo( mName, aRegionName,
                                                                         0, false );
    return marketInfo && marketInfo->getBoolean( "IsEnergyGood", false );
}

/*! \brief Return whether a good is a primary energy good.
* \param aRegionName Region name.
* \param aGoodName Good name.
* \return Whether the good is a primary energy price good.
*/
bool SGMInput::isPrimaryEnergyGood( const string& aRegionName ) const
{
    const IInfo* marketInfo = scenario->getMarketplace()->getMarketInfo( mName, aRegionName,
                                                                         0, false );
    return marketInfo && marketInfo->getBoolean( "IsPrimaryEnergyGood", false );
}

/*! \brief Return whether a good is a secondary energy good.
* \param aRegionName Region name.
* \return Whether the good is a secondary energy price good.
*/
bool SGMInput::isSecondaryEnergyGood( const string& aRegionName ) const
{
    const IInfo* marketInfo = scenario->getMarketplace()->getMarketInfo( mName, aRegionName,
                                                                         0, false );
    return marketInfo && marketInfo->getBoolean( "IsSecondaryEnergyGood", false );
}

/*! \brief For outputing SGM data to a flat csv File
 * 
 * \author Pralit Patel
 * \param period The period which we are outputing for
 */
void SGMInput::csvSGMOutputFile( ostream& aFile, const int period ) const {
    aFile << mName << ',';
    aFile.precision(0);
    aFile << mPhysicalDemand[ period ] << ',';
    aFile.precision(3);
    aFile << mPricePaid[ period ] << endl;
}

void SGMInput::doInterpolations( const int aYear, const int aPerviousYear,
                                 const int aNextYear, const IInput* aPreviousInput,
                                 const IInput* aNextInput )
{
    // TODO: have not even considered what this will mean for SGM
}

void SGMInput::accept( IVisitor* aVisitor, const int aPeriod ) const {
    aVisitor->startVisitSGMInput( this, aPeriod );
    aVisitor->endVisitSGMInput( this, aPeriod );
}
