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
 * \file cached_market.cpp
 * \ingroup Objects
 * \brief CachedMarket class source file.
 * \author Pralit Patel
 */

#include "util/base/include/definitions.h"
#include "marketplace/include/cached_market.h"

#include "marketplace/include/market.h"
#include "util/logger/include/ilogger.h"
#include "util/base/include/util.h"
#include "marketplace/include/marketplace.h"

using namespace std;

/*!
 * \brief Constructor which takes the parameters used to locate the given market.
 * \param aGoodName The good name used to locate aLocatedMarket.  Stored for debugging.
 * \param aGoodName The region name used to locate aLocatedMarket.  Stored for debugging.
 * \param aGoodName The period used to locate aLocatedMarket.  Stored for debugging.
 * \param aLocatedMarket A pointer to the actual market which was located.  Note that this
 *                       parameter can be null which indicates the market was not found.
 */
CachedMarket::CachedMarket( const string& aGoodName, const string& aRegionName, const int aPeriod,
                            Market* aLocatedMarket )
:
#if(!NDEBUG)
mGoodName( aGoodName ),
mRegionName( aRegionName ),
mPeriod( aPeriod ),
#endif
mCachedMarket( aLocatedMarket )
{
}

/*!
 * \brief Destructor
 */
CachedMarket::~CachedMarket() {
}

/*!
 * \brief Set the market price.
 * \details Mimics the behavior of Marketplace::setPrice.
 * \param aGoodName The good of the market.
 * \param aRegionName The region setting the price.
 * \param aValue The value to which to set price.
 * \param aPeriod The period in which to set the price.
 * \param aMustExist Whether it is an error for the market not to exist.
 * \see Marketplace::setPrice
 */
void CachedMarket::setPrice( const string& aGoodName, const string& aRegionName, const double aValue,
                             const int aPeriod, bool aMustExist )
{
    /*!
     * \invariant The given good name matches the one that was used when locating this market.
     */
    assert( aGoodName == mGoodName );
    
    /*!
     * \invariant The given region name matches the one that was used when locating this market.
     */
    assert( aRegionName == mRegionName );
    
    /*!
     * \invariant The given period matches the one that was used when locating this market.
     */
    assert( aPeriod == mPeriod );
    
    // Print a warning message if the new price is not a finite number.
    if ( !util::isValidNumber( aValue ) ) {
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::NOTICE );
        mainLog << "Error setting price in markeplace for: " << aGoodName << ", value: " << aValue << endl;
        return;
    }
    
    if ( mCachedMarket ) {
        mCachedMarket->setPrice( aValue );
    }
    else if( aMustExist ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::NOTICE );
        mainLog << "Cannot set price for market as it does not exist: " << aGoodName << " " 
            << aRegionName << endl;
    }
}

/*!
 * \brief Add to the supply for this market.
 * \details Mimics the behavior of Marketplace::addToSupply.
 * \param aGoodName The good of the market.
 * \param aRegionName The region adding supply.
 * \param aValue The supply value to add.
 * \param aPeriod The period in which to add the supply.
 * \param aMustExist Whether it is an error for the market not to exist.
 * \see Marketplace::addToSupply
 */
void CachedMarket::addToSupply( const string& aGoodName, const string& aRegionName, const double aValue,
                                const int aPeriod, bool aMustExist )
{
    /*!
     * \invariant The given good name matches the one that was used when locating this market.
     */
    assert( aGoodName == mGoodName );
    
    /*!
     * \invariant The given region name matches the one that was used when locating this market.
     */
    assert( aRegionName == mRegionName );
    
    /*!
     * \invariant The given period matches the one that was used when locating this market.
     */
    assert( aPeriod == mPeriod );
    
    // Print a warning message when adding infinity values to the supply.
    if ( !util::isValidNumber( aValue ) ) {
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::NOTICE );
        mainLog << "Error adding to supply in markeplace for: " << aGoodName << ", region: " << aRegionName << ", value: " << aValue << endl;
        return;
    }
    
    if ( mCachedMarket ) {
        mCachedMarket->addToSupply( aValue );
    }
    else if( aMustExist ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::NOTICE );
        mainLog << "Cannot add to supply for market as it does not exist: " << aGoodName << " " 
            << aRegionName << endl;
    }
}

/*!
 * \brief Add to the demand for this market.
 * \details Mimics the behavior of Marketplace::addToDemand.
 * \param aGoodName The good of the market.
 * \param aRegionName The region adding demand.
 * \param aValue The demand value to add.
 * \param aPeriod The period in which to add the demand.
 * \param aMustExist Whether it is an error for the market not to exist.
 * \see Marketplace::addToDemand
 */
void CachedMarket::addToDemand( const string& aGoodName, const string& aRegionName, const double aValue,
                                const int aPeriod, bool aMustExist )
{
    /*!
     * \invariant The given good name matches the one that was used when locating this market.
     */
    assert( aGoodName == mGoodName );
    
    /*!
     * \invariant The given region name matches the one that was used when locating this market.
     */
    assert( aRegionName == mRegionName );
    
    /*!
     * \invariant The given period matches the one that was used when locating this market.
     */
    assert( aPeriod == mPeriod );
    
    // Print a warning message when adding infinity values to the demand
    if ( !util::isValidNumber( aValue ) ) {
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::NOTICE );
        mainLog << "Error adding to demand in markeplace for: " << aGoodName << ", region: " << aRegionName << ", value: " << aValue << endl;
        return;
    }
    
    if ( mCachedMarket ) {
        mCachedMarket->addToDemand( aValue );
    }
    else if( aMustExist ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::NOTICE );
        mainLog << "Cannot add to demand for market as it does not exist: " << aGoodName << " " 
            << aRegionName << endl;
    }
}

/*!
 * \brief Get the market price.
 * \details Mimics the behavior of Marketplace::getPrice.
 * \param aGoodName The good of the market.
 * \param aRegionName The region getting the price.
 * \param aPeriod The period in which to get the price.
 * \param aMustExist Whether it is an error for the market not to exist.
 * \return The market price.
 * \see Marketplace::getPrice
 */  
double CachedMarket::getPrice( const string& aGoodName, const string& aRegionName, const int aPeriod,
                               bool aMustExist ) const {
    /*!
     * \invariant The given good name matches the one that was used when locating this market.
     */
    assert( aGoodName == mGoodName );
    
    /*!
     * \invariant The given region name matches the one that was used when locating this market.
     */
    assert( aRegionName == mRegionName );
    
    /*!
     * \invariant The given period matches the one that was used when locating this market.
     */
    assert( aPeriod == mPeriod );
    
    if( mCachedMarket ) {
        return mCachedMarket->getPrice();
    }
    
    if( aMustExist ) {
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::NOTICE );
        mainLog << "Called for price of non-existant market " << aGoodName << " in region " 
            << aRegionName << endl;
    }
    return Marketplace::NO_MARKET_PRICE;
}

/*!
 * \brief Get the market supply.
 * \details Mimics the behavior of Marketplace::getSupply.
 * \param aGoodName The good of the market.
 * \param aRegionName The region getting the supply.
 * \param aPeriod The period in which to get the supply.
 * \param aMustExist Whether it is an error for the market not to exist.
 * \return The market supply.
 * \see Marketplace::getSupply
 */
double CachedMarket::getSupply( const string& aGoodName, const string& aRegionName, const int aPeriod ) const {
    /*!
     * \invariant The given good name matches the one that was used when locating this market.
     */
    assert( aGoodName == mGoodName );
    
    /*!
     * \invariant The given region name matches the one that was used when locating this market.
     */
    assert( aRegionName == mRegionName );
    
    /*!
     * \invariant The given period matches the one that was used when locating this market.
     */
    assert( aPeriod == mPeriod );
    
    if ( mCachedMarket ) {
        return mCachedMarket->getSupply();
    }
    
    ILogger& mainLog = ILogger::getLogger( "main_log" );
    mainLog.setLevel( ILogger::NOTICE );
    mainLog << "Called for supply of non-existant market " << aGoodName << " in " << aRegionName << endl;
    return 0;
}

/*!
 * \brief Get the market demand.
 * \details Mimics the behavior of Marketplace::getDemand.
 * \param aGoodName The good of the market.
 * \param aRegionName The region getting the demand.
 * \param aPeriod The period in which to get the demand.
 * \param aMustExist Whether it is an error for the market not to exist.
 * \return The market demand.
 * \see Marketplace::getDemand
 */
double CachedMarket::getDemand(  const string& aGoodName, const string& aRegionName, const int aPeriod ) const {
    /*!
     * \invariant The given good name matches the one that was used when locating this market.
     */
    assert( aGoodName == mGoodName );
    
    /*!
     * \invariant The given region name matches the one that was used when locating this market.
     */
    assert( aRegionName == mRegionName );
    
    /*!
     * \invariant The given period matches the one that was used when locating this market.
     */
    assert( aPeriod == mPeriod );
    
    if ( mCachedMarket ) {
        return mCachedMarket->getDemand();
    }
    
    ILogger& mainLog = ILogger::getLogger( "main_log" );
    mainLog.setLevel( ILogger::NOTICE );
    mainLog << "Called for demand of non-existant market " << aGoodName << " in " << aRegionName << endl;
    return 0;
}

/*! \brief Get the information object for the specified market and period which
 *          can then be used to query for specific values.
 * \details Mimics the behavior of Marketplace::getMarketInfo.
 * \param aGoodName The good of the market for which to get the information object.
 * \param aRegionName The region used to find the market from which to get the
 *        information object.
 * \param aPeriod The period to fetch for which the information object. 
 * \param aMustExist Whether it is an error for the market not to exist.
 * \return A constant pointer to the market information object, null if the
 *         market does not exist.
 * \note This version of the function is required so that it can be called in
 *       constant functions. A second version is available which returns a
 *       mutable pointer.
 * \see Marketplace::getMarketInfo
 */
const IInfo* CachedMarket::getMarketInfo( const string& aGoodName, const string& aRegionName,
                                          const int aPeriod, const bool aMustExist ) const 
{
    /*!
     * \invariant The given good name matches the one that was used when locating this market.
     */
    assert( aGoodName == mGoodName );
    
    /*!
     * \invariant The given region name matches the one that was used when locating this market.
     */
    assert( aRegionName == mRegionName );
    
    /*!
     * \invariant The given period matches the one that was used when locating this market.
     */
    assert( aPeriod == mPeriod );

    const IInfo* info = 0;
    if ( mCachedMarket ) {
        info = mCachedMarket->getMarketInfo();
        /*! \invariant The market is required to return an information object
         *              that is non-null. 
         */
        assert( info );
    }
    
    // Report the error if requested.
    if( !info && aMustExist ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::NOTICE );
        mainLog << "Market info object cannot be returned because market "
            << aGoodName << " in " << aRegionName << " does not exist." << endl;
    }
    return info;
}

/*! \brief Get the information object for the specified market and period which
 *          can then be used to query, add, or modify specific values.
 * \details Mimics the behavior of Marketplace::getMarketInfo.
 * \param aGoodName The good of the market for which to get the information
 *        object.
 * \param aRegionName The region used to find the market from which to get the
 *        information object.
 * \param aPeriod The period to fetch for which the information object. 
 * \param aMustExist Whether it is an error for the market not to exist.
 * \return A mutable pointer to the market information object, null if the market
 *         does not exist.
 * \note This function returns a mutable pointer to the information object so it
 *       cannot be called from constant function.
 * \see Marketplace::getMarketInfo
 */
IInfo* CachedMarket::getMarketInfo( const string& aGoodName, const string& aRegionName,
                                    const int aPeriod, const bool aMustExist )
{
    /*!
     * \invariant The given good name matches the one that was used when locating this market.
     */
    assert( aGoodName == mGoodName );
    
    /*!
     * \invariant The given region name matches the one that was used when locating this market.
     */
    assert( aRegionName == mRegionName );
    
    /*!
     * \invariant The given period matches the one that was used when locating this market.
     */
    assert( aPeriod == mPeriod );
    
    IInfo* info = 0;
    if ( mCachedMarket ) {
        info = mCachedMarket->getMarketInfo();
        /*! \invariant The market is required to return an information object
         *              that is non-null. 
         */
        assert( info );
    }
    
    // Report the error if requested.
    if( !info && aMustExist ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::NOTICE );
        mainLog << "Market info object cannot be returned because market " 
            << aGoodName << " in " << aRegionName << " does not exist." << endl;
    }
    return info;
}
