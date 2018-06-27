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
* \file market_locator.cpp
* \ingroup Objects
* \brief MarketLocator class source file.
* \author Josh Lurz
*/

#include "util/base/include/definitions.h"
#include <cassert>
#include <string>

#include "marketplace/include/market_locator.h"
#include "util/base/include/hash_map.h"

#define PERFORM_TIMING 0
#if PERFORM_TIMING
#include <iostream>
#include "util/base/include/timer.h"
// Static variables used for timing. Static variables are automatically
// initialized to zero.
static double gTotalLookupTime;
static int gNumLookups;
#endif

using namespace std;



/*! \brief Constructor */
MarketLocator::MarketLocator(){
    mLastRegionLookup = 0;
    const unsigned int MARKET_REGION_LIST_SIZE = 71;
    mRegionList.reset( new RegionMarketList( MARKET_REGION_LIST_SIZE ) );
    mMarketList.reset( new RegionMarketList( MARKET_REGION_LIST_SIZE ) );
}

//! Destructor
MarketLocator::~MarketLocator(){
#if PERFORM_TIMING
    cout << "Total time spent in lookups: " << gTotalLookupTime << " seconds." << endl
         << "Average time per lookup: " << gTotalLookupTime / static_cast<double>( gNumLookups )
         << " seconds." << endl;
#endif
}

/*! \brief Add a market to the locator.
* \details This function adds a market's position to the MarketLocator.
* \param aMarket The name of the market area to which this market is being
*        added.
* \param aRegion The name of the region containing the new market.
* \param aGoodName The name of the Good this market number represents.
* \param aUniqueNumber A unique market index to use if the market does not
*        already exist.
* \return The market number that was either already existing or added.
*/
int MarketLocator::addMarket( const string& aMarket,
                              const string& aRegion,
                              const string& aGoodName,
                              const int aUniqueNumber )
{
    // Check if the market area exists in the market area list.
    RegionMarketList::iterator iter = mMarketList->find( aMarket );
    
    int goodNumber;
    // The market area does not exist. Create a new entry.
    if( iter == mMarketList->end() ){
        boost::shared_ptr<RegionOrMarketNode> newMarketNode( new RegionOrMarketNode( aMarket ) );
        // Add the node to the hashmap.
        mMarketList->insert( make_pair( aMarket, newMarketNode ) );

        // Add the item to it.
        goodNumber = newMarketNode->addGood( aGoodName, aUniqueNumber );
    }
    else {
        // The market area already exists. Add the item to it.
        goodNumber = iter->second->addGood( aGoodName, aUniqueNumber );
    }

    // Check if the region exists in the region list.
    iter = mRegionList->find( aRegion );

    // The region does not exist. Create a new entry.
    if( iter == mRegionList->end() ){
        boost::shared_ptr<RegionOrMarketNode> newRegionNode( new RegionOrMarketNode( aRegion ) );
        // Add the new region to the region list.
        mRegionList->insert( make_pair( aRegion, newRegionNode ) );
        
        // Add the item to the region list.
        newRegionNode->addGood( aGoodName, goodNumber );
    }
    else {
        // The region already exists. Add the item to it.
        iter->second->addGood( aGoodName, goodNumber );
    }

    // Return the good number used.
    return goodNumber;
}

/*! \brief Find the market number for a given region and good.
* \param aRegion Region for which to search.
* \param aGoodName Good for which to search.
* \return The market number or MARKET_NOT_FOUND if it is not present.
*/
int MarketLocator::getMarketNumber( const string& aRegion, const string& aGoodName ) const {
    // Compile in extra timing. Note that timing causes significant overhead, so
    // timed runs will take longer. The result is useful to compare across timed
    // runs, not vs non-timed runs.
#if PERFORM_TIMING
    Timer timer;
    timer.start();

    // Lookup the market in the marketList.
    int marketNumber = getMarketNumberInternal( aRegion, aGoodName );
    timer.stop();
    gTotalLookupTime += timer.getTimeDifference();
    ++gNumLookups;
    return marketNumber;
#else
    return getMarketNumberInternal( aRegion, aGoodName );
#endif
}

/*! \brief Internal calculation which determines the market number from a region
*          and good name.
* \details Performs the calculation which determines the market number from a
*          region and good name.
* \param aRegion Region for which to search.
* \param aGoodName Good for which to search.
* \return The number of the associated market or MARKET_NOT_FOUND if it does not
*         exist.
*/
int MarketLocator::getMarketNumberInternal( const string& aRegion,
                                            const string& aGoodName ) const
{
    // First check if the cached market number matched the region.
    const RegionOrMarketNode* region;
    if( mLastRegionLookup && mLastRegionLookup->getName() == aRegion ){
        region = mLastRegionLookup;
    }
    else {
        RegionMarketList::const_iterator iter = mRegionList->find( aRegion );
        // Check if the region was found.
        if( iter != mRegionList->end() ){
            /*! \invariant If the key is in the list the node must be non-null. */
            assert( iter->second.get() );

            // Cache the region found.
            mLastRegionLookup = region = iter->second.get();
        }
        // Search failed.
        else {
            region = 0;
        }
    }

    // If the region lookup succeeded search for the good, otherwise return
    // market not found.
    return region ? region->getMarketNumber( aGoodName ) : MARKET_NOT_FOUND;
}

//! Constructor
MarketLocator::RegionOrMarketNode::RegionOrMarketNode( const string& aName ):
mName( aName ){
    const unsigned int SECTOR_LIST_SIZE = 51;
    mSectorNodeList.reset( new SectorNodeList( SECTOR_LIST_SIZE ) );
}

//! Destructor
MarketLocator::RegionOrMarketNode::~RegionOrMarketNode(){
}

/*! \brief Add a Good to the RegionOrMarketNode.
* \param aGoodName Name of the Good to add.
* \param aUniqueNumber A unique market location to use if the good is added to
*        the list.
* \return aUniqueNumber if the good was added to the list, the market number if
*         it already existed.
*/
int MarketLocator::RegionOrMarketNode::addGood( const string& aGoodName,
                                                const int aUniqueNumber )
{
    // Check if it exists in the good list.
    SectorNodeList::iterator iter = mSectorNodeList->find( aGoodName );

    // Check if the good was found.
    if( iter != mSectorNodeList->end() ){
        /*! \invariant If the key is in the list the node must be non-null. */
        assert( iter->second.get() );

        // Return the good number.
        return iter->second->mNumber;
    }

    // The good node does not exist. Create a new one.
    boost::shared_ptr<GoodNode> newGoodNode( new GoodNode( aGoodName, aUniqueNumber ) );

    // Add the new node to the hashmap.
    mSectorNodeList->insert( make_pair( aGoodName, newGoodNode ) );

    // Return the new index.
    return aUniqueNumber;
}

/*! \brief Find a market number given a Good name.
* \param aGoodName The name of the Good to search for.
* \return The market number for the Good, MARKET_NOT_FOUND otherwise.
*/
int MarketLocator::RegionOrMarketNode::getMarketNumber( const string& aGoodName ) const {
    // Check if it exists in the good list.
    SectorNodeList::const_iterator iter = mSectorNodeList->find( aGoodName );
    if( iter != mSectorNodeList->end() ){
        return iter->second->mNumber;
    }

    // Return that the market was not found.
    return MARKET_NOT_FOUND;
}

//! Constructor
MarketLocator::GoodNode::GoodNode( const string& aName, const int aMarketNumber ):
mName( aName ),
mNumber( aMarketNumber )
{
}
