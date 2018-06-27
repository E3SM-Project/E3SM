#ifndef _CACHED_MARKET_H_
#define _CACHED_MARKET_H_
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
 * \file cached_market.h
 * \ingroup Objects
 * \brief The CachedMarket class header file.
 * \author Pralit Patel
 */

#include <string>

class Market;
class IInfo;

/*!
 * \brief A wrapper around a market which has been located through the marketplace.
 * \details This class along with Marketplace::locateMarket allows objects in the
 *          model to pay the relatively time intensive penalty of locating a market
 *          once and the subsequent operations are free.  The methods provided
 *          through this class mimic the behavior of the corresponding methods in
 *          Marketplace and thus could be used as a replacement for calls directly
 *          to the Marketplace.  It would be ideal to replace calls in objects
 *          for which there would be a large number of instance of or contain many
 *          calls to the Marketplace.
 *
 * \author Pralit Patel
 * \warning It is up to the user to ensure the cached market matches the intended
 *          market, i.e. the good name, region name, and period have not changed.
 *          To ensure this does not happen a user could run in debug mode to check
 *          asserts.
 */
class CachedMarket
{
public:
    CachedMarket( const std::string& aGoodName, const std::string& aRegionName, const int aPeriod, Market* aLocatedMarket );
    ~CachedMarket();

    void setPrice( const std::string& aGoodName, const std::string& aRegionName, const double aValue,
                   const int aPeriod, bool aMustExist = true );
    
    void addToSupply( const std::string& aGoodName, const std::string& aRegionName, const double aValue,
                      const int aPeriod, bool aMustExist = true );
    
    void addToDemand( const std::string& aGoodName, const std::string& aRegionName, const double aValue,
                      const int aPeriod, bool aMustExist = true );
    
    double getPrice( const std::string& aGoodName, const std::string& aRegionName, const int aPeriod,
                     bool aMustExist = true ) const;
    
    double getSupply( const std::string& aGoodName, const std::string& aRegionName,
                      const int aPeriod ) const;
    
    double getDemand( const std::string& aGoodName, const std::string& aRegionName,
                      const int aPeriod ) const;
    
    const IInfo* getMarketInfo( const std::string& aGoodName, const std::string& aRegionName,
                               const int aPeriod, const bool aMustExist ) const;
    
    IInfo* getMarketInfo( const std::string& aGoodName, const std::string& aRegionName,
                         const int aPeriod, const bool aMustExist );
private:
#if(!NDEBUG)
    //! The good name used when this market was located.  Used for debugging.
    const std::string mGoodName;
    
    //! The region name used when this market was located.  Used for debugging.
    const std::string mRegionName;
    
    //! The period used when this market was located.  Used for debugging.
    const int mPeriod;
#endif
    //! The actual market which is cached.
    Market* mCachedMarket;
};

#endif // _CACHED_MARKET_H_
