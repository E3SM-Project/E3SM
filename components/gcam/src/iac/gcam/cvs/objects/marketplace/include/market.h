#ifndef _MARKET_H_
#define _MARKET_H_
#if defined(_MSC_VER_)
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
* \file market.h
* \ingroup Objects
* \brief The Market class header file.
* \author Sonny Kim
*/

#include <vector>
#include <memory>
#include "marketplace/include/imarket_type.h"
#include "util/base/include/ivisitable.h"

class IInfo;
class Tabs;
class IVisitor;
namespace objects {
    class Atom;
}

/*!
 * \ingroup Objects
 * \brief A single market, or equation, in the model.
 * \details A Market conceptually represents the trade for a single good in an
 *          area in which there are no transportation costs. See Marketplace for
 *          an explanation of a market region vs. model region. The market has a
 *          price, supply and demand for the good. It also contains an object
 *          with additional information about the good, the market info. Market
 *          objects may also be used in cases where a true market is not
 *          required, but a solved equation. See the TrialValueMarket for an
 *          explanation of this feature.
 *
 *          The Market's functions are divided into two main areas:
 
 *          - There are setters and accessors which are called whenever a
 *            Marketplace function to set or get supply, demand, or price is
 *            called. These methods may be overridden by derived Market classes
 *            to add different behaviors to the markets. Because of this, calls
 *            to a getter may not return the Market variable of the same name.
 *            For example, getSupply does not necessarily return Market::supply.
 *            It will return whatever conceptually is the supply.
 *
 *          - There are functions to directly set and get the underlying supply,
 *            demand, and price variables. They are named
 *            (set|get)Raw(Supply|Demand|Price). These functions modify the left
 *            hand side, right hand side, and trial values of the equation the
 *            Market represents. These are only used by the solution mechanism,
 *            and cannot be overridden.
 *
 * \author Sonny Kim
 */

class Market: public IVisitable
{
    friend class XMLDBOutputter;
public:
    Market( const std::string& goodNameIn, const std::string& regionNameIn, const int periodIn );
    virtual ~Market();
    static std::auto_ptr<Market> createMarket( const IMarketType::Type aMarketType,
        const std::string& aGoodName, const std::string& aRegionName, const int aPeriod );
    void toDebugXML( const int period, std::ostream& out, Tabs* tabs ) const;
    static const std::string& getXMLNameStatic();
    void addRegion( const std::string& aRegion );
    const std::vector<const objects::Atom*>& getContainedRegions() const;

    virtual void initPrice() = 0;
    virtual void setPrice( const double priceIn ) = 0;
    void setRawPrice( const double priceIn );
    virtual void set_price_to_last_if_default( const double lastPrice ) = 0;
    virtual void set_price_to_last( const double lastPrice ) = 0;
    virtual double getPrice() const = 0;
    double getRawPrice() const;
    double getStoredRawPrice() const;

    virtual void nullDemand();
    void setRawDemand( const double value );
    virtual void addToDemand( const double demandIn ) = 0;
    virtual void removeFromRawDemand( const double demandIn );
    double getRawDemand() const;
    double getStoredRawDemand() const;
    virtual double getDemand() const = 0;

    virtual void nullSupply() = 0;
    double getRawSupply() const;
    double getStoredRawSupply() const;
    void setRawSupply( const double supplyIn );
    void removeFromRawSupply( const double supplyIn );
    virtual double getSupply() const = 0;
    virtual void addToSupply( const double supplyIn );
    const std::string& getName() const;
    const std::string& getRegionName() const;
    const std::string& getGoodName() const;
    const IInfo* getMarketInfo() const;
    IInfo* getMarketInfo();
    void storeInfo();
    void restoreInfo();
    void store_original_price();
    void restore_original_price();

    void setSolveMarket( const bool doSolve );
    virtual bool meetsSpecialSolutionCriteria() const = 0;
    virtual bool shouldSolve() const = 0;
    virtual bool shouldSolveNR() const = 0;
    bool isSolvable() const;

    /*!
    * \brief Return the type of the market as defined by the IMarketTypeEnum
    *        which is unique for each derived market class.
    * \return The type of the market.
    */
    virtual IMarketType::Type getType() const = 0;
    static const std::string& convert_type_to_string( const IMarketType::Type aType );

    void accept( IVisitor* aVisitor, const int aPeriod ) const;
protected:
    Market( const Market& aMarket );

    //! The name of the market.
    std::string mName;
    
    //! The good the market represents
    std::string good;
    
    //! The region of the market.
    std::string region;
    
    //! Whether to solve the market given other constraints are satisfied.
    bool solveMarket;
    
    //! The period the market is valid in.
    int period;
    
    //! The market price.
    double price;
    
    //! The stored market price.
    double storedPrice;
    
    //! The original market price.
    double original_price;
    
    //! The market demand.
    double demand;
    
    //! The stored demand.
    double storedDemand;
    
    //! The market supply.
    double supply;
    
    //! The stored supply.
    double storedSupply;
    
    //! Vector of atoms of all regions contained within this market.
    std::vector <const objects::Atom*> mContainedRegions;
    
    //! Object containing information related to the market.
    std::auto_ptr<IInfo> mMarketInfo;

        /*! \brief Add additional information to the debug xml stream for derived
    *          classes.
    * \details This method is inherited from by derived class if they which to
    *          add any additional information to the printout of the class.
    * \param out Output stream to print to.
    * \param tabs A tabs object responsible for printing the correct number of
    *        tabs. 
    */
    virtual void toDebugXMLDerived( std::ostream& out, Tabs* tabs ) const = 0;

};

#endif // _MARKET_H_
