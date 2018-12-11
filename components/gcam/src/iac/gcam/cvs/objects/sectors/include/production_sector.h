#ifndef _PRODUCTION_SECTOR_H_
#define _PRODUCTION_SECTOR_H_
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
* \file production_sector.h
* \ingroup Objects
* \brief The ProductionSector class header file.
* \author Sonny Kim
*/
#include <string>
#include <memory>
#include <map>
#include "sectors/include/sector.h"

class IInvestor;
class Demographic;
class NationalAccount;
class GDP;
class IInfo;
class DependencyFinder;

/*! 
 * \ingroup Objects
 * \brief This class represents a single production sector.
 * \details A ProductionSector basically consists of subsectors and an IInvestor to
 *          handle new investments.  It also contains other meta data such as wether
 *          it is an energy good, primary energy, etc however all of that info will
 *          be placed into the Sector Info object.  When operating a ProductionSector
 *          the follwoing steps will occur:
 *           - Calculate the price recieved and store it into a market info so that a
 *             technology will know how much to set aside for production taxes.
 *           - Operate old capital technologies.  Since old vintages already know how
 *             much investment they have they can operate it to calculate their supplies
 *             and demands.  Note that there are no ordering depencies which forces these
 *             to operate before new vintage technologies.
 *           - Calculate and distribute new investment to the technologies.  This task is
 *             delegated to the IInvestor object.  This must be done before new vintage
 *             technologies can operate.
 *           - Operate new vintage technologies so that they may add there demands into
 *             the marketplace.
 *
 *          <b>XML specification for ProductionSector</b>
 *          - XML name: \c productionSector
 *          - Contained by: Region
 *          - Parsing inherited from class: Sector
 *          - Attributes:
 *          - Elements:
 *              - \c Accelerator::getXMLNameStatic() ProductionSector::mInvestor
 *              - \c MarketBasedInvestor::getXMLNameStatic() ProductionSector::mInvestor
 *              - \c market-name ProductionSector::mMarketName
 *              - \c FixedPricePath ProductionSector::mIsFixedPrice
 *              - \c numeraire ProductionSector::mIsNumeraireSector
 *              - \c ghgEmissCoef ProductionSector::ghgEmissCoefMap
 *              - \c IsEnergyGood ProductionSector::mIsEnergyGood
 *              - \c IsPrimaryEnergyGood ProductionSector::mIsPrimaryEnergyGood
 *              - \c IsSecondaryEnergyGood ProductionSector::mIsSecondaryEnergyGood
 *              - \c sectorprice ProductionSector::mFixedPrices
 *
 * \author Sonny Kim
 */
class ProductionSector: public Sector
{

public:
    explicit ProductionSector ( const std::string& aRegionName );
    virtual ~ProductionSector();
    void calcFinalSupplyPrice( const GDP* aGDP, const int aPeriod ){};
    void supply( const GDP* aGDP, const int aPeriod ){};
    static const std::string& getXMLNameStatic();
    
    virtual void completeInit( const IInfo* aRegionInfo,
                               DependencyFinder* aDepFinder,
                               ILandAllocator* aLandAllocator );

    double getOutput( const int aPeriod ) const;
    
    virtual void initCalc( NationalAccount* aNationalAccount,
                           const Demographic* aDemographics,
                           const int aPeriod );

    virtual void operate( NationalAccount& aNationalAccount, const Demographic* aDemographics, const int aPeriod ); // Passing demographic here is not good.
    virtual void postCalc( const int aPeriod );
    virtual void accept( IVisitor* aVisitor, const int aPeriod ) const;

    virtual void dbOutput( const GDP* aGDP,
                           const IndirectEmissionsCalculator* aIndEmissCalc ) const {}
protected:
    std::map<std::string,double> ghgEmissCoefMap; //! Map of ghg name to emission coefficient
    void setMarket();

    virtual double getPrice( const GDP* aGDP,
                             const int aPeriod ) const;

    virtual const std::string& getXMLName() const;
    virtual bool XMLDerivedClassParse( const std::string& nodeName, const xercesc::DOMNode* curr ); 
    virtual void toInputXMLDerived( std::ostream& out, Tabs* tabs ) const;
    virtual void toDebugXMLDerived( const int period, std::ostream& out, Tabs* tabs ) const;
private:
    //! Vector of read-in prices for the production sector which will be used as
    //! fixed prices if mIsFixedPrice is true.
    std::vector<double> mFixedPrices;

    //! The market region into which the sector is exporting.
    std::string mMarketName;

    //! Whether this sector is on a fixed price path.
    bool mIsFixedPrice;

    //! Whether this is a numeraire sector.
    bool mIsNumeraireSector;

    //! If the sector has an energy product.
    bool mIsEnergyGood;

    //! If the sector has a primary energy product.
    bool mIsPrimaryEnergyGood;

    //! If the sector has a secondary energy product.
    bool mIsSecondaryEnergyGood;

    //! Object responsible for determining levels of investment for this sector
    //! in its various technologies. Different types of investment objects may
    //! be read in to change the investment behavior.
    std::auto_ptr<IInvestor> mInvestor;
    
    void calcInvestment( const Demographic* aDemographic, NationalAccount& aNationalAccount, const int aPeriod );
    void operateOldCapital( const Demographic* aDemographic, NationalAccount& aNationalAccount, const int aPeriod );
    void operateNewCapital( const Demographic* aDemographic, NationalAccount& aNationalAccount, const int aPeriod );
    void calcPriceReceived( const int period );
};

#endif // _PRODUCTION_SECTOR_H_
