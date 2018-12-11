#ifndef _PRODUCTION_TECHNOLOGY_H_
#define _PRODUCTION_TECHNOLOGY_H_
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
* \file production_technology.h
* \ingroup Objects
* \brief ProductionTechnology class header file.
* \author Pralit Patel
* \author Sonny Kim
*/

#include <string>
#include <xercesc/dom/DOMNode.hpp>
#include <memory>

#include "technologies/include/base_technology.h"
#include "functions/include/ifunction.h" // could remove if mTechChange was a pointer, but need a copyctor.
class Tabs;
class Demographic;
class MoreSectorInfo;
class IVisitor;
class IExpectedProfitRateCalculator;
class IShutdownDecider;

/*! 
 * \ingroup Objects
 * \brief A technology which accepts investment and produces output for SGM
 *        production sectors.
 * \details Represents a single vintage of an SGM production technology.  This technology
 *          is responsible calculating its own coefficients, passing its self forward,
 *          converting to an old vintage, calculating a technology cost, supply output
 *          and input demands based on investment, calculating emissions, and collecting
 *          taxes and various other accounting.
 *
 *          <b>XML specification for ProductionTechnology</b>
 *          - XML name: \c productionTechnology
 *          - Contained by: Subsector
 *          - Parsing inherited from class: BaseTechnology
 *          - Elements:
 *              - \c lifeTime TradeInput::lifeTime
 *                   The lifetime in years of this technology.
 *              - \c basePhysicalOutput TradeInput::mBasePhysicalOutput
 *                   Base year physical output of this technology.  This value is utilized
 *                   in conjunction with the total currency demand from the inputs to determine
 *                   the base year technology cost to use when initializing coefficients.  Note
 *                   that if this value is not provided we assume the cost is equal to the sector
 *                   price recieved.
 *              - \c fixed-investment TradeInput::mFixedInvestment
 *                   Fixed investment our output level of this technology.  Note that reading
 *                   this value implies this technology is a fixed technology and the investor will
 *                   not attempt to invest in this technology.  The technology cost will however
 *                   factor into the levelized cost of the sector.  Also note that even though
 *                   a fixed amount is read when the technology is an old vintage it is still
 *                   susceptible  to the shutdown coef.
 *              - \c technicalChangeHicks TradeInput::mTechChange.mHicksTechChange
 *                   Hicks neutral tech change applied in the initial investment year of this
 *                   technology.
 *              - \c technicalChangeEnergy TradeInput::mTechChange.mEnergyTechChange
 *                   Tech change for energy inputs applied in the initial investment year of this
 *                   technology.
 *              - \c technicalChangeMaterial TradeInput::mTechChange.mMaterialTechChange
 *                   Tech change for material inputs applied in the initial investment year of this
 *                   technology.
 *
 * \author Pralit Patel, Sonny Kim
 */
class ProductionTechnology : public BaseTechnology
{
    friend class SocialAccountingMatrix;
    friend class DemandComponentsTable;
    friend class SectorReport;
    friend class SGMGenTable;
    friend class SectorResults;
    friend class GovtResults;
    friend class XMLDBOutputter;
    // either these are friends or move isAvailable and isRetired to public
    friend class EnergyBalanceTable;
    friend class InputOutputTable;
public:
    ProductionTechnology();
    ProductionTechnology* clone() const;
    ~ProductionTechnology();

    virtual void copyParam( const BaseTechnology* baseTech,
                            const int aPeriod );

    void copyParamsInto( ProductionTechnology& prodTechIn,
                         const int aPeriod ) const;

    virtual void completeInit( const std::string& aRegionName,
                               const std::string& aSectorName,
                               const std::string& aSubsectorName );
    
    virtual void initCalc( const MoreSectorInfo* aMoreSectorInfo,
                           const std::string& aRegionName, 
                           const std::string& aSectorName,
                           NationalAccount& aNationalAccount, 
                           const Demographic* aDemographics,
                           const double aCapitalStock,
                           const int aPeriod );
    
    static const std::string& getXMLNameStatic();

    bool isNewInvestment( const int period ) const;

    virtual void operate( NationalAccount& aNationalAccount, const Demographic* aDemographic, 
        const MoreSectorInfo* aMoreSectorInfo, const std::string& aRegionName, const std::string& aSectorName, 
        const bool aIsNewVintageMode, const int aPeriod );

    virtual double setInvestment( const std::string& aRegionName, const double aAnnualInvestment,
                                  const double aTotalInvestment, const int aPeriod );

    virtual double getCapitalStock() const;
    virtual double getAnnualInvestment( const int aPeriod ) const;
    
    double getExpectedProfitRate( const NationalAccount& aNationalAccount,
                                  const std::string& aRegionName,
                                  const std::string& aSectorName,
                                  const IExpectedProfitRateCalculator* aExpProfitRateCalc,
                                  const double aInvestmentLogitExp,
                                  const bool aIsShareCalc,
                                  const bool aIsDistributing,
                                  const int aPeriod ) const;
    
    virtual double getCapitalOutputRatio( const IDistributor* aDistributor,
                                          const IExpectedProfitRateCalculator* aExpProfitRateCalc,
                                          const NationalAccount& aNationalAccount,
                                          const std::string& aRegionName,
                                          const std::string& aSectorName, 
                                          const int aPeriod ) const;

    double getFixedInvestment( const int aPeriod ) const;
    double distributeInvestment( const IDistributor* aDistributor, NationalAccount& aNationalAccount,
                                 const IExpectedProfitRateCalculator* aExpProfitRateCalc,
                                 const std::string& aRegionName, const std::string& aSectorName,
                                 const double aNewInvestment, const int aPeriod );
    void updateMarketplace( const std::string& sectorName, const std::string& regionName, 
        const int period );
    void csvSGMOutputFile( std::ostream& aFile, const int period ) const;
    void accept( IVisitor* aVisitor, const int aPeriod ) const;
    void setTypeHelper( TechnologyType* aTechType );
    void postCalc( const std::string& aRegionName, const std::string& aSectorName, const int aPeriod );
protected:
    virtual bool isCoefBased() const { return true; }
    virtual void calcEmissions( const std::string& aGoodName, const std::string& aRegionName, const int aPeriod );
    virtual const std::string& getXMLName() const;
    virtual bool XMLDerivedClassParse( const std::string& nodeName, const xercesc::DOMNode* curr );
    virtual void toInputXMLDerived( std::ostream& out, Tabs* tabs ) const;
    virtual void toDebugXMLDerived( const int period, std::ostream& out, Tabs* tabs ) const;
private:

    //! Costs stored by period for reporting.
    std::vector<double> mCostsReporting;

    //! Total profits stored by period.
    std::vector<double> mProfits;
    
    //! Structure containing types of technological change.
    TechChange mTechChange;
    
    //! Technology type containing this vintage. This may change.
    TechnologyType* mParentTechType;
    
    //! Expected profit rate stored for reporting.
    double mExpectedProfitRateReporting;
    
    /*! \brief An enum which describes the types of cached function return
    *          values for a specific year. 
    */
    enum CacheValues {
        //! The return value of the isAvailable function.
        AVAILABLE, 
        //! The return value of the isRetired function.
        RETIRED,
        //! The return value of the isNewInvestment function.
        NEW_INVESTMENT,
        //! A market for the last enum value, insert new enum values prior to this.
        END
    };
    //! A vector of cached results indexed by year for frequently called checks.
    std::vector<bool> mCachedValues;
    
    double mCapitalStock; //!< capital stock
    double indBusTax; //!< indirect business tax
    double mBasePhysicalOutput; //!< base year physical output for calculating conversion factor
    double mAnnualInvestment; //!< Annual investment
    double mFixedInvestment; //!< Quantity of fixed investment in the initial year for the technology.
    
    // All times must be in years because some periods are < 0.
    
    //! The period for which all cached values stored in mCachedValues are valid for.
    int mValidCachePeriod;

    int lifeTime; //!< nameplate lifetime of the technology
    int delayedInvestTime; //!< number of years between initial investment and first operation
    int maxLifeTime; //!< maximum allowable lifetime
    int retrofitLifeTime; //!< lifetime of the technology renovation
    int periodIniInvest; //!< number of periods until initial investment
    int periodInvestUnallowed; //!< period in which investment is no longer allowed

    
    void calcTaxes( NationalAccount& aNationalAccount, const MoreSectorInfo* aMoreSectorInfo, 
        const std::string& aRegionName, const std::string aSectorName, const int aPeriod );
    
    bool calcIsNewInvestment( const int aPeriod ) const;
    
    inline bool isAvailable( const int aPeriod ) const;
    bool calcIsAvailable( const int aPeriod ) const;
    
    inline bool isRetired( const int aPeriod ) const;
    bool calcIsRetired( const int aPeriod ) const;

    double calcShutdownCoef( const std::string& aRegionName,
                             const std::string& aSectorName,
                             const int aPeriod ) const;
};

/*! \brief Return whether a technology is available to go online.
* \param aPeriod The current period.
* \return Whether the technology has gone online.
* \author Josh Lurz
*/
bool ProductionTechnology::isAvailable( const int aPeriod ) const {
    // Check if we have a valid cached result for the period.
    if( aPeriod == mValidCachePeriod ){
        return mCachedValues[ AVAILABLE ];
    }
    // There is no cached result available for this period, so calculate it
    // dynamically.
    return calcIsAvailable( aPeriod );
}

/*! \brief Return whether a technology has been retired yet.
* \param aPeriod The current period.
* \return Whether the technology has been retired.
* \author Josh Lurz
*/
bool ProductionTechnology::isRetired( const int aPeriod ) const {
    // Check if we have a valid cached result for the period.
    if( aPeriod == mValidCachePeriod ){
        return mCachedValues[ RETIRED ];
    }
    // There is no cached result available for this period, so calculate it
    // dynamically.
    return calcIsRetired( aPeriod );
}

#endif // _PRODUCTION_TECHNOLOGY_H_
