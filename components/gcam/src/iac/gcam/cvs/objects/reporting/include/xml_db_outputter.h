#ifndef _XML_DB_OUTPUTTER_H_
#define _XML_DB_OUTPUTTER_H_
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
* \file xml_db_outputter.h
* \ingroup Objects
* \brief XMLDBOutputter class header file.
* \author Josh Lurz
*/

#if( __USE_XML_DB__ )
#include <iosfwd>
#include <sstream>
#include <stack>
#include <memory>
#include "util/base/include/default_visitor.h"

#include "dbxml/DbXml.hpp"
#include "dbxml/XmlInputStream.hpp"

class IndirectEmissionsCalculator;

/*! 
* \ingroup Objects
* \brief A visitor which writes model results to an XML database.
* \details
* \author Josh Lurz
*/

class XMLDBOutputter : public DefaultVisitor {
public:
    XMLDBOutputter();

    ~XMLDBOutputter();

    void finish() const;

    void startVisitScenario( const Scenario* aScenario, const int aPeriod );
    void endVisitScenario( const Scenario* aScenario, const int aPeriod );

    void startVisitOutputMetaData( const OutputMetaData* aOutputMetaData, const int aPeriod );
    void endVisitOutputMetaData( const OutputMetaData* aOutputMetaData, const int aPeriod );

    void startVisitWorld( const World* aWorld, const int aPeriod );
    void endVisitWorld( const World* aWorld, const int aPeriod );

    void startVisitRegion( const Region* aRegion, const int aPeriod );
    void endVisitRegion( const Region* aRegion, const int aPeriod );

    void startVisitRegionMiniCAM( const RegionMiniCAM* aRegionMiniCAM, const int aPeriod );
    void endVisitRegionMiniCAM( const RegionMiniCAM* aRegionMiniCAM, const int aPeriod );

    void startVisitRegionCGE( const RegionCGE* aRegionCGE, const int aPeriod );
    void endVisitRegionCGE( const RegionCGE* aRegionCGE, const int aPeriod );

    void startVisitResource( const AResource* aResource, const int aPeriod );
    void endVisitResource( const AResource* aResource, const int aPeriod );

    void startVisitSubResource( const SubResource* aSubResource, const int aPeriod );
    void endVisitSubResource( const SubResource* aSubResource, const int aPeriod );

    void startVisitGrade( const Grade* aGrade, const int aPeriod );
    void endVisitGrade( const Grade* aGrade, const int aPeriod );

    void startVisitSector( const Sector* aSector, const int aPeriod );
    void endVisitSector( const Sector* aSector, const int aPeriod );

    void startVisitSubsector( const Subsector* aSubsector, const int aPeriod );
    void endVisitSubsector( const Subsector* aSubsector, const int aPeriod );

    void startVisitEnergyFinalDemand( const EnergyFinalDemand* aEnergyFinalDemand, const int aPeriod );
    void endVisitEnergyFinalDemand( const EnergyFinalDemand* aEnergyFinalDemand, const int aPeriod );

    void startVisitBaseTechnology( const BaseTechnology* aBaseTech, const int aPeriod );
    void endVisitBaseTechnology( const BaseTechnology* aBaseTech, const int aPeriod );

    void startVisitTechnology( const Technology* aTechnology, const int aPeriod );
    void endVisitTechnology( const Technology* aTechnology, const int aPeriod );

    virtual void startVisitTranTechnology( const TranTechnology* aTranTechnology, const int aPeriod );
    virtual void endVisitTranTechnology( const TranTechnology* aTranTechnology, const int aPeriod );

    virtual void startVisitMiniCAMInput( const MiniCAMInput* aInput, const int aPeriod );
    virtual void endVisitMiniCAMInput( const MiniCAMInput* aInput, const int aPeriod );

    virtual void startVisitInput( const IInput* aInput, const int aPeriod );
    virtual void endVisitInput( const IInput* aInput, const int aPeriod );

    virtual void startVisitOutput( const IOutput* aOutput, const int aPeriod );
    virtual void endVisitOutput( const IOutput* aOutput, const int aPeriod );

    virtual void startVisitAgProductionTechnology( const AgProductionTechnology* AgProductionTechnology, const int aPeriod );
    virtual void endVisitAgProductionTechnology( const AgProductionTechnology* AgProductionTechnology, const int aPeriod );

    void startVisitGHG( const AGHG* aGHG, const int aPeriod );
    void endVisitGHG( const AGHG* aGHG, const int aPeriod );

    void startVisitMarketplace( const Marketplace* aMarketplace, const int aPeriod );
    void endVisitMarketplace( const Marketplace* aMarketplace, const int aPeriod );

    void startVisitMarket( const Market* aMarket, const int aPeriod );
    void endVisitMarket( const Market* aMarket, const int aPeriod );

    virtual void startVisitClimateModel( const IClimateModel* aClimateModel, const int aPeriod );
    virtual void endVisitClimateModel( const IClimateModel* aClimateModel, const int aPeriod );

    void startVisitDemographic( const Demographic* aDemographic, const int aPeriod );
    void endVisitDemographic( const Demographic* aDemographic, const int aPeriod );

    void startVisitPopulation( const Population* aPopulation, const int aPeriod );
    void endVisitPopulation( const Population* aPopulation, const int aPeriod );

    void startVisitPopulationMiniCAM( const PopulationMiniCAM* aPopulation, const int aPeriod );
    void endVisitPopulationMiniCAM( const PopulationMiniCAM* aPopulation, const int aPeriod );

    void startVisitPopulationSGMRate( const PopulationSGMRate* aPopulation, const int aPeriod );
    void endVisitPopulationSGMRate( const PopulationSGMRate* aPopulation, const int aPeriod );

    void startVisitPopulationSGMFixed( const PopulationSGMFixed* aPopulation, const int aPeriod );
    void endVisitPopulationSGMFixed( const PopulationSGMFixed* aPopulation, const int aPeriod );

    void startVisitAgeCohort( const AgeCohort* aAgeCohort, const int aPeriod );
    void endVisitAgeCohort( const AgeCohort* aAgeCohort, const int aPeriod );

    void startVisitGender( const Gender* aGender, const int aPeriod );
    void endVisitGender( const Gender* aGender, const int aPeriod );

    void startVisitGDP( const GDP* aGDP, const int aPeriod );
    void endVisitGDP( const GDP* aGDP, const int aPeriod );

    void startVisitCarbonCalc( const ICarbonCalc* aCarbonCalc, const int aPeriod );
    void endVisitCarbonCalc( const ICarbonCalc* aCarbonCalc, const int aPeriod );

    void startVisitLandNode( const LandNode* aLandNode, const int aPeriod );
    void endVisitLandNode( const LandNode* aLandNode, const int aPeriod );

    void startVisitLandLeaf( const LandLeaf* aLandLeaf, const int aPeriod );
    void endVisitLandLeaf( const LandLeaf* aLandLeaf, const int aPeriod );

    void startVisitExpenditure( const Expenditure* aExpenditure, const int aPeriod );
    void endVisitExpenditure( const Expenditure* aExpenditure, const int aPeriod );

    virtual void startVisitSGMInput( const SGMInput* aInput, const int aPeriod );
    virtual void endVisitSGMInput( const SGMInput* aInput, const int aPeriod );

    virtual void startVisitNodeInput( const NodeInput* aNodeInput, const int aPeriod );
    virtual void endVisitNodeInput( const NodeInput* aNodeInput, const int aPeriod );

    virtual void startVisitHouseholdConsumer( const HouseholdConsumer* aHouseholdConsumer, 
        const int aPeriod );
    virtual void endVisitHouseholdConsumer( const HouseholdConsumer* aHouseholdConsumer, 
        const int aPeriod );

    virtual void startVisitGovtConsumer( const GovtConsumer* aGovtConsumer, const int aPeriod );
    virtual void endVisitGovtConsumer( const GovtConsumer* aGovtConsumer, const int aPeriod );

    virtual void startVisitTradeConsumer( const TradeConsumer* aTradeConsumer, const int aPeriod );
    virtual void endVisitTradeConsumer( const TradeConsumer* aTradeConsumer, const int aPeriod );

    virtual void startVisitInvestConsumer( const InvestConsumer* aInvestConsumer, const int aPeriod );
    virtual void endVisitInvestConsumer( const InvestConsumer* aInvestConsumer, const int aPeriod );

    virtual void startVisitProductionTechnology( const ProductionTechnology* aProductionTechnology, 
        const int aPeriod );
    virtual void endVisitProductionTechnology( const ProductionTechnology* aProductionTechnology, 
        const int aPeriod );

    virtual void startVisitFactorSupply( const FactorSupply* aFactorySupply, const int aPeriod );
    virtual void endVisitFactorSupply( const FactorSupply* aFactorSupply, const int aPeriod );

    virtual void startVisitNationalAccount( const NationalAccount* aNationalAccount, const int aPeriod );
    virtual void endVisitNationalAccount( const NationalAccount* aNationalAccount, const int aPeriod );

    static bool appendData( const std::string& aData, const std::string& aLocation );
private:
    //! std::stringstream containing all the information until it is printed.
    mutable std::stringstream mBuffer;

    //! Current region name.
    std::string mCurrentRegion;

    //! Current sector name.
    std::string mCurrentSector;

    //! Current price unit.
    std::string mCurrentPriceUnit;

    //! Current output unit.
    std::string mCurrentOutputUnit;

    //! Current Input unit.
    std::string mCurrentInputUnit;

    //! Current market name.
    std::string mCurrentMarket;

    //! Current indirect emissions for the Technology. These are more easily
    //! calculated at the Technology but logically belong in the GHG output.
    objects::PeriodVector<double> mCurrIndirectEmissions;

    //! Tabs object.
    std::auto_ptr<Tabs> mTabs;

    //! Weak pointer to the current region's GDP object.
    const GDP* mGDP;

    //! Weak pointer to the current technology object.
    const Technology* mCurrentTechnology;
    
    //! A stack used to keep track of what needs to be written to the
    //! database.
    std::stack<std::stringstream*> mBufferStack;

    //! Indirect emissions calculator for the current region.
    std::auto_ptr<IndirectEmissionsCalculator> mIndirectEmissCalc;
    
    //! Flag to indicate if this was the first region processed
    bool mIsFirstRegion;

    /*! \brief Contains all objects necessary to operate on a container.
    * \details This struct defines the set of objects that must have the same
    *          lifetime so that the XML database outputter can operate on the
    *          container. The struct also ensures that the objects are deleted
    *          in the correct ordering to avoid accessing already deleted
    *          objects.
    * \note These objects must be in this order so destruction works correctly.
    */
    struct DBContainer {
        //! The database environment.
        DB_ENV* mDBEnvironment;

        //! The database manager.
        std::auto_ptr<DbXml::XmlManager> mManager;

        //! Wrapper around the XML container so the memory can be dynamically
        //! allocated.
        struct XMLContainerWrapper {
            XMLContainerWrapper( DbXml::XmlContainer aContainer );
            //! The XML container.
            DbXml::XmlContainer mContainer;
        };

        //! The wrapper around the XML container.
        std::auto_ptr<XMLContainerWrapper> mContainerWrapper;
        DBContainer();
        ~DBContainer();
    };
    static std::auto_ptr<DBContainer> createContainer();
    static const std::string createContainerName( const std::string& aScenarioName );
    bool appendBuffer( const std::string& aLocation );

    void writeItemToBuffer( const double aValue,
        const std::string& aName,
        std::ostream& out,
        const Tabs* tabs,
        const int aPeriod,
        const std::string& aUnit );

    void writeItem( const std::string& aName,
        const std::string& aUnit,
        const double aValue,
        const int aPeriod );

    void writeItemUsingYear( const std::string& aName,
        const std::string& aUnit,
        const double aValue,
        const int aYear );

    bool isTechnologyOperating( const int aPeriod );
    
    std::stringstream* popBufferStack();
    
    /*!
     * \brief Wrapper class around an stringstream to adapt the interface to the
     *        XmlInputStream interface used by the XML database.
     * \details Using this wrapper will allow us to import the XML string data
     *          generated by this visitor without having to create additional in
     *          memory copies.
     */
    class StringStreamXmlInputStream : public DbXml::XmlInputStream {
    public:
        StringStreamXmlInputStream( std::stringstream& aStream );
        virtual ~StringStreamXmlInputStream();
        
        // XmlInputStream methods
        virtual unsigned int curPos() const;
        virtual unsigned int readBytes( char *toFill, const unsigned int maxToRead );
    private:
        //! A weak pointer to the stringstream to wrap
        std::stringstream* mWrappedStringStream;
    };

};
#endif //  __USE_XML_DB__
#endif // _XML_DB_OUTPUTTER_H_
