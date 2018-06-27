#ifndef _IVISITOR_H_
#define _IVISITOR_H_
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
* \file ivisitor.h
* \ingroup Objects
* \brief IVisitor class header file.
* \author Josh Lurz
*/

#include <string>

class World;
class Region;
class RegionCGE;
class RegionMiniCAM;
class NationalAccount;
class Demographic;
class Sector;
class Subsector;
class BuildingDemandSubSector;
class BaseTechnology;
class Consumer;
class HouseholdConsumer;
class GovtConsumer;
class InvestConsumer;
class TradeConsumer;
class ProductionTechnology;
class DemandInput;
class ProductionInput;
class SGMInput;
class FactorSupply;
class ProductionSector;
class AResource;
class Technology;
class AFinalDemand;
class EnergyFinalDemand;
class Expenditure;
class AResource;
class Scenario;
class AGHG;
class OutputMetaData;
class Marketplace;
class Market;
class SubResource;
class SubRenewableResource;
class Grade;
class Population;
class PopulationMiniCAM;
class PopulationSGMFixed;
class PopulationSGMRate;
class AgeCohort;
class Gender;
class Male;
class Female;
class GDP;
class MiniCAMInput;
class IOutput;
class LandLeaf;
class UnmanagedLandLeaf;
class LandNode;
class ICarbonCalc;
class IClimateModel;
class IInput;
class TranTechnology;
class LandUseHistory;
class AgProductionTechnology;
class BuildingGenericDmdTechnology;
class NodeInput;

/*!
 * \brief An interface to a class which visits every node in the tree and
 *        optionally performs an operation on each.
 * \details Any class which implements the IVisitor interface may be passed to
 *          the accept method of any class that implements the Visitable
 *          interface. Once the Visitable class accepts the visitor, it will
 *          call startVisit for itself, call accept on all its children with the
 *          visitor, and then call endVisit with itself. The visitor may perform
 *          any desired work in the start and end visit methods, or it may
 *          choose to do nothing. The visitor will always perform a depth first
 *          traversal, which means it will visit all children of an item before
 *          visiting the next sibling of an item.
 *
 */
class IVisitor {
public:
    inline virtual ~IVisitor();
    virtual void finish() const = 0;

    virtual void startVisitScenario( const Scenario* aScenario, const int aPeriod ) = 0;
    virtual void endVisitScenario( const Scenario* aScenario, const int aPeriod ) = 0;
    virtual void startVisitWorld( const World* aWorld, const int aPeriod ) = 0;
    virtual void endVisitWorld( const World* aWorld, const int aPeriod ) = 0;

    virtual void startVisitRegion( const Region* aRegion, const int aPeriod ) = 0;
    virtual void endVisitRegion( const Region* aRegion, const int aPeriod ) = 0;

    virtual void startVisitRegionMiniCAM( const RegionMiniCAM* aRegion,
                                          const int aPeriod ) = 0;

    virtual void endVisitRegionMiniCAM( const RegionMiniCAM* aRegion,
                                        const int aPeriod ) = 0;

    virtual void startVisitRegionCGE( const RegionCGE* aRegionCGE, const int aPeriod ) = 0;
    virtual void endVisitRegionCGE( const RegionCGE* aRegionCGE, const int aPeriod ) = 0;
    
    virtual void startVisitDemographic( const Demographic* aDemographic, const int aPeriod ) = 0;
    virtual void endVisitDemographic( const Demographic* aDemographic, const int aPeriod ) = 0;

    virtual void startVisitPopulation( const Population* aPopulation, const int aPeriod ) = 0;
    virtual void endVisitPopulation( const Population* aPopulation, const int aPeriod ) = 0;

    virtual void startVisitPopulationMiniCAM( const PopulationMiniCAM* aPopulation, const int aPeriod ) = 0;
    virtual void endVisitPopulationMiniCAM( const PopulationMiniCAM* aPopulation, const int aPeriod ) = 0;

    virtual void startVisitPopulationSGMFixed( const PopulationSGMFixed* aPopulation, const int aPeriod ) = 0;
    virtual void endVisitPopulationSGMFixed( const PopulationSGMFixed* aPopulation, const int aPeriod ) = 0;

    virtual void startVisitPopulationSGMRate( const PopulationSGMRate* aPopulation, const int aPeriod ) = 0;
    virtual void endVisitPopulationSGMRate( const PopulationSGMRate* aPopulation, const int aPeriod ) = 0;

    virtual void startVisitAgeCohort( const AgeCohort* aAgeCohort, const int aPeriod ) = 0;
    virtual void endVisitAgeCohort( const AgeCohort* aAgeCohort, const int aPeriod ) = 0;
    
    virtual void startVisitGender( const Gender* aGender, const int aPeriod ) = 0;
    virtual void endVisitGender( const Gender* aGender, const int aPeriod ) = 0;

    virtual void startVisitFemale( const Female* aFemale, const int aPeriod ) = 0;
    virtual void endVisitFemale( const Female* aFemale, const int aPeriod ) = 0;

    virtual void startVisitMale( const Male* aMale, const int aPeriod ) = 0;
    virtual void endVisitMale( const Male* aMale, const int aPeriod ) = 0;

    virtual void startVisitResource( const AResource* aResource, const int aPeriod ) = 0;
    virtual void endVisitResource( const AResource* aResource, const int aPeriod ) = 0;
    
    virtual void startVisitSubResource( const SubResource* aSubResource, const int aPeriod ) = 0;
    virtual void endVisitSubResource( const SubResource* aSubResource, const int aPeriod ) = 0;

    virtual void startVisitSubRenewableResource( const SubRenewableResource* aSubResource, const int aPeriod ) = 0;
    virtual void endVisitSubRenewableResource( const SubRenewableResource* aSubResource, const int aPeriod ) = 0;

    virtual void startVisitGrade( const Grade* aGrade, const int aPeriod ) = 0;
    virtual void endVisitGrade( const Grade* aGrade, const int aPeriod ) = 0;

    virtual void startVisitSector( const Sector* aSector, const int aPeriod ) = 0;
    virtual void endVisitSector( const Sector* aSector, const int aPeriod ) = 0;
    
    virtual void startVisitFinalDemand( const AFinalDemand* aFinalDemand, const int aPeriod ) = 0;
    virtual void endVisitFinalDemand( const AFinalDemand* aFinalDemand, const int aPeriod ) = 0;

    virtual void startVisitEnergyFinalDemand( const EnergyFinalDemand* aEnergyFinalDemand, const int aPeriod ) = 0;
    virtual void endVisitEnergyFinalDemand( const EnergyFinalDemand* aEnergyFinalDemand, const int aPeriod ) = 0;

    virtual void startVisitProductionSector( const ProductionSector* aProdSector, const int aPeriod ) = 0;
    virtual void endVisitProductionSector( const ProductionSector* aProdSector, const int aPeriod ) = 0;
    
    virtual void startVisitSubsector( const Subsector* aSubsector, const int aPeriod ) = 0;
    virtual void endVisitSubsector( const Subsector* aSubsector, const int aPeriod ) = 0;
    
    virtual void startVisitBuildingDemandSubsector( const BuildingDemandSubSector* aSubsector,
                                                    const int aPeriod ) = 0;
    virtual void endVisitBuildingDemandSubsector( const BuildingDemandSubSector* aSubsector,
                                                  const int aPeriod ) = 0;

    virtual void startVisitTechnology( const Technology* aTechnology, const int aPeriod ) = 0;
    virtual void endVisitTechnology( const Technology* aTechnology, const int aPeriod ) = 0;

    virtual void startVisitAgProductionTechnology( const AgProductionTechnology* aTechnology, const int aPeriod ) = 0;
    virtual void endVisitAgProductionTechnology( const AgProductionTechnology* aTechnology, const int aPeriod ) = 0;

    virtual void startVisitBaseTechnology( const BaseTechnology* aBaseTechnology, const int aPeriod ) = 0;
    virtual void endVisitBaseTechnology( const BaseTechnology* aBaseTechnology, const int aPeriod ) = 0;

    virtual void startVisitConsumer( const Consumer* aConsumer, const int aPeriod ) = 0;
    virtual void endVisitConsumer( const Consumer* aConsumer, const int aPeriod ) = 0;

    virtual void startVisitHouseholdConsumer( const HouseholdConsumer* aHouseholdConsumer, 
        const int aPeriod ) = 0;
    virtual void endVisitHouseholdConsumer( const HouseholdConsumer* aHouseholdConsumer, 
        const int aPeriod ) = 0;

    virtual void startVisitGovtConsumer( const GovtConsumer* aGovtConsumer, const int aPeriod ) = 0;
    virtual void endVisitGovtConsumer( const GovtConsumer* aGovtConsumer, const int aPeriod ) = 0;

    virtual void startVisitInvestConsumer( const InvestConsumer* aInvestConsumer, const int aPeriod ) = 0;
    virtual void endVisitInvestConsumer( const InvestConsumer* aInvestConsumer, const int aPeriod ) = 0;

    virtual void startVisitTradeConsumer( const TradeConsumer* aTradeConsumer, const int aPeriod ) = 0;
    virtual void endVisitTradeConsumer( const TradeConsumer* aTradeConsumer, const int aPeriod ) = 0;

    virtual void startVisitProductionTechnology( const ProductionTechnology* aProductionTechnology, 
        const int aPeriod ) = 0;
    virtual void endVisitProductionTechnology( const ProductionTechnology* aProductionTechnology, 
        const int aPeriod ) = 0;

    virtual void startVisitFactorSupply( const FactorSupply* aFactorSupply, const int aPeriod ) = 0;
    virtual void endVisitFactorSupply( const FactorSupply* aFactorSupply, const int aPeriod ) = 0;

    virtual void startVisitNationalAccount( const NationalAccount* aNationalAccount, const int aPeriod ) = 0;
    virtual void endVisitNationalAccount( const NationalAccount* aNationalAccount, const int aPeriod ) = 0;

    virtual void startVisitInput( const IInput* aInput, const int aPeriod ) = 0;
    virtual void endVisitInput( const IInput* aInput, const int aPeriod ) = 0;

    virtual void startVisitProductionInput( const ProductionInput* aProductionInput, const int aPeriod ) = 0;
    virtual void endVisitProductionInput( const ProductionInput* aProductionInput, const int aPeriod ) = 0;

    virtual void startVisitDemandInput( const DemandInput* aDemandInput, const int aPeriod ) = 0;
    virtual void endVisitDemandInput( const DemandInput* aDemandInput, const int aPeriod ) = 0;

    virtual void startVisitExpenditure( const Expenditure* aExpenditure, const int aPeriod ) = 0;
    virtual void endVisitExpenditure( const Expenditure* aExpenditure, const int aPeriod ) = 0;

    virtual void startVisitMiniCAMInput( const MiniCAMInput* aInput, const int aPeriod ) = 0;
    virtual void endVisitMiniCAMInput( const MiniCAMInput* aInput, const int aPeriod ) = 0;

    virtual void startVisitSGMInput( const SGMInput* aSGMInput, const int aPeriod ) = 0;
    virtual void endVisitSGMInput( const SGMInput* aSGMInput, const int aPeriod ) = 0;
    
    virtual void startVisitOutput( const IOutput* aOutput, const int aPeriod ) = 0;
    virtual void endVisitOutput( const IOutput* aOutput, const int aPeriod ) = 0;

    virtual void startVisitGHG( const AGHG* aGHG, const int aPeriod ) = 0;
    virtual void endVisitGHG( const AGHG* aGHG, const int aPeriod ) = 0;
    
    virtual void startVisitOutputMetaData( const OutputMetaData* aOutputMetaData, const int aPeriod ) = 0;
    virtual void endVisitOutputMetaData( const OutputMetaData* aOutputMetaData, const int aPeriod ) = 0;
    
    virtual void startVisitMarketplace( const Marketplace* aMarketplace, const int aPeriod ) = 0;
    virtual void endVisitMarketplace( const Marketplace* aMarketplace, const int aPeriod ) = 0;
    
    virtual void startVisitMarket( const Market* aMarket, const int aPeriod ) = 0;
    virtual void endVisitMarket( const Market* aMarket, const int aPeriod ) = 0;

    virtual void startVisitGDP( const GDP* aGDP, const int aPeriod ) = 0;
    virtual void endVisitGDP( const GDP* aGDP, const int aPeriod ) = 0;
    
    virtual void startVisitLandNode( const LandNode* aLandNode, const int aPeriod ) = 0;
    virtual void endVisitLandNode( const LandNode* aLandNode, const int aPeriod ) = 0;

    virtual void startVisitLandLeaf( const LandLeaf* aLandLeaf, const int aPeriod ) = 0;
    virtual void endVisitLandLeaf( const LandLeaf* aLandLeaf, const int aPeriod ) = 0;
    
    virtual void startVisitLandUseHistory( const LandUseHistory* aLandUseHistory, const int aPeriod ) = 0;
    virtual void endVisitLandUseHistory( const LandUseHistory* aLandUseHistory, const int aPeriod ) = 0;

    virtual void startVisitCarbonCalc( const ICarbonCalc* aCarbonCalc, const int aPeriod ) = 0;
    virtual void endVisitCarbonCalc( const ICarbonCalc* aCarbonCalc, const int aPeriod ) = 0;

    virtual void startVisitClimateModel( const IClimateModel* aClimateModel, const int aPeriod ) = 0;
    virtual void endVisitClimateModel( const IClimateModel* aClimateModel, const int aPeriod ) = 0;

    virtual void startVisitNodeInput( const NodeInput* aNodeInput, const int aPeriod ) = 0;
    virtual void endVisitNodeInput( const NodeInput* aNodeInput, const int aPeriod ) = 0;

    // Following are the derived class accepts
    virtual void startVisitTranTechnology( const TranTechnology* aTranTechnology, const int aPeriod ) = 0;
    virtual void endVisitTranTechnology( const TranTechnology* aTranTechnology, const int aPeriod ) = 0;
    
    virtual void startVisitBuildingGenericDmdTechnology( const BuildingGenericDmdTechnology* aBuildingTechnology, const int aPeriod ) = 0;
    virtual void endVisitBuildingGenericDmdTechnology( const BuildingGenericDmdTechnology* aBuildingTechnology, const int aPeriod ) = 0;
};

IVisitor::~IVisitor(){
}
#endif // _IVISITOR_H_
