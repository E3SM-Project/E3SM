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

#ifndef _ENERGY_FINAL_DEMAND_H_
#define _ENERGY_FINAL_DEMAND_H_
#if defined(_MSC_VER)
#pragma once
#endif

/*! 
 * \file energy_final_demand.h
 * \ingroup Objects
 * \brief The EnergyFinalDemand class header file.
 * \author Josh Lurz
 */

#include <vector>
#include <xercesc/dom/DOMNode.hpp>
#include "sectors/include/afinal_demand.h"
#include "util/base/include/value.h"
#include "util/base/include/time_vector.h"

// Forward declarations
class GDP;
class Demographic;

/*! 
 * \ingroup Objects
 * \brief A class which represents a single end use of an energy product or
 *        service.
 * \details Energy final demands consume an energy derived good and are counted
 *          towards the total final energy of the region.
 */

class EnergyFinalDemand: public AFinalDemand
{
    friend class XMLDBOutputter;
    friend class EnergyBalanceTable; // TODO: currently only to get mServiceDemands

public:
    static const std::string& getXMLNameStatic();

    EnergyFinalDemand();

    virtual ~EnergyFinalDemand();

    virtual bool XMLParse( const xercesc::DOMNode* aNode );

    virtual void toInputXML( std::ostream& aOut,
                             Tabs* aTabs ) const;
    
    virtual void toDebugXML( const int aPeriod,
                             std::ostream& aOut,
                             Tabs* aTabs ) const;
    
    virtual const std::string& getName() const;
    
    virtual void completeInit( const std::string& aRegionName,
                               const IInfo* aRegionInfo );

    virtual void initCalc( const std::string& aRegionName,
                           const GDP* aGDP,
                           const Demographic* aDemographics,
                           const int aPeriod );

    virtual void setFinalDemand( const std::string& aRegionName,
                                 const Demographic* aDemographics,
                                 const GDP* aGDP,
                                 const int aPeriod );

    virtual double getWeightedEnergyPrice( const std::string& aRegionName,
                                           const int aPeriod ) const;

    virtual void tabulateFixedDemands( const std::string& aRegionName,
                                       const Demographic* aDemographic,
                                       const GDP* aGDP,
                                       const int aPeriod );

    virtual void csvOutputFile( const std::string& aRegionName ) const;

    virtual void dbOutput( const std::string& aRegionName ) const;

    virtual void accept( IVisitor* aVisitor, const int aPeriod ) const;
protected:

    // TODO: Move all these.
    class IDemandFunction {
    public:
        // TODO: Remove this function once construction is cleanly implemented.
        virtual bool isPerCapitaBased() const = 0;

        virtual double calcDemand( const Demographic* aDemographics,
                                   const GDP* aGDP,
                                   const double aPriceElasticity,
                                   const double aIncomeElasticity,
                                   const double aPriceRatio,
                                   const int aPeriod ) const = 0;
    };

    class PerCapitaGDPDemandFunction: public IDemandFunction {
    public:
        virtual bool isPerCapitaBased() const {
            return true;
        }

        virtual double calcDemand( const Demographic* aDemographics,
                                   const GDP* aGDP,
                                   const double aPriceElasticity,
                                   const double aIncomeElasticity,
                                   const double aPriceRatio,
                                   const int aPeriod ) const;
    };

    class TotalGDPDemandFunction: public IDemandFunction {
    public:
        virtual bool isPerCapitaBased() const {
            return false;
        }

        virtual double calcDemand( const Demographic* aDemographics,
                                   const GDP* aGDP,
                                   const double aPriceElasticity,
                                   const double aIncomeElasticity,
                                   const double aPriceRatio,
                                   const int aPeriod ) const;
    };

    class FinalEnergyConsumer {
    public:
        static const std::string& getXMLNameStatic();

        static double noCalibrationValue();

        FinalEnergyConsumer( const std::string& aFinalDemandName );

        bool XMLParse( const xercesc::DOMNode* aNode );

        void completeInit( const std::string& aRegionName,
                           const std::string& aFinalDemandName );

        double getCalibratedFinalEnergy( const int aPeriod ) const;

        void updateAEEI( const std::string& aRegionName,
                         const int aPeriod );

        double calcTechChange( const int aPeriod ) const;

        void toInputXML( std::ostream& aOut,
                         Tabs* aTabs ) const;
    
        void toDebugXML( const int aPeriod,
                         std::ostream& aOut,
                         Tabs* aTabs ) const;
    private:
        //! Name of the TFE market.
        std::string mTFEMarketName;

        //! Autonomous end-use energy intensity parameter.
        objects::PeriodVector<Value> mAEEI;

        //! Final energy to calibrate to.
        objects::PeriodVector<Value> mCalFinalEnergy;
    };
    
    //! Name of the final demand and the good it consumes.
    std::string mName;
    
    //! Total end-use sector service after technical change is applied.
    std::vector<double> mServiceDemands;

    //! Income elasticity 
    std::vector<Value> mIncomeElasticity;

    //! Price elasticity.
    std::vector<Value> mPriceElasticity;

    //! Service demand without technical change applied.
    std::vector<double> mPreTechChangeServiceDemand;

    //! Per capita service for each period to which to calibrate.
    std::vector<Value> mBaseService;

    //! Demand function used to calculate unscaled demand.
    std::auto_ptr<IDemandFunction> mDemandFunction;

    //! Object responsible for consuming final energy.
    std::auto_ptr<FinalEnergyConsumer> mFinalEnergyConsumer;
    
    virtual double calcFinalDemand( const std::string& aRegionName,
                                    const Demographic* aDemographics,
                                    const GDP* aGDP,
                                    const int aPeriod );

    virtual double calcMacroScaler( const std::string& aRegionName,
                                    const Demographic* aDemographics,
                                    const GDP* aGDP,
                                    const int aPeriod ) const;

    // Methods for deriving from EnergyFinalDemand.
    virtual bool XMLDerivedClassParse( const std::string& nodeName, const xercesc::DOMNode* curr );
    virtual void toInputXMLDerived( std::ostream& out, Tabs* tabs ) const;
    virtual void toDebugXMLDerived( const int period, std::ostream& out, Tabs* tabs ) const;
    virtual const std::string& getXMLName() const;
private:    
    void acceptDerived( IVisitor* aVisitor, const int aPeriod ) const;
};

#endif // _ENERGY_FINAL_DEMAND_H_

