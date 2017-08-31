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

#ifndef _NON_ENERGY_FINAL_DEMAND_H_
#define _NON_ENERGY_FINAL_DEMAND_H_
#if defined(_MSC_VER)
#pragma once
#endif

/*! 
 * \file non_energy_final_demand.h
 * \ingroup Objects
 * \brief The EnergyFinalDemand class header file.
 * \author Josh Lurz
 */

#include <vector>
#include <xercesc/dom/DOMNode.hpp>

#include "sectors/include/afinal_demand.h"

// Forward declarations
class GDP;

/*! 
 * \ingroup Objects
 * \brief A class which represents a single end use of a non-energy product or
 *        service.
 * \details Non-energy final demands consume a good which is not an energy good
 *          and are not counted towards the total final energy of the region.
 */
class NonEnergyFinalDemand: public AFinalDemand
{
public:
    static const std::string& getXMLNameStatic();

    NonEnergyFinalDemand();

    virtual ~NonEnergyFinalDemand();

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
                                       const int aPeriod ) const;

    virtual void scaleCalibratedValues( const std::string& aFuelName,
                                        const double aScaleValue,
                                        const int aPeriod );

    virtual void csvOutputFile( const std::string& aRegionName ) const;

    virtual void dbOutput( const std::string& aRegionName ) const;

    virtual void accept( IVisitor* aVisitor, const int aPeriod ) const;
protected:
    bool mIsPerCapitaBased; //!< demand equation based on per capita GNP
    
    std::string mName; //!< Name of the final demand

    std::vector<double> mServiceDemands; //!< total end-use sector service is applied.
    
    //! Income elasticity 
    // double mIncomeElasticity;
    std::vector<double> mIncomeElasticity; // TODO: Remove with reformulation.

    //! Price elasticity.
    // double mPriceElasticity;
    std::vector<double> mPriceElasticity; // TODO: Remove with reformulation.

    //! Base year service.
    double mBaseService;
 
    //! Autonomous end-use energy intensity parameter.
    std::vector<double> mAEEI;
    
    //! Calculated total technical change.
    std::vector<double> mTechnicalChange;

    virtual void calcTechChange( const int aPeriod );
    
    virtual double calcDemand( const std::string& aRegionName,
                               const Demographic* aDemographics,
                               const GDP* aGDP,
                               const int aPeriod ) const;
};

#endif // _ENERGY_FINAL_DEMAND_H_

