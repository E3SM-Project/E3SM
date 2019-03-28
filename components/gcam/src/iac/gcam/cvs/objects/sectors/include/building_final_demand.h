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

#ifndef _BUILDLING_DMD_SECTOR_H_
#define _BUILDLING_DMD_SECTOR_H_
#if defined(_MSC_VER)
#pragma once
#endif

/*! 
 * \file building_final_demand.h
 * \ingroup Objects
 * \brief BuildingFinalDemand class header file.
 * \author Steve Smith
 */

#include <vector>
#include <xercesc/dom/DOMNode.hpp>
#include "sectors/include/energy_final_demand.h"
#include "util/base/include/value.h"

// Forward declarations
class GDP;
class Demographics;
class IInfo;

/*! 
 * \ingroup Objects
 * \brief A class which defines a single building demand sector.
 * \details The building demand sector calculates the demand for building
 *          services in terms of square feet. The actual building service is
 *          supplied by a number of separate supply sectors. No energy is used
 *          directly by this sector or its technologies. The building model first
 *          calculates a per capita demand for the service, and then converts
 *          that to a total demand using the regional population. Service demand
 *          per capita is calculated as:
 *
 *          \f[ D_{per capita} = \gamma * \frac{i^{\alpha}}{1 + ({\frac{i}{Xc}})^\beta} \f]
 *
 *          where \n
 *          \f$D\f$ is the per capita demand for square footage. \n
 *          \f$\gamma\f$ is the calibrated base scaler. \n
 *          \f$i\f$ is the MER per capita GDP. \n
 *          \f$\alpha\f$ is the income elasticity. \n
 *          \f$\beta\f$ is the saturation elasticity. \n
 *          \f$Xc\f$ is the saturation control point. \n
 *
 *          <b>XML specification for BuildingFinalDemand</b>
 *          - XML name: \c building-demand-sector
 *          - Contained by: Region
 *          - Parsing inherited from class: EnergyFinalDemand
 *          - Attributes:
 *              - \c name EnergyFinalDemand::mName
 *          - Elements:
 *              - \c base-service BuildingFinalDemand::mBaseService
 *              - \c saturation-elasticity BuildingFinalDemand::mSaturationElasticity
 *              - \c saturation-point BuildingFinalDemand::mSaturationPoint
 *
 * \author Steve Smith, Josh Lurz
 */
class BuildingFinalDemand: public EnergyFinalDemand
{
public:
    static const std::string& getXMLNameStatic();

    BuildingFinalDemand();

    virtual ~BuildingFinalDemand();

    virtual void completeInit( const std::string& aRegionName,
                               const IInfo* aRegionInfo );

    virtual void initCalc( const std::string& aRegionName,
                           const GDP* aGDP,
                           const Demographic* aDemographics,
                           const int aPeriod );
                                       
protected:
    bool XMLDerivedClassParse( const std::string& nodeName, const xercesc::DOMNode* curr ); 
    void toInputXMLDerived( std::ostream& out, Tabs* tabs ) const;
    void toDebugXMLDerived( const int period, std::ostream& out, Tabs* tabs ) const;
    const std::string& getXMLName() const;
    virtual double calcMacroScaler( const std::string& aRegionName,
                                    const Demographic* aDemographics,
                                    const GDP* aGDP,
                                    const int aPeriod ) const;

    //! Elasticity for saturation.
    Value mSaturationElasticity;

    //! Saturation point of the demand function.
    Value mSaturationPoint;
};

#endif // _BUILDLING_DMD_SECTOR_H_
