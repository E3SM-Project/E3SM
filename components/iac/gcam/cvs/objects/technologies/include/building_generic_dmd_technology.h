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

#ifndef _BUILDING_GENERIC_TECHNOLOGY_H_
#define _BUILDING_GENERIC_TECHNOLOGY_H_
#if defined(_MSC_VER)
#pragma once
#endif

/*! 
 * \file building_generic_dmd_technology.h
 * \ingroup Objects
 * \brief BuildingGenericDmdTechnology header file.
 * \author Steve Smith, Josh Lurz
 */

#include <xercesc/dom/DOMNode.hpp>

#include "technologies/include/technology.h"

// Forward declarations
class GDP;
class IInfo;

/*! 
 * \ingroup Objects
 * \brief This building technology class calculates demand for building energy
 *        services.
 * \details Building demand technology objects, act differently than normal
 *          technology objects in that they each generate a demand for a specific
 *          building service (heating, cooling, lighting, etc.), which is then
 *          provided by a supply sector. These technologies do not consume fuels
 *          or generate GHG emissions. These come from the supply sectors.
 *
 *          <b>XML specification for BuildingGenericDmdTechnology</b>
 *          - XML name: \c building-demand-technology
 *          - Contained by: Subsector
 *          - Parsing inherited from class: Technology
 *          - Attributes:
 *              - \c name Technology::mName
 *          - Elements:
 *              - \c aveInsulation BuildingGenericDmdTechnology::mAveInsulation
 *              - \c floorToSurfaceArea BuildingGenericDmdTechnology::mFloorToSurfaceArea
 *
 * \author Steve Smith, Josh Lurz
 */

class BuildingGenericDmdTechnology : public Technology
{
    friend class CalibrateShareWeightVisitor;
public:
    static const std::string& getXMLNameStatic();
	static const std::string& getInternalGainsInfoName();

    BuildingGenericDmdTechnology( const std::string& aName, const int aYear );

    BuildingGenericDmdTechnology( const BuildingGenericDmdTechnology& aOther );

    virtual BuildingGenericDmdTechnology* clone() const;
    virtual ~BuildingGenericDmdTechnology();
    virtual const std::string& getXMLName() const;

    virtual void completeInit( const std::string& aRegionName,
                               const std::string& aSectorName,
                               const std::string& aSubsectorName,
                               DependencyFinder* aDepFinder,
                               const IInfo* aSubsectorInfo,
                               ILandAllocator* aLandAllocator );

    virtual void initCalc( const std::string& aRegionName,
                           const std::string& aSectorName,
                           const IInfo* aSubsectorInfo,
                           const Demographic* aDemographics,
                           PreviousPeriodInfo& aPrevPeriodInfo,
                           const int aPeriod );

    virtual void postCalc( const std::string& aRegionName,
                           const int aPeriod );	
	
    virtual void production( const std::string& aRegionName,
                             const std::string& aSectorName, 
		                     double aVariableDemand,
                             double aFixedOutputScaleFactor,
                             const GDP* aGDP,
                             const int aPeriod );
    
    virtual void calcCost( const std::string& aRegionName,
                           const std::string& aSectorName,
		                   const int aPeriod );
	
    virtual double calcShare( const std::string& aRegionName,
                              const std::string& aSectorName, 
		                      const GDP* aGDP,
                              const double aLogitExp,
                              const int aPeriod ) const;
    
    virtual void accept( IVisitor* aVisitor, const int aPeriod ) const;
    
    virtual void doInterpolations( const Technology* aPrevTech, const Technology* aNextTech );
protected:
    virtual bool XMLDerivedClassParse( const std::string& nodeName, const xercesc::DOMNode* curr );
    void toInputXMLDerived( std::ostream& out, Tabs* tabs ) const;
    void toDebugXMLDerived( const int period, std::ostream& out, Tabs* tabs ) const;
	
    double getTotalInputCost( const std::string& aRegionName,
                              const std::string& aSectorName,
							  const int aPeriod ) const;


    virtual const IInfo* getTechInfo() const;

    virtual const IFunction* getProductionFunction() const;

    void checkCoefficients( const std::string& aSectorName,
                            const int aPeriod ) const;

    //! Average insulation value (J/s-m^2) for this building type
    double mAveInsulation;

    //! Conversion from floor space to surface area for this building type
    double mFloorToSurfaceArea;

    //! Technology info object.
    std::auto_ptr<IInfo> mInfo;
};
#endif // _BUILDING_GENERIC_TECHNOLOGY_H_
