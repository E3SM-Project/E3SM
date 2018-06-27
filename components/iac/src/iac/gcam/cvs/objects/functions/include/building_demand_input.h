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

#ifndef _BUILDING_DEMAND_INPUT_H_
#define _BUILDING_DEMAND_INPUT_H_
#if defined(_MSC_VER)
#pragma once
#endif

/*! 
 * \file building_demand_input.h
 * \ingroup Objects
 * \brief BuildingDemandInput class header file.
 * \author Josh Lurz
 */

#include <string>
#include <xercesc/dom/DOMNode.hpp>
#include "functions/include/energy_input.h"
#include "util/base/include/value.h"
#include <vector>
#include <memory>

class Tabs;
class ICoefficient;

/*! 
 * \ingroup Objects
 * \brief Defines a single input to a building demand function.
 * \details Building demand inputs are used to determine the output of a
 *          building demand function. The building demand function determine the
 *          output of individual service demands based on an aggregate demand.
 *          The building input determines its output coefficient using various
 *          read in parameters such as heating degree days, average insulation
 *          of the building shell, etc. for floorspace. The building input
 *          coefficient must be calibrated in the base period.
 *
 *          <b>XML specification for BuildingDemandInput</b>
 *          - XML name: \c building-demand-input
 *          - Contained by: BuildingGenericDemandTechnology
 *          - Parsing inherited from class: MiniCAMInput
 *          - Attributes:
 *              - \c name MiniCAMInput::mName
 *          - Elements:
 *              - \c saturation BuildingDemandInput::mSaturation
 *              - \c type BuildingDemandInput::mType
 *              - \c coef-adjustment BuildingDemandInput::mCoefficientAdjustment
 *              - \c price-elasticity BuildingDemandInput::mPriceElasticity
 *          
 * \author Josh Lurz
 */
class BuildingDemandInput: public MiniCAMInput
{
    friend class InputFactory;
public:
    const static std::string& getXMLNameStatic();

    virtual const std::string& getXMLReportingName() const;

    virtual ~BuildingDemandInput();

    virtual BuildingDemandInput* clone() const;

    virtual void XMLParse( const xercesc::DOMNode* aNode );

    virtual bool isSameType( const std::string& aType ) const;

    virtual void toInputXML( std::ostream& aOut,
                             Tabs* aTabs ) const;

    virtual void toDebugXML( const int aPeriod,
                             std::ostream& aOut,
                             Tabs* aTabs ) const;

    virtual void copyParam( const IInput* aInput,
                            const int aPeriod );

    virtual void completeInit( const std::string& aRegionName,
                               const std::string& aSectorName,
                               const std::string& aSubsectorName,
                               const std::string& aTechName,
                               DependencyFinder* aDependencyFinder,
                               const IInfo* aTechInfo );

    virtual void initCalc( const std::string& aRegionName,
                           const std::string& aSectorName,
                           const bool aIsNewInvestmentPeriod,
                           const bool aIsTrade,
                           const int aPeriod );

    virtual  double getCO2EmissionsCoefficient( const std::string& aGHGName,
                                             const int aPeriod ) const;
    
    virtual double getPhysicalDemand( const int aPeriod ) const;
    
    virtual double getPrice( const std::string& aRegionName,
                             const int aPeriod ) const;

    virtual void setPrice( const std::string& aRegionName,
                           const double aPrice,
                           const int aPeriod );

    virtual void setPhysicalDemand( const double aPhysicalDemand,
                                    const std::string& aRegionName,
                                    const int aPeriod );

    virtual double getCoefficient( const int aPeriod ) const;

    virtual void setCoefficient( const double aCoefficient,
                                 const int aPeriod );

	virtual double getCalibrationQuantity( const int aPeriod ) const;

    virtual double getTechChange( const int aPeriod ) const;

    virtual bool hasTypeFlag( const int aTypeFlag ) const;
    
    virtual double getIncomeElasticity() const;

    virtual double getPriceElasticity() const;

    virtual void copyParamsInto( BuildingDemandInput& aInput,
                                 const int aPeriod ) const;
    
    virtual void doInterpolations( const int aYear, const int aPreviousYear,
                                   const int aNextYear, const IInput* aPreviousInput,
                                   const IInput* aNextInput );
protected:
    BuildingDemandInput();
    
    //! Physical Demand.
    std::vector<Value> mPhysicalDemand;

    /*!
     * \brief Describes possible types of building demand inputs.
     */
    enum BuildingInputType {
        //! Generic demand input.
        eGeneric,
        
        //! Heating demand input.
        eHeating,

        //! Cooling demand input.
        eCooling,

        //! End marker.
        eEnd
    };

    /*!
     * \brief Type of the building demand input.
     * \details This is stored as an enum but read in as a string. The type
     *          string must be one of:
     *          -Generic
     *          -Heating
     *          -Cooling
     */
    BuildingInputType mType;

    //! Structure which contains type dependent information.
    struct TypeInfo {
        //! The building type for which the typeinfo is valid.
        BuildingInputType mType;

        //! The string representation of the type.
        std::string mTypeString;
    };

    //! Penetration level for the service.
    Value mSaturation;

    //! Coefficient which includes average insulation, floor to surface area,
    //! and degree days.
    Value mCoefficient;

    //! Coeffient adjustment factor due to calibration.
    Value mCoefficientAdjustment;

    //! Average insulation value of the Technology.
    Value mAverageInsulation;

    //! Floor to surface area of the Technology.
    Value mFloorToSurfaceArea;

    //! Heating or cooling degree days depending on the type.
    Value mDegreeDays;

    //! Price elasticity.
    Value mPriceElasticity;

    void initializeParameters( const std::string& aSectorName,
                               const std::string& aSubsectorName,
                               const std::string& aTechName,
                               const IInfo* aTechInfo );

    double calcDemandAdjustmentFactor() const;

    static const BuildingDemandInput::TypeInfo& getTypeInfo( const BuildingInputType aType );

    static BuildingInputType convertStringToType( const std::string& aTypeString );

private:
    const static std::string XML_REPORTING_NAME; //!< tag name for reporting xml db 
};

#endif // _BUILDING_DEMAND_INPUT_H_
