#ifndef _ITECHNOLOGY_H_
#define _ITECHNOLOGY_H_
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
* \file itechnology.h
* \ingroup Objects
* \brief The technology interface header file.
* \author Pralit Patel
*/

#include <string>
#include <vector>
#include <map>
#include <xercesc/dom/DOMNode.hpp>
#include "util/base/include/istandard_component.h"
#include "util/base/include/value.h"

// Forward declaration
class AGHG;
class GDP;
class DependencyFinder;
class IInfo;
class ICalData;
class ILandAllocator;
class Demographic;
class IOutput;
class IInput;
class Technology;

/*!
* \brief A structure containing information about the same type Technology
*        in the previous period.
*/
struct PreviousPeriodInfo {
    //! The vector of inputs.
    const std::vector<IInput*>* mInputs;

    //! Cumulative Hicks-Neutral technical change information from the previous
    //! period.
    double mCumulativeHicksNeutralTechChange;
    
    //! If info was set by a technology.
    bool mIsFirstTech;
};

/*!
* \ingroup Objects
* \brief This class defines the interface for what it means to be a MiniCAM technology.
*
* This interface is abstract.  All MiniCAM technologies will implement these methods.
*
* \author Pralit Patel
*/
class ITechnology: public IParsedComponent
{
public:
    virtual ITechnology* clone() const = 0;

    inline virtual ~ITechnology() = 0;

    virtual void setYear( const int aNewYear ) = 0;

    virtual bool XMLParse( const xercesc::DOMNode* tempnode ) = 0;
    virtual void toInputXML( std::ostream& out, Tabs* tabs ) const = 0;
    virtual void toDebugXML( const int period, std::ostream& out, Tabs* tabs ) const = 0;
    
    virtual const std::string& getXMLName() const = 0;
    
    virtual void completeInit( const std::string& aRegionName,
                               const std::string& aSectorName,
                               const std::string& aSubsectorName,
                               DependencyFinder* aDepFinder,
                               const IInfo* aSubsectorIInfo,
                               ILandAllocator* aLandAllocator ) = 0;
    
    virtual void initCalc( const std::string& aRegionName,
                           const std::string& aSectorName,
                           const IInfo* aSubsectorInfo,
                           const Demographic* aDemographics,
                           PreviousPeriodInfo& aPrevPeriodInfo,
                           const int aPeriod ) = 0;
    
    virtual void postCalc( const std::string& aRegionName,
                           const int aPeriod ) = 0;

    virtual void production( const std::string& aRegionName,
                             const std::string& aSectorName, 
                             double aVariableDemand,
                             double aFixedOutputScaleFactor,
                             const GDP* aGDP,
                             const int aPeriod ) = 0;

    virtual double calcShare( const std::string& aRegionName,
                              const std::string& aSectorName, 
                              const GDP* aGDP,
                              const double aLogitExp,
                              const int aPeriod ) const = 0;
    
    virtual void calcCost( const std::string& aRegionName,
                           const std::string& aSectorName,
                           const int aPeriod ) = 0;

    virtual double getCost( const int aPeriod ) const = 0;
    
    virtual double getEnergyCost( const std::string& aRegionName,
                                  const std::string& aSectorName,
                                  const int aPeriod ) const = 0;


    virtual double getEnergyInput( const int aPeriod ) const = 0;

    virtual double getCalibrationOutput( const bool aHasRequiredInput,
                                         const std::string& aRequiredInput, 
                                         const int aPeriod ) const = 0;

    virtual bool hasCalibratedValue( const int aPeriod ) const = 0;

    virtual const std::map<std::string,double> getEmissions( const std::string& aGoodName,
                                                             const int aPeriod ) const = 0;

    virtual const std::map<std::string,double> getEmissionsByFuel( const std::string& aGoodName,
                                                                   const int aPeriod ) const = 0;

    virtual const std::string& getName() const = 0;

    virtual void setShareWeight( double shareWeightValue ) = 0;

    virtual bool isOutputFixed( const bool aHasRequiredInput,
                                const std::string& aRequiredInput, 
                                const int aPeriod ) const = 0;

    virtual bool isFixedOutputTechnology( const int aPeriod ) const = 0;

    virtual double getOutput( const int aPeriod ) const = 0;

    virtual double getTotalGHGCost( const std::string& aRegionName, const std::string& aSectorName, 
                            const int aPeriod ) const = 0;

    virtual Value getShareWeight() const = 0;
    virtual Value getParsedShareWeight() const = 0;
    virtual int getNumbGHGs()  const = 0;
    virtual void copyGHGParameters( const AGHG* prevGHG ) = 0;

    virtual const AGHG* getGHGPointer( const std::string& aGHGName ) const = 0;

    virtual const std::vector<std::string> getGHGNames() const = 0;
 
    virtual double getEmissionsByGas( const std::string& aGasName, const int aPeriod ) const = 0;

    virtual double getFixedOutput( const std::string& aRegionName,
                                   const std::string& aSectorName,
                                   const bool aHasRequiredInput,
                                   const std::string& aRequiredInput,
                                   const int aPeriod ) const = 0;
    
    virtual bool isAllCalibrated( const int aPeriod,
                          double aCalAccuracy,
                          const std::string& aRegionName,
                          const std::string& aSectorName,
                          const std::string& aSubsectorName,
                          const bool aPrintWarnings ) const = 0;

    virtual void accept( IVisitor* aVisitor, const int aPeriod ) const = 0;
    virtual void acceptDerived( IVisitor* aVisitor, const int aPeriod ) const = 0;

    virtual const std::map<std::string, double> getFuelMap( const int aPeriod ) const = 0;

    virtual bool isAvailable( const int aPeriod ) const = 0;
    
    virtual double calcFuelPrefElasticity( const int aPeriod ) const = 0;
    
    virtual void doInterpolations( const Technology* aPrevTech, const Technology* aNextTech ) = 0;

    protected:

    virtual double getTotalInputCost( const std::string& aRegionName,
                                    const std::string& aSectorName,
                                    const int aPeriod ) const = 0;
};

// Inline methods
ITechnology::~ITechnology(){
}

#endif // _ITECHNOLOGY_H_
