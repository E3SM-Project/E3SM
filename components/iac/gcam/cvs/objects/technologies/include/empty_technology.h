#ifndef _EMPTY_TECHNOLOGY_H_
#define _EMPTY_TECHNOLOGY_H_
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
 * \file empty_technology.h
 * \ingroup Objects
 * \brief The EmptyTechnology class header file.
 * \author Pralit Patel
 */
class IVisitor;

#include "technologies/include/itechnology.h"

/*!
 * \ingroup Objects
 * \brief An implementation of an ITechnology that does not do anything.
 * \details The idea here is this technology can be used as a fill in for years in which
 *          the real technology should not exist.  This will allow all calculations
 *          to remain the same without checking for a missing technology.  Note that
 *          there is only a single instance of this class and thus it can not be
 *          parsed nor should it be written in toInputXML.
 *
 * \author Pralit Patel
 */
class EmptyTechnology: public ITechnology
{
public:
    // Singleton class will only provide a getInstance method and no constructors
    static EmptyTechnology* getInstance();
    
    // ISimpleComponent methods
    virtual bool isSameType( const std::string& aType ) const;
    
    // ITechnology methods
    virtual ITechnology* clone() const;
    
    virtual void setYear( const int aNewYear );
    
    virtual bool XMLParse( const xercesc::DOMNode* tempnode );
    virtual void toInputXML( std::ostream& out, Tabs* tabs ) const;
    virtual void toDebugXML( const int period, std::ostream& out, Tabs* tabs ) const;
    
    virtual const std::string& getXMLName() const;
    
    virtual void completeInit( const std::string& aRegionName,
                              const std::string& aSectorName,
                              const std::string& aSubsectorName,
                              DependencyFinder* aDepFinder,
                              const IInfo* aSubsectorIInfo,
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
    
    virtual double calcShare( const std::string& aRegionName,
                             const std::string& aSectorName, 
                             const GDP* aGDP,
                             const double aLogitExp,
                             const int aPeriod ) const;
    
    virtual void calcCost( const std::string& aRegionName,
                          const std::string& aSectorName,
                          const int aPeriod );
    
    virtual double getCost( const int aPeriod ) const;
    
    virtual double getEnergyCost( const std::string& aRegionName,
                                 const std::string& aSectorName,
                                 const int aPeriod ) const;
    
    
    virtual double getEnergyInput( const int aPeriod ) const;
    
    virtual double getCalibrationOutput( const bool aHasRequiredInput,
                                        const std::string& aRequiredInput, 
                                        const int aPeriod ) const;
    
    virtual bool hasCalibratedValue( const int aPeriod ) const;
    
    virtual const std::map<std::string,double> getEmissions( const std::string& aGoodName,
                                                            const int aPeriod ) const;
    
    virtual const std::map<std::string,double> getEmissionsByFuel( const std::string& aGoodName,
                                                                  const int aPeriod ) const;
    
    virtual const std::string& getName() const;
    
    virtual void setShareWeight( double shareWeightValue );
    
    virtual bool isOutputFixed( const bool aHasRequiredInput,
                               const std::string& aRequiredInput, 
                               const int aPeriod ) const;
    
    virtual bool isFixedOutputTechnology( const int aPeriod ) const;
    
    virtual double getOutput( const int aPeriod ) const;
    
    virtual double getTotalGHGCost( const std::string& aRegionName, const std::string& aSectorName, 
                                   const int aPeriod ) const;
    
    virtual Value getShareWeight() const;
    virtual Value getParsedShareWeight() const;
    virtual int getNumbGHGs()  const;
    virtual void copyGHGParameters( const AGHG* prevGHG );
    
    virtual const AGHG* getGHGPointer( const std::string& aGHGName ) const;
    
    virtual const std::vector<std::string> getGHGNames() const;
    
    virtual double getEmissionsByGas( const std::string& aGasName, const int aPeriod ) const;
    
    virtual double getFixedOutput( const std::string& aRegionName,
                                  const std::string& aSectorName,
                                  const bool aHasRequiredInput,
                                  const std::string& aRequiredInput,
                                  const int aPeriod ) const;
    
    virtual bool isAllCalibrated( const int aPeriod,
                                 double aCalAccuracy,
                                 const std::string& aRegionName,
                                 const std::string& aSectorName,
                                 const std::string& aSubsectorName,
                                 const bool aPrintWarnings ) const;
    
    virtual void accept( IVisitor* aVisitor, const int aPeriod ) const;
    virtual void acceptDerived( IVisitor* aVisitor, const int aPeriod ) const;
    
    virtual const std::map<std::string, double> getFuelMap( const int aPeriod ) const;
    
    virtual bool isAvailable( const int aPeriod ) const;
    
    virtual double calcFuelPrefElasticity( const int aPeriod ) const;
    
    virtual void doInterpolations( const Technology* aPrevTech, const Technology* aNextTech );
    
protected:
    
    virtual double getTotalInputCost( const std::string& aRegionName,
                                     const std::string& aSectorName,
                                     const int aPeriod ) const;
    
private:
    // private constructors to enforce a single instance
    EmptyTechnology() {}
    ~EmptyTechnology() {}
    
    // intentionally undefined
    EmptyTechnology( const EmptyTechnology& aEmptyTechnology );
    EmptyTechnology& operator=( const EmptyTechnology& aEmptyTechnology );
};

#endif // _EMPTY_TECHNOLOGY_H_
