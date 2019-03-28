#ifndef _UNMANAGED_LAND_TECHNOLOGY_H_
#define _UNMANAGED_LAND_TECHNOLOGY_H_
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
* \file unmanaged_land_technology.h
* \ingroup Objects
* \brief The UnmanagedLandTechnology class header file.
* \author Steve Smith, Kate Calvin
*/  

#include <xercesc/dom/DOMNode.hpp>

#include "technologies/include/ag_production_technology.h"
#include "util/base/include/value.h"

// Forward declaration
class Tabs;
class DependencyFinder;

/*! 
* \ingroup Objects
* \brief Unmanaged land technology.
* \details 
* This class allows non-CO2 ghg emissions from land leaves. "Input" GHG objects will be driven
* by land area while "Output" GHG objects will be driven by land-area change
* (e.g. deforestation).  Note that this class has NO primary output and is here merely
* to allow calcultion of byproducts of land processes.  Values setin the primary
* output object are just for the sake of emissions calculations and will be removed
* immediately after those calculations are done.
* \todo Clean up this code and remove UnmanagedLandTechnology 
*       as a friend class in RenewableInput 
* \author Steve Smith, Kate Calvin
*/

class UnmanagedLandTechnology : public AgProductionTechnology {
public:
    UnmanagedLandTechnology( const std::string& aName,
                                const int aYear );

    ~UnmanagedLandTechnology();
    UnmanagedLandTechnology( const UnmanagedLandTechnology& aOther );
    static const std::string& getXMLNameStatic();
    const std::string& getXMLName() const;
    UnmanagedLandTechnology* clone() const;

    virtual void completeInit( const std::string& aRegionName,
                               const std::string& aSectorName,
                               const std::string& aSubsectorName,
                               DependencyFinder* aDepFinder,
                               const IInfo* aSubsectorIInfo,
                               ILandAllocator* aLandAllocator );
    
    virtual void initCalc( const std::string& aRegionName,
                           const std::string& aSectorName,
                           const IInfo* aSubsectorIInfo,
                           const Demographic* aDemographics,
                           PreviousPeriodInfo& aPrevPeriodInfo,
                           const int aPeriod );

    virtual double calcShare( const std::string& aRegionName,
                              const std::string& aSectorName,
                              const GDP* aGDP,
                              const double aLogitExp,
                              const int aPeriod ) const; 

    virtual void calcCost( const std::string& aRegionName,
                           const std::string& aSectorName,
                           const int aPeriod );

    virtual void production( const std::string& aRegionName,
                             const std::string& aSectorName, 
                             double aVariableDemand,
                             double aFixedOutputScaleFactor,
                             const GDP* aGDP,
                             const int aPeriod );

    virtual void calcEmissionsAndOutputs( const std::string& aRegionName,
                                  const double aPrimaryOutput,
                                  const GDP* aGDP,
                                  const int aPeriod );

protected:
    bool XMLDerivedClassParse( const std::string& nodeName, const xercesc::DOMNode* curr );

    void initializeInputLocations( const std::string& aRegionName,
                                   const std::string& aSectorName,
                                   const int aPeriod );

    virtual const std::string& getLandInputName( ) const;

private:
    //! Name of leaf or node to use as driver for this technology
    std::string mLandItemName;

    typedef std::vector<IInput*>::iterator InputIterator;

    //! Internal input that contains the amount of land.
    InputIterator mResourceInput;

    virtual void toInputXMLDerived( std::ostream& out, Tabs* tabs ) const;
    virtual void toDebugXMLDerived( const int period, std::ostream& out, Tabs* tabs ) const;
 };

#endif // _UNMANAGED_LAND_TECHNOLOGY_H_

