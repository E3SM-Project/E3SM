#ifndef _AG_PRODUCTION_TECHNOLOGY_H_
#define _AG_PRODUCTION_TECHNOLOGY_H_
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
* \file ag_production_technology.h
* \ingroup Objects
* \brief The AgProductionTechnology class header file.
* Main difference between agProductionTechnology class and default 
* technology class is that production is determined here bottom up
* by multiplying land qunatities by yields rather than any logit
* sharing among technologies.  Land sharing is done in the land allocator
*
* \author Marshall Wise, Kate Calvin
*/

#include <xercesc/dom/DOMNode.hpp>
#include "technologies/include/technology.h"

// Forward declaration
class Tabs;
class DependencyFinder;
class ILandAllocator;
class ALandAllocatorItem;

/*!
* \ingroup Objects
* \brief A technology which supplies ag products.
* \details TODO
* \author James Blackwood
*/

class AgProductionTechnology : public Technology {
    friend class XMLDBOutputter;
    friend class SetCarbonDensity;
public:
    AgProductionTechnology( const std::string& aName,
                              const int aYear );
    AgProductionTechnology( const AgProductionTechnology& aAgTech );
    ~AgProductionTechnology();
    static const std::string& getXMLNameStatic();
    AgProductionTechnology* clone() const;    
    
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

    virtual double calcShare( const std::string& aRegionName,
                              const std::string& aSectorName,
                              const GDP* aGDP,
                              const double aLogitExp,
                              const int aPeriod ) const; 
    
    virtual void production( const std::string& aRegionName,
                             const std::string& aSectorName, 
                             double aVariableDemand,
                             double aFixedOutputScaleFactor,
                             const GDP* aGDP,
                             const int aPeriod );

    virtual void calcCost( const std::string& aRegionName,
                           const std::string& aSectorName,
                           const int aPeriod );

    virtual double getNonEnergyCost( const std::string& aRegionName,
                                     const int aPeriod ) const;
    
    virtual void doInterpolations( const Technology* aPrevTech, const Technology* aNextTech );
    
    virtual Value getParsedShareWeight() const;
protected:

    double mNonLandVariableCost; //!< The non-land cost of producing a unit of product
    double mNonLandCostTechChange; //!< Annual percent reduction in non-land variableCost
    double mYield; //!< optional input of yield used only for sectors with no current production
	// but maybe also use for calibrated yield?

    double mAgProdChange;  //!< the technological change factor

    //! Measure of multiple cropping
    double mHarvestsPerYear;

    /*!
     * \brief Yield adjustment due to change in carbon densities.
     * \warning This change was made for iESM and this value is currently
     *          not linked to carbon densities but is mananged by the iESM
     *          visitor SetCarbonDensity.
     */
    double mYieldAdjForDensity;
 
    ILandAllocator* mLandAllocator;
    
    //! Weak pointer to the land leaf which corresponds to this technology
    //! used to save time finding it over and over
    ALandAllocatorItem* mProductLeaf;

    virtual void toInputXMLDerived( std::ostream& out, Tabs* tabs ) const;
    virtual void toDebugXMLDerived( const int period, std::ostream& out, Tabs* tabs ) const;
    virtual bool XMLDerivedClassParse( const std::string& nodeName, const xercesc::DOMNode* curr );
    virtual void acceptDerived( IVisitor* aVisitor, const int aPeriod ) const;    
    virtual const std::string& getXMLName() const;
    
    virtual double calcProfitRate( const std::string& aRegionName,
                                   const std::string& aProductName,
                                   const int aPeriod ) const;

    double calcSupply( const std::string& aRegionName,
                       const std::string& aProductName,
                       const int aPeriod ) const;

    void setCalYields( const std::string& aRegionName);
};

#endif // _AG_PRODUCTION_TECHNOLOGY_H_

