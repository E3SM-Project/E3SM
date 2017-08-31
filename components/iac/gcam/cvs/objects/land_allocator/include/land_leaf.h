#ifndef _LAND_LEAF_H_
#define _LAND_LEAF_H_
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
 * \file land_leaf.h
 * \ingroup Objects
 * \brief The LandLeaf class header file.
 * \author James Blackwood
 */

#include <xercesc/dom/DOMNode.hpp>
#include "land_allocator/include/aland_allocator_item.h"
#include "util/base/include/ivisitable.h"
#include "util/base/include/time_vector.h"
#include "util/base/include/value.h"

class Tabs;
class ICarbonCalc;
class LandNode;

/*!
 * \brief A LandLeaf is the leaf of a land allocation tree.
 * \details A leaf in the land allocator which represents the land used to
 *          produce a single crop. Land leaves can be separated into two
 *          categories, managed land leaves which are used by farming
 *          technologies, and unmanaged land leaves.
 *
 *          <b>XML specification for LandLeaf</b>
 *          - XML name: Not parsed
 *          - Contained by: LandNode
 *          - Parsing inherited from class: None
 *          - Attributes: None
 *          - Elements: None
 */
class LandLeaf : public ALandAllocatorItem {
    friend class XMLDBOutputter;
public:
    LandLeaf( const ALandAllocatorItem* aParent,
              const std::string& aName );

    virtual ~LandLeaf();

    static const std::string& getXMLNameStatic();

    // Tree Item methods.
    virtual size_t getNumChildren() const;

    virtual const ALandAllocatorItem* getChildAt( const size_t aIndex ) const;

    virtual ALandAllocatorItem* getChildAt( const size_t aIndex );

    virtual void completeInit( const std::string& aRegionName, 
                               const IInfo* aRegionInfo );
        
    virtual void initCalc( const std::string& aRegionName,
                           const int aPeriod );

    virtual void setInitShares( const std::string& aRegionName,
                                const double aLandAllocationAbove,
                                const int aPeriod );

    virtual void calculateProfitScalers( const std::string& aRegionName, 
                                const int aPeriod );

	virtual void setProfitRate( const std::string& aRegionName,
                                   const std::string& aProductName,
                                   const double aProfitRate,
                                   const int aPeriod );

    virtual void setCarbonPriceIncreaseRate( const double aCarbonPriceIncreaseRate, 
                                      const int aPeriod );

    virtual void setSoilTimeScale( const int aTimeScale );

    virtual double calcLandShares( const std::string& aRegionName,
                                   const double aLogitExpAbove,
                                   const int aPeriod );

    virtual void calcLandAllocation( const std::string& aRegionName,
                                     const double aLandAllocationAbove,
                                     const int aPeriod );

    virtual void calcLUCEmissions( const std::string& aRegionName,
                                   const int aPeriod, const int aEndYear );

    virtual double getLandAllocation( const std::string& aProductName,
                                      const int aPeriod ) const;

    virtual double getCalLandAllocation( const LandAllocationType aType,
                                         const int aPeriod ) const;
    virtual double getNewTechProfitScaler( const int aPeriod ) const;
        
    virtual double getLogitExponent( const int aPeriod ) const;

    virtual void setUnmanagedLandProfitRate( const std::string& aRegionName, 
                                             double aAverageProfitRate,
                                             const int aPeriod );

    virtual void calculateCalibrationProfitRate( const std::string& aRegionName, 
                                             double aAverageProfitRate,
                                             double aLogitExponentAbove,
                                             const int aPeriod );

    virtual void adjustProfitScalers( const std::string& aRegionName, 
                                const int aPeriod );

    virtual bool XMLParse( const xercesc::DOMNode* aNode );

    virtual void toInputXML( std::ostream& aOut,
                             Tabs* aTabs ) const;

    virtual void accept( IVisitor* aVisitor,
                         const int aPeriod ) const;

    virtual void acceptDerived( IVisitor* aVisitor,
                         const int aPeriod ) const;

    virtual ICarbonCalc* getCarbonContentCalc() const;
        
    virtual bool isManagedLandLeaf( )  const;

protected:
    //! Land allocated in 1000's of hectares
    objects::PeriodVector<Value> mLandAllocation;

    //! Carbon content and emissions calculator for the leaf.
    std::auto_ptr<ICarbonCalc> mCarbonContentCalc;

    //! Interest rate stored from the region info.
    Value mInterestRate;

    //! Minimum above ground carbon density (used for carbon subsidy and not emissions calculations)
    Value mMinAboveGroundCDensity;

    //! Minimum below ground carbon density (used for carbon subsidy and not emissions calculations)
    Value mMinBelowGroundCDensity;

    //! Expected rate of increase of the carbon price from the region info.
    objects::PeriodVector<Value> mCarbonPriceIncreaseRate;

    //! Container of historical land use.
    std::auto_ptr<LandUseHistory> mLandUseHistory;
    
    objects::PeriodVector<Value> mReadinLandAllocation;

    double getCarbonSubsidy( const std::string& aRegionName,
                           const int aPeriod ) const;

    virtual bool XMLDerivedClassParse( const std::string& aNodeName,
                                       const xercesc::DOMNode* aCurr );

    virtual void toDebugXMLDerived( const int aPeriod,
                                    std::ostream& aOut,
                                    Tabs* aTabs ) const;

    virtual const std::string& getXMLName() const;

    virtual void initLandUseHistory( const std::string& aRegionName );

};

#endif // _LAND_LEAF_H_
