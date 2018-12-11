#ifndef _LAND_ALLOCATOR_H_
#define _LAND_ALLOCATOR_H_
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
 * \file land_allocator.h
 * \ingroup Objects
 * \brief The LandAllocator class header file.
 * \author James Blackwood, Josh Lurz, Kate Calvin
 */
#include "land_allocator/include/iland_allocator.h"
#include "land_allocator/include/land_node.h"
#include "util/base/include/ivisitable.h"

class IInfo;

/*! 
 * \brief Root of a single land allocation tree.
 * \details The land allocator root contains the root of the land allocation
 *          tree and controls all access from the model into the land allocation
 *          system. This is accomplished by implementing the ILandAllocator
 *          interface, which is the only interface to which Regions have access.
 *          Many methods on this interface are implemented by directly calling
 *          the LandAllocatorNode functions.
 *
 *          <b>XML specification for LandAllocator</b>
 *          - XML name: \c LandAllocatorRoot
 *          - Contained by: Region
 *          - Parsing inherited from class: None
 *          - Attributes:
 *              - \c name ALandAllocatorItem::mName
 *          - Elements:
 *              - \c landAllocation ALandAllocatorItem::mLandAllocation
 */
class LandAllocator : public ILandAllocator,
                          public LandNode {
public:
    LandAllocator();
    virtual ~LandAllocator();
    static const std::string& getXMLNameStatic();

    // IParsable
    virtual bool XMLParse( const xercesc::DOMNode* aNode );
    
    // ILandAllocator methods.
    virtual void toDebugXML( const int aPeriod,
                             std::ostream& aOut,
                             Tabs* aTabs ) const;
    
    virtual void toInputXML( std::ostream& aOut,
                             Tabs* aTabs ) const;
    
    virtual double getLandAllocation( const std::string& aProductName,
                                      const int aPeriod ) const;

    virtual void setCarbonPriceIncreaseRate( const double aCarbonPriceIncreaseRate,
                                        const int aPeriod );

    virtual void setSoilTimeScale( const int aTimeScale );

    virtual void setProfitRate( const std::string& aRegionName,
                                   const std::string& aProductName,
                                   const double aProfitRate,
                                   const int aPeriod );
   
    virtual void calcFinalLandAllocation( const std::string& aRegionName, 
                                          const int aPeriod );

    virtual void completeInit( const std::string& aRegionName, 
                               const IInfo* aRegionInfo );
    
    virtual void initCalc( const std::string& aRegionName,
                           const int aPeriod );
                              
    virtual void postCalc( const std::string& aRegionName,
                           const int aPeriod );
    
    virtual void accept( IVisitor* aVisitor,
                         const int aPeriod ) const;

    // Land allocator node methods.
    virtual void setInitShares( const std::string& aRegionName,
                                const double aLandAllocationAbove,
                                const int aPeriod );

    virtual double calcLandShares( const std::string& aRegionName,
                                   const double aLogitExpAbove,
                                   const int aPeriod );

     virtual void calcLandAllocation( const std::string& aRegionName,
                                     const double aLandAllocationAbove,
                                     const int aYear );

    virtual void calcLUCEmissions( const std::string& aRegionName,
                                   const int aPeriod, const int aEndYear );
                              
    virtual ALandAllocatorItem* findProductLeaf( const std::string& aProductName );
protected:
    virtual const std::string& getXMLName() const;

    virtual bool XMLDerivedClassParse( const std::string& aNodeName,
                                       const xercesc::DOMNode* aCurr );

    virtual void toInputXMLDerived( std::ostream& aOutput,
                                    Tabs* aTabs ) const;
private:
    //! Land allocated in 1000's of hectares
    objects::PeriodVector<double> mLandAllocation;

    //! Rate at which carbon price is expected to increase
    objects::PeriodVector<double> mCarbonPriceIncreaseRate;

    //! Integer storing the soil time scale for a region
    int mSoilTimeScale;                              

    void calibrateLandAllocator( const std::string& aRegionName, const int aPeriod );

    void calculateProfitScalers( const std::string& aRegionName, 
                                const int aPeriod );

    void adjustProfitScalers( const std::string& aRegionName, 
                                const int aPeriod );

    void checkLandArea( const std::string& aRegionName, const int aPeriod );
};

#endif // _LAND_ALLOCATOR_H_
