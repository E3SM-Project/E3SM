#ifndef _LAND_NODE_H_
#define _LAND_NODE_H_
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
 * \file land_node.h
 * \ingroup Objects
 * \brief The LandNode class header file.
 * \author James Blackwood
 */

#include <vector>
#include <memory>
#include <xercesc/dom/DOMNode.hpp>
#include "land_allocator/include/aland_allocator_item.h"
#include "util/base/include/value.h"

// Forward declarations
class LandUseHistory;
class ICarbonCalc;
/*!
 * \brief A node in the land allocation tree.
 * \details A land allocator node represents a type of land available for
 *          producing crops. It does not itself produce anything, but may have
 *          leaves below it which produce products. The land node may also
 *          contain other land nodes, or land types. Land nodes are always
 *          read-in and are never created dynamically.
 *
 *          <b>XML specification for LandNode</b>
 *          - XML name: -c LandAllocatorNode
 *          - Contained by: LandAllocatorRoot
 *          - Parsing inherited from class: None
 *          - Attributes:
 *              - \c name ALandAllocatorItem::mName
 *          - Elements:
 *              - \c LandNode LandNode::children
 *              - \c UnmanagedLandLeaf LandNode::children
 *              - \c land-use-history LandNode::mLandUseHistory
 */
class LandNode : public ALandAllocatorItem {
public:
    explicit LandNode( const ALandAllocatorItem* aParent );

    virtual ~LandNode();

    // Tree Item methods.
    virtual size_t getNumChildren() const;

    virtual const ALandAllocatorItem* getChildAt( const size_t aIndex ) const;

    virtual ALandAllocatorItem* getChildAt( const size_t aIndex );

    static const std::string& getXMLNameStatic();
    
    virtual void completeInit( const std::string& aRegionName, 
                               const IInfo* aRegionInfo );
    
    virtual void initCalc( const std::string& aRegionName,
                           const int aPeriod );

    virtual void setInitShares( const std::string& aRegionName,
                                const double aLandAllocationAbove,
                                const int aPeriod );

    virtual void calculateProfitScalers( const std::string& aRegionName, 
                                const int aPeriod );

    virtual void adjustProfitScalers( const std::string& aRegionName, 
                                const int aPeriod );

    virtual void setProfitRate( const std::string& aRegionName,
                                   const std::string& aProductName,
                                   const double aProfitRate,
                                   const int aPeriod );
 
    virtual void setCarbonPriceIncreaseRate( const double aCarbonPriceIncreaseRate, 
                                      const int aPeriod );

    /*!
     * \brief Set the number of years needed to for soil carbons emissions/uptake
     * \details This method sets the soil time scale into the carbon calculator
     *          for each land leaf.
     * \param aTimeScale soil time scale (in years)
     * \author Kate Calvin
     */
    virtual void setSoilTimeScale( const int aTimeScale );

    virtual double calcLandShares( const std::string& aRegionName,
                                   const double aLogitExpAbove,
                                   const int aPeriod );

    virtual void calcLandAllocation( const std::string& aRegionName,
                                     const double aLandAllocationAbove,
                                     const int aPeriod );
    
    virtual void calcLUCEmissions( const std::string& aRegionName,
                                   const int aYear, const int aEndYear );
    
    virtual double getLandAllocation( const std::string& aProductName,
                                      const int aPeriod ) const;

    virtual double getCalLandAllocation( const LandAllocationType aType,
                                         const int aPeriod ) const;

    virtual LandUseHistory* getLandUseHistory();
        
    virtual double getLogitExponent( const int aPeriod ) const;

    virtual double getNewTechProfitScaler( const int aPeriod ) const;
    
    virtual void setUnmanagedLandProfitRate( const std::string& aRegionName, 
                                             double aAverageProfitRate,
                                             const int aPeriod );
           
    virtual bool isManagedLandLeaf( )  const;

    virtual void calculateCalibrationProfitRate( const std::string& aRegionName,
                                             double aAverageProfitRate,
                                             double aLogitExponentAbove,
                                             const int aPeriod );

    virtual void accept( IVisitor* aVisitor, 
                         const int aPeriod ) const;

    virtual void toInputXML( std::ostream& out, 
                             Tabs* tabs ) const;

    virtual bool XMLParse( const xercesc::DOMNode* aNode );

protected:
    virtual bool XMLDerivedClassParse( const std::string& nodeName, 
                                       const xercesc::DOMNode* curr );

    virtual void toDebugXMLDerived( const int period, 
                                    std::ostream& out, 
                                    Tabs* tabs ) const;

    virtual void toInputXMLDerived( std::ostream& aOutput, 
                                    Tabs* aTabs ) const;

    virtual const std::string& getXMLName() const;
    
    ALandAllocatorItem* findChild( const std::string& aName,
                                   const LandAllocatorItemType aType );
    
    const ALandAllocatorItem* findChild( const std::string& aName,
                                         const LandAllocatorItemType aType ) const;

    //! Land allocated -- used for conceptual roots
    objects::PeriodVector<double> mLandAllocation;
        
    //! Logit exponent -- should be positive since we are sharing on profit
    objects::PeriodVector<double> mLogitExponent;

    //! Share Profit scaler for new technologies in this node
    objects::PeriodVector<double> mNewTechProfitScaler;
    
    //! Numerator that determines share for new technologies IF the right profit conditions hold
	//! Share will equal ( mGhostShareNumerator / ( 1 + mGhostShareNumerator ) ) if and only if
	//! the profit of the new technology is equal to the profit of the dominant technology in 
	//! the base year, and all other profits stay the same.
    objects::PeriodVector<double> mGhostShareNumerator;
    
    //! Double storing the average price of land in a region or subregion
    double mUnManagedLandValue;

    //! Boolean indicating that scalers in this node should be adjusted for new technologies
    // TODO: we may want this boolean at the LandLeaf level, but it will need to be used in LandNode
    bool mAdjustScalersForNewTech;

    //! List of the children of this land node located below it in the land
    //! allocation tree.
    std::vector<ALandAllocatorItem*> mChildren;

    //! Container of historical land use.
    std::auto_ptr<LandUseHistory> mLandUseHistory;
};

#endif // _LAND_NODE_H_
