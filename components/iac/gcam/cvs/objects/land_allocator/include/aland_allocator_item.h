#ifndef _ALAND_ALLOCATOR_ITEM_H_
#define _ALAND_ALLOCATOR_ITEM_H_
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
 * \file aland_allocator_item.h
 * \ingroup Objects
 * \brief The ALandAllocatorItem class header file.
 * \author James Blackwood, Kate Calvin
 */

#include <vector>
#include <xercesc/dom/DOMNode.hpp>
#include <boost/noncopyable.hpp>

#include "containers/include/tree_item.h"
#include "util/base/include/ivisitable.h"
#include "util/base/include/iparsable.h"
#include "util/base/include/iround_trippable.h"
#include "util/base/include/time_vector.h"
#include "util/base/include/value.h"

// For LandUsageType enum.
#include "land_allocator/include/iland_allocator.h"

// Forward declarations
class IInfo;
class Tabs;
class LandUseHistory;
class LandNode;


/*!
* \brief An enum containing the possible types for items in the land allocator.
* \note This enum is defined outside the ALandAllocatorItem because
*       ALandAllocatorItem is a template class, so the enum would be a
*       dependent type.
*/
enum LandAllocatorItemType {
    /*!
    * \brief Node type.
    */
    eNode,

    /*!
    * \brief Leaf type.
    */
    eLeaf,

    /*!
    * \brief Any type, either node or leaf.
    */
    eAny
};

/*!
 * \brief A single item in the land allocator tree.
 * \details This is the abstract base class of all nodes and leaves in the land
 *          allocation tree, including the root. It inherits from the TreeItem
 *          class so that it can make use of the tree library functions.
 *
 *          <b>XML specification for ALandAllocatorItem</b>
 *          - XML name: Not parsed
 *          - Contained by: LandAllocator
 *          - Parsing inherited from class: None
 *          - Attributes: None
 *          - Elements: None
 */
class ALandAllocatorItem : public TreeItem<ALandAllocatorItem>,
                           public IVisitable,
                           public IParsable,
                           public IRoundTrippable,
                           private boost::noncopyable
{
    friend class XMLDBOutputter;
public:
    typedef TreeItem<ALandAllocatorItem> ParentTreeType;

    explicit ALandAllocatorItem( const ALandAllocatorItem* aParent,
                                 const LandAllocatorItemType aType );
   
    virtual ~ALandAllocatorItem();

    // TreeItem methods
    virtual size_t getNumChildren() const = 0;

    virtual const ALandAllocatorItem* getChildAt( const size_t aIndex ) const = 0;
    
    virtual ALandAllocatorItem* getChildAt( const size_t aIndex ) = 0;
    
    // IParsable
    virtual bool XMLParse( const xercesc::DOMNode* aNode ) = 0;
    
    // IRoundTrippable
    virtual void toInputXML( std::ostream& aOut,
                             Tabs* aTabs ) const = 0;
    
    // IVisitable
    virtual void accept( IVisitor* aVisitor,
                         const int aPeriod ) const = 0;

    /*!
     * \brief Write datamembers to datastream in XML format for debugging
     *        purposes.  
     * \param aPeriod Model time period
     * \param aOut Output file for debugging purposes in XML format
     * \param aTabs Tabs object used to track the number of tabs to print.
     */
    void toDebugXML( const int aPeriod,
                     std::ostream& aOut,
                     Tabs* aTabs ) const;

    /*!
     * \brief Returns the name.
     * \author James Blackwood
     * \return the name of this ALandAllocatorItem
     */
    const std::string& getName() const;

    /*!
     * \brief Complete the initialization of the ALandAllocatorItem.
     * \param aRegionName Region name.
     * \param aInfo Local info object.
     */
    virtual void completeInit( const std::string& aRegionName, 
                               const IInfo* aRegionInfo ) = 0;

    /*!
     * \brief Complete the initialization of the ALandAllocatorItem.
     * \param aRegionName Region name.
     * \param aInfo Local info object.
     */
    virtual void initCalc( const std::string& aRegionName, 
                           const int aPeriod ) {};

    /*!
     * \brief Get the land allocated to a given product.
     * \details Returns the land allocated to a given product in the
     *          given period.
     * \param aProductName Name of the product.
     * \return Land allocation.
     */
    virtual double getLandAllocation( const std::string& aProductName,
                                      const int aPeriod ) const = 0;

    /*!
     * \brief An enumeration of possible land allocation types.
     */
    enum LandAllocationType {
        // Managed land allocation.
        eManaged,

        //! Unmanaged land allocation.
        eUnmanaged,

        //! Any land allocation.
        eAnyLand
    };

    /*!
     * \brief Returns the total land allocated for calibration for this land type.
     * \details Returns all land allocated for this land allocator item of a
     *          given type, including all items below it that was read in for
     *          the calibration years.
     * \param aType The type of land allocation to return: unmanaged, managed,
     *        or either.
     * \param aPeriod Model period.
     * \author Steve Smith
     * \return The total calibration land allocated at or below the item.
     */
    virtual double getCalLandAllocation( const LandAllocationType aType,
                                         const int aPeriod ) const = 0;

    /*!
     * \brief Sets the initial shares and land allocation.
     * \param aRegionName Region name.
     * \param aLandAllocationAbove Land allocation of the node above this item.
     * \param aPeriod Period.
     * \author Kate Calvin
     */
    virtual void setInitShares( const std::string& aRegionName,
                                const double aLandAllocationAbove,
                                const int aPeriod ) = 0;

    /*!
     * \brief Calculates profit scalers
     * \param aRegionName Region name.
     * \param aPeriod Period.
     * \author Kate Calvin
     */
    virtual void calculateProfitScalers( const std::string& aRegionName,
                                const int aPeriod ) = 0;

    /*!
     * \brief Sets the profit rate for a given product.
     * \details Determines the appropriate land leaf and sets the profit rate
     *          for a given period.
     * \param aProductName Name of the product.
     * \param aProfitRate Profit rate of the product.
     * \param aPeriod Model period.
     * \author James Blackwood, Kate Calvin
     */
    virtual void setProfitRate( const std::string& aRegionName,
                                   const std::string& aProductName,
                                   const double aProfitRate,
                                   const int aPeriod ) = 0;

    /*!
     * \brief Set the rate at which the carbon price is expected to increase
     * \details This method sets expectations about the carbon price to be
     *          used in calculating the carbon subsidy on land. Setting the expected rate
     *          of increase of the carbon price to zero implies myopic decision making
     *          or an expectation of flat carbon prices. Setting this rate to a positive 
     *          number implies an expectation that the carbon price will rise exponentially
     *          at the rate specified.
     * \param aCarbonPriceIncreaseRate Expected rate of increase.
     * \param aPeriod Period.
     * \author Kate Calvin
     */
    virtual void setCarbonPriceIncreaseRate( const double aCarbonPriceIncreaseRate, 
                                      const int aPeriod ) = 0;

    /*!
     * \brief Set the number of years needed to for soil carbons emissions/uptake
     * \details This method sets the soil time scale into the carbon calculator
     *          for each land leaf.
     * \param aTimeScale soil time scale (in years)
     * \author Kate Calvin
     */
    virtual void setSoilTimeScale( const int aTimeScale ) = 0;

    /*!
     * \brief This method will calculate the share value for each leaf and node,
     *          and then normalize it.
     * \details If a child is a leaf then itwill call the calcLandShares 
     *          method in LandLeaf, where share is calculated. If a child 
     *          is a node there will be a recursive call to this method. 
     *          The second loop uses the sum of all the shares of the mChildren 
     *          vector and normalizes and overwrites share. Finally, share 
     *          is calculated for this node using the calculated profitRate 
     *          and the logit exponent from one level up. This method uses the
     *          modified logit from the energy system.
     * \param aRegionName Name of the containing region.
     * \param aLogitExpAbove the logit exponent value from the node above this level.
     * \param aPeriod Model period.
     * \return The unnormalized share.
     * \author Kate Calvin
     */
    virtual double calcLandShares( const std::string& aRegionName,
                                   const double aLogitExpAbove,
                                   const int aPeriod ) = 0;

    /*!
     * \brief Calculates the land allocation for all items in the land
     *        allocation tree.
     * \details Recursively calculates the landAllocation at each leaf and node
     *          using the shares. The land allocation is passed the value of 0
     *          at the root when this method is called, so the value in the land
     *          allocation variable at the root will not be changed and be
     *          passed down recursively.
     * \author Steve Smith, James Blackwood
     */
    virtual void calcLandAllocation( const std::string& aRegionName,
                                     const double aLandAllocationAbove,
                                     const int aPeriod ) = 0;

    /*!
     * \brief This calls the simple carbon calculator to calculate
     *        land-use change emissions.
     */
    virtual void calcLUCEmissions( const std::string& aRegionName,
                                   const int aPeriod, const int aEndYear ) {}

     /*!
     * \brief Set the profit rate of unmanaged land leafs
     * \details The profit rate of unmanaged land leafs is
     *          equal to the average rental rate on land in 
     *          a region/subregion.  This rate is read in to the land
     *          allocator and passed to all unmanaged land leafs
     * \param aRegionName Region name.
     * \param aAverageProfitRate Region's average profit rate.
     * \param aPeriod Model period
     * \author Kate Calvin
     */
    virtual void setUnmanagedLandProfitRate( const std::string& aRegionName, 
                                             double aAverageProfitRate,
                                             const int aPeriod ) = 0;

    /*!
     * \brief calculate the calibration profit rate
     * \details The calibration profit rate is
     *          determined by the shares and the avergae price of land 
     *          a region/subregion.  
     * \param aRegionName Region name.
     * \param aAverageProfitRate Region's average profit rate.
     * \param aPeriod Model period
     * \author Marshall Wise
     */
    virtual void calculateCalibrationProfitRate( const std::string& aRegionName, 
                                             double aAverageProfitRate,
                                             double aLogitExponentAbove,
                                             const int aPeriod ) = 0;

    /*!
     * \brief calculate the profit scaler adjustment factor
     * \details Profit scalers need to be adjusted when new technologies
     *          are added to a node.  This method calculates that adjustment
     * \param aRegionName Region name.
     * \param aPeriod Model period
     * \author Kate Calvin
     */
    virtual void adjustProfitScalers( const std::string& aRegionName, 
                                const int aPeriod ) = 0;

    /*!
     * \brief Set the share of this land item.
     * \param aShare Share of the land allocated to the parent.
     * \param aPeriod Period.
     * \author James Blackwood
     */
    void setShare( const double aShare,
                   const int aPeriod );
    
    /*!
     * \brief Set the profit scaler of this land item.
     * \param aProfitScaler Profit scaler of the item
     * \param aPeriod Period.
     * \author Kate Calvin
     */
    void setProfitScaler( const double aProfitScaler,
                         const int aPeriod );

    double getShare( const int aPeriod ) const;
        
    virtual double getLogitExponent( const int aPeriod ) const = 0;
 
    const ALandAllocatorItem* getParent() const;

    double getProfitRate( const int aPeriod ) const;

    double getScaledProfitRate( const int aPeriod ) const;

    double getProfitScaler( const int aPeriod ) const;

    virtual double getNewTechProfitScaler( const int aPeriod ) const = 0;

    double getShare( const double aPeriod ) const;

    bool isNewTech( const double aPeriod ) const;

    void setNewTechAdjustment( const double aAdjustment, const double aPeriod );

    LandAllocatorItemType getType() const;
            
    virtual bool isManagedLandLeaf( )  const = 0;

protected:
    virtual void toDebugXMLDerived( const int aPeriod,
                                    std::ostream& aOut,
                                    Tabs* aTabs ) const = 0;

    virtual const std::string& getXMLName() const = 0;

    //! Parent of this node
    const ALandAllocatorItem* mParent;

    /*!
     * \brief Share of parent's total land.
     * \details This is equal to the land allocated to this node divided by land
     *          allocated to node above. This is always the normalized share and
     *          so is always between zero and one inclusive.
     */
    objects::PeriodVector<double> mShare;
    
    //! Profit scaler 
    objects::PeriodVector<double> mProfitScaler;  

    //! Boolean indicating a node or leaf is new 
    objects::PeriodVector<bool> mIsNewTech;  

    //! Double that adjusts a profit scaler to account for the availability of new technologies
    objects::PeriodVector<double> mAdjustForNewTech;

    //! Calibration Profit or Calibration Land Rental Rate in dollars
    // This is the profit rate implied by the shares in the calibration data
    //It is not read in but computed as part of the calibration
    objects::PeriodVector<double> mCalibrationProfitRate;

    //! Land observed profit rate
    objects::PeriodVector<double> mProfitRate;
    
    //! Name of the land allocator item. This is the name of the product for
    //! leafs and name of the type of land for nodes.
    std::string mName;

    /*!
     * \brief Enum that stores the item's type.
     * \note This is stored to avoid a virtual function call.
     */
    LandAllocatorItemType mType;

    //! name of land expansion constraint cost curve
    std::string mLandExpansionCostName;
    bool mIsLandExpansionCost;
    
};

typedef std::unary_function<const ALandAllocatorItem*, bool> SearchPredicate;

/*!
 * \brief SearchPredicate that finds an item with the desired type and name.
 * \details This predicate should be passed to TreeItem's findItem method.
 *          It will be called on each item during the search.
 */
struct MatchesTypeAndName : public SearchPredicate {
    
    /*!
     * \brief Enum that stores the desired type.
     */
    LandAllocatorItemType mType;

    /*!
     * \brief String that stores the desired name.
     */
    const std::string& mName;

    /*!
     * \brief Constructor.
     * \param aName The desired name.
     * \param aType The desired type.
     */
    explicit MatchesTypeAndName( const std::string& aName, LandAllocatorItemType aType )
    : mType( aType ),
      mName( aName )
    {}

    /*!
     * \brief Operator() that returns whether this item is the desired type and name.
     * \param aItem The item to check.
     * \return Whether this item matches the desired type and name.
     */
    bool operator()( const ALandAllocatorItem* aItem ) const {
        return ( mType == eAny || aItem->getType() == mType ) 
               && ( mName == aItem->getName() );
    }
};

#endif // _ALAND_ALLOCATOR_ITEM_H_


