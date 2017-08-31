#ifndef _ILANDALLOCATOR_H_
#define _ILANDALLOCATOR_H_
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
 * \file iland_allocator.h
 * \ingroup Objects
 * \brief The ILandAllocator interface file.
 * \author Josh Lurz, Kate Calvin
 */
#include <xercesc/dom/DOMNode.hpp>
#include "util/base/include/ivisitable.h"
#include "util/base/include/iparsable.h"
#include "util/base/include/iround_trippable.h"

// Forward declarations
class Tabs;
class IInfo;
class ALandAllocatorItem;

/*!
 * \brief The interface to a land allocation system.
 * \details This interface represents a method for agricultural production
 *          technologies to interact with a system for distributing land between
 *          usages.
 */
class ILandAllocator : public IVisitable,
                       public IParsable,
                       public IRoundTrippable {
public:
    ILandAllocator();
    
    virtual ~ILandAllocator();

    virtual bool XMLParse( const xercesc::DOMNode* aNode ) = 0;
    
    virtual void toDebugXML( const int aPeriod, std::ostream& aOut, Tabs* aTabs ) const = 0;
    
    virtual void toInputXML( std::ostream& aOut, Tabs* aTabs ) const = 0;
    
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
     * \brief Sets the profit rate for a given product.
     * \details Determines the appropriate land leaf and sets the profit rate
     *          for a given period.
     * \param aRegionName Name of the containing region.
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
     * \brief Get the amount of land allocated for a type of land.
     * \param aProductName Product name.
     * \param aPeriod Model period.
     * \return The land allocated for the product.
     */
    virtual double getLandAllocation( const std::string& aProductName,
                                      const int aPeriod ) const = 0;
    
    /*!
     * \brief Calculate the final land allocation once all yields are known.
     * \param aRegionName Region name.
     * \param aPeriod Model period.
     */
    virtual void calcFinalLandAllocation( const std::string& aRegionName, 
                                          const int aPeriod ) = 0;
                           
    /*!
     * \brief Find a leaf with the given product name.
     * \details Test to improve performance by reducing the number
     *          of times the product need to be located.
     * \return The leaf node if found otherwise null.
     */
    virtual ALandAllocatorItem* findProductLeaf( const std::string& aProductName ) = 0;

    /*!
     * \brief Complete the initialization of the land allocator.
     * \param aRegionName Region name.
     * \param aRegionInfo Local info object.
     */
    virtual void completeInit( const std::string& aRegionName, 
                               const IInfo* aRegionInfo ) = 0;
    
    /*!
     * \brief Perform initializations that are needed before every model period
     * \param aRegionName Region name.
     * \param aPeriod Model period.
     */
    virtual void initCalc( const std::string& aRegionName, 
                           const int aPeriod ) = 0;
                           
    /*!
     * \brief Perform any calculations additional calculations after the model
     *        has finished attempting to solve the given period.
     * \param aRegionName Region name.
     * \param aPeriod Model period.
     */
    virtual void postCalc( const std::string& aRegionName, 
                           const int aPeriod ) = 0;

    virtual void accept( IVisitor* aVisitor, const int aPeriod ) const = 0;
};

// Inline function definitions.

//! Constructor
inline ILandAllocator::ILandAllocator(){
}

//! Destructor.
inline ILandAllocator::~ILandAllocator(){
}

#endif // _ILANDALLOCATOR_H_
