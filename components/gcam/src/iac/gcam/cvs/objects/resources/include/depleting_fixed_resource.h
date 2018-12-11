#ifndef _DEPLETING_FIXED_RESOURCE_H_
#define _DEPLETING_FIXED_RESOURCE_H_
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
 * \file depleting_fixed_resource.h
 * \ingroup Objects
 * \brief DepletingFixedResource header file.
 * \author Pralit Patel
 */
#include <xercesc/dom/DOMNode.hpp>
#include <vector>
#include "resources/include/aresource.h"
#include "util/base/include/value.h"

/*! 
 * \ingroup Objects
 * \brief A class which defines a fixed quantity of resource land which depletes 
 *        with resource extraction from the land.
 * \details The DepletingFixedResource allows for the depletion of a resource
 *          without defining any price behavoir or associtating any extraction costs.
 *          The rate of depletion should be calculated by taking the total quantity 
 *          of land and divide it by the total quantity of energy that would be available.
 *          Note that a rate of zero could be used to emulate a non-depleting resource.
 *          This resource will make that rate available through the "depletion-rate" in
 *          MarketInfo of this resource's market.  A technology would use the rate at the
 *          end of a period (postCalc) to determine the quantity of land that was depleted.
 *          That quantity will be added to the "depleted-resource" in the MarketInfo of this
 *          resource's market.  At the start of the subsequent period (initCalc) this resource
 *          will subtract the depleted quantity from the fixed amount.
 *
 *          <b>XML specification for DepletingFixedResource</b>
 *          - XML name: \c depleting-fixed-resource
 *          - Contained by: Region
 *          - Parsing inherited from class: None
 *          - Attributes: \c name AResource::mName
 *          - Elements:
 *              - \c market AResource::mMarket
 *              - \c fixed-resource DepletingFixedResource::mFixedResource
 *              - \c depletion-rate DepletingFixedResource::mDepletionRate
 *              - \c output-unit AResource::mOutputUnit
 *              - \c price-unit AResource::mPriceUnit
 *              - \c price DepletingFixedResource::mInitialPrice
 *
 * \author Pralit Patel
 */
class DepletingFixedResource: public AResource {
public:
    static const std::string& getXMLNameStatic();

    DepletingFixedResource();

    virtual ~DepletingFixedResource();

    virtual const std::string& getXMLName() const;

    virtual void XMLParse( const xercesc::DOMNode* aNode );

    virtual void toInputXML( std::ostream& aOut,
                             Tabs* aTabs ) const;

    virtual void toDebugXML( const int period,
                             std::ostream &out,
                             Tabs* tabs ) const;

    virtual const std::string& getName() const;

    virtual void completeInit( const std::string& aRegionName,
                               const IInfo* aRegionInfo );

    virtual void initCalc( const std::string& aRegionName,
                           const int aPeriod );

    virtual void postCalc( const std::string& aRegionName,
                           const int aPeriod );

    virtual void calcSupply( const std::string& aRegionName,
                             const GDP* aGDP,
                             const int aPeriod );

    virtual double getAnnualProd( const std::string& aRegionName,
                                  const int aPeriod ) const;

    virtual void dbOutput( const std::string& aRegionName );

    virtual void csvOutputFile( const std::string& aRegionName ); 

    virtual void setCalibratedSupplyInfo( const int aPeriod,
                                          const std::string& aRegionName );

    virtual void accept( IVisitor* aVisitor,
                         const int aPeriod ) const;

protected:

    //! Fixed resource quantity currently available.
    Value mFixedResource;

    //! The rate at which the resource depletes.
    Value mDepletionRate;

    //! The initial price for this resource.
    Value mInitialPrice;

    void setMarket( const std::string& aRegionName );
};

#endif // _DEPLETING_FIXED_RESOURCE_H_
