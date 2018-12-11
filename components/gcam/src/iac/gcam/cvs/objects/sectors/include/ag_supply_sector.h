#ifndef _AG_SUPPLY_SECTOR_H_
#define _AG_SUPPLY_SECTOR_H_
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
* \file ag_supply_sector.h
* \ingroup Objects
* \brief The AgSupplySector class header file.
* \author Marshall Wise, Kate Calvin
*/

#include <string>
#include <xercesc/dom/DOMNode.hpp>
#include "sectors/include/supply_sector.h"

// Forward declarations
class ILandAllocator;

/*!
 * \brief A sector which supplies ag products.
 * \details AgSupplySector differs from the standard supply sector class in one 
 *          critical area - Ag supply sectors are solved markets and determine
 *          their supply indepedently.
 *
 *          1. so methods that set supply equal to demand inside this class are 
 *             overwritten. Instead output is determined by the sum of contained
 *             AgTechnology objects that are based on land allocations and yield
 *
 *          2. the setMarket method is overwritten here to allow a marketname
 *             that can be multi-region (e.g., global)
 *
 *          3. the getPrice method is overwritten to get a price from the
 *             market rather than compute average cost like standard sector
 *
 *          notes by Marshall Wise, March 10, 2010
 */
class AgSupplySector : public SupplySector {
public:
    explicit AgSupplySector( std::string& aRegionName );
    virtual ~AgSupplySector();
    static const std::string& getXMLNameStatic();
    virtual void completeInit( const IInfo* aRegionInfo,
                               DependencyFinder* aDepFinder,
                               ILandAllocator* aLandAllocator );

    virtual void supply( const GDP* aGDP, const int aPeriod );
protected:
	virtual double getPrice( const GDP* aGDP,
                             const int aPeriod ) const;

    virtual bool XMLDerivedClassParse( const std::string& nodeName, const xercesc::DOMNode* curr );
    virtual void toInputXMLDerived( std::ostream& out, Tabs* tabs ) const;
    virtual void toDebugXMLDerived( const int period, std::ostream& out, Tabs* tabs ) const;
    virtual const std::string& getXMLName() const;
    virtual void setMarket();

    // TODO: Should this be a vector?
    double mCalPrice;

    //! Name of the market for this good.
    std::string mMarketName;
};

#endif // _AG_SUPPLY_SECTOR_H_
