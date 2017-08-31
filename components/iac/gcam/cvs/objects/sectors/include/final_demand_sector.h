#ifndef _FINAL_DEMAND_SECTOR_H_
#define _FINAL_DEMAND_SECTOR_H_
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
* \file final_demand_sector.h
* \ingroup Objects
* \brief FinalDemandSector class header file.
* \author Pralit Patel
* \author Sonny Kim
*/

#include <string>
#include <xercesc/dom/DOMNode.hpp>

#include "sectors/include/sector.h"

class Tabs;
class IInfo;
class DependencyFinder;

/*!
 * \brief A sector which calculates the final demands for a region using a
 *        series of consumers.
 */
class FinalDemandSector : public Sector
{
public:
	explicit FinalDemandSector( const std::string& aRegionName );
	virtual ~FinalDemandSector();
	virtual void calcFinalSupplyPrice( const GDP* aGDP, const int aPeriod ){};
	virtual void supply( const GDP* aGDP, const int aPeriod ){};
    virtual void operate( NationalAccount& aNationalAccount, const Demographic* aDemographic, 
        const int aPeriod );

	static const std::string& getXMLNameStatic();
    
    virtual void initCalc( NationalAccount* aNationalAccount,
                           const Demographic* aDemographics,
                           const int aPeriod );

    virtual void completeInit( const IInfo* aRegionInfo,
                               DependencyFinder* aDependencyFinder,
                               ILandAllocator* aLandAllocator );

    virtual void dbOutput( const GDP* aGDP,
                           const IndirectEmissionsCalculator* aIndEmissCalc ) const {}
protected:
    virtual double getOutput( const int aPeriod ) const { return 0; }
	
    virtual double getPrice( const GDP* aGDP,
                             const int aPeriod ) const;

    virtual void setMarket();
    virtual const std::string& getXMLName() const;
    virtual bool XMLDerivedClassParse( const std::string& nodeName, const xercesc::DOMNode* curr );
    virtual void toDebugXMLDerived( const int period, std::ostream& out, Tabs* tabs ) const;
	virtual void toInputXMLDerived( std::ostream& out, Tabs* tabs ) const;
};

#endif // _FINAL_DEMAND_SECTOR_H_

