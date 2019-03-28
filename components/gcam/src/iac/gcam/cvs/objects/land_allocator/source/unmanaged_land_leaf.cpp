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
 * \file unmanaged_land_leaf.cpp
 * \ingroup Objects
 * \brief UnmanagedLandLeaf class source file.
 * \author James Blackwood, Kate Calvin
 */

#include "util/base/include/definitions.h"
#include "land_allocator/include/unmanaged_land_leaf.h"
#include "land_allocator/include/land_use_history.h"
#include "util/base/include/xml_helper.h"
#include "ccarbon_model/include/land_carbon_densities.h"
#include "emissions/include/aghg.h"
#include "util/base/include/summary.h"
#include "util/base/include/ivisitor.h"
#include "emissions/include/ghg_factory.h"
#include "marketplace/include/marketplace.h"

using namespace std;
using namespace xercesc;

extern Scenario* scenario;

/*!
 * \brief Default constructor.
 * \param aParent Pointer to this leafs's parent.
 * \author James Blackwood
*/
UnmanagedLandLeaf::UnmanagedLandLeaf( const ALandAllocatorItem* aParent ):
// Default the name to the empty string. It will be read in during XML parsing.
LandLeaf( aParent, "" )
{
}

//! Destructor
UnmanagedLandLeaf::~UnmanagedLandLeaf() {
}

bool UnmanagedLandLeaf::XMLDerivedClassParse( const std::string& aNodeName,
                                              const xercesc::DOMNode* aCurr )
{
    return true;
}

const string& UnmanagedLandLeaf::getXMLName() const {
    return getXMLNameStatic();
}

/*! \brief Get the XML node name in static form for comparison when parsing XML.
*
* This public function accesses the private constant string, XML_NAME.
* This way the tag is always consistent for both read-in and output and can be easily changed.
* The "==" operator that is used when parsing, required this second function to return static.
* \note A function cannot be static and virtual.
* \author James Blackwood
* \return The constant XML_NAME as a static.
*/
const string& UnmanagedLandLeaf::getXMLNameStatic() {
    const static string XML_NAME = "UnmanagedLandLeaf";
    return XML_NAME;
}

/*!
* \brief Sets a the profit rate of a land leaf
* \details This method adjusts the profit rate of an unmanaged land leaf
*          to account for the carbon value of land if the ag subsidy is
*          is active and a carbon price exists. 
* \param aRegionName Region.
* \param aPeriod Period.
*/
void UnmanagedLandLeaf::setUnmanagedLandProfitRate( const string& aRegionName,
                                                    double aAverageProfitRate,
                                                    const int aPeriod )
{
    // Adjust profit rate for land expnasion costs if applicable
    double adjustedProfitRate = aAverageProfitRate;
    const Marketplace* marketplace = scenario->getMarketplace();

    if ( mIsLandExpansionCost ) {
        //subtract off expansion cost from profit rate
        double expansionCost = marketplace->getPrice( mLandExpansionCostName, aRegionName, aPeriod );
        adjustedProfitRate = adjustedProfitRate - expansionCost;
    }

    mProfitRate[ aPeriod ] = max( adjustedProfitRate + getCarbonSubsidy( aRegionName, aPeriod ), 0.0 );
}  


/*!
* \brief Returns the land allocation of this leaf if it is the specified type
* \param aType Land type.
* \param aPeriod Model Period
*/
double UnmanagedLandLeaf::getCalLandAllocation( const LandAllocationType aType,
                                                const int aPeriod ) const
{
    // Check if unmanaged land should be returned.
    if( aType == eAnyLand || aType == eUnmanaged ){
        return mReadinLandAllocation[ aPeriod ];
    }
    return 0;
}

/*!
* \brief Sets a the profit rate of a land leaf
* \param aRegionName Region.
* \param aProductName Name of land leaf
* \param aProfitRate Profit rate
* \param aPeriod Period.
*/
void UnmanagedLandLeaf::setProfitRate( const string& aRegionName,
                                 const string& aProductName,
                                 const double aProfitRate,
                                 const int aPeriod )
{
    // This shouldn't do anything for unmanaged land leafs
}


bool UnmanagedLandLeaf::isManagedLandLeaf( )  const 
{
    return false;
}
