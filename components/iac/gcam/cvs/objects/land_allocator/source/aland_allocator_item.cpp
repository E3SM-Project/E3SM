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
 * \file aland_allocator_item.cpp
 * \ingroup Objects
 * \brief ALandAllocatorItem class source file.
 * \author James Blackwood, Kate Calvin
 */

#include "util/base/include/definitions.h"
#include <xercesc/dom/DOMNodeList.hpp>
#include "util/base/include/xml_helper.h"
#include "land_allocator/include/aland_allocator_item.h"
#include "containers/include/scenario.h"

using namespace std;
using namespace xercesc;

extern Scenario* scenario;

/*!
 * \brief Constructor.
 * \param aParent Pointer to this item's parent.
 * \param aType Enum representing this nodes type.
 * \author James Blackwood, Kate Calvin
 */
ALandAllocatorItem::ALandAllocatorItem( const ALandAllocatorItem* aParent,
                                        const LandAllocatorItemType aType )
: mParent( aParent ),
  mType( aType ),
  mProfitScaler( -1.0 ), // this is so initialization can be checked.
  mAdjustForNewTech( 1.0 ),
  mIsNewTech( false ), 
  mProfitRate( 0 ), 
  mShare( -1.0 ), // this is so initialization can be checked.
  mIsLandExpansionCost( false )
{
}

//! Destructor
ALandAllocatorItem::~ALandAllocatorItem() {
}

void ALandAllocatorItem::setShare( const double aShare,
                                   const int aPeriod )
{
    assert( aShare >= 0 && aShare <= 1 );
    mShare[ aPeriod ] = aShare;
}

void ALandAllocatorItem::setProfitScaler( const double aProfitScaler,
                                         const int aPeriod )
{
    assert( aProfitScaler >= 0 );
    mProfitScaler[ aPeriod ] = aProfitScaler;
}

const string& ALandAllocatorItem::getName() const {
    return mName;
}

/*!
 * \brief Returns the parent of the item.
 * \return ALandAllocatorItem pointer to the parent of this item.
 */
const ALandAllocatorItem* ALandAllocatorItem::getParent() const {
    return mParent;
}

/*!
 * \brief Returns the profit rate for the specified period.
 * \param aPeriod The period to get the rate for.
 * \return double representing the profit rate of this item for the specified
 *         period.
 */
double ALandAllocatorItem::getProfitRate( const int aPeriod ) const {
    // Ensure that profit is positive
    if ( mProfitRate[ aPeriod ] < 0.0 ) {
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::ERROR );
        mainLog << "Profit is negative for leaf " << getName()
                << " in period " << aPeriod << endl;
        exit( -1 );
    }
    return mProfitRate[ aPeriod ];
}

/*!
 * \brief Returns the scaled profit rate for the specified period.
 * \param aPeriod The period to get the rate for.
 * \return double representing the profit rate of this item for the specified
 *         period.
 */
double ALandAllocatorItem::getScaledProfitRate( const int aPeriod ) const {
    // call to getprofitrate ensures that profit is positive
    double unScaledProfitRate = getProfitRate( aPeriod );
    return mProfitScaler[ aPeriod ] * unScaledProfitRate;
}
/*!
 * \brief Returns the share for the specified period.
 * \param aPeriod The period to get the rate for.
 * \return double representing the share of this item for the specified period.
 */
double ALandAllocatorItem::getShare( const int aPeriod ) const {
    return mShare[ aPeriod ];
}

double ALandAllocatorItem::getProfitScaler( const int aPeriod ) const {
    return mProfitScaler[ aPeriod ];
}

/*!
 * \brief Returns an enum representing the type of node (node/leaf).
 * \return Enum representing the type of this item.
 */
LandAllocatorItemType ALandAllocatorItem::getType() const {
    return mType;
}

/*!
 * \brief Returns an boolean indicating whether this is a new technology.
 */
bool ALandAllocatorItem::isNewTech( const double aPeriod ) const {
    return mIsNewTech[ aPeriod ];
}


void ALandAllocatorItem::setNewTechAdjustment( const double aAdjustment, const double aPeriod ) {
    mAdjustForNewTech[ aPeriod ] = aAdjustment;
}

void ALandAllocatorItem::toDebugXML( const int aPeriod, ostream& aOut, Tabs* aTabs ) const {
    XMLWriteOpeningTag ( getXMLName(), aOut, aTabs, mName );

    // write out basic data members
    XMLWriteElement( mProfitRate[ aPeriod ], "ProfitRate", aOut, aTabs );
    XMLWriteElement( mShare[ aPeriod ], "share", aOut, aTabs );
    XMLWriteElement( mProfitScaler[ aPeriod ], "profit-scaler", aOut, aTabs );
    XMLWriteElement( mAdjustForNewTech[ aPeriod ], "adjustment", aOut, aTabs );

    // Call derived class method
    toDebugXMLDerived( aPeriod, aOut, aTabs );

    XMLWriteClosingTag( getXMLName(), aOut, aTabs );
}

