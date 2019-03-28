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
 * \file cal_quantity_tabulator.cpp
 * \ingroup Objects
 * \brief The CalQuantityTabulator class source file.
 * \author Josh Lurz
 */

#include "util/base/include/definitions.h"
#include <cassert>
#include "sectors/include/cal_quantity_tabulator.h"
#include "technologies/include/technology.h"
#include "containers/include/region.h"
#include "technologies/include/ioutput.h"
#include "resources/include/resource.h"
#include "util/base/include/util.h"
#include "sectors/include/subsector.h"
#include "sectors/include/sector.h"

using namespace std;

/*!
 * \brief Constructor
 * \param aRegionName Name of the region if starting the visiting below the
 *        region level.
 */
CalQuantityTabulator::CalQuantityTabulator( const string& aRegionName )
:mCurrentRegionName( aRegionName ),
mCurrentOutput( 0 ),
mTechState( eUnknown ),
mSubsectorState( eUnknown ),
mShouldTabulateSector( false )
{
}

void CalQuantityTabulator::startVisitRegion( const Region* aRegion,
                                             const int aPeriod )
{
    mCurrentRegionName = aRegion->getName();
}

void CalQuantityTabulator::startVisitSector( const Sector* aSector,
                                             const int aPeriod )
{
    mCurrentSectorName = aSector->getName();

    // Check if the visitor should tabulate this sector.
    mShouldTabulateSector = mSectorType.empty() ||
                            mSectorType == aSector->mSectorType;
}

void CalQuantityTabulator::endVisitSector( const Sector* aSector,
                                           const int aPeriod )
{
    mCurrentSectorName.clear();
}

void CalQuantityTabulator::startVisitResource( const AResource* aResource,
                                               const int aPeriod )
{
    mCalSupplies[ aResource->getName() ].mAllFixed = false;
}

void CalQuantityTabulator::startVisitSubsector( const Subsector* aSubsector,
                                                const int aPeriod )
{
    // Check that the state is initialized.
    assert( mSubsectorState == eUnknown );
    assert( !mCurrentSectorName.empty() );

    // Check if the subsector has a zero share weight which implies a fixed
    // output of zero.
    if ( util::isEqual( aSubsector->getShareWeight( aPeriod ), 0.0 ) ){
        mSubsectorState = eFixed;
        mCurrentOutput = 0;
    }
    else {
        mSubsectorState = eVariable;
    }
}

void CalQuantityTabulator::endVisitSubsector( const Subsector* aSubsector,
                                             const int aPeriod )
{
    assert( mSubsectorState != eUnknown );
    mSubsectorState = eUnknown;
}

void CalQuantityTabulator::startVisitTechnology( const Technology* aTechnology,
                                                const int aPeriod )
{
    // Check that the subsector state is already set and the technology state is
    // not set yet.
    assert( mSubsectorState != eUnknown && mTechState == eUnknown );

    // If the subsector state is fixed or calibrated set the Technology to have
    // the same state.
    if( mSubsectorState != eVariable ){
        mTechState = mSubsectorState;
        // Even if the subsector is fixed to zero the technology will attempt to
        // produce its fixed output.
        if( mTechState == eFixed ){
            mCurrentOutput = aTechnology->getFixedOutput( mCurrentRegionName,
                                                          mCurrentSectorName,
                                                          false,
                                                          "",
                                                          aPeriod );

            // Don't allow the current output to be set to the -1 flag. This
            // could occur if the subsector had a zero shareweight and the
            // technology had variable output.
            mCurrentOutput = max( mCurrentOutput, 0.0 );
        }
    }

    // Output is fixed, it could be calibrated or fixed.
    else if( aTechnology->isOutputFixed( false, "", aPeriod ) ){
        // Check if the technology is calibrated.
        double calOutput = aTechnology->getCalibrationOutput( false, "", aPeriod );
        if( calOutput != -1 ) {
            mTechState = eCalibrated;
            // Add calibration outputs for all outputs.
            mCurrentOutput = calOutput;
        }
        else {
            mCurrentOutput = aTechnology->getFixedOutput( mCurrentRegionName,
                                                          mCurrentSectorName,
                                                          false,
                                                          "",
                                                          aPeriod );
            assert( mCurrentOutput != -1 );
            mTechState = eFixed;
        }
    }
    else {
        mTechState = eVariable;
        mCurrentOutput = 0;
    }
}

void CalQuantityTabulator::endVisitTechnology( const Technology* aTechnology,
                                               const int aPeriod )
{
    assert( mTechState != eUnknown );
    mTechState = eUnknown;
    mCurrentOutput = 0;
}

void CalQuantityTabulator::startVisitOutput( const IOutput* aOutput,
                                             const int aPeriod )
{
    if( !mShouldTabulateSector ){
        return;
    }

    // Make sure there is a stored region name.
    assert( !mCurrentRegionName.empty() );

    // Check that the current technology state is known.
    assert( mTechState != eUnknown );

    // Add calibrated and fixed supplies.
    typedef IOutput::OutputList::const_iterator COutputListIterator;
    if( mTechState == eCalibrated ){
        IOutput::OutputList outputList =
            aOutput->calcPhysicalOutput( mCurrentOutput,
                                         mCurrentRegionName,
                                         0,
                                         aPeriod );

        for( COutputListIterator i = outputList.begin(); i != outputList.end(); ++i ) {
            mCalSupplies[ i->first ].mCalQuantity += i->second;
        }
    }
    else if( mTechState == eFixed ){
        IOutput::OutputList outputList =
            aOutput->calcPhysicalOutput( mCurrentOutput,
                                         mCurrentRegionName,
                                         0,
                                         aPeriod );

        for( COutputListIterator i = outputList.begin(); i != outputList.end(); ++i ) {
            mCalSupplies[ i->first ].mFixedQuantity += i->second;
        }
    }
    // Set that at least one technology which produces the output good is not
    // all fixed.
    else {
        mCalSupplies[ aOutput->getName() ].mAllFixed = false;
    }
}

/*!
 * \brief Set the type of sector for which to tabulate calibrated output.
 * \details Instructs the object to only tabulate demands for sectors with a
 *          given type. If this function is not called the object will tabulate
 *          for all sector types. This function must be called before the
 *          visiting occurs.
 * \param aSectorType The sector type for which to tabulate calibrated output.
 */
void CalQuantityTabulator::setApplicableSectorType( const string& aSectorType ){
    mSectorType = aSectorType;
}

const CalQuantityTabulator::CalInfoMap& CalQuantityTabulator::getSupplyInfo() const {
    return mCalSupplies;
}
