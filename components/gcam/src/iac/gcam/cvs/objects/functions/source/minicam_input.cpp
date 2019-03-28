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
 * \file minicam_input.cpp
 * \ingroup Objects
 * \brief The MiniCAMInput class source file.
 * \author Josh Lurz
 */

#include "util/base/include/definitions.h"
#include <cassert>
#include "functions/include/minicam_input.h"
#include "util/base/include/ivisitor.h"

using namespace std;
using namespace xercesc;

//! Default Constructor
MiniCAMInput::MiniCAMInput()
{
}

//! Destructor
MiniCAMInput::~MiniCAMInput() {
}

const string& MiniCAMInput::getName() const {
    return mName;
}

double MiniCAMInput::getConversionFactor( const int aPeriod ) const {
    return 0; // check this.
}

double MiniCAMInput::getCurrencyDemand( const int aPeriod ) const {
    return 0;
}

// Return 0 for base class.
double MiniCAMInput::getCarbonContent( const int aPeriod ) const {
    return 0;
}

void MiniCAMInput::setCurrencyDemand( double aCurrencyDemand,
                                      const string& aRegionName,
                                      const int aPeriod )
{
    // MiniCAM cannot set currency demand directly.
}

double MiniCAMInput::getPricePaid( const string& aRegionName,
                                   const int aPeriod ) const
{
    // In MiniCAM, price, price paid, and price received are all equal.
    return getPrice( aRegionName, aPeriod );
}

void MiniCAMInput::setPricePaid( double aPricePaid, const int aPeriod ) {
    // MiniCAM cannot directly set price paid.
    assert( false );
}

double MiniCAMInput::getPriceReceived( const string& aRegionName,
                                       const int aPeriod ) const
{
    // In MiniCAM, price, price paid, and price received are all equal.
    return getPrice( aRegionName, aPeriod );
}

double MiniCAMInput::getPriceAdjustment() const {
    return 0;
}

void MiniCAMInput::csvSGMOutputFile( ostream& aFile,
                                     const int period ) const
{
    // MiniCAM should not be outputting a CSV file.
    assert( false );
}

void MiniCAMInput::doInterpolations( const int aYear, const int aPreviousYear,
                                     const int aNextYear, const IInput* aPreviousInput,
                                     const IInput* aNextInput )
{
    // most inputs will not need to do anything
}

void MiniCAMInput::accept( IVisitor* aVisitor, const int period ) const
{
    aVisitor->startVisitMiniCAMInput( this, period );
    aVisitor->endVisitMiniCAMInput( this, period );
}

void MiniCAMInput::copyParamsInto( ProductionInput& aProductionInput,
                                   const int aPeriod ) const
{
    // This should never be called.
    assert( false );
}

void MiniCAMInput::copyParamsInto( DemandInput& aDemandInput,
                                   const int aPeriod ) const
{
    // This should never be called.
    assert( false );
}

void MiniCAMInput::copyParamsInto( NodeInput& aNodeInput,
                                   const int aPeriod ) const
{
    // This should never be called.
    assert( false );
}

void MiniCAMInput::copyParamsInto( TradeInput& aDemandInput,
                                   const int aPeriod ) const
{
    // This should never be called.
    assert( false );
}
