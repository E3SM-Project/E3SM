#ifndef _CALIBRATE_RESOURCE_VISITOR_H_
#define _CALIBRATE_RESOURCE_VISITOR_H_
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
 * \file calibrate_resource_visitor.h
 * \ingroup Objects
 * \brief CalibrateResourceVisitor class header file.
 * \author Kate Calvin
 */

#include "util/base/include/default_visitor.h"
#include <string>

/*! 
 * \ingroup Objects
 * \brief A visitor which determines if a resource has calibrated values
 *        and will adjust price wedges algebraiclly to reproduce those values
 *        given the read in supply curves.
 * \author Kate Calvin
 * \warning This class never actually checks whether calibration is active.
 */
class CalibrateResourceVisitor : public DefaultVisitor {
public:

    CalibrateResourceVisitor( const std::string& aRegionName );

    virtual void startVisitResource( const AResource* aResource,
                                     const int aPeriod );

    virtual void endVisitResource( const AResource* aResource,
                                     const int aPeriod );

    virtual void startVisitSubResource( const SubResource* aSubResource,
                                     const int aPeriod );

    virtual void startVisitSubRenewableResource( const SubRenewableResource* aSubResource,
                                     const int aPeriod );

private:
    //! Name of the Region the for which we are calibrating
    std::string mCurrentRegionName;

    //! Name of the resource currently being tabulated.
    std::string mCurrentResourceName;
};

#endif // _CALIBRATE_RESOURCE_VISITOR_H_
