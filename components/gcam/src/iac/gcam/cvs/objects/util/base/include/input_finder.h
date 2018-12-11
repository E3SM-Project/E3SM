#ifndef _INPUT_FINDER_H_
#define _INPUT_FINDER_H_
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
* \file input_finder.h
* \ingroup Objects
* \brief InputFinder class header file.
* \author Josh Lurz
*/

#include "util/base/include/default_visitor.h"
#include <list>

class MiniCAMInput;

/*! 
* \ingroup Objects
* \brief A class which determines the list of inputs used for a portion of the
*        model tree.
* \details The input finder is a visitor which when passed through a portion of
*          the model tree will determine a unique list of inputs used. The input
*          finder offers a simple interface to access this list. The elements of
*          the list are guaranteed to be unique. No other contracts, such as
*          order, are specified for the list. There is also no information
*          stored about which technology used which input, or the quantity it
*          used. This allows the visitor to be used before demands are
*          calculated.
* \author Josh Lurz
*/

class InputFinder : public DefaultVisitor {
public:
    InputFinder();

    virtual void startVisitMiniCAMInput( const MiniCAMInput* aInput,
                                         const int aPeriod );

    virtual void startVisitFinalDemand( const AFinalDemand* aFinalDemand,
                                        const int aPeriod );

    virtual void finish() const;

    // Non-IVisitor interface method.
    const std::list<std::string>& getInputs() const;
private:
    //! List of inputs
    std::list<std::string> mInputs;
};

#endif // _INPUT_FINDER_H_
