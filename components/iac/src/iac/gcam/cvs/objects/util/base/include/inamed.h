#ifndef _INAMED_H_
#define _INAMED_H_
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
 * \file inamed.h
 * \ingroup Objects
 * \brief The INamed interface header file.
 * \author Josh Lurz
 */

#include <string>
/*! 
 * \ingroup Objects
 * \brief An interface that specifies a function to return the name of the
 *        object.
 * \details The getName function specified by this interface must return a
 *          constant string by reference. This means the string must be allocated
 *          permanently by the object, not created on the stack. This is to
 *          ensure that loops on objects calling this function will be fast. The
 *          identifier returned by getName must be unique within the current
 *          context of the object. This means that within any container each
 *          INamed object should return a different identifier. These identifiers
 *          do not have to be globally unique. Names should be non-null.
 */

class INamed {
public:
    //! Destructor.
    virtual inline ~INamed();

    /*!
     * \brief Get the name string from this object.
     * \return The name as a constant reference string.
     */
    virtual const std::string& getName() const = 0;
};

// Inline definitions.
INamed::~INamed(){
}

#endif // _INAMED_H_
