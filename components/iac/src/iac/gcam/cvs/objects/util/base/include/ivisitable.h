#ifndef _IVISITABLE_H_
#define _IVISITABLE_H_
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
 * \file ivisitable.h
 * \ingroup Objects
 * \brief IVisitable class header file.
 * \author Josh Lurz
 */
class IVisitor;

/*! \brief The IVisitable interface allows an object to be visited by an
 *          IVisitor.
 * \details The interface defines an ability to accept a IVisitor to the object.
 *          An object implementing this interface must define a single function
 *          accept which takes a visitor and a period as a parameters. The object
 *          must call the correct visit methods on the IVisitor with itself and
 *          the period as the parameters. This allows double-dispatch to occur,
 *          or more simply the function called is based on the type of both the
 *          object and the output container. The object must then pass the
 *          IVisitor object onto any of its children which can be visited, using
 *          the accept method.
 */

class IVisitable {
public:
    //! Virtual destructor so that instances of the interface may be deleted
    //! correctly through a pointer to the interface.
    inline virtual ~IVisitable();
    
    /*! \brief Accept a visitor to the object.
     * \details The accept method causes the object to inform the IVisitor that
     *          it is visiting the current object for the given period. The
     *          object must then pass the IVisitor onto any children which can be
     *          visited.
     * \param aVisitor Visitor to accept.
     * \param aPeriod Period
     */
    virtual void accept( IVisitor* aVisitor, const int aPeriod ) const = 0;
};

// Inline methods
IVisitable::~IVisitable(){
}

#endif // _IVISITABLE_H_
