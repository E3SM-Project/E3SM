#ifndef _IPARSABLE_H_
#define _IPARSABLE_H_
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
* \file iparsable.h  
* \ingroup Objects
* \brief Header file for the IParsable interface.
* \author Josh Lurz
*/
#include <xercesc/dom/DOMNode.hpp>
/*!
* \ingroup Objects
* \brief An interface to a class which can be parsed by the XMLParser.
* \details This interface represents a contract by a class that it implements
*          the XMLParse function. Classes that implement this interface are
*          classes that initiate an XML parse. Since they implement this
*          interface, these classes can be passed to the XMLHelper as a pointer
*          to the IParsable interface. This allows the XMLParser to then call
*          their parsing routine at the appropriate time, which still
*          controlling XML parsing. This class has no data members and so
*          classes that already inherit from another class may implement this
*          interface without multiple inheritance problems. 
* \author Josh Lurz
*/
class IParsable {
public:
	//! Virtual destructor so that instances of the interface may be deleted
    //! correctly through a pointer to the interface.
    inline virtual ~IParsable();
    
	/*!
     * \brief A function which parses XML starting at a given position in the
     *        DOM tree.
     * \param aNode The current node of a DOM tree.
     * \return Whether the parse completed successfully.
     */
    virtual bool XMLParse( const xercesc::DOMNode* aNode ) = 0;
};

// Inline function definitions.
IParsable::~IParsable(){
}
#endif // _IPARSABLE_H_
