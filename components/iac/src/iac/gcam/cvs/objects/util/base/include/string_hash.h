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

/* string_hash.h
 * Created: 02/28/2007
 * Version: 02/28/2007
 *
 * This software, which is provided in confidence, was prepared by employees
 * of Pacific Northwest National Laboratory operated by Battelle Memorial
 * Institute. Battelle has certain unperfected rights in the software
 * which should not be copied or otherwise disseminated outside your
 * organization without the express written authorization from Battelle.
 * All rights to the software are reserved by Battelle.   Battelle makes no
 * warranty, express or implied, and assumes no liability or responsibility
 * for the use of this software.
 */

#if !defined( __STRING_HASH_H )
#define __STRING_HASH_H     // prevent multiple includes

// include files ***********************************************************

#include <string>

// namespaces **************************************************************

namespace ObjECTS {

// string_hash *************************************************************

// VariableMap::hash *******************************************************

/*!
 * Return a numeric hash of the specified string
 * \param aString the string to hash
 * \return a numeric hash of the specified string
 * \remarks This is the default STL hash function from
 *          http://www.sgi.com/tech/stl/stl_hash_fun.h
 */
inline size_t string_hash( const std::string& aString )
{
   size_t       h   = 0;
   const char * ptr = aString.c_str();
   for ( ; *ptr; ++ptr )
   {
      h = 5 * h + *ptr;
   }

   return h;
}

} // namespace ObjECTS

#endif   // __STRING_HASH_H

// end of string_hash.h ****************************************************

