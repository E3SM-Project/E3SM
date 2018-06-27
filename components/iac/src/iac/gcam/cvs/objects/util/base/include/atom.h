#ifndef _ATOM_H_
#define _ATOM_H_
#if defined(_MSC_VER)
#pragma once
// Ignore warnings about finding the hash function through argument dependent
// lookup. This is the intended behavior.
#pragma warning( disable:4675 )
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
* \file atom.h  
* \ingroup util
* \brief Header file for the objects::Atom class.
* \author Josh Lurz
*/

#include <string>
#include <boost/noncopyable.hpp>

namespace objects {
    /*!
    * \ingroup util
    * \brief An immutable class which represents a single unique string.
    * \details Atoms are unique identifiers which may be located and used throughout
    *          the model. They allow for string sharing as all pointers to the same
    *          string point to the same Atom. It is also faster to compare two Atoms
    *          than two strings as all that is required is checking if the pointers
    *          are equal. Lastly, Atoms are optimized for hash-map storage as they
    *          precompute their hash code. This means searching for Atoms in a
    *          hashmap is nearly as fast as an array reference unless there is a
    *          collision. Atoms cannot be changed once they are created. The memory
    *          management of Atoms is handled by the AtomRegistry and occurs
    *          automatically.
    * \warning Atoms are automatically deallocated, they must always be allocated
    *          with new and should not be either explicitly deallocated or contained
    *          within an auto pointer.
    * \author Josh Lurz
    */
    class Atom : protected boost::noncopyable {
	public:
		Atom( const std::string& aUniqueID );
		~Atom();
		inline const std::string& getID() const;
		inline size_t getHashCode() const;
	protected:
		//! The unique string identifier of the atom.
		const std::string mUniqueID;

		//! The precomputed hash code of the atom.
		size_t mHashCode;
	};

	/*! \brief Atom hash function specialization.
	* \details This function allows the boost hash function to hook into the Atom's
	*          precomputed hash code so that calls to boost the hash function
	*          specialized on Atom pointers will call this function and use the
	*          precomputed value.
	*/
	inline size_t hash_value( const Atom* const& aAtom ){
		return aAtom->getHashCode();
	}


	// Inline function definitions.

	/*! \brief Get the unique identifier of the Atom.
	* \return The unique identifier of the Atom.
	*/
	const std::string& Atom::getID() const {
		return mUniqueID;
	}

	/*! \brief Get the precomputed hash code for the Atom.
	* \details Atoms precompute there hash code to optimize lookups. 
	* \return The hash code for the atom.
	*/
    size_t Atom::getHashCode() const {
		return mHashCode;
	}
}
#endif // _ATOM_H_
