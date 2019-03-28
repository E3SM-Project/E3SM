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
* \file atom_registry.cpp
* \ingroup Util
* \brief objects::AtomRegistry class source file.
* \author Josh Lurz
*/
#include "util/base/include/definitions.h"
#include <iostream>
#include <string>

#include "util/base/include/atom.h"
#include "util/base/include/atom_registry.h"
#include "util/base/include/util.h"
#include "util/base/include/hash_map.h"

using namespace std;

namespace objects {
	//! Private constructor to prevent creation of a second object.
	AtomRegistry::AtomRegistry() : mAtoms( new AtomMap( getInitialSize() ) ){
		mIsCurrentlyDeallocating = false;
	}

	//! Destructor
	AtomRegistry::~AtomRegistry() {
		// Set that deallocation has begun.
		mIsCurrentlyDeallocating = true;
	}

	/*! \brief Get a shared pointer to the instance of the AtomRegistry object.
	* \details If the static instance of the AtomRegistry class has not been
	*          created, getInstance() will create it. Otherwise getInstance()
	*          will return a shared pointer to the instance of the AtomRegistry
    *          class. This is the only available access to the registry.
	* \return A pointer to the single instance of the AtomRegistry class.
    */
	AtomRegistry* AtomRegistry::getInstance() {
		// Define a static instance variable. This will be created on the first
        // call of this function, and in every call after that will return a
        // reference to the same variable.
		static AtomRegistry instance;
		return &instance;
	}

	/*! \brief Find an atom the registry is responsible for by name.
	* \details Performs a search of the Atom lookup map to find an Atom with the
	*          requested name. The Atom is returned by constant pointer, so it
	*          cannot be modified or deallocated. If the search fails it returns
	*          null. This search should be very fast as the hashmap is initialized
	*          to a large size to avoid potential collisions.
	* \param aID The string identifier of the atom.
	* \return The atom with the ID aID, null if it is not found.
	*/
	const objects::Atom* AtomRegistry::findAtom( const string& aID ) const {
		AtomMap::const_iterator iter = mAtoms->find( aID );
		return ( iter != mAtoms->end() ) ? iter->second.get() : 0;
	}

	/*! \brief Register an atom with the Atom registry so that it can be fetched
	*          throughout the model and automatically deallocated.
	* \details This method registers an Atom with the registry. The atom list is
	*          first searched to see if the Atom already exists. If the Atom
	*          already exists, This function will delete the passed in atom,
    *          print a warning and return false. Otherwise the atom is added to
	*          the list.
	* \param aAtom An Atom to register.
	* \return Whether the atom was successfully registered.
    */
	bool AtomRegistry::registerAtom( Atom* aAtom ){
        /*! \pre Passed in atom is not null. */
		assert( aAtom );

		// Wrap the raw memory with a shared pointer so it will be deallocated if
		// this function fails.
		boost::shared_ptr<objects::Atom> atom( aAtom );

		// Search for the atom within the list of existing atoms.
		if( findAtom( atom->getID() ) ){
			// Using the output stream currently because the loggers may not
			// be created yet.
			cout << "Error: Attempting to register duplicate atom."
				 << "This may cause undesired behavior." << endl;
			return false;
		}

		// Add that atom to the list so it can be checked against for uniqueness and
		// deallocated.
		mAtoms->insert( make_pair( aAtom->getID(), atom ) );
		return true;
	}

	/*! \brief Return whether deallocation of atoms is currently occurring.
	* \details This function is only compiled if CHECK_DEALLOCATION is turned on. It
	*          allows Atoms to determine if the AtomRegistry is currently
	*          deallocating Atoms, which means it is permissible for the Atom
	*          destructor to be called. This is only for error checking.
	* \return Whether deallocation of atoms is currently occurring.
	*/
	bool AtomRegistry::isCurrentlyDeallocating() const {
		return mIsCurrentlyDeallocating;
	}

    /*! \brief Return the constant for the initial size of the map.
    * \return The initial size to use for the map.
    */
    unsigned int AtomRegistry::getInitialSize(){
	    const unsigned int INITIAL_MAP_SIZE = 103;
        return INITIAL_MAP_SIZE;
    }
}
