#ifndef _ATOM_REGISTRY_H_
#define _ATOM_REGISTRY_H_
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
* \file atom_registry.h
* \ingroup Util
* \brief The objects::AtomRegistry class header file.
* \author Josh Lurz
*/

#include <list>
#include <boost/shared_ptr.hpp>
#include <boost/noncopyable.hpp>
#include <memory>

// Forward declare the HashMap.
template <class T, class U> class HashMap;

/*!
* \brief The objects namespace.
* \note The AtomRegistry must be in a namespace so that it can interoperate
*       correctly with boost hashing code.
*/
namespace objects {
    /*! 
    * \ingroup util
    * \brief A registry that is responsible for tracking all Atoms created by the
    *        model.
    * \details The AtomRegistry is a class which tracks all Atoms within the model
    *          and ensures they are unique. It also is responsible for deallocating
    *          Atoms. The AtomRegistry is a singleton class which is only accessible
    *          through the static getInstance function. This returns the global
    *          instance of the registry. Registries cannot be allocated and the
    *          single instance is deallocated automatically when the model
    *          completes. The registry will deallocate all registered Atoms at this
    *          time. When Atoms are created they automatically register themselves
    *          with the registry. The registry ensures at this time that the atoms
    *          are unique. Non-unique atoms are not registered and will leak memory
    *          as they will never be deallocated. Atoms cannot be registered outside
    *          of their constructors. Registered atoms are kept for the entire
    *          lifetime of the model. They may be fetched using the findAtom
    *          function which searches the internal hashmap to find the requested
    *          Atom.
    * \author Josh Lurz
    */
    class AtomRegistry: boost::noncopyable {
        friend class Atom;
    public:
		~AtomRegistry();
		static AtomRegistry* getInstance();
		const Atom* findAtom( const std::string& aID ) const;
	private:
		AtomRegistry();
		bool registerAtom( Atom* aAtom );
		bool isCurrentlyDeallocating() const;
        static unsigned int getInitialSize();

		//! Boolean which tracks whether the atom registry is currently
		//! deallocating its atoms for error checking.
		bool mIsCurrentlyDeallocating;

		//! Typedef to simplify using the hashmap of atoms.
		typedef HashMap<const std::string, boost::shared_ptr<Atom> > AtomMap;

		/*! \brief A list of unique Atoms stored as a hashmap for quick
        *          searching.
		* \details The hashmap itself is an auto pointer so that the the hashmap
		*          can be forward declared to prevent excessive recompilation
        *          when the internal hashmap class changes. The Atoms are stored
		*          with a shared_ptr so that they are deallocated automatically.
		*          A standard auto_ptr is not used because that could cause an
		*          Atom to be unintentionally deallocated during a hashmap
		*          resize operation.
        */
		std::auto_ptr<AtomMap> mAtoms;
	};
}

#endif // _ATOM_REGISTRY_H_
