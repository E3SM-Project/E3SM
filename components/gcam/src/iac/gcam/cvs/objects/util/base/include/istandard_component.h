#ifndef _ISTANDARD_COMPONENT_H_
#define _ISTANDARD_COMPONENT_H_
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
 * \file istandard_component.h
 * \ingroup Objects
 * \brief IStandardComponent interface header file.
 * \author Josh Lurz
 */

#include <string>
#include <xercesc/dom/DOMNode.hpp>
#include "util/base/include/iparsable.h"
#include "util/base/include/iround_trippable.h"

class Tabs;

/*! 
 * \ingroup Objects
 * \brief Defines the interface to a standard component.
 * \details The Objects model contains a series of components to represent the
 *          parts of the economic model. These components all share several
 *          common function necessary for the operation of the model. This
 *          interface unifies these methods so that they are always defined in
 *          the exact same way.
 * \author Josh Lurz
 */
class ISimpleComponent { 
public:
	/*!
     * \brief Constructor.
	 * \details Inlined constructor to avoid compiler problems with abstract
     *          base classes. 
     */
    ISimpleComponent();

	/*!
     * \brief Virtual destructor so objects can be deleted through an interface
     *          pointer.
	 * \details Inlined destructor to avoid compiler problems with abstract base
     *          classes. 
     */
	virtual ~ISimpleComponent();

	/*!
     * \brief Creates an exact copy of the component.
	 * \return An exact copy of the component. 
     */
	virtual ISimpleComponent* clone() const = 0;

	/*!
     * \brief Returns whether the type of the object is the same as the passed
     *          in type.
	 * \param aType Type to check the object's type against.
	 * \return Whether the type of the object is the same as the passed in type.
     */
	virtual bool isSameType( const std::string& aType ) const = 0;
    
	/*!
     * \brief Get the name of the component.
	 * \return The name of the component.
	 */
	virtual const std::string& getName() const = 0;

	/*!
     * \brief Write data from this object in an XML format for debugging.
	 * \param aPeriod Period for which to write data.
	 * \param aOut Filestream to which to write.
	 * \param aTabs Object responsible for writing the correct number of tabs. 
     */
	virtual void toDebugXML( const int aPeriod,
                             std::ostream& aOut,
                             Tabs* aTabs ) const = 0;
};

// Inline function definitions.
inline ISimpleComponent::ISimpleComponent(){
}

inline ISimpleComponent::~ISimpleComponent(){
}

/*! 
 * \ingroup Objects
 * \brief Defines the interface to a standard component which is serialized to
 *        XML.
 * \details TODO
 * \todo Is there a better name? DataComponent, XMLComponent, etc?
 * \author Josh Lurz
 */
class IParsedComponent: public ISimpleComponent,
                        public IParsable,
                        public IRoundTrippable 
{ 
public:
	/*!
     * \brief Constructor.
	 * \details Inlined constructor to avoid compiler problems with abstract
     *          base classes. 
     */
    IParsedComponent();
};

// Inline function definitions.
inline IParsedComponent::IParsedComponent(){
}

#endif // _ISTANDARD_COMPONENT_H_
