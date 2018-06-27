#ifndef _EMISSIONS_DRIVER_FACTORY_H_
#define _EMISSIONS_DRIVER_FACTORY_H_
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

#include <memory>

/*!
 * \file emissions_driver_factory.h
 * \ingroup Objects
 * \brief EmissionsDriverFactory header file.
 * \author Jim Naslund
 */

class AEmissionsDriver;

/*! 
 * \ingroup Objects
 * \brief A factory which is used to create the various types of emissions drivers.
 * \details The factory encapsulates the creation of various types of emissions drivers.
 *          This simplifies adding new types and also minimizes recompilation.
 * \author Jim Naslund
 */
class EmissionsDriverFactory{

public:
    static std::auto_ptr<AEmissionsDriver> create( const std::string& aType );
    static bool isEmissionsDriverNode( const std::string& aNodeName );

};


#endif // _EMISSIONS_DRIVER_FACTORY_H_
