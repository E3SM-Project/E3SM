#ifndef _FUNCTION_MANAGER_H_
#define _FUNCTION_MANAGER_H_
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
* \file function_manager.h
* \ingroup Objects
* \brief FunctionManager class header file.
* \author Pralit Patel
* \author Josh Lurz
*/

#include <map>
#include <string>
#include <memory>

class IFunction;

/*! 
* \ingroup Objects
* \brief A static class used to access production and demand functions.
* \details The FunctionManager contains a mapping of function name to instances
*          of production and demand functions. Each production and demand
*          function is only instantiated once in the model. All technologies and
*          consumers contain a reference to these same functions. The
*          FunctionManager is responsible for instantiating these function
*          objects, distributing them to technologies and consumers, and
*          deallocating them at the end of the model run. The FunctionManager
*          cannot be directly instantiated as it is a static class, the only
*          access to the it is through the static getFunction method.
* \author Pralit Patel, Josh Lurz
* \todo Perhaps rename this to FunctionFactory to have more consistent naming.
*/
class FunctionManager {
public:
    static const IFunction* getFunction( const std::string& aFunctionName );
private:
    FunctionManager();
    ~FunctionManager();

    //! Maps the function's name to the pointer to the funtion
    std::map<std::string, IFunction*> mFunctions;
    typedef std::map<std::string, IFunction*>::iterator FunctionsIterator;
};

#endif // _FUNCTION_MANAGER_H_
