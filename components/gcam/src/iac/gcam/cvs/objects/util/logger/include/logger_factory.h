#ifndef _LOGGER_FACTORY_H_
#define _LOGGER_FACTORY_H_
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
* \file logger_factory.h
* \ingroup Objects
* \brief The LoggerFactory class header file.
* \author Josh Lurz
* \date $Date: 2007/01/11 23:52:34 $
* \version $Revision: 1.4.12.3 $
*/

#include <map>
#include <xercesc/dom/DOMNode.hpp>
#include "util/base/include/iparsable.h"

// Forward Declaration
class Logger;
class Tabs;

/*! 
* \ingroup Objects
* \brief This is a factory class which is used to instantiate create loggers.
* \author Josh Lurz
* \warning This class cannot be instantiated.
* \warning Loggers can only be created by the LoggerFactory.
*/

class LoggerFactory {
    friend class LoggerFactoryWrapper;
public:
    static Logger& getLogger( const std::string& aLogName );
    static void toDebugXML( std::ostream& aOut, Tabs* aTabs );
private:
    static std::map<std::string,Logger*> mLoggers; //!< Map of logger names to loggers.
    static void XMLParse( const xercesc::DOMNode* aRoot );
    static void cleanUp();
    //! Private undefined constructor to prevent creating a LoggerFactory.
    LoggerFactory();
    //! Private undefined copy constructor to prevent  copying a LoggerFactory.
    LoggerFactory( const LoggerFactory& );
    //! Private undefined assignment operator to prevent  copying a LoggerFactory.
    LoggerFactory& operator= ( const LoggerFactory& );
};

/*! 
* \ingroup Objects
* \brief This is a proxy or wrapper class which allows the IParsable functions to be translated into
* the static LoggerFactory calls. This is required because a static class cannot have virtual functions,
* nor can it inherit them. This class does not have any data members, or functions not defined by the
* IParsable class. This class also ensures that the LoggerFactory::cleanUp method is called, so this wrapper 
* class should not be destroyed until the end of the model.
* \todo LoggerFactory should be converted to a singleton instead of a static class so that this is not 
* necessary. 
* \author Josh Lurz
* \warning This class cannot be instantiated.
* \warning Loggers can only be created by the LoggerFactory.
*/
class LoggerFactoryWrapper: public IParsable {
public:
    ~LoggerFactoryWrapper() {
        LoggerFactory::cleanUp();
    }
    bool XMLParse( const xercesc::DOMNode* aRoot ){
        LoggerFactory::XMLParse( aRoot );
        return true;
    }
};

#endif // _LOGGER_FACTORY_H_

