#ifndef _XML_LOGGER_H_
#define _XML_LOGGER_H_
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
* \file xml_logger.h
* \ingroup Util
* \brief The XMLLogger class header file.
* \author Josh Lurz
* \date $Date: 2007/01/11 23:52:34 $
* \version $Revision: 1.4.2.3 $
*/

#include "util/logger/include/logger.h"

/*! 
* \ingroup Objects
* \brief This is a class which implements the Logger interface.. 
*
* This Logger prints message to a log file in an XML format as defined by the XMLLog schema..
* It does not support nesting or options to activate or deactivate sections of the log message.
* It does support the option to print absolute or relative pathnames.
*
* \author Josh Lurz
* \warning Since XMLLoggers can only be created by the LoggerFactory, public functions not in the Logger interface will be unusable.
*/

class XMLLogger: public Logger {
    friend class LoggerFactory;
public:
    void open( const char[] = 0 );
    void close();
    void logCompleteMessage( const std::string& aMessage );	

private:
    std::ofstream mLogFile; //!< The filestream to which data is written.
    XMLLogger( const std::string& loggerName ="" );
};

#endif // _XML_LOGGER_H_

