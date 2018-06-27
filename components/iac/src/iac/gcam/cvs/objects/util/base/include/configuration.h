#ifndef _CONFIGURATION_H_
#define _CONFIGURATION_H_
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
* \file configuration.h
* \ingroup Objects
* \brief The Configuration class header file.
* \author Josh Lurz
* \date $Date: 2007/01/11 00:13:33 $
* \version $Revision: 1.6.2.3 $
*/

#include <xercesc/dom/DOMNode.hpp>
#include <map>
#include <list>
#include <memory>
#include "util/base/include/iparsable.h"

class Tabs;

/*! 
* \ingroup Objects
* \brief This class is used as a container of configuration values which can be accessed throughout the program.
*
* The class is a singleton, so that only one can exist at any time during the program. It parses an XML file
* in the executable directory to instantiate its values. The values can be in one of 5 categories: files, strings,
* bools, ints, and doubles. To add a value it must only be added to the XML file in the appropriate section, the parser
* does not need to be changed. Then to access the variable, use the get method appropriate for the value's type.
*
* \author Josh Lurz
* \bug Bools are currently stored in XML as ints due to conversion problems.
* \warning The class is a singleton, so it may not be created with the constructor. Instead call getInstance to return a pointer to the instance.
* \warning The user must call delete on the object when they are finished with it.
*/

class Configuration: public IParsable {

public:
	static Configuration* getInstance();
	bool XMLParse( const xercesc::DOMNode* tempnode );
	void toDebugXML( std::ostream& out, Tabs* tabs ) const;
	const std::string& getFile( const std::string& key, const std::string& defaultValue = "" ) const;
	const std::string& getString( const std::string& key, const std::string& defaultValue = "" ) const;
	bool getBool( const std::string& key, const bool defaultValue = false ) const;
	int getInt( const std::string& key, const int defaultValue = 0 ) const;
	double getDouble( const std::string& key, const double defaultValue = 0 ) const;
    const std::list<std::string>& getScenarioComponents() const;
private:
    const std::string mLogFile; //!< The name of the log to use.
    static std::auto_ptr<Configuration> gInstance; //!< The static instance of the Configuration class.
	std::map<std::string, std::string> fileMap; //!< A map of the file names the program uses.
	std::map<std::string, std::string> stringMap; //!< A map of the strings the program uses.
	std::map<std::string, bool> boolMap; //!< A map of the bools the program uses.
	std::map<std::string, int> intMap;  //!< A map of the ints the program uses.
	std::map<std::string, double> doubleMap;  //!< A map of the doubles the program uses.
    std::list<std::string> scenarioComponents; //!< An ordered list of add-on files. 
	Configuration();

	//! Private undefined constructor to prevent a programmer from creating a second object.
	Configuration( const Configuration& );

	//! Private undefined assignment operator to prevent a programmer from creating a second object.
	Configuration& operator=( const Configuration& );


};

#endif // _CONFIGURATION_H_

