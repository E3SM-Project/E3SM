"""
Interface to the testreporter xml.  This class inherits from GenericXML.py

"""

from CIME.XML.standard_module_setup import *
from CIME.XML.generic_xml import GenericXML
from CIME.utils import expect,get_model

import urllib



class TestReporter(GenericXML):

    def __init__(self):
        """
        initialize an object
        """

        expect(get_model() == 'cesm', "testreport is only meant to populate the CESM test database." )
        self.root = None

        GenericXML.__init__(self)

    def setup_header(self, tagname,machine,compiler,mpilib,testroot,testtype,baseline):
        #
        # Create the XML header that the testdb is expecting to recieve
        #
        tlelem    = ET.Element("testrecord")
        elem      = ET.Element('tag_name')
        elem.text = tagname
        tlelem.append(elem)
        elem                = ET.Element('mach')
        elem.text           = machine
        tlelem.append(elem)
        elem                = ET.Element('compiler',attrib={"version":""})
        elem.text           = compiler
        tlelem.append(elem)
        elem                = ET.Element('mpilib',attrib={"version":""})
        elem.text           = mpilib
        tlelem.append(elem)
        elem                = ET.Element('testroot')
        elem.text           = testroot
        tlelem.append(elem)
        elem                = ET.Element('testtype')
        elem.text           = testtype
        tlelem.append(elem)
        elem   = ET.Element('baselinetag')
        elem.text   = baseline
        tlelem.append(elem)

        self.root=tlelem


    def add_result(self,test_name,test_status):
        #
        # Add a test result to the XML structure.
        #
        tlelem      = ET.Element('tests',attrib={"testname":test_name})
        elem=ET.Element('category',attrib={"name":"casestatus"})
        tlelem.append(elem)
        elem=ET.Element('category',attrib={"name":"comment"})
        elem.text= test_status['COMMENT']
        tlelem.append(elem)

        elem=ET.Element('category',attrib={"name":"compare"})
        elem.text= test_status['BASELINE']
        tlelem.append(elem)

        elem=ET.Element('category',attrib={"name":"memcomp"})
        elem.text= test_status['MEMCOMP']
        tlelem.append(elem)

        elem=ET.Element('category',attrib={"name":"memleak"})
        elem.text= test_status['MEMLEAK']
        tlelem.append(elem)

        elem=ET.Element('category',attrib={"name":"nlcomp"})
        elem.text= test_status['NLCOMP']
        tlelem.append(elem)

        elem=ET.Element('category',attrib={"name":"status"})
        elem.text= test_status['STATUS']
        tlelem.append(elem)

        elem=ET.Element('category',attrib={"name":"tputcomp"})
        elem.text= test_status['TPUTCOMP']
        tlelem.append(elem)

        self.root.append(tlelem)



    def push2testdb(self):
        #
        # Post test result XML to CESM test database
        #
        xmlstr = ET.tostring(self.root,method="xml",encoding="UTF-8")
        username=raw_input("Username:")
        os.system("stty -echo")
        password=raw_input("Password:")
        os.system("stty echo")
        params={'username':username,'password':password,'testXML':xmlstr}
        url="https://csegweb.cgd.ucar.edu/testdb/cgi-bin/processXMLtest.cgi"
        params = urllib.urlencode(params)
        f = urllib.urlopen(url, params)
        #
        # Print any messages from the post command
        #
        print(f.read())
        print(f.code)

