"""Interface to an expected failure xml file

Here is an example:

<?xml version= "1.0"?>

<expectedFails version="1.1">
  <test name="ERP_D_Ld10_P36x2.f10_f10_musgs.IHistClm50BgcCrop.cheyenne_intel.clm-ciso_decStart">
    <phase name="RUN">
      <status>FAIL</status>
      <issue>#404</issue>
    </phase>
    <phase name="COMPARE_base_rest">
      <status>PEND</status>
      <issue>#404</issue>
      <comment>Because of the RUN failure, this phase is listed as PEND</comment>
    </phase>
  </test>
  <test name="PFS_Ld20.f09_g17.I2000Clm50BgcCrop.cheyenne_intel">
    <phase name="GENERATE">
      <status>FAIL</status>
      <issue>ESMCI/cime#2917</issue>
    </phase>
    <phase name="BASELINE">
      <status>FAIL</status>
      <issue>ESMCI/cime#2917</issue>
    </phase>
  </test>
</expectedFails>

However, many of the above elements are optional, for human consumption only (i.e., not
parsed here). The only required elements are given by this example:

<?xml version= "1.0"?>

<expectedFails version="1.1">
  <test name="...">
    <phase name="...">
      <status>...</status>
    </phase>
  </test>
</expectedFails>
"""

from CIME.XML.standard_module_setup import *

from CIME.XML.generic_xml import GenericXML
from CIME.expected_fails import ExpectedFails

logger = logging.getLogger(__name__)

class ExpectedFailsFile(GenericXML):

    def __init__(self, infile):
        schema = os.path.join(get_cime_root(), "config", "xml_schemas", "expected_fails_file.xsd")
        GenericXML.__init__(self, infile, schema=schema)

    def get_expected_fails(self):
        """Returns a dictionary of ExpectedFails objects, where the keys are test names"""
        xfails = {}
        test_nodes = self.get_children("test")
        for tnode in test_nodes:
            test_name = self.attrib(tnode)["name"]
            phase_nodes = self.get_children("phase", root=tnode)
            for pnode in phase_nodes:
                phase_name = self.attrib(pnode)["name"]
                status_node = self.get_child("status", root=pnode)
                status = self.text(status_node)
                # issue and comment elements are not currently parsed
                if test_name not in xfails:
                    xfails[test_name] = ExpectedFails()
                xfails[test_name].add_failure(phase_name, status)

        return xfails
