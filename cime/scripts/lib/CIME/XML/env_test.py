"""
Interface to the env_test.xml file.  This class inherits from EnvBase
"""
from CIME.XML.standard_module_setup import *

from CIME.XML.env_base import EnvBase
from CIME.utils import convert_to_type

logger = logging.getLogger(__name__)

class EnvTest(EnvBase):
    # pylint: disable=unused-argument
    def __init__(self, case_root=None, infile="env_test.xml", components=None):
        """
        initialize an object interface to file env_test.xml in the case directory
        """
        EnvBase.__init__(self, case_root, infile)

    def add_test(self,testnode):
        self.root.append(testnode)
        self.write()

    def set_initial_values(self, case):
        """
        The values to initialize a test are defined in env_test.xml
        copy them to the appropriate case env files to initialize a test
        ignore fields set in the BUILD and RUN clauses, they are set in
        the appropriate build and run phases.
        """
        tnode = self.get_node("test")
        for child in tnode:
            if child.text is not None:
                logger.debug("Setting {} to {} for test".format(child.tag, child.text))
                if "$" in child.text:
                    case.set_value(child.tag,child.text,ignore_type=True)
                else:
                    item_type = case.get_type_info(child.tag)
                    value = convert_to_type(child.text,item_type,child.tag)
                    case.set_value(child.tag,value)
        case.flush()
        return

    def set_test_parameter(self, name, value):
        """
        If a node already exists update the value
        otherwise create a node and initialize it to value
        """
        case = self.get_value("TESTCASE")
        tnode = self.get_node("test",{"NAME":case})
        # This xpath statement gets immediate children only (ignoring
        # children of build and run blocks)
        idnode = self.get_optional_node(name, root=tnode, xpath=name)

        if idnode is None:
            newnode = ET.SubElement(tnode, name)
            newnode.text = value
        else:
            idnode.text = value

    def get_test_parameter(self, name):
        case = self.get_value("TESTCASE")
        tnode = self.get_node("test",{"NAME":case})
        # This xpath statement gets immediate children only (ignoring
        # children of build and run blocks)
        value = None
        idnode = self.get_optional_node(name, root=tnode, xpath=name)
        if idnode is not None:
            value = idnode.text
        return value

    def get_step_phase_cnt(self,step):
        bldnodes = self.get_nodes(step)
        cnt = 0
        for node in bldnodes:
            cnt = max(cnt, int(node.attrib["phase"]))
        return cnt

    def get_settings_for_phase(self, name, cnt):
        node = self.get_optional_node(name,attributes={"phase":cnt})
        settings = []
        if node is not None:
            for child in node:
                logger.debug ("Here child is {} with value {}".format(child.tag, child.text))
                settings.append((child.tag, child.text))

        return settings

    def run_phase_get_clone_name(self, phase):
        node = self.get_node("RUN",attributes={"phase":str(phase)})
        if "clone" in node.attrib:
            return node.attrib["clone"]
        return None

    def cleanupnode(self, node):
        '''
        keep the values component set
        '''
        fnode = node.find(".//file")
        node.remove(fnode)
        gnode = node.find(".//group")
        node.remove(gnode)
        dnode = node.find(".//default_value")
        if dnode is not None:
            node.remove(dnode)
        return node

    def set_value(self, vid, value, subgroup=None, ignore_type=False):
        """
        check if vid is in test section of file
        """
        newval = EnvBase.set_value(self, vid, value, subgroup, ignore_type)
        if newval is None:
            tnode = self.get_optional_node("test")
            if tnode is not None:
                newval = self.set_element_text(vid, value, root=tnode)
        return newval

    def get_value(self, vid, attribute=None, resolved=True, subgroup=None):
        value = EnvBase.get_value(self, vid, attribute, resolved, subgroup)
        if value is None:
            tnode = self.get_optional_node("test")
            if tnode is not None:
                value = self.get_element_text(vid, root=tnode)
        return value
