"""
Interface to the config_workflow.xml file.  This class inherits from GenericXML.py
"""

from CIME.XML.standard_module_setup import *
from CIME.XML.generic_xml import GenericXML
from CIME.XML.files import Files
from CIME.utils import expect

logger = logging.getLogger(__name__)

class Workflow(GenericXML):

    def __init__(self, infile=None, files=None):
        """
        initialize an object
        """
        if files is None:
            files = Files()
        if infile is None:
            infile = files.get_value("WORKFLOW_SPEC_FILE")
        expect(infile, "No workflow file defined in {}".format(files.filename))

        schema = files.get_schema("WORKFLOW_SPEC_FILE")

        GenericXML.__init__(self, infile, schema=schema)

        #Append the contents of $HOME/.cime/config_workflow.xml if it exists
        #This could cause problems if node matchs are repeated when only one is expected
        infile = os.path.join(os.environ.get("HOME"),".cime","config_workflow.xml")
        if os.path.exists(infile):
            GenericXML.read(self, infile)

    def get_workflow_jobs(self, machine, workflowid="default"):
        """
        Return a list of jobs with the first element the name of the script
        and the second a dict of qualifiers for the job
        """
        jobs = []
        bnodes = []
        findmore = True
        prepend = False
        while findmore:
            bnode = self.get_optional_child("workflow_jobs", attributes={"id":workflowid})
            expect(bnode,"No workflow {} found in file {}".format(workflowid, self.filename))
            if prepend:
                bnodes = [bnode] + bnodes
            else:
                bnodes.append(bnode)
            prepend = False
            workflow_attribs = self.attrib(bnode)
            if "prepend" in workflow_attribs:
                workflowid = workflow_attribs["prepend"]
                prepend = True
            elif "append" in workflow_attribs:
                workflowid = workflow_attribs["append"]
            else:
                findmore = False
        for bnode in bnodes:
            for jnode in self.get_children(root=bnode):
                if self.name(jnode) == "job":
                    name = self.get(jnode, "name")
                    jdict = {}
                    for child in self.get_children(root=jnode):
                        if self.name(child) == "runtime_parameters":
                            attrib = self.attrib(child)
                            if attrib and attrib == {'MACH' : machine}:
                                for rtchild in self.get_children(root=child):
                                    jdict[self.name(rtchild)] = self.text(rtchild)
                            elif not attrib:
                                for rtchild in self.get_children(root=child):
                                    if self.name(rtchild) not in jdict:
                                        jdict[self.name(rtchild)] = self.text(rtchild)

                        else:
                            jdict[self.name(child)] = self.text(child)

                jobs.append((name, jdict))

        return jobs
