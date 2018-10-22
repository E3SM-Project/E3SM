"""
Interface to the config_component.xml files.  This class inherits from EntryID.py
"""
from CIME.XML.standard_module_setup import *

from CIME.XML.entry_id import EntryID
from CIME.XML.files import Files
from CIME.utils import get_cime_root

logger = logging.getLogger(__name__)

class Component(EntryID):

    def __init__(self, infile, comp_class):
        """
        initialize a Component obect from the component xml file in infile
        associate the component class with comp_class if provided.
        """
        self._comp_class = comp_class
        if infile == 'testingonly':
            self.filename = infile
            return
        files = Files()
        schema = None
        EntryID.__init__(self, infile)
        schema = files.get_schema("CONFIG_{}_FILE".format(comp_class), attributes={"version":"{}".format(self.get_version())})

        if schema is not None:
            self.validate_xml_file(infile, schema)

    #pylint: disable=arguments-differ
    def get_value(self, name, attribute=None, resolved=False, subgroup=None):
        expect(subgroup is None, "This class does not support subgroups")
        return EntryID.get_value(self, name, attribute, resolved)

    def get_valid_model_components(self):
        """
        return a list of all possible valid generic (e.g. atm, clm, ...) model components
        from the entries in the model CONFIG_CPL_FILE
        """
        components = []
        comps_node = self.get_child("entry", {"id":"COMP_CLASSES"})
        comps = self.get_default_value(comps_node)
        components = comps.split(',')
        return components

    def _get_value_match(self, node, attributes=None, exact_match=False):
        """
        return the best match for the node <values> entries
        Note that a component object uses a different matching algorithm than an entryid object
        For a component object the _get_value_match used is below  and is not the one in entry_id.py
        """
        match_value = None
        match_max = 0
        match_count = 0
        match_values = []
        expect(not exact_match, " exact_match not implemented in this method")
        expect(node is not None," Empty node in _get_value_match")
        values = self.get_optional_child("values", root=node)
        if values is None:
            return

        # determine match_type if there is a tie
        # ASSUME a default of "last" if "match" attribute is not there
        match_type = self.get(values, "match", default="last")

        # use the default_value if present
        val_node = self.get_optional_child("default_value", root=node)
        if val_node is None:
            logger.debug("No default_value for {}".format(self.get(node, "id")))
            return val_node
        value = self.text(val_node)
        if value is not None and len(value) > 0 and value != "UNSET":
            match_values.append(value)

        for valnode in self.get_children("value", root=values):
            # loop through all the keys in valnode (value nodes) attributes
            for key,value in self.attrib(valnode).items():
                # determine if key is in attributes dictionary
                match_count = 0
                if attributes is not None and key in attributes:
                    if re.search(value, attributes[key]):
                        logger.debug("Value {} and key {} match with value {}".format(value, key, attributes[key]))
                        match_count += 1
                    else:
                        match_count = 0
                        break

            # a match is found
            if match_count > 0:
                # append the current result
                if self.get(values, "modifier") == "additive":
                    match_values.append(self.text(valnode))

                # replace the current result if it already contains the new value
                # otherwise append the current result
                elif self.get(values, "modifier") == "merge":
                    if self.text(valnode) in match_values:
                        del match_values[:]
                    match_values.append(self.text(valnode))

                else:
                    if match_type == "last":
                        # take the *last* best match
                        if match_count >= match_max:
                            del match_values[:]
                            match_max = match_count
                            match_value = self.text(valnode)
                    elif match_type == "first":
                        # take the *first* best match
                        if match_count > match_max:
                            del match_values[:]
                            match_max = match_count
                            match_value = self.text(valnode)
                    else:
                        expect(False, "match attribute can only have a value of 'last' or 'first'")

        if len(match_values) > 0:
            match_value = " ".join(match_values)

        return match_value

    #pylint: disable=arguments-differ
    def get_description(self, compsetname):
        if self.get_version() == 3.0:
            return self._get_description_v3(compsetname, self._comp_class)
        else:
            return self._get_description_v2(compsetname)

    def get_forcing_description(self, compsetname):
        if self.get_version() == 3.0:
            return self._get_description_v3(compsetname, 'forcing')
        else:
            return ""

    def _get_description_v3(self, compsetname, comp_class):
        """
        version 3 of the config_component.xml file has the description section at the top of the file
        the description field has one attribute 'modifier_mode' which has allowed values
        '*' 0 or more modifiers (default)
        '1' exactly 1 modifier
        '?' 0 or 1 modifiers
        '+' 1 or more modifiers

        modifiers are fields in the component section of the compsetname following the % symbol.

        The desc field can have an attribute which is the component class ('cpl', 'atm', 'lnd' etc)
        or it can have an attribute 'option' which provides descriptions of each optional modifier
        or (in the config_component_{model}.xml in the driver only) it can have the attribute 'forcing'

        component descriptions are matched to the compsetname using a set method
        """
        expect(comp_class is not None,"comp_class argument required for version3 files")
        comp_class = comp_class.lower()
        rootnode = self.get_child("description")
        desc = ""
        desc_nodes = self.get_children("desc", root=rootnode)

        modifier_mode = self.get(rootnode, 'modifier_mode')
        if modifier_mode is None:
            modifier_mode = '*'
        expect(modifier_mode in ('*','1','?','+'),
               "Invalid modifier_mode {} in file {}".format(modifier_mode, self.filename))
        optiondesc = {}
        if comp_class == "forcing":
            for node in desc_nodes:
                forcing = self.get(node, 'forcing')
                if forcing is not None and compsetname.startswith(forcing+'_'):
                    expect(len(desc)==0,
                           "Too many matches on forcing field {} in file {}".\
                               format(forcing, self.filename))
                    desc = self.text(node)
            if desc is None:
                desc = compsetname.split('_')[0]
            return desc


        # first pass just make a hash of the option descriptions
        for node in desc_nodes:
            option = self.get(node, 'option')
            if option is not None:
                optiondesc[option] = self.text(node)

        #second pass find a comp_class match
        desc = ""
        for node in desc_nodes:
            compdesc = self.get(node, comp_class)

            if compdesc is not None:
                opt_parts = [ x.rstrip("]") for x in compdesc.split("[%") ]
                parts = opt_parts.pop(0).split("%")
                reqset = set(parts)
                fullset = set(parts+opt_parts)
                match, complist =  self._get_description_match(compsetname, reqset, fullset, modifier_mode)
                if match:
                    desc = self.text(node)
                    for opt in complist:
                        if opt in optiondesc:
                            desc += optiondesc[opt]


        # cpl and esp components may not have a description
        if comp_class not in ['cpl','esp']:
            expect(len(desc) > 0,
                   "No description found for comp_class {} matching compsetname {} in file {}, expected match in {} % {}"\
                       .format(comp_class,compsetname, self.filename, list(reqset), list(opt_parts)))
        return desc

    def _get_description_match(self, compsetname, reqset, fullset, modifier_mode):
        """

        >>> obj = Component('testingonly', 'ATM')
        >>> obj._get_description_match("1850_DATM%CRU_FRED",set(["DATM"]), set(["DATM","CRU","HSI"]), "*")
        (True, ['DATM', 'CRU'])
        >>> obj._get_description_match("1850_DATM%FRED_Barn",set(["DATM"]), set(["DATM","CRU","HSI"]), "*")
        (False, None)
        >>> obj._get_description_match("1850_DATM_Barn",set(["DATM"]), set(["DATM","CRU","HSI"]), "?")
        (True, ['DATM'])
        >>> obj._get_description_match("1850_DATM_Barn",set(["DATM"]), set(["DATM","CRU","HSI"]), "1")
        Traceback (most recent call last):
        ...
        SystemExit: ERROR: Expected exactly one modifer found 0 in ['DATM']
        >>> obj._get_description_match("1850_DATM%CRU%HSI_Barn",set(["DATM"]), set(["DATM","CRU","HSI"]), "1")
        Traceback (most recent call last):
        ...
        SystemExit: ERROR: Expected exactly one modifer found 2 in ['DATM', 'CRU', 'HSI']
        >>> obj._get_description_match("1850_CAM50%WCCM%RCO2_Barn",set(["CAM50", "WCCM"]), set(["CAM50","WCCM","RCO2"]), "*")
        (True, ['CAM50', 'WCCM', 'RCO2'])

        # The following is not allowed because the required WCCM field is missing
        >>> obj._get_description_match("1850_CAM50%RCO2_Barn",set(["CAM50", "WCCM"]), set(["CAM50","WCCM","RCO2"]), "*")
        (False, None)
        >>> obj._get_description_match("1850_CAM50_Barn",set(["CAM50", "WCCM"]), set(["CAM50","WCCM","RCO2"]), "+")
        (False, None)
        >>> obj._get_description_match("1850_CAM50%WCCM_Barn",set(["CAM50", "WCCM"]), set(["CAM50","WCCM","RCO2"]), "+")
        (True, ['CAM50', 'WCCM'])
        """
        match = False
        comparts = compsetname.split('_')
        matchcomplist = None

        for comp in comparts:
            complist = comp.split('%')
            cset = set(complist)
            if cset == reqset or (cset > reqset and cset <= fullset):
                if modifier_mode == '1':
                    expect(len(complist) == 2,
                           "Expected exactly one modifer found {} in {}".format(len(complist)-1,complist))
                elif modifier_mode == '+':
                    expect(len(complist) >= 2,
                           "Expected one or more modifers found {} in {}".format(len(complist)-1, list(reqset)))
                elif modifier_mode == '?':
                    expect(len(complist) <= 2,
                           "Expected 0 or one modifers found {} in {}".format(len(complist)-1, complist))
                expect(not match,"Found multiple matches in file {} for {}".format(self.filename,comp))
                match = True
                matchcomplist = complist
                # found a match

        return match, matchcomplist

    def _get_description_v2(self, compsetname):
        rootnode = self.get_child("description")
        desc = ""
        desc_nodes = self.get_children("desc", root=rootnode)
        for node in desc_nodes:
            compsetmatch = self.get(node, "compset")
            if compsetmatch is not None and re.search(compsetmatch, compsetname):
                desc += self.text(node)

        return desc

    def print_values(self):
        """
        print values for help and description in target config_component.xml file
        """
        helpnode = self.get_child("help")
        helptext = self.text(helpnode)
        logger.info(" {}".format(helptext))
        entries = self.get_children("entry")
        for entry in entries:
            name = self.get(entry, "id")
            text = self.text(self.get_child("desc", root=entry))
            logger.info("   {:20s} : {}".format(name, text.encode('utf-8')))

    def return_values(self):
        """
        return a list of hashes from target config_component.xml file
        This routine is used by external tools in https://github.com/NCAR/CESM_xml2html
        """
        entry_dict = dict()
        items = list()
        helpnode = self.get_optional_child("help")
        if helpnode:
            helptext = self.text(helpnode)
        else:
            helptext = ''
        entries = self.get_children("entry")
        for entry in entries:
            item = dict()
            name = self.get(entry, "id")
            datatype = self.text(self.get_child("type", root=entry))
            valid_values = self.get_valid_values(name)
            default_value = self.get_default_value(node=entry)
            group = self.text(self.get_child("group", root=entry))
            filename = self.text(self.get_child("file", root=entry))
            text = self.text(self.get_child("desc", root=entry))
            item = {"name":name,
                    "datatype":datatype,
                    "valid_values":valid_values,
                    "value":default_value,
                    "group":group,
                    "filename":filename,
                    "desc":text.encode('utf-8')}
            items.append(item)
        entry_dict = {"items" : items}

        return helptext, entry_dict
