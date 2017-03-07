from CIME.XML.standard_module_setup import *
from CIME.BuildTools.macroconditiontree import MacroConditionTree

logger = logging.getLogger(__name__)

class PossibleValues(object):

    """Holds a list of settings for a single "Macros" variable.

    This helper class takes in variable settings and, for each one, decides
    whether to throw it out, add it to the list of values, or replace the
    existing list of values with the new, more specific setting.

    This class also performs ambiguity checking; if it is possible at build
    time for more than one setting to match the same variable, this is
    considered an error.

    Public attributes:
    name - The name of the variable.
    settings - The current list of possible initial settings for the
               variable.
    append_settings - A dictionary of lists of possible appending settings
                      for the variable, with the specificity of each list
                      as the associated dictionary key.
    depends - The current list of variables that this variable depends on
              to get its value.

    Public methods:
    add_setting
    ambiguity_check
    to_cond_trees
    """

    def __init__(self, name, setting, specificity, depends):
        """Construct a PossibleValues object.

        The name argument refers to the name of the variable. The other
        arguments are the same as for append_match.
        """
        self.name = name
        self.depends = depends
        # If this is an appending setting, its specificity can't cause it
        # to overwrite other settings, but we want to keep track of it.
        if setting.do_append:
            self.settings = []
            self.append_settings = {specificity: [setting]}
            self._specificity = 0
        else:
            self.settings = [setting]
            self.append_settings = {}
            self._specificity = specificity

    def add_setting(self, setting, specificity, depends):
        """Add a possible value for a variable.

        Arguments:
        setting - A ValueSetting to start the list.
        specificity - An integer representing how specific the setting is.
                      Only the initial settings with the highest
                      specificity and appending settings with at least that
                      specificity will actually be kept in the list. The
                      lowest allowed specificity is 0.
        depends - A set of variable names, specifying the variables that
                  have to be set before this setting can be used (e.g. if
                  SLIBS refers to NETCDF_PATH, then NETCDF_PATH has to be
                  set first).

        >>> a = ValueSetting('foo', False, dict(), [], [])
        >>> b = ValueSetting('bar', False, dict(), [], [])
        >>> vals = PossibleValues('var', a, 0, {'dep1'})
        >>> vals.add_setting(b, 1, {'dep2'})
        >>> a not in vals.settings and b in vals.settings
        True
        >>> 'dep1' not in vals.depends and 'dep2' in vals.depends
        True
        >>> vals.add_setting(a, 1, {'dep1'})
        >>> a in vals.settings and b in vals.settings
        True
        >>> 'dep1' in vals.depends and 'dep2' in vals.depends
        True
        """
        if setting.do_append:
            # Appending settings with at least the current level of
            # specificity should be kept.
            if specificity >= self._specificity:
                if specificity not in self.append_settings:
                    self.append_settings[specificity] = []
                self.append_settings[specificity].append(setting)
                self.depends |= depends
        else:
            # Add equally specific settings to the list.
            if specificity == self._specificity:
                self.settings.append(setting)
                self.depends |= depends
            # Replace the list if the setting is more specific.
            elif specificity > self._specificity:
                self.settings = [setting]
                self._specificity = specificity
                self.depends = depends
        # Do nothing if the setting is less specific.

    def ambiguity_check(self):
        """Check the current list of settings for ambiguity.

        This function raises an error if an ambiguity is found.
        """
        for i in range(len(self.settings)-1):
            for other in self.settings[i+1:]:
                expect(not self.settings[i].is_ambiguous_with(other),
                       "Variable "+self.name+" is set ambiguously in "
                       "config_build.xml. Check the file for these "
                       "conflicting settings: \n1: {}\n2: {}".format(
                           self.settings[i].conditions, other.conditions))

    def to_cond_trees(self):
        """Convert this object to a pair of MacroConditionTree objects.

        This represents the step where the list of possible values is
        frozen and we're ready to convert it into an actual text file. This
        object is checked for ambiguities before conversion.

        The return value is a tuple of two trees. The first contains all
        initial settings, and the second contains all appending settings.
        If either would be empty, None is returned instead.
        """
        self.ambiguity_check()
        if self.settings:
            normal_tree = MacroConditionTree(self.name, self.settings)
        else:
            normal_tree = None
        append_settings = []
        for specificity in self.append_settings:
            if specificity >= self._specificity:
                append_settings += self.append_settings[specificity]
        if append_settings:
            append_tree = MacroConditionTree(self.name, append_settings)
        else:
            append_tree = None
        return (normal_tree, append_tree)
