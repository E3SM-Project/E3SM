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

    Public methods:
    add_setting
    ambiguity_check
    dependencies
    to_cond_trees
    """

    def __init__(self, name, setting, specificity, depends):
        """Construct a PossibleValues object.

        The name argument refers to the name of the variable. The other
        arguments are the same as for append_match.
        """
        self.name = name
        # If this is an appending setting, its specificity can't cause it
        # to overwrite other settings.
        if setting.do_append:
            self.settings = []
            self.append_settings = [setting]
            self._specificities = []
            self._depends = []
            self._append_depends = depends
        else:
            self.settings = [setting]
            self.append_settings = []
            self._specificities = [specificity]
            self._depends = [depends]
            self._append_depends = set()

    def add_setting(self, setting, specificity, depends):
        """Add a possible value for a variable.

        Arguments:
        setting - A ValueSetting to start the list.
        specificity - An integer representing how specific the setting is.
                      Low-specificity settings that will never be used will be
                      dropped from the list. The lowest allowed specificity is
                      0.
        depends - A set of variable names, specifying the variables that
                  have to be set before this setting can be used (e.g. if
                  SLIBS refers to NETCDF_PATH, then NETCDF_PATH has to be
                  set first).

        >>> from CIME.BuildTools.valuesetting import ValueSetting
        >>> a = ValueSetting('foo', False, dict(), [], [])
        >>> b = ValueSetting('bar', False, dict(), [], [])
        >>> vals = PossibleValues('var', a, 0, {'dep1'})
        >>> vals.add_setting(b, 1, {'dep2'})
        >>> a not in vals.settings and b in vals.settings
        True
        >>> 'dep1' not in vals.dependencies() and 'dep2' in vals.dependencies()
        True
        >>> vals.add_setting(a, 1, {'dep1'})
        >>> a in vals.settings and b in vals.settings
        True
        >>> 'dep1' in vals.dependencies() and 'dep2' in vals.dependencies()
        True
        """
        if setting.do_append:
            self.append_settings.append(setting)
            self._append_depends |= depends
        else:
            mark_deletion = []
            for i in range(len(self.settings)):
                other_setting = self.settings[i]
                other_specificity = self._specificities[i]
                # Ignore this case if it's less specific than one we have.
                if other_specificity > specificity:
                    if other_setting.has_special_case(setting):
                        return
                # Override cases that are less specific than this one.
                elif other_specificity < specificity:
                    if setting.has_special_case(other_setting):
                        mark_deletion.append(i)
            mark_deletion.reverse()
            for i in mark_deletion:
                del self.settings[i]
                del self._specificities[i]
                del self._depends[i]
            self.settings.append(setting)
            self._specificities.append(specificity)
            self._depends.append(depends)

    def ambiguity_check(self):
        """Check the current list of settings for ambiguity.

        This function raises an error if an ambiguity is found.
        """
        for i in range(len(self.settings)-1):
            for j in range(i+1, len(self.settings)):
                if self._specificities[i] != self._specificities[j]:
                    continue
                other = self.settings[j]
                expect(not self.settings[i].is_ambiguous_with(other),
                       "Variable "+self.name+" is set ambiguously in "
                       "config_compilers.xml. Check the file for these "
                       "conflicting settings: \n1: {}\n2: {}".format(
                           self.settings[i].conditions, other.conditions))

    def dependencies(self):
        """Returns a set of names of variables needed to set this variable."""
        depends = self._append_depends.copy()
        for other in self._depends:
            depends |= other
        return depends

    def to_cond_trees(self):
        """Convert this object to a pair of MacroConditionTree objects.

        This represents the step where the list of possible values is
        frozen and we're ready to convert it into an actual text file. This
        object is checked for ambiguities before conversion.

        The return value is a tuple of two items. The first is a dict of
        condition trees containing all initial settings, with the specificities
        as the dictionary keys. The second is a single tree containing all
        appending settings. If the appending tree would be empty, None is
        returned instead.
        """
        self.ambiguity_check()
        # Get all values of specificity for which we need to make a tree.
        specificities = sorted(list(set(self._specificities)))
        # Build trees, starting from the least specific and working up.
        normal_trees = {}
        for specificity in specificities:
            settings_for_tree = [self.settings[i]
                                 for i in range(len(self.settings))
                                 if self._specificities[i] == specificity]
            normal_trees[specificity] = MacroConditionTree(self.name, settings_for_tree)
        if self.append_settings:
            append_tree = MacroConditionTree(self.name, self.append_settings)
        else:
            append_tree = None
        return (normal_trees, append_tree)
