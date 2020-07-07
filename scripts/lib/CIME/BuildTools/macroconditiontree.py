from CIME.XML.standard_module_setup import *
logger = logging.getLogger(__name__)

class MacroConditionTree(object): # pylint: disable=too-many-instance-attributes

    """Tree containing the various possible settings of a specific macro.

    Unlike the PossibleValues class, this class assumes that we have
    finished determining which settings could apply on a given machine. It
    also sorts the settings based on the conditions under which they take
    effect, in preparation for writing out the Macros file itself.

    Public methods:
    merge
    write_out
    """

    def __init__(self, name, settings):
        """Create a MacroConditionTree recursively.

        Arguments:
        name - Name of the variable.
        settings - A list of all settings for this variable.
        """
        # Search for any conditions controlling the number of settings.
        condition = None
        # Prefer the COMPILER attribute as the top level attribute, for
        # readability of the merged file.
        if any("COMPILER" in setting.conditions for setting in settings):
            condition = "COMPILER"
        else:
            # To make merging more effective, sort the conditions.
            all_conditions = []
            for setting in settings:
                all_conditions += setting.conditions.keys()
            if all_conditions:
                condition = sorted(all_conditions)[0]
        if condition is None:
            # If there are no conditions, we have reached a leaf.
            # We combine whatever settings are left; there should be at
            # most one non-appending setting, or an arbitrary number of
            # appending settings.
            self._is_leaf = True
            self._assignments = []
            self._set_up = []
            self._tear_down = []
            self._do_append = True
            for setting in settings:
                if not setting.do_append:
                    self._do_append = False
                    assert len(settings) == 1, \
                        "Internal error in macros: An ambiguity was " \
                        "found after the ambiguity check was complete, " \
                        "or there is a mixture of appending and initial " \
                        "settings in the condition tree."
                self._assignments.append((name, setting.value))
                self._set_up += setting.set_up
                self._tear_down += setting.tear_down
        else:
            # If a condition was found, partition the settings depending on
            # how they use it, and recursively create a tree for each
            # partition.
            self._is_leaf = False
            self._condition = condition
            partition = dict()
            for setting in settings:
                # If some of the settings don't use a condition, we use
                # None to represent that.
                cond_val = setting.conditions.pop(condition, None)
                if cond_val in partition:
                    partition[cond_val].append(setting)
                else:
                    partition[cond_val] = [setting]
            branches = dict()
            for cond_val in partition:
                branches[cond_val] = \
                            MacroConditionTree(name, partition[cond_val])
            self._branches = branches

    # pylint shouldn't concern itself with the way that we access other, since
    # it's actually a member of the same class.
    # pylint:disable=protected-access
    def merge(self, other):
        """Merge another tree with this one.

        This should be considered destructive to both trees. The only valid
        value is the one that's returned.
        """
        if self._is_leaf:
            if other._is_leaf:
                assert self._do_append == other._do_append, \
                    "Internal error in macros: Tried to merge an " \
                    "appending tree with a tree containing initial "\
                    "settings."
                # If both are leaves, just merge the values.
                self._assignments += other._assignments
                self._set_up += other._set_up
                self._tear_down += other._tear_down
                return self
            else:
                # If other is not a leaf, swap the arguments so that self
                # is the one that's not a leaf, handled below.
                return other.merge(self)
        else:
            # If self is not a leaf but other is, it should go in
            # self._branches[None]. The same goes for the case where the
            # conditions don't match, and self._condition is last
            # alphabetically.
            if other._is_leaf or self._condition > other._condition:
                if None in self._branches:
                    self._branches[None] = self._branches[None].merge(other)
                else:
                    self._branches[None] = other
                return self
            else:
                # If the other condition comes last alphabetically, swap
                # the order.
                if self._condition < other._condition:
                    return other.merge(self)
                # If neither is a leaf and their conditions match, merge
                # their sets of branches.
                for (cond_val, other_branch) in other._branches.items():
                    if cond_val in self._branches:
                        self._branches[cond_val] = \
                            self._branches[cond_val].merge(other_branch)
                    else:
                        self._branches[cond_val] = other_branch
                return self
    # pylint:enable=protected-access

    def write_out(self, writer):
        """Write tree to file.

        The writer argument is an object inheriting from MacroWriterBase.
        This function first writes out all the initial settings with
        appropriate conditionals, then the appending settings.
        """
        if self._is_leaf:
            for line in self._set_up:
                writer.write_line(line)
            for (name, value) in self._assignments:
                if self._do_append:
                    writer.append_variable(name, value)
                else:
                    writer.set_variable(name, value)
            for line in self._tear_down:
                writer.write_line(line)
        else:
            condition = self._condition
            # Take care of the settings that don't use this condition.
            if None in self._branches:
                self._branches[None].write_out(writer)
            # Now all the if statements.
            for cond_val in self._branches:
                if cond_val is None:
                    continue
                env_ref = writer.variable_string(condition)
                writer.start_ifeq(env_ref, cond_val)
                self._branches[cond_val].write_out(writer)
                writer.end_ifeq()

def merge_optional_trees(tree, big_tree):
    """Merge two MacroConditionTrees when one or both objects may be `None`."""
    if tree is not None:
        if big_tree is None:
            return tree
        else:
            return big_tree.merge(tree)
    else:
        return big_tree
