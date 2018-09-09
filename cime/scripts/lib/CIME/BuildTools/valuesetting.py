from CIME.XML.standard_module_setup import *

logger = logging.getLogger(__name__)

class ValueSetting(object):

    """Holds data about how a value can be assigned to a variable.

    Note that this class doesn't know or care *which* variable might be
    assigned in this way, only that there is a procedure to perform that
    operation

    Public attributes:
    value - The actual value that will be set.
    do_append - Boolean describing whether the value should be
                appended to the existing value of the variable rather
                than overwriting other settings.
    conditions - Dictionary containing the set of values that different
                 variables have to have to use this setting (e.g.
                 DEBUG="TRUE" might be a condition on a debug flag).
    set_up - List of any commands that have to be executed in the build
             system before this setting can occur.
    tear_down - List of any commands that should be executed to clean up
                after setting the variable.

    Public methods:
    is_ambiguous_with
    has_special_case
    """

    def __init__(self, value, do_append, conditions, set_up, tear_down, force_no_append=False): #  pylint: disable=too-many-arguments
        """Create a ValueSetting object by specifying all its data."""
        self.value = value
        self.do_append = do_append
        self.conditions = conditions
        self.set_up = set_up
        self.tear_down = tear_down
        self.force_no_append = force_no_append

    def is_ambiguous_with(self, other):
        """Check to see if this setting conflicts with another one.

        The purpose of this routine is to see if two settings can coexist
        in the same Macros file, or if doing so would raise an ambiguity
        about which one should be preferred over the other. Note that this
        is a symmetric relation (this function returns the same value if
        self and other are swapped).

        The rules to determine this are as follows:

        1) If one or both settings are appending to the value, there's no
           ambiguity, because both can cooperate to set the value.

        >>> a = ValueSetting('foo', True, dict(), [], [])
        >>> b = ValueSetting('bar', False, dict(), [], [])
        >>> a.is_ambiguous_with(b)
        False
        >>> b.is_ambiguous_with(a)
        False

        2) If the two settings have conflicting conditions, then there
           is no ambiguity because they can't both apply to the same
           build.

        >>> a = ValueSetting('foo', False, {"DEBUG": "TRUE"}, [], [])
        >>> b = ValueSetting('bar', False, {"DEBUG": "FALSE"}, [], [])
        >>> a.is_ambiguous_with(b)
        False

        3) If one setting is strictly more specific than the other, then
           there's no ambiguity, because we prefer the more specific
           setting whenever both apply to a build.

        >>> a = ValueSetting('foo', False, {"DEBUG": "TRUE"}, [], [])
        >>> b = ValueSetting('bar', False, {"DEBUG": "TRUE", "MPILIB": "mpich2"}, [], [])
        >>> a.is_ambiguous_with(b)
        False
        >>> b.is_ambiguous_with(a)
        False

        4) All other cases are considered ambiguous.

        >>> a = ValueSetting('foo', False, dict(), [], [])
        >>> b = ValueSetting('bar', False, dict(), [], [])
        >>> a.is_ambiguous_with(b)
        True
        >>> a = ValueSetting('foo', False, {"DEBUG": "TRUE"}, [], [])
        >>> b = ValueSetting('bar', False, {"MPILIB": "mpich2"}, [], [])
        >>> a.is_ambiguous_with(b)
        True
        """
        # Append check.
        if self.do_append or other.do_append:
            return False
        # Consistency check.
        for var_name in self.conditions:
            if var_name not in other.conditions:
                continue
            if self.conditions[var_name] != other.conditions[var_name]:
                return False
        # Specificity check.
        self_set = set(self.conditions.keys())
        other_set = set(other.conditions.keys())
        if self_set < other_set or other_set < self_set:
            return False
        # Any situation we couldn't resolve is ambiguous.
        return True

    def has_special_case(self, other):
        """Check to see if another setting is a special case of this one.

        The purpose of this routine is to see if one of the settings requires
        conditions that are a strict subset of another's conditions. This is
        used to check whether a setting can be thrown out entirely in favor of a
        more general, but machine-specific setting.

        >>> a = ValueSetting('foo', False, {"DEBUG": "TRUE"}, [], [])
        >>> b = ValueSetting('bar', False, {"DEBUG": "TRUE", "MPILIB": "mpich2"}, [], [])
        >>> c = ValueSetting('bar', False, {"DEBUG": "TRUE", "compile_threaded": "false"}, [], [])
        >>> d = ValueSetting('foo', False, {"DEBUG": "FALSE"}, [], [])
        >>> a.has_special_case(b)
        True
        >>> b.has_special_case(a)
        False
        >>> a.has_special_case(c)
        True
        >>> c.has_special_case(a)
        False
        >>> b.has_special_case(c)
        False
        >>> c.has_special_case(b)
        False
        >>> a.has_special_case(a)
        True
        >>> d.has_special_case(a)
        False
        >>> d.has_special_case(b)
        False
        """
        for var_name in self.conditions:
            if var_name not in other.conditions:
                return False
            if self.conditions[var_name] != other.conditions[var_name]:
                return False
        return True
