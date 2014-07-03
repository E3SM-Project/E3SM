"""Produce commands for environment module systems.

Public classes:
EnvSystemInterface - Abstract base class for interfaces.
NoModuleInterface - Can be constructed, but raises exception if used.
ModuleInterface - For standard module systems (including Lmod).
SoftEnvInterface - For SoftEnv, the ANL MCS environment system.

Public routines:
expand_env - Expand shell variables (like "${FOO}").
"""

import os
import re
import subprocess

__all__ = ("EnvSystemInterface", "NoModuleInterface", "ModuleInterface",
           "SoftEnvInterface", "expand_env")

class EnvSystemInterface(object):

    """Abstract base class for environment system interfaces.

    This class exists only to document the interface of subclasses, and to
    provide error messages if an unimplemented function is called by
    accident.

    Normally, none of these methods should actually be called. All of them
    raise NotImplementedError.

    Public methods:
    is_loaded
    purge
    load
    unload
    purge_str
    load_str
    unload_str
    """

    @classmethod
    def _raise_not_implemented(cls, method_name):
        raise NotImplementedError(
            cls.__name__ +" does not implement method "+method_name+".")

    def is_loaded(self, modname):
        """Method not implemented."""
        self._raise_not_implemented("is_loaded")

    def purge(self):
        """Method not implemented."""
        self._raise_not_implemented("purge")

    def load(self, modname):
        """Method not implemented."""
        self._raise_not_implemented("load")

    def unload(self, modname):
        """Method not implemented."""
        self._raise_not_implemented("unload")

    def purge_str(self):
        """Method not implemented."""
        self._raise_not_implemented("purge_str")

    def load_str(self, modname):
        """Method not implemented."""
        self._raise_not_implemented("load_str")

    def unload_str(self, modname):
        """Method not implemented."""
        self._raise_not_implemented("unload_str")


class NoModuleInterface(EnvSystemInterface):

    """Module interface for systems with no module system.

    The purpose of this class is to allow code to construct a module
    interface regardless of system, while still raising an error if the
    user attempts to actually interact with the non-existent interface.

    As a result, only purge_str is implemented.

    Public methods:
    purge_str - Returns null command.
    """

    def purge_str(self):
        """Returns ":", the null shell command."""
        return ":"


class ModuleInterface(EnvSystemInterface):

    """Module interface for systems with a typical module system.

    Class methods:
    python_init - Initialize python interface.

    Public methods:
    is_loaded - Tests if a module is currently loaded.
    purge - Purge all modules from environment.
    load - Load a module.
    unload - Unload a module.
    purge_str - Return command to purge modules.
    load_str - Return command to load the given module.
    unload_str  - Return command to unload the given module.
    """

    # Singleton declaring whether or not we have the module interface
    # started up yet.
    _python_initialized = False

    # Magic string to grep for in "module list" output to determine whether
    # or not a module is present.
    _not_loaded_string = "None found"

    @classmethod
    def python_init(cls, filename):
        """Initialize python to module system interface.

        This must be called before any commands that are intended to
        modify python's environment. The *_str methods can still be called
        without this. Only the first call has any effect.

        Arguments:
        filename - Python file to execute to load the interface.
        """
        if not cls._python_initialized:
            execfile(filename, globals())
            cls._python_initialized = True

    def is_loaded(self, modname):
        """"Whether the module is loaded in the current environment."""
        # This assertion message isn't quite true, but this behavior is
        # convenient for testing purposes, which is why is_loaded exists.
        assert self._python_initialized, \
            "Can't test modules without initializing the python interface!"
        process = subprocess.Popen("module list "+modname,
                                   shell=True,
                                   env=os.environ,
                                   stderr=subprocess.PIPE)
        stderr_output = process.communicate()[1]
        # Handle unexpected error in subprocess.
        if process.returncode not in (0,1):
            raise Exception("module list command failed with code: "+
                            str(process.returncode))
        # Hack! Assume that process.returncode == 1 means that no modules
        # are loaded, rather than some other possible error.
        return process.returncode != 1 and \
            re.search(self._not_loaded_string, stderr_output) is None

    def purge(self):
        """"Purge all modules from the current environment."""
        assert self._python_initialized, \
            "Can't purge modules without initializing the python interface!"
        module("purge")

    def load(self, modname):
        """"Load a module in the current environment."""
        assert self._python_initialized, \
            "Can't load modules without initializing the python interface!"
        module("load", modname)

    def unload(self, modname):
        """"Unload a module from the current environment."""
        assert self._python_initialized, \
            "Can't unload modules without initializing the python interface!"
        module("unload", modname)

    def purge_str(self):
        """Returns the shell command to purge modules as a string."""
        return "module purge"

    def load_str(self, modname):
        """Returns the shell command to load the module as a string."""
        return "module load "+modname

    def unload_str(self, modname):
        """Returns the shell command to unload the module as a string."""
        return "module unload "+modname


class SoftEnvInterface(EnvSystemInterface):

    """Module interface for systems with SoftEnv.

    Public methods:
    purge_str - Return command to reset SoftEnv.
    load_str - Return command to add the given keyword/macro.
    unload_str  - Return command to delete the given keyword/macro.
    """

    def purge_str(self):
        """Returns the shell command to reset SoftEnv as a string.

        Unfortunately, there's not a straightforward way to clear
        everything, so right now we just run "resoft". If it becomes
        necessary in the future, this may change to do something else (e.g.
        resoft using a custom file with only default settings).
        """
        return "resoft"

    def load_str(self, modname):
        """Returns the shell command to add the keyword as a string."""
        return "soft add "+modname

    def unload_str(self, modname):
        """Returns the shell command to delete the keyword as a string."""
        return "soft delete "+modname


# Regex that matches an environment variable reference, putting the name in
# the "name" group of the match.
_env_re = re.compile(
    """\$                      # Initial "$"
       (?P<brace>\{)?          # Open brace if present
       (?P<name>[A-Za-z0-9_]+) # Variable name
       (?(brace)\})            # Close brace if necessary""", re.X)


def expand_env(string, env):
    """Expand a shell string with given environment.

    Arguments:
    string - String to expand as if in a shell.
    env - dict specifying environment variables.

    This is very limited; right now it only handles variable substitution,
    and only with "${FOO}" or "$FOO" style syntax. The variable name must
    contain only alphanumeric characters and "_". There is no way to escape
    a "$".

    Expansion is recursive; the output of this function is its fixed point.
    """
    def expand_func(match):
        """Given an environment variable match, return replacement text."""
        var_name = match.group("name")
        if var_name in env:
            return env[var_name]
        else:
            # If there's no match in the environment, return original text.
            return match.group(0)
    old_string = ""
    new_string = string
    while new_string != old_string:
        old_string = new_string
        new_string = _env_re.sub(expand_func, new_string)

    return new_string
