# Add path to scream libs
import sys, os
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"))

from utils import ensure_yaml # pylint: disable=no-name-in-module
ensure_yaml()
import yaml

###############################################################################
# These are types that we use to differentiate lists of ints,bools,floats,strings
# We can use these types to tell YAML how to write them to file, which ultimately
# means to simply add the proper tag to the yaml file
###############################################################################
class Array(list):
    def __init__ (self, vals, t):
        super().__init__(t(v) for v in vals)
class Bools(Array):
    def __init__ (self,vals):
        Array.__init__(self,vals,bool)
class Ints(Array):
    def __init__ (self,vals):
        Array.__init__(self,vals,int)
class Floats(Array):
    def __init__ (self,vals):
        Array.__init__(self,vals,float)
class Strings(Array):
    def __init__ (self,vals):
        Array.__init__(self,vals,str)

###############################################################################
def make_array (vals,etype):
###############################################################################
    if etype=="bool" or etype=="logical":
        return Bools(vals)
    elif etype=="int" or etype=="integer":
        return Ints(vals)
    elif etype=="float" or etype=="real":
        return Floats(vals)
    elif etype=="string" or etype=="file":
        return Strings(vals)
    else:
        raise ValueError (f"Unsupported element type '{etype}' for arrays.")

###############################################################################
def array_constructor(loader: yaml.SafeLoader, node: yaml.nodes.SequenceNode) -> list:
###############################################################################
    entries = loader.construct_sequence(node)
    if node.tag=="!bools":
        return Bools(entries)
    elif node.tag=="!ints":
        return Ints(entries)
    elif node.tag=="!floats":
        return Floats(entries)
    elif node.tag=="!strings":
        return Strings(entries)
    else:
        raise ValueError(f"Invalid node tag={node.tag} for array constructor.")

###############################################################################
def array_representer(dumper,array) -> yaml.nodes.SequenceNode:
###############################################################################
    if isinstance(array,Bools):
        return dumper.represent_sequence('!bools',array)
    elif isinstance(array,Ints):
        return dumper.represent_sequence('!ints',array)
    elif isinstance(array,Floats):
        return dumper.represent_sequence('!floats',array)
    elif isinstance(array,Strings):
        return dumper.represent_sequence('!strings',array)
    else:
        raise ValueError (f"Unsupported array type: {type(array)}")

