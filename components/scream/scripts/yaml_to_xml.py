import os, sys, re
from collections import OrderedDict

from utils import ensure_yaml, expect
ensure_yaml()
import yaml

SWITCH_RE   = re.compile(r'<\s*([^:]+)\s*:')

###############################################################################
def ordered_load(item, Loader=yaml.SafeLoader, object_pairs_hook=OrderedDict):
###############################################################################
    """
    Copied from: https://stackoverflow.com/a/21912744
    Added ability to pass filename
    """
    class OrderedLoader(Loader):
        pass
    def construct_mapping(loader, node):
        loader.flatten_mapping(node)
        return object_pairs_hook(loader.construct_pairs(node))
    OrderedLoader.add_constructor(
        yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
        construct_mapping)

    if isinstance(item, str) and item.endswith(".yaml"):
        # Item is a filepath
        return yaml.load(open(item, "r"), OrderedLoader)
    else:
        return yaml.load(item, OrderedLoader)

###############################################################################
def do_switches(k, entry, indent):
###############################################################################
    m = SWITCH_RE.search(entry)
    while m:
        begin_switch_idx, end_cond_idx = m.span()
        expect(begin_switch_idx == 0, "Does not support nested switches")
        value = m.groups()[0].strip()
        value = value.lstrip("${").rstrip("}")
        opens = 1
        idx = end_cond_idx
        segments = [""]
        for idx in range(end_cond_idx, len(entry)):
            if entry[idx] == "<":
                opens += 1
            elif entry[idx] == ">" and entry[idx-1] != "=":
                opens -= 1
                if opens == 0:
                    break

            if entry[idx] == ":" and opens == 1:
                segments.append("")
            else:
                segments[-1] += entry[idx]

        expect(idx < len(entry), "Switch statement parse error in string '{}'".format(entry))
        end_switch_idx = idx + 1
        expect(end_switch_idx == len(entry), "Does not support nested switches")

        for segment in reversed(segments):
            components = [item.strip() for item in segment.split("=>", 1)]
            if len(components) == 1:
                print(f'{" "*indent}<{k}>{components[0]}</{k}>')
            else:
                print(f'{" "*indent}<{k} {value}="{components[0]}">{components[1]}</{k}>')

        m = None

###############################################################################
def process_leaf(k, v, indent):
###############################################################################
    if SWITCH_RE.match(str(v)):
        do_switches(k, v, indent)
    else:
        print(f'{" "*indent}<{k}>{v}</{k}>')

###############################################################################
def process_node(node, itr, indent):
###############################################################################
    for k, v in itr:

        if isinstance(v, dict):
            print(f'{" "*indent}<{k}>')
            process_node(v, v.items(), indent+2)
            print(f"{' '*indent}<\{k}>")

        elif isinstance(v, list):
            print(f'{" "*indent}<{k}>{", ".join([str(item) for item in v])}</{k}>')
        else:
            process_leaf(k, v, indent)

###############################################################################
def yaml_to_xml(filepath, indent):
###############################################################################
    expect(os.path.isfile(filepath), f"{filepath} does not exist")

    # Load scream inputs from yaml
    scream_input = ordered_load(filepath)

    process_node(scream_input, scream_input.items(), indent)
