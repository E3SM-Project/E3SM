#!/usr/bin/env python
"""
This script is used to setup individual test cases. Available test cases
can be see using the list_testcases.py script.

Specifically, this script parses XML files that define cases (steps in test
cases) and driver scripts, and generates directories and scripts to run each
step in the process of creating a test case.

This script requires a setup configuration file. Configuration files are
specific to each core. Template configuration files for each core can be seen
in this directory named 'general.config.{core}'. Each core may have different
requirements as far as what is required within a configuration file.
"""

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import sys
import os
import fnmatch
import argparse
import xml.etree.ElementTree as ET
import subprocess
from six.moves import configparser
import textwrap
import netCDF4
import shutil
import errno

try:
    from collections import defaultdict
except ImportError:
    from utils import defaultdict


# *** Namelist setup functions *** # {{{
def generate_namelist_files(config_file, case_path, configs):  # {{{
    config_tree = ET.parse(config_file)
    config_root = config_tree.getroot()

    # Iterate over all namelists to be generated
    for namelists in config_root.iter('namelist'):
        # Determine the name of the namelist that will be generated
        try:
            namelist_file = '{}/{}'.format(case_path, namelists.attrib['name'])
        except KeyError:
            print("ERROR: <namelist> tag is missing the 'name' attribute")
            print("Exiting...")
            sys.exit(1)

        try:
            namelist_mode = namelists.attrib['mode']
        except KeyError:
            print("ERROR: <namelist> tag is missing the 'mode' attribute.")
            print("Exiting...")
            sys.exit(1)

        if not configs.has_option("namelists", namelist_mode):
            print("Error. Configuration file '{}' requires paths for "
                  "streams and namelist files for '{}' mode.".format(
                      config_file, namelist_mode))
            print("Exiting...")
            sys.exit(1)

        template_namelist = configs.get("namelists", namelist_mode)

        # Ingest namelist template into a dictionary
        namelist_dict = defaultdict(lambda: defaultdict(list))
        ingest_namelist(template_namelist, namelist_dict)

        # Modify the dictionary to have the desired values
        configure_namelist(namelist_dict, namelists, configs)

        # Write the namelist using the template to determine the writing order.
        write_namelist(namelist_dict, namelist_file, template_namelist)
        del namelist_dict

    del config_root
    del config_tree
# }}}


def ingest_namelist(namelist_file, namelist_dict):  # {{{
    # Read the template file
    namelistfile = open(namelist_file, 'r')
    lines = namelistfile.readlines()

    record_name = 'NONE!!!'

    # Add each linke into the corresponding record / option entry in the
    # dictionary.
    for line in lines:
        if line.find('&') >= 0:
            record_name = line.strip().strip('&').strip('\n')
            namelist_dict[record_name] = defaultdict(list)
        elif line.find('=') >= 0:
            opt, val = line.strip().strip('\n').split('=')
            if record_name != "NONE!!!":
                namelist_dict[record_name][opt].append(val)
# }}}


def set_namelist_val(namelist_dict, option_name, option_val):  # {{{
    # Set the value of the namelist option.
    for record, opts in namelist_dict.items():
        for opt, val in opts.items():
            if opt.strip() == option_name:
                val[0] = option_val
# }}}


def configure_namelist(namelist_dict, namelist_tag, configs):  # {{{
    # Iterate over all children within the namelist tag.
    for child in namelist_tag:
        # Process <option> tags
        if child.tag == 'option':
            option_name = child.attrib['name']
            option_val = child.text
            set_namelist_val(namelist_dict, option_name, option_val)
        # Process <template> tags
        elif child.tag == 'template':
            apply_namelist_template(namelist_dict, child, configs)

# }}}


def apply_namelist_template(namelist_dict, template_tag, configs):  # {{{
    # Determine the template information, like it's path and the filename
    template_info = get_template_info(template_tag, configs)

    # Build the full filename for the template
    template_file = '{}/{}'.format(template_info['template_path'],
                                   template_info['template_file'])

    # Parse the template
    template_tree = ET.parse(template_file)
    template_root = template_tree.getroot()

    # Apply the template, by changing each option
    for child in template_root:
        if child.tag == 'namelist':
            for grandchild in child:
                if grandchild.tag == 'option':
                    option_name = grandchild.attrib['name']
                    option_val = grandchild.text
                    set_namelist_val(namelist_dict, option_name, option_val)
                elif grandchild.tag == 'template':
                    apply_namelist_template(namelist_dict, grandchild, configs)

    del template_root
    del template_tree
    del template_info
# }}}


def write_namelist(namelist_dict, outfilename, infilename):  # {{{
    # Write the namelist out, using the infilename as a template to determine
    # the writing order.
    in_namelist = open(infilename, 'r')
    lines = in_namelist.readlines()
    in_namelist.close()

    out_namelist = open(outfilename, 'w+')

    record_name = 'NONE!!!'

    for line in lines:
        if line.find('&') >= 0:
            if record_name != "NONE!!!":
                out_namelist.write('/\n')

            record_name = line.strip().strip('&').strip('\n')
            out_namelist.write(line)
        elif line.find('=') >= 0:
            opt, val = line.strip().strip('\n').split('=')
            if record_name != "NONE!!!":
                out_namelist.write('    {} = {}\n'.format(
                        opt.strip(),
                        namelist_dict[record_name][opt][0].strip()))

    if record_name != "NONE!!!":
        out_namelist.write('/\n')

    out_namelist.close()
# }}}
# }}}


# *** Streams setup functions *** # {{{

def generate_streams_files(config_file, case_path, configs):  # {{{
    config_tree = ET.parse(config_file)
    config_root = config_tree.getroot()

    # Iterate over all sterams files to be generated
    for streams in config_root:
        if streams.tag == "streams":
            # Determine the path to the template streams file
            streams_filename = '{}/{}'.format(case_path,
                                              streams.attrib['name'])

            try:
                streams_mode = streams.attrib['mode']
            except KeyError:
                print("ERROR: <streams> tag is missing the 'mode' attribute.")
                print("Exiting...")
                sys.exit(1)

            if not configs.has_option("streams", streams_mode):
                print("Error. Configuration file '{}' requires paths for "
                      "streams and namelist files for '{}' mode.".format(
                          config_file, streams_mode))
                print("Exiting...")
                sys.exit(1)

            template_streams = configs.get("streams", streams_mode)

            # Parse the template
            streams_tree = ET.parse(template_streams)
            streams_root = streams_tree.getroot()

            # Configure the new streams file, using the template as a starting
            # place.
            configure_streams_file(streams_root, streams, configs)

            # Write out the streams file
            write_streams_file(streams_root, config_file, streams_filename,
                               '{}'.format(case_path))

            del streams_root
            del streams_tree
# }}}


def flush_streams(streams, remove_mutable, remove_immutable):  # {{{
    if remove_mutable:
        # Remove all mutable streams from the template streams file
        for stream in streams.findall('stream'):
            streams.remove(stream)

    if remove_immutable:
        # Remove all immutable streams from the template streams file
        for stream in streams.findall('immutable_stream'):
            streams.remove(stream)
# }}}


def modify_stream_definition(streams_file, stream_conf):  # {{{
    # Determine the name of the stream to modify
    name_to_modify = stream_conf.attrib['name']

    found = False

    # Check if stream already exists:
    for stream in streams_file:
        if stream.tag == 'stream' or stream.tag == 'immutable_stream':
            name = stream.attrib['name']
            if name.strip() == name_to_modify.strip():
                if not found:
                    found = True
                    stream_to_modify = stream
                else:
                    print("ERROR: Stream {} found multiple times in "
                          "template. Exiting...".format(name.strip()))
                    sys.exit(1)

    # If not found, need to create it
    if not found:
        found = True
        stream_to_modify = ET.SubElement(streams_file, 'stream')
        stream_to_modify.set('name', name_to_modify)

    # Make all of the modifications from the config file
    for child in stream_conf:
        # Process attribute changes
        if child.tag == 'attribute':
            attr_name = child.attrib['name']
            attr_val = child.text
            stream_to_modify.set(attr_name, attr_val)
        # Process adding contents to the stream
        elif child.tag == 'add_contents':
            for member in child.findall('member'):
                member_name = member.attrib['name']
                member_type = member.attrib['type']
                sub_member = ET.SubElement(stream_to_modify, member_type)
                sub_member.set('name', member_name)
                if 'packages' in member.attrib.keys():
                    member_packages = member.attrib['packages']
                    sub_member.set('packages', member_packages)
        # Process removing contents from the stream
        elif child.tag == 'remove_contents':
            for member in child.findall('member'):
                member_name = member.attrib['name']
                for child in stream_to_modify.iter('*'):
                    try:
                        if child.attrib['name'] == member_name:
                            stream_to_modify.remove(child)
                    except KeyError:
                        print("   --- Tag: {} is missing a name "
                              "attribute".format(child.tag))

# }}}


def configure_streams_file(streams_file, streams_tag, configs):  # {{{
    keep_mode = streams_tag.attrib['keep']
    remove_immutable = False
    remove_mutable = False

    if keep_mode.strip() == 'immutable':
        remove_mutable = True
    if keep_mode.strip() == 'mutable':
        remove_immutable = True
    if keep_mode.strip() == 'none':
        remove_mutable = True
        remove_immutable = True

    # Flush requested streams
    flush_streams(streams_file, remove_mutable, remove_immutable)

    # Process all stream modifications
    for child in streams_tag:
        # Process all templates
        if child.tag == 'template':
            apply_stream_template(streams_file, child, configs)
        # Process stream definitions / modifications
        elif child.tag == 'stream':
            modify_stream_definition(streams_file, child)
# }}}


def apply_stream_template(streams_file, template_tag, configs):  # {{{
    # Determine template information, like path and filename
    template_info = get_template_info(template_tag, configs)

    # Build full path to template file
    template_file = '{}/{}'.format(template_info['template_path'],
                                   template_info['template_file'])

    # Parse the template
    template_tree = ET.parse(template_file)
    template_root = template_tree.getroot()

    # Apply the streams portion of the template to the streams file
    for child in template_root:
        if child.tag == 'streams':
            for grandchild in child:
                if grandchild.tag == 'stream':
                    modify_stream_definition(streams_file, grandchild)
                elif grandchild.tag == 'template':
                    apply_stream_template(streams_file, grandchild, configs)

    del template_tree
    del template_root
    del template_info
# }}}


def write_streams_file(streams, config_file, filename, init_path):  # {{{
    config_tree = ET.parse(config_file)
    config_root = config_tree.getroot()

    stream_file = open(filename, 'w')

    stream_file.write('<streams>\n')

    # Write out all immutable streams first
    for stream in streams.findall('immutable_stream'):
        stream_name = stream.attrib['name']

        stream_file.write('\n')
        stream_file.write('<immutable_stream name="{}"'.format(stream_name))
        # Process all attributes on the stream
        for attr, val in stream.attrib.items():
            if attr.strip() != 'name':
                stream_file.write('\n                  {}="{}"'.format(attr,
                                                                       val))

        stream_file.write('/>\n')

    # Write out all immutable streams
    for stream in streams.findall('stream'):
        stream_name = stream.attrib['name']

        stream_file.write('\n')
        stream_file.write('<stream name="{}"'.format(stream_name))

        # Process all attributes
        for attr, val in stream.attrib.items():
            if attr.strip() != 'name':
                stream_file.write('\n        {}="{}"'.format(attr, val))

        stream_file.write('>\n\n')

        # Write out all streams included in this stream
        for substream in stream.findall('stream'):
            substream_name = substream.attrib['name']
            if 'packages' in substream.attrib.keys():
                package_name = substream.attrib['packages']
                entry = '    <stream name="{}"'.format(substream_name) + \
                        ' packages="{}" '.format(package_name) + '/>\n'
            else:
                entry = '    <stream name="{}"'.format(substream_name) + '/>\n'
            stream_file.write(entry)

        # Write out all var_structs included in this stream
        for var_struct in stream.findall('var_struct'):
            var_struct_name = var_struct.attrib['name']
            if 'packages' in var_struct.attrib.keys():
                package_name = var_struct.attrib['packages']
                entry = '    <var_struct name="{}"'.format(var_struct_name) + \
                        ' packages="{}" '.format(package_name) + '/>\n'
            else:
                entry = '    <var_struct name="{}"'.format(var_struct_name) + \
                        '/>\n'
            stream_file.write(entry)

        # Write out all var_arrays included in this stream
        for var_array in stream.findall('var_array'):
            var_array_name = var_array.attrib['name']
            if 'packages' in var_array.attrib.keys():
                package_name = var_array.attrib['packages']
                entry = '    <var_array name="{}"'.format(var_array_name) + \
                        ' packages="{}" '.format(package_name) + '/>\n'
            else:
                entry = '    <var_array name="{}"'.format(var_array_name) + \
                        '/>\n'
            stream_file.write(entry)

        # Write out all vars included in this stream
        for var in stream.findall('var'):
            var_name = var.attrib['name']
            if 'packages' in var.attrib.keys():
                package_name = var.attrib['packages']
                entry = '    <var name="{}"'.format(var_name) + \
                        ' packages="{}" '.format(package_name) + '/>\n'
            else:
                entry = '    <var name="{}"'.format(var_name) + '/>\n'
            stream_file.write(entry)

        stream_file.write('</stream>\n')

    stream_file.write('\n')
    stream_file.write('</streams>\n')

    del config_tree
    del config_root
# }}}
# }}}


# *** Script Generation Functions *** # {{{
def generate_run_scripts(config_file, init_path, configs):  # {{{
    config_tree = ET.parse(config_file)
    config_root = config_tree.getroot()
    dev_null = open('/dev/null', 'r+')

    for run_script in config_root:
        # Process run_script
        if run_script.tag == 'run_script':
            # Determine the name of the script, and create the file
            script_name = run_script.attrib['name']
            script_path = "{}/{}".format(init_path, script_name)
            script = open(script_path, "w")

            # Write the script header
            script.write("#!/usr/bin/env python\n")
            script.write("\n")
            script.write("# This script was generated from "
                         "setup_testcases.py as part of a config file\n")
            script.write("\n")
            script.write('import sys\n')
            script.write('import os\n')
            script.write('import shutil\n')
            script.write('import glob\n')
            script.write("import subprocess\n\n\n")
            script.write("dev_null = open('/dev/null', 'w')\n")

            # Process each part of the run script
            for child in run_script:
                # Process each <step> tag
                if child.tag == 'step':
                    process_script_step(child, configs, '', script)
                # Process each <define_env_var> tag
                elif child.tag == 'define_env_var':
                    process_env_define_step(child, configs, '', script)
                elif child.tag == 'model_run':
                    process_model_run_step(child, configs, script)

            # Finish writing the script
            script.close()

            # Make the script executable
            subprocess.check_call(['chmod', 'a+x', '{}'.format(script_path)],
                                  stdout=dev_null, stderr=dev_null)

    dev_null.close()
    del config_tree
    del config_root
# }}}


def generate_driver_scripts(config_file, configs):  # {{{
    config_tree = ET.parse(config_file)
    config_root = config_tree.getroot()
    dev_null = open('/dev/null', 'r+')

    # init_path is where the driver script will live after it's generated.
    init_path = '{}/{}'.format(config.get('script_paths', 'work_dir'),
                               config.get('script_paths', 'config_path'))

    # Ensure we're in a <driver_script> tag
    if config_root.tag == 'driver_script':
        name = config_root.attrib['name']

        # Ensure work_dir exists before writing driver script there.
        if not os.path.exists(init_path):
            os.makedirs(init_path)

        link_load_compass_env(init_path, configs)

        # Create script file
        script = open('{}/{}'.format(init_path, name), 'w')

        # Write script header
        script.write('#!/usr/bin/env python\n')
        script.write('"""\n')
        script.write('This script was generated as part of a driver_script '
                     'file by the\nsetup_testcases.py script.\n')
        script.write('"""\n')
        script.write('import sys\n')
        script.write('import os\n')
        script.write('import shutil\n')
        script.write('import glob\n')
        script.write('import subprocess\n')
        script.write('import argparse\n')
        script.write('\n\n')
        script.write('# This script was generated by setup_testcases.py as '
                     'part of a driver_script\n'
                     '# file.\n')
        script.write("os.environ['PYTHONUNBUFFERED'] = '1'\n")
        script.write('parser = argparse.ArgumentParser(\n'
                     '        description=__doc__, '
                     'formatter_class=argparse.RawTextHelpFormatter)\n')

        case_dict = dict()
        for child in config_root:
            if child.tag == 'case':
                case_name = child.attrib['name']
                case_dict[case_name] = '1'
            if child.tag == 'template':
                print(" WARNING: use of templates outside of a case in a "
                      "driver_script is not supported!")
                print("          (name of template file is {})".format(
                    child.attrib['file']))
        for case_name in case_dict.keys():
            script.write('parser.add_argument("--no_{}", dest="no_{}",\n'
                         '                    help="If set, {} case will not '
                         'be run during "\n'
                         '                         "execution of this script'
                         '.",\n'
                         '                    action="store_true")\n'.format(
                             case_name, case_name, case_name))
            script.write('parser.add_argument("--finalize_{}", '
                         'dest="finalize_{}",\n'
                         '                    help="If set, {} case will have '
                         'symlinks replaced "\n'
                         '                         "with the files they point '
                         'to, this occurs after any "\n'
                         '                         "case runs that have been '
                         'requested.",\n'
                         '                    action="store_true")\n'.format(
                             case_name, case_name, case_name))

        script.write('\n')
        script.write('args = parser.parse_args()\n')
        script.write('base_path = os.getcwd()\n')
        script.write("dev_null = open('/dev/null', 'w')\n")
        script.write('error = False\n')
        script.write('\n')

        # Process children of driver_script
        for child in config_root:
            # Process each case, by changing into that directory, and
            # processing each step / define_env_var tag within it.
            if child.tag == 'case':
                case = child.attrib['name']
                script.write('if not args.no_{}:\n'.format(case))
                script.write('    os.chdir(base_path)\n')
                script.write('    os.chdir(' + "'{}')\n".format(case))
                # Process children of <case> tag
                for grandchild in child:
                    # Process <step> tags
                    if grandchild.tag == 'step':
                        process_script_step(grandchild, configs, '    ',
                                            script)
                    # Process <define_env_var> tags
                    elif grandchild.tag == 'define_env_var':
                        process_env_define_step(grandchild, configs, '    ',
                                                script)
            # Process <step> tags
            elif child.tag == 'step':
                script.write('os.chdir(base_path)\n')
                process_script_step(child, configs, '', script)
            # Process <compare_fields> tags
            elif child.tag == 'validation':
                script.write('os.chdir(base_path)\n')
                process_validation_step(child, configs, script)
            # Process <define_env_var> tags
            elif child.tag == 'define_env_var':
                script.write('os.chdir(base_path)\n')
                process_env_define_step(child, configs, '', script)

        # Write script footer, that ensures a 1 is returned if the script
        # encountered an error. This happens before finalizing a case
        # directory.
        script.write('if error:\n')
        script.write('    sys.exit(1)\n')

        for case_name in case_dict.keys():
            script.write('if args.finalize_{}:\n'.format(case_name))
            script.write('    old_dir = os.getcwd()\n')
            script.write('    os.chdir("{}")\n'.format(case_name))
            script.write('    file_list = glob.glob("*")\n')
            script.write('    for file in file_list:\n')
            script.write('        if os.path.islink(file):\n')
            script.write('            link_path = os.readlink(file)\n')
            script.write('            os.unlink(file)\n')
            script.write('            shutil.copyfile(link_path, file)\n')
            script.write('    os.chdir(old_dir)\n')
            script.write('\n')

        script.write('sys.exit(0)\n')
        script.close()
        del case_dict

        # Make script executable
        subprocess.check_call(['chmod', 'a+x',
                               '{}/{}'.format(init_path, name)],
                              stdout=dev_null, stderr=dev_null)
# }}}


def process_env_define_step(var_tag, configs, indentation, script_file):  # {{{
    try:
        var_name = var_tag.attrib['name']
    except KeyError:
        print("ERROR: <define_env_var> tag is missing 'name' attribte")
        print('Exiting...')
        sys.exit(1)

    try:
        var_val = var_tag.attrib['value']
    except KeyError:
        print("ERROR: <define_env_var> tag is missing 'value' attribute")
        print('Exiting...')
        sys.exit(1)

    # Write line to define the environment variable
    script_file.write("{}os.environ['{}'] = '{}'\n".format(indentation,
                                                           var_name,
                                                           var_val))
# }}}


def process_script_step(step, configs, indentation, script_file):  # {{{
    # Determine step attributes.
    if 'executable_name' in step.attrib.keys() and 'executable' in \
            step.attrib.keys():
        print("ERROR: <step> tag has both an 'executable' and "
              "'executable_name' attribute. Only one is allowed per step.")
        print("Exiting...")
        sys.exit(1)

    try:
        quiet_val = step.attrib['quiet']
        if quiet_val == "true":
            quiet = True
        else:
            quiet = False
    except KeyError:
        quiet = False

    try:
        step_pre_message = step.attrib['pre_message']
        write_pre_message = True
    except KeyError:
        write_pre_message = False

    try:
        step_post_message = step.attrib['post_message']
        write_post_message = True
    except KeyError:
        write_post_message = False

    try:
        executable_name = step.attrib['executable_name']
        executable = configs.get('executables', executable_name)
    except KeyError:
        executable = step.attrib['executable']

    # Write step header
    script_file.write("\n")

    # If a pre_message attribute was supplied, write it before adding the
    # command.
    if write_pre_message:
        script_file.write('{}print("{}")\n'.format(indentation,
                                                   step_pre_message))

    script_file.write("{}# Run command is:\n".format(indentation))

    command_args = [executable]
    # Process step arguments
    for argument in step:
        if argument.tag == 'argument':
            flag = argument.attrib['flag']
            val = argument.text

            if flag.strip() != "":
                command_args.append(flag)

            if val is not None:
                command_args.append(val)

    # Build comment and command bases
    comment = wrap_subprocess_comment(command_args, indentation)
    command = wrap_subprocess_command(command_args, indentation, quiet)

    # Write the comment, and the command. Also, ensure the command has the same
    # environment as the calling script.
    script_file.write("{}\n".format(comment))
    script_file.write("{}\n".format(command))

    # If a post_message attribute was supplied, write it after the command.
    if write_post_message:
        script_file.write('{}print("{}")\n'.format(indentation,
                                                   step_post_message))

# }}}


def process_validation_step(validation_tag, configs, script):  # {{{
    for child in validation_tag:
        if child.tag == 'compare_fields':
            process_compare_fields_step(child, configs, script)
        if child.tag == 'compare_timers':
            process_compare_timers_step(child, configs, script)
# }}}


# *** Field Comparison Functions *** ##{{{
def process_compare_fields_step(compare_tag, configs, script):  # {{{
    missing_file1 = False
    missing_file2 = False
    # Determine comparison attributes
    try:
        file1 = compare_tag.attrib['file1']
    except KeyError:
        missing_file1 = True

    try:
        file2 = compare_tag.attrib['file2']
    except KeyError:
        missing_file2 = True

    if missing_file1 and missing_file2:
        print("ERROR: <compare_fields> tag is missing both 'file1' and "
              "'file2' tags. At least one is required.")
        print("Exiting...")
        sys.exit(1)

    baseline_root = configs.get('script_paths', 'baseline_dir')

    if baseline_root != 'NONE':
        baseline_root = '{}/{}'.format(baseline_root,
                                       configs.get('script_paths', 'test_dir'))

    for child in compare_tag:
        # Process field comparisons
        if child.tag == 'field':
            if not (missing_file1 or missing_file2):
                process_field_definition(child, configs, script, file1, file2,
                                         False)

            if not missing_file1 and baseline_root != 'NONE':
                process_field_definition(child, configs, script, file1,
                                         '{}/{}'.format(baseline_root, file1),
                                         True)

            if not missing_file2 and baseline_root != 'NONE':
                process_field_definition(child, configs, script, file2,
                                         '{}/{}'.format(baseline_root, file2),
                                         True)
        # Process field comparison template
        elif child.tag == 'template':
            apply_compare_fields_template(child, compare_tag, configs, script)
# }}}


def apply_compare_fields_template(template_tag, compare_tag, configs, script):
    # {{{
    missing_file1 = False
    missing_file2 = False
    # Determine comparison attributes
    try:
        file1 = compare_tag.attrib['file1']
    except KeyError:
        missing_file1 = True

    try:
        file2 = compare_tag.attrib['file2']
    except KeyError:
        missing_file2 = True

    if missing_file1 and missing_file2:
        print("ERROR: <compare_fields> tag is missing both 'file1' and "
              "'file2' tags. At least one is required.")
        print("Exiting...")
        sys.exit(1)

    # Build the path to the baselines
    baseline_root = configs.get('script_paths', 'baseline_dir')
    if baseline_root != 'NONE':
        baseline_root = '{}/{}'.format(baseline_root,
                                       configs.get('script_paths', 'test_dir'))

    # Determine template information, like path and filename
    template_info = get_template_info(template_tag, configs)

    template_file = '{}/{}'.format(template_info['template_path'],
                                   template_info['template_file'])

    # Parse the template
    template_tree = ET.parse(template_file)
    template_root = template_tree.getroot()

    # Find a child tag that is validation->compare_fields->field, and add each
    # field
    for validation in template_root:
        if validation.tag == 'validation':
            for compare_fields in validation:
                if compare_fields.tag == 'compare_fields':
                    for field in compare_fields:
                        if field.tag == 'field':
                            if not (missing_file1 or missing_file2):
                                process_field_definition(field, configs,
                                                         script, file1, file2,
                                                         False)

                            if not missing_file1 and baseline_root != 'NONE':
                                process_field_definition(
                                    field, configs, script, file1,
                                    '{}/{}'.format(baseline_root, file1), True)

                            if not missing_file2 and baseline_root != 'NONE':
                                process_field_definition(
                                    field, configs, script, file2,
                                    '{}/{}'.format(baseline_root, file2), True)
                        elif field.tag == 'template':
                            apply_compare_fields_template(field, compare_tag,
                                                          configs, script)

    del template_root
    del template_tree
    del template_info
# }}}


def process_field_definition(field_tag, configs, script, file1, file2,
                             baseline_comp):  # {{{
    # Build the path to the comparison script.
    compare_executable = '{}/compare_fields.py'.format(
        configs.get('script_paths', 'utility_scripts'))

    field_name = field_tag.attrib['name']

    # Build the base command to compare the fields
    command_args = [compare_executable, '-q', '-1', file1, '-2', file2, '-v',
                    field_name]

    # Determine norm thresholds
    if baseline_comp:
        command_args.extend(['--l1', '0.0', '--l2', '0.0', '--linf', '0.0'])
    else:
        if 'l1_norm' in field_tag.attrib.keys():
            command_args.extend(['--l1', field_tag.attrib['l1_norm']])
        if 'l2_norm' in field_tag.attrib.keys():
            command_args.extend(['--l2', field_tag.attrib['l2_norm']])
        if 'linf_norm' in field_tag.attrib.keys():
            command_args.extend(['--linf', field_tag.attrib['linf_norm']])

    command = wrap_subprocess_command(command_args, indentation='    ',
                                      quiet=False)

    # Write the pass/fail logic.
    script.write('try:\n')
    script.write('{}\n'.format(command))
    script.write("    print(' ** PASS Comparison of {} between {} and\\n'\n"
                 "          '    {}')\n".format(field_name, file1, file2))
    script.write('except subprocess.CalledProcessError:\n')
    script.write("    print(' ** FAIL Comparison of {} between {} and\\n'\n"
                 "          '    {}')\n".format(field_name, file1, file2))
    script.write('    error = True\n')
# }}}
# }}}


# *** Timer Comparison Functions *** # {{{
def process_compare_timers_step(compare_tag, configs, script):  # {{{
    baseline_root = configs.get('script_paths', 'baseline_dir')
    baseline_root = '{}/{}'.format(baseline_root,
                                   configs.get('script_paths', 'test_dir'))

    missing_rundir1 = True
    missing_rundir2 = True

    try:
        rundir1 = compare_tag.attrib['rundir1']
        missing_rundir1 = False
    except KeyError:
        missing_rundir1 = True

    try:
        rundir2 = compare_tag.attrib['rundir2']
        missing_rundir2 = False
    except KeyError:
        missing_rundir2 = True

    for child in compare_tag:
        if child.tag == 'timer':
            if not (missing_rundir1 or missing_rundir2):
                process_timer_definition(child, configs, script, rundir1,
                                         rundir2)

            if not missing_rundir1:
                process_timer_definition(
                    child, configs, script,
                    '{}/{}'.format(baseline_root, rundir1), rundir1)

            if not missing_rundir2:
                process_timer_definition(
                    child, configs, script,
                    '{}/{}'.format(baseline_root, rundir2), rundir2)
        elif child.tag == 'template':
            apply_compare_timers_template(child, compare_tag, configs, script)

# }}}


def apply_compare_timers_template(template_tag, compare_tag, configs, script):
    # {{{
    # Build the path to the baselines
    baseline_root = configs.get('script_paths', 'baseline_dir')
    baseline_root = '{}/{}'.format(baseline_root, configs.get('script_paths',
                                                              'test_dir'))

    missing_rundir1 = True
    missing_rundir2 = True

    try:
        rundir1 = compare_tag.attrib['rundir1']
        missing_rundir1 = False
    except KeyError:
        missing_rundir1 = True

    try:
        rundir2 = compare_tag.attrib['rundir2']
        missing_rundir2 = False
    except KeyError:
        missing_rundir2 = True

    # Get the template information and build the template file
    template_info = get_template_info(template_tag, configs)
    template_file = '{}/{}'.format(template_info['template_path'],
                                   template_info['template_file'])

    # Parse template file
    template_tree = ET.parse(template_file)
    template_root = template_tree.getroot()

    for validation in template_root:
        if validation.tag == 'validation':
            for compare_timers in validation:
                if compare_timers.tag == 'compare_timers':
                    for timer in compare_timers:
                        if timer.tag == 'timer':
                            if not (missing_rundir1 or missing_rundir2):
                                process_timer_definition(
                                    timer, configs, script, rundir1, rundir2)

                            if not missing_rundir1:
                                process_timer_definition(
                                    timer, configs, script,
                                    '{}/{}'.format(baseline_root, rundir1),
                                    rundir1)

                            if not missing_rundir2:
                                process_timer_definition(
                                    timer, configs, script,
                                    '{}/{}'.format(baseline_root, rundir2),
                                    rundir2)
                        elif timer.tag == 'template':
                            apply_compare_timers_template(timer, compare_tag,
                                                          configs, script)
    del template_root
    del template_tree
    del template_info

# }}}


def process_timer_definition(timer_tag, configs, script, basedir, compdir):
    # {{{
    compare_script = '{}/compare_timers.py'.format(
        configs.get('script_paths', 'utility_scripts'))

    try:
        timer_name = timer_tag.attrib['name']
    except KeyError:
        print("ERROR: <timer> tag is missing the 'name' attribute.")
        print("Exiting...")
        sys.exit(1)

    command = 'subprocess.check_call(["{}", "-b", "{}", "-c", "{}", "-t", ' \
              '"{}"])'.format(compare_script, compdir, basedir, timer_name)

    script.write('\n')
    script.write('if os.path.exists("{}") and \\\n'
                 '        os.path.exists("{}"):\n'.format(compdir, basedir))
    script.write('    try:\n')
    script.write('        {}\n'.format(command))
    script.write("        print(' ** PASS Comparison of timer {} between {} "
                 "and\\n'\n"
                 "              '    {}')\n".format(timer_name, compdir,
                                                    basedir))
    script.write('    except subprocess.CalledProcessError:\n')
    script.write("        print(' ** FAIL Comparison of timer {} between {} "
                 "and\\n'\n"
                 "              '    {}')\n".format(timer_name, compdir,
                                                    basedir))
    script.write("        error = True\n")
# }}}
# }}}


def process_model_run_step(model_run_tag, configs, script):  # {{{
    run_definition_file = configs.get('script_input_arguments',
                                      'model_runtime')
    run_config_tree = ET.parse(run_definition_file)
    run_config_root = run_config_tree.getroot()

    dev_null = open('/dev/null', 'r+')

    try:
        executable_name = model_run_tag.attrib['executable']
    except KeyError:
        executable_name = 'model'

    script.write('print("\\n")\n')
    script.write('print("     *****************************")\n')
    script.write('print("     ** Starting model run step **")\n')
    script.write('print("     *****************************")\n')
    script.write('print("\\n")\n')

    # Process each part of the run script
    for child in run_config_root:
        # Process each <step> tag
        if child.tag == 'step':
            # Setup child step, and it's attributes to be correct for the
            # process_script_step function
            for grandchild in child:
                if grandchild.tag == 'argument':
                    arg_text = grandchild.text

                    if arg_text == 'model':
                        executable_full_path = config.get('executables',
                                                          executable_name)
                        executable_parts = executable_full_path.split('/')
                        executable_link = \
                            executable_parts[len(executable_parts) - 1]
                        link_path = '{}/{}/{}'.format(
                            config.get('script_paths', 'work_dir'),
                            config.get('script_paths', 'case_dir'),
                            executable_link)
                        subprocess.check_call(
                            ['ln', '-sf',
                             config.get('executables', executable_name),
                             link_path],
                            stdout=dev_null, stderr=dev_null)
                        grandchild.text = './{}'.format(executable_link)
                    elif arg_text.find('attr_') >= 0:
                        attr_array = arg_text.split('_')
                        try:
                            grandchild.text = \
                                model_run_tag.attrib[attr_array[1]]
                        except KeyError:
                            print(" <step> tag defined within a <model_run> "
                                  "tag requires attribute '{}', but it is "
                                  "not defined.".format(attr_array[1]))
                            print(" Exiting...")
                            sys.exit(1)

            # Process the resulting element, instead of the original step.
            process_script_step(child, configs, '', script)
        # Process each <define_env_var> tag
        elif child.tag == 'define_env_var':
            if child.attrib['value'].find('attr_') >= 0:
                attr_array = child.attrib['value'].split('_')
                try:
                    child.attrib['value'] = model_run_tag.attrib[attr_array[1]]
                except KeyError:
                    print(" <define_env_var> tag defined within a "
                          "<model_run> tag requires attribute '{}', but it "
                          "is not defined.".format(attr_array[1]))
                    print(" Exiting...")
                    sys.exit(1)

            process_env_define_step(child, configs, '', script)

    script.write('print("\\n")\n')
    script.write('print("     *****************************")\n')
    script.write('print("     ** Finished model run step **")\n')
    script.write('print("     *****************************")\n')
    script.write('print("\\n")\n')
    dev_null.close()
# }}}


def wrap_subprocess_comment(command_args, indentation):  # {{{
    # Build comment and command bases
    comment = textwrap.fill(' '.join(command_args), width=79,
                            initial_indent="{}# ".format(indentation),
                            subsequent_indent="{}# ".format(indentation),
                            break_on_hyphens=False, break_long_words=False)
    return comment
# }}}


def wrap_subprocess_command(command_args, indentation, quiet):  # {{{
    # Setup command redirection
    if quiet:
        redirect = ", stdout=dev_null, stderr=None"
    else:
        redirect = ""

    prefix = "{}subprocess.check_call(".format(indentation)
    command = textwrap.wrap("'{}'".format("', '".join(command_args)), width=79,
                            initial_indent="{}[".format(prefix),
                            subsequent_indent=' ' * (len(prefix)+1),
                            break_on_hyphens=False, break_long_words=False)

    last_line = command.pop()
    command.extend(textwrap.wrap(
            "{}]{})".format(last_line, redirect),
            width=80, subsequent_indent=' ' * (len(prefix)),
            break_on_hyphens=False, break_long_words=False))
    command = '\n'.join(command)
    return command
# }}}
# }}}


# *** General Utility Functions *** #{{{
def add_links(config_file, configs):  # {{{
    config_tree = ET.parse(config_file)
    config_root = config_tree.getroot()

    case = config_root.attrib['case']

    dev_null = open('/dev/null', 'r+')

    # Determine the path for the case directory
    test_path = '{}/{}'.format(configs.get('script_paths', 'test_dir'), case)
    base_path = '{}/{}'.format(configs.get('script_paths', 'work_dir'),
                               test_path)

    # Process all children tags
    for child in config_root:
        # Process an <add_link> tag
        if child.tag == 'add_link':
            source_file = get_source_file(child, configs)

            dest = child.attrib['dest']
            old_cwd = os.getcwd()
            os.chdir(base_path)

            subprocess.check_call(['ln', '-sfn', '{}'.format(source_file),
                                   '{}'.format(dest)],
                                  stdout=dev_null, stderr=dev_null)
            os.chdir(old_cwd)
        # Process an <add_executable> tag
        elif child.tag == 'add_executable':
            source_attr = child.attrib['source']
            dest = child.attrib['dest']
            if not configs.has_option("executables", source_attr):
                raise ValueError('Configuration {} requires a definition of '
                      '{}.'.format(config_file, source_attr))
            source = configs.get("executables", source_attr)

            subprocess.check_call(['ln', '-sf', '{}'.format(source),
                                   '{}/{}'.format(base_path, dest)],
                                  stdout=dev_null, stderr=dev_null)

    dev_null.close()
# }}}


def get_source_file(child, config):  # {{{
    try:
        source = child.attrib['source']
    except KeyError:
        raise KeyError("{} tag missing a 'source' attribute.".format(
            child.tag))

    try:
        source_path_name = child.attrib['source_path']
    except KeyError:
        return source

    keyword_path = any(substring in source_path_name for
                       substring in ['work_', 'script_'])

    if keyword_path:
        source_arr = source_path_name.split('_')
        base_name = source_arr[0]
        if base_name == 'work':
            file_base_path = 'work_dir'
        elif base_name == 'script':
            file_base_path = 'script_path'
        else:
            raise ValueError('Unexpected source prefix {} in {} tag'.format(
                base_name, child.tag))

        subname = '{}_{}'.format(source_arr[1], source_arr[2])
        if subname not in ['core_dir', 'configuration_dir',
                           'resolution_dir', 'test_dir', 'case_dir']:
            raise ValueError('Unexpected source suffix {} in {} tag'.format(
                subname, child.tag))

        source_path = '{}/{}'.format(
            config.get('script_paths', file_base_path),
            config.get('script_paths', subname))
    else:
        if config.has_option('paths', source_path_name):
            source_path = config.get('paths', source_path_name)
        else:
            if not config.has_option('script_paths', source_path_name):
                raise ValueError('Undefined source_path on {} tag: {}'.format(
                    child.tag, source_path_name))
            source_path = config.get('script_paths',  source_path_name)

    source_file = '{}/{}'.format(source_path, source)
    return source_file
# }}}


def copy_files(config_file, config):  # {{{
    config_tree = ET.parse(config_file)
    config_root = config_tree.getroot()

    case = config_root.attrib['case']

    # Determine the path for the case directory
    test_path = '{}/{}'.format(config.get('script_paths', 'test_dir'), case)
    base_path = '{}/{}'.format(config.get('script_paths', 'work_dir'),
                               test_path)

    # Process all children tags
    for child in config_root:
        # Process an <copy_file> tag
        if child.tag == 'copy_file':
            source = get_source_file(child, config)

            dest = '{}/{}'.format(base_path, child.attrib['dest'])

            shutil.copy(source, dest)
# }}}


def make_case_dir(config_file, base_path):  # {{{
    config_tree = ET.parse(config_file)
    config_root = config_tree.getroot()

    case_name = config_root.attrib['case']

    # Build the case directory, if it doesn't already exist
    if not os.path.exists('{}/{}'.format(base_path, case_name)):
        os.makedirs('{}/{}'.format(base_path, case_name))

    del config_root
    del config_tree

    return case_name
# }}}


def get_defined_files(config_file, init_path, configs):  # {{{
    config_tree = ET.parse(config_file)
    config_root = config_tree.getroot()
    dev_null = open('/dev/null', 'w')

    for get_file in config_root:
        # Process <get_file> tag
        if get_file.tag == 'get_file':
            # Determine dest_path
            try:
                dest_path_name = get_file.attrib['dest_path']
            except KeyError:
                print(" get_file tag is missing the 'dest_path' attribute.")
                print(" Exiting...")
                sys.exit(1)

            # Determine file_name
            try:
                file_name = get_file.attrib['file_name']
            except KeyError:
                print(" get_file tag is missing a 'file_name' attribute.")
                print(" Exiting...")
                sys.exit(1)

            # Build out the dest path
            keyword_path = False
            if dest_path_name.find('work_') >= 0:
                keyword_path = True
            elif dest_path_name.find('script_') >= 0:
                keyword_path = True
            else:
                if configs.has_option('paths', dest_path_name):
                    dest_path = '{}'.format(
                        configs.get('paths', dest_path_name))
                else:
                    print(" Path '{}' is not defined in the config file, "
                          "but is required to get a file.".format(
                              dest_path_name))
                    print(" Exiting...")
                    sys.exit(1)

            if keyword_path:
                dest_arr = dest_path_name.split('_')
                base_name = dest_arr[0]
                subname = '{}_{}'.format(dest_arr[1], dest_arr[2])

                if base_name == 'work':
                    base_path = 'work_dir'
                elif base_name == 'script':
                    base_path = 'script_path'

                if subname in {'core_dir', 'configuration_dir',
                               'resolution_dir', 'test_dir', 'case_dir'}:
                    dest_path = '{}/{}'.format(
                        configs.get('script_paths', base_path),
                        configs.get('script_paths', subname))
                else:
                    print(" Path '{}' is not defined.".format(dest_path_name))
                    print(" Exiting...")
                    sys.exit(1)

            # if the dest_path doesn't exist, create it
            if not os.path.exists(dest_path):
                os.makedirs(dest_path)

            # If the file doesn't exist in dest_path, process it's mirrors
            if not os.path.exists('{}/{}'.format(dest_path, file_name)):
                file_found = False
                if not file_found:
                    for mirror in get_file:
                        # Process each mirror
                        if mirror.tag == 'mirror':
                            # Determine the protocol for the mirror
                            try:
                                protocol = mirror.attrib['protocol']
                            except KeyError:
                                print("Mirror is missing the 'protocol' "
                                      "attribute.")
                                print("Exiting...")
                                sys.exit(1)

                            # Process a wget mirror
                            if protocol == 'wget':
                                if config.get('script_input_arguments',
                                              'no_download') == 'no':
                                    try:
                                        path = '{}/{}'.format(
                                            mirror.attrib['url'], file_name)
                                    except KeyError:
                                        print(" Mirror with protocol 'wget' "
                                              "is missing a 'url' attribute")
                                        print(" Exiting...")
                                        sys.exit(1)

                                    try:
                                        subprocess.check_call(
                                            ['wget', '-q', '{}'.format(path)],
                                            stdout=dev_null, stderr=dev_null)
                                    except subprocess.CalledProcessError:
                                        print(" Wget failed....")

                                    try:
                                        subprocess.check_call(
                                            ['mv', '{}'.format(file_name),
                                             '{}/{}'.format(dest_path,
                                                            file_name)],
                                            stdout=dev_null, stderr=dev_null)
                                    except subprocess.CalledProcessError:
                                        print("  -- Web mirror attempt failed."
                                              " Trying other mirrors...")

                        # If the file exists in dest_path, determine if it
                        # should be validated
                        if os.path.exists('{}/{}'.format(dest_path,
                                          file_name)):
                            try:
                                expected_hash = get_file.attrib['hash']
                                validate_file = True
                            except KeyError:
                                validate_file = False

                            # Validate the file, if requested (i.e. file had a
                            # hash attribute)
                            if validate_file:
                                if os.path.exists('{}/{}'.format(dest_path,
                                                                 file_name)):
                                    nc = netCDF4.Dataset('{}/{}'.format(
                                        dest_path, file_name), 'r')
                                    try:
                                        file_hash = nc.file_id
                                    except AttributeError:
                                        nc.close()
                                        print(" Downloaded file '{}' does "
                                              "not have a 'file_id' "
                                              "attribute.".format(file_name))
                                        print(" Deleting file and exiting...")
                                        subprocess.check_call(
                                            ['rm', '-f',
                                             '{}/{}'.format(dest_path,
                                                            file_name)],
                                            stdout=dev_null, stderr=dev_null)
                                        sys.exit(1)

                                    nc.close()

                                    if file_hash.strip() != \
                                            expected_hash.strip():
                                        print("*** ERROR: Base mesh has "
                                              "hash of '{}' which does not "
                                              "match expected hash of "
                                              "'{}'.".format(file_hash,
                                                             expected_hash))
                                        print(" Deleting file...")
                                        subprocess.check_call(
                                            ['rm', '-f',
                                             '{}/{}'.format(dest_path,
                                                            file_name)],
                                            stdout=dev_null, stderr=dev_null)

                    # IF validation valied, exit.
                    if not os.path.exists('{}/{}'.format(dest_path,
                                                         file_name)):
                        print(" Failed to acquire required file '{}'.".format(
                            file_name))
                        print(" Exiting...")
                        sys.exit(1)

    del config_tree
    del config_root
    dev_null.close()
# }}}


def get_template_info(template, configs):  # {{{
    # Determine template attributes
    try:
        template_file = template.attrib['file']
    except KeyError:
        print("ERROR: <template> tag is missing 'file' attribute.")
        print('Exiting...')
        sys.exit(1)

    try:
        template_path_base = template.attrib['path_base']
    except KeyError:
        print("ERROR: <template> tag is missing 'path_base' attribute.")
        print('Exiting...')
        sys.exit(1)

    try:
        template_path_modifier = '/{}'.format(template.attrib['path'])
    except KeyError:
        template_path_modifier = ''

    # Determine template path from it's path_base attributes
    if template_path_base.find('work_') >= 0:
        keyword_path = True
    elif template_path_base.find('script_') >= 0:
        keyword_path = True

    if keyword_path:
        template_arr = template_path_base.split('_')
        base = template_arr[0]
        if base == 'work':
            base_path = 'work_dir'
        elif base == 'script':
            base_path = 'script_path'

        sub_path = '{}_{}'.format(template_arr[1], template_arr[2])

        if sub_path in {'core_dir', 'configuration_dir', 'resolution_dir',
                        'test_dir'}:
            template_path = '{}/{}{}'.format(
                configs.get('script_paths', base_path),
                configs.get('script_paths', sub_path), template_path_modifier)
        else:
            print(" Template path_base '{}' does not exist".format(
                template_path_base))
            print(" Exiting...")
            sys.exit(1)

    info = {}
    info['template_file'] = template_file
    info['template_path'] = template_path

    return info
# }}}


def get_config_file_type(config_file):  # {{{
    config_tree = ET.parse(config_file)
    config_root = config_tree.getroot()

    # Determine file type
    file_type = config_root.tag

    del config_root
    del config_tree

    return file_type
# }}}


def get_case_name(config_file):  # {{{
    config_tree = ET.parse(config_file)
    config_root = config_tree.getroot()

    # Determine file type
    if config_root.tag == 'config':
        name = config_root.attrib['case']

    del config_root
    del config_tree

    return name
# }}}

def link_load_compass_env(init_path, configs):  # {{{

    if configs.getboolean('conda', 'link_load_compass'):
        target = '{}/{}/load_compass_env.sh'.format(
            configs.get('script_paths', 'script_path'),
            configs.get('script_paths', 'core_dir'))

        link_name = '{}/load_compass_env.sh'.format(init_path)
        try:
            os.symlink(target, link_name)
        except OSError as e:
            if e.errno == errno.EEXIST:
                os.remove(link_name)
                os.symlink(target, link_name)
            else:
                raise e
# }}}
# }}}


if __name__ == "__main__":
    # Define and process input arguments
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-o", "--core", dest="core",
                        help="Core that contains configurations",
                        metavar="CORE")
    parser.add_argument("-c", "--configuration", dest="configuration",
                        help="Configuration to setup", metavar="CONFIG")
    parser.add_argument("-r", "--resolution", dest="resolution",
                        help="Resolution of configuration to setup",
                        metavar="RES")
    parser.add_argument("-t", "--test", dest="test",
                        help="Test name within a resolution to setup",
                        metavar="TEST")
    parser.add_argument("-n", "--case_number", dest="case_num",
                        help="Case number to setup, as listed from "
                             "list_testcases.py. Can be a comma delimited "
                             "list of case numbers.", metavar="NUM")
    parser.add_argument("-f", "--config_file", dest="config_file",
                        help="Configuration file for test case setup",
                        metavar="FILE")
    parser.add_argument("-m", "--model_runtime", dest="model_runtime",
                        help="Definition of how to build model run commands "
                             "on this machine", metavar="FILE")
    parser.add_argument("-b", "--baseline_dir", dest="baseline_dir",
                        help="Location of baseslines that can be compared to",
                        metavar="PATH")
    parser.add_argument("-q", "--quiet", dest="quiet",
                        help="If set, script will not write a command_history "
                             "file", action="store_true")
    parser.add_argument("--no_download", dest="no_download",
                        help="If set, script will not auto-download base_mesh "
                             "files", action="store_true")
    parser.add_argument("--work_dir", dest="work_dir",
                        help="If set, script will create case directories in "
                             "work_dir rather than the current directory.",
                        metavar="PATH")
    parser.add_argument("--link_load_compass", dest="link_load_compass",
                        action="store_true",
                        help="If set, a link to <core>/load_compass_env.sh is "
                             "included with each test case")

    args = parser.parse_args()

    if not args.config_file:
        print("WARNING: No configuration file specified. Using the default "
              "of 'local.config'")
        args.config_file = 'local.config'

    if not os.path.exists(args.config_file):
        parser.error(" Configuration file '{}' does not exist. Please create "
                     "and setup before running again.".format(
                         args.config_file))

    if not args.case_num and not (args.core and args.configuration and
                                  args.resolution and args.test):
        print('Must be run with either the --case_number argument, or the '
              'core, configuration, resolution, and test arguments.')
        parser.error(' Invalid configuration. Exiting...')

    if args.case_num and args.core and args.configuration and args.resoltuion \
            and args.test:
        print('Can only be configured with either --case_number (-n) or '
              '--core (-o), --configuration (-c), --resolution (-r), and '
              '--test (-t).')
        parser.error(' Invalid configuration. Too many options used. '
                     'Exiting...')

    if args.case_num:
        use_case_list = True
        case_list = args.case_num.split(',')
    else:
        use_case_list = False
        case_list = list()
        case_list.append(0)

    if sys.version_info >= (3, 2):
        config = configparser.ConfigParser()
    else:
        config = configparser.SafeConfigParser()
    config.read(args.config_file)

    if not args.no_download:
        args.no_download = False

    if not args.work_dir:
        args.work_dir = os.getcwd()

    # Add configuation information to the config object.
    # This allows passing config around with all of the config options needed
    # to build paths, and determine options.
    config.add_section('script_input_arguments')
    config.add_section('script_paths')

    if not use_case_list:
        config.set('script_input_arguments', 'core', args.core)
        config.set('script_input_arguments', 'configuration',
                   args.configuration)
        config.set('script_input_arguments', 'resolution', args.resolution)
        config.set('script_input_arguments', 'test', args.test)

    if args.baseline_dir:
        config.set('script_paths', 'baseline_dir', args.baseline_dir)
    else:
        config.set('script_paths', 'baseline_dir', 'NONE')

    if args.no_download:
        config.set('script_input_arguments', 'no_download', 'yes')
    else:
        config.set('script_input_arguments', 'no_download', 'no')

    config.set('script_paths', 'script_path',
               os.path.dirname(os.path.realpath(__file__)))
    config.set('script_paths', 'work_dir', os.path.abspath(args.work_dir))
    config.set('script_paths', 'utility_scripts',
               '{}/utility_scripts'.format(config.get('script_paths',
                                                      'script_path')))

    if not args.model_runtime:
        config.set('script_input_arguments', 'model_runtime',
                   '{}/runtime_definitions/mpirun.xml'.format(
                       config.get('script_paths', 'script_path')))
        print(' WARNING: No runtime definition selected. Using the default '
              'of {}'.format(config.get('script_input_arguments',
                                        'model_runtime')))
    else:
        config.set('script_input_arguments', 'model_runtime',
                   args.model_runtime)

    if not config.has_section('conda'):
        config.add_section('conda')

    if not config.has_option('conda', 'link_load_compass'):
        config.set('conda', 'link_load_compass', 'False')

    if args.link_load_compass:
        config.set('conda', 'link_load_compass', 'True')

    # Build variables for history output
    old_dir = os.getcwd()
    os.chdir(config.get('script_paths', 'script_path'))
    git_version = subprocess.check_output(
        ['git', 'describe', '--tags', '--dirty']).decode('utf-8')
    git_version = git_version.strip('\n')
    calling_command = ""
    for arg in sys.argv:
        calling_command = "{}{} ".format(calling_command, arg)
    os.chdir(old_dir)

    # Iterate over all cases in the case_list.
    # There is only one if the (-o, -c, -r) options were used in place of (-n)
    for case_num in case_list:

        # If we're using a case_list, determine the core, configuration, and
        # resolution for the current test case.
        if use_case_list:
            core_configuration = subprocess.check_output(
                ['{}/list_testcases.py'.format(config.get('script_paths',
                                                          'script_path')),
                 '-n', '{:d}'.format(int(case_num))]).decode('utf-8')
            config_options = core_configuration.strip('\n').split(' ')
            config.set('script_input_arguments', 'core', config_options[1])
            config.set('script_input_arguments', 'configuration',
                       config_options[3])
            config.set('script_input_arguments', 'resolution',
                       config_options[5])
            config.set('script_input_arguments', 'test', config_options[7])

        # Setup each xml file in the configuration directory:
        test_path = '{}/{}/{}/{}'.format(
            config.get('script_input_arguments', 'core'),
            config.get('script_input_arguments', 'configuration'),
            config.get('script_input_arguments', 'resolution'),
            config.get('script_input_arguments', 'test'))
        work_dir = '{}/{}'.format(config.get('script_paths', 'work_dir'),
                                  test_path)

        # Set paths to core, configuration, resolution, and case for use in
        # functions
        config.set('script_paths', 'core_dir',
                   config.get('script_input_arguments', 'core'))
        config.set('script_paths', 'configuration_dir',
                   '{}/{}'.format(config.get('script_paths', 'core_dir'),
                                  config.get('script_input_arguments',
                                             'configuration')))
        config.set('script_paths', 'resolution_dir',
                   '{}/{}'.format(config.get('script_paths',
                                             'configuration_dir'),
                                  config.get('script_input_arguments',
                                             'resolution')))
        config.set('script_paths', 'test_dir', test_path)
        config.set('script_paths', 'config_path', test_path)

        # Only write history if we did something...
        write_history = False

        # Loop over all files in test_path that have the .xml extension.
        for file in os.listdir('{}'.format(test_path)):
            if fnmatch.fnmatch(file, '*.xml'):
                # Build full file name
                config_file = '{}/{}'.format(test_path, file)

                # Determine the type of the config file
                # Could be config, driver_script, template, etc.
                # The type is defined by the parent tag.
                config_type = get_config_file_type(config_file)

                # Process config files
                if config_type == 'config':
                    write_history = True
                    # Ensure the case directory exists
                    case_dir = make_case_dir(config_file, work_dir)
                    case_name = get_case_name(config_file)

                    # Set case_dir path for function calls
                    config.set('script_paths', 'case_dir',
                               '{}/{}'.format(config.get('script_paths',
                                                         'test_dir'),
                                              case_name))

                    case_path = '{}/{}'.format(work_dir, case_dir)

                    # Generate all namelists for this case
                    generate_namelist_files(config_file, case_path, config)

                    # Generate all streams files for this case
                    generate_streams_files(config_file, case_path, config)

                    # Ensure required files exist for this case
                    get_defined_files(config_file, '{}'.format(case_path),
                                      config)

                    # Process all links for this case
                    add_links(config_file, config)

                    copy_files(config_file, config)

                    # Generate run scripts for this case.
                    generate_run_scripts(config_file, '{}'.format(case_path),
                                         config)

                    print(" -- Set up case: {}/{}".format(work_dir, case_dir))
                # Process driver scripts
                elif config_type == 'driver_script':
                    write_history = True

                    # Generate driver scripts.
                    generate_driver_scripts(config_file, config)
                    print(" -- Set up driver script in {}".format(work_dir))

    # Write the history of this command to the command_history file, for
    # provenance.
    if write_history and not args.quiet:
        history_file_path = '{}/command_history'.format(
            config.get('script_paths', 'work_dir'))
        if os.path.exists(history_file_path):
            history_file = open(history_file_path, 'a')
            history_file.write('\n')
        else:
            history_file = open(history_file_path, 'w')

        history_file.write('**************************************************'
                           '*********************\n')
        history_file.write('git_version: {}\n'.format(git_version))
        history_file.write('command: {}\n'.format(calling_command))
        history_file.write('setup the following cases:\n')
        if use_case_list:
            for case_num in case_list:
                core_configuration = subprocess.check_output(
                    ['./list_testcases.py', '-n',
                     '{:d}'.format(int(case_num))]).decode('utf-8')
                config_options = core_configuration.strip('\n').split(' ')
                history_file.write('\n')
                history_file.write('    core: {}\n'.format(config_options[1]))
                history_file.write('    configuration: {}\n'.format(
                    config_options[3]))
                history_file.write('    resolution: {}\n'.format(
                    config_options[5]))
                history_file.write('    test: {}\n'.format(config_options[7]))
        else:
            history_file.write('core: {}\n'.format(
                config.get('script_input_arguments', 'core')))
            history_file.write('configuration: {}\n'.format(
                config.get('script_input_arguments', 'configuration')))
            history_file.write('resolution: {}\n'.format(
                config.get('script_input_arguments', 'resolution')))
            history_file.write('test: {}\n'.format(
                config.get('script_input_arguments', 'test')))

        history_file.write('**************************************************'
                           '*********************\n')
        history_file.close()

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
