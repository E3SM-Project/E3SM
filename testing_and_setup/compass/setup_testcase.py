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

import sys, os, glob, shutil, numpy, math
import fnmatch
import argparse
import xml.etree.ElementTree as ET
import subprocess
import ConfigParser

import netCDF4

try:
	from collections import defaultdict
except ImportError:
	from utils import defaultdict

# *** Namelist setup functions *** #{{{
def generate_namelist_files(config_file, case_path, configs):#{{{
	config_tree = ET.parse(config_file)
	config_root = config_tree.getroot()

	# Iterate over all namelists to be generated
	for namelists in config_root.iter('namelist'):
		# Determine the name of the namelist that will be generated
		try:
			namelist_file = '%s/%s'%(case_path, namelists.attrib['name'])
		except:
			print "ERROR: <namelist> tag is missing the 'name' attribute"
			print "Exiting..."
			sys.exit(1)


		try:
			namelist_mode = namelists.attrib['mode']
		except:
			print "ERROR: <namelist> tag is missing the 'mode' attribute."
			print "Exiting..."
			sys.exit(1)

		if not configs.has_option("namelists", namelist_mode):
			print "Error. Configuration file '%s' requires paths for streams and namelist files for '%s' mode."%(config_file, namelist_mode)
			print "Exiting..."
			sys.exit(1)

		template_namelist = configs.get("namelists", namelist_mode)

		# Ingest namelist template into a dictionary
		namelist_dict = defaultdict(lambda : defaultdict(list))
		ingest_namelist(template_namelist, namelist_dict)

		# Modify the dictionary to have the desired values
		configure_namelist(namelist_dict, namelists, configs)

		# Write the namelist using the template to determine the writing order.
		write_namelist(namelist_dict, namelist_file, template_namelist)
		del namelist_dict

	del config_root
	del config_tree
#}}}

def ingest_namelist(namelist_file, namelist_dict):#{{{
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
#}}}

def set_namelist_val(namelist_dict, option_name, option_val):#{{{
	# Set the value of the namelist option.
	for record, opts in namelist_dict.items():
		for opt, val in opts.items():
			if opt.strip() == option_name:
				val[0] = option_val
#}}}

def configure_namelist(namelist_dict, namelist_tag, configs):#{{{
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

#}}}

def apply_namelist_template(namelist_dict, template_tag, configs):#{{{
	# Determine the template information, like it's path and the filename
	template_info = get_template_info(template_tag, configs)

	# Build the full filename for the template
	template_file = '%s/%s'%(template_info['template_path'], template_info['template_file'])

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
#}}}

def write_namelist(namelist_dict, outfilename, infilename):#{{{
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
			out_namelist.write(line);
		elif line.find('=') >= 0:
			opt, val = line.strip().strip('\n').split('=')
			if record_name != "NONE!!!":
				out_namelist.write('    %s = %s\n'%(opt.strip(), namelist_dict[record_name][opt][0].strip()))

	if record_name != "NONE!!!":
		out_namelist.write('/\n')

	out_namelist.close()
#}}}
#}}}

# *** Streams setup functions *** #{{{

def generate_streams_files(config_file, case_path, configs):#{{{
	config_tree = ET.parse(config_file)
	config_root = config_tree.getroot()


	# Iterate over all sterams files to be generated
	for streams in config_root:
		if streams.tag == "streams":
			# Determine the path to the template streams file
			streams_filename = '%s/%s'%(case_path, streams.attrib['name'])

			try:
				streams_mode = streams.attrib['mode']
			except:
				print "ERROR: <streams> tag is missing the 'mode' attribute."
				print "Exiting..."
				sys.exit(1)

			if not configs.has_option("streams", streams_mode):
				print "Error. Configuration file '%s' requires paths for streams and namelist files for '%s' mode."%(config_file, streams_mode)
				print "Exiting..."
				sys.exit(1)

			template_streams = configs.get("streams", streams_mode)

			# Parse the template
			streams_tree = ET.parse(template_streams)
			streams_root = streams_tree.getroot()

			# Configure the new streams file, using the template as a starting place.
			configure_streams_file(streams_root, streams, configs)

			# Write out the streams file
			write_streams_file(streams_root, config_file, streams_filename, '%s'%(case_path))

			del streams_root
			del streams_tree
#}}}

def flush_streams(streams, remove_mutable, remove_immutable):#{{{
	if remove_mutable:
		# Remove all mutable streams from the template streams file
		for stream in streams.findall('stream'):
			streams.remove(stream)

	if remove_immutable:
		# Remove all immutable streams from the template streams file
		for stream in streams.findall('immutable_stream'):
			streams.remove(stream)
#}}}

def modify_stream_definition(streams_file, stream_conf):#{{{
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
					print "ERROR: Stream %s found multiple times in template. Exiting..."%(name.strip())
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
					except:
						print "   --- Tag: %s is missing a name attribute"%(child.tag)

#}}}

def configure_streams_file(streams_file, streams_tag, configs):#{{{
	keep_mode = streams_tag.attrib['keep']
	remove_immutable = False
	remove_mutable = False

	if ( keep_mode.strip() == 'immutable' ):
		remove_mutable = True
	if ( keep_mode.strip() == 'mutable' ):
		remove_immutable = True
	if ( keep_mode.strip() == 'none' ):
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
#}}}

def apply_stream_template(streams_file, template_tag, configs):#{{{
	# Determine template information, like path and filename
	template_info = get_template_info(template_tag, configs)

	# Build full path to template file
	template_file = '%s/%s'%(template_info['template_path'], template_info['template_file'])

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
#}}}

def write_streams_file(streams, config_file, filename, init_path):#{{{
	config_tree = ET.parse(config_file)
	config_root = config_tree.getroot()

	stream_file = open(filename, 'w')

	stream_file.write('<streams>\n')

	# Write out all immutable streams first
	for stream in streams.findall('immutable_stream'):
		stream_name = stream.attrib['name']

		stream_file.write('\n')
		stream_file.write('<immutable_stream name="%s"'%(stream_name))
		# Process all attributes on the stream
		for attr, val in stream.attrib.items():
			if ( not attr.strip() == 'name' ):
				stream_file.write('\n                  %s="%s"'%(attr, val))

		stream_file.write('/>\n')

	# Write out all immutable streams
	for stream in streams.findall('stream'):
		stream_name = stream.attrib['name']

		stream_file.write('\n')
		stream_file.write('<stream name="%s"'%(stream_name))

		# Process all attributes
		for attr, val in stream.attrib.items():
			if ( not attr.strip() == 'name' ):
				stream_file.write('\n        %s="%s"'%(attr, val))

		stream_file.write('>\n\n')

		# Write out all streams included in this stream
		for substream in stream.findall('stream'):
			substream_name = substream.attrib['name']
			if 'packages' in substream.attrib.keys():
				package_name = substream.attrib['packages']
				entry = '\t<stream name="%s"'%(substream_name) + ' packages="%s" '%(package_name) +'/>\n'
                        else:
				entry = '\t<stream name="%s"'%(substream_name) +'/>\n'
			stream_file.write(entry)

		# Write out all var_structs included in this stream
		for var_struct in stream.findall('var_struct'):
			var_struct_name = var_struct.attrib['name']
			if 'packages' in var_struct.attrib.keys():
				package_name = var_struct.attrib['packages']
				entry = '\t<var_struct name="%s"'%(var_struct_name) + ' packages="%s" '%(package_name) +'/>\n'
                        else:
				entry = '\t<var_struct name="%s"'%(var_struct_name) +'/>\n'
			stream_file.write(entry)

		# Write out all var_arrays included in this stream
		for var_array in stream.findall('var_array'):
			var_array_name = var_array.attrib['name']
			if 'packages' in var_array.attrib.keys():
				package_name = var_array.attrib['packages']
				entry = '\t<var_array name="%s"'%(var_array_name) + ' packages="%s" '%(package_name) +'/>\n'
                        else:
				entry = '\t<var_array name="%s"'%(var_array_name) +'/>\n'
			stream_file.write(entry)

		# Write out all vars included in this stream
		for var in stream.findall('var'):
			var_name = var.attrib['name']
			if 'packages' in var.attrib.keys():
				package_name = var.attrib['packages']
				entry = '\t<var name="%s"'%(var_name) + ' packages="%s" '%(package_name) +'/>\n'
                        else:
				entry = '\t<var name="%s"'%(var_name) +'/>\n'
			stream_file.write(entry)

		stream_file.write('</stream>\n')

	stream_file.write('\n')
	stream_file.write('</streams>\n')

	del config_tree
	del config_root
#}}}
#}}}

# *** Script Generation Functions *** #{{{
def generate_run_scripts(config_file, init_path, configs):#{{{
	config_tree = ET.parse(config_file)
	config_root = config_tree.getroot()
	dev_null = open('/dev/null', 'r+')

	for run_script in config_root:
		# Process run_script
		if run_script.tag == 'run_script':
			# Determine the name of the script, and create the file
			script_name = run_script.attrib['name']
			script_path = "%s/%s"%(init_path, script_name)
			script = open(script_path, "w")

			# Write the script header
			script.write("#!/usr/bin/env python\n")
			script.write("\n")
			script.write("### This script was generated from setup_testcases.py as part of a config file\n")
			script.write("\n")
			script.write("import sys, os, fnmatch, resource\n")
			script.write("import xml.etree.ElementTree as ET\n")
			script.write("import subprocess\n")
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
			subprocess.check_call(['chmod', 'a+x', '%s'%(script_path)], stdout=dev_null, stderr=dev_null, env=os.environ.copy())


	dev_null.close()
	del config_tree
	del config_root
#}}}

def generate_driver_scripts(config_file, configs):#{{{
	config_tree = ET.parse(config_file)
	config_root = config_tree.getroot()
	dev_null = open('/dev/null', 'r+')

	# init_path is where the driver script will live after it's generated.
	init_path = '%s/%s'%(config.get('script_paths', 'work_dir'), config.get('script_paths', 'config_path'))

	# Ensure we're in a <driver_script> tag
	if config_root.tag == 'driver_script':
		name = config_root.attrib['name']

		# Ensure work_dir exists before writing driver script there.
		if not os.path.exists(init_path):
			os.makedirs(init_path)

		# Create script file
		script = open('%s/%s'%(init_path, name), 'w')

		# Write script header
		script.write('#!/usr/bin/env python\n')
		script.write('"""\n')
		script.write('This script was generated as part of a driver_script file by the setup_testcases.py script.\n')
		script.write('"""\n')
		script.write('import sys, os, shutil, glob, subprocess\n')
		script.write("import xml.etree.ElementTree as ET\n")
		script.write('import argparse\n')
		script.write('\n')
		script.write('## This script was generated by setup_testcases.py as part of a driver_script file.\n')
		script.write("os.environ['PYTHONUNBUFFERED'] = '1'\n")
		script.write('parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)\n')

		case_dict = dict()
		for child in config_root:
			if child.tag == 'case':
				case_name = child.attrib['name']
				case_dict[case_name] = '1'

		for case_name in case_dict.keys():
			script.write('parser.add_argument("--no_%s", dest="no_%s", help="If set, %s case will not be run during execution of this script.", action="store_true")\n'%(case_name, case_name, case_name))
			script.write('parser.add_argument("--finalize_%s", dest="finalize_%s", help="If set, %s case will have symlinks replaced with the files they point to, this occurs after any case runs that have been requested..", action="store_true")\n'%(case_name, case_name, case_name))

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
				script.write('if not args.no_%s:\n'%(case))
				script.write('\tos.chdir(base_path)\n')
				script.write('\tos.chdir(' + "'%s')\n"%(case))
				# Process children of <case> tag
				for grandchild in child:
					# Process <step> tags
					if grandchild.tag == 'step':
						process_script_step(grandchild, configs, '\t', script)
					# Process <define_env_var> tags
					elif grandchild.tag == 'define_env_var':
						process_env_define_step(grandchild, configs, '\t', script)
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

		# Write script footer, that ensures a 1 is returned if the script encountered an error. This happens before finalizing a case directory.
		script.write('if error:\n')
		script.write('\tsys.exit(1)\n')

		for case_name in case_dict.keys():
			script.write('if args.finalize_%s:\n'%(case_name))
			script.write('    old_dir = os.getcwd()\n')
			script.write('    os.chdir("%s")\n'%(case_name))
			script.write('    file_list = glob.glob("*")\n')
			script.write('    for file in file_list:\n')
			script.write('       if os.path.islink(file):\n')
			script.write('          link_path = os.readlink(file)\n')
			script.write('          os.unlink(file)\n')
			script.write('          shutil.copyfile(link_path, file)\n')
			#script.write('          subprocess.check_call(["cp", "--remove-destination", link_path, file], stdout=dev_null, stderr=dev_null)\n')
			script.write('    os.chdir(old_dir)\n')
			script.write('\n')


		script.write('sys.exit(0)\n')
		script.close()
		del case_dict

		# Make script executable
		subprocess.check_call(['chmod', 'a+x', '%s/%s'%(init_path, name)], stdout=dev_null, stderr=dev_null, env=os.environ.copy())
#}}}

def process_env_define_step(var_tag, configs, indentation, script_file):#{{{
	try:
		var_name = var_tag.attrib['name']
	except:
		print "ERROR: <define_env_var> tag is missing 'name' attribte"
		print 'Exiting...'
		sys.exit(1)
	
	try:
		var_val = var_tag.attrib['value']
	except:
		print "ERROR: <define_env_var> tag is missing 'value' attribute"
		print 'Exiting...'
		sys.exit(1)
	
	# Write line to define the environment variable
	script_file.write("%sos.environ['%s'] = '%s'\n"%(indentation, var_name, var_val))
#}}}

def process_script_step(step, configs, indentation, script_file):#{{{
	# Determine step attributes.
	if 'executable_name' in step.attrib.keys() and 'executable' in step.attrib.keys():
		print "ERROR: <step> tag has both an 'executable' and 'executable_name' attribute. Only one is allowed per step."
		print "Exiting..."
		sys.exit(1)

	try:
		quiet_val = step.attrib['quiet']
		if quiet_val == "true":
			quiet = True
		else:
			quiet = False
	except:
		quiet = False

	try:
		step_pre_message = step.attrib['pre_message']
		write_pre_message = True
	except:
		write_pre_message = False

	try:
		step_post_message = step.attrib['post_message']
		write_post_message = True
	except:
		write_post_message = False

	try:
		executable_name = step.attrib['executable_name']
		executable = configs.get('executables', executable_name)
	except:
		executable = step.attrib['executable']

	# Write step header
	script_file.write("\n")
	script_file.write("# Run command is:\n")

	# Build comment and command bases
	comment = "# %s "%(executable)
	command = "subprocess.check_call(['%s'"%(executable)

	# Process step arguments
	for argument in step:
		if argument.tag == 'argument':
			flag = argument.attrib['flag']
			val = argument.text

			if not flag.strip() == "":
				comment = "%s %s"%(comment, flag)
				command = "%s, '%s'"%(command, flag)

			comment = "%s %s"%(comment, val)
			command = "%s, '%s'"%(command, val)

	# If a pre_message attribute was supplied, write it before adding the command.
	if write_pre_message:
		script_file.write('%sprint "%s"\n'%(indentation, step_pre_message))

	# Setup command redirection
	if quiet:
		command = "%s], stdout=dev_null, stderr=dev_null"%(command)
	else:
		command = "%s]"%(command)

	# Write the comment, and the command. Also, ensure the command has the same
	# environment as the calling script.
	script_file.write("%s%s\n"%(indentation, comment))
	script_file.write("%s%s, env=os.environ.copy())\n"%(indentation, command))

	# If a post_message attribute was supplied, write it after the command.
	if write_post_message:
		script_file.write('%sprint "%s"\n'%(indentation, step_post_message))

	script_file.write("\n");
#}}}

def process_validation_step(validation_tag, configs, script):#{{{
	for child in validation_tag:
		if child.tag == 'compare_fields':
			process_compare_fields_step(child, configs, script)
		if child.tag == 'compare_timers':
			process_compare_timers_step(child, configs, script)
#}}}

# *** Field Comparison Functions *** ##{{{
def process_compare_fields_step(compare_tag, configs, script):#{{{
	missing_file1 = False
	missing_file2 = False
	# Determine comparision attributes
	try:
		file1 = compare_tag.attrib['file1']
	except:
		missing_file1 = True

	try:
		file2 = compare_tag.attrib['file2']
	except:
		missing_file2 = True

	if missing_file1 and missing_file2:
		print "ERROR: <compare_fields> tag is missing both 'file1' and 'file2' tags. At least one is required."
		print "Exiting..."
		sys.exit(1)
	
	baseline_root = configs.get('script_paths', 'baseline_dir')

	if not baseline_root == 'NONE':
		baseline_root = '%s/%s'%(baseline_root, configs.get('script_paths', 'test_dir'))
	
	for child in compare_tag:
		# Process field comparisions
		if child.tag == 'field':
			if not (missing_file1 or missing_file2):
				process_field_definition(child, configs, script, file1, file2, False)

			if not missing_file1 and not baseline_root == 'NONE':
				process_field_definition(child, configs, script, file1, '%s/%s'%(baseline_root, file1), True)

			if not missing_file2 and not baseline_root == 'NONE':
				process_field_definition(child, configs, script, file2, '%s/%s'%(baseline_root, file2), True)
		# Process field comparision template
		elif child.tag == 'template':
			apply_compare_fields_template(child, compare_tag, configs, script)
#}}}

def apply_compare_fields_template(template_tag, compare_tag, configs, script):#{{{
	missing_file1 = False
	missing_file2 = False
	# Determine comparision attributes
	try:
		file1 = compare_tag.attrib['file1']
	except:
		missing_file1 = True

	try:
		file2 = compare_tag.attrib['file2']
	except:
		missing_file2 = True

	if missing_file1 and missing_file2:
		print "ERROR: <compare_fields> tag is missing both 'file1' and 'file2' tags. At least one is required."
		print "Exiting..."
		sys.exit(1)

	# Build the path to the baselines
	baseline_root = configs.get('script_paths', 'baseline_dir')
	if not baseline_root == 'NONE':
		baseline_root = '%s/%s'%(baseline_root, configs.get('script_paths', 'test_dir'))

	# Determine template information, like path and filename
	template_info = get_template_info(template_tag, configs)

	template_file = '%s/%s'%(template_info['template_path'], template_info['template_file'])

	# Parse the template
	template_tree = ET.parse(template_file)
	template_root = template_tree.getroot()

	# Find a child tag that is validation->compare_fields->field, and add each field
	for validation in template_root:
		if validation.tag == 'validation':
			for compare_fields in validation:
				if compare_fields.tag == 'compare_fields':
					for field in compare_fields:
						if field.tag == 'field':
							if not (missing_file1 or missing_file2):
								process_field_definition(field, configs, script, file1, file2, False)

							if not missing_file1 and not baseline_root == 'NONE':
								process_field_definition(field, configs, script, file1, '%s/%s'%(baseline_root, file1), True)

							if not missing_file2 and not baseline_root == 'NONE':
								process_field_definition(field, configs, script, file2, '%s/%s'%(baseline_root, file2), True)
						elif field.tag == 'template':
							apply_compare_fields_template(field, compare_tag, configs, script)

	del template_root
	del template_tree
	del template_info
#}}}

def process_field_definition(field_tag, configs, script, file1, file2, baseline_comp):#{{{
	# Build the path to the comparision script.
	compare_executable = '%s/compare_fields.py'%(configs.get('script_paths', 'utility_scripts'))

	# Build the base command to compare the fields
	base_command = "subprocess.check_call(['%s', '-q', '-1', '%s', '-2', '%s'"%(compare_executable, file1, file2)

	field_name = field_tag.attrib['name']

	command = "%s, '-v', '%s'"%(base_command, field_name)
	
	# Determine norm thresholds
	if not baseline_comp:
		if 'l1_norm' in field_tag.attrib.keys():
			command = "%s, '--l1', '%s'"%(command, field_tag.attrib['l1_norm'])
		if 'l2_norm' in field_tag.attrib.keys():
			command = "%s, '--l2', '%s'"%(command, field_tag.attrib['l2_norm'])
		if 'linf_norm' in field_tag.attrib.keys():
			command = "%s, '--linf', '%s'"%(command, field_tag.attrib['linf_norm'])
	else:
		command = "%s, '--l1', '0.0', '--l2', '0.0', '--linf', '0.0'"%(command)

	# Ensure the comparision script has the same environment as the calling script.
	command = '%s], env=os.environ.copy())'%(command)

	# Write the pass/fail logic.
	script.write('try:\n')
	script.write('\t%s\n'%(command))
	script.write("\tprint ' ** PASS Comparison of %s between %s and %s'\n"%(field_name, file1, file2))
	script.write('except:\n')
	script.write("\tprint ' ** FAIL Comparison of %s between %s and %s'\n"%(field_name, file1, file2))
	script.write('\terror = True\n')
#}}}
#}}}

# *** Timer Comparison Functions *** ##{{{
def process_compare_timers_step(compare_tag, configs, script):#{{{
	baseline_root = configs.get('script_paths', 'baseline_dir')
	baseline_root = '%s/%s'%(baseline_root, configs.get('script_paths', 'test_dir'))

	missing_rundir1 = True
	missing_rundir2 = True

	try:
		rundir1 = compare_tag.attrib['rundir1']
		missing_rundir1 = False
	except:
		missing_rundir1 = True

	try:
		rundir2 = compare_tag.attrib['rundir2']
		missing_rundir2 = False
	except:
		missing_rundir2 = True

	for child in compare_tag:
		if child.tag == 'timer':
			try:
				child_name = child.attrib['name']
			except:
				print "ERROR: <timer> tag is missing the 'name' attribute."
				print "Exiting..."
				sys.exit(1)

			if not (missing_rundir1 or missing_rundir2):
				process_timer_definition(child, configs, script, rundir1, rundir2)

			if not missing_rundir1:
				process_timer_definition(child, configs, script, '%s/%s'%(baseline_root, rundir1), rundir1)

			if not missing_rundir2:
				process_timer_definition(child, configs, script, '%s/%s'%(baseline_root, rundir2), rundir2)
		elif child.tag == 'template':
			apply_compare_timers_template(child, compare_tag, configs, script)

#}}}

def apply_compare_timers_template(template_tag, compare_tag, configs, script):#{{{
	# Build the path to the baselines
	baseline_root = configs.get('script_paths', 'baseline_dir')
	baseline_root = '%s/%s'%(baseline_root, configs.get('script_paths', 'test_dir'))

	missing_rundir1 = True
	missing_rundir2 = True

	try:
		rundir1 = compare_tag.attrib['rundir1']
		missing_rundir1 = False
	except:
		missing_rundir1 = True

	try:
		rundir2 = compare_tag.attrib['rundir2']
		missing_rundir2 = False
	except:
		missing_rundir2 = True

	# Get the template information and build the template file
	template_info = get_template_info(template_tag, configs)
	template_file = '%s/%s'%(template_info['template_path'], template_info['template_file'])
	
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
								process_timer_definition(timer, configs, script, rundir1, rundir2)

							if not missing_rundir1:
								process_timer_definition(timer, configs, script, '%s/%s'%(baseline_root, rundir1), rundir1)

							if not missing_rundir2:
								process_timer_definition(timer, configs, script, '%s/%s'%(baseline_root, rundir2), rundir2)
						elif timer.tag == 'template':
							apply_compare_timers_template(timer, compare_tag, configs, script)
	del template_root
	del template_tree
	del template_info

#}}}

def process_timer_definition(timer_tag, configs, script, basedir, compdir):#{{{
	compare_script = '%s/compare_timers.py'%(configs.get('script_paths', 'utility_scripts'))

	try:
		timer_name = timer_tag.attrib['name']
	except:
		print "ERROR: <timer> tag is missing the 'name' attribute."
		print "Exiting..."
		sys.exit(1)

	command = 'subprocess.check_call(["%s", "-b", "%s", "-c", "%s", "-t", "%s"], env=os.environ.copy())'%(compare_script, basedir, compdir, timer_name)

	script.write('\n')
	script.write('if os.path.exists("%s") and os.path.exists("%s"):\n'%(basedir, compdir))
	script.write('\ttry:\n')
	script.write('\t\t%s\n'%(command))
	script.write("\t\tprint ' ** PASS Comparision of timer %s between %s and %s'\n"%(timer_name, basedir, compdir))
	script.write('\texcept:\n')
	script.write("\t\tprint ' ** FAIL Comparision of timer %s between %s and %s'\n"%(timer_name, basedir, compdir))
	script.write("\t\terror = True\n")
#}}}
#}}}

def process_model_run_step(model_run_tag, configs, script):#{{{
	run_definition_file = configs.get('script_input_arguments', 'model_runtime')
	run_config_tree = ET.parse(run_definition_file)
	run_config_root = run_config_tree.getroot()

	dev_null = open('/dev/null', 'r+')

	try:
		executable_name = model_run_tag.attrib['executable']
	except:
		executable_name = 'model'

	script.write('print "\\n"\n')
	script.write('print "     *****************************"\n')
	script.write('print "     ** Starting model run step **"\n')
	script.write('print "     *****************************"\n')
	script.write('print "\\n"\n')

	# Process each part of the run script
	for child in run_config_root:
		# Process each <step> tag
		if child.tag == 'step':
			# Setup child step, and it's attributes to be correct for the process_script_step function
			for grandchild in child:
				if grandchild.tag == 'argument':
					arg_text = grandchild.text

					if arg_text == 'model':
						executable_full_path = config.get('executables', executable_name)
						executable_parts = executable_full_path.split('/')
						executable_link = executable_parts[ len(executable_parts) - 1]
						link_path = '%s/%s/%s'%(config.get('script_paths', 'work_dir'), config.get('script_paths', 'case_dir'), executable_link)
						subprocess.check_call(['ln', '-sf', config.get('executables', executable_name), link_path], stdout=dev_null, stderr=dev_null)
						grandchild.text = executable_link
					elif arg_text.find('attr_') >= 0:
						attr_array = arg_text.split('_')
						try:
							grandchild.text = model_run_tag.attrib[attr_array[1]]
						except:
							print " <step> tag defined within a <model_run> tag requires attribute '%s', but it is not defined."%(attr_array[1])
							print " Exiting..."
							sys.exit(1)

			# Process the resulting element, instead of the original step.
			process_script_step(child, configs, '', script)
		# Process each <define_env_var> tag
		elif child.tag == 'define_env_var':
			if child.attrib['value'].find('attr_') >= 0:
				attr_array = child.attrib['value'].split('_')
				try:
					child.attrib['value'] = model_run_tag.attrib[attr_array[1]]
				except:
					print " <define_env_var> tag defined within a <model_run> tag requires attribute '%s', but it is not defined."%(attr_array[1])
					print " Exiting..."
					sys.exit(1)
				
			process_env_define_step(child, configs, '', script)

	script.write('print "\\n"\n')
	script.write('print "     *****************************"\n')
	script.write('print "     ** Finished model run step **"\n')
	script.write('print "     *****************************"\n')
	script.write('print "\\n"\n')
	dev_null.close()
#}}}
#}}}

# *** General Utility Functions *** #{{{
def add_links(config_file, configs):#{{{
	config_tree = ET.parse(config_file)
	config_root = config_tree.getroot()

	case = config_root.attrib['case']

	dev_null = open('/dev/null', 'r+')

	# Determine the path for the case directory
	test_path = '%s/%s'%(configs.get('script_paths', 'test_dir'), case)
	base_path = '%s/%s'%(configs.get('script_paths', 'work_dir'), test_path)

	# Process all children tags
	for child in config_root:
		# Process an <add_link> tag
		if child.tag == 'add_link':
			try:
				source = child.attrib['source']
			except:
				print " add_link tag missing a 'source' attribute."
				print " Exiting..."
				sys.exit(1)

			try:
				source_path_name = child.attrib['source_path']

				keyword_path = False
				if source_path_name.find('work_') >= 0:
					keyword_path = True
				elif source_path_name.find('script_') >= 0:
					keyword_path = True

				if not keyword_path:
					try:
						source_path = configs.get('paths', source_path_name)
					except:
						source_path = 'NONE'

					if source_path == 'NONE':
						try:
							source_path = configs.get('script_paths', source_path_name)
						except:
							source_path = 'NONE'

					if source_path == 'NONE':
						print "ERROR: source_path on <add_link> tag is '%s' which is not defined"%(source_path_name)
						print "Exiting..."
						sys.exit(1)

				else:
					source_arr = source_path_name.split('_')
					base_name = source_arr[0]
					subname = '%s_%s'%(source_arr[1], source_arr[2])

					if base_name == 'work':
						file_base_path = 'work_dir'
					elif base_name == 'script':
						file_base_path = 'script_path'

					if subname in {'core_dir', 'configuration_dir', 'resolution_dir', 'test_dir', 'case_dir'}:
						source_path = '%s/%s'%(configs.get('script_paths', file_base_path), configs.get('script_paths', subname))

				source_file = '%s/%s'%(source_path, source)
			except:
				source_file = '%s'%(source)

			dest = child.attrib['dest']
			old_cwd = os.getcwd()
			os.chdir(base_path)

			subprocess.check_call(['ln', '-sfn', '%s'%(source_file), '%s'%(dest)], stdout=dev_null, stderr=dev_null, env=os.environ.copy())
			os.chdir(old_cwd)
			del source
			del dest
		# Process an <add_executable> tag
		elif child.tag == 'add_executable':
			source_attr = child.attrib['source']
			dest = child.attrib['dest']
			if not configs.has_option("executables", source_attr):
				print 'ERROR: Configuration %s requires a definition of %s.'%(config_file, source_attr)
				sys.exit(1)
			else:
				source = configs.get("executables", source_attr)

			subprocess.check_call(['ln', '-sf', '%s'%(source), '%s/%s'%(base_path, dest)], stdout=dev_null, stderr=dev_null, env=os.environ.copy())
			del source_attr
			del source
			del dest

	del config_tree
	del config_root
	dev_null.close()
#}}}

def make_case_dir(config_file, base_path):#{{{
	config_tree = ET.parse(config_file)
	config_root = config_tree.getroot()

	case_name = config_root.attrib['case']

	# Build the case directory, if it doesn't already exist
	if not os.path.exists('%s/%s'%(base_path, case_name)):
		os.makedirs('%s/%s'%(base_path, case_name))

	del config_root
	del config_tree

	return case_name
#}}}

def get_defined_files(config_file, init_path, configs):#{{{
	config_tree = ET.parse(config_file)
	config_root = config_tree.getroot()
	dev_null = open('/dev/null', 'w')

	for get_file in config_root:
		# Process <get_file> tag
		if get_file.tag == 'get_file':
			# Determine dest_path
			try:
				dest_path_name = get_file.attrib['dest_path']
			except:
				print " get_file tag is missing the 'dest_path' attribute."
				print " Exiting..."
				sys.exit(1)

			# Determine file_name
			try:
				file_name = get_file.attrib['file_name']
			except:
				print " get_file tag is missing a 'file_name' attribute."
				print " Exiting..."
				sys.exit(1)

			# Build out the dest path
			keyword_path = False
			if dest_path_name.find('work_') >= 0:
				keyword_path = True
			elif dest_path_name.find('script_') >= 0:
				keyword_path = True
			else:
				try:
					dest_path = '%s'%(configs.get('paths', dest_path_name))
				except:
					print " Path '%s' is not defined in the config file, but is required to get a file."%(dest_path_name)
					print " Exiting..."
					sys.exit(1)


			if keyword_path:
				dest_arr = dest_path_name.split('_')
				base_name = dest_arr[0]
				subname = '%s_%s'%(dest_arr[1], dest_arr[2])

				if base_name == 'work':
					base_path = 'work_dir'
				elif base_name == 'script':
					base_path = 'script_path'

				if subname in {'core_dir', 'configuration_dir', 'resolution_dir', 'test_dir', 'case_dir'}:
					dest_path = '%s/%s'%(configs.get('script_paths', base_path), configs.get('script_paths', subname))
				else:
					print " Path '%s' is not defined."%(dest_path_name)
					print " Exiting..."
					sys.exit(1)

			# if the dest_path doesn't exist, create it
			if not os.path.exists(dest_path):
				os.makedirs(dest_path)

			# If the file doesn't exist in dest_path, process it's mirrors
			if not os.path.exists('%s/%s'%(dest_path, file_name)):
				file_found = False
				if not file_found:
					for mirror in get_file:
						# Process each mirror
						if mirror.tag == 'mirror':
							# Determine the protocol for the mirror
							try:
								protocol = mirror.attrib['protocol']
							except:
								print "Mirror is missing the 'protocol' attribute."
								print "Exiting..."
								sys.exit(1)

							# Process a wget mirror
							if protocol == 'wget':
								if config.get('script_input_arguments', 'no_download') == 'no':
									try:
										path = '%s/%s'%( mirror.attrib['url'], file_name)
									except:
										print " Mirror with protocol 'wget' is missing a 'url' attribute"
										print " Exiting..."
										sys.exit(1)

									try:
										subprocess.check_call(['wget', '-q', '%s'%(path)], stdout=dev_null, stderr=dev_null, env=os.environ.copy())
									except:
										print " Wget failed...."

									try:
										subprocess.check_call(['mv', '%s'%(file_name), '%s/%s'%(dest_path, file_name)], stdout=dev_null, stderr=dev_null, env=os.environ.copy())
									except:
										print "  -- Web mirror attempt failed. Trying other mirrors..."


						# If the file exists in dest_path, determine if it should be validated
						if os.path.exists('%s/%s'%(dest_path, file_name)):
							try:
								expected_hash = get_file.attrib['hash']
								validate_file = True
							except:
								validate_file = False

							# Validate the file, if requested (i.e. file had a hash attribute)
							if validate_file:
								if os.path.exists('%s/%s'%(dest_path, file_name)):
									nc = netCDF4.Dataset('%s/%s'%(dest_path, file_name), 'r')
									try:
										file_hash = nc.file_id
									except:
										nc.close()
										print " Downloaded file '%s' does not have a 'file_id' attribute."%(file_name)
										print " Deleting file and exiting..."
										subprocess.check_call(['rm', '-f', '%s/%s'%(dest_path, file_name)], stdout=dev_null, stderr=dev_null, env=os.environ.copy())
										sys.exit(1)

									nc.close()

									if not file_hash.strip() == expected_hash.strip():
										print "*** ERROR: Base mesh has hash of '%s' which does not match expected hash of '%s'."%(file_hash, expected_hash)
										print " Deleting file..."
										subprocess.check_call(['rm', '-f', '%s/%s'%(dest_path, file_name)], stdout=dev_null, stderr=dev_null, env=os.environ.copy())

					# IF validation valied, exit.
					if not os.path.exists('%s/%s'%(dest_path, file_name)):
						print " Failed to acquire required file '%s'."%(file_name)
						print " Exiting..."
						sys.exit(1)

	del config_tree
	del config_root
	dev_null.close()
#}}}

def get_template_info(template, configs):#{{{
	# Determine template attributes
	try:
		template_file = template.attrib['file']
	except:
		print "ERROR: <template> tag is missing 'file' attribute."
		print 'Exiting...'
		sys.exit(1)
	
	try:
		template_path_base = template.attrib['path_base']
	except:
		print "ERROR: <template> tag is missing 'path_base' attribute."
		print 'Exiting...'
		sys.exit(1)

	try:
		template_path_modifier = '/%s'%(template.attrib['path'])
	except:
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

		sub_path = '%s_%s'%(template_arr[1], template_arr[2])

		if sub_path in {'core_dir', 'configuration_dir', 'resolution_dir', 'test_dir'}:
			template_path = '%s/%s%s'%(configs.get('script_paths', base_path), configs.get('script_paths', sub_path), template_path_modifier)
		else:
			print " Template path_base '%s' does not exist"%(template_path_base)
			print " Exiting..."
			sys.exit(1)

	info = {}
	info['template_file'] = template_file
	info['template_path'] = template_path

	return info
#}}}

def get_config_file_type(config_file):#{{{
	config_tree = ET.parse(config_file)
	config_root = config_tree.getroot()

	# Determine file type
	file_type = config_root.tag

	del config_root
	del config_tree

	return file_type
#}}}

def get_case_name(config_file):#{{{
	config_tree = ET.parse(config_file)
	config_root = config_tree.getroot()

	# Determine file type
	if config_root.tag == 'config':
		name = config_root.attrib['case']

	del config_root
	del config_tree

	return name
#}}}
#}}}

if __name__ == "__main__":
	# Define and process input arguments
	parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument("-o", "--core", dest="core", help="Core that contains configurations", metavar="CORE")
	parser.add_argument("-c", "--configuration", dest="configuration", help="Configuration to setup", metavar="CONFIG")
	parser.add_argument("-r", "--resolution", dest="resolution", help="Resolution of configuration to setup", metavar="RES")
	parser.add_argument("-t", "--test", dest="test", help="Test name within a resolution to setup", metavar="TEST")
	parser.add_argument("-n", "--case_number", dest="case_num", help="Case number to setup, as listed from list_testcases.py. Can be a comma delimited list of case numbers.", metavar="NUM")
	parser.add_argument("-f", "--config_file", dest="config_file", help="Configuration file for test case setup", metavar="FILE")
	parser.add_argument("-m", "--model_runtime", dest="model_runtime", help="Definition of how to build model run commands on this machine", metavar="FILE")
	parser.add_argument("-b", "--baseline_dir", dest="baseline_dir", help="Location of baseslines that can be compared to", metavar="PATH")
	parser.add_argument("-q", "--quiet", dest="quiet", help="If set, script will not write a command_history file", action="store_true")
	parser.add_argument("--no_download", dest="no_download", help="If set, script will not auto-download base_mesh files", action="store_true")
	parser.add_argument("--work_dir", dest="work_dir", help="If set, script will create case directories in work_dir rather than the current directory.", metavar="PATH")

	args = parser.parse_args()

	if not args.config_file:
		print "WARNING: No configuration file specified. Using the default of 'local.config'"
		args.config_file = 'local.config'

	if not os.path.exists(args.config_file):
		parser.error(" Configuration file '%s' does not exist. Please create and setup before running again."%(args.config_file))

	if not args.case_num and not ( args.core and args.configuration and args.resolution and args.test):
		print 'Must be run with either the --case_number argument, or the core, configuration, resolution, and test arguments.'
		parser.error(' Invalid configuration. Exiting...')

	if args.case_num and args.core and args.configuration and args.resoltuion and args.test:
		print 'Can only be configured with either --case_number (-n) or --core (-o), --configuration (-c), --resolution (-r), and --test (-t).'
		parser.error(' Invalid configuration. Too many options used. Exiting...')

	if args.case_num:
		use_case_list = True
		case_list = args.case_num.split(',')
	else:
		use_case_list = False
		case_list = list()
		case_list.append(0)

	config = ConfigParser.SafeConfigParser()
	config.read(args.config_file)

	if not args.no_download:
		args.no_download = False

	if not args.work_dir:
		args.work_dir = os.getcwd()

	# Add configuation information to the config object.
	# This allows passing config around with all of the config options needed to
	# build paths, and determine options.
	config.add_section('script_input_arguments')
	config.add_section('script_paths')

	if not use_case_list:
		config.set('script_input_arguments', 'core', args.core)
		config.set('script_input_arguments', 'configuration', args.configuration)
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

	config.set('script_paths', 'script_path', os.path.dirname(os.path.realpath(__file__)))
	config.set('script_paths', 'work_dir', os.path.abspath(args.work_dir) )
	config.set('script_paths', 'utility_scripts', '%s/utility_scripts'%(config.get('script_paths', 'script_path')))

	if not args.model_runtime:
		config.set('script_input_arguments', 'model_runtime', '%s/runtime_definitions/mpirun.xml'%(config.get('script_paths', 'script_path')))
		print ' WARNING: No runtime definition selected. Using the default of %s'%(config.get('script_input_arguments', 'model_runtime'))
	else:
		config.set('script_input_arguments', 'model_runtime', args.model_runtime)

	# Build variables for history output
	old_dir = os.getcwd()
	os.chdir( config.get('script_paths', 'script_path' ))
	git_version = subprocess.check_output(['git', 'describe', '--tags', '--dirty'])
	git_version = git_version.strip('\n')
	calling_command = ""
	for arg in sys.argv:
		calling_command = "%s%s "%(calling_command, arg)
	os.chdir(old_dir)

	# Iterate over all cases in the case_list.
	# There is only one if the (-o, -c, -r) options were used in place of (-n)
	for case_num in case_list:

		# If we're using a case_list, determine the core, configuration, and
		# resolution for the current test case.
		if use_case_list:
			core_configuration = subprocess.check_output(['%s/list_testcases.py'%(config.get('script_paths', 'script_path')), '-n', '%d'%(int(case_num))])
			config_options = core_configuration.strip('\n').split(' ')
			config.set('script_input_arguments', 'core', config_options[1])
			config.set('script_input_arguments', 'configuration', config_options[3])
			config.set('script_input_arguments', 'resolution', config_options[5])
			config.set('script_input_arguments', 'test', config_options[7])

		# Setup each xml file in the configuration directory:
		test_path = '%s/%s/%s/%s'%(config.get('script_input_arguments', 'core'), config.get('script_input_arguments', 'configuration'), config.get('script_input_arguments', 'resolution'), config.get('script_input_arguments', 'test'))
		work_dir = '%s/%s'%(config.get('script_paths', 'work_dir'), test_path)

		# Set paths to core, configuration, resolution, and case for use in functions
		config.set('script_paths', 'core_dir', config.get('script_input_arguments', 'core'))
		config.set('script_paths', 'configuration_dir', '%s/%s'%(config.get('script_paths', 'core_dir'), config.get('script_input_arguments', 'configuration')))
		config.set('script_paths', 'resolution_dir', '%s/%s'%(config.get('script_paths', 'configuration_dir'), config.get('script_input_arguments', 'resolution')))
		config.set('script_paths', 'test_dir', test_path)
		config.set('script_paths', 'config_path', test_path)

		# Only write history if we did something...
		write_history = False

		# Loop over all files in test_path that have the .xml extension.
		for file in os.listdir('%s'%(test_path)):
			if fnmatch.fnmatch(file, '*.xml'):
				# Build full file name
				config_file = '%s/%s'%(test_path, file)

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
					config.set('script_paths', 'case_dir', '%s/%s'%(config.get('script_paths', 'test_dir'), case_name))

					case_path = '%s/%s'%(work_dir, case_dir)

					# Generate all namelists for this case
					generate_namelist_files(config_file, case_path, config)

					# Generate all streams files for this case
					generate_streams_files(config_file, case_path, config)

					# Ensure required files exist for this case
					get_defined_files(config_file, '%s'%(case_path), config)

					# Process all links for this case
					add_links(config_file, config)

					# Generate run scripts for this case.
					generate_run_scripts(config_file, '%s'%(case_path), config)

					print " -- Set up case: %s/%s"%(work_dir, case_dir)
				# Process driver scripts
				elif config_type == 'driver_script':
					write_history = True

					# Generate driver scripts.
					generate_driver_scripts(config_file, config)
					print " -- Set up driver script in %s"%(work_dir)

	# Write the history of this command to the command_history file, for
	# provenance.
	if write_history and not args.quiet:
		history_file_path = '%s/command_history'%(config.get('script_paths', 'work_dir'))
		if os.path.exists(history_file_path):
			history_file = open(history_file_path, 'a')
			history_file.write('\n')
		else:
			history_file = open(history_file_path, 'w')

		history_file.write('***********************************************************************\n')
		history_file.write('git_version: %s\n'%(git_version))
		history_file.write('command: %s\n'%(calling_command))
		history_file.write('setup the following cases:\n')
		if use_case_list:
			for case_num in case_list:
				core_configuration = subprocess.check_output(['./list_testcases.py', '-n', '%d'%(int(case_num))])
				config_options = core_configuration.strip('\n').split(' ')
				history_file.write('\n')
				history_file.write('\tcore: %s\n'%(config_options[1]))
				history_file.write('\tconfiguration: %s\n'%(config_options[3]))
				history_file.write('\tresolution: %s\n'%(config_options[5]))
				history_file.write('\ttest: %s\n'%(config_options[7]))
		else:
			history_file.write('core: %s\n'%(config.get('script_input_arguments', 'core')))
			history_file.write('configuration: %s\n'%(config.get('script_input_arguments', 'configuration')))
			history_file.write('resolution: %s\n'%(config.get('script_input_arguments', 'resolution')))
			history_file.write('test: %s\n'%(config.get('script_input_arguments', 'test')))

		history_file.write('***********************************************************************\n')
		history_file.close()

