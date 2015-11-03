#!/usr/bin/env python

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
def generate_namelist_files(config_file, case_path, template_namelist, templates_path):#{{{
	config_tree = ET.parse(config_file)
	config_root = config_tree.getroot()

	for namelists in config_root.iter('namelist'):
		# -- Generate namelist
		namelist_file = '%s/%s'%(case_path, namelists.attrib['name'])

		namelist_dict = defaultdict(lambda : defaultdict(list))
		ingest_namelist(template_namelist, namelist_dict)
		configure_namelist(namelist_dict, templates_path, namelists)
		write_namelist(namelist_dict, namelist_file, template_namelist)
		del namelist_dict

	del config_root
	del config_tree
#}}}

def ingest_namelist(namelist_file, namelist_dict):#{{{
	namelistfile = open(namelist_file, 'r+')
	lines = namelistfile.readlines()

	record_name = 'NONE!!!'

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
	for record, opts in namelist_dict.items():
		for opt, val in opts.items():
			if opt.strip() == option_name:
				val[0] = option_val
#}}}

def apply_namelist_template(namelist_dict, template_name, template_path):#{{{
	template_file = '%s/%s.xml'%(template_path, template_name)
	template_tree = ET.parse(template_file)
	template_root = template_tree.getroot()

	for namelist in template_root.findall('namelist'):
		for option in namelist.findall('option'):
			option_name = option.attrib['name']
			option_val = option.text
			set_namelist_val(namelist_dict, option_name, option_val)
#}}}

def configure_namelist(namelist_dict, template_path, namelist_tag):#{{{
	for child in namelist_tag.iter('*'):
		if child.tag == 'option':
			option_name = child.attrib['name']
			option_val = child.text
			set_namelist_val(namelist_dict, option_name, option_val)
		elif child.tag == 'template':
			apply_namelist_template(namelist_dict, child.attrib['name'], template_path)
	
#}}}

def write_namelist(namelist_dict, outfilename, infilename):#{{{
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

#	for record, opts in namelist_dict.items():
#		out_namelist.write('&%s\n'%(record))
#
#		for opt, val in opts.items():
#			out_namelist.write('    %s = %s\n'%(opt.strip(), val[0].strip()))
#
#		out_namelist.write('/\n\n')
	out_namelist.close()
#}}}
#}}}

# *** Streams setup functions *** #{{{

def generate_streams_files(config_file, case_path, template_streams, templates_path):#{{{
	config_tree = ET.parse(config_file)
	config_root = config_tree.getroot()

	for streams in config_root.iter('streams'):
		streams_filename = '%s/%s'%(case_path, streams.attrib['name'])

		streams_tree = ET.parse(template_streams)
		streams_root = streams_tree.getroot()

		configure_streams_file(streams_root, templates_path, streams)
		write_streams_file(streams_root, config_file, streams_filename, '%s'%(case_path))

		del streams_root
		del streams_tree
#}}}

def flush_all_streams(streams):#{{{
	for stream in streams.findall('stream'):
		streams.remove(stream)

	for stream in streams.findall('immutable_stream'):
		streams.remove(stream)
#}}}

def flush_mutable_streams(streams):#{{{
	for stream in streams.findall('stream'):
		name = stream.attrib['name']
		streams.remove(stream)
#}}}

def flush_immutable_streams(streams):#{{{
	for stream in streams.findall('immutable_stream'):
		streams.remove(stream)
#}}}

def modify_stream_definition(streams_file, stream_conf):#{{{
	name_to_modify = stream_conf.attrib['name']

	found = False

	# Check if stream already exists:
	for stream in streams_file.findall('immutable_stream'):
		name = stream.attrib['name']
		if name.strip() == name_to_modify.strip():
			if not found:
				found = True
				stream_to_modify = stream
			else:
				print "ERROR: Stream %s found multiple times in template. Exiting..."%(name.strip())
				quit(1)

	for stream in streams_file.findall('stream'):
		name = stream.attrib['name']
		if name.strip() == name_to_modify.strip():
			if not found:
				found = True
				stream_to_modify = stream
			else:
				print "ERROR: Stream %s found multiple times in template. Exiting..."%(name.strip())
				quit(1)

	# If not found, need to create it
	if not found:
		found = True
		stream_to_modify = ET.SubElement(streams_file, 'stream')
		stream_to_modify.set('name', name_to_modify)

	for child in stream_conf.iter('*'):
		if child.tag == 'attribute':
			attr_name = child.attrib['name']
			attr_val = child.text
			stream_to_modify.set(attr_name, attr_val)
		elif child.tag == 'add_contents':
			for member in child.findall('member'):
				member_name = member.attrib['name']
				member_type = member.attrib['type']
				sub_member = ET.SubElement(stream_to_modify, member_type)
				sub_member.set('name', member_name)
				if 'packages' in member.attrib.keys():
					member_packages = member.attrib['packages']
					sub_member.set('packages', member_packages)
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

def apply_stream_template(streams_file, template_name, template_path):#{{{
	template_file = '%s/%s.xml'%(template_path, template_name)

	template_tree = ET.parse(template_file)
	template_root = template_tree.getroot()

	for streams in template_root.findall('streams'):
		for stream in streams.findall('stream'):
			modify_stream_definition(streams_file, stream)

	del template_tree
	del template_root
#}}}

def configure_streams_file(streams_file, template_path, streams_tag):#{{{
	keep_mode = streams_tag.attrib['keep']

	if ( keep_mode.strip() == 'immutable' ):
		flush_mutable_streams(streams_file)
	if ( keep_mode.strip() == 'none' ):
		flush_all_streams(streams_file)
	if ( keep_mode.strip() == 'mutable' ):
		flush_immutable_streams(streams_file)

	for template in streams_tag.findall('template'):
		template_name = template.attrib['name']
		apply_stream_template(streams_file, template_name, template_path)

	for stream in streams_tag.findall('stream'):
		modify_stream_definition(streams_file, stream)
#}}}

def write_streams_file(streams, config_file, filename, init_path):#{{{
	config_tree = ET.parse(config_file)
	config_root = config_tree.getroot()

	stream_file = open(filename, 'w')

	stream_file.write('<streams>\n')

	for stream in streams.findall('immutable_stream'):
		stream_name = stream.attrib['name']

		stream_file.write('\n')
		stream_file.write('<immutable_stream name="%s"'%(stream_name))
		for attr, val in stream.attrib.items():
			if ( not attr.strip() == 'name' ):
				stream_file.write('\n                  %s="%s"'%(attr, val))

		stream_file.write('/>\n')

	for stream in streams.findall('stream'):
		stream_name = stream.attrib['name']

		stream_file.write('\n')
		stream_file.write('<stream name="%s"'%(stream_name))

		for attr, val in stream.attrib.items():
			if ( not attr.strip() == 'name' ):
				stream_file.write('\n        %s="%s"'%(attr, val))

		stream_file.write('>\n\n')

		for substream in stream.findall('stream'):
			substream_name = substream.attrib['name']
			if 'packages' in substream.attrib.keys():
				package_name = substream.attrib['packages']
				entry = '\t<stream name="%s"'%(substream_name) + ' packages="%s" '%(package_name) +'/>\n'
                        else:
				entry = '\t<stream name="%s"'%(substream_name) +'/>\n'
			stream_file.write(entry)

		for var_struct in stream.findall('var_struct'):
			var_struct_name = var_struct.attrib['name']
			if 'packages' in var_struct.attrib.keys():
				package_name = var_struct.attrib['packages']
				entry = '\t<var_struct name="%s"'%(var_struct_name) + ' packages="%s" '%(package_name) +'/>\n'
                        else:
				entry = '\t<var_struct name="%s"'%(var_struct_name) +'/>\n'
			stream_file.write(entry)

		for var_array in stream.findall('var_array'):
			var_array_name = var_array.attrib['name']
			if 'packages' in var_array.attrib.keys():
				package_name = var_array.attrib['packages']
				entry = '\t<var_array name="%s"'%(var_array_name) + ' packages="%s" '%(package_name) +'/>\n'
                        else:
				entry = '\t<var_array name="%s"'%(var_array_name) +'/>\n'
			stream_file.write(entry)

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

# *** General Utility Functions *** #{{{
def add_links(config_file, args, configs):#{{{
	config_tree = ET.parse(config_file)
	config_root = config_tree.getroot()

	case = config_root.attrib['case']

	dev_null = open('/dev/null', 'r+')

	base_path = '%s/%s/%s/%s'%(args.core, args.configuration, args.resolution, case)

	for link in config_root.findall('add_link'):
		source = link.attrib['source']
		dest = link.attrib['dest']
		old_cwd = os.getcwd()
		os.chdir(base_path)
		subprocess.check_call(['ln', '-sf', '%s'%(source), '%s'%(dest)], stdout=dev_null, stderr=dev_null)
		os.chdir(old_cwd)
		del source
		del dest


	for executable in config_root.findall('add_executable'):
		source_attr = executable.attrib['source']
		dest = executable.attrib['dest']
		if not configs.has_option("executables", source_attr):
			print 'ERROR: Configuration %s requires a definition of %s.'%(config_file, source_attr)
			quit(1)
		else:
			source = configs.get("executables", source_attr)

		subprocess.check_call(['ln', '-sf', '%s'%(source), '%s/%s'%(base_path, dest)], stdout=dev_null, stderr=dev_null)
		del source_attr
		del source
		del dest

	del config_tree
	del config_root
	dev_null.close()
#}}}

def get_case_mode(config_file):#{{{
	config_tree = ET.parse(config_file)
	config_root = config_tree.getroot()

	mode = config_root.attrib['mode']

	del config_tree
	del config_root

	return mode
#}}}

def make_case_dir(config_file, base_path):#{{{
	config_tree = ET.parse(config_file)
	config_root = config_tree.getroot()

	case_name = config_root.attrib['case']

	if not os.path.exists('%s/%s'%(base_path, case_name)):
		os.makedirs('%s/%s'%(base_path, case_name))

	del config_root
	del config_tree

	return case_name
#}}}

def get_defined_files(config_file, init_path, args, configs):#{{{
	config_tree = ET.parse(config_file)
	config_root = config_tree.getroot()
	dev_null = open('/dev/null', 'w')

	for get_file in config_root.findall('get_file'):
		dest_name = get_file.attrib['dest']

		for mirror in get_file.findall('mirror'):
			try:
				protocol = mirror.attrib['protocol']
			except:
				print "Mirror is missing the 'protocol' attribute."
				print "Exiting..."
				quit(1)

			if protocol == 'local':

				try:
					path = '%s/%s'%( configs.get('paths', mirror.attrib['path_name']), mirror.attrib['file_name'])
				except:
					print " Mirror with protocol 'local' is missing one or more attribute of 'path_name' or 'file_name'"
					print " Or paths.%s doesn't exist in the configuration file."%(mirror.attrib['path_name'])
					print " Exiting..."
					quit(1)

				try:
					subprocess.check_call(['cp', '-f', '%s'%(path), '%s/%s'%(init_path, dest_name)], stdout=dev_null, stderr=dev_null)
				except:
					print "  -- Local mirror attempt failed. Trying other mirrors..."

			elif protocol == 'wget':
				if not args.no_download:
					try:
						name = mirror.attrib['file_name']
						path = '%s/%s'%( mirror.attrib['url'], mirror.attrib['file_name'])
					except:
						print " Mirror with protocol 'wget' is missing one or more attribute of 'url' or 'file_name'"
						print " Exiting..."
						quit(1)

					try:
						subprocess.check_call(['wget', '-q', '%s'%(path)], stdout=dev_null, stderr=dev_null)
						subprocess.check_call(['mv', '%s'%(name), '%s/%s'%(init_path, dest_name)], stdout=dev_null, stderr=dev_null)
					except:
						print "  -- Web mirror attempt failed. Trying other mirrors..."

	
		try:
			expected_hash = get_file.attrib['hash']
			validate_file = True
		except:
			validate_file = False

		if validate_file:
			if os.path.exists('%s/base_mesh.nc'%(init_path)):
				nc = netCDF4.Dataset('%s/base_mesh.nc'%(init_path), 'r')
				try:
					file_hash = nc.file_id
				except:
					nc.close()
					print " base_mesh does not have file_id attribute. Exiting..."
					quit(1)

				nc.close()

				if not file_hash.strip() == expected_hash.strip():
					print "*** ERROR: Base mesh has hash of '%s' which does not match expected hash of '%s'."%(file_hash, expected_hash)

	del config_tree
	del config_root
	dev_null.close()
#}}}

def generate_run_scripts(config_file, init_path):#{{{
	config_tree = ET.parse(config_file)
	config_root = config_tree.getroot()
	dev_null = open('/dev/null', 'r+')

	for run_script in config_root.findall('run_script'):
		script_name = run_script.attrib['name']
		script_path = "%s/%s"%(init_path, script_name)
		script = open(script_path, "w")

		script.write("#!/usr/bin/env python\n")
		script.write("import subprocess\n")

		for step in run_script.findall('step'):
			executable = step.attrib['executable']

			script.write("\n")
			script.write("# Run command is:\n")

			comment = "# %s "%(executable)
			command = "subprocess.check_call(['%s'"%(executable)
			
			for argument in step.findall('argument'):
				flag = argument.attrib['flag']
				val = argument.text

				if not flag.strip() == "":
					comment = "%s %s"%(comment, flag)
					command = "%s, '%s'"%(command, flag)

				comment = "%s %s"%(comment, val)
				command = "%s, '%s'"%(command, val)

			command = "%s])"%(command)
			script.write("%s\n"%(comment))
			script.write("%s\n"%(command))

			script.write("\n");

		script.close()

		subprocess.check_call(['chmod', 'a+x', '%s'%(script_path)], stdout=dev_null, stderr=dev_null)


	dev_null.close()
	del config_tree
	del config_root
#}}}
#}}}

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-o", "--core", dest="core", help="Core that contains configurations", metavar="CORE", required=True)
parser.add_argument("-c", "--configuration", dest="configuration", help="Configuration to setup", metavar="CONFIG", required=True)
parser.add_argument("-r", "--resolution", dest="resolution", help="Resolution of configuration to setup", metavar="RES", required=True)
parser.add_argument("-f", "--config_file", dest="config_file", help="Configuration file for test case setup", metavar="FILE", required=True)
parser.add_argument("--no_download", dest="no_download", help="If set, script will not auto-download base_mesh files", action="store_true")
parser.add_argument("--copy_instead_of_symlink", dest="nolink_copy", help="If set, script will replace symlinks with copies of files.", action="store_true")

args = parser.parse_args()

config = ConfigParser.SafeConfigParser()
config.read(args.config_file)

if not args.no_download:
	args.no_download = False

if not args.nolink_copy:
	args.nolink_copy = False

# Setup each xml file in the configuration directory:
template_path = '%s/templates/%s'%(os.getcwd(), args.core)
base_path = '%s/%s/%s'%(args.core, args.configuration, args.resolution)
for file in os.listdir('%s'%(base_path)):
	if fnmatch.fnmatch(file, '*.xml'):
		config_file = '%s/%s'%(base_path, file)

		case_dir = make_case_dir(config_file, base_path)
		case_mode = get_case_mode(config_file)

		case_path = '%s/%s'%(base_path, case_dir)

		if not config.has_option("streams", case_mode) or not config.has_option("namelists", case_mode):
			print "Error. Configuration file %s is requires paths for streams and namelist files for %s mode."%(args.config_file, case_mode)
			quit(1)
		else:
			case_streams_template = config.get("streams", case_mode)
			case_namelist_template = config.get("namelists", case_mode)
			
		generate_namelist_files(config_file, case_path, case_namelist_template, template_path)

		generate_streams_files(config_file, case_path, case_streams_template, template_path)

		add_links(config_file, args, config)

		get_defined_files(config_file, '%s'%(case_path), args, config)

		generate_run_scripts(config_file, '%s'%(case_path))

		print " -- Set up case: %s/%s"%(base_path, case_dir)

