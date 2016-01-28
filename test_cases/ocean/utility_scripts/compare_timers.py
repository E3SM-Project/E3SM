#!/usr/bin/env python

import sys, os, fnmatch
import argparse
import re

def find_timer_value(timer_name, directory):#{{{
	# Build a regular expression for any two characters with a space between them.
	regex = re.compile(r'(\S) (\S)')

	sub_timer_name = timer_name.replace(' ', '_')

	timer = 0.0
	timer_found = False
	for file in os.listdir(directory):
		if not timer_found:
			process_file = False
			# Compare files written using built in MPAS timers
			if fnmatch.fnmatch(file, "log.*.out"):
				timer_line_size = 6
				name_index = 0
				total_index = 1
				process_file = True
			# Compare files written using GPTL timers
			elif fnmatch.fnmatch(file, "timing.*"):
				timer_line_size = 6
				name_index = 0
				total_index = 3
				process_file = True

			if process_file:
				stats_file = open('%s/%s'%(directory, file), "r")
				for block in iter(lambda: stats_file.readline(), ""):
					new_block = regex.sub(r"\1_\2", block[2:])
					new_block_arr = new_block.split()
					if len(new_block_arr) >= timer_line_size:
						if sub_timer_name.find(new_block_arr[name_index]) >= 0:
							try:
								timer = timer + float(new_block_arr[total_index])
								timer_found = True
							except ValueError:
								timer = timer
						del new_block
					del new_block_arr

	return timer_found, timer
#}}}

# Define and process input arguments
parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-b', '--base_directory', dest="base_directory", help="Directory with the baseline timer information.", required=True)
parser.add_argument('-c', '--comparison_directory', dest="comparison_directory", help="Directory with the comparison timer information.", required=True)
parser.add_argument('-t', '--timer', dest="timer", help="Name of the timer to compare", required=True)
parser.add_argument('-s', '--speedup', dest="speedup", help="If set, only speedup will be printed. This is useful when making speedup plots.", action="store_true")

args = parser.parse_args()

timer1_found, timer1 = find_timer_value(args.timer, args.base_directory)
timer2_found, timer2 = find_timer_value(args.timer, args.comparison_directory)

if timer1_found and timer2_found:
	try:
		speedup = timer1 / timer2
	except:
		speedup = 1.0

	percent = (timer2 - timer1) / timer1

	if not args.speedup:
		print "Comparing timer %s:"%(args.timer)
		print "             Base: %lf"%(timer1)
		print "          Compare: %lf"%(timer2)
		print "   Percent Change: %lf%%"%(percent*100)
		print "          Speedup: %lf"%(speedup)
	else:
		print "%lf"%(speedup)
