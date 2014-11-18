# Formats json file with spacegroup symbol and corresponding symmetry
# operations.

# Input : formated_groups.dat:
#          1  FM-3M
#          2  P 1
#          ....
#          where first element is line number, second space group symbol
#
#
# Uses cctbx toolkit:
#  http://cctbx.sourceforge.net/
# 
# Based on: cctbx_sources/cctbx/examples/getting_started.py 
# in cctbx toolkit file structure
# 
# The complete set of spacegroup symbols can be found in:
#  cctbx_sources/cctbx/sgtbx/symbols.cpp
# in cctbx toolkit file structure
#
# To run the script, cctbx toolkit must be installed:
# 1) source /build/setpaths.sh  in cctbx toolkit (after its compiled)
# 2) libtbx.python format_Spacegroup_json.py
#
# Jakub Kaminski, UCLA, 2014
#

from cctbx import crystal
import json

file = open("formated_groups.dat",'r')

jsonData = {}

while 1:
	line = file.readline()
	if line == "": break
	line = line.split()
	symbol = line[1]
	group = crystal.symmetry(space_group_symbol=symbol)
	operations = []
	for operation in group.space_group():
		operations.append(operation.as_xyz())
	jsonData[symbol] = operations

with open("spacegroups.json",'w') as outfile:
	json.dump(jsonData, outfile, indent = 4)

