#!/usr/bin/python

# Creates complete list with space group symbols, including version 
# with different origins and axes
# The input file is taken from cctbx toolkit from directory:
#  cctbx_sources/cctbx/sgtbx/symbols.cpp
#
# The notation symbol:{R,H} for different axes
# 	       symbol:{1,2} for different origins
#
# The default name of the group (withot :{R,H,1,2} extension) is also 
# preserved, and it points to :H for different axes and :2 for different 
# origins

input = "cctbx2dat-orgins.in"

file = open(input,'r')
counter = 1
while 1:
	line = file.readline()
	if line == "":break
	line = line.split(',')
	# HM symbol
	symbol = line [0]
	symbol = symbol.split()
	# check if the group has two different axes or orgins
	mult = line[1]
	mult = mult.split("\\0")
	mult = mult[:-1]

	newSymbol = ""
	for i in symbol:
		newSymbol = newSymbol + i
	
	print "%5i   %s"%(counter,newSymbol)
	counter += 1

	if len(mult) >1: # two different orgins of two axes
		elem = mult[0].strip()
		if (elem[0] == "R") or (elem[0:2] == "-R"):
			# Two axes
			newSymbolH = newSymbol+":H" # hexagonal (default)
			newSymbolR = newSymbol+":R" # rhombohedral
			print "%5i   %s"%(counter,newSymbolH)
			counter += 1
			print "%5i   %s"%(counter,newSymbolR)
			counter += 1
		else:
			# Two different origins
			newSymbol1 = newSymbol+":1" # origin 1
			newSymbol2 = newSymbol+":2" # origin 2 (default)
			print "%5i   %s"%(counter,newSymbol1)
			counter += 1
			print "%5i   %s"%(counter,newSymbol2)
			counter += 1

		

