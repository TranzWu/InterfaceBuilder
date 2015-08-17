#!/usr/bin/python

import numpy as np
import fractions
import math as m
import os
import decimal as dc
import sys
import json

#
# Interface Builder 
# Author: Jakub Kaminski, UCLA 2014
#

if len(sys.argv) <=1:
	print "Usage:\n%s optionsFile"%sys.argv[0]
	exit()

inputFile = sys.argv[1]

def readInput(inputFile):

	file = open(inputFile,'r')
	line = file.readline()
	line = line.split()
	subCIF = line[1]

	line = file.readline()
	line = line.split()
	subMillerString = line[1]

	line = file.readline()
	line = line.split()
	depCIF = line[1]

	line = file.readline()
	line = line.split()
	depMillerString = line[1]
	
	line = file.readline()
	line = line.split()
	maxArea = float(line[1])
	
	line = file.readline()
	line = line.split()
	areaThres = float(line[1])

	line = file.readline()
	line = line.split()
	vecThres = float(line[1])

	line = file.readline()
	line = line.split()
	angleThres = float(line[1])

	line = file.readline()
	line = line.split()
	if line[1] == "None":
		capAtmS = None
	else:
		capAtmS = line[1]

	line = file.readline()
	line = line.split()
	if line[1]  == "None":
		capAtmD = None
	else:
		capAtmD = line[1]

	line = file.readline()
	line = line.split()
	fparam = float(line[1])

	line = file.readline()
	line = line.split()
	nLS = float(line[1])

	line = file.readline()
	line = line.split()
	nLD = float(line[1])

	line = file.readline()
	line = line.split()
	nConf = int(line[1])

	line = file.readline()
	line = line.split()
	subAtRad = float(line[1])

	line = file.readline()
	line = line.split()
	depAtRad = float(line[1])

	line = file.readline()
	line = line.split()
	tmp = line[1]
	if tmp == "True":
		skipStep1 = True
	elif tmp == "False":
		skipStep1 = False
	else:
		print "Value for skipStep1 not recognized"
		exit()
	# Read Poisson ratio
	line = file.readline()
	line = line.split()
	tmp = line[1]
	if tmp == "True":
		poissonRatio = True
	elif tmp == "False":
		poissonRatio = False
	else:
		print "Value for poissonRatio not recognized"
		exit()

	# Read if we want to split Substrate
	line = file.readline()
	line = line.split()
	tmp = line[1]
	if tmp == "True":
		sandwich = True
	elif tmp == "False":
		sandwich = False
	else:
		print "Value for sandwich not recognized"
		exit()

	# Read nVac
	line = file.readline()
	line = line.split()
	nVac= int(line[1])
	if nVac <= 1: nVac = 1



	return subCIF, subMillerString,\
	       depCIF, depMillerString,\
	       maxArea, areaThres, vecThres, angleThres,\
	       capAtmS, capAtmD, fparam, nLS, nLD, nConf, subAtRad, depAtRad,\
	       skipStep1, poissonRatio, sandwich, nVac

def getMillerFromString(millerString):

	out = np.array((0,0,0))
	for item in range(len(millerString)):
		out[item] = millerString[item]
	
	return out

def createMillerList(maxMillerInd):
        # create a list that contains strings for all the possible Miller indices
        # from "0 0 1" to "maxMillerInd maxMillerInd maxMillerInd"

        MillerList = list() 

	# If the length of maxMillerInd is 3 (i.e. 100, 111, 212) 
	# it means we are dealing with the indyvidual Miller orienation,
	# hence don't generate it, just return it value packed in the list
	if len(maxMillerInd) >= 3:
		MillerList.append(maxMillerInd)
		return MillerList
	
	# Generate all possible combinations up to specified Miller index 
        x1 = 0
        x2 = 0
        x3 = 0
	maxMillerInd = int(maxMillerInd)
        while x1 <= maxMillerInd:
                x2 = 0
                while x2 <= maxMillerInd:
                        x3 = 0
                        while x3 <= maxMillerInd:
                                MillerList.append(str(x1)+str(x2)+str(x3))
                                x3 = x3+1
                        x2 = x2+1
                x1 = x1+1

        MillerList.remove("000") # we cannot have "0 0 0"

        return MillerList

def uqLabels(labels,types):
        # Find unique labels in the list
        # Return dictionary with key=integer number, value=labels
        tmp = []
	# Get all the unique atoms
        for i in labels:
                if i not in tmp:
                        tmp.append(i)

	# Get all the atoms that are in dictionary to the list
	vals = types.values()

	nextIdx = 0
	if types != {}:
		maxTyp = max(types)
		nextIdx = maxTyp + 1

	# Put the atoms to the dictionary
	for i in tmp:
		if i not in vals:
			types[nextIdx] = i
			nextIdx += 1

        return types

def ReadCIF(filename,atomTypes):
	# Reads information from CIF file:
	# - ICSD database code for the material
	# - unit cell lengths and angles
	# - evaluates xyz coordinates based on fractional coordinates and 
	# symmetry position of unique atoms
	# Returns:
	# - ICSD database code for the material
	# - unit cell lengths and angles
	# - list with atom labels
	# - array with xyz coordinates

	file = open(filename,'r')

	lines = file.readlines()
	file.close()

	# Strip whitespaces and line endings
	for i in range(len(lines)):
		lines[i]=lines[i].strip()

	for line in lines:
		if line.find("data") >= 0:
			MatID = line

		elif line.find("_cell_length_a") >= 0:
			line = line.split()
			cellA =  CheckParentheses(line[1])

		elif line.find("_cell_length_b") >= 0:
			line = line.split()
			cellB =  CheckParentheses(line[1])

		elif line.find("_cell_length_c") >= 0:
			line = line.split()
			cellC =  CheckParentheses(line[1])

		elif line.find("_cell_angle_alpha") >= 0:
			line = line.split()
			alpha = CheckParentheses(line[1])

		elif line.find("_cell_angle_beta") >= 0:
			line = line.split()
			beta = CheckParentheses(line[1])

		elif line.find("_cell_angle_gamma") >= 0:
			line = line.split()
			gamma = CheckParentheses(line[1])
		elif line.find("_symmetry_space_group_name_H-M") >= 0:
			SGsymbol = extractSpaceGroup(line)

	f2cmat = FracToCart(cellA,cellB,cellC,alpha,beta,gamma)
	print "SPACE GROUP: ",SGsymbol

	# READ ATOM SYMMETRY POSITION CORRESPONDING TO SPACE GROUP
	sgInput = open('spacegroups.json','r')
	sgData = json.load(sgInput)
	try:
		atompos = sgData[SGsymbol]
	except KeyError:
		print "Unrecognized spacegroup symbol read from CIF file %s"%SGsymbol
		exit()

	# READ ATOM FRACTIONAL COORDINATES
	# Find index of the the beginning of atom fractional coordinates
	# Find either _atom_site_label or _atom_site_symbol and see what is first
	# Set the ib to the fist one
	try:
		ibL = lines.index("_atom_site_label")
	except ValueError:
		ibL = 100000 # Just put big unresonable value

	try:
		ibS = lines.index("_atom_site_type_symbol")
	except ValueError:
		ibS = 100000 # Just put big unresonable value

	# Find where is the begining of atom definition loop in CIF file.
	if ibL < ibS:
		ib = ibL
	else:
		ib = ibS

	allLabelsRead = False
	ibSave = ib
	while not allLabelsRead:
		ib += 1
		if lines[ib][0] == "_":
			if lines[ib][0:18] == "_atom_site_fract_x":
				xPos = ib - ibSave
			if lines[ib][0:18] == "_atom_site_fract_y":
				yPos = ib - ibSave
			if lines[ib][0:18] == "_atom_site_fract_z":
				zPos = ib - ibSave
		else:
			allLabelsRead = True

	# Get atom positions
	# We don't know how many atoms ther is the structure
	# Read them unitl keyworld "loop", or line begining with "_"
	# or line begining with "#" or end of file is encountered

	frapos = []
	while 1:
		try:
			line = lines[ib]
		except IndexError:
			break # end of file reached 

		# If file finishes with empty line
		try:
			test = line[0]
		except IndexError:
			break

		if (line[0:4] != "loop") and (line[0] != "_") \
		   and (line[0] != "#"):
			frapos.append(line)
			ib += 1
		else:
			break

	# Now, knowing fractional coordinates and the symmetry positions, 
	# build fractional xyz coordinates.
	# First we will put fractional coordinates into array fracoord

	fracoord = np.zeros((len(frapos),3))
	atomlabels = []
	ii = 0 #index
	for atom in frapos:
		atom = atom.split()
		symbol = atom[0]
		# Strip the digits from the symbol name:
		if symbol[-1].isdigit():
			symbol = symbol[0:-1]

		x = float(CheckParentheses(atom[xPos]))
		y = float(CheckParentheses(atom[yPos]))
		z = float(CheckParentheses(atom[zPos]))

		fracoord[ii,0] = x
		fracoord[ii,1] = y
		fracoord[ii,2] = z
		atomlabels.append(symbol)
		ii = ii+1
	# Create array with fractional coordinates based on
	# symmetry positions for each atom type
#	positions = np.zeros((len(fracoord)*len(atompos),3))

	iatlab = 0 # index of symmetry unique atom
	atoms = []# list to be fileld with atom indexes after symmetry operation
	for atom in fracoord:
		postmp = []
		for pos in atompos: #For each element in atompos, convert it 
			            # to fraction, multiply by fractional 
				    # coordinate and save to positions array
#			pos = pos.split()
			pos = pos.split(",") # For group
			coordtmp = []
#			for elem in pos[1:]:
			for elem in pos:     # For group
				# Read symmetry operation and constuct 
				# fractional coordiante

				s,o = ReadSym(elem)

				coord = s[0]*o[0]*atom[0] + \
					s[1]*o[1]*atom[1] + \
					s[2]*o[2]*atom[2] + \
					s[3]*o[3]             

				# Normalise to [0,1)
				if coord < 0 : coord = 1+coord 
				if coord > 1 : coord = coord-1
# WARNING
				# Round to 4 digits to avoid rounding error
# END WARNING
				coord = round(coord,3)
				coordtmp.append(coord)
			#Find duplicates	
			if coordtmp not in postmp:
				postmp.append(coordtmp)
				atoms.append(atomlabels[iatlab])

		# If first atom or only one atom, create position array
		if iatlab == 0:
			positions = np.array(postmp)				
		else: # add 2nd and other atoms to position array
			tmp = np.array(postmp)
			positions = np.vstack([positions,tmp])
		
	# Find and add 1 to every 0 in postmp. Add new coordiante to positions
		for elem in postmp:
			if elem[0] == 0:
				newrow = np.array(elem)
				newrow[0] += 1
				diff = CmpRows(positions,newrow)
				if diff:
					positions = \
						np.vstack([positions,newrow])
					atoms.append(atomlabels[iatlab])
			if elem[1] == 0:
				newrow = np.array(elem)
				newrow[1] += 1
				diff = CmpRows(positions,newrow)
				if diff:
					positions = \
						np.vstack([positions,newrow])
					atoms.append(atomlabels[iatlab])
			if elem[2] == 0:
				newrow = np.array(elem)
				newrow[2] += 1
				diff = CmpRows(positions,newrow)
				if diff:
					positions = \
						np.vstack([positions,newrow])
					atoms.append(atomlabels[iatlab])
			if elem[0] == 0 and elem[1] == 0:
				newrow = np.array(elem)
				newrow[0] += 1
				newrow[1] += 1
				diff = CmpRows(positions,newrow)
				if diff:
					positions = \
						np.vstack([positions,newrow])
					atoms.append(atomlabels[iatlab])
			if elem[0] == 0 and elem[2] == 0:
				newrow = np.array(elem)
				newrow[0] += 1
				newrow[2] += 1
				diff = CmpRows(positions,newrow)
				if diff:
					positions = \
						np.vstack([positions,newrow])
					atoms.append(atomlabels[iatlab])
			if elem[1] == 0 and elem[2] == 0:
				newrow = np.array(elem)
				newrow[1] += 1
				newrow[2] += 1
				diff = CmpRows(positions,newrow)
				if diff:
					positions = \
						np.vstack([positions,newrow])
					atoms.append(atomlabels[iatlab])
			if elem[1] == 0 and elem[1] == 0 and elem[2] == 0:
				newrow = np.array(elem)
				newrow[0] += 1
				newrow[1] += 1
				newrow[2] += 1
				diff = CmpRows(positions,newrow)
				if diff:
					positions = \
						np.vstack([positions,newrow])
				atoms.append(atomlabels[iatlab])
		iatlab += 1 # next atom

        #Convert atoms list to numpy array with atom types. Atom types will be hold in dictonary
        # Get the dictionary
        atomTypes = uqLabels(atoms,atomTypes)
        # Convert atoms list to array
        ii = 0
        atomstmp = np.zeros(len(atoms))
        for label in atoms:
                for uqLabIdx,uqLab in atomTypes.iteritems():
                        if label == uqLab:
                                atomstmp[ii] = uqLabIdx
                ii+=1
        atoms = atomstmp

	positions = np.dot(positions,f2cmat)
	#TODO: check if this is always true.
	# This is done for the cases when if CIF file the coordinates of first atom are different from 0.0 0.0 0.0
	# Such cases prevents proper rotations on def plane, where the desired plane doesnt have z=0
	shift = positions[0].copy()
	positions -= shift
	#end TODO
	return MatID,f2cmat,atoms,positions,atomTypes

def CmpRows(mat1,vector):

	diff = True
	for i in mat1:
		if list(i) == list(vector):
			diff = False
	return diff

def FracToCart(a,b,c,alpha,beta,gamma):
	# Calculate transformation matrix to transform fractional coordinates
	# to cartesian

	# Degrees to radians

	alpha = alpha*m.pi/180
	beta = beta*m.pi/180
	gamma = gamma*m.pi/180

	cosa = m.cos(alpha)
	cosb = m.cos(beta)
	cosg = m.cos(gamma)

	cos2a = cosa**2
	cos2b = cosb**2
	cos2g = cosg**2

	sing = m.sin(gamma)

	vol = m.sqrt(1-cos2a-cos2b-cos2g + (2*cosa*cosb*cosg))

	mat = np.ones((3,3))

	mat[0,0] = a
	mat[0,1] = b*cosg
	mat[0,2] = c*cosb

	mat[1,0] = 0.0
	mat[1,1] = b*sing
	mat[1,2] = c*(cosa - (cosb*cosg))/sing

	mat[2,0] = 0.0
	mat[2,1] = 0.0
	mat[2,2] = c*vol/sing

	mat = CleanMatElements(mat)

	# Transpose transformation matrix, so the vectors are aligned in
	# columns, not rows, i.e
	# mat = [a1,a2,a3
	#        b1,b2,b3
	#        c1,c2,c3]
	mat = mat.T
	
	return mat

def FracToCart2(a,b,c,alpha,beta,gamma):
	# Calculate transformation matrix to transform fractional coordinates
	# to cartesian

	# Degrees to radians

	alpha = alpha*m.pi/180
	beta = beta*m.pi/180
	gamma = gamma*m.pi/180

	cosa = m.cos(alpha)
	cosb = m.cos(beta)
	cosg = m.cos(gamma)

	cos2a = cosa**2
	cos2b = cosb**2
	cos2g = cosg**2

	sinb = m.sin(beta)
	sing = m.sin(gamma)

	vol = m.sqrt(1-cos2a-cos2b-cos2g + (2*cosa*cosb*cosg))

	mat = np.ones((3,3))

	mat[0,0] = a
	mat[0,1] = b*cosg
	mat[0,2] = c*cosb

	mat[1,0] = 0.0
	mat[1,1] = b*sing
	mat[1,2] = c*(cosa - cosb*cosg)/sing

	mat[2,0] = 0.0
	mat[2,1] = 0.0
	mat[2,2] = m.sqrt((c**2)-cos2b-(c*(cosa-cosb*cosg)/sing)**2)

	mat = CleanMatElements(mat)
	
	return mat



def CleanMatElements(mat):

	# Clean numeric mess in numpy array, by replacing really small 
	# elements with 0
	# For instance:
	# input
	# [  0.00000000e+00   5.43053000e+00   3.32524059e-16] 
	# output
	# [  0.0   5.43053   0.0]

	small =  abs(mat) < 1e-12
	mat[small] = 0.0

	return mat

def ReadSym(string):
	# Reads the _symmetry_equiv_pos_as_xyz operation in CIF 
	# Assumes that general format is:
	# 'x+y+z+int1/int2'
	# Example: for 'x-y+1/2' it can be written in form of lists:
	# operation=[1,1,0,0.5]  x=1, y=1, z=0, frac=0.5
	# signs=[1,-1,1,1]   +x -y +z +fraction

	operation = [0,0,0,0]
	# operation[0]  -  1 if there is x symmetry operation, 0 otherwise
	# operation[1]  -  1 if there is y symmetry operation, 0 otherwise
	# operation[2]  -  1 if there is z symmetry operation, 0 otherwise
	# operation[3]  -  1 if there is fraction, 0 otherwise
	
	signs = [1,1,1,1] 
	# signs[0] - 'x' sign, -1 means -x
	# signs[1] - 'y' sign, -1 means -y
	# signs[2] - 'z' sign, -1 means -z
	# signs[3] - fraction sign, -1 means -fraction
	
	# Find x,y,z 
	idx = string.find("x")
	if idx > -1: 
		operation[0] = 1
		# Find sign before x
		if idx != 0: #check in the case fraction is first element
			if string[idx-1] == '-':
				signs[0] = -1

	idx = string.find("y")
	if idx > -1: 
		operation[1] = 1
		# Find sign before x
		if idx != 0: #check in the case fraction is first element
			if string[idx-1] == '-':
				signs[1] = -1

	idx = string.find("z")
	if idx > -1: 
		operation[2] = 1
		# Find sign before x
		if idx != 0: #check in the case fraction is first element
			if string[idx-1] == '-':
				signs[2] = -1


	# Find fraction
	idx = string.find("/")
	if idx > -1:
		f1idx = idx-1
		f2idx = idx+1
		frac = float(string[f1idx])/float(string[f2idx])
		operation[3] = frac
		
		# Find sign before fraction
		if f1idx != 0: #check in the case fraction is first element
			if string[f1idx-1] == '-':
				signs[3] = -1

	return signs,operation

def ReadSymFrac(string):
	# This is specyfic routine for CIF file, to transform symmetry 
	# equivalent position to actual true numbers.
	# Works for the cases when the symmetry operation is of the forms:
	# "x", "-x", "x+1/2", "-z+3/4"...

	# There can be only two operation op1*x+op2*ratio, where 
	# op1 and op2 = "+" or "-". Define op list as two "+" operations
	# 
	# Input:
	# - string - symmetry operation as from CIF file, eg. "x+1/2"
	# Output:
	# - coordid - number for coordinate in symmetry operation: 
	#             'x'=0, 'y'=1, 'z'=2
	# - op - signs in symmetry operation, e.g:
	#        "x+1/2" : op=[1,1]
	#        "-x+1/2": op=[-1,1]
	# 	 "-x-1/2": op=[-1,-1] 
	#        "-x"    : op=[-1,1] ....
	# - digits - the digits in the fraction, e.g.:
	#        "1/4" : digits=[1,4]
	#
	# Author: Jakub Kaminski, UCLA, 04/2013


	op = [1,1]

	# There can be two digits defining ratio, e.g. 1/4. 
	# Define digits list with two "0,1" for the start

	digits = [0,1]

	# Find what is the coordiante
	coorindex = string.find("x")
	coordid = 0
	if coorindex == -1:
		coorindex = string.find("y")
		coordid = 1
		if coorindex == -1:
			coorindex = string.find("z")
			coordid = 2
	
	# Check if the coordinate is not negative 
	# and mark it by -1 in op variable
	if coorindex > 0: # just in case to be sure we are not accessing beyond
		          # string length
		optype = string[coorindex-1]
		if optype == "-":
			op[0] = -1
	
	# Check if we add or substract ratio.
	# If we add, op[1]=1, if we substract op[1]=-1
	if coorindex <len(string): # just in case to be sure we are not accesing
		                   # beyond string length 
		optype = string[coorindex+1]
		if optype == "-":
			op[1] = -1
	

	# Now read the ratio
	# Find digit by digit and save them to digit index
	digitindex = 0
	for i in string:
		if i.isdigit():
			digits[digitindex] = float(i)
			digitindex = digitindex + 1
	
	return coordid, op, digits


def CheckParentheses(input):
	# In some CIF file the cell dimensions and anglesare given with 
	# standard deviation in final digits. It is usally the last 
	# number given in parentheses. 

	# CheckParentheses checks is string read from CIF has paretheses and 
	# disregards them.
	#
	# Variables:
	# input is as string
	# Returns float
	#
	# Author: Jakub Kaminski, UCLA, 04/2013

	i = input.find("(")

	if i> -1: # found
		input = input[:i]

	return float(input)

def extractSpaceGroup(symbol):
	# Extract spacegroup symbol from CIF file line
	symbol = symbol.strip("_symmetry_space_group_name_H-M")
	symbol = symbol.strip() # strip any leading spaces
	symbol = symbol.strip("'") # strip quote signs
	symbol = symbol.strip('"') # strip double quotes (just in case)

        # In case symbol contains "_", i.e "I4_1/amd"
        if symbol.find("_") > 0:
                tmp = ""
                for i in symbol:
                        if i != "_":
                                tmp += i
                symbol = tmp

	# In the case the symbol is given with white space, i.e. "F M -3 M"
	symbol = symbol.split()
	symbolOut = ""
	for i in symbol:
		symbolOut += i
	symbolOut = symbolOut.capitalize()

	# In ISCD "S" or ":S" extension on the end of the symbol is equivalent
	# to :1 notation for origin choice for space groups with two origins,
	# .i.e "F d -3 m S"
	# Our notoations follows :1 :2 extension for different origins and 
	# :R :H for different axes, so convert everything to this notation

	if symbolOut[-1] == "s": 
		symbolOut = symbolOut[:-1] +":1"
	
	if symbolOut[-2:] == ":s": 
		symbolOut = symbolOut[:-2] +":1"
	
	# Intel is using follwing convetion in their MadeA software to
	# label group with different origins
	# FD-3MO1 - origin 1; FD-3MO2 - origin 2

	if symbolOut[-2:] == "o1":
		symbolOut = symbolOut[:-2] +":1"

	if symbolOut[-2:] == "o2":
		symbolOut = symbolOut[:-2] +":2"

	# TODO: What are other convetions?

	return symbolOut

def unique2(a,labels):
	# Remove duplicate elements from array keeping it order
	# Not the fastest but will do for now
	# 3.5s to sort 13122 elements array

	dups = set()

	newlist = []
	newlabels = []
	ii = 0

	for elem in a:
		if str(elem) not in dups:
			newlist.append(elem)
			newlabels.append(labels[ii])
			dups.add(str(elem))
		ii += 1

	return np.array(newlist), newlabels

def fneigh(point,setp):

	# Find of nearest neighbours of the "point" in the "setp" of points
	# 
	# Return array of the shape, where
	# - 1st column is the distance of "point" to point p in "setp"
	# - 2nd column is the index of point p in "setp"

	neigh = np.zeros((len(setp),2))
	idx = 0
	for p in setp:
		d = np.linalg.norm(point-p)
		neigh[idx,0] = d
		neigh[idx,1] = idx
		idx += 1

	# Sort the distances according to the distance
	neigh = neigh[np.argsort(neigh[:,0])]

	neigh = CleanMatElements(neigh)
	
	return neigh

class Surface:
	def __init__(self,transM,positions,atoms,atomTyp,Midx):
		self.A = transM[0] # Vector A of the unit cell
		self.B = transM[1] # Vector B of the unit cell
		self.C = transM[2] # Vector C of the unit cell
		self.positions = positions # coordinates of material
		self.atoms = atoms # list with atom labels
		self.Midx = Midx # Miller indices of the plane
		#Define plane vectors and normal
		self.u = np.array((0,0,0))
		self.v = np.array((0,0,0))
		self.n = np.array((0,0,0))
		self.origin = np.array((0,0,0)) 
		self.a = np.array((0,0,0)) # primitive vector
		self.b = np.array((0,0,0)) # primitive vector

		self.unitcell = self.positions.copy() # save unit cell 
		self.unitcellatoms = atoms # save unit cell atom labels
		self.atomTyp = atomTyp
		self.primarea = 0 # area of the primitive cell
		self.phiab = 0 # angle between primtive vectors
		self.norma = 0 # norm of primitive vector a
		self.normb = 0 # norm of primitive vector b

		self.planepos = np.array([[]])#empty array for surface positions
		self.planeatms = [] # empty list for plane atom labels

		self.planeposblk = np.array([[]])# empty array for surface 
			         		 # positions and atoms below
		self.planeatmsblk = [] # empty list for plane atom labels
		                       # and atome below
		self.avecs = np.array((0,0,0)) # anticpatory vectors
		self.nneigh = 0 # number of nearst neigbours
		self.exists = False

	def __removeperiodic3D(self,pos,vec1,vec2,vec3):

		# Find the atoms in the superlattice that are
		# translations of the other atoms 

		# Rude and not effecient solution 
		# Try more through array interception

		r = range(2)
		uniqlist=[]
		poscheck = []
		for x in r:
			for y in r:
				for z in r:
					if x != 0 or y!= 0 or z!=0:
						if poscheck == []:
							poscheck = pos+x*vec1+y*vec2+z*vec3
						else:
							tmp = pos + x*vec1 + y*vec2+z*vec3
							poscheck = np.vstack\
							          ([poscheck,tmp])

		# Find indicies of unique elements in pos array that 
		# are not there is poscheck array
		ii = 0
		for i in pos:
			uq = True
			for j in poscheck:
				# Smaller accuracy required for z-axis
				# so round it to 5
				# Gives problems otherwise
				if round(i[0],8) == round(j[0],8) and \
				   round(i[1],8) == round(j[1],8) and \
				   round(i[2],8) == round(j[2],8):
					   uq = False
			if uq:
				if ii not in uniqlist:
					uniqlist.append(ii)
			ii += 1
		return uniqlist

        def bulkNEW(self,ncells):
                # 
                # OLD ROUTINE TO CONSTRUCT BULK MATERIAL USING 
                # __UNIQUE ROUTINE TO FIND DUPLICATE ELEMENTS.
                # THE PROBLEM WAS THAT IT DID NOT KEEP THE ORDER 
                # OF THE POSITION ARRAY
                # Consruct bulk surface from "ncells x unit cell"
                #
                # The results is saved in the positions variable

                r = range(-ncells,ncells+1)

                # Initial cooridnates and labels 
                posout = self.unitcell.copy()

                # Make unit cell periodic
                perIDX = self.__removeperiodic3D(self.unitcell,self.A,self.B,self.C)
                unitcellPer = self.unitcell.copy()[perIDX]
                unitcellatomsPer = self.unitcellatoms.copy()[perIDX]

                nElems = len(r)
                nAtoms = len(unitcellPer)
                nElems = nElems * nElems * nElems * nAtoms
                posout = np.zeros((nElems,3))
                newlabels = np.zeros((nElems))
                posout[0:nAtoms] = unitcellPer
                newlabels[0:nAtoms] = unitcellatomsPer

                if ncells> 2:
                        middle = len(r)/2
                        closeR = r[middle-2:middle+3]
                else:
                        closeR = r
                # Generate the middle of the bulk first. This is usefull in the later parts of the code to 
                # limit the size of the bulk to small number of atoms, for instance to look for nearest 
                # neighbors 
                i = nAtoms
                for x in closeR:
                        for y in closeR:
                                for z in closeR:
                                        if x != 0 or y!= 0 or z != 0:
                                                posout[i:nAtoms+i] = unitcellPer \
                                                           + x*self.A \
                                                           + y*self.B \
                                                           + z*self.C
                                                newlabels[i:nAtoms+i] = unitcellatomsPer
                                                i += nAtoms

                idxMiddle = i
                # Continue with the rest of bulk
                if ncells >2:
                        for x in r:
                                for y in r:
                                        for z in r:
                                                if x not in closeR or y not in closeR or z not in closeR:
                                                        posout[i:nAtoms+i] = unitcellPer \
                                                                   + x*self.A \
                                                                   + y*self.B \
                                                                   + z*self.C
                                                        newlabels[i:nAtoms+i] = unitcellatomsPer
                                                        i += nAtoms


#                ii = 0
#                for i in posout:
#                       #print newlabels[ii],i[0],i[1],i[2]
#                       print self.atomTyp[newlabels[ii]],i[0],i[1],i[2]
#                       ii+=1
                self.positions = posout
                self.atoms = newlabels
                self.positionsSmall = posout[0:idxMiddle]
                self.atomsSmall = newlabels[0:idxMiddle]

	def bulk(self,ncells):
		# 
		# OLD ROUTINE TO CONSTRUCT BULK MATERIAL USING 
		# __UNIQUE ROUTINE TO FIND DUPLICATE ELEMENTS.
		# THE PROBLEM WAS THAT IT DID NOT KEEP THE ORDER 
		# OF THE POSITION ARRAY
		# Consruct bulk surface from "ncells x unit cell"
		#
		# The results is saved in the positions variable

		r = range(-ncells,ncells+1)
		
		# Initial cooridnates and labels 
		posout = self.unitcell.copy()
		newlabels = self.__addatoms([])

		for x in r:
			for y in r:
				for z in r:
					if x != 0 or y!= 0 or z != 0:
						newpos = self.unitcell \
						           + x*self.A \
							   + y*self.B \
						           + z*self.C
						posout = \
						      np.vstack([posout,newpos])
						newlabels = \
						      self.__addatoms(newlabels)


		# Remove all duplicates

#		self.positions, self.atoms = self.__unique(posout,newlabels)
		posout = CleanMatElements(posout)
		self.positions, self.atoms = unique2(posout,newlabels)

	def bulkNAIVE(self,ncells):
		#
		# Consruct bulk surface from "ncells x unit cell"
		#
		# The results is saved in the positions variable

		r = range(-ncells,ncells+1)
		
		# Initial cooridnates and labels 
		posout = self.unitcell.copy()
		newlabels = self.__addatoms([])


		for x in r:
			for y in r:
				for z in r:
					if x != 0 or y!= 0 or z != 0:
						newpos = self.unitcell \
						           + x*self.A \
							   + y*self.B \
						           + z*self.C
						posout,newlabels = \
							self.__uniqueNAIVE\
							 (newpos,posout,\
							 newlabels)

		self.positions = posout
		self.atoms = newlabels
		
	def __addatoms(self, newlabels):
		# Add unit cell atoms labels to atom label list
		for i in range(len(self.unitcellatoms)):
			newlabels.append(self.unitcellatoms[i])
		return newlabels

	def __addatomsNEW(self, nl,base):
		# Add unit cell atoms labels to atom label list
		for i in range(len(base)):
			nl.append(base[i])
		return nl

	def __uniqueNAIVE(self,newmat,oldmat,labels):

		# VERY NAIVE WAY OF FINDING DUPLICATE ROWS IN ARRAY
		# MEMORY AND CPU TIME INEFFECIENT
		# HAS THE ADVANTAGE OF KEEPING THE ORDER
		# WILL WORK FOR NOW, BUT FIND BETTER WAY ASAP
		# 42s to sort 13122 elements array

		out = oldmat.copy()
		#duplicate = False
		ii = 0
		count= 0
		for elem1 in newmat:
			duplicate = False
			for elem2 in oldmat:
				if (elem1[0] == elem2[0]) and \
				   (elem1[1] == elem2[1]) and \
				   (elem1[2] == elem2[2]):
					duplicate = True
					break # break the for loop

			if not duplicate:
				count += 1
				out = np.vstack([out,elem1])
				labels.append(self.unitcellatoms[ii])
			ii += 1
			

		return out,labels

	def __unique(self,a,labels):

		# Check the numpy array with coordiantes for duplicate entries
		# Remove duplicates, remove lables correspoding to duplicated 
		# coordinates. 
		# The returned coordiantes are sorted in different order than 
		# orginals
		#
		# Subroutine take from:
		# http://stackoverflow.com/questions/8560440/removing-duplicate-columns-and-rows-from-a-numpy-2d-array

		# 0.5s to sort 13122 elements array

		order = np.lexsort(a.T)
		a = a[order]
		diff = np.diff(a, axis=0)
		ui = np.ones(len(a), 'bool')
		ui[1:] = (diff != 0).any(axis=1) 

		ii = 0
		newlabels = []
		for idx in order:
			if ui[ii] != False:
				newlabels.append(labels[idx])
			ii += 1

		return a[ui], newlabels

	def construct(self):
		h = self.Midx[0]
		k = self.Midx[1]
		l = self.Midx[2]

		if h != 0 and k == 0 and l == 0:
			# (100) surface
			self.u = self.B.copy()
			self.v = self.C.copy()

			# We take abs(h) as (-100) will be same as (100)
			self.origin = (1.0/abs(h)) * self.A

		elif h == 0 and k != 0 and l == 0:
			# (010) surfae
			self.u = self.A.copy()
			self.v = self.C.copy()
			
			# We take abs(k) as (0-10) will be same as (010)
			self.origin = (1.0/abs(k)) * self.B

		elif h == 0 and k == 0 and l != 0:
			# (001) surface
			self.u = self.A.copy()
			self.v = self.B.copy()

			# We take abs(l) as (00-1) will be same as (001)
			self.origin = (1.0/abs(l)) * self.C

		elif h != 0 and k != 0 and l == 0:
			# (hk0) surface

			# By symmetry (-h,-k,0) is equvalent to (h,k,0)
			if h < 0 and k < 0:
				h *= -1
				k *= -1
			
			self.u = (1./h)*self.A - (1./k)*self.B
			self.v = self.C.copy()

			if h > 0 and k > 0: #110
#				self.origin = (1./k)*self.B
				self.origin = (1./k)*self.B

			if h < 0 and k > 0: #(-110)
				self.origin = (1./k)*self.B + self.A

			if h > 0 and k <-1: #(1-10)
				self.origin = (1./abs(k))*self.B


		elif h != 0 and k == 0 and l != 0:
			# (h0l) surface

			# By symmetry (-h,-k,0) is equvalent to (h,k,0)
			if h < 0 and l < 0:
				h *= -1
				l *= -1

			self.u = (1./h)*self.A - (1./l)*self.C
			self.v = self.B.copy()

			if h > 0 and l > 0: #(101)
				self.origin = (1./l)*self.C

			if h < 0 and l > 0: #(-101)
				self.origin = (1./l)*self.C + self.A

			if  h > 0 and l < -1:
				self.origin = (1./abs(l))*self.C


		elif h == 0 and k != 0 and l!= 0:
			# (0kl) surface

			# By symmetry (-h,-k,0) is equvalent to (h,k,0)
			if k < 0 and l < 0:
				k *= -1
				l *= -1

			self.u = (1./k)*self.B - (1./l)*self.C
			self.v = self.A.copy()

			if k > 0 and l > 0: #(011)
				self.origin = (1./l)*self.C

			if k < 0 and l > 0: #(0-11)
				self.origin = (1./abs(l))*self.C + self.B

			if k > 0 and l < -1: #(01-1)
				self.origin = (1./abs(l))*self.C


		elif h != 0 and k != 0 and l != 0:
			# (hkl) surface

#			self.u = (1.0/gcd_hk) * (k*self.A - h*self.B)
#			self.v = (1.0/gcd_hl) * (l*self.A - h*self.C)

			# The equivalent planes
			# (-1-1-1)=(111)
			# (1-1-1)=(-111)
			# (-1-11)=(11-1)
			# (-11-1)=(1-11)
			# Find all equvalent planes
			counter = 0
			for idx in self.Midx:
				if idx < 0:
					counter += 1

			if counter >= 2:
				h *= -1
				k *= -1
				l *= -1

			gcd_hk = fractions.gcd(abs(h),abs(k))
			gcd_hl = fractions.gcd(abs(h),abs(l))

#			self.u = (1.0/gcd_hk) * (h*self.A - k*self.B)
#			self.v = (1.0/gcd_hl) * (h*self.A - l*self.C)


			self.u =  (1./h)*self.A - (1./k)*self.B
			self.v =  (1./l)*self.C - (1./k)*self.B

			if h > 0 and k > 0 and l > 0:
				self.origin = (1./k)*self.B

			if h < 0 and k > 0 and l > 0:
				self.origin = self.A + (1./k)*self.B

			if h > 0 and k < -1 and l > 0:
				self.origin = (1./abs(k))*self.B

			if h > 0 and k > 0 and l < 0:
				self.origin = self.C + (1./k)*self.B


	def plane(self):
		
		# Cut the plane from the postitions
		
		# Define normal to the plane
		self.n = normal(self.u,self.v)

		# Shift the orgin of the unit cell so it 
		# corresponds to the plane vectors
		self.positions = self.positions - self.origin

		# Rotate the surface to xy plane
		# Algorithm:
		# 1) Find angle of normal of the plane to the z axis
		# 2) Define vecor of rotation as cross product between
		#    normal and z axis
		# 3) Define rotation matrix for this angle and vector
		# 4) Rotate the plane

		# Find the angle
		normN = np.linalg.norm(self.n) # norm of plane normal

		z = np.array((0,0,1))
		normZ = 1 # norm of Z
		
		cosphi = np.dot(self.n,z)/(normN*normZ)
		phi = m.acos(cosphi)

		phi = phi * 180/m.pi
		
		# Create diagonal rotation matrix 
		R = np.diag((1.0,1.0,1.0))

		if phi != 0: # angle 0 gives problems with NaN in R matrix,
			     # but we dont need it anyways in this case
			phi = 180 - phi # Rotation is counterclockwise, 
			                #so define angle as 180 - phi
			phi = phi * m.pi/180

			cosphi = m.cos(phi)
			sinphi = m.sin(phi)

			rv = np.cross(self.n, z) # rotation vector
			normrv = np.linalg.norm(rv)
			rv = rv/normrv
			
			# now define rotation matrix
#			R = np.zeros((3,3))
			R[0,0] = cosphi + (rv[0]**2)*(1-cosphi)
			R[0,1] = rv[0]*rv[1]*(1-cosphi) - rv[2]*sinphi
			R[0,2] = rv[0]*rv[2]*(1-cosphi) + rv[1]*sinphi
	
			R[1,0] = rv[1]*rv[0]*(1-cosphi) + rv[2]*sinphi
			R[1,1] = cosphi + (rv[1]**2)*(1-cosphi)
			R[1,2] = rv[1]*rv[2]*(1-cosphi)-rv[0]*sinphi
	
			R[2,0] = rv[2]*rv[0]*(1-cosphi) - rv[1]*sinphi
			R[2,1] = rv[2]*rv[1]*(1-cosphi) + rv[0]*sinphi
			R[2,2] = cosphi + (rv[2]**2)*(1-cosphi)

		# Rotate normal
		self.n = np.dot(self.n,R)
		self.n = CleanMatElements(self.n)

		# Rotate plane vectors
		self.u = np.dot(self.u,R)
		self.v = np.dot(self.v,R)
		self.u = CleanMatElements(self.u)
		self.v = CleanMatElements(self.v)
		# Rotate positions of atoms
		posrot = np.dot(self.positions,R)
		posrot = CleanMatElements(posrot)

		# list of indexes of atom belonging to the plane
		idxlist = []
		# list of indexes of atom belonging to the plane and below
		idxlistblk = []

                idxlistBOOL = posrot[:,2] == 0
                idxlist = np.where(idxlistBOOL == 1)[0]
                idxlistblkBOOL = posrot[:,2] < 0
                idxlistblk = np.where(idxlistblkBOOL == 1)[0]

		# Check if plane exists:
		if len(idxlist) != 0:
			self.exists = True

		if self.exists:
			# Crate anticipatory vectors
			# Find all atoms lying above the plane
			idxabove = posrot[:,2]>0
			posabove = posrot[idxabove]
	
			# Labels of atoms on surface only
			self.planeatms = self.atoms[idxlist]

			uqatoms,uqidx,uqidxbulk = self.__finduqatomsCENTER(self.planeatms,\
					                   posrot,idxlist)
				
			# Find anticipatory vectors for each unique atom on the plane
			avecslist = []
			for idx in uqidxbulk:
				avecsitem = self.__anticipatoryvecs\
   					     (posrot[idx],posabove)
				avecslist.append(avecsitem)
	
			self.avecs = avecslist[0]
			# Create a dictionary avecall to hold all anticipatory
			# vectors assiciated to given atom type. This will be 
			# usefull for scoring function
			self.avecsall = {}
			for i in range(len(uqatoms)):
				lab = uqatoms[i]
				ii = avecslist[i]
				self.avecsall[lab]=ii


		#	self.avecs = self.__anticipatoryvecs(posrot[idxlist[0]],\
		#			                     posabove)
	
			# Find the nearest neighbours of the atom in bulk
	#		self.nneigh = self.__anticipatoryvecs(posrot[0],\
#		   			    posrot[1:],neighonly=True)
	
			# Find nearest unique atoms in the whole bulk 
			# This is needed to find nearest neigbours of each atom type
			uqatoms,uqidx,uqidxbulk = self.__finduqatomsCENTER(self.atomsSmall,\
					                   self.positionsSmall,\
							   idxlist=range(len(self.positionsSmall)))
	
			neighlist = []
                        for i in uqidxbulk:
                                # Construct array with poistion without atom i
                                postmpa = self.positionsSmall[:i]
                                postmpb = self.positionsSmall[i+1:]
                                postmp  = np.concatenate((postmpa,postmpb))
                                nneighat = self.__anticipatoryvecs(self.positionsSmall[i],\
                                                    postmp,neighonly=True)
                                neighlist.append(nneighat)

			# Create dictionary that as a key will have label 
			# of unique atom and as value, number of its nearest 
			# neighbors
			self.uqneigh = {}
			for i in range(len(uqatoms)):
				lab = uqatoms[i]
				ii = neighlist[i]
				self.uqneigh[lab]=ii
			self.nneigh = neighlist[0]
			
			# Revert to orginal origin of atom coordinates
			self.positions = self.positions + self.origin
	
			# Create surface coordinates
			self.planepos = posrot[idxlist]

			# Create surace coordiantes including atoms below it
			# Make sure that surface atoms are first in the array
			#idxlistblk = idxlist + idxlistblk
			idxlistblk = np.concatenate((idxlist,idxlistblk))
			self.planeposblk = posrot[idxlistblk]
	
			# Labels of atoms on surface and below
                        self.planeatmsblk = self.atoms[idxlistblk]
	
			# Translate plane so coordinate z=0	
#			if abs(self.planepos[:,2]).all > 0:
#				self.planepos[:,2] -= self.planepos[:,2]

		#end if self.exists

	def __finduqatoms(self,labels,pos,idxlist):
		# Routine to find unique types of atoms from set of
		# atom lables
		# Return atom labels, and index of the atom on the plane
		# and index of atoms in the bulk structure

		uql = [] # unique labels
		uqi = [] # index of representative atom on surface
		uqibulk = [] # index of representative atom in bulk structure

		postmp = pos[idxlist]	
		for i in range(len(postmp)):
			atom = labels[i]
			if atom not in uql:
				uql.append(atom)
				uqi.append(i)

		for i in uqi:
			uqibulk.append(idxlist[i])

		return uql,uqi,uqibulk

	def __finduqatomsCENTER(self,labels,pos,idxlist):
		# Routine to find unique types of atoms from set of
		# atom lables
		# Return atom labels, and index of the atom on the plane
		# and index of atoms in the bulk structureq
		#
		# Modified version of __finduqatoms to take atom in the center 
		# as the base to look

		uql = [] # unique labels
		uqi = [] # index of representative atom on surface
		uqibulk = [] # index of representative atom in bulk structure

		postmp = pos[idxlist]	

                # Find centroid 
                cX = sum(postmp[:,0])/len(postmp[:,0])
                cY = sum(postmp[:,1])/len(postmp[:,1])
                cZ = sum(postmp[:,2])/len(postmp[:,2])
                cXYZ = np.array((cX,cY,cZ))
                nnXYZ = fneigh(cXYZ,postmp)
                originIdx = int(nnXYZ[0][1])

		for i in range(originIdx,len(postmp)):
			atom = labels[i]
			if atom not in uql:
				uql.append(atom)
				uqi.append(i)

		for i in uqi:
			uqibulk.append(idxlist[i])

		return uql,uqi,uqibulk

	def __anticipatoryvecs(self,atom,bulk,neighonly=False):
		# Routine to find anticipatory vectors for atom
		# by looking through it nearest neighbours 

		# Find nearest neighbours of the atom
		# Construct array with the distance and indices of atoms in 
		# bulk
		# Introcude threshold to find number neartest neighbours
		# There are cases where atoms has N nearest neighbours, but
		# but they have only slighlty different distances, for instance
		# in the case of a-quarts-SiO2, Si has four NN, two has distance
		# 1.62A and two 1.63A. Use threshold to find them
		Nthresh = 0.1

		neigh = np.zeros((len(bulk),2))
		idx = 0
		for atm in bulk:
			d = np.linalg.norm(atom-atm)
			neigh[idx,0]=d
			neigh[idx,1]=idx
			idx += 1

		# Sort the distances according to the distance
		neigh = neigh[np.argsort(neigh[:,0])]

		# To avoid small differences in floating point comparisons
		# substract first elements from all others and round the 
		# result to 12 decimal place
		shift = neigh[0,0].copy()
#		neigh[:,0] -= neigh[0,0] 
		neigh[:,0] -= shift

		neigh = CleanMatElements(neigh)
		
		# Select only atoms with the shortest distance
		#idx = neigh[:,0] == neigh[0,0]
		#idx = neigh[:,0] == 0.0
		idx = neigh[:,0] <= Nthresh
		neigh = neigh[idx]

		if neighonly:
			return len(neigh)

		# Find the coordinates of nearest neighbours in bulk
		idxs = neigh[:,1].astype(int)
		avecs = bulk[idxs]
		
		# Shift avecs to origin of coordinate system
		avecs = avecs - atom
		avecs = CleanMatElements(avecs)

		# If there is more than one anitcipatory vector,
		# see if they are quivalent. If yes, remove the duplicates
		# The vectors are  assumed to be equivalent if the angle
		# between vectors and z-axis is the same
		
#		if len(avecs) > 1:
#			z = np.array((0,0,1))
#			nz = 1
#			ang = []
#			idx = []
#			ii = 0
#			for v in avecs:
#				nv = np.linalg.norm(v)
#				cosphi = np.dot(v,z)/(nv*nz)
#				if cosphi not in ang:
#					ang.append(cosphi)
#					idx.append(ii)
#				ii += 1
#			avecs = avecs[idx]

		return avecs

	def initpvec(self):
		# Initialize primitive vecotrs
		# Algorithm:
		# Take first atom in the surface, calculate its nearest 
		# neighbours, and define primitive vectors as the 
		# vectors pointing to the nearest two atoms. 

		# Shift coordinates so that 1st atom is in origin
		orgsave = self.planepos[0]
		self.planepos -= self.planepos[0]

		# Find all the nearest neighbours from atom 0
		distmat = np.sum(np.abs(self.planepos)**2,axis=1)**(1./2)
		# Sort distmat. 
		idx = distmat.argsort()

		# Find two non-linear vectors. They will be initial 
		# primitive vectors
		# It is most likely to be fist two smallest vectors,
		# but check for linearity just in case

		linear = True
		i = 1
		while i < len(distmat):
			j = i + 1
			while j < len(distmat):
				self.a = self.planepos[idx[i]]
				self.b = self.planepos[idx[j]]

				# In the case of crystals containig different
				# atoms, the primitive vecors need to be defined
				# between atoms of the same type.
				if self.planeatms[idx[i]] == self.planeatms[0]\
				and self.planeatms[idx[j]] == self.planeatms[0]:

					check = np.cross(self.a,self.b)
					# Clean numeric noise
					check = CleanMatElements(check)
					#TMP
					self.norma=np.linalg.norm(self.a)
					self.normb=np.linalg.norm(self.b)
					print "NORM PRIM ",self.norma,self.normb
					cosphi = np.dot(self.a,self.b)/(self.norma*self.normb)
					if np.linalg.norm(check) != 0:
						linear = False

				if not linear: break
				j += 1
			i += 1
			if not linear: break

                if linear:
                        print "Only linear primitive vectors found"
                        print "Something is wrong. Exiting."
                        exit()

		# Reduce a and b
		self.a, self.b = reduction(self.a, self.b)

		print "REDUCED VECTORS"
		print "A", self.a, np.linalg.norm(self.a)
		print "B", self.b, np.linalg.norm(self.b)
		
		# Shift coordinates back to orginals positions 
		self.planepos += orgsave

	def initpvecNEW(self):
		# Initialize primitive vecotrs
		# Algorithm:
		# Take first atom in the surface, calculate its nearest 
		# neighbours, and define primitive vectors as the 
		# vectors pointing to the nearest two atoms. 

		# Algortihm used:
		# Check two conditions:
		# 1) If lattice constant for given plane can be expresed 
		#    in term of intiger multiplications of primitive vector:
		#        n*|a| = |u| , n=1,2,3,....
		# 2) If length of the scalar projection of the primitive 
		#    vector on the lattice vector is equal to 0.5*(|u|**2)
		#    This can be shown from the properties of dot product:
		#   dot(a,u) == 0.5*(|u|**2) when a_u = |a|*cos(phi) == 0.5*|u|
		#
		#   If any of those conditions is met, the pair of vectors 
		#   are the primitive vectors for this lattice

		# Primitive vectors not found yet
		self.exists = False

		# Shift coordinates so that 1st atom is in origin
		orgsave = self.planepos[0].copy()
		self.planepos -= orgsave

		# Find all the nearest neighbours from atom 0
		distmat = np.sum(np.abs(self.planepos)**2,axis=1)**(1./2)
		# Sort distmat. 
		idx = distmat.argsort()


		# Find two non-linear vectors. They will be initial 
		# primitive vectors
		# It is most likely to be fist two smallest vectors,
		# but check for linearity just in case

		nu = np.linalg.norm(self.u)
		nv = np.linalg.norm(self.v)

		primitive = False

		# Search for 1st primitive vector
		i = 1
		while i < len(distmat):
			# In the case of crystals containig different
			# atoms, the primitive vecors need to be defined
			# between atoms of the same type.
			if self.planeatms[idx[i]] == self.planeatms[0]:
				self.a = self.planepos[idx[i]].copy()
	
				na = np.linalg.norm(self.a)
	
				# 1st condtition
				if nu >= na:
					ma = nu%na
				else:
					ma = na%nu
	
				if round(ma,6) == 0.0: primitive = True
	
				# 2nd condition
				dau = np.dot(self.a,self.u)
	
				if round(abs(dau),5) == round((nu**2)/2,5): primitive = True

				if primitive: break
			i += 1

		# Search for 2nd primitive vector
		linear = True
		primitive = False
		i = 1
		while i < len(distmat):
			# In the case of crystals containig different
			# atoms, the primitive vecors need to be defined
			# between atoms of the same type.
			if self.planeatms[idx[i]] == self.planeatms[0]:
				self.b = self.planepos[idx[i]].copy()
	
				nb = np.linalg.norm(self.b)
	
				# 1st condtition
				if nv >= nb:
					mb = nv%nb
				else:
					mb = nb%nv
	
				if round(mb,6) == 0.0 : primitive = True
	
				# 2nd condition
				dbv = np.dot(self.b,self.v)
				
				if round(abs(dbv),5) == round((nv**2)/2,5): primitive = True
	
				# Check if vector b is not linear with vector a
				if primitive:
					check = np.cross(self.a,self.b)
					# Clean numeric noise
					check = CleanMatElements(check)
					self.norma=np.linalg.norm(self.a)
					self.normb=np.linalg.norm(self.b)
					cosphi = np.dot(self.a,self.b)/\
						       (self.norma*self.normb)
					if np.linalg.norm(check) != 0:
						linear = False
	
				if primitive and (not linear): break
			i += 1

		if primitive:
			self.exists = True
			# Reduce a and b
			self.a, self.b = reduction(self.a, self.b)
#		else:
#			print "COULDN'T FIND PRIMITIVE VECTORS"
#			exit()
#
#                print "REDUCED VECTORS"
#                print "A", self.a, np.linalg.norm(self.a)
#                print "B", self.b, np.linalg.norm(self.b)
		
		# Shift coordinates back to orginals positions 
		self.planepos += orgsave

	def initpvecNEW2(self):
		# Initialize primitive vecotrs
		# Algorithm:
		# Take first atom in the surface, calculate its nearest 
		# neighbours, and define primitive vectors as the 
		# vectors pointing to the nearest two atoms. 

		# Algortihm used:
		# Check two conditions:
		# 1) If lattice constant for given plane can be expresed 
		#    in term of intiger multiplications of primitive vector:
		#        n*|a| = |u| , n=1,2,3,....
		# 2) If length of the scalar projection of the primitive 
		#    vector on the lattice vector is equal to 0.5*(|u|**2)
		#    This can be shown from the properties of dot product:
		#   dot(a,u) == 0.5*(|u|**2) when a_u = |a|*cos(phi) == 0.5*|u|
		#
		#   If any of those conditions is met, the pair of vectors 
		#   are the primitive vectors for this lattice

		# Primitive vectors not found yet
		self.exists = False

		# Shift coordinates so that 1st atom is in origin
		orgsave = self.planepos[0].copy()
		self.planepos -= orgsave

		# Find all the nearest neighbours from atom 0
		distmat = np.sum(np.abs(self.planepos)**2,axis=1)**(1./2)
		# Sort distmat. 
		idx = distmat.argsort()

		# Find two non-linear vectors. They will be initial 
		# primitive vectors
		# It is most likely to be fist two smallest vectors,
		# but check for linearity just in case

		nu = np.linalg.norm(self.u)
		nv = np.linalg.norm(self.v)

		primitive = False

		# Search for 1st primitive vector
		i = 1
		while i < len(distmat):
			# In the case of crystals containig different
			# atoms, the primitive vecors need to be defined
			# between atoms of the same type.
			if self.planeatms[idx[i]] == self.planeatms[0]:
				self.a = self.planepos[idx[i]].copy()
	
				na = np.linalg.norm(self.a)
	
				# 1st condtition
				if nu >= na:
					ma = nu%na
				else:
					ma = na%nu
	
				if round(ma,6) == 0.0: primitive = True
	
		#		# 2nd condition
		#		dau = np.dot(self.a,self.u)
		#		print "2nd condition:",round(abs(dau),5),round((nu**2)/2,5)
		#		if round(abs(dau),5) == round((nu**2)/2,5): primitive = True

				# 3rd condition u and a are colinear
				check = np.cross(self.a,self.u)
				# Clean numeric noise
				check = CleanMatElements(check)
				cosphi = np.dot(self.a,self.u)/\
					       (na*nu)
				if np.linalg.norm(check) == 0:
					primitive = True


				if primitive: break
			i += 1
		# Search for 2nd primitive vector
		linear = True
		primitive = False
		i = 1
		while i < len(distmat):
			# In the case of crystals containig different
			# atoms, the primitive vecors need to be defined
			# between atoms of the same type.
			if self.planeatms[idx[i]] == self.planeatms[0]:
				self.b = self.planepos[idx[i]].copy()
	
				nb = np.linalg.norm(self.b)
	
				# 1st condtition
				if nv >= nb:
					mb = nv%nb
				else:
					mb = nb%nv
	
				if round(mb,6) == 0.0 : primitive = True
	
				# 2nd condition
		#		dbv = np.dot(self.b,self.v)
		#		
		#		if round(abs(dbv),5) == round((nv**2)/2,5): primitive = True
		#
				# 3rd condition u and a are colinear
				check = np.cross(self.b,self.v)
				# Clean numeric noise
				check = CleanMatElements(check)
				cosphi = np.dot(self.b,self.v)/\
					       (nb*nv)
				if np.linalg.norm(check) == 0:
					primitive = True

				# Check if vector b is not linear with vector a
				if primitive:
					check = np.cross(self.a,self.b)
					# Clean numeric noise
					check = CleanMatElements(check)
					self.norma=np.linalg.norm(self.a)
					self.normb=np.linalg.norm(self.b)
					cosphi = np.dot(self.a,self.b)/\
						       (self.norma*self.normb)
					if np.linalg.norm(check) != 0:
						linear = False
	
				if primitive and (not linear): break

			i += 1
		
		if primitive:
			self.exists = True

			# Reduce a and b
			self.a, self.b = reduction(self.a, self.b)

		# Shift coordinates back to orginals positions 
		self.planepos += orgsave

	def primitivecell(self):
		# Calculate the norm of primitive vectors a,b , 
		# angle between them and the area of primitive cell

		self.norma=np.linalg.norm(self.a)
		self.normb=np.linalg.norm(self.b)

		cosphi = np.dot(self.a,self.b)/(self.norma*self.normb)
		phi = m.acos(cosphi)
		self.phiab = phi * 180/m.pi
		sinphi = m.sin(phi)

		self.primarea = self.norma * self.normb * sinphi

def normal(vec1,vec2):
	# Calculate normal to the plane given by vec1 and vec2
		result = np.cross(vec1,vec2)
		return result

def rotmat2d(phi):
	# Calculate clockwise rotation matrix around z-axis 
	# 
	# phi - input angle, in radians

	mat = np.zeros((2,2))

	mat[0,0] =  m.cos(phi)
	mat[0,1] = -1*m.sin(phi)
	mat[1,0] =  m.sin(phi)
	mat[1,1] =  m.cos(phi)

	return mat


#CONTINUED FRACTIONS CODE		
# source:
# http://tech.groups.yahoo.com/group/tuning-math/message/14958
def contfrac(a): ### continued fraction expansion
	terms=[]
	count=0
	b=1
	while ((b != 0) and (count < 24)): ### 24 times, emprical accuracy
		### limit
		terms.append(m.floor(a/(b+0.0)))
		a,b = b, a % b
		count = count + 1
	return terms

def ra(x): ### 'rational approximation', or convergent
	numerators=[0,1]
	denominators=[1,0]
	expansion=contfrac(x) ### call the contfrac function
	for num in expansion: ### [-1] and [-2] index 'previous'
		### and 'previous-previous'
		num = int(num)
		numerators.append((num*numerators[-1])+numerators[-2])
		denominators.append((num*denominators[-1])+denominators[-2])
#	for index in range(len(numerators)):
#		print "%i  %i" % (numerators[index], denominators[index])
	return numerators[2:], denominators[2:]
# END OF CONTINUED FRACTIONS CODE

def FindAllDivisors(x):
	# Algorithm taken from
	# http://stackoverflow.com/questions/12421969/finding-all-divisors-of-a-number-optimization
	divList = []
	y = 1
	while y <= m.sqrt(x):
		if x % y == 0:
			divList.append(y)
			# If x is a square of y, don't add it to div list to
			# avoid doubling the same divisor
			if (int(x/y)) != y:
				divList.append(int(x / y))
		y += 1
	divList.sort()
	return divList


def AreaMatch(A1, A2, threshold, maxarea):

	# A1 - area of primitive cell of material1
	# A2 - area of primitive cell of material2
	# threshold - threshold for are misfit (in %)
	# maxarea - maximum size of p*A1 or q*A2 (in Angstrom)

	# Find p and q, such that:
	# p*A1 = q*A2 => p/q = A2/A1

	ratio = A2/A1

	# Use continued fraction algorithm
	p,q = ra(ratio)

	# Prepare output lists	
	pout = []
	qout = []

	# Filter the results
	# Convert threshold value from % to decimal
	threshold = threshold /100.0
	for i in range(len(p)):

		area1 = p[i] * A1
		area2 = q[i] * A2

		ratio = area1/area2

		misfit = 1-ratio

		if abs(misfit) < threshold and area1 < maxarea \
				 and area2 < maxarea:
			pout.append(p[i])
			qout.append(q[i])

	return pout, qout

def AreaMatch2(area1,area2,maxarea,threshold):

	# area1 - area of primitive cell of material1
	# area2 - area of primitive cell of material2
	# threshold - threshold for are misfit (in %)
	# maxarea - maximum size of p*A1 or q*A2 (in Angstrom)
	# 
	# Generate all the p and q, where
	#
	#   pA1 = qA2   (A1 - area1, A2 - area2)
	#
	# that give pA1 and qA2 that  is smaller than maxarea. 
	# Return only those  p and q that pA1/qA <= threshold
		
	p = []
	q = []
	
	i = 1 
	threshold = threshold/100.0
	while i*area1 <maxarea:
		j = 1
		while j*area2 < maxarea:
			r1 = i*area1
			r2 = j*area2
			pq = (1-(r1/r2))
			if abs(pq) <= threshold:
				p.append(i)
				q.append(j)
		#		print i,j,abs(pq)
			j += 1
		i += 1
		
	return p,q

def reduction(a,b):
	# Primitive vectors reduction algortihm
	# Refeference:
	# Zur, McGil, J. Appl. Phys. 55, 378 (1984)

	reduced = False

	while not reduced:

		dot = np.dot(a,b)

		if dot < 0:
			b = -1 * b

		norma = np.linalg.norm(a)
		normb = np.linalg.norm(b)

		if norma <= normb:
			ab = a+b
			normab = np.linalg.norm(ab)

			if normb <= normab:
				amb = a - b
				normamb = np.linalg.norm(amb)

				if normb <=normamb:
					reduced = True
				else:
					b = b - a
			else:
				b = b + a
		else: 
			tmp = a
			a = b
			b = tmp
	
	return a,b



def Superlattice(a,b,n):

	# Find all the possible superlattices formed from primitive cell
	# (a,b) multiplied by n
	# a - primitive cell vector a
	# b - primitive cell vector b
	# n - primitive cell multiplicator

	# Based on  Zur, McGil, J. Appl. Phys. 55, 378 (1984) eqs. (2.3)-(2.6)

	# Find the divisors on n

	divlist = FindAllDivisors(n)

	# Construct transformation (2.3) from Zur-McGil
	results = []
	for m in divlist:
		i = n/m # Eq. (2.4)
		for j in range(m):
			tmat = np.zeros((2,2))
			tmat[0,0] = i
			tmat[0,1] = j
			tmat[1,1] = m

			vec = np.array((a[0:2],b[0:2]))

			result = np.dot(tmat,vec)

			# Reduce pair of superlattice vectors
			reda, redb = reduction(result[0], result[1])

			# store the result as vector in 3D (with z = 0)
			tmp1 = np.array((0.0,0.0,0.0))
			tmp2 = np.array((0.0,0.0,0.0))

			tmp1[0:2] = tmp1[0:2] + reda
			tmp2[0:2] = tmp2[0:2] + redb

			pair = np.array((tmp1,tmp2))

			results.append(pair)
	
	return results

def SuperlatticeMatch(vecset1, vecset2, vecthresh, anglethresh):
	# Check if two supperlattices matches by comapring their 
	# vectors and angle between them

	# vecset1 - [2,3] array with set of vectors of lattice1
	# vecset2 - [2,3] array with set of vectors of lattice2
	# vecthresh - threshold for vector difference (in %)
	# anglethresh - threshold for angle diffrence (in %)

	fit = False

	# Convert thresholds to decimals
	vecthresh = vecthresh / 100.0
	anglethresh = anglethresh / 100.0

	u1 = vecset1[0]
	v1 = vecset1[1]

	u2 = vecset2[0]
	v2 = vecset2[1]

	normu1 = np.linalg.norm(u1)
	normv1 = np.linalg.norm(v1)

	normu2 = np.linalg.norm(u2)
	normv2 = np.linalg.norm(v2)

	cosphi = np.dot(u1,v1)/(normu1*normv1)
	sinphi = m.sqrt(1-cosphi**2)
	phiuv1 = m.acos(cosphi)*180/m.pi

	cosphi = np.dot(u2,v2)/(normu2*normv2)
	sinphi = m.sqrt(1-cosphi**2)
	phiuv2 = m.acos(cosphi)*180/m.pi

	umisfit = abs(1 - normu1/normu2)
	vmisfit = abs(1 - normv1/normv2)

	anglemisfit = abs(1 - phiuv1/phiuv2)

	if umisfit <  vecthresh and  vmisfit < vecthresh \
			and anglemisfit <  anglethresh:
		fit = True
	
	return fit, normu1, normv1, phiuv1, normu2, normv2, phiuv2


def SuperlatticeParams(vecset):

	u = vecset[0]
	v = vecset[1]


	normu = np.linalg.norm(u)
	normv = np.linalg.norm(v)


	cosphi = np.dot(u,v)/(normu*normv)
	sinphi = m.sqrt(1-cosphi**2)
	phiuv = m.acos(cosphi)*180/m.pi
	area = normu*normv*sinphi
	
	return normu, normv, phiuv, area

def CheckListDuplicates(inlist, elem, in1, in2, in3, in4, in5, in6):

	# Check if floating point elements of list elem with indicies 
	# in1, in2, in3 are already in the in the list of listst inlist 

	duplicate = False
	for i in inlist:
		e1 = abs(i[in1]-elem[in1])
		e2 = abs(i[in2]-elem[in2])
		e3 = abs(i[in3]-elem[in3])
		e4 = abs(i[in4]-elem[in4])
		e5 = abs(i[in5]-elem[in5])
		e6 = abs(i[in6]-elem[in6])
		if e1 < 1e-9 and e2 <1e-9 and e3 < 1e-9 \
		and e4 < 1e-9 and e5 <1e-9 and e6 < 1e-9:
			duplicate = True
	return duplicate


class Interface:
	def __init__(self,vecDep,vecSub,Deposit,Substrate,vecpair,bfrac,wbond,\
	    	     atomicRadius, nLS, nLD, capAtD, capAtS, sandwich, genHD=True,genHS=True,\
		     poissonRatio=True):

		# Deposit surface and atoms below it
		posDepBlk = Deposit.planeposblk.copy()
		# Substrate surface and atoms below it
		posSubBlk = Substrate.planeposblk.copy() 
		# Deposit atom labels
		DepAtomsBlk = Deposit.planeatmsblk
		# Substrate atom labels
		SubAtomsBlk = Substrate.planeatmsblk
		# Deposit and Substrate atom types
		self.atomTyp = Deposit.atomTyp
		# Number of nearest neigbors of the atom in the bulk of Depoit
		DepNneigh = Deposit.nneigh
		#Number of nearest neigbors of the atom in the bulk of Substrate
		SubNneigh = Substrate.nneigh	

		DepNneighuq = Deposit.uqneigh
		SubNneighuq = Substrate.uqneigh

		self.vecDep = vecDep.copy() # vectors of defining deposit surface
		self.vecSub = vecSub.copy() # vectors defining substrate surface
		self.Depavecs = Deposit.avecs.copy() # Deposit anticipatory vectors
		self.Subavecs = Substrate.avecs.copy() # Substrate anticipatory vectors	
		# Dictionary with vectors for all the atom types
		self.Depavecsall = Deposit.avecsall.copy()
		self.Subavecsall = Substrate.avecsall.copy()
		# wbond - parameter for scoring function


		# Prepare xyz output for calculations

		periodic = True
		#periodic = False
		nlayersS = nLS# 0 for scoring func - gives two same surfaces for
		nlayersD = nLD
#					SiO2
#		nlayers = 14# 0 for scoring func

		# Create Periodic Deposit and Substrate superlattices

		DSurfPos,DSurfAtm,DepavecsR,vecDepR,DH =\
						 self.__CreateSurface\
				                (posDepBlk,DepAtomsBlk,\
				                 self.vecDep[0],self.vecDep[1],\
						 self.Depavecs,nlayersD,\
						 DepNneighuq,periodic,\
						 genHD,genHS)

#		ii=0
#		for i in DSurfPos:
#			print DSurfAtm[ii],i[0],i[1],i[2]
#			ii+=1
#		print "SubNneighuq",SubNneighuq
		SSurfPos,SSurfAtm,SubavecsR,vecSubR,SH =\
		 	 			self.__CreateSurface\
				                (posSubBlk,SubAtomsBlk,\
				                 self.vecSub[0],self.vecSub[1],\
						 self.Subavecs,nlayersS,\
						 SubNneighuq,periodic,\
						 genHD,genHS)

		# Scale deposit so it matches substrate exactely
		DSurfPosScale = self.__scale(DSurfPos,\
   				             vecDepR[0],vecDepR[1],\
					     vecSubR[0],vecSubR[1],poissonRatio)
		self.IfacePosSC,self.IfaceAtmSC,self.idxDep,self.idxSub,alignvec \
					     = self.__alignsurfaces\
				              (DepavecsR,SubavecsR,\
					      DSurfPosScale,SSurfPos,\
					      DSurfAtm,SSurfAtm,bfrac,\
					      vecpair,vecSubR)
		# Make vectors globaly available:
		self.IfaceVecs = [vecSubR[0],vecSubR[1],alignvec]

		# If vectors are passing through rectangular interface, 
		# translate those atoms so PBC conditions can be met.

	   	# Define new indieces: idxDepH and idxSubH, which contain
		# the number list of indieces of Depost and Hydrogens and 
		# Substrate and Hydrogens
		# idxSub and idxDep contain indices of only Deposit and only
		# Substrate
		if genHD or genHS:
			self.IfacePosSC,self.IfaceAtmSC,self.idxDep,\
					self.idxSub,self.idxDepH,self.idxSubH\
					= self.__addHydrogens(self.IfacePosSC,\
					self.IfaceAtmSC,self.idxDep,\
					self.idxSub,DH,SH,genHD,genHS)
		else:
			self.idxDepH = self.idxDep
			self.idxSubH = self.idxSub


		if capAtD or capAtS:
			if sandwich: capAtD = False
			self.IfaceAtmSC = self.__addCapAtoms(self.IfacePosSC,\
					self.IfaceAtmSC, capAtD, capAtS) 

		botSub = min(self.IfacePosSC[:,2])
		topDep = max(self.IfacePosSC[:,2])
       		topSub = max(self.IfacePosSC[self.idxSubH][:,2]) # Substrate layer in the interface
       		botDep = min(self.IfacePosSC[self.idxDepH][:,2]) # Deposit layer in the interface
	        self.SDdist = botDep - topSub

		if sandwich:
			moveDist = topDep - botSub + self.SDdist
			newSub = self.IfacePosSC[self.idxSubH].copy()
			newSubLabels = self.IfaceAtmSC[self.idxSubH].copy()
			newSub[:,2] -= topSub
			newSub[:,2] *= -1
			newSub[:,2] += moveDist
			newIdxSub = np.arange(1,len(newSub)+1)
			newIdxSub += self.idxSubH[-1]
			self.idxSubH = np.concatenate((self.idxSubH,newIdxSub))
			self.IfacePosSC = np.concatenate((self.IfacePosSC, newSub))
			self.IfaceAtmSC = np.concatenate((self.IfaceAtmSC, newSubLabels))
			self.IfaceVecs[2][-1] += topSub - botSub + self.SDdist

		self.IfacePosSC = self.__checkedge(self.IfacePosSC,self.IfaceVecs)
		self.IfacePosSC = self.__checkedge(self.IfacePosSC,self.IfaceVecs)
		self.IfacePosSC = self.__checkedge(self.IfacePosSC,self.IfaceVecs)
#		# 3 was orginal
		self.IfacePosSC = self.__checkedge(self.IfacePosSC,self.IfaceVecs)
		self.IfacePosSC = self.__checkedge(self.IfacePosSC,self.IfaceVecs)
		# Some more
		self.IfacePosSC = self.__checkedge(self.IfacePosSC,self.IfaceVecs)
		self.IfacePosSC = self.__checkedge(self.IfacePosSC,self.IfaceVecs)
		self.IfacePosSC = self.__checkedge(self.IfacePosSC,self.IfaceVecs)

		### BEGIN FOR SCORING FUNC
		#print "STARTING GENERATING STRUCTURES FOR SCORING FUNCTION"
		print "Generating structures for scoring function"
		nlayersS = 0# 0 for scoring func
		nlayersD = 0# 0 for scoring func
		nlayersS2 = 5
		nlayersD2 = 5

		DSurfPos,DSurfAtm,DepavecsR,vecDepR,DH =\
						 self.__CreateSurface\
				                (posDepBlk,DepAtomsBlk,\
				                 self.vecDep[0],self.vecDep[1],\
						 self.Depavecs,nlayersD,\
						 DepNneighuq,periodic,\
						 genHD,genHS)

		# Deposit with few layers on top
		DSurfPosL,DSurfAtmL,DepavecsR,vecDepR,DH =\
						 self.__CreateSurface\
				                (posDepBlk,DepAtomsBlk,\
				                 self.vecDep[0],self.vecDep[1],\
						 self.Depavecs,nlayersD2,\
						 DepNneighuq,periodic,\
						 genHD,genHS)



#		ii=0
#		for i in DSurfPos:
#			print DSurfAtm[ii],i[0],i[1],i[2]
#			ii+=1
#		exit()
		SSurfPos,SSurfAtm,SubavecsR,vecSubR,SH =\
		 	 			self.__CreateSurface\
				                (posSubBlk,SubAtomsBlk,\
				                 self.vecSub[0],self.vecSub[1],\
						 self.Subavecs,nlayersS,\
						 SubNneighuq,periodic,\
						 genHD,genHS)

		# Substrate with few layers on top
		SSurfPosL,SSurfAtmL,SubavecsR,vecSubR,SH =\
		 	 			self.__CreateSurface\
				                (posSubBlk,SubAtomsBlk,\
				                 self.vecSub[0],self.vecSub[1],\
						 self.Subavecs,nlayersS2,\
						 SubNneighuq,periodic,\
						 genHD,genHS)


#		ii=0
#		for i in SSurfPos:
#			print SSurfAtm[ii],i[0],i[1],i[2]
#			ii+=1
#		exit()

		# Scale deposit so it matches substrate exactely
		DSurfPosScale = self.__scale(DSurfPos,\
   				             vecDepR[0],vecDepR[1],\
					     vecSubR[0],vecSubR[1],poissonRatio)

		DSurfPosScaleL = self.__scale(DSurfPosL,\
   				             vecDepR[0],vecDepR[1],\
					     vecSubR[0],vecSubR[1],poissonRatio)

		# extend scaled deposit, use Sub vecs as they are the same
#		DSurfPosScaleEx,DSurfAtmExSC = self.__extendsurface\
#				                  (DSurfPosScale,\
#				DSurfAtm,vecSubR[0],vecSubR[1],periodic)

		DSurfPosScaleExL,DSurfAtmExSCL = self.__extendsurface\
				                  (DSurfPosScaleL,\
				DSurfAtmL,vecSubR[0],vecSubR[1],periodic)
		# Same for substrate
#		SSurfPosEx,SSurfAtmEx = self.__extendsurface\
#				                  (SSurfPos,\
#				SSurfAtm,vecSubR[0],vecSubR[1],periodic)

		SSurfPosExL,SSurfAtmExL = self.__extendsurface\
				                  (SSurfPosL,\
				SSurfAtmL,vecSubR[0],vecSubR[1],periodic)

		# substrate on extended deposit
#		self.IfacePosExSC,self.IfaceAtmExSC,self.idxDepExD,\
#						     self.idxSubExD,alignvec \
#					     = self.__alignsurfaces\
#				              (DepavecsR,SubavecsR,\
#					      DSurfPosScaleEx,SSurfPos,\
#					      DSurfAtmExSC,SSurfAtm,bfrac,\
#					      vecpair,vecSubR)

		self.IfacePosExSC,self.IfaceAtmExSC,self.idxDepExD,\
						     self.idxSubExD,alignvec \
					     = self.__alignsurfaces\
				              (DepavecsR,SubavecsR,\
					      DSurfPosScaleExL,SSurfPos,\
					      DSurfAtmExSCL,SSurfAtm,bfrac,\
					      vecpair,vecSubR)

		# Scaled Deposit on extented substrate
#		self.IfacePosSExSC,self.IfaceAtmSExSC,self.idxDepExS,\
#						     self.idxSubExS,alignvec \
#					     = self.__alignsurfaces\
#				              (DepavecsR,SubavecsR,\
#					      DSurfPosScale,SSurfPosEx,\
#					      DSurfAtm,SSurfAtmEx,bfrac,\
#					      vecpair,vecSubR)

		self.IfacePosSExSC,self.IfaceAtmSExSC,self.idxDepExS,\
						     self.idxSubExS,alignvec \
					     = self.__alignsurfaces\
				              (DepavecsR,SubavecsR,\
					      DSurfPosScale,SSurfPosExL,\
					      DSurfAtm,SSurfAtmExL,bfrac,\
					      vecpair,vecSubR)

#		ii=0
#		for i in self.IfacePosExSC:
#			print self.IfaceAtmExSC[ii],i[0],i[1],i[2]
#			ii+=1
		### END FOR SCORING FUNC 

#		print "SC1 DEPOSIT"
		SubVecX = np.linalg.norm(vecSubR[0])
		SubVecY = np.linalg.norm(vecSubR[1])
		DepVecX = np.linalg.norm(vecDepR[0])
		DepVecY = np.linalg.norm(vecDepR[1])

	        cosphi = np.dot(vecSubR[0],vecSubR[1])/(SubVecX*SubVecY)
	        phiSub = m.acos(cosphi)*180/m.pi

         	cosphi = np.dot(vecDepR[0],vecDepR[1])/(DepVecX*DepVecY)
       	        phiDep = m.acos(cosphi)*180/m.pi

		misfitX = abs(DepVecX/SubVecX)
		misfitY = abs(DepVecY/SubVecY)
		misfitPhi = (phiDep/phiSub)

		if misfitX >1:
			misfitX = 1-abs(misfitX-1)

		if misfitY >1:
			misfitY = 1-abs(1-misfitY)

		if misfitPhi >1:
			misfitPhi = 1-abs(misfitPhi-1)

		self.sc1 = self.SCORE7(self.IfacePosSExSC,self.idxDepExS,\
				       self.idxSubExS,self.Depavecs,\
				       self.Subavecs,self.IfaceVecs,\
				       self.IfaceAtmSExSC,self.Depavecsall,\
				       wbond,atomicRadius,misfitX,misfitY,\
				       misfitPhi)
#		self.sc1 = 0
#		print "SCORE OUTPUT SC1" ,self.sc1
#		print "SC2 SUBSTRATE"
		self.sc2 = self.SCORE7(self.IfacePosExSC,self.idxSubExD,\
				       self.idxDepExD,self.Subavecs,\
				       self.Depavecs,self.IfaceVecs,\
				       self.IfaceAtmExSC,self.Subavecsall,\
				       wbond,atomicRadius,misfitX,misfitY,\
				       misfitPhi)
#		self.sc2 = 0
#		print "SCORE OUTPUT SC2", self.sc2

	def __CreateSurface(self,posin,atoms,vec1,vec2,avecsin,nlayers,\
			    bulkNN,periodic=True,genHD=True,genHS=True):
		# Create surface given reduced superlattice vectors 
		# and plane coordinates coordinates
		# The output surface is rotated so that 
		# vec1 point x direction
		# Return
		# - positions and labels of surface
		# - rotated anticipatory vectors
		# - rotated vectors
		
		# list of indexes of atom belonging to the plane
		#idxlist = []
		#planeatms = []

		pos = posin.copy()
		avecs = avecsin.copy()

		# Limit only to nlayers of atoms
		iidx = abs(pos[:,2]) <= nlayers
		pos = pos[iidx]

		# Generate atom labels for limited positions
#		tmplabels = []
#		for i in range(len(iidx)):
#			if iidx[i]:
#				tmplabels.append(atoms[i])
		tmplabels = atoms[iidx]

		#TODO: this is not used anymore since we have addCapAtom
		if genHD or genHS:
			Hpos = self.__genHydrogen(pos,tmplabels,bulkNN)
		else:
			Hpos = 0
		idx = 0

                # Orient the atoms in such a way that the coorindate origin is 
                # in the middle of the top plane
                # Find positions of only top plane
                idxPlane = abs(pos[:,2]) == 0.0
                posPlane = pos[idxPlane]
                # Find centroid of this plane
                cX = sum(posPlane[:,0])/len(posPlane[:,0])
                cY = sum(posPlane[:,1])/len(posPlane[:,1])
                # The z-variable of centroid is constant
                cZ = posPlane[0,2]
                cXYZ = np.array((cX,cY,cZ))
                nnXYZ = fneigh(cXYZ,posPlane)
                originIdx = int(nnXYZ[0][1])

		# Find the all the atoms that are inside the area 
		# designated by vectors vec1 and vec2
		# We will use barycentric coordinates to do this
		# by defining two triangles with following corners:
		# 1) (0,0,0), vec1, vec2
		# 2) vec1+vec2, vec1, vec2
	
		# Shift the coordiante system so the origin is on first atom
		# We need this for calculations in barycentric coordinates

#		pos -= pos[0] # Does not work for large number of atoms 
		#shift = pos[0].copy()
		shift = posPlane[originIdx].copy()
		pos -= shift

		idxlist = np.zeros(len(pos),bool)

		for atom in pos:

			# Lower triangle
			# p1 = [0,0]
			# p2 = vec1
			# p3 = vec2
			p1 = np.array((0.0,0.0))

			alpha1,beta1,gamma1 = self.__barycentric(atom,p1,\
						              vec1,vec2)

			a1 = alpha1 >=0 and alpha1 <= 1
			b1 = beta1  >=0 and beta1  <= 1
			c1 = gamma1 >=0 and gamma1 <= 1

			# Upper triangle
			# p1 = vec1+vec2
			# p2 = vec1
			# p3 = vec2

			p1 = vec1+vec2

			alpha2,beta2,gamma2 = self.__barycentric(atom,p1,\
						              vec1,vec2)


			a2 = alpha2 >= 0 and alpha2 <= 1
			b2 = beta2  >= 0 and beta2  <= 1
			c2 = gamma2 >= 0 and gamma2 <= 1

			if (a1 and b1 and c1) or (a2 and b2 and c2):
			#	idxlist.append(idx)
				idxlist[idx] = True
			
			idx += 1
		# Create the surface
		planepos = pos[idxlist]

		# and corresponind atom labels
#		for i in idxlist:
#			planeatms.append(tmplabels[i])
		planeatms = tmplabels[idxlist]

#		ii=0
#		for i in planepos:
#			print planeatms[ii],i[0],i[1],i[2]
#			ii+=1

		# Rotate the plane in such a way that vec1 is aligned 
		# with x axis
		# Find the angle
		vecx = np.array((1.0,0.0,0.0))
		normvec1 = np.linalg.norm(vec1)
		normx = np.linalg.norm(vecx)
		cosphi = np.dot(vec1,vecx)/(normvec1*normx)
		phi = m.acos(cosphi)
		
		# If vec1 is in 1st and 2nd quater of coordinate system
		# define the angle such that  rotation is clockwise to x axis
		if vec1[1] > 0:
	        	phi =  phi * 180/m.pi
		       	phi = 360 - phi
		        phi = phi * m.pi/180
		
		# Find the rotation matrix
		rotmat = rotmat2d(phi)

		# Rotate plane positions so vector u is aligned with x
		tmp = np.dot(rotmat,planepos[:,:2].T)
		planepos[:,:2] = tmp.T
		
		planepos = CleanMatElements(planepos)

		# Rotate vectors
		vec1r = np.array((0.0,0.0,0.0))
		vec2r = np.array((0.0,0.0,0.0))
		vec1r[:2] = np.dot(rotmat,vec1[:2])
		vec2r[:2] = np.dot(rotmat,vec2[:2])
		vec1r = CleanMatElements(vec1r)
		vec2r = CleanMatElements(vec2r)
		# If vec2 is pointing toward 3th and/or 4th quater of coordinate
		# system, flip the coordinates so they are in 1st and 2st 
		# quater 
		if vec2r[1] <0:
			planepos[:,1] *= -1
			vec2r[1] *= -1

		# Rotate anticipatory vectors to together with the plane

		for i in range(len(avecs)):
			tmp = avecs[i][:2]
			tmp = np.dot(rotmat,tmp)
			avecs[i][:2] = tmp

		if periodic: 
		# Find only unique atoms
			uniqueidx = self.__removeperiodic(planepos,vec1r,vec2r)
			planepos = planepos[uniqueidx]
#			tmpatm = []
#	
	#		for i in uniqueidx:
	#			tmpatm.append(planeatms[i])
	#		planeatms = tmpatm
			planeatms = planeatms[uniqueidx]

		return planepos, planeatms, avecs, [vec1r,vec2r], Hpos


	def __barycentric(self,p,p1,p2,p3):

		# Converstion to barycentric coordinates
		# Source:
		# http://en.wikipedia.org/wiki/Barycentric_coordinate_system_%28mathematics%29
		# TODO:
		# BE CAREFUL ABOUT THE ROUNDING!
		alpha = ((p2[1] - p3[1])*(p[0] - p3[0]) +\
			 (p3[0] - p2[0])*(p[1] - p3[1])) /\
			((p2[1] - p3[1])*(p1[0] - p3[0]) + \
			 (p3[0] - p2[0])*(p1[1] - p3[1]))

		beta =  ((p3[1] - p1[1])*(p[0] - p3[0])  +\
			 (p1[0] - p3[0])*(p[1] - p3[1])) /\
			((p2[1] - p3[1])*(p1[0] - p3[0]) + \
			 (p3[0] - p2[0])*(p1[1] - p3[1]))

		# Clean numeric noise - round to 6th place
		alpha = round(alpha,6)
		beta = round(beta,6)

		gamma = 1.0 - alpha - beta
		gamma = round(gamma,6)
		
		return alpha, beta, gamma

	def __removeperiodic(self,pos,vec1,vec2):

		# Find the atoms in the superlattice that are
		# translations of the other atoms 

		# Rude and not effecient solution 
		# Try more through array interception

		r = range(2)
		uniqlist=[]
		poscheck = []
		for x in r:
			for y in r:
				if x != 0 or y!= 0:
					if poscheck == []:
						poscheck = pos+x*vec1+y*vec2
					else:
						tmp = pos + x*vec1 + y*vec2
						poscheck = np.vstack\
						          ([poscheck,tmp])

		# Find indicies of unique elements in pos array that 
		# are not there is poscheck array
		ii = 0
		for i in pos:
			uq = True
			for j in poscheck:
				# Smaller accuracy required for z-axis
				# so round it to 5
				# Gives problems otherwise
				if round(i[0],8) == round(j[0],8) and \
				   round(i[1],8) == round(j[1],8) and \
				   round(i[2],8) == round(j[2],8):
					   uq = False
			if uq:
				if ii not in uniqlist:
					uniqlist.append(ii)
			ii += 1
		return uniqlist

	def __genHydrogen(self,pos,labels,refNN):

		# Generate the positions of the hydrogens on
		# the bottom of the position array. It can be either deposit
		# or substrate. Generate N-hydrogens, where N-is the number
		# of dangling bonds. 

		# pos - position of atoms after the cut along the planes
		# refNN -  number of nearest neighbeours of atoms in the 
		#          given material in bulk. Generated in "plane" 
		#          routine from "class Surface"

		# Introcude threshold to find number neartest neighbours
		# There are cases where atoms has N nearest neighbours, but
		# but they have only slighlty different distances, for instance
		# in the case of a-quartz-SiO2, Si has four NN, two has distance
		# 1.62A and two 1.63A. Use threshold to find them
		Nthresh = 0.1

		# Find the atoms on the bottom plane. At this point,
		# they have always negative z-components
		bottom = min(pos[:,2])
		# Add threshold, so one can avoid floating point comparrison, 
		# and it is highly unlikely that any other atom in the 
		# next layer on solid will have z-component only 0.05 
		# higher
		thresh = 0.7
		bottom += thresh
		# Find the indices of these atoms
		idx = pos[:,2] <= bottom
		postmp = pos[idx]

		# Generate atom labels for limited positions
		#tmplabels = []
		#for i in range(len(idx)):
		#	if idx[i]:
		#		tmplabels.append(labels[i])
		tmplabels = labels[idx]
		print tmplabels


		# Find point in the middle of the bottom plane by
		# finding centorid 

		cX = sum(postmp[:,0])/len(postmp[:,0])
		cY = sum(postmp[:,1])/len(postmp[:,1])
		# The z-variable of centroid is constant
		cZ = postmp[0,2]
		cXYZ = np.array((cX,cY,cZ))

		# Find the nearest neighbor of centroid. I will be our
		# reference atom.
		# If there are more than one atom types in the system,
		# find them too
			
		nnXYZ = fneigh(cXYZ,postmp)

		cidx = nnXYZ[0,1]
		cidx = int(cidx)
		catomlabel = tmplabels[cidx]
		uqlabels = []
		uqidx = []
		uqlabels.append(catomlabel) # store unique lables
		uqidx.append(cidx) # store unique indices

		# Find 2nd atom, if exists
		for i in nnXYZ[1:,1]:
			i = int(i)
			if tmplabels[i] not in uqlabels:
				uqlabels.append(tmplabels[i])
				uqidx.append(i)

		# Find the nearest neighbors of catoms in the whole set
		Hpostotal={}
		for ii in range(len(uqidx)):
			idx = uqidx[ii]
			catom = postmp[idx] # atom positions
			catomlab = uqlabels[ii] # atom labels
			nn = fneigh(catom,pos)
			shift = nn[1,0].copy()
			nn[:,0] -= shift
			nn = CleanMatElements(nn)
#		
			# Select only atoms with the shortest distance
#			idx = nn[:,0] == 0.0
			# take absolute value, as fist element in nn
			# is the reference atom itself, hence it will give
			# negative value of the distance
			idx = abs(nn[:,0]) <= Nthresh
			print nn[idx]
			nnlist = nn[idx]

			# Number of hyrdogens will be equal to number of 
			# number of nearest nieghbors in the bulk-number of 
			# nearest neighbours on the bottom surface

			print "catomlab",catomlab
			bulkNN = refNN[catomlab]
			print "refNN",refNN
			print "bulkNN",bulkNN
			nH = bulkNN - len(nnlist)

			# Generete list with hydrogen coordinate, that will
			# be later attached to each atom on the bottom surface

			# Method to generate hydrogens poistions:
			#  - put hydrogens on the bottom of the cone, where 
			# top of the
			#    cone is attached to the surface atom,
			#  - height of the cone is Hz=1.41
			#  - radius of the bottom = 1
			#  If there is only one hydrogen - put it along the 
			#  height of the cone
			#  If there is more than one hydrogen, put in on the 
			# circumferecne of the cone 360/nH apart of each other
			# It is achieved by doing nH rotations of the initial
			# point poistioned on the circumference of the cone 
			# bottom

#			Hz = -0.5
#			Hz = -1.70
#			Hz = -2.0
			Hz  = -1.11 # then bond length = 1.5
			Hzs = -1.49 # bond length when there is only one hydrogen
			print "NH",nH
			if nH <= 0: Hpos = []
			if nH == 1:
#				Hpos = np.array((0.0, 0.0, Hz))
				Hpos = np.zeros((1,3))
#				# TMP for SiO2
#				Hpos[0,0] = 1.0 
#				# END TMP for SiO2
				Hpos[0,2] = Hzs
			elif nH>1:
				Hpos = np.zeros((nH,3))
				Hpos[0] = np.array((1.0, 0.0, Hz))
	
				phi = 360/nH
			        phi = phi * m.pi/180
				rotmat = rotmat2d(phi)
				for i in range(nH)[1:]:
					Hpos[i][:2] = np.dot(rotmat,\
							     Hpos[i-1][:2])
					Hpos[i][2] = Hz
			Hpostotal[catomlab] = Hpos

		return Hpostotal

	def __getHydrogenPos(self,pos,labels,Hpos,top=False):

		# Add hydrogens to the each atom on the bottom of the slab

		# Find indicies of atom at the bottom. For explenation,
		# see comment to "def __genHydrogen"

		# if top=False - add hydrogens to the bottom of the slab, 
		#                i.e. substrate
		#
		# if top=True - add hydrogens to the top of the slab, 
		#               i.e. deposit

		# Return numpy array with hydrogen positions

		# Add threshold, so one can avoid floating point comparrison, 
		# and it is highly unlikely that any other atom in the 
		# next layer on solid will have z-component only 0.05 
		# higher
		thresh = 0.7

		if top: # Deposit
			plane = round(max(pos[:,2]),6)
			plane -= thresh
			idx = pos[:,2] >= plane
			for i in Hpos:
				Hpos[i] *= -1
		else:   # Substrate
			plane = round(min(pos[:,2]),6)
			plane += thresh
			idx = pos[:,2] <= plane

		postmp = pos[idx]

		# Generate atom labels for limited positions
		tmplabels = []
		for i in range(len(idx)):
			if idx[i]:
				tmplabels.append(labels[i])

		# Tmp list to hold H positions
		Hxyztmp = []
		
		# For each atom in the layer, check its atom label
		# and find appropirate H vectors in Hpos dictionary.
		# Than for each vector found in Hpos, generate hydrogen 
		# positions. 
		for i in range(len(postmp)):
			pvec = postmp[i] 
			atom = tmplabels[i]
			Hvec = Hpos[atom] # find H vectors for given atom
			
			#For each vector in the set, get H postions
			for j in Hvec:
				Hxyztmp.append(j+pvec)

		Hxyz = np.array(Hxyztmp)
		Hlabels = ["H"]*len(Hxyz)

		return Hxyz,Hlabels

	def __addHydrogens(self,posin,labels,idxD,idxS,HposD,HposS,genHD,genHS):

		# Add hydrogens to the bottom of substrate and top of deposit
		# 
		# Hydrogen coordinates are added using "__getHydrogenPos" 
		# routine, where "HposD" and "HposS" are arrays of hydrogens
		# atom positions for a single atom of the surface obtained
		# from  "__genHydrogen" subroutine.

		# Add hydrogens in such away that hydrogens connected to
		# Deposit are after Deposit posistions in the position array,
		# and same for Substrate

		pos = posin.copy()
		
		if genHD:
			DHxyz,DHlab = self.__getHydrogenPos\
				      (pos,labels,HposD,top=True)
		else:
			DHxyz = []
			DHlab = []
	
		if genHS:
			SHxyz,SHlab = self.__getHydrogenPos\
				      (pos,labels,HposS,top=False)
  	        else:
			SHxyz = []
			SHlab = []


		natD = len(idxD)
		natS = len(idxS)
		natDH = len(DHxyz)
		natSH = len(SHxyz)
		nattot = natD + natS + natDH + natSH
		posout = np.zeros((nattot,3))

#		pos = np.vstack([pos,DHxyz])
#		pos = np.vstack([pos,SHxyz])
		#idxDep = range(0,len(posDep))
		#idxSub = range(len(posDep),len(posout))
	
		# Deposit + Deposit hydrogens positions
		posout[0:natD] = pos[idxD]
		if genHD and DHxyz != []:
			posout[natD:natD+natDH] = DHxyz

		# Substrate + Substrate hydrogens positions
		posout[natD+natDH:natD+natDH+natS] = pos[idxS]
		if genHS and SHxyz != []:
			posout[natD+natDH+natS:] = SHxyz
		
		# Update indices of atoms of Deposit and Substrate to include
		# hydrogens
		idxDH = range(natD+natDH)
		idxSH = range(natD+natDH,len(posout))

		idxD = range(natD)
		idxS = range(natD+natDH,natD+natDH+natS)

#		labels = labels + DHlab + SHlab

		labelsD = labels[:natD]
		labelsS = labels[natD:]

		labels = labelsD + DHlab + labelsS + SHlab 
		# We have added hydrogens below the substrate. 
		# Shift everything up, so hydrogen atom are at z=0
		m = min(posout[:,2])
		posout[:,2] -= m
		
		# Shift the position for additional 2 angstroms, so 
		# in the case when hydrogens position are optimized, they
		# dont get shiffted to the top of tho box
		posout[:,2] += 2

		ii = 0

		#for i in pos:
	#		print labels[ii],i[0],i[1],i[2]
#			ii += 1

		return posout,labels,idxD,idxS,idxDH,idxSH 

	def __extendsurface(self,pos,labels,vec1,vec2,periodic=True):
		# Extend the surface, by adding copies of the surface
		# in xy plane

#		r = range(-1,2)
		r = range(-2,3)

	#	lx = np.linalg.norm(vec1)
	#	ly = np.linalg.norm(vec2)

		posout = pos.copy()
		newlabels = []
		for i in labels:
			newlabels.append(i)
#		print "VEC2", vec2
		for x in r:
			for y in r:
				if x != 0 or y!= 0:
#					newpos = np.zeros((len(pos),3))
	#				newpos[:,0] = pos[:,0] + x*lx 
#					newpos[:,1] = pos[:,1] + y*ly 
					newpos = pos + x*vec1 + y*vec2
					posout = \
					      np.vstack([posout,newpos])
					for i in labels:
						newlabels.append(i)

		posout = CleanMatElements(posout)

		# Remove all duplicates
		if not periodic:
			posout,newlabels = unique2(posout,newlabels)

		return posout,newlabels

	def __alignsurfaces(self,avecsDep,avecsSub,posDep,posSub,labDep,\
			    labSub,bfrac,vecpair,vecSubR):

		# Pick appropriate pair of anticipatory 
		# vectors from substrate and deposit 
		avecsSidx = vecpair[1]
		avecsDidx = vecpair[2]
		if (avecsSidx >= 0) and (avecsDidx >= 0): 
			avecS = avecsSub[avecsSidx]
			avecD = avecsDep[avecsDidx]
		else:
			# if indices are negative, it means alignment 
			# along vertical pointing vectors
			avecS = np.array([0,0,np.linalg.norm(avecsSub[0])])
			avecD = np.array([0,0,np.linalg.norm(avecsDep[0])])
		
		# Define resulting anticipatory vector
		vector = avecS+avecD
		vector = avecS

		# Define spacing between layers as average of substrate and 
		# deposit spacing
#		vector[2] /= 2 enable if avecS+avecD
#		vector[0] /= 2 enable if avecS+avecD
#		vector[1] /= 2 enable if avecS+avecD
		vn = np.linalg.norm(vector)

		# Calculate cosine between vector and z-axis
#		nv = np.linalg.norm(vector)
#		zvec = np.array((0.0,0.0,vector[2]))
#		nz = np.linalg.norm(zvec)
#		cosphi = np.dot(vector,zvec)/(nv*nz)
#		bfrac = 1.802/cosphi
#		print "BFRAC", cosphi,bfrac
		#For the straight vector, bfrac = 0.9
#		if (avecsSidx < 0) and (avecsDidx < 0): 
#			bfrac = 0.9


		# Our system are oriented in the 1st quarter of coordinate system
		# Force anticipatory vector to point there to
		if vector[0] <0: vector[0] *= -1
		if vector[1] <0: vector[1] *= -1

		# start TEST covalent radius displacement
		covrad = 2.35# For Si
		lv = np.linalg.norm(vector)
#		bfrac = covrad/lv
		# end TEST 
		#bfrac = 0.9
		bfrac = bfrac/lv
		vector = vector*bfrac
		#vector[2]=1.8024400000000007
		#vector[2]=bfrac

		
		# Now, move the deposit according the resulting anticipatory 
		# vector

		# Reverse positions so they are put on top of the 
		# substrate face down
		#posDep = self.DepSurfPos
		posDeptmp = posDep.copy()
		posDeptmp[:,2] = posDep[:,2] * -1
		# Rorate deposit by 90 deg along z-axis (so its not just mirror
		# image)
		angle = 180 * m.pi/180
		rotmat = rotmat2d(angle)
		tmp = np.dot(rotmat,posDeptmp[:,:2].T)
		posDeptmp[:,:2] = tmp.T
		posDeptmp = posDeptmp + vecSubR[0] + vecSubR[1]
		# Translate according to vector
		posDeptmp = posDeptmp + vector 

		# Save Substrate + Deposit positions
		posout = np.zeros((len(posDep)+len(posSub),3))
		
		#Deposit first
		posout[0:len(posDep)] = posDeptmp
		posout[len(posDep):] = posSub

		idxDep = np.arange(0,len(posDep))
		idxSub = np.arange(len(posDep),len(posout))

		#Create labels for the interface
		#labels = labDep + labSub
		labels = np.concatenate((labDep,labSub))

		# Shift the coordinates so the corner of the box is 0,0,0
		shift = min(posout[:,2])
		posout[:,2] -= shift
		
		# Return also vector which will be used as z-vector for the
		# direction vector of the interface
		# Find z-variable of top of the interface
		# TODO:
		# We are taking advantage that interface is parallel to x-y
		# plane, and we can find end of the vetor by taking 
		# biggest value of z-coordiante in the interface and scaling
		# anticpatory vector accordingly
		# More elegant solution: find the intercept of the line one
		# can draw aling anticipatory vector with the plane designated
		# by top atoms in the interface
		top = max(posout[:,2])
		ratio = top/vector[2]
		vector *= ratio

		topSub = max(posout[idxSub][:,2])
		botDep = min(posout[idxDep][:,2])
		spacing = botDep - topSub

		vector[-1] += spacing

		return posout, labels, idxDep, idxSub, np.array(vector)


	def __scale(self,pos,vec1D,vec2D,vec1S,vec2S,poissonRatio):
		
		# Scale the coordinates such that deposit
		# is always matching substrate


		posout = pos.copy()

		n1D = np.linalg.norm(vec1D)
		n2D = np.linalg.norm(vec2D)

		n1S = np.linalg.norm(vec1S)
		n2S = np.linalg.norm(vec2S)

		cosphi = np.dot(vec1D,vec2D)/(n1D*n2D)
		phiD = m.acos(cosphi)*180/m.pi

		cosphi = np.dot(vec1S,vec2S)/(n1S*n2S)
		phiS = m.acos(cosphi)*180/m.pi

		vec2Dr = vec2D.copy()

		# If angle between the vectors is different
		# scale it and the postitions accordingly
		# by rotating vec2D and saving it to vec2Dr
		if round(phiS,4) != round(phiD,4):
			delta = phiS - phiD
			# Back to radians
			delta = delta*m.pi/180
			rotmat = rotmat2d(delta)
			tmp = np.dot(rotmat,posout[:,:2].T)
			posout[:,:2] = tmp.T
			posout = CleanMatElements(posout)

			vec2Dr[:2] = np.dot(rotmat,vec2D[:2])
			vec2Dr = CleanMatElements(vec2Dr)

		vec1scale = np.array((1.0,1.0,1.0))
		vec2scale = np.array((1.0,1.0,1.0))

		# Scale vector 1
		if not (vec1D == vec1S).all():
			xs = vec1S[0]/vec1D[0]
			if vec1D[1] == 0.0:
				ys = 1.0
			else:
				ys = vec1S[1]/vec1D[1]
			vec1scale = np.array((xs,ys,1.0))

		# Scale vector 2
		if not (vec2Dr == vec2S).all():
			if vec2Dr[0] == 0:
				xs = 1.0
			else:
				xs = vec2S[0]/vec2Dr[0]
			ys = vec2S[1]/vec2Dr[1]
			vec2scale = np.array((xs,ys,1.0))

		# Scale the positons:
		posout = posout * vec1scale * vec2scale

		#poissonRatio = True
		if poissonRatio:
			m1 = n1D/n1S # misfit in x
			m2 = n2D/n2S # misfit in y
			ma = phiD/phiS # misfit in angle

			mz = m1 * m2 * ma
			vec3scale = np.array((1.0, 1.0, mz))
			posout = posout * vec3scale

		return posout

	def __checkedge(self,posinp,vecs):

		# In the case z-vector is not pointing aling 0,0,z direction,
		# move the atoms from the side of the interface that lay
		# left of (x,y,Z) to the right side on the interface

		# Algoritm:
		# 1) Find normal to the plane defined by the vecors x and z
		# 2) For each atom in the interface, calculate dot product with
		#    normal. If >0 then point lay above the plane, and needs to
		#    be translated.
		# 3) Repatet step 1 and 2 for plane defined by y and z

		pos = posinp.copy()

		xvec = vecs[0]
		yvec = vecs[1]
		zvec = vecs[2]

		# x-z plane
		# Find normal
		n = np.cross(xvec,zvec)
		n = CleanMatElements(n)
		dd = np.dot(xvec,zvec)
		dd = round(dd,12)
		# For each atom check the dot product with normal
		# If >0, need to move atom (check righ-hand-rule for sign)

		result = np.dot(pos,n)
		result = CleanMatElements(result)
		idx = result > 0.0
		# Translate
		pos[idx] = pos[idx] + yvec

		# Same for y-z plane
		# Find normal
		n = np.cross(yvec,zvec)
		n = CleanMatElements(n)
		dd = np.dot(yvec,zvec)
		dd = round(dd,12)
		# For each atom check the dot product with normal
		# If <0, need to move atoma (check righ-hand-rule for sign)

		result = np.dot(pos,n)
		result = CleanMatElements(result)
		idx = result < 0.0
		# Translate
		pos[idx] = pos[idx] + xvec


		#for i in pos[idx]:
#		for i in pos:
#			print "Al", i[0],i[1],i[2]

		return pos

	def __addCapAtoms(self, pos, labels, capAtD, capAtS):

		thresh = 0.1

		top = max(pos[:,2]) # Deposit surface
		bottom = min(pos[:,2]) # Substrate surface

		idxtop = pos[:,2] >= top - 6*thresh
		idxbot = pos[:,2] <= bottom + thresh

		# Substitute labels in Deposit
		if capAtD:
			# Check if capAtD in atomTyp
			maxType = max(self.atomTyp)
			types = self.atomTyp.values()
			if capAtD not in types:
				self.atomTyp[maxType+1] = capAtD
				atIdx = maxType+1
			else:
				atIdx = types.index(capAtD)

			for i in range(len(idxtop)):
				atom = idxtop[i]
				if atom:
				#	labels[i] = capAtD
					labels[i] = atIdx

		# Substitute labels in Substrate
		if capAtS:
			# Check if capAtD in atomTyp
			maxType = max(self.atomTyp)
			types = self.atomTyp.values()
			if capAtS not in types:
				self.atomTyp[maxType+1] = capAtS
				atIdx = maxType+1
			else:
				atIdx = types.index(capAtS)

			for i in range(len(idxbot)):
				atom = idxbot[i]
				if atom:
					labels[i] = atIdx

		return labels

	def SCORE4(self,pos,idx1,idx2,avecs1,avecs2,latticevecs):
		# NUMBER OF BONDS PER UNIT ARES
		#
		# WEIGHTED BY BONDLENGTH
		# Q= sum(cos(delta)[1:Nav] - sum(cos(delta))[Nav+1,:]
		#
		# ORDERED BY ANGLE

		u = latticevecs[0]
		v = latticevecs[1]

		# Cosine Parameter
		f = 5.0

		z = np.array((0.0,0.0,1.0))
		nz = 1 # norm

		nb = len(avecs1)
		natms = len(idx1)
		dist = (np.linalg.norm(avecs1[0])+np.linalg.norm(avecs2[0]))/2
		# MAKE THE REFERENCE DISTANCE THE LENGTH OF ALIGNEMENT VECTOR

		#calculate angle of anticipatory vector to z axis
		na = np.linalg.norm(avecs1[0])
		cosphi = np.dot(avecs1[0],z)/(na*nz)
		phiref = m.acos(cosphi)
		phiref = phiref * 180/m.pi

		blist = []
		for i in idx1: # iterate surface atoms of material1
			bcount = []
			angles = []
			for j in idx2: # iterate atoms of material2
					d = np.linalg.norm(pos[i]-pos[j])
					delta = abs(d - dist)
					# Calculate the angle
					vec = pos[i]-pos[j]
					if vec[2] <0:vec[2]*=-1
					nvec = np.linalg.norm(vec)
					cosphi = np.dot(vec,z)/(nvec*nz)
					phi = m.acos(cosphi)
					phi = phi * 180/m.pi

					cut = m.acos(0)/f # where cos ends
					w = m.cos(f*delta)
					if delta > cut: w = 0
					bcount.append(w)
					# Store the difference of the 
					# bond angle and anticipatory
					# vector angle
					angles.append(abs(phi-phiref))
			# Sorting operations
			bcount = np.array(bcount)
			angles = np.array(angles)
#			# Get the indicies after sorting
#			idx = bcount.argsort()
#			bcount = bcount[idx]
#			bcount = bcount[::-1] #reverse array
#			angles = angles[idx]
#			angles = angles[::-1]
#			
#			# Sort the bonds according to the angle bond-anticpatory
#			# vector angle difference
#			bidx = bcount > 0
#			bcounttmp = bcount[bidx] # only bonds >0
#			anglestmp = angles[bidx] # angles for bonds >0
#
#			idx = anglestmp.argsort()
#			bcounttmp = bcounttmp[idx]
#			anglestmp = anglestmp[idx]
#
#			Q1 = sum(bcounttmp[0:nb])
#			Q2 = sum(bcounttmp[nb:])
#
			# Other way - sort angles and bond later
			idx = angles.argsort()
			angles = angles[idx]
			bcount = bcount[idx]
			Q1 = sum(bcount[0:nb])
			Q2 = sum(bcount[nb:])
		
			Q = Q1-Q2
			if Q < 0: Q = 0

			blist.append(Q)

		nbonds = sum(blist)
		
		#CALCULATE UNIT AREA
		normu = np.linalg.norm(u)
		normv = np.linalg.norm(v)
		cosphi = np.dot(u,v)/(normu*normv)
		sinphi = m.sqrt(1-cosphi**2)
		area = normu*normv*sinphi

		score = nbonds/area
		#Reference
		scoreref = (natms*nb)/area
		score = score/scoreref

		return score 


	def SCORE5(self,pos,idx1,idx2,avecs1,avecs2,latticevecs):
		# NUMBER OF BONDS PER UNIT ARES
		#
		# WEIGHTED BY BONDLENGTH
		# Q= sum(cos(delta)[1:Nav] - sum(cos(delta))[Nav+1,:]
		#
		# WEIGHTED BY ANGLE
		#
		# ORDERED BY BOND

		u = latticevecs[0]
		v = latticevecs[1]

		# Cosine Parameter
		f = 5.0

		z = np.array((0.0,0.0,1.0))
		nz = 1 # norm

		nb = len(avecs1)
		natms = len(idx1)
		dist = (np.linalg.norm(avecs1[0])+np.linalg.norm(avecs2[0]))/2
		# MAKE THE REFERENCE DISTANCE THE LENGTH OF ALIGNEMENT VECTOR

		#calculate angle of anticipatory vector to z axis
		na = np.linalg.norm(avecs1[0])
		cosphi = np.dot(avecs1[0],z)/(na*nz)
		phiref = m.acos(cosphi)
		phiref = phiref * 180/m.pi

		blist = []
		for i in idx1: # iterate surface atoms of material1
			bcount = []
			angles = []
			for j in idx2: # iterate atoms of material2
					d = np.linalg.norm(pos[i]-pos[j])
					delta = abs(d - dist)
					# Calculate the angle
					vec = pos[i]-pos[j]
					if vec[2] <0:vec[2]*=-1
					nvec = np.linalg.norm(vec)
					cosphi = np.dot(vec,z)/(nvec*nz)
					phi = m.acos(cosphi)
					phi = phi * 180/m.pi

					deltaphi = abs(phi-phiref)
					# To radians
					deltaphi *= (m.pi/180)
					cut = m.acos(0)/f # where cos ends

					w = m.cos(f*delta) # bond legth weight
					wa = m.cos(deltaphi) # angle weight 
					if delta > cut: w = 0
					bcount.append(w*wa)

			# Sorting	
			bcount.sort(reverse=True)
			Q1 = sum(bcount[0:nb])
			Q2 = sum(bcount[nb:])
		
			Q = Q1-Q2
			if Q < 0: Q = 0

			blist.append(Q)

		nbonds = sum(blist)
		
		#CALCULATE UNIT AREA
		normu = np.linalg.norm(u)
		normv = np.linalg.norm(v)
		cosphi = np.dot(u,v)/(normu*normv)
		sinphi = m.sqrt(1-cosphi**2)
		area = normu*normv*sinphi

		score = nbonds/area
		#Reference
		scoreref = (natms*nb)/area
		score = score/scoreref

		return score 

	def SCORE6(self,pos,idx1,idx2,avecs1,avecs2,latticevecs):
		# NUMBER OF BONDS PER UNIT ARES
		#
		# WEIGHTED BY BONDLENGTH
		# Q= sum(cos(delta)[1:Nav] - sum(cos(delta))[Nav+1,:]
		#
		# WEIGHTED BY ANGLE
		#
		# ORDERED BY ANGLE

		u = latticevecs[0]
		v = latticevecs[1]

		# Cosine Parameter
		f = 0.5
		fa = 14.0

		z = np.array((0.0,0.0,1.0))
		nz = 1 # norm

		nb = len(avecs1)
		natms = len(idx1)
		dist = (np.linalg.norm(avecs1[0])+np.linalg.norm(avecs2[0]))/2
		# MAKE THE REFERENCE DISTANCE THE LENGTH OF ALIGNEMENT VECTOR

		#calculate angle of anticipatory vector to z axis
		na = np.linalg.norm(avecs1[0])
		cosphi = np.dot(avecs1[0],z)/(na*nz)
		phiref = m.acos(cosphi)
		phiref = phiref * 180/m.pi

		blist = []
		for i in idx1: # iterate surface atoms of material1
			bcount = []
			angles = []
			bdelta = []
			for j in idx2: # iterate atoms of material2
					d = np.linalg.norm(pos[i]-pos[j])
					delta = abs(d - dist)
					# Calculate the angle
					vec = pos[i]-pos[j]
					if vec[2] <0:vec[2]*=-1
					nvec = np.linalg.norm(vec)
					cosphi = np.dot(vec,z)/(nvec*nz)
					phi = m.acos(cosphi)
					phi = phi * 180/m.pi

					deltaphi = abs(phi-phiref)
					# To radians
					deltaphi *= (m.pi/180)
					cut = m.acos(0)/f # where cos ends
					cuta = m.acos(0)/fa # where angle cos end

					w = m.cos(f*delta) # bond legth weight
					wa = m.cos(fa*deltaphi) # angle weight 
					if delta > cut: w = 0
					if deltaphi > cuta: wa = 0
					bcount.append(w*wa)
					bdelta.append(delta)
					angles.append(deltaphi*180/m.pi)

			bcount = np.array(bcount)
			angles = np.array(angles)
			bdelta= np.array(bdelta)

			idx = angles.argsort()
			angles = angles[idx]
			bcount = bcount[idx]
			bdelta = bdelta[idx] #For debugging only
			Q1 = sum(bcount[0:nb])
			Q2 = sum(bcount[nb:])

#			# Sorting	
#			bcount.sort(reverse=True)
#			Q1 = sum(bcount[0:nb])
#			Q2 = sum(bcount[nb:])

			Q = Q1-Q2
			if Q < 0: Q = 0
			Q = Q1

			print bcount[0:5]
			print angles[0:5]
			print bdelta[0:5]
			print Q
			print 

			blist.append(Q)

		nbonds = sum(blist)
		print "NBONDS",nbonds
		
		#CALCULATE UNIT AREA
		normu = np.linalg.norm(u)
		normv = np.linalg.norm(v)
		cosphi = np.dot(u,v)/(normu*normv)
		sinphi = m.sqrt(1-cosphi**2)
		area = normu*normv*sinphi

		score = nbonds/area
		#Reference
		scoreref = (natms*nb)/area
		score = score/scoreref

		print "SCORE",score

		return score 

	def SCORE7(self,pos,idx1,idx2,avecs1,avecs2,latticevecs,\
			labels,avecsall1,f,atomicRadius,mx,my,mphi):
		# NUMBER OF BONDS PER UNIT ARES
		#
		# WEIGHTED BY BONDLENGTH
		# Q= sum(cos(delta)[1:Nav] 
		#
		# ORDERED BY DELTA BOND LENGTH

		u = latticevecs[0]
		v = latticevecs[1]
		
		#Cosine parameter
#		f = 3.0
		f=float(f)

		nb = len(avecs1)
		#nb = 1
		natms = len(idx1)
		dist = (np.linalg.norm(avecs1[0])+np.linalg.norm(avecs2[0]))/2
	##	print "DIST IN SCORE7", dist
		dist = 2.20
		dist = atomicRadius
		# MAKE THE REFERENCE DISTANCE THE LENGTH OF ALIGNEMENT VECTOR

		blist = []
		for i in idx1: # iterate surface atoms of material1
			bcount = []
			bdist = []
			ilist = []
			atom = labels[i]
			if atom in avecsall1:
				nb = len(avecsall1[atom])
				for j in idx2: # iterate atoms of material2
						d = np.linalg.norm(pos[i]-pos[j])
						delta = abs(d - dist)
						cut = m.acos(0)/f # where cos ends
						w = m.cos(f*delta)
						bdist.append(d)
						if delta >cut: w = 0
						bcount.append(w)
				# Sort bcount
				#bcount.sort(reverse=True)
				bdist = np.array(bdist)
				bcount = np.array(bcount)
				idx=bdist.argsort()
				bcount=bcount[idx]
				Q1 = sum(bcount[0:nb])
				Q2 = sum(bcount[nb:])
				Q = Q1-Q2
				if Q < 0: Q = 0

#			print bcount[0:5]
#			print bdist[0:5]
#			print Q
#			print 

				blist.append(Q)

		nbonds = sum(blist)
		
		#CALCULATE UNIT AREA
		normu = np.linalg.norm(u)
		normv = np.linalg.norm(v)
		cosphi = np.dot(u,v)/(normu*normv)
		sinphi = m.sqrt(1-cosphi**2)
		area = normu*normv*sinphi

		score = nbonds/area
		#Reference
		scoreref2 = (natms*nb)/area
#		score = score/scoreref
		
		#Reference
		scoreref = 0 
		for i in idx1:
			atom = labels[i]
			if atom in avecsall1:
				nb = len(avecsall1[atom])
			else:
				nb = 0
			scoreref += nb
		scoreref /= area
	#	print "SCORE, SCORE2",score,scoreref
		score = score/scoreref

		score = score*mx*my*mphi

		return score

	def SCORE7metal(self,pos,idx1,idx2,avecs1,avecs2,latticevecs,f):
		# NUMBER OF BONDS PER UNIT ARES
		#
		# WEIGHTED BY BONDLENGTH
		# Q= sum(cos(delta)[1:Nav] 
		#
		# ORDERED BY DELTA BOND LENGTH
		#
		# FIXSED NUMBER OF BONDS
		#

		u = latticevecs[0]
		v = latticevecs[1]
		
		#Cosine parameter
#		f = 3.0
		f=float(f)

		nb = len(avecs1)
		#nb = 4
		natms = len(idx1)
		dist = (np.linalg.norm(avecs1[0])+np.linalg.norm(avecs2[0]))/2
##		print "DIST IN SCORE7", dist
		dist=2.35
		# MAKE THE REFERENCE DISTANCE THE LENGTH OF ALIGNEMENT VECTOR

		blist = []
		for i in idx1: # iterate surface atoms of material1
			bcount = []
			bdist = []
			ilist = []
			for j in idx2: # iterate atoms of material2
					d = np.linalg.norm(pos[i]-pos[j])
					delta = abs(d - dist)
					cut = m.acos(0)/f # where cos ends
					w = m.cos(f*delta)
					bdist.append(d)
					if delta >cut: w = 0
					bcount.append(w)
			# Sort bcount
			#bcount.sort(reverse=True)
			bdist = np.array(bdist)
			bcount = np.array(bcount)
			idx=bdist.argsort()
			bcount=bcount[idx]
			Q1 = sum(bcount[0:nb])
			if bcount[0] == 0.0:Q1=0
			Q2 = sum(bcount[nb:])
			Q = Q1#-Q2
			if Q < 0: Q = 0

			blist.append(Q)

		nbonds = sum(blist)
		
		#CALCULATE UNIT AREA
		normu = np.linalg.norm(u)
		normv = np.linalg.norm(v)
		cosphi = np.dot(u,v)/(normu*normv)
		sinphi = m.sqrt(1-cosphi**2)
		area = normu*normv*sinphi

		score = nbonds/area
		#Reference
		scoreref = (natms*nb)/area
		score = score/scoreref

		return score

def mkOutput(iface,confno,Dname,Dface,Sname,Sface,vecpairidx,nLS,nLD,mfit,\
				     score,alignNo,iterNo,nVac,sandwich,\
			             writeGEN=False,\
                                     writeAIMS=False,writeGULP=False,\
				     cstrlxD=0,cstrlxS=0,\
				     cstrlxH=False,\
				     bondmod=False): 
	#
	# Write to XYZ file with coordinates for given interface
	# Optional create DFTB+ GEN file format using external program
	#
	# cstrlx - optional create AIMS file with contrained atoms on the 
	# edge of deposit and substrate
	# and with constrained hydrogen atoms:
	# cstrlxHD - hydrogens on deposit
	# cstrlxHS - hydrogens on substrate
	# 

	norma = np.linalg.norm(iface.IfaceVecs[0])
	normb = np.linalg.norm(iface.IfaceVecs[1])

	cosphi = np.dot(iface.IfaceVecs[0],iface.IfaceVecs[1])/(norma*normb)
	phi = m.acos(cosphi)
	sinphi = m.sin(phi)

	area = norma * normb * sinphi

	# Find distance between substrate and deposit in straight line
	subPos = iface.IfacePosSC[iface.idxSubH]
	depPos = iface.IfacePosSC[iface.idxDepH]
	topSub = max(subPos[:,2]) # Substrate layer in the interface
	botDep = min(depPos[:,2]) # Deposit layer in the interface
	topDep = max(depPos[:,2]) # Top of deposit
	SDdist = botDep - topSub

#	if aligNo == 0 and confno == 0:
	if iterNo == 1:
		summaryFile = open('summary.txt','w')
		summaryFileCSV = open('summary.csv','w')
		summaryFile.write("%12s | %12s | %12s | %12s | %12s | %12s | %12s | %12s | %12s | %12s | %12s | %12s | %12s\n"\
			     %("Substrate","Sub. Orient.", "Deposit", "Dep. Orient.",\
			     "Alignment","Conf. No",\
		     	     "area","x - stress","y-stress","angle-stress","area-stress",\
		   	     "score","S-D distance"))

		summaryFileCSV.write("%12s , %12s , %12s , %12s , %12s , %12s , %12s , %12s , %12s , %12s , %12s , %12s , %12s\n"\
			     %("Substrate","Sub. Orient.", "Deposit", "Dep. Orient.",\
			     "Alignment","Conf. No",\
		     	     "area","x - stress","y-stress","angle-stress","area-stress",\
		   	     "score","S-D distance"))


	else:
		summaryFile = open('summary.txt','a')
		summaryFileCSV = open('summary.csv','a')

	summaryFile.write("%12s   %12s   %12s   %12s   %12i   %12i   %12.2f   %12.2f   %12.2f   %12.2f   %12.2f   %12.2f   %12.2f\n"\
			%(Sname,Sface,Dname,Dface,aligNo,confno,\
			  area,mfit[0],mfit[1],mfit[2],mfit[3],score,iface.SDdist))

	summaryFileCSV.write("%12s , %12s , %12s , %12s , %12i , %12i , %12.2f , %12.2f , %12.2f , %12.2f , %12.2f , %12.2f , %12.2f\n"\
			%(Sname,"'"+Sface+"'",Dname,"'"+Dface+"'",aligNo,confno,\
			  area,mfit[0],mfit[1],mfit[2],mfit[3],score,iface.SDdist))

	summaryFile.close()
	summaryFileCSV.close()
		       
	
	#dirname = "%s%s-%s%s/%i"%(Dname,Dface,Sname,Sface,vecpairidx)
	dirnameConf = "%s%s-%s%s-%i"%(Dname,Dface,Sname,Sface,confno)
	if bondmod:
		dirnameConf = "%s%s-%s%s-%2.1f"%(Dname,Dface,Sname,Sface,confno)
	dirname = "%s%s-%s%s/%i/%s"%(Dname,Dface,Sname,Sface,vecpairidx,dirnameConf)
	os.system("mkdir -p %s"%dirname)

	# Write which vectors we are using for alignement
	ftmp = open("%s/vectors.txt"%dirname,'w')
	ftmp.write("%5i   %5i\n"%(vecpair[1],vecpair[2]))
	ftmp.close()
	
	if bondmod:
		filename  = "%s/%s%s-%s%s-%2.1f.xyz"%(dirname,Dname,Dface,\
			                   Sname,Sface,confno)
		filenameS = "%s/%s%s-%s%s-%2.1f-S.xyz"%(dirname,Dname,Dface,\
			                   Sname,Sface,confno)
		filenameD = "%s/%s%s-%s%s-%2.1f-D.xyz"%(dirname,Dname,Dface,\
			                   Sname,Sface,confno)
		if writeAIMS:
			filenameAIMS = "%s/%s%s-%s%s-%2.1f.in"\
				%(dirname,Dname,Dface,Sname,Sface,confno)
			filenameSAIMS = "%s/%s%s-%s%s-%2.1f-S.in"\
				%(dirname,Dname,Dface,Sname,Sface,confno)
			filenameDAIMS = "%s/%s%s-%s%s-%2.1f-D.in"\
				%(dirname,Dname,Dface,Sname,Sface,confno)

		if writeGULP:
			filenameGULP = "%s/%s%s-%s%s-%2.1f.gin"\
				%(dirname,Dname,Dface,Sname,Sface,confno)
			filenameGULPS = "%s/%s%s-%s%s-%2.1f-S.gin"\
				%(dirname,Dname,Dface,Sname,Sface,confno)
			filenameGULPD = "%s/%s%s-%s%s-%2.1f-D.gin"\
				%(dirname,Dname,Dface,Sname,Sface,confno)

		# Index file with number of atoms in 
		# Interface,Deposit and Substate, and the area
		filenameIDX  = "%s/%s%s-%s%s-%2.1f.idx"%(dirname,Dname,Dface,\
			                   Sname,Sface,confno)

	else:
		filename  = "%s/%s%s-%s%s-%i.xyz"%(dirname,Dname,Dface,\
			                   Sname,Sface,confno)
		filenameS = "%s/%s%s-%s%s-%i-S.xyz"%(dirname,Dname,Dface,\
			                   Sname,Sface,confno)
		filenameD = "%s/%s%s-%s%s-%i-D.xyz"%(dirname,Dname,Dface,\
			                   Sname,Sface,confno)
		if writeAIMS:
			filenameAIMS  = "%s/%s%s-%s%s-%i.in"\
				%(dirname,Dname,Dface,Sname,Sface,confno)
			filenameSAIMS = "%s/%s%s-%s%s-%i-S.in"\
				%(dirname,Dname,Dface,Sname,Sface,confno)
			filenameDAIMS = "%s/%s%s-%s%s-%i-D.in"\
				%(dirname,Dname,Dface,Sname,Sface,confno)
		if writeGULP:
			filenameGULP  = "%s/%s%s-%s%s-%i.gin"\
				%(dirname,Dname,Dface,Sname,Sface,confno)
			filenameGULPS = "%s/%s%s-%s%s-%i-S.gin"\
				%(dirname,Dname,Dface,Sname,Sface,confno)
			filenameGULPD = "%s/%s%s-%s%s-%i-D.gin"\
				%(dirname,Dname,Dface,Sname,Sface,confno)

		# Index file with number of atoms in 
		# Interface,Deposit and Substate
		filenameIDX  = "%s/%s%s-%s%s-%i.idx"%(dirname,Dname,Dface,\
			                   Sname,Sface,confno)

	file  = open(filename,'w')
	fileS = open(filenameS,'w')
	fileD = open(filenameD,'w')
	if writeAIMS:
		fileAIMS  = open(filenameAIMS,'w')
		fileSAIMS = open(filenameSAIMS,'w')
		fileDAIMS = open(filenameDAIMS,'w')

	if writeGULP:
		fileGULP  = open(filenameGULP,'w')
		fileGULPS = open(filenameGULPS,'w')
		fileGULPD = open(filenameGULPD,'w')
	
	fileIDX = open(filenameIDX,'w')

	natoms  = len(iface.IfacePosSC)
	# in the case there is no hydrogens added - iface.idxSubH and 
	# iface.idxDepH are set to iface.idxSub and iface.idxDep in
	# class Interface
	natomsS = len(iface.idxSubH) 
	natomsD = len(iface.idxDepH)

	# For isolated deposit, we want to have it shifted to the bottom 
	# of the box
	Dshift = min(iface.IfacePosSC[iface.idxDepH,2])

	file.write("%i\n"%natoms)
	file.write("\n")

	fileS.write("%i\n"%natomsS)
	fileS.write("\n")

	fileD.write("%i\n"%natomsD)
	fileD.write("\n")

	# If outputting AIMS, write the lattice vectors
	if writeAIMS:
		print "Generating AIMS .in file"
		fileAIMS.write("lattice_vector %12.6f  %12.6f  %12.6f\n"%\
		             (iface.IfaceVecs[0][0],\
			     iface.IfaceVecs[0][1],iface.IfaceVecs[0][2]))
		fileAIMS.write("lattice_vector %12.6f  %12.6f  %12.6f\n"%\
		             (iface.IfaceVecs[1][0],\
			      iface.IfaceVecs[1][1],iface.IfaceVecs[1][2]))
		fileAIMS.write("lattice_vector %12.6f  %12.6f  %12.6f\n"%\
		             (iface.IfaceVecs[2][0]*nVac,\
			      iface.IfaceVecs[2][1]*nVac,\
			      iface.IfaceVecs[2][2]*nVac))

		fileSAIMS.write("lattice_vector %12.6f  %12.6f  %12.6f\n"%\
		             (iface.IfaceVecs[0][0],\
			     iface.IfaceVecs[0][1],iface.IfaceVecs[0][2]))
		fileSAIMS.write("lattice_vector %12.6f  %12.6f  %12.6f\n"%\
		             (iface.IfaceVecs[1][0],\
			      iface.IfaceVecs[1][1],iface.IfaceVecs[1][2]))
		fileSAIMS.write("lattice_vector %12.6f  %12.6f  %12.6f\n"%\
		             (iface.IfaceVecs[2][0]*nVac,\
			      iface.IfaceVecs[2][1]*nVac,\
			      iface.IfaceVecs[2][2]*nVac))

		fileDAIMS.write("lattice_vector %12.6f  %12.6f  %12.6f\n"%\
		             (iface.IfaceVecs[0][0],\
			     iface.IfaceVecs[0][1],iface.IfaceVecs[0][2]))
		fileDAIMS.write("lattice_vector %12.6f  %12.6f  %12.6f\n"%\
		             (iface.IfaceVecs[1][0],\
			      iface.IfaceVecs[1][1],iface.IfaceVecs[1][2]))
		fileDAIMS.write("lattice_vector %12.6f  %12.6f  %12.6f\n"%\
		             (iface.IfaceVecs[2][0]*nVac,\
			      iface.IfaceVecs[2][1]*nVac,\
			      iface.IfaceVecs[2][2]*nVac))
	
	# Outputting GULP
	if writeGULP:
		print "Generating GULP .gin file"
		fileGULP.write("opti\n")
		fileGULP.write("maxcyc opt 20000\n")
		fileGULP.write("dump\n")
		fileGULP.write("vectors\n")
		fileGULP.write("%12.6f  %12.6f  %12.6f\n"%\
                             (iface.IfaceVecs[0][0],\
                              iface.IfaceVecs[0][1],iface.IfaceVecs[0][2]))
		fileGULP.write("%12.6f  %12.6f  %12.6f\n"%\
                             (iface.IfaceVecs[1][0],\
                              iface.IfaceVecs[1][1],iface.IfaceVecs[1][2]))
		fileGULP.write("%12.6f  %12.6f  %12.6f\n"%\
		             (iface.IfaceVecs[2][0]*nVac,\
			      iface.IfaceVecs[2][1]*nVac,\
			      iface.IfaceVecs[2][2]*nVac))
		fileGULP.write("0 0 0 0 0 0\n")
		fileGULP.write("cartesian\n")

		fileGULPS.write("opti\n")
		fileGULPS.write("maxcyc opt 20000\n")
		fileGULPS.write("dump\n")
		fileGULPS.write("vectors\n")
		fileGULPS.write("%12.6f  %12.6f  %12.6f\n"%\
                             (iface.IfaceVecs[0][0],\
                              iface.IfaceVecs[0][1],iface.IfaceVecs[0][2]))
		fileGULPS.write("%12.6f  %12.6f  %12.6f\n"%\
                             (iface.IfaceVecs[1][0],\
                              iface.IfaceVecs[1][1],iface.IfaceVecs[1][2]))
		fileGULPS.write("%12.6f  %12.6f  %12.6f\n"%\
		             (iface.IfaceVecs[2][0]*nVac,\
			      iface.IfaceVecs[2][1]*nVac,\
			      iface.IfaceVecs[2][2]*nVac))
		fileGULPS.write("0 0 0 0 0 0\n")
		fileGULPS.write("cartesian\n")

		fileGULPD.write("opti\n")
		fileGULPD.write("maxcyc opt 20000\n")
		fileGULPD.write("dump\n")
		fileGULPD.write("vectors\n")
		fileGULPD.write("%12.6f  %12.6f  %12.6f\n"%\
                             (iface.IfaceVecs[0][0],\
                              iface.IfaceVecs[0][1],iface.IfaceVecs[0][2]))
		fileGULPD.write("%12.6f  %12.6f  %12.6f\n"%\
                             (iface.IfaceVecs[1][0],\
                              iface.IfaceVecs[1][1],iface.IfaceVecs[1][2]))
		fileGULPD.write("%12.6f  %12.6f  %12.6f\n"%\
		             (iface.IfaceVecs[2][0]*nVac,\
			      iface.IfaceVecs[2][1]*nVac,\
			      iface.IfaceVecs[2][2]*nVac))
		fileGULPD.write("0 0 0 0 0 0\n")
		fileGULPD.write("cartesian\n")
	#Find top and bottom layer
	threshold = 0.2 # define threshol to add to top and bottom for 
	                # float comparisons
	frozenidx=[]
	if (cstrlxS):
#		order  = np.sort(iface.IfacePosSC[iface.idxSub,2])
#		order = np.unique(order)
#		bot  = round(order[cstrlxS-1],6) # bottom atom (substrate)
#		bot += threshold
		bot = nLS/3.0
		if sandwich:
			topS = topDep + iface.SDdist + 2*bot #1/3 of the top part of substarte in sandwich
			botS = botDep - iface.SDdist - 2*bot #1/3 of the bottom part of substrate
			topI = iface.IfacePosSC[:,2] > topS
			botI = iface.IfacePosSC[:,2] < botS

			for i in range(len(topI)):
				if topI[i]:
					frozenidx.append(i)

		else:
			botI = iface.IfacePosSC[:,2] < bot

		for i in range(len(botI)):
			if botI[i]:
				frozenidx.append(i)

	if (cstrlxD):
#		order  = np.sort(iface.IfacePosSC[iface.idxDep,2])
#		order = np.unique(order)
#		top = round(order[(-1*cstrlxD)-1],6) # top atom (deposit)
#		top -= threshold
		top = nLS +iface.SDdist + (2*(nLD/3.0))
		if sandwich:
			dFrac = nLD/3.0
			mTop = topDep - dFrac
			mBot = botDep + dFrac
			mFrozenI=(iface.IfacePosSC[:,2] <mTop) & (iface.IfacePosSC[:,2] >mBot)

			for i in range(len(mFrozenI)):
				if mFrozenI[i]:
					frozenidx.append(i)
		else: 
			topI = iface.IfacePosSC[:,2] > top

			for i in range(len(topI)):
				if topI[i]:
					frozenidx.append(i)
	frozenidx=np.array(frozenidx)

	# Freeze the middle atoms
	# Freeze 1/3 of the atom in S and D
	# S: 1/3 * nL :  2/3 *nL
	# D: nL + 1/3 * nL : nl+2/3 * nL

	mBottomS = nLS/3.0 #middle later Bottom of Substrate
	mTopS    = 2*(nLS/3.0)

	mBottomD = nLD + (nLD/3.0)
	mTopD    = nLD + (2*(nLD/3.0))

	for i in range(len(iface.IfacePosSC)):
		atom = iface.atomTyp[iface.IfaceAtmSC[i]]
		x = iface.IfacePosSC[i][0]
		y = iface.IfacePosSC[i][1]
		z = iface.IfacePosSC[i][2]
		file.write("%-4s   %12.6f   %12.6f   %12.6f\n"%(atom,x,y,z))
		if writeAIMS:
			fileAIMS.write("atom  %12.6f   %12.6f   %12.6f  %4s\n"%\
				      (x,y,z,atom))
			# Constrain deposit hydrogen
			if (cstrlxH) and atom == "H":
				#fileAIMS.write("constrain_relaxation .true.\n")
				fileAIMS.write("constrain_relaxation x\n")
				fileAIMS.write("constrain_relaxation y\n")

			if i in frozenidx:
				fileAIMS.write\
				("constrain_relaxation x\n")
				fileAIMS.write\
				("constrain_relaxation y\n")

#			# Constrain middle layers
#			if ((round(z,6) >= mBottomS ) and (round(z,6) <= mTopS) ):
#				fileAIMS.write\
#				("constrain_relaxation x\n")
#				fileAIMS.write\
#				("constrain_relaxation y\n")
#			if ((round(z,6) >= mBottomD ) and (round(z,6) <= mTopD) ):
#				fileAIMS.write\
#				("constrain_relaxation x\n")
#				fileAIMS.write\
#				("constrain_relaxation y\n")

		if writeGULP:
			if (cstrlxH) and atom == "H":
				fileGULP.write(("%4s    core   %12.6f   %12.6f"
				"%12.6f 0.0 1.0 0.0 0 0 1\n")%(atom,x,y,z))
			elif (((cstrlxS) or (cstrlxD)) and \
			    ((round(z,6) <= bot) or (round(z,6) >= top) )):
				fileGULP.write(("%4s    core   %12.6f   %12.6f"
				"%12.6f 0.0 1.0 0.0 0 0 1\n")%(atom,x,y,z))
			else:
				fileGULP.write(("%4s    core   %12.6f   %12.6f"
				"%12.6f 0.0 1.0 0.0 1 1 1\n")%(atom,x,y,z))


		# Write isolated Substrate and Deposit
		# Only XYZ and AIMS file are written. Since 
		# Since GULP is used only for optimization, 
		# the isolated not opimized monomers are not 
		# very usefull.

		if i in iface.idxSubH: # Write isolated Substrate
			fileS.write("%-4s   %12.6f   %12.6f   %12.6f\n"%\
				   (atom,x,y,z))
			if writeAIMS:
				fileSAIMS.write\
				     ("atom  %12.6f   %12.6f   %12.6f  %4s\n"%\
			             (x,y,z,atom))

				if (cstrlxH) and atom == "H":
#					fileSAIMS.write\
#					    ("constrain_relaxation .true.\n")
					fileSAIMS.write\
					    ("constrain_relaxation x\n")
					fileSAIMS.write\
					    ("constrain_relaxation y\n")
				if (cstrlxS) and (round(z,6) <= bot):
#					fileSAIMS.write\
#					    ("constrain_relaxation .true.\n")
					fileSAIMS.write\
					    ("constrain_relaxation x\n")
					fileSAIMS.write\
					    ("constrain_relaxation y\n")
			if writeGULP:
				if (cstrlxH) and atom == "H":
					fileGULPS.write\
					(("%4s    core   %12.6f   %12.6f"
					"%12.6f 0.0 1.0 0.0 0 0 1\n")\
							%(atom,x,y,z))

				elif (cstrlxS) and (round(z,6) <= bot):
					fileGULPS.write\
					(("%4s    core   %12.6f   %12.6f"
					"%12.6f 0.0 1.0 0.0 0 0 1\n")\
							%(atom,x,y,z))
				else:
					fileGULPS.write\
					(("%4s    core   %12.6f   %12.6f"
					"%12.6f 0.0 1.0 0.0 1 1 1\n")\
							%(atom,x,y,z))


		if i in iface.idxDepH: # Write isolated Deposit
			fileD.write("%-4s   %12.6f   %12.6f   %12.6f\n"%\
				   (atom,x,y,z-Dshift))
			if writeAIMS:
				fileDAIMS.write\
				     ("atom  %12.6f   %12.6f   %12.6f  %4s\n"%\
			             (x,y,z,atom))

				if (cstrlxH) and atom == "H":
#					fileDAIMS.write\
#					    ("constrain_relaxation .true.\n")
					fileDAIMS.write\
					    ("constrain_relaxation x\n")
					fileDAIMS.write\
					    ("constrain_relaxation y\n")
				if (cstrlxD) and (round(z,6) >= top):
#					fileDAIMS.write\
#					    ("constrain_relaxation .true.\n")
					fileDAIMS.write\
					    ("constrain_relaxation x\n")
					fileDAIMS.write\
					    ("constrain_relaxation y\n")

			if writeGULP:
				if (cstrlxH) and atom == "H":
					fileGULPD.write\
					(("%4s    core   %12.6f   %12.6f"
					"%12.6f 0.0 1.0 0.0 0 0 1\n")\
							%(atom,x,y,z))

				elif (cstrlxD) and (round(z,6) >= top):
					fileGULPD.write\
					(("%4s    core   %12.6f   %12.6f"
					"%12.6f 0.0 1.0 0.0 0 0 1\n")\
							%(atom,x,y,z))
				else:
					fileGULPD.write\
					(("%4s    core   %12.6f   %12.6f"
					"%12.6f 0.0 1.0 0.0 1 1 1\n")\
							%(atom,x,y,z))

	# Finish GULP file - write potentials. 
	if writeGULP:
		fileGULP.write("\n")
		fileGULP.write("species\n")
		fileGULP.write("Si core Si\n")
		fileGULP.write("O core O\n")
		fileGULP.write("H core H\n")
		fileGULP.write("library reaxff_general.lib\n")

		fileGULPS.write("\n")
		fileGULPS.write("species\n")
		fileGULPS.write("Si core Si\n")
		fileGULPS.write("O core O\n")
		fileGULPS.write("H core H\n")
		fileGULPS.write("library reaxff_general.lib\n")

		fileGULPD.write("\n")
		fileGULPD.write("species\n")
		fileGULPD.write("Si core Si\n")
		fileGULPD.write("O core O\n")
		fileGULPD.write("H core H\n")
		fileGULPD.write("library reaxff_general.lib\n")

	#Write Interface, Deposit, Substarte atoms numbers, and area

	fileIDX.write("%i:%i:%i\n"%(len(iface.IfacePosSC), len(iface.idxDepH),\
			            len(iface.idxSubH)))
	text="Outside idx:   "+"".join("%5i"% t for t in frozenidx)+"\n"
	fileIDX.write(text)
	fileIDX.write("Area: %12.6f\n"%area)
	fileIDX.write("Misfit: %12.6f %12.6f, Angle: %12.6f\n"%(mfit[0],mfit[1],mfit[2]))

	file.close()	
	fileS.close()	
	fileD.close()	
	fileAIMS.close()
	fileSAIMS.close()
	fileDAIMS.close()
	if writeGULP:
		fileGULP.close()
		fileGULPS.close()
		fileGULPD.close()
	fileIDX.close()
			
	if writeGEN:
		print "Generating DFTB+ .gen file"
		# Create DFTB+ GEN file using external xyz2gen tool. 
		# Must be installed and be on the PATH
		
		# Decide how much vacumm to add on top of the interface
		# The amount of vacumm is multiplication of the 
		# interface z-vector
		
		latfile = open("lattice.tmp",'w')
		latfile.write("%12.6f  %12.6f  %12.6f\n"%\
		             (iface.IfaceVecs[0][0],\
			     iface.IfaceVecs[0][1],iface.IfaceVecs[0][2]))
		latfile.write("%12.6f  %12.6f  %12.6f\n"%\
		             (iface.IfaceVecs[1][0],\
			      iface.IfaceVecs[1][1],iface.IfaceVecs[1][2]))
		latfile.write("%12.6f  %12.6f  %12.6f\n"%\
		             (iface.IfaceVecs[2][0]*vac,\
			      iface.IfaceVecs[2][1]*vac,\
			      iface.IfaceVecs[2][2]*vac))
		latfile.close()
		os.system("xyz2gen -l lattice.tmp %s"%filename)
		os.system("xyz2gen -l lattice.tmp %s"%filenameS)
		os.system("xyz2gen -l lattice.tmp %s"%filenameD)

		os.system("rm -f lattice.tmp")

	return

def RecipLattice(a1,a2,a3):

	c23 = np.cross(a2,a3)
	c31 = np.cross(a3,a1)
	c12 = np.cross(a1,a2)

	b1 = 2*m.pi*c23/np.dot(a1,c23)
	b2 = 2*m.pi*c31/np.dot(a1,c23)
	b3 = 2*m.pi*c12/np.dot(a1,c23)

	return [b1,b2,b3]



#########################################
#					#
#					#
#          Start program 		#
#					#
#					#
#########################################

# Set empty global dictionary to hold atom types in both systems
atomTyp = {}
# Read input
subCIF, subMiller, \
depCIF, depMiller, \
maxArea, areaThres, vecThres, angleThres,\
capAtmS, capAtmD, fparam, nLS, nLD, nConf, subAtRad, depAtRad,\
skipStep1, poissonRatio, sandwich, nVac = readInput(inputFile)
print 
print "Substrate CIF..."
i,transM,atoms,positions,atomTyp = ReadCIF(subCIF,atomTyp)
print "Deposit CIF..."
i1,transM1,atoms1,positions1,atomTyp = ReadCIF(depCIF,atomTyp)

subMillerList = createMillerList(subMiller)
depMillerList = createMillerList(depMiller)

# Create list for to store failed results
failedResults = []
# Lists to store when Deposit/Substrate does not have given orientations
depPlaneNE = []
subPlaneNE = []
# Lists to store when looking for primitive vectors failed
depPrimNE  = []
subPrimNE  = []

# Create big bulk of Substarte and Deposit.
# It will be re-used each time an interface needs to be created
print 
print "Constructing bulk strucutre for Substrate..."
subBigBulk = Surface(transM,positions,atoms,atomTyp,np.array((0,0,0)))
subBigBulk.bulkNEW(16)
print "Constructing bulk strucutre for Deposit..."
depBigBulk = Surface(transM1,positions1,atoms1,atomTyp,np.array((0,0,0)))
depBigBulk.bulkNEW(8)

# Iterate over Miller indieces 
iterNo = 0
for subMillerString in subMillerList:
	for depMillerString in depMillerList:
		newIndex = True
		#MATERIAL 1
		print "Constructing bulk Substrate"
		Miller = getMillerFromString(subMillerString)
		# Set size of bulk to be twice the max Miller index
		nbulk = max(Miller)*2

		Sub=Surface(transM,positions,atoms,atomTyp,Miller)
		Sub.bulkNEW(nbulk)

		Sub.construct()
		Sub.plane()
		if not Sub.exists:
			print "Given plane for Substrate does not exists!"
			subPlaneNE.append("%s-%s"%(subMillerString,depMillerString))
			continue
			#exit()
		Sub.initpvecNEW2()
		if not Sub.exists:
			print "Couldn't find primitive vectors for Substrate"
			subPrimNE.append("%s-%s"%(subMillerString,depMillerString))
			continue
			#exit()
		Sub.primitivecell()

		print
		#MATERIAL 2
		print "Constructing bulk Deposit"
		Miller1 = getMillerFromString(depMillerString)
		# Set size of bulk to be twice the max Miller index
		nbulk = max(Miller1)*2

		Dep = Surface(transM1,positions1,atoms1,atomTyp,Miller1)
		Dep.bulkNEW(nbulk)
		Dep.construct()
		Dep.plane()
		if not Dep.exists:
			print "Given plane for Deposit does not exists!"
			depPlaneNE.append("%s-%s"%(subMillerString,depMillerString))
			continue
			#exit()
		Dep.initpvecNEW2()
		if not Dep.exists:
			print "Couldn't find primitive vectors for Deposit"
			depPrimNE.append("%s-%s"%(subMillerString,depMillerString))
			continue
			#exit()
		Dep.primitivecell()

		print 
		print 
		print 
		print "Matching superlattices" 
		print

		# ASSUME THAT P IS FOR DEPOSIT, Q FOR SUBSTRATE
		p,q=AreaMatch2(Dep.primarea,Sub.primarea, maxArea, areaThres)
		print "List of p:"
		print p
		print "List of q:"
		print q

		fitcount = 0 
		totalcomb = 0 
		resultlist=[] 
		vecsDeposit = []
		vecsSubstrate = []

		for i in range(len(p)):
			vecsD = Superlattice(Dep.a,Dep.b,p[i])
			vecsS = Superlattice(Sub.a,Sub.b,q[i])
			for ii in vecsD:
				for jj in vecsS: 
					fit,u1,v1,p1,u2,v2,p2 = SuperlatticeMatch(ii,jj,\
								 vecThres, angleThres)
					if fit:
		#				u1,v1,p1,ar1 = SuperlatticeParams(ii)
		#				u2,v2,p2,ar2 = SuperlatticeParams(jj)
						result=[p[i],q[i],u1,v1,p1,u2,v2,p2]

						if resultlist == []:
							resultlist.append(result)
							vecsSubstrate.append(jj)
							vecsDeposit.append(ii)
						else:
							duplicate = CheckListDuplicates\
								(resultlist,result,\
								 2,3,4,5,6,7)
							if not duplicate:
								resultlist.append(result)
								vecsSubstrate.append(jj)
								vecsDeposit.append(ii)

						fitcount += 1
					totalcomb += 1

		print "TOTAL NUMBER OF COMBINATIONS FOR SET OF p and q PAIRS: %i"%(totalcomb)
		print "TOTAL NUMBER OF FITTED STRUCTURES FOUND: %i"%(fitcount)
		print "NUMBER OF UNIQUE MATCHES FOUND: %i"%(len(resultlist))

		#Didn't find any match
		if len(resultlist) == 0: 
			print 
			print "********************************************************" 
			print "Couldn't find any matching superlattices." 
			print "Materials : %s(%s) - %s(%s) "%(subCIF[0:-4],subMillerString, 
							      depCIF[0:-4],depMillerString) 
			print "Try increasing the max. area or loosening the thresholds" 
			print "********************************************************" 
			print
			failedResults.append("%s-%s"%(subMillerString,depMillerString))
			continue
#			exit()

		print "RESULTS:"
		print "                 Vecx[A]         Vecy[A]    Angle[deg]  Area[A^2]"
		counter = 0
		nResults = len(resultlist)
		misfitList =[]
		for i in resultlist:
			
			ar1 = i[2]*i[3]*m.sin(i[4]*m.pi/180)
			ar2 = i[5]*i[6]*m.sin(i[7]*m.pi/180)

			mu    = (1-(i[2]/i[5]))*100  # Deposit/Substrate
			mv    = (1-(i[3]/i[6]))*100  # Deposit/Substrate
			mang  = (1-(i[4]/i[7]))*100
			marea = (1-(ar1/ar2))*100
			print "Configuration no. %i"%counter
			print "p = %i, q = %i"%(i[0],i[1])
			print "Deposit  : %12.2f   %12.2f   %10.1f   %10.4f"\
					%(i[2],i[3],i[4],ar1)
			print "Substrate: %12.2f   %12.2f   %10.1f   %10.4f"\
					%(i[5],i[6],i[7],ar2)
			print "Misfits  : %11.2f%%   %11.2f%%   %9.2f%%   %9.2f%%"\
					%(mu,mv,mang,marea)
			counter += 1
			misfitList.append([mu,mv,mang,marea])
		
		print "CONSTRUCTING INTERFACE"
		#Construct big planes for Substrate and Deposit
		#nbulk2 = 8
		#Sub.bulk(nbulk2)
		print 
		Sub.positions = subBigBulk.positions.copy()
		Sub.atoms = subBigBulk.atoms.copy()
	        Sub.positionsSmall = subBigBulk.positionsSmall.copy()
      	        Sub.atomsSmall = subBigBulk.atomsSmall.copy()
		Sub.construct()
		Sub.plane()


		#ii=0
		#for i in Sub.planepos:
		#	print Sub.planeatms[ii],i[0],i[1],i[2]
		#	ii+=1


#		Dep.bulk(nbulk2)
		Dep.positions = depBigBulk.positions.copy()
		Dep.atoms = depBigBulk.atoms
	        Dep.positionsSmall = depBigBulk.positionsSmall.copy()
      	        Dep.atomsSmall = depBigBulk.atomsSmall.copy()
		Dep.construct()
		Dep.plane()

		#ii=0
		#for i in Dep.planeposblk:
		#	print Dep.planeatmsblk[ii],i[0],i[1],i[2]
		#	ii+=1


		# Sort the Substrate-Deposit alignemnt by checking coolinarity of
		# pairs of anticipatory vectors.
		avecsSub = Sub.avecs.copy()
		avecsDep = Dep.avecs.copy()
		ii = 0
		# If there is only ona pair of anticipatory vectors,
		# it means those are vectors pointing straight up, eg. [0,0,z].
		# Indicate it, by giving them negative labels.
		if (len(avecsDep) != 1) or (len(avecsSub) != 1):
			resmat = np.zeros(((len(avecsDep)*len(avecsSub)),3))
			for i in range(len(avecsSub)):
				for j in range(len(avecsDep)):
					prod = np.dot(avecsSub[i],avecsDep[j])
					nS = np.linalg.norm(avecsSub[i])
					nD = np.linalg.norm(avecsDep[j])
					cosphi = prod/(nS*nD)
					cosphi = round(cosphi,4)
					resmat[ii,0] = m.acos(cosphi)
					resmat[ii,1] = i  
					resmat[ii,2] = j 
					ii += 1
			resmat = CleanMatElements(resmat)
			# Sort the results and pick the smallest dot product
			resmat = resmat[np.argsort(resmat[:,0])]
			# Create additional alignment - where vectors are pointing straight up
			resmat = np.vstack([np.array((0.0,-1.0,-1.0)),resmat])

		else:
			resmat = np.array([[0.0,-1.0,-1.0]])

		print
		print "****************************"
		print "Found %i possible alignments"%len(resmat) 
		print "****************************"
		print

		#bondl=[0.9,1.6,1.9,1.4,2.0]
		#resmat=np.array([[0.0,-1.0,-1.0]])
		#resmat=np.array([[0.0,1.0,1.0]])
		#resmat=np.array([resmat[2]])
		scoreDirName = "SCORE-%s-%s"%(subMillerString,depMillerString)
		os.system("mkdir -p %s/DIST"%scoreDirName)

		# if user asked for more results than there actually is, limit it.
		if nConf > nResults: nConf = nResults

		# Initialize variables
		vecpairidx = 0   # counter for anticicipatory vector pairs

		confnor = range(nConf) # number of configurations to consider
		#fparam = 1         # f-parameter for scoring function
		genHD = False       # generate hydrogens on Deposit?
		genHS = False      # generate hydrogens on Substrate?
		#capAtmS = "Cl"
		#capAtmD = "H"
		atomicRadius = (subAtRad + depAtRad)/2 

		bond = np.arange(atomicRadius-0.5,atomicRadius+1.0,0.1) # range of bond lengths to scan for search of 
					      # optimal Substrate-Deposit separation
#		bond = np.arange(0.5,3.0,0.1) 

		wflist = [fparam]  # list of f-parameters for Scoring function
		#skipStep1 = True # option to skip Step 1. You must specify bondlist 
		#                   # by hand below

		bondlist = []      # list for Substrate-Deposit optimal sepration length
		if skipStep1:
			for i in resmat:
				bondlist.append(atomicRadius*2)
#			bondlist = [0.9,2.2,1.3,1.4,2.2]

		# Start generating structures
		# Step 1: Optimal distance between Substrate and Deposit 
		# Find the optimal distance between Substrate and Deposit using 
		# scoring function. We do this by scanning the Substrate and Deposit 
		# distance of first configurations for every anticipatory vectors pair 
		# for the range of distances given in "bond" variable. We choose the 
		# distance corresponding to maximal score and save it in "bondlist"

		if not skipStep1:
			print 
			print "FLAG skipStep1 IS NOT SET."
			print "USING SCORING FUNCTION TO DETERMINE OPTIMAL DISTANCE"
			print "BETWEEN SLABS. THIS MIGHT TAKE A WHILE."
			print
			for vecpair in resmat:
				print
				print "***************************************"
				print "Working on Alignment no. %i out of %i"\
					%(vecpairidx+1,len(resmat))
				print "***************************************"
				print
				bscoreS=[]
				bscoreD=[]
				bscoreSD=[]
				fS=open("%s/DIST/OUTPUT-S-%i.txt"%(scoreDirName,vecpairidx),'w')
				fD=open("%s/DIST/OUTPUT-D-%i.txt"%(scoreDirName,vecpairidx),'w')
				fSD=open("%s/DIST/OUTPUT-SD-%i.txt"%(scoreDirName,vecpairidx),'w')
				for i in bond:
					iface = Interface(vecsDeposit[0],vecsSubstrate[0],\
						     Dep,Sub,vecpair,i,fparam,atomicRadius,False,\
						     False,False,False,False,False,False)
					s1=iface.sc1
					s2=iface.sc2	
					sca=(s1+s2)/2
			
					bscoreD.append(s1)
					bscoreS.append(s2)
					bscoreSD.append(sca)
					# Write to score to files
					fD.write("%8.4f   %12.6f\n"%(i,s1))
					fS.write("%8.4f   %12.6f\n"%(i,s2))
					fSD.write("%8.4f   %12.6f\n"%(i,sca))

				#Find max score for given configuration
				maxSD = max(bscoreSD)
				scoreidx = bscoreSD.index(maxSD)
				minblength = bond[scoreidx]
				bondlist.append(minblength)
				fD.close()
				fS.close()
				fSD.close()

				vecpairidx += 1
		print "BOND LIST",bondlist
		# Step 2: Generte configurations
		# For every anticipatory vector pair in "resmat" generate number of structures
		# given in "confnor", setting the spacing between them to one from "bondlist" 
		# found in step 1.

		vecpairidx = 0 
		for vecpair in resmat:
			aligNo = vecpairidx
			print
			print "***************************************"
			print "Working on Alignment no. %i out of %i"%(aligNo+1, len(resmat))
			print "***************************************"
			print
			fS=open("%s/OUTPUT-S-%i.txt"%(scoreDirName,vecpairidx),'w')
			fD=open("%s/OUTPUT-D-%i.txt"%(scoreDirName,vecpairidx),'w')
			fSD=open("%s/OUTPUT-SD-%i.txt"%(scoreDirName,vecpairidx),'w')
			scresultsS=[]
			scresultsD=[]
			scresultsSD=[]

			for wf in wflist: # Old option to be able to scan the entire space of 
					  # f-parameters in scoring function in one run.
				resultstmpS=[]
				resultstmpD=[]
				resultstmpSD=[]
				for confno in confnor:
					iterNo += 1 # Mark iteration no to be used in mkOutput to create summary.txt file
					print
					print "Configuration no. %5i"%confno

					i=bondlist[vecpairidx]
					iface = Interface(vecsDeposit[confno],\
							vecsSubstrate[confno],Dep,Sub,\
							vecpair,i,wf,atomicRadius, nLS, nLD,\
							capAtmD, capAtmS,sandwich,\
							genHD,genHS,poissonRatio)
					s1=iface.sc1
					s2=iface.sc2	
					sca=(s1+s2)/2
			
					resultstmpD.append(s1)
					resultstmpS.append(s2)
					resultstmpSD.append(sca)

					i=confno
					mfit = misfitList[i]
					mkOutput(iface,i,depCIF[0:-4],depMillerString,\
							subCIF[0:-4],subMillerString,\
							vecpairidx,nLS,nLD,mfit,sca,\
							aligNo,iterNo,nVac,sandwich,\
							writeGEN=False,\
							writeAIMS=True,\
							writeGULP=False,\
							cstrlxD=3,cstrlxS=2,\
							cstrlxH=False,\
							bondmod=False)
				scresultsS.append(resultstmpS)
				scresultsD.append(resultstmpD)
				scresultsSD.append(resultstmpSD)
			vecpairidx += 1

			#Prepeare OUTPUT HEADER
			text = "".join("%10i" %i for i in wflist)
			fS.write("     %s\n"%text)
			fD.write("     %s\n"%text)
			fSD.write("     %s\n"%text)
			#Write to file
			########
			##WHEN CHANGING BOND LENGTH
			#########
		#	for i in range(len(bond)):
		#		# Subsrtate
		#		text = ""
		#		text = "".join("%10.3f" % t[i] for t in scresultsS)
		#		fS.write("%2.1f  %s\n"%(bond[i],text))
		#		# Deposit
		#		text = ""
		#		text = "".join("%10.3f" % t[i] for t in scresultsD)
		#		fD.write("%2.1f  %s\n"%(bond[i],text))
				# Deposit+Substrate
		#		text = ""
		#		text = "".join("%10.3f" % t[i] for t in scresultsSD)
		#		fSD.write("%2.1f  %s\n"%(bond[i],text))
		#
			########
			##WHEN CHANGING CONFIGURATION
			#########
			for i in range(len(confnor)):
				# Subsrtate
				text = ""
				text = "".join("%10.3f" % t[i] for t in scresultsS)
				fS.write("%2.1f  %s\n"%(confnor[i],text))
				# Deposit
				text = ""
				text = "".join("%10.3f" % t[i] for t in scresultsD)
				fD.write("%2.1f  %s\n"%(confnor[i],text))
				# Deposit+Substrate
				text = ""
				text = "".join("%10.3f" % t[i] for t in scresultsSD)
				fSD.write("%2.1f  %s\n"%(confnor[i],text))

			fS.close()
			fD.close()
			fSD.close()
#List failed results
fileF = open("FAILED_RESULTS.txt",'w')
if len(failedResults) >0:
	fileF.write("The following orientations coudn't be found:\n")
	fileF.write("Substrate: %s  -  Deposit: %s\n"%(subCIF[0:-4],depCIF[0:-4]))
	for i in failedResults:
		fileF.write("%s\n"%i)

	fileF.write("\nTry increasing max area, or loosening the thresholds\n\n")

if len(subPlaneNE) >0:
	fileF.write("Following planes does not exist for the Substrate:\n")
	fileF.write("Substrate: %s\n"%subCIF[0:-4])
	for i in subPlaneNE:
		fileF.write("%s\n"%i)
	fileF.write("\n\n")

if len(depPlaneNE) >0:
	fileF.write("Following planes does not exist for the Deposit:\n")
	fileF.write("Substrate: %s\n"%depCIF[0:-4])
	for i in depPlaneNE:
		fileF.write("%s\n"%i)
	fileF.write("\n\n")

if len(subPrimNE) >0:
	fileF.write("Couldn't find primitive vectors for Substrate planes:\n")
	fileF.write("Substrate: %s\n"%subCIF[0:-4])
	for i in subPrimNE:
		fileF.write("%s\n"%i)
	fileF.write("\n\n")

if len(depPrimNE) >0:
	fileF.write("Couldn't find primitive vectors for Deposit planes:\n")
	fileF.write("Substrate: %s\n"%depCIF[0:-4])
	for i in depPrimNE:
		fileF.write("%s\n"%i)
	fileF.write("\n\n")

fileF.close()	


for i in range(len(misfitList)):
	n = ((misfitList[i][0]**2) + (misfitList[i][1]**2))**0.5
	print "%i   %12.6f"%(i,n)

