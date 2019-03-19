#!/usr/bin/python

input = "cctbx2dat.in"

file = open(input,'r')
counter = 1
while 1:
	line = file.readline()
	if line == "":break
	line = line.split()
	new=""
	for i in line:
		new = new + i
	new = new[:-2]
	print "%5i   %s"%(counter,new)
	counter += 1
	

