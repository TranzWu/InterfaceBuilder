#!/usr/bin/python

import json

input = open('spacegroups.json','r')

data = json.load(input)

input.close()

for i in  data["Fm-3m"]:
	print i
