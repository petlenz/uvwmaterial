import os
import glob

directions = ['../include/uvwmaterial']

for direction in directions:
	list_of_files = glob.glob(direction + '/*.h')
	for name in list_of_files:
		print(name)
		data_file = []
		with open(name, 'r') as file:
			data_file = file.readlines()
		data_out = []
		data_out.append('/***************************************************************************\n')
		data_out.append('* Copyright (c) Peter Lenz                                                 *\n')
		data_out.append('*                                                                          *\n')
		data_out.append('* Distributed under the terms of the BSD 3-Clause License.                 *\n')
		data_out.append('*                                                                          *\n')
		data_out.append('* The full license is in the file LICENSE, distributed with this software. *\n')
		data_out.append('****************************************************************************/\n')
		for line in data_file:
			data_out.append(line)
		with open(name, 'w') as file:
			file.writelines(data_out)
