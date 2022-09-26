import os
import glob


examples = ['linear_elasticity']
direction = '../examples'


material_file = []
material_file.append('.. Copyright (c) 2022, Peter Lenz\n\n')
material_file.append('   Distributed under the terms of the BSD 3-Clause License.\n\n')
material_file.append('   The full license is in the file LICENSE, distributed with this software.\n\n\n')
material_file.append('=========\n')
material_file.append('Materials\n')
material_file.append('=========\n')
material_file.append('\n\n')

for name in next(os.walk(direction))[1]:
	print(name)
	Name = name.replace("_", " ")
	material_file.append(Name.title() + '\n')
	material_file.append(len(Name)*'='+ '\n')
	material_file.append('\n')

	for example in next(os.walk(direction + '/' + name))[1]:
		print(example)
		Example = example.replace("_", " ")
		material_file.append(Example.title() + '\n')
		material_file.append(len(Example)*'-'+ '\n')
		material_file.append('\n')
		main_file = []
		main_path = direction + "/" + name + "/" + example + "/" + "main.cpp"
		
		if os.path.exists(main_path):
			print(main_path)
			with open(main_path, 'r') as f:
				main_file = f.readlines()
			material_file.append(".. code::\n\n")
			for line in main_file:
				material_file.append("    " + line)
			material_file.append("\n\n")
				


with open("../docs/source/api/materials.rst", 'w') as file:
	file.writelines(material_file)
