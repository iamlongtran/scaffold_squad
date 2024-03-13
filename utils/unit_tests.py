import os,glob,sys
import parsers


def pdb_hetatmlines(file): #check to see if pdb has hetatm lines
	f = open(file, 'r')
	lines = f.readlines()
	lines = [line for line in lines if 'HETATM' in line]
	l = len(lines)
	if l == 0:
		raise ValueError("File does not contain HETATM lines")
	else:
		pass
	if l > 1:
	    raise ValueError("File cannot have more than one ligand")
	else:
		pass


def pdb_containsmetal(file, metal):
	f = open(file, 'r')
	lines = f.readlines()
	lines = [line for line in lines if 'HETATM' in line]
	assert len(lines) == 1
	line = lines[0]
	if not metal in line:
		raise ValueError(f"File {pdb} does not contain metal atom {metal}")
	else:
		pass


def pdb_has_metalsite(pdb, metal_id, atom_type): 
	parsed = parsers.parse_pdb(pdb)
	parsed = parsed['info_het']
	n = 0
	for i in parsed:
		idx = i['idx']
		if int(idx) == int(metal_id):
			if atom_type in i['atom_type']:
				n+=1
	if n == 0:
		raise ValueError(f"Metal atom {atom_type} {metal_id} not in pdb")
	else:
		pass
	if n > 1:
		raise ValueError(f"Error: PDB contains duplicates of metal atom!")
	else:
		pass
	if n == 1:
		return True

def check_inputs(pdb, metal_id, atom_type):
	assert pdb_hetatmlines(pdb) == True
	assert pdb_containsmetal(pdb, atom_type) == True
	assert pdb_has_metalsite(pdb, metal_id, atom_type) == True
	print("Inputs accepted!")
