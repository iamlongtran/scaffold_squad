import os,glob,sys
import parsers


def pdb_hetatmlines(file, silent=False): #check to see if pdb has hetatm lines
	f = open(file, 'r')
	lines = f.readlines()
	lines = [line for line in lines if 'HETATM' in line]
	l = len(lines)
	if l == 0:
		if silent == False:
			raise ValueError("File does not contain HETATM lines")
		else:
			return False
	else:
		return True
    

def pdb_containsmetal(file, metal, silent=False):
	f = open(file, 'r')
	lines = f.readlines()
	lines = [line for line in lines if 'HETATM' in line]
	n = 0
	for line in lines:
		if metal in line:
			n+=1
	if n == 0:
		if silent == False:
			raise ValueError(f"Atom {metal} not present in pdb file")
		else:
			return False
	else:
		return True



def pdb_has_metalsite(pdb, metal_id, atom_type, silent=False): 
	parsed = parsers.parse_pdb(pdb)
	parsed = parsed['info_het']
	n = 0
	for i in parsed:
		idx = i['idx']
		if int(idx) == int(metal_id):
			if atom_type in i['name']:
				n+=1
	if n == 0:
		if silent == False:
			raise ValueError(f"Metal atom {atom_type} {metal_id} not in pdb")
		else:
			return False
	
	else:
		return True

def check_inputs(pdb, metal_id, atom_type, silent=False):
	if silent == False:
		assert pdb_hetatmlines(pdb) == True
		assert pdb_has_metalsite(pdb, metal_id, atom_type) == True
		print("Inputs accepted!")
	else:
		het = pdb_hetatmlines(pdb, silent=True)
		site = pdb_has_metalsite(pdb, metal_id, atom_type, silent=True)
		if (het == False or site == False):
			return False
		else:
			return True
