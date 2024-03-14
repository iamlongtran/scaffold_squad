import argparser
import os,sys,glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
script_dir = os.path.abspath('') #only in jupyter demo, different in .py
parent_dir = os.path.join(os.path.dirname(script_dir), '')
sys.path.append(parent_dir) #only in jupyter demo, different in .py
from utils import util
from tqdm import tqdm
import joblib
from utils import unit_tests

parser = ArgumentParser()

parser.add_arg('--pdb', help='pdb file with metal motif for prediction')
parser.add_arg('--metal_id', help='atom index number of metal in pdb file')
parser.add_arg('--metal_type', help='type of metal for binding prediction (e.g ZN)')
parser.add_arg('--model_wt', help='path to model weight', default=None)

args = parser.parse_args()

def main(args):
	
	pdb = args.pdb
	idx = args.metal_id
	atom_type = args.metal_type
	model = args.model

	unit_tests.check_inputs(pdb, idx, atom_type) #verify pdb is in correct format

	name = pdb.split('/')[-1]
	name = name.split('.pdb')[0]

	if model == None:
		if atom_type == 'CA':
			model = joblib.load(f'{parent_dir}/weight3/one_class_svm_CA_nu_0.1_roc_0.72.pkl')
		if atom_type == 'NA':
			model = joblib.load(f'{parent_dir}/weight3/')
		if atom_type == 'K':
			model = joblib.load(f'{parent_dir}/weight3/')
		if atom_type == 'ZN':
			model = joblib.load(f'{parent_dir}/weight3/')
	else:
		model = joblib.load(model)

	pdb_data = util.preprocess_tip_atom(pdb, metal=atom_type)
	pdb_data = np.reshape(pdb_data, (1,-1))
	pred = model.predict(pdb_data)
	X = ['no result\n']
	if int(pred[0]) == 1:
		print(f'Scaffold predicted to bind {atom_type}!')
		X = [f'prediction for {atom_type} binding of scaffold {name}: 1']
	if int(pred[0]) == -1:
		print(f'Scaffold predicted NOT TO bind {atom_type}!')
		X = [f'prediction for {atom_type} binding of scaffold {name}:-1']
	f = open(f'{name}_{atom_type}_prediction.txt', 'w')
	f.writelines(X)
	f.close()

if __name__ == "__main__":
	main(args)



