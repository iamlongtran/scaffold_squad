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

parser.add_arg('--csv', help='csv file with pdbs for prediction, accompanying metal atom index numbers, and atom type. Columns should be labeled "pdb", "metal_idx", and "metal"')
parser.add_arg('--metal_type', help='type of metal for binding prediction (e.g ZN)')
parser.add_arg('--model_ZN', help='path to model weight for ZN', default=None)
parser.add_arg('--model_CA', help='path to model weight for CA', default=None)
parser.add_arg('--model_NA', help='path to model weight for NA', default=None)
parser.add_arg('--model_K', help='path to model weight for K', default=None)
parser.add_arg('--model_MG', help='path to model weight for MG', default=None)
parser.add_arg('--output', help='name of results csv file to output')

args = parser.parse_args()

def main(args):
	
	csv = args.csv
	atom_type = args.metal_type
	model_ZN = args.model_ZN
	model_CA = args.model_CA
	model_NA = args.model_NA
	model_K = args.model_K
	model_MG = args.model_MG
	output = args.output

	unit_tests.check_inputs(pdb, idx, atom_type) #verify pdb is in correct format

	name = pdb.split('/')[-1]
	name = name.split('.pdb')[0]

	if model_CA == None:	
		model_CA = joblib.load(f'{parent_dir}/weight4/svm_svc_CA_roc_0.81.pkl')
	else:
		model_CA = joblib.load(model_CA)
	
	if model_NA == None:
		model_NA = joblib.load(f'{parent_dir}/weight4/svm_svc_NA_roc_0.871.pkl')
	else:
		model_NA = joblib.load(model_NA)
	
	if model_K == None:
		model_K = joblib.load(f'{parent_dir}/weight4/svm_svc_K_roc_0.929.pkl')
	else:
		model_K = joblib.load(model_K)
		
	if model_ZN == None:
		model_ZN = joblib.load(f'{parent_dir}/weight4/svm_svc_ZN_roc_0.873.pkl')
	else:
		model_ZN = joblib.load(model_ZN)

	if model_MG == None:
		model_MG = joblib.load(f'{parent_dir}/weight4/svm_svc_MG_roc_0.962.pkl')
	else:
		model_MG = joblib.load(model_MG)
	



	df = pd.read_csv(csv)
	l = df.shape[0]
	results = []

	for i in range(l):

		pdb = df.iloc[i]['pdb']
		idx = df.iloc[i]['metal_idx']
		metal = df.iloc[i]['metal']
		tf = unit_tests.check_inputs(pdb, idx, metal, silent=True)
		if tf == True:


			pdb_data = util.preprocess_tip_atom(pdb, metal=metal)
			pdb_data = np.reshape(pdb_data, (1,-1))
			if metal == 'ZN':
				pred = model_ZN.predict(pdb_data)
				pred = int(pred[0])
				results.append([pdb, idx, metal, pred])
			if metal == 'CA':
				pred = model_CA.predict(pdb_data)
				pred = int(pred[0])
				results.append([pdb, idx, metal, pred])
			if metal == 'K':
				pred = model_K.predict(pdb_data)
				pred = int(pred[0])
				results.append([pdb, idx, metal, pred])
			if metal == 'NA':
				pred = model_NA.predict(pdb_data)
				pred = int(pred[0])
				results.append([pdb, idx, metal, pred])
			if metal == 'MG':
				pred = model_MG.predict(pdb_data)
				pred = int(pred[0])
				results.append([pdb, idx, metal, pred])

		
		else:
			results.append([pdb, idx, metal, 'INPUT ERROR! CHECK YOUR PDB!'])

		dfout = pandas.DataFrame(results, columns=['pdb', 'metal_idx', 'metal', 'prediction'])
		dfout.to_csv(output)



if __name__ == "__main__":
	main(args)