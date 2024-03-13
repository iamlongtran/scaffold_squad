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


test_pdb = '/home/lhtran/class/scaffold_squad/datasets/rotamer_motif/7MCI_1.pdb'
pdb_data = util.preprocess_tip_atom(test_pdb, metal='CA')#, metal_id=15859))
print(pdb_data)
#print(util.get_motif(test_pdb, metal='CA'))#, metal_id=15859)

