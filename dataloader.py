import os
import sys
import glob
import numpy as np
import pandas as pd
from utils import util
from urllib.request import urlretrieve
import requests
import validators

def load_raw_data(intermetaldb_csv):
    assert os.path.exists(intermetaldb_csv), f"File not found: {intermetaldb_csv}"
    df = pd.read_csv(intermetaldb_csv)
    #modifying raw data
    df = df[(df['Is homomer'] == 1) & ~(df['Amino acids or nucleotide residues names'].str.contains('DG')) & ~(df['Amino acids or nucleotide residues names'].str.contains('DT')) & ~(df['Amino acids or nucleotide residues names'].str.contains('DC')) & ~(df['Amino acids or nucleotide residues names'].str.contains('DA'))] #exclude non-protein macromolecule data
    df['pdb_id'] = df['Id'].str.split('-', expand=True)[0]
    df['metal'] = df['Id'].str.split('-', expand=True)[1]
    df['metal_binding_id'] = df['Id'].str.split('-', expand=True)[2]
    return df

def retrieve_dataset(df, rcsb_link='https://files.rcsb.org/download/', log_dir ='./', dataset_dir='./datasets/rcsb_pdbs/'):
    logs = {'pdb_id':[], 'status':[]}
    for pdb_id in (df['pdb_id'].unique()):
        logs['pdb_id'].append(pdb_id)
        url = f'{rcsb_link}{pdb_id}.pdb'
        file = f'{dataset_dir}{pdb_id}.pdb'
        r = requests.get(url)
        if '404' in str(r):
            logs['status'].append(0) #failed
            continue
        urlretrieve(url,file)
        logs['status'].append(1) #success
    logs = pd.DataFrame(logs)
    logs.to_csv(f'{log_dir}/retrieve_dataset.log', index=False)
    return None

def analyze_dataset(df, logs, out_dir='./'):
    logs = pd.read_csv('datasets.log')
    pdb_ok = logs[logs['status']==1]['pdb_id']
    df_filter = df[df['pdb_id'].isin(pdb_ok)]
    df_filter.to_csv(f'{out_dir}filtered_intermetaldb.csv', index=False)
    return None


def serch_metal_site(df):
    for pdb_id in df['pdb_id'].unique():
        if pdb_id != id_debug and debug:
            continue
        sample_df = df[df['pdb_id'] == pdb_id]
        row = sample_df[sample_df['Representative'] == 1].head(1)
        if len(row) == 0:
            row = sample_df.sort_values(by='metal_binding_id').head(1)
        hotspots_raw = row['Amino acids or nucleotide residues names'].values[0].split('.')
        hotspots_raw.remove('')
        hotspots = []
        for h in hotspots_raw:
            if h in util.THREE_TO_ONE.keys():
                hotspots.append(h)
        metal_name = row['metal'].values[0]
        metal_binding_id = row['metal_binding_id'].values[0]
        no_bound_residues = row['Number of bound amino acids or nucleotide residues'].values[0]
        motif_residue = []

        if hotspots==[]:
            continue #sample more cause oversampling, or non-canonical resi
        metal_infos = util.get_ligand_num(f'datasets/rcsb_pdbs/{pdb_id}.pdb', metal_name, multi=True)
        for metal_info in metal_infos:
            metal_num, metal_chain, metal_name, metal_xyz = metal_info
            motif = util.get_ligand_contacts(f'datasets/rcsb_pdbs/{pdb_id}.pdb', hetatm_num=metal_num, hetatm=metal_name,hetatm_chain=metal_chain, cutoff=6.0,hotspot=hotspots) # motif_resi_num, motif_chain, motif_resi_id
            motif_residue = [item[-1] for item in motif]
            if len(motif) >=  len(hotspots):# and :
                if set(hotspots) == set(motif_residue):
                    break
                elif len(set(motif_residue)) == len(hotspots):
                    break
    return metal_info, motif


def extract_metal_motif(pdb_id, metal_info, motif, dataset_dir='./datasets/rcsb_pdbs/', motif_dir='./datasets/rotamer_motif/'):
    metal_num, metal_chain, metal_name, metal_xyz = metal_info
    with open(f'{dataset_dir}{pdb_id}.pdb', 'r') as f:
        w_lines = []
        for line in f.readlines():
            if line.startswith('HETATM'):
                if line[17:21].strip() == metal_name and line[22:26].strip() == metal_num and line[21:22].strip() == metal_chain:
                    w_lines.append(line[:30] + '{:8.3f}'.format(0.0) + '{:8.3f}'.format(0.0) + '{:8.3f}'.format(0.0) + line[54:])
            if line.startswith('ATOM'):
                for m in motif:
                    if line[22:26].strip()== m[0] and line[21:22].strip() == m[1] and line[17:21].strip()== m[2]:
                        line_xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                        w_lines.append(line[:30] + '{:8.3f}'.format(line_xyz[0]-metal_xyz[0]) + '{:8.3f}'.format(line_xyz[1]-metal_xyz[1]) + '{:8.3f}'.format(line_xyz[2]-metal_xyz[2]) + line[54:])

    with open(f'{motif_dir}{pdb_id}_{metal_binding_id}.pdb', 'w') as f1:
        f1.writelines(w_lines)
