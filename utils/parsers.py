import numpy as np
import scipy
import scipy.spatial
import string
import sys,os,re
import random
script_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(script_dir)
import util
#from . import util


num2aa=[
    'ALA','ARG','ASN','ASP','CYS',
    'GLN','GLU','GLY','HIS','ILE',
    'LEU','LYS','MET','PHE','PRO',
    'SER','THR','TRP','TYR','VAL',
    'UNK','MAS',
    ' DA',' DC',' DG',' DT', ' DX',
    ' RA',' RC',' RG',' RU', ' RX',
    'HIS_D', # only used for cart_bonded
    'Al', 'As', 'Au', 'B',
    'Be', 'Br', 'C', 'Ca', 'Cl',
    'Co', 'Cr', 'Cu', 'F', 'Fe',
    'Hg', 'I', 'Ir', 'K', 'Li', 'Mg',
    'Mn', 'Mo', 'N', 'Ni', 'O',
    'Os', 'P', 'Pb', 'Pd', 'Pr',
    'Pt', 'Re', 'Rh', 'Ru', 'S',
    'Sb', 'Se', 'Si', 'Sn', 'Tb',
    'Te', 'U', 'W', 'V', 'Y', 'Zn',
    'ATM'
]
aa2num= {x:i for i,x in enumerate(num2aa)}
aa2num['MEN'] = 20
aa2long=[
    (" N  "," CA "," C  "," O  "," CB ",  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ","3HB ",  None,  None,  None,  None,  None,  None,  None,  None), #0  ala
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD "," NE "," CZ "," NH1"," NH2",  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ","1HG ","2HG ","1HD ","2HD "," HE ","1HH1","2HH1","1HH2","2HH2"), #1  arg
    (" N  "," CA "," C  "," O  "," CB "," CG "," OD1"," ND2",  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ","1HD2","2HD2",  None,  None,  None,  None,  None,  None,  None), #2  asn
    (" N  "," CA "," C  "," O  "," CB "," CG "," OD1"," OD2",  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ",  None,  None,  None,  None,  None,  None,  None,  None,  None), #3  asp
    (" N  "," CA "," C  "," O  "," CB "," SG ",  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB "," HG ",  None,  None,  None,  None,  None,  None,  None,  None), #4  cys
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD "," OE1"," NE2",  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ","1HG ","2HG ","1HE2","2HE2",  None,  None,  None,  None,  None), #5  gln
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD "," OE1"," OE2",  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ","1HG ","2HG ",  None,  None,  None,  None,  None,  None,  None), #6  glu
    (" N  "," CA "," C  "," O  ",  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None," H  ","1HA ","2HA ",  None,  None,  None,  None,  None,  None,  None,  None,  None,  None), #7  gly
    (" N  "," CA "," C  "," O  "," CB "," CG "," ND1"," CD2"," CE1"," NE2",  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ","2HD ","1HE ","2HE ",  None,  None,  None,  None,  None,  None), #8  his
    (" N  "," CA "," C  "," O  "," CB "," CG1"," CG2"," CD1",  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None," H  "," HA "," HB ","1HG2","2HG2","3HG2","1HG1","2HG1","1HD1","2HD1","3HD1",  None,  None), #9  ile
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD1"," CD2",  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB "," HG ","1HD1","2HD1","3HD1","1HD2","2HD2","3HD2",  None,  None), #10 leu
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD "," CE "," NZ ",  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ","1HG ","2HG ","1HD ","2HD ","1HE ","2HE ","1HZ ","2HZ ","3HZ "), #11 lys
    (" N  "," CA "," C  "," O  "," CB "," CG "," SD "," CE ",  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ","1HG ","2HG ","1HE ","2HE ","3HE ",  None,  None,  None,  None), #12 met
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD1"," CD2"," CE1"," CE2"," CZ ",  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ","1HD ","2HD ","1HE ","2HE "," HZ ",  None,  None,  None,  None), #13 phe
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD ",  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None," HA ","1HB ","2HB ","1HG ","2HG ","1HD ","2HD ",  None,  None,  None,  None,  None,  None), #14 pro
    (" N  "," CA "," C  "," O  "," CB "," OG ",  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None," H  "," HG "," HA ","1HB ","2HB ",  None,  None,  None,  None,  None,  None,  None,  None), #15 ser
    (" N  "," CA "," C  "," O  "," CB "," OG1"," CG2",  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None," H  "," HG1"," HA "," HB ","1HG2","2HG2","3HG2",  None,  None,  None,  None,  None,  None), #16 thr
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD1"," CD2"," NE1"," CE2"," CE3"," CZ2"," CZ3"," CH2",  None,  None,  None,  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ","1HD ","1HE "," HZ2"," HH2"," HZ3"," HE3",  None,  None,  None), #17 trp
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD1"," CD2"," CE1"," CE2"," CZ "," OH ",  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ","1HD ","1HE ","2HE ","2HD "," HH ",  None,  None,  None,  None), #18 tyr
    (" N  "," CA "," C  "," O  "," CB "," CG1"," CG2",  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None," H  "," HA "," HB ","1HG1","2HG1","3HG1","1HG2","2HG2","3HG2",  None,  None,  None,  None), #19 val
    (" N  "," CA "," C  "," O  "," CB ",  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ","3HB ",  None,  None,  None,  None,  None,  None,  None,  None), #20 unk
    (" N  "," CA "," C  "," O  "," CB ",  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ","3HB ",  None,  None,  None,  None,  None,  None,  None,  None), #21 mask
    (" OP1"," P  "," OP2"," O5'"," C5'"," C4'"," O4'"," C3'"," O3'"," C2'"," C1'"," N9 "," C4 "," N3 "," C2 "," N1 "," C6 "," C5 "," N7 "," C8 "," N6 ",  None,  None,"H5''"," H5'"," H4'"," H3'","H2''"," H2'"," H1'"," H2 "," H61"," H62"," H8 ",  None,  None), #22  DA
    (" OP1"," P  "," OP2"," O5'"," C5'"," C4'"," O4'"," C3'"," O3'"," C2'"," C1'"," N1 "," C2 "," O2 "," N3 "," C4 "," N4 "," C5 "," C6 ",  None,  None,  None,  None,"H5''"," H5'"," H4'"," H3'","H2''"," H2'"," H1'"," H42"," H41"," H5 "," H6 ",  None,  None), #23  DC
    (" OP1"," P  "," OP2"," O5'"," C5'"," C4'"," O4'"," C3'"," O3'"," C2'"," C1'"," N9 "," C4 "," N3 "," C2 "," N1 "," C6 "," C5 "," N7 "," C8 "," N2 "," O6 ",  None,"H5''"," H5'"," H4'"," H3'","H2''"," H2'"," H1'"," H1 "," H22"," H21"," H8 ",  None,  None), #24  DG
    (" OP1"," P  "," OP2"," O5'"," C5'"," C4'"," O4'"," C3'"," O3'"," C2'"," C1'"," N1 "," C2 "," O2 "," N3 "," C4 "," O4 "," C5 "," C7 "," C6 ",  None,  None,  None,"H5''"," H5'"," H4'"," H3'","H2''"," H2'"," H1'"," H3 "," H71"," H72"," H73"," H6 ",  None), #25  DT
    (" OP1"," P  "," OP2"," O5'"," C5'"," C4'"," O4'"," C3'"," O3'"," C2'"," C1'",  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,"H5''"," H5'"," H4'"," H3'","H2''"," H2'"," H1'",  None,  None,  None,  None,  None,  None), #26  DX (unk DNA)
    (" OP1"," P  "," OP2"," O5'"," C5'"," C4'"," O4'"," C3'"," O3'"," C1'"," C2'"," O2'"," N1 "," C2 "," N3 "," C4 "," C5 "," C6 "," N6 "," N7 "," C8 "," N9 ",  None," H5'","H5''"," H4'"," H3'"," H2'","HO2'"," H1'"," H2 "," H61"," H62"," H8 ",  None,  None), #27   A
    (" OP1"," P  "," OP2"," O5'"," C5'"," C4'"," O4'"," C3'"," O3'"," C1'"," C2'"," O2'"," N1 "," C2 "," O2 "," N3 "," C4 "," N4 "," C5 "," C6 ",  None,  None,  None," H5'","H5''"," H4'"," H3'"," H2'","HO2'"," H1'"," H42"," H41"," H5 "," H6 ",  None,  None), #28   C
    (" OP1"," P  "," OP2"," O5'"," C5'"," C4'"," O4'"," C3'"," O3'"," C1'"," C2'"," O2'"," N1 "," C2 "," N2 "," N3 "," C4 "," C5 "," C6 "," O6 "," N7 "," C8 "," N9 "," H5'","H5''"," H4'"," H3'"," H2'","HO2'"," H1'"," H1 "," H22"," H21"," H8 ",  None,  None), #29   G
    (" OP1"," P  "," OP2"," O5'"," C5'"," C4'"," O4'"," C3'"," O3'"," C1'"," C2'"," O2'"," N1 "," C2 "," O2 "," N3 "," C4 "," O4 "," C5 "," C6 ",  None,  None,  None," H5'","H5''"," H4'"," H3'"," H2'","HO2'"," H1'"," H3 "," H5 "," H6 ",  None,  None,  None), #30   U
    (" OP1"," P  "," OP2"," O5'"," C5'"," C4'"," O4'"," C3'"," O3'"," C1'"," C2'"," O2'",  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None," H5'","H5''"," H4'"," H3'"," H2'","HO2'"," H1'",  None,  None,  None,  None,  None,  None), #31  RX (unk RNA)
    (" N  "," CA "," C  "," O  "," CB "," CG "," NE2"," CD2"," CE1"," ND1",  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ","2HD ","1HE ","1HD ",  None,  None,  None,  None,  None,  None), #-1 his_d
]

NTOTAL = 36
CHAIN_GAP = 200
to1letter = {
    "ALA":'A', "ARG":'R', "ASN":'N', "ASP":'D', "CYS":'C',
    "GLN":'Q', "GLU":'E', "GLY":'G', "HIS":'H', "ILE":'I',
    "LEU":'L', "LYS":'K', "MET":'M', "PHE":'F', "PRO":'P',
    "SER":'S', "THR":'T', "TRP":'W', "TYR":'Y', "VAL":'V' }


def parse_a3m(filename):
    '''read A3M and convert letters into integers in the 0..20 range,
    also keep track of insertions
    '''

    # read A3M file line by line
    lab,seq = [],[] # labels and sequences
    for line in open(filename, "r"):
        if line[0] == '>':
            lab.append(line.split()[0][1:])
            seq.append("")
        else:
            seq[-1] += line.rstrip()

    # parse sequences
    msa,ins = [],[]
    table = str.maketrans(dict.fromkeys(string.ascii_lowercase))
    nrow,ncol = len(seq),len(seq[0])

    for seqi in seq:

        # remove lowercase letters and append to MSA
        msa.append(seqi.translate(table))

        # 0 - match or gap; 1 - insertion
        a = np.array([0 if c.isupper() or c=='-' else 1 for c in seqi])
        i = np.zeros((ncol))

        if np.sum(a) > 0:
            # positions of insertions
            pos = np.where(a==1)[0]

            # shift by occurrence
            a = pos - np.arange(pos.shape[0])

            # position of insertions in the cleaned sequence
            # and their length
            pos,num = np.unique(a, return_counts=True)
            i[pos[pos<ncol]] = num[pos<ncol]

        # append to the matrix of insetions
        ins.append(i)

    # convert letters into numbers
    alphabet = np.array(list("ARNDCQEGHILKMFPSTWYV-"), dtype='|S1').view(np.uint8)
    msa = np.array([list(s) for s in msa], dtype='|S1').view(np.uint8)
    for i in range(alphabet.shape[0]):
        msa[msa == alphabet[i]] = i

    # treat all unknown characters as gaps
    msa[msa > 20] = 20

    ins = np.array(ins, dtype=np.uint8)

    return {"msa":msa, "labels":lab, "insertions":ins}


def parse_pdb(filename, **kwargs):
    '''extract xyz coords for all heavy atoms'''
    lines = open(filename,'r').readlines()
    return parse_pdb_lines(lines, **kwargs)

def parse_pdb_lines(lines, parse_hetatom=False, ignore_het_h=True):
    # indices of residues observed in the structure
    res = [(l[22:26],l[17:20]) for l in lines if l[:4]=="ATOM" and l[12:16].strip()=="CA"]
    seq = [util.aa2num[r[1]] if r[1] in util.aa2num.keys() else 20 for r in res]
    pdb_idx = [( l[21:22].strip(), int(l[22:26].strip()) ) for l in lines if l[:4]=="ATOM" and l[12:16].strip()=="CA"]  # chain letter, res num


    # 4 BB + up to 10 SC atoms
    xyz = np.full((len(res), 14, 3), np.nan, dtype=np.float32)
    for l in lines:
        if l[:4] != "ATOM":
            continue
        chain, resNo, atom, aa = l[21:22], int(l[22:26]), ' '+l[12:16].strip().ljust(3), l[17:20]
        idx = pdb_idx.index((chain,resNo))
        for i_atm, tgtatm in enumerate(util.aa2long[util.aa2num[aa]]):
            if tgtatm is not None and tgtatm.strip() == atom.strip(): # ignore whitespace
                xyz[idx,i_atm,:] = [float(l[30:38]), float(l[38:46]), float(l[46:54])]
                break

    # save atom mask
    mask = np.logical_not(np.isnan(xyz[...,0]))
    xyz[np.isnan(xyz[...,0])] = 0.0

    # remove duplicated (chain, resi)
    new_idx = []
    i_unique = []
    for i,idx in enumerate(pdb_idx):
        if idx not in new_idx:
            new_idx.append(idx)
            i_unique.append(i)
            
    pdb_idx = new_idx
    xyz = xyz[i_unique]
    mask = mask[i_unique]
    seq = np.array(seq)[i_unique]

    out = {'xyz':xyz, # cartesian coordinates, [Lx14]
            'mask':mask, # mask showing which atoms are present in the PDB file, [Lx14]
            'idx':np.array([i[1] for i in pdb_idx]), # residue numbers in the PDB file, [L]
            'seq':np.array(seq), # amino acid sequence, [L]
            'pdb_idx': pdb_idx,  # list of (chain letter, residue number) in the pdb file, [L]
           }

    # heteroatoms (ligands, etc)
    if parse_hetatom:
        xyz_het, info_het = [], []
        for l in lines:
            if l[:6]=='HETATM' and not (ignore_het_h and l[77]=='H'):
                info_het.append(dict(
                    idx=int(l[7:11]),
                    atom_id=l[12:16],
                    atom_type=l[77],
                    name=l[16:20]
                ))
                xyz_het.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])

        out['xyz_het'] = np.array(xyz_het)
        out['info_het'] = info_het

    return out

def parse_pdb_multi(filename, parse_hetatom=False, ignore_het_h=True, reverse=False):

    lines = open(filename,'r').readlines()

    lines_groups = []
    lines_curr = []
    for line in lines:
        if 'ENDMDL' in line:
            lines_groups.append(lines_curr)
            lines_curr = []
        else:
            lines_curr.append(line)


    all_pdb_idx = []
    all_idx = []
    all_xyz = []
    all_mask = []
    all_seq = []
    all_xyz_het = []
    all_info_het = []

    for lines in lines_groups:
        # indices of residues observed in the structure
        res = [(l[22:26],l[17:20]) for l in lines if l[:4]=="ATOM" and l[12:16].strip()=="CA"]
        seq = [util.aa2num[r[1]] if r[1] in util.aa2num.keys() else 20 for r in res]
        pdb_idx = [( l[21:22].strip(), int(l[22:26].strip()) ) for l in lines if l[:4]=="ATOM" and l[12:16].strip()=="CA"]  # chain letter, res num

        # 4 BB + up to 10 SC atoms
        xyz = np.full((len(res), 14, 3), np.nan, dtype=np.float32)
        for l in lines:
            if l[:4] != "ATOM":
                continue
            chain, resNo, atom, aa = l[21:22], int(l[22:26]), ' '+l[12:16].strip().ljust(3), l[17:20]
            idx = pdb_idx.index((chain,resNo))
            # for i_atm, tgtatm in enumerate(util.aa2long[util.aa2num[aa]]):
            for i_atm, tgtatm in enumerate(util.aa2long[util.aa2num[aa]][:14]): # Nate's proposed change
                if tgtatm is not None and tgtatm.strip() == atom.strip(): # ignore whitespace
                    xyz[idx,i_atm,:] = [float(l[30:38]), float(l[38:46]), float(l[46:54])]
                    break

        # save atom mask
        mask = np.logical_not(np.isnan(xyz[...,0]))
        xyz[np.isnan(xyz[...,0])] = 0.0

        # remove duplicated (chain, resi)
        new_idx = []
        i_unique = []
        for i,idx in enumerate(pdb_idx):
            if idx not in new_idx:
                new_idx.append(idx)
                i_unique.append(i)

        pdb_idx = new_idx
        xyz = xyz[i_unique]
        mask = mask[i_unique]
        seq = np.array(seq)[i_unique]

        all_pdb_idx.append(pdb_idx)
        all_xyz.append(xyz)
        all_idx.append(np.array([i[1] for i in pdb_idx]))
        all_mask.append(mask)
        all_seq.append(np.array(seq))

        # heteroatoms (ligands, etc)
        if parse_hetatom:
            xyz_het, info_het = [], []
            for l in lines:
                if l[:6]=='HETATM' and not (ignore_het_h and l[77]=='H'):
                    info_het.append(dict(
                        idx=int(l[7:11]),
                        atom_id=l[12:16],
                        atom_type=l[77],
                        name=l[16:20]
                    ))
                    xyz_het.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])

            all_xyz_het.append(np.array(xyz_het))
            all_info_het.append(info_het)

    if reverse:
        all_pdb_idx = all_pdb_idx[::-1]
        all_idx = all_idx[::-1]
        all_xyz = all_xyz[::-1]
        all_mask = all_mask[::-1]
        all_seq = all_seq[::-1]
        all_xyz_het = all_xyz_het[::-1]
        all_info_het = all_info_het[::-1]

    out = {'xyz':np.stack(all_xyz), # cartesian coordinates, [Lx14]
            'mask': np.stack(all_mask), # mask showing which atoms are present in the PDB file, [Lx14]
            'idx': np.stack(all_idx), # residue numbers in the PDB file, [L]
            'seq': np.stack(all_seq), # amino acid sequence, [L]
            'pdb_idx': all_pdb_idx,  # list of (chain letter, residue number) in the pdb file, [L]
           }

    if parse_hetatom:
        out['xyz_het']: np.stack(all_xyz_het)
        out['info_het']: all_info_het

    return out
def parse_pdb_rf2aa(filename, seq=True, lig=True): #false/false
    lines = open(filename,'r').readlines()
    if seq:
        if lig:
            return parse_pdb_lines_w_seq_w_lig(lines)
        return parse_pdb_lines_w_seq(lines)
    return parse_pdb_lines(lines)

def parse_pdb_lines_w_seq_w_lig(lines):

    # indices of residues observed in the structure
    res = [(l[21:22].strip(), l[22:26],l[17:20]) for l in lines if l[:4]=="ATOM" and l[12:16].strip()=="CA"] # (chain letter, res num, aa)
    pdb_idx_s = [(r[0], int(r[1])) for r in res]
    idx_s = [int(r[1]) for r in res]
    seq = [aa2num[r[2]] if r[2] in aa2num.keys() else 20 for r in res]

    # 4 BB + up to 10 SC atoms
    xyz = np.full((len(idx_s), NTOTAL, 3), np.nan, dtype=np.float32)
    for l in lines:
        if l[:4] != "ATOM":
            continue
        chain, resNo, atom, aa = l[21:22].strip(), int(l[22:26]), l[12:16], l[17:20]
        idx = pdb_idx_s.index((chain,resNo))
        for i_atm, tgtatm in enumerate(aa2long[aa2num[aa]]):
            if tgtatm == atom:
                xyz[idx,i_atm,:] = [float(l[30:38]), float(l[38:46]), float(l[46:54])]
                break

    # parse ligand atoms
    offset = 0
    if len(idx_s) > 0:
        offset = max(idx_s)
    res_lig = [l[12:16].strip() for l in lines if l[:6]=="HETATM"]
    res_lig = [(i+offset+CHAIN_GAP,l) for i,l in enumerate(res_lig)]
    idx_s_lig = [int(r[0]) for r in res_lig]
    seq_lig = [aa2num[r[1]] if r[1] in aa2num.keys() else 20 for r in res_lig]

    xyz_s_lig = [[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in lines if l[:6]=='HETATM']

    if len(xyz_s_lig)>0:
        xyz_lig = np.full((len(idx_s_lig), NTOTAL, 3), np.nan, dtype=np.float32)
        xyz_lig[:,1,:] = np.array(xyz_s_lig)
        xyz = np.concatenate([xyz, xyz_lig],axis=0)
        idx_s = idx_s + idx_s_lig
        seq = seq + seq_lig

    # save atom mask
    mask = np.logical_not(np.isnan(xyz[...,0]))
    xyz[np.isnan(xyz[...,0])] = 0.0

    return xyz,mask,np.array(idx_s), np.array(seq)

def parse_fasta(filename):
    '''
    Return dict of name: seq
    '''
    out = {}
    with open(filename, 'r') as f_in:
        while True:
            name = f_in.readline().strip()[1:]
            seq = f_in.readline().strip()
            if not name: break

            out[name] = seq

    return out
