repo for cheme 546 project

Metal Binding Prediction for Protein Scaffolds

### environment
conda env create -f environment.yml

### dataset curation
initial pdb was curateed from https://intermetaldb.biotech.uni.wroc.pl/search_metal_site/
pdb was download from rcsb.com
extract motif on curated pdbs

### demo
demo jupyternote book is in demo/

load model: model loaded from weight4/(model_metal).pkl
example_data = preprocess pdb by utl.preprocess_tip_atom(example.pdb)
model.predict(example_data)


### EXAMPLE USAGE:
# for predicting a single Zn binding site

python run_scaffold_pred.py --pdb PATH_TO_PDB --metal_id ATOM_ID_NUM --metal_type ZN
outputs: a txt file PDB_ZN.txt with prediction

# for predicting binding ability of multiple scaffolds - make csv of pdb paths, pdb atom index number of metal, and metal type

# csv columns must be labeled "pdb", "metal_idx", and "metal"

python run_scaffold_pred_batch.py --csv PATH_TO_CSV --output OUTPUT_CSV_PATH
outputs: csv file with pdb ids and binding prediction




### references
util function is built on https://github.com/RosettaCommons/RFdiffusion
original periodic tabale data sourced from https://gist.github.com/GoodmanSciences/c2dd862cd38f21b0ad36b8f96b4bf1ee#file-periodic-table-of-elements-csv-L9 and https://gist.github.com/speters33w/e5b1246d6859f29c4f02a299714d4c20
metal binding database sourced from https://intermetaldb.biotech.uni.wroc.pl/search_metal_site/https://arxiv.org/pdf/2101.03064.pdf
One class SVM pseudocodes 
https://openaccess.thecvf.com/content_ICCV_2019/papers/Wang_GODS_Generalized_One-Class_Discriminative_Subspaces_for_Anomaly_Detection_ICCV_2019_paper.pdf
https://github.com/iqiukp/SVDD-Python
