repo for cheme 546 project
scaffold generation for ion binding sites
Division of labor of tasks:
Jack- Train GNN model
Long- Aid with training, optimize learning rate and hyperparameters
Peik - Extract and normalize data for training and input
Annika- 
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

### references
util function is built on https://github.com/RosettaCommons/RFdiffusion
original periodic tabale data sourced from https://gist.github.com/GoodmanSciences/c2dd862cd38f21b0ad36b8f96b4bf1ee#file-periodic-table-of-elements-csv-L9 and https://gist.github.com/speters33w/e5b1246d6859f29c4f02a299714d4c20
metal binding database sourced from https://intermetaldb.biotech.uni.wroc.pl/search_metal_site/https://arxiv.org/pdf/2101.03064.pdf
One class SVM pseudocodes 
https://openaccess.thecvf.com/content_ICCV_2019/papers/Wang_GODS_Generalized_One-Class_Discriminative_Subspaces_for_Anomaly_Detection_ICCV_2019_paper.pdf
https://github.com/iqiukp/SVDD-Python