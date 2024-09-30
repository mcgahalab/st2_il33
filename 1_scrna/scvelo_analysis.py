# Workaround for ImportError: cannot import name 'Iterable' from 'collections' (/cluster/home/quever/miniconda3/envs/scvelo/lib/python3.10/collections/__init__.py)
# https://stackoverflow.com/questions/72032032/importerror-cannot-import-name-iterable-from-collections-in-python
import collections.abc
#hyper needs the four following aliases to be done manually.
collections.Iterable = collections.abc.Iterable
collections.Mapping = collections.abc.Mapping
collections.MutableSet = collections.abc.MutableSet
collections.MutableMapping = collections.abc.MutableMapping

import os
import numpy as np
import pandas as pd
import scvelo as scv
import scanpy as sc
import cellrank as cr
from cellrank.kernels import VelocityKernel, ConnectivityKernel, CytoTRACEKernel
from cellrank.estimators import GPCCA
from pathlib import Path
import matplotlib.backends.backend_pdf as backend_pdf
from matplotlib import pyplot as plt
import matplotlib as mp
import re

import sys
sys.path.append("/cluster/home/quever/git/scvelo_helper")
import custom.analysis as st2
import custom.scvelo as customvelo
print("Loaded")
####################
#### Parameters ####
# Setup:
model = 'Tumor'   # LN or Tumor or LN_Tumor
# term = 'wt'  # ko, wt, or merged
batch='B1' # B1, B2
label_id = 'treg_anno' #'cell_type' # 'treg_anno'
sample_id = 'orig.ident'
basis_id = 'umap_orig' # umap umap_orig
celltype = 'cd8' # '.treg' # 'tregs' 'cd8'

outdir = Path("/cluster/home/quever/xfer/")
indir = Path("/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/scrna_tdln_tumor_7d/results/scVelo/", celltype)
loomdir = Path("/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/scrna_tdln_tumor_7d/data/velocyto/loom_files")
# indir = Path("/cluster/projects/mcgahalab/data/mcgahalab/siran_ifnblockade/results/velocity")
# loomdir = Path("/cluster/projects/mcgahalab/data/mcgahalab/siran_ifnblockade/data/velocyto/loom_files")
# wtdata = scv.read(indir / ("seu_velo." + model + "." + batch + ".WT.h5ad"))
# kodata = scv.read(indir / ("seu_velo." + model + "." + batch + ".KO.h5ad"))
# mergeddata = scv.read(indir / ("seu_velo." + model + "." + batch + ".all.h5ad"))
# datadict = {
#     'wt' : wtdata,
#     'ko' : kodata,
#     'merged' : mergeddata
# }

#adata = scv.read("seu_velo.h5ad")
scv.settings.verbosity = 3
scv.settings.presenter_view = True
scv.settings.set_figure_params('scvelo')
cr.settings.verbosity = 2


##  read in seurat RNA counts data
os.chdir(indir)
model_celltype = model if celltype == '' else (model + "." + celltype)
adata = st2.read_in_seurat_flatfile(
    counts=model_celltype + '.counts.mtx',
    metadata=model_celltype + '.metadata.csv',
    pca=model_celltype + '.pca.csv',
    features=model_celltype + '.genes.csv',
    h5ad=model_celltype + '.h5ad')

##  read in velocyto spliced/unspliced counts data
# adata = sc.read_h5ad('LN.treg.h5ad')
barcodes = adata.obs.index.to_list()
samples = [re.sub("_[ACGT]*$", "", a) for a in barcodes]
samples = pd.Series(samples).unique()

## read in velocyto loom filees and merge matrices into the original adata object
ldata = st2.read_velocyto_files(loomdir, samples)
scv.utils.clean_obs_names(adata)
scv.utils.clean_obs_names(ldata)
adata = scv.utils.merge(adata, ldata)

# plot umap to check

# adata.obs['ln_tumor'] = [re.sub("_.*", "", a) for a in adata.obs[label_id]]
# sc.pl.umap(adata, color='ln_tumor', frameon=False, legend_loc='on data',
#     title='', save=('_lntumor.' + model + '.pdf'))

# pre-process
adatas={}
# for subset in ['Tumor_KO', 'Tumor_WT', 'LN_KO', 'LN_WT']:
for subset in ['Tumor_KO_7d', 'Tumor_WT_7d', 'Tumor_KO_Un', 'Tumor_WT_Un']:
# for subset in ['all']:
    print(subset)
    if subset == 'all':
        model_sub = model
        adata2 = adata
    else:
        model_sub = (model + ".sub" + subset)
        adata2 = adata[ [subset in i for i in adata.obs['orig.ident']] ]
    adata2 = st2.preprocess_velo(adata2, model_sub, indir, label_id=label_id, mode='stochastic', issubset=True)
    st2.plot_velocity(adata2, model_sub, label_id, basis='umap')
    adatas[subset]=adata2

id1='Tumor_WT_7d'
id2='Tumor_KO_7d'
outf = "./figures/" + "diff." + id1 + "." + id2 + ".pdf"
adatasub = {key: adatas[key] for key in [id1, id2]}
diff_velocity_graph(adatasub, 'umap', label_id, legend_loc='best',
    legend_fontsize=6, legend_fontweight='normal',
    refkey=id1, altkey=id2, save=outf)
