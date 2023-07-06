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



####################
#### Parameters ####
# Setup:
model = 'Tumor'   # LN or Tumor
term = 'wt'  # ko, wt, or merged
batch='B1' # B1, B2
label_id = 'manual_anno'
sample_id = 'orig.ident'
basis_id = 'umap_orig' # umap umap_orig


outdir = Path("/cluster/home/quever/xfer/")
indir = Path("/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/scrna_tdln_tumor_7d/results/scVelo/")
wtdata = scv.read(indir / ("seu_velo." + model + "." + batch + ".WT.h5ad"))
kodata = scv.read(indir / ("seu_velo." + model + "." + batch + ".KO.h5ad"))
mergeddata = scv.read(indir / ("seu_velo." + model + "." + batch + ".all.h5ad"))
datadict = {
    'wt' : wtdata,
    'ko' : kodata,
    'merged' : mergeddata
}

#adata = scv.read("seu_velo.h5ad")
scv.settings.verbosity = 3
scv.settings.presenter_view = True
scv.settings.set_figure_params('scvelo')
cr.settings.verbosity = 2
os.chdir("/cluster/home/quever")

adatas = {}
for term in ['ko', 'wt']:
    print(term)
    adata = preprocess_velo(model, batch, term, indir, datadict)
    os.chdir("/cluster/home/quever")
    adata = velocity_analysis(adata, outdir, model, batch, term, sample_id, basis_id, label_id)
    adatas[term] = adata

outf = Path(outdir, model + "_" + batch + ".pdf")
diff_velocity_graph(adatas, basis_id, label_id, legend_loc='best', legend_fontsize=6, legend_fontweight='normal')
plt.savefig(outf, format="pdf")
plt.show()


########################################
#### Preprocessing and loading data ####
# Preprocess the data
def preprocess_velo(model, batch, term, indir, datadict):
    processed_f = "seu_velo." + model + "." + batch + "." + term + ".processed.h5ad"
    processed_f_path = indir / processed_f
    adata = datadict[term]
    if os.path.exists(processed_f_path):
        print("loading...")
        adata = scv.read(processed_f_path)
    else:
        print("processing...")
        if True:
            scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=3000)
            sc.tl.pca(adata)
            sc.pp.neighbors(adata, n_pcs=30, n_neighbors=30)
            scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
            scv.tl.recover_dynamics(adata) # learn full transcriptional dynamics of splicing kinetics
            scv.tl.velocity(adata, mode='dynamical')
            # var_names = ['Ly6c2', 'Ly6g', 'Itgam', 'H2-Ab1']
            #scv.tl.differential_kinetic_test(adata, groupby='seurat_clusters') # var_names=var_names,
            #scv.tl.velocity(adata, diff_kinetics=True)
            scv.tl.velocity_graph(adata)
            scv.tl.velocity_confidence(adata)
            scv.tl.latent_time(adata)
            
            # Saving error _index fix: https://github.com/theislab/scvelo/issues/255
            adata.__dict__['_raw'].__dict__['_var'] = adata.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
            # del adata.raw # dotplot error: f"Could not find keys '{not_found}' in columns of `adata.{dim}` or in"
            # del(adata.var['_index']) #if still error while saving
            adata.write(filename=processed_f_path, compression='gzip')
    return adata

adata = preprocess_velo(model, batch, term, indir, datadict)

############
#### QC ####
def velo_qc(adata, outdir, model, term, label_id):
    # Kinetic Rate parameters
    df = adata.var
    df = df[(df['fit_likelihood'] > .1) & df['velocity_genes'] == True]
    scv.get_df(adata, 'fit*', dropna=True).head()
    
    kwargs = dict(xscale='log', fontsize=16)
    with scv.GridSpec(ncols=3) as pl:
        pl.hist(df['fit_alpha'], xlabel='transcription rate', **kwargs)
        pl.hist(df['fit_beta'] * df['fit_scaling'], xlabel='splicing rate', xticks=[.1, .4, 1], **kwargs)
        pl.hist(df['fit_gamma'], xlabel='degradation rate', xticks=[.1, .4, 1], **kwargs)
    
    outf = outdir / (model + "." + term + ".kinetics.pdf")
    plt.savefig(outf, format="pdf", bbox_inches="tight")
    plt.show()
    
    # Proportions of spliced to unspliced
    scv.pl.proportions(adata, groupby=label_id)
    outf = outdir / (model + "." + term + ".proportions.pdf")
    plt.savefig(outf, format="pdf", bbox_inches="tight")
    plt.show()
    
    print("xfer " + model + "." + term + ".kinetics.pdf")
    print("xfer " + model + "." + term + ".proportions.pdf")
    
    return adata

adata = velo_qc(adata, outdir, model, term, label_id)


###########################
#### Velocity Analysis ####
def velocity_analysis(adata, outdir, model, batch, term, sample_id, basis_id, label_id):
    #---------- Overall velocity and latent-time analysis
    # Velocity_1: Streamline velocity projection
    outf = outdir / (model + "." + batch + "." + term + ".scvelo.vel_stream.png")
    scv.pl.velocity_embedding_stream(adata, basis=basis_id, color=label_id)
    plt.savefig(outf, format="png", bbox_inches="tight")
    plt.show()
    
    
    # Velocity_2: Cellular and Gridlines level velocity projection
    outf = outdir / (model + "." + batch + "." + term + ".scvelo.vel_projection.pdf")
    pdf = backend_pdf.PdfPages(outf)
    # ---- Gridline velocity
    fig = scv.pl.velocity_embedding_grid(adata, basis=basis_id, color=label_id,
        arrow_length=3, arrow_size=2, dpi=120,
        legend_loc='best', legend_fontsize=6, legend_fontweight='normal')
    pdf.savefig( fig )
    fig = scv.pl.velocity_embedding_grid(adata, basis=basis_id, color=sample_id,
        arrow_length=3, arrow_size=2, dpi=120,
        legend_loc='best', legend_fontsize=6, legend_fontweight='normal')
    pdf.savefig( fig )
    # ---- Cellular velocity
    fig = scv.pl.velocity_embedding(adata, basis=basis_id, color=label_id,
        arrow_length=3, arrow_size=2, dpi=120,
        legend_loc='right margin', legend_fontsize=4, legend_fontweight='normal')
    pdf.savefig( fig )
    fig = scv.pl.velocity_embedding(adata, basis=basis_id, color=sample_id,
        arrow_length=3, arrow_size=2, dpi=120,
        legend_loc='right margin', legend_fontsize=4, legend_fontweight='normal')
    pdf.savefig( fig )
    # ---- Velocity confidence
    keys = 'velocity_length', 'velocity_confidence'
    fig = scv.pl.scatter(adata, basis=basis_id, c=keys, cmap='coolwarm', perc=[5, 95])
    pdf.savefig( fig )
    # ---- Splicing efficiency: fraction of unspliced counts [https://github.com/theislab/scvelo/issues/173]
    counts_s = scv.utils.sum_var(adata.layers['spliced'])
    counts_u = scv.utils.sum_var(adata.layers['unspliced'])
    fractions_u = counts_u / (counts_s + counts_u)
    fig = scv.pl.scatter(adata, color=fractions_u, basis=basis_id, smooth=True)
    pdf.savefig( fig )
    # ---- Pseudotime projections
    fig = scv.pl.scatter(adata,  basis=basis_id , color="latent_time", color_map="gnuplot")
    pdf.savefig( fig )
    pdf.close()
    
    #---------- Top and selected velocity driver genes
    # ---- Top 100 genes for velocity drivers
    all_genes = adata.var["fit_likelihood"].sort_values(ascending=False).index
    top_genes = all_genes[:100]
    scv.pl.heatmap(adata, var_names=top_genes, sortby="latent_time", col_color=label_id,
        n_convolve=100, font_scale=0.5, figsize=(8,15))
    outf = outdir / (model + "." + batch + "." + term + ".scvelo.topheatmap.png")
    plt.savefig(outf, format="png", bbox_inches="tight")
    plt.show()
    
    outf = outdir / (model + "." + batch + "." + term + ".scvelo.gene_velocities.pdf")
    pdf = backend_pdf.PdfPages(outf)
    # ---- Scatterplot of top genes velocities
    fig = scv.pl.scatter(adata, basis=top_genes[:15], ncols=5, color=label_id, frameon=False)
    pdf.savefig( fig )
    ## ---- Scatterplot of selected genes gene_velocities
    #fig = scv.pl.scatter(adata, ['Ahr', 'Ly6c2', 'Ly6g'], ncols=2, color=label_id)
    #pdf.savefig( fig )
    ## ---- Cluster-specific top-likelihood genes
    #scv.tl.rank_dynamical_genes(adata, groupby=label_id)
    #df = scv.get_df(adata, 'rank_dynamical_genes/names')
    ##df.head(5)
    #for cluster in ['Treg_LTlike', 'Treg_NLT', 'Treg_NLTlike', 'Treg_Stat1', 'Treg_effector', 'Treg_lymphoid']:
    #    fig = scv.pl.scatter(adata, df[cluster][:5], ylabel=cluster, color=label_id, frameon=False)
    #    pdf.savefig( fig )
    pdf.close()
    print("xfer " + model + "." + batch + "." + term + ".scvelo.topheatmap.png")
    print("xfer " + model + "." + batch + "." + term + ".scvelo.vel_stream.png")
    print("xfer " + model + "." + batch + "." + term + ".scvelo.vel_projection.pdf")
    print("xfer " + model + "." + batch + "." + term + ".scvelo.gene_velocities.pdf")
    return adata

adata = velocity_analysis(adata, outdir, model, batch, term, sample_id, basis_id, label_id)

for term in ['ko', 'wt']:
    print(term)
    adata = preprocess_velo(model, batch, term, indir, datadict)
    os.chdir("/cluster/home/quever")
    adata = velocity_analysis(adata, outdir, model, batch, term, sample_id, basis_id, label_id)

###########################
#### CellRank Analysis ####
adata,g = cellrank_analysis(adata, outdir, model, term, label_id, sample_id, basis_id)

def cellrank_analysis(adata, outdir, model, term, label_id, sample_id, basis_id):
    vk = VelocityKernel(adata)
    vk.compute_transition_matrix()
    ck = ConnectivityKernel(adata).compute_transition_matrix()
    combined_kernel = 0.8 * vk + 0.2 * ck
    
    g = GPCCA(combined_kernel)
    g.compute_schur(n_components=20, method='brandts')
    g.compute_macrostates(n_states=6, cluster_key=label_id)
    
    g.predict_initial_states()
    g.predict_terminal_states()
    g.compute_absorption_times(use_petsc = False)
    print("xfer " + model + "." + term + ".cellrank.pdf")
    return(adata, g)
    
outf = outdir / (model + "." + term + ".cellrank.pdf")
pdf = backend_pdf.PdfPages(outf)
fig = g.plot_spectrum()
pdf.savefig(fig)
fig = g.plot_macrostates(basis=basis_id)
pdf.savefig(fig)
fig = g.plot_macrostates(same_plot=False, basis=basis_id)
pdf.savefig(fig)
fig = g.plot_macrostates(discrete=True, basis=basis_id)
pdf.savefig(fig)
pdf.close()
    
    
    #alpha_drivers = g.compute_lineage_drivers(lineages="Treg_effector_3", return_drivers=True)
    #alpha_drivers.sort_values(by="Treg_effector_3_corr", ascending=False)
    
    #outf = outdir / (model + "." term + ".cellrank2.pdf")
    #pdf = backend_pdf.PdfPages(outf)
    #fig = g.plot_lineage_drivers("Treg_effector_3", n_genes=5, basis=basis_id)
    #pdf.savefig(fig)
    
    #fig = cr.pl.circular_projection(adata, keys="functional_cluster", legend_loc="right")
    #pdf.savefig(fig)
    #pdf.close()
    
    ctk = CytoTRACEKernel(adata).compute_cytotrace()
    ctk.compute_transition_matrix(threshold_scheme="soft", nu=0.5)
    ctk.plot_projection(basis=basis_id)
    
    outf = outdir / (model + "." + term + ".cellrank3.pdf")
    pdf = backend_pdf.PdfPages(outf)
    fig = scv.pl.scatter(
        adata,
        c=["ct_pseudotime", sample_id],
        basis=basis_id,
        legend_loc="right",
        color_map="gnuplot2",
    )
    pdf.savefig(fig)
    
    fig = scv.pl.scatter(
        adata, basis=basis_id, c=["term_states_fwd", 'init_states_fwd'], legend_loc="right"
    )
    pdf.savefig(fig)
    
    fig = scv.pl.scatter(
        adata, basis=basis_id, c=["seurat_clusters", label_id], legend_loc="right"
    )
    pdf.savefig(fig)
    
    fig = ctk.plot_random_walks(
        n_sims=15,
        start_ixs={"seurat_clusters": "3"},
        basis=basis_id,
        color=sample_id,
        legend_loc="right",
        seed=1,
    )
    pdf.savefig(fig)
    
    fig = scv.pl.velocity_embedding_stream(
        adata, color="ct_pseudotime", basis=basis_id, legend_loc="right"
    )
    pdf.savefig(fig)
    pdf.close()
    
    g = GPCCA(ctk)
    g.compute_schur(n_components=20, method='brandts')
    g.compute_macrostates(n_states=5, cluster_key=label_id)
    
    outf = outdir / (model + "." + term + ".cellrank4.pdf")
    pdf = backend_pdf.PdfPages(outf)
    fig=g.plot_spectrum(real_only=True)
    pdf.savefig(fig)
    fig = g.plot_macrostates(
        discrete=True, legend_loc="right", size=100, basis=basis_id, which='all'
    )
    pdf.savefig(fig)
    pdf.close()
