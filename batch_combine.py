#!/home/pranav/mambaforge/bin/python

# In[ ]: Import
import os
import scanpy as sc
import scvi
import seaborn as sns
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix

import logging

scvi.settings.verbosity = logging.ERROR

# In[]: Get ribosome genes

ribo_url = "http://software.broadinstitute.org/gsea/msigdb/download_geneset.jsp?geneSetName=KEGG_RIBOSOME&fileType=txt"

ribo_genes = pd.read_table(ribo_url, skiprows=2, header=None)
ribo_genes


# In[]: Define functions
def pp(csv_path, i):
    adata = sc.read_csv(csv_path).T

    adata.X = csr_matrix(adata.X)

    sc.pp.filter_genes(adata, min_cells=10)
    sc.pp.highly_variable_genes(
        adata, n_top_genes=2000, subset=True, flavor="seurat_v3"
    )
    scvi.model.SCVI.setup_anndata(adata)
    vae = scvi.model.SCVI(adata)
    vae.train()
    solo = scvi.external.SOLO.from_scvi_model(vae)
    solo.train()
    df = solo.predict()
    df["prediction"] = solo.predict(soft=False)
    df.index = df.index.map(lambda x: x[:-2])
    df["dif"] = df.doublet - df.singlet
    doublets = df[(df.prediction == "doublet") & (df.dif > 1)]

    adata = sc.read_csv(csv_path).T
    adata.obs["Sample"] = csv_path.split("_")[
        2
    ]  #'raw_counts/GSM5226574_C51ctr_raw_counts.csv'

    adata.obs["doublet"] = adata.obs.index.isin(doublets.index)
    adata = adata[~adata.obs.doublet]

    sc.pp.filter_cells(
        adata, min_genes=200
    )  # get rid of cells with fewer than 200 genes
    # sc.pp.filter_genes(adata, min_cells=3) #get rid of genes that are found in fewer than 3 cells
    adata.var["mt"] = adata.var_names.str.startswith(
        "mt-"
    )  # annotate the group of mitochondrial genes as 'mt'
    adata.var["ribo"] = adata.var_names.isin(ribo_genes[0].values)
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo"], percent_top=None, log1p=False, inplace=True
    )
    upper_lim = np.quantile(adata.obs.n_genes_by_counts.values, 0.98)
    adata = adata[adata.obs.n_genes_by_counts < upper_lim]
    adata = adata[adata.obs.pct_counts_mt < 20]
    adata = adata[adata.obs.pct_counts_ribo < 2]

    adata.write_h5ad("./output/{i}_data.h5ad")

    return adata


# In[]: Run Loop

folder_dir = "/home/pranav/work/pranavmishra90/research-reference/nextgen-analysis/lung-atlas-tutorial"


# In[]:

files_list = []

for root, directories, files in os.walk(folder_dir):
    for name in files:
        if name.endswith(("csv")):
            files_list.append(os.path.join(root, name))

files_list = sorted(files_list)

print(files_list)

# In[]:

combined_adata = []
i = 1

for sample in files_list:
    pp(sample, i)

    combined_adata.append(pp(sample, i))
    adata = sc.concat(combined_adata)
    adata.write_h5ad(f'./output/{i}_combined.h5ad')

    i = i + 1

# In[]: Combine all and save

combined_adata.write_h5ad(f'./output/all_combined.h5ad')