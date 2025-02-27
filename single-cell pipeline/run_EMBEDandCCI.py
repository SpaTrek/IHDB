import anndata as ad
import scanpy as sc
import celltypist
import pandas as pd
import numpy as np
import os
import liana as li
import copy
from liana.method import cellphonedb
human_datasets=[]#list of human sample IDs

resource = pd.read_table("/pathway/to/mouse_lr_pair.txt")
resource = resource.iloc[:,1:3]
resource.columns = pd.Index(["ligand","receptor"])

path = "/path/to/all/preprocessed/h5ad/folder/"
path_csv = "/path/to/all/output/csv/folder/"
li = list(os.walk(path))[0][2]

ot = list(os.walk(path_csv))[0][2]
for n in li:
    
    adata = sc.read_h5ad(path+n)
    
    # remove subtypes with few cells
    adata.obs["subtype_revised"] = copy.deepcopy(adata.obs["subtype"])
    ids = pd.read_csv(path_csv+n[:10]+"_hv.csv")
    app = pd.read_csv(path_csv+n[:10]+"_app.csv")
    marker = pd.read_csv(path_csv+n[:10]+"_marker.csv")
    for i in adata.obs["subtype_revised"].value_counts().index[adata.obs["subtype_revised"].value_counts()<=10].to_list():
        if i not in list(set(sum([j.split(" | ") for j in list(ids["cell_type"].unique())],[]))):
            print(i)

        if "." in i:
            adata.obs["subtype_revised"][adata.obs["subtype_revised"]==i] = i.split(".")[0]+".0"
            adata.obs["subtype_revised"] = adata.obs["subtype_revised"].cat.remove_unused_categories()
            genes_rem = ids.loc[(ids["cell_type"].str.startswith(i)),]["gene"].to_list()
            ids = ids.loc[~(ids["cell_type"].str.startswith(i)),]
            for k in ids.index:
                if i in ids.loc[k,"cell_type"]:
                    ids.loc[k,"cell_type"] = ids.loc[k,"cell_type"].replace(" | "+i,"")
            app = app.loc[:,~app.columns.isin(genes_rem)]
            marker = marker.loc[:,~marker.columns.isin(genes_rem)]
        else:
            continue
    
    app.to_csv(path_csv+n[:10]+"_app.csv", index = False, lineterminator = "\r\n")
    marker.to_csv(path_csv+n[:10]+"_marker.csv", index = False, lineterminator = "\r\n")
    ids.to_csv(path_csv+n[:10]+"_hv.csv", index = False, lineterminator = "\r\n")
    
    sc.tl.draw_graph(adata, layout = "fa")
    sc.tl.tsne(adata)
    sc.tl.draw_graph(adata, layout = "fr")
    sc.tl.draw_graph(adata, layout = "drl")
    sc.tl.draw_graph(adata, layout = "kk")
    sc.tl.diffmap(adata, n_comps=6)

    basic = pd.DataFrame(np.concatenate(
            (np.array(adata.obs_names)[:,None],    # cell barcode
            np.array(adata.obs["celltypist_cell_label"])[:,None],    # leiden group numbers
            np.multiply(adata.obsm["X_tsne"],100).astype(int),    # tsne (x,y) => (100x,100y)
            np.multiply(adata.obsm["X_umap"],100).astype(int),    # umap (x,y) => (100x,100y)
            np.multiply(adata.obsm["X_diffmap"],100000).astype(int),    # diffusion components
            np.array(adata.obs["subtype_revised"])[:,None],    # celltypist cell type
            np.multiply(adata.obsm["X_draw_graph_fa"],100).astype(int),    # umap (x,y) => (100x,100y)
            np.multiply(adata.obsm["X_draw_graph_fr"],100).astype(int),    # umap (x,y) => (100x,100y)
            np.multiply(adata.obsm["X_draw_graph_drl"],100).astype(int),    # umap (x,y) => (100x,100y)
            np.multiply(adata.obsm["X_draw_graph_kk"],100).astype(int),    # umap (x,y) => (100x,100y)
            #np.round(log1p_deg50.T,decimals=2),    # log1p data
            #np.round(cell_embeddings,decimals=3),    # programs expression level
            #np.multiply(adata.obsm["X_diffmap"],100000).astype(int),    # diffusion components
            #####np.round(size_factors,decimals=4)[:,None],    #scran's size factors
            #np.round(np.array(adata.obs["celltypist_conf_score"]),decimals=3)[:,None],    # celltypist confidence score
            #np.round(kernel_density, decimals=2),    # Gaussion kernel estimation
            ), axis=1))
    basic.columns = pd.Series(["index","celltype","snex","sney","umapx","umapy","diff_1","diff_2","diff_3","diff_4","diff_5","diff_6","subtype","fax","fay","frx","fry","drlx","drly","kkx","kky",])
    basic.to_csv(path_csv+n[:10]+"_id.csv", index = False, lineterminator = "\r\n")

    if n[:10] in human_datasets:
        cellphonedb(adata, use_raw=False, groupby='subtype_revised', expr_prop=0.1, verbose=False, key_added='cci_minor', resource_name="celltalkdb")
        cellphonedb(adata, use_raw=False, groupby='celltypist_cell_label', expr_prop=0.1, verbose=False, key_added='cci_major', resource_name="celltalkdb")
    else:
        cellphonedb(adata, use_raw=False, groupby='subtype_revised', expr_prop=0.1, verbose=False, key_added='cci_minor', resource=resource)
        cellphonedb(adata, use_raw=False, groupby='celltypist_cell_label', expr_prop=0.1, verbose=False, key_added='cci_major', resource=resource)




    save_cci_minor = adata.uns["cci_minor"].loc[(adata.uns["cci_minor"]["cellphone_pvals"]<0.05),][1:500]
    save_cci_minor = save_cci_minor.round({'ligand_means': 2, 'ligand_props': 2, 'receptor_means':2, 'receptor_props':2, 'lr_means':2})
    save_cci_minor.to_csv(path_csv+n[:10]+"_cci.csv", index = True, lineterminator = "\r\n")

    save_cci_major = adata.uns["cci_major"].loc[(adata.uns["cci_major"]["cellphone_pvals"]<0.05),][1:500]
    save_cci_major = save_cci_major.round({'ligand_means': 2, 'ligand_props': 2, 'receptor_means':2, 'receptor_props':2, 'lr_means':2})
    save_cci_major.to_csv(path_csv+n[:10]+"_cci_major.csv", index = True, lineterminator = "\r\n")
    adata.write(path+n)