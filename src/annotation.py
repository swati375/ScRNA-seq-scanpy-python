import scanpy as sc
import numpy as np

## add more markers from publications as required for your data. Some common ones are listed below.

def get_marker_dict():
    return {
        "B cells": ["MS4A1", "CD79A", "CD79B"],
        "T cells": ["CD3D", "CD3E", "CD2"],
        "NK cells": ["GNLY", "NKG7", "KLRD1"],
        "Monocytes": ["LYZ", "CD14", "FCGR3A"],
        "Epithelial": ["EPCAM", "KRT18", "KRT8"],
        "Endothelial": ["PECAM1", "VWF", "KDR"],
        "Mast cells": ["TPSAB1", "KIT", "CPA3"]
    }

def annotate_clusters(adata):

    marker_dict = get_marker_dict()
    cluster_labels = {}

    for cluster in adata.obs["clusters"].unique():
        
        # subset cluster
        cells = adata[adata.obs["clusters"] == cluster]

        scores = {}

        for cell_type, markers in marker_dict.items():
            valid_genes = [g for g in markers if g in adata.var_names]

            if len(valid_genes) == 0:
                scores[cell_type] = 0
                continue

            # mean expression across markers
            expr = cells[:, valid_genes].X

            if hasattr(expr, "toarray"):
                expr = expr.toarray()

            scores[cell_type] = np.mean(expr)

        # assign best matching cell type
        best_type = max(scores, key=scores.get)
        cluster_labels[cluster] = best_type

    return cluster_labels

def annotate_clusters_deg(adata, n_top=10):

    marker_dict = get_marker_dict()
    result = adata.uns["rank_genes_groups"]
    groups = result["names"].dtype.names

    cluster_labels = {}

    for group in groups:

        # top DE genes for this cluster
        top_genes = list(result["names"][group][:n_top])

        scores = {}

        for cell_type, markers in marker_dict.items():
            overlap = len(set(top_genes) & set(markers))
            scores[cell_type] = overlap

        # assign best match
        best_type = max(scores, key=scores.get)

        # fallback if no overlap
        if scores[best_type] == 0:
            best_type = "Unknown"

        cluster_labels[group] = best_type

    return cluster_labels



