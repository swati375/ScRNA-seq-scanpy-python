import pandas as pd
import scanpy as sc
from src.deg import build_deg_table,filter_deg_table

def get_marker_dict():
    return {
        "B cells": ["MS4A1", "CD79A", "CD79B"],
        "T cells": ["CD3D", "CD3E", "CD2"],
        "NK cells": ["GNLY", "NKG7", "KLRD1"],
        "Monocytes": ["LYZ", "CD14", "FCGR3A"],
        "Epithelial": ["EPCAM", "KRT18", "KRT8"],
        "Endothelial": ["PECAM1", "VWF", "KDR"],
        "Mast cells": ["TPSAB1", "KIT", "CPA3"],
        "DC": ["FCER1A","CST3"],
        "Platelet":["PPBP"]
    }


# -------------------------------
# DEG-based Annotation
# -------------------------------
def annotate_clusters_deg_filtered(
    adata,
    logfc_thresh=1.0,
    pval_thresh=0.05
):

    marker_dict = get_marker_dict()

    # Build DEG table
    deg_df = build_deg_table(adata)

    # Filter DEG
    filtered_deg = filter_deg_table(
        deg_df,
        logfc_thresh=logfc_thresh,
        pval_thresh=pval_thresh
    )

    cluster_labels = {}

    for cluster in filtered_deg["cluster"].unique():

        cluster_genes = filtered_deg[
            filtered_deg["cluster"] == cluster
        ]["gene"].tolist()

        scores = {}

        for cell_type, markers in marker_dict.items():

            # keep only markers present in dataset
            valid_markers = [g for g in markers if g in adata.var_names]

            overlap = len(set(cluster_genes) & set(valid_markers))
            scores[cell_type] = overlap

        best_type = max(scores, key=scores.get)

        if scores[best_type] == 0:
            best_type = "Unknown"

        cluster_labels[cluster] = best_type

    return cluster_labels, deg_df, filtered_deg