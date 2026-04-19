import scanpy as sc
import pandas as pd


def compute_deg(adata, groupby="clusters"):
    
    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        method="wilcoxon"
    )

    return adata

def get_top_genes(adata, n_genes=5):
    
    result = adata.uns["rank_genes_groups"]
    groups = result["names"].dtype.names

    top_genes = {
        group: result["names"][group][:n_genes]
        for group in groups
    }

    return top_genes

# def save_deg_results(adata, filename="results/deg_results.csv"):
    
#     result = adata.uns["rank_genes_groups"]
#     groups = result["names"].dtype.names

#     df = pd.DataFrame({
#         group: result["names"][group]
#         for group in groups
#     })

#     df.to_csv(filename)


# Step 1: Build DEG DataFrame
# -------------------------------
def build_deg_table(adata):

    result = adata.uns["rank_genes_groups"]
    groups = result["names"].dtype.names

    dfs = []

    for group in groups:
        df = pd.DataFrame({
            "gene": result["names"][group],
            "logfc": result["logfoldchanges"][group],
            "pval_adj": result["pvals_adj"][group],
            "cluster": group
        })
        dfs.append(df)

    deg_df = pd.concat(dfs)

    return deg_df


# -------------------------------
# Step 2: Filter DEG
# -------------------------------
def filter_deg_table(deg_df, logfc_thresh=1.0, pval_thresh=0.05):

    filtered = deg_df[
        (deg_df["logfc"] > logfc_thresh) &
        (deg_df["pval_adj"] < pval_thresh)
    ]

    return filtered