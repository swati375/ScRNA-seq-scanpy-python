import scanpy as sc

def run_dimensionality_reduction(adata, n_pcs=30):
    
    # Highly variable genes (important!)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    # adata = adata[:, adata.var.highly_variable]
    # Create separate object, keeping only highly variable genes. This is used for dim. reduction but for deg all genes used. 
    adata_hvg = adata[:, adata.var.highly_variable].copy()

    # Scale
    sc.pp.scale(adata_hvg, max_value=10)

    # PCA
    sc.tl.pca(adata_hvg, svd_solver="arpack")

    return adata, adata_hvg


def run_neighbors(adata_hvg, n_neighbors=10, n_pcs=30):
    sc.pp.neighbors(adata_hvg, n_neighbors=n_neighbors, n_pcs=n_pcs)
    return adata_hvg


def run_clustering(adata_hvg, resolution=0.5):
    sc.tl.leiden(adata_hvg, resolution=resolution, key_added="clusters",flavor="igraph",n_iterations=4)
    return adata_hvg


def run_umap(adata_hvg):
    sc.tl.umap(adata_hvg)
    return adata_hvg