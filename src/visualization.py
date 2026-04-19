import scanpy as sc

def plot_violin(adata, genes, groupby="clusters",save=None):
    
    sc.pl.violin(
        adata,
        keys=genes,
        groupby=groupby,
        stripplot=False
    )
    if save:
        import matplotlib.pyplot as plt
        plt.savefig(save, bbox_inches="tight")
        plt.close()

def plot_umap(adata, color="clusters",save=None):
    sc.pl.umap(adata, color=color)
    if save:
        import matplotlib.pyplot as plt
        plt.savefig(save, bbox_inches="tight")
        plt.close()
