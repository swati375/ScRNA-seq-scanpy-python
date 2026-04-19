from src.preprocessing import load_data, preprocess,visualize
from src.preprocessing import add_cell_cycle_scores
from src.clustering import run_neighbors,run_dimensionality_reduction,run_clustering,run_umap
from src.visualization import plot_umap,plot_violin
from src.deg import compute_deg
from src.deg import get_top_genes
# from src.annotation import annotate_clusters, annotate_clusters_deg
from src.annotation_new import annotate_clusters_deg_filtered
import pandas as pd

DATA_PATH = "data/filtered_gene_bc_matrices/hg19/"

def main():
    print("Loading data...")
    adata = load_data(DATA_PATH)

    print("Before filtering:", adata.shape)
    print("Preprocessing...")
    adata = preprocess(adata)
    print("After filtering:", adata.shape)

    # visualize(adata)
    
    adata=add_cell_cycle_scores(adata)
    print("Added cell cycle scores")
    # print(adata.obs[["S_score", "G2M_score", "phase"]].head())

    ## to regress out s/g2m score
    # sc.pp.regress_out(adata, ["S_score", "G2M_score"])

    # # Save processed data
    # adata.write("data/pbmc3k_processed.h5ad")

    ##clustering
    adata,adata_hvg = run_dimensionality_reduction(adata)
    adata_hvg = run_neighbors(adata_hvg)# run_neighbors(adata, n_neighbors=10, n_pcs=30)
    #leiden clustering
    adata_hvg = run_clustering(adata_hvg)# run_clustering(adata, resolution=0.5)
    adata_hvg = run_umap(adata_hvg)
    

    print("Clustering done")
    adata.obs["clusters"]=adata_hvg.obs["clusters"]
    adata.obsm["X_umap"]=adata_hvg.obsm["X_umap"]
    # plot_umap(adata,"clusters", save="figures/umap_clusters.png")

    print("\nCluster sizes:")
    print(adata.obs["clusters"].value_counts())

    adata.write("results/pbmc3k_processed.h5ad")

    ##DEG analysis
    adata=compute_deg(adata)
    top_genes = get_top_genes(adata)
    print("Top marker genes per cluster:")
    for cluster, genes in top_genes.items():
        print(cluster, ":", genes)

    # plot some genes
    plot_violin(adata, ["MS4A1", "CD74", "LYZ","CD3D"],
    save="figures/violin_markers.png")

    ##cluster annotation
    #------------------------------------------------------------------
    #----------------------------annotation.py--------------------------
    #------------------------------------------------------------------
    # ##Annotation a. Marker based
    # cluster_labels = annotate_clusters(adata)

    # print("\nCluster Annotations by finding set genes:")
    # for cl, label in cluster_labels.items():
    #     print(f"Cluster {cl}: {label}")

    # adata.obs["cell_type"] = adata.obs["clusters"].map(cluster_labels)
    # plot_umap(adata,"cell_type", save="figures/umap_celltype.png")

    # ##Annotation b. DEG based and then checked with markers 
    # cluster_labels = annotate_clusters_deg(adata)
    # print("\nDEG-based annotations:")
    # for cl, label in cluster_labels.items():
    #     print(f"Cluster {cl}: {label}")

    # adata.obs["cell_type_deg"] = adata.obs["clusters"].map(cluster_labels)
    # plot_umap(adata,"cell_type_deg", save="figures/umap_celltype_deg.png")
    # adata.write("results/pbmc3k_processed_annotated.h5ad")

    # df = adata.obs[["clusters", "cell_type", "cell_type_deg"]]
    # df.to_csv("results/cell_annotations.csv")

#------------------------------------------------------------------
#----------------------------annotation_new.py--------------------------
#------------------------------------------------------------------
    cluster_labels_deg, deg_df, filtered_deg = annotate_clusters_deg_filtered(adata,1.5)

    # print
    print("\nDEG-based annotations:")
    for cl, label in cluster_labels_deg.items():
        print(f"Cluster {cl}: {label}")

    # map to adata
    adata.obs["cell_type_deg"] = adata.obs["clusters"].map(cluster_labels_deg)

    # save DEG tables
    deg_df.to_csv("results/deg_full.csv", index=False)
    filtered_deg.to_csv("results/deg_filtered.csv", index=False)

    # plot
    plot_umap(adata, "cell_type_deg", save="figures/umap_celltype_deg.png")
    adata.write("results/pbmc3k_processed_annotated.h5ad")


if __name__ == "__main__":
    main()