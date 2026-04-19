import scanpy as sc

def load_data(data_path):
    adata = sc.read_10x_mtx(
        data_path,
        var_names="gene_symbols",
        cache=True
    )
    return adata


def preprocess(adata, min_genes=200, min_cells=3, mt_threshold=10):
    # Identify mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith(("MT-","mt-"))

    # Compute QC metrics
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

    # Filter cells
    adata = adata[adata.obs.n_genes_by_counts > min_genes, :]
    adata = adata[adata.obs.pct_counts_mt < mt_threshold, :]

    # Filter genes
    sc.pp.filter_genes(adata, min_cells=min_cells)

    # Normalize
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    return adata

def load_human_cc_genes():## Tirosh et. al 2016
    s_genes = [
        "MCM5","PCNA","TYMS","FEN1","MCM2","MCM4","RRM1","UNG","GINS2","MCM6",
        "CDCA7","DTL","PRIM1","UHRF1","HELLS","RFC2","RPA2","NASP","RAD51AP1"
    ]

    g2m_genes = [
        "HMGB2","CDK1","NUSAP1","UBE2C","BIRC5","TPX2","TOP2A","NDC80",
        "CKS2","NUF2","CKS1B","MKI67","TMPO","CENPF","TACC3"
    ]

    return s_genes, g2m_genes

def add_cell_cycle_scores(adata):

    s_genes, g2m_genes = load_human_cc_genes()

    # cc_genes = sc.queries.cell_cycle_genes()
    # s_genes = [g for g in cc_genes["S"] if g in adata.var_names]
    # g2m_genes = [g for g in cc_genes["G2/M"] if g in adata.var_names]

    s_genes = [g for g in s_genes if g in adata.var_names]
    g2m_genes = [g for g in g2m_genes if g in adata.var_names]

    sc.tl.score_genes_cell_cycle(
        adata,
        s_genes=s_genes,
        g2m_genes=g2m_genes
    )

    return adata    


def visualize(adata):
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True
    )

    