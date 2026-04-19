🧬 Single-cell RNA-seq Analysis Toolkit (Python)

A modular and reproducible single-cell RNA-seq analysis pipeline implemented in Python using Scanpy. This set of codes is a starting point for your single-cell data obtained from the sequencer after running on cellranger (once you obtained your gene expression matrix)


⚙️ Installation
conda create -n scrna-toolkit python=3.10
conda activate scrna-toolkit
pip install -r requirements.txt

🧪 Dataset

Tested using: PBMC 3k dataset from 10x Genomics
Can be used with your in-house generated data by replacing with your own data in the data/ folder

🚀 Functionality provides by the module:
Preprocessing and QC (including mitochondrial filtering)
Cell cycle scoring
Highly variable gene selection
PCA, neighbor graph construction
Leiden clustering
UMAP visualization
Differential gene expression (DEG) analysis
Marker-based and DEG-based cluster annotation
Export of results and plots

🧠 Methods
Clustering: Leiden algorithm
DEG: Wilcoxon rank-sum test
Annotation:
Marker-based scoring
DEG-based overlap with canonical markers



▶️ Run Pipeline
python main.py
📊 Outputs
results/deg_filtered.csv → high-confidence marker genes
results/cell_annotations.csv → cluster annotations
figures/umap_*.png → visualizations

📌 Upcoming
Reference-based annotation
Integration with spatial transcriptomics
User-friendly app interface