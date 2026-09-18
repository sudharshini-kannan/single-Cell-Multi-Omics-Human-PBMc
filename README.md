# 🧬 Single-Cell Multi-Omics Analysis: Human PBMC RNA & ATAC

A hands-on bioinformatics project using real **10x Genomics PBMC multiome data** to explore single-cell RNA expression and ATAC peak annotations.

## 🎯 Objectives

* Analyze single-cell RNA-seq data using **Scanpy**
* Explore ATAC peak annotations and genomic regions
* Compare RNA genes with ATAC-associated genes
* Perform QC, normalization, PCA, UMAP and clustering
* Identify marker genes and perform tentative cell-type annotation
* Build a lightweight, reproducible Python workflow

## 🧪 Dataset

**PBMC from a healthy donor with granulocytes removed (3k)** from 10x Genomics.

Current downloaded data:

* RNA expression matrix
* ATAC peak coordinates
* ATAC peak-to-gene annotations

**RNA:** 2,711 cells × 36,601 features
**Filtered RNA:** 2,645 cells
**Highly variable genes:** 2,000
**Leiden clusters:** 10

> ⚠️ The current dataset does not contain the cell × ATAC peak accessibility matrix. Therefore, this project currently focuses on **single-cell RNA analysis and gene-level RNA–ATAC exploration**, rather than full cell-level multiome integration.

## 🧬 RNA-seq Workflow

```text
RNA Matrix
   ↓
Quality Control
   ↓
Filtering
   ↓
Normalization
   ↓
Highly Variable Genes
   ↓
PCA
   ↓
UMAP
   ↓
Leiden Clustering
   ↓
Marker Genes
   ↓
Cell-Type Annotation
```

### QC filtering

```text
n_genes_by_counts > 200
n_genes_by_counts < 6000
pct_counts_mt < 20
```

**2,645 cells retained (97.57%)**

## 🔗 RNA–ATAC Exploration

RNA genes were compared with genes associated with ATAC peak annotations.

| Category              | Result |
| --------------------- | -----: |
| RNA genes             | 36,591 |
| ATAC-associated genes | 27,844 |
| Shared genes          | 27,844 |

ATAC annotations included **promoter, distal and intergenic** regions.

## 📊 Outputs

### Figures

```text
figures/
├── QC plots
├── PCA / UMAP plots
├── cluster plots
├── marker gene plots
└── cell-type UMAP
```

### Results

```text
results/
├── highly_variable_genes.csv
├── rna_atac_shared_genes.csv
├── rna_cell_type_annotations.csv
├── rna_clusters.csv
├── rna_marker_genes.csv
├── rna_pca_coordinates.csv
├── rna_qc_metrics.csv
└── rna_umap_coordinates.csv
```

## 📁 Project Structure

```text
Single-Cell-Multiomics/
├── data/raw/
├── figures/
├── results/
├── src/
│   ├── real_multiome.py
│   └── rna_analysis.py
├── requirements.txt
├── test_installation.py
└── README.md
```

## 🛠️ Tools

**Python · Scanpy · AnnData · Muon · MuData · pandas · NumPy · Bash · Linux · WSL · Git · GitHub**

## 🚀 Run

```bash
git clone https://github.com/sudharshini-kannan/single-Cell-Multi-Omics-Human-PBMc.git
cd single-Cell-Multi-Omics-Human-PBMc

python3 -m venv multiomics_env
source multiomics_env/bin/activate
pip install -r requirements.txt

python src/real_multiome.py
python src/rna_analysis.py
```

## 🔮 Future Work

* Cell-level ATAC analysis
* ATAC QC and peak filtering
* TF-IDF + LSI
* MuData construction
* RNA–ATAC integration
* Regulatory and promoter analysis

## 👩‍💻 Author

**Sudharshini Kannan**
*Molecular Biology Researcher | Bioinformatics & Genomics Enthusiast*

🔗 https://github.com/sudharshini-kannan
