import scanpy as sc
import pandas as pd

# Load RNA data
print("Loading RNA data...")

rna = sc.read_10x_h5(
    "data/raw/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix.h5"
)

# Make gene names unique
rna.var_names_make_unique()

print("\n===== RNA DATA =====")
print(rna)
print("Cells:", rna.n_obs)
print("Genes:", rna.n_vars)

# Calculate quality-control metrics
# Identify mitochondrial genes
rna.var["mt"] = rna.var_names.str.startswith("MT-")

# Calculate QC metrics
sc.pp.calculate_qc_metrics(
    rna,
    qc_vars=["mt"],
    inplace=True
)

print("\n===== QC METRICS =====")
print(
    rna.obs[
        [
            "total_counts",
            "n_genes_by_counts",
            "pct_counts_mt"
        ]
    ].describe()
)

# Save QC information
rna.obs[
    [
        "total_counts",
        "n_genes_by_counts",
        "pct_counts_mt"
    ]
].to_csv(
    "results/rna_qc_metrics.csv"
)

print("\nQC metrics saved to results/rna_qc_metrics.csv")

# Plot QC distributions
sc.pl.violin(
    rna,
    ["n_genes_by_counts", "total_counts"],
    jitter=0.4,
    multi_panel=True,
    show=False
)

sc.pl.scatter(
    rna,
    x="total_counts",
    y="n_genes_by_counts",
    show=False
)

# Save plots
sc.pl.violin(
    rna,
    ["n_genes_by_counts", "total_counts"],
    jitter=0.4,
    multi_panel=True,
    save="_qc_violin.png"
)

sc.pl.scatter(
    rna,
    x="total_counts",
    y="n_genes_by_counts",
    save="_qc_scatter.png"
)

print("\n===== QC METRICS =====")

print(
    rna.obs[
        [
            "total_counts",
            "n_genes_by_counts"
        ]
    ].describe()
)

# Save basic QC information
rna.obs[
    [
        "total_counts",
        "n_genes_by_counts"
    ]
].to_csv(
    "results/rna_qc_metrics.csv"
)

print("\nQC metrics saved to results/rna_qc_metrics.csv")

# Filter low-quality cells
print("\n===== FILTERING CELLS =====")

before_filtering = rna.n_obs

rna = rna[
    (rna.obs["n_genes_by_counts"] > 200)
    & (rna.obs["n_genes_by_counts"] < 6000)
    & (rna.obs["pct_counts_mt"] < 20)
].copy()

after_filtering = rna.n_obs

print("Cells before filtering:", before_filtering)
print("Cells after filtering:", after_filtering)
print("Cells removed:", before_filtering - after_filtering)

print(
    "Percentage retained:",
    round((after_filtering / before_filtering) * 100, 2),
    "%"
)

# Preserve raw counts before normalization
rna.layers["counts"] = rna.X.copy()

# Normalize total counts per cell
print("\n===== NORMALIZATION =====")

sc.pp.normalize_total(
    rna,
    target_sum=1e4
)

# Log-transform normalized values
sc.pp.log1p(rna)

# Preserve log-normalized expression for marker analysis
rna.layers["lognorm"] = rna.X.copy()

print("Normalization complete.")
print("Target counts per cell: 10,000")

# Identify highly variable genes
print("\n===== HIGHLY VARIABLE GENES =====")

sc.pp.highly_variable_genes(
    rna,
    n_top_genes=2000,
    flavor="seurat"
)

print(
    "Highly variable genes:",
    rna.var["highly_variable"].sum()
)

# Save highly variable gene information
rna.var[
    ["highly_variable"]
].to_csv(
    "results/highly_variable_genes.csv"
)

print(
    "HVG information saved to results/highly_variable_genes.csv"
)

# PCA using highly variable genes
print("\n===== PCA =====")

sc.pp.scale(
    rna,
    max_value=10
)

sc.tl.pca(
    rna,
    n_comps=30,
    use_highly_variable=True
)

print("PCA completed.")
print(
    "Number of PCs:",
    rna.obsm["X_pca"].shape[1]
)

# Save PCA coordinates
pd.DataFrame(
    rna.obsm["X_pca"],
    index=rna.obs_names
).to_csv(
    "results/rna_pca_coordinates.csv"
)

print(
    "PCA coordinates saved to results/rna_pca_coordinates.csv"
)

# PCA variance explained
print("\n===== PCA VARIANCE =====")

variance_ratio = rna.uns["pca"]["variance_ratio"]

for i, variance in enumerate(
    variance_ratio,
    start=1
):
    print(
        f"PC{i}: {variance * 100:.2f}%"
    )

print(
    "\nTotal variance explained by 30 PCs:",
    round(
        variance_ratio.sum() * 100,
        2
    ),
    "%"
)

# Build neighborhood graph using PCA
print("\n===== NEIGHBOR GRAPH =====")

sc.pp.neighbors(
    rna,
    n_neighbors=15,
    n_pcs=20
)

print("Neighbor graph completed.")

# Compute UMAP
print("\n===== UMAP =====")

sc.tl.umap(
    rna,
    random_state=42
)

print("UMAP completed.")

# Save UMAP coordinates
pd.DataFrame(
    rna.obsm["X_umap"],
    index=rna.obs_names,
    columns=["UMAP1", "UMAP2"]
).to_csv(
    "results/rna_umap_coordinates.csv"
)

print(
    "UMAP coordinates saved to results/rna_umap_coordinates.csv"
)

# Save UMAP plot
sc.pl.umap(
    rna,
    show=False,
    save="_rna.png"
)

# Save colored UMAP plots
print("\n===== COLORED UMAP PLOTS =====")

# UMAP colored by number of detected genes
sc.pl.umap(
    rna,
    color="n_genes_by_counts",
    show=False,
    save="_n_genes.png"
)

# UMAP colored by mitochondrial percentage
sc.pl.umap(
    rna,
    color="pct_counts_mt",
    show=False,
    save="_mitochondrial.png"
)

print("Colored UMAP plots saved.")

# Cell clustering using Leiden algorithm
print("\n===== CELL CLUSTERING =====")

sc.tl.leiden(
    rna,
    resolution=0.5,
    random_state=42
)

print("Clustering completed.")

print("\nCluster counts:")

print(
    rna.obs[
        "leiden"
    ].value_counts().sort_index()
)

# Save cluster assignments
rna.obs[
    ["leiden"]
].to_csv(
    "results/rna_clusters.csv"
)

# UMAP colored by clusters
sc.pl.umap(
    rna,
    color="leiden",
    palette="tab10",
    legend_loc="on data",
    show=False,
    save="_clusters_colored.png"
)

print(
    "Colored cluster UMAP saved to "
    "figures/umap_clusters_colored.png"
)

# Find marker genes for each cluster
print("\n===== MARKER GENE ANALYSIS =====")

sc.tl.rank_genes_groups(
    rna,
    groupby="leiden",
    mask_var="highly_variable",
    layer="lognorm",
    method="wilcoxon",
    n_genes=20
)

print("Marker gene analysis completed.")

# Save marker genes
markers = rna.uns[
    "rank_genes_groups"
]

marker_df = pd.DataFrame({
    "cluster": markers[
        "names"
    ].dtype.names
})

marker_table = pd.DataFrame(
    markers["names"]
)

marker_table.to_csv(
    "results/rna_marker_genes.csv"
)

print(
    "Marker genes saved to "
    "results/rna_marker_genes.csv"
)

# Plot top marker genes
sc.pl.rank_genes_groups(
    rna,
    n_genes=10,
    sharey=False,
    show=False,
    save="_marker_genes.png"
)

print("Marker gene plot saved.")

# Cell-type annotation
print("\n===== CELL-TYPE ANNOTATION =====")

cell_type_annotations = {
    "0": "Classical Monocytes",
    "1": "Naive/CD4 T cells",
    "2": "Naive T cells",
    "3": "CD8 T cells",
    "4": "T cells",
    "5": "Cytotoxic lymphocytes",
    "6": "B cells",
    "7": "NK cells",
    "8": "Antigen-presenting cells",
    "9": "Rare/Unclear"
}

rna.obs["cell_type"] = (
    rna.obs["leiden"]
    .map(cell_type_annotations)
)

print("\nCell-type counts:")

print(
    rna.obs["cell_type"]
    .value_counts()
)

# Save cell-type annotations
rna.obs[
    ["leiden", "cell_type"]
].to_csv(
    "results/rna_cell_type_annotations.csv"
)

print(
    "\nCell-type annotations saved to "
    "results/rna_cell_type_annotations.csv"
)

# Plot UMAP with cell-type annotations
sc.pl.umap(
    rna,
    color="cell_type",
    legend_loc="on data",
    show=False,
    save="_cell_types.png"
)

print(
    "Cell-type annotated UMAP saved."
)