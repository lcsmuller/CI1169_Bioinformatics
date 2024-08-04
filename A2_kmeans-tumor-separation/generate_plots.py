import scanpy as sc
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
import anndata as ad
import matplotlib.pyplot as plt
import infercnvpy as cnv
import pandas as pd

SAMPLES_DIR = './samples'
SAMPLES = ['PT1', 'PT2', 'PT3', 'LN1', 'LN2', 'Li1', 'Li2']
RESULTS_DIR = './results'
GENE_POSITIONS_FILE = './gene_positions.csv'
CELL_MARKERS_FILE = './cell_markers.csv'


def load_and_preprocess_data(samples_dir, samples):
    adatas = []
    for sample in samples:
        adata = sc.read_10x_mtx(f'{samples_dir}/{sample}', var_names='gene_symbols', cache=True)
        adata.obs_names_make_unique()

        # Calculate n_counts and percent.mt
        adata.obs['n_counts'] = adata.X.sum(axis=1).A1
        adata.obs['percent_mt'] = (adata[:, adata.var_names.str.startswith('MT-')].X.sum(axis=1).A1 / adata.obs['n_counts']) * 100

        # Filtering low-quality cells based on the provided thresholds
        sc.pp.filter_cells(adata, min_genes=500)
        adata = adata[adata.obs['n_counts'] > 1000, :]
        adata = adata[adata.obs['n_counts'] < 20000, :]
        adata = adata[adata.obs['percent_mt'] < 10, :]

        sc.pp.filter_genes(adata, min_cells=3)
        adata.raw = adata  # Store raw data
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        adatas.append(adata)

    # Concatenate all AnnData objects
    adata = ad.concat(adatas, label='sample', keys=samples, index_unique='-')
    return adata, adatas


def annotate_clusters(adata, markers_file):
    markers = pd.read_csv(markers_file)
    for cell_type in markers['cell_type'].unique():
        cell_type_genes = markers[markers['cell_type'] == cell_type]['gene'].values
        existing_genes = [gene for gene in cell_type_genes if gene in adata.var_names]
        if existing_genes:
            sc.tl.score_genes(adata, existing_genes, score_name=cell_type)
            adata.obs.loc[adata.obs[cell_type] > 0, 'cell_type'] = cell_type
        else:
            print(f"WARNING: genes for {cell_type} are not in var_names and ignored.")
    return adata[adata.obs['cell_type'].notna()]


def load_gene_positions(gene_positions_file, adata):
    gene_positions = pd.read_csv(gene_positions_file)
    gene_positions = gene_positions.drop_duplicates(subset='gene_symbol')
    adata.var = adata.var.merge(gene_positions, left_index=True, right_on='gene_symbol', how='left')
    adata.var.set_index('gene_symbol', inplace=True)

    if not {'chromosome', 'start', 'end'}.issubset(adata.var.columns):
        raise ValueError("Gene position information is incomplete. Ensure 'chromosome', 'start', and 'end' columns are present.")
    return adata


def perform_cnv_analysis(adata):
    reference_cells = adata[adata.obs['cell_type'].isin(['Myeloid'])].obs.index.tolist()
    adata.obs['cell_kind'] = 'tumor'
    adata.obs.loc[reference_cells, 'cell_kind'] = 'normal'

    cnv.tl.infercnv(
        adata,
        reference_key='cell_kind',
        reference_cat='normal',
        window_size=200,
        step=100
    )

    adata.obs['cnv_score'] = adata.obsm['X_cnv'].mean(axis=1)

    cnv_scores = adata.obs['cnv_score'].values.reshape(-1, 1)
    kmeans_cnv = KMeans(n_clusters=2, random_state=0).fit(cnv_scores)
    adata.obs['cnv_cluster'] = kmeans_cnv.labels_

    cluster_means = adata.obs.groupby('cnv_cluster')['cnv_score'].mean()
    malignant_cluster = cluster_means.idxmax()

    adata.obs['malign'] = adata.obs['cnv_cluster'].apply(lambda x: 'Malignant' if x == malignant_cluster else 'Non-malignant')
    return adata


def plot_results(adata, results_dir):
    sc.pl.umap(adata, color='sample', title='UMAP colored by sample', show=False)
    plt.savefig(f'{results_dir}/umap_samples.png', bbox_inches='tight')

    sc.pl.umap(adata, color='kmeans', title='UMAP colored by K-means clusters', show=False, palette=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b'])
    plt.savefig(f'{results_dir}/umap_clusters.png', bbox_inches='tight')

    sc.pl.umap(adata, color='cell_type', title='UMAP colored by cell types', show=False)
    plt.savefig(f'{results_dir}/umap_celltypes.png', bbox_inches='tight')

    epithelial_cells = adata[adata.obs['cell_type'] == 'Epithelial'].copy()
    sc.pl.umap(epithelial_cells, color='malign', title='UMAP colored by Malignant and Non-malignant cells (Epithelial)', show=False)
    plt.savefig(f'{results_dir}/umap_epithelial_malign_nonmalign.png', bbox_inches='tight')


def basic_analysis(adata):
    celltype_counts = adata.obs['cell_type'].value_counts()
    print("Distribution of cell types:")
    print(celltype_counts)

    malign_counts = adata.obs.groupby(['cell_type', 'malign'], observed=True).size().unstack(fill_value=0)
    malign_counts['Total'] = malign_counts.sum(axis=1)
    malign_counts['Percent Malignant'] = 100 * malign_counts['Malignant'] / malign_counts['Total']
    print("\nMalignant and Non-malignant cell counts in each K-means cluster:")
    print(malign_counts)


if __name__ == '__main__':
    adata, adatas = load_and_preprocess_data(SAMPLES_DIR, SAMPLES)

    adata = annotate_clusters(adata, CELL_MARKERS_FILE)

    adata_raw = ad.concat([a.raw.to_adata() for a in adatas], label='sample', keys=SAMPLES, index_unique='-')
    sc.pp.highly_variable_genes(adata_raw, n_top_genes=2000, flavor='seurat_v3')
    adata = adata[:, adata_raw.var.highly_variable].copy()

    sc.pp.scale(adata, max_value=10)

    adata = load_gene_positions(GENE_POSITIONS_FILE, adata)

    sc.pp.pca(adata, n_comps=50)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=20)
    sc.tl.umap(adata, min_dist=0.3, spread=1.0)

    X = adata.X
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    kmeans = KMeans(n_clusters=19)
    kmeans.fit(X_scaled)
    clusters = kmeans.labels_
    adata.obs['kmeans'] = clusters + 1

    adata = perform_cnv_analysis(adata)
    plot_results(adata, RESULTS_DIR)
    basic_analysis(adata)
