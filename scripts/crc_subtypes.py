#!/usr/bin/env python
# coding: utf-8

# In[1]:


import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
import glob
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.colors as mcolors
# Inital setting for plot size
from matplotlib import rcParams
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')
def get_top_genes(adata, n_genes=25 ):
    # Get the dictionary of ranked genes
    ranked_genes = adata.uns["rank_genes_groups"]['names']
    # Create a DataFrame to store the top genes for each cluster
    top_genes_df = pd.DataFrame({cluster: ranked_genes[cluster][:n_genes]
                                 for cluster in ranked_genes.dtype.names})

    return top_genes_df


# In[2]:


adata = sc.read('/mnt/beegfs/userdata/d_papakonstantinou/crc/data/integrated_crc.h5ad')
adata.obs['sex'] = adata.obs['sex'].replace({22:'F',2:"F",1:'M'})


# In[3]:

print(adata.obs.columns)
meta = adata.obs[['orig.ident','cohort','age','sex','clinical_trials','grade_of_differentiation','tumor_location','histological_variants','location_of_metastasis','tumor_type','histological_type']].drop_duplicates()
print(len(meta))
meta['age'].plot(kind='hist', bins=20, alpha=0.7, color='blue', edgecolor='black', linewidth=1.2)
plt.title('Age', fontsize=15)
plt.ylabel('Patients', fontsize=12)
plt.grid(False)
plt.show()
meta['tumor_location'].value_counts().plot(kind='bar', alpha=0.8, color='blue', edgecolor='black', linewidth=1.2)
plt.title('Location', fontsize=15)
plt.xlabel('')
plt.xticks(rotation= 45,ha='right', fontsize =9 )
plt.grid(False)
plt.show()

meta['grade_of_differentiation'].value_counts().plot(kind='bar', alpha=0.8, color='blue', edgecolor='black', linewidth=1.2)
plt.title('differentiation', fontsize=15)
plt.xlabel('')
plt.grid(False)
plt.show()

meta['histological_variants'].value_counts().plot(kind='bar', alpha=0.8, color='blue', edgecolor='black', linewidth=1.2)
plt.title('variants', fontsize=15)
plt.xlabel('')
plt.xticks(rotation= 45,ha='right', fontsize =9 )
plt.grid(False)
plt.show()

meta['location_of_metastasis'].value_counts().plot(kind='bar', alpha=0.8, color='blue', edgecolor='black', linewidth=1.2)
plt.title('metastasis', fontsize=15)
plt.xlabel('')
plt.grid(False)
plt.show()

meta['sex'].value_counts().plot(kind='bar', alpha=0.8, color='blue', edgecolor='black', linewidth=1.2)
plt.title('sex', fontsize=15)
plt.xlabel('')
plt.grid(False)
plt.show()


# In[4]:


adata


# In[5]:


ncols = 2
nrows = 1
figsize = 5
wspace = 0.5
fig, axs = plt.subplots(
    nrows=1,
    ncols=2,
    figsize=(ncols * figsize + figsize * wspace * (ncols - 1), nrows * figsize),
)
plt.subplots_adjust(wspace=wspace)
# This produces two Axes objects in a single Figure
print("axes:", axs)
sc.pl.umap(adata,color="cell_type",legend_loc="on data",frameon=False,legend_fontsize=8,legend_fontoutline=2,ax=axs[0], show=False)
sc.pl.umap(adata,color="leiden",legend_loc="on data",frameon=False,legend_fontsize=10,legend_fontoutline=2,ax=axs[1])


# In[6]:


def freq_(x,a):
    return round((x*100)/a,2)

def stacked_barplot(ax,metadata,title):
    cell_counts = metadata.groupby(["cell_type"]).apply(len)
    freqs = cell_counts.apply(freq_,a= sum(cell_counts))
    bottom = np.zeros(1)
    width = 1
    for type_, freq in freqs.items():
        p = ax.bar('freqs',freq, width=0.1, label=type_, bottom=bottom)
        bottom += freq

    ax.set_title(title,fontsize= 10, y=0.95,  loc='center', va='baseline')

    return ax, freqs

fig, axs = plt.subplots(1, 5, figsize=(10, 10))
axs[0],_ = stacked_barplot(axs[0],adata.obs,"total")
axs[0].axis('off')

for i, variant in enumerate(['Mucinous','Adenoma-like (mucinous 5%)','Mixed (NOS+Muc)','NOS']):
    sub = adata.obs[adata.obs['histological_variants'] == variant]
    variant = 'mucinous 5%' if variant == 'Adenoma-like (mucinous 5%)' else variant
    axs[i +1], freqs = stacked_barplot(axs[i+1],sub,variant)
    axs[i +1].axis('off')

order = ['T cell', 'Monocyte', 'Fibroblast', 'Epithelial', 'Endothelial', 'B cell','Plasma cell']
handles, labels = plt.gca().get_legend_handles_labels()
ordered_handles = [handles[labels.index(item)] for item in order]
ordered_labels = [item for item in order]
plt.legend(ordered_handles, ordered_labels,bbox_to_anchor=(1.05, 1.0))
plt.grid(False)
plt.axis('off')
plt.tight_layout()
plt.show()


# In[7]:


subset = adata[adata.obs['cell_type'].isin(['Monocyte']), :].copy()

#sc.pl.highest_expr_genes(subset, n_top=20)
subset.uns['log1p']["base"] = None
#print(subset.X)
#sc.pp.normalize_total(subset, target_sum=1e4)
#sc.pp.highly_variable_genes(subset, min_mean=0.0125, max_mean=3, min_disp=0.5)
#sc.pp.highly_variable_genes(subset, n_bins = 20)
#sc.pl.highly_variable_genes(subset)
#subset = subset[:, subset.var.highly_variable]
subset.raw = subset
#sc.pp.regress_out(subset, ["total_counts", "pct_counts_mt"])
#sc.pp.scale(subset, max_value=10)

sc.tl.pca(subset, svd_solver="arpack")
sc.pl.pca_variance_ratio(subset, log=True)
sc.pp.neighbors(subset, n_neighbors=10, n_pcs=50)
sc.tl.umap(subset)
sc.external.pp.bbknn(subset, batch_key='orig.ident',neighbors_within_batch =1)
sc.tl.umap(subset)
sc.tl.leiden(subset,resolution = 0.5,key_added = 'clusters')
sc.pl.umap(subset, color=['histological_variants', 'TREM2'])

sc.tl.rank_genes_groups(subset, groupby="clusters", method="wilcoxon")

#sc.pl.heatmap(subset, marker_genes_dict, groupby="clusters" ,cmap="RdBu_r", dendrogram=True, swap_axes=True, figsize=(11, 4))
sc.pl.rank_genes_groups(subset,groupby = 'clusters', n_genes=25, sharey=False)


# In[8]:


plt.rcParams['font.size'] = 6
sc.pl.rank_genes_groups_heatmap(subset, n_genes=10, use_raw=False, swap_axes=True, show_gene_labels=True, vmin=-3, vmax=3, cmap="bwr",  figsize=(10, 7), show=False)
#sc.pl.rank_genes_groups_heatmap(subset, n_genes=10, use_raw=True, swap_axes=True, show_gene_labels=False, vmin=-3, vmax=3, cmap="bwr")


# In[9]:


top_genes = get_top_genes(subset, n_genes=25)
for i in top_genes.columns:
    gene_markers =top_genes[i].to_list()  # Replace with your actual gene list
    print(i,':', ', '.join(gene_markers))
    # put results to http://xteam.xbio.top/ACT/ResultAction.action?jobID=20240628063549UAWQVFE7YRSQ15


cluster_annotations = {
    "0":"Classical monocyte",
    "1":"Tumor-associated macrophage",
    "2":"Unknown",
    "3":"Dendritic cell",
    "4":"Plasmacytoid dendritic cell"
}
subset.obs["cell_type"] = subset.obs["clusters"].map(cluster_annotations).astype("category")

sc.settings.set_figure_params(dpi=80, facecolor='white')
#plt.axes('off')
sc.pl.umap(subset, color='cell_type',legend_loc="on data",frameon=False,legend_fontsize=8,legend_fontoutline=2,)
sc.pl.umap(subset, color=['TREM2','histological_variants'],cmap="YlGn")

classical_monocyte =  ['TIMP1', 'S100A8', 'IL1B', 'S100A9', 'SOD2', 'BCL2A1', 'PLAUR', 'VCAN', 'EREG', 'FCN1', 'SLC11A1', 'SERPINA1', 'G0S2', 'ETS2', 'CD300E', 'SLC25A37', 'IL1RN', 'OLR1', 'VEGFA', 'PPIF', 'TREM1', 'C15orf48', 'PTGS2', 'ACSL1', 'SAT1']
TAM = ['C1QA', 'C1QB', 'CTSD', 'C1QC', 'CTSB', 'CTSZ', 'APOE', 'TMEM176B', 'RNASE1', 'SELENOP', 'GRN', 'APOC1', 'GPNMB', 'CAPG', 'CD14', 'SLC40A1', 'SLC7A8', 'INF2', 'DSC2', 'GCNT1', 'MUC5B', 'NHSL1', 'HLA-DRA', 'EIF5AL1', 'TPTE2P5']
Unknown = ['PLAAT2', 'IGHV1-3', 'IGLL5', 'STBD1', 'AL139383.1', 'GAS1RR', 'HCG27', 'AC138028.4', 'ZNF608', 'IGLV2-8', 'CD79A', 'IGHV7-4-1', 'UCHL1', 'PLEKHN1', 'IGLV4-69', 'IGLV2-23', 'HTRA2', 'PAEP', 'IGLV7-43', 'STIM2-AS1', 'AC245060.4', 'CCNI2', 'IGHV3-72', 'EEF2KMT', 'FAM149A']
Dendritic= ['HLA-DPB1', 'HLA-DPA1', 'HLA-DQB1', 'HLA-DRB1', 'NAPSB', 'HLA-DQA1', 'HLA-DRA', 'CST7', 'CST3', 'LGALS2', 'LTB', 'TUBA1B', 'CCND1', 'PADI2', 'RMI2', 'AL109615.3', 'C12orf75', 'ADGRG6', 'DTL', 'CCSER1', 'NR2C2AP', 'VRK1', 'IRF6', 'LINC02195', 'AMOT']
PlasmacytoidDendritic = ['GZMB', 'PLAC8', 'IRF4', 'CLIC3', 'BCL11A', 'TSPAN13', 'IL3RA', 'C12orf75', 'SOX4', 'NAPSB', 'IRF8', 'ALOX5AP', 'SERPINF1', 'SCT', 'SLC15A4', 'MZB1', 'AREG', 'PPP1R14A', 'LTB', 'TCL1A', 'TPM2', 'SMPD3', 'STMN1', 'SIDT1', 'SMIM20']
subset.obs['classical_monocyte'] = np.mean([subset.X[:, subset.var_names == gene].sum(axis=1) for gene in classical_monocyte], axis=0)
subset.obs['TAM'] = np.mean([subset.X[:, subset.var_names == gene].sum(axis=1) for gene in TAM], axis=0)
subset.obs['Unknown'] = np.mean([subset.X[:, subset.var_names == gene].sum(axis=1) for gene in Unknown], axis=0)
subset.obs['Dendritic'] = np.mean([subset.X[:, subset.var_names == gene].sum(axis=1) for gene in Dendritic], axis=0)
subset.obs['PlasmacytoidDendritic'] = np.mean([subset.X[:, subset.var_names == gene].sum(axis=1) for gene in PlasmacytoidDendritic], axis=0)
# Plot the UMAP with aggregate expression
sc.pl.umap(subset, color=['classical_monocyte','TAM','Unknown','Dendritic','PlasmacytoidDendritic'],
           legend_loc="on data",
           frameon=False,
           ncols=5, cmap="YlGn"
          )


# In[10]:


marker_genes_dict = {
    "TREM2":["TREM2","APOC1","APOE","SPP1","FABP5","LGALS3","CD9","LIPA"],
    "HES1":["HES1","C1QA","C1QB","C1QC","CCL3","CCL4","HBEGF","APOE","LIPA","CCL18"],
    "IL4I1":["IL4I1","IDO1","CXCL9","CXCL10","CXCL11","CD40","CCL8","LAMP3"],
    "RGS1":["RGS1",'HLA-DRB1','HLA-DPB1'],
    "IL1B":['VEGFA','CCL20','IL1RN','CXCL2','CXCL3','EREG','PLAUR','CCL4','TNFAIP3','CD44'],
    "DC2/3":['FCER1A','CLEC10A','CD1C','CD1E','CD74'],
    "C1Q":['C1Q','C1QA','C1QB','C1QB','C1QC','HLA-DRB1'],
    'STMN1':['STMN1','TUBB','TOP2A','HMGB2','CDK1','MKI67'],
    'S100A8':['S100A8','S100A9','S100A12','VCAN','S100A6','FOS','LYZ','HMGB2','KLF6']
}
non_existent = ['C1Q', 'CD44', 'CD74', 'FOS', 'HMGB2', 'KLF6', 'RGS1', 'S100A6', 'TNFAIP3', 'TUBB']

for key in marker_genes_dict:
    marker_genes_dict[key] = [item for item in marker_genes_dict[key] if item not in non_existent]
    subset.obs[key.lower()] = np.mean([subset.X[:, subset.var_names == gene].sum(axis=1) for gene in marker_genes_dict[key]], axis=0)

sc.pl.umap(subset, color=[x.lower() for x in marker_genes_dict.keys()],
           legend_loc="on data",
           frameon=False,
           ncols=5, cmap="YlGn"
          )


# In[69]:


normalize = lambda x: (x - x.min()) / (x.max() - x.min())
import numpy as np
from sklearn.metrics import mutual_info_score
import warnings
from scipy.special import kl_div
from sklearn.metrics.pairwise import cosine_similarity
from scipy.spatial.distance import hamming
def kl_divergence_binary(p, q):
    # Convert to numpy arrays
    p = np.array(p)
    q = np.array(q)

    # Check if inputs are binary
    if not set(p).issubset({0, 1}) or not set(q).issubset({0, 1}):
        raise ValueError("Inputs must be binary (0s and 1s only)")

    # Calculate probabilities
    p_1 = np.mean(p)
    p_0 = 1 - p_1
    q_1 = np.mean(q)
    q_0 = 1 - q_1

    # Avoid log(0) by adding a small epsilon
    epsilon = 1e-10

    # Calculate KL divergence
    kl = p_0 * np.log((p_0 + epsilon) / (q_0 + epsilon)) + p_1 * np.log((p_1 + epsilon) / (q_1 + epsilon))

    return kl

def mutual_information(x, y):
    x = np.array(x).astype(int)
    y = np.array(y).astype(int)

    if not set(x).issubset({0, 1}) or not set(y).issubset({0, 1}):
        raise ValueError("Inputs must be binary (0s and 1s only)")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mi = mutual_info_score(x, y)
    return mi

def entropy(x):
    p1 = np.mean(x)
    p0 = 1 - p1
    return -p0 * np.log2(p0) - p1 * np.log2(p1) if 0 < p0 < 1 else 0

def cross_entropy_loss(y_true, y_pred):

    y_true = np.array(y_true)
    y_pred = np.array(y_pred)

    # Check input validity
    if not set(y_true).issubset({0, 1}):
        raise ValueError("y_true must contain only binary values (0 or 1)")
    if np.any((y_pred < 0) | (y_pred > 1)):
        raise ValueError("y_pred must contain values between 0 and 1")

    # Avoid log(0) by clipping predictions
    epsilon = 1e-15
    y_pred = np.clip(y_pred, epsilon, 1 - epsilon)

    # Calculate cross-entropy loss
    loss = -np.mean(y_true * np.log(y_pred) + (1 - y_true) * np.log(1 - y_pred))

    return loss

def hamming_distance(x, y):
    """Calculate the Hamming distance between two binary vectors."""
    x, y = np.array(x), np.array(y)
    if not (set(x).issubset({0, 1}) and set(y).issubset({0, 1})):
        raise ValueError("Inputs must be binary (0s and 1s only)")
    return hamming(x, y) * len(x)

def jaccard_similarity(x, y):
    """Calculate the Jaccard similarity between two binary vectors."""
    x, y = np.array(x), np.array(y)
    if not (set(x).issubset({0, 1}) and set(y).issubset({0, 1})):
        raise ValueError("Inputs must be binary (0s and 1s only)")
    intersection = np.logical_and(x, y)
    union = np.logical_or(x, y)
    return np.sum(intersection) / np.sum(union)

def cosine_similarity_binary(x, y):
    """Calculate the cosine similarity between two binary vectors."""
    x, y = np.array(x).reshape(1, -1), np.array(y).reshape(1, -1)
    if not (set(x[0]).issubset({0, 1}) and set(y[0]).issubset({0, 1})):
        raise ValueError("Inputs must be binary (0s and 1s only)")
    return cosine_similarity(x, y)[0][0]

def phi_coefficient(x, y):
    """Calculate the phi coefficient (correlation) between two binary vectors."""
    x, y = np.array(x), np.array(y)
    if not (set(x).issubset({0, 1}) and set(y).issubset({0, 1})):
        raise ValueError("Inputs must be binary (0s and 1s only)")
    n = len(x)
    sum_x = np.sum(x)
    sum_y = np.sum(y)
    sum_xy = np.sum(x * y)

    numerator = n * sum_xy - sum_x * sum_y
    denominator = np.sqrt((n * sum_x - sum_x**2) * (n * sum_y - sum_y**2))

    if denominator == 0:
        return 0  # To handle the case where one or both vectors are all 0s or all 1s
    return numerator / denominator
res =  []
for histo in set(subset.obs['histological_variants']):
    histo_cells = [1  if histo == row['histological_variants'] else 0 for index,row in subset.obs.iterrows()]
    for name in marker_genes_dict.keys():
        array = subset.obs[name.lower()]
        binarized_expression = [1 if x > 0.5 else 0 for x in normalize(array)]
        mi = mutual_information(histo_cells, binarized_expression)
        h_x = entropy(histo_cells)
        h_y = entropy(binarized_expression)
        kl = kl_divergence_binary(histo_cells, binarized_expression)
        ce = cross_entropy_loss(histo_cells, binarized_expression)
        pc = phi_coefficient(histo_cells,binarized_expression)
        cs = cosine_similarity_binary(histo_cells,binarized_expression)
        hd = hamming_distance(histo_cells,binarized_expression)
        js = jaccard_similarity(histo_cells,binarized_expression)
        d = {"Method":"%s - %s"%(histo,name),"Kullback-Leibler Divergence":kl,"Mutual Information":mi / np.sqrt(h_x * h_y),"Cross-Entropy":ce,'Phi coefficient':pc,'Cosine Similarity':cs,'Hamming Distance':hd,'Jaccard Similarity':js}
        res.append(d)

df = pd.DataFrame(res)
display(df)


# In[74]:


import numpy as np
from scipy import stats
import pandas as pd

def kruskal_wallis_test(values, categories):
    """
    Perform Kruskal-Wallis H-test on continuous values grouped by categories.

    Args:
    values (list or np.array): Continuous values
    categories (list or np.array): Corresponding categories for each value

    Returns:
    dict: A dictionary containing the test statistic, p-value, and a summary DataFrame
    """
    # Convert inputs to numpy arrays
    values = np.array(values)
    categories = np.array(categories)

    # Perform Kruskal-Wallis H-test
    h_statistic, p_value = stats.kruskal(*[values[categories == cat] for cat in np.unique(categories)])

    # Create a summary DataFrame
    summary = pd.DataFrame({
        'Category': np.unique(categories),
        'Count': [np.sum(categories == cat) for cat in np.unique(categories)],
        'Mean Rank': [np.mean(stats.rankdata(values)[categories == cat]) for cat in np.unique(categories)]
    })

    return {
        'H-statistic': h_statistic,
        'p-value': p_value,
        'summary': summary
    }

def create_boxplot(values, categories, test_result):
    """
    Create an enhanced boxplot with whiskers for each category, including p-value and other statistics.

    Args:
    values (list or np.array): Continuous values
    categories (list or np.array): Corresponding categories for each value
    test_result (dict): Result from kruskal_wallis_test function
    """
    plt.figure(figsize=(12, 8))
    sns.set_style("whitegrid")
    sns.set_palette("deep")

    # Create the boxplot
    ax = sns.boxplot(x=categories, y=values, width=0.6)

    # Add strip plot for individual data points
    sns.stripplot(x=categories, y=values, color=".3", size=3, alpha=0.4)

    # Customize the plot
    plt.title('Distribution of Values by Category', fontsize=16, fontweight='bold')
    plt.xlabel('Category', fontsize=12)
    plt.ylabel('Values', fontsize=12)

    # Add p-value annotation
    p_value = test_result['p-value']
    p_value_text = f"p = {p_value:.4f}" if p_value >= 0.0001 else "p < 0.0001"
    plt.text(0.95, 0.95, p_value_text, transform=ax.transAxes,
             verticalalignment='top', horizontalalignment='right',
             fontsize=12, fontweight='bold',
             bbox=dict(facecolor='white', edgecolor='black', alpha=0.8))

    # Add mean values on the plot
    for i, cat in enumerate(np.unique(categories)):
        mean_val = np.mean(values[categories == cat])
        plt.text(i, plt.ylim()[1], f'Mean: {mean_val:.2f}',
                 horizontalalignment='center', verticalalignment='bottom')

    # Add Kruskal-Wallis test statistic
    plt.text(0.05, 0.95, f"H-statistic: {test_result['H-statistic']:.2f}",
             transform=ax.transAxes, verticalalignment='top',
             fontsize=10, bbox=dict(facecolor='white', edgecolor='black', alpha=0.8))

    plt.tight_layout()
    plt.show()



# Example usage
np.random.seed(42)  # for reproducibility

# Generate sample data
categories = np.repeat(['A', 'B', 'C'], [30, 30, 30])
values_A = np.random.normal(loc=5, scale=1, size=30)
values_B = np.random.normal(loc=6, scale=1, size=30)
values_C = np.random.normal(loc=5.5, scale=1, size=30)
values = np.concatenate([values_A, values_B, values_C])

# Perform Kruskal-Wallis test
result = kruskal_wallis_test(values, categories)

# Print results
print(f"Kruskal-Wallis H-statistic: {result['H-statistic']:.4f}")
print(f"p-value: {result['p-value']:.4f}")
print("\nSummary:")
print(result['summary'])

# Interpret the result
alpha = 0.05  # significance level
if result['p-value'] < alpha:
    print(f"\nReject the null hypothesis (p-value < {alpha})")
    print("There are significant differences between the groups.")
else:
    print(f"\nFail to reject the null hypothesis (p-value >= {alpha})")
    print("There is not enough evidence to conclude significant differences between the groups.")

# Create and display the boxplot
create_boxplot(values, categories,result)


# In[12]:


sc.tl.rank_genes_groups(subset, groupby="histological_variants", method="t-test")


# In[13]:


histo_variants = {
'Adenoma_like' : ['HSPA1A', 'LYZ', 'BAG3', 'ZFAND2A', 'DNAJA4', 'SAT1', 'MTRNR2L8', 'THBS1', 'MTCO3P12', 'WDR61', 'HSPA6', 'HLA-DRB1', 'AD000090.1', 'PLIN2', 'MTND2P28', 'DNAJB4', 'AOAH', 'FGL2', 'RNF144B', 'TENT5A', 'MS4A6A', 'MGAT4EP', 'RNASE6', 'AP000648.4', 'C1orf162'],
'Mixed' : ['HLA-DRB5', 'AGR2', 'MUC2', 'REG4', 'TFF3', 'CNNM2', 'HLA-DQA1', 'ZNF302', 'HLA-DRB1', 'TFF1', 'SPARC', 'HLA-DRA', 'MMP9', 'FN1', 'GPNMB', 'PIGR', 'CTSZ', 'FPR3', 'MRC1', 'CAPG', 'LYZ', 'TMEM176B', 'COL1A2', 'HLA-DPA1', 'COL1A1'],
'Mucinous' : ['MTRNR2L1', 'SLC11A1', 'CD163', 'CCL20', 'TIMP1', 'IGHV3-33', 'THBS1', 'RNASE1', 'ACSL1', 'REG1A', 'CLEC5A', 'C1orf162', 'CALB2', 'PLAUR', 'LUCAT1', 'TREM1', 'TGFBI', 'ALOX5AP', 'FPR1', 'SPRR2A', 'IGHV3-64D', 'EREG', 'IGHV1-69', 'SAT1', 'IL1R2'],
'NOS' : ['G0S2', 'OLFM4', 'MTCO1P12', 'NLRP3', 'TIMP1', 'BCL2A1', 'MTND4P12', 'PPIF', 'TNIP3', 'RPL10P9', 'ETS2', 'AC016739.1', 'ACOD1', 'RNF19B', 'BASP1', 'APOBEC3A', 'PLEK', 'SAT1', 'IGHG4', 'MTCO2P22', 'MMP19', 'S100A12', 'OASL', 'CCRL2', 'OLR1']
}

for key in histo_variants:
    subset.obs[key.lower()] = np.mean([subset.X[:, subset.var_names == gene].sum(axis=1) for gene in histo_variants[key]], axis=0)

sc.pl.umap(subset, color=[x.lower() for x in histo_variants.keys()],
           legend_loc="on data",
           frameon=False,
           ncols=5, cmap="YlGn"
          )


# In[14]:


sc.pl.rank_genes_groups_dotplot(
    subset,n_genes=10,values_to_plot="logfoldchanges",
    min_logfoldchange=3, vmax=7,vmin=-7,
    cmap="bwr", #groups= ['NOS','Mucinous']
)

sc.pl.dotplot(subset, marker_genes_dict, groupby="histological_variants",  cmap="bwr")


# In[15]:


top_genes = get_top_genes(subset, n_genes=25)
for i in top_genes.columns:
    gene_markers =top_genes[i].to_list()  # Replace with your actual gene list
    print(i,':', ', '.join(gene_markers))


# In[16]:


print(len(subset.var_names))
print(len(subset.obs))
print(subset.X.shape)
subset.obs = subset.obs.assign(numeric_index = lambda x: range(len(x)))
NOS_indices = subset.obs.loc[subset.obs['histological_variants'] == 'NOS'].numeric_index.to_list()
mixed_indices = subset.obs.loc[subset.obs['histological_variants'] == 'Mixed (NOS+Muc)'].numeric_index.to_list()
adenoma_like_indices = subset.obs.loc[subset.obs['histological_variants'] == 'Adenoma-like (mucinous 5%)'].numeric_index.to_list()
mucinous_indices = subset.obs.loc[subset.obs['histological_variants'] == 'Mucinous'].numeric_index.to_list()


# In[17]:


from scipy.stats import wilcoxon,mannwhitneyu,kruskal,f_oneway,ttest_ind
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt


def deg_function(mat1,mat2):
    pvalues = []
    logfolds = []
    for i in range(len(mat1)):
        statistic, pvalue = mannwhitneyu(mat1[i,:], mat2[i,:])
        logfold = np.log2(np.mean(mat1[i,:])/np.mean(mat2[i,:]))
        logfolds.append(logfold)
        pvalues.append(pvalue)

    reject, adj_pvalues, _, _ = multipletests(pvalues, alpha=0.05, method='fdr_bh')
    return list(zip(logfolds,adj_pvalues))

def plot_de(coords,title):
    xs = [x[0] for x in coords]
    ys = [ -np.log10(x[1]) for x in coords]
    plt.figure(figsize=(10, 5))
    plt.scatter(xs, ys, s=2)
    for i, (x,y) in enumerate(zip(xs,ys)):
        label = genes[i]
        if label == 'TREM2':
            plt.annotate(label, # Text to display
            (x, y), # Position of the text (x, y)
            textcoords="offset points", # How to position the text
            xytext=(0, 10), # Offset of the text from the point (x, y)
            ha='center',
            fontsize = 10,color = 'r') # Horizontal alignment

        if (x >1 or x < -1) and y > 10 :
            plt.annotate(label, # Text to display
            (x, y), # Position of the text (x, y)
            textcoords="offset points", # How to position the text
            xytext=(0, 10), # Offset of the text from the point (x, y)
            ha='center',
            fontsize = 6) # Horizontal alignment

    plt.axhline(y=-np.log10(0.01), color='grey', linestyle='--')
    plt.axvline(x=1, color='grey', linestyle='--', label='Line at x=1')
    plt.axvline(x=-1, color='grey', linestyle='--', label='Line at x=-1')
    plt.xlim(left=-5, right=5) # Set the range from 1 to 5
    plt.xlabel(title)
    plt.ylabel('Adjusted p-value')
    plt.show()

genes = subset.var_names.tolist()
adenoma_like_subset = subset.X[adenoma_like_indices].transpose()
nos_subset = subset.X[NOS_indices].transpose()
mixed_subset = subset.X[mixed_indices].transpose()
mucinous_subset = subset.X[mucinous_indices].transpose()

coords = deg_function(nos_subset,adenoma_like_subset)
plot_de(coords,"adenoma-like <- log2(FC) -> NOS")

coords = deg_function(nos_subset,mixed_subset)
plot_de(coords,"Mixed <- log2(FC) -> NOS")

coords = deg_function(nos_subset,mucinous_subset)
plot_de(coords,"Mucinous <- log2(FC) -> NOS")


# In[18]:


# Retrieve the results for the treated vs. control comparison
results = subset.uns["rank_genes_groups"]['params']
print(results)

