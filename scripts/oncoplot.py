import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from matplotlib.colors import ListedColormap
import matplotlib.patches as patches
import numpy.ma as ma
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import textwrap
from matplotlib.cm import viridis
from scipy.spatial.distance import squareform
import uuid

numeric_categories = lambda color_dict:  {k: int(i) for i, k in enumerate(color_dict.keys())}

def create_colormap(metadata, column, color_mapping):
    column_nan_filled = metadata[column].fillna("nan")
    numeric_mapping = column_nan_filled.map(numeric_categories(color_mapping))
    colormap = ListedColormap([color_mapping[k] for k in color_mapping.keys()])
    return numeric_mapping, colormap

# Matching dissimilarity function
def matching_dissim(a, b):
    return np.sum(a != b)

# Compute distance matrix and clustering
def compute_distance_and_clustering(df):
    num_samples = len(df)
    distance_matrix = np.zeros((num_samples, num_samples))
    for i in range(num_samples):
        for j in range(num_samples):
            distance_matrix[i, j] = matching_dissim(df.iloc[i].values, df.iloc[j].values)

    # Assuming distance_matrix is your full distance matrix
    condensed_distance_matrix = squareform(distance_matrix)
    linkage_matrix = linkage(condensed_distance_matrix, method="ward")
    #linkage_matrix = linkage(distance_matrix, method="median")
    dendro = dendrogram(linkage_matrix, labels=df.index, no_plot=True)
    return dendro['ivl']

def plot_annotation_legend(annotations):
    n_features = len(annotations)

    # Determine if the last feature needs to be split into two columns
    last_feature = list(annotations.keys())[-1]
    last_feature_categories = len(annotations[last_feature])

    # Create a new annotations dict for the modified layout
    modified_annotations = {}
    for i, (feature, mapping) in enumerate(annotations.items()):
        if i < n_features - 1:
            # Keep all features except the last one as is
            modified_annotations[feature] = mapping
        else:
            # For the last feature, split it into two parts if needed
            categories = list(mapping.keys())
            colors = list(mapping.values())

            # Calculate split point (approximately half)
            split_point = len(categories) // 2

            # Create two new features from the split
            modified_annotations[f"{feature} (1/2)"] = {cat: color for cat, color in zip(categories[:split_point], colors[:split_point])}
            modified_annotations[f"{feature} (2/2)"] = {cat: color for cat, color in zip(categories[split_point:], colors[split_point:])}

    # Adjust the plotting code for the modified layout
    n_modified_features = len(modified_annotations)
    fig, axes = plt.subplots(nrows=1, ncols=n_modified_features, figsize=(12 + 4, 5))
    fig.tight_layout(pad=1)

    for i, (feature, mapping) in enumerate(modified_annotations.items()):
        ax = axes[i]

        # Handle title differently for the split feature
        if "(1/2)" in feature or "(2/2)" in feature:
            base_feature = feature.split(" (")[0]
            part_indicator = feature.split(" ")[-1]

            # Only show the main title on the first part
            if "(1/2)" in feature:
                wrapped_feature = textwrap.fill(base_feature, width=30)
                ax.set_title(wrapped_feature, fontsize=12, fontweight='bold', loc='left')
            else:
                # For the second part, we can either leave it blank or show a continuation indicator
                ax.set_title("", fontsize=10, loc='left')
        else:
            wrapped_feature = textwrap.fill(feature, width=30)
            ax.set_title(wrapped_feature, fontsize=12, fontweight='bold', loc='left')

        categories = list(mapping.keys())
        colors = list(mapping.values())

        # Draw color patches and labels
        ax.set_aspect('equal')
        for j, (cat, color) in enumerate(zip(categories, colors)):
            # Adjust fontsize based on feature name (similar to your original code)
            base_feature_name = feature.split(" (")[0] if "(" in feature else feature
            fontsize = 9 if base_feature_name == 'organs' else 10

            ax.add_patch(
                plt.Rectangle(
                    (0, j), 1.5, 0.7,  # Width, Height
                    facecolor=color,
                    edgecolor='black',
                    linewidth=1
                )
            )
            ax.text(1.7, j + 0.35, str(cat), va='center', fontsize=fontsize)

        ax.set_xlim(0, 4)
        ax.set_ylim(0, len(categories))
        ax.invert_yaxis()
        ax.axis('off')

    plt.savefig("legend.svg")
    plt.show()
# Load and preprocess data
df = pd.read_excel('base oncoplots.xlsx')
df['unique_id'] = [str(uuid.uuid4()) for _ in range(len(df))]  # Generate UUIDs for each row
df.set_index('unique_id', inplace=True)

mutations = df.iloc[:,30:132]
metadata = pd.concat([df.iloc[:, :30], df.iloc[:, 132:]], axis=1)

mutations = mutations.fillna("WT").replace({
    0: "WT", 1: "Mutated", '1*': "Multiple Mutations", '1**': "Multiple Mutations", 3: "Deletion", 4: "Amplification"
    }).T

mutation_counts = (mutations != "WT").sum(axis=1)
mutation_frequency = (mutation_counts * 100) / mutations.shape[1]
number_of_genes = 25
top_genes = mutation_frequency.sort_values(ascending=False).head(number_of_genes).index
patients = compute_distance_and_clustering(mutations.T)
metadata = metadata.loc[patients]

mutation_colors = {
        "WT": "#cccfc4",
        "Mutated": "#ff00ff",
        "Multiple Mutations": "#ff99cc",
        "Deletion": "#008000",
        "Amplification": "#ff0000"
        }

ajcc_colors = {'II': '#268bd2', 'III': '#cb4b16', 'IV': '#6c71c4', 'I': '#859900', 'nan': '#586e75'}
organs_colors={'Colon': '#ff6f61', 'Rectum': '#f06358', 'Anal canal': '#dc5a52', 'Small intestine': '#c55651', 'Appendix': '#a95754', 'Gastric': '#85605f', 'Pancreas': '#666', 'Ovary': '#d33682', 'Testis': '#8d3bc3', 'Lung': '#3bf2e7', 'Skin': '#f4a261', 'Hamartoma': '#68a1b3', 'Teratoma': '#66a894', 'Urachal': '#649c6d', 'Head and neck': '#758e66', 'Breast': '#7e8066'}
msi_ihc_colors = { 'dMMR': '#4b8cba', 'pMMR': '#cb4b16', 'nan': '#839496'}
tmb_status_colors =  {'Low': '#268bd2', 'High': '#dc322f', 'Intermediate': '#859900', 'nan': '#657b83' }
msi_ngs_colors = { 'MSI_dMMR': '#1f77b4', 'MSS_pMMR': '#ac4d25', 'nan': '#586e75' }

# Apply the function to each column
ajcc_numeric, ajcc_cmap = create_colormap(metadata, "AJCC 9th ", ajcc_colors)
organ_numeric, organ_cmap = create_colormap(metadata, "organs", organs_colors)
msi_ihc_numeric, msi_ihc_cmap = create_colormap(metadata, "MSI status (IHC)", msi_ihc_colors)
tmb_status_numeric, tmb_status_cmap = create_colormap(metadata, "TMB status", tmb_status_colors)
msi_ngs_numeric, msi_ngs_cmap = create_colormap(metadata, "MSI status (NGS)", msi_ngs_colors)

annotation_bars = np.vstack([metadata['TMB'],ajcc_numeric, organ_numeric, msi_ihc_numeric, msi_ngs_numeric])

annotation_colors = {"Stage": ajcc_colors,"MSI (IHC)":msi_ihc_colors,"MSI (NGS)":msi_ngs_colors, "organs":organs_colors}
plot_annotation_legend(annotation_colors)

# Use the keys of annotations to determine annotation columns
number_of_annotations = 5

bar_height = 0.5
heatmap_width = 16
heatmap_height = 8
barplot_width = 1.5
legend_width = 1

fig_height = heatmap_width + number_of_annotations  * bar_height
fig_width = heatmap_width + barplot_width + legend_width
fig = plt.figure(figsize=(fig_width, fig_height))

gs = fig.add_gridspec(1 + number_of_annotations, 3,
                      width_ratios=[heatmap_width, barplot_width, legend_width],
                      height_ratios=[1.5] + [bar_height] * (number_of_annotations -1) + [heatmap_height])

ax0 = fig.add_subplot(gs[0, 0])
sns.barplot(annotation_bars[0], ax=ax0,width =1 ,linewidth=0,color="skyblue")
ax0.text(-1, 20, 'TMB', va='center', ha='right', fontsize=11)
ax0.spines['top'].set_visible(False)
ax0.spines['right'].set_visible(True)
ax0.spines['bottom'].set_visible(True)
ax0.spines['left'].set_visible(False)
ax0.yaxis.tick_right()
ax0.set_xticklabels([])
ax0.set_xlabel('')
ax0.set_ylim(bottom=-0.5)
ax0.axhline(y=20, color='black', linestyle='--', linewidth=1)
ax0.axhline(y=0, color='black', linewidth=1)

ax1 = fig.add_subplot(gs[1, 0])
sns.heatmap([annotation_bars[1]], cmap=ajcc_cmap, cbar=False, ax=ax1,xticklabels=False,yticklabels=False)
ax1.text(-1, 0.5, 'Stage', va='center', ha='right', fontsize=11)

ax2 = fig.add_subplot(gs[2, 0])
sns.heatmap([annotation_bars[2]], cmap=organ_cmap, cbar=False, ax=ax2,xticklabels=False,yticklabels=False)
ax2.text(-1, 0.5, 'Organs', va='center', ha='right', fontsize=11)


ax3 = fig.add_subplot(gs[3, 0])
sns.heatmap([annotation_bars[3]], cmap=msi_ihc_cmap, cbar=False, ax=ax3,xticklabels=False,yticklabels=False)
ax3.text(-1, 0.5, 'MSI (IHC)', va='center', ha='right', fontsize=11)

ax4 = fig.add_subplot(gs[4, 0])
sns.heatmap([annotation_bars[4]], cmap=msi_ngs_cmap, cbar=False, ax=ax4,xticklabels=False,yticklabels=False)
ax4.text(-1, 0.5, 'MSI (NGS)', va='center', ha='right', fontsize=11)
#################
##   Heatmap   ##
#################
value_to_int = {j: i for i, j in enumerate(mutation_colors.keys())}
cmap = ListedColormap([mutation_colors[cat] for cat in mutation_colors.keys()])
# Ensure the patients are in the correct order for the heatmap and annotations
mutations_top = mutations.loc[top_genes, patients]

ax5 = fig.add_subplot(gs[5, 0])
sns.heatmap(mutations_top.replace(value_to_int), cmap=cmap, linewidths=0.3, linecolor="white", cbar=False, ax=ax5,xticklabels= False)
ax5.set_ylabel("Genes")
ax5.set_xlabel("")

ax6 = fig.add_subplot(gs[5, 1])
ax6.barh(np.arange(len(top_genes)), mutation_frequency[top_genes], color="skyblue", height=0.9, align="edge")
ax6.set_xlabel("Frequency (%)", fontsize=12)
ax6.set_yticks([])
ax6.set_ylim(ax5.get_ylim())
ax6.set_xticks([20])
ax6.grid(axis="x", linestyle="--", alpha=0.7)
ax6.spines['top'].set_visible(False)
ax6.spines['right'].set_visible(False)
ax6.spines['bottom'].set_visible(False)
ax6.spines['left'].set_visible(False)
ax6.set_position([0.7, ax5.get_position().y0, 0.1, ax5.get_position().height])

ax7 = fig.add_subplot(gs[5, 2])
num_categories = len(mutation_colors)
spacing = 0.03
start_y = 0.9
rect_width = 0.1  # Reduce rectangle width
rect_height = spacing * 0.4
for i, (category, color) in enumerate(mutation_colors.items()):
    y_position = start_y - (i * spacing)
    rect = patches.Rectangle((0.1, y_position - rect_height/2), rect_width, rect_height, facecolor=color, edgecolor='black', linewidth=0.5)
    ax7.add_patch(rect)
    ax7.text(0.1 + rect_width + 0.02, y_position, category, va='center', fontsize=8)  # Reduce font size
ax7.set_xlim(0, 1)
ax7.set_ylim(0, 1)
ax7.axis('off')

fig.suptitle("Top Mutated Genes", fontsize=16)
plt.subplots_adjust(wspace=0, hspace=0)
plt.savefig("oncoplot.svg")
plt.show()
