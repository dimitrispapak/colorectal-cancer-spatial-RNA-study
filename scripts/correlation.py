import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

def plot_correlation_boxplots(df, categorical_columns, numerical_columns, test='anova',
                            color_dicts=None, title="", file=None):
    """
    Plot boxplots only for significant correlations (p-value < 0.05) with custom colors.

    Parameters:
    - df: DataFrame containing the data
    - categorical_columns: List of categorical column names
    - numerical_columns: List of numerical column names
    - test: Statistical test to use ('anova', 'kruskal', 'mannwhitney', 'chi2', or 'spearman')
    - color_dicts: List of dictionaries (one per categorical column) mapping categories to colors
    - file: Optional file path to save the figure
    """
    # Set a style and context for the plots
    sns.set(style="whitegrid", context="notebook")

    # Validate color_dicts
    if color_dicts is not None:
        if len(color_dicts) != len(categorical_columns):
            raise ValueError("color_dicts must have the same length as categorical_columns")

    # First pass: calculate all p-values to determine which plots to show
    significant_plots = []
    for i, cat_col in enumerate(categorical_columns):
        row_plots = []
        for j, num_col in enumerate(numerical_columns):
            # Perform statistical test
            if test == 'anova':
                groups = [group for _, group in df.groupby(cat_col)[num_col]]
                _, p_value = stats.f_oneway(*groups)
            elif test == 'kruskal':
                groups = [group for _, group in df.groupby(cat_col)[num_col]]
                _, p_value = stats.kruskal(*groups)
            elif test == 'mannwhitney':
                unique_categories = df[cat_col].unique()
                if len(unique_categories) != 2:
                    continue  # Skip if not exactly two groups
                group1 = df[df[cat_col] == unique_categories[0]][num_col]
                group2 = df[df[cat_col] == unique_categories[1]][num_col]
                _, p_value = stats.mannwhitneyu(group1, group2)
            elif test == 'chi2':
                contingency_table = pd.crosstab(df[cat_col], df[num_col])
                _, p_value, _, _ = stats.chi2_contingency(contingency_table)
            elif test == 'spearman':
                _, p_value = stats.spearmanr(df[cat_col].astype('category').cat.codes, df[num_col])
            else:
                raise ValueError("Invalid test type. Choose from 'anova', 'kruskal', 'mannwhitney', 'chi2', or 'spearman'")

            if p_value < 0.05:
                row_plots.append((cat_col, num_col, p_value))

        if row_plots:  # Only add rows that have at least one significant plot
            significant_plots.append((i, row_plots))  # Store original index for color_dicts reference

    # If no significant plots, show message and return
    if not significant_plots:
        print("No significant correlations found (p-value < 0.05)")
        return

    # Calculate figure size (4 inches per column)
    max_cols = max(len(row_plots) for (_, row_plots) in significant_plots)
    fig_width = 4 * max_cols
    fig_height = 5 * len(significant_plots)

    fig, axes = plt.subplots(len(significant_plots), max_cols,
                             figsize=(fig_width, fig_height),
                             constrained_layout=True)

    # Ensure axes is a 2D array for consistent indexing
    if len(significant_plots) == 1:
        axes = [axes]
    if max_cols == 1:
        axes = [[ax] for ax in axes]

    # Set super title for the plot
    fig.suptitle(title, fontsize=16)

    # Plot only significant correlations
    for row_idx, (original_idx, row_plots) in enumerate(significant_plots):
        for col_idx, (cat_col, num_col, p_value) in enumerate(row_plots):
            ax = axes[row_idx][col_idx]

            # Get unique categories and their order
            categories = df[cat_col].unique()
            category_order = sorted(categories)

            # Create palette based on color_dicts if provided
            if color_dicts is not None:
                color_dict = color_dicts[original_idx]  # Use the original index
                palette = [color_dict.get(cat, "#999999") for cat in category_order]
            else:
                # Default palette
                palette = sns.color_palette("Set2", n_colors=len(category_order))

            # Create an enhanced boxplot with specified colors
            sns.boxplot(x=cat_col, y=num_col, data=df, ax=ax,
                       palette=palette, showfliers=False,
                       hue=cat_col, legend=False, order=category_order)

            # Overlay swarmplot with BLACK points
            sns.swarmplot(x=cat_col, y=num_col, data=df, ax=ax,
                         color="black", size=3, alpha=0.7,
                         order=category_order)

            # Get test name for display
            test_name = {
                'anova': 'ANOVA',
                'kruskal': 'Kruskal-Wallis',
                'mannwhitney': 'Mann-Whitney U',
                'chi2': 'Chi-square',
                'spearman': 'Spearman'
            }[test]

            # Add p-value annotation
            ax.text(0.05, 0.95, f'{test_name} p-value: {p_value:.4f}', transform=ax.transAxes,
                    verticalalignment='top', fontsize=10, bbox=dict(boxstyle="round", facecolor="white", alpha=0.5))

            # Set labels and title
            ax.set_xlabel(cat_col, fontsize=12)
            ax.set_ylabel(num_col, fontsize=12)
            ax.set_title(f'{num_col} vs {cat_col}', fontsize=14)

            # Rotate x-axis labels for readability if necessary
            if len(category_order) > 3:
                plt.setp(ax.get_xticklabels(), rotation=45, ha='right')

        # Hide any unused subplots in this row
        for col_idx in range(len(row_plots), max_cols):
            axes[row_idx][col_idx].axis('off')

    # Save the figure if a file path is provided
    if file:
        plt.savefig(file, bbox_inches='tight', dpi=300)
    plt.show()

def plot_numerical_correlations(df, correlations,numeric_columns,color_dict, title = "", file=None):
    # Set a style and context for the plots
    sns.set(style="whitegrid", context="notebook")

    num_plots = len(correlations) * 2
    nrows = 2
    ncols = (num_plots + 1) // 2
    fig, axes = plt.subplots(nrows, ncols, figsize=(len(correlations) * 4 , 8), constrained_layout=True)

    # Set super title for the plot
    fig.suptitle(title, fontsize=16)


    for idx, cell_type in enumerate(correlations):
        print("Correlation of {}".format(cell_type))
        ax0 = axes[0,idx]
        ax1 = axes[1,idx]

        # Calculate correlation
        correlation1 = df[cell_type].corr(df[numeric_columns[0]])

        # Create an enhanced scatter plot with regression line
        sns.regplot(x=cell_type, y=numeric_columns[0], data=df, ax=ax0,
                    scatter_kws={'alpha': 0.5, 'color': color_dict[cell_type]},
                    line_kws={'color': 'red'})

        # Add correlation coefficient as annotation
        ax0.legend([f'Correlation: {correlation1:.4f}'], loc='best', fontsize=10)
        ax0.set_xlabel("")
        ax0.set_ylabel( f'{numeric_columns[0]} - proximity to tumor', fontsize=10)
        ax0.set_title(cell_type, fontsize=8)

        correlation2 = df[cell_type].corr(df[numeric_columns[1]])
        # Create an enhanced scatter plot with regression line
        sns.regplot(x=cell_type, y=numeric_columns[1], data=df, ax=ax1,
                    scatter_kws={'alpha': 0.5, 'color': color_dict[cell_type]},
                    line_kws={'color': 'red'})

        # Add correlation coefficient as annotation
        ax1.legend([f'Correlation: {correlation2:.4f}'], loc='best', fontsize=10)
        ax1.set_xlabel("")
        ax1.set_ylabel('{numeric_columns[1]}- proximity to leukocytes', fontsize=10)


    # Turn off unused axes
    for idx in range(len(correlations), len(axes)):
        axes[0,idx].axis('off')
        axes[1,idx].axis('off')

    if file:
        plt.savefig(file, bbox_inches='tight', dpi=300)
    plt.show()
