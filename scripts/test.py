from colour import Color  # Install using `pip install colours`

# Define organ groups
organ_groups = {
    'Digestive': ['Colon', 'Rectum', 'Anal canal', 'Small intestine', 'Appendix', 'Gastric', 'Pancreas'],
    'Reproductive': ['Ovary', 'Testis'],
    'Respiratory': ['Lung'],
    'Skin': ['Skin'],
    'Misc': ['Hamartoma', 'Teratoma', 'Urachal', 'Head and neck', 'Breast']
}

# Assign base colors to each group
base_colors = {
    'Digestive': '#ff6f61',    # Warm red
    'Reproductive': '#d33682', # Pink
    'Respiratory': '#3bf2e7',  # Teal
    'Skin': '#f4a261',         # Orange
    'Misc': '#68a1b3'          # Neutral gray
}

# Generate a color gradient for each group
num_shades = 7  # Adjust based on group size
color_maps = {group: list(Color(base).range_to(Color("#666666"), num_shades)) for group, base in base_colors.items()}

# Assign colors to organs
organs_colors = {}
for group, organs in organ_groups.items():
    for i, organ in enumerate(organs):
        organs_colors[organ] = color_maps[group][i].hex

# Print the new color mapping
print(organs_colors)
import matplotlib.pyplot as plt
# Plotting
fig, ax = plt.subplots(figsize=(6, 8))
ax.set_xlim(0, 1)
ax.set_ylim(0, len(organs_colors))
ax.set_xticks([])
ax.set_yticks(range(len(organs_colors)))
ax.set_yticklabels(organs_colors.keys(), fontsize=12)

# Add colored rectangles
for i, (organ, color) in enumerate(organs_colors.items()):
    ax.add_patch(plt.Rectangle((0, i - 0.4), 1, 0.8, color=color))

# Remove axes
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)

plt.show()

