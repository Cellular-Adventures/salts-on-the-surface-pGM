import yaml
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from pathlib import Path

plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.size'] = 10
fontfam = 'sans-serif'
plt.rcParams['font.family'] = fontfam
plt.rcParams[f'font.{fontfam}'] = ['helvetica', 'arial', *plt.rcParams[f'font.{fontfam}']]

plot_properties_yaml = "06_figures/plot_properties.yaml"
# Load plot properties
with open(plot_properties_yaml) as prop_yaml:
    plot_props = yaml.safe_load(prop_yaml)

salt_colors = plot_props["salt_colors"]
salt_labels = plot_props["salt_labels"]
markers = plot_props["markers"]
positions = plot_props["positions"]
r_column = plot_props["d_column"] / 2

# Load experiment database
exp_h5 = pd.HDFStore("01_data/02_processed/per-experiment-data.h5")
exp_df = exp_h5['data']
exp_h5.close()

salts = exp_df["Salt"].unique()

png_path = Path("06_figures/png")
svg_path = Path("06_figures/svg")
eps_path = Path("06_figures/eps")

# Plot comparison between pp, fp (Schweitzer) and fp (augmented) as parity
fig, ax = plt.subplots(1, 1, figsize=(4.8, 3), layout="constrained")
# fig.set_constrained_layout_pads(wspace=0.01)
hu_min = min(exp_df["Holdup PP"].min(),
             exp_df["Holdup FP (Schweitzer mean)"].min(),
             exp_df["Holdup FP (Augmented mean)"].min()) * 0.9
hu_max = max(exp_df["Holdup PP"].max(),
             exp_df["Holdup FP (Schweitzer mean)"].max(),
             exp_df["Holdup FP (Augmented mean)"].max()) * 1.1
ax.plot([hu_min, hu_max],
        [hu_min, hu_max],
        '-k')
ax.plot([hu_min, hu_max], [hu_min * 0.8, hu_max * 0.8], '--k')
ax.plot([hu_min, hu_max], [hu_min * 1.2, hu_max * 1.2], '--k')
for salt in salts[::-1]:
    salt_boi = exp_df.loc[exp_df["Salt"] == salt]
    line_schw = ax.plot(
        salt_boi["Holdup PP"],
        salt_boi["Holdup FP (Schweitzer mean)"],
        's',
        color=salt_colors[salt],
        markeredgewidth=0.5,
        markeredgecolor='black',
        label="FP Only"
    )
ax.set_aspect('equal')
ax.set_title("Parity plot of global holdup")
ax.set_xlabel("Pressure probe holdup [-]")
ax.set_ylabel("Fiber probe holdup [-]")
ax.set_ylim(hu_min, hu_max)
ax.set_xlim(hu_min, hu_max)
for salt in salts[::-1]:
    salt_boi = exp_df.loc[exp_df["Salt"] == salt]
    line_aug = ax.plot(
        salt_boi["Holdup PP"],
        salt_boi["Holdup FP (Augmented mean)"],
        'o',
        color=salt_colors[salt],
        markeredgewidth=0.5,
        markeredgecolor='black',
        label="Combined"
    )
schw_aug_legend = ax.legend(handles=[line_schw[0], line_aug[0]])
ax.add_artist(schw_aug_legend)
fig_legend_elements = []
for salt in salts:
    fig_legend_elements.append(Patch(
        facecolor=salt_colors[salt],
        edgecolor='k',
        linewidth=0.5,
        label=salt_labels[salt]
    ))
fig.legend(
    handles=fig_legend_elements,
    loc='outside right center',
    # bbox_to_anchor=(1.02, 0.5),
    # borderaxespad=1
)
fig.savefig(png_path / "Fig05_holdup_parity.png", dpi=900)
fig.savefig(svg_path / "Fig05_holdup_parity.svg")
fig.savefig(eps_path / "Fig05_holdup_parity.eps")
