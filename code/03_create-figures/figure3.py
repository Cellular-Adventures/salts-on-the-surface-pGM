import yaml
import pandas as pd
import matplotlib.pyplot as plt
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
exp_h5 = pd.HDFStore("data/02_processed/per-experiment-data.h5")
exp_df = exp_h5['data']
exp_h5.close()

png_path = Path("06_figures/png")
svg_path = Path("06_figures/svg")
eps_path = Path("06_figures/eps")

# For all salts - 3x3 figure, created plot-by-plot
salts = exp_df["Salt"].unique()
row_indexes = ["dmedian", "Holdup FP (Augmented mean)", "Interfacial area (m2/m3) - Q1-Q99"]
row_labels = ["$d_{B,median}$ [$mm$]", "Gas fraction [-]", "Interfacial area\n[$m^2/m^3$]"]
column_indexes = ["Concentration (mol/l)", "Ionic strength", "pGM"]
column_labels = ["Concentration [mol/l]", "Ionic strength [mol/l]", "$p_{GM}$ [Pa]"]
column_titles = ["Concentration", "Ionic strength", "GM pressure"]

fig, axs = plt.subplots(3, 3, figsize=(16/2.5, 11/2.5), layout="constrained", 
                        sharex='col', sharey='row')
rows = axs.shape[0]
columns = axs.shape[1]
for i in range(columns):
    c_index = column_indexes[i]
    c_label = column_labels[i]
    c_title = column_titles[i]

    for j in range(rows):
        r_index = row_indexes[j]
        r_label = row_labels[j]
        ax = axs[j, i]
        ax.grid(which='both', color="#eeeeee")
        if i == 0:
            ax.set_ylabel(r_label)
        if j == rows - 1:
            ax.set_xlabel(c_label)
            ax.set_xscale('log')
        if j == 0:
            ax.set_title(c_title)

        for salt in salts:
            salt_boi = exp_df.loc[exp_df["Salt"] == salt]
            if j == 0:
                # d32 values are recorded at 5 positions for every concentration
                for c in range(len(salt_boi[c_index])):
                    if len(salt_boi[r_index].iloc[c][0]) != 5:
                        # This one only has 4 points. Doesn't work nicely
                        # Found for NaCl, 0.037 M
                        # TODO Fix data for these ones. Data is available
                        continue
                    ax.plot(
                        salt_boi[c_index].iloc[c],
                        salt_boi[r_index].iloc[c][0][3],
                        color=salt_colors[salt],
                        marker=markers[salt],
                        markeredgewidth=0.4,
                        markeredgecolor='k',
                        linestyle='none',
                        label=salt_labels[salt]
                    )
            else:
                ax.plot(
                    salt_boi[c_index],
                    salt_boi[r_index],
                    color=salt_colors[salt],
                    marker=markers[salt],
                    markeredgewidth=0.4,
                    markeredgecolor='k',
                    linestyle='none',
                    label=salt_labels[salt]
                )

fig.legend(
    ax.get_legend_handles_labels()[0][1:],
    ax.get_legend_handles_labels()[1][1:],
    loc="outside lower center",
    ncols=len(ax.get_legend_handles_labels()[0][1:])
)

# save figure
fig.savefig(png_path / "Fig03_Salts_over_c_proxies.png", dpi=900)
fig.savefig(svg_path / "Fig03_Salts_over_c_proxies.svg")
fig.savefig(eps_path / "Fig03_Salts_over_c_proxies.eps")
