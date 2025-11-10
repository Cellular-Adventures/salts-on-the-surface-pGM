import yaml
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.size'] = 10
fontfam = 'sans-serif'
plt.rcParams['font.family'] = fontfam
plt.rcParams[f'font.{fontfam}'] = ['helvetica', 'arial', *plt.rcParams[f'font.{fontfam}']]

plot_properties_yaml = "results/figures/plot_properties.yaml"
# Load plot properties
with open(plot_properties_yaml) as prop_yaml:
    plot_props = yaml.safe_load(prop_yaml)

positions = plot_props["positions"]
r_column = plot_props["d_column"] / 2

# Load experiment database
exp_h5 = pd.HDFStore("data/02_processed/per-experiment-data.h5")
exp_df = exp_h5['data']
exp_h5.close()

png_path = Path("results/figures/png")
svg_path = Path("results/figures/svg")
eps_path = Path("results/figures/eps")

# Left: BSD of position 5 with vertical line for median.
# Right: Median over positions
water_df = exp_df.loc[exp_df["Salt"] == "Water"]
# Manually found the first index where position == 5
bsd_pos_5 = water_df["BSD"].values[0][4] / 1000
median = np.median(bsd_pos_5)
# Select the first measurement at that position
fig, (hist_ax, r_ax) = plt.subplots(1, 2, figsize=(5.7, 2.5), layout="constrained")
# Plot a histogram of bubble sizes
hist_ax.hist(bsd_pos_5, bins=50, facecolor='gray', edgecolor='black', label="BSD")
hist_ax.plot([median, median], hist_ax.get_ylim(), '--', label="Median")
hist_ax.set_ylabel("Count [-]")
hist_ax.set_xlabel("Chord length [mm]")
hist_ax.set_title(f"BSD at center\n({len(bsd_pos_5)} bubbles)")
hist_ax.legend()

# Position 5 was done twice
# pos_indexes = [0:8]

bsd_all_pos = water_df["BSD"].values[0][0:9] / 1000
median_all_pos = [np.median(bsd) for bsd in bsd_all_pos]
pos = water_df["FP Positions"].values[0][0:9]

r_ax.plot(pos, median_all_pos, 'dk')
r_ax.set_ylabel("$d_{B,median}$ [mm]")
r_ax.set_xlabel("Normalized radius [-]")
r_ax.set_title("Radial bubble size profile")
r_ax.set_xlim((-1, 1))
r_ax.set_ylim((0, 4))
# fix axes: y 0-3000; x -1 - 1

fig.savefig(png_path / "Fig02_BSD_and_radial.png", dpi=900)
fig.savefig(svg_path / "Fig02_BSD_and_radial.svg")
fig.savefig(eps_path / "Fig02_BSD_and_radial.eps")
