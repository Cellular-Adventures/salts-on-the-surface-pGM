import yaml
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
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

bubble_h5 = pd.HDFStore("data/02_processed/consolidated_measurement_data.h5")
bubbles_df = bubble_h5['data']
bubble_h5.close()

flow = 100

water_data = bubbles_df.loc[
    (bubbles_df["Concentration (mol/l)"] == 0.0)
    & (bubbles_df["Ignore"] == False)   # Linters complain about `== False`. It's needed here.
    & (bubbles_df["Gas flow (l/min)"] == flow)
    & (bubbles_df["FP position"] == 5)
]
# df extracted is still going to contain multiple measurements.
# Need to combine if that's the case.
water_bubbles = [[], []]
for i, row in water_data.iterrows():
    for size, velo in zip(row['Valid bubbles'], row['FP velocity']):
        water_bubbles[0].append(size)
        water_bubbles[1].append(velo)


Na2SO4_data = bubbles_df.loc[
    (bubbles_df["Salt"] == "Na2SO4")
    & (bubbles_df["Concentration (mol/l)"] == 0.5)
    & (bubbles_df["Ignore"] == False)   # Linters complain about `== False`. It's needed here.
    & (bubbles_df["Gas flow (l/min)"] == flow)
    & (bubbles_df["FP position"] == 5)
]
Na2SO4_bubbles = [[], []]
for i, row in Na2SO4_data.iterrows():
    for size, velo in zip(row['Valid bubbles'], row['FP velocity']):
        Na2SO4_bubbles[0].append(size)
        Na2SO4_bubbles[1].append(velo)

fig, ax = plt.subplots(1, 1, figsize=(16/2.5, 11/2.5), layout="constrained")
ax.set_xlabel("Bubble size ($\\mu m$)")
ax.set_ylabel("Bubble velocity ($m/s$)")
ax.scatter(water_bubbles[0], water_bubbles[1], alpha=0.3, label="Tap water")
ax.scatter(Na2SO4_bubbles[0], Na2SO4_bubbles[1], alpha=0.3, label=f"0.5 M {plot_props['salt_labels']['Na2SO4']}")


def scatter_hist(x, y, ax, ax_histx, ax_histy, scatter_label, color, align):
    # no labels
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    # the scatter plot:
    ax.scatter(x, y, label=scatter_label, alpha=0.5, c=color, edgecolors='grey')

    # now determine nice limits by hand:
    binwidth_x = 1000
    xmax = np.max(np.abs(x))
    xlim = (int(xmax / binwidth_x) + 1) * binwidth_x
    xbins = np.arange(-xlim, xlim + binwidth_x, binwidth_x)

    binwidth_y = 0.1
    ymax = np.max(np.abs(y))
    ylim = (int(ymax / binwidth_y) + 1) * binwidth_y
    ybins = np.arange(-ylim, ylim + binwidth_y, binwidth_y)

    ax_histx.hist(x, bins=xbins, density=True, rwidth=0.5, align=align, color=color)
    ax_histy.hist(y, bins=ybins, density=True, rwidth=0.5, align=align, orientation='horizontal', color=color)


fig, axs = plt.subplot_mosaic([['histx', '.'],
                               ['scatter', 'histy']],
                              figsize=(16/2.5, 11/2.5),
                              width_ratios=(4, 1), height_ratios=(1, 4),
                              layout='constrained')
scatter_hist(
    np.array(water_bubbles[0]),
    np.array(water_bubbles[1]),
    axs['scatter'],
    axs['histx'],
    axs['histy'],
    scatter_label="Tap water",
    color=plot_props['salt_colors']['Water'],
    align='right'
)
scatter_hist(
    np.array(Na2SO4_bubbles[0]),
    np.array(Na2SO4_bubbles[1]),
    axs['scatter'],
    axs['histx'],
    axs['histy'],
    scatter_label="0.5 M Na2SO4",
    color=plot_props['salt_colors']['Na2SO4'],
    align='mid')

axs['histx'].set_xlim(axs['scatter'].get_xlim())
axs['histy'].set_ylim(axs['scatter'].get_ylim())
axs['scatter'].set_ylabel("Bubble velocity ($m/s$)")
axs['scatter'].set_xlabel("Bubble size ($\\mu m$)")

res_water = stats.linregress(water_bubbles[0], water_bubbles[1])
res_salt = stats.linregress(Na2SO4_bubbles[0], Na2SO4_bubbles[1])

xx = np.linspace(min(water_bubbles[0]), max(water_bubbles[0]), 1000)

axs['scatter'].plot(xx, res_water.intercept + res_water.slope * xx,
                    c=plot_props['salt_colors']['Water'],
                    linewidth=3)
axs['scatter'].plot(xx, res_water.intercept - res_water.intercept_stderr + (res_water.slope - res_water.stderr) * xx,
                    c=plot_props['salt_colors']['Water'],
                    linestyle='--',
                    linewidth=3)
axs['scatter'].plot(xx, res_water.intercept + res_water.intercept_stderr + (res_water.slope + res_water.stderr) * xx,
                    c=plot_props['salt_colors']['Water'],
                    linestyle='--',
                    linewidth=3)

axs['scatter'].plot(xx, res_salt.intercept + res_salt.slope * xx,
                    c=plot_props['salt_colors']['Na2SO4'],
                    linewidth=3)
axs['scatter'].plot(xx, res_salt.intercept - res_salt.intercept_stderr + (res_salt.slope - res_salt.stderr) * xx,
                    c=plot_props['salt_colors']['Na2SO4'],
                    linestyle=':',
                    linewidth=3)
axs['scatter'].plot(xx, res_salt.intercept + res_salt.intercept_stderr + (res_salt.slope + res_salt.stderr) * xx,
                    c=plot_props['salt_colors']['Na2SO4'],
                    linestyle=':',
                    linewidth=3)
axs['scatter'].legend()

png_path = Path("results/figures/png")
svg_path = Path("results/figures/svg")
eps_path = Path("results/figures/eps")

fig.savefig(png_path / "FigC7_size_velocity.png", dpi=900)
fig.savefig(svg_path / "FigC7_size_velocity.svg")
fig.savefig(eps_path / "FigC7_size_velocity.eps")

plt.show()


