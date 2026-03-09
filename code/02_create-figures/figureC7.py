import yaml
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path

# Stats stuff to turn standard error into a 95% confidence interval.
from scipy.stats import t

tinv = lambda p, df: abs(t.ppf(p/2, df))

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

bubbles = {}
# df extracted is still going to contain multiple measurements.
# Need to combine if that's the case.
bubbles['Water'] = [[], []]
for i, row in water_data.iterrows():
    for size, velo in zip(row['Valid bubbles'], row['FP velocity']):
        bubbles['Water'][0].append(size / 1_000_000)
        bubbles['Water'][1].append(velo)


Na2SO4_data = bubbles_df.loc[
    (bubbles_df["Salt"] == "Na2SO4")
    & (bubbles_df["Concentration (mol/l)"] == 0.5)
    & (bubbles_df["Ignore"] == False)   # Linters complain about `== False`. It's needed here.
    & (bubbles_df["Gas flow (l/min)"] == flow)
    & (bubbles_df["FP position"] == 5)
]
bubbles['Na2SO4'] = [[], []]
for i, row in Na2SO4_data.iterrows():
    for size, velo in zip(row['Valid bubbles'], row['FP velocity']):
        bubbles['Na2SO4'][0].append(size / 1_000_000)
        bubbles['Na2SO4'][1].append(velo)

NaCl_data = bubbles_df.loc[
    (bubbles_df["Salt"] == "NaCl")
    & (bubbles_df["Concentration (mol/l)"] == 0.5)
    & (bubbles_df["Ignore"] == False)   # Linters complain about `== False`. It's needed here.
    & (bubbles_df["Gas flow (l/min)"] == flow)
    & (bubbles_df["FP position"] == 5)
]
bubbles['NaCl'] = [[], []]
for i, row in NaCl_data.iterrows():
    for size, velo in zip(row['Valid bubbles'], row['FP velocity']):
        bubbles['NaCl'][0].append(size / 1_000_000)
        bubbles['NaCl'][1].append(velo)

NH4Cl_data = bubbles_df.loc[
    (bubbles_df["Salt"] == "NH4Cl")
    & (bubbles_df["Concentration (mol/l)"] == 0.5)
    & (bubbles_df["Ignore"] == False)   # Linters complain about `== False`. It's needed here.
    & (bubbles_df["Gas flow (l/min)"] == flow)
    & (bubbles_df["FP position"] == 5)
]
bubbles['NH4Cl'] = [[], []]
for i, row in NH4Cl_data.iterrows():
    for size, velo in zip(row['Valid bubbles'], row['FP velocity']):
        bubbles['NH4Cl'][0].append(size / 1_000_000)
        bubbles['NH4Cl'][1].append(velo)

NH42SO4_data = bubbles_df.loc[
    (bubbles_df["Salt"] == "NH42SO4")
    & (bubbles_df["Concentration (mol/l)"] == 0.5)
    & (bubbles_df["Ignore"] == False)   # Linters complain about `== False`. It's needed here.
    & (bubbles_df["Gas flow (l/min)"] == flow)
    & (bubbles_df["FP position"] == 5)
]
bubbles['NH42SO4'] = [[], []]
for i, row in NH42SO4_data.iterrows():
    for size, velo in zip(row['Valid bubbles'], row['FP velocity']):
        bubbles['NH42SO4'][0].append(size / 1_000_000)
        bubbles['NH42SO4'][1].append(velo)

NaAc_data = bubbles_df.loc[
    (bubbles_df["Salt"] == "NaAc")
    & (bubbles_df["Concentration (mol/l)"] == 0.92)
    & (bubbles_df["Ignore"] == False)   # Linters complain about `== False`. It's needed here.
    & (bubbles_df["Gas flow (l/min)"] == flow)
    & (bubbles_df["FP position"] == 5)
]
bubbles['NaAc'] = [[], []]
for i, row in NaAc_data.iterrows():
    for size, velo in zip(row['Valid bubbles'], row['FP velocity']):
        bubbles['NaAc'][0].append(size / 1_000_000)
        bubbles['NaAc'][1].append(velo)

fig, ax = plt.subplots(1, 1, figsize=(16/2.5, 11/2.5), layout="constrained")
ax.set_xlabel("Bubble size (m)")
ax.set_ylabel("Bubble velocity (m/s)")
ax.scatter(bubbles['Water'][0], bubbles['Water'][1], alpha=0.3, label="Tap water")
ax.scatter(bubbles['Na2SO4'][0], bubbles['Na2SO4'][1], alpha=0.3, label=f"0.5 M {plot_props['salt_labels']['Na2SO4']}")


def scatter_hist(x, y, ax, ax_histx, ax_histy, scatter_label, color, align):
    # no labels
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    # the scatter plot:
    ax.scatter(x, y, label=scatter_label, alpha=0.5, c=color, edgecolors='grey')

    # now determine nice limits by hand:
    binwidth_x = 1e-3
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
    np.array(bubbles['Water'][0]),
    np.array(bubbles['Water'][1]),
    axs['scatter'],
    axs['histx'],
    axs['histy'],
    scatter_label=f"Tap water ({len(bubbles['Water'][0])} bubbles)",
    color=plot_props['salt_colors']['Water'],
    align='right'
)
scatter_hist(
    np.array(bubbles['Na2SO4'][0]),
    np.array(bubbles['Na2SO4'][1]),
    axs['scatter'],
    axs['histx'],
    axs['histy'],
    scatter_label=f"0.5 M Na2SO4 ({len(bubbles['Na2SO4'][0])} bubbles)",
    color=plot_props['salt_colors']['Na2SO4'],
    align='mid')

axs['histx'].set_xlim(axs['scatter'].get_xlim())
axs['histy'].set_ylim(axs['scatter'].get_ylim())

axs['histx'].set_ylabel('PD $(1/m)$')
axs['histy'].set_xlabel('PD $(s/m)$')

axs['scatter'].set_ylabel("Bubble velocity $(m/s)$")
axs['scatter'].set_xlabel("Bubble size $(m)$")

res_water = stats.linregress(bubbles['Water'][0], bubbles['Water'][1])
res_salt = stats.linregress(bubbles['Na2SO4'][0], bubbles['Na2SO4'][1])

xx = np.linspace(min(bubbles['Water'][0]), max(bubbles['Water'][0]), 1000)

ts = tinv(0.05, len(bubbles['Water'][0]) - 2)
axs['scatter'].plot(xx, res_water.intercept + res_water.slope * xx,
                    c=plot_props['salt_colors']['Water'],
                    linewidth=3)
axs['scatter'].plot(xx, res_water.intercept - ts*res_water.intercept_stderr + (res_water.slope - ts*res_water.stderr) * xx,
                    c=plot_props['salt_colors']['Water'],
                    linestyle='--',
                    linewidth=3)
axs['scatter'].plot(xx, res_water.intercept + ts*res_water.intercept_stderr + (res_water.slope + ts*res_water.stderr) * xx,
                    c=plot_props['salt_colors']['Water'],
                    linestyle='--',
                    linewidth=3)

ts = tinv(0.05, len(bubbles['Na2SO4'][0]) - 2)
axs['scatter'].plot(xx, res_salt.intercept + res_salt.slope * xx,
                    c=plot_props['salt_colors']['Na2SO4'],
                    linewidth=3)
axs['scatter'].plot(xx, res_salt.intercept - ts*res_salt.intercept_stderr + (res_salt.slope - ts*res_salt.stderr) * xx,
                    c=plot_props['salt_colors']['Na2SO4'],
                    linestyle='--',
                    linewidth=3)
axs['scatter'].plot(xx, res_salt.intercept + ts*res_salt.intercept_stderr + (res_salt.slope + ts*res_salt.stderr) * xx,
                    c=plot_props['salt_colors']['Na2SO4'],
                    linestyle='--',
                    linewidth=3)
axs['scatter'].legend()

png_path = Path("results/figures/png")
svg_path = Path("results/figures/svg")
eps_path = Path("results/figures/eps")

fig.savefig(png_path / "FigC7_size_velocity.png", dpi=900)
fig.savefig(svg_path / "FigC7_size_velocity.svg")
fig.savefig(eps_path / "FigC7_size_velocity.eps")

res = {}
for salt in bubbles.keys():
    res[salt] = stats.linregress(bubbles[salt][0], bubbles[salt][1])
    ts = tinv(0.05, len(bubbles[salt][0]) - 2)
    print(f"{salt}: {res[salt].slope:.3f} +/- {ts*res[salt].stderr:.3f} * x + "
          f"{res[salt].intercept:.3f} +/- {ts*res[salt].intercept_stderr:.3f}")

plt.show()
