import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt

root_data_folder = Path('data/02_processed/01_Literature-Data')

d = {}
d_keys = ["KCl", "NaCl", "MgSO4", "NaClO4", "HCl", "HClO4", "NaCl+HClO4",
          "NaNO3+HCl", "NaAc+KClO3", "NaClO3+HClO4", "HCl+HClO4", "HCl+HNO3"]
for d_key in d_keys:
    d[d_key] = pd.read_csv(root_data_folder / f"LitData_{d_key}.csv").to_dict(orient="list")

legend_keys1 = [
    "$\\mathrm{KCl}$",
    "$\\mathrm{NaCl}$",
    "$\\mathrm{MgSO_4}$",
    "$\\mathrm{NaClO_4}$",
    "$\\mathrm{HCl}$",
    "$\\mathrm{HClO_4}$",
    
    "$\\mathrm{NaCl+HClO_4}$",
    "$\\mathrm{NaNO_3+HCl}$",
    "$\\mathrm{NaAc+KClO_3}$",
    "$\\mathrm{NaClO_3+HClO_4}$",
    "$\\mathrm{HCl+HClO_4}$",
    "$\\mathrm{HCl+HNO_3}$"
]

legend_keys2 = [
    "$\\mathrm{(\\alpha\\alpha)}$",
    "$\\mathrm{(\\alpha\\alpha)}$",
    "$\\mathrm{(\\alpha\\alpha)}$",
    "$\\mathrm{(\\alpha\\beta)}$",
    "$\\mathrm{(\\beta\\alpha)}$",
    "$\\mathrm{(\\beta\\beta)}$",

    "$\\mathrm{(\\alpha\\alpha + \\beta\\beta)}$",
    "$\\mathrm{(\\alpha\\alpha + \\beta\\alpha)}$",
    "$\\mathrm{(\\alpha\\beta + \\alpha\\beta)}$",
    "$\\mathrm{(\\alpha\\beta + \\beta\\beta)}$",
    "$\\mathrm{(\\beta\\alpha + \\beta\\beta)}$",
    "$\\mathrm{(\\beta\\alpha + \\beta\\alpha)}$"
]

markers = [
    "o",
    "o",
    "o",
    "^",
    "^",
    "o",
    "o",
    "s",
    "^",
    "s",
    "s",
    "^",
]

legend_keys = legend_keys1 + legend_keys2

plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.size'] = 7
fontfam = 'sans-serif'
plt.rcParams['font.family'] = fontfam
plt.rcParams[f'font.{fontfam}'] = ['helvetica', 'arial', *plt.rcParams[f'font.{fontfam}']]

color_keys = [
    '#a6cee3',
    '#1f78b4',
    '#b2df8a',
    '#33a02c',
    '#fb9a99',
    '#e31a1c',
    '#fdbf6f',
    '#ff7f00',
    '#cab2d6',
    '#6a3d9a',
    '#ffff99',
    '#b15928']

order_keys = np.array(12*[10])
plot_keys = ["KCl", "NaCl", "MgSO4", "NaClO4", "HCl", "HClO4",
             "NaCl+HClO4", "NaNO3+HCl", "NaAc+KClO3", "NaClO3+HClO4", "HCl+HClO4", "HCl+HNO3"]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=[12/2.54, 2.5], layout="constrained")
fig.suptitle('Salt concentration proxies for coalescence')

ax1.set_title(r'Ionic strength')
ax1.set_xlabel('I [M]')
ax1.grid(which='both', color="#eeeeee")
ax1.text(-.05, 1.05, "a", fontsize=9, fontweight='bold', transform=ax1.transAxes)

for i, key in enumerate(plot_keys):
    # print(d[key]["ionS"])
    ax1.scatter(
        d[key]["ionS"],
        d[key]["D"],
        zorder=order_keys[i],
        label=legend_keys[i],
        marker=markers[i],
        color=color_keys[i],
        s=30,
        edgecolors='black',
        lw=0.5)

ax2.set_title(r'Gibbs-Marangoni pressure')
ax2.grid(which='both', color="#eeeeee")
ax2.text(-.05, 1.05, "b", fontsize=9, fontweight='bold', transform=ax2.transAxes)

for i, key in enumerate(plot_keys):
    ax2.scatter(
        d[key]["pGM"],
        d[key]["D"],
        zorder=order_keys[i],
        marker=markers[i],
        color=color_keys[i],
        s=30,
        edgecolors='black',
        lw=0.5)
for i, key in enumerate(plot_keys):
    ax2.plot(np.nan, np.nan, '.', ms=0)
ax2.set_xlabel('$\\mathrm{p_{GM}\\ [Pa]}$')
ax2.legend(
    legend_keys,
    fontsize='x-small',
    ncols=2,
    # markerscale=0.5,
    columnspacing=-2,
    labelspacing=0,
    # prop={'family': 'DejaVu Sans', 'size': 'x-small'},
    loc='lower left')

for ax in fig.get_axes():
    ax.set_xscale('log')
#     ax.set_xlim(0.015,1.5)
#     ax.set_ylim(0.015,1.5)
    # ax.tick_params(labelsize=12)
    ax.set_ylabel('$\\mathrm{Coalescence\\ [\\%]}$')
    ax.label_outer()
plt.savefig('results/figures/png/Fig01_LiteratureCoalescence.png', bbox_inches='tight', dpi=900)
plt.savefig('results/figures/svg/Fig01_LiteratureCoalescence.svg')
plt.savefig('results/figures/eps/Fig01_LiteratureCoalescence.eps')

plt.show()
