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

salts = exp_df["Salt"].unique()

png_path = Path("06_figures/png")
svg_path = Path("06_figures/svg")
eps_path = Path("06_figures/eps")

water_df = exp_df.loc[exp_df["Salt"] == "Water"]
# plot a comparison of interfacial area with some literature
water_ia = water_df["Interfacial area (m2/m3) - Q1-Q99"][0]
# Xue 2008 - https://doi.org/10.1002/AIC.11386
# looking at data at 8 cm/s
Xue_radial_ia = np.array([
    [-0.0024630541871920597, 289.03903903903904],
    [0.29556650246305427, 283.033033033033],
    [-0.3054187192118226, 285.2852852852853],
    [-0.5985221674876847, 232.73273273273273],
    [0.6009852216748768, 223.7237237237237],
    [0.8990147783251232, 243.24324324324323],
    [-0.9039408866995073, 189.93993993993993],
    [1, 0]
])

Xue_radial_ia[:, 0] = np.abs(Xue_radial_ia[:, 0])
Xue_radial_pos = Xue_radial_ia[:, 0]
Xue_radial_ia = Xue_radial_ia[:, 1]
Xue_2008_ia_calc = 2 * np.trapezoid(Xue_radial_ia * Xue_radial_pos, Xue_radial_pos)
print(Xue_2008_ia_calc)

Xue_2008_ia = Xue_2008_ia_calc
Xue_2008_ia_pm = np.nan
Xue_2008_ugs = 8

# Schumpe 1982 - https://doi.org/10.1021/I200019A028
# Looking at data at 5 cm/s. A dotted line on a graph with solid lines.
Schumpe_1982_ia = 150
Schumpe_1982_ia_pm = np.nan
Schumpe_1982_ugs = 5

# Besagni 2017 - https://doi.org/10.1016/J.CES.2017.03.043
# Data at 2 cm/s. Data shows some weird behaviour at this (or the previous) 
# datapoint.
Besagni_2017_ia = 90       # Converted to m2/m3
Besagni_2017_ia_pm = np.nan
Besagni_2017_ugs = 2

tick_labels = ["This work", "Xue, 2008", "Schumpe, 1982", "Besagni, 2017"]
fig, ax = plt.subplots(1, 1, figsize=(4, 4), layout="constrained")
ax.bar(
    tick_labels,
    [water_ia, Xue_2008_ia, Schumpe_1982_ia, Besagni_2017_ia],
    yerr=[np.nan, Xue_2008_ia_pm, Schumpe_1982_ia_pm, Besagni_2017_ia_pm]
)
ax.set_ylabel("Interfacial area ($m^2/m^3$)")
ax.set_xticks(range(4), tick_labels, rotation=45, ha='right')

fig.savefig(png_path / "FigB1_Literature_comparison.png", dpi=900)
fig.savefig(svg_path / "FigB1_Literature_comparison.svg")
fig.savefig(eps_path / "FigB1_Literature_comparison.eps")
