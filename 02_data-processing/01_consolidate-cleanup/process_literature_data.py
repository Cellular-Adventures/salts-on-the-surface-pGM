import pandas as pd
from pathlib import Path

import warnings
warnings.filterwarnings('ignore')

import juliacall
from juliacall import Main as jl

input_data_folder = Path("01_data/01_raw/02_unpacked/LitData")
output_data_folder = Path("01_data/02_processed/01_Literature-Data")

df_NaCl = pd.read_csv(input_data_folder / R'1. single salts\NaCl.csv', header=None)
df_MgSO4 = pd.read_csv(input_data_folder / R'1. single salts\MgSO4.csv', header=None)
df_KCl = pd.read_csv(input_data_folder / R'1. single salts\KCl.csv', header=None)
df_HCl = pd.read_csv(input_data_folder / R'1. single salts\HCl.csv', header=None)
df_HClO4 = pd.read_csv(input_data_folder / R'1. single salts\HClO4.csv', header=None)
df_NaClO4 = pd.read_csv(input_data_folder / R'1. single salts\NaClO4.csv', header=None)
C_NaCl, D_NaCl = df_NaCl[0], df_NaCl[1]
C_MgSO4, D_MgSO4 = df_MgSO4[0], df_MgSO4[1]
C_KCl, D_KCl = df_KCl[0], df_KCl[1]
C_HCl, D_HCl = df_HCl[0], df_HCl[1]
C_HClO4, D_HClO4 = df_HClO4[0], df_HClO4[1]
C_NaClO4, D_NaClO4 = df_NaClO4[0], df_NaClO4[1]

df_HCl_HClO4 = pd.read_csv(input_data_folder / R'2. salt mixtures\HCl_HClO4.csv', header=None)
df_HCl_HNO3 = pd.read_csv(input_data_folder / R'2. salt mixtures\HCl_HNO3.csv', header=None)
df_NaCl_HClO4 = pd.read_csv(input_data_folder / R'2. salt mixtures\NaCl_HClO4.csv', header=None)
df_NaClO3_HClO4 = pd.read_csv(input_data_folder / R'2. salt mixtures\NaClO3_HClO4.csv', header=None)
df_NaNO3_HCl = pd.read_csv(input_data_folder / R'2. salt mixtures\NaNO3_HCl.csv', header=None)
df_NaAc_KClO3 = pd.read_csv(input_data_folder / R'2. salt mixtures\NaOAc_KClO3.csv', header=None)

C_HCl_HClO4, D_HCl_HClO4 = df_HCl_HClO4[0], df_HCl_HClO4[1]
C_HCl_HNO3, D_HCl_HNO3 = df_HCl_HNO3[0], df_HCl_HNO3[1]
C_NaCl_HClO4, D_NaCl_HClO4 = df_NaCl_HClO4[0], df_NaCl_HClO4[1]
C_NaClO3_HClO4, D_NaClO3_HClO4 = df_NaClO3_HClO4[0], df_NaClO3_HClO4[1]
C_NaNO3_HCl, D_NaNO3_HCl = df_NaNO3_HCl[0], df_NaNO3_HCl[1]
C_NaAc_KClO3, D_NaAc_KClO3 = df_NaAc_KClO3[0], df_NaAc_KClO3[1]
jl.seval("using Logging")
jl.seval("Logging.disable_logging(Logging.Warn)")
jl.include("04_pGM-model/code/main.jl")

# Hydrated radii
ah_SO4 = 3.79e-10
ah_Ac = 3.22e-10
ah_Cl = 2.00e-10
ah_NO3 = 1.98e-10
ah_ClO3 = 2.16e-10
ah_ClO4 = 2.83e-10
ah_NH4 = 2.30e-10
ah_K = 2.40e-10
ah_Na = 2.50e-10
ah_H = 1.97e-10
ah_Ca = 4.12e-10
ah_Mg = 4.28e-10

d = {}

# --- aa ---
print('NaCl')
pGM_NaCl, ionS_NaCl = [], []
for c in C_NaCl:
    ion_list = [[c, +1, ah_Na, "alpha"], [c, -1, ah_Cl, "alpha"]]
#     pGM = jl.main(juliacall.convert(jl.Vector[jl.Vector[jl.Float64]], ion_list), print_flag=True)
    pGM, ST, success = jl.main(ion_list, print_flag=False)
    pGM_NaCl.append(pGM)
    ionS_NaCl.append(c)
d["NaCl"] = {"pGM": pGM_NaCl, "ionS": ionS_NaCl, "D": D_NaCl}

print('MgSO4')
pGM_MgSO4, ionS_MgSO4 = [], []
for c in C_MgSO4:
    ion_list = [[c, +2, ah_Mg, "alpha"], [c, -2, ah_SO4, "alpha"]]
    pGM, ST, success = jl.main(ion_list, print_flag=False)
    pGM_MgSO4.append(pGM)
    ionS_MgSO4.append(4*c)
d["MgSO4"] = {"pGM": pGM_MgSO4, "ionS": ionS_MgSO4, "D": D_MgSO4}

print('KCl')
pGM_KCl, ionS_KCl = [], []
for c in C_KCl:
    ion_list = [[c, +1, ah_K, "alpha"], [c, -1, ah_Cl, "alpha"]]
    pGM, ST, success = jl.main(ion_list, print_flag=False)
    pGM_KCl.append(pGM)
    ionS_KCl.append(c)
d["KCl"] = {"pGM": pGM_KCl, "ionS": ionS_KCl, "D": D_KCl}

# --- ab/ba ---
print('HCl')
pGM_HCl, ionS_HCl = [], []
for c in C_HCl:
    ion_list = [[c, +1, ah_H, "proton"], [c, -1, ah_Cl, "alpha"]]
    pGM, ST, success = jl.main(ion_list, print_flag=False)
    pGM_HCl.append(pGM)
    ionS_HCl.append(c)
d["HCl"] = {"pGM": pGM_HCl, "ionS": ionS_HCl, "D": D_HCl}

print('HClO4')
pGM_HClO4, ionS_HClO4 = [], []
for c in C_HClO4:
    ion_list = [[c, +1, ah_H, "proton"], [c, -1, ah_ClO4, "beta"]]
    pGM, ST, success = jl.main(ion_list, print_flag=False)
    pGM_HClO4.append(pGM)
    ionS_HClO4.append(c)
d["HClO4"] = {"pGM": pGM_HClO4, "ionS": ionS_HClO4, "D": D_HClO4}

print('NaClO4')
pGM_NaClO4, ionS_NaClO4 = [], []
for c in C_NaClO4:
    ion_list = [[c, +1, ah_Na, "alpha"], [c, -1, ah_ClO4, "beta"]]
    pGM, ST, success = jl.main(ion_list, print_flag=False)
    pGM_NaClO4.append(pGM)
    ionS_NaClO4.append(c)
d["NaClO4"] = {"pGM": pGM_NaClO4, "ionS": ionS_NaClO4, "D": D_NaClO4}

# --- Mixtures ---
print('HCl+HClO4')
pGM_HCl_HClO4, ionS_HCl_HClO4 = [], []
for c in C_HCl_HClO4:
    ion_list = [[c, +1, ah_H, "proton"], [c, -1, ah_Cl, "alpha"],
                [c, +1, ah_H, "proton"], [c, -1, ah_ClO4, "beta"]]
    pGM, ST, success = jl.main(ion_list, print_flag=False)
    pGM_HCl_HClO4.append(pGM)
    ionS_HCl_HClO4.append(2*c)
d["HCl+HClO4"] = {"pGM": pGM_HCl_HClO4, "ionS": ionS_HCl_HClO4, "D": D_HCl_HClO4}

print('HCl+HNO3')
pGM_HCl_HNO3, ionS_HCl_HNO3 = [], []
for c in C_HCl_HNO3:
    ion_list = [[c, +1, ah_H, "proton"], [c, -1, ah_Cl, "alpha"],
                [c, +1, ah_H, "proton"], [c, -1, ah_NO3, "alpha"]]
    pGM, ST, success = jl.main(ion_list, print_flag=False)
    pGM_HCl_HNO3.append(pGM)
    ionS_HCl_HNO3.append(2*c)
d["HCl+HNO3"] = {"pGM": pGM_HCl_HNO3, "ionS": ionS_HCl_HNO3, "D": D_HCl_HNO3}

print('NaCl+HClO4')
pGM_NaCl_HClO4, ionS_NaCl_HClO4, ST_NaCl_HClO4 = [], [], []
for c in C_NaCl_HClO4:
    ion_list = [[c, +1, ah_Na, "alpha"], [c, -1, ah_Cl, "alpha"],
                [c, +1, ah_H, "proton"], [c, -1, ah_ClO4, "beta"]]
    pGM, ST, success = jl.main(ion_list, print_flag=False)
    pGM_NaCl_HClO4.append(pGM)
    ST_NaCl_HClO4.append(ST)
    ionS_NaCl_HClO4.append(2*c)
d["NaCl+HClO4"] = {"pGM": pGM_NaCl_HClO4, "ionS": ionS_NaCl_HClO4, "D": D_NaCl_HClO4, "ST": ST_NaCl_HClO4}

print('NaClO3+HClO4')
pGM_NaClO3_HClO4, ionS_NaClO3_HClO4 = [], []
for c in C_NaClO3_HClO4:
    ion_list = [[c, +1, ah_Na, "alpha"], [c, -1, ah_ClO3, "beta"],
                [c, +1, ah_H, "proton"], [c, -1, ah_ClO4, "beta"]]
    pGM, ST, success = jl.main(ion_list, print_flag=False)
    pGM_NaClO3_HClO4.append(pGM)
    ionS_NaClO3_HClO4.append(2*c)
d["NaClO3+HClO4"] = {"pGM": pGM_NaClO3_HClO4, "ionS": ionS_NaClO3_HClO4, "D": D_NaClO3_HClO4}

print('NaNO3+HCl')
pGM_NaNO3_HCl, ionS_NaNO3_HCl = [], []
for c in C_NaNO3_HCl:
    ion_list = [[c, +1, ah_Na, "alpha"], [c, -1, ah_NO3, "alpha"],
                [c, +1, ah_H, "proton"], [c, -1, ah_Cl, "alpha"]]
    pGM, ST, success = jl.main(ion_list, print_flag=False)
    pGM_NaNO3_HCl.append(pGM)
    ionS_NaNO3_HCl.append(2*c)
d["NaNO3+HCl"] = {"pGM": pGM_NaNO3_HCl, "ionS": ionS_NaNO3_HCl, "D": D_NaNO3_HCl}

print('NaAc+KClO3')
pGM_NaAc_KClO3, ionS_NaAc_KClO3 = [], []
for c in C_NaAc_KClO3:
    ion_list = [[c, +1, ah_Na, "alpha"], [c, -1, ah_Ac, "beta"],
                [c, +1, ah_K, "alpha"], [c, -1, ah_ClO3, "beta"]]
    pGM, ST, success = jl.main(ion_list, print_flag=False)
    pGM_NaAc_KClO3.append(pGM)
    ionS_NaAc_KClO3.append(2*c)
d["NaAc+KClO3"] = {"pGM": pGM_NaAc_KClO3, "ionS": ionS_NaAc_KClO3, "D": D_NaAc_KClO3}

print('end')

for d_key in d.keys():
    d_df = pd.DataFrame.from_dict(d[d_key])
    d_df.to_csv(output_data_folder / f"LitData_{d_key}.csv")
