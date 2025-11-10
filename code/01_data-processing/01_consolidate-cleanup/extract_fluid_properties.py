import warnings
import numpy as np
import pandas as pd
import juliacall
from juliacall import Main as jl

from helpers.density_interpolator import density


def ionic_strength(conc, salt):
    if salt in ["NaCl", "NH4Cl", "NaAc"]:
        c_I = conc
    elif salt in ["Na2SO4", "NH42SO4"]:
        c_cat = 2 * conc
        c_SO4 = conc
        c_I = .5 * ((c_cat * 1**2) + (c_SO4 * 2**2))
    else:
        warnings.warn(f"Unknown salt used: {salt}. Assuming ionic strength equal to concentration")
        c_I = conc
    return c_I


def molar_to_grams_per_liter(c_mol, salt):
    global MW
    if salt == "Water":
        return 0.0
    c_grams = c_mol * MW[salt]
    return c_grams


def ion_concentrations(conc, salt):
    """Convert molar concentration to concentration per ion

    Args:
        conc (float): Molar concentration of the salt
        salt (str): String representation of the salt

    Raises:
        ValueError: Unknown salt string provided

    Returns:
        dict: Dictionary with a concentration entry for each ion
    """
    c_ion = {
        'Na': 0.0,
        'NH4': 0.0,
        'Cl': 0.0,
        'SO4': 0.0,
        'Ac': 0.0,
    }
    if salt == "NaCl":
        c_ion['Na'] += conc
        c_ion['Cl'] += conc
    elif salt == "NH4Cl":
        c_ion['NH4'] += conc
        c_ion['Cl'] += conc
    elif salt == "Na2SO4":
        c_ion['Na'] += 2 * conc
        c_ion['SO4'] += conc
    elif salt == "NH42SO4":
        c_ion['NH4'] += 2 * conc
        c_ion['SO4'] += conc
    elif salt == "NaAc":
        c_ion['Na'] += conc
        c_ion['Ac'] += conc
    elif salt == "Water":
        pass
    else:
        raise ValueError(f"Did not recognize salt: {salt}")
    return c_ion


def conc_to_pGM(c_ion, salt, c_mol):
    if salt == "Water":
        return 0
    if c_mol == 0.0:
        return 0
    ion_concentrations = [c_ion['Na'], c_ion['Cl'], c_ion['NH4'], c_ion['SO4'], c_ion['Ac']]
    ion_charges = {
        'Na': +1,
        'Cl': -1,
        'NH4': +1,
        'SO4': -2,
        'Ac': -1
    }
    hydrated_radii = {
        'Na': 2.5e-10,
        'Cl': 2.0e-10,
        'NH4': 2.5e-10,
        'SO4': 3.79e-10,
        'Ac': 3.22e-10
    }
    ion_types = {
        'Na': "alpha",
        'Cl': "alpha",
        'NH4': "alpha",
        'SO4': "alpha",
        'Ac': "beta"
    }
    ion_list = np.empty(len(ion_concentrations), dtype=tuple)
    for i, ion in enumerate(c_ion.keys()):
        ion_list[i] = (c_ion[ion], ion_charges[ion], hydrated_radii[ion], ion_types[ion])
    total_charge = 0.0
    for ion in ion_list:
        total_charge += ion[0] * ion[1]
    if not np.isclose(total_charge, 0, atol=1e-5):
        raise ValueError(f"Check your input. Sum of charges is not 0 but {total_charge}")
    p_GM, _, success = jl.main(
        juliacall.convert(jl.Vector[jl.Tuple[jl.Float64, jl.Float64, jl.Float64, jl.String]], ion_list),
        print_flag=False
    )
    if not success:
        p_GM = None

    return p_GM


MW = {
    "NaCl": 58.4,
    "NH4Cl": 53.49,
    "Na2SO4": 142.04,
    "NH42SO4": 132.14,
    "NaAc": 82.0,
}

if __name__ == "__main__":
    # Load bubbles dataframe
    input_df = pd.read_csv(
        'data/02_processed/consolidated-experiment-logs.csv',
        dtype={
            "Exp. No.": int,
            "Salt": str,
            "Concentration (mol/l)": float,
            "FP position": int,
            "Gas flow (l/min)": float,
            "Date": str,
            "PP Timestamp": str,
            "FP Timestamp": str,
            "Comment": str,
            "Ignore": bool
        }
    )

    output_df = input_df.copy()
    print("Adding fluid properties")

    output_df['Concentration (g/l)'] = pd.Series(dtype='float', index=output_df.index)
    output_df['Ionic strength'] = pd.Series(dtype='float', index=output_df.index)
    output_df['p_GM'] = pd.Series(dtype='float', index=output_df.index)
    output_df['Na (mol/l)'] = pd.Series(dtype='float', index=output_df.index)
    output_df['NH4 (mol/l)'] = pd.Series(dtype='float', index=output_df.index)
    output_df['Cl (mol/l)'] = pd.Series(dtype='float', index=output_df.index)
    output_df['SO4 (mol/l)'] = pd.Series(dtype='float', index=output_df.index)
    output_df['Ac (mol/l)'] = pd.Series(dtype='float', index=output_df.index)
    output_df['Density'] = pd.Series(dtype='float', index=output_df.index)
    output_df['Viscosity'] = pd.Series(dtype='float', index=output_df.index)
    output_df['Surface tension'] = pd.Series(dtype='float', index=output_df.index)

    jl.println("Printing from Julia!")
    jl.include("code/03_pGM-model/code/main.jl")
    for i, exp in output_df.iterrows():
        c_mol = exp["Concentration (mol/l)"]
        salt = exp["Salt"]

        # ionic strength
        c_I = ionic_strength(c_mol, salt)
        output_df.at[i, "Ionic strength"] = c_I

        # Density
        rho = density(salt, c_mol)
        output_df.at[i, "Density"] = rho

        # concentration as g/l
        c_gl = molar_to_grams_per_liter(c_mol, salt)
        output_df.at[i, "Concentration (g/l)"] = c_gl

        # Individual ion concentrations
        c_ion = ion_concentrations(c_mol, salt)
        output_df.at[i, "Na (mol/l)"] = c_ion['Na']
        output_df.at[i, "NH4 (mol/l)"] = c_ion['NH4']
        output_df.at[i, "Cl (mol/l)"] = c_ion['Cl']
        output_df.at[i, "Ac (mol/l)"] = c_ion['Ac']
        output_df.at[i, "SO4 (mol/l)"] = c_ion['SO4']

        # Gibbs-Marangoni pressure
        p_GM = conc_to_pGM(c_ion, salt, c_mol)
        output_df.at[i, "p_GM"] = p_GM
        if p_GM is None:
            print(f"No pGM for experiment {i}")
        print(f"\rDone fluid properties of experiment {i}", end="", flush=True)
    print("\n")

    output_df.to_csv("data/02_processed/exp-with-fluid-properties.csv")
