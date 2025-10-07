import numpy as np

# Dictionary to store properties for each compound
# Source: https://chemistry.mdma.ch/hiveboard/rhodium/pdf/chemical-data/prop_aq.pdf
# Which, in turn, is based on Soehnel, O., and Novotny, P., Densities of Aqueous 
# Solutions of Inorganic Substances, Elsevier, Amsterdam, 1985.
# Columns:
# Mass% Mass of solute divided by total mass of solution, expressed as percent.
# m     Molality (moles of solute per kg of water).
# c     Molarity (moles of solute per liter of solution).
# rho   Density of solution in g/cm 3.
# n     Index of refraction, relative to air, at lambda 589 nm
# Delta Freezing point depression in °C relative to pure water.
# eta   Absolute (dynamic) viscosity in mPa s (equal to centipoise, cP)

properties = {
    "NaAc": [
        # mass%,  m,     c,    rho,      n,delta,   eta
        [0.0, 0.000, 0.000, 0.9982, 1.3330, 0.00, 1.002],
        [0.5, 0.061, 0.061, 1.0008, 1.3337, 0.22, 1.021],
        [1.0, 0.123, 0.122, 1.0034, 1.3344, 0.43, 1.040],
        [2.0, 0.249, 0.246, 1.0085, 1.3358, 0.88, 1.080],
        [3.0, 0.377, 0.371, 1.0135, 1.3372, 1.34, 1.124],
        [4.0, 0.508, 0.497, 1.0184, 1.3386, 1.82, 1.171],
        [5.0, 0.642, 0.624, 1.0234, 1.3400, 2.32, 1.222],
        [6.0, 0.778, 0.752, 1.0283, 1.3414, 2.85, 1.278],
        [7.0, 0.918, 0.882, 1.0334, 1.3428, 3.40, 1.337],
        [8.0, 1.060, 1.013, 1.0386, 1.3442, 3.98, 1.401],
        [9.0, 1.206, 1.145, 1.0440, 1.3456, 4.57, 1.468],
    ],
    "NaCl": [
        [0.0, 0.000, 0.000, 0.9982, 1.3330, 0.00, 1.002],
        [0.5, 0.086, 0.086, 1.0018, 1.3339, 0.30, 1.011],
        [1.0, 0.173, 0.172, 1.0053, 1.3347, 0.59, 1.020],
        [2.0, 0.349, 0.346, 1.0125, 1.3365, 1.19, 1.036],
        [3.0, 0.529, 0.523, 1.0196, 1.3383, 1.79, 1.052],
        [4.0, 0.713, 0.703, 1.0268, 1.3400, 2.41, 1.068],
        [5.0, 0.901, 0.885, 1.0340, 1.3418, 3.05, 1.085],
        [6.0, 1.092, 1.069, 1.0413, 1.3435, 3.70, 1.104],
    ],
    "Na2SO4": [
        [0.0, 0.000, 0.000, 0.9982, 1.3330, 0.00, 1.002],
        [0.5, 0.035, 0.035, 1.0027, 1.3338, 0.17, 1.013],
        [1.0, 0.071, 0.071, 1.0071, 1.3345, 0.32, 1.026],
        [2.0, 0.144, 0.143, 1.0161, 1.3360, 0.61, 1.058],
        [3.0, 0.218, 0.217, 1.0252, 1.3376, 0.87, 1.091],
        [4.0, 0.293, 0.291, 1.0343, 1.3391, 1.13, 1.126],
        [5.0, 0.371, 0.367, 1.0436, 1.3406, 1.36, 1.163],
        [6.0, 0.449, 0.445, 1.0526, 1.3420, 1.56, 1.202],
        [7.0, 0.530, 0.523, 1.0619, 1.3435, 0.00, 1.244],
        [8.0, 0.612, 0.603, 1.0713, 1.3449, 0.00, 1.289],
    ],
    "NH4Cl": [
        [0.0, 0.000, 0.000, 0.9982, 1.3330, 0.00, 1.002],
        [0.5, 0.094, 0.093, 0.9998, 1.3340, 0.32, 0.999],
        [1.0, 0.189, 0.187, 1.0014, 1.3349, 0.64, 0.996],
        [2.0, 0.382, 0.376, 1.0045, 1.3369, 1.27, 0.992],
        [3.0, 0.578, 0.565, 1.0076, 1.3388, 1.91, 0.988],
        [4.0, 0.779, 0.756, 1.0107, 1.3407, 2.57, 0.985],
    ],
    "NH42SO4": [
        [0.0, 0.000, 0.000, 0.9982, 1.3330, 0.00, 1.002],
        [0.5, 0.038, 0.038, 1.0012, 1.3338, 0.17, 1.008],
        [1.0, 0.076, 0.076, 1.0042, 1.3346, 0.33, 1.014],
        [2.0, 0.154, 0.153, 1.0101, 1.3363, 0.63, 1.027],
        [3.0, 0.234, 0.231, 1.0160, 1.3379, 0.92, 1.041],
        [4.0, 0.315, 0.309, 1.0220, 1.3395, 1.21, 1.057],
        [5.0, 0.398, 0.389, 1.0279, 1.3411, 1.49, 1.073],
        [6.0, 0.483, 0.469, 1.0338, 1.3428, 1.77, 1.090],
        [7.0, 0.570, 0.551, 1.0397, 1.3444, 2.05, 1.108],
    ],
    "Water": [
        [0.0, 0.000, 0.000, 0.9982, 1.3330, 0.00, 1.002],
    ],
}


def density(compound, concentration):
    """
    General function to calculate density for a given compound and concentration.

    Args:
        compound (str): The key for the compound in the properties dictionary.
        concentration (float): The concentration (mol/l) for which to calculate density.

    Returns:
        float: The interpolated density value in kg/m3.
    """
    compound_data = np.array(properties[compound])
    rho_gcm3 = np.interp(concentration, compound_data[:, 2], compound_data[:, 3])
    rho_kgm3 = rho_gcm3 * 1000
    return rho_kgm3
