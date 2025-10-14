import yaml
import numpy as np
import pandas as pd
import scipy.optimize as scopt


plot_properties_yaml = "06_figures/plot_properties.yaml"
# Load plot properties
with open(plot_properties_yaml) as prop_yaml:
    plot_props = yaml.safe_load(prop_yaml)

positions = plot_props["positions"]
r_column = plot_props["d_column"] / 2


def mean_holdup_schweitzer(holdups, positions):
    """
    Calculate the mean holdup using the Schweitzer method.

    Args:
        holdups (array-like): Holdup measurements.
        positions (array-like): Corresponding positions.

    Returns:
        float: The calculated mean holdup.
    """
    a = ((np.sum(holdups * positions**2) - (np.sum(holdups) * np.sum(positions**2) / len(holdups)))
         / (np.sum(positions**4) - (np.sum(positions**2)**2 / len(holdups))))
    b = (np.sum(holdups) - a * np.sum(positions**2)) / len(holdups)

    mean_holdup = a/2 + b

    return mean_holdup


def radial_holdup_schweitzer_normalized(r, a, b, c):
    """
    Calculate the normalized radial holdup using Schweitzer's model.

    Args:
        r (float): Normalized radius.
        a (float): Parameter a.
        b (float): Parameter b.
        c (float): Parameter c.

    Returns:
        float: The normalized radial holdup.
    """
    return a * (r**6 - 1) + b * (r**4 - 1) + c * (r**2 - 1)


def radial_holdup_schweitzer(r, mean_hu, a, b, c):
    """
    Calculate the radial holdup using Schweitzer's model.

    Args:
        r (float): Normalized radius.
        mean_hu (float): Mean holdup.
        a (float): Parameter a.
        b (float): Parameter b.
        c (float): Parameter c.

    Returns:
        float: The radial holdup.
    """
    return mean_hu * radial_holdup_schweitzer_normalized(r, a, b, c)


def schweitzer_obj(x, r, hu_meas, mean_hu):
    """
    Objective function for Schweitzer's model optimization.

    Args:
        x (array-like): Parameters [a, b, c].
        r (array-like): Normalized radii.
        hu_meas (array-like): Measured holdups.
        mean_hu (float): Mean holdup.

    Returns:
        float: The sum of squared errors.
    """
    # Extract variables
    a, b, c = x

    # Calculate radial holdup
    hu_pred = radial_holdup_schweitzer(r, mean_hu, a, b, c)

    # Sum squared error
    err_FP = np.sum((hu_meas - hu_pred)**2)
    return err_FP


def schweitzer_obj_aug(x, r, hu_meas, pp_hu):
    """
    Augmented objective function for Schweitzer's model optimization.

    Args:
        x (array-like): Parameters [a, b, c, mean_hu].
        r (array-like): Normalized radii.
        hu_meas (array-like): Measured holdups.
        pp_hu (float): Pressure probe holdup.

    Returns:
        float: The sum of squared errors (FP and PP combined).
    """
    # Extract variables
    a, b, c, mean_hu = x

    # Calculate radial holdup
    hu_pred = radial_holdup_schweitzer(r, mean_hu, a, b, c)

    # Sum squared error
    err_FP = np.sum((hu_meas - hu_pred)**2)

    # Compare mean holdup to pressure probe measurement
    r_vec = np.linspace(0, 1, 100)
    hu_vec = radial_holdup_schweitzer(r_vec, mean_hu, a, b, c)
    hu_mean_pred = 2 * np.trapezoid(hu_vec * r_vec, r_vec)
    # Squared error scaled to number of FP measurements
    err_PP = len(hu_meas) * (hu_mean_pred - pp_hu)**2

    err = err_FP + err_PP
    return err


def schweitzer_constraint(x):
    """
    Constraint function to ensure the normalized radial holdup integrates to 1.

    Args:
        x (array-like): Parameters [a, b, c] or [a, b, c, mean_hu].

    Returns:
        float: The constraint value (should be 0 for valid parameters).
    """
    if len(x) == 3:
        a, b, c = x
    elif len(x) == 4:
        a, b, c, _ = x
    r_vec = np.linspace(0, 1, 100)

    hu_vec = radial_holdup_schweitzer_normalized(r_vec, a, b, c)
    return 2 * np.trapezoid(hu_vec * r_vec, r_vec) - 1


def fit_and_find_mean(positions, holdups):
    """
    Fit Schweitzer's model to holdup data and find the mean holdup.

    Args:
        positions (array-like): Normalized positions.
        holdups (array-like): Measured holdups.

    Returns:
        tuple: Mean holdup and optimized parameters [a, b, c].
    """
    x0 = [-1.7889, 1.228, -0.939]   # Parameters from Schweitzer, 2001
    positions = np.array(positions)
    holdups = np.array(holdups)
    mean_holdup = mean_holdup_schweitzer(holdups, positions)
    min_res = scopt.minimize(
        lambda x: schweitzer_obj(x, positions, holdups, mean_holdup),
        x0,
        constraints={
            'type': 'eq',
            'fun': schweitzer_constraint
        },
        method="SLSQP"
    )
    x_opt = min_res.x

    if not min_res.success:
        print(min_res.message)
        input("Errors in holdup curve fitting...")
        raise ValueError

    return mean_holdup, x_opt


def fit_and_find_mean_aug(positions, holdups, pp_holdup, x_fit):
    """
    Fit Schweitzer's model with pressure probe augmentation.

    Args:
        positions (array-like): Normalized positions.
        holdups (array-like): Measured holdups.
        pp_holdup (float): Pressure probe holdup.
        x_fit (array-like): Initial parameters [a, b, c].

    Returns:
        tuple: Augmented mean holdup and optimized parameters [a, b, c].
    """
    x0 = [*x_fit, mean_holdup_schweitzer(holdups, positions)]
    positions = np.array(positions)
    holdups = np.array(holdups)
    min_res = scopt.minimize(
        lambda x: schweitzer_obj_aug(x, positions, holdups, pp_holdup),
        x0,
        constraints={
            'type': 'eq',
            'fun': schweitzer_constraint
        },
        method="SLSQP"
    )
    x_opt = min_res.x
    mean_holdup = x_opt[-1]

    if not min_res.success:
        print(min_res.message)
        input("Errors in holdup curve fitting...")

    return mean_holdup, x_opt[:-1]


def process_holdup_data(output_df, boi, salt, conc, flow, i):
    """
    Process holdup data for a single experiment and update the output DataFrame.

    Args:
        output_df (pd.DataFrame): The output DataFrame to update.
        boi (pd.DataFrame): Subset of the bubbles DataFrame for the experiment.
        salt (str): The salt used in the experiment.
        conc (float): The concentration (mol/l).
        flow (float): The gas flow rate (l/min).
        i (int): The index for the output DataFrame.

    Returns:
        pd.DataFrame: The updated output DataFrame.
    """
    # Add salt, concentration, flowrate to dataframe
    output_df.at[i, "Salt"] = salt
    output_df.at[i, "Concentration (mol/l)"] = conc
    output_df.at[i, "Flowrate (l/min)"] = flow
    output_df.at[i, "Ionic strength"] = boi["Ionic strength"].mean()
    output_df.at[i, "pGM"] = boi["p_GM"].mean()

    # Calculate holdup PP mean and std and add to dataframe
    PP_holdup = boi["PP holdup"].mean()
    output_df.at[i, "Holdup PP"] = PP_holdup
    output_df.at[i, "Holdup PP std"] = boi["PP holdup"].std()

    # Concatenate FP positions for holdup measurements, add to dataframe
    FP_positions = boi["FP position"].apply(
        lambda p: (positions[p] - r_column) / r_column).values
    FP_holdups = boi["FP holdup"].values
    output_df.at[i, "FP Positions"] = FP_positions
    # Concatenate FP holdup measurements, add to dataframe
    output_df.at[i, "Holdup FP (meas)"] = FP_holdups

    # Calculate FP mean holdup following Schweitzer, add value and parameters to dataframe
    FP_mean, x_opt = fit_and_find_mean(FP_positions, FP_holdups)
    output_df.at[i, "Holdup FP (Schweitzer mean)"] = FP_mean
    output_df.at[i, "x_opt (Schweitzer)"] = x_opt

    # Calculate FP mean holdup with PP augmentation, add value and parameters to dataframe
    FP_mean_aug, x_opt_aug = fit_and_find_mean_aug(FP_positions, FP_holdups, PP_holdup, x_opt)
    output_df.at[i, "Holdup FP (Augmented mean)"] = FP_mean_aug
    output_df.at[i, "x_opt (Augmented)"] = x_opt_aug

    return output_df


def bound_quantiles(data, QL, QR):
    """Crop data to quantiles QL and QR

    Args:
        data (array-like): Measurement points
        QL (int): Left quantile as integer
        QR (int): Right quantile as integer

    Returns:
        np.ndarray: Data cropped to quantiles
    """
    data = np.array(data)
    data.sort()
    # get first quantile
    Q1 = np.quantile(data, QL/100)
    # find corresponding index
    i_Q1 = np.where(data > Q1)[0][0]

    Q3 = np.quantile(data, QR/100)
    i_Q3 = np.where(data <= Q3)[0][-1]
    data = data[i_Q1:i_Q3]

    return data


def calc_d23(bsd):
    return 6 * np.sum(np.power(bsd, 2)) / np.sum(np.power(bsd, 3))


def process_a_data(output_df, boi, i):
    # Concatenate FP BSD measurements, add to dataframe
    BSD_list = boi["Valid bubbles"].values
    output_df.at[i, "BSD"] = BSD_list

    mean_hu = output_df.at[i, "Holdup FP (Augmented mean)"]
    x_opt = output_df.at[i, "x_opt (Augmented)"]
    # Create function for holdup
    hu_fun = lambda r: radial_holdup_schweitzer(r, mean_hu, *x_opt)

    # if multiple measurements are done at a single position, concatenate
    pos_BSDs = []
    unique_pos = []
    for pos in np.unique(output_df.at[i, "FP Positions"]):
        pos_loc = output_df.at[i, "FP Positions"] == pos
        if np.sum(pos_loc) > 1:
            positional_BSD = np.concatenate(output_df.at[i, "BSD"][pos_loc])
        else:
            positional_BSD = output_df.at[i, "BSD"][pos_loc][0]
        pos_BSDs.append(positional_BSD)
        unique_pos.append(pos)

    d23_meas = np.zeros(len(pos_BSDs))
    d_B_median = np.zeros(len(pos_BSDs))
    for d23_i in range(len(d23_meas)):
        d_B = pos_BSDs[d23_i]
        d23_meas[d23_i] = calc_d23(d_B)
        d_B_median[d23_i] = np.median(d_B)

    output_df.at[i, "d23"] = [d23_meas]
    output_df.at[i, "d32"] = [np.power(d23_meas, -1)]
    output_df.at[i, "dmedian"] = [d_B_median]
    # Create function for d23
    d23_fun = lambda r: d_B_step(r, d23_meas, unique_pos)

    # Multiply for range of relative positions
    r_vec = np.linspace(-1, 1, 100)

    a_vec = hu_fun(r_vec) * d23_fun(r_vec)

    # Calculate interfacial area for right part (-pi/2 < theta < pi/2)
    a_right = np.pi * np.trapezoid(a_vec[r_vec > 0] * r_vec[r_vec > 0], r_vec[r_vec > 0])
    # Calculate interfacial area for left part (pi/2 < theta -pi/2)
    a_left = - np.pi * np.trapezoid(a_vec[r_vec < 0] * r_vec[r_vec < 0], r_vec[r_vec < 0])
    a_all = a_right + a_left

    output_df.at[i, "Interfacial area (m2/m3)"] = a_all * 1e6

    # if multiple measurements are done at a single position, concatenate
    pos_BSDs = []
    unique_pos = []
    for pos in np.unique(output_df.at[i, "FP Positions"]):
        pos_loc = output_df.at[i, "FP Positions"] == pos
        if np.sum(pos_loc) > 1:
            positional_BSD = np.concatenate(output_df.at[i, "BSD"][pos_loc])
        else:
            positional_BSD = output_df.at[i, "BSD"][pos_loc][0]
        positional_BSD = bound_quantiles(positional_BSD, 1, 99)
        pos_BSDs.append(positional_BSD)
        unique_pos.append(pos)

    output_df.at[i, "BSD - Q1-Q99"] = np.array(pos_BSDs, dtype=object)

    d23_meas = np.zeros(len(pos_BSDs))
    for d23_i in range(len(d23_meas)):
        d_B = pos_BSDs[d23_i]
        d23_meas[d23_i] = calc_d23(d_B)

    output_df.at[i, "d23 - Q1-Q99"] = [d23_meas]
    output_df.at[i, "d32 - Q1-Q99"] = [np.power(d23_meas, -1)]
    # Create function for d23
    d23_fun = lambda r: d_B_step(r, d23_meas, unique_pos)

    # Multiply for range of relative positions
    r_vec = np.linspace(-1, 1, 100)

    a_vec = hu_fun(r_vec) * d23_fun(r_vec)

    # Calculate interfacial area for right part (-pi/2 < theta < pi/2)
    a_right = np.pi * np.trapezoid(a_vec[r_vec > 0] * r_vec[r_vec > 0], r_vec[r_vec > 0])
    # Calculate interfacial area for left part (pi/2 < theta -pi/2)
    a_left = - np.pi * np.trapezoid(a_vec[r_vec < 0] * r_vec[r_vec < 0], r_vec[r_vec < 0])
    a_all = a_right + a_left

    output_df.at[i, "Interfacial area (m2/m3) - Q1-Q99"] = a_all * 1e6

    return output_df


def d_B_step(r_R, d_B, rel_pos):
    """return d_B for rel_pos closest to r_R

    Args:
        r_R (float or np.ndarray): Relative position(s) to get d_B for
        d_B (np.ndarray): Bubble sizes corresponding to rel_pos
        rel_pos (np.ndarray): Relative positions corresponding to BSD measurements

    Returns:
        float or np.ndarray: Bubble sizes corresponding to r_R
    """
    assert len(d_B) == len(rel_pos), "len(d_B) =/= len(rel_pos)"
    d_B = np.array(d_B)
    rel_pos = np.array(rel_pos)
    crosses = rel_pos[:-1] + np.abs(rel_pos[1:] - rel_pos[:-1]) / 2
    crosses = [-1, *crosses, 1]
    # find i where crosses[i-1] < r_R
    i = np.searchsorted(crosses, r_R)
    # corresponds to d_median[i-1], but avoid negative values
    d_i = np.where(i-1 < 0, 0, i-1)
    return d_B[d_i]


def interfacial_area_d23(d23, hu):
    return d23 * hu


if __name__ == "__main__":
    # Load data
    bubble_h5 = pd.HDFStore("data/02_processed/consolidated_measurement_data.h5")
    bubbles_df = bubble_h5['data']
    bubble_h5.close()

    # Generate database with per-experiment row (experiment defined as one concentration / flowrate)
    # conc | flowrate | holdup (PP) | holdup p1 (FP) | ... | holdup p9 (FP) | BSD p1 (FP) | ... | BSD p9 (FP) | BSD full (FP) |
    output_df = pd.DataFrame()
    output_df["Salt"] = pd.Series(index=output_df.index, dtype=str)
    output_df["Concentration (mol/l)"] = pd.Series(index=output_df.index, dtype=float)
    output_df["Flowrate (l/min)"] = pd.Series(index=output_df.index, dtype=float)
    output_df["Ionic strength"] = pd.Series(index=output_df.index, dtype=float)
    output_df["pGM"] = pd.Series(index=output_df.index, dtype=float)
    output_df["Holdup PP"] = pd.Series(index=output_df.index, dtype=float)
    output_df["Holdup PP std"] = pd.Series(index=output_df.index, dtype=float)
    # list of the FP positions for each measurement
    output_df["FP Positions"] = pd.Series(index=output_df.index, dtype=object)
    # list of holdup measurement values
    output_df["Holdup FP (meas)"] = pd.Series(index=output_df.index, dtype=object)
    # List of lists of bubble size values for each measurement
    output_df["BSD"] = pd.Series(index=output_df.index, dtype=object)
    output_df["BSD - Q1-Q99"] = pd.Series(index=output_df.index, dtype=object)
    output_df["d32"] = pd.Series(index=output_df.index, dtype=object)
    output_df["d32 - Q1-Q99"] = pd.Series(index=output_df.index, dtype=object)
    output_df["d23"] = pd.Series(index=output_df.index, dtype=object)
    output_df["d23 - Q1-Q99"] = pd.Series(index=output_df.index, dtype=object)
    output_df["dmedian"] = pd.Series(index=output_df.index, dtype=object)
    output_df["Interfacial area (m2/m3)"] = pd.Series(index=output_df.index, dtype=float)
    output_df["Interfacial area (m2/m3) - Q1-Q99"] = pd.Series(index=output_df.index, dtype=float)
    # Mean FP holdup calculated following Schweitzer, 2001
    output_df["Holdup FP (Schweitzer mean)"] = pd.Series(index=output_df.index, dtype=float)
    # Mean FP holdup combining method of Schweitzer with pressure probe data
    output_df["Holdup FP (Augmented mean)"] = pd.Series(index=output_df.index, dtype=float)
    # Optimized parameters for the two methods
    output_df["x_opt (Schweitzer)"] = pd.Series(index=output_df.index, dtype=object)   # contains (a, b, c)
    output_df["x_opt (Augmented)"] = pd.Series(index=output_df.index, dtype=object)    # contains (a, b, c)
    # Go through bubbles dataframe
    # First, find all water experiments
    # where concentration is 0.0 and ignore is False
    i = 0
    conc = 0.0
    ionic_strength = 0.0
    p_GM = 0.0
    flow = 100
    salt = "Water"
    boi = bubbles_df.loc[
        (bubbles_df["Concentration (mol/l)"] == 0.0)
        & (bubbles_df["Ignore"] == False)   # Linters complain about `== False`. It's needed here.
        & (bubbles_df["Gas flow (l/min)"] == flow)
    ]
    output_df = process_holdup_data(output_df, boi, salt, conc, flow, i)
    output_df = process_a_data(output_df, boi, i)
    i += 1

    salts = ["NaCl", "NH4Cl", "NaAc", "NH42SO4", "Na2SO4"]
    # For each salt
    for salt in salts:
        salt_boi = bubbles_df.loc[
            (bubbles_df["Salt"] == salt)
            & (bubbles_df["Ignore"] == False)
            & (bubbles_df["Gas flow (l/min)"] == flow)
        ]
        # find unique concentrations > 0.0
        salt_boi = salt_boi.loc[salt_boi["Concentration (mol/l)"] > 0.0]
        unique_c = salt_boi["Concentration (mol/l)"].unique()
        for conc in unique_c:
            boi = salt_boi.loc[salt_boi["Concentration (mol/l)"] == conc]
            output_df = process_holdup_data(output_df, boi, salt, conc, flow, i)
            output_df = process_a_data(output_df, boi, i)
            i += 1

    h5Store = pd.HDFStore('data/02_processed/per-experiment-data.h5', 'w')

    h5Store['data'] = output_df
    h5Store.close()
