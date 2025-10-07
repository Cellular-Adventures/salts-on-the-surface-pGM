import csv
import yaml
import warnings
import numpy as np
import xml.etree.ElementTree as ET
import pandas as pd
from nptdms import TdmsFile
from pathlib import Path
# Reads raw data (voltage) from U:drive 
# saves summary (bubbles / holdup) in data/preprocessed

root_path = Path("01_data/01_raw/02_unpacked")
print(f"\nRoot: {root_path}")


def get_valid_bubbles(evt_path):
    valid_bubbles = []
    with open(evt_path) as evt_file:
        evt_reader = csv.DictReader(evt_file, delimiter="\t")
        for event in evt_reader:
            if int(event['Valid']):
                # convert to float (first change decimal to point)
                size_float = float(event['Size'].replace(',', '.'))
                # Super tiny bubble events - I don't trust them if based on entry
                if size_float < 10 and event['VeloOut'] == "-1":
                    continue
                valid_bubbles.append(size_float)
                if size_float < 10:
                    print(f"Tiny bubble: {size_float:.2} micrometer")
    return np.array(valid_bubbles)


def bin_bubbles(bubbles):
    # bins should always be equal size.
    bin_size = 100
    if any(bubbles <= 0):
        raise ValueError("Negative or 'zero' bubble sizes found...")
    # bins should never cut off the largest / smallest bubbles
    val_max = int(max(max(bubbles), 50_000))
    if val_max > 50_000:
        print("Big boi found!")
        print(max(bubbles))
    val_min = 0

    bin_edges = np.arange(val_min, val_max+2*bin_size, bin_size)
    binned_bubbles = np.histogram(bubbles, bins=bin_edges, density=False)
    return binned_bubbles


def get_holdup(evtlog):
    logfile = ET.parse(evtlog)
    root = logfile.getroot()
    if 'processingMeanVoidFraction' not in root.keys():
        warnings.warn(f"No void fraction found in {evtlog}")
        return -1
    holdup = float(root.attrib['processingMeanVoidFraction'])
    return holdup


def get_holdup_signal(stream, log):
    # stream is a tab-separated file
    # load as pd dataframe and immediately convert to np
    bubbles_df = pd.read_csv(stream, delimiter='\t', header=0)
    values = np.array([bubbles_df['Arrival'], bubbles_df['Duration']])

    # Get values from log file.
    logfile = ET.parse(log)
    root = logfile.getroot()
    frequency = float(root[0].attrib['streamingFrequency'])
    duration = float(root[0].attrib['streamingDuration'])
    params = {
        'f': frequency,     # frequency in Hz
        't': duration,      # duration in s
    }

    return values, params


def process_FP(fp_folder, input):
    evt_file = list(fp_folder.glob(f"*{input['FP Timestamp']}.evt"))

    if len(evt_file) == 0:
        raise AssertionError(f"Found no .evt files with the right timestamp: "
                             f"{input['Date']} - {input['FP Timestamp']}")
    if len(evt_file) > 1:
        raise AssertionError(f"Found multiple .evt files with the right timestamp: "
                             f"{input['Date']} - {input['FP Timestamp']}")
    evt_file = evt_file[0]

    valid_bubbles = get_valid_bubbles(evt_file)

    evtlog_file = list(fp_folder.glob(f"*{input['FP Timestamp']}.evtlog"))[0]
    # some files don't have holdup. Skip those. Figure way to indicate.
    FP_holdup = get_holdup(evtlog_file)
    if FP_holdup == -1:
        FP_holdup_signal = (-1, {'f': -1, 't': -1})
    else:
        hu_stream_file = list(fp_folder.glob(f"*{input['FP Timestamp']}_stream.evt"))[0]
        hu_log_file = list(fp_folder.glob(f"*{input['FP Timestamp']}.binlog"))[0]
        FP_holdup_signal = get_holdup_signal(hu_stream_file, hu_log_file)

    return valid_bubbles, FP_holdup, FP_holdup_signal


def process_PP(pp_folder, input: pd.Series):
    max_points = 120_000
    if input['PP Timestamp'] in ["91203", "091203"]:
        # something went wrong in this one...
        return np.nan

    pp_path = list(pp_folder.glob(f"*{input['PP Timestamp']}*.tdms"))[0]
    probes = list(CALIBRATION.keys())
    with TdmsFile.open(pp_path) as pp_file:
        group = pp_file.groups()[0]
        # Store voltage for air measurements
        voltages = np.zeros(2)
        pressures = np.zeros(2)
        for i, probe in enumerate(probes):
            # Read in a maximum of 120 000 points (2 minutes)
            # Some measurements ran for too long, this avoids the tails with gas flow changes
            voltages[i] = group[probe][:max_points].mean()
            min_voltage = CALIBRATION[probe]['range']['min']
            max_voltage = CALIBRATION[probe]['range']['max']
            probe_a = CALIBRATION[probe]['linear']['a']
            probe_b = CALIBRATION[probe]['linear']['b']
            if voltages[i] < min_voltage:
                warnings.warn(f"\nPP measurement {input['PP Timestamp']} dips " +
                              "below minimum calibrated voltage\n")
            if voltages[i] > max_voltage:
                warnings.warn(f"\nPP measurement {input['PP Timestamp']} exceeds " +
                              "maximum calibrated voltage\n")
            pressures[i] = voltages[i] * probe_a + probe_b

    if "Density" not in input:
        rho_L = 1000
    elif np.isnan(input['Density']):
        rho_L = 1000
    else:
        rho_L = input["Density"]

    d_H = abs(CALIBRATION[probes[0]]['height']
              - CALIBRATION[probes[1]]['height']) / 1000
    holdup = 1 - (pressures[0] - pressures[1]) / (9.81 * rho_L * d_H)

    return holdup


def pre_process(input_line):
    # date
    date = input_line.loc['Date']

    # open corresponding folder
    date_folder = root_path.glob(f'{date}*')
    date_folder = list(date_folder)     # We don't want fancy generator properties

    assert len(date_folder) == 1, "Found multiple folders for the same date"
    date_folder = date_folder[0]

    valid_bubbles, FP_holdup, FP_holdup_signal = process_FP(
        date_folder / "FP",
        input_line)

    PP_holdup = process_PP(date_folder / "PP", input_line)

    return valid_bubbles, FP_holdup, FP_holdup_signal, PP_holdup


MW = {
    "NaCl": 58.4,
    "NH4Cl": 53.49,
    "Na2SO4": 142.04,
    "NH42SO4": 132.14,
    "NaAc": 82.0,
}


with open("01_data/02_processed/probe-characterization.yaml") as cal_yaml:
    CALIBRATION = yaml.safe_load(cal_yaml)

input_df = pd.read_csv(
    "01_data/02_processed/exp-with-fluid-properties.csv",
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
        "Ignore": bool,
        "Concentration (g/l)": float,
        "Ionic strength": float,
        "p_GM": float,
        "Na (mol/l)": float,
        "NH4 (mol/l)": float,
        "Cl (mol/l)": float,
        "SO4 (mol/l)": float,
        "Ac (mol/l)": float,
        "Density": float,
        "Viscosity": float,
        "Surface tension": float,
    }
)

output_df = input_df.copy()
output_df['Valid bubbles'] = pd.Series(dtype='object', index=output_df.index)
output_df['FP holdup'] = pd.Series(dtype='float', index=output_df.index)
output_df['FP holdup bubbles'] = pd.Series(dtype='object', index=output_df.index)
output_df['FP holdup frequency'] = pd.Series(dtype='float', index=output_df.index)
output_df['FP holdup duration'] = pd.Series(dtype='float', index=output_df.index)
output_df['PP holdup'] = pd.Series(dtype='float', index=output_df.index)

for i, row in input_df.iterrows():
    if row["Ignore"]:
        continue
    print(f"\rReading {row['Date']} - {row['FP Timestamp']}...   ", end='')
    valid_bubbles, FP_holdup, FP_holdup_signal, PP_holdup = pre_process(row)
    FP_holdup_bubbles = FP_holdup_signal[0]
    FP_holdup_frequency = FP_holdup_signal[1]['f']
    FP_holdup_duration = FP_holdup_signal[1]['t']
    if len(valid_bubbles) == 0:
        print(f"Empty bubbles for {row['Date']} - FP {row['FP Timestamp']}. Alarm!")
        input()
    output_df.at[i, 'Valid bubbles'] = valid_bubbles
    output_df.at[i, 'FP holdup'] = FP_holdup
    output_df.at[i, 'FP holdup bubbles'] = FP_holdup_bubbles
    output_df.at[i, 'FP holdup frequency'] = FP_holdup_frequency
    output_df.at[i, 'FP holdup duration'] = FP_holdup_duration
    output_df.at[i, 'PP holdup'] = PP_holdup

print("\rRead all files.                      \n")
h5Store = pd.HDFStore('01_data/02_processed/consolidated_measurement_data.h5', 'w')

h5Store['data'] = output_df
h5Store.close()
