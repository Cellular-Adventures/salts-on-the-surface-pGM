# Pressure probe calibration
import yaml
import numpy as np
import numpy.linalg as linalg
from nptdms import TdmsFile
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator

g = 9.81    # [m/s2]
density_water = 998.19      # [kg/m3]
input_file = "01_data/01_raw/2503_PP_setup.yaml"
output_file = "01_data/02_processed/probe-characterization.yaml"

# import info from setup.yaml
with open(input_file) as setup_file:
    settings = yaml.safe_load(setup_file)

database = {}
channel_names = {}
probe_characterization = {}

# create placeholder for each channel
for calibration in settings['calibrations']:
    for probe in calibration['probes']:
        if probe['channel'] not in database.keys():
            # add channel to placeholders
            database[probe['channel']] = {
                'voltage': [],
                'height': [],
                'hs_pressure': [],
                'mean_voltage': [],
                'mean_height': [],
                'mean_press': []
            }
            # add channel/name combo to dictionary
            channel_names[probe['channel']] = probe['name']

for calibration in settings['calibrations']:
    # Create database of probe height per channel
    probe_height = {}
    for probe in calibration['probes']:
        probe_height[probe['channel']] = probe['height']

    for measurement in calibration['measurements']:
        # Construct filename
        tdms_filename = (
            calibration['dir'] +
            calibration['file']['prefix'] +
            measurement['timestamp'] +
            calibration['file']['postfix']
        )

        with TdmsFile.open(tdms_filename) as tdms_file:
            # There is only 1 group of measurements per file (tdms shenanigans)
            group = tdms_file.groups()[0]
            # Loop through probes for easier calculation of relative height
            for probe in calibration['probes']:
                rel_height = max(0, measurement['waterlevel'] - probe_height[probe['channel']])
                voltage = list(group[probe['channel']][0:27000])
                database[probe['channel']]['voltage'] += voltage
                database[probe['channel']]['height'] += [rel_height] * len(voltage)
                database[probe['channel']]['mean_voltage'] += [np.mean(voltage)]
                database[probe['channel']]['mean_height'] += [rel_height]

for probe in database.keys():
    database[probe]['hs_pressure'] = list(np.array(database[probe]['height']) / 1000 * g * density_water)
    database[probe]['mean_press'] = list(np.array(database[probe]['mean_height']) / 1000 * g * density_water)

    x_meas = database[probe]['voltage']

    A = np.vstack((x_meas, np.ones(len(x_meas)))).T

    y_meas = database[probe]['hs_pressure']

    sol = linalg.lstsq(A, y_meas, rcond=None)
    beta = sol[0]
    residual = sol[1]
    # print(beta.shape)
    # print(residual.shape)

    print('\nbeta = ')
    print(beta)

    x_plot = np.linspace(min(x_meas), max(x_meas), 1000)

    # statistics
    R2 = 1 - residual / sum((y_meas - np.mean(y_meas))**2)

    probe_characterization[probe] = {
        'linear': {
            'a': float(beta[0]),
            'b': float(beta[1])
        },
        'R2': float(R2[0]),
        'range': {
            'min': float(min(x_meas)),
            'max': float(max(x_meas))
        },
        'height': probe_height[probe]
    }

    # figure of regression
    plt.figure(figsize=(4, 2))
    plt.plot(x_plot, beta[0] * x_plot + beta[1], label=f"{beta[0]:.2f} x + {beta[1]:.2f}, $R^2$ = {R2[0]:.3f}")
    plt.plot(database[probe]['mean_voltage'], database[probe]['mean_press'], "o", label="Mean")
    plt.plot(x_meas, y_meas, "x", alpha=0.01)
    plt.legend()
    plt.xlabel("Probe output (V)")
    plt.ylabel("$p_{hydrostatic}$ (Pa)")
    plt.title("Regression of voltage vs. pressure")
    # plt.ylim([0, 1500])
    plt.tight_layout()
    plt.savefig(f"06_figures/PP_calibration/regression_{channel_names[probe]}.png", dpi=600)
    
    # extrapolation to full measurement range
    x_plot_extra = np.linspace(min(x_meas), 10, 1000)

    fig, ax = plt.subplots(figsize=(4, 2))
    ax.plot(x_meas, y_meas, "x")
    ax.plot(x_plot_extra, beta[0] * x_plot_extra + beta[1], label=f"{beta[0]:.2f} x + {beta[1]:.2f}")
    # plt.ylim([0, 10_000])
    # plt.xlim([0, 10])
    plt.grid(True, 'major')
    plt.grid(True, "minor", linestyle=":", color="silver")
    # ax.yaxis.set_major_locator(MultipleLocator(5000))
    # ax.yaxis.set_minor_locator(MultipleLocator(1000))
    # ax.xaxis.set_minor_locator(MultipleLocator(1))
    plt.xlabel("Probe output (V)")
    plt.ylabel("Gauge pressure (Pa)")
    plt.title("Extrapolated 0-10 V")
    plt.tight_layout()
    plt.savefig(f"06_figures/PP_calibration/extrapolation_{channel_names[probe]}.png", dpi=600)

with open(output_file, 'w+') as output_yaml:
    yaml.dump(probe_characterization, output_yaml)

plt.show()
