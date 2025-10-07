import pandas as pd
from pathlib import Path
from datetime import datetime, timedelta


def correct_input(fp, pp, field):
    return input(
        f"Found conflicting {field} for FP {fp['timestamp']} with PP {pp['timestamp']}\n"
        f"FP {field}: {fp[field]}\n"
        f"PP {field}: {pp[field]}\n"
        f"Please enter value to use (-1 to skip entry): "
    )


def extract_field(fp, pp, field, requested_type):
    if fp[field] == pp[field]:
        field_value = fp[field]
    else:
        field_value = correct_input(fp, pp, field)
    if requested_type is float:
        return float(field_value)
    elif requested_type is str:
        return str(field_value)
    return field_value


def compare_entries(fp, pp):
    # flowrate
    flowrate = extract_field(fp, pp, 'flowrate', float)

    # substance
    substance = extract_field(fp, pp, 'substance', str)

    # concentration
    concentration = extract_field(fp, pp, 'concentration', float)

    return flowrate, substance, concentration


def check_inclusion(fp, pp):
    if isinstance(fp['comments'], str) or isinstance(pp['comments'], str):
        combo_comment = f"FP: {fp['comments']}.\nPP: {pp['comments']}"
        if "SKIP" in combo_comment:
            flag = True
            print(f"Found 'SKIP' in comments. Skipping {fp['timestamp']}.")
        else:
            flag = not bool(int(input(f"\nFound comments:\n"
                                    f"{combo_comment}\nInclude in analyses? (1 for yes, 0 for no) ")))
            print("")
    else:
        combo_comment = ''
        flag = False
    return combo_comment, flag


if __name__ == "__main__":
    root_folder = Path("01_data/01_raw/02_unpacked")

    output_file = "01_data/02_processed/consolidated-experiment-logs.csv"

    # Dates of the folders with bubble measurements
    meas_dates = [
        20250321,
        20250324,
        20250326,
        20250327,
        20250328,
        20250331,
        20250401,
        20250403,
        20250404,
        20250407,
        20250408,
        20250409,
        20250410,
        20250411,
        20250414,
        20250417,
        20250422,
    ]

    # load input.csv into a pd dataframe
    if Path(output_file).exists():
        # Load existing input_automated.csv
        output = pd.read_csv(output_file, index_col=0,
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
                             })
    else:
        # Empty dataframe
        output = pd.DataFrame()
        output["Exp. No."] = pd.Series(dtype=int, index=output.index)
        output["Salt"] = pd.Series(dtype=str, index=output.index)
        output["Concentration (mol/l)"] = pd.Series(dtype=float, index=output.index)
        output["FP position"] = pd.Series(dtype=int, index=output.index)
        output["Gas flow (l/min)"] = pd.Series(dtype=float, index=output.index)
        output["Date"] = pd.Series(dtype=str, index=output.index)
        output["PP Timestamp"] = pd.Series(dtype=str, index=output.index)
        output["FP Timestamp"] = pd.Series(dtype=str, index=output.index)
        output["Datetime"] = pd.Series(dtype=str, index=output.index)
        output["Comment"] = pd.Series(dtype=str, index=output.index)
        output["Ignore"] = pd.Series(dtype=bool, index=output.index)
    n_exp = 0
    new_info = False
    for meas_date in meas_dates:
        # Convert meas_date to datetime object
        date = datetime.strptime(str(meas_date), "%Y%m%d")
        print(f"Analyzing experiments on {date.strftime('%d/%m/%Y')}\n")

        day_folder = sorted(root_folder.glob(f"{str(meas_date)}*"))
        if len(day_folder) > 1:
            raise ValueError("Multiple folders for the same day found."
                            "I wasn't meant tO DEAL WITH THIS!")
        # load fp_logs
        fp_path = day_folder[0] / "FP" / (date.strftime("%m%d") + "_FP Log.csv")
        fp_log = pd.read_csv(fp_path)
        # load pp_logs
        pp_path = day_folder[0] / "PP" / (date.strftime("%m%d") + "_PP Log.csv")
        pp_log = pd.read_csv(pp_path)
        for i, fp_exp in fp_log.iterrows():
            # Sometimes there's more experiments in fp log than in pp log. Avoid errors.
            if i >= len(pp_log):
                i = len(pp_log) - 1

            # Check for existence of input entry with matching n_exp and fp timestamp
            if not output[
                (output["Exp. No."] == n_exp) &
                (output["FP Timestamp"] == str(fp_exp["timestamp"]))
            ].empty:
                print(f"Skipping existing entry for Exp. No. {n_exp} and "
                      f"FP Timestamp {fp_exp['timestamp']}")
                n_exp += 1
                continue
            else:
                new_info = True

            # Find pp log within 30 seconds of fp_log stamp
            pp_timestamp = pp_log.iloc[i]['timestamp']
            pp_exp = pp_log.loc[pp_log["timestamp"] == pp_timestamp]

            fp_datetime = datetime.strptime(str(fp_exp['timestamp']), '%H%M%S')
            pp_datetime = datetime.strptime(str(pp_timestamp), '%H%M%S')
            # Compare based on timedelta, to deal well with post-minute differences 
            if abs(fp_datetime - pp_datetime) > timedelta(seconds=30):
                pp_exp = pp_log.loc[abs(
                    pp_log["timestamp"].apply(lambda ts: datetime.strptime(str(ts), '%H%M%S'))
                    - fp_datetime)
                    < timedelta(seconds=30)]
            if pp_exp.empty:
                pp_timestamp = int(input(
                    f"Could not find a matching PP exp for FP "
                    f"{fp_exp['timestamp']}\nComment: '"
                    f"{fp_exp['comments']}'\nProvide timestamp for PP: "))
                print("")
                pp_exp = pp_log.loc[pp_log["timestamp"] == pp_timestamp]
            else:
                pp_timestamp = pp_exp["timestamp"].iloc[0]
            pp_exp = pp_exp.iloc[0]
            fp_time = fp_datetime.time()
            pp_time = datetime.strptime(str(pp_timestamp), '%H%M%S').time()

            # Check if flowrate, substance and concentration match
            flowrate, substance, concentration = compare_entries(fp_exp, pp_exp)
            # User indicates pass with -1 value
            if flowrate == -1 or substance == -1 or concentration == -1:
                pass
            # check for inclusion
            combo_comment, ignore_flag = check_inclusion(fp_exp, pp_exp)
            series_index = len(output)
            # Add entry to output
            output.at[series_index, "Exp. No."] = int(n_exp)
            output.at[series_index, "Salt"] = str(substance)
            output.at[series_index, "Concentration (mol/l)"] = float(concentration)
            output.at[series_index, "FP position"] = int(fp_exp['position'])
            output.at[series_index, "Gas flow (l/min)"] = float(flowrate)
            output.at[series_index, "Date"] = date.strftime("%Y%m%d")
            output.at[series_index, "PP Timestamp"] = str(pp_timestamp)
            output.at[series_index, "FP Timestamp"] = str(fp_exp['timestamp'])
            output.at[series_index, "Datetime"] = datetime.combine(date, fp_time).isoformat()
            output.at[series_index, "Comment"] = str(combo_comment)
            output.at[series_index, "Ignore"] = bool(ignore_flag)
            n_exp += 1

    # Write updated input dataframe to input.csv. Overwrite old one.
    # input_df = pd.DataFrame(input_dict)

    print(output.head())

    if new_info:
        output.to_csv(output_file, lineterminator='\r\n')
        print("Expanded input.\n")
    else:
        print("No new info. No new file.")

    print("All done.")
