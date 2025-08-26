

# Imports converters
import converters.behavior_to_nwb
import converters.nwb_saving
import converters.general_to_nwb
import converters.Initiation_nwb
import converters.acquisition_to_nwb
import converters.intervals_to_nwb

# Imports modules
from typing import List, Optional, Tuple, Set
from pynwb import NWBHDF5IO, validate
from pathlib import Path
from tqdm import tqdm
import pandas as pd
import numpy as np
import importlib
import argparse
import os
import gc



############################################################
# Function that converts PL data to NWB format
############################################################

def convert_data_to_nwb_pl(output_folder,Folder_sessions_info,Folder_general_info,mouses_name = None ) :
    """
        Convert PL session data into validated NWB files.

        The function:
        - Finds matching processed/raw files for given mice.
        - Builds session DataFrames and Creates NWB files with metadata, LFP signals, intervals, and behavior.
        - Saves and validates them.

        Parameters
        ----------
        output_folder : str
            Destination folder for NWB files.
        Folder_sessions_info : str
            Path to session information.
        Folder_general_info : str
            Path to general information.
        mouses_name : list of str, optional
            Mouse names to include (default: predefined PL200–PL225).

        Returns
        -------
        None
            NWB files are written to disk; errors are printed.
        """
    # Reload the converters module
    importlib.reload(converters.Initiation_nwb)
    importlib.reload(converters.Initiation_nwb)
    importlib.reload(converters.acquisition_to_nwb)
    importlib.reload(converters.general_to_nwb)
    importlib.reload(converters.intervals_to_nwb)
    importlib.reload(converters.behavior_to_nwb)
    importlib.reload(converters.nwb_saving)

    # Default mouse list if none provided
    if mouses_name is None:
        mouses_name = [
            "PL200", "PL201", "PL202", "PL203", "PL204", "PL205", "PL206", "PL207", "PL208",
            "PL209", "PL210", "PL211", "PL212", "PL213", "PL214", "PL215", "PL216", "PL217",
            "PL218", "PL219", "PL220", "PL221", "PL222", "PL223", "PL224", "PL225"
        ]


    print("**************************************************************************")
    print("-_-_-_-_-_-_-_-_-_-_-_-_-_-_- NWB conversion _-_-_-_-_-_-_-_-_-_-_-_-_-_-_")

    # Find pairs of processed/raw data files
    pairs = converters.Initiation_nwb.find_pl_pairs(Folder_general_info, Folder_sessions_info, mouse_names=mouses_name)
    failures: List[Tuple[str, str]] = []   # (session id, error message)
    seen_mice = set()            # keep track of which mice were processed
    print("Converting data to NWB format for mouse: ", mouses_name)
    pbar = tqdm(pairs, unit="file", desc="Processing...")
    for General_path, sessions_path in pbar:
        session_id = Path(General_path).name
        pbar.set_postfix_str("Loading sessions for:" + session_id)

        # Build a small DataFrame for this session only
        try:
            sessions_df_for_this_pair = converters.Initiation_nwb.files_to_dataframe(
                General_path, sessions_path, pd.DataFrame(columns=["Mouse Name"])
            )
        except Exception as e:
            failures.append((session_id, f"files_to_dataframe: {e}"))
            gc.collect()
            continue

        # 2) Filter by requested mice
        if "Mouse Name" not in sessions_df_for_this_pair.columns:
            failures.append((session_id, "Missing 'Mouse Name' column in session DataFrame"))
            gc.collect()
            continue

        sessions_df_for_this_pair = sessions_df_for_this_pair[sessions_df_for_this_pair["Mouse Name"].isin(mouses_name)]
        if sessions_df_for_this_pair.empty:
            # Skip if no relevant mouse in this session
            gc.collect()
            continue

        # 3) Process each row (one session may produce multiple rows)
        for _, row in sessions_df_for_this_pair.iterrows():
            mouse = str(row.get("Mouse Name", "Unknown"))
            pbar.set_postfix_str("Mouse :" + mouse + " Session :" + row.get("Session", "Unknown"))
            seen_mice.add(mouse)

            # Determine behavior type
            behavior = str(row.get("Behavior Type", ""))
            if behavior == "Detection Task":
                rewarded = True
            elif behavior == "Neutral Exposition":
                rewarded = False
            else:
                failures.append((str(row.get("Session", session_id)),
                                 f"Unknown behavior type: {behavior!r}"))
                continue

            try:
                # Creating configs and  NWB files
                output_path, _ = converters.Initiation_nwb.files_to_config(
                    subject_info=row, output_folder=output_folder
                )
                nwb_file = converters.Initiation_nwb.create_nwb_file_an(config_file=output_path)

                # Add General metadata 
                signal_LFP, regions = converters.acquisition_to_nwb.extract_lfp_signal(csv_data_row=row)

                electrode_region, labels = converters.general_to_nwb.add_general_container(nwb_file=nwb_file, regions=regions)

                # Add acquisition
                converters.acquisition_to_nwb.add_acquisitions(
                    nwb_file=nwb_file,
                    lfp_array=signal_LFP,
                    electrode_region_all=electrode_region,
                    channel_labels=labels,
                )

                # Add Intervals 
                converters.intervals_to_nwb.add_intervals_container(
                    nwb_file=nwb_file, csv_data_row=row, Rewarded=rewarded
                )

                # Add Behavior
                converters.behavior_to_nwb.add_behavior_container(
                    nwb_file=nwb_file, csv_data_row=row, Rewarded=rewarded
                )

                # Save and validate
                subfolder = "Detection Task" if rewarded else "Neutral Exposition"
                output_folder_save = os.path.join(output_folder, subfolder)
                os.makedirs(output_folder_save, exist_ok=True)

                nwb_path = converters.nwb_saving.save_nwb_file(
                    nwb_file=nwb_file, output_folder=output_folder_save
                )

                with NWBHDF5IO(nwb_path, "r") as io:
                    nwb_errors = validate(io=io)

                if nwb_errors:
                    # Delete invalid file
                    try:
                        os.remove(nwb_path)
                    except OSError:
                        pass
                    raise RuntimeError(
                        "NWB validation failed: " + "; ".join(map(str, nwb_errors))
                    )

                # Remove .yaml config file
                if output_path and os.path.exists(output_path):
                    try:
                        os.remove(output_path)
                    except OSError:
                        pass

            except Exception as e:
                failures.append((str(row.get("Session", session_id)), str(e)))
                continue
            finally:
                gc.collect()

        del sessions_df_for_this_pair
        gc.collect()

    # Report missing mice
    missing = [m for m in mouses_name if m not in seen_mice]
    pbar.close()
    if failures or missing:
        if missing:
            print(f"⚠️ No sessions found for: {missing}")
        if failures:
            print("⚠️ Conversion errors:")
            for ident, err in failures:
                print(f"    - {ident}: {err}")

    # Final cleanup: remove orphan .yaml files in root output folder
    for f in Path(output_folder).glob("*.yaml"):
        try:
            f.unlink()
        except OSError:
            pass

    print("**************************************************************************")



#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_ MAIN _-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert data to NWB format for PL sessions")

    parser.add_argument("output_folder", type=str,
                        help="Path to the folder where the NWB file will be saved")

    parser.add_argument("Folder_sessions_info", type=str,
                        help="Path to the folder containing session information")

    parser.add_argument("Folder_general_info", type=str,
                        help="Path to the folder containing general information")

    parser.add_argument("--mouses_name", nargs='+', default=None,
                        help="Mouse name(s) to process (e.g., PL200 PL201)")

    args = parser.parse_args()

    convert_data_to_nwb_pl(
        output_folder=args.output_folder,
        Folder_sessions_info=args.Folder_sessions_info,
        Folder_general_info=args.Folder_general_info,
        mouses_name=args.mouses_name
    )
