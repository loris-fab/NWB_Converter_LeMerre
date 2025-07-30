"""_summary_
"""
import importlib
import os
import h5py
import pandas as pd
import numpy as np
from pynwb import NWBHDF5IO, validate
import argparse
import converters.behavior_to_nwb
import converters.nwb_saving
import converters.general_to_nwb
import converters.Initiation_nwb
import converters.acquisition_to_nwb
import converters.analysis_to_nwb
import converters.intervals_to_nwb


############################################################
# Functions for converting data to NWB format for AN sessions
#############################################################


def convert_data_to_nwb_pl(csv_file, 
                           output_folder, 
                           mouses_name = None):

    csv_data = pd.read_csv(csv_file, sep=";")
    csv_data.columns = csv_data.columns.str.strip() 

    if mouses_name is None:
        mouses_name = csv_data["Mouse Name"].unique().tolist()
    else:
        missing = [name for name in mouses_name if name not in csv_data["Mouse Name"].unique()]
        if missing:
            raise ValueError(f"Mouse name(s) not found in csv file: {missing}")

    
    csv_data = csv_data[csv_data["Mouse Name"].isin(mouses_name)]
    all_sessions = csv_data["Session"]
    
    print("**************************************************************************")
    print("-_-_-_-_-_-_-_-_-_-_-_-_-_-_- NWB conversion _-_-_-_-_-_-_-_-_-_-_-_-_-_-_")
    print("Converting data to NWB format for mouse:", list(all_sessions))
    i = 0
    for index, csv_data_row in csv_data.iterrows():
        i += 1
        print("📃 Creating configs for NWB conversion :") if i == 1 else None
        importlib.reload(converters.Initiation_nwb)
        output_path, config_file = converters.Initiation_nwb.files_to_config(csv_data_row=csv_data_row, output_folder=output_folder)

        print("📑 Created NWB files :") if i == 1 else None
        importlib.reload(converters.general_to_nwb)
        nwb_file = converters.Initiation_nwb.create_nwb_file_an(config_file=output_path) 
                                  
        print("     o 📌 Add general metadata") if i == 1 else None
        importlib.reload(converters.acquisition_to_nwb)
        signal_LFP, regions, EMG, EEG = converters.acquisition_to_nwb.extract_lfp_signal(csv_data_row=csv_data_row)
        electrode_table_region, unique_values = converters.general_to_nwb.add_general_container(nwb_file=nwb_file, csv_data_row=csv_data_row, regions=regions)
        print("         - Subject metadata")
        print("         - Session metadata")
        print("         - Device metadata")
        print("         - Extracellular electrophysiology metadata")

     
        print("     o 📶 Add acquisition container") if i == 1 else None
        #converters.acquisition_to_nwb.add_acquisitions_3series(nwb_file, lfp_array=signal_LFP, electrode_region_all=electrode_table_region, channel_labels=unique_values, emg=EMG, eeg=EEG)


        if False:  
            print("     o ⏸️ Add intervall container")
            importlib.reload(converters.intervals_to_nwb)
            if Rewarded:
                converters.intervals_to_nwb.add_intervals_container_Rewarded(nwb_file=nwb_file, data=data, mat_file=mat_file)
            #else:
                #converters.intervals_to_nwb.add_intervals_container_NonRewarded(nwb_file=nwb_file, data=data, mat_file=mat_file)

            print("     o ⚙️ Add processing container")
            importlib.reload(converters.behavior_to_nwb)
            importlib.reload(converters.analysis_to_nwb)
            if Rewarded:
                print("         - Behavior data")
                converters.behavior_to_nwb.add_behavior_container_Rewarded(nwb_file=nwb_file, data=data, config=config_file)
            else:
                print("         - Behavior data")
                converters.behavior_to_nwb.add_behavior_container_NonRewarded(nwb_file=nwb_file, data=data, config_file=config_file)

            print("         - No ephys data for AN sessions")
            print("         - Analysis complementary information")
            converters.analysis_to_nwb.add_analysis_container(nwb_file=nwb_file, Rewarded=Rewarded, psth_window=psth_window, psth_bin=psth_bin) # almost same for rewarded and non-rewarded sessions
        else :
            importlib.reload(converters.nwb_saving)
            nwb_path = converters.nwb_saving.save_nwb_file(nwb_file=nwb_file, output_folder=output_folder) # same for rewarded and non-rewarded sessions

            print(" ") if i == 1 else None
            print(f"🔎 Validating NWB file before saving...") if i == 1 else None
            with NWBHDF5IO(nwb_path, 'r') as io:
                errors = validate(io=io)

            if not errors:
                print(f"     o ✅ File {nwb_path} is valid, no errors detected and saved successfully.")
            else:
                print(f"     o ❌ NWB file {nwb_path} is invalid, deleting file..")
                os.remove(nwb_path)
                for err in errors:
                    print("         -", err)

            # Delete .yaml config file 
            if os.path.exists(output_path):
                os.remove(output_path)
            
            # Stop after processing 2 sessions for testing purposes
            if i == 1:
                break
    print("No forget to delete the testing purpose")
    print("**************************************************************************")

"""
#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_ MAIN _-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert data to NWB format for AN sessions")
    parser.add_argument("mat_file", type=str, help="Path to the .mat file containing the data")
    parser.add_argument("output_folder", type=str, help="Path to the folder where the NWB file will be saved")
    parser.add_argument("--psth_window", type=float, nargs=2, default=(-0.2, 0.5), help="PSTH and LFP window in seconds for analysis")
    parser.add_argument("--psth_bin", type=float, default=0.010, help="PSTH and LFP bin size in seconds for analysis")

    args = parser.parse_args()

    convert_data_to_nwb_an(
        mat_file=args.mat_file,
        output_folder=args.output_folder,
        psth_window=tuple(args.psth_window),
        psth_bin=args.psth_bin
    )
    
"""