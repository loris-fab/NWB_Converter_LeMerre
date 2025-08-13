"""_summary_
"""
import importlib
import os
import h5py
import pandas as pd
import numpy as np
from tqdm import tqdm
from pynwb import NWBHDF5IO, validate
import argparse
import converters.behavior_to_nwb
import converters.nwb_saving
import converters.general_to_nwb
import converters.Initiation_nwb
import converters.acquisition_to_nwb
import converters.analysis_to_nwb
import converters.intervals_to_nwb
from pathlib import Path
import gc

############################################################
# Functions for converting data to NWB format for AN sessions
#############################################################


def convert_data_to_nwb_pl(output_folder,Folder_sessions_info, Folder_general_info = "/Volumes/Petersen-Lab/publications/2018/2018_LeMerre_Neuron/2018_LeMerre_Neuron_data/processed_data",mouses_name=None, ):


    # Find pairs of processed and raw data files 
    importlib.reload(converters.Initiation_nwb)
    pairs = converters.Initiation_nwb.find_pl_pairs(Folder_general_info, Folder_sessions_info, mouse_names=mouses_name)
    #print("pairs:", pairs)
    csv_data = pd.DataFrame(columns=['Mouse Name'])
    csv_data.columns = csv_data.columns.str.strip() 

    print("**************************************************************************")
    print("-_-_-_-_-_-_-_-_-_-_-_-_-_-_- NWB conversion _-_-_-_-_-_-_-_-_-_-_-_-_-_-_")
    # Loop through each pair of processed and raw data files
    for pair in tqdm(pairs, desc="Loading data ..."):
        PL, PLLA = pair
        csv_data = converters.Initiation_nwb.files_to_dataframe(PL, PLLA, csv_data)
    gc.collect()
    #display("Dataframe after loading:", csv_data)

    if mouses_name is None:
        mouses_name = csv_data["Mouse Name"].unique().tolist()
    else:
        missing = [name for name in mouses_name if name not in csv_data["Mouse Name"].unique()]
        if missing:
            raise ValueError(f"Mouse name(s) not found in csv file: {missing}")

    
    csv_data = csv_data[csv_data["Mouse Name"].isin(mouses_name)]
    all_sessions = csv_data["Session"]
    

    print("Converting data to NWB format for mouse:", list(np.unique(csv_data["Mouse Name"])))
    failures = []  # (mouse_name, err_msg)
    bar = tqdm(total=len(csv_data), desc="Processing ")
    for _, csv_data_row in csv_data.iterrows():
        bar.set_postfix_str(str(csv_data_row["Mouse Name"])) 
        bar.update(1)
        try:
            if csv_data_row["Behavior Type"] == "Detection Task":
                Rewarded = True
            elif csv_data_row["Behavior Type"] == "Neutral Exposition":
                Rewarded = False
            else :
                raise ValueError(f"Unknown behavior type: {csv_data_row['Behavior Type']}")

            # Creating configs for NWB conversion
            importlib.reload(converters.Initiation_nwb)
            output_path, config_file = converters.Initiation_nwb.files_to_config(csv_data_row=csv_data_row, output_folder=output_folder)  #same between Rewarded and NonRewarded sessions

            # 📑 Created NWB files
            importlib.reload(converters.general_to_nwb)
            nwb_file = converters.Initiation_nwb.create_nwb_file_an(config_file=output_path)  #same between Rewarded and NonRewarded sessions

            # o 📌 Add general metadata
            importlib.reload(converters.acquisition_to_nwb)
            signal_LFP, regions= converters.acquisition_to_nwb.extract_lfp_signal(csv_data_row=csv_data_row)
            electrode_table_region, labels = converters.general_to_nwb.add_general_container(nwb_file=nwb_file, regions=regions)  #same between Rewarded and NonRewarded sessions

            # o 📶 Add acquisition container
            converters.acquisition_to_nwb.add_acquisitions_3series(nwb_file=nwb_file, lfp_array=signal_LFP, electrode_region_all=electrode_table_region, channel_labels=labels)  #same between Rewarded and NonRewarded sessions

            # o ⏸️ Add intervall container
            importlib.reload(converters.intervals_to_nwb)
            converters.intervals_to_nwb.add_intervals_container(nwb_file=nwb_file,csv_data_row=csv_data_row, Rewarded=Rewarded)

            # o ⚙️ Add behavior container
            importlib.reload(converters.behavior_to_nwb)
            converters.behavior_to_nwb.add_behavior_container(nwb_file=nwb_file,csv_data_row=csv_data_row, Rewarded=Rewarded)

            # 🔎 Validating NWB file and saving...
            importlib.reload(converters.nwb_saving)
            if Rewarded:
                output_folder_save = os.path.join(output_folder, "Detection Task")
            else:
                output_folder_save = os.path.join(output_folder, "Neutral Exposition")
            os.makedirs(output_folder_save, exist_ok=True)
            nwb_path = converters.nwb_saving.save_nwb_file(nwb_file=nwb_file, output_folder=output_folder_save) #same between Rewarded and NonRewarded sessions

            with NWBHDF5IO(nwb_path, 'r') as io:
                nwb_errors = validate(io=io)

            if nwb_errors:
                os.remove(nwb_path)
                raise RuntimeError("NWB validation failed: " + "; ".join(map(str, nwb_errors)))

            # Delete .yaml config file 
            if os.path.exists(output_path):
                os.remove(output_path)
        except Exception as e:
            failures.append((csv_data_row["Session"], str(e)))
            continue
        finally:
            bar.update(1)
        gc.collect()
    if len(failures) > 0:
            print(f"⚠️ Conversion completed with errors for {len(failures)} files")
            for i, (mouse_name, error) in enumerate(failures):
                print(f"    - {mouse_name}: {error}")

    bar.close()
    for f in Path(output_folder).glob("*.yaml"):  
        f.unlink()

    print("**************************************************************************")


#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_ MAIN _-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert data to NWB format for AN sessions")
    parser.add_argument("csv_file", type=str, help="Path to the .csv file containing the data")
    parser.add_argument("output_folder", type=str, help="Path to the folder where the NWB file will be saved")
    parser.add_argument("--mouses_name", nargs='+', default=None, help="Mouse name(s) to process (e.g., PL200 PL201)")

    args = parser.parse_args()

    convert_data_to_nwb_pl(
        csv_file=args.csv_file,
        output_folder=args.output_folder,
        mouses_name=args.mouses_name
    )
