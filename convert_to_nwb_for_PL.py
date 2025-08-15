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
import converters.intervals_to_nwb
from pathlib import Path
import gc

Folder_general_info_server= "/Volumes/Petersen-Lab/publications/2018/2018_LeMerre_Neuron/2018_LeMerre_Neuron_data/processed_data"
Folder_sessions_info_server = "/Volumes/Petersen-Lab/analysis/Sylvain_Crochet/DATA_REPOSITORY/LeMerre_mPFC_2018/Chronic_LFPs_Preprocessed"

############################################################
# Functions for converting data to NWB format for AN sessions
#############################################################


def convert_data_to_nwb_pl(output_folder,Folder_sessions_info = Folder_general_info_server, Folder_general_info = Folder_general_info_server, mouses_name=None, ):

    if mouses_name == None: 
        mouses_name = ["PL200", "PL201", "PL202", "PL203", "PL204", "PL205", "PL206", "PL207", "PL208", "PL209", "PL210", "PL211", "PL212", "PL213","PL214", "PL215", "PL216", "PL217", "PL218", "PL219", "PL220","PL221", "PL222", "PL223", "PL224", "PL225"]

    # Find pairs of processed and raw data files 
    importlib.reload(converters.Initiation_nwb)
    pairs = converters.Initiation_nwb.find_pl_pairs(Folder_general_info, Folder_sessions_info, mouse_names=mouses_name)


    print("**************************************************************************")
    print("-_-_-_-_-_-_-_-_-_-_-_-_-_-_- NWB conversion _-_-_-_-_-_-_-_-_-_-_-_-_-_-_")
    # Loop through each pair of processed and raw data files
    csv_data = pd.DataFrame(columns=['Mouse Name'])
    failures = []  # (mouse_name, err_msg)

    pbar = tqdm(pairs, unit="file")
    
    for PL, PLLA in pbar:
    #try:
        pbar.set_description(f"Loading data {Path(PL).name} ...")
        csv_data= converters.Initiation_nwb.files_to_dataframe(PL, PLLA, csv_data)
        gc.collect()
    #except Exception as e:
    #    failures.append((Path(PL).name, str(e)))

    for mouse, err_msg in failures:
        csv_data = converters.Initiation_nwb.remove_rows_of_df(csv_data, mouse, "Mouse Name")


    csv_data = csv_data[csv_data["Mouse Name"].isin(mouses_name)]

    # Check for missing mouse names
    missing = [name for name in mouses_name if name not in csv_data["Mouse Name"].unique().tolist()]


    print("Converting data to NWB format for mouse:", list(np.unique(csv_data["Mouse Name"])))
    bar = tqdm(total=len(csv_data), desc="Processing ")
    for _, csv_data_row in csv_data.iterrows():
        bar.set_postfix_str(str(csv_data_row["Mouse Name"])) 
        bar.update(1)
    #try:
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
    #except Exception as e:
    #    failures.append((csv_data_row["Session"], str(e)))
    #    continue
    #finally:
        bar.update(1)
    gc.collect()
#if len(failures) > 0:
    #    print(f"⚠️ Conversion completed except for : {missing} because of the following errors:")
    #    for i, (mouse_name, error) in enumerate(failures):
    #        print(f"    - {mouse_name}: {error}")

    bar.close()
    for f in Path(output_folder).glob("*.yaml"):  
        f.unlink()

    print("**************************************************************************")


#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_ MAIN _-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert data to NWB format for AN sessions")
    parser.add_argument("output_folder", type=str, help="Path to the folder where the NWB file will be saved")
    parser.add_argument("--Folder_sessions_info", type=str, default=Folder_sessions_info_server, help="Path to the folder containing session information")
    parser.add_argument("--Folder_general_info", type=str, default=Folder_general_info_server, help="Path to the folder containing general information")
    parser.add_argument("--mouses_name", nargs='+', default=None, help="Mouse name(s) to process (e.g., PL200 PL201)")

    args = parser.parse_args()

    convert_data_to_nwb_pl(
        output_folder=args.output_folder,
        Folder_sessions_info=args.Folder_sessions_info,
        Folder_general_info=args.Folder_general_info,
        mouses_name=args.mouses_name
    )
