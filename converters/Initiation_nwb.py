from datetime import datetime, timedelta
import pandas as pd
import numpy as np
import os
import yaml
from dateutil.tz import tzlocal
from pynwb import NWBFile
from pynwb.file import Subject
from scipy.io import loadmat
import re
from pathlib import Path
from typing import List, Tuple, Optional, Sequence
import mat73 as mat73

#############################################################################
# Function that creates the nwb file object using all metadata
#############################################################################

def create_nwb_file_an(config_file):
    """
    Create an NWBFile from a YAML configuration.

    Args:
        config_file (str): Absolute path to a YAML file containing
            "subject_metadata" and "session_metadata" sections.

    Returns:
        pynwb.file.NWBFile: The in-memory NWB file object.
        None: If the YAML cannot be read or required fields are missing.

    Raises:
        ValueError: If NWB file creation fails after parsing the config.
    """

    try:
        with open(config_file, 'r', encoding='utf8') as stream:
            config = yaml.safe_load(stream)
    except:
        print("Issue while reading the config file.")
        return

    subject_data_yaml = config['subject_metadata']
    session_data_yaml = config['session_metadata']

    # Subject info
    keys_kwargs_subject = ['age', 'age__reference', 'description', 'genotype', 'sex', 'species', 'subject_id',
                           'weight', 'date_of_birth', 'strain']
    kwargs_subject = dict()
    for key in keys_kwargs_subject:
        kwargs_subject[key] = subject_data_yaml.get(key)
        if kwargs_subject[key] is not None:
            kwargs_subject[key] = str(kwargs_subject[key])
    if 'date_of_birth' in kwargs_subject and kwargs_subject['date_of_birth'] != "na":
        date_of_birth = datetime.strptime(kwargs_subject['date_of_birth'], '%m/%d/%Y')
        date_of_birth = date_of_birth.replace(tzinfo=tzlocal())
        kwargs_subject['date_of_birth'] = date_of_birth
    else:
        kwargs_subject.pop('date_of_birth', None)
    subject = Subject(**kwargs_subject)

    # Session info
    keys_kwargs_nwb_file = ['session_description', 'identifier', 'session_id', 'session_start_time',
                            'experimenter', 'experiment_description', 'institution', 'keywords',
                            'notes', 'pharmacology', 'protocol', 'related_publications',
                            'source_script', 'source_script_file_name',  'surgery', 'virus',
                            'stimulus_notes', 'slices', 'lab']

    kwargs_nwb_file = dict()
    for key in keys_kwargs_nwb_file:
        kwargs_nwb_file[key] = session_data_yaml.get(key)
        if kwargs_nwb_file[key] is not None:
            if not isinstance(kwargs_nwb_file[key], list):
                kwargs_nwb_file[key] = str(kwargs_nwb_file[key])
    if 'session_description' not in kwargs_nwb_file:
        print('session_description is needed in the config file.')
        return
    if 'identifier' not in kwargs_nwb_file:
        print('identifier is needed in the config file.')
        return
    if 'session_start_time' not in kwargs_nwb_file:
        print('session_start_time is needed in the config file.')
        return
    else:

        if kwargs_nwb_file['session_start_time'][-2:] == '60': # handle non leap seconds
            kwargs_nwb_file['session_start_time'] = kwargs_nwb_file['session_start_time'][0:-2] + '59'
        try:
            session_start_time = datetime.strptime(kwargs_nwb_file['session_start_time'], '%Y%m%d %H%M%S')
        except ValueError as e:
            if "day is out of range for month" in str(e):
                y = int(kwargs_nwb_file['session_start_time'][0:4])
                m = int(kwargs_nwb_file['session_start_time'][4:6])
                d = int(kwargs_nwb_file['session_start_time'][6:8])
                H = int(kwargs_nwb_file['session_start_time'][9:11])
                M = int(kwargs_nwb_file['session_start_time'][11:13])
                S = int(kwargs_nwb_file['session_start_time'][13:15])

                # Construire une date valide en partant du premier jour du mois
                base_date = datetime(y, m, 1, H, M, S) + timedelta(days=d - 1)
                # Puis soustraire un jour
                session_start_time = base_date - timedelta(days=1)
            else:
                raise

        session_start_time = session_start_time.replace(tzinfo=tzlocal())
        kwargs_nwb_file['session_start_time'] = session_start_time
    if 'session_id' not in kwargs_nwb_file:
        kwargs_nwb_file['session_id'] = kwargs_nwb_file['identifier']


    # Create NWB file object
    kwargs_nwb_file['subject'] = subject
    kwargs_nwb_file['file_create_date'] = datetime.now(tzlocal())

    nwb_file = NWBFile(**kwargs_nwb_file)
    if nwb_file is None:
        raise ValueError("❌ NWB file creation failed. Please check the provided metadata in the config file.")

    return nwb_file


#############################################################################
# Function that creates the config file for the NWB conversion
#############################################################################


def files_to_config(subject_info,output_folder="data"):
    """
    Build a session/subject NWB config from one CSV row and save it as YAML.

    Args:
        subject_info (pandas.Series or Mapping): Row with fields such as
            "Mouse Name", "Session", "Session Date (yyymmdd)", "Start Time (hhmmss)",
            "Behavior Type", etc.
        output_folder (str or pathlib.Path): Folder where the YAML file is saved.

    Returns:
        tuple[str, dict]: (output_path to the written YAML file, in-memory config dict).
    """

    ##  Session metadata extraction 

    ### Experiment_description

    experiment_description = {
    'wh_reward': 1 if str(subject_info.get("Behavior Type", "Unknown").strip()) == "Detection Task" else 0,
    'reward_proba': 1 if str(subject_info.get("Behavior Type", "Unknown").strip()) == "Detection Task" else 0,
    'wh_stim_amps': '5',
    'behavioral_type': "Whisker rewarded (WR+)" if str(subject_info.get("Behavior Type", "Unknown").strip()) == "Detection Task" else "Whisker non-rewarded (WR-)",
    "Session_Counter": float(subject_info.get("counter", 0)),
    "session type":"ephys session",
    'licence': str(subject_info.get("licence", "")).strip(),
    "Software and Algorithms": "Labview, Klusta, MATLAB R2015b",
    'Ambient noise': '80 dB',

    }

    # ---------- Subject and session metadata text ----------

    ### Mouse name, Session name, Related Publications, Experimenter, Session_id, identifier, institution, keywords, Session start time
    mouse = subject_info['Mouse Name']
    session_name = subject_info['Session']
    related_publications = 'Le Merre P, Esmaeili V, Charrière E, Galan K, Salin PA, Petersen CCH, Crochet S. Reward-Based Learning Drives Rapid Sensory Signals in Medial Prefrontal Cortex and Dorsal Hippocampus Necessary for Goal-Directed Behavior. Neuron. 2018 Jan 3;97(1):83-91.e5. doi: 10.1016/j.neuron.2017.11.031. Epub 2017 Dec 14. PMID: 29249287; PMCID: PMC5766832.'
    experimenter = "Pierre Le Merre"
    session_id = subject_info["Session"].strip() 
    identifier = session_id + "_" + str(subject_info["Start Time (hhmmss)"])
    keywords = ["neurophysiology", "behaviour", "mouse", "electrophysiology"] 
    session_start_time = str(subject_info["Session Date (yyymmdd)"])+" " + str(subject_info["Start Time (hhmmss)"])

    ### Birth date and age calculation
    if subject_info["Birth date"] != "Unknown":
        birth_date = pd.to_datetime(subject_info["Birth date"], dayfirst=True).strftime('%m/%d/%Y')
    else:
        birth_date = 'na'

    ### Age
    age = subject_info.get("Mouse Age (d)", "na")
    try:
        age = float(age)
    except Exception:
        age = 'na'

    ### Genotype 
    genotype = str(subject_info.get("mutations", "WT")).strip()
    
    ### weight
    weight = subject_info.get("Weight Session", "na")
    try:
        weight = float(weight)
    except Exception:
        weight = 'na'


    # ---------- Behavior-dependent text ----------
    behavior_type = str(subject_info.get("Behavior Type", "Unknown").strip())
    if behavior_type == "Detection Task":
        session_description = "ephys Whisker Rewarded (WR+) mouse: the mouse was trained to lick within 1 s following a whisker stimulus (go trials) but not in the absence of the stimulus (no-go trials). Chronic multisite LFP recordings were performed using high-impedance tungsten electrodes (FHC, 10-12 MOhms, Shaft 0.075 mm, Catalog numb: UEWSCGSELNND)"
        stimulus_notes = "Whisker stimulation was applied to the C2 region to evoke sensory responses."
    elif behavior_type == "Neutral Exposition":
        session_description = "ephys Whisker non-Rewarded (WR-) mouse: the mouse was rewarded by licking randomly (Free-licking) and exposed to brief whisker stimuli that did not predict reward availability. Chronic multisite LFP recordings were performed using high-impedance tungsten electrodes (FHC, 10-12 MOhms, Shaft 0.075 mm, Catalog numb: UEWSCGSELNND)"
        stimulus_notes = "Whisker stimulation was applied to the C2 region to evoke sensory responses."

    else:
        raise ValueError(f"Unknown behavior type: {behavior_type}")
    
    # ---------- Build config dict ----------
    config = {
        'session_metadata': {
            'experiment_description' : experiment_description,
            'experimenter': experimenter,
            'identifier': identifier,
            'institution': "Ecole Polytechnique Federale de Lausanne",
            'keywords': keywords,
            'lab' : "Laboratory of Sensory Processing",
            'notes': 'Combination of chronic multisite LFP recordings from 5 cortical areas with nuchal EMG recording across multiple sessions as mice learned a whisker sensory detection task reported by licking.' if behavior_type == "Detection Task" else "Combination of chronic multisite LFP recordings from 5 cortical areas with nuchal EMG recording across multiple sessions as mice were exposed to a whisker stimulus that did not predict reward availability.",
            'pharmacology': 'na',
            'protocol': 'na',
            'related_publications': related_publications,
            'session_description': session_description,
            'session_id': session_id,
            'session_start_time': session_start_time,
            'slices': "na", 
            'source_script': 'na',
            'source_script_file_name': 'na',
            'stimulus_notes': stimulus_notes,
            'surgery': 'na',
            'virus': 'na',

        },
        'subject_metadata': {
            'age': age,
            'age__reference': 'birth',
            'date_of_birth': birth_date,
            'description': mouse,
            'genotype': genotype,
            'sex': subject_info.get("Sex_bin", "").upper().strip(),
            'species': "Mus musculus",
            'strain': subject_info.get("strain", "").strip(),
            'subject_id': mouse,
            'weight': weight,

        },
    }

    # save config
    output_path = os.path.join(output_folder, f"{session_name}_config.yaml")
    with open(output_path, 'w') as f:
        yaml.dump(config, f, default_flow_style=False)
    return output_path, config


#############################################################################
# Function that creates the row pandas DataFrame for the NWB conversion
#############################################################################

def remove_rows_of_df(df, info_to_remove, col):
    """
    Return a copy of the DataFrame without rows matching a given mouse name.

    Args:
        df (pandas.DataFrame): Input table.
        mouse_name (str): Exact value to remove from the "Mouse Name" column.

    Returns:
        pandas.DataFrame: Filtered DataFrame.
    """
    return df[df[col] != info_to_remove]

def remove_nwb_files(folder_path):
    """
    Delete all .nwb files in a folder (best-effort, non-recursive).

    Args:
        folder_path (str or pathlib.Path): Directory to scan.

    Returns:
        None

    Notes:
        Prints an error message for files that cannot be removed.
    """
    for filename in os.listdir(folder_path):
        if filename.endswith('.nwb'):
            file_path = os.path.join(folder_path, filename)
            try:
                os.remove(file_path)
            except Exception as e:
                print(f"Erreur lors de la suppression de {file_path}: {e}")

def find_pl_pairs(
    path_folder1,
    path_folder2,
    mouse_names = None,
) :
    """
    Match PL mouse folders with the corresponding WR group folder.

    This scans `path_folder1` for subfolders whose names start with "PL2"
    (optionally filtered by `mouse_names`). For each mouse name, it searches
    in `path_folder2/"WR+ mice"` first, then `path_folder2/"WR- mice"`, for
    any file whose name starts with that mouse name. It returns the pair:
    (absolute path to the PL folder, absolute path to the matched WR folder).

    Args:
        path_folder1 (str): Root directory containing PL mouse folders.
        path_folder2 (str): Root directory containing "WR+ mice" / "WR- mice".
        mouse_names (Sequence[str], optional): Subset of mouse names to keep.

    Returns:
        list[tuple[str, str]]: List of (pl_dir, wr_dir) absolute paths.

    Raises:
        FileNotFoundError: If no file starting with the mouse name is found
            in either "WR+ mice" or "WR- mice".
    """

    p1, p2 = Path(path_folder1), Path(path_folder2)
    pairs: List[Tuple[str, str]] = []
    name_filter = set(mouse_names) if mouse_names else None

    for d1 in sorted(p1.iterdir()):
        if not d1.is_dir():
            continue
        name = d1.name 

        # filtre
        if name_filter is not None:
            if name not in name_filter or not name.startswith("PL2"):
                continue
        else:
            if not name.startswith("PL2"):
                continue

        match = None
        # 1) WR+ mice
        wr = p2 / "WR+ mice"
        if wr.is_dir():
            for p in wr.glob(f"{name}*"):
                if p.is_file():
                    match = p
                    break

        # 2) WR- mice si pas trouvé
        if match is None:
            wr = p2 / "WR- mice"
            if wr.is_dir():
                for p in wr.glob(f"{name}*"):
                    if p.is_file():
                        match = p
                        break

        # 3) rien trouvé -> erreur
        if match is None:
            raise FileNotFoundError(
                f"No file starting with '{name}' in '{p2 / 'WR+ mice'}' or '{p2 / 'WR- mice'}'."
            )

        pairs.append((str(d1.resolve()), str(wr.resolve())))
    return pairs

def files_to_dataframe(PL, PLALL, dataframe):
    """
    Append session rows to a DataFrame from PL general info and per-session .mat files.

    It reads one general .mat in `PL` (SciPy loadmat) and iterates over session
    files in `PLALL` (Matlab v7.3 via mat73). It extracts metadata (dates,
    behavior type, stim/catch onsets, lick flags/times, signals) and appends a
    row per session to `dataframe`.

    Args:
        PL (str or pathlib.Path): Folder with the general .mat file.
        PLALL (str or pathlib.Path): Folder with per-session .mat files for that mouse.
        dataframe (pandas.DataFrame): Existing table to extend.

    Returns:
        pandas.DataFrame: The input DataFrame with appended session rows.
    """
    csv_data = dataframe
    csv_data.columns = csv_data.columns.str.strip()

    # 1) Locate and read the general .mat file from PL

    General_data = next((os.path.join(PL, f) for f in os.listdir(PL) 
                         if (f.endswith('.mat') and f.startswith("PL2"))), 
                         None)

    if General_data:
        General_data = loadmat(General_data)
    else:
        raise FileNotFoundError("No .mat file found in the specified directory.")
    
    # 2) Subject-level metadata (from general .mat)
    # Unique indices that define each session in general data (used to align dates/times)
    _, id_session_indices = np.unique(General_data["LFP_Data"][0][0][4], return_index=True)

    # load mouse_name, strain, sex 
    # Cast to str for consistent downstream usage
    mouse_name = str(General_data["LFP_Data"][0][0][0][0][0][0])
    strain = str(General_data["LFP_Data"][0][0][1][0][0][0])
    sex = str(General_data["LFP_Data"][0][0][2][0][0][0])

    # Load birth date
    raw_birth = General_data["LFP_Data"][0][0][3][0][0][0][0]
    birth_date = "Unknown" if raw_birth is None or str(raw_birth).strip() in ["", "[]"] else str(raw_birth)

    # 3) Iterate over session files in PLALL (MAT v7.3 via mat73)
    #    Files are sorted by the trailing D<number> token (e.g., "...D5.mat").

    def extract_number(file_name):
        """
        Extract the integer that follows 'D' in filenames ending with '.mat'.
        Example: 'PL200_xxx_D12.mat' -> 12.
        Returns -1 for non-.mat files (so they sink in sorting).
        Raises ValueError if a .mat file lacks a 'D<digits>' pattern.
        """
        if file_name.endswith('.mat') :
            match = re.search(r'D(\d+)', file_name)
            if match:
                return int(match.group(1))
            raise ValueError(f"No number found in a .mat file: {file_name}")
        return -1
    
    i = -1 # Index `i` aligns session order from general data with per-session files
    
    for  file_name in sorted(os.listdir(PLALL), key=extract_number):
        file_path = os.path.join(PLALL, file_name)
        # Only process good files:
        if os.path.isfile(file_path) and file_name.endswith('.mat') and file_name.startswith(mouse_name):
            pli = mat73.loadmat(file_path)
            i += 1
            # --------------------- Dates (dd.mm.yyyy & yyyymmdd) ---------------------
            # General_data["LFP_Data"][0][0][5][idx] = [yy, mm, dd, hh, mm, ss] (nested)
            ## dd
            dd = str(General_data["LFP_Data"][0][0][5][id_session_indices[i]][2][0][0])
            if len(dd) == 1:
                dd = "0" + dd
            ## mm
            mm = str(General_data["LFP_Data"][0][0][5][id_session_indices[i]][1][0][0])
            if len(mm) == 1:
                mm = "0" + mm
            ## yy
            yy = str(General_data["LFP_Data"][0][0][5][id_session_indices[i]][0][0][0])

            start_date = dd + "." + mm + "." + yy
            start_date_2 = yy + mm + dd
            End_date = start_date

            session = mouse_name + "_" + start_date_2

            # --------------------- Start time (hhmmss) ---------------------
            ## hh
            hh = str(General_data["LFP_Data"][0][0][5][id_session_indices[i]][3][0][0])
            if len(hh) == 1:
                hh = "0" + hh
            ## mm
            mm = str(General_data["LFP_Data"][0][0][5][id_session_indices[i]][4][0][0])
            if len(mm) == 1:
                mm = "0" + mm
            ## ss
            ss = str(General_data["LFP_Data"][0][0][5][id_session_indices[i]][5][0][0])
            if len(str(int(float(ss)))) == 1:
                ss = "0" + str(int(float(ss)))
            else:
                ss = str(int(float(ss)))
                
            start_time = hh + mm +  ss

            # --------------------- Behavior type ---------------------
            behaviortype = str(General_data["LFP_Data"][0][0][6][id_session_indices[i]][0][0])
            if behaviortype == "DT":
                behaviortype = "Detection Task"
            elif behaviortype == "X":
                behaviortype = "Neutral Exposition"
            else:
                raise ValueError(f"Unknown behavior type: {behaviortype}")


            # --------------------- Events & signals (per-session .mat) ---------------------

            # Stim_times
            stim_onset= np.asarray(pli["Performance"]["data"].T[3])

            #reward_onset
            if behaviortype == "Neutral Exposition":
                reward_onset = np.asarray(pli["Valve_times"]["data"]/2000)

            # trial_onset
            Trial_onset = np.asarray(pli["Performance"]["data"].T[0])

            #stim_indice
            stim_indices = np.asarray(pli["Performance"]["data"].T[1]).flatten()
            #stim_amp
            stim_amp = np.asarray(pli["Performance"]["data"].T[2]).flatten()

            # Catch_times
            catch_onset = [np.nan if el in Trial_onset else el for el in stim_onset]
            
            #Response_data
            response_data = np.asarray(pli["Performance"]["data"].T[8]).flatten()

            # Lick flags and lick times
            lickflag = np.asarray(pli["Performance"]["data"].T[6]).flatten().astype(int)
            lick_time = np.asarray(pli["Performance"]["data"].T[7]).flatten()

            # signals: use NaN if key is missing (keep original behavior)
            EMG  = np.asarray(pli["EMG"]["data"]).flatten()  if "EMG"  in pli.keys() else np.nan
            PtA  = np.asarray(pli["PtA"]["data"]).flatten()  if "PtA"  in pli.keys() else np.nan
            dCA1 = np.asarray(pli["dCA1"]["data"]).flatten() if "dCA1" in pli.keys() else np.nan
            mPFC = np.asarray(pli["mPFC"]["data"]).flatten() if "mPFC" in pli.keys() else np.nan
            wM1  = np.asarray(pli["wM1"]["data"]).flatten()  if "wM1"  in pli.keys() else np.nan
            wS1  = np.asarray(pli["wS1"]["data"]).flatten()  if "wS1"  in pli.keys() else np.nan
            wS2  = np.asarray(pli["wS2"]["data"]).flatten()  if "wS2"  in pli.keys() else np.nan
            antM1 = np.asarray(pli["antM1"]["data"]).flatten() if "antM1" in pli.keys() else np.nan
            EEG   = np.asarray(pli["EEG"]["data"]).flatten()   if "EEG"   in pli.keys() else np.nan


            # Create a new row for the session
            new_row = {
                "Mouse Name": mouse_name,
                "User (user_userName)": "PL",
                "Ear tag": "Unknown",
                "Start date (dd.mm.yy)": start_date,
                "End date": End_date,
                "Sex_bin": sex,
                "strain": strain,
                "mutations": "WT",
                "Birth date": birth_date,
                "licence": "DR2013-47",
                "DG": "",
                "ExpEnd": "",
                "counter": float(pli["Session_Info"]["Session_Counter"]),
                "Created on": "Unknown",
                "Session": session,
                "Session Date (yyymmdd)": start_date_2,
                "Start Time (hhmmss)": start_time,
                "Behavior Type": behaviortype,
                "Session Type": "Whisker Rewarded",
                "Mouse Age (d)": "Unknown",
                "Weight of Reference": "Unknown",
                "Weight Session": "Unknown",
                "Trial_onset" : ';'.join(map(str, Trial_onset)),
                "stim_indices": ';'.join(map(str, stim_indices)),
                "stim_onset": ';'.join(map(str, stim_onset)),
                "catch_onset": ';'.join(map(str, catch_onset)),
                "stim_amp": ';'.join(map(str, stim_amp)),
                "lickflag": ';'.join(map(str, lickflag)),
                "lick_time": ';'.join(map(str, lick_time)),
                "reward_onset": ';'.join(map(str, reward_onset)) if behaviortype == "Neutral Exposition" else np.nan,
                "PiezoLickSignal": ';'.join(map(str, np.asarray(pli["lick"]["data"]))),
                "response_data": ';'.join(map(str, response_data)) ,
                "EMG": ';'.join(map(str, EMG)) if not (isinstance(EMG, float) and np.isnan(EMG)) else np.nan,
                "PtA": ';'.join(map(str, PtA)) if not (isinstance(PtA, float) and np.isnan(PtA)) else np.nan,
                "dCA1": ';'.join(map(str, dCA1)) if not (isinstance(dCA1, float) and np.isnan(dCA1)) else np.nan,
                "mPFC": ';'.join(map(str, mPFC)) if not (isinstance(mPFC, float) and np.isnan(mPFC)) else np.nan,
                "wM1": ';'.join(map(str, wM1)) if not (isinstance(wM1, float) and np.isnan(wM1)) else np.nan,
                "wS1": ';'.join(map(str, wS1)) if not (isinstance(wS1, float) and np.isnan(wS1)) else np.nan,
                "wS2": ';'.join(map(str, wS2)) if not (isinstance(wS2, float) and np.isnan(wS2)) else np.nan,
                "antM1": ';'.join(map(str, antM1)) if not (isinstance(antM1, float) and np.isnan(antM1)) else np.nan,
                "EEG": ';'.join(map(str, EEG)) if not (isinstance(EEG, float) and np.isnan(EEG)) else np.nan,
            }

            # Append the new row to the DataFrame
            csv_data = pd.concat([csv_data, pd.DataFrame([new_row])], ignore_index=True)

    return csv_data
