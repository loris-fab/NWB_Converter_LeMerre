

import numpy as np
import h5py
from pynwb.epoch import TimeIntervals   

###############################################################
# Functions for converting intervals to NWB format for AN sessions
###############################################################

def add_intervals_container(nwb_file,csv_data_row,Rewarded) -> None:
    """
    Populate `nwb_file.trials` from one CSV-like row by parsing trial-level fields
    and adding one NWB trial per entry.

    Args:
        nwb_file (pynwb.file.NWBFile): Target NWB file to which trials are added.
        csv_data_row (pandas.Series | Mapping): Row containing at least the keys:
            "Trial_onset", "stim_indices", "stim_amp", "stim_onset",
            "lickflag", "response_data", "lick_time".
            Values are expected as semicolon-separated strings.
        Rewarded (bool): If True, sets `reward_available=1` for all trials (else 0).

    Returns:
        None

    """

    # --- Extract trial data ---
    trial_onsets = list(map(float, csv_data_row["Trial_onset"].split(";")))
    stim_indices = np.asarray(list(map(float, csv_data_row["stim_indices"].split(";"))))
    stim_amp = list(map(float, csv_data_row["stim_amp"].split(";")))
    stim_onset = list(map(float, csv_data_row["stim_onset"].split(";")))
    catch_onset = np.where(stim_indices == 0, trial_onsets, np.nan)
    lick_flag = list(map(int, csv_data_row["lickflag"].split(";")))
    response_data_type = list(map(float, csv_data_row["response_data"].split(";")))
    response_window = 2
    n_trials = len(trial_onsets)
    lick_time = list(map(float, csv_data_row["lick_time"].split(";")))

    # --- Define new trial columns ---
    new_columns = {
        'trial_type': 'stimulus Whisker vs no stimulus trial',
        'whisker_stim': '1 if whisker stimulus delivered, else 0',
        'perf': '0= whisker miss; 1= whisker hit ; 2= correct rejection ; 3= false alarm',
        "whisker_stim_amplitude": "Amplitude of the whisker stimulation between 0 and 5",
        "whisker_stim_duration": "Duration of the whisker stimulation (ms)",
        "whisker_stim_time": "trial start time for whisker_stim=1 else NaN",
        "no_stim" : "1 if no stimulus delivered, else 0",
        "no_stim_time": "trial start time for no_stim=1 (catch trial) else NaN",
        "reward_available": "1 if reward is available, else 0",
        "response_window_start_time": "Start of response window",
        "response_window_stop_time": "Stop of response window",
        "lick_flag": "1 if lick detected, else 0",
        "lick_time": "Within response window lick time. Absolute time (s) relative to session start time"
    }

    # --- Add columns before inserting trials ---
    if nwb_file.trials is None:
        # This creates an empty trial table
        for col, desc in new_columns.items():
            nwb_file.add_trial_column(name=col, description=desc)

    else:
        # Add only missing columns if table already exists
        for col, desc in new_columns.items():
            if col not in nwb_file.trials.colnames:
                nwb_file.add_trial_column(name=col, description=desc)

    # --- Add trials ---
    for i in range(n_trials):
        nwb_file.add_trial(
            start_time=float(trial_onsets[i]),
            stop_time=float(trial_onsets[i]) + response_window,
            trial_type='whisker_trial' if stim_indices[i] else 'no_whisker_trial',
            whisker_stim=stim_indices[i],
            perf=response_data_type[i],
            whisker_stim_time=stim_onset[i],
            whisker_stim_amplitude=stim_amp[i],
            whisker_stim_duration=str("1 (ms)"),
            no_stim= 1 if stim_indices[i] == 0 else 0,
            no_stim_time=catch_onset[i],
            reward_available=1 if Rewarded else 0,
            response_window_start_time=float(trial_onsets[i]) + 0.05,
            response_window_stop_time=float(trial_onsets[i]) + 1,
            lick_flag=lick_flag[i],
            lick_time=lick_time[i]
        )

