

import numpy as np
import h5py
from pynwb.epoch import TimeIntervals   

###############################################################
# Functions for converting intervals to NWB format for AN sessions
###############################################################

def add_intervals_container_Rewarded(nwb_file, info_trials) -> None:
    """
    Add detailed trial information to the NWBFile for a rewarded whisker detection task.
    """

    # --- Extract trial data ---
    trial_onsets = np.asarray(info_trials[0])
    stim_indices = np.asarray(info_trials[1]).flatten().astype(bool)
    response_data_type = np.asarray(info_trials[2]).flatten()
    window_trial = float(info_trials[3])
    n_trials = len(trial_onsets)


    # --- Define new trial columns ---
    new_columns = {
        'trial_type': 'stimulus Whisker vs no stimulus trial',
        'whisker_stim': '1 if whisker stimulus delivered, else 0',
        #'response_window_start_time': 'Start of response window',
        'ResponseType': 'Trial outcome label (Hit, Miss, etc.)',
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
            stop_time=float(trial_onsets[i]) + window_trial,
            trial_type='whisker_trial' if stim_indices[i] else 'no_whisker_trial',
            whisker_stim=int(stim_indices[i]),
            #response_window_start_time=float(reaction_abs[i]),
            ResponseType=response_data_type[i],
            #lick_time=lick_time_per_trial[i]
        )

