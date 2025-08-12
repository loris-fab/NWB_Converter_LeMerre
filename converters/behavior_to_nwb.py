import numpy as np
from pynwb.base import TimeSeries
from pynwb.behavior import BehavioralEvents, BehavioralTimeSeries
from pynwb import TimeSeries
import pandas as pd


################################################################
# Functions for adding behavior container to NWB file
################################################################

def add_behavior_container_Rewarded(nwb_file,csv_data_row):
    """
    Adds a 'behavior' container to the NWB file from the loaded .mat data.

    :param nwb_file: existing NWB file
    :param csv_data_row: a single row from the CSV file containing behavior data
    :return: None
    """
    response_window = 2


    # 1. Created behavior processing module
    bhv_module = nwb_file.create_processing_module('behavior', 'contains behavioral processed data')

    ###############################################
    ### Add behavioral events                    ###
    ###############################################


    behavior_events = BehavioralEvents(name='BehavioralEvents')
    bhv_module.add_data_interface(behavior_events)


    # --- TRIAL ONSETS ---
    trial_onsets = list(map(float, csv_data_row["Trial_onset"].split(";")))
    ts_trial = TimeSeries(
        name='TrialOnsets',
        data=np.ones_like(trial_onsets),
        unit='n.a.',
        timestamps=trial_onsets,
        description='Timestamps marking the onset of each trial.',
        comments='time start of each trial',
        rate = None,
    )
    behavior_events.add_timeseries(ts_trial)


    # --- STIMULATION FLAGS (stim et flag) ---    
    stim_tms = list(map(float, csv_data_row["stim_onset"].split(";")))  
    stim_data = list()
    for stim in trial_onsets:
        if stim in stim_tms:
            stim_data.append(1)
        elif stim not in stim_tms:
            stim_data.append(0)

    ts_stim_flags = TimeSeries(
        name='StimFlags',
        data=stim_data,
        timestamps=trial_onsets,
        unit='n.a.',
        description='Timestamps marking the whisker stimulation for each trial',
        comments='Whisker stimulation : 0 = no stimulus (Catch trial), 1 = deflection of the C2 whisker (stim trial).',
        rate = None,
    )
    behavior_events.add_timeseries(ts_stim_flags)
    
    
    # --- Responses_times ---
    reaction_times = np.array(list(map(float, csv_data_row["Responses_times"].split(";"))))
    Responses_tms_per_trial = []
    Responses_tms_per_trial_tms = []
    response_data = []
    for i, t0 in enumerate(trial_onsets):
        t1 = t0 + response_window
        indices = np.where((reaction_times >= t0) & (reaction_times < t1))[0]
        if len(indices) > 0:
            Responses_tms_per_trial.append(1)
            Responses_tms_per_trial_tms.append(reaction_times[indices[0]])
            if stim_data[i] == 1:
                response_data.append(1)
            else:
                response_data.append(3)
        else:
            Responses_tms_per_trial.append(0)
            Responses_tms_per_trial_tms.append(t0) # Change to np.nan if you want to indicate no response
            if stim_data[i] == 1:
                response_data.append(0)
            else:
                response_data.append(2)
    if (len(Responses_tms_per_trial) != len(Responses_tms_per_trial_tms)) and  (len(Responses_tms_per_trial) != len(trial_onsets)):
        raise ValueError("Mismatch between number of trials and response timestamps.")

    ts_reaction = TimeSeries(
        name='ReactionTimes',
        data=np.array(Responses_tms_per_trial),
        timestamps=np.array(Responses_tms_per_trial_tms),
        unit='n.a.',
        description = "Timestamps of reaction events defined as a lick occurring after trial onset.",
        comments = "Encoded as 1 at time of reaction, 0 if no reaction occurred with the corresponding trial timestamp.",
    )
    behavior_events.add_timeseries(ts_reaction)


    # ---- "ResponseType" ------

    response_labels_ts = TimeSeries(
        name='ResponseType',
        data=response_data,
        unit='code',
        timestamps=np.array(Responses_tms_per_trial_tms),
        description = "Response type for each trial",
        comments='Integer-encoded trial responses: 0 = MISS, 1 = HIT, 2 = CR (Correct Rejection). We havent defined FA (False Alarm) in this task, but it could be added as 3 if needed.',

    )

    behavior_events.add_timeseries(response_labels_ts)



    #########################################################
    ### Add continuous traces  ###
    #########################################################
    bts = bhv_module.data_interfaces.get('BehavioralTimeSeries')
    if bts is None:
        bts = BehavioralTimeSeries(name='BehavioralTimeSeries')
        bhv_module.add(bts)

    if pd.notna(csv_data_row["EMG"]):
        EMG = list(map(float, csv_data_row["EMG"].split(";")))
    else:
        EMG = None

    # ---------- EMG ----------
    RATE = 2000.0
    UNIT = "V"

    if EMG is not None :
        es_emg = TimeSeries(
            name="ElectricalSeries_EMG",
            data=EMG,
            starting_time=0.0,
            rate=RATE,
            unit=UNIT,
            description="EMG recorded differentially from 2 electrodes, resulting in a single EMG signal",
            comments = "2000 Hz, in V."
        )
        bts.add_timeseries(es_emg)

    trial_onsets = list(map(float, csv_data_row["Trial_onset"].split(";")))
    stim_tms = list(map(float, csv_data_row["stim_onset"].split(";")))  
    stim_data = list()
    for stim in trial_onsets:
        if stim in stim_tms:
            stim_data.append(1)
        elif stim not in stim_tms:
            stim_data.append(0)
    

    return trial_onsets, stim_data , response_data, response_window



def add_behavior_container_NonRewarded(nwb_file,csv_data_row):
    """
    Adds a 'behavior' container to the NWB file from the loaded .mat data.

    :param nwb_file: existing NWB file
    :param csv_data_row: a single row from the CSV file containing behavior data
    :return: None
    """

    # 1. Created behavior processing module
    bhv_module = nwb_file.create_processing_module('behavior', 'contains behavioral processed data')

    ###############################################
    ### Add behavioral events                    ###
    ###############################################


    behavior_events = BehavioralEvents(name='BehavioralEvents')
    bhv_module.add_data_interface(behavior_events)

    
    # --- TRIAL ONSETS ---
    trial_onsets = list(map(float, csv_data_row["Trial_onset"].split(";")))
    ts_trial = TimeSeries(
        name='TrialOnsets',
        data=np.ones_like(trial_onsets),
        unit='n.a.',
        timestamps=trial_onsets,
        description='Timestamps marking the onset of each trial.',
        comments='Encoded as 1 at each trial onset timestamp & the trial duration is 1 seconds.',
        rate = None,
    )
    behavior_events.add_timeseries(ts_trial)
    

    # --- STIMULATION ---    
    stim_tms = list(map(float, csv_data_row["stim_onset"].split(";")))  

    ts_stim_flags = TimeSeries(
        name='StimFlags',
        data=np.ones_like(stim_tms),
        timestamps=stim_tms,
        unit='n.a.',
        description='Timestamps marking the whisker stimulation for each trial',
        comments='Whisker stimulation :1 = deflection of the C2 whisker (stim trial).',
        rate = None,
    )
    behavior_events.add_timeseries(ts_stim_flags)
    
    
    # --- Valve_onset ---
    reaction_times = list(map(float, csv_data_row["Responses_times"].split(";"))) 

    ts_reaction = TimeSeries(
        name='ValveOnsets',
        data=np.ones_like(reaction_times),
        timestamps=reaction_times,
        unit='n.a.',
        description = "Timestamps marking the onset of the valve activation.",
        comments = "Encoded as 1 at each valve activation timestamp. The whisker stimulus was not correlated to the delivery of the reward, therefore, no association between the stimulus and the delivery of the reward could be made.",
    )
    behavior_events.add_timeseries(ts_reaction)

    #########################################################
    ### Add continuous traces  ###
    #########################################################
    bts = bhv_module.data_interfaces.get('BehavioralTimeSeries')
    if bts is None:
        bts = BehavioralTimeSeries(name='BehavioralTimeSeries')
        bhv_module.add(bts)

    if pd.notna(csv_data_row["EMG"]):
        EMG = list(map(float, csv_data_row["EMG"].split(";")))
    else:
        EMG = None

    # ---------- EMG ----------
    RATE = 2000.0
    UNIT = "V"

    if EMG is not None :
        es_emg = TimeSeries(
            name="ElectricalSeries_EMG",
            data=EMG,
            starting_time=0.0,
            rate=RATE,
            unit=UNIT,
            description="EMG recorded differentially from 2 electrodes, resulting in a single EMG signal",
            comments = "2000 Hz, in V."
        )
        bts.add_timeseries(es_emg)

    return None
