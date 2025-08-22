import numpy as np
from pynwb.base import TimeSeries
from pynwb.behavior import BehavioralEvents, BehavioralTimeSeries
from pynwb import TimeSeries
import pandas as pd


################################################################
# Functions for adding behavior container to NWB file
################################################################

def add_behavior_container(nwb_file,csv_data_row,Rewarded):
    """
    Adds a 'behavior' container to the NWB file from the loaded .mat data.

   Args:
       nwb_file (pynwb.file.NWBFile): Target NWB file to which behavior data is added.
       csv_data_row (pandas.Series | Mapping): Row containing behavior data fields.
       Rewarded (bool): if the mouse has a rewarded task.
       
    return: None
    """
    # --- Extract behavior data ---
    trial_onsets = list(map(float, csv_data_row["Trial_onset"].split(";")))
    stim_amp = np.asarray(list(map(float, csv_data_row["stim_amp"].split(";"))))
    stim_onset = np.asarray(list(map(float, csv_data_row["stim_onset"].split(";"))))
    response_data_type = np.asarray(list(map(float, csv_data_row["response_data"].split(";"))))
    lick_time = np.asarray(list(map(float, csv_data_row["lick_time"].split(";"))))
    PiezoLickSignal = np.asarray(list(map(float, csv_data_row["PiezoLickSignal"].split(";"))))
    if not Rewarded:
        reward_onset = np.asarray(list(map(float, csv_data_row["reward_onset"].split(";"))))

    # 1. Created behavior processing module
    bhv_module = nwb_file.create_processing_module('behavior', 'contains behavioral processed data')

    ###############################################
    ### Add behavioral events                    ###
    ###############################################


    behavior_events = BehavioralEvents(name='BehavioralEvents')
    bhv_module.add_data_interface(behavior_events)


    # --- TRIAL ONSETS ---
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

    
    # --- ReactionTimes  ---
    lick_time1 = lick_time[lick_time > 0]
    ts_reaction = TimeSeries(
        name='ReactionTimes',
        data=np.ones_like(lick_time1),
        unit='n.a.',
        timestamps=lick_time1,
        description='Timestamps of response-time defined as lick-onset occurring after trial onset.',
        comments='reaction time from PiezoLickSignal',
        rate = None,
    )
    behavior_events.add_timeseries(ts_reaction)


    # --- STIMULATION FLAGS (stim et flag) ---    
    ts_stim_flags = TimeSeries(
        name='StimFlags',
        data=stim_amp,
        timestamps=trial_onsets,
        unit='code',
        description='Timestamps marking the amplitude of whisker stimulation for each trial',
        comments='Whisker stimulation amplitudes are encoded as integers: 0 = no stimulus (Catch trial), 1 = deflection of the C2 whisker, and higher values indicate increasing stimulation amplitudes.',
        rate = None,
    )
    behavior_events.add_timeseries(ts_stim_flags)
    

    # ---- "ResponseType" ------

    response_labels_ts = TimeSeries(
        name='ResponseType',
        data=response_data_type,
        unit='code',
        timestamps=trial_onsets,
        description = "Response type for each trial",
        comments='trial responses: 0 = MISS, 1 = HIT, 2 = CR (Correct Rejection), 3 = FA (False Alarm), 4 = Unlabeled (no assigned response).',

    )

    behavior_events.add_timeseries(response_labels_ts)

    # ---- "Whisker_hit_trial" ------
    ts_whisker_hit = TimeSeries(
        name='whisker_hit_trial',
        data=(response_data_type == 1).astype(int), 
        unit='n.a.',
        timestamps=trial_onsets,
        description='Timestamps for whisker_hit_trial',
        comments='time of each whisker_hit_trial event.',
        rate=None,
    )
    behavior_events.add_timeseries(ts_whisker_hit)

    # --- whisker_miss_trial ----
    ts_whisker_miss = TimeSeries(
        name='whisker_miss_trial',
        data=(response_data_type == 0).astype(int), 
        unit='n.a.',
        timestamps=trial_onsets,
        description='Timestamps for whisker_miss_trial',
        comments='time of each whisker_miss_trial event.',
        rate=None,
    )
    behavior_events.add_timeseries(ts_whisker_miss)

    # ---- correct_rejection_trial ----
    ts_correct_rejection = TimeSeries(
        name='correct_rejection_trial',
        data=(response_data_type == 2).astype(int),  
        unit='n.a.',
        timestamps=trial_onsets,
        description='Timestamps for correct_rejection_trial',
        comments='time of each correct_rejection_trial event.',
        rate=None,
    )
    behavior_events.add_timeseries(ts_correct_rejection)

    # ---- false_alarm_trial ----
    ts_false_alarm = TimeSeries(
        name='false_alarm_trial',
        data=(response_data_type == 3).astype(int),  
        unit='n.a.',
        timestamps=trial_onsets,
        description='Timestamps for false_alarm_trial',
        comments='time of each false_alarm_trial event.',
        rate=None,
    )
    behavior_events.add_timeseries(ts_false_alarm)

    if not Rewarded:
        # --- reward_onset ---
        ts_reward_onset = TimeSeries(
            name='reward_onset',
            data=np.ones_like(reward_onset),
            timestamps=reward_onset,
            unit='n.a.',
            description = "Timestamps for reward-times",
            comments = "time of each reward delivery event.",
        )
        behavior_events.add_timeseries(ts_reward_onset)
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

    es_PiezoLickSignal = TimeSeries(
    name="ElectricalSeries_PiezoLickSignal",
    data=PiezoLickSignal,
    starting_time=0.0,
    rate=RATE,
    unit=UNIT,
    description="Lick signal over time (V, Sampling rate = 2000 Hz)",
    comments="PiezoLickSignal is the continuous electrical signal recorded from the piezo film attached to the water spout to detect when the mouse contacts the water spout with its tongue."
    )
    bts.add_timeseries(es_PiezoLickSignal)

    return None

