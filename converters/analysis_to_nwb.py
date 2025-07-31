# converters/analysis_to_nwb.py
import numpy as np
from pynwb import ProcessingModule
from pynwb.base import TimeSeries
import matplotlib.pyplot as plt

#######################################################
# Functions for converting analysis data to NWB format
#######################################################


def add_analysis_container(
    *,
    nwb_file, 
    psth_window,         # seconds around stimulus
    psth_bin,       # 10-ms bins
    rewarded
):
    """
    Populate ``nwb_file.analysis`` with:
      • a PSTH (all units, stimulus-aligned)
      • a mean stimulus-locked LFP trace per channel

    Parameters
    ----------
    nwb_file : pynwb.NWBFile
    data : dict
        Already-loaded MATLAB 'Data' struct (keys → ndarrays / HDF refs)
    Rewarded : bool
        True if the session is rewarded, False otherwise
    psth_window : tuple
        (start, stop) window around each stimulus for the PSTH (seconds)
    psth_bin : float
        Bin width for PSTH (seconds)
    """

    ####################################################
    ###  Make / get the "analysis" processing module ###
    ####################################################
    
    ana_mod = nwb_file.create_processing_module(
            name="analysis",
            description="Secondary analyses: mean LFP and global LFP")

    start_w, stop_w = psth_window
    if rewarded:
        TrialOnsets_All = np.array(nwb_file.processing['behavior'].data_interfaces['BehavioralEvents'].time_series["TrialOnsets"].timestamps[:])
        stim_indice = nwb_file.processing['behavior'].data_interfaces['BehavioralEvents'].time_series["StimFlags"].data[:]
        stim_times = TrialOnsets_All[np.array(stim_indice) >= 1]
    else:
        stim_times = nwb_file.processing['behavior'].data_interfaces['BehavioralEvents'].time_series["StimFlags"].timestamps[:]


    ###################
    ###  Mean LFP   ###
    ###################

    # fetch raw LFP acquisition
    lfp_acq = nwb_file.acquisition["ElectricalSeries_LFP"]
    lfp_rate = lfp_acq.rate
    lfp_data = lfp_acq.data 

    # determine sample indices for the window
    idx_pre  = int(np.floor(start_w * lfp_rate))
    idx_post = int(np.ceil(stop_w  * lfp_rate))
    span = idx_post - idx_pre

    mean_lfp = np.zeros((lfp_data.shape[1], span), dtype=np.float32)

    for s in stim_times:
        i_start = int(np.round((s - lfp_acq.starting_time) * lfp_rate)) + idx_pre
        i_end   = i_start + span
        # skip if window exceeds bounds
        if i_start < 0 or i_end > lfp_data.shape[0]:
            continue
        mean_lfp += lfp_data[i_start:i_end, :].T   # channels × time

    mean_lfp /= len(stim_times)

    lfp_times = np.linspace(start_w, stop_w, span, endpoint=False)

    mean_lfp_ts = TimeSeries(
        name="MeanLFP_all_electrodes",
        data=mean_lfp.T,
        unit="volts",
        timestamps=lfp_times,
        description="Stimulus-aligned average LFP (channels × time)",
        comments="Averaged across all trials; same window and alignment as PSTH."
    )
    ana_mod.add_data_interface(mean_lfp_ts)

    # Compute mean across electrodes (1D signal : time only)
    mean_lfp_all_channels = mean_lfp.mean(axis=0)  # shape: (time,)

    # Create new TimeSeries for the global mean LFP
    mean_lfp_global_ts = TimeSeries(
        name="MeanLFP_global",
        data=mean_lfp_all_channels,
        unit="volts",
        timestamps=lfp_times,
        description="Average LFP across all electrodes, aligned to stimulus",
        comments="Mean of MeanLFP_all_electrodes across channels and averaged across trials"
    )

    # Add it to the analysis module
    ana_mod.add_data_interface(mean_lfp_global_ts)


    # plot the mean LFP
    plt.figure(figsize=(8, 4))
    for ch in range(min(1, mean_lfp.shape[0])):
        plt.plot(lfp_times, mean_lfp[ch, :], label=f"Channel {ch+1}")
    plt.axvline(x=0, color='red', linestyle='--', label="Stimulus")
    plt.xlabel("Time (s)")
    plt.ylabel("LFP (V)")
    plt.title("Mean LFP - First Channel")
    plt.legend()
    plt.tight_layout()
    #plt.show()
    #plt.savefig("data/analysis/mean_lfp_example.png")
    #plt.close()

    # plot the global mean LFP
    plt.figure(figsize=(8, 4))
    plt.plot(lfp_times, mean_lfp_all_channels, label="Global Mean LFP")
    plt.axvline(x=0, color='red', linestyle='--', label="Stimulus")
    plt.xlabel("Time (s)")
    plt.ylabel("LFP (V)")
    plt.title("Global Mean LFP")
    plt.legend()
    plt.tight_layout()
    #plt.show()
    #plt.savefig("data/analysis/mean_lfp_global.png")
    #plt.close()


    return None