from pynwb.ecephys import ElectricalSeries
import h5py
import numpy as np
import warnings
import os
import pandas as pd


#####################################################
# Functions that handle LFP acquisition in NWB files
#####################################################

def add_acquisitions_3series(nwb_file, lfp_array, electrode_region_all, channel_labels, emg=None, eeg=None):
    """
    Creates 3 separate ElectricalSeries (LFP, EMG, EEG) if the data/electrodes exist.
    - rate = 2000.0
    - unit = "V"

    Parameters
    ----------
    nwb_file : pynwb.NWBFile
        The NWB file object to which the acquisitions will be added.
    lfp_array : np.ndarray
        Array of shape (T, n_lfp) containing LFP data.
    electrode_region_all : DynamicTableRegion
        DynamicTableRegion containing all created channels.
    channel_labels : list[str]
        List of channel labels (returned by add_general_container).
    emg : np.ndarray or None, optional
        Array of shape (T,) or (T,2) containing EMG data, or None.
    eeg : np.ndarray or None, optional
        Array of shape (T,) or (T,2) containing EEG data, or None.
    """


    RATE = 2000.0
    UNIT = "V"
    LFP_LABELS = ["PtA", "dCA1", "mPFC", "wM1", "wS1", "wS2", "antM1"]


    def _region_indices(dtr):
        # indices globaux des électrodes
        try:
            return list(dtr.data[:])
        except Exception:
            try:
                print("second")
                return list(dtr.region.data[:])
            except Exception:
                print("third")
                return list(dtr.region)

    # mapping label -> index global d'électrode
    idx_all = _region_indices(electrode_region_all)
    label_to_idx = dict(zip(channel_labels, idx_all))

    # ---------- LFP ----------
    lfp_labels = [lab for lab in channel_labels if lab in LFP_LABELS][:lfp_array.shape[1]]
    if lfp_array.shape[1] != len(lfp_labels):
        raise ValueError(f"Inconsistency in LFP: data has {lfp_array.shape[1]} columns, found {len(lfp_labels)} LFP labels.")
    
    lfp_idx = [label_to_idx[lab] for lab in lfp_labels]
    lfp_region = nwb_file.create_electrode_table_region(lfp_idx, "LFP electrodes")
    es_lfp = ElectricalSeries(
        name="ElectricalSeries_LFP",
        data=lfp_array,
        electrodes=lfp_region,
        starting_time=0.0,
        rate=RATE,
        description="LFPs recorded from multiple electrodes",
        comments = "2000 Hz, in V."
    )
    nwb_file.add_acquisition(es_lfp)

    # ---------- EMG ----------
    es_emg = None
    emg_labels = [lab for lab in channel_labels if lab.startswith("EMG")]
    if emg is not None and len(emg_labels) > 0:
        emg_idx = [label_to_idx[lab] for lab in emg_labels]
        emg_region = nwb_file.create_electrode_table_region(emg_idx, "EMG electrodes")
        es_emg = ElectricalSeries(
            name="ElectricalSeries_EMG",
            data=emg,
            electrodes=emg_region,
            starting_time=0.0,
            rate=RATE,
            description="EMG recorded differentially from 2 electrodes, resulting in a single EMG signal",
            comments = "2000 Hz, in V."
        )
        nwb_file.add_acquisition(es_emg)

    # ---------- EEG  ----------
    es_eeg = None
    eeg_labels = [lab for lab in channel_labels if lab.startswith("EEG")]
    if eeg is not None and len(eeg_labels) > 0:
        eeg_idx = [label_to_idx[lab] for lab in eeg_labels]
        eeg_region = nwb_file.create_electrode_table_region(eeg_idx, "EEG electrodes")
        es_eeg = ElectricalSeries(
            name="ElectricalSeries_EEG",
            data=eeg,
            electrodes=eeg_region,
            starting_time=0.0,
            rate=RATE,
            description="EEG recorded differentially from 2 electrodes, resulting in a single EMG signal",
            comments = "2000 Hz, in V."
        )
        nwb_file.add_acquisition(es_eeg)

    return {"lfp": es_lfp, "emg": es_emg, "eeg": es_eeg}




def extract_lfp_signal(csv_data_row):
    """
    Extract and merge LFP signals from multiple shanks into one array of shape (T, n_channels).

    Parameters
    ----------
    csv_data_row : pd.Series
        Row from the CSV file containing session data eg LFPs

    Returns
    -------
    np.ndarray
        Array of shape (n_timepoints, n_channels)
    list
        Names of LFP regions (columns)
    """


    All_LFP = ["PtA", "dCA1", "mPFC", "wM1", "wS1", "wS2", "antM1"]
    LFPs = csv_data_row[All_LFP].dropna()

    if pd.notna(csv_data_row["EMG"]):
        EMG = list(map(float, csv_data_row["EMG"].split(";")))
    else:
        EMG = None
    
    if pd.notna(csv_data_row["EEG"]):
        EEG = list(map(float, csv_data_row["EEG"].split(";")))
    else:
        EEG = None

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        parsed_LFPs = LFPs.apply(lambda x: list(map(float, x.split(";"))))

    lfp_matrix = np.vstack(parsed_LFPs.values)

    return lfp_matrix.T, list(parsed_LFPs.index) , EMG , EEG # lfp.T to have shape (n_timepoints, n_channels) and regions as columns



def remove_nwb_files(folder_path):
    for filename in os.listdir(folder_path):
        if filename.endswith('.nwb'):
            file_path = os.path.join(folder_path, filename)
            try:
                os.remove(file_path)
            except Exception as e:
                print(f"Erreur lors de la suppression de {file_path}: {e}")