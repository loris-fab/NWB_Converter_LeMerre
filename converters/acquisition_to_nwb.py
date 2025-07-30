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
    Crée 3 ElectricalSeries séparées (LFP, EMG, EEG) si les données/électrodes existent.
    - rate = 2000.0
    - unit = "V"

    Paramètres
    ----------
    nwb_file : pynwb.NWBFile
    lfp_array : (T, n_lfp)
    electrode_region_all : DynamicTableRegion (tous les canaux créés)
    channel_labels : list[str]  (renvoyé par add_general_container)
    emg : (T,) ou (T,2) ou None
    eeg : (T,) ou (T,2) ou None
    """


    RATE = 2000.0
    UNIT = "V"
    LFP_LABELS = ["PtA", "dCA1", "mPFC", "wM1", "wS1", "wS2", "antM1"]

    # --- util ---
    def _to_2d(a):
        a = np.asarray(a, dtype=float)
        return a[:, None] if a.ndim == 1 else a

    def _region_indices(dtr):
        # indices globaux des électrodes
        try:
            return list(dtr.data[:])
        except Exception:
            try:
                return list(dtr.region.data[:])
            except Exception:
                return list(dtr.region)

    # mapping label -> index global d'électrode
    idx_all = _region_indices(electrode_region_all)
    label_to_idx = dict(zip(channel_labels, idx_all))

    # ---------- LFP ----------
    lfp_array = _to_2d(lfp_array)
    lfp_labels = [lab for lab in channel_labels if lab in LFP_LABELS][:lfp_array.shape[1]]
    if lfp_array.shape[1] != len(lfp_labels):
        raise ValueError(
            f"Incohérence LFP: data a {lfp_array.shape[1]} colonnes, labels LFP trouvés {len(lfp_labels)}."
        )
    lfp_idx = [label_to_idx[lab] for lab in lfp_labels]
    lfp_region = nwb_file.create_electrode_table_region(lfp_idx, "LFP electrodes")
    es_lfp = ElectricalSeries(
        name="ElectricalSeries_LFP",
        data=lfp_array,
        electrodes=lfp_region,
        starting_time=0.0,
        rate=RATE,
        unit=UNIT,
        description="Traces LFP (brut), 2000 Hz, en V."
    )
    nwb_file.add_acquisition(es_lfp)

    # ---------- EMG (optionnel) ----------
    es_emg = None
    emg_labels = [lab for lab in channel_labels if lab.startswith("EMG")]
    if emg is not None and len(emg_labels) > 0:
        emg = _to_2d(emg)
        # Ajuste au nombre d'électrodes disponibles (1 ou 2 typiquement)
        if emg.shape[1] > len(emg_labels):
            raise ValueError(f"EMG a {emg.shape[1]} colonnes mais {len(emg_labels)} électrodes EMG définies.")
        emg_labels = emg_labels[:emg.shape[1]]
        emg_idx = [label_to_idx[lab] for lab in emg_labels]
        emg_region = nwb_file.create_electrode_table_region(emg_idx, "EMG electrodes")
        es_emg = ElectricalSeries(
            name="ElectricalSeries_EMG",
            data=emg,
            electrodes=emg_region,
            starting_time=0.0,
            rate=RATE,
            unit=UNIT,
            description="Traces EMG (brut), 2000 Hz, en V."
        )
        nwb_file.add_acquisition(es_emg)

    # ---------- EEG (optionnel) ----------
    es_eeg = None
    eeg_labels = [lab for lab in channel_labels if lab.startswith("EEG")]
    if eeg is not None and len(eeg_labels) > 0:
        eeg = _to_2d(eeg)
        if eeg.shape[1] > len(eeg_labels):
            raise ValueError(f"EEG a {eeg.shape[1]} colonnes mais {len(eeg_labels)} électrodes EEG définies.")
        eeg_labels = eeg_labels[:eeg.shape[1]]
        eeg_idx = [label_to_idx[lab] for lab in eeg_labels]
        eeg_region = nwb_file.create_electrode_table_region(eeg_idx, "EEG electrodes")
        es_eeg = ElectricalSeries(
            name="ElectricalSeries_EEG",
            data=eeg,
            electrodes=eeg_region,
            starting_time=0.0,
            rate=RATE,
            unit=UNIT,
            description="Traces EEG (brut), 2000 Hz, en V."
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