from pynwb.ecephys import ElectricalSeries
import numpy as np
import warnings
import os



#####################################################
# Functions that handle LFP acquisition in NWB files
#####################################################

def add_acquisitions_3series(nwb_file, lfp_array, electrode_region_all, channel_labels):
    """
    Add LFP data as a single ElectricalSeries acquisition to an NWB file.

    This function builds an electrodes region from the provided mapping and writes
    one ElectricalSeries named "ElectricalSeries_LFP" sampled at 2000 Hz.

    Args:
        nwb_file (pynwb.file.NWBFile): Target NWB file object (opened for writing).
        lfp_array (np.ndarray): LFP data of shape (n_timepoints, n_channels), in Volts.
        electrode_region_all: Electrode table region/selection containing global electrode indices.
        channel_labels (Sequence[str]): Per-channel labels; must include a subset of
            {"PtA","dCA1","mPFC","wM1","wS1","wS2","antM1"} matching the columns in `lfp_array`.

    Returns:
        None
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
        description="LFPs recorded from multiple electrodes, bandpass filtering 0.1-1000 Hz",
        comments = "sampling rate 2000 Hz, in V."
    )
    nwb_file.add_acquisition(es_lfp)



    return None



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

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        parsed_LFPs = LFPs.apply(lambda x: list(map(float, x.split(";"))))

    lfp_matrix = np.vstack(parsed_LFPs.values)

    return lfp_matrix.T, list(parsed_LFPs.index)  # lfp.T to have shape (n_timepoints, n_channels) and regions as columns



def remove_nwb_files(folder_path):
    """
    Delete all `.nwb` files in the given directory (non-recursive, best effort).

    Args:
        folder_path (str | pathlib.Path): Directory to scan.

    Returns:
        None

    Side Effects:
        Removes files on disk. Prints an error message if a file cannot be deleted.
    """
    for filename in os.listdir(folder_path):
        if filename.endswith('.nwb'):
            file_path = os.path.join(folder_path, filename)
            try:
                os.remove(file_path)
            except Exception as e:
                print(f"Erreur lors de la suppression de {file_path}: {e}")