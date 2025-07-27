from pynwb.ecephys import ElectricalSeries
import h5py
import numpy as np
import warnings

#####################################################
# Functions that handle LFP acquisition in NWB files
#####################################################

def add_lfp_acquisition(nwb_file, signal_array, electrode_region):
    """
    Add LFP signal to NWB file acquisition as an ElectricalSeries.

    Parameters
    ----------
    nwb_file : pynwb.NWBFile
    signal_array : np.ndarray, shape (n_timepoints, n_channels)
        The LFP signal data.
    electrode_region : DynamicTableRegion
    """
    sampling_rate = float(2000)

    e_series = ElectricalSeries(
        name="ElectricalSeries_LFP",                             
        data=signal_array,
        electrodes=electrode_region,
        starting_time=0.0,
        rate=sampling_rate,
        description=f"Raw acquisition traces: Local Field Potential from {signal_array.shape[1]} electrodes"
    )
    nwb_file.add_acquisition(e_series)


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

    All_LFP = ["EMG", "PtA", "dCA1", "mPFC", "wM1", "wS1", "wS2", "antM1"]
    LFPs = csv_data_row[All_LFP].dropna()

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        parsed_LFPs = LFPs.apply(lambda x: list(map(float, x.split(";"))))

    lfp_matrix = np.vstack(parsed_LFPs.values)

    return lfp_matrix.T, list(parsed_LFPs.index)

