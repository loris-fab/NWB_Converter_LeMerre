
from uuid import uuid4
import numpy as np
import h5py
from pynwb.device import Device
import numpy as np
import pandas as pd




def add_general_container(nwb_file, csv_data_row, regions):
    """
    Crée devices, groupes d’électrodes et électrodes (LFP/EMG/EEG),
    puis renvoie (electrode_table_region, channel_labels) dans l’ordre des colonnes de données.

    Paramètres
    ----------
    nwb_file : pynwb.NWBFile
    csv_data_row : pd.Series
    regions : list[str]
        LFP regions présentes (p.ex. issues de extract_lfp_signal)

    Retour
    ------
    electrode_table_region : pynwb.core.DynamicTableRegion
    channel_labels : list[str]
        Étiquettes de canaux dans l’ordre (ex : ["PtA","dCA1","...","EMG1","EMG2","EEG1","EEG2"])
    """

    # ##############################################################
    # 1. Add Device 
    # ##############################################################

    def get_or_create_device(name, manufacturer, description):
        if name in nwb_file.devices:
            return nwb_file.devices[name]
        dev = Device(name=name, description=description, manufacturer=manufacturer)
        nwb_file.add_device(dev)
        return dev

    ampli_lfp = get_or_create_device(
        "Amplifier LFP",
        "A-M Systems",
        "Extracellular amplifier for LFP recording - Model 3000 AC/DC Differential Amplifier (custom modified).",
    )
    electrode_device = get_or_create_device(
        "Tungsten Microelectrodes",
        "FHC",
        "High-impedance sharp tungsten microelectrodes (10–12 MΩ), 75 μm shaft diameter; stereotaxically implanted individually. Model UEWSCGSELNND.",
    )
    eeg_device = get_or_create_device(
        "EEG Device",
        "N/A",
        "Surface electrodes for EEG recording (contralateral hemisphere).",
    )
    digitizer = get_or_create_device(
        "Digitizer",
        "LDS Nicolet",
        "Signals digitized and recorded at 2 kHz on Vision XP.",
    )
    _ = (ampli_lfp, digitizer)  # référencés pour traçabilité


    # ##############################################################
    # 2. Create Electrode Group and Add electrodes to the NWB file
    # ##############################################################

    def get_or_create_group(name, description, location, device):
        if name in nwb_file.electrode_groups:
            return nwb_file.electrode_groups[name]
        return nwb_file.create_electrode_group(
            name=name, description=description, location=location, device=device
        )

    lfp_group = get_or_create_group(
        "LFP",
        "LFP via sharp tungsten microelectrodes (~10 MΩ). Ref: cerebellar silver wire; band-pass 0.1–1000 Hz.",
        "Multiple (PtA, dCA1, mPFC, wM1, wS1, wS2, antM1) using interaural coordinates (Paxinos and Franklin 2008)",
        electrode_device,
    )
    emg_group = get_or_create_group(
        "EMG",
        "Differential EMG from neck muscles; band-pass 10–20000 Hz.",
        "Neck muscles",
        electrode_device,
    )
    eeg_group = get_or_create_group(
        "EEG",
        "Differential EEG from parietal and frontal; band-pass 0.1–1000 Hz.",
        "Cortex dura (parietal & frontal) using interaural coordinates (Paxinos and Franklin 2008)",
        eeg_device,
    )


    MM = 1 #because the coordinates are in mm

    AREA_COORDS_MM = {
        "wS1": (1.95, 3.5, 0.5),
        "wS2": (2.1, 4.2, 0.5),
        "wM1": (4.8, 1.0, 0.4),
        "PtA": (1.85, 1.6, 0.5),
        "mPFC": (5.8, 0.3, 1.85),
        "dCA1": (1.3, 2.0, 1.3),
        "antM1": (np.nan, np.nan, np.nan),
    }
    EEG_POS_MM = {
        "parietal": (2.0, 1.5, 0.0),
        "frontal":  (5.3, 1.5, 0.0),
    }

    def _last_row_index():
        """Index de la dernière ligne de la table electrodes (0 based)."""
        try:
            return len(nwb_file.electrodes.id.data) - 1
        except Exception:
            # fallback selon versions de pynwb
            return nwb_file.electrodes.table.length - 1

    electrode_indices = []
    channel_labels = []

    # ---------- LFP ----------
    ALL_LFP = ["PtA", "dCA1", "mPFC", "wM1", "wS1", "wS2", "antM1"]
    lfp_regions_present = [r for r in regions if r in ALL_LFP]
    if len(lfp_regions_present) > 7:
        lfp_regions_present = lfp_regions_present[:7]  # LFP:7 max

    for r in lfp_regions_present:
        ap, lat, depth = AREA_COORDS_MM.get(r, (np.nan, np.nan, np.nan))
        nwb_file.add_electrode(
            x=np.nan if np.isnan(ap) else ap * MM,
            y=np.nan if np.isnan(lat) else lat * MM,
            z=np.nan if np.isnan(depth) else depth * MM,

            imp=10e6,  # Ohms (~10 MΩ)
            location=r,
            filtering="0.1–1000 Hz band-pass",
            group=lfp_group,
            reference="Cerebellar silver wire",
            group_name="LFP",
        )
        electrode_indices.append(_last_row_index())
        channel_labels.append(r)

    # ---------- EMG (2 si présent) ----------
    has_emg = pd.notna(csv_data_row.get("EMG", None))
    if has_emg:
        for i in range(2):
            nwb_file.add_electrode(
                x=np.nan, y=np.nan, z=np.nan,
                imp=np.nan,
                location="Neck muscles",
                filtering="10–20000 Hz band-pass",
                group=emg_group,
                reference="Differential (neck 1 ↔ neck 2)",
                group_name="EMG",
            )
            electrode_indices.append(_last_row_index())
            channel_labels.append(f"EMG{i+1}")

    # ---------- EEG (2 si présent) ----------
    has_eeg = pd.notna(csv_data_row.get("EEG", None))
    if has_eeg:
        for name in ("parietal", "frontal"):
            ap, lat, depth = EEG_POS_MM[name]
            nwb_file.add_electrode(
                x=ap * MM, y=lat * MM, z=depth * MM,
                imp=np.nan,
                location=f"EEG {name} (contralateral; dura)",
                filtering="0.1–1000 Hz band-pass",
                group=eeg_group,
                reference="Differential (parietal ↔ frontal)",
                group_name="EEG",
            )
            electrode_indices.append(_last_row_index())
            channel_labels.append("EEG1" if name == "parietal" else "EEG2")

    # ##############################################################
    # 3. Return region 
    # ##############################################################

    electrode_table_region = nwb_file.create_electrode_table_region(
        region=electrode_indices,
        description="Electrodes used for LFP/EMG/EEG in this session (order matches data columns).",
    )

    return electrode_table_region, channel_labels
