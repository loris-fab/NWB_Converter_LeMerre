
from uuid import uuid4
import numpy as np
import h5py
from pynwb.device import Device
import numpy as np
import pandas as pd




def add_general_container(nwb_file, regions):
    """
    Create general NWB structures (devices, electrode group) and register electrodes
    for the provided LFP regions, then return an electrode table region aligned with
    data columns.

    Args:
        nwb_file (pynwb.file.NWBFile): Target NWB file to populate.
        regions (Sequence[str]): LFP region labels to register (subset of
            ["PtA", "dCA1", "mPFC", "wM1", "wS1", "wS2", "antM1"]).

    Returns:
        tuple[pynwb.core.DynamicTableRegion, list[str]]:
            - electrode_table_region: DynamicTableRegion pointing to the added electrodes
              (order matches the returned `channel_labels` and expected data columns).
            - channel_labels: List of region labels actually added, in the same order.
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
        " Differential extracellular amplifier for LFP recording - Model 3000 AC/DC Differential Amplifier (custom modified).",
    )
    electrode_device = get_or_create_device(
        "Tungsten Microelectrodes",
        "FHC",
        "High-impedance sharp tungsten microelectrodes (10–12 MΩ), 75 μm shaft diameter; stereotaxically implanted individually. Model UEWSCGSELNND.",
    )

    digitizer = get_or_create_device(
        "Digitizer",
        "LDS Nicolet",
        "Signals digitized and recorded at 2 kHz on Vision XP.",
    )
    _ = (ampli_lfp, digitizer)  


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
        "Local Field Potential (LFP) recorded via sharp tungsten microelectrodes (~10 MΩ) : Reference: cerebellar silver wire; band-pass filtering  0.1–1000 Hz.",
        "2-6 areas from (PtA, dCA1, mPFC, wM1, antM1, wS1, wS2,) using interaural stereotaxic coordinates (Paxinos and Franklin 2008)",
        electrode_device,
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

    def _last_row_index():
        """Index de la dernière ligne de la table electrodes (0 based)."""
        try:
            return len(nwb_file.electrodes.id.data) - 1
        except Exception:
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
        )
        electrode_indices.append(_last_row_index())
        channel_labels.append(r)


    # ##############################################################
    # 3. Return region 
    # ##############################################################

    electrode_table_region = nwb_file.create_electrode_table_region(
        region=electrode_indices,
        description="Electrodes used for LFP in this session (order matches data columns).",
    )

    return electrode_table_region, channel_labels
