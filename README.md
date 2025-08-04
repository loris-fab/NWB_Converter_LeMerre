
# 🧠 NWB Conversion Pipeline for PL Sessions

This project provides a conversion pipeline for behavioral and electrophysiological data from the article: **Pierre Le Merre et al., 2028, Cell Reports** (only chronic LFP data) into the standard **Neurodata Without Borders (NWB)** format.

## 📚 Reference

Pierre Le Merre et al., *Reward-Based Learning Drives Rapid Sensory
Signals in Medial Prefrontal Cortex and Dorsal
Hippocampus Necessary for Goal-Directed Behavior*, Cell Reports, 2028.  
👉 [https](https://pubmed.ncbi.nlm.nih.gov/29249287/)



## ⚙️ Features

- Reads `.csv` files containing raw data
- Converts to NWB structure including:
  - General metadata (subject, session…)
  - Time intervals (e.g., trials)
  - Behavioral data (licks, rewards…)
  - Optional analysis containers (e.g., LFPmean)
- Validates the NWB file after conversion



## 📁 Project Structure

```
NWB\_converter\_LeMerre
│
├── converters/
│   ├── acquisition\_to\_nwb.py
│   ├── analysis\_to\_nwb.py
│   ├── behavior\_to\_nwb.py
│   ├── general\_to\_nwb.py
│   ├── intervals\_to\_nwb.py
│   ├── Initiation\_nwb.py
│   └── nwb\_saving.py
├── requirement.txt
├── README.md
├── convert_to_nwb_for_PL.py  ← Main conversion script
````

---

## 🚀 Usage

Create environment :
```bash
conda create -n nwb_env python=3.9

conda activate nwb_env
```

Install dependencies with:

```bash
pip install -r requirement.txt
```


## 🧩 How to use
Run the following command in the terminal, replacing `csv_path` with the`.csv` file I created, located in the data storage of the LSENS laboratory at EPFL, and `output_folder` with the directory where you want the NWB file to be saved. `--mouses_name` lets you specify one or more mouse names to process, separated by spaces (e.g., `--mouses_name PL200 PL201`).


```bash
python convert_to_nwb_for_PL.py csv_path output_folder --mouses_name PL200 PL201 (...)
```
*Options:*
* `--mouses_name` : Name of the mouse/session to convert (default: all sessions)
* `--psth_window`: time window for PSTH (default: -0.2 0.5 seconds)

 for exemple :

```bash
python convert_to_nwb_for_PL.py Subject_Session_Selection.csv data --mouses_name PL200
```


If everything runs correctly, you should see an output similar to this:

```bash
**************************************************************************
-_-_-_-_-_-_-_-_-_-_-_-_-_-_- NWB conversion _-_-_-_-_-_-_-_-_-_-_-_-_-_-_
Converting data to NWB format for mouse: ['PL200_20140619', 'PL200_20140620', 'PL200_20140621', 'PL200_20140622', 'PL200_20140623', 'PL200_20140624', 'PL200_20140625', 'PL200_20140626']
📃 Creating configs for NWB conversion :
📑 Created NWB files :
     o 📌 Add general metadata
         - Subject metadata
         - Session metadata
         - Device metadata
         - Extracellular electrophysiology metadata
     o 📶 Add acquisition container
     o ⚙️ Add processing container
         - Behavior data
         - No ephys data for AN sessions
         - Analysis complementary information
             > Added LFP_mean_across_all_units to analysis module
             > Added global_LFP to analysis module
     o ⏸️ Add intervall container

🔎 Validating NWB file before saving...
     o ✅ File data/PL200_20140619_142055.nwb is valid, no errors detected and saved successfully.
     o ✅ File data/PL200_20140620_153120.nwb is valid, no errors detected and saved successfully.
     o ✅ File data/PL200_20140621_161645.nwb is valid, no errors detected and saved successfully.
     o ✅ File data/PL200_20140622_160902.nwb is valid, no errors detected and saved successfully.
     o ✅ File data/PL200_20140623_142533.nwb is valid, no errors detected and saved successfully.
     o ✅ File data/PL200_20140624_150428.nwb is valid, no errors detected and saved successfully.
     o ✅ File data/PL200_20140625_142129.nwb is valid, no errors detected and saved successfully.
     o ✅ File data/PL200_20140626_142654.nwb is valid, no errors detected and saved successfully.

**************************************************************************
```

## ✍️ Author

Project developed as part of a student project focused on organizing and converting neuroscience data from the above-mentioned publication.
Main code by **@loris-fab**

For any questions related to the code, please contact: loris.fabbro@epfl.ch


---


