
# 🧠 NWB Conversion Pipeline for PL Sessions

This project provides a conversion pipeline for behavioral and electrophysiological data from the article: **Pierre Le Merre et al., 2028, Cell Reports** (only chronic LFP data) into the standard **Neurodata Without Borders (NWB)** format.

## 📚 Reference

Pierre Le Merre et al., *Reward-Based Learning Drives Rapid Sensory
Signals in Medial Prefrontal Cortex and Dorsal
Hippocampus Necessary for Goal-Directed Behavior*, Cell Reports, 2028.  
👉 [DOI](https://pubmed.ncbi.nlm.nih.gov/29249287/)



## ⚙️ Features

- Reads `.mat` files containing raw data
- Converts to NWB structure including:
  - General metadata (subject, session…)
  - Time intervals (e.g., trials)
  - Behavioral data (licks, rewards…)
- Validates the NWB file after conversion



## 📁 Project Structure

```
NWB\_converter\_LeMerre
│
├── converters/
│   ├── acquisition\_to\_nwb.py
│   ├── behavior\_to\_nwb.py
│   ├── general\_to\_nwb.py
│   ├── intervals\_to\_nwb.py
│   ├── Initiation\_nwb.py
│   └── nwb\_saving.py
├── README.md
├── convert_to_nwb_for_PL.py  ← Main conversion script
````

## 💻 Work Environment

Follow the **< Environment setup >** instructions provided in [LSENS-Lab-Immersion repository](https://github.com/loris-fab/LSENS-Lab-Immersion.git), and include the link to it.


## 🧩 How to use

Whether you run the pipeline from the **terminal** or from **Jupyter**, it is essential to ensure that you are using the correct environment. If you are working in *Visual Studio Code*, follow the **< Verification >** steps in the [LSENS-Lab-Immersion repository](https://github.com/loris-fab/LSENS-Lab-Immersion.git) to confirm that you are using the right environment either in the terminal when executing the pipeline there, or in Jupyter when running it from notebooks. Once confirmed, you can proceed with the instructions further down to run the pipeline.

Now, please find below the key information

1. `output_folder` → directory where you want the NWB file to be saved
2. `Folder_sessions_info` → path to the directory containing session information
3. `Folder_general_info` → path to the directory containing general information
4. `mouses_name` → lets you specify one or more mouse names to process.

### Commande in the terminal
Run the following command in the terminal, replacing the arguments :

```bash
python convert_to_nwb_for_PL.py output_folder Folder_sessions_info Folder_general_info --mouses_name PL200 PL201 (...)
```

*Options:*
* `--mouses_name` : Name(s) of the mouse/session(s) to convert (default: all mice), separated by spaces (e.g., `--mouses_name PL200 PL201`).

*for exemple for window:* 
```bash
python convert_to_nwb_for_PL.py \
"//sv-nas1.rcp.epfl.ch/Petersen-Lab/z_LSENS/Share/Loris_Fabbro/PL/NWB_files" \
"//sv-nas1.rcp.epfl.ch/Petersen-Lab/analysis/Sylvain_Crochet/DATA_REPOSITORY/LeMerre_mPFC_2018/Chronic_LFPs_Preprocessed" \
"//sv-nas1.rcp.epfl.ch/Petersen-Lab/publications/2018/2018_LeMerre_Neuron/2018_LeMerre_Neuron_data/processed_data" \
--mouses_name PL200
```

### Run inside a Jupyter Notebook

You can also call the conversion function directly in a Jupyter Notebook without using the command line.
Simply import the function `convert_data_to_nwb_pl` from your script and call it with the proper arguments:

*for exemple for window:* 
```python
import importlib
import convert_to_nwb_for_PL


Folder_general_info_server= "//sv-nas1.rcp.epfl.ch/Petersen-Lab/publications/2018/2018_LeMerre_Neuron/2018_LeMerre_Neuron_data/processed_data"
Folder_sessions_info_server = "//sv-nas1.rcp.epfl.ch/Petersen-Lab/analysis/Sylvain_Crochet/DATA_REPOSITORY/LeMerre_mPFC_2018/Chronic_LFPs_Preprocessed"
output_folder_serveur = "//sv-nas1.rcp.epfl.ch/Petersen-Lab/z_LSENS/Share/Loris_Fabbro/PL/NWB_files"
importlib.reload(convert_to_nwb_for_PL)
nwb_path = convert_to_nwb_for_PL.convert_data_to_nwb_pl(output_folder= output_folder_serveur, Folder_sessions_info=Folder_sessions_info_server, Folder_general_info=Folder_general_info_server, mouses_name=["PL200"])
```

*Options:*
* `mouses_name` : Name(s) of the mouse/session(s) to convert (default: all mice).Use a Python list for multiple mice (e.g. ["PL200", "PL201"]).

### Outcome
If everything runs correctly, you should see an output similar to this:

```bash
**************************************************************************
-_-_-_-_-_-_-_-_-_-_-_-_-_-_- NWB conversion _-_-_-_-_-_-_-_-_-_-_-_-_-_-_
Converting data to NWB format for mouse:  ['PL200']
Conversion to NWB is finished: 100%|██████████| 1/1 [05:43<00:00, 343.61s/file]
**************************************************************************
```

For each session, an `.nwb` file named with the session identifier is written under `Detection Task/` or `Neutral Exposition/` inside the chosen output directory. The script prints a summary and performs NWB validation before final save.

## ✍️ Author

Project developed as part of a student project focused on organizing and converting neuroscience data from the above-mentioned publication.
Main code by **@loris-fab**

For any questions related to the code, please contact: loris.fabbro@epfl.ch


---


