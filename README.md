
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
├── README.md
├── convert_to_nwb_for_PL.py  ← Main conversion script
````

---

## 🚀 Usage

Follow the environment setup instructions provided in [LSENS-Lab-Immersion repository](https://github.com/loris-fab/LSENS-Lab-Immersion.git), and include the link to it.

## 🧩 How to use
Run the following command in the terminal, replacing `csv_path` with the`.csv` file I created, located in the data storage of the LSENS laboratory at EPFL, and `output_folder` with the directory where you want the NWB file to be saved. `--mouses_name` lets you specify one or more mouse names to process, separated by spaces (e.g., `--mouses_name PL200 PL201`).


```bash
python convert_to_nwb_for_PL.py output_folder --mouses_name PL200 PL201 (...)
```
*Options:*
* `--mouses_name` : Name of the mouse/session to convert (default: all sessions)


If everything runs correctly, you should see an output similar to this:

```bash
**************************************************************************
-_-_-_-_-_-_-_-_-_-_-_-_-_-_- NWB conversion _-_-_-_-_-_-_-_-_-_-_-_-_-_-_
Loading data PL202 ...: 100%|██████████| 1/1 [01:05<00:00, 65.93s/file]
Converting data to NWB format for mouse: ['PL202']
Processing : 16it [00:43,  2.70s/it, PL202]                      
**************************************************************************
```

## ✍️ Author

Project developed as part of a student project focused on organizing and converting neuroscience data from the above-mentioned publication.
Main code by **@loris-fab**

For any questions related to the code, please contact: loris.fabbro@epfl.ch


---


