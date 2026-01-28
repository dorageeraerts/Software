__all__ = ['get_orientation', 'create_database']


from glob import glob
import pandas as pd
from muraves_lib import run_manager
from muraves_lib import file_handler
from tqdm import tqdm
import logging
import numpy as np
import json
import os

def get_orientation(run_number: int, color: str) -> str:
    """
    Docstring for get_orientation
    
    :param run_number: number of the run
    :type run_number: int
    :param color: color of the hodoscope
    :type color: str
    :return: orientation of the hodoscope VESUVIUS or FRE_SKY or unknown
    :rtype: str

    Determines the orientation of the hodoscope based on the run number and color.
    The orientation is determined using predefined run number ranges for each color.
    These numbers are hardcoded based on comunication with the operations team and on Yanwen's presentation: 
    https://agenda.infn.it/event/47877/contributions/272281/attachments/139076/209525/DataAnalysis%202.pdf , Slide 7.
    """
    
    orientation = "unknown"

    if color == "BLU":
        if run_number in np.arange(2800, 8819, 1) : #8819 is excluded
            orientation = "vesuvius"
        elif run_number >= 8819:
            orientation = "fre_sky"
    elif color == "NERO":
        if run_number in np.arange(2500, 5601, 1) or run_number >=7894: #5601 is excluded
            orientation = "vesuvius"
        elif run_number in np.arange(5636, 7894, 1): #7894 is excluded
            orientation = "fre_sky" 
    elif color == "ROSSO":
        if run_number in np.arange(4000, 10061, 1) or run_number >= 13458: #10061 is excluded
            orientation = "vesuvius"
        elif run_number in np.arange(10061, 13458, 1): #13458 is excluded
            orientation = "fre_sky"
    
    return orientation

# save info in a dataframe alice function
def create_database(raw_gz_path, output_path='./', color_list = ['BLU', 'NERO', 'ROSSO'], suffix = '' ):

    rows=[]
    for color in color_list:
        #List the files to parse
        path = f'{raw_gz_path}/{color}/'
        file_list = glob(f"{path}/SLOWCONTROL*{suffix}*.gz")
        logging.info(f"Number of files to decompress and parse:{len(file_list)}")

        #Loop over the files
        for file in tqdm(file_list): #tqdm for progress bar
            decompressed_file, err = file_handler.decompress(file)
            if err is None:
                run_filename = decompressed_file.split('run')[-1]
                slowcontrol_data = file_handler.parse_slow_control(decompressed_file, target_run=run_filename)
                timestamp = slowcontrol_data[0]['timestamp'] 
                run = slowcontrol_data[0]['run']+1 
                day = slowcontrol_data[0]['day'] 
                hour = slowcontrol_data[0]['hour']
                working_point = slowcontrol_data[0]['wp']
                temperature = slowcontrol_data[0]['temperature'] 
                humidity = slowcontrol_data[0]['humidity']
                trigger_rate = slowcontrol_data[0]['tr'] 
                accidental_rate = slowcontrol_data[0]['ar']
                orientation = run_manager.get_orientation(run, color)
                rows.append({
                    "hodoscope": color,
                    "run": run,
                    "timestamp": timestamp,
                    "day": day,
                    "hour": hour,
                    "wp":working_point,
                    "temperature": temperature,
                    "humidity": humidity,
                    "trigger_rate": trigger_rate,
                    "accidental_rate": accidental_rate,
                    "orientation": orientation,
                })
            else:
                print(f"Error decompressing file: {err}")

    df = pd.DataFrame(rows)
    #df.to_csv(f"{output_path}/run_index.csv", index=True)
    df.to_pickle(f"{output_path}/data_scan.pkl")
    logging.debug(df.head())

    return df

def create_database_chunk(
    raw_gz_path,
    output_path='./',
    color_list=['BLU', 'NERO', 'ROSSO'],
    suffix='',
    chunk_size=2000,
    state_file='state.json'
):
    # Carica lo stato precedente, se esiste
    if os.path.exists(state_file):
        with open(state_file, 'r') as f:
            state = json.load(f)
    else:
        state = {"processed_files": [], "last_chunk": 0}

    rows = []
    processed_files = set(state["processed_files"])

    for color in color_list:
        path = f'{raw_gz_path}/{color}/'
        file_list = glob(f"{path}/SLOWCONTROL*{suffix}*.gz")
        file_list = [f for f in file_list if f not in processed_files]

        logging.info(f"Number of files to process for {color}: {len(file_list)}")

        for file in tqdm(file_list):
            decompressed_file, err = file_handler.decompress(file)
            if err is None:
                run_filename = decompressed_file.split('run')[-1]
                slowcontrol_data = file_handler.parse_slow_control(decompressed_file, target_run=run_filename)
                rows.append({
                    "hodoscope": color,
                    "run": slowcontrol_data[0]['run'] + 1,
                    "timestamp": slowcontrol_data[0]['timestamp'],
                    "day": slowcontrol_data[0]['day'],
                    "hour": slowcontrol_data[0]['hour'],
                    "wp": slowcontrol_data[0]['wp'],
                    "temperature": slowcontrol_data[0]['temperature'],
                    "humidity": slowcontrol_data[0]['humidity'],
                    "trigger_rate": slowcontrol_data[0]['tr'],
                    "accidental_rate": slowcontrol_data[0]['ar'],
                    "orientation": run_manager.get_orientation(slowcontrol_data[0]['run'] + 1, color),
                })
                processed_files.add(file)
            else:
                logging.error(f"Error decompressing file: {err}")

            # Salva ogni chunk_size file
            if len(processed_files) % chunk_size == 0:
                df = pd.DataFrame(rows)
                df.to_pickle(f"{output_path}/data_scan_chunk_{state['last_chunk'] + 1}.pkl")
                state["last_chunk"] += 1
                state["processed_files"] = list(processed_files)
                with open(state_file, 'w') as f:
                    json.dump(state, f)
                logging.info(f"Saved chunk {state['last_chunk']} with {len(processed_files)} files processed.")

    # Salva l'ultimo chunk se necessario
    if len(rows) > 0:
        df = pd.DataFrame(rows)
        df.to_pickle(f"{output_path}/data_scan_chunk_{state['last_chunk'] + 1}.pkl")
        state["last_chunk"] += 1
        state["processed_files"] = list(processed_files)
        with open(state_file, 'w') as f:
            json.dump(state, f)

    # Unisci tutti i chunk in un unico file (opzionale)
    all_chunks = []
    for i in range(1, state["last_chunk"] + 1):
        chunk_df = pd.read_pickle(f"{output_path}/data_scan_chunk_{i}.pkl")
        all_chunks.append(chunk_df)
    final_df = pd.concat(all_chunks, ignore_index=True)
    final_df.to_pickle(f"{output_path}/data_scan_final.pkl")

    return final_df
