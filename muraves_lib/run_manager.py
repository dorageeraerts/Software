__all__ = ['get_orientation', 'create_database', 'update_database', 'read_summary_file', 'is_good_run']


from glob import glob
import pandas as pd
from muraves_lib import run_manager
from muraves_lib import file_handler
from tqdm import tqdm
import logging
import numpy as np
import json
import os
from typing import NamedTuple, List
import signal
from pathlib import Path
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
def create_database(raw_gz_path, output_path='./', output_filename = "data_scan", color_list = ['BLU', 'NERO', 'ROSSO'], suffix = '' ):

    rows=[]
    for color in color_list:
        #List the files to parse
        path = f'{raw_gz_path}/{color}/'
        file_list = glob(f"{path}/SLOWCONTROL*{suffix}*.gz")
        logging.info(f"Number of files to decompress and parse:{len(file_list)}")

        #Loop over the files
        for file in tqdm(file_list, position=0): #tqdm for progress bar
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
                dac10 = slowcontrol_data[0]['DAC10']
                or32_counts = slowcontrol_data[0]['OR32_counts']
                len_slowcontrol = slowcontrol_data[0]['length']
                orientation = run_manager.get_orientation(run, color)
                rows.append({
                    "color": color,
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
                    "DAC10": dac10,
                    "OR32_counts": or32_counts,
                    "slowcontrol_length": len_slowcontrol

                })
            else:
                logging.error(f"Error decompressing file: {err}")

    df = pd.DataFrame(rows)
    #df.to_csv(f"{output_path}/run_index.csv", index=True)
    df.to_pickle(f"{output_path}/{output_filename}.pkl")
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


def update_database(new_df, existing_df_path):
    # Check that 'run' column exists in new_df
    if 'run' not in new_df.columns:
        raise ValueError("'run' column not found in new dataframe")
    
    try:
        
        existing_df = pd.read_pickle(existing_df_path)
        
        # Check that 'run' column exists in existing_df
        if 'run' not in existing_df.columns:
            raise ValueError("'run' column not found in existing dataframe")
        
        # Identify new columns in new_df that don't exist in existing_df
        new_columns = [col for col in new_df.columns if col not in existing_df.columns]
        
        # Merge: keep all rows from existing_df and update/add from new_df
        updated_df = existing_df.merge(new_df, on='run', how='left', suffixes=('', '_new'))
        
        # For columns that exist in both, update with new values (preferring non-NaN)
        for col in new_df.columns:
            if col != 'run' and f'{col}_new' in updated_df.columns:
                updated_df[col] = updated_df[f'{col}_new'].fillna(updated_df[col])
                updated_df.drop(columns=[f'{col}_new'], inplace=True)
        
        # Add rows from new_df that don't exist in existing_df
        new_runs = new_df[~new_df['run'].isin(existing_df['run'])]
        if len(new_runs) > 0:
            logging.warning(f"New run found! Adding {len(new_runs)} new runs to the database.")
            updated_df = pd.concat([updated_df, new_runs], ignore_index=True)
        
        updated_df.to_pickle(existing_df_path)
        logging.info(f"Database updated with {len(new_df)} entries. New columns added: {new_columns}")
    except FileNotFoundError:
        logging.error(f"No pre-existing database has been found: {existing_df_path}. Impossible to update.")

class FileTimeoutError(Exception):
    pass

def _timeout_handler(signum, frame):
    raise FileTimeoutError("File read timed out")

class RunValidation(NamedTuple):
    """Result of run validation."""
    is_run_ok: bool
    is_parsing_failed: bool
    has_less_events: bool
    has_mismatches: bool
    has_unrecoverable_mismatches: bool
    mismatch_counter: int
    unrecoverable_mismatch_counter: int
    bit_flip_counter: int
    bit_missing_counter: int
    missing_run: bool
    mismatch_blockindex: List[int]
    unrecoverable_mismatches_details: dict
    type: str
    
    


def read_summary_file(filename, timeout_seconds=5) -> RunValidation:
    #print(f"Reading summary file: {filename}")

    try:
        with open(filename, "r", encoding="utf-8") as f:
            data = json.load(f)

    except (json.JSONDecodeError, OSError, FileTimeoutError) as e:
        print(f"File corrupted or unreadable: {filename} ({e})")

        return RunValidation(
            type=None,
            is_run_ok=None,
            is_parsing_failed=None,
            has_less_events=None,
            has_mismatches=None,
            has_unrecoverable_mismatches=None,
            mismatch_counter=None,
            unrecoverable_mismatch_counter = None,
            bit_flip_counter=None,
            bit_missing_counter=None,
            missing_run=None,
            mismatch_blockindex=None,
            unrecoverable_mismatches_details=None
        )

    # ---- Se arrivi qui il file è valido ----
    type = Path(filename).stem.split("_")[0] 
    assert type in ["PIEDISTALLI", "ADC"], f"Unexpected file type: {type} in file {filename}"

    missing_run = len(data.get("subruns", [])) == 0
    is_run_ok = data.get("status") == "ok"
    qc = data.get("quality_checks", {})

    has_less_events = True if (qc.get("has_less_nevents") == True) else False
    #import pdb;pdb.set_trace()
    has_mismatches = True if (qc.get("has_mismatches") == True) else False
    has_unrecoverable_mismatches = True if (qc.get("has_unrecoverable_errors") == True) else False
    mismatch_counter = qc.get("mismatch_counter", None)
    unrecoverable_mismatch_counter = qc.get("unrecoverable_mismatch_counter", None)
    mpl=qc.get("mismatches_per_line", None)
    mismatch_linenumbers = mpl.get("line_number", None) if mpl else None
    mismatch_eventnumbers = mpl.get("event_number", None) if mpl else None
    mismatch_blockindex = mpl.get("block_idx", None) if mpl else None
    unrecoverable_mismatches_details = qc.get("unrecoverable_mismatches_details", None)
    #mismatch_subevent_expected = mpl.get("subevent_number_expected", None) if mpl else None
    #mismatch_subevent_found = mpl.get("subevent_number_found", None) if mpl else None
    bit_flips = mpl.get("bit_flip", None) if mpl else None
    bit_missing = mpl.get("bit_missing", None) if mpl else None
    bit_flip_counter = sum(bit_flips) if bit_flips else None
    bit_missing_counter = sum(bit_missing) if bit_missing else None


    #bit_mismatch_test_passed = missing_run==False and has_less_events and ((has_mismatches == False) | ( (has_mismatches == True) & (has_unrecoverable_mismatches == False) ))
    is_parsing_failed = True if has_less_events==False and missing_run==False and is_run_ok==False else False

    return RunValidation(
        type=type,
        is_run_ok = is_run_ok,
        is_parsing_failed=is_parsing_failed,
        has_less_events=has_less_events,
        has_mismatches=has_mismatches,
        has_unrecoverable_mismatches=has_unrecoverable_mismatches,
        missing_run=missing_run,
        mismatch_counter=mismatch_counter,
        unrecoverable_mismatch_counter=unrecoverable_mismatch_counter,
        bit_flip_counter=bit_flip_counter,
        bit_missing_counter=bit_missing_counter,
        mismatch_blockindex=mismatch_blockindex,
        unrecoverable_mismatches_details=unrecoverable_mismatches_details,

    )

def is_good_run(filename) :

    try:
        with open(filename, "r", encoding="utf-8") as f:
            data = json.load(f)
        is_good = True if data.get("status")== "ok" else False
        file = "found" if len(data.get("subruns"))!= 0 else "not found"
    except FileNotFoundError:
        file =  "not found"
        is_good = "N/A"

    return file, is_good