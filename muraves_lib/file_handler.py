__all__ = ['decompress', 'parse_slow_control', 'parse_log_file', 'parse_conteggi']

import logging
import os
import gzip
import shutil
import tempfile
from pathlib import Path
import json

def read_run_info_from_json(filename, subdict='slowcontrol') -> dict:
    """
    Reads a JSON file containing run information and extracts:
    - timestamp
    - working point (wp)
    - temperature
    - humidity
    - trigger rate
    - accidental rate
    """

    # ---- Read JSON file ----
    with open(filename, "r") as f:
        data = json.load(f)

    # ---- Extract slow control values ----
    slowcontrol = data.get(subdict, [{}])
    sc = slowcontrol[0] if slowcontrol else {}

    wp = sc.get("wp")
    temperature = sc.get("temperature")
    humidity = sc.get("humidity")    # may be missing

    # ---- Extract log values ----
    log = data.get("log", {})
    timestamp = log.get("timestamp")
    trigger_rate = log.get("trigger_rate")
    accidental_rate = log.get("accidental_rate")

    return {
        "timestamp": timestamp,
        "working_point": wp,
        "temperature": temperature,
        "humidity": humidity,
        "trigger_rate": trigger_rate,
        "accidental_rate": accidental_rate
    }

def decompress(filename, destination_path = "/workspace/tmp/DECOMPRESSED/") -> str:
    gz_path = Path(filename)
    gz_filename = gz_path.stem # removes .gz extension
    out_path_root = Path(destination_path)
    os.makedirs(out_path_root, exist_ok=True)
    out_filename = out_path_root / gz_filename
    #out_path = str(out_filename.parent)
    #os.makedirs(out_path, exist_ok=True)


    ctrl = 0
    if not gz_path.exists():
        logging.error(f"{gz_path} does not exist.")   
        ctrl=1     
    if out_filename.exists():
        logging.info(f"Uncompressed file {out_filename} already exists. Skipping.")
        #sys.exit(0)
        ctrl=1
    err = None
    if ctrl==0:
        try:
            # create temp file in the same directory to ensure atomic rename
            with tempfile.NamedTemporaryFile(delete=False, dir=out_filename.parent) as tmp_file:
                tmp_name = tmp_file.name
                with gzip.open(gz_path, 'rb') as f_in:
                    shutil.copyfileobj(f_in, tmp_file)
    
            # atomic move to final destination
            Path(tmp_name).rename(out_filename)
            logging.debug(f"Successfully uncompressed {gz_path} → {out_filename}")
        except Exception as e:
            # if something goes wrong, remove the temp file
            if 'tmp_name' in locals() and Path(tmp_name).exists():
                Path(tmp_name).unlink()
            err = f'Decompression failed'
            logging.error(f"{err} for {gz_path}: {e}")
            file_path = Path(out_filename)
            file_path.touch(exist_ok=True)
    return str(out_filename), err


def parse_slow_control(filename, target_run)-> dict:
    """
    Reads a one-line slow-control TSV file.
    Requirements:
      - File has exactly one meaningful line
      - Column 0 must match target_run
      - Must have at least 44 columns

    Returns:
      {"temperature": float, "wp": float, "tr": float} 
    or None if validation fails.
    """

    with open(filename, "r") as f:

        return_list = []
        logging.debug("Parsing SLOWCONTROL file.")
    
        for line in f:
            slowcontrol_dict = None
        #line = f.readline().strip()
            cols = line.strip().split("\t")
            # Check column count
            if len(cols) < 44:
                logging.warning(
                    f"SLOWCONTROL file has only {len(cols)} columns, expected at least 44: {filename}"
                )
                #slowcontrol_dict =  None

            # Validate run number
            try:
                run = int(float(cols[0]))
                if run != target_run:
                    logging.info(
                    f"SLOWCONTROL file, expected run number mismatch: file has run {run}, expected {target_run}: {filename}"
                )
            except ValueError:
                logging.error("Run not found")
                run = None
            try:
                ts = int(cols[1])
                try:
                    from datetime import datetime
                    dt = datetime.fromtimestamp(ts / 1000.0)
                    dt_hour = dt.replace(second=0, microsecond=0)
                    ts_day = dt_hour.strftime("%Y-%m-%d")
                    ts_hour = dt_hour.strftime("%H:%M")
                except Exception as e:
                    logging.error(f"SLOWCONTROL file, error converting timestamp: {e}")
                    ts_day = None
                    ts_hour = None
            except ValueError:
                logging.error("SLOWCONTROL file, could not retrieve effective time stamp.")
                ts = None
            try:
                temperature = float(cols[3])
            except ValueError:
                logging.warning("SLOWCONTROL file, could not retrieve effective temperature.")
                temperature = None
            try:
                humidity = float(cols[4])
            except ValueError:
                logging.warning("SLOWCONTROL file, could not retrieve humidity (%).")
                humidity = None
            try:
                wp = float(cols[5])
            except ValueError:
                logging.warning("SLOWCONTROL file, could not retrieve working point temperature.")
                wp = None
            try:
                tr = float(cols[43])
            except:
                logging.warning("SLOWCONTROL file, could not retrieve tr.")
                tr = None
            try:
                ar = float(cols[44])
            except:
                logging.warning("SLOWCONTROL file, could not retrieve accidental rate.")
                ar = None

            slowcontrol_dict = {
                "run": run,
                "timestamp": ts,
                "day": ts_day,
                "hour": ts_hour,
                "temperature": temperature,
                "humidity": humidity,
                "wp": wp,
                "tr": tr,
                "ar": ar,
            }
            return_list.append(slowcontrol_dict)
    return return_list


def parse_log_file(file_path, target_run)-> dict:
    """
    Reads a LOG file and extracts:
      - timestamp  (col 1)
      - trigger_rate (col 43)
      - accidental_rate (col 44)
      - OR32 counts (cols 45–60 : 16 values)
    The run (col 0) must match target_run.
    
    Returns a dictionary with these values (run excluded).
    """
    logging.debug("Parsing LOG file.")
    data = {}
    or32_counts = {}
    run = None

    with open(file_path, 'r') as f:
        raw_lines = f.readlines()

    for i, raw_line in enumerate(raw_lines):
        line = raw_line.strip()
        if line.startswith("Current run"):
           run = int(line.split(":")[1].strip())
        elif line.startswith("Timestamp"):
            data["timestamp"] = int(line.split(":")[1].strip())
        elif line.startswith("Trigger rate"):
            data["trigger_rate"] = float(line.split(":")[1].strip())
        elif line.startswith("Accidental rate"):
            data["accidental_rate"] = float(line.split(":")[1].strip())
        elif line.startswith("OR32 Counts"):
            # OR32 counts start after this line
            for count_line in raw_lines[i+1:]:
                cl = count_line.strip()
                if not cl:
                    continue
                
                # Stop if we reach next section (non-numeric lines)
                parts = cl.split()
                if len(parts) != 2:
                    break
                
                key, value = parts
                or32_counts[int(key)] = float(value)

            break  # all remaining lines are OR32 counts

    data["OR32_counts"] = or32_counts
    if int(run)!=int(target_run):
        logging.error(f"Run mismatch: file has run {run}, expected {target_run}.")
        data = None

    return data

def parse_conteggi(file_path) -> list:
    data = []
    logging.debug("Parsing CONTEGGI file.")
    with open(file_path, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 36:
                continue  # skip incomplete lines
            
            board_idx = int(parts[0])
            channels = [int(x) for x in parts[1:33]]  # 32 channels
            total_or = int(parts[33])
            flag = int(parts[34])
            timestamp = int(parts[35])

            data.append( {
                "board": board_idx,
                "channels": channels,
                "total_or": total_or,
                "flag": flag,
                "timestamp": timestamp
            })

    return data






