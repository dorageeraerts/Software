__all__ = ['decompress', 'parse_slow_control', 'parse_log_file', 'parse_conteggi', 'temp_file_manager', 'copy_on_system', 'temp_to_output']

import logging
import os
import gzip
import shutil
import tempfile
from pathlib import Path
import json
from contextlib import contextmanager
import time

logger = logging.getLogger(__name__)

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

def _atomic_replace(src: Path, dst: Path, retries: int = 3, backoff: float = 0.1):
    """Move *src* to *dst* atomically, retrying on transient failures.

    The source and destination **must live on the same filesystem**.  A
    simple ``os.replace`` is used (it always overwrites the target) but
    we wrap it in a retry loop so that a busy network filesystem has a
    chance to recover when many clients are renaming concurrently.  A
    directory lock can be added later if even retries are not enough.
    """
    for attempt in range(1, retries + 1):
        try:
            os.replace(src, dst)
            return
        except OSError as e:
            if attempt == retries:
                raise
            time.sleep(backoff * attempt)
        src.unlink(missing_ok=True)
        raise RuntimeError(f"Failed to replace {src} with {dst} after {retries} attempts")


def decompress(filename) -> str:

    gz_path = Path(filename)
    
    if not gz_path.exists():
        logger.error(f"{gz_path} does not exist.")   
        return None, f"File not found: {gz_path}"

    try:
        # create temp file in the same directory to ensure atomic rename
        with tempfile.NamedTemporaryFile(dir= "/tmp",delete=False) as tmp_file:
            tmp_path = Path(tmp_file.name)
            # copy decompressed content to temp file
            with gzip.open(gz_path, 'rb') as f_in:
                shutil.copyfileobj(f_in, tmp_file)
        
        return tmp_path, None
    except Exception as e:
        # if something goes wrong, remove the temp file
        if tmp_path and tmp_path.exists():
            tmp_path.unlink()

        return None, f"Decompression failed for {gz_path}: {e}"



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
        logger.debug("Parsing SLOWCONTROL file.")
    
        for line in f:
            slowcontrol_dict = None
        #line = f.readline().strip()
            cols = line.strip().split("\t")
            # Check column count
            sc_lenght_ctrl = 0
            slowcontrol_lenght = len(cols)
            if slowcontrol_lenght < 60:
                sc_lenght_ctrl = 1
                #logging.warning(
                #    f"SLOWCONTROL file has only {len(cols)} columns, expected 60: {filename}"
                #)
            
            # Validate run number
            try:
                run = int(float(cols[0]))
                if run != target_run:
                    logger.info(
                    f"SLOWCONTROL file, expected run number mismatch: file has run {run}, expected {target_run}: {filename}"
                )
            except ValueError:
                logger.error("Run not found")
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
                    logger.error(f"SLOWCONTROL file, error converting timestamp: {e}")
                    ts_day = None
                    ts_hour = None
            except ValueError:
                logging.error("SLOWCONTROL file, could not retrieve effective time stamp.")
                ts = None
            try:
                temperature = float(cols[3])
            except ValueError:
                logger.warning("SLOWCONTROL file, could not retrieve effective temperature.")
                temperature = None
            try:
                humidity = float(cols[4])
            except ValueError:
                logger.warning("SLOWCONTROL file, could not retrieve humidity (%).")
                humidity = None
            try:
                wp = float(cols[5])
            except ValueError:
                logger.warning("SLOWCONTROL file, could not retrieve working point temperature.")
                wp = None

            if sc_lenght_ctrl ==1:
                tr = None
            else: tr = float(cols[43])

            if sc_lenght_ctrl ==1:
                ar = None
            else: ar = float(cols[44])
 
            if sc_lenght_ctrl ==1:
                    DAC10 = None
            else:
                DAC10  = []
                for i in range(6, 22):
                    DAC10.append(float(cols[i]))

            if sc_lenght_ctrl ==1:
                OR32_counts = None
            else:
                OR32_counts = []
                for i in range(45, 61):
                    OR32_counts.append(float(cols[i]))

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
                "DAC10": DAC10,
                "OR32_counts": OR32_counts,
                "length": len(cols),

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
    logger.debug("Parsing LOG file.")
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
        logger.error(f"Run mismatch: file has run {run}, expected {target_run}.")
        data = None

    return data

def parse_conteggi(file_path) -> list:
    data = []
    logger.debug("Parsing CONTEGGI file.")
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






@contextmanager
def temp_file_manager(delete_on_failure=True):
    """Context manager for safely creating a temporary file."""
    tmp_path = None
    try:
        with tempfile.NamedTemporaryFile(mode='w', dir="/tmp", delete=False) as tf:
            tmp_path = Path(tf.name)
            yield tmp_path  # give control back to the block
    except Exception as e:
        if delete_on_failure and tmp_path and tmp_path.exists():
            tmp_path.unlink()
        logger.error(f"Error creating temp file {tmp_path}: {e}")
        open(tmp_path, "w").close()
        

@contextmanager
def temp_to_output(output_path: Path, success_message= ""):
    """
    Context manager that:
    - creates a temp file in scratch_dir
    - yields its Path for writing
    - copies it safely to output_path on success
    - ensures output file exists even on failure
    - cleans scratch
    """

    output_path = Path(output_path)
    #scratch_dir = Path(scratch_dir)

    scratch_tmp = None
    final_tmp = output_path.with_suffix(output_path.suffix + ".tmp")

    try:
        # create scratch temp file
        with tempfile.NamedTemporaryFile(
            #mode="w", commented out because specific for txt file. 
            dir="/tmp",
            delete=False
        ) as tf:
            scratch_tmp = Path(tf.name)

        # give control back to user code
        yield scratch_tmp
        logger.debug(f"Successfully filled {scratch_tmp}, copying into {final_tmp}.")

    except Exception as e:
        logging.error(f"Failed filling {scratch_tmp}: {e}. Creating empty output at {final_tmp}.")
        # Ensure Snakemake sees the file
        scratch_tmp.touch(exist_ok=True)
        #raise  # propagate if you want job to fail
    finally:
        #import pdb; pdb.set_trace()
        # Always create output
        shutil.copy2(scratch_tmp, final_tmp)
        _atomic_replace(final_tmp, output_path, retries=10)
        # Always cleanup
        scratch_tmp.unlink(missing_ok=True)
        if success_message=="":
            logger.info(f"Successfully renamed temporary file into its final destination: {output_path} and temporary file deleted.")
        else:
            logger.info(success_message + f"\nSuccessfully renamed file into its final destination: {output_path} and temporary file deleted.")
        #final_tmp.unlink(missing_ok=True)


def copy_on_system(temp_file, final_file):
    try:
        tmp_onsystem = final_file.with_suffix(final_file.suffix + ".tmp")
        shutil.copy2(temp_file, tmp_onsystem)
        _atomic_replace(tmp_onsystem, final_file, retries=10)
        #os.replace(tmp_onsystem, final_file)
    except Exception as e:
        logger.error(f"Failed moving {temp_file} to {final_file}: {e}")
        open(final_file, "w").close()
    finally:
        #tmp_onsystem.unlink(missing_ok=True)
        temp_file.unlink(missing_ok=True)