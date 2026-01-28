
from pathlib import Path
import argparse
from glob import glob
import parserer
import os
import logging
import json
import time
import re
from muraves_lib import file_handler
import tempfile
import shutil

parser = argparse.ArgumentParser(description='Input configuration for cutting script')
parser.add_argument("-t", "--type", dest="type", required=True,
                    help="Choose pedestal or slave data: [ADC | PIEDISTALLI]")
parser.add_argument("-r", "--run", dest="run", required=True,
                    help="Choose the run number")
parser.add_argument("-i", "--input_optional_files", dest="input_optional_files", nargs="+", required=True,
                    help="Input slowcontrol, log and conteggi files")
parser.add_argument("-l", "--log_on_console", dest="log_on_console", required=True,
                    help="If true logs are printed on terminal, if False they are printed on a file.")
parser.add_argument("-j", "--status_file", dest="status_file", required=True,
                    help="Status file for this rule: json file for summury analysis")
parser.add_argument("-o", "--output_filename", dest="output_filename", required=True,
                    help="Name of the output")
parser.add_argument("-e", "--n_events", dest="n_events", type=int, default=10000,
                    help="Number of events in the file (defalut = 10 000)")
parser.add_argument("-n", "--rows_to_combine", dest="rows_to_combine", type=int, default=16,
                    help="Number of raw in the file to be combined under a unique events (default = 16, number of slave board) ")
args = parser.parse_args()


type = args.type
run = args.run
input_slowcontrol = Path(args.input_optional_files[0])
input_log = Path(args.input_optional_files[1])
input_conteggi = Path(args.input_optional_files[2])
rawfile_path = str(input_log.parent)
output_filename = args.output_filename
file_path = Path(output_filename)
events_number = args.n_events
rows_to_combine = args.rows_to_combine
log_on_console= args.log_on_console
status_file = args.status_file

decompress_path = str(file_path.parent).replace("PARSED", "DECOMPRESSED")

# setup logging
if log_on_console == 'True':
    logging.basicConfig(
        #filename=log_file,
        level=logging.DEBUG,
        format="%(asctime)s [%(levelname)s] %(message)s"
    )
else:
    log_file = output_filename.replace(".txt", ".log")
    Path(log_file).unlink(missing_ok=True) # se esiste lo svuota, altrimenti più tardi viene creato
    log_dir = Path(log_file).parent
    log_dir.mkdir(exist_ok=True)
    logging.basicConfig(
        filename=log_file,
        level=logging.DEBUG,
        format="%(asctime)s [%(levelname)s] %(message)s"
    )

    
#print("Log filename", log_file)
# setup status dictionary
result = {
    "run": run,
    "subruns": [],
    "status": "ok",
    "runtime": None
}


if type=="ADC":
    n_evts = 40000 
    type_raw = "slave" 
else:
    n_evts = 50000
    type_raw = "ped"
filename_root = f"{rawfile_path}/{type_raw}Data_evts{n_evts}_run{run}" 
logging.debug(f"Searching for subruns of {filename_root}")
subrun_list_notsorted = glob(f"{filename_root}_*.gz")
if len(subrun_list_notsorted) == 0:
    error = f"No subrun files found for run {run} and type {type} in path {rawfile_path}"
    logging.error(error)
    file_path.touch(exist_ok=True)
    result["status"] = "failed"
    result["error"] = error
    with open(status_file, 'w') as f:
        json.dump(result, f)
    exit(0)
# sorting by subrun
subrun_list = sorted(subrun_list_notsorted, key=lambda f: int(re.findall(r'\d+', f)[-1]))
logging.info(f"List of files to decompress and parse:{subrun_list}")

if type=="ADC" and len(subrun_list)!=4:
    logging.warning(f"Unconsistent number of subrun for type ADC. Number of subrun found: {len(subrun_list)}")
if type=="PIEDISTALLI" and len(subrun_list)!=5:
    logging.warning(f"Unconsistent number of subrun for type PIEDISTALLI. Number of subrun found: {len(subrun_list)}")


start_global = time.time()
logging.debug(f"Starting... [{start_global}]")
temp_filelist = []
for subrun in subrun_list:  
    start_subrun = time.time()
    subrun_dict = {
    "id": None,
    "status": "ok",
    "error": [],
    "runtime": None
    }
    # decompress the file in a temporary folder which is not in /data (alias /pnfs)
    decompressed_filename, error = file_handler.decompress(subrun, destination_path=decompress_path)
    subrun_dict["id"] = decompressed_filename.split("_")[-1]
    if error is not None:
        subrun_dict["status"] = "failed"
        subrun_dict["error"].append(error)
        result["status"] = "failed"
    logging.info(f"parsing  {decompressed_filename} ")
    ctrl = 0
    try:
        if "NERO" in rawfile_path.split("/"):
            logging.info("Parsing data from NERO")
            #ctrl = parserer.parser_nero(decompressed_filename, output_filename, events_number, rows_to_combine )
            ctrl, tmp_filename = parserer.parser_nero_new(decompressed_filename, 'placeholder', events_number, rows_to_combine )
            logging.info(f"tmp_filename: {tmp_filename}, ctrl: {ctrl}")
        elif "ROSSO" in rawfile_path.split("/"):
            logging.info("Parsing data from ROSSO")
            ctrl = parserer.parser_rosso(decompressed_filename, output_filename, events_number, rows_to_combine )
        elif "BLU" in rawfile_path.split("/"): 
            logging.info("Parsing data from BLU")
            ctrl = parserer.parser_blu(decompressed_filename, output_filename, events_number, rows_to_combine )
        else:
            NameError("Data path must contain the color of the hodoscope. Retry.")
            exit(1) 
        temp_filelist.append(tmp_filename)
  
    except:
        error = f"Parsing failed unexpectedly. ctrl = {ctrl}"
        logging.error(error)
        file_path.touch(exist_ok=True)
        subrun_dict["status"] = "failed"
        subrun_dict["error"].append(error)
        result["status"] = "failed"
    finally:
        os.remove(decompressed_filename)
        subrun_dict["runtime"] = round(time.time() - start_subrun, 2)
        result["subruns"].append(subrun_dict)
        if ctrl == 1:
            error = "Parsing failed with controller = 1."
            logging.error(error)
            file_path.touch(exist_ok=True)
            subrun_dict["status"] = "failed"
            subrun_dict["error"].append(error)
            result["status"] = "failed"

# unisci tutti i temp_files in un file finale
try:
    # create temp file in the same directory to ensure atomic rename
    with tempfile.NamedTemporaryFile(mode='w', delete=False, dir=file_path.parent) as tf:
        tmp_name = tf.name
        for f in temp_filelist:
            with open(f) as fin:
                tf.writelines(fin)
            os.remove(f)
    # atomic move to final destination
    Path(tmp_name).replace(file_path)
    logging.debug(f"Successfully created total parsed file: {output_filename}")
except Exception as e:
    # if something goes wrong, remove the temp file
    if 'tmp_name' in locals() and Path(tmp_name).exists():
        Path(tmp_name).unlink()
    err = f'Failing creating total parsed file'
    logging.error(f"{err} for {output_filename}: {e}")
    file_path = Path(output_filename)
    file_path.touch(exist_ok=True)


end_global = time.time()
logging.debug(f"End [{end_global}]")
result["runtime"] = round(end_global- start_global, 2)

# Quality control on Parsed file:
# X if status is failed -> create an empty file
# V Relative ev. number must be the same 16 times in the same row
# 

try:
    rel_ev_num_check_errors = parserer.check_event_number_spacing(output_filename, 39, 0)
    if len(rel_ev_num_check_errors) > 0:
        logging.error[f"Parser file didn't satisfy minimal check on relative event number and reported the following mismatches: \n {rel_ev_num_check_errors}"]
except:
    if file_path.exists():
        file_path.unlink()     # delete the file
        file_path.touch() 


# Read SLOWCONTROL and store relevant info
if not input_slowcontrol.is_file():
    logging.warning(f"Slow control file not found: {input_slowcontrol}. Slowcontrol dictionary will be empty.")
else:
    decompressed_slowfile, error = file_handler.decompress(input_slowcontrol, decompress_path)
    slowcontrol_dict = file_handler.parse_slow_control(decompressed_slowfile, run)
    result["slowcontrol"] = slowcontrol_dict
    #os.remove(decompressed_slowfile)

# Read LOG and store relevand info
if not input_log.is_file():
    logging.warning(f"Input log file not found: {input_log}. Log dictionary will be empty.")
else:
    decompressed_logfile, error = file_handler.decompress(input_log, decompress_path)
    log_dict = file_handler.parse_log_file(decompressed_logfile, run)
    result['log'] = log_dict

# Read CONTEGGI
if not input_conteggi.is_file():
    logging.warning(f"Input conteggi file not found: {input_conteggi}. Log dictionary will be empty.")
else:
    decompressed_conteggifile, _ = file_handler.decompress(input_conteggi, decompress_path)
    conteggi_dict = file_handler.parse_conteggi(decompressed_conteggifile)
    result['conteggi'] = conteggi_dict
    
# NB: Cannot delete decompressed LOG, CONTEGGI, SLOW because when running in parallel could rise problem as the file is the same for ADC and PIEDISTALLI
# Save dictionary
with open(status_file, 'w') as f:
    json.dump(result, f)


