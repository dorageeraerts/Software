
from pathlib import Path
import argparse
from glob import glob

import os
import logging
import json
import time
import re

import tempfile
import shutil
from multiprocessing import Pool

logger = logging.getLogger(__name__)


def process_run(type, run, rawfile_path, output_filename, n_events, rows_to_combine):
    """
    Process a single run with its subruns.
    
    Args:
        type: Data type (ADC or PIEDISTALLI)
        run: Run number
        input_optional_files: List of [slowcontrol, log, conteggi] file paths
        output_filename: Output file path
        n_events: Number of events to process
        rows_to_combine: Number of rows to combine as a single event
        log_on_console: Whether to log to console (True/False as string)
        status_file: Path to status file (json)
        num_threads: Number of threads to use for processing (default=1)
    
    Returns:
        result dict with status information
    """
    
    output_filename = Path(output_filename)
    json_filename = output_filename.with_suffix('.json')
    #input_slowcontrol = Path(input_optional_files[0])
    #input_log = Path(input_optional_files[1])
    #input_conteggi = Path(input_optional_files[2])
    #decompress_path = str(file_path.parent).replace("PARSED", "DECOMPRESSED")

    # setup status dictionary
    result = {
        "run": run,
        "subruns": [],
        "status": "ok",
        "quality_checks": {
            "has_less_nevents": None,
            "has_mismatches": None,
            "mismatch_counter": None,
            "has_unrecoverable_errors": None,
            "unrecoverable_mismatch_counter": None,
            "errors": None
        },
        "runtime": None
    }

    if type == "ADC":
        n_evts = 40000 
        type_raw = "slave" 
    else:
        n_evts = 50000
        type_raw = "ped"
    
    raw_filename_root = f"{rawfile_path}/{type_raw}Data_evts{n_evts}_run{run}" 
    logger.debug(f"Searching for subruns of {raw_filename_root}")
    subrun_list_notsorted = glob(f"{raw_filename_root}_*.gz")
    
    if len(subrun_list_notsorted) == 0:
        error = f"No subrun files found for run {run} and type {type} in path {rawfile_path}"
        logger.error(error)
        output_filename.touch(exist_ok=True)
        result["status"] = "failed"
        result["error"] = error
        with open(json_filename, 'w') as f:
            json.dump(result, f)
        return result
        exit(0)
    
    # sorting by subrun
    subrun_list = sorted(subrun_list_notsorted, key=lambda f: int(re.findall(r'\d+', f)[-1]))
    logger.info(f"List of files to decompress and parse:{subrun_list}")

    if type == "ADC" and len(subrun_list) != 4:
        logger.warning(f"Unconsistent number of subrun for type ADC. Number of subrun found: {len(subrun_list)}")
    if type == "PIEDISTALLI" and len(subrun_list) != 5:
        logger.warning(f"Unconsistent number of subrun for type PIEDISTALLI. Number of subrun found: {len(subrun_list)}")

    start_global = time.time()
    logger.debug(f"Starting... [{start_global}]")
    temp_filelist = []
    for subrun in subrun_list:  
        subrun_number = str(Path(subrun).stem).split("_")[-1].replace('sr', '')
        start_subrun = time.time()
        subrun_dict = {
            "id": 'sr' + str(subrun_number),
            "decompression": None,
            "parsing": "ok",
            "error": [],
            "runtime": None
        }
        # decompress the file in a temporary folder
        print(f"Decompressing file: {subrun}")
        decompressed_filename, error = file_handler.decompress(subrun)
        if error is not None:
            subrun_dict["parsing"] = "failed"
            subrun_dict["error"].append(error)
            result["status"] = "failed"
            subrun_dict["decompression"] = "failed"
            pass
        else:
            subrun_dict["decompression"] = "ok"
        logger.info(f"Subrun: {subrun_number}, Decompressed filename: {decompressed_filename} ")

        ctrl = 0
        error = None
        try:
            if "NERO" in rawfile_path.split("/"):
                logger.info("Parsing data from NERO")
                ctrl, tmp_filename = parserer.parser_new(decompressed_filename, 40000 if type == "ADC" else 50000, run, subrun_number, n_events, rows_to_combine, color="NERO")
                logger.info(f"tmp_filename: {tmp_filename}, ctrl: {ctrl}")
            elif "ROSSO" in rawfile_path.split("/"):
                logger.info("Parsing data from ROSSO")
                ctrl, tmp_filename = parserer.parser_new(decompressed_filename, 40000 if type == "ADC" else 50000, run, subrun_number, n_events, rows_to_combine, color="ROSSO")
                logger.info(f"tmp_filename: {tmp_filename}, ctrl: {ctrl}")
            elif "BLU" in rawfile_path.split("/"): 
                logger.info("Parsing data from BLU")
                ctrl, tmp_filename = parserer.parser_new(decompressed_filename, 40000 if type == "ADC" else 50000, run, subrun_number, n_events, rows_to_combine, color="BLU")
                logger.info(f"tmp_filename: {tmp_filename}, ctrl: {ctrl}")
            else:
                NameError("Data path must contain the color of the hodoscope. Retry.")
                exit(1) 
            temp_filelist.append(tmp_filename)
        except:
            error = f"Parsing failed unexpectedly, ctrl = {ctrl}. Decompressed file: {decompressed_filename}."
            logger.error(error)
            subrun_dict["parsing"] = "failed"
            subrun_dict["error"].append(error)
            result["status"] = "failed"
        finally:
            print(f"Subrun {subrun_number} processed in {round(time.time() - start_subrun, 2)} seconds with ctrl = {ctrl}")
            if decompressed_filename and os.path.exists(decompressed_filename):
                os.remove(decompressed_filename)
                print(f"Removing temporary file: {decompressed_filename}")
            subrun_dict["runtime"] = round(time.time() - start_subrun, 2)
            result["subruns"].append(subrun_dict)
            if ctrl == 1:
                error = "Parsing failed with controller = 1."
                logger.error(error)
                subrun_dict["parsing"] = "failed"
                subrun_dict["error"].append(error)
                result["status"] = "failed"

    # merge all temp_files (of each subrun) in one final file
    logger.info(f"Merging temporary files: {temp_filelist} into final file: {output_filename}")
    with file_handler.temp_to_output(output_filename) as tmp_path:
        with open(tmp_path, "w") as tf:
            for f in temp_filelist:
                with open(f) as fin:
                    tf.writelines(fin)
                os.remove(f)
        # perform quality control on the merged file
        with open(tmp_path, "r") as tf:
            # check that number of lines is 40 000 for ADC and 50 000 for PIEDISTALLI   
            try:
                n_lines = len(tf.readlines())
                if type == "ADC":
                    assert n_lines == 40000, f"Number of lines in the parsed file is not 40 000 but {n_lines}"
                else:
                    assert n_lines == 50000, f"Number of lines in the parsed file is not 50 000 but {n_lines}"
                tf.seek(0)  # reset file pointer to the beginning of the file for the next check
            except AssertionError as e:
                logger.error(str(e))
                result["status"] = "failed"
                result["quality_checks"]["has_less_nevents"] = True
                raise # error is handles by temp_to_output context manager that will ensure that an empty file is written in output and the job is not failed.

                
            check_result = parserer.check_subevent_number_spacing(tf, line_to_print=105)
            result["quality_checks"]["has_mismatches"] = check_result.has_mismatches
            result["quality_checks"]["mismatch_counter"] = check_result.mismatches_counter
            result["quality_checks"]["has_unrecoverable_errors"] = check_result.has_unrecoverable_errors
            result["quality_checks"]["unrecoverable_mismatch_counter"] = check_result.unrecoverable_mismatches_counter
            result["quality_checks"]["errors"] = check_result.errors
            result["quality_checks"]["mismatches_per_line"] = check_result.mismatches_per_line
            result["quality_checks"]["unrecoverable_mismatches_details"] = check_result.unrecoverable_mismatches_details
            if check_result.has_unrecoverable_errors == True:
                logger.error(f"Parser file didn't satisfy minimal check on relative event number and reported the following mismatches: \n {check_result.errors}")
                result["status"] = "failed"
                #raise # error is handles by temp_to_output context manager that will ensure that an empty file is written in output and the job is not failed.
            elif check_result.has_unrecoverable_errors == False and check_result.has_mismatches == True :
                logger.warning(f"Parser file didn't satisfy check on relative event number but the errors are recoverable.")
            else:
                logger.info("Parser file satisfied check on relative event number.")


    end_global = time.time()
    logger.debug(f"End [{end_global}]")
    result["runtime"] = round(end_global - start_global, 2)

    # Save dictionary
    with file_handler.temp_to_output(json_filename) as tmp:
        with open(tmp, "w") as f:
            json.dump(result, f)
    


    return result


if __name__ == "__main__":

    
    parser = argparse.ArgumentParser(description='Input configuration for cutting script')
    parser.add_argument("-t", "--type", dest="type", required=True,
                        help="Choose pedestal or slave data: [ADC | PIEDISTALLI]")
    parser.add_argument("-r", "--run_list", dest="run_list", required=True, nargs= "+",
                        help="List of runs to process (e.g. [2500 2501 2502] )")
    #parser.add_argument("-i", "--input_optional_files", dest="input_optional_files", nargs="+", required=True,
    #                    help="Input slowcontrol, log and conteggi files")
    parser.add_argument("-i", "--raw_data_path", dest="raw_data_path", required=True,
                        help="Path to raw data files")
    parser.add_argument("-l", "--log_on_console", dest="log_on_console", required=True,
                        help="If true logs are printed on terminal, if False they are printed on a file.")
    parser.add_argument("-ow", "--overwrite_outputs", dest="overwrite_outputs", required=True, 
                        help="If true, existing output files will be overwritten. If False, existing output files will be kept and the corresponding runs will be skipped.")
    parser.add_argument("-v", "--verbose", dest="verbose", required=False, default="info",
                        help="Logging level: debug/info/warning/error/critical (default = info)")
    #parser.add_argument("-j", "--status_file", dest="status_file", required=True,
    #                    help="Status file for this rule: json file for summury analysis")
    parser.add_argument("-o", "--output_filename", dest="output_filename", required=True,
                        help="Name of the output")
    parser.add_argument("-e", "--n_events", dest="n_events", type=int, default=10000,
                        help="Number of events in the file (defalut = 10 000)")
    parser.add_argument("-n", "--rows_to_combine", dest="rows_to_combine", type=int, default=16,
                        help="Number of raw in the file to be combined under a unique events (default = 16, number of slave board) ")
    parser.add_argument("-th", "--num_threads", dest="num_threads", type=int, default=1,
                        help="Number of threads/cores to use for processing (default = 1)")
    args = parser.parse_args()

    log_on_console = args.log_on_console
    overwrite_outputs = args.overwrite_outputs
    verbose = args.verbose
    batch_idx = int(args.output_filename.split("_batch")[-1].split(".")[0])  # Extract batch index from output filename
    output_filename_root = str(Path(args.output_filename).name).split("_batch")[0]
    output_filename_path = Path(args.output_filename).parent

    # setup logging
    if log_on_console == 'True':
        tmp_log = None
        logging.basicConfig(
            level=getattr(logging, verbose.upper()),
            format="%(asctime)s [%(levelname)s] %(message)s"
        )
    else:
        log_file = "logs/PARSED/" + Path(args.output_filename).with_suffix(".log").name
        Path(log_file).parent.mkdir(parents=True, exist_ok=True)
        open(log_file, "w").close()
        logging.basicConfig(
            filename=log_file,
            level=getattr(logging, verbose.upper()),
            format="%(asctime)s [%(levelname)s] %(message)s"
        )
    # logging of file_handler are not printed for some reasons.
    from muraves_lib import file_handler
    import parserer
    # Create a wrapper function to pass run-specific output filename
    def process_run_with_output(run):
        run_output_filename = output_filename_path / Path(f"{output_filename_root}_run{run}.txt")

        if overwrite_outputs == 'False' and (run_output_filename.exists() and run_output_filename.with_suffix('.json').exists()):
            logger.info(f"Output file {run_output_filename} already exists and overwrite_outputs is set to False. Skipping run {run}.")
        elif overwrite_outputs == 'True' or (not run_output_filename.exists() or not run_output_filename.with_suffix('.json').exists()):
                logger.info(f"Processing run {run}...")
                process_run(
                    type=args.type,
                    run=run,
                    rawfile_path=args.raw_data_path,
                    output_filename=str(run_output_filename),
                    n_events=args.n_events,
                    rows_to_combine=args.rows_to_combine,
                )
        else:
            logger.error(f"Invalid value for overwrite_outputs: {overwrite_outputs}. It should be 'True' or 'False'.")

        return run_output_filename

    # Processa tutti i run del batch in parallelo usando il numero di threads fornito da Snakemake
    run_output_filename_list = Pool(args.num_threads).map(process_run_with_output, args.run_list)

    with file_handler.temp_to_output(args.output_filename) as tmp:
        with open(tmp, "a") as f:
            for filename in run_output_filename_list:
                f.write(f"{filename}\n")
    print(f'Batch {batch_idx} completato con {args.num_threads} thread: {args.run_list}')
