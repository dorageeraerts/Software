import argparse as argp
import logging
from muraves_lib import file_handler
from pathlib import Path

import uproot
import numpy as np

from multiprocessing import Pool
logger = logging.getLogger(__name__)


def parse_strips(strip_str):
    """Convertit une chaîne comme '_15_14_' en liste d'entiers [15, 14]."""
    if strip_str.startswith('_') and strip_str.endswith('_'):
        return [int(x) for x in strip_str[1:-1].split('_') if x]
    return []

def parse_adc_line(line, nBoards, nInfoBoard, nChannels):
    """Parse une ligne ADC en extrayant strips et canaux."""
    data = line.split()
    boards = []
    for n in range(nBoards):
        start = n * nInfoBoard
        try:
            board = {
                "boardNumber": int(data[start + 2]),
                "temperature": int(data[start + 36]),
                "channels": [int(x) for x in data[start + 3 : start + 3 + nChannels]],
                "strips": parse_strips(data[start + 39])  # Champ FastOr $ex: '_15_14_'$
            }
            boards.append(board)
        except (IndexError, ValueError):
            continue  # Ignore les boards incomplets/malformés
    return boards

def make_tree(output_filename, datalines, nBoards, nInfoBoard, nChannels):
    """Crée un TTree avec uproot, incluant les strips déclenchés."""
    total_entries = len(datalines) * nBoards
    scheda = np.empty(total_entries, dtype=np.int32)
    temperature = np.empty_like(scheda)
    channels = np.empty((total_entries, nChannels), dtype=np.int32)
    strips = []  # Liste de listes (non vectorisable directement)

    idx = 0
    for line in datalines:
        boards = parse_adc_line(line, nBoards, nInfoBoard, nChannels)
        for board in boards:
            scheda[idx] = board["boardNumber"]
            temperature[idx] = board["temperature"]
            channels[idx] = board["channels"]
            strips.append(board["strips"])  # Stocke la liste des strips
            idx += 1

    # Redimensionne les tableaux
    scheda = scheda[:idx]
    temperature = temperature[:idx]
    channels = channels[:idx]

    # Convertit les strips en tableau de type "object" pour uproot
    strips_array = np.empty(len(strips), dtype=object)
    strips_array[:] = strips

    with uproot.recreate(output_filename) as f:
        f.mktree('Tree', {
            "scheda": scheda,
            "temperature": temperature,
            **{f"adc_{i}": channels[:, i] for i in range(nChannels)},
            "strips": strips_array  # Champ supplémentaire pour les strips
        })
    




if __name__ == "__main__":

    Description = ' This code takes as input the data from ADC and unpacks them to obtain a track collection '
    parser = argp.ArgumentParser(description = Description)
    parser.add_argument('-i', '--input_filename', dest="input_filename", required = True, help = 'This files contains the list of ADC files.')
    parser.add_argument("-l", "--log_on_console", dest="log_on_console", required=True,
                        help="If true logs are printed on terminal, if False they are printed on a file.")
    parser.add_argument("-ow", "--overwrite_outputs", dest="overwrite_outputs", required=True, 
                        help="If true, existing output files will be overwritten. If False, existing output files will be kept and the corresponding runs will be skipped.")
    parser.add_argument("-v", "--verbose", dest="verbose", required=False, default="info",
                        help="Logging level: debug/info/warning/error/critical (default = info)")
    parser.add_argument("-th", "--num_threads", dest="num_threads", type=int, default=1,
                        help="Number of threads/cores to use for processing (default = 1)")
    parser.add_argument("-info", "--info_board", dest="info_board", nargs="+", required=True,
                        help="Input nBoards, nInfoBoard and nChannels integer values")
    parser.add_argument("-o", "--output_filename", dest="output_filename", required=True,
                        help="Name of the output")
    args = parser.parse_args()
    
    input_filename = args.input_filename
    nBoards = int(args.info_board[0])
    nInfoBoard = int(args.info_board[1])
    nChannels= int(args.info_board[2])
    output_filename = args.output_filename
    batch_idx = int(args.output_filename.split("_batch")[-1].split(".")[0])
    overwrite_outputs = args.overwrite_outputs

    # setup logging per batch
    if args.log_on_console == 'True':
        logging.basicConfig(
            level=getattr(logging, args.verbose.upper()),
            format="%(asctime)s [%(levelname)s] %(message)s"
        )
    else:
        log_file = "logs/PRERECONSTRUCTED/" + Path(args.output_filename).with_suffix(".log").name
        Path(log_file).parent.mkdir(parents=True, exist_ok=True)
        open(log_file, "w").close()
        logging.basicConfig(
            filename=log_file,
            level=getattr(logging, args.verbose.upper()),
            format="%(asctime)s [%(levelname)s] %(message)s"
        )


    with open(input_filename, "r") as f_input:
        file_list = [line.strip() for line in f_input]

    def process_run(parsed_file):
        output_filename_per_run = Path(parsed_file.replace("PARSED", "PRERECONSTRUCTED").replace(".txt", ".root"))
        runnumber = str(output_filename_per_run.stem).split("run")[-1]
        datalines = []
        if overwrite_outputs == 'False' and output_filename_per_run.exists():
            logger.info(f"Output file {output_filename_per_run} already exists and overwrite_outputs is set to False. Skipping run.")
        elif overwrite_outputs == 'True' or not output_filename_per_run.exists():
            logger.info(f"Processing run {runnumber}...")
            try:
                with open(parsed_file, "r") as f_parsed:
                    datalines = f_parsed.read().splitlines()
            except:
                logger.error("File do not exist! It should always exist as the parsing stage. If it doesn't probably there was an error while parsing this batch of runs. ")          
            if len(datalines)==0:
                logger.info(f"Parsed file {parsed_file} is empty. Parsing failed or the run is missing. Creating an empty ouput for snakemake.")
                with file_handler.temp_to_output(output_filename_per_run) as tmp_path: 
                    pass
            else:
                with file_handler.temp_to_output(output_filename_per_run) as tmp_path:           
                    make_tree(tmp_path, datalines, nBoards, nInfoBoard, nChannels)
        return output_filename_per_run, runnumber
    

    results = Pool(args.num_threads).map(process_run, file_list)
    output_filename_per_run_list, run_list = zip(*results)
        

    with file_handler.temp_to_output(args.output_filename) as tmp:
        with open(tmp, "a") as f:
            for filename in output_filename_per_run_list:
                f.write(f"{filename}\n")
    print(f'Batch {batch_idx} completato con {args.num_threads} thread. \nRuns: {run_list}')
    #print(f"time of writing: {round(end_time - start_time, 2)}")
    


