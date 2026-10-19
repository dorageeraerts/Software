from ROOT import TH1F
from ROOT import TFile
from ROOT import TSpectrum
from ROOT import TCanvas
import ROOT
ROOT.gROOT.SetBatch(True)

#import sys
import matplotlib.pyplot as plt
import glob
import numpy as  np
import math as mt
#import subprocess
import argparse as argp
import logging
import uproot 
import pandas as pd
import time
from array import array
from pathlib import Path
from muraves_lib import file_handler
from multiprocessing import Pool
logger = logging.getLogger(__name__)

def FillVector(input_file, nBoards):

    print("Reading the input tree...")
    start_time = time.time()
    with uproot.open(input_file) as f:
        tree = f[f.keys()[0]]
        df = tree.arrays(library="pd")

    # Prepare output structure
    boardALL = [[] for _ in range(nBoards)]

    # For each board index
    for i in range(nBoards):

        # Select rows where scheda == i
        df_i = df[df["scheda"] == i]

        # Extract the 32 ADC columns
        adc_cols = [f"adc_{k}" for k in range(32)]

        # Convert each row to a list of 32 values
        boardALL[i] = df_i[adc_cols].values.tolist()

    print("End.")
    print("Run Time: ", time.time() - start_time )
    return boardALL


def analyse_pedestal_file(tmp_filename, inputfile, run,  nBoards, nInfoBoard, nChannels):


    print("Started the script.")
    #output_path = Path(output_filename).parent / Path(run)
    ###########################################################################################

    AllRuns = []

    ################################### Fill Arrays  ###########################################

    boardALL = FillVector(inputfile, nBoards)
    boardT = []
    for nBoard in range(nBoards):
        boardT.append(np.array(boardALL[nBoard]).transpose())
    AllRuns.append(boardT)
    print("LEN(boardALL):", len(boardALL))

    for b in range(nBoards):
        print(f"Board {b}: entries = {len(boardALL[b])}")
    del boardALL
    del boardT
    #import pdb; pdb.set_trace()
    outfile_final = ROOT.TFile(str(tmp_filename), "recreate")
    tree_final = ROOT.TTree("PedTree", "Global pedestal analysis results")
    run_arr     = array('i', [0])
    board_arr   = array('i', [0])
    channel_arr = array('i', [0])
    npeaks_arr  = array('i', [0])
    pedestal_arr = array('d', [0])
    onephe_arr   = array('d', [0])
    second_peak = array('d', [0])

    tree_final.Branch("run", run_arr, "run/I")
    tree_final.Branch("board", board_arr, "board/I")
    tree_final.Branch("channel", channel_arr, "channel/I")
    tree_final.Branch("npeaks", npeaks_arr, "npeaks/I")
    tree_final.Branch("pedestal", pedestal_arr, "pedestal/D")
    tree_final.Branch("onephe", onephe_arr, "onephe/D")
    tree_final.Branch("second_peak", second_peak, "second_peak/D")

    OnePhes = []
    PedValues = []

    for board in range(nBoards):
        print('Board: ',board)
        OnePhes_singleBoard = []
        PedValues_singleBoard = []
        for i in range(0,8):
            for j in range(0,4):

                OnePhes_channel = []
                PedValues_channel = []
                for r in  range(len(AllRuns)):
                    channel = 4*i + j
                    ###### Start filling tree 
                    run_arr[0]     = int(run)
                    board_arr[0]   = board
                    channel_arr[0] = channel

                    ################### Peacks searching ###########################
                    histo = TH1F("Ped_run_"+run+"_board_"+str(board)+"_ch_"+str(4*i+j),"PED distribution",300,1450,1750)
                    histo2 = TH1F("Ped_run_"+run+"_board_"+str(board)+"_ch_"+str(4*i+j)+"_logscale","PED distribution",300,1450,1750)
                    #### Using log scale to enanhance peacks 
                    #print(f"Run={run}, board={board}, channel={4*i+j}")
                    #print(f"Shape of AllRuns[run][board]: {np.shape(AllRuns[run][board])}")

                    ################## Histograms filling #####################################
                    for index in range(len(AllRuns[r][board][4*i+j])):
                        histo.Fill(int(AllRuns[r][board][4*i+j][index]))

                    for bin in range(histo.GetXaxis().GetNbins()):
                        if histo.GetBinContent(bin)!=0:
                            histo2.SetBinContent(bin,mt.log10(histo.GetBinContent(bin)))

                    ##################### Searching peacks with TSpectrum ##################
                    s = TSpectrum()
                    nfoundPeacks = s.Search(histo2,0,"",0.00001)

                    #histo2.Write()
                    xPeacks= s.GetPositionX()
                    xPositions=[]
                    for x in range(nfoundPeacks):
                        xPositions.append(xPeacks[x])
                    xPositions.sort()
                    yPeacks = []

                    yPeacks_real = []
                    for ip in range(nfoundPeacks):
                        bin_index = histo.FindBin(xPeacks[ip])
                        yPeacks_real.append(histo.GetBinContent(bin_index))

                    ############################# Evaluate 1phe ############################
                    for  h in range(nfoundPeacks):
                        yPeacks.append(histo2.GetBinContent(histo.FindBin(xPositions[h])))
                    MaxIndex  = yPeacks.index(max(yPeacks))
                    SecondMaxPeacks = MaxIndex+1
                    try:
                        Onephe =  xPositions[SecondMaxPeacks]-xPositions[MaxIndex]
                        if Onephe <= 20 and  ( xPositions[SecondMaxPeacks+1] -  xPositions[SecondMaxPeacks]) <15 :
                            Onephe =  xPositions[SecondMaxPeacks+1]-xPositions[MaxIndex]
                        OnePhe_num = Onephe
                    except:
                        Onephe = "no 2nd peack"
                        OnePhe_num = -1

                    OnePhes_channel.append(Onephe)
                    PedValues_channel.append(xPositions[MaxIndex])

                    ###### Keep filling tree
                    pedestal_arr[0] = xPositions[MaxIndex]
                    onephe_arr[0]   = OnePhe_num
                    second_peak[0] = xPositions[MaxIndex] + OnePhe_num if OnePhe_num != -1 else -1

                    tree_final.Fill()
                    ###############################################################

                ############### Storing values  ########################################
                OnePhes_singleBoard.append(OnePhes_channel)
                PedValues_singleBoard.append(PedValues_channel)
                del  OnePhes_channel
                del PedValues_channel
        OnePhes.append(OnePhes_singleBoard)
        PedValues.append(PedValues_singleBoard)
    ################# 1phe  writings  ##############################
    
    outfile_final.Write()
    outfile_final.Close()
    print("Final ROOT file saved!")
    return  OnePhes, PedValues, AllRuns


def write_results(nBoard, filename, nChannels, OnePhes, PedValues, AllRuns):
    print("Computing fallsback values...")
    file=open(filename,"w")
    file.write("ch \t ped \t 1pe \n")
    for ch in range(nChannels):
        for r in range(len(AllRuns)):
            if  OnePhes[nBoard][ch][r] == "no 2nd peack":
                 file.write(str(ch)+"\t"+str(int(PedValues[nBoard][ch][r]))+ "\t" + str(1000)+ "\t")
            else:
    #               if int(OnePhes[nBoard][ch][r]) > 38 or int(OnePhes[nBoard][ch][r]) < 25: #nero
                if  nBoard !=5 and (int(OnePhes[nBoard][ch][r]) > 38 or int(OnePhes[nBoard][ch][r]) < 20):  
                    if ch!=0 and OnePhes[nBoard][ch-1][r] != "no 2nd peack":
                        file.write(str(ch)+"\t"+str(int(PedValues[nBoard][ch][r]))+ "\t" + str(int(OnePhes[nBoard][ch-1][r]))+ "\t 1" )
                        OnePhes[nBoard][ch][r] = OnePhes[nBoard][ch-1][r]
                    else:
                        file.write(str(ch)+"\t"+str(int(PedValues[nBoard][ch][r]))+ "\t" + str(28)+ "\t 0")
                        OnePhes[nBoard][ch][r] = 30
                elif nBoard==5 and (int(OnePhes[nBoard][ch][r]) > 28 or int(OnePhes[nBoard][ch][r]) < 10):
                    if ch!=0 and OnePhes[nBoard][ch-1][r] != "no 2nd peack":
                        file.write(str(ch)+"\t"+str(int(PedValues[nBoard][ch][r]))+ "\t" + str(int(OnePhes[nBoard][ch-1][r]))+ "\t 1")
                        OnePhes[nBoard][ch][r] = OnePhes[nBoard][ch-1][r]
                    else:
                        file.write(str(ch)+"\t"+str(int(PedValues[nBoard][ch][r]))+ "\t" + str(18)+ "\t")
                        OnePhes[nBoard][ch][r] = 18
                else:
                    file.write(str(ch)+"\t"+str(int(PedValues[nBoard][ch][r]))+ "\t" + str(int(OnePhes[nBoard][ch][r]))+ "\t 0")
        file.write("\n")


if __name__ == "__main__":

    Description = ' This code takes as input the pedestal run and computes the pedestal position and 1phe conversion.'
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

    ##################################

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
        log_file = "logs/PEDESTAL/" + Path(args.output_filename).with_suffix(".log").name
        Path(log_file).parent.mkdir(parents=True, exist_ok=True)
        open(log_file, "w").close()
        logging.basicConfig(
            filename=log_file,
            level=getattr(logging, args.verbose.upper()),
            format="%(asctime)s [%(levelname)s] %(message)s"
        )

    with open(input_filename, "r") as f_input:
        file_list = [line.strip() for line in f_input]

    def process_run(prereco_root_file):
        runnumber = str(Path(prereco_root_file).stem).split("run")[-1]
        root_output_filename_per_run = str(Path(prereco_root_file.replace("PRERECONSTRUCTED", "PEDESTAL")).parent) + f"/{runnumber}"+ "/pedestal_analysis.root" 
        Path(root_output_filename_per_run).parent.mkdir(parents=True, exist_ok=True)
        outputs_list_per_run = [str(Path(root_output_filename_per_run).parent) + f"/pedestal_{n}" for n in range(nBoards)]
        should_process = (
            overwrite_outputs == 'True'
            or not Path(root_output_filename_per_run).exists()
            or any(not Path(p).exists() for p in outputs_list_per_run)
        )
        if not should_process:
            logger.info(f"Output files already exists and overwrite_outputs is set to False. Skipping run.")
        else:
            logger.info(f"Analysing file {prereco_root_file}...")
            try:
                with file_handler.temp_to_output(root_output_filename_per_run, f"Successfully analysed pedestal inputfile {prereco_root_file}") as tmp_path:                            
                     OnePhes, PedValues, AllRuns = analyse_pedestal_file(tmp_path, prereco_root_file, runnumber, nBoards, nInfoBoard, nChannels)
            except:
                print("[ERROR] Something went wrong while analysing!")
                raise
            try:
                for n, file in enumerate(outputs_list_per_run):
                    with file_handler.temp_to_output(file) as tmp_file:
                        write_results(n, tmp_file, nChannels,  OnePhes, PedValues, AllRuns)
            except:
                print("[ERROR] Something went wrong while writing results!")
                raise 
        return root_output_filename_per_run, runnumber
    

    results = Pool(args.num_threads).map(process_run, file_list)
    output_filename_per_run_list, run_list = zip(*results)
        

    with file_handler.temp_to_output(args.output_filename) as tmp:
        with open(tmp, "a") as f:
            for filename in output_filename_per_run_list:
                f.write(f"{filename}\n")
    print(f'Batch {batch_idx} completato con {args.num_threads} thread. \nRuns: {run_list}')





