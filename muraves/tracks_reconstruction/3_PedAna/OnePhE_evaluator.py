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

Description = ' This code takes as input the pedestal run and computes the pedestal position and 1phe conversion.'
parser = argp.ArgumentParser(description = Description)
parser.add_argument('-data', '--ped_file', dest = "ped_file", required = True, help = 'This file cointains pedestal data')
parser.add_argument('-r','--run', dest = "run", required=True, help = 'run number')
#parser.add_argument("-info", "--info", dest = "info", nargs="+", required=True,
#                    help="Input number of runs and delta around the central run")
parser.add_argument("-o", "--output_list", dest="output_list", nargs = '+', required=True,
                    help="List of outputs")
args = parser.parse_args()

##################################

run = args.run
inputfile = args.ped_file
output_list = args.output_list
nBoards = len(output_list)
nChannels = 32
print("Started the script.")
root_path = str(Path(output_list[0]).parent)


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

"""
### Useless code for now ###
# It just saves the TSpectrum output in a ROOT file without any processing, this is anyway done below.
# The only difference is that here also the normal histograms shows TSpectrum peaks.
# Kept for reference.
             
# ===============================
# Salvataggio ROOT senza elaborazioni
# ===============================
if not Path("tspectrum_output.root").exists():
    print("Saving ROOT file...")
    outfile = ROOT.TFile(f"tspectrum_output.root", "RECREATE")
    canvas_dir = outfile.mkdir("canvases")
    tree = ROOT.TTree("PeakTree", "TSpectrum results")

    # Branches (come nel tuo codice)
    run_arr     = array('i', [0])
    board_arr   = array('i', [0])
    channel_arr = array('i', [0])
    npeaks_arr  = array('i', [0])
    pedestal_arr = array('d', [0])
    onephe_err   = array('d', [0])
    xpeaks = ROOT.std.vector('double')()
    ypeaks = ROOT.std.vector('double')()

    tree.Branch("run", run_arr, "run/I")
    tree.Branch("board", board_arr, "board/I")
    tree.Branch("channel", channel_arr, "channel/I")
    tree.Branch("npeaks", npeaks_arr, "npeaks/I")
    tree.Branch("xpeaks", xpeaks)
    tree.Branch("ypeaks", ypeaks)
    tree.Branch("pedestal", pedestal_arr, "pedestal/D")
    tree.Branch("onephe", onephe_err, "onephe/D")

    # Loop su board e canali
    for board in range(nBoards):
        for i in range(8):
            for j in range(4):
                histo = TH1F(f"h_board{board}_ch{4*i+j}", "PED distribution", 300, 1450, 1750)
                for val in AllRuns[0][board][4*i+j]:
                    histo.Fill(val)

                histo2 = TH1F(f"h_log_board{board}_ch{4*i+j}", "PED distribution log", 300, 1450, 1750)
                for bin in range(histo.GetXaxis().GetNbins()):
                    if histo.GetBinContent(bin) != 0:
                        histo2.SetBinContent(bin, mt.log10(histo.GetBinContent(bin)))

                s = TSpectrum()
                nfoundPeacks = s.Search(histo2, 0, "goff", 0.00001)

                xPeacks = s.GetPositionX()
                yPeacks = s.GetPositionY()
                xpeaks.clear()
                ypeaks.clear()
                for ip in range(nfoundPeacks):
                    xpeaks.push_back(xPeacks[ip])
                    ypeaks.push_back(yPeacks[ip])

                # Fill the tree
                run_arr[0] = int(args.run)
                board_arr[0] = board
                channel_arr[0] = 4*i + j
                npeaks_arr[0] = nfoundPeacks
                pedestal_arr[0] = float(xPeacks[0]) if nfoundPeacks > 0 else -1
                onephe_err[0] = -1
                tree.Fill()

                # =============================
                # Canvas grafici
                # =============================
                canvas_dir.cd()

                # Canvas normale
                c1 = TCanvas(f"c_board{board}_ch{4*i+j}_norm", "", 800, 600)
                histo.Draw("hist")
                if nfoundPeacks > 0:
                    pm1 = ROOT.TPolyMarker()
                    for ip in range(nfoundPeacks):
                        bin_index = histo.FindBin(xPeacks[ip])
                        y_value = histo.GetBinContent(bin_index)
                        pm1.SetPoint(ip, xPeacks[ip], y_value)
                    pm1.SetMarkerStyle(23)
                    pm1.SetMarkerColor(ROOT.kRed)
                    pm1.SetMarkerSize(0.9)
                    histo.GetListOfFunctions().Add(pm1)
                c1.Write()
                c1.Close()

                # Canvas log
                c2 = TCanvas(f"c_board{board}_ch{4*i+j}_log", "", 800, 600)
                histo2.Draw("hist")
                if nfoundPeacks > 0:
                    pm2 = ROOT.TPolyMarker()
                    for ip in range(nfoundPeacks):
                        bin_index = histo.FindBin(xPeacks[ip])
                        y_value = histo.GetBinContent(bin_index)
                        pm2.SetPoint(ip, xPeacks[ip], mt.log10(y_value))
                    pm2.SetMarkerStyle(23)
                    pm2.SetMarkerColor(ROOT.kRed)
                    pm2.SetMarkerSize(0.9)
                    histo2.GetListOfFunctions().Add(pm2)
                c2.Write()
                c2.Close()

    outfile.Write()
    outfile.Close()
    print("ROOT file saved!")

"""
    


outfile_final = ROOT.TFile(f"{root_path}/pedestal_output.root", "RECREATE")
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

print("Computing fallsback values...")
for nBoard, filename in enumerate(output_list):
    file=open(filename,"w")
    file.write("ch \t ped \t 1pe \n")
    for ch in range(nChannels):
        for r in range(len(AllRuns)):
            if  OnePhes[nBoard][ch][r] == "no 2nd peack":
                 file.write(str(ch)+"\t"+str(int(PedValues[nBoard][ch][r]))+ "\t" + str(1000)+ "\t")
            else:
#                if int(OnePhes[nBoard][ch][r]) > 38 or int(OnePhes[nBoard][ch][r]) < 25: #nero
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

