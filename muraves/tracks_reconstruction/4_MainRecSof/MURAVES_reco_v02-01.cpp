#include <stdio.h>
#include <Python.h>
#include "pyhelper.hpp"
#include <string>
#include <sstream>
#include <stdlib.h>
#include <regex>
#include <iostream>
#include <fstream>
#include <time.h>
#include <ctime>
#include <numeric>
#include <algorithm>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TDatime.h"
#include "TMath.h"
#include "TObject.h"

////// MURAVES-FUNCTIONS LIBRARIES

#include "EvaluateAngularCoordinates.h"
#include "ClusterLists.h"
#include "ReadEvent.h"
#include "Tracking.h"
using namespace std;

int main(int argc, char* argv[]) {

  /////// welcome message ////////
  cout << " ~~~~~~~  Welcome to the MURAVES reconstruction 2.0 ~~~~~~~~  " << endl;
  cout << "                      .-----. "<< endl;
  cout << "             .----. .'       ' "<< endl;
  cout << "            '      V           '  "<< endl;
  cout << "          '                      ' "<< endl;
  cout << "        '                          '   " << endl;
  cout << "      '                              ' " << endl;
  cout << "       _  _        _   _        _  _  " << endl;
  cout << "      |  V | |  | |_| |_| \\  / |_ |_  " << endl;
  cout << "      |    | |__| | \\ | |  \\/  |_  _| " << endl;
  cout  << endl;

  cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  " << endl;

  // Measure the total time of the code 
  clock_t start = clock();
  
  // TELESCOPE  AND RUN SELECTION //
  char color[10];
  strcpy(color,argv[1]);
  int run;
  stringstream color_ss(argv[1]);
  string color_st;
  color_ss>>color_st;
  stringstream run_ss(argv[2]);
  run_ss >> run;
  char run_string[10];
  sprintf(run_string, "%d", run);
  /////////////////////////////////////

  int useSingleRunPED=1;
  
  /////// GEOMETRICAL PARAMETERS ///////
  double Z_add[4], X_pos[4];
  if(color_st == "ROSSO") {
    double Z_add_ROSSO[]= {0.292, 0.251, 0.207,0.0};
    double X_pos_ROSSO[] = {-0.265,0,0.262,1.475}; 
    for(int g=0; g<4; g++) {
      Z_add[g] = Z_add_ROSSO[g];
      X_pos[g] = X_pos_ROSSO[g];
    }
  }
  if(color_st == "NERO") {
    double Z_add_NERO[]= {0.293,0.251,0.210,0.0};
    //    double Z_add_NERO[]= {0.293,0.253,0.210,0.0}; PreCorrection z2 = z2 - 0.002
    double X_pos_NERO[] = {-0.26,0,0.262,1.492};
    for(int g=0; g<4; g++) {
      Z_add[g] = Z_add_NERO[g];
      X_pos[g] = X_pos_NERO[g];
    }
  }
  if(color_st == "BLU") {
    double Z_add_BLU[]= {0.2712,0.2312,0.1892,0.0};
    double X_pos_BLU[] = {-0.26,0,0.262,1.492};
    for(int g=0; g<4; g++) {
      Z_add[g] = Z_add_BLU[g];
      X_pos[g] = X_pos_BLU[g];
    }
  }

  double Y_add[] = {0.,0.,0.};
  
  /// Spatial resolutions 
  double sigma_z = 0.0040;
  double sigma_y = 0.0035;
  ///////////////////////////////////////////////////////

  //CLUSTERING PARAMETERS
  const double s1=6.; // cluster strips energy 
  const double s2 = 10.; // single-strip energy
  const double s3 = 2.; // adiacent strips of a single-strip cluster
  //////////////////////////////////

  //TRACKING PARAMETERS
  double proximity_cut_xz = 5*sigma_z; // point-trackcandidate distance requirement
  double proximity_cut_xy = 5*sigma_y; // point-trackcandidate distance requirement
  
  // CONSTANT PARAMETERS
  const int nBoards = 16;
  const int nInfo = 168;
  const int nChInfo = 5;
  const int nChannels = 32;
  const double FirstStripPos = -0.528;
  const double AdiacentStripsDistance = 0.0165;
  /////////////////////////////////

  ///////////// OUTPUTS //////////////////
  /// Mini-Tree
  char  MiniRunTreeName[200];
  char SavingPath_miniTree[100];
  strcpy(SavingPath_miniTree,"/workspace/test/RECONSTRUCTED/");
  strncat(SavingPath_miniTree,color,10);
  strncat(SavingPath_miniTree,"/",1);
  //strncat(SavingPath_miniTree,"/MINI_TREESLowerEnTh/",45);
  //strncat(SavingPath_miniTree,"/MINI_TREES/",35);
  strcpy(MiniRunTreeName,SavingPath_miniTree);
  strncat(MiniRunTreeName,"MURAVES_miniRunTree_run",30);
  strncat(MiniRunTreeName,run_string,10);
  //strncat(MiniRunTreeName,"_",2);
  //strncat(MiniRunTreeName,color,10);
  strncat(MiniRunTreeName,".root",10);
  
  // Analysis Tree
  char ROOTfileName[200];
  //char SavingPath_AnalysisTree[100];
  //strcpy(SavingPath_AnalysisTree,"/workspace/test/RECONSTRUCTED/");
  //strncat(SavingPath_AnalysisTree,color,10);
  //  strncat(SavingPath_AnalysisTree,"/ANALYZED_RECOv2/ClusteringWithTrMask/",50);
  //strncat(SavingPath_AnalysisTree,"/ANALYZED_RECOv2/",45);
  //if(useSingleRunPED==0) strncat(SavingPath_AnalysisTree,"GLOB_PED/",12);
  //cout << SavingPath_AnalysisTree << endl;
  strcpy(ROOTfileName,SavingPath_miniTree);
  strncat(ROOTfileName,	"MURAVES_AnalyzedData_run",30);
  strncat(ROOTfileName,run_string,10);
  //strncat(ROOTfileName,"_",2);
  //strncat(ROOTfileName,color,10);
  strncat(ROOTfileName,".root",10);

  
  // PRINT INFORMATIONS
  cout << "ANALYZING RUN : " << run << " OF "<< color_st << " DETECTOR " << endl;
  cout << "Clustering parameters: " << endl;
  cout << "Single strip min energy: " << s1 << endl;
  cout << "Single strip cluster min energy: " << s2 << endl;
  cout << "Single strip cluster adiacent strips  min energy: " << s3 << endl;

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  // FIND THE RUN ADC FILE ///////////////////////////////////////////////////////////////////////////////
  string line;
  char  ADC_File_name[150];
  char PATH_ADC[100];
  char FileName[50];

  strcpy(FileName,"ADC_run");
  strncat(FileName,run_string,10);
  strncat(FileName, ".txt", 5);
  strcpy(PATH_ADC,"/workspace/test/PARSED/");
  strncat(PATH_ADC,color,10);
  strncat(PATH_ADC, "/v0/", 6); // should be changed to wildcard version
  strcpy(ADC_File_name,PATH_ADC);
  //strncat(ADC_File_name,"/PARSED/",30);
  strncat(ADC_File_name,FileName,50);

  cout<< ADC_File_name <<endl;
  // USE PYTHON GLOB.GLOB TO FIND THE TAIL OF THE FILE NAME /////
  Py_Initialize();
  // Add the current directory to sys.path
  PyRun_SimpleString("import sys");
  PyRun_SimpleString("import os");
  PyRun_SimpleString("sys.path.append(os.getcwd())");
  CPyInstance hInstance;
  CPyObject pName_searchFile = PyUnicode_FromString("SearchFileName");
  CPyObject Module_SearchAFile = PyImport_Import(pName_searchFile);
  if (!Module_SearchAFile) {
    PyErr_Print(); // This prints the Python error message
    std::cerr << "Failed to import module SearchFileName\n";
    return 1;
  }
   
  CPyObject SearchFile_func = PyObject_GetAttrString(Module_SearchAFile, "Search_File");
  if (!SearchFile_func) {
    PyErr_Print(); // this will show the actual Python exception
    std::cerr << "SearchFile_func not found\n";
    return 1; // or handle differently
  } 
  PyObject *ADC_name_toPy = PyTuple_New(1);
  PyTuple_SetItem(ADC_name_toPy, 0, PyUnicode_FromString(ADC_File_name));
  CPyObject CompleteFileName;
  CompleteFileName =  PyObject_CallObject(SearchFile_func,ADC_name_toPy);
  if (!CompleteFileName) {
    PyErr_Print();  // this will show the actual Python exception
    std::cerr << "Error: call to SearchFile_func failed" << std::endl;
    return 1; // or handle differently
  }

  if (!PyUnicode_Check(CompleteFileName.getObject())) {
    std::cerr << "Search_File did not return a string\n";
    return 1;
}

  Py_ssize_t size;
  const char *Complete_ADCfile_name = PyUnicode_AsUTF8AndSize(CompleteFileName.getObject(), &size); // ----> COMPLETE NAME 
  cout << "ADC file: " ;

  puts(Complete_ADCfile_name); 
  Py_Finalize();
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////

  ///////////////// SLOW CONTROL INFORMATIONS: Working Point, Trigger Rate, External Temperature ///////////
  char SlowControlFileName[100];
  strcpy(SlowControlFileName,"/workspace/test/RAW_GZ/");
  strncat(SlowControlFileName,color,5);
  strncat(SlowControlFileName,"/SLOWCONTROL_run",30);
  strncat(SlowControlFileName, run_string, 10);
  //strncat(SlowControlFileName,".txt",5);
 
  ifstream SlowControl(SlowControlFileName);
  string slow_line;
  char c_slow_line[1000]; 
  char *char_run, *char_tr, *char_wp, *char_temperature;
  double SC_run,SC_tr,SC_temperature,SC_wp;
  int SC_run_int;
  char *slow_ptr;
  vector<char*>  SLOW_DATAdata_splitted;

  ////////////////////////// READ SLOW CONTROL LINE /////////////////////////////////////////
  while (getline(SlowControl, slow_line))
    {
      SLOW_DATAdata_splitted.clear();
      memset(c_slow_line, 0, sizeof slow_line);
      strcpy (c_slow_line, slow_line.c_str());
      slow_ptr = strtok(c_slow_line,"\t");
      while(slow_ptr!=NULL) {
	
	  SLOW_DATAdata_splitted.push_back(slow_ptr);
	  slow_ptr = strtok(NULL,"\t");
	}
	SC_run = strtod(SLOW_DATAdata_splitted.at(0),&char_run); 
	SC_run_int = (int) SC_run;
	//////////// WHEN THE CORRESPONDING RUN IS FOUND //////////////////////
	if(SC_run == run) {
	SC_tr = (double) strtod(SLOW_DATAdata_splitted.at(43),&char_tr); 
	SC_temperature = (double) strtod(SLOW_DATAdata_splitted.at(3),&char_temperature); 
	SC_wp = strtod(SLOW_DATAdata_splitted.at(5),&char_wp);

	  break;
	}
      }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////// PREPARING THE TREE ////////////////////////////////////////////  
  // 1) CREATE THE ROOT FILE + TTREE 
  TTree *tree = new TTree("AnalyzedData","A tree to store MURAVES analyzed data");

  // 2) INSTANTIATE THE USEFUL VARIABLES 
  Int_t Board, channel, strip, ADCValue, Nmodule, Nplane, Nstation;
  Double_t Deposit, OnePhe, timeExp, temperature;

  // CLUSTERING VARIABLES //
  Int_t Nclusters_Z1, Nclusters_Z2, Nclusters_Z3, Nclusters_Z4, Nclusters_Y1, Nclusters_Y2, Nclusters_Y3, Nclusters_Y4;
  vector<int> ClusterSize_Z1, ClusterSize_Z2, ClusterSize_Z3, ClusterSize_Z4, ClusterSize_Y1, ClusterSize_Y2, ClusterSize_Y3, ClusterSize_Y4;
  vector<double> ClusterEnergy_Z1, ClusterEnergy_Z2, ClusterEnergy_Z3, ClusterEnergy_Z4, ClusterEnergy_Y1, ClusterEnergy_Y2, ClusterEnergy_Y3, ClusterEnergy_Y4;
  vector<double> ClusterPosition_Z1, ClusterPosition_Z2, ClusterPosition_Z3, ClusterPosition_Z4, ClusterPosition_Y1, ClusterPosition_Y2, ClusterPosition_Y3, ClusterPosition_Y4;
  vector <double> ClusterZ1_Texp,ClusterZ2_Texp,ClusterZ3_Texp,ClusterZ4_Texp,ClusterY1_Texp,ClusterY2_Texp,ClusterY3_Texp,ClusterY4_Texp;
  vector<vector<double>> StripsEnergy_Z1,StripsEnergy_Z2,StripsEnergy_Z3, StripsEnergy_Z4,StripsEnergy_Y1,StripsEnergy_Y2,StripsEnergy_Y3, StripsEnergy_Y4;
  vector<vector<double>> StripsPosition_Z1,StripsPosition_Z2,StripsPosition_Z3, StripsPosition_Z4,StripsPosition_Y1,StripsPosition_Y2,StripsPosition_Y3, StripsPosition_Y4;
  vector<vector<double>> StripsID_Z1,StripsID_Z2,StripsID_Z3, StripsID_Z4,StripsID_Y1,StripsID_Y2,StripsID_Y3, StripsID_Y4;

  // TRACKING VARIABLES //
  Int_t Ntracks_3p_xz, Ntracks_3p_xy;
  vector<int> Ntracks_4p_xy, Ntracks_4p_xz,Plane4th_isIntercepted_xz, Plane4th_isIntercepted_xy; 
  vector<double> Intercept_3p_xz, Slope_3p_xz, chiSquare_3p_xz,Intercept_3p_xy, Slope_3p_xy, chiSquare_3p_xy, TrackCluster_z1_Position,TrackCluster_z2_Position,TrackCluster_z3_Position,TrackCluster_y1_Position,TrackCluster_y2_Position,TrackCluster_y3_Position;
  vector<double> Residue_Track3p_z1,Residue_Track3p_z2,Residue_Track3p_z3,Residue_Track3p_y1,Residue_Track3p_y2,Residue_Track3p_y3;
  vector<int> TrackCluster_z1_index,TrackCluster_z2_index,TrackCluster_z3_index,TrackCluster_z4_index, TrackCluster_y1_index,TrackCluster_y2_index,TrackCluster_y3_index,TrackCluster_y4_index;
  vector<double> TrackCluster_z1_Texp,TrackCluster_z2_Texp,TrackCluster_z3_Texp,TrackCluster_z4_Texp, TrackCluster_y1_Texp,TrackCluster_y2_Texp,TrackCluster_y3_Texp,TrackCluster_y4_Texp;  
  vector<double> TrackCluster_z1_Size,TrackCluster_z2_Size,TrackCluster_z3_Size, TrackCluster_y1_Size,TrackCluster_y2_Size,TrackCluster_y3_Size;
  vector<double> TrackEnergy_3p_xy, TrackEnergy_3p_xz,Phi_3p,Elevation_3p, clusters4_Texp_single,ExpectedPosition_OnPlane4th_xy,ExpectedPosition_OnPlane4th_xz;
  vector<vector<double>> TrackCluster_z4_Position, TrackCluster_y4_Position, Intercept_4p_xz, Slope_4p_xz, chiSquare_4p_xz,Intercept_4p_xy, Slope_4p_xy, chiSquare_4p_xy,displacement_p4_xz, displacement_p4_xy, residue_c1_p4_xz,residue_c2_p4_xz,residue_c3_p4_xz,residue_c4_p4_xz, residue_c1_p4_xy,residue_c2_p4_xy,residue_c3_p4_xy,residue_c4_p4_xy, TrackEnergy_4p_xy,TrackEnergy_4p_xz;
  vector<vector<double>> cluster_c4_index_xy, cluster_c4_index_xz, ScatteringAngle_xy,ScatteringAngle_xz,ScatteringAngle_3D,Phi_4p,Elevation_4p, cluster_c4_Texp_xy, cluster_c4_Texp_xz;
  double phi3p, elev3p, phi4p, elev4p, Dx,Dz,Dy;
  vector<int> isInTrack_3p_clZ1, isInTrack_3p_clZ2,isInTrack_3p_clZ3,isInTrack_3p_clZ4,isInTrack_3p_clY1, isInTrack_3p_clY2,isInTrack_3p_clY3,isInTrack_3p_clY4;
  vector<int> isInTrack_4p_clZ1, isInTrack_4p_clZ2,isInTrack_4p_clZ3,isInTrack_4p_clZ4,isInTrack_4p_clY1, isInTrack_4p_clY2,isInTrack_4p_clY3,isInTrack_4p_clY4;
  vector<double> Tracks3p_ExpectedRes_z2, Tracks3p_ExpectedRes_y2;
  ///// BEST TRACK VARIABLES //////
  double BestChi_xy, BestEnergy_xy,BestChi_xz, BestEnergy_xz; 
  int BestChi_xy_index,BestEnergy_xy_index,BestChi_xz_index,BestEnergy_xz_index;
  double theta_3p,phi_3p, theta_4p,phi_4p;
  vector<int> Track_3p_to_4p_index_xy,Track_3p_to_4p_index_xz;
  vector<int> isBestTrack_3p_xy,isBestTrack_3p_xz,isBestTrack_4p_xy_single,isBestTrack_4p_xz_single;
  vector<vector<double>> isBestTrack_4p_xy,isBestTrack_4p_xz;
  int bestTracks_4p_index_xz=0, bestTracks_4p_index_xy=0,index_4p_xz,index_4p_xy;
  vector<int>::iterator track_4p_index_xy,track_4p_index_xz;
  int BestTrack_3p_xy_index, BestTrack_3p_xz_index, Track_3p_of_4p_index_xy,Track_3p_of_4p_index_xz,Track_4p_index_xy,Track_4p_index_xz;
  double BestTrack_3p_Energy_xy,BestTrack_3p_Energy_xz, BestTrack_4p_Energy_xy,BestTrack_4p_Energy_xz;
  double BestTrack_3p_ChiSquare_xy,BestTrack_3p_ChiSquare_xz, BestTrack_4p_ChiSquare_xy,BestTrack_4p_ChiSquare_xz;
  double BestTracks_ScatteringAngle_xy, BestTracks_ScatteringAngle_xz;
  double Residue_BestTracks3p_z1,Residue_BestTracks3p_z2,Residue_BestTracks3p_z3,Residue_BestTracks3p_y1,Residue_BestTracks3p_y2,Residue_BestTracks3p_y3;
  double Residue_BestTracks4p_z1,Residue_BestTracks4p_z2,Residue_BestTracks4p_z3,Residue_BestTracks4p_z4,Residue_BestTracks4p_y1,Residue_BestTracks4p_y2,Residue_BestTracks4p_y3,Residue_BestTracks4p_y4;
  int Best_track_4p_isTexpNULL_xz, Best_track_4p_isTexpNULL_xy;

  //// TRIGGER MASK //////
  vector<vector<double>> TriggerMaskStrips, TriggerMaskChannels;
  vector<int> TriggerMaskSize;
    
  // OTHER VARIABLES (MINI TREE) ///
  double RunDuration;
  Long_t timestamp;
  time_t Time; 
  TDatime datime;
  string view;
  double WorkingPoint,Temperature, TriggerRate;
  int W_P; 
  double EnergyCut_clusterStrip=s1, EnergyCut_singleStrip=s2, EnergyCut_additionalStrip=s3;
  // CLUSTERS BRANCHES //
  tree->Branch("Run",&run);
  tree->Branch("Nclusters_Z1",&Nclusters_Z1,"Nclusters_Z1/I");
  tree->Branch("Nclusters_Z2",&Nclusters_Z2,"Nclusters_Z2/I");
  tree->Branch("Nclusters_Z3",&Nclusters_Z3,"Nclusters_Z3/I");
  tree->Branch("Nclusters_Z4",&Nclusters_Z4,"Nclusters_Z4/I");
  tree->Branch("Nclusters_Y1",&Nclusters_Y1,"Nclusters_Y1/I");
  tree->Branch("Nclusters_Y2",&Nclusters_Y2,"Nclusters_Y2/I");
  tree->Branch("Nclusters_Y3",&Nclusters_Y3,"Nclusters_Y3/I");
  tree->Branch("Nclusters_Y4",&Nclusters_Y4,"Nclusters_Y4/I");

  tree->Branch("ClusterSize_Z1",&ClusterSize_Z1);
  tree->Branch("ClusterSize_Z2",&ClusterSize_Z2);
  tree->Branch("ClusterSize_Z3",&ClusterSize_Z3);
  tree->Branch("ClusterSize_Z4",&ClusterSize_Z4);
  tree->Branch("ClusterZ1_Texp",&ClusterZ1_Texp);
  tree->Branch("ClusterZ2_Texp",&ClusterZ2_Texp);
  tree->Branch("ClusterZ3_Texp",&ClusterZ3_Texp);
  tree->Branch("ClusterZ4_Texp",&ClusterZ4_Texp);
  tree->Branch("ClusterEnergy_Z1",&ClusterEnergy_Z1);
  tree->Branch("ClusterEnergy_Z2",&ClusterEnergy_Z2);
  tree->Branch("ClusterEnergy_Z3",&ClusterEnergy_Z3);
  tree->Branch("ClusterEnergy_Z4",&ClusterEnergy_Z4);
  tree->Branch("ClusterPosition_Z1",&ClusterPosition_Z1);
  tree->Branch("ClusterPosition_Z2",&ClusterPosition_Z2);
  tree->Branch("ClusterPosition_Z3",&ClusterPosition_Z3);
  tree->Branch("ClusterPosition_Z4",&ClusterPosition_Z4),

  tree->Branch("ClusterSize_Y1",&ClusterSize_Y1);
  tree->Branch("ClusterSize_Y2",&ClusterSize_Y2);
  tree->Branch("ClusterSize_Y3",&ClusterSize_Y3);
  tree->Branch("ClusterSize_Y4",&ClusterSize_Y4);
  tree->Branch("ClusterY1_Texp",&ClusterY1_Texp);
  tree->Branch("ClusterY2_Texp",&ClusterY2_Texp);
  tree->Branch("ClusterY3_Texp",&ClusterY3_Texp);
  tree->Branch("ClusterY4_Texp",&ClusterY4_Texp);  
  tree->Branch("ClusterEnergy_Y1",&ClusterEnergy_Y1);
  tree->Branch("ClusterEnergy_Y2",&ClusterEnergy_Y2);
  tree->Branch("ClusterEnergy_Y3",&ClusterEnergy_Y3);
  tree->Branch("ClusterEnergy_Y4",&ClusterEnergy_Y4);
  tree->Branch("ClusterPosition_Y1",&ClusterPosition_Y1);
  tree->Branch("ClusterPosition_Y2",&ClusterPosition_Y2);
  tree->Branch("ClusterPosition_Y3",&ClusterPosition_Y3);
  tree->Branch("ClusterPosition_Y4",&ClusterPosition_Y4);

  tree->Branch("StripsEnergy_Z1", &StripsEnergy_Z1);
  tree->Branch("StripsEnergy_Z2", &StripsEnergy_Z2);
  tree->Branch("StripsEnergy_Z3", &StripsEnergy_Z3);
  tree->Branch("StripsEnergy_Z4", &StripsEnergy_Z4);
  tree->Branch("StripsEnergy_Y1", &StripsEnergy_Y1);
  tree->Branch("StripsEnergy_Y2", &StripsEnergy_Y2);
  tree->Branch("StripsEnergy_Y3", &StripsEnergy_Y3);
  tree->Branch("StripsEnergy_Y4", &StripsEnergy_Y4);
  tree->Branch("StripsPosition_Z1",&StripsPosition_Z1);
  tree->Branch("StripsPosition_Z2",&StripsPosition_Z2);
  tree->Branch("StripsPosition_Z3",&StripsPosition_Z3);
  tree->Branch("StripsPosition_Z4",&StripsPosition_Z4);
  tree->Branch("StripsPosition_Y1",&StripsPosition_Y1);
  tree->Branch("StripsPosition_Y2",&StripsPosition_Y2);
  tree->Branch("StripsPosition_Y3",&StripsPosition_Y3);
  tree->Branch("StripsPosition_Y4",&StripsPosition_Y4);
  tree->Branch("StripsID_Z1",&StripsID_Z1);
  tree->Branch("StripsID_Z2",&StripsID_Z2);
  tree->Branch("StripsID_Z3",&StripsID_Z3);
  tree->Branch("StripsID_Z4",&StripsID_Z4);
  tree->Branch("StripsID_Y1",&StripsID_Y1);
  tree->Branch("StripsID_Y2",&StripsID_Y2);
  tree->Branch("StripsID_Y3",&StripsID_Y3);
  tree->Branch("StripsID_Y4",&StripsID_Y4);



  // TRACKING BRANCHES
  tree->Branch("Ntracks_3p_xz",&Ntracks_3p_xz,"Ntracks_3p_xz/I");
  tree->Branch("Ntracks_3p_xy",&Ntracks_3p_xy,"Ntracks_3p_xy/I");
  tree->Branch("Intercept_3p_xz",&Intercept_3p_xz);
  tree->Branch("Slope_3p_xz",&Slope_3p_xz);
  tree->Branch("chiSquare_3p_xz",&chiSquare_3p_xz);
  tree->Branch("Intercept_3p_xy",&Intercept_3p_xy);
  tree->Branch("Slope_3p_xy",&Slope_3p_xy);
  tree->Branch("chiSquare_3p_xy",&chiSquare_3p_xy);
  tree->Branch("TrackCluster_z1_index",&TrackCluster_z1_index);
  tree->Branch("TrackCluster_z2_index",&TrackCluster_z2_index);
  tree->Branch("TrackCluster_z3_index",&TrackCluster_z3_index);
  tree->Branch("TrackCluster_y1_index",&TrackCluster_y1_index);
  tree->Branch("TrackCluster_y2_index",&TrackCluster_y2_index);
  tree->Branch("TrackCluster_y3_index",&TrackCluster_y3_index);

  /*
  tree->Branch("TrackCluster_z1_Position",&TrackCluster_z1_Position);
  tree->Branch("TrackCluster_z2_Position",&TrackCluster_z2_Position);
  tree->Branch("TrackCluster_z3_Position",&TrackCluster_z3_Position);
  tree->Branch("TrackCluster_y1_Position",&TrackCluster_y1_Position);
  tree->Branch("TrackCluster_y2_Position",&TrackCluster_y2_Position);
  tree->Branch("TrackCluster_y3_Position",&TrackCluster_y3_Position);
  tree->Branch("TrackCluster_z1_Texp",&TrackCluster_z1_Texp);
  tree->Branch("TrackCluster_z2_Texp",&TrackCluster_z2_Texp);
  tree->Branch("TrackCluster_z3_Texp",&TrackCluster_z3_Texp);
  tree->Branch("TrackCluster_y1_Texp",&TrackCluster_y1_Texp);
  tree->Branch("TrackCluster_y2_Texp",&TrackCluster_y2_Texp);
  tree->Branch("TrackCluster_y3_Texp",&TrackCluster_y3_Texp);
  tree->Branch("TrackCluster_z1_Size",&TrackCluster_z1_Size);
  tree->Branch("TrackCluster_z2_Size",&TrackCluster_z2_Size);
  tree->Branch("TrackCluster_z3_Size",&TrackCluster_z3_Size);
  tree->Branch("TrackCluster_y1_Size",&TrackCluster_y1_Size);
  tree->Branch("TrackCluster_y2_Size",&TrackCluster_y2_Size);
  tree->Branch("TrackCluster_y3_Size",&TrackCluster_y3_Size);
  */  
  tree->Branch("Residue_Track3p_z1",&Residue_Track3p_z1);
  tree->Branch("Residue_Track3p_z2",&Residue_Track3p_z2);
  tree->Branch("Residue_Track3p_z3",&Residue_Track3p_z3);
  tree->Branch("Residue_Track3p_y1",&Residue_Track3p_y1);
  tree->Branch("Residue_Track3p_y2",&Residue_Track3p_y2);
  tree->Branch("Residue_Track3p_y3",&Residue_Track3p_y3);
  tree->Branch("TrackEnergy_3p_xy",&TrackEnergy_3p_xy);
  tree->Branch("TrackEnergy_3p_xz",&TrackEnergy_3p_xz);
  tree->Branch("Plane4th_isIntercepted_xz",&Plane4th_isIntercepted_xz);
  tree->Branch("Plane4th_isIntercepted_xy",&Plane4th_isIntercepted_xy);  
  tree->Branch("BestChi_xy_index", &BestChi_xy_index);
  tree->Branch("BestEnergy_xy_index", &BestEnergy_xy_index);
  tree->Branch("BestChi_xz_index", &BestChi_xz_index);
  tree->Branch("BestEnergy_xz_index", &BestEnergy_xz_index);
  tree->Branch("BestChi_xy", &BestChi_xy);
  tree->Branch("BestEnergy_xy", &BestEnergy_xy);
  tree->Branch("BestChi_xz", &BestChi_xz);
  tree->Branch("BestEnergy_xz", &BestEnergy_xz);
  tree->Branch("ExpectedPosition_OnPlane4th_xy",&ExpectedPosition_OnPlane4th_xy);
  tree->Branch("ExpectedPosition_OnPlane4th_xz",&ExpectedPosition_OnPlane4th_xz);
  //TRACKING BRANCHES ---> 4 PLANES
  //tree->Branch("Ntracks_xz_4p",&Ntracks_4p_xz);
  //tree->Branch("Ntracks_xy_4p",&Ntracks_4p_xy);
  //tree->Branch("TrackCluster_z4_Position",&TrackCluster_z4_Position);
  //tree->Branch("TrackCluster_y4_Position",&TrackCluster_y4_Position);
  //tree->Branch("cluster_c4_Texp_xz",&cluster_c4_Texp_xz);
  //tree->Branch("cluster_c4_Texp_xy",&cluster_c4_Texp_xy);
  /*  tree->Branch("residue_c1_p4_xy",&residue_c1_p4_xy);
  tree->Branch("residue_c2_p4_xy",&residue_c2_p4_xy);
  tree->Branch("residue_c3_p4_xy",&residue_c3_p4_xy);
  tree->Branch("residue_c4_p4_xy",&residue_c4_p4_xy);
  tree->Branch("residue_c1_p4_xz",&residue_c1_p4_xz);
  tree->Branch("residue_c2_p4_xz",&residue_c2_p4_xz);
  tree->Branch("residue_c3_p4_xz",&residue_c3_p4_xz);
  tree->Branch("residue_c4_p4_xz",&residue_c4_p4_xz);
  tree->Branch("Intercept_4p_xz",&Intercept_4p_xz);
  */
  tree->Branch("Slope_4p_xz",&Slope_4p_xz);
  tree->Branch("chiSquare_4p_xz",&chiSquare_4p_xz);
  tree->Branch("displacement_p4_xz",&displacement_p4_xz);


  tree->Branch("Intercept_4p_xy",&Intercept_4p_xy);
  tree->Branch("Slope_4p_xy",&Slope_4p_xy);
  tree->Branch("chiSquare_4p_xy",&chiSquare_4p_xy);
  tree->Branch("displacement_p4_xy",&displacement_p4_xy);
  tree->Branch("cluster_c4_index_xz",&cluster_c4_index_xz);
  tree->Branch("cluster_c4_index_xy",&cluster_c4_index_xy);
  tree->Branch("ScatteringAngle_xy",&ScatteringAngle_xy);
  tree->Branch("ScatteringAngle_xz",&ScatteringAngle_xz);

  //CLUSTERS OF TRACKS
  tree->Branch("isInTrack_3p_clZ1", &isInTrack_3p_clZ1);
  tree->Branch("isInTrack_3p_clZ2", &isInTrack_3p_clZ2);
  tree->Branch("isInTrack_3p_clZ3", &isInTrack_3p_clZ3);
  tree->Branch("isInTrack_3p_clZ4", &isInTrack_3p_clZ4);
  tree->Branch("isInTrack_3p_clY1", &isInTrack_3p_clY1);
  tree->Branch("isInTrack_3p_clY2", &isInTrack_3p_clY2);
  tree->Branch("isInTrack_3p_clY3", &isInTrack_3p_clY3);
  tree->Branch("isInTrack_3p_clY4", &isInTrack_3p_clY4);
  tree->Branch("isInTrack_4p_clZ1", &isInTrack_4p_clZ1);
  tree->Branch("isInTrack_4p_clZ2", &isInTrack_4p_clZ2);
  tree->Branch("isInTrack_4p_clZ3", &isInTrack_4p_clZ3);
  tree->Branch("isInTrack_4p_clZ4", &isInTrack_4p_clZ4);
  tree->Branch("isInTrack_4p_clY1", &isInTrack_4p_clY1);
  tree->Branch("isInTrack_4p_clY2", &isInTrack_4p_clY2);
  tree->Branch("isInTrack_4p_clY3", &isInTrack_4p_clY3);
  tree->Branch("isInTrack_4p_clY4", &isInTrack_4p_clY4);

  /////  BEST TRACKS  ////////
  tree->Branch("Theta_3p", &theta_3p);
  tree->Branch("Theta_4p", &theta_4p);
  tree->Branch("Phi_3p", &phi_3p);
  tree->Branch("Phi_4p", &phi_4p);
  tree->Branch("isBestTrack_3p_xy", &isBestTrack_3p_xy);
  tree->Branch("isBestTrack_3p_xz", &isBestTrack_3p_xz);
   tree->Branch("BestTrack_3p_xy_index",&BestTrack_3p_xy_index);
  tree->Branch("BestTrack_3p_xz_index",&BestTrack_3p_xz_index);
  tree->Branch("Track_3p_of_4p_index_xy",&Track_3p_of_4p_index_xy);
  tree->Branch("Track_3p_of_4p_index_xz",&Track_3p_of_4p_index_xz);
  tree->Branch("Track_4p_index_xy",&Track_4p_index_xy);
  tree->Branch("Track_4p_index_xz",&Track_4p_index_xz);
  /*
  tree->Branch("BestTrack_3p_Energy_xy",&BestTrack_3p_Energy_xy);
  tree->Branch("BestTrack_3p_Energy_xz",&BestTrack_3p_Energy_xz);
  tree->Branch("BestTrack_4p_Energy_xy",&BestTrack_4p_Energy_xy);
  tree->Branch("BestTrack_4p_Energy_xz",&BestTrack_4p_Energy_xz);
  tree->Branch("Residue_BestTracks3p_z1",&Residue_BestTracks3p_z1);
  tree->Branch("Residue_BestTracks3p_z2",&Residue_BestTracks3p_z2);
  tree->Branch("Residue_BestTracks3p_z3",&Residue_BestTracks3p_z3);
  tree->Branch("Residue_BestTracks3p_y1",&Residue_BestTracks3p_y1);
  tree->Branch("Residue_BestTracks3p_y2",&Residue_BestTracks3p_y2);
  tree->Branch("Residue_BestTracks3p_y3",&Residue_BestTracks3p_y3);
  tree->Branch("Tracks3p_ExpectedRes_z2",&Tracks3p_ExpectedRes_z2);
  tree->Branch("Tracks3p_ExpectedRes_y2",&Tracks3p_ExpectedRes_y2);
  tree->Branch("Residue_BestTracks4p_z1",&Residue_BestTracks4p_z1);
  tree->Branch("Residue_BestTracks4p_z2",&Residue_BestTracks4p_z2);
  tree->Branch("Residue_BestTracks4p_z3",&Residue_BestTracks4p_z3);
  tree->Branch("Residue_BestTracks4p_z4",&Residue_BestTracks4p_z4);
  tree->Branch("Residue_BestTracks4p_y1",&Residue_BestTracks4p_y1);
  tree->Branch("Residue_BestTracks4p_y2",&Residue_BestTracks4p_y2);
  tree->Branch("Residue_BestTracks4p_y3",&Residue_BestTracks4p_y3);
  tree->Branch("Residue_BestTracks4p_y4",&Residue_BestTracks4p_y4);
  */
  tree->Branch("BestTrack_3p_ChiSquare_xy",&BestTrack_3p_ChiSquare_xy);
  tree->Branch("BestTrack_3p_ChiSquare_xz",&BestTrack_3p_ChiSquare_xz);
  tree->Branch("BestTrack_4p_ChiSquare_xy",&BestTrack_4p_ChiSquare_xy);
  tree->Branch("BestTrack_4p_ChiSquare_xz",&BestTrack_4p_ChiSquare_xz);

  tree->Branch("BestTracks_ScatteringAngle_xy",&BestTracks_ScatteringAngle_xy);
  tree->Branch("BestTracks_ScatteringAngle_xz",&BestTracks_ScatteringAngle_xz);
  tree->Branch("Best_track_4p_isTexpNULL_xz",&Best_track_4p_isTexpNULL_xz);
  tree->Branch("Best_track_4p_isTexpNULL_xy",&Best_track_4p_isTexpNULL_xy);
  
  // OTHER BRACHES //
  tree->Branch("datime",&datime);
  tree->Branch("WorkingPoint",&WorkingPoint);
  tree->Branch("Temperature",&Temperature);
  tree->Branch("TriggerRate",&TriggerRate);


  // TRIGGER MASK /////////
  tree->Branch("TriggerMaskChannels",&TriggerMaskChannels);
  tree->Branch("TriggerMaskStrips",&TriggerMaskStrips);
  tree->Branch("TriggerMaskSize",&TriggerMaskSize);
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////
  // SLOW CONTROL DATA ///
  TriggerRate = SC_tr;
  Temperature = SC_temperature;
  WorkingPoint = (int) SC_wp;
  cout << "Working Point: " << WorkingPoint << " Temperature: " << Temperature << " Trigger Rate: " << TriggerRate << endl;
  W_P = (int)WorkingPoint;
  //FIST AND LAST LINES TO EVALUATE THE TIME  
  char first_line[10000], last_line[10000];
  ////////////////////////////////////////////////// VARIABLES NEEDED TO THE CLUSTERING ///////////////////////////////////////////////////
  // Vectors of module deposits 
  vector <double> Deposits_p1x1, Deposits_p1x2, Deposits_p2x1, Deposits_p2x2, Deposits_p3x1, Deposits_p3x2, Deposits_p4x1, Deposits_p4x2;
  vector <double> Deposits_p1y1, Deposits_p1y2, Deposits_p2y1, Deposits_p2y2, Deposits_p3y1, Deposits_p3y2, Deposits_p4y1, Deposits_p4y2;
  vector <double> Deposits_p1x, Deposits_p2x, Deposits_p3x,Deposits_p4x;
  vector <double> Deposits_p1y, Deposits_p2y, Deposits_p3y,Deposits_p4y;
  // Clustering outputs 
  ClusterCollection resuts_p1x, resuts_p2x,resuts_p3x,resuts_p4x, resuts_p1y,resuts_p2y,resuts_p3y,resuts_p4y;
  struct_Event EventInfo; 
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////// VARIABLES NEEDED TO THE TRACKING ///////////////////////////////////////////////////
  vector<double> tracks_3p_intercepts, tracks_3p_slopes, tracks_3p_chiSquares;
  TracksCollection tracks_xz,tracks_xy;
  int BestTrack_xy_ind, BestTrack_xz_ind;
  vector<double> Coordinates_3p, Coordinates_4p;
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  // CONFIGURATION OF THE CHANNELS ---> READ SPIROC/HYBRID MAP
  ifstream Spiroc_config;
  int st, ch;
  vector <int>  strips, channels, sorted_ch,sorted_strips;
  string line_c;  
  Spiroc_config.open("/workspace/Software/tracks_reconstruction/AncillaryFiles/spiroc-hybrid-map.cfg");
  int l=0;
  /// //READING SPIROC/HYRID MAP FILE //////////////////////////////////////////////////////////////
  if(Spiroc_config.is_open()) {
    while (getline(Spiroc_config, line_c))
      {
	istringstream iss(line_c);
	if (l>0) {
	  iss >> st >> ch;
	  strips.push_back(st);
	  channels.push_back(ch);
	}
	  l++;
      }
  }
  Spiroc_config.close();
  //////////////////////////////////////////////////////////////////////////////////////////////////////
    
  // SORT CHANNELS AND STRIP IN ORDER OF STRIPS
  vector<int> strip_indices(strips.size());
  iota(strip_indices.begin(), strip_indices.end(), 0);
  sort(strip_indices.begin(), strip_indices.end(),
           [&](int A, int B) -> bool {
                return strips[A] < strips[B];
            });
  for(int ind=0; ind<strip_indices.size(); ind++) {
    sorted_strips.push_back(strips.at(strip_indices.at(ind)));
    sorted_ch.push_back(channels.at(strip_indices.at(ind)));
   

  }
    ///////////////////////////////////////////////////
  
  //////READ AND STORE ONEPHE AND PEDESTAL
  ///////// PED CONFIGURATION INPUTS /////////////////
  char ConfigPATH[100];
  char ConfigPed_file[100];

  ///// USE SINGLE_RUN PED:
  if(useSingleRunPED==1) {
    strcpy(ConfigPATH,"/workspace/test/PEDESTAL/");
    strcpy(ConfigPed_file,ConfigPATH);
    //char PED_FOLDER[15];
    //strcpy(PED_FOLDER,"ped");
    strncat(ConfigPed_file,color,6);
    strncat(ConfigPed_file, "/", 1);
    strncat(ConfigPed_file, run_string, 6);
    strncat(ConfigPed_file, "/", 1);
    cout<<ConfigPed_file<<endl;
    //char Run_PART[50];
    //strcpy(Run_PART,PED_FOLDER);
    //strncat(Run_PART,"/ped_run",20);
    //strncat(Run_PART,run_string,20);
    //strncat(ConfigPed_file,Run_PART,50);
    //strncat(ConfigPed_file, PATH_ADC, 100);
    //strncat(ConfigPed_file, "PEDESTAL/", 30);
  }
  else {
    strcpy(ConfigPATH,"config/");
    strcpy(ConfigPed_file,ConfigPATH);
    strncat(ConfigPed_file,color,6);
    strncat(ConfigPed_file,"/ped_WP",7);
    char workPoint_char[3];
    sprintf(workPoint_char,"%d",W_P);
    strncat(ConfigPed_file,workPoint_char,3);
    
  }
    char final_string[30];
    strcpy(final_string,"pedestal_"); 
    strncat(ConfigPed_file,final_string,30);
    cout << ConfigPed_file << endl;
  ifstream FilePED;
  double ped, OnePh;
  int n_ch,Iscopy;
  vector <double> Peds, sorted_Peds;
  vector<vector<double>> boards_Peds;
  vector<double> OnePhes,sorted_OnePhes,sorted_Is1phe_copy,Is1phe_copy;
  vector<vector<double>> boards_OnePhes,boards_Is1phe_copy;
  char Ped_fileName[100];
  char Nboard_string[5];

  for(int board_n=0; board_n<nBoards; board_n++) {
    sprintf(Nboard_string, "%d", board_n);  
    strcpy(Ped_fileName,ConfigPed_file);
    strncat(Ped_fileName,Nboard_string,5);
    strncat(Ped_fileName,".cfg",5);
    FilePED.open(Ped_fileName);
    l=0;
    if(FilePED.is_open()) {
      while (getline(FilePED, line_c))
	{
	  istringstream iss(line_c);
	  if (l>0) {
	    
	    //iss >> n_ch >> ped >> OnePh >> Iscopy;
	    iss >> n_ch >> ped >> OnePh;

	    Peds.push_back(ped);
	    OnePhes.push_back(OnePh);
	    Is1phe_copy.push_back(Iscopy);
	  }
          l++;
	}
      FilePED.close();
    }
    //////////// SORT FOLLOWING THE STRIP ORDER   /////////////////
    for(int ind=0; ind<sorted_ch.size(); ind++) {
      sorted_Peds.push_back(Peds.at( sorted_ch.at(ind)));
      sorted_OnePhes.push_back(OnePhes.at(sorted_ch.at(ind)));
      sorted_Is1phe_copy.push_back(Is1phe_copy.at(sorted_ch.at(ind)));
      
     
    }
    boards_Peds.push_back(sorted_Peds);
    boards_OnePhes.push_back(sorted_OnePhes);
    boards_Is1phe_copy.push_back(sorted_Is1phe_copy);
    boards_Is1phe_copy;
    sorted_Is1phe_copy.clear();
    Is1phe_copy.clear();
    sorted_OnePhes.clear();
    sorted_Peds.clear();
    Peds.clear();
    OnePhes.clear();
  }
  ///////////////////////////////////////////////////////////////////////  

  TFile *fileROOT = new TFile(ROOTfileName,"recreate");

  // READ TELESCOPE CONFIGURATION  ----> CORRESPONDACE BOARD - N STATION - VIEW //////// ///////////////////////////////
  char  telescopeConfig_name[10000];
  strcpy(telescopeConfig_name,"/workspace/Software/tracks_reconstruction/AncillaryFiles/telescope");
  strncat(telescopeConfig_name,color,10);
  strncat(telescopeConfig_name,".cfg",5);
  
  ifstream file_telescope;
  file_telescope.open(telescopeConfig_name);
  string tel_line;
  char *tel_ptr;
  char c_line[100];
  string view_teel;
  int n_Station;
  vector<int> nStations;
  vector<string> singleboard_tel_info, Views;
  vector<vector<string>> allboards_tel_info;
  int t=0;

  while (getline(file_telescope, tel_line)) {
    singleboard_tel_info.clear();
    strcpy (c_line, tel_line.c_str());
    if(t>0) {
	tel_ptr = strtok(c_line,"\t");
	while(tel_ptr!=NULL) {
	  string st_tel(tel_ptr);
	  singleboard_tel_info.push_back(st_tel);
	  tel_ptr = strtok(NULL,"\t");
	}
	stringstream intValue(singleboard_tel_info.at(2));
	intValue >> n_Station;
	nStations.push_back(n_Station);
	view_teel = singleboard_tel_info.at(4);
	Views.push_back(view_teel);
    }

    t++;
  }
  //////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////                               //////////////////////////////////
  /////////////////////////////  /// READ THE ADC FILE ////  //////////////////////////////////
  ////////////////////////////                              //////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////
  
  ifstream ADCfile;
  string event;
  char c_event[10000];
  ADCfile.open(Complete_ADCfile_name);

  // CALCULATE SOME OUTPUT INFORMATIONS
  int Nevents_withAtrack3p=0, Nevents_withAtrack4p=0, n4p_xy=0,n4p_xz=0, Nevents_withX1_cl=0, Nevents_withX2_cl=0,Nevents_withX3_cl=0,Nevents_withX4_cl=0,Nevents_withY1_cl=0, Nevents_withY2_cl=0,Nevents_withY3_cl=0,Nevents_withY4_cl=0, Nevents_withAtrack3p_xz=0, Nevents_withAtrack3p_xy=0; 
  int N_tracks_xy_3p=0, N_tracks_xz_3p=0, N_tracks_3p=0;
  int N_tracks_xy_4p=0, N_tracks_xz_4p=0, N_tracks_4p=0; 
  // Calculate Deposits ADC - PED / OnePhe -> NEEDED VARIABLES ////
  vector <vector<double>> AllBoards_deposits, AllBoards_ADC;
  vector<double> ADCcounts, ADCDeposits;
  double ADCcount,Depos;
  vector<double> Peds_board;
  vector<double> OnePhe_board;

  ////////////////////////////////////////////////////////////////


  //// TIME EXPANSION  ////
  double Texp_1x1, Texp_1x2, Texp_2x1, Texp_2x2, Texp_3x1, Texp_3x2, Texp_4x1, Texp_4x2;
  double Texp_1y1, Texp_1y2, Texp_2y1, Texp_2y2, Texp_3y1, Texp_3y2, Texp_4y1, Texp_4y2;

  // COUNT THE DIGIT OF THE RUN (NEEDED TO ASSIGN THE DATE)
  int count =0;
  int  temp = run;  

  while(temp != 0) {
    count++;
    temp /= 10;
  }
  int letters=0;

    while(color[letters] != '\0') {  
    letters++;
  }

  char dateTime[100];
  char year[5], month[5], day[5], hour[5], min[5], sec[5];
  strncpy(dateTime,Complete_ADCfile_name+38+letters+count,20);
  strncpy(year,dateTime,4);
  strncpy(month,dateTime+4,2);
  strncpy(day,dateTime+6,2);
  strncpy(hour,dateTime+9,2);
  strncpy(min,dateTime+11,2);
  strncpy(sec,dateTime+13,2);
  int y, mt,d, h,m,s;

  stringstream ss_year(year);
  stringstream ss_month(month);
  stringstream ss_day(day);
  stringstream ss_hour(hour);
  stringstream ss_min(min);
  stringstream ss_sec(sec);
  ss_year >> y;
  ss_month >> mt;
  ss_day >> d;
  ss_hour>> h;
  ss_min >> m;
  ss_sec >> s;
  datime.Set(y,mt,d,h,m,s);
  datime.Print();
  int ev =0;
int  NGooDTracks3p=0;
int  NGooDTracks4p=0;
  if (ADCfile.is_open()) {
    while (getline(ADCfile, event))
      {

	ev++;
	// LOOP OVER THE EVENTS
	memset(c_event, 0, sizeof c_event);
	strcpy (c_event, event.c_str());
	if(ev==1) strcpy(first_line, event.c_str());
	if(event!='\n')  	  strcpy(last_line, event.c_str());
	
	//// READ DATA RESULTS AS STRUCT 
	EventInfo = ReadEvent(c_event, sorted_ch);
	AllBoards_ADC = EventInfo.boards;
	timestamp = EventInfo.timeStamp;

	//// TriggerMask //////
	TriggerMaskSize = EventInfo.TrMask_size;
	TriggerMaskStrips = EventInfo.TrMask_strips;
	TriggerMaskChannels =EventInfo.TrMask_channels;
	  
	// CALCULATE ENERGY DEPOSITS
	for(int b=0; b<nBoards; b++) {
	  view = Views.at(b);
	  Nstation = nStations.at(b);
	  ADCcounts = AllBoards_ADC.at(b);
	  Peds_board = boards_Peds.at(b);
	  OnePhe_board = boards_OnePhes.at(b);
	  
	  ////////////////////// TIME EXPANSIONS ///////////////////////////
	  if (Nstation == 1) {
	    if(view == "x1") Texp_1x1 = EventInfo.timeExp.at(b);
	    if(view == "x2") Texp_1x2 = EventInfo.timeExp.at(b);
	    if(view == "y1") Texp_1y1 = EventInfo.timeExp.at(b);
	    if(view == "y2") Texp_1y2 = EventInfo.timeExp.at(b);
	  }
	  if (Nstation == 2) {
	    if(view == "x1") Texp_2x1 = EventInfo.timeExp.at(b);
	    if(view == "x2") Texp_2x2 = EventInfo.timeExp.at(b);
	    if(view == "y1") Texp_2y1 = EventInfo.timeExp.at(b);
	    if(view == "y2") Texp_2y2 = EventInfo.timeExp.at(b);
	  }
	  if (Nstation == 3) {
	    if(view == "x1") Texp_3x1 = EventInfo.timeExp.at(b);
	    if(view == "x2") Texp_3x2 = EventInfo.timeExp.at(b);
	    if(view == "y1") Texp_3y1 = EventInfo.timeExp.at(b);
	    if(view == "y2") Texp_3y2 = EventInfo.timeExp.at(b);
	  }
	  if (Nstation == 4) {
	    if(view == "x1") Texp_4x1 = EventInfo.timeExp.at(b);
	    if(view == "x2") Texp_4x2 = EventInfo.timeExp.at(b);
	    if(view == "y1") Texp_4y1 = EventInfo.timeExp.at(b);
	    if(view == "y2") Texp_4y2 = EventInfo.timeExp.at(b);
	  }
	  ////////////////////////////////////////////////////////////

	  for(int ADC_ch=0; ADC_ch<nChannels; ADC_ch++) {

	    ADCcount = ADCcounts.at(ADC_ch);
	    Deposit = (ADCcount-Peds_board.at(ADC_ch))/OnePhe_board.at(ADC_ch);

	    //Fill vectors of modules
	    if (Nstation == 1) {
	      if(view == "x1")  Deposits_p1x1.push_back(Deposit);
	      if(view == "x2")  Deposits_p1x2.push_back(Deposit);
	      if(view == "y1")  Deposits_p1y1.push_back(Deposit);
	      if(view == "y2")  Deposits_p1y2.push_back(Deposit);
	    }

	    if (Nstation == 2) {
	      if(view == "x1")  Deposits_p2x1.push_back(Deposit);
	      if(view == "x2")  Deposits_p2x2.push_back(Deposit);
	      if(view == "y1")  Deposits_p2y1.push_back(Deposit);
	      if(view == "y2")  Deposits_p2y2.push_back(Deposit);
	    }

	    if (Nstation == 3) {
	      if(view == "x1")  Deposits_p3x1.push_back(Deposit);
	      if(view == "x2")  Deposits_p3x2.push_back(Deposit);
	      if(view == "y1")  Deposits_p3y1.push_back(Deposit);
	      if(view == "y2")  Deposits_p3y2.push_back(Deposit);
	    }
	    if (Nstation == 4) {
	      if(view == "x1")  Deposits_p4x1.push_back(Deposit);
	      if(view == "x2")  Deposits_p4x2.push_back(Deposit);
	      if(view == "y1")  Deposits_p4y1.push_back(Deposit);
	      if(view == "y2")  Deposits_p4y2.push_back(Deposit);
	    }
	    
	  }
	}

	// Concatenate vectors to form view-plane arrays of deposits (64 strips) 
	Deposits_p1x.insert(Deposits_p1x.end(), Deposits_p1x1.begin(), Deposits_p1x1.end());
	Deposits_p1x.insert(Deposits_p1x.end(), Deposits_p1x2.begin(), Deposits_p1x2.end());
	Deposits_p2x.insert(Deposits_p2x.end(), Deposits_p2x1.begin(), Deposits_p2x1.end());
	Deposits_p2x.insert(Deposits_p2x.end(), Deposits_p2x2.begin(), Deposits_p2x2.end());
	Deposits_p3x.insert(Deposits_p3x.end(), Deposits_p3x1.begin(), Deposits_p3x1.end());
	Deposits_p3x.insert(Deposits_p3x.end(), Deposits_p3x2.begin(), Deposits_p3x2.end());
	Deposits_p4x.insert(Deposits_p4x.end(), Deposits_p4x1.begin(), Deposits_p4x1.end());
	Deposits_p4x.insert(Deposits_p4x.end(), Deposits_p4x2.begin(), Deposits_p4x2.end());

	Deposits_p1y.insert(Deposits_p1y.end(), Deposits_p1y2.begin(), Deposits_p1y2.end());
	Deposits_p1y.insert(Deposits_p1y.end(), Deposits_p1y1.begin(), Deposits_p1y1.end());
	Deposits_p2y.insert(Deposits_p2y.end(), Deposits_p2y2.begin(), Deposits_p2y2.end());
	Deposits_p2y.insert(Deposits_p2y.end(), Deposits_p2y1.begin(), Deposits_p2y1.end());
	Deposits_p3y.insert(Deposits_p3y.end(), Deposits_p3y2.begin(), Deposits_p3y2.end());
	Deposits_p3y.insert(Deposits_p3y.end(), Deposits_p3y1.begin(), Deposits_p3y1.end());
	Deposits_p4y.insert(Deposits_p4y.end(), Deposits_p4y2.begin(), Deposits_p4y2.end());
	Deposits_p4y.insert(Deposits_p4y.end(), Deposits_p4y1.begin(), Deposits_p4y1.end());
	
	/////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////            ///////////////////////////////////////
	///////////////////////////////// CLUSTERING ///////////////////////////////////////
	////////////////////////////////            ///////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////

	
	resuts_p1x = CreateClusterList(Deposits_p1x,s1,s2,s3,Texp_1x1,Texp_1x2,TriggerMaskStrips.at(0),TriggerMaskStrips.at(1));
	resuts_p2x = CreateClusterList(Deposits_p2x,s1,s2,s3,Texp_2x1,Texp_2x2,TriggerMaskStrips.at(4),TriggerMaskStrips.at(5));
	resuts_p3x = CreateClusterList(Deposits_p3x,s1,s2,s3,Texp_3x1,Texp_3x2,TriggerMaskStrips.at(8),TriggerMaskStrips.at(9));
	resuts_p4x = CreateClusterList(Deposits_p4x,s1,s2,s3,Texp_4x1,Texp_4x2,TriggerMaskStrips.at(12),TriggerMaskStrips.at(13));

	resuts_p1y = CreateClusterList(Deposits_p1y,s1,s2,s3,Texp_1y2,Texp_1y1,TriggerMaskStrips.at(3),TriggerMaskStrips.at(2));
	resuts_p2y = CreateClusterList(Deposits_p2y,s1,s2,s3,Texp_2y2,Texp_2y1,TriggerMaskStrips.at(7),TriggerMaskStrips.at(6));
	resuts_p3y = CreateClusterList(Deposits_p3y,s1,s2,s3,Texp_3y2,Texp_3y1,TriggerMaskStrips.at(11),TriggerMaskStrips.at(10));
	resuts_p4y = CreateClusterList(Deposits_p4y,s1,s2,s3,Texp_4y2,Texp_4y1,TriggerMaskStrips.at(15),TriggerMaskStrips.at(14));

	ClusterZ1_Texp = resuts_p1x.TimeExpansions;
	ClusterZ2_Texp = resuts_p2x.TimeExpansions;
	ClusterZ3_Texp = resuts_p3x.TimeExpansions;
	ClusterZ4_Texp = resuts_p4x.TimeExpansions;
	
	ClusterY1_Texp = resuts_p1y.TimeExpansions;
	ClusterY2_Texp = resuts_p2y.TimeExpansions;
	ClusterY3_Texp = resuts_p3y.TimeExpansions;
	ClusterY4_Texp = resuts_p4y.TimeExpansions;
	
	Nclusters_Z1 = resuts_p1x.ClustersEnergy.size();
	Nclusters_Z2 = resuts_p2x.ClustersEnergy.size();
	Nclusters_Z3 = resuts_p3x.ClustersEnergy.size();
	Nclusters_Z4 = resuts_p4x.ClustersEnergy.size();
	
	Nclusters_Y1 = resuts_p1y.ClustersEnergy.size();
	Nclusters_Y2 = resuts_p2y.ClustersEnergy.size();
	Nclusters_Y3 = resuts_p3y.ClustersEnergy.size();
	Nclusters_Y4 = resuts_p4y.ClustersEnergy.size();

	if(Nclusters_Z1>0) Nevents_withX1_cl++;
	if(Nclusters_Z2>0) Nevents_withX2_cl++;
	if(Nclusters_Z3>0) Nevents_withX3_cl++;
	if(Nclusters_Z4>0) Nevents_withX4_cl++;

	if(Nclusters_Y1>0) Nevents_withY1_cl++;
	if(Nclusters_Y2>0) Nevents_withY2_cl++;
	if(Nclusters_Y3>0) Nevents_withY3_cl++;
	if(Nclusters_Y4>0) Nevents_withY4_cl++;
	
	ClusterEnergy_Z1 = resuts_p1x.ClustersEnergy;
	ClusterPosition_Z1 =resuts_p1x.ClustersPositions;
	ClusterSize_Z1 = resuts_p1x.ClustersSize;

	ClusterEnergy_Z2 = resuts_p2x.ClustersEnergy;
	ClusterPosition_Z2 =resuts_p2x.ClustersPositions;
	ClusterSize_Z2 = resuts_p2x.ClustersSize;

	ClusterEnergy_Z3 = resuts_p3x.ClustersEnergy;
	ClusterPosition_Z3 =resuts_p3x.ClustersPositions;
	ClusterSize_Z3 = resuts_p3x.ClustersSize;

	ClusterEnergy_Z4 = resuts_p4x.ClustersEnergy;
	ClusterPosition_Z4 =resuts_p4x.ClustersPositions;
	ClusterSize_Z4 = resuts_p4x.ClustersSize;

	ClusterEnergy_Y1 = resuts_p1y.ClustersEnergy;
	ClusterPosition_Y1 =resuts_p1y.ClustersPositions;
	ClusterSize_Y1 = resuts_p1y.ClustersSize;

	ClusterEnergy_Y2 = resuts_p2y.ClustersEnergy;
	ClusterPosition_Y2 =resuts_p2y.ClustersPositions;
	ClusterSize_Y2 = resuts_p2y.ClustersSize;

	ClusterEnergy_Y3 = resuts_p3y.ClustersEnergy;
	ClusterPosition_Y3 =resuts_p3y.ClustersPositions;
	ClusterSize_Y3 = resuts_p3y.ClustersSize;

	ClusterEnergy_Y4 = resuts_p4y.ClustersEnergy;
	ClusterPosition_Y4 =resuts_p4y.ClustersPositions;
	ClusterSize_Y4 = resuts_p4y.ClustersSize;

	StripsEnergy_Z1 = resuts_p1x.StripsEnergy;
	StripsEnergy_Z2 = resuts_p2x.StripsEnergy;
	StripsEnergy_Z3 = resuts_p3x.StripsEnergy;
	StripsEnergy_Z4 = resuts_p4x.StripsEnergy;

	StripsEnergy_Y1 = resuts_p1y.StripsEnergy;
	StripsEnergy_Y2 = resuts_p2y.StripsEnergy;
	StripsEnergy_Y3 = resuts_p3y.StripsEnergy;
	StripsEnergy_Y4 = resuts_p4y.StripsEnergy;

	StripsPosition_Z1 = resuts_p1x.StripsPositions;
	StripsPosition_Z2 = resuts_p2x.StripsPositions;
	StripsPosition_Z3 = resuts_p3x.StripsPositions;
	StripsPosition_Z4 = resuts_p4x.StripsPositions;

	StripsPosition_Y1 = resuts_p1y.StripsPositions;
	StripsPosition_Y2 = resuts_p2y.StripsPositions;
	StripsPosition_Y3 = resuts_p3y.StripsPositions;
	StripsPosition_Y4 = resuts_p4y.StripsPositions;

	StripsID_Z1 = resuts_p1x.StripsID;
	StripsID_Z2 = resuts_p2x.StripsID;
	StripsID_Z3 = resuts_p3x.StripsID;
	StripsID_Z4 = resuts_p4x.StripsID;

	StripsID_Y1 = resuts_p1y.StripsID;
	StripsID_Y2 = resuts_p2y.StripsID;
	StripsID_Y3 = resuts_p3y.StripsID;
	StripsID_Y4 = resuts_p4y.StripsID;
	
	//////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////            ///////////////////////////////////////
	/////////////////////////////////  TRACKING  ////////////////////////////////////////
	////////////////////////////////            ////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	
	tracks_xz = MakeTracks(ClusterPosition_Z1,ClusterPosition_Z2,ClusterPosition_Z3,ClusterPosition_Z4,ClusterEnergy_Z1,ClusterEnergy_Z2,ClusterEnergy_Z3,ClusterEnergy_Z4,ClusterZ1_Texp, ClusterZ2_Texp, ClusterZ3_Texp, ClusterZ4_Texp,proximity_cut_xz,X_pos,Z_add,sigma_z);
	tracks_xy = MakeTracks(ClusterPosition_Y1,ClusterPosition_Y2,ClusterPosition_Y3,ClusterPosition_Y4,ClusterEnergy_Y1,ClusterEnergy_Y2,ClusterEnergy_Y3,ClusterEnergy_Y4,ClusterY1_Texp, ClusterY2_Texp, ClusterY3_Texp,ClusterY4_Texp,proximity_cut_xy,X_pos,Y_add,sigma_y);
	// 3 PLANES TRACKS ----> VARIABLES //////////////
	Intercept_3p_xz = tracks_xz.intercepts_3p;
	Slope_3p_xz =tracks_xz.slopes_3p;
	chiSquare_3p_xz=tracks_xz.chiSquares_3p;
	TrackCluster_z1_Position=tracks_xz.position_c1;
	TrackCluster_z2_Position=tracks_xz.position_c2;
	TrackCluster_z3_Position=tracks_xz.position_c3;
	TrackCluster_z1_index =tracks_xz.cluster_index_1;
	TrackCluster_z2_index =tracks_xz.cluster_index_2;
	TrackCluster_z3_index =tracks_xz.cluster_index_3;
	Intercept_3p_xy = tracks_xy.intercepts_3p;
	Slope_3p_xy =tracks_xy.slopes_3p;
	Ntracks_3p_xz = Intercept_3p_xz.size();
	Ntracks_3p_xy = Intercept_3p_xy.size();
	Tracks3p_ExpectedRes_z2 = tracks_xz.Track_3p_ExpectedRes_p2;
	Tracks3p_ExpectedRes_y2 = tracks_xy.Track_3p_ExpectedRes_p2;
	// CALCULATE THE TOTAL NUMBER OF TRACKS IN THE RUN
	N_tracks_xy_3p+=Ntracks_3p_xy;
	N_tracks_xz_3p+=Ntracks_3p_xz;
	N_tracks_3p+=(Ntracks_3p_xy*Ntracks_3p_xz);
	//////////////////////////////////////////////////


	//////////////// INFO ABOUT CLUSTERS OF THE TRACKS ////////////////////////////////////////
	///////////////// Track Clusters time expansions /////////////////////////////////////////
	TrackCluster_z1_Texp.clear();
	TrackCluster_z2_Texp.clear();
	TrackCluster_z3_Texp.clear();

	for (int tr=0; tr<Ntracks_3p_xz; tr++) {
	  
	  TrackCluster_z1_Texp.push_back(ClusterZ1_Texp.at(tracks_xz.cluster_index_1.at(tr)));
	  TrackCluster_z2_Texp.push_back(ClusterZ2_Texp.at(tracks_xz.cluster_index_2.at(tr)));
	  TrackCluster_z3_Texp.push_back(ClusterZ3_Texp.at(tracks_xz.cluster_index_3.at(tr)));
	}

	TrackCluster_y1_Texp.clear();
        TrackCluster_y2_Texp.clear();
        TrackCluster_y3_Texp.clear();
	  
	for (int tr=0; tr<Ntracks_3p_xy; tr++) {
	  
	  TrackCluster_y1_Texp.push_back(ClusterY1_Texp.at(tracks_xy.cluster_index_1.at(tr)));
	  TrackCluster_y2_Texp.push_back(ClusterY2_Texp.at(tracks_xy.cluster_index_2.at(tr)));
	  TrackCluster_y3_Texp.push_back(ClusterY3_Texp.at(tracks_xy.cluster_index_3.at(tr)));
	}

	/////////////////////////////////////////////////////////////////////////////////////////

	///////////////// Track Clusters Size  ////////////////////////////////////////////////
	TrackCluster_z1_Size.clear();
	TrackCluster_z2_Size.clear();
	TrackCluster_z3_Size.clear();

	for (int tr=0; tr<Ntracks_3p_xz; tr++) {
	  
	  TrackCluster_z1_Size.push_back(ClusterSize_Z1.at(tracks_xz.cluster_index_1.at(tr)));
	  TrackCluster_z2_Size.push_back(ClusterSize_Z2.at(tracks_xz.cluster_index_2.at(tr)));
	  TrackCluster_z3_Size.push_back(ClusterSize_Z3.at(tracks_xz.cluster_index_3.at(tr)));
	}

	TrackCluster_y1_Size.clear();
        TrackCluster_y2_Size.clear();
        TrackCluster_y3_Size.clear();
	  
	for (int tr=0; tr<Ntracks_3p_xy; tr++) {
	  
	  TrackCluster_y1_Size.push_back(ClusterSize_Y1.at(tracks_xy.cluster_index_1.at(tr)));
	  TrackCluster_y2_Size.push_back(ClusterSize_Y2.at(tracks_xy.cluster_index_2.at(tr)));
	  TrackCluster_y3_Size.push_back(ClusterSize_Y3.at(tracks_xy.cluster_index_3.at(tr)));
	}

	/////////////////////////////////////////////////////////////////////////////////////////

	
	chiSquare_3p_xy=tracks_xy.chiSquares_3p;
	TrackCluster_y1_Position=tracks_xy.position_c1;
	TrackCluster_y2_Position=tracks_xy.position_c2;
	TrackCluster_y3_Position=tracks_xy.position_c3;
	TrackCluster_y1_index =tracks_xy.cluster_index_1;
	TrackCluster_y2_index =tracks_xy.cluster_index_2;
	TrackCluster_y3_index =tracks_xy.cluster_index_3;
	Residue_Track3p_z1 = tracks_xz.residue_c1;
	Residue_Track3p_z2 = tracks_xz.residue_c2;
	Residue_Track3p_z3 = tracks_xz.residue_c3;
	Residue_Track3p_y1 = tracks_xy.residue_c1;
	Residue_Track3p_y2 = tracks_xy.residue_c2;
	Residue_Track3p_y3 = tracks_xy.residue_c3;
	Ntracks_4p_xy = tracks_xy.Ntracks_4p;
	Ntracks_4p_xz = tracks_xz.Ntracks_4p;
	TrackEnergy_3p_xy = tracks_xy.TrackEnergy_3p;
	TrackEnergy_3p_xz = tracks_xz.TrackEnergy_3p;
	BestChi_xy =  tracks_xy.BestChi;
	BestEnergy_xy =  tracks_xy.BestEnergy;
	BestChi_xz =  tracks_xz.BestChi;
	BestEnergy_xz =  tracks_xz.BestEnergy;
	BestChi_xy_index =  tracks_xy.BestChi_index;
	BestEnergy_xy_index =  tracks_xy.BestEnergy_index;
	BestChi_xz_index =  tracks_xz.BestChi_index;
	BestEnergy_xz_index =  tracks_xz.BestEnergy_index;
	ExpectedPosition_OnPlane4th_xy = tracks_xy.ExpectedPosition_OnPlane4th;
	ExpectedPosition_OnPlane4th_xz = tracks_xz.ExpectedPosition_OnPlane4th;
	Track_3p_to_4p_index_xy = tracks_xy.Track_3p_to_4p_index;
	Track_3p_to_4p_index_xz = tracks_xz.Track_3p_to_4p_index;
	
	/////////////////////////// 3D BEST TRACK ////////////////////////////////////////
	if(Ntracks_3p_xz>0 && Ntracks_3p_xy>0) {
	  
	  if(BestChi_xy_index==BestEnergy_xy_index) BestTrack_xy_ind = BestEnergy_xy_index;
	  else {
	    //// CRITERIA TO CHOOSE THE BEST TRACK BETWEEN BEST CHI AND BEST ENERGY ///////
	    BestTrack_xy_ind =BestChi_xy_index; 
	  }	
	  if(BestChi_xz_index==BestEnergy_xz_index) BestTrack_xz_ind = BestEnergy_xz_index;
	  else {
	    //// CRITERIA TO CHOOSE THE BEST TRACK BETWEEN BEST CHI AND BEST ENERGY ///////
	    BestTrack_xz_ind =BestChi_xz_index; 
	  }
	  if(BestTrack_xz_ind>=0 && BestTrack_xy_ind >=0) {
	    Coordinates_3p = TrackAngularCoordinates(Slope_3p_xy.at(BestTrack_xy_ind),Slope_3p_xz.at(BestTrack_xz_ind),X_pos[0],X_pos[2]);
	    theta_3p = Coordinates_3p.at(0);
	    phi_3p =  Coordinates_3p.at(1);
	    NGooDTracks3p++;
	    
	  BestTrack_3p_xy_index= BestTrack_xy_ind;
	  BestTrack_3p_xz_index= BestTrack_xz_ind;
	  BestTrack_3p_Energy_xy = TrackEnergy_3p_xy.at(BestTrack_xy_ind);
	  BestTrack_3p_Energy_xz = TrackEnergy_3p_xz.at(BestTrack_xz_ind);
	  BestTrack_3p_ChiSquare_xy = chiSquare_3p_xy.at(BestTrack_xy_ind);
	  BestTrack_3p_ChiSquare_xz = chiSquare_3p_xz.at(BestTrack_xz_ind);
	  Residue_BestTracks3p_y1 = Residue_Track3p_y1.at(BestTrack_xy_ind);
	  Residue_BestTracks3p_y2= Residue_Track3p_y2.at(BestTrack_xy_ind);
	  Residue_BestTracks3p_y3= Residue_Track3p_y3.at(BestTrack_xy_ind);
	  Residue_BestTracks3p_z1 = Residue_Track3p_z1.at(BestTrack_xz_ind);
	  Residue_BestTracks3p_z2= Residue_Track3p_z2.at(BestTrack_xz_ind);
	  Residue_BestTracks3p_z3= Residue_Track3p_z3.at(BestTrack_xz_ind);
	  for(int nT_xz=0; nT_xz<Ntracks_3p_xz;nT_xz++) {
	    if(nT_xz==BestTrack_xz_ind) isBestTrack_3p_xz.push_back(1);
	    else isBestTrack_3p_xz.push_back(0);
	  }
	  for(int nT_xy=0; nT_xy<Ntracks_3p_xy;nT_xy++) {
	      if(nT_xy ==BestTrack_xy_ind) isBestTrack_3p_xy.push_back(1);
	      else isBestTrack_3p_xy.push_back(0);
	    }
	  }
	  else {
	    for(int nT_xz=0; nT_xz<Ntracks_3p_xz;nT_xz++)  isBestTrack_3p_xz.push_back(0);
	    for(int nT_xy=0; nT_xy<Ntracks_3p_xy;nT_xy++)  isBestTrack_3p_xy.push_back(0);
	    BestTrack_3p_ChiSquare_xy=-1;
	    BestTrack_3p_ChiSquare_xz=-1;
	    Residue_BestTracks3p_y1=-1;
	    Residue_BestTracks3p_y2=-1;
	    Residue_BestTracks3p_y3=-1;
	    Residue_BestTracks3p_z1=-1;
	    Residue_BestTracks3p_z2=-1;
	    Residue_BestTracks3p_z3=-1;
	    BestTrack_3p_xy_index=-1;
	    BestTrack_3p_xz_index=-1;
	    theta_3p=-1;
	    phi_3p=-1;
	    BestTrack_3p_Energy_xz=-1;
	    BestTrack_3p_Energy_xy=-1;
	    
	  }
	}
	else {
	  BestTrack_3p_ChiSquare_xy=-1;
	  BestTrack_3p_ChiSquare_xz=-1;
	  Residue_BestTracks3p_y1=-1;
	  Residue_BestTracks3p_y2=-1;
	  Residue_BestTracks3p_y3=-1;
	  Residue_BestTracks3p_z1=-1;
	  Residue_BestTracks3p_z2=-1;
	  Residue_BestTracks3p_z3=-1;
	  BestTrack_3p_xy_index=-1;
	  BestTrack_3p_xz_index=-1;
	  theta_3p=-1;
	  phi_3p=-1;
	  BestTrack_3p_Energy_xz=-1;
	  BestTrack_3p_Energy_xy=-1;
	  
	}
	// CLUSTER IS IN TRACK
	isInTrack_3p_clZ1 = tracks_xz.IsInTrack_clusters1;
	isInTrack_3p_clZ2 = tracks_xz.IsInTrack_clusters2;
	isInTrack_3p_clZ3 = tracks_xz.IsInTrack_clusters3;
	isInTrack_3p_clY1 = tracks_xy.IsInTrack_clusters1;
	isInTrack_3p_clY2 = tracks_xy.IsInTrack_clusters2;
	isInTrack_3p_clY3 = tracks_xy.IsInTrack_clusters3;

	isInTrack_4p_clZ1 = tracks_xz.IsInTrack_clusters1_4p;
	isInTrack_4p_clZ2 = tracks_xz.IsInTrack_clusters2_4p;
	isInTrack_4p_clZ3 = tracks_xz.IsInTrack_clusters3_4p;
	isInTrack_4p_clY1 = tracks_xy.IsInTrack_clusters1_4p;
	isInTrack_4p_clY2 = tracks_xy.IsInTrack_clusters2_4p;
	isInTrack_4p_clY3 = tracks_xy.IsInTrack_clusters3_4p;
	
	/////////////////////////////////////////////////

	// 4 PLANES TRACKING ----> VARIABLES ///////////
	Plane4th_isIntercepted_xz=tracks_xz.Plane4th_isIntercepted;
	Plane4th_isIntercepted_xy=tracks_xy.Plane4th_isIntercepted;
	TrackCluster_z4_Position=tracks_xz.positions_c4;
	TrackCluster_y4_Position=tracks_xy.positions_c4;
	Intercept_4p_xz =tracks_xz.intercept_4p;
	Slope_4p_xz = tracks_xz.slope_4p;	
	chiSquare_4p_xz =tracks_xz.chiSquares_4p;
	displacement_p4_xz =tracks_xz.displacement_p4;
	residue_c1_p4_xz = tracks_xz.residue_c1_p4;
	residue_c2_p4_xz = tracks_xz.residue_c2_p4;
	residue_c3_p4_xz = tracks_xz.residue_c3_p4;
	residue_c4_p4_xz = tracks_xz.residue_c4_p4;
	cluster_c4_index_xz = tracks_xz.cluster_indices_4;

	cluster_c4_Texp_xz.clear();
	cluster_c4_Texp_xy.clear();
	// Tiime expansions of 4 planes tracks (cluster of the 4th plane)
	for(int k=0; k<cluster_c4_index_xz.size(); k++) {
	  clusters4_Texp_single.clear();
	  for(int j=0; j<cluster_c4_index_xz.at(k).size(); j++) {
	    clusters4_Texp_single.push_back(ClusterZ4_Texp.at(cluster_c4_index_xz.at(k).at(j)));
	  }
	  cluster_c4_Texp_xz.push_back(clusters4_Texp_single);
	}
	cluster_c4_index_xy = tracks_xy.cluster_indices_4;

	for(int k=0; k<tracks_xy.cluster_indices_4.size(); k++) {
	  clusters4_Texp_single.clear();
	  for(int j=0; j<tracks_xy.cluster_indices_4.at(k).size(); j++) {
	    clusters4_Texp_single.push_back(ClusterY4_Texp.at(tracks_xy.cluster_indices_4.at(k).at(j)));
	  }
	  cluster_c4_Texp_xy.push_back(clusters4_Texp_single);
	}
	/////////////////////////////////////////////////////////////////////
       
	Intercept_4p_xy =tracks_xy.intercept_4p;
	Slope_4p_xy = tracks_xy.slope_4p;	
	chiSquare_4p_xy =tracks_xy.chiSquares_4p;
	displacement_p4_xy =tracks_xy.displacement_p4;
	residue_c1_p4_xy = tracks_xy.residue_c1_p4;
	residue_c2_p4_xy = tracks_xy.residue_c2_p4;
	residue_c3_p4_xy = tracks_xy.residue_c3_p4;
	residue_c4_p4_xy = tracks_xy.residue_c4_p4;
	TrackEnergy_4p_xy = tracks_xy.TrackEnergy_4p;
	TrackEnergy_4p_xz = tracks_xz.TrackEnergy_4p;
	ScatteringAngle_xy =tracks_xy.ScatteringAngles;
	ScatteringAngle_xz =tracks_xz.ScatteringAngles;

	////////////7////////////////////////////////////////////////////////////////////
	////////////////////////////// BEST TRACK 4 PLANES /////////////////////////////
	
	if(Track_3p_to_4p_index_xz.size()>0 &&  Track_3p_to_4p_index_xy.size()>0) {

	  track_4p_index_xy = find(Track_3p_to_4p_index_xy.begin(), Track_3p_to_4p_index_xy.end(),BestTrack_xy_ind);
	  track_4p_index_xz = find(Track_3p_to_4p_index_xz.begin(), Track_3p_to_4p_index_xz.end(),BestTrack_xz_ind);
	  if(track_4p_index_xy!= Track_3p_to_4p_index_xy.end() && track_4p_index_xz!= Track_3p_to_4p_index_xz.end()) {
	    index_4p_xz = std::distance(Track_3p_to_4p_index_xz.begin(), track_4p_index_xz);
	    index_4p_xy = std::distance(Track_3p_to_4p_index_xy.begin(), track_4p_index_xy);
	    Track_3p_of_4p_index_xy=index_4p_xy;
	    Track_3p_of_4p_index_xz=index_4p_xz;
	    
	    for(int k =0; k<Intercept_4p_xy.at(index_4p_xy).size(); k++) {
	      if(cluster_c4_Texp_xy.at(index_4p_xy).at(k) >0 )
		{
		  Best_track_4p_isTexpNULL_xy=0;
		  bestTracks_4p_index_xy  = k;
		  break;
		}else {
		Best_track_4p_isTexpNULL_xy=1;
		bestTracks_4p_index_xy =0;
		continue;
	      }
	    }
	    Track_4p_index_xy=bestTracks_4p_index_xy;
	    for(int k =0; k<Intercept_4p_xz.at(index_4p_xz).size(); k++) {
	      if(cluster_c4_Texp_xz.at(index_4p_xz).at(k) >0 )
		{
		  Best_track_4p_isTexpNULL_xz=0;
		  bestTracks_4p_index_xz  = k;
		  break;
		}else {
		bestTracks_4p_index_xz =0;
		Best_track_4p_isTexpNULL_xz=1;
		continue;
	      }
	      
	    }
	    Residue_BestTracks4p_y1=residue_c1_p4_xy.at(index_4p_xy).at(bestTracks_4p_index_xy);
	    Residue_BestTracks4p_y2=residue_c2_p4_xy.at(index_4p_xy).at(bestTracks_4p_index_xy);
	    Residue_BestTracks4p_y3=residue_c3_p4_xy.at(index_4p_xy).at(bestTracks_4p_index_xy);
	    Residue_BestTracks4p_y4=residue_c4_p4_xy.at(index_4p_xy).at(bestTracks_4p_index_xy);
	    
	    Residue_BestTracks4p_z1=residue_c1_p4_xz.at(index_4p_xz).at(bestTracks_4p_index_xz);
	    Residue_BestTracks4p_z2=residue_c2_p4_xz.at(index_4p_xz).at(bestTracks_4p_index_xz);
	    Residue_BestTracks4p_z3=residue_c3_p4_xz.at(index_4p_xz).at(bestTracks_4p_index_xz);
	    Residue_BestTracks4p_z4=residue_c4_p4_xz.at(index_4p_xz).at(bestTracks_4p_index_xz);
	    BestTracks_ScatteringAngle_xy = ScatteringAngle_xy.at(index_4p_xy).at(bestTracks_4p_index_xy);
	    BestTracks_ScatteringAngle_xz = ScatteringAngle_xz.at(index_4p_xz).at(bestTracks_4p_index_xz);
	    BestTrack_4p_ChiSquare_xy = chiSquare_4p_xy.at(index_4p_xy).at(bestTracks_4p_index_xy);
	    BestTrack_4p_ChiSquare_xz = chiSquare_4p_xz.at(index_4p_xz).at(bestTracks_4p_index_xz);
	    Track_4p_index_xz=bestTracks_4p_index_xz;
	    Coordinates_4p = TrackAngularCoordinates(Slope_4p_xy.at(index_4p_xy).at(bestTracks_4p_index_xy),Slope_4p_xz.at(index_4p_xz).at(bestTracks_4p_index_xz),X_pos[0],X_pos[3]);
	    theta_4p = Coordinates_4p.at(0);
	    phi_4p =  Coordinates_4p.at(1);
	    NGooDTracks4p++;
	  }else {
          Residue_BestTracks4p_y1=-1;
          Residue_BestTracks4p_y2=-1;
          Residue_BestTracks4p_y3=-1;
          Residue_BestTracks4p_y4=-1;
          Residue_BestTracks4p_z1=-1;
          Residue_BestTracks4p_z2=-1;
          Residue_BestTracks4p_z3=-1;
          Residue_BestTracks4p_z4=-1;

          BestTrack_4p_ChiSquare_xy=-1;
          BestTrack_4p_ChiSquare_xz=-1;
          BestTracks_ScatteringAngle_xy=-1;
          BestTracks_ScatteringAngle_xz=-1;
          Best_track_4p_isTexpNULL_xz=-1;
          Best_track_4p_isTexpNULL_xy=-1;
          Track_4p_index_xz=-1;
          Track_4p_index_xy=-1;
          Track_3p_of_4p_index_xy=-1;
          Track_3p_of_4p_index_xz=-1;
          theta_4p =-1;
          phi_4p=-1;
	  }

	}else {
	  Residue_BestTracks4p_y1=-1;
	  Residue_BestTracks4p_y2=-1;
	  Residue_BestTracks4p_y3=-1;
	  Residue_BestTracks4p_y4=-1;
	  Residue_BestTracks4p_z1=-1;
	  Residue_BestTracks4p_z2=-1;
	  Residue_BestTracks4p_z3=-1;
	  Residue_BestTracks4p_z4=-1;
	  
	  BestTrack_4p_ChiSquare_xy=-1;
	  BestTrack_4p_ChiSquare_xz=-1;
	  BestTracks_ScatteringAngle_xy=-1;
	  BestTracks_ScatteringAngle_xz=-1;
	  Best_track_4p_isTexpNULL_xz=-1;
	  Best_track_4p_isTexpNULL_xy=-1;
	  Track_4p_index_xz=-1;
	  Track_4p_index_xy=-1;
	  Track_3p_of_4p_index_xy=-1;
	  Track_3p_of_4p_index_xz=-1;
	  theta_4p =-1;
	  phi_4p=-1;
	}
	  n4p_xy=0;
	  n4p_xz=0;
	if(Ntracks_3p_xz>0 && Ntracks_3p_xy>0) Nevents_withAtrack3p+=1;
	if(Ntracks_3p_xz>0) Nevents_withAtrack3p_xz++;
	if(Ntracks_3p_xy>0) Nevents_withAtrack3p_xy++;

	for(int n=0; n<Ntracks_4p_xy.size(); n++) {
	  if(Ntracks_4p_xy.at(n) >0 ) n4p_xy++;
	}

	for(int n=0; n<Ntracks_4p_xz.size(); n++) {
	  if(Ntracks_4p_xz.at(n) >0 ) n4p_xz++;
	}
	if(n4p_xz>0 && n4p_xy>0) {
	  Nevents_withAtrack4p+=1;
	  N_tracks_xy_4p+=n4p_xy;
	  N_tracks_xz_4p+=n4p_xz;
	  N_tracks_4p+=(n4p_xz*n4p_xy);
	}
	tree->Fill();
	isBestTrack_3p_xy.clear();
	isBestTrack_3p_xz.clear();
	ClusterPosition_Z1.clear();
	ClusterEnergy_Z1.clear();
	ClusterPosition_Z2.clear();
	ClusterEnergy_Z2.clear();
	ClusterPosition_Z3.clear();
	ClusterEnergy_Z3.clear();
	ClusterPosition_Z4.clear();
	ClusterEnergy_Z4.clear();

	ClusterPosition_Y1.clear();
	ClusterEnergy_Y1.clear();
	ClusterPosition_Y2.clear();
	ClusterEnergy_Y2.clear();
	ClusterPosition_Y3.clear();
	ClusterEnergy_Y3.clear();
	ClusterPosition_Y4.clear();
	ClusterEnergy_Y4.clear();

	ClusterSize_Z1.clear();
	ClusterSize_Z2.clear();
	ClusterSize_Z3.clear();
	ClusterSize_Z4.clear();
	
	ClusterSize_Y1.clear();
	ClusterSize_Y2.clear();
	ClusterSize_Y3.clear();
	ClusterSize_Y4.clear();
	
	Deposits_p1x1.clear();
	Deposits_p1x2.clear();
	Deposits_p2x1.clear();
	Deposits_p2x2.clear();
	Deposits_p3x1.clear();
	Deposits_p3x2.clear();
	Deposits_p4x1.clear();
	Deposits_p4x2.clear();
	Deposits_p1y1.clear();
	Deposits_p1y2.clear();
	Deposits_p2y1.clear();
	Deposits_p2y2.clear();
	Deposits_p3y1.clear();
	Deposits_p3y2.clear();
	Deposits_p4y1.clear();
	Deposits_p4y2.clear();

	Deposits_p1x.clear();
	Deposits_p2x.clear();
	Deposits_p3x.clear();
	Deposits_p4x.clear();

	Deposits_p1y.clear();
	Deposits_p2y.clear();
	Deposits_p3y.clear();
	Deposits_p4y.clear(); 
      
      }
  

// CALCULATE THE DURATION OF THE RUN ///////////////////////////////////////////////
  //CPyObject RunDurationModuleName =PyUnicode_FromString("RunDuration");
  //CPyObject RunDurationModule =PyImport_Import(RunDurationModuleName);
  //CPyObject RunnTime_function =PyObject_GetAttrString(RunDurationModule, "runTime");
  //PyObject *RunTime_args =PyTuple_New(2);
  //PyTuple_SetItem(RunTime_args, 0, PyUnicode_FromString(first_line));
  //PyTuple_SetItem(RunTime_args, 1, PyUnicode_FromString(last_line));
  //RunDuration =PyFloat_AsDouble(PyObject_CallObject(RunnTime_function,RunTime_args));
  //cout << "Run Duration: "<<  RunDuration<<" minutes"<<endl;
  ////////////////////////////////////////////////////////////////////////////////////

  cout << "N events with at least a 3p track: " << Nevents_withAtrack3p << endl;
  cout << "N events with at least a 4p track: " << Nevents_withAtrack4p << endl;
  cout << "N events with at least a 3p track on the plane xy: " << Nevents_withAtrack3p_xy << endl;
  cout << "N events with at least a 3p track on the plane xz: " << Nevents_withAtrack3p_xz << endl;
  cout << "Events with at least 1 cluster in Z1: " << Nevents_withX1_cl << endl;
  cout << "Events with at least 1 cluster in Z2: " << Nevents_withX2_cl << endl;
  cout << "Events with at least 1 cluster in Z3: " << Nevents_withX3_cl << endl;
  cout << "Events with at least 1 cluster in Z4: " << Nevents_withX4_cl << endl;
  cout << "Events with at least 1 cluster in Y1: " << Nevents_withY1_cl << endl;
  cout << "Events with at least 1 cluster in Y2: " << Nevents_withY2_cl << endl;
  cout << "Events with at least 1 cluster in Y3: " << Nevents_withY3_cl << endl;
  cout << "Events with at least 1 cluster in Y4: " << Nevents_withY4_cl << endl;
  TBranch *b_RunDur =   tree->Branch("RunDuration",&RunDuration);
  TBranch *b_Nevents =   tree->Branch("Nevents",&ev);

  for(int entr=0; entr<ev; entr++) {
  b_RunDur->Fill();
  b_Nevents->Fill();
  }


/////////// PARAMETERS TO CONTROL THE ANALYSIS PERFORMANCIES --> FILLED IN THE MINITREE //////////
  
  double Mean_clusterEnergy_Z1,Mean_clusterEnergy_Z2,Mean_clusterEnergy_Z3,Mean_clusterEnergy_Z4;
  double RMS_clusterEnergy_Z1,RMS_clusterEnergy_Z2,RMS_clusterEnergy_Z3,RMS_clusterEnergy_Z4;
  TH1F *h_EnergyZ1 = new TH1F("h_EnergyZ1","",1000,0,1000);
  TH1F *h_EnergyZ2 = new TH1F("h_EnergyZ2","",1000,0,1000);
  TH1F *h_EnergyZ3 = new TH1F("h_EnergyZ3","",1000,0,1000);
  TH1F *h_EnergyZ4 = new TH1F("h_EnergyZ4","",1000,0,1000);
  tree->Draw("ClusterEnergy_Z1>>h_EnergyZ1","","goff");
  Mean_clusterEnergy_Z1= h_EnergyZ1->GetMean();
  RMS_clusterEnergy_Z1= h_EnergyZ1->GetRMS();
  h_EnergyZ1->Delete();
  tree->Draw("ClusterEnergy_Z2>>h_EnergyZ2","","goff");
  Mean_clusterEnergy_Z2= h_EnergyZ2->GetMean();
  RMS_clusterEnergy_Z2= h_EnergyZ2->GetRMS();
  h_EnergyZ2->Delete();
  tree->Draw("ClusterEnergy_Z3>>h_EnergyZ3","","goff");
  Mean_clusterEnergy_Z3= h_EnergyZ3->GetMean();
  RMS_clusterEnergy_Z3= h_EnergyZ3->GetRMS();
  h_EnergyZ3->Delete();
  tree->Draw("ClusterEnergy_Z4>>h_EnergyZ4","","goff");
  Mean_clusterEnergy_Z4= h_EnergyZ4->GetMean();
  RMS_clusterEnergy_Z4= h_EnergyZ4->GetRMS();
  h_EnergyZ4->Delete();

  double Mean_clusterEnergy_Y1,Mean_clusterEnergy_Y2,Mean_clusterEnergy_Y3,Mean_clusterEnergy_Y4;
  double RMS_clusterEnergy_Y1,RMS_clusterEnergy_Y2,RMS_clusterEnergy_Y3,RMS_clusterEnergy_Y4;
  TH1F *h_EnergyY1 = new TH1F("h_EnergyY1","",1000,0,1000);
  TH1F *h_EnergyY2 = new TH1F("h_EnergyY2","",1000,0,1000);
  TH1F *h_EnergyY3 = new TH1F("h_EnergyY3","",1000,0,1000);
  TH1F *h_EnergyY4 = new TH1F("h_EnergyY4","",1000,0,1000);
  tree->Draw("ClusterEnergy_Y1>>h_EnergyY1","","goff");
  Mean_clusterEnergy_Y1= h_EnergyY1->GetMean();
  RMS_clusterEnergy_Y1= h_EnergyY1->GetRMS();
  h_EnergyY1->Delete();
  tree->Draw("ClusterEnergy_Y2>>h_EnergyY2","","goff");
  Mean_clusterEnergy_Y2= h_EnergyY2->GetMean();
  RMS_clusterEnergy_Y2= h_EnergyY2->GetRMS();
  h_EnergyY2->Delete();
  tree->Draw("ClusterEnergy_Y3>>h_EnergyY3","","goff");
  Mean_clusterEnergy_Y3= h_EnergyY3->GetMean();
  RMS_clusterEnergy_Y3= h_EnergyY3->GetRMS();
  h_EnergyY3->Delete();
  tree->Draw("ClusterEnergy_Y4>>h_EnergyY4","","goff");
  Mean_clusterEnergy_Y4= h_EnergyY4->GetMean();
  RMS_clusterEnergy_Y4= h_EnergyY4->GetRMS();
  h_EnergyY4->Delete();

  double Mean_NclustersZ1, RMS_NclustersZ1,  Mean_NclustersZ2, RMS_NclustersZ2, Mean_NclustersZ3, RMS_NclustersZ3, Mean_NclustersZ4, RMS_NclustersZ4;
  double Mean_NclustersY1, RMS_NclustersY1,  Mean_NclustersY2, RMS_NclustersY2, Mean_NclustersY3, RMS_NclustersY3, Mean_NclustersY4, RMS_NclustersY4;

  TH1F *h_NclustersZ1 = new TH1F("h_NclustersZ1","",10,0,10);
  tree->Draw("Nclusters_Z1>>h_NclustersZ1","","goff");
  Mean_NclustersZ1=h_NclustersZ1->GetMean();
  RMS_NclustersZ1=h_NclustersZ1->GetRMS();
  h_NclustersZ1->Delete();

  TH1F *h_NclustersZ2 = new TH1F("h_NclustersZ2","",10,0,10);
  tree->Draw("Nclusters_Z2>>h_NclustersZ2","","goff");
  Mean_NclustersZ2=h_NclustersZ2->GetMean();
  RMS_NclustersZ2=h_NclustersZ2->GetRMS();
  h_NclustersZ2->Delete();

  TH1F *h_NclustersZ3 = new TH1F("h_NclustersZ3","",10,0,10);
  tree->Draw("Nclusters_Z3>>h_NclustersZ3","","goff");
  Mean_NclustersZ3=h_NclustersZ3->GetMean();
  RMS_NclustersZ3=h_NclustersZ3->GetRMS();
  h_NclustersZ3->Delete();

  TH1F *h_NclustersZ4 = new TH1F("h_NclustersZ4","",10,0,10);
  tree->Draw("Nclusters_Z4>>h_NclustersZ4","","goff");
  Mean_NclustersZ4=h_NclustersZ4->GetMean();
  RMS_NclustersZ4=h_NclustersZ4->GetRMS();
  h_NclustersZ4->Delete();

  TH1F *h_NclustersY1 = new TH1F("h_NclustersY1","",10,0,10);
  tree->Draw("Nclusters_Y1>>h_NclustersY1","","goff");
  Mean_NclustersY1=h_NclustersY1->GetMean();
  RMS_NclustersY1=h_NclustersY1->GetRMS();
  h_NclustersY1->Delete();

  TH1F *h_NclustersY2 = new TH1F("h_NclustersY2","",10,0,10);
  tree->Draw("Nclusters_Y2>>h_NclustersY2","","goff");
  Mean_NclustersY2=h_NclustersY2->GetMean();
  RMS_NclustersY2=h_NclustersY2->GetRMS();
  h_NclustersY2->Delete();

  TH1F *h_NclustersY3 = new TH1F("h_NclustersY3","",10,0,10);
  tree->Draw("Nclusters_Y3>>h_NclustersY3","","goff");
  Mean_NclustersY3=h_NclustersY3->GetMean();
  RMS_NclustersY3=h_NclustersY3->GetRMS();
  h_NclustersY3->Delete();

  TH1F *h_NclustersY4 = new TH1F("h_NclustersY4","",10,0,10);
  tree->Draw("Nclusters_Y4>>h_NclustersY4","","goff");
  Mean_NclustersY4=h_NclustersY4->GetMean();
  RMS_NclustersY4=h_NclustersY4->GetRMS();
  h_NclustersY4->Delete();
  
  tree->Write();
  fileROOT->Close();

  double Perc_withAtrack3p, Perc_withAtrack4p,Perc_withAtrack3p_xy,Perc_withAtrack3p_xz, Perc_withX1_cl,Perc_withX2_cl, Perc_withX3_cl,Perc_withX4_cl,Perc_withY1_cl,Perc_withY2_cl, Perc_withY3_cl,Perc_withY4_cl, Perc_Ntracks_xy_3p,Perc_Ntracks_xz_3p,Perc_Ntracks_3p,Perc_Ntracks_xy_4p,Perc_Ntracks_xz_4p,Perc_Ntracks_4p;
  Perc_withAtrack3p=(double)Nevents_withAtrack3p;
  Perc_withAtrack4p=(double)Nevents_withAtrack4p;
  Perc_withAtrack3p_xy=(double)Nevents_withAtrack3p_xy;
  Perc_withAtrack3p_xz=(double)Nevents_withAtrack3p_xz;
  Perc_withX1_cl=(double)Nevents_withX1_cl;
  Perc_withX2_cl=(double)Nevents_withX2_cl;
  Perc_withX3_cl=(double)Nevents_withX3_cl;
  Perc_withX4_cl=(double)Nevents_withX4_cl;
  Perc_withY1_cl=(double)Nevents_withY1_cl;
  Perc_withY2_cl=(double)Nevents_withY2_cl;
  Perc_withY3_cl=(double)Nevents_withY3_cl;
  Perc_withY4_cl=(double)Nevents_withY4_cl;
  Perc_Ntracks_xy_3p =(double) N_tracks_xy_3p;
  Perc_Ntracks_xz_3p =(double) N_tracks_xz_3p;
  Perc_Ntracks_3p= (double) N_tracks_3p;
  Perc_Ntracks_xy_4p =(double) N_tracks_xy_4p;
  Perc_Ntracks_xz_4p =(double) N_tracks_xz_4p;
  Perc_Ntracks_4p= (double) N_tracks_4p;
  
  // MINI TREE WITH GENERAL RUN INFO //////////////////////////////////////////////////
  TFile *file_OnePHe = new TFile(MiniRunTreeName,"recreate");
  TTree *tree_OnePhe = new TTree("Run_info","A tree containing some run general info");
  tree_OnePhe->Branch("Nevents",&ev);
  tree_OnePhe->Branch("Run",&run);
  tree_OnePhe->Branch("boards_OnePhes",&boards_OnePhes);
  tree_OnePhe->Branch("WorkingPoint",&WorkingPoint);
  tree_OnePhe->Branch("Temperature",&Temperature);
  tree_OnePhe->Branch("TriggerRate",&TriggerRate);
  tree_OnePhe->Branch("Nev_withAtrack3p",&Perc_withAtrack3p);
  tree_OnePhe->Branch("Nev_withAtrack3p_xy",&Perc_withAtrack3p_xy);
  tree_OnePhe->Branch("Nev_withAtrack3p_xz",&Perc_withAtrack3p_xz);
  tree_OnePhe->Branch("Nev_withAtrack4p",&Perc_withAtrack4p);
  tree_OnePhe->Branch("Ntracks_xy_3p", &Perc_Ntracks_xy_3p);
  tree_OnePhe->Branch("Ntracks_xz_3p", &Perc_Ntracks_xz_3p);
  tree_OnePhe->Branch("Ntracks_3p", &Perc_Ntracks_3p);
  tree_OnePhe->Branch("Ntracks_xy_4p", &Perc_Ntracks_xy_4p);
  tree_OnePhe->Branch("Ntracks_xz_4p", &Perc_Ntracks_xz_4p);
  tree_OnePhe->Branch("Ntracks_4p", &Perc_Ntracks_4p);
  tree_OnePhe->Branch("Nev_withAcluster_Z1",&Perc_withX1_cl);
  tree_OnePhe->Branch("Nev_withAcluster_Z2",&Perc_withX2_cl);
  tree_OnePhe->Branch("Nev_withAcluster_Z3",&Perc_withX3_cl);
  tree_OnePhe->Branch("Nev_withAcluster_Z4",&Perc_withX4_cl);
  tree_OnePhe->Branch("Nev_withAcluster_Y1",&Perc_withY1_cl);
  tree_OnePhe->Branch("Nev_withAcluster_Y2",&Perc_withY2_cl);
  tree_OnePhe->Branch("Nev_withAcluster_Y3",&Perc_withY3_cl);
  tree_OnePhe->Branch("Nev_withAcluster_Y4",&Perc_withY4_cl);
  tree_OnePhe->Branch("EnergyCut_clusterStrip", &EnergyCut_clusterStrip);
  tree_OnePhe->Branch("EnergyCut_singleStrip", &EnergyCut_singleStrip);
  tree_OnePhe->Branch("EnergyCut_additionalStrip", &EnergyCut_additionalStrip);
  tree_OnePhe->Branch("proximity_cut_xz", &proximity_cut_xz);
  tree_OnePhe->Branch("proximity_cut_xy", &proximity_cut_xy);
  tree_OnePhe->Branch("RunDuration",&RunDuration);
  tree_OnePhe->Branch("Mean_clusterEnergy_Z1", &Mean_clusterEnergy_Z1);
  tree_OnePhe->Branch("Mean_clusterEnergy_Z2", &Mean_clusterEnergy_Z2);
  tree_OnePhe->Branch("Mean_clusterEnergy_Z3", &Mean_clusterEnergy_Z3);
  tree_OnePhe->Branch("Mean_clusterEnergy_Z4", &Mean_clusterEnergy_Z4);
  tree_OnePhe->Branch("RMS_clusterEnergy_Z1", &RMS_clusterEnergy_Z1);
  tree_OnePhe->Branch("RMS_clusterEnergy_Z2", &RMS_clusterEnergy_Z2);
  tree_OnePhe->Branch("RMS_clusterEnergy_Z3", &RMS_clusterEnergy_Z3);
  tree_OnePhe->Branch("RMS_clusterEnergy_Z4", &RMS_clusterEnergy_Z4);
  tree_OnePhe->Branch("Mean_clusterEnergy_Y1", &Mean_clusterEnergy_Y1);
  tree_OnePhe->Branch("Mean_clusterEnergy_Y2", &Mean_clusterEnergy_Y2);
  tree_OnePhe->Branch("Mean_clusterEnergy_Y3", &Mean_clusterEnergy_Y3);
  tree_OnePhe->Branch("Mean_clusterEnergy_Y4", &Mean_clusterEnergy_Y4);
  tree_OnePhe->Branch("RMS_clusterEnergy_Y1", &RMS_clusterEnergy_Y1);
  tree_OnePhe->Branch("RMS_clusterEnergy_Y2", &RMS_clusterEnergy_Y2);
  tree_OnePhe->Branch("RMS_clusterEnergy_Y3", &RMS_clusterEnergy_Y3);
  tree_OnePhe->Branch("RMS_clusterEnergy_Y4", &RMS_clusterEnergy_Y4);
  tree_OnePhe->Branch("Mean_Nclusters_Z1",&Mean_NclustersZ1);
  tree_OnePhe->Branch("RMS_Nclusters_Z1",&RMS_NclustersZ1);
  tree_OnePhe->Branch("Mean_Nclusters_Z2",&Mean_NclustersZ2);
  tree_OnePhe->Branch("RMS_Nclusters_Z2",&RMS_NclustersZ2);
  tree_OnePhe->Branch("Mean_Nclusters_Z3",&Mean_NclustersZ3);
  tree_OnePhe->Branch("RMS_Nclusters_Z3",&RMS_NclustersZ3);
  tree_OnePhe->Branch("Mean_Nclusters_Z4",&Mean_NclustersZ4);
  tree_OnePhe->Branch("RMS_Nclusters_Z4",&RMS_NclustersZ4);
  tree_OnePhe->Branch("Mean_Nclusters_Y1",&Mean_NclustersY1);
  tree_OnePhe->Branch("RMS_Nclusters_Y1",&RMS_NclustersY1);
  tree_OnePhe->Branch("Mean_Nclusters_Y2",&Mean_NclustersY2);
  tree_OnePhe->Branch("RMS_Nclusters_Y2",&RMS_NclustersY2);
  tree_OnePhe->Branch("Mean_Nclusters_Y3",&Mean_NclustersY3);
  tree_OnePhe->Branch("RMS_Nclusters_Y3",&RMS_NclustersY3);
  tree_OnePhe->Branch("Mean_Nclusters_Y4",&Mean_NclustersY4);
  tree_OnePhe->Branch("RMS_Nclusters_Y4",&RMS_NclustersY4);
  tree_OnePhe->Branch("datime",&datime);
  tree_OnePhe->Branch("IsOnePhe_aCopy",&boards_Is1phe_copy);
  tree_OnePhe->Branch("NGoodTracks3p",&NGooDTracks3p);
  tree_OnePhe->Branch("NGoodTracks4p",&NGooDTracks4p);
  tree_OnePhe->Fill();
  tree_OnePhe->Write();
  file_OnePHe->Close();
  

  ADCfile.close();
  }
  //FINAL TIME OF EXECUTION --> PRINT
  clock_t end =clock();
  double elapsed = (double) (end - start)/CLOCKS_PER_SEC;
  cout << "Number of events: " << ev <<endl;
  cout << "RUN " << run << " COMPLETED. Find your data in ---> "<< ROOTfileName << endl;
  cout << " Find your  minitree in ---> "<<MiniRunTreeName<< endl;
  //printf("\nEXECUTION TIME: %.1f seconds \n", (double)elapsed);


}
