#ifndef __SORTDATA_HH
#define __SORTDATA_HH

#include <stdlib.h>
#include <fstream>
#include <iomanip>
#include <vector>

#include "TH2.h"
#include "TH1.h"
#include "TFile.h"
#include "TList.h"
#include "TTree.h"
#include "TChain.h"
#include "TTigress.h"
#include "TGenericDetector.h"
#include "TRF.h"
#include "TFragment.h"
#include "TCutG.h"
#include "TString.h"

#include "ReadConfig.hh"

// *****************************************************************
//                                                                 *
// Honestly this whole thing should be a class instead of all this *
// global scope nonsense. But that's for another day...            *
// -ASC, 16 Nov 2018                                               *
//                                                                 *
// *****************************************************************

const Int_t S1K = 1024;
const Int_t S2K = 2048;
const Int_t S4K = 4096;
const Int_t S8K = 8192;
const Int_t S16K = 16384;
const Int_t S32K = 32768;
const Int_t S65K = 65536;
const Int_t NTIG = 64;
const Int_t NTRF = 32;
const Int_t NRING = 7;
const Int_t MAXCHAN = 48;

std::vector<std::string> store;

// config file data and associated/derived quantities
ReadConfig* readConfig;
TString inafile, inffile, calfile, outfile;
Bool_t fragSort;
Int_t refCh; // reference channel for time difference
std::vector<Double_t> timeGate, bgLow, bgHigh;
std::vector<Int_t> dE1, dE2;
Double_t beta;
Double_t bgLowScale, bgHighScale;
TString ArGateFileName, KGateFileName;

int ReadGateFile(const char*, std::vector< std::vector<Double_t> >&);
// end config 

TList *trfRawEnergyList, 
  *trfRawTimeList, 
  *chRawEnergyTimeList, 
  *tigRawEnergyList, 
  *tigRawEnergySumList, 
  *tigRawTimeList, 
  *trfGatedEnergyList, 
  *trfGatedTimeList, 
  *tigGatedEnergyList, 
  *tigGatedEnergyABSumList, 
  *tigGatedEnergyABRingList, 
  *tigGatedEnergyABDopSumList, 
  *tigGatedEnergyABRingDopSumList, 
  *tigGatedEnergyBkgSumList, 
  *tigGatedTimeList, 
  *rawList, 
  *gatedList;

// raw histograms
// TRIFIC
TH1I* trific_multiplicity;
TH1I* trific_idmin;
TH1I* trific_idmax;
TH2I* trific_tigress_multiplicity;
TH2D* trific_bragg;
TH2D* trific_pid;
TH2I* trific_kvalue;
TH1D* trific_energy[NTRF];
TH1D* trific_time[NTRF];
TH2D* channel_vs_time[MAXCHAN];
TH1D* trific_grrgrX_tdiff[NTRF];
// TIGRESS
TH1I* tigress_multiplicity;
TH2I* tigress_kvalue;
// Energy
TH1D* tigress_energy_sum;
// Addback energy
TH1D* tigress_energyAB[NTIG];
TH1D* tigress_energyAB_sum;
// Doppler-corrected addback energy
TH1D* tigress_energyABDop_sum;
// Time
TH1D* tigress_time[NTIG];
// TIGRESS-TRIFIC time difference
TH1D* tigress_trificgr_tdiff;
TH1D* tigress_trificgr_PIDAr_tdiff;
TH1D* tigress_trificgr_PIDK_tdiff;

// gated histograms
// TRIFIC
TH1I* gtrific_multiplicity;
TH1I* gtrific_idmin;
TH1I* gtrific_idmax;
TH2D* gtrific_bragg;
TH2D* gtrific_bragg_PIDAr;
TH2D* gtrific_bragg_PIDK;
TH2D* gtrific_pid;
TH2D* gtrific_pid_tcut;
TH1D* gtrific_energy[NTRF];
TH1D* gtrific_time[NTRF];
// TIGRESS
TH1I* gtigress_multiplicity;
// Addback energy
TH1D* gtigress_energyAB[NTIG];
TH1D* gtigress_energyAB_sum;
TH1D* gtigress_energyAB_bgLow_sum;
TH1D* gtigress_energyAB_bgHigh_sum;
TH1D* gtigress_energyABSub_sum;
TH1D* gtigress_energyABRing[NRING];
// Doppler-corrected addback energy
TH1D* gtigress_energyABDop_sum;
TH1D* gtigress_energyABDop_bgLow_sum;
TH1D* gtigress_energyABDop_bgHigh_sum;
TH1D* gtigress_energyABDopSub_sum;
// Doppler-corrected addback PID-gated energy
TH1D* gtigress_energyABDop_PIDAr_sum;
TH1D* gtigress_energyABDop_bgLow_PIDAr_sum;
TH1D* gtigress_energyABDop_bgHigh_PIDAr_sum;
TH1D* gtigress_energyABDopSub_PIDAr_sum;
TH1D* gtigress_energyABRingDop_PIDAr_sum[NRING];

TH1D* gtigress_energyABDop_PIDArCoulex_sum;
TH1D* gtigress_energyABDop_bgLow_PIDArCoulex_sum;
TH1D* gtigress_energyABDop_bgHigh_PIDArCoulex_sum;
TH1D* gtigress_energyABDopSub_PIDArCoulex_sum;
TH1D* gtigress_energyABRingDop_PIDArCoulex_sum[NRING];

TH1D* gtigress_energyABDop_PIDK_sum;
TH1D* gtigress_energyABDop_bgLow_PIDK_sum;
TH1D* gtigress_energyABDop_bgHigh_PIDK_sum;
TH1D* gtigress_energyABDopSub_PIDK_sum;
TH1D* gtigress_energyABRingDop_PIDK_sum[NRING];

// Grid-by-grid energy gates
TH1D* gtigress_energyABDop_EPIDAr_sum;
TH1D* gtigress_energyABDop_bgLow_EPIDAr_sum;
TH1D* gtigress_energyABDop_bgHigh_EPIDAr_sum;
TH1D* gtigress_energyABDopSub_EPIDAr_sum;

TH1D* gtigress_energyABDop_EPIDK_sum;
TH1D* gtigress_energyABDop_bgLow_EPIDK_sum;
TH1D* gtigress_energyABDop_bgHigh_EPIDK_sum;
TH1D* gtigress_energyABDopSub_EPIDK_sum;

// Time
TH1D* gtigress_time[NTIG];
// TIGRESS-TRIFIC time difference
TH1D* gtigress_trificgr_tdiff;
TH1D* gtigress_trificgr_PIDAr_tdiff;
TH1D* gtigress_trificgr_PIDK_tdiff;
// TIGRESS-TRIFIC 2D coincidence plots
TH2D* gtrific_mult_vs_tigress_energyABDop_sum;

// some other TRIFIC-TIGRESS energy plots
TH1D* gtrificEsum;
TH2D* gtrificEsum_vs_tigressE;
TH2D* gtrificEsum_vs_tigressE_PIDAr;
TH2D* gtrificEsum_vs_tigressE_PIDK;


// Grid-by-grid energy gates
std::vector< std::vector<Double_t> > ArGates, KGates;


void InitializeHistograms();
Int_t SortData();
void Print();

#endif
