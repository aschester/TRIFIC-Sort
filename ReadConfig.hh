#ifndef __READCONFIG_HH
#define __READCONFIG_HH

#include <iostream>
#include <string>
#include <vector>

#include "TEnv.h"
#include "TString.h"

class ReadConfig {
 public:

  ReadConfig(const char*, const Int_t);
  ~ReadConfig(){};

  // stuff for accessing class members
  TString GetAnalysisFileName() { return fAnalysisFileName; };
  TString GetFragmentFileName() { return fFragmentFileName; };
  TString GetCalibrationFileName() { return fCalibrationFileName; };
  TString GetOutputFileName() { return fOutputFileName; };
  Bool_t GetFragSortFlag() { return (fFragSort==1) ? true : false; };
  Int_t GetTimeRefCh() { return fTimeRefCh; };
  std::vector<Double_t> GetTimeGateLimits();
  std::vector<Double_t> GetBGLowTimeGateLimits();
  std::vector<Double_t> GetBGHighTimeGateLimits();
  Double_t GetDopplerBeta() { return fDopplerBeta; };
  TString GetPid1FileName() { return fPid1Name; }; // grid-by-grid energy gates
  TString GetPid2FileName() { return fPid2Name; }; // grid-by-grid energy gates
  std::vector<Int_t> GetDeltaE1Range();
  std::vector<Int_t> GetDeltaE2Range();


 private:

  Int_t fVerbosityLevel;

  // file I/O
  TString fAnalysisFileName;
  TString fFragmentFileName;
  TString fCalibrationFileName;
  TString fOutputFileName;

  // Fragment tree
  Int_t fFragSort; // TEnv reads boolean but returns int (sigh)  

  // Timing
  Int_t fTimeRefCh;
  Double_t fTimeLimLow;
  Double_t fTimeLimHigh;
  Double_t fBGLowLimLow;
  Double_t fBGLowLimHigh;
  Double_t fBGHighLimLow;
  Double_t fBGHighLimHigh;

  // Doppler correction
  Double_t fDopplerBeta;

  // PID gates
  TString fPid1Name;
  TString fPid2Name;

  // DeltaE grids
  Int_t fdE1low,fdE1high,fdE2low,fdE2high;

};

#endif
