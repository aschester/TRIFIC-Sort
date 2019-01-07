#include "ReadConfig.hh"

using namespace std;

ReadConfig::ReadConfig(const char* filename, const Int_t verbosityLevel)
  : fVerbosityLevel(verbosityLevel) {

  TEnv env;
  env.ReadFile(filename,kEnvLocal);
  // env.Print();

  // file I/O
  fAnalysisFileName = env.GetValue("ANALYSIS.FILE","analysis.root");
  fFragmentFileName = env.GetValue("FRAGMENT.FILE","fragment.root");
  fCalibrationFileName = env.GetValue("CALIBRATION.FILE","CalFiles/TRIFIC_Cal.cal");
  fOutputFileName = env.GetValue("OUTPUT.FILE","Histograms.root");

  // Do you want to sort the fragment tree?
  fFragSort = env.GetValue("SORT.FRAGMENT.TREE",0); // don't sort by default

  // Timing
  fTimeRefCh = env.GetValue("TRIFIC.TIME.REF.CHANNEL",1);
  fTimeLimLow = env.GetValue("TIME.COINC.LIMIT.LOW",-125.);
  fTimeLimHigh = env.GetValue("TIME.COINC.LIMIT.HIGH",125.);
  fBGLowLimLow = env.GetValue("TIME.BG.LOW.LIMIT.LOW",-1000.);
  fBGLowLimHigh = env.GetValue("TIME.BG.LOW.LIMIT.HIGH",-250.);
  fBGHighLimLow = env.GetValue("TIME.BG.HIGH.LIMIT.LOW",250.);
  fBGHighLimHigh = env.GetValue("TIME.BG.HIGH.LIMIT.HIGH",1000.);

  // Doppler shift
  fDopplerBeta = env.GetValue("DOPPLER.BETA",0.1);

  // Gates
  fPid1Name = env.GetValue("PID1.GATE.FILE","CalFiles/ArGates.dat");
  fPid2Name = env.GetValue("PID2.GATE.FILE","CalFiles/KGates.dat");

  // Delta E grid ranges
  fdE1low = env.GetValue("DELTA.E1.GRID.LOW",1);
  fdE1high = env.GetValue("DELTA.E1.GRID.HIGH",2);
  fdE2low = env.GetValue("DELTA.E2.GRID.LOW",5);
  fdE2high = env.GetValue("DELTA.E2.GRID.HIGH",6);


  if(fVerbosityLevel >= 1) {
    cout << "----> Config file settings <---------------------------" << endl;
    cout << endl;

    cout << "----- File I/O ----------------------------------------" << endl;
    cout << "Analysis tree read from:\t" << fAnalysisFileName << endl;
    cout << "Fragment tree read from:\t" << fFragmentFileName << endl;
    cout << "Calibraion read from:\t\t" << fCalibrationFileName << endl;
    cout << "Histograms saved to:\t\t" << fOutputFileName << endl;

    cout << "----- Fragment tree sort -------------------------------" << endl;
    if(fFragSort== 1) cout << "Sorting from fragment tree: YES" << endl;
    if(fFragSort == 0) cout << "Sorting from fragment tree: NO" << endl;
    else cout << "Fragment tree sorting condition is undefined!" << endl;


    cout << "----- Timing -------------------------------------------" << endl;
    cout << "TRIFIC channel\t" << fTimeRefCh << "\tused as time reference." << endl;
    cout << "TIGRESS-TRIFIC time coincidence window is\t" << fTimeLimLow << "\tto\t" << fTimeLimHigh << "\t[ns]." << endl;
    cout << "TIGRESS-TRIFIC low background window is \t" << fBGLowLimLow << "\tto\t" << fBGLowLimHigh << "\t[ns]." << endl;
    cout << "TIGRESS-TRIFIC high background window is\t" << fBGHighLimLow << "\tto\t" << fBGHighLimHigh << "\t[ns]." << endl;

    cout << "----- Doppler shift ------------------------------------" << endl;
    cout << "Beta for Doppler correction is\t" << fDopplerBeta << "\t[v/c]." << endl;

    cout << "----- PID gates  ---------------------------------------" << endl;
    cout << "PID 1 gate file name is\t" << fPid1Name << endl;
    cout << "PID 2 gate file name is\t" << fPid2Name << endl;

    cout << "----- Delta E Grid ranges  -----------------------------" << endl;
    cout << "dE1 range is\t" << fdE1low << "\tto\t" << fdE1high << "." << endl;
    cout << "dE2 range is\t" << fdE2low << "\tto\t" << fdE2high << "." << endl;

    cout << endl;
    cout << "----> End reading config file <-------------------------" << endl;

  }

 if(fVerbosityLevel == 2)
    env.Print();

}

vector<Double_t> ReadConfig::GetTimeGateLimits() {
  vector<Double_t> lim;
  lim.push_back(fTimeLimLow);
  lim.push_back(fTimeLimHigh);
  return lim;
}

vector<Double_t> ReadConfig::GetBGLowTimeGateLimits() {
  vector<Double_t> lim;
  lim.push_back(fBGLowLimLow);
  lim.push_back(fBGLowLimHigh);
  return lim;
}

vector<Double_t> ReadConfig::GetBGHighTimeGateLimits() {
  vector<Double_t> lim;
  lim.push_back(fBGHighLimLow);
  lim.push_back(fBGHighLimHigh);
  return lim;
}

vector<Int_t> ReadConfig::GetDeltaE1Range() {
  vector<Int_t> lim;
  lim.push_back(fdE1low);
  lim.push_back(fdE1high);
  return lim;
}

vector<Int_t> ReadConfig::GetDeltaE2Range() {
  vector<Int_t> lim;
  lim.push_back(fdE2low);
  lim.push_back(fdE2high);
  return lim;
}
