#include "SortData.hh"

using namespace std;

Int_t SortData() {
  TFile* inputafile = new TFile(inafile.Data(),"READ");
  if(!inputafile->IsOpen()) {
    cout << "Opening analysis file failed, aborting" << endl;
    return 1;
  } else {
    cout << "File " << inafile << " opened successfully..." << endl;
  }
  TChain* AnalysisTree = (TChain*)inputafile->Get("AnalysisTree");

  TFile* inputffile = new TFile(inffile.Data(),"READ");
  if(!inputffile->IsOpen()){
    cout << "Opening fragment file failed, aborting" << endl;
    return 1;
  } else {
    cout << "File " << inffile << " opened successfully..." << endl;
  }
  TChain* FragmentTree = (TChain*)inputffile->Get("FragmentTree");

  cout << "Reading calibration file: " << calfile << endl;;
  TChannel::ReadCalFile(calfile.Data());

  TTigress* tigress = nullptr;
  TGenericDetector* trific = nullptr;
  if(AnalysisTree->FindBranch("TTigress")) {
    AnalysisTree->SetBranchAddress("TTigress",&tigress);
  } else {
    cout << "Branch 'TTigress' not found! TTigress variable is NULL pointer!" << endl;
  }
  if(AnalysisTree->FindBranch("TGenericDetector")) {
    AnalysisTree->SetBranchAddress("TGenericDetector",&trific);
  } else {
    cout << "Branch 'TGenericDetector' not found! TGenericDetector variable is NULL pointer" << endl;
  }

  TFragment* frag = nullptr;
  if(FragmentTree->FindBranch("TFragment")) {
    FragmentTree->SetBranchAddress("TFragment",&frag);
  } else {
    cout << "Branch 'TFragment' not found! TFragment vairable is NULL pointer!" << endl;
  }
  
  Long_t naentries = AnalysisTree->GetEntries();

  Long_t nfentries = FragmentTree->GetEntries();

  cout << "Number of entries in fragment tree: " << nfentries << endl;
  cout << "Number of entries in analysis tree: " << naentries << endl;
  cout << "Fragments per event: " << 
    static_cast<Double_t>(nfentries)/static_cast<Double_t>(naentries) << endl;

  // useful stuff for later
  // time gate window sizes
  // TCuts for PID gating
   TCutG *CutAr = new TCutG("CutAr",18);
   CutAr->SetPoint(0,1308.66,1583.7);
   CutAr->SetPoint(1,1365.47,1627.81);
   CutAr->SetPoint(2,1412.94,1630.75);
   CutAr->SetPoint(3,1443.14,1592.52);
   CutAr->SetPoint(4,1482.7,1486.65);
   CutAr->SetPoint(5,1515.78,1239.62);
   CutAr->SetPoint(6,1533.04,1069.05);
   CutAr->SetPoint(7,1527.28,992.592);
   CutAr->SetPoint(8,1507.15,989.651);
   CutAr->SetPoint(9,1473.35,1177.86);
   CutAr->SetPoint(10,1430.2,1354.31);
   CutAr->SetPoint(11,1395.68,1419.01);
   CutAr->SetPoint(12,1311.54,1377.84);
   CutAr->SetPoint(13,1273.42,1374.9);
   CutAr->SetPoint(14,1249.69,1413.13);
   CutAr->SetPoint(15,1250.41,1471.95);
   CutAr->SetPoint(16,1274.86,1545.47);
   CutAr->SetPoint(17,1308.66,1583.7);
   TCutG *CutK = new TCutG("CutK",13);
   CutK->SetPoint(0,1446.02,1636.63);
   CutK->SetPoint(1,1476.94,1539.58);
   CutK->SetPoint(2,1515.78,1392.54);
   CutK->SetPoint(3,1548.86,1174.92);
   CutK->SetPoint(4,1594.16,801.438);
   CutK->SetPoint(5,1598.48,319.144);
   CutK->SetPoint(6,1668.24,313.262);
   CutK->SetPoint(7,1697.72,442.658);
   CutK->SetPoint(8,1684.78,924.953);
   CutK->SetPoint(9,1622.21,1433.71);
   CutK->SetPoint(10,1539.51,1692.51);
   CutK->SetPoint(11,1452.49,1686.63);
   CutK->SetPoint(12,1446.02,1636.63);
   
  cout << "Starting sort..." << endl;

  // // Sorting from fragment trees
  // // TRIFIC
  if(fragSort) {
    FragmentTree->GetEntry(0); 
    Long_t tstart = frag->GetMidasTimeStamp();
    for(Long_t entry=0; entry<nfentries; entry++) {
      //for(Long_t entry=0; entry<10000; entry++) {
      FragmentTree->GetEntry(entry);
      Int_t Id = frag->GetChannelNumber();    
      Double_t E = frag->GetCharge();
      Long_t t = frag->GetMidasTimeStamp()-tstart; // in s
      channel_vs_time[Id]->Fill(t,E);
      if(entry%10000 == 0) {
	cout << setiosflags(ios::fixed) << "Fragment tree: entry " << entry << " of " 
	     << nfentries << ", " << 100 * entry/nfentries << "% complete" << "\r" << flush;
      }
    }
    cout << "Fragement tree: entry " << nfentries << " of " << nfentries << ", 100% Complete!\n";
  }
  
  // Sorting from analysis trees
  Int_t nAr = 0, nK = 0, nArK = 0, nNone = 0, nUnk = 0; // number of time-gated events in Ar only, K only, Ar and K, ~Ar and ~K, unknown event type
  for(Long_t entry=0; entry<naentries; entry++) {

    // cout << "---------- New event ---------- " << endl;

    AnalysisTree->GetEntry(entry);

    // Flags and global event data
    Int_t IdMin = 100, IdMax = -1, trfMult = -1;
    Double_t tref = -S1K; // negative if not set
    Double_t deltaE1 = 0.0, deltaE2 = 0.0;
    Double_t Esum = 0.0;
    Bool_t goodMultFlag = false, refGridFlag = false; // multiplicity and reference grid flags
    Bool_t timeGatedFlag = false; // time coincidence flag
    Bool_t bgLowFlag = false, bgHighFlag = false; // background flags
    Bool_t ArEnergyGateFlag = true, KEnergyGateFlag = true; // Good energy gating condition

    if(trific && tigress)
      trific_tigress_multiplicity->Fill(trific->GetMultiplicity(),tigress->GetMultiplicity());

    // TRIFIC
    if(trific) {
      trfMult = trific->GetMultiplicity(); // built event multiplicity
      trific_multiplicity->Fill(trfMult);      

      for(Int_t i=0; i<trfMult; i++) {
	TGRSIDetectorHit* trfh = trific->GetHit(i);
	Int_t Id = trfh->GetDetector();

	// find first and last grid hit
	if(Id < IdMin) IdMin = Id;
	if(Id > IdMax) IdMax = Id;

	Double_t E = trfh->GetEnergy(); // charge for uncalibrated

	// reasonable energy
	if(E > 5.) { 
	  if(Id == refCh) { 
	    refGridFlag = true; // do we have reference grid in the event
	    tref = trfh->GetTime(); // in ns
	  }
	}
      }

      // "Good multiplicity" every grid from the 
      // first to the last in the event registers a hit
      if(trfMult == IdMax-IdMin+1)
	goodMultFlag = true;

      trific_idmin->Fill(IdMin);
      trific_idmax->Fill(IdMax);

      // only want events where grids are not dropped
      if(goodMultFlag) { 
	for(Int_t i=0; i<trfMult; i++) {
	  TGRSIDetectorHit* trfh = trific->GetHit(i);
	  Int_t Id = trfh->GetDetector();
	  Double_t E = trfh->GetEnergy(); // charge for uncalibrated
	  Int_t k = static_cast<Int_t>(trfh->GetKValue());
	  Double_t t = trfh->GetTime(); // in ns
	  
	  // energy is something reasonable
	  if(E > 5.) { 

	    Esum += E;

	    // Argon grid energy gates
	    if( (ArGates[Id][0] > 0) && (ArGates[Id][1] > 0) ) {
	      if( (E > ArGates[Id][0]) && (E < ArGates[Id][1]) )
		ArEnergyGateFlag = ArEnergyGateFlag && true;
	      else
		ArEnergyGateFlag = ArEnergyGateFlag && false;
	    }    
	    // Potassium grid energy gates
	    if( (KGates[Id][0] > 0) && (KGates[Id][1] > 0) ) {
	      if( (E > KGates[Id][0]) && (E < KGates[Id][1]) )
		KEnergyGateFlag = KEnergyGateFlag && true;
	      else
		KEnergyGateFlag = KEnergyGateFlag && false;
	    }    
	    // cout << trfId << " " << trfE << " " << ArGates[trfId][0] << " " << ArGates[trfId][1] << " " << ArEnergyGateFlag << " " 
	    // 	 << KGates[trfId][0] << " " << KGates[trfId][1] << " " << KEnergyGateFlag << endl;
	  	    
	    // Get energies for deltaE1-deltaE2 PID
	    // Reaction data: 1-2, 5-6
	    // 5.1 MeV rate test: 1-3, 9-11
	    // 6.75 MeV rate test: 1-5, 14-18
	    // Calibration data: 1-4, 17-20
	    if(Id >= dE1[0] && Id <= dE1[1]) deltaE1 += E; 
	    else if(Id >= dE2[0] && Id <= dE2[1]) deltaE2 += E;    
	    
	    trific_bragg->Fill(static_cast<Double_t>(Id),E);
	    trific_energy[Id]->Fill(E);
	    trific_time[Id]->Fill(t/1.0E9); // contract time to s
	    trific_kvalue->Fill(Id,k);

	    if(refGridFlag) {
	      trific_grrgrX_tdiff[Id]->Fill(tref-t); // absolute time diff is ns
	    }
	  }
	}
	// range condition
	if(IdMax > 8)
	  ArEnergyGateFlag = ArEnergyGateFlag && true;
	else 
	  ArEnergyGateFlag = ArEnergyGateFlag && false;
	if( (IdMax > 5) && (IdMax < 10) )
	  KEnergyGateFlag = KEnergyGateFlag && true;
	else 
	  KEnergyGateFlag = KEnergyGateFlag && false;

        // getc(stdin);	
	trific_pid->Fill(deltaE1,deltaE2); // fill PID after looping over multiplicity
      }
    } 
    
    // TIGRESS and TRIFIC coincidences
    // Require good TRIFIC multiplicity and presence of reference grid
    if(goodMultFlag && refGridFlag && tigress) {
      Int_t tigMult = tigress->GetMultiplicity();
      // TIGRESS no addback
      for(Int_t i=0; i<tigMult; i++) {
	TTigressHit* tigh = tigress->GetTigressHit(i);
	Int_t tigId = tigh->GetArrayNumber();
	Double_t tigE = tigh->GetEnergy(); // in keV	  
	Int_t tigk = static_cast<Int_t>(tigh->GetKValue());
	if(tigE > 10.) {
	  tigress_kvalue->Fill(tigId,tigk);
	  tigress_energy_sum->Fill(tigE);
	}
      }
      // End TIGRESS no addback
      
      // TIGRESS addback
      // Event multiplicity, time, time gates, etc performed with addback event info
      tigMult = tigress->GetAddbackMultiplicity(); // reset mult to addback mult
      tigress_multiplicity->Fill(tigMult); // addback multiplicity
      if(tigMult > 0) gtrific_pid->Fill(deltaE1,deltaE2); // in coinc. with TIGRESS
      for(Int_t i=0; i<tigMult; i++) {
	TTigressHit* tigh = tigress->GetAddbackHit(i);
	Int_t tigId = tigh->GetArrayNumber();
	Int_t tigRing = tigh->GetRing();
	Double_t tigE = tigh->GetEnergy(); // in keV	  
	Double_t tigE_dop = tigh->GetDoppler(beta); // in keV
	Double_t tigt = tigh->GetTime(); // in ns
	Double_t tdiff = S65K; // some big number	   
	
	// Reasonable energy check
	if(tigE > 10.) { 
	  tigress_energyAB[tigId]->Fill(tigE);
	  tigress_time[tigId]->Fill(tigt/1.0E9); // contract time to s
	  tigress_energyAB_sum->Fill(tigE);
	  tigress_energyABDop_sum->Fill(tigE_dop);		  
	}
	
	// get coincidence with TRIFIC
	if(trific) {
	  tdiff = tigt-tref; // TIGRESS-reference grid time difference in ns
	  tigress_trificgr_tdiff->Fill(tdiff);
	   if(CutAr->IsInside(deltaE1,deltaE2)) 
	  // if(ArEnergyGateFlag && !KEnergyGateFlag)
	  // if(ArEnergyGateFlag)
	    tigress_trificgr_PIDAr_tdiff->Fill(tdiff);
	   if(CutK->IsInside(deltaE1,deltaE2)) 
	  // if(KEnergyGateFlag && !ArEnergyGateFlag)
	  // if(KEnergyGateFlag)
	    tigress_trificgr_PIDK_tdiff->Fill(tdiff);	  
	  
	  // check if time is in 1D time gate and flag background events
	  if(tdiff > timeGate[0] && tdiff < timeGate[1]) { 
	    timeGatedFlag = true;
	    gtigress_trificgr_tdiff->Fill(tdiff);
	    if(CutAr->IsInside(deltaE1,deltaE2)) 
	    // if(ArEnergyGateFlag && !KEnergyGateFlag)
	    // if(ArEnergyGateFlag)
	      gtigress_trificgr_PIDAr_tdiff->Fill(tdiff);
	    if(CutK->IsInside(deltaE1,deltaE2)) 
	    // if(KEnergyGateFlag && !ArEnergyGateFlag)
	    // if(KEnergyGateFlag)
	      gtigress_trificgr_PIDK_tdiff->Fill(tdiff);
	  }
	  else if(tdiff > bgLow[0]  && tdiff < bgLow[1])  bgLowFlag = true;
	  else if(tdiff > bgHigh[0] && tdiff < bgHigh[1]) bgHighFlag = true; 
	}
	// End TRIFIC coincidence analysis

	// if good time write out TIGRESS in coinc. with TRIFIC reference grid
	if(timeGatedFlag) { 
	  if(tigE > 10.) {
	    gtigress_multiplicity->Fill(tigMult);
	    gtigress_energyAB[tigId]->Fill(tigE);
	    gtigress_energyAB_sum->Fill(tigE);
	    gtigress_energyABRing[tigRing]->Fill(tigE);
	    gtigress_time[tigId]->Fill(tigt/1.0E9); // contracted time	    
	    gtigress_energyABDop_sum->Fill(tigE_dop);
	    if(CutAr->IsInside(deltaE1,deltaE2)) {
	      gtigress_energyABDop_PIDAr_sum->Fill(tigE_dop);
	      // if(trfMult >= 11) gtigress_energyABDop_PIDAr_sum->Fill(tigE_dop);
	      // else if(trfMult >= 8 && trfMult <=10) gtigress_energyABDop_PIDArCoulex_sum->Fill(tigE_dop);
	      gtigress_energyABRingDop_PIDAr_sum[tigRing]->Fill(tigE_dop);
	      gtrificEsum_vs_tigressE_PIDAr->Fill(Esum,tigE_dop);
	      
	    }
	    if(CutK->IsInside(deltaE1,deltaE2)) {
	      gtigress_energyABDop_PIDK_sum->Fill(tigE_dop);
	      gtigress_energyABRingDop_PIDK_sum[tigRing]->Fill(tigE_dop);
	      gtrificEsum_vs_tigressE_PIDK->Fill(Esum,tigE_dop);
	    }

	    // grid energy gates
	    // if(ArEnergyGateFlag && !KEnergyGateFlag)
	    if(ArEnergyGateFlag)
	      gtigress_energyABDop_EPIDAr_sum->Fill(tigE_dop);
	    // if(KEnergyGateFlag && !ArEnergyGateFlag)
	    if(KEnergyGateFlag)
	      gtigress_energyABDop_EPIDK_sum->Fill(tigE_dop);

	    // time-gated with good energy
	    if(ArEnergyGateFlag && !KEnergyGateFlag)
	      nAr++;
	    if(KEnergyGateFlag && !ArEnergyGateFlag)
	      nK++;
	    if(ArEnergyGateFlag && KEnergyGateFlag)
	      nArK++;
	    if(!ArEnergyGateFlag && !KEnergyGateFlag)
	      nNone++;
	    else
	      nUnk++;
	  }	  

	  // TIGRESS-TRIFIC 2D coincidence spectra
	  gtrific_mult_vs_tigress_energyABDop_sum->Fill(trfMult,tigE_dop);
	  gtrificEsum_vs_tigressE->Fill(Esum,tigE_dop);
	  // End TIGRESS-TRIFIC 2D coincidence spectra
	}
	
	// fill background spectra
	if(bgLowFlag) {
	  if(tigE > 10.) {
	    gtigress_energyAB_bgLow_sum->Fill(tigE); 
	    gtigress_energyABDop_bgLow_sum->Fill(tigE_dop); 
	    if(CutAr->IsInside(deltaE1,deltaE2)) {
	      gtigress_energyABDop_bgLow_PIDAr_sum->Fill(tigE_dop);
	      // if(trfMult >= 11) gtigress_energyABDop_bgLow_PIDAr_sum->Fill(tigE_dop);
	      // else if(trfMult >= 8 && trfMult <= 10) gtigress_energyABDop_bgLow_PIDArCoulex_sum->Fill(tigE_dop);
	    }
	    if(CutK->IsInside(deltaE1,deltaE2))
	      gtigress_energyABDop_bgLow_PIDK_sum->Fill(tigE_dop);

	    // grid energy gates
	    // if(ArEnergyGateFlag && !KEnergyGateFlag)
	    if(ArEnergyGateFlag)
	      gtigress_energyABDop_bgLow_EPIDAr_sum->Fill(tigE_dop);
	    // if(KEnergyGateFlag && !ArEnergyGateFlag)
	    if(KEnergyGateFlag)
	      gtigress_energyABDop_bgLow_EPIDK_sum->Fill(tigE_dop);

	  }
	}
	if(bgHighFlag) {
	  if(tigE > 10.) {
	    gtigress_energyAB_bgHigh_sum->Fill(tigE); 
	    gtigress_energyABDop_bgHigh_sum->Fill(tigE_dop); 
	    if(CutAr->IsInside(deltaE1,deltaE2)) {
	      gtigress_energyABDop_bgHigh_PIDAr_sum->Fill(tigE_dop);
	      // if(trfMult >= 11) gtigress_energyABDop_bgHigh_PIDAr_sum->Fill(tigE_dop);
	      // else if(trfMult >= 8 && trfMult <= 10) gtigress_energyABDop_bgHigh_PIDArCoulex_sum->Fill(tigE_dop);
	    }
	    if(CutK->IsInside(deltaE1,deltaE2))
	      gtigress_energyABDop_bgHigh_PIDK_sum->Fill(tigE_dop);
	  }

	    // grid energy gates
	    // if(ArEnergyGateFlag && !KEnergyGateFlag)
	    if(ArEnergyGateFlag)
	      gtigress_energyABDop_bgHigh_EPIDAr_sum->Fill(tigE_dop);
	    // if(KEnergyGateFlag && !ArEnergyGateFlag)
	    if(KEnergyGateFlag)
	      gtigress_energyABDop_bgHigh_EPIDK_sum->Fill(tigE_dop);
	}
      }
      // End TIGRESS addback

      // now do time coincidence for TRIFIC data
      if(timeGatedFlag) { 
	// in principle all this global stuff can go earlier
	// but lets keep it together in one place
	gtrific_multiplicity->Fill(trfMult);
	gtrific_idmin->Fill(IdMin);
	gtrific_idmax->Fill(IdMax);
	gtrific_pid_tcut->Fill(deltaE1,deltaE2); // in good time coinc. with TIGRESS
	gtrificEsum->Fill(Esum);
	// have to loop through TRIFIC events again for individual grid data :(
	for(Int_t i=0; i<trfMult; i++) {
	  TGRSIDetectorHit* trfh = trific->GetHit(i);
	  Double_t trfE = trfh->GetEnergy();
	  Int_t trfId = trfh->GetDetector();
	  Double_t trft = trfh->GetTime(); // in ns
	  // check for reasonable energy
	  if(trfE > 5.) { 

	    gtrific_bragg->Fill(static_cast<Double_t>(trfId),trfE);

	    // PID-gated Bragg curves
	    if(CutAr->IsInside(deltaE1,deltaE2))
	      gtrific_bragg_PIDAr->Fill(static_cast<Double_t>(trfId),trfE);
	    if(CutK->IsInside(deltaE1,deltaE2))
	      gtrific_bragg_PIDK->Fill(static_cast<Double_t>(trfId),trfE);
	    
	    gtrific_energy[trfId]->Fill(trfE);
	    gtrific_time[trfId]->Fill(trft/1.0E9); // contracted time
	  }
	}
	// getc(stdin);
      }
      // End time-gated TRIFIC
    }
    // End TIGRESS-TRIFIC coincidences
    if(entry%10000 == 0) {
      cout << setiosflags(ios::fixed) << "Analysis tree: entry " << entry << " of " 
	   << naentries << ", " << 100 * entry/naentries << "% complete" << "\r" << flush;  
    }
  }
  cout << "Analysis Tree: entry " << naentries << " of " << naentries << ", 100% Complete!\n";

  cout << "PID gate 1 integral: " << CutAr->IntegralHist(gtrific_pid_tcut) << "\tPID gate 2 integral: " << CutK->IntegralHist(gtrific_pid_tcut) 
       << " TOTAL: " << gtrific_pid_tcut->Integral() << endl;
  cout << "Ar gridE gate: " << nAr << " K gridE gate: " << nK << " Ar&&K gridE gate: " << nArK << " ~Ar&&~K gridE gate: " << nNone 
       << " unknown: " << nUnk << " TOTAL: " << nAr+nK+nArK+nNone+nUnk << endl;
  
  // do background subtraction
  gtigress_energyABSub_sum->Add(gtigress_energyAB_sum, 1.0);
  gtigress_energyABSub_sum->Add(gtigress_energyAB_bgLow_sum, bgLowScale);
  gtigress_energyABSub_sum->Add(gtigress_energyAB_bgHigh_sum, bgHighScale);

  gtigress_energyABDopSub_sum->Add(gtigress_energyABDop_sum, 1.0);
  gtigress_energyABDopSub_sum->Add(gtigress_energyABDop_bgLow_sum, bgLowScale);
  gtigress_energyABDopSub_sum->Add(gtigress_energyABDop_bgHigh_sum, bgHighScale);

  gtigress_energyABDopSub_PIDAr_sum->Add(gtigress_energyABDop_PIDAr_sum, 1.0);
  gtigress_energyABDopSub_PIDAr_sum->Add(gtigress_energyABDop_bgLow_PIDAr_sum, bgLowScale);
  gtigress_energyABDopSub_PIDAr_sum->Add(gtigress_energyABDop_bgHigh_PIDAr_sum, bgHighScale);

  gtigress_energyABDopSub_PIDK_sum->Add(gtigress_energyABDop_PIDK_sum, 1.0);
  gtigress_energyABDopSub_PIDK_sum->Add(gtigress_energyABDop_bgLow_PIDK_sum, bgLowScale);
  gtigress_energyABDopSub_PIDK_sum->Add(gtigress_energyABDop_bgHigh_PIDK_sum, bgHighScale);

  gtigress_energyABDopSub_EPIDAr_sum->Add(gtigress_energyABDop_EPIDAr_sum, 1.0);
  gtigress_energyABDopSub_EPIDAr_sum->Add(gtigress_energyABDop_bgLow_EPIDAr_sum, bgLowScale);
  gtigress_energyABDopSub_EPIDAr_sum->Add(gtigress_energyABDop_bgHigh_EPIDAr_sum, bgHighScale);

  gtigress_energyABDopSub_EPIDK_sum->Add(gtigress_energyABDop_EPIDK_sum, 1.0);
  gtigress_energyABDopSub_EPIDK_sum->Add(gtigress_energyABDop_bgLow_EPIDK_sum, bgLowScale);
  gtigress_energyABDopSub_EPIDK_sum->Add(gtigress_energyABDop_bgHigh_EPIDK_sum, bgHighScale);
  // End background subtraction

  cout << "Writing histograms to " << outfile << "..." << endl;
  TFile *myfile = new TFile(outfile.Data(),"RECREATE");
  myfile->cd();
  TDirectory *rawhist = myfile->mkdir("RawHistograms");
  rawhist->cd();
  rawList->Write();
  TDirectory *chEnergyTimeDir = rawhist->mkdir("ChannelEnergyVsTime");
  chEnergyTimeDir->cd();
  chRawEnergyTimeList->Write();
  rawhist->cd();
  TDirectory *trfEnergyDir =rawhist->mkdir("TrificEnergy");
  trfEnergyDir->cd();
  trfRawEnergyList->Write();
  rawhist->cd();
  TDirectory *trfTimeDir = rawhist->mkdir("TrificTime");
  trfTimeDir->cd();
  trfRawTimeList->Write();
  rawhist->cd();
  TDirectory *tigEnergyDir = rawhist->mkdir("TigressEnergy");
  tigEnergyDir->cd();
  tigRawEnergyList->Write();
  rawhist->cd();
  TDirectory *tigEnergySumDir = rawhist->mkdir("TigressEnergySum");
  tigEnergySumDir->cd();
  tigRawEnergySumList->Write();
  rawhist->cd();
  TDirectory *tigTimeDir = rawhist->mkdir("TigressTime");
  tigTimeDir->cd();
  tigRawTimeList->Write();
  rawhist->cd();
  myfile->cd();
  TDirectory *gatedhist = myfile->mkdir("GatedHistograms");
  gatedhist->cd();
  gatedList->Write();
  TDirectory *trfGatedEnergyDir = gatedhist->mkdir("TrificEnergy");
  trfGatedEnergyDir->cd();
  trfGatedEnergyList->Write();
  gatedhist->cd();
  TDirectory *trfGatedTimeDir = gatedhist->mkdir("TrificTime");
  trfGatedTimeDir->cd();
  trfGatedTimeList->Write();
  gatedhist->cd();
  TDirectory *tigGatedEnergyDir = gatedhist->mkdir("TigressEnergy");
  tigGatedEnergyDir->cd();
  tigGatedEnergyList->Write();
  gatedhist->cd();
  TDirectory *tigGatedEnergyABRingDir = gatedhist->mkdir("TigressAddbackEnergyRing");
  tigGatedEnergyABRingDir->cd();
  tigGatedEnergyABRingList->Write();
  gatedhist->cd();
  TDirectory *tigGatedEnergyABSumDir = gatedhist->mkdir("TigressAddbackEnergySum");
  tigGatedEnergyABSumDir->cd();
  tigGatedEnergyABSumList->Write();
  gatedhist->cd();
  TDirectory *tigGatedEnergyABDopSumDir = gatedhist->mkdir("TigressAddbackDopplerCorrectedEnergySum");
  tigGatedEnergyABDopSumDir->cd();
  tigGatedEnergyABDopSumList->Write();
  gatedhist->cd();
  TDirectory *tigGatedEnergyABRingDopSumDir = gatedhist->mkdir("TigressAddbackRingDopplerCorrectedEnergySum");
  tigGatedEnergyABRingDopSumDir->cd();
  tigGatedEnergyABRingDopSumList->Write();
  gatedhist->cd();
  TDirectory *tigGatedEnergyBkgSumDir = gatedhist->mkdir("TigressAddbackEnergyBackgroundSum");
  tigGatedEnergyBkgSumDir->cd();
  tigGatedEnergyBkgSumList->Write();
  gatedhist->cd();
  TDirectory *tigGatedTimeDir = gatedhist->mkdir("TigressTime");
  tigGatedTimeDir->cd();
  tigGatedTimeList->Write();
  gatedhist->cd();
  myfile->cd();
  myfile->Close();
  cout << "Histograms written, sorting complete!" << endl;

  cout << "Closing analysis and fragment files..." << endl;
  inputafile->Close();
  inputffile->Close();

  return 0;
}

int main(Int_t argc, Char_t** argv) 
{
  cout << "Starting sort code..." << endl;

  if(argc != 2) {
    cout << "Usage: SortData ConfigFileName" << endl;
    exit(1);
  }

  // Read in config file information
  const char* confName = argv[1];
  const Int_t verbosityLevel = 0;
  readConfig = new ReadConfig(confName,verbosityLevel);
  inafile = readConfig->GetAnalysisFileName();;
  inffile = readConfig->GetFragmentFileName();
  calfile = readConfig->GetCalibrationFileName();
  outfile = readConfig->GetOutputFileName();
  fragSort = readConfig->GetFragSortFlag();
  refCh = readConfig->GetTimeRefCh();
  timeGate = readConfig->GetTimeGateLimits();
  bgLow = readConfig->GetBGLowTimeGateLimits();
  bgHigh = readConfig->GetBGHighTimeGateLimits();
  bgLowScale = -0.5*(timeGate[1]-timeGate[0])/(bgLow[1]-bgLow[0]);
  bgHighScale = -0.5*(timeGate[1]-timeGate[0])/(bgHigh[1]-bgHigh[0]);
  beta = readConfig->GetDopplerBeta();
  ArGateFileName = readConfig->GetPid1FileName();
  KGateFileName = readConfig->GetPid2FileName();
  dE1 = readConfig->GetDeltaE1Range();
  dE2 = readConfig->GetDeltaE2Range();

  if(ReadGateFile(ArGateFileName,ArGates)) {
    cout << "Gate file " << ArGateFileName << " not read properly, exiting..." << endl;
    return 1;
  }
  if(ReadGateFile(KGateFileName,KGates)) {
    cout << "Gate file " << KGateFileName << " not read properly, exiting..." << endl;
    return 1;
  }

  // // output filenames, gate limits, etc.
  Print();

  // done reading config file

  InitializeHistograms();

  if(!SortData())
    cout << "Data sort SUCCESS!" << endl;
  else
    {
      cout << "ERROR encountered during data sort!" << endl;
      return 1;
    }

  delete readConfig;

  return 0;

}

void InitializeHistograms() 
{
  cout << "Initializing histograms for sort..." << endl;
  trfRawEnergyList = new TList(); 
  trfRawTimeList = new TList(); 
  chRawEnergyTimeList = new TList(); 
  tigRawEnergyList = new TList(); 
  tigRawEnergySumList = new TList(); 
  tigRawTimeList = new TList(); 
  trfGatedEnergyList = new TList(); 
  trfGatedTimeList = new TList(); 
  tigGatedEnergyList = new TList(); 
  tigGatedEnergyABSumList = new TList(); 
  tigGatedEnergyABRingList = new TList(); 
  tigGatedEnergyABDopSumList = new TList(); 
  tigGatedEnergyABRingDopSumList = new TList(); 
  tigGatedEnergyBkgSumList = new TList(); 
  tigGatedTimeList = new TList(); 
  rawList = new TList(); 
  gatedList = new TList();
  
  // raw histograms
  // TRIFIC
  trific_multiplicity = new TH1I("Trific multiplicity","Trific multiplicity",NTRF,0,NTRF); 
  rawList->Add(trific_multiplicity);
  trific_idmin = new TH1I("Trific minimum grid","Trific minimum grid",NTRF,0,NTRF); 
  rawList->Add(trific_idmin);
  trific_idmax = new TH1I("Trific maximum grid","Trific maximum grid",NTRF,0,NTRF); 
  rawList->Add(trific_idmax);
  trific_tigress_multiplicity = new TH2I("Trific-Tigress multiplicity","Trific-Tigress multiplicity",NTRF,0,NTRF,NTIG,0,NTIG); 
  rawList->Add(trific_tigress_multiplicity);
  trific_bragg = new TH2D("Trific Bragg curve","Trific Bragg curve",NTRF,0,NTRF,S8K,0,S16K); 
  rawList->Add(trific_bragg);
  trific_pid = new TH2D("Trific PID","Trific PID",S8K,0,S16K,S8K,0,S16K); 
  rawList->Add(trific_pid);
  trific_kvalue = new TH2I("Trific kValue","Trific kValue",NTRF,0,NTRF,S1K,0,S1K); 
  rawList->Add(trific_kvalue);

  // Detector spectra
  for(Int_t i=0; i<NTRF; i++) {
    Char_t hname[132];
    sprintf(hname,"Trific %d energy",i);
    trific_energy[i] = new TH1D(hname,hname,S16K,0,S16K); 
    trfRawEnergyList->Add(trific_energy[i]);
    sprintf(hname,"Trific %d time",i); 
    trific_time[i] = new TH1D(hname,hname,S8K,0,S8K); 
    trfRawTimeList->Add(trific_time[i]);
    sprintf(hname,"Trific grid %d-grid %d tdiff",refCh,i);
    trific_grrgrX_tdiff[i] = new TH1D(hname,hname,S16K,-S1K,S1K); 
    trfRawTimeList->Add(trific_grrgrX_tdiff[i]);
  }
  for(Int_t i=0; i<MAXCHAN; i++) {
    Char_t hname[132];
    sprintf(hname,"Channel %d energy vs MIDAS time",i);
    channel_vs_time[i] = new TH2D(hname,hname,S4K,0,S4K,S1K,0,S1K); 
    chRawEnergyTimeList->Add(channel_vs_time[i]);
  }

  // TIGRESS
  tigress_multiplicity = new TH1I("Tigress multiplicity","Tigress multiplicity",NTIG,0,NTIG); 
  rawList->Add(tigress_multiplicity);
  tigress_energy_sum = new TH1D("Tigress energy sum","Tigress energy sum",S8K,0,S16K); 
  tigRawEnergySumList->Add(tigress_energy_sum);
  tigress_energyAB_sum = new TH1D("Tigress addback energy sum","Tigress addback energy sum",S8K,0,S16K); 
  tigRawEnergySumList->Add(tigress_energyAB_sum);
  tigress_energyABDop_sum = new TH1D("Doppler-corrected Tigress addback energy sum","Doppler-corrected Tigress addback energy sum",S8K,0,S16K); 
  tigRawEnergySumList->Add(tigress_energyABDop_sum);
  Char_t refname[132];
  sprintf(refname,"Tigress-Trific grid %d tdiff",refCh);
  tigress_trificgr_tdiff = new TH1D(refname,refname,S16K,-S8K,S8K); 
  rawList->Add(tigress_trificgr_tdiff);
  sprintf(refname,"PIDAr-gated Tigress-Trific grid %d tdiff",refCh);
  tigress_trificgr_PIDAr_tdiff = new TH1D(refname,refname,S16K,-S8K,S8K); 
  rawList->Add(tigress_trificgr_PIDAr_tdiff);
  sprintf(refname,"PIDK-gated Tigress-Trific grid %d tdiff",refCh);
  tigress_trificgr_PIDK_tdiff = new TH1D(refname,refname,S16K,-S8K,S8K); 
  rawList->Add(tigress_trificgr_PIDK_tdiff);
  tigress_kvalue = new TH2I("Tigress kValue","Tigress kValue",NTIG,0,NTIG,S1K,0,S1K); 
  rawList->Add(tigress_kvalue);

  // Detector spectra
  for(Int_t i=0; i<NTIG; i++) {
    Char_t hname[132]; 
    sprintf(hname,"Tigress %d energy",i);
    tigress_energyAB[i] = new TH1D(hname,hname,S16K,0,S16K);  
    tigRawEnergyList->Add(tigress_energyAB[i]);
    sprintf(hname,"Tigress %d time",i);
    tigress_time[i] = new TH1D(hname,hname,S8K,0,S8K);  
    tigRawTimeList->Add(tigress_time[i]);
  }
  
  // gated histograms
  // TRIFIC
  gtrific_multiplicity = new TH1I("Gated Trific multiplicity","Gated Trific multiplicity",NTRF,0,NTRF); 
  gatedList->Add(gtrific_multiplicity);
  gtrific_idmin = new TH1I("Gated Trific minimum grid","Trific minimum grid",NTRF,0,NTRF); 
  gatedList->Add(gtrific_idmin);
  gtrific_idmax = new TH1I("Gated Trific maximum grid","Trific maximum grid",NTRF,0,NTRF); 
  gatedList->Add(gtrific_idmax);
  gtrific_bragg = new TH2D("Gated Trific bragg curve","Gated Trific bragg curve",NTRF,0,NTRF,S8K,0,S16K); 
  gatedList->Add(gtrific_bragg);
  gtrific_bragg_PIDAr = new TH2D("PIDAr-Gated Trific bragg curve","Gated Trific bragg curve",NTRF,0,NTRF,S8K,0,S16K); 
  gatedList->Add(gtrific_bragg_PIDAr);
  gtrific_bragg_PIDK = new TH2D("PIDK-Gated Trific bragg curve","Gated Trific bragg curve",NTRF,0,NTRF,S8K,0,S16K); 
  gatedList->Add(gtrific_bragg_PIDK);
  gtrific_pid = new TH2D("Trific PID coinc. with Tigress","Gated Trific PID coinc. with Tigress",S4K,0,S16K,S4K,0,S16K); 
  gatedList->Add(gtrific_pid);
  gtrific_pid_tcut = new TH2D("Gated Trific PID coinc. with Tigress","Gated Trific PID coinc. with Tigress",S4K,0,S16K,S4K,0,S16K); 
  gatedList->Add(gtrific_pid_tcut);

  // Detector spectra
  for(Int_t i=0; i<NTRF; i++) {
    Char_t hname[132];
    sprintf(hname,"Gated Trific %d energy",i);
    gtrific_energy[i] = new TH1D(hname,hname,S8K,0,S16K); 
    trfGatedEnergyList->Add(gtrific_energy[i]);
    sprintf(hname,"Gated Trific %d time",i);
    gtrific_time[i] = new TH1D(hname,hname,S8K,0,S8K); 
    trfGatedTimeList->Add(gtrific_time[i]);
  }

  // TIGRESS
  gtigress_multiplicity = new TH1I("Gated Tigress multiplicity","Gated Tigress multiplicity",NTIG,0,NTIG); 
  gatedList->Add(gtigress_multiplicity);

  // Addback energy
  gtigress_energyAB_sum = new TH1D("Gated Tigress addback energy sum","Gated Tigress addback energy sum",S8K,0,S16K); 
  tigGatedEnergyABSumList->Add(gtigress_energyAB_sum);
  gtigress_energyAB_bgLow_sum = new TH1D("Tigress addback energy bkg low sum","Tigress addback energy bkg low sum",S8K,0,S16K); 
  tigGatedEnergyBkgSumList->Add(gtigress_energyAB_bgLow_sum);
  gtigress_energyAB_bgHigh_sum = new TH1D("Tigress addback energy bkg high sum","Tigress addback energy bkg high sum",S8K,0,S16K); 
  tigGatedEnergyBkgSumList->Add(gtigress_energyAB_bgHigh_sum);
  gtigress_energyABSub_sum = new TH1D("Gated background-subtracted Tigress addback energy sum",
				      "Gated background-subtracted Tigress addback energy sum",S8K,0,S16K); 
  tigGatedEnergyABSumList->Add(gtigress_energyABSub_sum);

  // Doppler-corrected addback energy
  gtigress_energyABDop_sum = new TH1D("Gated Doppler-corrected Tigress addback energy sum",
				      "Gated Doppler-corrected Tigress addback energy sum",S8K,0,S16K); 
  tigGatedEnergyABDopSumList->Add(gtigress_energyABDop_sum);
  
  gtigress_energyABDop_bgLow_sum = new TH1D("Doppler-corrected Tigress addback energy bkg low sum",
					    "Doppler-corrected Tigress addback energy bkg low sum",S8K,0,S16K); 
  tigGatedEnergyBkgSumList->Add(gtigress_energyABDop_bgLow_sum);
  gtigress_energyABDop_bgHigh_sum = new TH1D("Doppler-corrected Tigress addback energy bkg high sum",
					     "Doppler-corrected Tigress addback energy bkg high sum",S8K,0,S16K); 
  tigGatedEnergyBkgSumList->Add(gtigress_energyABDop_bgHigh_sum);
  
  gtigress_energyABDopSub_sum = new TH1D("Gated Doppler-corrected background-subtracted Tigress addback energy sum",
					 "Gated Doppler-corrected background-subtracted Tigress addback energy sum",S8K,0,S16K); 
  tigGatedEnergyABDopSumList->Add(gtigress_energyABDopSub_sum);

  // PID gated
  gtigress_energyABDop_PIDAr_sum = new TH1D("PIDAr-gated Doppler-corrected Tigress addback energy sum",
					   "PIDAr-gated Gated Doppler-corrected Tigress addback energy sum",S8K,0,S16K); 
  tigGatedEnergyABDopSumList->Add(gtigress_energyABDop_PIDAr_sum);
  
  gtigress_energyABDop_bgLow_PIDAr_sum = new TH1D("PIDAr-gated Doppler-corrected Tigress addback energy bkg low sum",
							 "PIDAr-gated Doppler-corrected Tigress addback energy bkg low sum",
						  S8K,0,S16K); 
  tigGatedEnergyBkgSumList->Add(gtigress_energyABDop_bgLow_PIDAr_sum);
  gtigress_energyABDop_bgHigh_PIDAr_sum = new TH1D("PIDAr-gated Doppler-corrected Tigress addback energy bkg high sum",
						  "PIDAr-gated Doppler-corrected Tigress addback energy bkg high sum",S8K,0,S16K); 
  tigGatedEnergyBkgSumList->Add(gtigress_energyABDop_bgHigh_PIDAr_sum);
  
  gtigress_energyABDopSub_PIDAr_sum = new TH1D("PIDAr-gated Gated Doppler-corrected background-subtracted Tigress addback energy sum",
					      "PIDAr-gated Gated Doppler-corrected background-subtracted Tigress addback energy sum",
					       S8K,0,S16K); 
  tigGatedEnergyABDopSumList->Add(gtigress_energyABDopSub_PIDAr_sum);

  // PID K
  gtigress_energyABDop_PIDK_sum = new TH1D("PIDK-gated Doppler-corrected Tigress addback energy sum",
					   "PIDK-gated Gated Doppler-corrected Tigress addback energy sum",S8K,0,S16K); 
  tigGatedEnergyABDopSumList->Add(gtigress_energyABDop_PIDK_sum);
  
  gtigress_energyABDop_bgLow_PIDK_sum = new TH1D("PIDK-gated Doppler-corrected Tigress addback energy bkg low sum",
						 "PIDK-gated Doppler-corrected Tigress addback energy bkg low sum",S8K,0,S16K); 
  tigGatedEnergyBkgSumList->Add(gtigress_energyABDop_bgLow_PIDK_sum);

  gtigress_energyABDop_bgHigh_PIDK_sum = new TH1D("PIDK-gated Doppler-corrected Tigress addback energy bkg high sum",
						  "PIDK-gated Doppler-corrected Tigress addback energy bkg high sum",S8K,0,S16K); 
  tigGatedEnergyBkgSumList->Add(gtigress_energyABDop_bgHigh_PIDK_sum);
  
  gtigress_energyABDopSub_PIDK_sum = new TH1D("PIDK-gated Gated Doppler-corrected background-subtracted Tigress addback energy sum",
					      "PIDK-gated Gated Doppler-corrected background-subtracted Tigress addback energy sum",
					      S8K,0,S16K); 
  tigGatedEnergyABDopSumList->Add(gtigress_energyABDopSub_PIDK_sum);

  // EPID Ar
  gtigress_energyABDop_EPIDAr_sum = new TH1D("EPIDAr-gated Doppler-corrected Tigress addback energy sum",
					   "EPIDAr-gated Gated Doppler-corrected Tigress addback energy sum",S8K,0,S16K); 
  tigGatedEnergyABDopSumList->Add(gtigress_energyABDop_EPIDAr_sum);
  
  gtigress_energyABDop_bgLow_EPIDAr_sum = new TH1D("EPIDAr-gated Doppler-corrected Tigress addback energy bkg low sum",
						 "EPIDAr-gated Doppler-corrected Tigress addback energy bkg low sum",S8K,0,S16K); 
  tigGatedEnergyBkgSumList->Add(gtigress_energyABDop_bgLow_EPIDAr_sum);

  gtigress_energyABDop_bgHigh_EPIDAr_sum = new TH1D("EPIDAr-gated Doppler-corrected Tigress addback energy bkg high sum",
						  "EPIDAr-gated Doppler-corrected Tigress addback energy bkg high sum",S8K,0,S16K); 
  tigGatedEnergyBkgSumList->Add(gtigress_energyABDop_bgHigh_EPIDAr_sum);
  
  gtigress_energyABDopSub_EPIDAr_sum = new TH1D("EPIDAr-gated Gated Doppler-corrected background-subtracted Tigress addback energy sum",
					      "EPIDAr-gated Gated Doppler-corrected background-subtracted Tigress addback energy sum",
						S8K,0,S16K); 
  tigGatedEnergyABDopSumList->Add(gtigress_energyABDopSub_EPIDAr_sum);

  // EPID K
  gtigress_energyABDop_EPIDK_sum = new TH1D("EPIDK-gated Doppler-corrected Tigress addback energy sum",
					   "EPIDK-gated Gated Doppler-corrected Tigress addback energy sum",S8K,0,S16K); 
  tigGatedEnergyABDopSumList->Add(gtigress_energyABDop_EPIDK_sum);
  
  gtigress_energyABDop_bgLow_EPIDK_sum = new TH1D("EPIDK-gated Doppler-corrected Tigress addback energy bkg low sum",
						 "EPIDK-gated Doppler-corrected Tigress addback energy bkg low sum",S8K,0,S16K); 
  tigGatedEnergyBkgSumList->Add(gtigress_energyABDop_bgLow_EPIDK_sum);

  gtigress_energyABDop_bgHigh_EPIDK_sum = new TH1D("EPIDK-gated Doppler-corrected Tigress addback energy bkg high sum",
						  "EPIDK-gated Doppler-corrected Tigress addback energy bkg high sum",S8K,0,S16K); 
  tigGatedEnergyBkgSumList->Add(gtigress_energyABDop_bgHigh_EPIDK_sum);
  
  gtigress_energyABDopSub_EPIDK_sum = new TH1D("EPIDK-gated Gated Doppler-corrected background-subtracted Tigress addback energy sum",
					      "EPIDK-gated Gated Doppler-corrected background-subtracted Tigress addback energy sum",
					       S8K,0,S16K); 
  tigGatedEnergyABDopSumList->Add(gtigress_energyABDopSub_EPIDK_sum);


  // TIGRESS ring spectra
  for(Int_t i=0; i<NRING; i++) {
    Char_t hname[132];
    sprintf(hname,"Tigress ring %d energy sum",i);
    gtigress_energyABRing[i] = new TH1D(hname,hname,S8K,0,S16K);
    tigGatedEnergyABRingList->Add(gtigress_energyABRing[i]);
    sprintf(hname,"PIDAr-gated Gated Doppler-corrected Tigress addback ring %d energy sum",i);
    gtigress_energyABRingDop_PIDAr_sum[i] = new TH1D(hname,hname,S8K,0,S16K); 
    tigGatedEnergyABRingDopSumList->Add(gtigress_energyABRingDop_PIDAr_sum[i]);
    sprintf(hname,"PIDK-gated Gated Doppler-corrected Tigress addback ring %d energy sum",i);
    gtigress_energyABRingDop_PIDK_sum[i] = new TH1D(hname,hname,S8K,0,S16K); 
    tigGatedEnergyABRingDopSumList->Add(gtigress_energyABRingDop_PIDK_sum[i]);
  }
  
  // TIGRESS-TRIFIC time difference
  sprintf(refname,"Time-gated Tigress-Trific grid %d tdiff",refCh);
  gtigress_trificgr_tdiff = new TH1D(refname,refname,S16K,-S8K,S8K); 
  gatedList->Add(gtigress_trificgr_tdiff);
  sprintf(refname,"PIDAr-gated time-gated Tigress-Trific grid %d tdiff",refCh);
  gtigress_trificgr_PIDAr_tdiff = new TH1D(refname,refname,S16K,-S8K,S8K); 
  gatedList->Add(gtigress_trificgr_PIDAr_tdiff);
  sprintf(refname,"PIDK-gated time-gated Tigress-Trific grid %d tdiff",refCh);
  gtigress_trificgr_PIDK_tdiff = new TH1D(refname,refname,S16K,-S8K,S8K); 
  gatedList->Add(gtigress_trificgr_PIDK_tdiff);

  // Detector spectra
  for(Int_t i=0; i<NTIG; i++) {
    Char_t hname[132];
    sprintf(hname,"Gated Tigress addback %d energy",i);
    gtigress_energyAB[i] = new TH1D(hname,hname,S8K,0,S16K); 
    tigGatedEnergyList->Add(gtigress_energyAB[i]);
    sprintf(hname,"Gated Tigress %d time",i);
    gtigress_time[i] = new TH1D(hname,hname,S8K,0,S8K); 
    tigGatedTimeList->Add(gtigress_time[i]);
  }

  // TIGRESS-TRIFIC 2D coincidence spectra
  gtrific_mult_vs_tigress_energyABDop_sum = new TH2D("Trific multiplicity vs Tigress energy",
						     "Trific multiplicity vs Tigress energy",NTRF,0,NTRF,S8K,0,S16K);
  gatedList->Add(gtrific_mult_vs_tigress_energyABDop_sum);

  // Additional TRIFIC and TIGRESS-TRIFIC things
  gtrificEsum = new TH1D("Trific Esum","Trific Esum",S2K,0,S16K);
  gatedList->Add(gtrificEsum);
  gtrificEsum_vs_tigressE = new TH2D("Trific Esum vs Tigress E","Trific Esum vs Tigress E",S1K,0,S16K,S2K,0,S16K);
  gatedList->Add(gtrificEsum_vs_tigressE);
  gtrificEsum_vs_tigressE_PIDAr = new TH2D("Ar-gated Trific Esum vs Tigress E","Ar-gated Trific Esum vs Tigress E",S1K,0,S16K,S2K,0,S16K);
  gatedList->Add(gtrificEsum_vs_tigressE_PIDAr);
  gtrificEsum_vs_tigressE_PIDK = new TH2D("K-gated Trific Esum vs Tigress E","K-gated Trific Esum vs Tigress E",S1K,0,S16K,S2K,0,S16K);
  gatedList->Add(gtrificEsum_vs_tigressE_PIDK);
    
  cout << "Done initializing histograms, continuing sort..." << endl;

}

int ReadGateFile(const char* gfname, vector< vector<Double_t> > &gates) {

  // gates[detector#][2]
  // 0: low limit
  // 1: high limit
  gates.clear();

  ifstream gateFile(gfname);
  if(!gateFile.is_open()) {
    cout << "Opening gate file " << gfname << " failed, aborting!" << endl;
    return 1;
  } else {
    // cout << "Gate file " << gfname << " opened successfully..." << endl;
    string line;

    // Read header
    Int_t nline = 0;
    while (nline < 3) {
      getline(gateFile,line);
      nline++;
    }
    
    Double_t det, ELow, EHigh;
    while( gateFile >> det >> ELow >> EHigh ) {
      vector<Double_t> elim;

      // hacky way to rescale
      Double_t offset = 0.0;
      if( (ELow > 0) && (EHigh > 0) && (det <= 4) ) {
	Double_t range = EHigh-ELow;
	Double_t fact = 0.60; // range modification factor
	offset = 0.5*range*(1.0-fact);      
      }

      elim.push_back(ELow+offset);
      elim.push_back(EHigh-offset);
      gates.push_back(elim);
    }

    // // print out the gates
    // cout << "Energy gates defined in " << gfname << " are:" << endl;
    // cout << "Energy\tELow\tEHigh" << endl;
    // for(ULong_t i=0; i<gates.size(); i++) {
    //   cout << i << "\t" << gates[i][0] << "\t" << gates[i][1] << endl;
    // }    
    
    gateFile.close();
    return 0;
  }
  return 1;
}

void Print() {

  cout << "Analysis file: " << inafile << endl; 
  cout << "Fragment file: " << inffile << endl; 
  cout << "Calibration file: " << calfile << endl;
  cout << "Output file: " << outfile << endl;
  cout << "TRIFIC time reference channel: " << refCh << endl;
  cout << "Doppler correction beta: " << beta << endl;
  if(fragSort == true) cout << "Sorting from fragment tree: YES" << endl;
  else cout << "Sorting from fragment tree: NO" << endl;
  cout << "Time coincidence gate [ns]: ";
  for(vector<Double_t>::iterator it = timeGate.begin(); it!=timeGate.end(); ++it)
    cout << *it << " ";
  cout << endl;
  cout << "Low background gate [ns]: ";
  for(vector<Double_t>::iterator it = bgLow.begin(); it!=bgLow.end(); ++it)
    cout << *it << " ";
  cout << endl;
  cout << "High background gate [ns]: ";
  for(vector<Double_t>::iterator it = bgHigh.begin(); it!=bgHigh.end(); ++it)
    cout << *it << " ";
  cout << endl;

  cout << "Background subtraction scaling factors are: " << endl;
  cout << "----> Low : " << bgLowScale << endl;
  cout << "----> High: " << bgHighScale << endl;

  cout << "Argon PID file: " << ArGateFileName << endl; 
  // print out the gates
  cout << "Ar energy gates defined in " << ArGateFileName << " are:" << endl;
  cout << "Energy\tELow\tEHigh" << endl;
  for(ULong_t i=0; i<ArGates.size(); i++) {
    cout << i << "\t" << ArGates[i][0] << "\t" << ArGates[i][1] << endl;
  } 
  cout << "Potassium PID file: " << KGateFileName << endl; 
  // print out the gates
  cout << "K energy gates defined in " << KGateFileName << " are:" << endl;
  cout << "Energy\tELow\tEHigh" << endl;
  for(ULong_t i=0; i<KGates.size(); i++) {
    cout << i << "\t" << KGates[i][0] << "\t" << KGates[i][1] << endl;
  } 

  cout << "DeltaE1 region is: ";
  for(vector<Int_t>::iterator it = dE1.begin(); it!=dE1.end(); ++it)
    cout << *it << " ";
  cout << endl;
  cout << "DeltaE2 region is: ";
  for(vector<Int_t>::iterator it = dE2.begin(); it!=dE2.end(); ++it)
    cout << *it << " ";
  cout << endl;

  // getc(stdin);
  
}
