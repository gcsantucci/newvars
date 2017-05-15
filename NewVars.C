#define NewVars_cxx
#include "NewVars.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TVector3.h>
#include <fstream>

using namespace std;


int GetLoc(float mom, std::vector <double> &de, std::vector <double> &dx){
  float diff = 0;
  float diff_min = 100000;
  int loc = 0;
  for (int i=0; i < de.size(); i++){
    diff = mom - de.at(i);
    if (diff < diff_min && diff > 0) loc = i;
  }
  return loc;
}

float GetRange(float mom, float x1, float x2, float y1, float y2){
  float range = 0;
  range = y1 + ((y2-y1)/(x2-x1))*(mom - x1);
  return range;
}

void Get_dedx(std::vector <double> &de, std::vector <double> &dx){ 
  std::ifstream dedx;
  dedx.open("/Users/santucci/Dropbox/PhD/SK/fiTQun_analysis/Knu_muGamma/codes/mu_dedx_table.txt", std::ifstream::in);
  string line;
  if (dedx.is_open()){
    while ( getline (dedx,line) ){
      std::istringstream iss(line); // access line as a stream
      double a, b;
      iss >> a >> b;
      de.push_back(a);
      dx.push_back(b);
    }
    dedx.close();
  }
}

void NewVars::Loop(TString outputfile, int option){

  std::vector<double> de;
  std::vector<double> dx;
  Get_dedx(de, dx);
  int loc = 0; 
  float range = 0; 

  float temp = 0;

  TFile *outfile = new TFile(outputfile,"RECREATE");

  int foundmuon = 0;
  int foundelec = 0;
  int foundgamma = 0;
  int muoncounter = 0;
  int gammacounter = 0;
  int eleccounter = 0;

  bool notMuGamma = false;

  TVector3 truemomK;
  TVector3 truemomMu;
  TVector3 truemomEl;
  TVector3 truemomGamma;
  TVector3 recmomMu;
  TVector3 recmomEl;
  TVector3 recmomGamma;
  TVector3 onermomMu;
  TVector3 onermomEl;
  TVector3 trueMuElVec;
  TVector3 finalMuElVec;
  TVector3 recMuElVec;
  TVector3 onerMuElVec;

  // New Tree: treePDK
  TTree *treePDK;

  // variables in the tree        
  int loose = 0;
  int loose2 = 0;
  int loose3 = 0;
  int presel = 0;

  int MuGamma = 0;
  int MuGamma2 = 0;
  int EventID = 0;

  int FC;
  int FV;
  int FCFV;
  float FQdwall;
  float MCdwall;

  int nmichel;
  float tmichel[5];

  // Vertices Positions and Diferences
  float MCposMu[3];
  float MCfinalMu[3];
  float MCposGamma[3];
  float MCposEl[3];
  float FQposMu[3];
  float FQfinalMu[3];
  float SRfinalMu[3];
  float FQposGamma[3];
  float FQposEl[3];

  float MCrange;
  float FQrange;
  float SRrange;

  float MCposMuGamma;
  float FQposMuGamma;
  float MCposMuEl;
  float FQposMuEl;
  float MCfinalMuEl;
  float FQfinalMuEl;
  float SRposMuEl;
  float SRfinalMuEl;
  float MCimpact;
  float FQimpact;
  float SRimpact;

  // Momenta
  float MCmomK[3];
  float MCmomMu[3];
  float MCmomEl[3];
  float MCmomGamma[3];
  float FQmomMu[3];
  float FQmomEl[3];
  float FQmomGamma[3];

  float MCptotK;
  float MCptotMu;
  float MCptotEl;
  float MCptotGamma;
  float FQptotMu;
  float FQptotEl;
  float FQptotGamma;

  float MCtimeMu; float MCtimeGamma; float MCtimeEl;
  float FQtimeMu; float FQtimeGamma;
  float MCDeltaT; float FQDeltaT;

  //angle between the muon direction and the vector connecting Mu and Electron Vertices  
  float MCangMuElVec; float FQangMuElVec;
  float MCangMuElVec2; float FQangMuElVec2;

  //angle between TRUE and RECONS. directions for MU, Gamma and Electron            
  float angMu; float angEl; float angGamma;

  //Reconstructed angle between directions for Mu-Gamma, El-Gamma and Mu-Electron  
  float FQangMuEl; float FQangMuGamma;  float FQangElGamma;
  //True angle between directions for Mu-Gamma, El-Gamma and Mu-Electron                 
  float MCangMuEl;  float MCangMuGamma; float MCangElGamma;

  float MCthetaK; float MCthetaMu; float MCthetaEl; float MCthetaGamma;
  float FQthetaMu; float FQthetaEl; float FQthetaGamma;

  float MCphiK; float MCphiMu; float MCphiEl; float MCphiGamma;
  float FQphiMu; float FQphiEl; float FQphiGamma;

  float CCangMuGamma;
  int ccqe;
  int cc1p;
  int ccmp;
  int other;

  treePDK = new TTree("treePDK", "new Variables for PDK analysis");

  treePDK->Branch("EventID", &EventID, "EventID/I");

  treePDK->Branch("loose", &loose, "loose/I");
  treePDK->Branch("loose2", &loose2, "loose2/I");
  treePDK->Branch("loose3", &loose3, "loose3/I");
  treePDK->Branch("presel", &presel, "presel/I");
  treePDK->Branch("MuGamma", &MuGamma, "MuGamma/I");
  treePDK->Branch("MuGamma2", &MuGamma2, "MuGamma2/I");
  treePDK->Branch("foundgamma", &foundgamma, "foundgamma/I");
  treePDK->Branch("ccqe", &ccqe, "ccqe/I");
  treePDK->Branch("cc1p", &cc1p, "cc1p/I");
  treePDK->Branch("ccmp", &ccmp, "ccmp/I");
  treePDK->Branch("other", &other, "other/I");

  treePDK->Branch("FC", &FC, "FC/I");
  treePDK->Branch("FV", &FV, "FV/I");
  treePDK->Branch("FCFV", &FCFV, "FCFV/I");
  treePDK->Branch("FQdwall", &FQdwall, "FQdwall/F");
  treePDK->Branch("MCdwall", &MCdwall, "MCdwall/F");

  treePDK->Branch("nmichel", &nmichel, "nmichel/I");
  treePDK->Branch("tmichel", tmichel, "tmichel[5]/F");

  treePDK->Branch("MCposMu", MCposMu, "MCposMu[3]/F");
  treePDK->Branch("MCfinalMu", MCfinalMu, "MCfinalMu[3]/F");
  treePDK->Branch("MCposGamma", MCposGamma, "MCposGamma[3]/F");
  treePDK->Branch("MCposEl", MCposEl, "MCposEl[3]/F");
  treePDK->Branch("FQposMu", FQposMu, "FQposMu[3]/F");
  treePDK->Branch("FQfinalMu", FQfinalMu, "FQfinalMu[3]/F");
  treePDK->Branch("SRfinalMu", SRfinalMu, "SRfinalMu[3]/F");
  treePDK->Branch("FQposGamma", FQposGamma, "FQposGamma[3]/F");
  treePDK->Branch("FQposEl", FQposEl, "FQposEl[3]/F");

  treePDK->Branch("MCrange", &MCrange, "MCrange/F");
  treePDK->Branch("FQrange", &FQrange, "FQrange/F");
  treePDK->Branch("SRrange", &SRrange, "SRrange/F");

  treePDK->Branch("MCposMuGamma", &MCposMuGamma, "MCposMuGamma/F");
  treePDK->Branch("FQposMuGamma", &FQposMuGamma, "FQposMuGamma/F");
  treePDK->Branch("MCposMuEl", &MCposMuEl, "MCposMuEl/F");
  treePDK->Branch("FQposMuEl", &FQposMuEl, "FQposMuEl/F");
  treePDK->Branch("MCfinalMuEl", &MCfinalMuEl, "MCfinalMuEl/F");
  treePDK->Branch("FQfinalMuEl", &FQfinalMuEl, "FQfinalMuEl/F");
  treePDK->Branch("SRposMuEl", &SRposMuEl, "SRposMuEl/F");
  treePDK->Branch("SRfinalMuEl", &SRfinalMuEl, "SRfinalMuEl/F");
  treePDK->Branch("MCimpact", &MCimpact, "MCimpact/F");
  treePDK->Branch("FQimpact", &FQimpact, "FQimpact/F");
  treePDK->Branch("SRimpact", &SRimpact, "SRimpact/F");

  treePDK->Branch("MCDeltaT", &MCDeltaT, "MCDeltaT/F");
  treePDK->Branch("FQDeltaT", &FQDeltaT, "FQDeltaT/F");

  treePDK->Branch("MCtimeMu", &MCtimeMu, "MCtimeMu/F");
  treePDK->Branch("MCtimeEl", &MCtimeEl, "MCtimeEl/F");
  treePDK->Branch("MCtimeGamma", &MCtimeGamma, "MCtimeGamma/F");
  treePDK->Branch("FQtimeMu", &FQtimeMu, "FQtimeMu/F");
  treePDK->Branch("FQtimeGamma", &FQtimeGamma, "FQtimeGamma/F");

  treePDK->Branch("angMu", &angMu, "angMu/F");
  treePDK->Branch("angEl", &angEl, "angEl/F");
  treePDK->Branch("angGamma", &angGamma, "angGamma/F");
  treePDK->Branch("MCangMuEl", &MCangMuEl, "MCangMuEl/F");
  treePDK->Branch("MCangMuGamma", &MCangMuGamma, "MCangMuGamma/F");
  treePDK->Branch("MCangElGamma", &MCangElGamma, "MCangElGamma/F");
  treePDK->Branch("FQangMuEl", &FQangMuEl, "FQangMuEl/F");
  treePDK->Branch("FQangMuGamma", &FQangMuGamma, "FQangMuGamma/F");
  treePDK->Branch("FQangElGamma", &FQangElGamma, "FQangElGamma/F");
  treePDK->Branch("MCangMuElVec", &MCangMuElVec, "MCangMuElVec/F");
  treePDK->Branch("FQangMuElVec", &FQangMuElVec, "FQangMuElVec/F");
  treePDK->Branch("MCangMuElVec2", &MCangMuElVec2, "MCangMuElVec2/F");
  treePDK->Branch("FQangMuElVec2", &FQangMuElVec2, "FQangMuElVec2/F");

  treePDK->Branch("MCthetaK", &MCthetaK, "MCthetaK/F");
  treePDK->Branch("MCthetaMu", &MCthetaMu, "MCthetaMu/F");
  treePDK->Branch("MCthetaEl", &MCthetaEl, "MCthetaEl/F");
  treePDK->Branch("MCthetaGamma", &MCthetaGamma, "MCthetaGamma/F");
  treePDK->Branch("MCphiK", &MCphiK, "MCphiK/F");
  treePDK->Branch("MCphiMu", &MCphiMu, "MCphiMu/F");
  treePDK->Branch("MCphiEl", &MCphiEl, "MCphiEl/F");
  treePDK->Branch("MCphiGamma", &MCphiGamma, "MCphiGamma/F");
  treePDK->Branch("FQthetaMu", &FQthetaMu, "FQthetaMu/F");
  treePDK->Branch("FQthetaEl", &FQthetaEl, "FQthetaEl/F");
  treePDK->Branch("FQthetaGamma", &FQthetaGamma, "FQthetaGamma/F");
  treePDK->Branch("FQphiMu", &FQphiMu, "FQphiMu/F");
  treePDK->Branch("FQphiEl", &FQphiEl, "FQphiEl/F");
  treePDK->Branch("FQphiGamma", &FQphiGamma, "FQphiGamma/F");
  treePDK->Branch("CCangMuGamma", &CCangMuGamma, "CCangMuGamma/F");

  treePDK->Branch("MCmomK", MCmomK, "MCmomK[3]/F");
  treePDK->Branch("MCmomMu", MCmomMu, "MCmomMu[3]/F");
  treePDK->Branch("MCmomEl", MCmomEl, "MCmomEl[3]/F");
  treePDK->Branch("MCmomGamma", MCmomGamma, "MCmomGamma[3]/F");
  treePDK->Branch("MCptotK", &MCptotK, "MCptotK/F");
  treePDK->Branch("MCptotMu", &MCptotMu, "MCptotMu/F");
  treePDK->Branch("MCptotEl", &MCptotEl, "MCptotEl/F");
  treePDK->Branch("MCptotGamma", &MCptotGamma, "MCptotGamma/F");

  treePDK->Branch("FQmomMu", FQmomMu, "FQmomMu[3]/F");
  treePDK->Branch("FQmomEl", FQmomEl, "FQmomEl[3]/F");
  treePDK->Branch("FQmomGamma", FQmomGamma, "FQmomGamma[3]/F");
  treePDK->Branch("FQptotMu", &FQptotMu, "FQptotMu/F");
  treePDK->Branch("FQptotEl", &FQptotEl, "FQptotEl/F");
  treePDK->Branch("FQptotGamma", &FQptotGamma, "FQptotGamma/F");
  
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if( jentry % 50000 == 0 ) std::cout << "Event = " << jentry << std::endl;

    //initialize variables each entry
    foundmuon = 0;
    foundelec = 0;
    foundgamma = 0;
    notMuGamma = false;
    MuGamma = 0;
    MuGamma2 = 0;
    EventID = jentry;

    nmichel = 0;
    for (int i=0; i<5; i++) tmichel[i] = 0;

    MCdwall = wallv;
    MCptotK = 0;
    MCptotMu = 1000;
    MCthetaMu = 10;
    MCphiMu = 10;
    MCtimeMu = 10000;
    MCptotEl = 2000;
    MCthetaEl = 20;
    MCphiEl = 20;
    MCtimeEl = 20000;
    MCptotGamma = 500;
    MCthetaGamma = 10;
    MCphiGamma = 10;
    MCtimeGamma = -10;
    MCDeltaT = -500;
    MCimpact = -10;
    for (int i = 0; i < 3; i++){
      MCmomK[i] = 5000;
      MCmomMu[i] = 1000;
      MCposMu[i] = 20000;
      MCfinalMu[i] = 20000;
      MCmomEl[i] = 2000;
      MCposEl[i] = 30000;
      MCmomGamma[i] = 500;
      MCposGamma[i] = 10000;
    }

    MCposMuGamma = 0;
    MCposMuEl = -10;
    MCfinalMuEl = -10;

    MCangMuEl = 10;
    MCangMuGamma = 10;
    MCangElGamma = 10;
    MCangMuElVec = 10;
    MCangMuElVec2 = 10;

    MCrange = -5;

    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    if (option == 1){ // PDK MC 
      for (int ipart = 0; ipart < nscndprt; ipart++){
	if ( iprntprt[ipart] == 321 && iprtscnd[ipart] != -13 && iprtscnd[ipart] != 12 )
	  notMuGamma = true;
	if (iprtscnd[ipart] == -13 && iprntprt[ipart] == 321){
	  muoncounter++;
	  foundmuon = 1;
	  for (int i = 0; i < 3; i++){
	    MCmomMu[i] = pscnd[ipart][i];
	    MCposMu[i] = vtxscnd[ipart][i];
	  }
	  truemomMu.SetX(MCmomMu[0]);
	  truemomMu.SetY(MCmomMu[1]);
	  truemomMu.SetZ(MCmomMu[2]);
	  MCptotMu = truemomMu.Mag();
	  MCthetaMu = truemomMu.Theta();
	  MCphiMu = truemomMu.Phi();
	  MCtimeMu = tscnd[ipart];

	  loc = GetLoc(MCptotMu, de, dx);
	  range = GetRange(MCptotMu, de.at(loc), de.at(loc+1), dx.at(loc), dx.at(loc+1));
	  MCrange = range;
	  MCfinalMu[0] = MCposMu[0] + range*sin(MCthetaMu)*cos(MCphiMu);
	  MCfinalMu[1] = MCposMu[1] + range*sin(MCthetaMu)*sin(MCphiMu);
	  MCfinalMu[2] = MCposMu[2] + range*cos(MCthetaMu);
	} //look for muon

        else if (iprtscnd[ipart] == -11 && iprntprt[ipart] == -13){
	  eleccounter++;
	  foundelec = 1;
	  for (int i = 0; i < 3; i++){
	    MCmomEl[i] = pscnd[ipart][i];
	    MCposEl[i] = vtxscnd[ipart][i];
	  }
	  truemomEl.SetX(MCmomEl[0]);
	  truemomEl.SetY(MCmomEl[1]);
	  truemomEl.SetZ(MCmomEl[2]);
	  MCptotEl = truemomEl.Mag();
	  MCthetaEl = truemomEl.Theta();
	  MCphiEl = truemomEl.Phi();
	  MCtimeEl = tscnd[ipart];
	} //look for Michel-electron

	if (foundmuon && foundelec){
	  trueMuElVec.SetX(MCposMu[0] - MCposEl[0]);
	  trueMuElVec.SetY(MCposMu[1] - MCposEl[1]);
          trueMuElVec.SetZ(MCposMu[2] - MCposEl[2]);
	  finalMuElVec.SetX(MCfinalMu[0] - MCposEl[0]);
	  finalMuElVec.SetY(MCfinalMu[1] - MCposEl[1]);
          finalMuElVec.SetZ(MCfinalMu[2] - MCposEl[2]);
	  MCposMuEl = trueMuElVec.Mag();
	  MCfinalMuEl = finalMuElVec.Mag();
	  temp = 0;
	  MCimpact = 0;
	  for (int i = 0; i < 3; i++)
	    temp += (MCposMu[i] - MCposEl[i])*(MCmomMu[i]/MCptotMu);
	  for (int i = 0; i < 3; i++)
	    MCimpact += pow((MCposMu[i] - MCposEl[i]) - (temp*MCmomMu[i]/MCptotMu),2);
	  MCimpact = sqrt(MCimpact);
	}
      } //loop over secondary particles                         

      for (int iprim = 0; iprim < npar; iprim++){
	int pid = ipv[iprim];
	if ( pid == 1 ){
	  gammacounter++;
	  foundgamma = 1;
	  for (int i = 0; i < 3; i++){
	    MCposGamma[i] = posv[i];
	    MCmomGamma[i] = pmomv[iprim]*dirv[iprim][i];
	  }
	  truemomGamma.SetX(MCmomGamma[0]);
	  truemomGamma.SetY(MCmomGamma[1]);
	  truemomGamma.SetZ(MCmomGamma[2]);
	  MCptotGamma = truemomGamma.Mag();
	  MCthetaGamma = truemomGamma.Theta();
	  MCphiGamma = truemomGamma.Phi();
	  MCtimeGamma = 0;
	} //look for Gamma
	if ( pid == 11 ){
          for (int i = 0; i < 3; i++) MCmomK[i] = pmomv[iprim]*dirv[iprim][i];
	  truemomK.SetX(MCmomK[0]);
          truemomK.SetY(MCmomK[1]);
          truemomK.SetZ(MCmomK[2]);
          MCptotK = truemomK.Mag();
	  MCthetaK = truemomK.Theta();
	  MCphiK = truemomK.Phi();
        } //look for Kaon+
      }//loop over all primary particles              

      if(foundmuon && foundgamma && !notMuGamma){
	MuGamma = 1;
	MCDeltaT = MCtimeMu - MCtimeGamma;
	MCposMuGamma = sqrt( pow(MCposMu[0]-MCposGamma[0],2)+pow(MCposMu[1]-MCposGamma[1],2)+pow(MCposMu[2]-MCposGamma[2],2) );
	if (MCptotGamma < 11) MuGamma2 = 1;
	MCangMuEl = truemomMu.Angle(truemomEl);
	MCangMuGamma = truemomMu.Angle(truemomGamma);
	MCangElGamma = truemomEl.Angle(truemomGamma);
	MCangMuElVec = trueMuElVec.Angle(truemomMu);
	MCangMuElVec2 = trueMuElVec.Angle(truemomEl);
      } //look for MuGamma signal
    }//if option== 1 PDK                 

    else if (option == 0){ //atm nu MC
      ccqe = 0;
      cc1p = 0;
      ccmp = 0;
      other = 0;
      if (abs(mode)==1 || abs(mode)==2) ccqe = 1;
      else if (abs(mode)==11 || abs(mode)==12 || abs(mode)==13) cc1p = 1;
      else if(abs(mode)==21) ccmp = 1;
      else other = 1;
      for (int iprim = 0; iprim < npar; iprim++){
	int pid = ipv[iprim];
	if ( pid == 1 ){
	  foundgamma = 1;
	  for (int i = 0; i < 3; i++){
	    MCposGamma[i] = posv[i];
	    MCmomGamma[i] = pmomv[iprim]*dirv[iprim][i];
	  }
	  truemomGamma.SetX(MCmomGamma[0]);
	  truemomGamma.SetY(MCmomGamma[1]);
	  truemomGamma.SetZ(MCmomGamma[2]);
	  MCptotGamma = truemomGamma.Mag();
	  MCthetaGamma = truemomGamma.Theta();
	  MCphiGamma = truemomGamma.Phi();
	  MCtimeGamma = 0;
	} //look for Gamma   
	if ( pid == 5 || pid == 6 ){
	  foundmuon = 1;
	  for (int i = 0; i < 3; i++){
	    MCposMu[i] = posv[i];
	    MCmomMu[i] = pmomv[iprim]*dirv[iprim][i];
	  }
	  truemomMu.SetX(MCmomMu[0]);
	  truemomMu.SetY(MCmomMu[1]);
	  truemomMu.SetZ(MCmomMu[2]);
	  MCptotMu = truemomMu.Mag();
	  MCthetaMu = truemomMu.Theta();
	  MCphiMu = truemomMu.Phi();
	  MCtimeMu = 0;
	} //look for Mu   
	if (foundmuon && foundgamma)
	  MCangMuGamma = truemomMu.Angle(truemomGamma);
      } // loop over primaries
      for (int ipart = 0; ipart < nscndprt; ipart++){
	if ( abs(iprntprt[ipart]) == 13 && abs(iprtscnd[ipart]) == 11 && nmichel < 5){
          tmichel[nmichel] = tscnd[ipart];
          nmichel++;
	}
      }
    } // else option == 0 atm MC

    //set FQ reconstructed values:          
    FQptotMu = fqpmgmom1[0];
    FQptotEl = fq1rmom[1][1];
    FQptotGamma= fqpmgmom2[0];
    FQtimeMu = fqpmgt01[0];
    FQtimeGamma = fqpmgt02[0];
    FQDeltaT = FQtimeMu-FQtimeGamma;

    for (int i = 0; i < 3; i++){
      FQmomMu[i] = fqpmgmom1[0]*fqpmgdir1[0][i];
      FQmomEl[i] = fq1rmom[1][1]*fq1rdir[1][1][i];
      FQmomGamma[i] = fqpmgmom2[0]*fqpmgdir2[0][i];
      FQposMu[i] = fqpmgpos1[0][i];
      FQposGamma[i] = fqpmgpos2[0][i];
      FQposEl[i] = fq1rpos[1][1][i];
    }

    recmomMu.SetX(FQmomMu[0]);
    recmomMu.SetY(FQmomMu[1]);
    recmomMu.SetZ(FQmomMu[2]);
    recmomEl.SetX(FQmomEl[0]);
    recmomEl.SetY(FQmomEl[1]);
    recmomEl.SetZ(FQmomEl[2]);
    recmomGamma.SetX(FQmomGamma[0]);
    recmomGamma.SetY(FQmomGamma[1]);
    recmomGamma.SetZ(FQmomGamma[2]);
    recMuElVec.SetX(fqpmgpos1[0][0] - fq1rpos[1][1][0]);
    recMuElVec.SetY(fqpmgpos1[0][1] - fq1rpos[1][1][1]);
    recMuElVec.SetZ(fqpmgpos1[0][2] - fq1rpos[1][1][2]);
    onerMuElVec.SetX(fq1rpos[0][2][0] - fq1rpos[1][1][0]);
    onerMuElVec.SetY(fq1rpos[0][2][1] - fq1rpos[1][1][1]);
    onerMuElVec.SetZ(fq1rpos[0][2][2] - fq1rpos[1][1][2]);

    FQthetaMu = recmomMu.Theta();
    FQthetaEl = recmomEl.Theta();
    FQthetaGamma = recmomGamma.Theta();
    FQphiMu = recmomMu.Phi();
    FQphiEl = recmomEl.Phi();
    FQphiGamma = recmomGamma.Phi();

    loc = GetLoc(FQptotMu, de, dx);
    range = GetRange(FQptotMu, de.at(loc), de.at(loc+1), dx.at(loc), dx.at(loc+1));
    FQrange = range;
    for (int i=0; i<3;i++)
      FQfinalMu[i] = FQposMu[i] + range*fqpmgdir1[0][i];    
    loc = GetLoc(fq1rmom[0][2], de, dx);
    range = GetRange(fq1rmom[0][2], de.at(loc), de.at(loc+1), dx.at(loc), dx.at(loc+1));
    SRrange = range;
    for (int i=0; i<3;i++)
      SRfinalMu[i] = fq1rpos[0][2][i] + range*fq1rdir[0][2][i];

    FQposMuGamma = sqrt( pow(FQposMu[0]-FQposGamma[0],2)+pow(FQposMu[1]-FQposGamma[1],2)+pow(FQposMu[2]-FQposGamma[2],2) );
    FQposMuEl = recMuElVec.Mag();
    SRposMuEl = onerMuElVec.Mag();
    FQfinalMuEl = sqrt( pow(FQfinalMu[0]-FQposEl[0], 2) + pow(FQfinalMu[1]-FQposEl[1], 2) + pow(FQfinalMu[2]-FQposEl[2], 2) );
    SRfinalMuEl = sqrt(pow(SRfinalMu[0]-FQposEl[0], 2) + pow(SRfinalMu[1]-FQposEl[1], 2) + pow(SRfinalMu[2]-FQposEl[2], 2) );

    temp = 0;
    FQimpact = 0;
    for (int i = 0; i < 3; i++)
      temp += (FQposMu[i] - FQposEl[i])*fqpmgdir1[0][i];
    for (int i = 0; i < 3; i++)
      FQimpact += pow((FQposMu[i] - FQposEl[i]) - (temp*fqpmgdir1[0][i]), 2);
    FQimpact = sqrt(FQimpact);

    temp = 0;
    SRimpact = 0;
    for (int i = 0; i < 3; i++)
      temp += (fq1rpos[0][2][i] - fq1rpos[1][1][i])*fq1rdir[0][2][i];
    for (int i = 0; i < 3; i++)
      SRimpact += pow((fq1rpos[0][2][i] - fq1rpos[1][1][i]) - (temp*fq1rdir[0][2][i]),2);
    SRimpact = sqrt(SRimpact);

    angMu = truemomMu.Angle(recmomMu);
    angEl = truemomEl.Angle(recmomEl);
    angGamma = truemomGamma.Angle(recmomGamma);
    FQangMuEl = recmomMu.Angle(recmomEl);
    FQangMuGamma = recmomMu.Angle(recmomGamma);
    FQangElGamma = recmomEl.Angle(recmomGamma);
    FQangMuElVec = recMuElVec.Angle(recmomMu);
    FQangMuElVec2 = recMuElVec.Angle(recmomEl);

    if ( 1810 - fq1rpos[0][2][2] <= 1810 + fq1rpos[0][2][2] ){
      if ( 1810 - fq1rpos[0][2][2] <= 1690 - sqrt( pow(fq1rpos[0][2][0],2) + pow(fq1rpos[0][2][1],2)) )
	FQdwall = 1810 - fq1rpos[0][2][2];
      else
	FQdwall = 1690 - sqrt( pow(fq1rpos[0][2][0],2) + pow(fq1rpos[0][2][1],2));
    }
    else{
      if ( 1810 + fq1rpos[0][2][2] <= 1690 - sqrt( pow(fq1rpos[0][2][0],2) + pow(fq1rpos[0][2][1],2)) )
	FQdwall = 1810 + fq1rpos[0][2][2];
      else
	FQdwall= 1690 - sqrt( pow(fq1rpos[0][2][0],2) + pow(fq1rpos[0][2][1],2));
    }

    FC = 0;
    FV = 0;
    FCFV = 0;
    if(nhitac < 16)
      FC = 1;
    if(FQdwall > 200)
      FV = 1;
    if (FC==1 && FV==1)
      FCFV = 1;


    //pre-selection cuts:
    if ( abs(mode)==1 && (ipv[2]==5||ipv[2]==6) && pmomv[2]>170 && pmomv[2]<310 )
      presel = 1;
    else if ( abs(mode)==2 && (ipv[3]==5||ipv[3]==6) && pmomv[3]>50 && pmomv[3]<310 )    
      presel = 1;
    else if ( ((abs(mode)==11)||(abs(mode)==12)||(abs(mode)==13)||(abs(mode)==21)) && (ipv[2]==5||ipv[2]==6) && pmomv[2] < 310 )
      presel = 1;
    else if ( (abs(mode)==16)||(abs(mode)==17)||(abs(mode)>21) )
      presel = 1;
    else
      presel = 0;

    if (presel){
      if (potot > 300 && potot < 1600)
	presel = 1;
      else
	presel = 0;
    }
      

    loose2 = 0;
    loose3 = 0;
    //loose-cuts
    if (FCFV && fqmrnring[0]==1 && fqnse==2 && fq1rnll[0][1]>fq1rnll[0][2] && fq1rmom[0][2]>210 && fq1rmom[0][2]<270 && FQptotGamma>2 && FQptotGamma<12 && FQDeltaT>3 ) loose2 = 1;
    if ( FCFV && fqmrnring[0]==1 && fqnse==3 && fq1rnll[0][1]-fq1rnll[0][2]<100 && fq1rnll[1][1]-fq1rnll[1][2]>-50 && fq1rnll[2][1]-fq1rnll[2][2]<30 && fq1rmom[0][1]<20 && fq1rmom[1][2]>210 && fq1rmom[1][2]<270 && fq1rt0[1][2]-fq1rt0[0][1]>20 && fq1rt0[1][2]-fq1rt0[0][1]<120) loose3 = 1;
    if (loose2 == 1 || loose3 == 1) loose = 1;
    else
      loose = 0;

    treePDK->Fill();
      
  } //for loop over all entries (jentry)

  outfile->Write();
  outfile->Close();

} // Loop function
