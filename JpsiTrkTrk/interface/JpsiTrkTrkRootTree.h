#ifndef BHadron_JpsiTrkTrk_JpsiTrkTrkRootTree_h
#define BHadron_JpsiTrkTrk_JpsiTrkTrkRootTree_h

#include <string>
#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TChain.h>
#include <vector>

class JpsiTrkTrkRootTree{
  public: 
	
	  JpsiTrkTrkRootTree();
	
  
	  ~JpsiTrkTrkRootTree();
  
		
  void resetEntries();
  void writeFile();
  void createTree(const std::string filename);
  void readTree(const std::string filename);
  void readTree(std::vector<std::string> filenames);
  void getAngles(const double aa, const double bb, const double cc, const double dd);
  void getVtx(const double aa, const double bb, const double cc, const double dd, const double ee, const double ff, const double gg, const double hh, const double ii);
  void getDeDx(const double f1, const double f2, const int f3);
  void fill();
  void setBranchAddresses();
  
  public: 
  int Mu1Charge_;

  TFile* bFile_;
  TTree* bTree_;

  int    runNumber_;
  int    PUinteraction_;
  unsigned int eventNumber_;
  int    lumiSection_;

  int    isPV_;
  int    isBS_;
  int    NVertices_;
  double BSx_ ;
  double BSy_ ;
  double BSz_ ;
  double BSdx_ ;
  double BSdy_ ;
  double BSdz_ ;
  double BSdydz_;
  double BSdxdz_; 
  double PVx_ ;
  double PVy_ ;
  double PVz_ ;
  double PVerrx_ ;
  double PVerry_ ;
  double PVerrz_ ;
  
  
  int    Mu1SoftID_;
  int    Mu2SoftID_; 
  int    triggerbit_HLTDimuon4JpsiDisplaced_;
  int    triggerbit_HLTDimuon4JpsiNoVertexing_;
  int    triggerbit_HLTDimuon4JpsiTrkTrkDisplaced_;
  int    triggerbit_HLTDimuon4JpsiTrkBc_;
  int    MuonMultiplicity_ ;
  int    iPassedCutIdentBd_;
  int    iPassedCutIdent_;
  int    ihaveajpsi_;
  int    JpsiNumberOfCandidates_;
  double JpsiM_alone_;
  double JpsiPhi_alone_;
  double JpsiEta_alone_;
  double JpsiPt_alone_;
  double Mu1Pt_beffit_;
  double Mu1Pz_beffit_;
  double Mu1Eta_beffit_;
  double Mu1Phi_beffit_;
  double Mu2Pt_beffit_;
  double Mu2Pz_beffit_;
  double Mu2Eta_beffit_;
  double Mu2Phi_beffit_;
  int    JpsiMuMuDCA_beffit_;
  double JpsiCosDeltaAlpha_;
  double JpsiLxySigma_;
  double JpsiLxy_;
  double JpsiLxyOverPt_;
  double JpsiTrigVtxProb_;
  int    JpsiMuonCat1_;
  int    JpsiMuonCat2_;
  int    JpsiMuon1Cat_alone_;
  int    JpsiMuon2Cat_alone_;
  int    BdJpsiMuon1Cat_alone_;
  int    BdJpsiMuon2Cat_alone_;
  int    JpsiMu1nPixHits_alone_;
  int    JpsiMu2nPixHits_alone_;

  double K1Pt_beffit_;
  double K1Pz_beffit_;
  double K1Eta_beffit_;
  double K1Phi_beffit_;
  double PiPt_beffit_;
  double PiPz_beffit_;
  double PiEta_beffit_;
  double PiPhi_beffit_;
  int    Bd0NumberOfCandidatesAfterFit_;
  double Bd0FitM_;
  double Bd0FitEta_;
  double Bd0FitPt_;
  double Bd0FitPz_;
  double Bd0FitPhi_;



};
#endif
