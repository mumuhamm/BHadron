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

};
#endif
