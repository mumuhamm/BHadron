#include "/eos/user/a/almuhamm/01.MuonTech/WorkArea_BPhysics/CMSSW_14_0_3/src/BHadron/JpsiTrkTrk/interface/JpsiTrkTrkRootTree.h"
#include "BHadron/JpsiTrkTrk/interface/JpsiTrkTrkRootTree.h"
#include <string>
#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TChain.h>
#include <vector>
using namespace std;



JpsiTrkTrkRootTree::JpsiTrkTrkRootTree()
{
  resetEntries();
  bTree_ = 0;
  bFile_ = 0;
}

void JpsiTrkTrkRootTree::createTree(const std::string filename)
{
  bFile_ = new TFile (filename.c_str(), "RECREATE" );
  int bufsize = 256000;
  // create tree structure
  bTree_ = new TTree("bTree","bTree",bufsize);
  bTree_->Branch( "runNumber"   , &runNumber_   ,         "runNumber/I"   );
  bTree_->Branch( "eventNumber"   , &eventNumber_   ,       "eventNumber/i"   );
  bTree_->Branch( "lumiSection"   , &lumiSection_   ,       "lumiSection/I"   );
  bTree_->Branch( "PUinteraction"   , &PUinteraction_   ,         "PUinteraction/I"   );
  bTree_->Branch( "isPV"   , &isPV_   ,         "isPV/I"   );
  bTree_->Branch( "isBS"   , &isBS_   ,         "isBS/I"   );
  bTree_->Branch(  "NVertices"             , &NVertices_,                "NVertices/I");
  bTree_->Branch(  "BSx"				  , &BSx_,                              "BSx/D");
  bTree_->Branch(  "BSy"				  , &BSy_,                              "BSy/D");
  bTree_->Branch(  "BSz"				  , &BSz_,                              "BSz/D");
  bTree_->Branch(  "BSdx"                           , &BSdx_,                              "BSdx/D");
  bTree_->Branch(  "BSdy"                           , &BSdy_,                              "BSdy/D");
  bTree_->Branch(  "BSdz"                           , &BSdz_,                              "BSdz/D");
  bTree_->Branch(  "BSdxdz"                           , &BSdxdz_,                              "BSdxdz/D");
  bTree_->Branch(  "BSdydz"                           , &BSdydz_,                              "BSdydz/D");
  bTree_->Branch(  "PVx"				  , &PVx_,                              "PVx/D");
  bTree_->Branch(  "PVy"				  , &PVy_,                              "PVy/D");
  bTree_->Branch(  "PVz"				  , &PVz_,                              "PVz/D");
  bTree_->Branch(  "PVerrx"			  , &PVerrx_,                           "PVerrx/D");
  bTree_->Branch(  "PVerry"			  , &PVerry_,                           "PVerry/D");
  bTree_->Branch(  "PVerrz"			  , &PVerrz_,                           "PVerrz/D");
  bTree_->Branch("triggerbit_HLTDimuon4JpsiDisplaced", &triggerbit_HLTDimuon4JpsiDisplaced_,"triggerbit_HLTDimuon4JpsiDisplaced_/I"); 
  bTree_->Branch("triggerbit_HLTDimuon4JpsiNoVertexing",&triggerbit_HLTDimuon4JpsiNoVertexing_,"triggerbit_HLTDimuon4JpsiNoVertexing_/I");
  bTree_->Branch("triggerbit_HLTDimuon4JpsiTrkTrkDisplaced",&triggerbit_HLTDimuon4JpsiTrkTrkDisplaced_,"triggerbit_HLTDimuon4JpsiTrkTrkDisplaced_/I");
  bTree_->Branch("triggerbit_HLTDimuon4JpsiTrkBc",&triggerbit_HLTDimuon4JpsiTrkBc_,"triggerbit_HLTDimuon4JpsiTrkBc_/I");
  bTree_->Branch("MuonMultiplicity", &MuonMultiplicity_, "MuonMultiplicity/I");
  bTree_->Branch(  "iPassedCutIdent"		  , &iPassedCutIdent_,                  "iPassedCutIdent/I");
  bTree_->Branch(  "iPassedCutIdentBd"		  , &iPassedCutIdentBd_,                "iPassedCutIdentBd/I");  
  bTree_->Branch(  "ihaveajpsi"             , &ihaveajpsi_,                "ihaveajpsi/I");
  bTree_->Branch(  "JpsiNumberOfCandidates"          , &JpsiNumberOfCandidates_,            "JpsiNumberOfCandidates/I");
  bTree_->Branch(  "JpsiM_alone"			  , &JpsiM_alone_,                      "JpsiM_alone/D");
  bTree_->Branch(  "JpsiPhi_alone"		  , &JpsiPhi_alone_,                    "JpsiPhi_alone/D");
  bTree_->Branch(  "JpsiEta_alone"		  , &JpsiEta_alone_,                    "JpsiEta_alone/D");
  bTree_->Branch(  "JpsiPt_alone"		  , &JpsiPt_alone_,                     "JpsiPt_alone/D");
  bTree_->Branch(  "Mu1Pt_beffit"			  , &Mu1Pt_beffit_,                       "Mu1Pt_beffit/D");
  bTree_->Branch(  "Mu1Pz_beffit"			  , &Mu1Pz_beffit_,                       "Mu1Pz_beffit/D");
  bTree_->Branch(  "Mu1Eta_beffit"			  , &Mu1Eta_beffit_,                      "Mu1Eta_beffit/D");
  bTree_->Branch(  "Mu1Phi_beffit"			  , &Mu1Phi_beffit_,                      "Mu1Phi_beffit/D");
  bTree_->Branch(  "Mu2Pt_beffit"			  , &Mu2Pt_beffit_,                       "Mu2Pt_beffit/D");
  bTree_->Branch(  "Mu2Pz_beffit"			  , &Mu2Pz_beffit_,                       "Mu2Pz_beffit/D");
  bTree_->Branch(  "Mu2Eta_beffit"			  , &Mu2Eta_beffit_,                      "Mu2Eta_beffit/D");
  bTree_->Branch(  "Mu2Phi_beffit"			  , &Mu2Phi_beffit_,                      "Mu2Phi_beffit/D");
  bTree_->Branch(  "JpsiCosDeltaAlpha"         , &JpsiCosDeltaAlpha_, "JpsiCosDeltaAlpha/D");
  bTree_->Branch(  "JpsiLxySigma"         , &JpsiLxySigma_, "JpsiLxySigma/D");
  bTree_->Branch(  "JpsiLxy"         , &JpsiLxy_, "JpsiLxy/D");
  bTree_->Branch(  "JpsiLxyOverPt"         , &JpsiLxyOverPt_, "JpsiLxyOverPt/D");
  bTree_->Branch(  "JpsiTrigVtxProb"         , &JpsiTrigVtxProb_,"JpsiTrigVtxProb/D"); 
  bTree_->Branch(  "JpsiMuon1Cat_alone"		  , &JpsiMuon1Cat_alone_,               "JpsiMuon1Cat_alone/I");
  bTree_->Branch(  "JpsiMuon2Cat_alone"		  , &JpsiMuon2Cat_alone_,               "JpsiMuon2Cat_alone/I");
  bTree_->Branch(  "BdJpsiMuon1Cat_alone"		  , &BdJpsiMuon1Cat_alone_,               "BdJpsiMuon1Cat_alone/I");
  bTree_->Branch(  "BdJpsiMuon2Cat_alone"		  , &BdJpsiMuon2Cat_alone_,               "BdJpsiMuon2Cat_alone/I");
  bTree_->Branch(  "JpsiMuonCat1"		  , &JpsiMuonCat1_ , "JpsiMuonCat1/I"  );
  bTree_->Branch(  "JpsiMuonCat2"		  , &JpsiMuonCat2_ , "JpsiMuonCat2/I"  ); 
  bTree_->Branch(  "JpsiMu1nPixHits_alone"         , &JpsiMu1nPixHits_alone_,               "JpsiMu1nPixHits_alone/I");
  bTree_->Branch(  "JpsiMu2nPixHits_alone"         , &JpsiMu2nPixHits_alone_,               "JpsiMu2nPixHits_alone/I");

  bTree_->Branch(  "K1Pt_beffit"			  , &K1Pt_beffit_,                       "K1Pt_beffit/D");
  bTree_->Branch(  "K1Pz_beffit"			  , &K1Pz_beffit_,                       "K1Pz_beffit/D");
  bTree_->Branch(  "K1Eta_beffit"			  , &K1Eta_beffit_,                      "K1Eta_beffit/D");
  bTree_->Branch(  "K1Phi_beffit"			  , &K1Phi_beffit_,                      "K1Phi_beffit/D");

  bTree_->Branch(  "PiPt_beffit"			  , &PiPt_beffit_,                       "PiPt_beffit/D");
  bTree_->Branch(  "PiPz_beffit"			  , &PiPz_beffit_,                       "PiPz_beffit/D");
  bTree_->Branch(  "PiEta_beffit"			  , &PiEta_beffit_,                      "PiEta_beffit/D");
  bTree_->Branch(  "PiPhi_beffit"			  , &PiPhi_beffit_,                      "PiPhi_beffit/D");
  bTree_->Branch(  "Bd0NumberOfCandidatesAfterFit"          , &Bd0NumberOfCandidatesAfterFit_,            "Bd0NumberOfCandidatesAfterFit/I");
  bTree_->Branch(  "Bd0FitM"			  , &Bd0FitM_,                           "Bd0FitM/D");
  bTree_->Branch(  "Bd0FitEta"			  , &Bd0FitEta_,                         "Bd0FitEta/D");
  bTree_->Branch(  "Bd0FitPt"			  , &Bd0FitPt_,                          "Bd0FitPt/D");
  bTree_->Branch(  "Bd0FitPz"			  , &Bd0FitPz_,                          "Bd0FitPz/D");
  bTree_->Branch(  "Bd0FitPhi"			  , &Bd0FitPhi_,                         "Bd0FitPhi/D");
}

JpsiTrkTrkRootTree::~JpsiTrkTrkRootTree()
{}
void JpsiTrkTrkRootTree::writeFile()
{
  bFile_->Write();
  bFile_->Close();

}
void JpsiTrkTrkRootTree::resetEntries()
{       
	runNumber_ =  -9999999;
        eventNumber_ =  -9999999;
        lumiSection_ =  -9999999;
        PUinteraction_ = -9999999;
	NVertices_ = -9999999;
	isPV_ = -9999999;
        isBS_ = -9999999;
	BSx_ = -9999999;
        BSy_ = -9999999;
        BSz_ = -9999999;
        BSdx_ = -9999999;
        BSdy_ = -9999999;
        BSdz_ = -9999999;
        BSdxdz_ = -9999999;
        BSdydz_ = -9999999;
        PVx_ = -9999999;
        PVy_ = -9999999;
        PVz_ = -9999999;
        PVerrx_ = -9999999;
        PVerry_ = -9999999;
        PVerrz_ = -9999999;
	triggerbit_HLTDimuon4JpsiDisplaced_ = -9999999;
	triggerbit_HLTDimuon4JpsiNoVertexing_ = -9999999;
	triggerbit_HLTDimuon4JpsiTrkTrkDisplaced_ = -9999999;
	triggerbit_HLTDimuon4JpsiTrkBc_ = -9999999;
        MuonMultiplicity_ = -9999999;
	iPassedCutIdent_ = -9999999;
        iPassedCutIdentBd_ = -9999999;
	Mu1Charge_ = -999999;
        ihaveajpsi_=  -9999999;
	JpsiNumberOfCandidates_ = 0;
	JpsiM_alone_ = -9999999;
        JpsiPhi_alone_ = -9999999;
        JpsiEta_alone_ = -9999999;
        JpsiPt_alone_ = -9999999;
        Mu1Pt_beffit_ = -9999999;
        Mu1Pz_beffit_ = -9999999;
        Mu1Eta_beffit_ = -9999999;
        Mu1Phi_beffit_ = -9999999;
        Mu2Pt_beffit_ = -9999999;
        Mu2Pz_beffit_ = -9999999;
        Mu2Eta_beffit_ = -9999999;
        Mu2Phi_beffit_ = -9999999;
	JpsiCosDeltaAlpha_ = -9999999;
        JpsiLxySigma_ = -9999999;
        JpsiLxy_ = -9999999;
	JpsiLxyOverPt_= -9999999;
        JpsiTrigVtxProb_ = -9999999;
	JpsiMuonCat1_ = -9999999;
        JpsiMuonCat2_ = -9999999;
        JpsiMuon1Cat_alone_ = -9999999;
        JpsiMuon2Cat_alone_ = -9999999;
        BdJpsiMuon1Cat_alone_ = -9999999;
        BdJpsiMuon2Cat_alone_ = -9999999;
 	JpsiMu1nPixHits_alone_ = -9999999;
        JpsiMu2nPixHits_alone_ = -9999999;
        K1Pt_beffit_ = -9999999;
        K1Pz_beffit_ = -9999999;
        K1Eta_beffit_ = -9999999;
        K1Phi_beffit_ = -9999999;
        PiPt_beffit_ = -9999999;
        PiPz_beffit_ = -9999999;
        PiEta_beffit_ = -9999999;
        PiPhi_beffit_ = -9999999;
        Bd0NumberOfCandidatesAfterFit_ =  0;
        Bd0FitM_ = -9999999;
	Bd0FitEta_ = -9999999;
        Bd0FitPt_ = -9999999;
        Bd0FitPz_ = -9999999;
        Bd0FitPhi_ = -9999999;

}
void JpsiTrkTrkRootTree::getDeDx(const double f1, const double f2, const int f3)
{  
}
void JpsiTrkTrkRootTree::getVtx(const double aa, const double bb, const double cc, const double dd, const double ee, const double ff,
                                 const double gg, const double hh, const double ii)
{
	
  BSx_ = aa;
  BSy_ = bb;
  BSz_ = cc;
  PVx_ = dd;
  PVy_ = ee;
  PVz_ = ff;
  PVerrx_ = gg;
  PVerry_ = hh;
  PVerrz_ = ii;

}
void JpsiTrkTrkRootTree::getAngles(const double aa, const double bb, const double cc, const double dd)
{
 
}
void JpsiTrkTrkRootTree::fill()
{
  bTree_->Fill();
}
void JpsiTrkTrkRootTree::readTree(const std::string filename)
{

  
  bFile_ = new TFile (filename.c_str(), "READ" );
  // create tree structure
  bTree_ =  (TTree*) bFile_->Get("bTree");
  setBranchAddresses();
}
void JpsiTrkTrkRootTree::readTree(std::vector<std::string> filenames){
  TChain * myChain = new TChain("bTree");
  //  for(int i=0;i<filenames.size();++i) {
    //    std::string filenome = filenames[i] ;
  //    myChain->Add(filenome,-1);
  //  }

  for(std::vector<std::string>::iterator it = filenames.begin(); it != filenames.end(); it++){
    myChain->Add( (*it).c_str());
  }

  bTree_ = myChain;
  setBranchAddresses();
}

void JpsiTrkTrkRootTree::setBranchAddresses(){

  bTree_->SetBranchAddress(  "runNumber"             , &runNumber_  );
  bTree_->SetBranchAddress(  "eventNumber"             , &eventNumber_  );
  bTree_->SetBranchAddress(  "lumiSection"             , &lumiSection_  );
  bTree_->SetBranchAddress(  "PUinteraction"             , &PUinteraction_  );
  bTree_->SetBranchAddress(  "isPV"                               , &isPV_  );
  bTree_->SetBranchAddress(  "isBS"                               , &isBS_  );
  bTree_->SetBranchAddress(  "NVertices"             , &NVertices_  ); 
  bTree_->SetBranchAddress(  "BSx"				  , &BSx_  );
  bTree_->SetBranchAddress(  "BSy"				  , &BSy_  );
  bTree_->SetBranchAddress(  "BSz"				  , &BSz_  );
  bTree_->SetBranchAddress(  "BSdx"                                , &BSdx_  );
  bTree_->SetBranchAddress(  "BSdy"                                , &BSdy_  );
  bTree_->SetBranchAddress(  "BSdz"                                , &BSdz_  );
  bTree_->SetBranchAddress("PVx"				  , &PVx_  );
  bTree_->SetBranchAddress("PVy"				  , &PVy_  );
  bTree_->SetBranchAddress("PVz"				  , &PVz_  );
  bTree_->SetBranchAddress("PVerrx"			  , &PVerrx_  );
  bTree_->SetBranchAddress("PVerry"			  , &PVerry_  );
  bTree_->SetBranchAddress("PVerrz"			  , &PVerrz_  );
  bTree_->SetBranchAddress("triggerbit_HLTDimuon4JpsiDisplaced",&triggerbit_HLTDimuon4JpsiDisplaced_);
  bTree_->SetBranchAddress("triggerbit_HLTDimuon4JpsiNoVertexing",&triggerbit_HLTDimuon4JpsiNoVertexing_);
  bTree_->SetBranchAddress("triggerbit_HLTDimuon4JpsiTrkTrkDisplaced_",&triggerbit_HLTDimuon4JpsiTrkTrkDisplaced_);
  bTree_->SetBranchAddress("triggerbit_HLTDimuon4JpsiTrkBc_",&triggerbit_HLTDimuon4JpsiTrkBc_);
  bTree_->SetBranchAddress("iPassedCutIdent"		  , &iPassedCutIdent_  );
  bTree_->SetBranchAddress("iPassedCutIdentBd"		  , &iPassedCutIdentBd_  );
  bTree_->SetBranchAddress("ihaveajpsi"             , &ihaveajpsi_);
  bTree_->SetBranchAddress(  "JpsiNumberOfCandidates"        , &JpsiNumberOfCandidates_);
  bTree_->SetBranchAddress(  "JpsiM_alone"			  , &JpsiM_alone_  );
  bTree_->SetBranchAddress(  "JpsiPhi_alone"		  , &JpsiPhi_alone_  );
  bTree_->SetBranchAddress(  "JpsiEta_alone"		  , &JpsiEta_alone_  );
  bTree_->SetBranchAddress(  "JpsiPt_alone"		  , &JpsiPt_alone_  );
  bTree_->SetBranchAddress(  "Mu1Pt_beffit"			  , &Mu1Pt_beffit_  );
  bTree_->SetBranchAddress(  "Mu1Pz_beffit"			  , &Mu1Pz_beffit_  );
  bTree_->SetBranchAddress(  "Mu1Eta_beffit"			  , &Mu1Eta_beffit_  );
  bTree_->SetBranchAddress(  "Mu1Phi_beffit"			  , &Mu1Phi_beffit_  );
  bTree_->SetBranchAddress(  "Mu2Pt_beffit"			  , &Mu2Pt_beffit_  );
  bTree_->SetBranchAddress(  "Mu2Pz_beffit"			  , &Mu2Pz_beffit_  );
  bTree_->SetBranchAddress(  "Mu2Eta_beffit"			  , &Mu2Eta_beffit_  );
  bTree_->SetBranchAddress(  "Mu2Phi_beffit"			  , &Mu2Phi_beffit_  );
  bTree_->SetBranchAddress(  "JpsiCosDeltaAlpha"         , &JpsiCosDeltaAlpha_);
  bTree_->SetBranchAddress(  "JpsiLxySigma"         , &JpsiLxySigma_);
  bTree_->SetBranchAddress(  "JpsiLxy"         , &JpsiLxy_);
  bTree_->SetBranchAddress(  "JpsiLxyOverPt"         , &JpsiLxyOverPt_);
  bTree_->SetBranchAddress(  "MuonMultiplicity",&MuonMultiplicity_);
  bTree_->SetBranchAddress(  "JpsiTrigVtxProb"         , &JpsiTrigVtxProb_); 
  bTree_->SetBranchAddress(  "JpsiMuon1Cat_alone"		  , &JpsiMuon1Cat_alone_  );
  bTree_->SetBranchAddress(  "JpsiMuon2Cat_alone"		  , &JpsiMuon2Cat_alone_  );
  bTree_->SetBranchAddress(  "BdJpsiMuon1Cat_alone"		  , &BdJpsiMuon1Cat_alone_  );
  bTree_->SetBranchAddress(  "BdJpsiMuon2Cat_alone"		  , &BdJpsiMuon2Cat_alone_  ); 
  bTree_->SetBranchAddress(  "JpsiMuonCat1"		  , &JpsiMuonCat1_  );
  bTree_->SetBranchAddress(  "JpsiMuonCat2"		  , &JpsiMuonCat2_  );
  bTree_->SetBranchAddress(  "JpsiMu1nPixHits_alone"               , &JpsiMu1nPixHits_alone_  );
  bTree_->SetBranchAddress(  "JpsiMu2nPixHits_alone"               , &JpsiMu2nPixHits_alone_  );
  bTree_->SetBranchAddress(  "K1Pt_beffit"			  , &K1Pt_beffit_  );
  bTree_->SetBranchAddress(  "K1Pz_beffit"			  , &K1Pz_beffit_  );
  bTree_->SetBranchAddress(  "K1Eta_beffit"			  , &K1Eta_beffit_  );
  bTree_->SetBranchAddress(  "K1Phi_beffit"			  , &K1Phi_beffit_  );
  bTree_->SetBranchAddress(  "PiPt_beffit"			  , &PiPt_beffit_  );
  bTree_->SetBranchAddress(  "PiPz_beffit"			  , &PiPz_beffit_  );
  bTree_->SetBranchAddress(  "PiEta_beffit"			  , &PiEta_beffit_  );
  bTree_->SetBranchAddress(  "PiPhi_beffit"			  , &PiPhi_beffit_  );
  bTree_->SetBranchAddress(  "Bd0NumberOfCandidatesAfterFit"        , &Bd0NumberOfCandidatesAfterFit_);
  bTree_->SetBranchAddress(  "Bd0FitM"			  , &Bd0FitM_  );
  bTree_->SetBranchAddress(  "Bd0FitEta"			  , &Bd0FitEta_  );
  bTree_->SetBranchAddress(  "Bd0FitPt"			  , &Bd0FitPt_  );
  bTree_->SetBranchAddress(  "Bd0FitPz"			  , &Bd0FitPz_  );
  bTree_->SetBranchAddress(  "Bd0FitPhi"			  , &Bd0FitPhi_  );
}


