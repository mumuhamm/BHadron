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
	Mu1Charge_ = -999999;

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
  bTree_->SetBranchAddress(  "PVx"				  , &PVx_  );
  bTree_->SetBranchAddress(  "PVy"				  , &PVy_  );
  bTree_->SetBranchAddress(  "PVz"				  , &PVz_  );
  bTree_->SetBranchAddress(  "PVerrx"			  , &PVerrx_  );
  bTree_->SetBranchAddress(  "PVerry"			  , &PVerry_  );
  bTree_->SetBranchAddress(  "PVerrz"			  , &PVerrz_  );
}


