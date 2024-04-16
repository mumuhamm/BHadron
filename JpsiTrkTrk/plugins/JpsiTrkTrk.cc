// -*- C++ -*-
// Package:    BHadron/JpsiTrkTrk
// Class:      JpsiTrkTrk
//class JpsiTrkTrk JpsiTrkTrk.cc BHadron/JpsiTrkTrk/plugins/JpsiTrkTrk.cc

// Original Author:  Alibordi Muhammad
//         Created:  Sat, 13 Apr 2024 13:11:49 GMT



#include <memory>
#include "BHadron/JpsiTrkTrk/interface/JpsiTrkTrkRootTree.h"
#include "BHadron/JpsiTrkTrk/interface/KinematicFitInterface.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "BHadron/JpsiTrkTrk/interface/JpsiTrkTrk.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include <cstddef>
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"
#include "FWCore/Utilities/interface/Exception.h"
#include <cfloat>
#include <string>
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "DataFormats/PatCandidates/interface/IsolatedTrack.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/Registry.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/Association.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/TriggerNamesService.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
//#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include <RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h>
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidateFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonTrackLinks.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/CaloMuon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Candidate/interface/OverlapChecker.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TVector.h"
#include "TLorentzRotation.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicVertex.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticleFwd.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/TrackReco/interface/DeDxHit.h"
#include "RecoTracker/DeDx/interface/DeDxTools.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/ProjectedSiStripRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexProducerAlgorithm.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/MuonSeed/interface/L2MuonTrajectorySeed.h"
#include "DataFormats/MuonSeed/interface/L2MuonTrajectorySeedCollection.h"
#include "DataFormats/MuonSeed/interface/L3MuonTrajectorySeed.h"
#include "DataFormats/MuonSeed/interface/L3MuonTrajectorySeedCollection.h"
#include "MuonAnalysis/MuonAssociators/interface/PropagateToMuon.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "MuonAnalysis/MuonAssociators/interface/L1MuonMatcherAlgo.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/JetCorrFactors.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include <iostream>
#include <TMath.h>
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/PatCandidates/interface/IsolatedTrack.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
using namespace reco;
using namespace edm;
using namespace std;
using namespace pat;





JpsiTrkTrk::JpsiTrkTrk(const edm::ParameterSet& iConfig):
	theConfig_(iConfig),
	nominalJpsiMass( 3.096916 ),
        nominalPhiMass(1.019 ),
        nominalElectronMass(0.00051099893),
        nominalMuonMass(0.1056583),
        nominalKaonMass(0.493677),
        nominalPionMass(0.139570),
        nominalKstarMass(0.892),
        nominalBplusMass(5.2792)
{
	isMCstudy_ = iConfig.getParameter<bool>("isMCstudy");
	genParticlesLabel                 = iConfig.getParameter<InputTag>("genParticlesLabel");
        genParticlesTok                   = consumes<edm::View<reco::GenParticle>>(genParticlesLabel);
        MuonTag                           = iConfig.getParameter<edm::InputTag>("MuonTag");
        MuonTagTok                        = consumes<edm::View<pat::Muon>>(MuonTag);
        JetCollection                     = iConfig.getParameter<edm::InputTag>("JetCollection");
        JetCollectionTok                  = consumes<edm::View<pat::Jet>>(JetCollection);
        PUInfo                            = iConfig.getParameter<InputTag>("PUInfo");
        PUInfoTok                         = consumes<edm::View<PileupSummaryInfo>>(PUInfo);
        vertexBeamSpot                    = iConfig.getParameter<edm::InputTag>("vertexBeamSpot");
        vertexBeamSpotTok                 = consumes<reco::BeamSpot>(vertexBeamSpot);
        primaryvertex                     = iConfig.getParameter<edm::InputTag>("primaryvertex");
        primaryvertexTok                  = consumes<edm::View<reco::Vertex>>(primaryvertex);
        triggerresults                    = iConfig.getParameter<edm::InputTag>("triggerresults");
        triggerresultsTok                 = consumes<edm::TriggerResults>(triggerresults);
        ElectronTag                       = iConfig.getParameter<edm::InputTag>("ElectronTag");
        ElectronTagTok                    = consumes<edm::View<pat::Electron>>(ElectronTag);
        track                             = iConfig.getParameter<edm::InputTag>("track");
        trackLabelK                       = consumes<edm::View<pat::PackedCandidate>>(track);
        isotrack                          = iConfig.getParameter<edm::InputTag>("isotrack");
        isotrackTok                       = consumes<edm::View<pat::IsolatedTrack>>(isotrack);
        trackBuilderToken                 = esConsumes(edm::ESInputTag("", "TransientTrackBuilder"));
//=============================================================================================
  StoreDeDxInfo_ = iConfig.getParameter<bool>("StoreDeDxInfo");
  JpsiMassWindowBeforeFit_ = iConfig.getParameter<double>("JpsiMassWindowBeforeFit");

  BsLowerMassCutBeforeFit_  = iConfig.getParameter<double>("BsLowerMassCutBeforeFit");
  BsUpperMassCutBeforeFit_  = iConfig.getParameter<double>("BsUpperMassCutBeforeFit");
  BsLowerMassCutAfterFit_  = iConfig.getParameter<double>("BsLowerMassCutAfterFit");
  BsUpperMassCutAfterFit_  = iConfig.getParameter<double>("BsUpperMassCutAfterFit");

  JpsiMassWindowAfterFit_ = iConfig.getParameter<double>("JpsiMassWindowAfterFit");
  JpsiPtCut_ =  iConfig.getParameter<double>("JpsiPtCut");
  KaonTrackPtCut_ = iConfig.getParameter<double>("KaonTrackPtCut");//https://arxiv.org/pdf/1307.2782.pdf
  PionTrackPtCut_ = iConfig.getParameter<double>("PionTrackPtCut");
  BdKaonTrackPtCut_ = iConfig.getParameter<double>("BdKaonTrackPtCut");
  PhiMassWindowBeforeFit_ = iConfig.getParameter<double>("PhiMassWindowBeforeFit");
  PhiMassWindowAfterFit_ = iConfig.getParameter<double>("PhiMassWindowAfterFit");

  KstarMassWindowBeforeFit_ = iConfig.getParameter<double>("KstarMassWindowBeforeFit");
  KstarMassWindowAfterFit_ = iConfig.getParameter<double>("KstarMassWindowAfterFit");
  BdLowerMassCutBeforeFit_ = iConfig.getParameter<double>("BdLowerMassCutBeforeFit");
  BdUpperMassCutBeforeFit_ = iConfig.getParameter<double>("BdUpperMassCutBeforeFit");

  BdLowerMassCutAfterFit_ = iConfig.getParameter<double>("BdLowerMassCutAfterFit");
  BdUpperMassCutAfterFit_ = iConfig.getParameter<double>("BdUpperMassCutAfterFit");

  BdPDGMass_ = iConfig.getParameter<double>("BdPDGMass");
  BpPDGMass_ = iConfig.getParameter<double>("BpPDGMass");
  BsPDGMass_ = iConfig.getParameter<double>("BsPDGMass");
  
  outputFile_ = iConfig.getUntrackedParameter<std::string>("outputFile");
  event_counter_ = 0;
  elecounter_    = 0;
  muoncounter_   = 0;
  jetcounter_    = 0;
  tagmucounter_ =0;


  verbose_                = iConfig.getParameter<bool>("verbose");
  TestVerbose_            = iConfig.getParameter<bool>("TestVerbose");
  edm::LogInfo("RecoVertex/JpsiTrkTrk")<< "Initializing b-hadron with displaced track analyser  - Output file: " << outputFile_ <<"\n";
}

JpsiTrkTrk::~JpsiTrkTrk() {}
void JpsiTrkTrk::beginJob() {
  bRootTree_ = new JpsiTrkTrkRootTree();
  bRootTree_->createTree(outputFile_);
}
void JpsiTrkTrk::endJob() {
  bRootTree_->writeFile();
  delete bRootTree_;
  cout << "Total number of Events          : " << event_counter_ << endl;
  cout << "Total number of Tagged muons    : " << muoncounter_   << endl;
  cout << "Total number of Tagged electrons: " << elecounter_    << endl;
  cout << "Total number of Tagged jets     : " << jetcounter_    << endl;
  cout << "Max amount of Tag muons		    : " << tagmucounter_ << endl;
}


void JpsiTrkTrk::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  event_counter_++;
  bRootTree_->resetEntries();
 //===========================================================TLorentzVectors of B mesons after the kin fit
  TLorentzVector TheBp;
  TLorentzVector TheBs;
  TLorentzVector TheBd;
 //========================================================Create objects needed for building the B candidate
  pat::CompositeCandidate Bd0Cand_best;
  TrackRef trk1Ref_best;
  TrackRef trk2Ref_best;
  TrackRef trkMu1Ref_best;
  TrackRef trkMu2Ref_best;
  RefCountedKinematicParticle bd0_best;
  pat::Muon mu1_best;
  pat::Muon mu2_best;
 //=========================================================Define Primary vertices(PVs)
  Vertex PVvtxCosTheta;
  Vertex BpPVvtxCosTheta;
  Vertex BdPVvtxCosTheta;
  int Bd0PVVtxInd=0;
  vector<TrackRef> BpTrkRefs;
  vector<TrackRef> BsTrkRefs;
  vector<TrackRef> BdTrkRefs;
  if(isMCstudy_){
  //===================================================Creat handle for PU information
    
    edm:: Handle<edm::View<PileupSummaryInfo> > PUinfo;
    iEvent.getByToken( PUInfoTok, PUinfo);
    edm::View<PileupSummaryInfo>::const_iterator PVI;
    int numInteraction = 0;
    for(PVI = PUinfo->begin(); PVI != PUinfo->end(); ++PVI)
    {
      if (PVI->getBunchCrossing()==0)
        numInteraction += PVI->getPU_NumInteractions();
    }
    bRootTree_->PUinteraction_ = numInteraction; /// SaveToTree -> total number of interactions (with PU)
    //cout<<numInteraction<<"\n";
   }
   //==========================================================Get Primary Vertices
   //==================================================Handle and assign beam spot parameters
   edm::Handle<reco::BeamSpot> vertexBeamSpot ;
   iEvent.getByToken(vertexBeamSpotTok,vertexBeamSpot);

         BSx = vertexBeamSpot->x0(); 
         BSy = vertexBeamSpot->y0();
         BSz = vertexBeamSpot->z0(); 
         BSdxdz = vertexBeamSpot->dxdz();
         BSdydz = vertexBeamSpot->dydz(); 
         
	 edm::Handle<edm::View<reco::Vertex> > recVtxs;
         iEvent.getByToken(primaryvertexTok, recVtxs);

	 bRootTree_->NVertices_ = recVtxs->size();
	 for(size_t iVtx = 0; iVtx < recVtxs->size(); ++ iVtx)
 {
           VtxIndex = iVtx;
         //  cout<<iVtx<<"\n";
         // const Vertex &RecVtx = (*recVtxs)[iVtx];
         /*
         //cout<<vtx<<"\n";
           double PtSumVertex = 0.;
           for (reco::Vertex::trackRef_iterator trackvertex = vtx.tracks_begin(); trackvertex != vtx.tracks_end(); trackvertex++)
           {//  cout<<trackvertex<<"\n";
          cout<<"pt: "<<(**trackvertex).pt()<<"\n";
          const reco::Track & VtxTrack = *(trackvertex->get());
          //const pat::PackedCandidate & VtxTrack = *(trackvertex->get());
          PtSumVertex = PtSumVertex + abs(VtxTrack.pt());
          // cout<<"ptsumvertex: "<<PtSumVertex<<"\n";
           }
           if(PtSumVertex > MinPtVertex)
          {int VtxIndex = -99;
           VtxIndex = iVtx;
           MinPtVertex = PtSumVertex;
           cout<<"minvertex"<<MinPtVertex<<"\n";
           }*/
 }//-----------------------------------------------------VERTEX LOOP ENDS
        const Vertex &RecVtx = (*recVtxs)[VtxIndex];
	if(VtxIndex!=-99) /// Use PV
    {
      bRootTree_->isPV_ = 1;
      PVx = RecVtx.x();
      PVy= RecVtx.y();
      PVz= RecVtx.z();
      PVerrx=RecVtx.xError();
      PVerry=RecVtx.yError();
      PVerrz=RecVtx.zError();
    }
  else {  /// Use BS
            bRootTree_->isBS_ = 1;
            PVx=BSx;
            PVy=BSy;
            PVz=BSz;
            PVerrx=BSdx;
            PVerry=BSdy;
            PVerrz=BSdz;
            std:: cout<<PVerrz<<"\n";
       }
             bRootTree_->getVtx(BSx,BSy,BSz,PVx,PVy,PVz,PVerrx,PVerry,PVerrz);
             bRootTree_->BSdx_ = BSdx;
             bRootTree_->BSdy_ = BSdy;
             bRootTree_->BSdz_ = BSdz;
             bRootTree_->BSdxdz_ = vertexBeamSpot->dxdz();
             bRootTree_->BSdydz_ = vertexBeamSpot->dydz();

            if(verbose_ == true){
                         std::cout<<"BeamSpot   (x,y,z) = ("<< BSx << ", " << BSy << ", "<< BSz << ")\n";
                         std::cout<<"PrimaryVtx (x,y,z) = ("<< PVx <<" , " << PVy << ", "<< PVz << ")\n";
                     }

        //========================================================================================================Run/Event/Lumiblock
         bRootTree_->runNumber_ = iEvent.id().run();
        //cout<<iEvent.id().run()<<"\n";
        bRootTree_->eventNumber_ = (unsigned int)iEvent.id().event();
        bRootTree_->lumiSection_ = iEvent.luminosityBlock();
        
	edm::Handle<edm::TriggerResults> hltresults;
        iEvent.getByToken(triggerresultsTok, hltresults);
	//std::cout<<"hlt size:  "<<hltresults->size()<<"\n";
  

	const  edm::TriggerNames & triggerNames_ = iEvent.triggerNames(*hltresults);
        int ntrigs = hltresults->size();
	for (int itrig = 0; itrig != ntrigs; ++itrig)
        {
                                                  TString trigName = triggerNames_.triggerName(itrig);
                                                  //std::cout<<"triggernames:"<<itrig<<"::"<<trigName<<"\n";

                if (trigName=="triggerbit_HLTDimuon4JpsiDisplaced_")      bRootTree_->triggerbit_HLTDimuon4JpsiDisplaced_       = hltresults->accept(itrig);
                if (trigName=="triggerbit_HLTDimuon4JpsiNoVertexing_")    bRootTree_->triggerbit_HLTDimuon4JpsiNoVertexing_     = hltresults->accept(itrig);
                if (trigName=="triggerbit_HLTDimuon4JpsiTrkTrkDisplaced_")bRootTree_->triggerbit_HLTDimuon4JpsiTrkTrkDisplaced_ = hltresults->accept(itrig);
		if (trigName=="triggerbit_HLTDimuon4JpsiTrkBc_")          bRootTree_->triggerbit_HLTDimuon4JpsiTrkBc_ = hltresults->accept(itrig);

                string str = (string) trigName  ;

                   /*  if (str.compare(0,18,"HLT_DoubleMu4_3_Jpsi_Displaced_v12") == 0)
                     {
                        cout << "trigger found! " << str << endl;
                        cout <<"dimuon4_jpsidis  " << bsRootTree_->triggerbit_HLTmu4TkDis_  << endl;
                        cout <<"dimuon4_jpsitktkdis  " << bsRootTree_->triggerbit_HLTmu4TkTkDis_ << endl;
                     }*/

         }//Trigger loop ends
 
	 edm::Handle< View<pat::Muon> > allmuons;
         iEvent.getByToken(MuonTagTok, allmuons);
         //if(allmuons->size()>0)std::cout<<"muonmultiplicity"<<allmuons->size()<<"\n";
         bRootTree_->MuonMultiplicity_ = allmuons->size();
         //==================================================================================================================Muon LOOP
//-------------------------------------------------------First Muon
      for(size_t i=0; i < allmuons->size(); ++i)
          {
                    const pat::Muon & mu1 = (*allmuons)[i];
                    if (mu1.innerTrack().isNull()){continue;}
                    //cout<<"mu1pt: "<<mu1.pt()<<"\n";
                    if(verbose_ == true)
                        {
                                std::cout<<"Got one muon "<<mu1.pt()<<std::endl;
                        }

//------------------------------------------------------Loop over 2nd muon
                   for (size_t j=i+1; j < allmuons->size(); ++j)
                        {
                              const pat::Muon & mu2 = (*allmuons)[j];
                              if (mu2.innerTrack().isNull()){continue;}	   
                              //  cout<<"mu2pt: "<<mu2.pt()<<"\n";
                              if(verbose_ == true)
                                  {
                                         std::cout<<"Got the second muon "<<mu2.pt()<<std::endl;
                                   }
               
           	   bRootTree_->Mu1SoftID_ = mu1.passed(reco::Muon::SoftCutBasedId);
                   bRootTree_->Mu2SoftID_ = mu2.passed(reco::Muon::SoftCutBasedId);
                   if(!mu1.isGlobalMuon() && !mu1.isTrackerMuon()) continue; // skip if mu1 is not GLB or TRK (PASS IF LooseMuId)
                   if(!mu2.isGlobalMuon() && !mu2.isTrackerMuon()) continue; // skip if mu2 is not GLB or TRK (PASS IF LooseMuId)
                   if(verbose_==true) {
                                          std::cout << "******mu1.isGlobalMuon() == "<<mu1.isGlobalMuon() << "\n";
                                          std::cout << "      mu1.isTrackerMuon()== "<<mu1.isTrackerMuon()<< "\n";
                                          std::cout << "      mu2.isGlobalMuon() == "<<mu2.isGlobalMuon() << "\n";
                                          std::cout << "      mu2.isTrackerMuon()== "<<mu2.isTrackerMuon()<< "\n";
                                       }
		   bRootTree_->ihaveajpsi_=1;
		    if(mu1.charge()==mu2.charge()) continue; // Skip iif mu1 e mu2 have the same charge
                          if(verbose_==true) {
                                                  std::cout<<"******MUONS HAVE OPPOSITE CHARGE: mu1.charge() = " <<mu1.charge()<<" , mu2.charge() = "<<mu2.charge()<<std::endl;
                                             }

                                                if(bRootTree_->iPassedCutIdent_   < 1 )   bRootTree_->iPassedCutIdent_ = 1 ;
                                                if(bRootTree_->iPassedCutIdentBd_   < 1 ) bRootTree_->iPassedCutIdentBd_ = 1 ;

                    pat::CompositeCandidate Jpsi;
                    Jpsi.addDaughter(mu1);
                    Jpsi.addDaughter(mu2);
                    AddFourMomenta addP4;
                    addP4.set(Jpsi);
		    if(verbose_==true)
                                          {
                                                 std::cout<<"******Di-Muon Mass= " <<Jpsi.mass()<<std::endl;
                                           }
                    if ( abs(Jpsi.mass() - nominalJpsiMass ) > JpsiMassWindowBeforeFit_ ) continue; // skip if mu1-mu2 combination mass is far from JPsi
                    if ( Jpsi.pt() < JpsiPtCut_) continue;                                          // skip if mu1-mu2 combination pt is less than JPsi Pt cut

                    //Jpsi window validity

                   if(bRootTree_->iPassedCutIdent_   < 2 )   bRootTree_->iPassedCutIdent_   = 2 ;
                   if(bRootTree_->iPassedCutIdentBd_   < 2 ) bRootTree_->iPassedCutIdentBd_ = 2 ;

                  //edm::ESHandle<TransientTrackBuilder> theB;
                  //iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
		  const auto& trackBuilder = iSetup.getData(trackBuilderToken);
                  TrackRef trkMu1Ref = mu1.get<TrackRef>();
                  TrackRef trkMu2Ref = mu2.get<TrackRef>();
                  TrackRef muonTrkP = mu1.track();
                  TrackRef muonTrkM = mu2.track();
                  vector<TransientTrack> trk_all;
                  TransientTrack mu1TT=trackBuilder.build(*trkMu1Ref);
                  TransientTrack mu2TT=trackBuilder.build(*trkMu2Ref);
                  trk_all.push_back(mu1TT);
                  trk_all.push_back(mu2TT);
                  KalmanVertexFitter kvf(true);
                  TransientVertex tv = kvf.vertex(trk_all);

    
                  if(verbose_==true) 
                  {
                              std::cout<<"****** MUONS HAVE VALID VERTEX FIT"<< std::endl;
                  }



                 // vertex validity
                 if(bRootTree_->iPassedCutIdent_   < 3 )   bRootTree_->iPassedCutIdent_   = 3 ;     
                 if(bRootTree_->iPassedCutIdentBd_   < 3 ) bRootTree_->iPassedCutIdentBd_ = 3 ;
		 
                 if (mu1.charge()==1) {
                 float mupPhi = atan(mu1.py()/mu1.px());
                 if ( mu1.px() < 0 && mu1.py() < 0 ) mupPhi -= TMath::Pi();
                 if ( mu1.px() < 0 && mu1.py() > 0 ) mupPhi += TMath::Pi();
                 float mumPhi = atan(mu2.py()/mu2.px());
                 if ( mu2.px() < 0 && mu2.py() < 0 ) mumPhi -= TMath::Pi();
                 if ( mu2.px() < 0 && mu2.py() > 0 ) mumPhi += TMath::Pi();
                 //std::cout<<"mumphi"<<mumPhi<<"\n";
                 if ( (mupPhi - mumPhi)>0 ) isCowboy=1; }
                 else {
                   float mupPhi = atan(mu2.py()/mu2.px());
                   if ( mu2.px() < 0 && mu2.py() < 0 ) mupPhi -= TMath::Pi();
                   if ( mu2.px() < 0 && mu2.py() > 0 ) mupPhi += TMath::Pi();
                   float mumPhi = atan(mu1.py()/mu1.px());
                   if ( mu1.px() < 0 && mu1.py() < 0 ) mumPhi -= TMath::Pi();
                   if ( mu1.px() < 0 && mu1.py() > 0 ) mumPhi += TMath::Pi();
                   if ( (mupPhi - mumPhi)>0 ) isCowboy=1;
                  }
                 //==========================================================================ceate vertex for mu1-mu2 combination
                 Vertex vertex = tv;
                 //calculate variable in the cloasest way to trigger
                 double  vtxProb_Jpsi = TMath::Prob(vertex.chi2(),(int)vertex.ndof());
                 math::XYZVector      pperp(mu1.px() + mu2.px(), mu1.py() + mu2.py(), 0.);
                 reco::Vertex::Point  vpoint=vertex.position();
                 //std::cout<<"vpoint"<<vpoint<<"\n";

		 GlobalPoint secondaryVertex (vpoint.x(), vpoint.y(), vpoint.z());
                 GlobalPoint displacementFromBeamspot( -1*((BSx -  secondaryVertex.x()) +  (secondaryVertex.z() - BSz) * BSdxdz),-1*((BSy - secondaryVertex.y())+  (secondaryVertex.z() - BSz) * BSdydz), 0);
                reco::Vertex::Point vperp(displacementFromBeamspot.x(),displacementFromBeamspot.y(),0.);
                double CosAlpha = vperp.Dot(pperp)/(vperp.R()*pperp.R());
                //std::cout<<"opening angle: "<<CosAlpha<<"\n";
                 TrajectoryStateClosestToPoint mu1TS = mu1TT.impactPointTSCP();
                 TrajectoryStateClosestToPoint mu2TS = mu2TT.impactPointTSCP();
                 if (mu1TS.isValid() && mu2TS.isValid()) {

                 ClosestApproachInRPhi cApp;
                 cApp.calculate(mu1TS.theState(), mu2TS.theState());
                 MuonsDCA=cApp.distance();
               }
                 double max_Dr1=fabs( (- (mu1.vx()-BSx) * mu1.py() + (mu1.vy()-BSy) * mu1.px() ) / mu1.pt() );
                 double max_Dr2=fabs( (- (mu2.vx()-BSx) * mu2.py() + (mu2.vy()-BSy) * mu2.px() ) / mu2.pt() );
                 // std::cout<<"maxDR1: "<<max_Dr1<<"\n";
//========================================================================================================================muon overlaping remover
                if ( muon::overlap(mu1,mu2,1,1,true) ) continue; /// Skip the mu-mu combination if the two muons overlap
                reco::Vertex::Error verr = vertex.error();
                //std::cout<<"Verror: "<<verr<<"\n";
//==================================================================================================================Translate to global error
                GlobalError err(verr.At(0,0), verr.At(1,0), verr.At(1,1), verr.At(2,0), verr.At(2,1), verr.At(2,2) );
                float lxy = displacementFromBeamspot.perp();
                float lxyerr = sqrt(err.rerr(displacementFromBeamspot));
                //std::cout<<"LXYerror: "<<lxyerr<<"\t"<<"LXY: "<<lxy<<"\n";
                //Requires a filter match : bRootTree_->JpsiNumberOfCandidates_++;

		bRootTree_->JpsiM_alone_ = Jpsi.mass();
		bRootTree_->JpsiPhi_alone_ = Jpsi.phi();
                bRootTree_->JpsiEta_alone_ = Jpsi.eta();
                bRootTree_->JpsiPt_alone_ = Jpsi.pt();
                if( mu1.charge() == -1 )
       { 

                      bRootTree_->Mu1Pt_beffit_   = mu1.pt();
                      bRootTree_->Mu1Pz_beffit_   = mu1.pz();
                      bRootTree_->Mu1Eta_beffit_  = mu1.eta();
                      bRootTree_->Mu1Phi_beffit_  = mu1.phi();
                     // std::cout<<"mu1phi:"<<mu1.phi();
                      bRootTree_->Mu2Phi_beffit_  = mu2.phi();
                      bRootTree_->Mu2Pt_beffit_  = mu2.pt();
                      bRootTree_->Mu2Pz_beffit_   = mu2.pz();
                      bRootTree_->Mu2Eta_beffit_  = mu2.eta();
       }

      else 
      {
                            bRootTree_->Mu2Pt_beffit_   = mu1.pt();
                            bRootTree_->Mu2Pz_beffit_   = mu1.pz();
                            bRootTree_->Mu2Eta_beffit_  = mu1.eta();
                            bRootTree_->Mu2Phi_beffit_  = mu1.phi();
                            //std::cout<<"mu1phi"<<mu1.phi();
                            bRootTree_->Mu1Phi_beffit_  = mu2.phi();
                            bRootTree_->Mu1Pt_beffit_  = mu2.pt();
                            bRootTree_->Mu1Pz_beffit_   = mu2.pz();
                            bRootTree_->Mu1Eta_beffit_  = mu2.eta();
                     
      }

                bRootTree_->JpsiMuMuDCA_beffit_  = MuonsDCA;
 		bRootTree_->JpsiCosDeltaAlpha_  = CosAlpha;
                // std::cout<<"cosalpha"<<CosAlpha<<"\n";
 		bRootTree_->JpsiLxySigma_   = lxyerr;
		bRootTree_->JpsiLxy_  = lxy;
		bRootTree_->JpsiLxyOverPt_ =	( displacementFromBeamspot.x() * Jpsi.px() + 	displacementFromBeamspot.y() * Jpsi.py() ) /	( Jpsi.pt()* Jpsi.pt()) ;
                //std::cout<<"jpsilxyoverpt: "<<( displacementFromBeamspot.x() * Jpsi.px() +    displacementFromBeamspot.y() * Jpsi.py() ) /    ( Jpsi.pt()* Jpsi.pt())<<"\n";
                if( tv.isValid() )
                {
	          	 bRootTree_->JpsiTrigVtxProb_  =  vtxProb_Jpsi;
                        // std::cout<<"prob:"<<vtxProb_Jpsi;
		}
   
		
		if (mu1.isTrackerMuon() && !mu1.isGlobalMuon())bRootTree_->JpsiMuon1Cat_alone_ = 1;
                else if (!mu1.isTrackerMuon() && mu1.isGlobalMuon())bRootTree_->JpsiMuon1Cat_alone_ = 2;
                else if (mu1.isTrackerMuon() && mu1.isGlobalMuon())bRootTree_->JpsiMuon1Cat_alone_ = 3;
                else if (!mu1.isTrackerMuon() && !mu1.isGlobalMuon())bRootTree_->JpsiMuon1Cat_alone_ = 4;
                if (mu1.isPFMuon())       bRootTree_->JpsiMuonCat1_ = 1;

                if (mu2.isTrackerMuon() && !mu2.isGlobalMuon())       bRootTree_->JpsiMuon2Cat_alone_ = 1;
                else if (!mu2.isTrackerMuon() && mu2.isGlobalMuon())  bRootTree_->JpsiMuon2Cat_alone_ = 2;
                else if (mu2.isTrackerMuon() && mu2.isGlobalMuon())   bRootTree_->JpsiMuon2Cat_alone_ = 3;
                else if (!mu2.isTrackerMuon() && !mu2.isGlobalMuon()) bRootTree_->JpsiMuon2Cat_alone_ = 4;
                if (mu2.isPFMuon())       bRootTree_->JpsiMuonCat2_ = 1;

		int pixhits1 = 0;
                const reco::HitPattern& pp1 = trkMu1Ref.get()->hitPattern();
                for (int iter=0; iter<pp1.numberOfAllHits(reco::HitPattern::TRACK_HITS); iter++) {
			uint32_t hit = pp1.getHitPattern(reco::HitPattern::TRACK_HITS,iter);
			if (pp1.validHitFilter(hit) && pp1.pixelBarrelHitFilter(hit)) pixhits1++;
			if (pp1.validHitFilter(hit) && pp1.pixelEndcapHitFilter(hit)) pixhits1++;
                  }
                bRootTree_->JpsiMu1nPixHits_alone_   = pixhits1;

                int pixhits2 = 0;
                const reco::HitPattern& pp2 = trkMu2Ref.get()->hitPattern();
                for (int iter=0; iter<pp2.numberOfAllHits(reco::HitPattern::TRACK_HITS); iter++)
                 {
			uint32_t hit = pp2.getHitPattern(reco::HitPattern::TRACK_HITS,iter);
                        //std::cout<<"hits"<<hit<<"\n";
			if (pp2.validHitFilter(hit) && pp2.pixelBarrelHitFilter(hit)) pixhits2++;
			if (pp2.validHitFilter(hit) && pp2.pixelEndcapHitFilter(hit)) pixhits2++;
                 }

                bRootTree_->JpsiMu2nPixHits_alone_   = pixhits2;

		edm::Handle<View<pat::PackedCandidate>> allTracks;
                iEvent.getByToken(trackLabelK, allTracks);
		 for (size_t k=0; k< allTracks->size(); ++k)
                  {
                           const pat::PackedCandidate & track1 = (*allTracks)[k];
                             if (!track1.hasTrackDetails())continue;

                             if (track1.charge()<0)continue;
                             if (track1.pt() < KaonTrackPtCut_) continue;
                             if (track1.numberOfHits() < 5)continue;
                             if(!track1.trackHighPurity()) continue;
			     DeltaRKaonJpsi = deltaR(Jpsi.eta(), Jpsi.phi(), track1.eta(), track1.phi());
                             if (DeltaRKaonJpsi > 2.2) continue;
                             const reco::Track &  rtrk1 = (*allTracks)[k].pseudoTrack();
                             if (rtrk1.charge()<0) continue;
                             TransientTrack KPTT = trackBuilder.build(&rtrk1);
                             TrajectoryStateClosestToPoint KPTS = KPTT.impactPointTSCP();
                             if(!KPTS.isValid())continue;
                             if (!track1.clone()->hasTrackDetails())continue;
                             pat::PackedCandidate *track11 = track1.clone();


			     for (size_t l=k+1; l< allTracks->size(); ++l)
                           { 
		             const pat::PackedCandidate & track2 = (*allTracks)[l];
	                     if ( !track2.hasTrackDetails() )continue;
                             if (track2.charge()>0) continue;
                             if (track2.pt() < PionTrackPtCut_) continue;
                             if ( track2.numberOfHits()<5) continue;
                             if(!track2.trackHighPurity()) continue;
			     DeltaRPionJpsi = deltaR(Jpsi.eta(), Jpsi.phi(), track2.eta(), track2.phi());
			     if (DeltaRPionJpsi >2.2) continue;
                             const reco::Track &  rtrk2 = (*allTracks)[l].pseudoTrack();
                             if (rtrk2.charge()>0) continue;
                             TransientTrack PiMTT = trackBuilder.build(&rtrk2);
                             TrajectoryStateClosestToPoint PiMTS = PiMTT.impactPointTSCP();
                             if(!PiMTS.isValid())continue;
			     if (KPTS.isValid() && PiMTS.isValid()) {
				     ClosestApproachInRPhi cAppK;
                                     cAppK.calculate(KPTS.theState(), PiMTS.theState());
                                     KPiDCA=cAppK.distance();
			     }
			     if(KPiDCA > 0.5)continue;
			     if (!track2.clone()->hasTrackDetails())continue;
                             pat::PackedCandidate *track22 = track2.clone();
                             //bRootTree_->KaonsDCA_   = KaonsDCA;
                             if(bRootTree_->iPassedCutIdent_ < 4 ) bRootTree_->iPassedCutIdent_ = 4 ;
                             pat::CompositeCandidate KStarCand;
			     track11->setMass(kaonmass);
                             KStarCand.addDaughter(*track11);
			     track22->setMass(pionmass);
                             KStarCand.addDaughter(*track22);
			     AddFourMomenta ad;
                             ad.set(KStarCand);
			     if (abs(KStarCand.mass()- nominalKstarMass) > KstarMassWindowBeforeFit_) continue;
                             if(bRootTree_->iPassedCutIdent_   < 5 ) bRootTree_->iPassedCutIdent_ = 5 ;
                             pat::CompositeCandidate Bd0Cand;
       	  		       Bd0Cand.addDaughter(mu1);
       	  		       Bd0Cand.addDaughter(mu2);
       	  	               Bd0Cand.addDaughter(*track11);
       	  		       Bd0Cand.addDaughter(*track22);
       	  		       AddFourMomenta add4mom;
       	  		       add4mom.set(Bd0Cand);
                               //std::cout<<"mass Bd0 : "<<BCand.mass()<<"\n";
       	  		       if (Bd0Cand.mass() < BdLowerMassCutBeforeFit_ || Bd0Cand.mass() > BdUpperMassCutBeforeFit_) continue;
                               if(bRootTree_->iPassedCutIdent_   < 8 ) bRootTree_->iPassedCutIdent_ = 8 ;


			      vector<TransientTrack> t_tracks;
       	  	              t_tracks.push_back(trackBuilder.build(*trkMu1Ref));
       	  		      t_tracks.push_back(trackBuilder.build(*trkMu2Ref));
       	  		      t_tracks.push_back(trackBuilder.build(&rtrk1));
       	  		      t_tracks.push_back(trackBuilder.build(&rtrk2));
			      if (!trkMu1Ref.isNonnull() || !trkMu2Ref.isNonnull() )continue;
			      if(bRootTree_->iPassedCutIdent_   < 9 ) bRootTree_->iPassedCutIdent_ = 9 ;
			       vector<TransientTrack> kstar_tracks;
                               kstar_tracks.push_back(trackBuilder.build(&rtrk1));//trk1Ref
                               kstar_tracks.push_back(trackBuilder.build(&rtrk2));//trk2Ref
                               KalmanVertexFitter kvfkstar;
                               TransientVertex tvkstar = kvfkstar.vertex(kstar_tracks);
                               if (!tvkstar.isValid()) continue;
			       Vertex vertexkstar = tvkstar;
                               double vtxProb_Kstar = TMath::Prob(vertexkstar.chi2(),(int)vertexkstar.ndof());

			       KalmanVertexFitter kvfbs;
                               TransientVertex kvfbd0vertex = kvfbs.vertex(t_tracks);
                               Vertex vertexbd0kalman = kvfbd0vertex;
                               if (!kvfbd0vertex.isValid()) continue;
                               GlobalError gigi=kvfbd0vertex.positionError();
			       bRootTree_->K1Pt_beffit_   = track1.pt();
	 	 	       bRootTree_->K1Pz_beffit_   = track1.pz();
	  	         	bRootTree_->K1Eta_beffit_  = track1.eta();
	  		        bRootTree_->K1Phi_beffit_  = track1.phi();
	 	 		bRootTree_->PiPt_beffit_   = track2.pt();
	  			bRootTree_->PiPz_beffit_   = track2.pz();
	  			bRootTree_->PiEta_beffit_  = track2.eta();
	  			bRootTree_->PiPhi_beffit_  = track2.phi();
				KinematicFitInterface Kfitter;
				bool fitSuccess = Kfitter.doFit(t_tracks, nominalMuonMass,  nominalKaonMass, nominalPionMass);//Actual Fit for BD candidate vertex
                                if(fitSuccess != 1) continue;
				if(bRootTree_->iPassedCutIdent_   < 10 ) bRootTree_->iPassedCutIdent_ = 10 ;
                                double vtxprob_Bd0 = TMath::Prob(vertexbd0kalman.chi2(),(int)vertexbd0kalman.ndof());
				RefCountedKinematicParticle bd0 = Kfitter.getParticle();
	  		        RefCountedKinematicVertex bVertex = Kfitter.getVertex();
	  		        AlgebraicVector7 b_par = bd0->currentState().kinematicParameters().vector();
	  		        //AlgebraicSymMatrix77 bd0_er = bd0->currentState().kinematicParametersError().matrix();
                                AlgebraicMatrix33 BVError(bVertex->error().matrix());
	  		        double fittedBd0Mass = b_par[6];
				if(!bVertex->vertexIsValid()) continue;
				TMatrix cova(2,2);
                                cova.IsSymmetric();
                                cova(0,0)=gigi.cxx();
                                cova(1,1)=gigi.cyy();
                                cova(0,1)=gigi.cyx();
                                cova(1,0)=gigi.cyx();
				if(abs(Jpsi.mass() - nominalJpsiMass) < JpsiMassWindowAfterFit_  
				&& Jpsi.pt() > JpsiPtCut_   
				&& abs(KStarCand.mass()- nominalKstarMass) > KstarMassWindowBeforeFit_    
				&&  fittedBd0Mass > BdLowerMassCutAfterFit_  
				&&  fittedBd0Mass < BdUpperMassCutAfterFit_ ){
                                bRootTree_->Bd0NumberOfCandidatesAfterFit_++;
				}
				if(vtxprob_Bd0 > minVtxP)
                               {
				       reco::Track trk1Ref_best = rtrk1;//trk1Ref;
                              reco::Track trk2Ref_best = rtrk2;//trk2Ref;
                              trkMu1Ref_best=trkMu1Ref;
                              trkMu2Ref_best=trkMu2Ref;
                              mu1_best=mu1;
                              mu2_best=mu2;



                               TrackRef mu1trkref = mu1.get<TrackRef>();
                               TrackRef mu2trkref = mu2.get<TrackRef>();
                               TrackRef K1trkRef = track1.get<TrackRef>();
                               TrackRef PitrkRef = track2.get<TrackRef>();
			       BdTrkRefs.clear();
                               BdTrkRefs.push_back(mu1trkref);
                               BdTrkRefs.push_back(mu2trkref);
                               BdTrkRefs.push_back(K1trkRef);
                               BdTrkRefs.push_back(PitrkRef);

		               if (abs(Jpsi.mass() - nominalJpsiMass) > JpsiMassWindowAfterFit_ || Jpsi.pt() < JpsiPtCut_) continue;
                               if (abs(KStarCand.mass()- nominalKstarMass) > KstarMassWindowBeforeFit_)continue;
		               if (fittedBd0Mass < BdLowerMassCutAfterFit_ || fittedBd0Mass > BdUpperMassCutAfterFit_) continue;		
			       minVtxP = vtxprob_Bd0;
                               Bd0Cand_best = Bd0Cand;
			       GlobalVector Bd0vec(b_par[3], b_par[4], b_par[5]); // the fitted momentum vector of the Bs
			       bRootTree_->Bd0FitM_ = fittedBd0Mass;
                               bRootTree_->Bd0FitEta_ = Bd0vec.eta();
                               bRootTree_->Bd0FitPt_  = Bd0vec.perp();
                               bRootTree_->Bd0FitPz_  = Bd0vec.z();
                               bRootTree_->Bd0FitPhi_ = Bd0vec.phi();
			       }
			   }
		  }

			}}	      
              bRootTree_->fill();

}



// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void JpsiTrkTrk::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //edm::ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks", edm::InputTag("ctfWithMaterialTracks"));
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(JpsiTrkTrk);
