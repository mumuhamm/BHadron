// -*- C++ -*-
// Package:    BHadron/JpsiTrkTrk
// Class:      JpsiTrkTrk
//class JpsiTrkTrk JpsiTrkTrk.cc BHadron/JpsiTrkTrk/plugins/JpsiTrkTrk.cc

// Original Author:  Alibordi Muhammad
//         Created:  Sat, 13 Apr 2024 13:11:49 GMT



#include <memory>
#include "BHadron/JpsiTrkTrk/interface/JpsiTrkTrkRootTree.h"
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
//#include "PhysicsTools/CandUtils/interface/AddFourMomenta.h"
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

//=============================================================================================
  StoreDeDxInfo_ = iConfig.getParameter<bool>("StoreDeDxInfo");
  JpsiMassWindowBeforeFit_ = iConfig.getParameter<double>("JpsiMassWindowBeforeFit");

  BsLowerMassCutBeforeFit_  = iConfig.getParameter<double>("BsLowerMassCutBeforeFit");
  BsUpperMassCutBeforeFit_  = iConfig.getParameter<double>("BsUpperMassCutBeforeFit");
  BsLowerMassCutAfterFit_  = iConfig.getParameter<double>("BsLowerMassCutAfterFit");
  BsUpperMassCutAfterFit_  = iConfig.getParameter<double>("BsUpperMassCutAfterFit");

  JpsiMassWindowAfterFit_ = iConfig.getParameter<double>("JpsiMassWindowAfterFit");
  JpsiPtCut_ =  iConfig.getParameter<double>("JpsiPtCut");
  KaonTrackPtCut_ = iConfig.getParameter<double>("KaonTrackPtCut");
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
  pat::CompositeCandidate BCand_best;
  TrackRef trk1Ref_best;
  TrackRef trk2Ref_best;
  TrackRef trkMu1Ref_best;
  TrackRef trkMu2Ref_best;
  RefCountedKinematicParticle bs_best;
  pat::Muon mu1_best;
  pat::Muon mu2_best;
 //=========================================================Define Primary vertices(PVs)
  Vertex PVvtxCosTheta;
  Vertex BpPVvtxCosTheta;
  Vertex BdPVvtxCosTheta;
  int BsPVVtxInd=0;
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
 int    VtxIndex    = -99;
//====================================================Beam Spot
  double BSx         = -9999999.;
  double BSy         = -9999999.;
  double BSz         = -9999999.;
  double BSdx        = -9999999.;
  double BSdy        = -9999999.;
  double BSdz        = -9999999.;
  double BSdxdz      = -9999999.;
  double BSdydz      = -9999999.;
  double BSsigmaZ    = -9999999.;
  double BSdsigmaZ   = -9999999.;
//========================================================PV
  double PVx         = -9999999.;
  double PVy         = -9999999.;
  double PVz         = -9999999.;
  double PVerrx      = -9999999.;
  double PVerry      = -9999999.;
  double PVerrz      = -9999999.;
   TLorentzVector kaontrack1, kaontrack2;
//==================================================Handle and assign beam spot parameters
   edm::Handle<reco::BeamSpot> vertexBeamSpot ;
   iEvent.getByToken(vertexBeamSpotTok,vertexBeamSpot);

         BSx = vertexBeamSpot->x0(); 
         BSy = vertexBeamSpot->y0();
         BSz = vertexBeamSpot->z0(); 
         //BSdxdz = vertexBeamSpot->dxdz();
         //BSdydz = vertexBeamSpot->dydz(); 
         
	 edm::Handle<edm::View<reco::Vertex> > recVtxs;
         iEvent.getByToken(primaryvertexTok, recVtxs);

	 //bRootTree_->NVertices_ = recVtxs->size();
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
