#ifndef BHadron_JpsiTrkTrk_JpsiTrkTrk_h
#define BHadron_JpsiTrkTrk_JpsiTrkTrk_h

#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "MuonAnalysis/MuonAssociators/interface/PropagateToMuon.h"
#include "MuonAnalysis/MuonAssociators/interface/L1MuonMatcherAlgo.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicVertex.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelTopologyBuilder.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"
#include <TFile.h>
#include <TH1F.h>
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidateFwd.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include "DataFormats/PatCandidates/interface/IsolatedTrack.h"
#include "TLorentzVector.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "BHadron/JpsiTrkTrk/interface/JpsiTrkTrkRootTree.h"


class JpsiTrkTrk : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit JpsiTrkTrk(const edm::ParameterSet&);
  ~JpsiTrkTrk() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;
   void fillMCInfo( edm::Handle<edm::View<reco::GenParticle>>& genParticles);
  void setFitParKK(RefCountedKinematicTree& myTree);
  void setFitParHyp1(RefCountedKinematicTree& myTree);
  void setFitParHyp2(RefCountedKinematicTree& myTree);
private:
  bool MCmatching(const reco::Candidate & track1,  edm::Handle<edm::View<reco::GenParticle>>& genParticles,
                  int &K1mcId, int &K1momId, int &K1gmomId,
                  int condMom, int condGMom);
  bool MCmatchingBplusK(const reco::Candidate & track1, edm::Handle<edm::View<reco::GenParticle>>& genParticles,
                        int &K1mcId, int &K1momId,
                        int condMom);
  bool MCmatchingJpsi(const reco::Candidate& track1,  edm::Handle<edm::View<reco::GenParticle>>& genParticles, int &K1mcId, int &K1momId,	int condMom);

     reco::Vertex reVertex(const edm::Event& theEvent,  const edm::EventSetup& iSetup, edm::Handle<reco::BeamSpot>& vertexBeamSpot, reco::Vertex RecVtx, pat::Muon mu1, pat::Muon mu2, reco::Track & rtrk1, reco::Track & rtrk2);
     //reco::Vertex reVertex(const edm::EventSetup& ,edm::Handle<reco::BeamSpot> , reco::Vertex , pat::Muon , pat::Muon , reco::TrackRef , reco::TrackRef );
     //reco::Vertex reVertex2(const edm::EventSetup& ,edm::Handle<reco::BeamSpot>, reco::Vertex , pat::Muon , pat::Muon , reco::TrackRef, double &vtxCL);
     reco::Vertex reVertex2(const edm::EventSetup& iSetup, edm::Handle<reco::BeamSpot>& vertexBeamSpot, reco::Vertex RecVtx, pat::Muon mu1, pat::Muon mu2, reco::TrackRef trk3, double &vtxCL);

     double CalculateCtErrvertex(const edm::EventSetup& iSetup, reco::Vertex PV, TransientVertex SV, TMatrix SVcovariance2D ,GlobalVector Ptvec, double Bmass);
     double CalculateCtErrbeamspot(const edm::EventSetup& iSetup, edm::Handle<reco::BeamSpot>& PV, TransientVertex SV, TMatrix SVcovariance2D ,GlobalVector Ptvec, double Bmass); // for Beam spot
 
     // double CalculateCtErr(const edm::EventSetup& iSetup, edm::Handle<reco::BeamSpot>& PV, TransientVertex SV, TMatrix SVcovariance2D ,GlobalVector Ptvec, double Bmass);
     // double CalculateCtErr(const edm::EventSetup& iSetup, reco::Vertex PV, TransientVertex SV, TMatrix SVcovariance2D ,GlobalVector Ptvec, double Bmass); // for Beam spot 
  void RecursivelyPrintMother( const reco::Candidate & genp );

  GlobalVector flightDirection(const reco::Vertex &pv, reco::Vertex &sv);


  edm::ParameterSet theConfig_;

  bool selGlobalMuon(const pat::Muon aMuon, const math::XYZPoint RefVtx);
  bool selTrackerMuon(const pat::Muon aMuon, const math::XYZPoint RefVtx);
  bool selTightMuon(const pat::Muon aMuon, const reco::Vertex RecVtx);

  /// MuonChargeCone - TrackBased
  double MuonChargeCone(const edm::Event& theEvent, reco::TrackRef muTrackRef, const double Dr, const double KExp, bool IncludingMuon);
  double MuonChargeConeWrtPV(const edm::Event& theEvent, reco::TrackRef muTrackRef, reco::Vertex PVtx, const double Dr, const double KExp, bool IncludingMuon);
  
  /// LeptonChargeCone - PFcandidateBased
  double LeptonChargeCone(const reco::PFCandidateCollection & PFCand, const reco::PFCandidate theLept, const double Dr, const double KExp, bool IncludingLepton);
   /// MuonMC
  short int FindMuonMCCode(const pat::Muon theMu, edm::Handle<edm::View<reco::GenParticle>>& genParticles);
  short int FindMuonMCSimpleCode(const pat::Muon theMu, edm::Handle<edm::View<reco::GenParticle>>& genParticles);
  int FindMuonAncestor(const pat::Muon theMu, edm::Handle<edm::View<reco::GenParticle>>& genParticles);
   /// GenericMC
  int LookForMotherStringId(reco::GenParticle theGenP);
  short int LookForMotherString(reco::GenParticle theGenP);
  const TrackerGeometry* m_tracker;
  bool isMCstudy_;
  std::string outputFile_; // output file
  JpsiTrkTrkRootTree *bRootTree_;

  edm::InputTag genParticlesLabel;
  edm::EDGetTokenT<edm::View<reco::GenParticle>> genParticlesTok;
  edm::InputTag MuonTag;
  edm::EDGetTokenT<edm::View<pat::Muon>> MuonTagTok;
  edm::InputTag JetCollection;
  edm::EDGetTokenT<edm::View<pat::Jet>> JetCollectionTok;
  edm::InputTag PUInfo;
  edm::EDGetTokenT<edm::View<PileupSummaryInfo>> PUInfoTok;
  edm::InputTag vertexBeamSpot;
  edm::EDGetTokenT<reco::BeamSpot> vertexBeamSpotTok;
  edm::InputTag primaryvertex;
  edm::EDGetTokenT<edm::View<reco::Vertex>> primaryvertexTok;
  edm::InputTag triggerresults;
  edm::EDGetTokenT<edm::TriggerResults> triggerresultsTok;
  edm::InputTag ElectronTag;
  edm::EDGetTokenT<edm::View<pat::Electron>> ElectronTagTok;
  edm::InputTag track;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate>> trackLabelK;
  edm::InputTag isotrack;
  edm::EDGetTokenT<edm::View<pat::IsolatedTrack>> isotrackTok;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> trackBuilderToken;
  bool StoreDeDxInfo_;
  bool verbose_;
  bool TestVerbose_;



  
  const double nominalJpsiMass;
  const double nominalPhiMass;
  const double nominalElectronMass;
  const double nominalMuonMass;
  const double nominalKaonMass;
  const double nominalPionMass;
  const double nominalKstarMass;
  const double nominalBplusMass;

  double JpsiMassWindowBeforeFit_;
  double JpsiMassWindowAfterFit_;
  double JpsiPtCut_;
  double KaonTrackPtCut_;
  double PionTrackPtCut_;
  double BdKaonTrackPtCut_;
  double PhiMassWindowAfterFit_;
  double PhiMassWindowBeforeFit_;
  double BsLowerMassCutBeforeFit_;
  double BsUpperMassCutBeforeFit_;
  double BsLowerMassCutAfterFit_ ;
  double BsUpperMassCutAfterFit_ ;
  double KstarMassWindowBeforeFit_;
  double KstarMassWindowAfterFit_;
  double BdLowerMassCutBeforeFit_;
  double BdUpperMassCutBeforeFit_;
  double BdLowerMassCutAfterFit_;
  double BdUpperMassCutAfterFit_;

  double BsPDGMass_;
  double BdPDGMass_;
  double BpPDGMass_;
  
   int Mu1Truth;


  int match[15][10];
  int match2[15][10];
  int matching[15][10];
  int L1_mu_size0;
  int L1_mu_size1;
  int L1_mu_size;
  int L1_mu_size2;

  unsigned int tagmucounter_;
  unsigned int event_counter_;
  unsigned int elecounter_;
  unsigned int muoncounter_;
  unsigned int jetcounter_;
  
  int    isCowboy       = 0;
  double MuonsDCA       = 999;
  int    VtxIndex       = -99;
  double kaonmass       = 0.493677;
  double pionmass       = 0.139570; 
  double DeltaRKaonJpsi = 999;
  double DeltaRPionJpsi = 999;
  double KPiDCA         = 999;
  double minVtxP = -99.; 
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
};
#endif 

