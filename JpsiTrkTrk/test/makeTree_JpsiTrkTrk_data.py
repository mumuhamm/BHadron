import FWCore.ParameterSet.Config as cms
process = cms.Process("MUMU")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi') 
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("FWCore.MessageLogger.MessageLogger_cfi")

from PhysicsTools.PatAlgos.tools.coreTools import *
process.load("PhysicsTools.PatAlgos.patSequences_cff")
from PhysicsTools.PatAlgos.tools.pfTools import *

process.maxEvents = cms.untracked.PSet(  input = cms.untracked.int32(-1) )
process.source = cms.Source("PoolSource",
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                            skipEvents = cms.untracked.uint32(0),
                            fileNames = cms.untracked.vstring(
'root://xrootd-cms.infn.it//store/data/Run2023D/ParkingDoubleMuonLowMass0/MINIAOD/22Sep2023_v1-v1/2550000/0419eec5-0ae4-4732-8f06-6d72dd25a149.root',
#'root://cms-xrd-global.cern.ch//store/mc/RunIIAutumn18MiniAOD/BcToBsPi_JpsiPhiPi_MuMuKKPi_JpsiPhiFilter_TuneCP5_13TeV-bcvegpy2-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/4E457795-DA2E-7340-9F04-9FD229FE5465.root',
#'root://cms-xrd-global.cern.ch//store/mc/RunIIFall17MiniAODv2/BsToJpsiPhi_BMuonFilter_DGamma0_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/FEE2F092-6C44-E811-BAAE-AC1F6B1AF002.root',
#'root://cms-xrd-global.cern.ch//store/data/Run2017C/Charmonium/MINIAOD/12Sep2017-v1/70000/B6F7B207-48A7-E711-9229-FA163EC5EF3C.root',
 )
)

#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
#from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '130X_dataRun3_HLT_v2','')
"""
process.MessageLogger = cms.Service(
                                    "MessageLogger",
                                    #destinations = cms.untracked.vstring('BuToJPsiK_log.txt'),#not checked output
                                    #destinations = cms.untracked.vstring('BsMCtest_log.txt'),#need to be check output
                                    #destinations = cms.untracked.vstring('BdToJPsiKstar_log.txt'),#Alpakaproblem output
                                    default = cms.untracked.PSet( reportEvery = cms.untracked.int32(1000) ) #Marcins way 
)
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
"""


#--PatOverlap, mu/ele--#
process.load("PhysicsTools.PatAlgos.cleaningLayer1.genericTrackCleaner_cfi")
process.cleanPatTracks.checkOverlaps.muons.requireNoOverlaps     = cms.bool(False)
process.cleanPatTracks.checkOverlaps.electrons.requireNoOverlaps = cms.bool(False)

#--Pat Matching --#
#MUON MC-MATCHING VALUES FROM BsMuMu MUON-ID STUDIES
process.load("PhysicsTools.PatAlgos.mcMatchLayer0.muonMatch_cfi")
process.muonMatch.matched = cms.InputTag("prunedGenParticles")
process.muonMatch.maxDeltaR = cms.double(0.12)
process.muonMatch.maxDPtRel = cms.double(0.3)
process.muonMatch.checkCharge = cms.bool(True)
process.muonMatch.resolveAmbiguities = cms.bool(True)
process.muonMatch.resolveByMatchQuality = cms.bool(True)

#--Do we need electron exclusion, I am sure , for student project it does not matter, but for future this requires a bit careful treatment--#
#-- PAT MC MATCHING ele --#
#ELECTRON MC-MATCHING VALUES FROM: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATMCMatching
# object        electron        photon  muon    tau to jet      jet
# maxDPtRel     0.5             1.0     0.5     3.0             3.0
# maxDeltaR     0.5             0.2     0.5     0.1             0.4
process.load("PhysicsTools.PatAlgos.mcMatchLayer0.electronMatch_cfi")
process.electronMatch.matched = cms.InputTag("prunedGenParticles")
process.electronMatch.maxDeltaR = cms.double(0.5)
process.electronMatch.maxDPtRel = cms.double(0.5)
process.electronMatch.checkCharge = cms.bool(True)
process.electronMatch.resolveAmbiguities = cms.bool(True)
process.electronMatch.resolveByMatchQuality = cms.bool(True)

"""
#--PATTracks--Marcin had some question on this part, Now I cant remember--#
process.allKTracks = cms.EDProducer("ConcreteChargedCandidateProducer",
                                    src = cms.InputTag("packedPFCandidates"),
                                    particleType = cms.string('K+')
                                    )
process.allPiTracks = cms.EDProducer("ConcreteChargedCandidateProducer",
                                    src = cms.InputTag("packedPFCandidates"),
                                    particleType = cms.string('pi+')
                                    )
process.kTracks = cms.EDFilter("CandViewRefSelector",
                               src = cms.InputTag("packedPFCandidates"),
                               cut = cms.string("pt > 0.3 & abs(eta) < 2.5")
                               #cut = cms.string("pt > 0.59 & abs(eta) < 2.5")
                               )
process.piTracks = cms.EDFilter("CandViewRefSelector",
                               src = cms.InputTag("packedPFCandidates"),
                               cut = cms.string("pt > 0.3 & abs(eta) < 2.5")
                               #cut = cms.string("pt > 0.59 & abs(eta) < 2.5")
                               )
"""
#-- ANALYZER TAGS AND PARAMETERS --#
process.bVertexAnalysis = cms.EDAnalyzer("JpsiTrkTrk",
                                          isMCstudy                     = cms.bool(False),
                                          genParticlesLabel             = cms.InputTag("prunedGenParticles"),
                                          MuonTag                       = cms.InputTag("slimmedMuons"),
                                          JetCollection                 = cms.InputTag("slimmedJets"),
                                          PUInfo                        = cms.InputTag("slimmedAddPileupInfo"),
                                          vertexBeamSpot                = cms.InputTag("offlineBeamSpot"),
                                          primaryvertex                 = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                          triggerresults                = cms.InputTag("TriggerResults",'',"HLT"),
                                          ElectronTag                   = cms.InputTag("slimmedElectrons"),
                                          track                         = cms.InputTag("packedPFCandidates"),
                                          isotrack                      = cms.InputTag("isolatedTracks"),
                                          StoreDeDxInfo                 = cms.bool( False ),
                                          JpsiMassWindowBeforeFit       = cms.double(0.31), #leave this selection looser than the trigger one for the efficiency calculation
                                          JpsiMassWindowAfterFit        = cms.double(0.150),
                                          MuonPtCut                     = cms.double(4),
                                          JpsiPtCut                     = cms.double(7),
                                          KaonTrackPtCut                = cms.double(0.7),
                                          PionTrackPtCut                = cms.double(0.5),
                                          BdKaonTrackPtCut              = cms.double(0.6),
                                          PhiMassWindowBeforeFit        = cms.double(0.03),
                                          PhiMassWindowAfterFit         = cms.double(0.02),
                                          BsLowerMassCutBeforeFit       = cms.double(4.5),
                                          BsUpperMassCutBeforeFit       = cms.double(6),
                                          BsLowerMassCutAfterFit        = cms.double(5),
                                          BsUpperMassCutAfterFit        = cms.double(6),
                                          KstarMassWindowBeforeFit      = cms.double(0.2),
                                          KstarMassWindowAfterFit       = cms.double(0.15),
                                          BdLowerMassCutBeforeFit       = cms.double(4.5),
                                          BdUpperMassCutBeforeFit       = cms.double(6),
                                          BdLowerMassCutAfterFit        = cms.double(4.9),
                                          BdUpperMassCutAfterFit        = cms.double(5.7),
                                          verbose                       = cms.bool( False ),
                                          TestVerbose                   = cms.bool( False ),
                                          BsPDGMass = cms.double(5.3699),
                                          BdPDGMass = cms.double(5.2794),
                                          BpPDGMass = cms.double(5.2790),
                                          #outputFile                   = cms.untracked.string("BuToJPsiK_AllOktest.root"),
                                          #outputFile                   = cms.untracked.string("BuToJPsiKCt_Dimuon8Trig.root"),
                                          outputFile                   = cms.untracked.string("BdToJpsiKstar_data.root"),#BsJpsiPhi_dG0.root
                                          #outputFile                   = cms.untracked.string("BdToJPsiKstar.root"),
)

# New (easier) Onia2mumu trigger matching
#
#    # Make PAT Muons

process.load("MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff")
import PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi
process.patMuonsWithoutTrigger = PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi.patMuons.clone()


from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import addMCinfo, changeRecoMuonInput, useL1MatchingWindowForSinglets, changeTriggerProcessName, switchOffAmbiguityResolution
useL1MatchingWindowForSinglets(process)
#changeTriggerProcessName(process, "REDIGI36X")
switchOffAmbiguityResolution(process) # Switch off ambiguity resolution: allow multiple reco muons to match to the same trigger muon
process.muonMatchHLTL3.maxDeltaR = 0.1
process.muonMatchHLTL3.maxDPtRel = 10.0
process.muonMatchHLTCtfTrack.maxDeltaR = 0.1
process.muonMatchHLTCtfTrack.maxDPtRel = 10.0
process.muonMatchHLTTrackMu.maxDeltaR = 0.1
process.muonMatchHLTTrackMu.maxDPtRel = 10.0
"""
### ==== Apply some final selection (none by default) ====
process.patMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTrigger"),
    cut = cms.string("p>2 && abs(eta)<2.4"),
)

#apply the scraping event filter here
process.noScraping= cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False), ## Or 'True' to get some per-event info
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.2)
)

process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                           minimumNDOF = cms.uint32(4),
                                           maxAbsZ = cms.double(15),
                                           maxd0 = cms.double(2)
                                           )
"""
# can I do a replace of patMuons with the sequence that includes the trigger matching?
#process.patDefaultSequence.replace(process.patMuons,process.patMuonsWithoutTrigger * process.patTriggerMatching * process.patMuons)
#process.vertex = cms.Path(process.inclusiveVertexing * process.inclusiveMergedVertices * process.selectedVertices * process.bcandidates)
#process.pat = cms.Path( process.patDefaultSequence )
#process.pat = cms.Path(process.patDefaultSequence)
#print(process.pat)

#process.ntup = cms.Path(process.allPiTracks * process.allKTracks * process.kTracks * process.piTracks * process.bVertexAnalysis )
process.ntup = cms.Path(process.bVertexAnalysis )
#process.filter = cms.Path(process.noScraping)
process.schedule = cms.Schedule(process.ntup)

