import FWCore.ParameterSet.Config as cms

process = cms.Process("CLEANJETS")

#######################################
# Loading all processes and functions
#######################################
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')
process.GlobalTag.globaltag = cms.string('76X_dataRun2_v15')
process.load("Tools.CleanJets.cleanjets_cfi")


#######################################
# Declaring Input and configurations
#######################################
process.source = cms.Source("PoolSource",
         fileNames = cms.untracked.vstring('root://eoscms//eos/cms/store/user/mshi/gg2H2aa2mumutautau_19GeV_step2/gg2H2aa2mumutautau_19GeV_step2_NUM.root')
)


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )


process.CleanJets.outFileName = cms.string('/afs/cern.ch/user/k/ktos/GroupDir/CMSSW_7_6_3/src/Tools/CleanJets/BSUB/DIRNAME/CleanJets_Plots_NUM.root')

#######################################
# HPS Tau Reconstruction alterations 
#######################################
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
process.ak4PFJetTracksAssociatorAtVertex.jets = cms.InputTag('CleanJets', 'ak4PFJetsNoMu', 'CLEANJETS')
process.recoTauAK4PFJets08Region.src = cms.InputTag('CleanJets', 'ak4PFJetsNoMu', 'CLEANJETS')

process.ak4PFJetsLegacyHPSPiZeros.jetSrc = cms.InputTag('CleanJets', 'ak4PFJetsNoMu', 'CLEANJETS')
process.ak4PFJetsRecoTauChargedHadrons.jetSrc = cms.InputTag('CleanJets', 'ak4PFJetsNoMu', 'CLEANJETS')
process.combinatoricRecoTaus.jetSrc = cms.InputTag('CleanJets', 'ak4PFJetsNoMu', 'CLEANJETS')

#######################################
# Configuring Output
#######################################
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('file:DIRNAME_NUM.root'),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('p')),
)

process.p = cms.Path(process.muonsRef*
                     process.CleanJets*
                     process.PFTau)
process.e = cms.EndPath(process.out)

