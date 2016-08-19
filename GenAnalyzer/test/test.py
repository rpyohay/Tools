### PDG IDs ###

A_PDGID = 36
H_PDGID = 25
MU_PDGID = 13
TAU_PDGID = 15
ANY_PDGID = 0

### Tau decay types ###

TAU_HAD = 0
TAU_MU = 1
TAU_E = 2
TAU_ALL = 3

### Tau hadronic decay types ###

TAU_ALL_HAD = -1
TAU_1PRONG_0NEUTRAL = 0
TAU_1PRONG_1NEUTRAL = 1
TAU_1PRONG_2NEUTRAL = 2
TAU_1PRONG_3NEUTRAL = 3
TAU_1PRONG_NNEUTRAL = 4
TAU_2PRONG_0NEUTRAL = 5
TAU_2PRONG_1NEUTRAL = 6
TAU_2PRONG_2NEUTRAL = 7
TAU_2PRONG_3NEUTRAL = 8
TAU_2PRONG_NNEUTRAL = 9
TAU_3PRONG_0NEUTRAL = 10
TAU_3PRONG_1NEUTRAL = 11
TAU_3PRONG_2NEUTRAL = 12
TAU_3PRONG_3NEUTRAL = 13
TAU_3PRONG_NNEUTRAL = 14
TAU_RARE = 15

### No consideration of pT rank ###

ANY_PT_RANK = -1

### Process declaration ###

import FWCore.ParameterSet.Config as cms
from subprocess import *

process = cms.Process("TEST")

### Verbosity ###

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)
process.MessageLogger.categories.append('TriggerSummaryAnalyzerAOD')
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
    )

### Input ###

import FWCore.Utilities.FileUtils as FileUtils
mylist = FileUtils.loadListFromFile('NMSSM_ggH_a9_H1125_MG5aMCatNLO_2mu2tau_GEN-SIM.txt')
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )
process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(*mylist),
    skipEvents = cms.untracked.uint32(0)
    )

### Select gen tau decays from a-->tautau ###

APSet = cms.PSet(momPDGID = cms.vint32(A_PDGID),
                 chargedHadronPTMin = cms.double(0.0), #should always be 0.0
                 neutralHadronPTMin = cms.double(0.0), #should always be 0.0
                 chargedLeptonPTMin = cms.double(0.0), #should always be 0.0
                 totalPTMin = cms.double(0.0))         #should always be 0.0

process.genTauATauTauSelector = cms.EDFilter(
    'GenObjectProducer',
    genParticleTag = cms.InputTag('genParticles'),
    absMatchPDGIDs = cms.vuint32(TAU_PDGID),     #choose a gen tau...
    sisterAbsMatchPDGID = cms.uint32(TAU_PDGID), #...whose sister is another gen tau...
    genTauDecayIDPSet = APSet,                   #...and whose mother is a pseudoscalar a
    primaryTauDecayType = cms.uint32(TAU_ALL),   #primary tau decay mode is anything...
    sisterTauDecayType = cms.uint32(TAU_ALL),    #...sister tau decay mode is anything
    primaryTauPTRank = cms.int32(ANY_PT_RANK),   #should always be ANY_PT_RANK
    primaryTauHadronicDecayType = cms.int32(TAU_ALL_HAD), #choose TAU_ALL_HAD when the tau decay 
                                                          #type is hadronic and you want any 
                                                          #hadronic mode
    sisterHadronicDecayType = cms.int32(TAU_ALL_HAD),     #choose TAU_ALL_HAD when the tau decay 
                                                          #type is non-hadronic
    primaryTauAbsEtaMax = cms.double(-1.0), #no cut on |eta|
    primaryTauPTMin = cms.double(-1.0),     #no cut on pT
    countSister = cms.bool(True),           #put both taus in the output collection
    applyPTCuts = cms.bool(False),          #should always be False
    countKShort = cms.bool(False),          #should always be False
    minNumGenObjectsToPassFilter = cms.uint32(0), #EDFilter only returns true if >=0 taus are 
                                                  #found
    makeAllCollections = cms.bool(False) #should always be False
    )

### Select gen muons from a-->mumu ###

process.genMuAMuMuSelector = cms.EDFilter(
    'GenObjectProducer',
    genParticleTag = cms.InputTag('genParticles'),
    absMatchPDGIDs = cms.vuint32(MU_PDGID),      #choose a gen muon...
    sisterAbsMatchPDGID = cms.uint32(MU_PDGID),  #...whose sister is another gen muon...
    genTauDecayIDPSet = APSet,                   #...and whose mother is a pseudoscalar a
    primaryTauDecayType = cms.uint32(TAU_ALL),
    sisterTauDecayType = cms.uint32(TAU_ALL),
    primaryTauPTRank = cms.int32(ANY_PT_RANK),   #should always be ANY_PT_RANK
    primaryTauHadronicDecayType = cms.int32(TAU_ALL_HAD), #choose TAU_ALL_HAD when the primary 
                                                          #particle is not a tau
    sisterHadronicDecayType = cms.int32(TAU_ALL_HAD),     #choose TAU_ALL_HAD when the primary 
                                                          #particle is not a tau
    primaryTauAbsEtaMax = cms.double(-1.0), #no cut on |eta|
    primaryTauPTMin = cms.double(-1.0),     #no cut on pT
    countSister = cms.bool(True),           #put both muons in the output collection
    applyPTCuts = cms.bool(False),          #should always be False
    countKShort = cms.bool(False),          #should always be False
    minNumGenObjectsToPassFilter = cms.uint32(0), #EDFilter only returns true if >=0 muons are 
                                                  #found
    makeAllCollections = cms.bool(False) #should always be False
    )

### Select gen pseudoscalars from H-->aa ###

HPSet = APSet.clone()
HPSet.momPDGID = cms.vint32(H_PDGID)

process.genAHAASelector = cms.EDFilter(
    'GenObjectProducer',
    genParticleTag = cms.InputTag('genParticles'),
    absMatchPDGIDs = cms.vuint32(A_PDGID),       #choose a gen pseudoscalar...
    sisterAbsMatchPDGID = cms.uint32(A_PDGID),   #...whose sister is another gen pseudoscalar...
    genTauDecayIDPSet = HPSet,                   #...and whose mother is a scalar H
    primaryTauDecayType = cms.uint32(TAU_ALL),
    sisterTauDecayType = cms.uint32(TAU_ALL),
    primaryTauPTRank = cms.int32(ANY_PT_RANK),   #should always be ANY_PT_RANK
    primaryTauHadronicDecayType = cms.int32(TAU_ALL_HAD), #choose TAU_ALL_HAD when the primary 
                                                          #particle is not a tau
    sisterHadronicDecayType = cms.int32(TAU_ALL_HAD),     #choose TAU_ALL_HAD when the primary 
                                                          #particle is not a tau
    primaryTauAbsEtaMax = cms.double(-1.0), #no cut on |eta|
    primaryTauPTMin = cms.double(-1.0),     #no cut on pT
    countSister = cms.bool(True),           #put both pseudoscalars in the output collection
    applyPTCuts = cms.bool(False),          #should always be False
    countKShort = cms.bool(False),          #should always be False
    minNumGenObjectsToPassFilter = cms.uint32(0), #EDFilter only returns true if >=0 pseudoscalars 
                                                  #are found
    makeAllCollections = cms.bool(False) #should always be False
    )

### Analyze gen H-->aa-->2mu2tau events ###

process.genObjectAnalyzer = cms.EDAnalyzer(
    'GenObjectAnalyzer',
    genParticleTag = cms.InputTag('genParticles'),
    selectedGenTauTag = cms.InputTag('genTauATauTauSelector'),
    selectedGenMuonTag = cms.InputTag('genMuAMuMuSelector'),
    selectedGenPseudoscalarTag = cms.InputTag('genAHAASelector'),
    triggerBitsTag = cms.InputTag("TriggerResults", "", "HLT"),
    triggerEventTag = cms.InputTag("hltTriggerSummaryAOD", "", "HLT"),
    genTauDecayIDPSet = APSet
    )

process.TFileService = cms.Service(
    'TFileService',
    fileName = cms.string('test_MG5aMCatNLO_H125_GEN-SIM.root')
)

### Path ###

process.path = cms.Path(
    process.genTauATauTauSelector*
    process.genMuAMuMuSelector*
    process.genAHAASelector*
    process.genObjectAnalyzer
    )
