import FWCore.ParameterSet.Config as cms


#########################################################
# this will produce a ref to the original muon collection
#########################################################
muonsRef = cms.EDFilter('MuonRefSelector',
                   src = cms.InputTag('muons'),
                   cut = cms.string('pt > 0.0'),
                   filter = cms.bool(True)
)

#############################
# Clean Jets Definition
##############################
CleanJets = cms.EDProducer(
    'CleanJets',
    jetSrc = cms.InputTag("ak4PFJets"),
    muonSrc = cms.InputTag("muonsRef"),
    baseMuonSrc = cms.InputTag("muons"),
    PFCandSrc = cms.InputTag("pfIsolatedMuonsEI"),
    genParticleTag = cms.InputTag("genParticles")
)
