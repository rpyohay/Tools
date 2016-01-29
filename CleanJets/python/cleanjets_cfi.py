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
    PFCandSrc = cms.InputTag("pfIsolatedMuonsEI"),
    outFileName = cms.string('file:/afs/cern.ch/user/k/ktos/NMSSM_Analysis/CMSSW_7_4_12_patch4/src/Tools/CleanJets/BSUB/DIRNAME/CleanJets_Plots.root')
)
