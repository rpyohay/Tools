# tools

## GenObjectProducer

This module produces a reco::GenParticleRefVector (a vector of reco::GenParticleRef objects) of gen particles chosen by PDG ID, mother, sister, and, if applicable, tau decay type.  Here are some example Python configuration snippets.  In all CMSSW cfg files, the following macros should be defined in order to use the GenObjectProducer module.

```python
### PDG IDs ###

A_PDGID = 36
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
```

Select both gen muons from a-->mumu:

```python
### Choose mother (pseudoscalar a) ###

AMuMuPSet = cms.PSet(momPDGID = cms.vint32(A_PDGID),
                     chargedHadronPTMin = cms.double(0.0), #should always be 0.0
                     neutralHadronPTMin = cms.double(0.0), #should always be 0.0
                     chargedLeptonPTMin = cms.double(0.0), #should always be 0.0
                     totalPTMin = cms.double(0.0)) #should always be 0.0

process.genAMuSelector = cms.EDFilter(
    'GenObjectProducer',
    genParticleTag = cms.InputTag('genParticles'),
    absMatchPDGIDs = cms.vuint32(MU_PDGID),     #choose a gen muon...
    sisterAbsMatchPDGID = cms.uint32(MU_PDGID), #...whose sister is another gen muon...
    genTauDecayIDPSet = AMuMuPSet,              #...and whose mother is a pseudoscalar a
    primaryTauDecayType = cms.uint32(TAU_ALL),  #choose TAU_ALL when the gen particle is not a tau
    sisterTauDecayType = cms.uint32(TAU_ALL),   #choose TAU_ALL when the gen particle sister is not a tau
    primaryTauPTRank = cms.int32(ANY_PT_RANK),  #should always be ANY_PT_RANK
    primaryTauHadronicDecayType = cms.int32(TAU_ALL_HAD), #choose TAU_ALL_HAD when the gen particle is not a tau
    sisterHadronicDecayType = cms.int32(TAU_ALL_HAD),     #choose TAU_ALL_HAD when the gen particle sister is not a tau
    primaryTauAbsEtaMax = cms.double(-1.0), #no cut on gen particle |eta|
    primaryTauPTMin = cms.double(-1.0),     #no cut on gen particle pT
    countSister = cms.bool(True),           #True if you want to put both muons in the output collection, False if just one
    applyPTCuts = cms.bool(False),          #should always be False
    countKShort = cms.bool(False),          #should always be False
    minNumGenObjectsToPassFilter = cms.uint32(2), #EDFilter only returns true if >=2 muons are found
    makeAllCollections = cms.bool(False) #should always be False
    )
```

Select only 1 gen muon from a-->mumu:

```python
### Choose mother (pseudoscalar a) ###

AMuMuPSet = cms.PSet(momPDGID = cms.vint32(A_PDGID),
                     chargedHadronPTMin = cms.double(0.0), #should always be 0.0
                     neutralHadronPTMin = cms.double(0.0), #should always be 0.0
                     chargedLeptonPTMin = cms.double(0.0), #should always be 0.0
                     totalPTMin = cms.double(0.0)) #should always be 0.0

process.genAMuSelector = cms.EDFilter(
    'GenObjectProducer',
    genParticleTag = cms.InputTag('genParticles'),
    absMatchPDGIDs = cms.vuint32(MU_PDGID),     #choose a gen muon...
    sisterAbsMatchPDGID = cms.uint32(MU_PDGID), #...whose sister is another gen muon...
    genTauDecayIDPSet = AMuMuPSet,              #...and whose mother is a pseudoscalar a
    primaryTauDecayType = cms.uint32(TAU_ALL),  #choose TAU_ALL when the gen particle is not a tau
    sisterTauDecayType = cms.uint32(TAU_ALL),   #choose TAU_ALL when the gen particle sister is not a tau
    primaryTauPTRank = cms.int32(ANY_PT_RANK),  #should always be ANY_PT_RANK
    primaryTauHadronicDecayType = cms.int32(TAU_ALL_HAD), #choose TAU_ALL_HAD when the gen particle is not a tau
    sisterHadronicDecayType = cms.int32(TAU_ALL_HAD),     #choose TAU_ALL_HAD when the gen particle sister is not a tau
    primaryTauAbsEtaMax = cms.double(-1.0), #no cut on gen particle |eta|
    primaryTauPTMin = cms.double(-1.0),     #no cut on gen particle pT
    countSister = cms.bool(False),          #True if you want to put both muons in the output collection, False if just one
    applyPTCuts = cms.bool(False),          #should always be False
    countKShort = cms.bool(False),          #should always be False
    minNumGenObjectsToPassFilter = cms.uint32(1), #EDFilter only returns true if >=1 muon is found
    makeAllCollections = cms.bool(False) #should always be False
    )
```

Select only 1 gen muon from a-->mumu with pT > 5 GeV and |eta| < 2.1:

```python
### Choose mother (pseudoscalar a) ###

AMuMuPSet = cms.PSet(momPDGID = cms.vint32(A_PDGID),
                     chargedHadronPTMin = cms.double(0.0), #should always be 0.0
                     neutralHadronPTMin = cms.double(0.0), #should always be 0.0
                     chargedLeptonPTMin = cms.double(0.0), #should always be 0.0
                     totalPTMin = cms.double(0.0)) #should always be 0.0

process.genAMuSelector = cms.EDFilter(
    'GenObjectProducer',
    genParticleTag = cms.InputTag('genParticles'),
    absMatchPDGIDs = cms.vuint32(MU_PDGID),     #choose a gen muon...
    sisterAbsMatchPDGID = cms.uint32(MU_PDGID), #...whose sister is another gen muon...
    genTauDecayIDPSet = AMuMuPSet,              #...and whose mother is a pseudoscalar a
    primaryTauDecayType = cms.uint32(TAU_ALL),  #choose TAU_ALL when the gen particle is not a tau
    sisterTauDecayType = cms.uint32(TAU_ALL),   #choose TAU_ALL when the gen particle sister is not a tau
    primaryTauPTRank = cms.int32(ANY_PT_RANK),  #should always be ANY_PT_RANK
    primaryTauHadronicDecayType = cms.int32(TAU_ALL_HAD), #choose TAU_ALL_HAD when the gen particle is not a tau
    sisterHadronicDecayType = cms.int32(TAU_ALL_HAD),     #choose TAU_ALL_HAD when the gen particle sister is not a tau
    primaryTauAbsEtaMax = cms.double(2.1),  #|eta| < 2.1 on primary muon
    primaryTauPTMin = cms.double(5.0),      #pT > 5 GeV on primary muon
    countSister = cms.bool(False),          #True if you want to put both muons in the output collection, False if just one
    applyPTCuts = cms.bool(False),          #should always be False
    countKShort = cms.bool(False),          #should always be False
    minNumGenObjectsToPassFilter = cms.uint32(1), #EDFilter only returns true if >=1 muon is found satisfying pT and |eta| cuts
    makeAllCollections = cms.bool(False) #should always be False
    )
```

Select gen muon from a-->tau(-->mu)tau(-->had) with pT > 5 GeV and |eta| < 2.1:

```python
### Choose mother (pseudoscalar a) ###

ATauTauPSet = cms.PSet(momPDGID = cms.vint32(A_PDGID),
                       chargedHadronPTMin = cms.double(0.0), #should always be 0.0
                       neutralHadronPTMin = cms.double(0.0), #should always be 0.0
                       chargedLeptonPTMin = cms.double(0.0), #should always be 0.0
                       totalPTMin = cms.double(0.0)) #should always be 0.0

process.genATauMuSelector = cms.EDFilter(
    'GenObjectProducer',
    genParticleTag = cms.InputTag('genParticles'),
    absMatchPDGIDs = cms.vuint32(TAU_PDGID),     #choose a gen tau...
    sisterAbsMatchPDGID = cms.uint32(TAU_PDGID), #...whose sister is another gen tau...
    genTauDecayIDPSet = ATauTauPSet,             #...and whose mother is a pseudoscalar a
    primaryTauDecayType = cms.uint32(TAU_MU),    #primary tau decay mode is mu...
    sisterTauDecayType = cms.uint32(TAU_HAD),    #...sister tau decay mode is hadronic
    primaryTauPTRank = cms.int32(ANY_PT_RANK),  #should always be ANY_PT_RANK
    primaryTauHadronicDecayType = cms.int32(TAU_ALL_HAD), #choose TAU_ALL_HAD when the tau decay type is non-hadronic
    sisterHadronicDecayType = cms.int32(TAU_ALL_HAD),     #choose TAU_ALL_HAD when the tau decay type is hadronic and you want any hadronic mode
    primaryTauAbsEtaMax = cms.double(2.1),  #|eta| < 2.1 on muon from tau-->mu
    primaryTauPTMin = cms.double(5.0),      #pT > 5 GeV on muon from tau-->mu
    countSister = cms.bool(False),          #only put the muon from tau-->mu in the output collection (i.e. the object in the output collection will be the gen particle with |PDG ID| = 13 and status = 1 that is decayed from the tau
    applyPTCuts = cms.bool(False),          #should always be False
    countKShort = cms.bool(False),          #should always be False
    minNumGenObjectsToPassFilter = cms.uint32(1), #EDFilter only returns true if >=1 tau-->mu is found satisfying pT, |eta|, and decay mode cuts
    makeAllCollections = cms.bool(False) #should always be False
    )
```

Select gen 1-prong tau from a-->tau(-->1-prong)tau(-->3-prong):

```python
### Choose mother (pseudoscalar a) ###

ATauTauPSet = cms.PSet(momPDGID = cms.vint32(A_PDGID),
                       chargedHadronPTMin = cms.double(0.0), #should always be 0.0
                       neutralHadronPTMin = cms.double(0.0), #should always be 0.0
                       chargedLeptonPTMin = cms.double(0.0), #should always be 0.0
                       totalPTMin = cms.double(0.0)) #should always be 0.0

process.genATau1ProngTau3ProngSelector = cms.EDFilter(
    'GenObjectProducer',
    genParticleTag = cms.InputTag('genParticles'),
    absMatchPDGIDs = cms.vuint32(TAU_PDGID),     #choose a gen tau...
    sisterAbsMatchPDGID = cms.uint32(TAU_PDGID), #...whose sister is another gen tau...
    genTauDecayIDPSet = ATauTauPSet,             #...and whose mother is a pseudoscalar a
    primaryTauDecayType = cms.uint32(TAU_HAD),   #primary tau decay mode is hadronic...
    sisterTauDecayType = cms.uint32(TAU_HAD),    #...sister tau decay mode is hadronic
    primaryTauPTRank = cms.int32(ANY_PT_RANK),  #should always be ANY_PT_RANK
    primaryTauHadronicDecayType = cms.int32(TAU_1PRONG_0NEUTRAL), #primary tau hadronic decay mode is 1-prong...
    sisterHadronicDecayType = cms.int32(TAU_3PRONG_0NEUTRAL),     #...sister tau hadronic decay mode is 3-prong
    primaryTauAbsEtaMax = cms.double(-1.0),
    primaryTauPTMin = cms.double(-1.0),
    countSister = cms.bool(False),          #only put the 1-prong tau in the output collection (i.e. the object in the output collection will be the gen particle with |PDG ID| = 15 and status = 2 whose decay products are the tau neutrino and charged hadron
    applyPTCuts = cms.bool(False),          #should always be False
    countKShort = cms.bool(False),          #should always be False
    minNumGenObjectsToPassFilter = cms.uint32(1), #EDFilter only returns true if >=1 tau-->1-prong is found with sister tau-->3-prong
    makeAllCollections = cms.bool(False) #should always be False
    )
```

## GenMatchedRecoObjectProducer

Once the collection of gen particles to match is produced via the GenObjectProducer and a set of reco objects is produced via, for example, MuonRefSelector, GenMatchedRecoObjectProducer is run to select only those reco objects that match the specified gen particles.  Let's say the gen particles to match are produced by a module called genATauMuSelector configured as above, and a reco::MuonRefVector of reco muons with pT > 5 GeV and |eta| < 2.1 are produced according to the following snippet:

```python
process.recoMuonSelector = cms.EDFilter('MuonRefSelector',
                                        src = cms.InputTag('muons'),
                                        cut = cms.string('(pt > 5.0) && (abs(eta) < 2.1)'),
                                        filter = cms.bool(True)
                                        )
```

The gen matching can be done via the following snippet:

```python
process.genATauMuMatchedRecoMuonSelector = cms.EDFilter(
    'GenMatchedMuonProducer',
    genParticleTag = cms.InputTag('genParticles'),
    selectedGenParticleTag = cms.InputTag('genATauMuSelector'), #must be a reco::GenParticleRefVector
    recoObjTag = cms.InputTag('recoMuonSelector'),              #must be a reco::MuonRefVector
    baseRecoObjTag = cms.InputTag('muons'),
    genTauDecayIDPSet = ATauTauPSet,      #need to know the pseudoscalar a mother
    applyPTCuts = cms.bool(False),        #should always be false
    countKShort = cms.bool(False),        #should always be false
    dR = cms.double(0.1),                 #dR criteria for matching
    minNumGenObjectsToPassFilter = cms.uint32(1) #EDFilter returns true if >=1 gen-matched reco muon is found
    )
```

The CMSSW Path object would then be:

```python
process.p = cms.Path(
    process.genATauMuSelector*
    process.recoMuonSelector*
    process.genATauMuMatchedRecoMuonSelector*
    process.myGenMatchedRecoMuonAnalyzer
    )
```

Note that GenMatchedRecoObjectProducer is a template class, with (so far) 3 implementations defined:

* GenMatchedMuonProducer: recoObjTag must refer to a reco::MuonRefVector
* GenMatchedJetProducer: recoObjTag must refer to a reco::PFJetRefVector
* GenMatchedTauProducer: recoObjTag must refer to a reco::PFTauRefVector

## Hadronic tau decays

The algorithm coded in GenTauDecayID recursively counts the number of status 1 charged pions, status 1 charged kaons, and status 2 neutral pions in the decay of a tau to determine the hadronic decay mode and the visible 4-vector.  So, all intermediate neutrals are decayed to pions and kaons before performing the final count.

Klongs and Kshorts, which are status 1 in Pythia 8, are ignored in the determination of decay mode.  Decays via these particles account for 2% of tau decays and would be not be reconstructed as a valid gen tau decay in this scheme.  The HPS tau reconstruction algorithm is not designed to reconstruct decays via Klongs and Kshorts.