// -*- C++ -*-
//
// Package:    TauMatchedRecoObjectProducer
// Class:      TauMatchedRecoObjectProducer
// 
/**\class TauMatchedRecoObjectProducer TauMatchedRecoObjectProducer.cc 
   Tools/TauMatchedRecoObjectProducer/src/TauMatchedRecoObjectProducer.cc

Description: produce a collection of reco objects matched to gen boosted di-tau objects

Implementation:

*/
//
// Original Author:  Rachel Yohay,512 1-010,+41227670495,
//         Created:  Thu Aug 23 11:23:58 CEST 2012
// $Id: TauMatchedRecoObjectProducer.cc,v 1.3 2012/09/25 11:44:34 yohay Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Tools/Common/interface/GenTauDecayID.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "Tools/Common/interface/Common.h"
#include "DataFormats/MuonReco/interface/Muon.h"

//
// class declaration
//

class TauMatchedRecoObjectProducer : public edm::EDFilter {
public:
  explicit TauMatchedRecoObjectProducer(const edm::ParameterSet&);
  ~TauMatchedRecoObjectProducer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
      
  virtual bool beginRun(edm::Run&, edm::EventSetup const&);
  virtual bool endRun(edm::Run&, edm::EventSetup const&);
  virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
  virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

  // ----------member data ---------------------------

  //input tag for base gen particle collection
  edm::InputTag genParticleTag_;

  /*input tag for gen particle collection to match
    count on the user to pass in a collection that will not lead to the same reco object being 
    matched to multiple different gen objects
    for example, if the input object is a boosted di-tau pair, only 1 member of the pair should be 
    in the input collection*/
  edm::InputTag selectedGenParticleTag_;

  //input tag for reco object collection
  edm::InputTag recoObjTag_;

  //input tag for base reco object collection
  edm::InputTag baseRecoObjTag_;

  //set of parameters for GenTauDecayID class
  edm::ParameterSet genTauDecayIDPSet_;

  //flag indicating whether pT cuts should be applied in determining valid gen objects
  bool applyPTCuts_;

  //flag indicating whether KShorts should be counted as neutral hadrons
  bool countKShort_;

  //dR matching cut
  double dR_;

  //If wanting to make the ValueMaps for CleanJets Analysis
  bool ifTauColl_;

  //minimum number of gen objects passing cuts that must be found for event to pass filter
  unsigned int minNumGenObjectsToPassFilter_;
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
TauMatchedRecoObjectProducer::TauMatchedRecoObjectProducer(const edm::ParameterSet& iConfig) :
  genParticleTag_(iConfig.getParameter<edm::InputTag>("genParticleTag")),
  selectedGenParticleTag_(iConfig.getParameter<edm::InputTag>("selectedGenParticleTag")),
  recoObjTag_(iConfig.getParameter<edm::InputTag>("recoObjTag")),
  baseRecoObjTag_(iConfig.getParameter<edm::InputTag>("baseRecoObjTag")),
  genTauDecayIDPSet_(iConfig.getParameter<edm::ParameterSet>("genTauDecayIDPSet")),
  applyPTCuts_(iConfig.getParameter<bool>("applyPTCuts")),
  countKShort_(iConfig.getParameter<bool>("countKShort")),
  dR_(iConfig.getParameter<double>("dR")),
  ifTauColl_(iConfig.getParameter<bool>("ifTauColl")),
  minNumGenObjectsToPassFilter_(iConfig.getParameter<unsigned int>("minNumGenObjectsToPassFilter"))
{
  //register your products
  produces<edm::RefVector<std::vector<reco::PFTau> > >();
  if (ifTauColl_)
  {
    produces<std::vector<reco::PFTau> >( "valMapAccessers" );
    produces<edm::ValueMap<int> >( "TauMatchedMap" );
    produces<edm::ValueMap<int> >( "TauDecayModeMap" );
    produces<edm::ValueMap<double> >( "TauVisiblePtMap" );
  }//if ifTauColl_

  //now do what ever other initialization is needed
  
}


TauMatchedRecoObjectProducer::~TauMatchedRecoObjectProducer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
bool TauMatchedRecoObjectProducer::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  std::cout << "<-------------------New TauMatchedRecoObjectProducer--------------------------->" << std::endl;
  //get base gen particles
  edm::Handle<reco::GenParticleCollection> pGenParticles;
  iEvent.getByLabel(genParticleTag_, pGenParticles);

  //get selected gen particles
  edm::Handle<reco::GenParticleRefVector> pSelectedGenParticles;
  iEvent.getByLabel(selectedGenParticleTag_, pSelectedGenParticles);

  //get reco object collection
  edm::Handle<edm::RefVector<std::vector<reco::PFTau> > > pRecoObjs;
  iEvent.getByLabel(recoObjTag_, pRecoObjs);

  //get base reco object collection
  edm::Handle<std::vector<reco::PFTau> > pBaseRecoObjs;
  iEvent.getByLabel(baseRecoObjTag_, pBaseRecoObjs);

  //fill STL container of pointers to reco objects
  std::vector<reco::PFTau*> recoObjPtrs;
  for (edm::RefVector<std::vector<reco::PFTau> >::const_iterator iRecoObj = pRecoObjs->begin(); iRecoObj != pRecoObjs->end(); ++iRecoObj) 
    recoObjPtrs.push_back(const_cast<reco::PFTau*>(iRecoObj->get()));
  //fill STL container of selected gen objects
  std::vector<GenTauDecayID> selectedGenObjs;
  for (unsigned int iGenParticle = 0; iGenParticle < pSelectedGenParticles->size(); ++iGenParticle) 
  {
    try 
    {
      const reco::GenParticleRef& genParticleRef = pSelectedGenParticles->at(iGenParticle);
      GenTauDecayID tauDecay(genTauDecayIDPSet_, pGenParticles, genParticleRef.key()); /*this now assumes that if iGenParticle refers to a tau, it's a status 2 tau*/

      /*if matching to gen taus: scalar --> status 23 tau --> status 2 tau --> status 1/2 tau decay products ==> the gen particle associated to the GenTauDecayID object is the status 2 tau*/
      if (((fabs(genParticleRef->pdgId()) == GenTauDecayID::TAUPDGID) && tauDecay.isStatus2DecayProduct()) || ((fabs(genParticleRef->pdgId()) != GenTauDecayID::TAUPDGID) && tauDecay.hasRightMother()))
	selectedGenObjs.push_back(tauDecay);
      /* ^^ if matching to gen non-taus: status whatever mother --> status whatever non-gen tau ==> the gen particle associated to the GenTauDecayID object should be the gen non-tau itself*/
    }//try
    catch (std::string& ex) { throw cms::Exception("TauMatchedRecoObjectProducer") << ex; }
  }//for iGenParticle

  /*find the decay type of each gen object if it's a tau (needed to get the visible 
    4-vector), otherwise use the gen particle's 4-vector*/
  for (std::vector<GenTauDecayID>::iterator iGenObj = selectedGenObjs.begin(); iGenObj != selectedGenObjs.end(); ++iGenObj) 
  {
    try { iGenObj->tauDecayType(applyPTCuts_, countKShort_); }
    catch (std::string& ex) { throw cms::Exception("GenObjectProducer") << ex; }
  }//for iGenObj

  //declare pointer to output collection to produce
  std::auto_ptr<edm::RefVector<std::vector<reco::PFTau> > > genMatchedRecoObjs(new edm::RefVector<std::vector<reco::PFTau> >);
  std::auto_ptr<std::vector<reco::PFTau> > SetOfTaus ( new std::vector<reco::PFTau> ); //Only Needed for ifTauColl_ 
  std::vector<int>   genMatchedRecoObjs_Matched; //Only Needed for ifTauColl_ 
  std::vector<int>    genMatchedRecoObjs_DecayMode; //Only Needed for ifTauColl_
  std::vector<double> genMatchedRecoObjs_VisiblePt; //Only Needed for ifTauColl_
  if (ifTauColl_)
  {
    for (std::vector<reco::PFTau>::const_iterator iRecoObj = pBaseRecoObjs->begin(); iRecoObj != pBaseRecoObjs->end(); ++iRecoObj) 
    {
/*      int charge = iRecoObj->charge();
      reco::LeafCandidate::LorentzVector p4 = iRecoObj->p4();
      reco::Candidate::Point point = iRecoObj->point();*/
      SetOfTaus->push_back(*iRecoObj);
      genMatchedRecoObjs_Matched.push_back(0);
      genMatchedRecoObjs_DecayMode.push_back(-1);
      genMatchedRecoObjs_VisiblePt.push_back(-1.0);
    }//for iRecoObj
  }//if ifTauColl_

  // //debug
  std::vector<edm::Ref<std::vector<reco::PFTau> > > recoObjsToSave;

  //loop over selected gen particles
  for (std::vector<GenTauDecayID>::iterator iGenObj = selectedGenObjs.begin(); iGenObj != selectedGenObjs.end(); ++iGenObj) 
  {
    //make a dummy LeafCandidate and ref out of the visible 4-vector
    reco::LeafCandidate::LorentzVector visibleGenP4 = iGenObj->getVisibleTauP4();
    std::vector<reco::LeafCandidate>   visibleGenParticle(1, reco::LeafCandidate(0.0, visibleGenP4));
    edm::Ref<std::vector<reco::LeafCandidate> > visibleGenParticleRef(&visibleGenParticle, 0);

    //find the nearest reco object to the gen particle
    int nearestRecoObjKey = -1;
    const reco::PFTau* nearestRecoObj = Common::nearestObject(visibleGenParticleRef, recoObjPtrs, nearestRecoObjKey);

    //if nearest reco object is within dR_ of the gen object, save
    if ((nearestRecoObj != NULL) && (reco::deltaR(*nearestRecoObj, *visibleGenParticleRef) < dR_)) 
    {
      genMatchedRecoObjs->push_back(edm::Ref<std::vector<reco::PFTau> >(pBaseRecoObjs, pRecoObjs->at(nearestRecoObjKey).key()));
//      std::cout << "\tFound GenRecoMatch" << std::endl;
      if (ifTauColl_)
      {
	int iter = 0;
	edm::Ref<std::vector<reco::PFTau> > RecoTau = edm::Ref<std::vector<reco::PFTau> >(pBaseRecoObjs, pRecoObjs->at(nearestRecoObjKey).key());
        for (std::vector<reco::PFTau>::const_iterator iRecoObj = SetOfTaus->begin(); iRecoObj != SetOfTaus->end(); ++iRecoObj)
        {
          double dPhi = reco::deltaPhi(iRecoObj->phi(), RecoTau->phi() ), dEta = iRecoObj->eta() - RecoTau->eta(), dRTaus = sqrt( dPhi*dPhi + dEta*dEta );
//	  std::cout << "\tRecoTau->pt= " << RecoTau->pt() << "  \tiRecoObj->pt= " << iRecoObj->pt() << "  \tdRTaus= " << dRTaus << std::endl;
          if (dRTaus < .1 && abs(iRecoObj->pt() - RecoTau->pt() ) < 1)
          {
std::cout << "\t\titer= " << iter << "  iRecoObj - SetOfTaus->begin()= " <<  iRecoObj - SetOfTaus->begin() << std::endl;
std::cout << "\t\t<-------FOUND IN SetOfTaus------------>" << std::endl;
            genMatchedRecoObjs_Matched[iter] = 1;
std::cout << "check1" << std::endl;
	    std::pair<reco::PFTau::hadronicDecayMode, GenTauDecayID::DecayType> pair = iGenObj->tauDecayType(applyPTCuts_, countKShort_);
std::cout << "check2" << std::endl;
            genMatchedRecoObjs_DecayMode[iter] = std::get<1>(pair);
std::cout << "check2" << std::endl;
            genMatchedRecoObjs_VisiblePt[iter] = visibleGenP4.Pt();
          }//dRTaus
	  iter++;
        }//for iRecoObj
      }//if ifTauColl_
    }//if nearestRecoObj
  }//for iGenObj
  //flag indicating whether right number of gen-matched reco objects were found
  bool foundGenMatchedRecoObject = genMatchedRecoObjs->size() >= minNumGenObjectsToPassFilter_;

  iEvent.put(genMatchedRecoObjs); //this function frees the auto_ptr argument

  if (ifTauColl_)
  {
    const edm::OrphanHandle<std::vector<reco::PFTau> > SetOfRecoTausProd = iEvent.put( SetOfTaus, "valMapAccessers"  );

    std::auto_ptr<edm::ValueMap<int> > TauMatchedMap(new edm::ValueMap<int>());
    edm::ValueMap<int>::Filler int1Filler(*TauMatchedMap);
    int1Filler.insert(SetOfRecoTausProd, genMatchedRecoObjs_Matched.begin(), genMatchedRecoObjs_Matched.end() );
    int1Filler.fill();
    iEvent.put(TauMatchedMap, "TauMatchedMap" );

    std::auto_ptr<edm::ValueMap<int> > TauDecayModeMap(new edm::ValueMap<int>());
    edm::ValueMap<int>::Filler int2Filler(*TauDecayModeMap);
    int2Filler.insert(SetOfRecoTausProd, genMatchedRecoObjs_DecayMode.begin(), genMatchedRecoObjs_DecayMode.end());
    int2Filler.fill();
    iEvent.put(TauDecayModeMap, "TauDecayModeMap" );

    std::auto_ptr<edm::ValueMap<double> > TauVisiblePtMap(new edm::ValueMap<double>());
    edm::ValueMap<double>::Filler doubleFiller(*TauVisiblePtMap);
    doubleFiller.insert(SetOfRecoTausProd, genMatchedRecoObjs_VisiblePt.begin(), genMatchedRecoObjs_VisiblePt.end());
    doubleFiller.fill();
    iEvent.put(TauVisiblePtMap, "TauVisiblePtMap" );

  }//ifTauColl_
  std::cout << "<-------------------End TauMatchedRecoObjectProducer--------------------------->" << std::endl;
  //stop processing if no gen-matched objects were found
  return foundGenMatchedRecoObject;
}

// ------------ method called once each job just before starting event loop  ------------
void TauMatchedRecoObjectProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void TauMatchedRecoObjectProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool TauMatchedRecoObjectProducer::beginRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool TauMatchedRecoObjectProducer::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool TauMatchedRecoObjectProducer::beginLuminosityBlock(edm::LuminosityBlock&, 
							   edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool TauMatchedRecoObjectProducer::endLuminosityBlock(edm::LuminosityBlock&, 
							 edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  -----------
void 
TauMatchedRecoObjectProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TauMatchedRecoObjectProducer);
