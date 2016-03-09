// -*- C++ -*-
//
// Package:    GenMatchedRecoObjectProducer
// Class:      GenMatchedRecoObjectProducer
// 
/**\class GenMatchedRecoObjectProducer GenMatchedRecoObjectProducer.cc 
   Tools/GenMatchedRecoObjectProducer/src/GenMatchedRecoObjectProducer.cc

Description: produce a collection of reco objects matched to gen boosted di-tau objects

Implementation:

*/
//
// Original Author:  Rachel Yohay,512 1-010,+41227670495,
//         Created:  Thu Aug 23 11:23:58 CEST 2012
// $Id: GenMatchedRecoObjectProducer.cc,v 1.3 2012/09/25 11:44:34 yohay Exp $
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

template<class T>
class GenMatchedRecoObjectProducer : public edm::EDFilter {
public:
  explicit GenMatchedRecoObjectProducer(const edm::ParameterSet&);
  ~GenMatchedRecoObjectProducer();

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
  edm::EDGetTokenT<reco::GenParticleCollection> genParticleTag_;

  /*input tag for gen particle collection to match
    count on the user to pass in a collection that will not lead to the same reco object being 
    matched to multiple different gen objects
    for example, if the input object is a boosted di-tau pair, only 1 member of the pair should be 
    in the input collection*/
  edm::EDGetTokenT<reco::GenParticleRefVector> selectedGenParticleTag_;

  //input tag for reco object collection
  edm::EDGetTokenT<edm::RefVector<std::vector<T> > > recoObjTag_;

  //input tag for base reco object collection
  edm::EDGetTokenT<std::vector<T> > baseRecoObjTag_;

  //set of parameters for GenTauDecayID class
  edm::ParameterSet genTauDecayIDPSet_;

  //flag indicating whether pT cuts should be applied in determining valid gen objects
  bool applyPTCuts_;

  //flag indicating whether KShorts should be counted as neutral hadrons
  bool countKShort_;

  //dR matching cut
  double dR_;

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
template<class T>
GenMatchedRecoObjectProducer<T>::GenMatchedRecoObjectProducer(const edm::ParameterSet& iConfig) :
  genParticleTag_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticleTag"))),
  selectedGenParticleTag_(consumes<reco::GenParticleRefVector>(iConfig.getParameter<edm::InputTag>("selectedGenParticleTag"))),
  recoObjTag_(consumes<edm::RefVector<std::vector<T> > >(iConfig.getParameter<edm::InputTag>("recoObjTag"))),
  baseRecoObjTag_(consumes<std::vector<T> >(iConfig.getParameter<edm::InputTag>("baseRecoObjTag"))),
  genTauDecayIDPSet_(iConfig.getParameter<edm::ParameterSet>("genTauDecayIDPSet")),
  applyPTCuts_(iConfig.getParameter<bool>("applyPTCuts")),
  countKShort_(iConfig.getParameter<bool>("countKShort")),
  dR_(iConfig.getParameter<double>("dR")),
  minNumGenObjectsToPassFilter_
  (iConfig.getParameter<unsigned int>("minNumGenObjectsToPassFilter"))
{
  //register your products
  produces<edm::RefVector<std::vector<T> > >();

  //now do what ever other initialization is needed
  
}


template<class T>
GenMatchedRecoObjectProducer<T>::~GenMatchedRecoObjectProducer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
template<class T>
bool GenMatchedRecoObjectProducer<T>::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //get base gen particles
  edm::Handle<reco::GenParticleCollection> pGenParticles;
  iEvent.getByToken(genParticleTag_, pGenParticles);

  //get selected gen particles
  edm::Handle<reco::GenParticleRefVector> pSelectedGenParticles;
  iEvent.getByToken(selectedGenParticleTag_, pSelectedGenParticles);

  //get reco object collection
  edm::Handle<edm::RefVector<std::vector<T> > > pRecoObjs;
  iEvent.getByToken(recoObjTag_, pRecoObjs);

  //get base reco object collection
  edm::Handle<std::vector<T> > pBaseRecoObjs;
  iEvent.getByToken(baseRecoObjTag_, pBaseRecoObjs);

  //fill STL container of pointers to reco objects
  std::vector<T*> recoObjPtrs;
  for (typename edm::RefVector<std::vector<T> >::const_iterator iRecoObj = pRecoObjs->begin(); 
       iRecoObj != pRecoObjs->end(); ++iRecoObj) {
    recoObjPtrs.push_back(const_cast<T*>(iRecoObj->get()));
  }

  //fill STL container of selected gen objects
  std::vector<GenTauDecayID> selectedGenObjs;
  for (unsigned int iGenParticle = 0; iGenParticle < pSelectedGenParticles->size(); 
       ++iGenParticle) {
    try {
      const reco::GenParticleRef& genParticleRef = pSelectedGenParticles->at(iGenParticle);
      GenTauDecayID tauDecay(genTauDecayIDPSet_, pGenParticles, 
			     genParticleRef.key()); /*this now assumes that if iGenParticle refers 
						      to a tau, it's a status 2 tau*/

      /*if matching to gen taus: scalar --> status 23 tau --> status 2 tau --> 
	status 1/2 tau decay products ==> the gen particle associated to the GenTauDecayID object 
	is the status 2 tau*/
      if (((fabs(genParticleRef->pdgId()) == GenTauDecayID::TAUPDGID) && 
	   tauDecay.isStatus2DecayProduct()) || 

	  /*if matching to gen non-taus: status whatever mother --> status whatever non-gen tau ==> 
	    the gen particle associated to the GenTauDecayID object should be the gen non-tau 
	    itself*/
	  ((fabs(genParticleRef->pdgId()) != GenTauDecayID::TAUPDGID) && 
	   tauDecay.hasRightMother())) selectedGenObjs.push_back(tauDecay);
    }
    catch (std::string& ex) { throw cms::Exception("GenMatchedRecoObjectProducer") << ex; }
  }

  /*find the decay type of each gen object if it's a tau (needed to get the visible 
    4-vector), otherwise use the gen particle's 4-vector*/
  for (std::vector<GenTauDecayID>::iterator iGenObj = selectedGenObjs.begin(); 
       iGenObj != selectedGenObjs.end(); ++iGenObj) {
    try { iGenObj->tauDecayType(applyPTCuts_, countKShort_); }
    catch (std::string& ex) { throw cms::Exception("GenObjectProducer") << ex; }
  }

  //declare pointer to output collection to produce
  std::auto_ptr<edm::RefVector<std::vector<T> > > 
    genMatchedRecoObjs(new edm::RefVector<std::vector<T> >);

  // //debug
  std::vector<edm::Ref<std::vector<T> > > recoObjsToSave;

  //loop over selected gen particles
  for (std::vector<GenTauDecayID>::iterator iGenObj = selectedGenObjs.begin(); 
       iGenObj != selectedGenObjs.end(); ++iGenObj) {

    //make a dummy LeafCandidate and ref out of the visible 4-vector
    reco::LeafCandidate::LorentzVector visibleGenP4 = iGenObj->getVisibleTauP4();
    std::vector<reco::LeafCandidate> 
      visibleGenParticle(1, reco::LeafCandidate(0.0, visibleGenP4));
    edm::Ref<std::vector<reco::LeafCandidate> > visibleGenParticleRef(&visibleGenParticle, 0);

    //find the nearest reco object to the gen particle
    int nearestRecoObjKey = -1;
    const T* nearestRecoObj = 
      Common::nearestObject(visibleGenParticleRef, recoObjPtrs, nearestRecoObjKey);

    //if nearest reco object is within dR_ of the gen object, save
    if ((nearestRecoObj != NULL) && (reco::deltaR(*nearestRecoObj, *visibleGenParticleRef) < dR_)) {
	genMatchedRecoObjs->
	  push_back(edm::Ref<std::vector<T> >(pBaseRecoObjs, 
					      pRecoObjs->at(nearestRecoObjKey).key()));
    }
  }

  //flag indicating whether right number of gen-matched reco objects were found
  bool foundGenMatchedRecoObject = genMatchedRecoObjs->size() >= minNumGenObjectsToPassFilter_;

  iEvent.put(genMatchedRecoObjs); //this function frees the auto_ptr argument

  //stop processing if no gen-matched objects were found
  return foundGenMatchedRecoObject;
}

// ------------ method called once each job just before starting event loop  ------------
template<class T>
void GenMatchedRecoObjectProducer<T>::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
template<class T>
void GenMatchedRecoObjectProducer<T>::endJob() {
}

// ------------ method called when starting to processes a run  ------------
template<class T>
bool GenMatchedRecoObjectProducer<T>::beginRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a run  ------------
template<class T>
bool GenMatchedRecoObjectProducer<T>::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
template<class T>
bool GenMatchedRecoObjectProducer<T>::beginLuminosityBlock(edm::LuminosityBlock&, 
							   edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
template<class T>
bool GenMatchedRecoObjectProducer<T>::endLuminosityBlock(edm::LuminosityBlock&, 
							 edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  -----------
template<class T>
void 
GenMatchedRecoObjectProducer<T>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
typedef GenMatchedRecoObjectProducer<reco::Muon> GenMatchedMuonProducer;
typedef GenMatchedRecoObjectProducer<reco::PFJet> GenMatchedJetProducer;
typedef GenMatchedRecoObjectProducer<reco::PFTau> GenMatchedTauProducer;
DEFINE_FWK_MODULE(GenMatchedMuonProducer);
DEFINE_FWK_MODULE(GenMatchedJetProducer);
DEFINE_FWK_MODULE(GenMatchedTauProducer);
