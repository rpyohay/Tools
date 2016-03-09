// -*- C++ -*-
//
// Package:    GenMatchedRecoObjectProducer
// Class:      CustomTauSelector
// 
/**\class CustomTauSelector CustomTauSelector.cc 
   BoostedTauAnalysis/GenMatchedRecoObjectProducer/src/CustomTauSelector.cc

 Description: create a collection of custom selected hadronic taus to put in the event, and stop 
 processing if no taus are selected

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Rachel Yohay,512 1-010,+41227670495,
//         Created:  Fri Aug 24 17:10:12 CEST 2012
// $Id: CustomTauSelector.cc,v 1.6 2012/12/12 16:02:05 yohay Exp $
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
#include "Tools/Common/interface/Common.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TFile.h"
//
// class declaration
//
using namespace edm;
template<class T>
class CustomTauSelector : public edm::EDFilter {
public:
  explicit CustomTauSelector(const edm::ParameterSet&);
  ~CustomTauSelector();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob();
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob();
  virtual bool beginRun(edm::Run&, edm::EventSetup const&);
  virtual bool endRun(edm::Run&, edm::EventSetup const&);
  virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
  virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

  // ----------member data ---------------------------

  TFile* out_;
  std::string outFileName_;

  //input tag for reco tau collection
  edm::EDGetTokenT<reco::PFTauRefVector> tauTag_;

  //input tag for base tau collection
  edm::EDGetTokenT<reco::PFTauCollection> baseTauTag_;

  //input tag for tau isolation energy
  edm::EDGetTokenT<reco::PFTauDiscriminator> tauHadIsoTag_;

  //input tag for clean jet collection
  edm::EDGetTokenT<reco::PFJetCollection> jetTag_;

  //input tag for map jet muon removal decisions
  edm::EDGetTokenT<edm::ValueMap<bool> > muonRemovalDecisionTag_;

  //input tag for overlap candidate collection
  edm::EDGetTokenT<edm::RefVector<std::vector<T> > > overlapCandTag_;

  //vector of input tags, 1 for each discriminator the tau should pass
  //std::vector<edm::EDGetTokenT<reco::PFTauDiscriminator> > tauDiscriminatorTags_;
  //std::vector<edm::InputTag> tauDiscriminatorTags_;
  typedef edm::EDGetTokenT<reco::PFTauDiscriminator> PFTauDiscriminatorToken;
  std::vector<PFTauDiscriminatorToken> tauDiscriminatorTags_; 

  //flag indicating whether the selected taus should pass or fail the discriminator
  bool passDiscriminator_;

  //pT cut
  double pTMin_;

  //|eta| cut
  double etaMax_;

  //maximum isolation cut
  double isoMax_;

  //overlap candidate matching cut
  double dR_;

  //minimum number of objects that must be found to pass the filter
  unsigned int minNumObjsToPassFilter_;
  //  std::map<std::string, TH1D*> histos1D_;
  TH1F *TestDRValMap_;
  TH1F *TestDR_;
  TH1F *TauMuBranchingRatio_;
  TH1F *TauPt_;
  TH1F *NPassing_;

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
CustomTauSelector<T>::CustomTauSelector(const edm::ParameterSet& iConfig) :
  outFileName_(iConfig.getParameter<std::string>("outFileName")),
  tauTag_(iConfig.existsAs<edm::InputTag>("tauTag") ? 
	  consumes<reco::PFTauRefVector>(iConfig.getParameter<edm::InputTag>("tauTag")) : edm::EDGetTokenT<reco::PFTauRefVector>()),
  baseTauTag_(consumes<reco::PFTauCollection>(iConfig.getParameter<edm::InputTag>("baseTauTag"))),
  tauHadIsoTag_(consumes<reco::PFTauDiscriminator>(iConfig.getParameter<edm::InputTag>("tauHadIsoTag"))),
  jetTag_(iConfig.existsAs<edm::InputTag>("jetTag") ? 
	  consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("jetTag")) : edm::EDGetTokenT<reco::PFJetCollection>()),
  muonRemovalDecisionTag_(iConfig.existsAs<edm::InputTag>("muonRemovalDecisionTag") ? 
			  consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("muonRemovalDecisionTag")) : 
			  edm::EDGetTokenT<edm::ValueMap<bool> >()),
  overlapCandTag_(iConfig.existsAs<edm::InputTag>("overlapCandTag") ? 
		  consumes<edm::RefVector<std::vector<T> > >(iConfig.getParameter<edm::InputTag>("overlapCandTag")) : edm::EDGetTokenT<edm::RefVector<std::vector<T> > >()),
//  tauDiscriminatorTags_(consumes<std::vector<edm::EDGetTokenT<reco::PFTauDiscriminator> > >(iConfig.getParameter<std::vector<edm::InputTag> >("tauDiscriminatorTags"))),
//  tauDiscriminatorTags_(iConfig.getParameter<std::vector<edm::InputTag> >("tauDiscriminatorTags")),
  passDiscriminator_(iConfig.getParameter<bool>("passDiscriminator")),
  pTMin_(iConfig.getParameter<double>("pTMin")),
  etaMax_(iConfig.getParameter<double>("etaMax")),
  isoMax_(iConfig.getParameter<double>("isoMax")),
  dR_(iConfig.getParameter<double>("dR")),
  minNumObjsToPassFilter_(iConfig.getParameter<unsigned int>("minNumObjsToPassFilter"))
//  histos1D_()
{
  typedef std::vector<edm::InputTag> vInputTag;
  vInputTag srcTauDiscriminatorTags_ = iConfig.getParameter<vInputTag>("tauDiscriminatorTags");
  for ( vInputTag::const_iterator it = srcTauDiscriminatorTags_.begin(); it != srcTauDiscriminatorTags_.end(); ++it )
      tauDiscriminatorTags_.push_back(/*ccollector.*/consumes<reco::PFTauDiscriminator>(*it));




  if (((jetTag_.isUninitialized()) && !(muonRemovalDecisionTag_.isUninitialized())) || 
      (!(jetTag_.isUninitialized()) && (muonRemovalDecisionTag_.isUninitialized()))) {
    std::cerr << "Warning: only one of jetTag or muonRemovalDecisionTag was supplied.  No ";
    std::cerr << "decision on tau seed jet will be made.\n";
  }
  if ((overlapCandTag_.isUninitialized()) && !(jetTag_.isUninitialized()) && 
      !(muonRemovalDecisionTag_.isUninitialized())) {
    std::cerr << "Warning: both jetTag and muonRemovalDecisionTag were supplied, but not ";
    std::cerr << "overlapCandTag.  Overlap with overlap candidates will not be checked.\n";
  }
  produces<reco::PFTauRefVector>();
}

template<class T>
CustomTauSelector<T>::~CustomTauSelector()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
template<class T>
bool CustomTauSelector<T>::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
std::cout << "<---------------------In CustomTauSelector-------------------->" << std::endl;
  //create pointer to output collection
  std::auto_ptr<reco::PFTauRefVector> tauColl(new reco::PFTauRefVector);

  //get taus
  edm::Handle<reco::PFTauRefVector> pTaus;
  if (tauTag_.isUninitialized()) {}
  else iEvent.getByToken(tauTag_, pTaus);

  //get base tau collection
  edm::Handle<reco::PFTauCollection> pBaseTaus;
  iEvent.getByToken(baseTauTag_, pBaseTaus);

  //get hadronic tau deltaBeta-corrected isolation
  edm::Handle<reco::PFTauDiscriminator> pTauHadIso;
  iEvent.getByToken(tauHadIsoTag_, pTauHadIso);

  //get tau discriminators
  std::vector<edm::Handle<reco::PFTauDiscriminator> > pTauDiscriminators(tauDiscriminatorTags_.size(), edm::Handle<reco::PFTauDiscriminator>());
  for (std::vector<edm::EDGetTokenT<reco::PFTauDiscriminator> >::const_iterator iTag = tauDiscriminatorTags_.begin();  iTag != tauDiscriminatorTags_.end(); ++iTag) 
    iEvent.getByToken(*iTag, pTauDiscriminators[iTag - tauDiscriminatorTags_.begin()]);

/*  std::vector<edm::Handle<reco::PFTauDiscriminator> >  pTauDiscriminators(tauDiscriminatorTags_.size(), edm::Handle<reco::PFTauDiscriminator>());
  for (std::vector<edm::InputTag>::const_iterator iTag = tauDiscriminatorTags_.begin(); iTag != tauDiscriminatorTags_.end(); ++iTag) 
    iEvent.getByLabel(*iTag, pTauDiscriminators[iTag - tauDiscriminatorTags_.begin()]);
  
  for ( std::vector<PFTauDiscriminatorToken>::const_iterator it = tauDiscriminatorTags_.begin(); it != tauDiscriminatorTags_.end(); ++it ) 
  {
    edm::Handle<reco::PFTauDiscriminator> particlesNotToBeFiltered;
    iEvent.getByToken(*it, particlesNotToBeFiltered);
  }//for
*/

  //get jet collection
  edm::Handle<reco::PFJetCollection> pJets;
  if (jetTag_.isUninitialized()) {}
  else iEvent.getByToken(jetTag_, pJets);

  //get map of jet muon removal decisions
  edm::Handle<edm::ValueMap<bool> > pMuonRemovalDecisions;
  if (muonRemovalDecisionTag_.isUninitialized()) {}
  else iEvent.getByToken(muonRemovalDecisionTag_, pMuonRemovalDecisions);

  //get overlap candidates
  edm::Handle<edm::RefVector<std::vector<T> > > pOverlapCands;
  if (overlapCandTag_.isUninitialized()) {}
  else iEvent.getByToken(overlapCandTag_, pOverlapCands);

  //fill STL container of pointers to overlap candidates
  std::vector<T*> overlapCandPtrs;
  if (pOverlapCands.isValid()) {
    for (typename edm::RefVector<std::vector<T> >::const_iterator iOverlapCand = 
	   pOverlapCands->begin(); iOverlapCand != pOverlapCands->end(); 
	 ++iOverlapCand) { overlapCandPtrs.push_back(const_cast<T*>(iOverlapCand->get())); }
  }

//   //debug
//   std::cerr << "Jets " << pJets.isValid() << std::endl;
 //  std::cerr << "Value map " << pMuonRemovalDecisions.isValid() << std::endl;
  // std::cerr << "Taus " << pTaus.isValid() << std::endl;
  // std::cerr << "Base taus " << pBaseTaus.isValid() << std::endl;
  // std::cout<<"pBaseTaus.size()=="<<pBaseTaus->size()<<std::endl;

  //fill STL container with taus passing specified discriminators in specified eta and pT range
  std::vector<reco::PFTauRef> taus = pTaus.isValid() ? 
    Common::getRecoTaus(pTaus, pBaseTaus, pTauDiscriminators, pTauHadIso, pTMin_, etaMax_, 
			passDiscriminator_, isoMax_) : 
    Common::getRecoTaus(pBaseTaus, pTauDiscriminators, pTauHadIso, pTMin_, etaMax_, 
			passDiscriminator_, isoMax_);
 
  //loop over selected taus
  unsigned int nPassingTaus = 0;
  for (std::vector<reco::PFTauRef>::const_iterator iTau = taus.begin(); iTau != taus.end(); ++iTau) 
  {
    if (fabs((*iTau)->charge() ) == 1)
    {
      //find the nearest overlap candidate to the tau
      int nearestMuonIndex = -1;
      const reco::Candidate* nearestMuon = Common::nearestObject(*iTau, overlapCandPtrs, nearestMuonIndex);
      if( (*pMuonRemovalDecisions)[(*iTau)->jetRef()]&& pOverlapCands.isValid())
      {
        TestDRValMap_->Fill(reco::deltaR(**iTau, *nearestMuon));
        std::cout<<"nearestMuon!=NULL is true or not=="<<(nearestMuon!=NULL)<<std::endl;
      }
   
      //if tau doesn't overlap with overlap candidate (or no overlap checking requested)...
      if(pOverlapCands.isValid())   {
        TestDR_->Fill(reco::deltaR(**iTau,*nearestMuon));
        std::cout<<"nearestMuon!=NULL is true or not=="<<(nearestMuon!=NULL)<<std::endl;
        std::cout<<"(reco::deltaR(**iTau, *nearestMuon) > dR_)=="<<(reco::deltaR(**iTau, *nearestMuon) > dR_)<<std::endl;
      }
     
      if (!(pOverlapCands.isValid()) || 
  	((nearestMuon != NULL) && (reco::deltaR(**iTau, *nearestMuon) > dR_))) {
        /*...if jet collection and muon removal decision map exist, fill output collection if tau is 
  	matched to jet tagged for muon removal*/
        if (pJets.isValid() && pMuonRemovalDecisions.isValid()) {
  //        std::cout<<"((*iTau)->jetRef())"<< (*pMuonRemovalDecisions)[(*iTau)->jetRef()]<<std::endl;
  	if ((*pMuonRemovalDecisions)[(*iTau)->jetRef()]) {
  	  TauPt_->Fill((*iTau)->pt());
            tauColl->push_back(*iTau);
  	  ++nPassingTaus;
  	}
        }
  
        /*...if jet collection and muon removal decision map do not exist, assume no selection on 
  	tau seed jet is desired and fill output collection*/
        else {
  	tauColl->push_back(*iTau);
  	++nPassingTaus;
        }//pjets.isValid()==0||pMuonRemovalDecisions.isValid()==0
      }//outside if(!(pOver... but already store info of nPassingTaus
      if(!pOverlapCands.isValid()){
        std::cout<<"inside loop of iTau, if no overlap tag is provided"<<std::endl;
        std::cout<<"((*iTau)->jetRef())"<< (*pMuonRemovalDecisions)[(*iTau)->jetRef()]<<std::endl;
        std::cout<<"pJets.isValid()"<<pJets.isValid()<<std::endl;
        std::cout<<"pMuonRemovalDecisions.isValid()"<<pMuonRemovalDecisions.isValid()<<std::endl;
        std::cout<<"nPassingTaus"<<nPassingTaus<<std::endl;
        TauMuBranchingRatio_->Fill((*pMuonRemovalDecisions)[(*iTau)->jetRef()]);
      }
    }//if iTaucharge
  }
  iEvent.put(tauColl);

  //if not enough taus passing cuts were found in this event, stop processing
  NPassing_->Fill((nPassingTaus >= minNumObjsToPassFilter_));
  return (nPassingTaus >= minNumObjsToPassFilter_);
}

// ------------ method called once each job just before starting event loop  ------------
template<class T>
void CustomTauSelector<T>::beginJob()
{
//edm::Service< TFileService > fileService; 
//thistos1D_["testDRValMap"]=fileService->make<TH1D>("testDRValMap","testDRValMap", 23, 0,7.0);
//histos1D_["testDR"]=fileService->make<TH1D>("test","test",23,0,7.0);
//histos1D_["tau_mu branching ratio"]=fileService->make<TH1D>("tau_mu branching ratio", "tau_mu branching rattio",2,0,2);
//histos1D_["TauPt"]=fileService->make<TH1D>("TauPt", "TauPt",100,0,300);
//histos1D_["nPassing"]=fileService->make<TH1D>("nPassing","nPassing",2,0,2);
out_ = new TFile(outFileName_.c_str(), "RECREATE");
TestDRValMap_ = new TH1F("TestDRValMap", ";dR valMap;", 23, 0, 7.0);
TestDR_ = new TH1F("TestDR", ";dR;", 23, 0, 7.0);
TauMuBranchingRatio_ = new TH1F("tau_mu branching ratio", "tau_mu branching rattio",2,0,2);
TauPt_ = new TH1F("TauPt", "TauPt",100,0,300);
NPassing_ = new TH1F("nPassing","nPassing",2,0,2);

}

// ------------ method called once each job just after ending the event loop  ------------
template<class T>
void CustomTauSelector<T>::endJob() {
  out_->cd();

  TCanvas TestDRValMapCanvas_("TestDRValMapCanvas","",600,600);
  TCanvas TestDRCanvas_("TestDRCanvas","",600,600);
  TCanvas TauMuBranchingRatioCanvas_("TauMuBranchingRatioCanvas","",600,600);
  TCanvas TauPtCanvas_("TauPtCanvas","",600,600);
  TCanvas NPassingCanvas_("NPassingCanvas","",600,600);

  Common::draw1DHistograms(TestDRValMapCanvas_, TestDRValMap_);
  Common::draw1DHistograms(TestDRCanvas_, TestDR_);
  Common::draw1DHistograms(TauMuBranchingRatioCanvas_, TauMuBranchingRatio_);
  Common::draw1DHistograms(TauPtCanvas_, TauPt_);
  Common::draw1DHistograms(NPassingCanvas_, NPassing_);

  TestDRValMapCanvas_.Write();
  TestDRCanvas_.Write();
  TauMuBranchingRatioCanvas_.Write();
  TauPtCanvas_.Write();
  NPassingCanvas_.Write();

  out_->Write();
  out_->Close();


}

// ------------ method called when starting to processes a run  ------------
template<class T>
bool CustomTauSelector<T>::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}
 
// ------------ method called when ending the processing of a run  ------------
template<class T>
bool CustomTauSelector<T>::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
template<class T>
bool CustomTauSelector<T>::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
template<class T>
bool CustomTauSelector<T>::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
template<class T>
void CustomTauSelector<T>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
typedef CustomTauSelector<reco::Muon> CustomTauSepFromMuonSelector;
typedef CustomTauSelector<reco::Photon> CustomTauSepFromPhotonSelector;
DEFINE_FWK_MODULE(CustomTauSepFromMuonSelector);
DEFINE_FWK_MODULE(CustomTauSepFromPhotonSelector);
