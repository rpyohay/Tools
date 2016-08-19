// -*- C++ -*-
//
// Package:    Tools/GenMatchedRecoObjectProducer
// Class:      GenObjectAnalyzer
// 
/**\class GenObjectAnalyzer GenObjectAnalyzer.cc Tools/GenMatchedRecoObjectProducer/src/GenObjectAnalyzer.cc

 Description: debug tau gen tools and MC generation features

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Rachel Yohay
//         Created:  Tue, 12 Apr 2016 11:20:53 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "Tools/Common/interface/GenTauDecayID.h"

#include "TH1D.h"
#include "TH2D.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "Tools/Common/interface/Common.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class GenObjectAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit GenObjectAnalyzer(const edm::ParameterSet&);
      ~GenObjectAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------

  edm::EDGetTokenT<reco::GenParticleCollection> genParticleTag_;
  edm::EDGetTokenT<reco::GenParticleRefVector> selectedGenTauTag_;
  edm::EDGetTokenT<reco::GenParticleRefVector> selectedGenMuonTag_;
  edm::EDGetTokenT<reco::GenParticleRefVector> selectedGenPseudoscalarTag_;

  edm::ParameterSet genTauDecayIDPSet_;

  edm::EDGetTokenT<edm::TriggerResults> triggerBitsTag_;
  edm::EDGetTokenT<trigger::TriggerEvent> triggerEventTag_;

  std::map<std::string, TH1D*> hists1D_;
  std::map<std::string, TH2D*> hists2D_;

  unsigned int nTriggerableEvts_;

  double diTauDRMax_;
  double diMuonDRMax_;
  double mu1PTMax_;
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
GenObjectAnalyzer::GenObjectAnalyzer(const edm::ParameterSet& iConfig) :

  genParticleTag_(consumes<reco::GenParticleCollection>
		  (iConfig.getParameter<edm::InputTag>("genParticleTag"))),
  selectedGenTauTag_(consumes<reco::GenParticleRefVector>
		     (iConfig.getParameter<edm::InputTag>("selectedGenTauTag"))),
  selectedGenMuonTag_(consumes<reco::GenParticleRefVector>
		     (iConfig.getParameter<edm::InputTag>("selectedGenMuonTag"))),
  selectedGenPseudoscalarTag_(consumes<reco::GenParticleRefVector>
			      (iConfig.getParameter<edm::InputTag>("selectedGenPseudoscalarTag"))),

  genTauDecayIDPSet_(iConfig.getParameter<edm::ParameterSet>("genTauDecayIDPSet")),

  triggerBitsTag_(consumes<edm::TriggerResults>
		  (iConfig.getParameter<edm::InputTag>("triggerBitsTag"))),
  triggerEventTag_(consumes<trigger::TriggerEvent>
		   (iConfig.getParameter<edm::InputTag>("triggerEventTag"))),

  nTriggerableEvts_(0),

  diTauDRMax_(0),
  diMuonDRMax_(0),
  mu1PTMax_(0)
{
   //now do what ever initialization is needed
   usesResource("TFileService");

}


GenObjectAnalyzer::~GenObjectAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
GenObjectAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // Get collections from event //

  edm::Handle<reco::GenParticleCollection> pGenParticles;
  iEvent.getByToken(genParticleTag_, pGenParticles);

  edm::Handle<reco::GenParticleRefVector> pSelectedGenTaus;
  iEvent.getByToken(selectedGenTauTag_, pSelectedGenTaus);

  edm::Handle<reco::GenParticleRefVector> pSelectedGenMuons;
  iEvent.getByToken(selectedGenMuonTag_, pSelectedGenMuons);

  edm::Handle<reco::GenParticleRefVector> pSelectedGenPseudoscalars;
  iEvent.getByToken(selectedGenPseudoscalarTag_, pSelectedGenPseudoscalars);

  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBitsTag_, triggerBits);

  edm::Handle<trigger::TriggerEvent> triggerEvent;
  iEvent.getByToken(triggerEventTag_, triggerEvent);

  // Analyze decay properties and kinematics //

  reco::LeafCandidate::LorentzVector fourBodyP4;
  reco::LeafCandidate::LorentzVector diTauP4;
  reco::LeafCandidate::LorentzVector diMuonP4;
  std::vector<reco::LeafCandidate> taus;
  // std::vector<reco::GenParticleRef> taus;
  unsigned int nLeptons = 0;
  unsigned int nTaus = 0;
  unsigned int nMuons = 0;

  for (reco::GenParticleRefVector::const_iterator iSelectedGenTau = pSelectedGenTaus->begin(); 
       iSelectedGenTau != pSelectedGenTaus->end(); ++iSelectedGenTau) {
    try {

      GenTauDecayID tauDecay(genTauDecayIDPSet_, pGenParticles, iSelectedGenTau->key());
      const std::pair<reco::PFTau::hadronicDecayMode, GenTauDecayID::DecayType> decayType = 
	tauDecay.tauDecayType(false, false);

      hists1D_["genDecayMode"]->Fill(decayType.first);

      reco::LeafCandidate::LorentzVector visibleGenP4 = tauDecay.getVisibleTauP4();
      fourBodyP4+=visibleGenP4;
      diTauP4+=visibleGenP4;
      taus.push_back(reco::LeafCandidate(0.0, visibleGenP4));

    }
    catch (std::string& ex) { throw cms::Exception("GenObjectAnalyzer") << ex; }

    // fourBodyP4+=(*iSelectedGenTau)->p4();
    // diTauP4+=(*iSelectedGenTau)->p4();
    // taus.push_back(*iSelectedGenTau);
    ++nLeptons;
    ++nTaus;
  }

  std::vector<reco::GenParticleRef> muons;
  std::vector<reco::GenParticle*> muonPtrs;
  double mu1PT = 0.0;
  double mu1Eta = 0.0;
  double mu1Phi = 0.0;

  for (reco::GenParticleRefVector::const_iterator iSelectedGenMuon = pSelectedGenMuons->begin(); 
       iSelectedGenMuon != pSelectedGenMuons->end(); ++iSelectedGenMuon) {

    fourBodyP4+=(*iSelectedGenMuon)->p4();
    diMuonP4+=(*iSelectedGenMuon)->p4();
    muons.push_back(*iSelectedGenMuon);
    muonPtrs.push_back(const_cast<reco::GenParticle*>(iSelectedGenMuon->get()));
    ++nLeptons;
    ++nMuons;

    if ((*iSelectedGenMuon)->pt() > mu1PT) {
      mu1PT = (*iSelectedGenMuon)->pt();
      mu1Eta = (*iSelectedGenMuon)->eta();
      mu1Phi = (*iSelectedGenMuon)->phi();
    }

  }

  for (reco::GenParticleRefVector::const_iterator iSelectedGenPseudoscalar = 
	 pSelectedGenPseudoscalars->begin(); iSelectedGenPseudoscalar != 
	 pSelectedGenPseudoscalars->end(); ++iSelectedGenPseudoscalar) {

    hists1D_["genPseudoscalarEta"]->Fill((*iSelectedGenPseudoscalar)->eta());

    if (nTaus == 0) {
      for (reco::GenParticleRefVector::const_iterator iDaughter = 
	     (*iSelectedGenPseudoscalar)->daughterRefVector().begin(); iDaughter != 
	     (*iSelectedGenPseudoscalar)->daughterRefVector().end(); ++iDaughter) {
	std::cout << "Daughter PDG ID: " << (*iDaughter)->pdgId() << std::endl;
      }

    }
  }

  if (nLeptons == 4) hists1D_["gen4BodyMass"]->Fill(fourBodyP4.M());
  if (nTaus == 2) {
    hists1D_["genDiTauMass"]->Fill(diTauP4.M());
    // hists1D_["genDiTauDRBeforeHLT"]->Fill(std::sqrt(deltaR2(*(taus[0]), *(taus[1]))));
    const double diTauDR = std::sqrt(deltaR2(taus[0], taus[1]));
    hists1D_["genDiTauDRBeforeHLT"]->Fill(diTauDR);
    if (diTauDR > diTauDRMax_) diTauDRMax_ = diTauDR;
  }
  if (nMuons == 2) {
    hists1D_["genDiMuonMass"]->Fill(diMuonP4.M());
    const double diMuonDR = std::sqrt(deltaR2(*(muons[0]), *(muons[1])));
    hists1D_["genDiMuonDRBeforeHLT"]->Fill(diMuonDR);
    if (diMuonDR > diMuonDRMax_) diMuonDRMax_ = diMuonDR;
  }
  hists1D_["nLeptons"]->Fill(nLeptons);
  hists1D_["nTaus"]->Fill(nTaus);
  hists1D_["nMuons"]->Fill(nMuons);

  hists1D_["genMu1PT"]->Fill(mu1PT);
  hists1D_["genMu1Eta"]->Fill(mu1Eta);
  hists2D_["genMu1AbsEtaVsPT"]->Fill(mu1PT, fabs(mu1Eta));

  if (mu1PT > mu1PTMax_) mu1PTMax_ = mu1PT;

  for (reco::GenParticleCollection::const_iterator iGenParticle = pGenParticles->begin(); 
       iGenParticle != pGenParticles->end(); ++iGenParticle) {

    if (iGenParticle->pdgId() == GenTauDecayID::HPDGID) {
      hists1D_["genHiggsEta"]->Fill(iGenParticle->eta());
    }

  }

  // Analyze trigger decision //

  bool passMu45Eta2p1 = false;
  bool passMu17Mu8SameSignDZ = false;
  bool passTripleMu12105 = false;
  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  for (unsigned int iBit = 0; iBit < triggerBits->size(); ++iBit) {
    if (names.triggerName(iBit).find("HLT_Mu45_eta2p1") != std::string::npos) {
      passMu45Eta2p1 = triggerBits->accept(iBit);
    }
    if (names.triggerName(iBit).find("HLT_Mu17_Mu8_SameSign_DZ") != std::string::npos) {
      passMu17Mu8SameSignDZ = triggerBits->accept(iBit);
    }
    if (names.triggerName(iBit).find("HLT_TripleMu_12_10_5") != std::string::npos) {
      passTripleMu12105 = triggerBits->accept(iBit);
    }
  }
  if (passMu45Eta2p1) hists1D_["HLTDecisionByPath"]->Fill(0);
  else if (passMu17Mu8SameSignDZ) hists1D_["HLTDecisionByPath"]->Fill(1);
  else if (passTripleMu12105) hists1D_["HLTDecisionByPath"]->Fill(2);
  else hists1D_["HLTDecisionByPath"]->Fill(3);

  if ((mu1PT > 45.0/*GeV*/) && (fabs(mu1Eta) < 2.1)) {
    ++nTriggerableEvts_;

    if (nTaus == 2) {
      // hists1D_["genDiTauDRAfterHLT"]->Fill(std::sqrt(deltaR2(*(taus[0]), *(taus[1]))));
      hists1D_["genDiTauDRAfterHLT"]->Fill(std::sqrt(deltaR2(taus[0], taus[1])));
    }
    if (nMuons == 2) {
      hists1D_["genDiMuonDRAfterHLT"]->Fill(std::sqrt(deltaR2(*(muons[0]), *(muons[1]))));
    }
    bool HLTDecision = false;
    const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
    for (unsigned int iBit = 0; iBit < triggerBits->size(); ++iBit) {
      if (names.triggerName(iBit).find("HLT_Mu45_eta2p1") != std::string::npos) {
    	HLTDecision = triggerBits->accept(iBit);
    	hists1D_["HLTDecision"]->Fill(HLTDecision);
      }
    }

    double mu2PT = 0.0;
    double mu2Eta = 0.0;
    double mu2Phi = 0.0;
    if (mu1PT == muons[0]->pt()) {
      mu2PT = muons[1]->pt();
      mu2Eta = muons[1]->eta();
      mu2Phi = muons[1]->phi();
    }
    else {
      mu2PT = muons[0]->pt();
      mu2Eta = muons[0]->eta();
      mu2Phi = muons[0]->phi();
    }
    std::cerr << "\nmu1\n";
    std::cerr << "pT (GeV): " << mu1PT << " | eta: " << mu1Eta << " | phi: " << mu1Phi << std::endl;
    std::cerr << "mu2\n";
    std::cerr << "pT (GeV): " << mu2PT << " | eta: " << mu2Eta << " | phi: " << mu2Phi << std::endl;
    const trigger::TriggerObjectCollection& triggerObjects(triggerEvent->getObjects());
    const trigger::size_type 
      iFilter(triggerEvent->filterIndex
    	      (edm::InputTag("hltL3fL1sMu16orMu25L1f0L2f10QL3Filtered45e2p1Q", "", "HLT")));
    // const trigger::size_type 
    //   iFilter(triggerEvent->filterIndex
    // 	      (edm::InputTag("hltL2fL1sMu16orMu25L1f0L2Filtered10Q", "", "HLT")));
    unsigned int nMu1MatchedTrgObjs = 0;
    unsigned int nMu2MatchedTrgObjs = 0;
    if (iFilter < triggerEvent->sizeFilters()) {
      const trigger::Keys& HLTKeys(triggerEvent->filterKeys(iFilter));
      std::cerr << "HLT_Mu45_eta2p1 trigger objects\n";
      for (trigger::Keys::const_iterator iHLTKey = HLTKeys.begin(); iHLTKey != HLTKeys.end(); 
    	   ++iHLTKey) {
    	const trigger::TriggerObject& triggerObject(triggerObjects[*iHLTKey]);
    	std::cerr << "ID: " << triggerObject.id() << " | pT (GeV): " << triggerObject.pt();
    	std::cerr << " | eta: " << triggerObject.eta() << " | phi: " << triggerObject.phi();
    	std::cerr << std::endl;
    	int nearestGenMuKey = -1;
    	const reco::Particle triggerObjectParticle(triggerObject.particle());
    	const reco::GenParticle* nearestGenMu = 
    	  Common::nearestObject(triggerObjectParticle, muonPtrs, nearestGenMuKey);
    	const double trgObjNearestGenMuonDR = 
    	  std::sqrt(deltaR2(triggerObjectParticle, *nearestGenMu));
    	hists1D_["trgObjNearestGenMuonDR"]->Fill(trgObjNearestGenMuonDR);
    	if (trgObjNearestGenMuonDR < 0.1) {
    	  if (nearestGenMu->pt() == mu1PT) ++nMu1MatchedTrgObjs;
    	  if (nearestGenMu->pt() == mu2PT) ++nMu2MatchedTrgObjs;
    	}
      }
      std::cerr << std::endl;
    }
    else std::cerr << "Event did not pass HLT_Mu45_eta2p1\n\n";
    if (HLTDecision) {
      if ((nMu1MatchedTrgObjs > 0) && (nMu2MatchedTrgObjs > 0)) {
    	hists1D_["HLTDecisionByGenMuon"]->Fill(0);
      }
      if ((nMu1MatchedTrgObjs > 0) && (nMu2MatchedTrgObjs == 0)) {
    	hists1D_["HLTDecisionByGenMuon"]->Fill(1);
      }
      if ((nMu1MatchedTrgObjs == 0) && (nMu2MatchedTrgObjs > 0)) {
    	hists1D_["HLTDecisionByGenMuon"]->Fill(2);
      }
      if ((nMu1MatchedTrgObjs == 0) && (nMu2MatchedTrgObjs == 0)) {
    	hists1D_["HLTDecisionByGenMuon"]->Fill(3);
      }
    }
    else hists1D_["HLTDecisionByGenMuon"]->Fill(4);

  }

}


// ------------ method called once each job just before starting event loop  ------------
void 
GenObjectAnalyzer::beginJob()
{
  edm::Service<TFileService> fileService;

  hists1D_["genDecayMode"] = 
    fileService->make<TH1D>("genDecayMode", "", 
			    reco::PFTau::kThreeProngNPiZero - reco::PFTau::kOneProng0PiZero + 1, 
			    reco::PFTau::kOneProng0PiZero - 0.5, 
			    reco::PFTau::kThreeProngNPiZero + 0.5);
  TAxis* xAxis = hists1D_["genDecayMode"]->GetXaxis();
  xAxis->SetBinLabel(reco::PFTau::kOneProng0PiZero + 1, "1-prong + 0-#pi^{0}");
  xAxis->SetBinLabel(reco::PFTau::kOneProng1PiZero + 1, "1-prong + 1-#pi^{0}");
  xAxis->SetBinLabel(reco::PFTau::kOneProng2PiZero + 1, "1-prong + 2-#pi^{0}");
  xAxis->SetBinLabel(reco::PFTau::kOneProng3PiZero + 1, "1-prong + 3-#pi^{0}");
  xAxis->SetBinLabel(reco::PFTau::kOneProngNPiZero + 1, "1-prong + N-#pi^{0}");
  xAxis->SetBinLabel(reco::PFTau::kTwoProng0PiZero + 1, "2-prong + 0-#pi^{0}");
  xAxis->SetBinLabel(reco::PFTau::kTwoProng1PiZero + 1, "2-prong + 1-#pi^{0}");
  xAxis->SetBinLabel(reco::PFTau::kTwoProng2PiZero + 1, "2-prong + 2-#pi^{0}");
  xAxis->SetBinLabel(reco::PFTau::kTwoProng3PiZero + 1, "2-prong + 3-#pi^{0}");
  xAxis->SetBinLabel(reco::PFTau::kTwoProngNPiZero + 1, "2-prong + N-#pi^{0}");
  xAxis->SetBinLabel(reco::PFTau::kThreeProng0PiZero + 1, "3-prong + 0-#pi^{0}");
  xAxis->SetBinLabel(reco::PFTau::kThreeProng1PiZero + 1, "3-prong + 1-#pi^{0}");
  xAxis->SetBinLabel(reco::PFTau::kThreeProng2PiZero + 1, "3-prong + 2-#pi^{0}");
  xAxis->SetBinLabel(reco::PFTau::kThreeProng3PiZero + 1, "3-prong + 3-#pi^{0}");
  xAxis->SetBinLabel(reco::PFTau::kThreeProngNPiZero + 1, "3-prong + N-#pi^{0}");

  hists1D_["gen4BodyMass"] = 
    fileService->make<TH1D>("gen4BodyMass", ";m_{#mu#mu#tau#tau} (GeV);", 34, 0.0, 750.0);
  hists1D_["genDiTauMass"] = 
    fileService->make<TH1D>("genDiTauMass", ";m_{#tau#tau} (GeV);", 20, 0.0, 10.0);
  hists1D_["genDiMuonMass"] = 
    fileService->make<TH1D>("genDiMuonMass", ";m_{#mu#mu} (GeV);", 20, 8.0, 10.0);
  hists1D_["genDiTauDRBeforeHLT"] = 
    fileService->make<TH1D>("genDiTauDRBeforeHLT", ";#DeltaR(#tau,#tau);", 35, 0.0, 3.5);
  hists1D_["genDiMuonDRBeforeHLT"] = 
    fileService->make<TH1D>("genDiMuonDRBeforeHLT", ";#DeltaR(#mu,#mu);", 35, 0.0, 3.5);
  hists1D_["nLeptons"] = 
    fileService->make<TH1D>("nLeptons", ";No. leptons from H decay;", 5, -0.5, 4.5);
  hists1D_["nTaus"] = fileService->make<TH1D>("nTaus", ";No. taus from H decay;", 3, -0.5, 2.5);
  hists1D_["nMuons"] = fileService->make<TH1D>("nMuons", ";No. muons from H decay;", 3, -0.5, 2.5);
  hists1D_["genPseudoscalarEta"] = 
    fileService->make<TH1D>("genPseudoscalarEta", ";a #eta;", 20, -5.0, 5.0);
  hists1D_["genHiggsEta"] = 
    fileService->make<TH1D>("genHiggsEta", ";a #eta;", 20, -5.0, 5.0);

  hists1D_["genMu1PT"] = 
    fileService->make<TH1D>("genMu1PT", ";#mu_{1} p_{T} (GeV);", 28, 0.0, 560.0);
  hists1D_["genMu1Eta"] = 
    fileService->make<TH1D>("genMu1Eta", ";#mu_{1} #eta;", 20, -5.0, 5.0);
  hists2D_["genMu1AbsEtaVsPT"] = 
    fileService->make<TH2D>("genMu1AbsEtaVsPT", ";#mu_{1} p_{T} (GeV);#mu_{1} |#eta|", 
			    28, 0.0, 560.0, 20, -5.0, 5.0);

  hists1D_["HLTDecisionByPath"] = 
    fileService->make<TH1D>("HLTDecisionByPath", "", 4, -0.5, 3.5);
  xAxis = hists1D_["HLTDecisionByPath"]->GetXaxis();
  xAxis->SetBinLabel(1, "(1)"/*"HLT_Mu45_eta2p1"*/);
  xAxis->SetBinLabel(2, "(2) && !(1)"/*"HLT_Mu17_Mu8_SameSign_DZ && !HLT_Mu45_eta2p1"*/);
  xAxis->SetBinLabel(3, "(3) && !(1) && !(2)"/*"HLT_TripleMu_12_10_5 && !HLT_Mu45_eta2p1 && !HLT_Mu17_Mu8_SameSign_DZ"*/);
  xAxis->SetBinLabel(4, "!(1) && !(2) && !(3)"/*"!HLT_TripleMu_12_10_5 && !HLT_Mu45_eta2p1 && !HLT_Mu17_Mu8_SameSign_DZ"*/);
  hists1D_["genDiTauDRAfterHLT"] = 
    fileService->make<TH1D>("genDiTauDRAfterHLT", ";#DeltaR(#tau,#tau);", 35, 0.0, 3.5);
  hists1D_["genDiMuonDRAfterHLT"] = 
    fileService->make<TH1D>("genDiMuonDRAfterHLT", ";#DeltaR(#mu,#mu);", 35, 0.0, 3.5);
  hists1D_["HLTDecision"] = fileService->make<TH1D>("HLTDecision", "", 2, -0.5, 1.5);
  xAxis = hists1D_["HLTDecision"]->GetXaxis();
  xAxis->SetBinLabel(1, "Fail");
  xAxis->SetBinLabel(2, "Pass");
  hists1D_["trgObjNearestGenMuonDR"] = 
    fileService->make<TH1D>("trgObjNearestGenMuonDR", 
			    ";#DeltaR(trigger object, nearest gen muon);", 25, 0.0, 0.005);
  hists1D_["HLTDecisionByGenMuon"] = 
    fileService->make<TH1D>("HLTDecisionByGenMuon", "", 4, -0.5, 3.5);
  xAxis = hists1D_["HLTDecisionByGenMuon"]->GetXaxis();
  xAxis->SetBinLabel(1, "Both");
  xAxis->SetBinLabel(2, "#mu_{1} only");
  xAxis->SetBinLabel(3, "#mu_{2} only");
  xAxis->SetBinLabel(4, "Neither");

}

// ------------ method called once each job just after ending the event loop  ------------
void 
GenObjectAnalyzer::endJob() 
{
  std::cout << std::endl << "Number of triggerable events: " << nTriggerableEvts_ << std::endl;
  std::cout << std::endl << "Max. di-tau dR: " << diTauDRMax_ << std::endl;
  std::cout << std::endl << "Max. di-muon dR: " << diMuonDRMax_ << std::endl;
  std::cout << std::endl << "Max. mu1 pT (GeV): " << mu1PTMax_ << std::endl;
  std::cout << std::endl;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenObjectAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("genParticleTag", edm::InputTag("genParticles"));
  desc.add<edm::InputTag>("selectedGenTauTag", edm::InputTag("genTauHadATauMuTauHadSelector"));
  desc.add<edm::InputTag>("selectedGenMuonTag", edm::InputTag("genAMuMuSelector"));
  desc.add<edm::InputTag>("selectedGenPseudoscalarTag", edm::InputTag("genAHAASelector"));
  desc.add<edm::InputTag>("triggerBitsTag", edm::InputTag("TriggerResults", "", "HLT"));
  desc.add<edm::InputTag>("triggerEventTag", edm::InputTag("hltTriggerSummaryAOD", "", "HLT"));

  edm::ParameterSetDescription ATauTauPSet;
  ATauTauPSet.add<std::vector<int> >("momPDGID", std::vector<int>(1, GenTauDecayID::APDGID));
  ATauTauPSet.add<double>("chargedHadronPTMin", 0.0);
  ATauTauPSet.add<double>("neutralHadronPTMin", 0.0);
  ATauTauPSet.add<double>("chargedLeptonPTMin", 0.0);
  ATauTauPSet.add<double>("totalPTMin", 0.0);
  desc.add<edm::ParameterSetDescription>("genTauDecayIDPSet", ATauTauPSet);

  descriptions.add("genObjectAnalyzer", desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenObjectAnalyzer);
