// -*- C++ -*-
//
// Package:    Tools/GenMatchedRecoObjectProducer
// Class:      GenObjectAnalyzer
// 
/**\class GenObjectAnalyzer GenObjectAnalyzer.cc Tools/GenMatchedRecoObjectProducer/src/GenObjectAnalyzer.cc

 Description: debug tau gen tools

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

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "Tools/Common/interface/GenTauDecayID.h"

#include "TH1D.h"
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
  edm::EDGetTokenT<reco::GenParticleRefVector> selectedGenParticleTag_;

  edm::ParameterSet genTauDecayIDPSet_;

  std::map<std::string, TH1D*> hists_;
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
  selectedGenParticleTag_(consumes<reco::GenParticleRefVector>
			  (iConfig.getParameter<edm::InputTag>("selectedGenParticleTag"))),

  genTauDecayIDPSet_(iConfig.getParameter<edm::ParameterSet>("genTauDecayIDPSet"))
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
  edm::Handle<reco::GenParticleCollection> pGenParticles;
  iEvent.getByToken(genParticleTag_, pGenParticles);

  edm::Handle<reco::GenParticleRefVector> pSelectedGenParticles;
  iEvent.getByToken(selectedGenParticleTag_, pSelectedGenParticles);

  for (reco::GenParticleRefVector::const_iterator iSelectedGenParticle = 
	 pSelectedGenParticles->begin(); iSelectedGenParticle != pSelectedGenParticles->end(); 
       ++iSelectedGenParticle) {
    try {

      GenTauDecayID tauDecay(genTauDecayIDPSet_, pGenParticles, iSelectedGenParticle->key());
      const std::pair<reco::PFTau::hadronicDecayMode, GenTauDecayID::DecayType> decayType = 
	tauDecay.tauDecayType(false, false);

      hists_["genDecayMode"]->Fill(decayType.first);

    }
    catch (std::string& ex) { throw cms::Exception("GenObjectAnalyzer") << ex; }
  }
}


// ------------ method called once each job just before starting event loop  ------------
void 
GenObjectAnalyzer::beginJob()
{
  edm::Service<TFileService> fileService;

  hists_["genDecayMode"] = 
    fileService->make<TH1D>("genDecayMode", "", 
			    reco::PFTau::kThreeProngNPiZero - reco::PFTau::kOneProng0PiZero + 1, 
			    reco::PFTau::kOneProng0PiZero - 0.5, 
			    reco::PFTau::kThreeProngNPiZero + 0.5);
  TAxis* xAxis = hists_["genDecayMode"]->GetXaxis();
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
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GenObjectAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenObjectAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("genParticleTag", edm::InputTag("genParticles"));
  desc.add<edm::InputTag>("selectedGenParticleTag", edm::InputTag("genTauHadATauMuTauHadSelector"));

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
