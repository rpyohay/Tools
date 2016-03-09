// -*- C++ -*-
//
// Package:    CleanJets
// Class:      CleanJets
// 
/**\class CleanJets CleanJets.cc BoostedTauAnalysis/CleanJets/src/CleanJets.cc

 Description: Removes PF muons from PFJet candidates and reconstructs the jets
              Matches those PF muons to muons from a --> tau --> muon decays
              Studies the kinematic properties of the discarded muons
	      Associates those muons to the jets from which they were removed

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Francesca Ricci-Tam,6 R-025,+41227672274,
//         Created:  Fri Aug 31 13:01:48 CEST 2012
// $Id: CleanJets.cc,v 1.6 2012/12/06 17:44:58 yohay Exp $
//
//


// system include files
#include <memory>
#include <string>
#include <sstream>
#include <fstream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "Tools/Common/interface/Common.h"
#include "Tools/Common/interface/GenTauDecayID.h"
#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/TauReco/interface/PFTauTagInfo.h"
#include "DataFormats/TauReco/interface/PFTauDecayMode.h"
#include "DataFormats/TauReco/interface/PFTauDecayModeAssociation.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TLorentzVector.h"
#include "TGraphAsymmErrors.h"
#include "TMath.h"

using namespace edm;
using namespace reco;
using namespace std;


//
// class declaration
//

class CleanJets : public edm::EDProducer {
   public:
      explicit CleanJets(const edm::ParameterSet&);
      ~CleanJets();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------

      // source of the jets to be cleaned of muons
      edm::EDGetTokenT<PFJetCollection> jetSrc_;

      // source of muons that, if found within jet, should be removed
      edm::EDGetTokenT<MuonRefVector> muonSrc_;

      // source of PF candidates
      edm::EDGetTokenT<PFCandidateCollection> PFCandSrc_;

      //input, output
      TFile* out_;
      std::string outFileName_;
      edm::ParameterSet* cfg_;
      bool debug_;

      //Histograms
      TH1F* NMu_;
      TH1F* MuEnergy_;
      TH2F* DiffJetEvsMuonE_;
      TH2F* JetEtavsMuEta_;
      TH2F* JetConstBeforevsAfter_;
      TH1F* DiffConstituents_;
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
CleanJets::CleanJets(const edm::ParameterSet& iConfig):
  jetSrc_(consumes<PFJetCollection>(iConfig.getParameter<edm::InputTag>("jetSrc"))),
  muonSrc_(consumes<MuonRefVector>(iConfig.getParameter<edm::InputTag>("muonSrc"))),
  PFCandSrc_(consumes<PFCandidateCollection>(iConfig.getParameter<edm::InputTag>("PFCandSrc"))),
  outFileName_(iConfig.getParameter<std::string>("outFileName"))
{
  cfg_ = const_cast<edm::ParameterSet*>(&iConfig);
  debug_ = true; //set to true if you want to draw validation histograms

  //register your products
  produces<PFJetCollection>( "ak4PFJetsNoMu" );
  produces<edm::ValueMap<bool> >("valMap" );
  produces<edm::ValueMap<MuonRefVector> >( "muonValMap" );
  produces<edm::ValueMap<PFJetRef> >( "jetValMap" );
  produces<PFCandidateCollection>();

}


CleanJets::~CleanJets()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
CleanJets::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   Handle<PFJetCollection> PFJets;
   iEvent.getByToken(jetSrc_, PFJets);
   auto_ptr<PFJetCollection> SetOfJets( new PFJetCollection );

   Handle<MuonRefVector> muons;
   iEvent.getByToken(muonSrc_, muons);

   Handle<PFCandidateCollection> PFCands;
   iEvent.getByToken(PFCandSrc_, PFCands);
   auto_ptr<PFCandidateCollection> PFCandsExcludingSoftMuons(new PFCandidateCollection);

   //fill an STL container with muon ref keys
   std::vector<unsigned int> muonRefKeys;
   if (muons.isValid()) 
   {
     for (MuonRefVector::const_iterator iMuon = muons->begin(); iMuon != muons->end(); ++iMuon) 
       muonRefKeys.push_back(iMuon->key());
   }//if muons.isValid()

   //vector of bools holding the signal muon tag decision for each jet
   std::vector<bool> muonTagDecisions;

   //map between new jet and refs to muons in the original collection that were removed
   std::vector<MuonRefVector> removedMuonMap;

   //map between new jet and ref to original jet
   std::vector<PFJetRef> oldJets;

   std::vector<reco::PFJet> pfjetVector;
   for (PFJetCollection::const_iterator j = PFJets->begin(); j != PFJets->end(); ++j)
     pfjetVector.push_back(*j);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////// NOW DO THE JET-CLEANING /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   int PFCandMuons = 0;
   
   for(reco::PFJetCollection::const_iterator iJet = PFJets->begin(); iJet != PFJets->end(); ++iJet)
   { // loop over jets
       std::vector<reco::PFCandidatePtr> JetPFCands = iJet->getPFConstituents();
       reco::PFJet::Specific specs = iJet->getSpecific();
       double jetEBefore = iJet->energy(), muEta = 1000, muE = -10;
       unsigned int nConstBefore = JetPFCands.size();
       math::XYZTLorentzVector pfmomentum;
       std::vector<edm::Ptr<Candidate> > jetConstituents;
       jetConstituents.clear();

       //flag indicating whether >=0 muons were tagged for removal
       bool taggedMuonForRemoval = false;

       //vector of removed muons for this jet
       MuonRefVector removedMuons;

       for(std::vector<edm::Ptr<reco::PFCandidate> >::iterator i = JetPFCands.begin(); i != JetPFCands.end(); ++i)
       { // loop over PF candidates
	 reco::PFCandidate pfCand = *i;
	 
	 /* Is the PF Candidate a muon? */
	 if (pfCand.particleId() == 3) //Reference: https://cmssdt.cern.ch/SDT/doxygen/CMSSW_7_1_17/doc/html/d8/d17/PFCandidate_8h_source.html
	 { // if it's a muon
	   // get the ref to the corresponding muon
	   // and count one more PF muon
	   reco::MuonRef theRecoMuon = pfCand.muonRef();
	   PFCandMuons += 1;

           //does this muon pass the desired muon ID?
 	   std::vector<unsigned int>::const_iterator iSoftMuon = std::find(muonRefKeys.begin(), muonRefKeys.end(), theRecoMuon.key());
	
           /*if we're not requiring gen matching but instead looking for muons in jets that pass the desired muon ID...*/	       
           if (iSoftMuon != muonRefKeys.end()) 
    	   {
             muE = pfCand.p4().e();
             muEta = pfCand.p4().eta();
             specs.mMuonEnergy -= pfCand.p4().e();
             specs.mMuonMultiplicity -= 1;
             specs.mChargedMuEnergy -= pfCand.p4().e();
             specs.mChargedMultiplicity -= 1;
             //save tag decision for this muon
             taggedMuonForRemoval = true;

             /*add this muon ref to the vector of removed muons for this jet
                iSoftMuon - muonRefKeys.begin() is the index into muonRefKeys of the soft muon
                since muonRefKeys was filled in order of muons, it is also the index into 
                muons of the soft muon*/
	     removedMuons.push_back(muons->at(iSoftMuon - muonRefKeys.begin()));
	   }//if iSoftMuon
	   else
	   {//this muon doesn't pass the soft ID, so keep it in the jet
	     pfmomentum += pfCand.p4(); // total p4()
	     jetConstituents.push_back((*i));
	   }//else
	 }//if pfCand.particleId() == 3
	 else // if it's not a muon
	 { // get p4 and constituents
	   pfmomentum += pfCand.p4(); // total p4()
 	   jetConstituents.push_back((*i));
	 } //get p4 and constituents
       } // loop over PF candidates

       ////// Build a new jet without the muon /////////////

       PFJet muonfreePFJet(pfmomentum, specs, jetConstituents);
       SetOfJets->push_back( muonfreePFJet );

       //if at least 1 muon was tagged for removal, save a positive muon tag decision for this jet
       muonTagDecisions.push_back(taggedMuonForRemoval);

       //save the ref vector of removed muons
       removedMuonMap.push_back(removedMuons);

       //ref to this (old) jet
       oldJets.push_back(PFJetRef(PFJets, iJet - PFJets->begin()));
      
       NMu_->Fill(PFCandMuons );
       if (taggedMuonForRemoval)
       {
         JetConstBeforevsAfter_->Fill(nConstBefore, jetConstituents.size());
         DiffJetEvsMuonE_->Fill(jetEBefore - muonfreePFJet.energy(), muE );
         JetEtavsMuEta_->Fill(iJet->eta(), muEta );
	 MuEnergy_->Fill(muE );
         //std::cout << "\tnConst: Before= " << nConstBefore << " \tjetConstituents.size= " << jetConstituents.size() << std::endl;
         //std::cout << "\tJetE  : Before= " << jetEBefore << " \tmuonfreePFJet= " << muonfreePFJet.energy() << " \tMuon E= " << muE << std::endl;
         //std::cout << "\tMuEta= " << muEta << " \tetiJet->eta()= " << iJet->eta() << std::endl;
         DiffConstituents_->Fill( nConstBefore - jetConstituents.size());
       }

   } // loop over jets
   
   //fill an STL container of keys of removed muons
   std::vector<unsigned int> removedMuRefKeys;
   for (std::vector<MuonRefVector>::const_iterator iJet = removedMuonMap.begin(); iJet != removedMuonMap.end(); ++iJet) 
   {
     for (MuonRefVector::const_iterator iRemovedMu = iJet->begin(); iRemovedMu != iJet->end(); ++iRemovedMu) 
       removedMuRefKeys.push_back(iRemovedMu->key()); 
   }//for iJet

   /*build a collection of PF candidates excluding soft muons
     we will still tag the jet as signal-like by the presence of a soft muon IN the jet, but this 
     ensures that such jets also cannot have the soft muon enter the isolation candidate 
     collection
     right now only remove muons that are inside jets; later expand to muons within dR = X.X of 
     the jets*/
   for (PFCandidateCollection::const_iterator iPFCand = PFCands->begin(); iPFCand != PFCands->end(); ++iPFCand) 
   {
     MuonRef removedMuRef = iPFCand->muonRef();
     if ((removedMuRef.isNonnull() && (std::find(removedMuRefKeys.begin(), removedMuRefKeys.end(), removedMuRef.key()) == removedMuRefKeys.end())) || removedMuRef.isNull()) 
       PFCandsExcludingSoftMuons->push_back(*iPFCand);
   }//for iPFCand

   const OrphanHandle<PFJetCollection> cleanedJetsRefProd = iEvent.put( SetOfJets, "ak4PFJetsNoMu"  );

   //fill the value map of muon tag decision for each cleaned jet
   std::auto_ptr<edm::ValueMap<bool> > valMap(new edm::ValueMap<bool>());
   edm::ValueMap<bool>::Filler filler(*valMap);
   filler.insert(cleanedJetsRefProd, muonTagDecisions.begin(), muonTagDecisions.end());
   filler.fill();
   iEvent.put(valMap, "valMap" );

   //fill the value map of removed muon refs for each cleaned jet
   std::auto_ptr<edm::ValueMap<MuonRefVector> > muonValMap(new edm::ValueMap<MuonRefVector>());
   edm::ValueMap<MuonRefVector>::Filler muonFiller(*muonValMap);
   muonFiller.insert(cleanedJetsRefProd, removedMuonMap.begin(), removedMuonMap.end());
   muonFiller.fill();
   iEvent.put(muonValMap, "muonValMap" );

   //fill the value map of old jet refs for each cleaned jet
   std::auto_ptr<edm::ValueMap<PFJetRef> > jetValMap(new edm::ValueMap<PFJetRef>());
   edm::ValueMap<PFJetRef>::Filler jetFiller(*jetValMap);
   jetFiller.insert(cleanedJetsRefProd, oldJets.begin(), oldJets.end());
   jetFiller.fill();
   iEvent.put(jetValMap, "jetValMap" );

   //put the soft-muon-free PF cands into the event
   iEvent.put(PFCandsExcludingSoftMuons);

}

// ------------ method called once each job just before starting event loop  ------------
void 
CleanJets::beginJob()
{
  if (debug_) out_ = new TFile(outFileName_.c_str(), "RECREATE");

  NMu_= new TH1F("NMu", ";Number of muons;", 7, -.5, 6.5);
  MuEnergy_= new TH1F("MuEnergy", ";Muon p_{T} (GeV);", 50, -.5, 50);
  DiffJetEvsMuonE_ = new TH2F("DiffJetEvsMuonE", 
			      ";E_{uncleaned jet} - E_{cleaned jet} (GeV);E_{#mu} (GeV)", 
			      50, 0, 50, 50, 0, 50);
  JetEtavsMuEta_ = new TH2F("JetEtavsMuEta", ";#eta_{jet};#eta_{#mu}", 20, -5, 5, 20, -5, 5);
  JetConstBeforevsAfter_ = 
    new TH2F("JetConstBeforevsAfter", 
	     ";No. uncleaned jet constituents;No. cleaned jet constituents", 20, 0, 20, 20, 0, 20);
  DiffConstituents_= new TH1F("DiffConstituents", ";N_{const}(uncleaned) - N_{const}(cleaned);", 
			      8, -1.5, 6.5);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
CleanJets::endJob()
{

  if (debug_) {

    out_->cd();

    TCanvas NMuCanvas_("NMuCanvas","",600,600);
    TCanvas MuEnergyCanvas_("MuEnergyCanvas","",600,600);
    TCanvas DiffJetEvsMuonECanvas_("DiffJetEvsMuonECanvas","",600,600);
    TCanvas JetEtavsMuEtaCanvas_("JetEtavsMuEtaCanvas","",600,600);
    TCanvas JetConstBeforevsAfterCanvas_("JetConstBeforevsAfterCanvas","",600,600);
    TCanvas DiffConstituentsCanvas_("DiffConstituentsCanvas","",600,600);

    Common::draw1DHistograms(NMuCanvas_, NMu_);
    Common::draw1DHistograms(MuEnergyCanvas_, MuEnergy_);
    Common::draw1DHistograms(DiffConstituentsCanvas_, DiffConstituents_);

    Common::draw2DHistograms(DiffJetEvsMuonECanvas_, DiffJetEvsMuonE_);
    Common::draw2DHistograms(JetEtavsMuEtaCanvas_, JetEtavsMuEta_);
    Common::draw2DHistograms(JetConstBeforevsAfterCanvas_, JetConstBeforevsAfter_);

    NMuCanvas_.Write();
    MuEnergyCanvas_.Write();
    DiffJetEvsMuonECanvas_.Write();
    JetEtavsMuEtaCanvas_.Write();
    JetConstBeforevsAfterCanvas_.Write();
    DiffConstituentsCanvas_.Write();
    out_->Write();
    out_->Close();

  }

}

// ------------ method called when starting to processes a run  ------------
void 
CleanJets::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
CleanJets::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
CleanJets::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
CleanJets::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
CleanJets::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(CleanJets);
