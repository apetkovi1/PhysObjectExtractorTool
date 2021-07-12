// -*- C++ -*-
//
// Package:    MatchingAnalyzer
// Class:      MatchingAnalyzer
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

//classes to extract Muon information
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

//classes to extract Photon information
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"

//classes to extract Electron information
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

//classes to extract Tau information
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"

//classes to extract Generator particle information
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

//class to extract separation
#include "DataFormats/Math/interface/deltaR.h"

//classes to save data
#include "TTree.h"
#include "TFile.h"
#include<vector>

//matching functions

template <typename T>
void subtractInvisible(T g, reco::Candidate::LorentzVector& p4) {
  auto daughters = (*g).daughterRefVector();
  for (auto d = daughters.begin(); d != daughters.end(); d++) {
    const auto pdgId = (*d)->pdgId();
    if (std::abs(pdgId) == 12 || std::abs(pdgId) == 14 ||
        std::abs(pdgId) == 16 || std::abs(pdgId) == 18) {
      p4 = p4 - (*d)->p4();
    }
    subtractInvisible(*d, p4);
  }
}

template <typename T>
int findBestVisibleMatch(T& gens, reco::Candidate::LorentzVector& p4) {
  float minDeltaR = 999.0;
  int idx = -1;
  for (auto g = gens.begin(); g != gens.end(); g++) {
    auto tmp_p4 = g->p4();
    subtractInvisible(g, tmp_p4);
    const auto tmp = deltaR(tmp_p4, p4);
    if (tmp < minDeltaR) {
      minDeltaR = tmp;
      idx = g - gens.begin();
    }
  }
  return idx;
}

//
// class declaration
//

class MatchingAnalyzer : public edm::EDAnalyzer {
   public:
      explicit MatchingAnalyzer(const edm::ParameterSet&);
      ~MatchingAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------

      TTree *mtree;
      std::vector<int> GenPart_status;
      std::vector<float> GenPart_pt;
      std::vector<float> GenPart_eta;
      std::vector<float> GenPart_mass;
      std::vector<int> GenPart_pdgId;
      std::vector<float> GenPart_phi; 
    
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

MatchingAnalyzer::MatchingAnalyzer(const edm::ParameterSet& iConfig)

{
//now do what ever initialization is needed
	 
	edm::Service<TFileService> fs;
	mtree = fs->make<TTree>("Events", "Events");

	mtree->Branch("GenPart_pt",&GenPart_pt);
    mtree->GetBranch("GenPart_pt")->SetTitle("generator particle transverse momentum");
    mtree->Branch("GenPart_eta",&GenPart_eta);
    mtree->GetBranch("GenPart_eta")->SetTitle("generator particle pseudorapidity");
    mtree->Branch("GenPart_mass",&GenPart_mass);
    mtree->GetBranch("GenPart_mass")->SetTitle("generator particle mass");
    mtree->Branch("GenPart_phi",&GenPart_phi);
    mtree->GetBranch("GenPart_phi")->SetTitle("generator particle azimuthal angle of momentum vector");
    mtree->Branch("GenPart_pdgId",&GenPart_pdgId);
    mtree->GetBranch("GenPart_pdgId")->SetTitle("generator particle PDG id");
    mtree->Branch("GenPart_status",&GenPart_status);
    mtree->GetBranch("GenPart_status")->SetTitle("Particle status. 1=stable");
}

MatchingAnalyzer::~MatchingAnalyzer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called for each event  ------------
void
MatchingAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;
   using namespace reco;

   Handle<GsfElectronCollection> electrons;
   iEvent.getByLabel(InputTag("gsfElectrons"), electrons);

   Handle<MuonCollection> muons;
   iEvent.getByLabel(InputTag("muons"), muons);

   Handle<PFTauCollection> taus;
   iEvent.getByLabel(InputTag("hpsPFTauProducer"), taus);

   Handle<GenParticleCollection> gens;
   iEvent.getByLabel(InputTag("genParticles"), gens);

   Handle<PhotonCollection> photons;
   iEvent.getByLabel(InputTag("photons"), photons);

   GenPart_pt.clear();
   GenPart_eta.clear();
   GenPart_mass.clear();
   GenPart_pdgId.clear();
   GenPart_phi.clear();
   GenPart_status.clear();
   
   std::vector<GenParticle> interestingGenParticles;

   //find interesting gen particles (electrons, muons, taus, photons)
    
    for (auto it = gens->begin(); it != gens->end(); it++) {
      const auto status = it->status();
      const auto pdgId = std::abs(it->pdgId());
      if (status == 1 && pdgId == 13) // muon
      { 
        interestingGenParticles.emplace_back(*it);
      }
      if (status == 1 && pdgId == 11) // electron
      { 
        interestingGenParticles.emplace_back(*it);
      }
      
      if (status == 1 && pdgId == 22) // photon
      { 
        interestingGenParticles.emplace_back(*it);
      }
      
      if (status == 2 && pdgId == 15) // tau
      { 
        interestingGenParticles.emplace_back(*it);
      }
    }

   //match generator particles with electrons
   for (auto p = electrons->begin(); p != electrons->end(); p++) {
      auto p4 = p->p4();
      auto idx = findBestVisibleMatch(interestingGenParticles, p4);
      if (idx != -1) {
        auto g = interestingGenParticles.begin() + idx;
        GenPart_pt.push_back(g->pt());
        GenPart_eta.push_back(g->eta());
        GenPart_mass.push_back(g->mass());
        GenPart_pdgId.push_back(g->pdgId());
        GenPart_phi.push_back(g->phi());
        GenPart_status.push_back(g->status());
      }
   
  }

   //match generator particles with muons
   for (auto p = muons->begin(); p != muons->end(); p++) {
      auto p4 = p->p4();
      auto idx = findBestVisibleMatch(interestingGenParticles, p4);
      if (idx != -1) {
        auto g = interestingGenParticles.begin() + idx;
        GenPart_pt.push_back(g->pt());
        GenPart_eta.push_back(g->eta());
        GenPart_mass.push_back(g->mass());
        GenPart_pdgId.push_back(g->pdgId());
        GenPart_phi.push_back(g->phi());
        GenPart_status.push_back(g->status());
      }
   
  }
  
   //match generator particles with taus
  for (auto p = taus->begin(); p != taus->end(); p++) {
      auto p4 = p->p4();
      auto idx = findBestVisibleMatch(interestingGenParticles, p4);
      if (idx != -1) {
        auto g = interestingGenParticles.begin() + idx;
        GenPart_pt.push_back(g->pt());
        GenPart_eta.push_back(g->eta());
        GenPart_mass.push_back(g->mass());
        GenPart_pdgId.push_back(g->pdgId());
        GenPart_phi.push_back(g->phi());
        GenPart_status.push_back(g->status());
      }
   
  }

  //match generator particles with photons
  for (auto p = photons->begin(); p != photons->end(); p++) {
      auto p4 = p->p4();
      auto idx = findBestVisibleMatch(interestingGenParticles, p4);
      if (idx != -1) {
        auto g = interestingGenParticles.begin() + idx;
        GenPart_pt.push_back(g->pt());
        GenPart_eta.push_back(g->eta());
        GenPart_mass.push_back(g->mass());
        GenPart_pdgId.push_back(g->pdgId());
        GenPart_phi.push_back(g->phi());
        GenPart_status.push_back(g->status());
      }
   
  }
	
  mtree->Fill();
  return;

}

// ------------ method called once each job just before starting event loop  ------------
void
MatchingAnalyzer::beginJob()
{}

// ------------ method called once each job just after ending the event loop  ------------
void
MatchingAnalyzer::endJob()
{}

// ------------ method called when starting to processes a run  ------------
void
MatchingAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{}

// ------------ method called when ending the processing of a run  ------------
void
MatchingAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{}
// ------------ method called when starting to processes a luminosity block  ------------
void
MatchingAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{}

// ------------ method called when ending the processing of a luminosity block  ------------
void
MatchingAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MatchingAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
 //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MatchingAnalyzer);
