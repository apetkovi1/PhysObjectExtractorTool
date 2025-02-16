#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <vector>
#include <TChain.h>
#include "math.h"
#include "DataFormats/Math/interface/deltaR.h"

//Compile with:
//g++ -o test MatchingAnalysis.cc $(root-config --cflags --libs)


void BestMatch(float object_eta, float object_phi, int object_pdgId,Long64_t event){

  using namespace std;
   
  TFile* infile = new TFile("myoutput.root", "READ");
  
  TTree* tevent = (TTree*)infile->Get("myevents/Events");
  TTree* tgenparticles = (TTree*)infile->Get("mygenparticle/Events");

  float minDeltaR = 999.0;
  float r;
  Int_t numGenPart=0;
  Int_t j_matched;
  vector<float>*  GenPart_eta=0;
  vector<float>*  GenPart_phi=0;
  vector<int> *GenPart_status=0;
  vector<float> *GenPart_pt=0;
  vector<float> *GenPart_mass=0;
  vector<int> *GenPart_pdgId=0;
    
  tgenparticles->SetBranchAddress("numGenPart",&numGenPart);
  tgenparticles->SetBranchAddress("GenPart_eta",&GenPart_eta);
  tgenparticles->SetBranchAddress("GenPart_phi",&GenPart_phi);
  tgenparticles->SetBranchAddress("GenPart_status",&GenPart_status);
  tgenparticles->SetBranchAddress("GenPart_pdgId",&GenPart_pdgId);
  tgenparticles->SetBranchAddress("GenPart_mass",&GenPart_mass);
  tgenparticles->SetBranchAddress("GenPart_pt",&GenPart_pt);
   
  tevent->AddFriend(tgenparticles);
   
  tevent->GetEntry(event);
  for (Int_t j=0; j<numGenPart;++j)
    {
       r=deltaR(GenPart_eta->at(j),GenPart_phi->at(j),object_eta,object_phi);
       if (r < minDeltaR && abs(GenPart_pdgId->at(j))==object_pdgId) 
       {
        minDeltaR = r;
        j_matched=j;
       }
    }    
  
  cout<<"Matched generator particle has index "<<j_matched+1<<" and the following properties:"<<endl;
  cout<<"pdg ID="<<GenPart_pdgId->at(j_matched)<<endl;
  cout<<"mass="<<GenPart_mass->at(j_matched)<<" "<<"Gev"<<endl;
  cout<<"pt="<<GenPart_pt->at(j_matched)<<" "<<"Gev"<<endl;
  cout<<"eta="<<GenPart_eta->at(j_matched)<<endl;
  cout<<"status="<<GenPart_status->at(j_matched)<<endl;
  cout<<"phi="<<GenPart_phi->at(j_matched)<<endl<<endl;
  delete infile;
  
}


int main (){
  
  using namespace std;
  
  std::vector<float>* electron_eta=0;
  std::vector<float>* electron_phi=0;
  Int_t numelectron=0;
  
  std::vector<float>* muon_eta=0;
  std::vector<float>* muon_phi=0;
  Int_t nummuon=0;
  
  std::vector<float>* tau_eta=0;
  std::vector<float>* tau_phi=0;
  Int_t numtau=0;

  std::vector<float>* photon_eta=0;
  std::vector<float>* photon_phi=0;
  Int_t numphoton=0;

  TFile* infile = new TFile("myoutput.root", "READ");
  TTree* tevent = (TTree*)infile->Get("myevents/Events");
  TTree* telectrons = (TTree*)infile->Get("myelectrons/Events");
  TTree* tmuons = (TTree*)infile->Get("mymuons/Events");
  TTree* ttaus = (TTree*)infile->Get("mytaus/Events");
  TTree* tphotons = (TTree*)infile->Get("myphotons/Events");

  telectrons->SetBranchAddress("electron_eta",&electron_eta);
  telectrons->SetBranchAddress("electron_phi",&electron_phi);
  telectrons->SetBranchAddress("numberelectron",&numelectron);

  tmuons->SetBranchAddress("muon_eta",&muon_eta);
  tmuons->SetBranchAddress("muon_phi",&muon_phi);
  tmuons->SetBranchAddress("numbermuon",&nummuon);

  ttaus->SetBranchAddress("tau_eta",&tau_eta);
  ttaus->SetBranchAddress("tau_phi",&tau_phi);
  ttaus->SetBranchAddress("numbertau",&numtau);

  tphotons->SetBranchAddress("photon_eta",&photon_eta);
  tphotons->SetBranchAddress("photon_phi",&photon_phi);
  tphotons->SetBranchAddress("numberphoton",&numphoton);

  tevent->AddFriend(telectrons);
  Long64_t nentries = tevent->GetEntries(); 
  for (Long64_t i=0;i<nentries;i++) 
  {  
    tevent->GetEntry(i);
    for (Int_t j=0; j<numelectron;++j)
    {
      cout<<"Finding best match for "<<j+1<<". electron"<<" in "<<i+1<<". event"<<endl;
      BestMatch(electron_eta->at(j),electron_phi->at(j),11,i);
    }  
  }

  tevent->AddFriend(tmuons);
  nentries = tevent->GetEntries();
  for (Long64_t i=0;i<nentries;i++) 
  {  
    tevent->GetEntry(i);
    for (Int_t j=0; j<nummuon;++j)
    {
      cout<<"Finding best match for "<<j+1<<". muon"<<" in "<<i+1<<". event"<<endl;
      BestMatch(muon_eta->at(j),muon_phi->at(j),13,i);
    }    
  }
  
  tevent->AddFriend(ttaus);
  nentries = tevent->GetEntries(); 
  for (Long64_t i=0;i<nentries;i++) 
  {  
    tevent->GetEntry(i);
    for (Int_t j=0; j<numtau;++j)
    {
      cout<<"Finding best match for "<<j+1<<". tau"<<" in "<<i+1<<". event"<<endl;
      BestMatch(tau_eta->at(j),tau_phi->at(j),15,i);
    }  
  }

  tevent->AddFriend(tphotons);
  nentries = tevent->GetEntries(); 
  for (Long64_t i=0;i<nentries;i++) 
  {  
    tevent->GetEntry(i);
    for (Int_t j=0; j<numphoton;++j)
    {
      cout<<"Finding best match for "<<j+1<<". photon"<<" in "<<i+1<<". event"<<endl;
      BestMatch(photon_eta->at(j),photon_phi->at(j),22,i);
    }  
  }
}
