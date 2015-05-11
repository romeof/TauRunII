// -*- C++ -*-
//
// Package:    JetTauFakeRateDisc/JetTauFakeRateDisc
// Class:      JetTauFakeRateDisc
// 
/**\class JetTauFakeRateDisc JetTauFakeRateDisc.cc JetTauFakeRateDisc/JetTauFakeRateDisc/plugins/JetTauFakeRateDisc.cc

 Description:
 This class is intended for studying tau properties in order to reduce jet -> tau fake rate

 Implementation:
 The implementation follows some rules. 
 If you modify part of the code, you are kindly invited to be coherent with the style it is written.
 For instance, pay attention to  
  the way you write comments
  the indentation
  the use of methods and functions in the main code
*/
//
// Original Author:  Francesco Romeo
//         Created:  Wed, 01 Apr 2015 15:26:06 GMT
//
//
/////
//   Headers
/////
//system and events
#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TStopwatch.h"
//Gen info
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
//Trigger
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
//Vertex
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
//Candidate
#include "DataFormats/Candidate/interface/CandidateFwd.h"
//Muon
//Electron
//Tau
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
//#include "DataFormats/TauReco/interface/PFTauDecayMode.h"//How Tau decays
#include "DataFormats/TauReco/interface/PFTau.h"
#include "RecoTauTag/RecoTau/interface/HPSPFRecoTauAlgorithm.h"
//Photon
//Jet and Met
#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
//Track builder infos
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexUpdator.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexTrackCompatibilityEstimator.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexSmoother.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
//Math
#include "DataFormats/Math/interface/deltaR.h"
//Store info
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
//User defined classes
#include "TauRunII/JetTauFakeRateDisc/interface/TreeVariables.h"
#include "TauRunII/JetTauFakeRateDisc/interface/HelpingFunctions.h"
#include "TauRunII/JetTauFakeRateDisc/interface/IPvars.h"
/////
//   Namespace
/////
using namespace edm;
using namespace reco;
using namespace std;
/////
//   Class declaration
/////
class JetTauFakeRateDisc : public edm::EDAnalyzer {
 public:
 explicit JetTauFakeRateDisc(const edm::ParameterSet&);
 ~JetTauFakeRateDisc();
 static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
 private:
 virtual void beginJob() override;
 virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
 virtual void endJob() override;
 //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
 //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
 //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
 //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
 //Values from python/ConfFile_cfg.py
 const double tau_min_pt;
 const double tau_max_eta;
 const bool tau_gentauh_match;
 const bool tau_genjet_match;
 // ----------member data ---------------------------
 //Watch time and cpu for the analysis
 TStopwatch* stopwatch;
 int num_evt_tot;
 //Tree
 CTree *tree;
 const edm::Service<TFileService> fs;
};
/////
//   Constructors and destructor
/////
JetTauFakeRateDisc::JetTauFakeRateDisc(const edm::ParameterSet& iConfig):
 //Values from python/ConfFile_cfg.py
 tau_min_pt(iConfig.getUntrackedParameter<double>("tau_min_pt")),
 tau_max_eta(iConfig.getUntrackedParameter<double>("tau_max_eta")),
 tau_gentauh_match(iConfig.getParameter<bool>("tau_gentauh_match")),
 tau_genjet_match(iConfig.getParameter<bool>("tau_genjet_match")),
 //TTree 
 tree(new CTree(fs->make<TTree>("tree", "tree")))
{
 //Now do what ever initialization is needed
 tree->make_branches();
 stopwatch = new TStopwatch();
}
JetTauFakeRateDisc::~JetTauFakeRateDisc()
{
 // do anything here that needs to be done at desctruction time
 // (e.g. close files, deallocate resources etc.)
 delete stopwatch;
}
/////
//   Member functions
/////
// ------------ method called for each event  ------------
void JetTauFakeRateDisc::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
 num_evt_tot++;
 //Vertex collection
 Handle<vector<Vertex> > vertices;
 iEvent.getByLabel("offlinePrimaryVertices", vertices);
 if(vertices->empty()) return; // skip the event if no PV found
 const reco::Vertex &pv = vertices->front();
 //TransientTrackBuilder
 edm::ESHandle<TransientTrackBuilder> ttrkbuilder;
 iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",ttrkbuilder); 
 //Tau collection
 Handle<vector<PFTau> > PFTaus;
 iEvent.getByLabel("hpsPFTauProducer", PFTaus);
 //Order in Pt (Not too elegant sort, improve it)  
 vector<pair<double,int> > TauPtPos;
 for(int taupos=0; taupos<int(PFTaus->size()); taupos++) TauPtPos.push_back(make_pair((*PFTaus)[taupos].pt(),taupos)); 
 sort(TauPtPos.rbegin(),TauPtPos.rend());
 //Access tau info
 for(int taupos=0; taupos<int(PFTaus->size()); taupos++){
  const PFTau& pftau = (*PFTaus)[TauPtPos[taupos].second];
  //Minimum tau requirements
  if(!(pftau.pt()>tau_min_pt && fabs(pftau.eta())<tau_max_eta && pftau.decayMode()!=-1)) continue;
  //Matching with GenTauh (posgenmatchedcand = -1 if there is no matching)
  int posgenmatchedcand = -1;
  if(tau_gentauh_match) posgenmatchedcand = MatchRecoTauhGenTauh(pftau,iEvent);
  if(tau_genjet_match)  posgenmatchedcand = MatchRecoTauhGenJet(pftau,iEvent);
  //const PFJetRef & pfjetref = pftau.jetRef();
  if(posgenmatchedcand==-1) continue;
  //Now get the info you need for the study
  tree->loop_initialize();
  //Refit PV removing pftau ch hads
  //TransientVertex refpv = refitted_vertex(pv,*ttrkbuilder);
  TransientVertex unbpv = unbiased_vertex(pv,*ttrkbuilder,pftau);
  if(!unbpv.isValid()) continue;
  //Kinematic
  tree->tau_pt  = pftau.pt();
  tree->tau_eta = pftau.eta();
  tree->tau_phi = pftau.phi();
  tree->tau_en  = pftau.energy();
  //Charge
  tree->tau_ch  = pftau.charge();
  //Tau IP
  tree->tau_vtxdz   = TMath::Abs(pftau.vertex().z()-pv.position().z());
  tree->tau_vtxdxy  = sqrt(pow(pftau.vertex().x()-pv.position().x(),2)+pow(pftau.vertex().y()-pv.position().y(),2));
  tree->tau_vtxdxyz = sqrt(pow(pftau.vertex().x()-pv.position().x(),2)+pow(pftau.vertex().y()-pv.position().y(),2)+pow(pftau.vertex().z()-pv.position().z(),2));
  //Info of tau constituents (for tau ch had)
  GlobalVector pftaugv(pftau.px(), pftau.py(), pftau.pz());
  const vector<reco::PFCandidatePtr>& sigpfchhadcands = pftau.signalPFChargedHadrCands();
  for(uint p=0; p<sigpfchhadcands.size(); p++){
   PFCandidatePtr cand = sigpfchhadcands[p];
   const Track* candtrk = cand->bestTrack();
   if(!candtrk) continue;
   if(!is_goodtrk(candtrk,pv)) continue;
   //Kinematic
   tree->pftauchhads_pt[p]  = cand->pt();
   tree->pftauchhads_eta[p] = cand->eta();
   tree->pftauchhads_phi[p] = cand->phi();
   tree->pftauchhads_en[p]  = cand->energy();
   //Charge
   tree->pftauchhads_ch[p]  = cand->charge();
   //IP and signed IP
   double pftauchhads_IP3D_val = 0;
   double pftauchhads_IP3D_err = 0;
   double pftauchhads_IP3D_sig = 0;
   double pftauchhads_IP2D_val = 0;
   double pftauchhads_IP2D_err = 0;
   double pftauchhads_IP2D_sig = 0;
   double pftauchhads_sIP3D_val = 0;
   double pftauchhads_sIP3D_err = 0;
   double pftauchhads_sIP3D_sig = 0;
   double pftauchhads_sIP2D_val = 0;
   double pftauchhads_sIP2D_err = 0;
   double pftauchhads_sIP2D_sig = 0;
   IP3D2D( candtrk,*ttrkbuilder,unbpv,pftaugv, pftauchhads_IP3D_val,pftauchhads_IP3D_err,pftauchhads_IP3D_sig,  pftauchhads_IP2D_val,pftauchhads_IP2D_err,pftauchhads_IP2D_sig, pftauchhads_sIP3D_val,pftauchhads_sIP3D_err,pftauchhads_sIP3D_sig,  pftauchhads_sIP2D_val,pftauchhads_sIP2D_err,pftauchhads_sIP2D_sig );
   tree->pftauchhads_IP3D_val[p] = pftauchhads_IP3D_val;
   tree->pftauchhads_IP3D_err[p] = pftauchhads_IP3D_err;
   tree->pftauchhads_IP3D_sig[p] = pftauchhads_IP3D_sig;
   tree->pftauchhads_IP2D_val[p] = pftauchhads_IP2D_val;
   tree->pftauchhads_IP2D_err[p] = pftauchhads_IP2D_err;
   tree->pftauchhads_IP2D_sig[p] = pftauchhads_IP2D_sig;
   tree->pftauchhads_sIP3D_val[p] = pftauchhads_sIP3D_val;
   tree->pftauchhads_sIP3D_err[p] = pftauchhads_sIP3D_err;
   tree->pftauchhads_sIP3D_sig[p] = pftauchhads_sIP3D_sig;
   tree->pftauchhads_sIP2D_val[p] = pftauchhads_sIP2D_val;
   tree->pftauchhads_sIP2D_err[p] = pftauchhads_sIP2D_err;
   tree->pftauchhads_sIP2D_sig[p] = pftauchhads_sIP2D_sig;
   double pftauchhads_IP1D_val = 0;
   double pftauchhads_IP1D_err = 0;
   double pftauchhads_IP1D_sig = 0;
   double pftauchhads_sIP1D_val = 0;
   double pftauchhads_sIP1D_err = 0;
   double pftauchhads_sIP1D_sig = 0;
   IP1D( candtrk,*ttrkbuilder,unbpv,pftaugv, pftauchhads_IP1D_val,pftauchhads_IP1D_err,pftauchhads_IP1D_sig, pftauchhads_sIP1D_val,pftauchhads_sIP1D_err,pftauchhads_sIP1D_sig );
   tree->pftauchhads_IP1D_val[p]  = pftauchhads_IP1D_val;
   tree->pftauchhads_IP1D_err[p]  = pftauchhads_IP1D_err;
   tree->pftauchhads_IP1D_sig[p]  = pftauchhads_IP1D_sig;
   tree->pftauchhads_sIP1D_val[p] = pftauchhads_sIP1D_val;
   tree->pftauchhads_sIP1D_err[p] = pftauchhads_sIP1D_err;
   tree->pftauchhads_sIP1D_sig[p] = pftauchhads_sIP1D_sig;
  } 
  tree->pftauchhads_numv = sigpfchhadcands.size(); 
  tree->pftauchhads_num  = sigpfchhadcands.size(); 
  //Fill tree
  tree->tree->Fill();
  break; //One candidate per event is ok
 }
}
// ------------ method called once each job just before starting event loop  ------------
void JetTauFakeRateDisc::beginJob(){
 stopwatch->Start();
 num_evt_tot = 0;
}
// ------------ method called once each job just after ending the event loop  ------------
void JetTauFakeRateDisc::endJob(){
 stopwatch->Stop();
 cout<<endl;
 cout<<"Rapid job summary "<<endl;
 cout<<num_evt_tot<<" events analysed in "<<stopwatch->RealTime()<<" seconds"<<endl;
 cout<<endl;
}
// ------------ method called when starting to processes a run  ------------
/*
void JetTauFakeRateDisc::beginRun(edm::Run const&, edm::EventSetup const&){
}
*/
// ------------ method called when ending the processing of a run  ------------
/*
void JetTauFakeRateDisc::endRun(edm::Run const&, edm::EventSetup const&){
}
*/
// ------------ method called when starting to processes a luminosity block  ------------
/*
void JetTauFakeRateDisc::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&){
}
*/
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void JetTauFakeRateDisc::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&){
}
*/
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void JetTauFakeRateDisc::fillDescriptions(edm::ConfigurationDescriptions& descriptions){
 //The following says we do not know what parameters are allowed so do no validation
 // Please change this to state exactly what you do use, even if it is no parameters
 edm::ParameterSetDescription desc;
 desc.setUnknown();
 descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(JetTauFakeRateDisc);
//Checks
  //cout<<pftau.pt()<<setw(15)<<pftau.eta()<<setw(15)<<pftau.decayMode()<<setw(15)<<posgenmatchedcand<<endl;  
   //cout<<p<<setw(15)<<cand->pt()<<setw(15)<<cand->eta()<<setw(15)<<cand->phi()<<setw(15)<<cand->charge()<<endl;
   //cout<<p<<setw(15)<<candtrk->pt()<<setw(15)<<candtrk->eta()<<setw(15)<<candtrk->phi()<<setw(15)<<candtrk->charge()<<endl;
   //cout<<"Lead val"<<endl;
   //cout<<p<<setw(15)<<pftau.leadPFChargedHadrCand()->pt()<<setw(15)<<cand->pt()<<endl;
   //cout<<pftauchhads_sIP2D_val<<setw(15)<<candtrk->dxy(pv.position())<<endl;
   //cout<<"pftauchhads kin"<<endl;
   //cout<<p<<setw(15)<<cand->pt()<<setw(15)<<cand->eta()<<setw(15)<<cand->phi()<<setw(15)<<cand->energy()<<setw(15)<<cand->charge()<<setw(15)<<endl;
   //cout<<"pftauchhads IP"<<endl;
   //cout<<p<<setw(15)<<pftauchhads_IP3D_val<<setw(15)<<pftauchhads_IP3D_err<<setw(15)<<pftauchhads_IP3D_sig<<setw(15)<<endl;
   //cout<<p<<setw(15)<<pftauchhads_IP2D_val<<setw(15)<<pftauchhads_IP2D_err<<setw(15)<<pftauchhads_IP2D_sig<<setw(15)<<endl;
   //cout<<p<<setw(15)<<pftauchhads_sIP3D_val<<setw(15)<<pftauchhads_sIP3D_err<<setw(15)<<pftauchhads_sIP3D_sig<<setw(15)<<endl;
   //cout<<p<<setw(15)<<pftauchhads_sIP2D_val<<setw(15)<<pftauchhads_sIP2D_err<<setw(15)<<pftauchhads_sIP2D_sig<<setw(15)<<endl;
