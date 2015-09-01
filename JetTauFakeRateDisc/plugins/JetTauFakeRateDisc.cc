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
#include "DataFormats/TauReco/interface/PFTauTransverseImpactParameterAssociation.h"
#include "RecoTauTag/RecoTau/interface/HPSPFRecoTauAlgorithm.h"
//Photon
//Jet and Met
#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
//Track and Vertex infos
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
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "RecoBTag/BTagTools/interface/SignedTransverseImpactParameter.h"
//Math
#include "DataFormats/Math/interface/deltaR.h"
//Store info
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
//User defined classes
#include "TauRunII/JetTauFakeRateDisc/interface/TreeVariables.h"
#include "TauRunII/JetTauFakeRateDisc/interface/HelpingFunctions.h"
#include "TauRunII/JetTauFakeRateDisc/interface/VertexRefit.h"
#include "TauRunII/JetTauFakeRateDisc/interface/TauVSTauJet.h"
#include "TauRunII/JetTauFakeRateDisc/interface/Collimation.h"
#include "TauRunII/JetTauFakeRateDisc/interface/LifeTime.h"
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
 const bool   tau_gentauh_match;
 const bool   tau_genjet_match;
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
 /////
 //   Get info by Label
 /////
 //Vertices
 Handle<vector<Vertex> > vertices;
 iEvent.getByLabel("offlinePrimaryVertices", vertices);
 //TransientTrackBuilder
 edm::ESHandle<TransientTrackBuilder> ttrkbuilder;
 iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",ttrkbuilder); 
 //GenPart Collection
 Handle<vector<GenParticle> > genParts;
 iEvent.getByLabel("genParticles", genParts);
 //GenJet Collection
 edm::Handle<vector<GenJet> > genjets;
 iEvent.getByLabel("ak4GenJets", genjets);
 //Tau collection
 Handle<vector<PFTau> > PFTaus;
 iEvent.getByLabel("hpsPFTauProducer", PFTaus);
 //For Tau lifetime info 
 typedef edm::AssociationVector<reco::PFTauRefProd, std::vector<reco::PFTauTransverseImpactParameterRef> > PFTauTIPAssociationByRef;
 edm::Handle<PFTauTIPAssociationByRef> recTauLifetimeInfos;
 iEvent.getByLabel("hpsPFTauTransverseImpactParameters", recTauLifetimeInfos);
 //For tau discriminators
 edm::Handle<reco::PFTauDiscriminator> discriminator;
 iEvent.getByLabel("hpsPFTauDiscriminationByDecayModeFindingNewDMs", discriminator);
 /////
 //   Access vertex information
 ///// 
 if(vertices->empty()) return; // skip the event if no PV found
 //Initialize tree
 tree->loop_initialize();
 tree->numrecovtcs = vertices->size();
 const reco::Vertex &pv = vertices->front();
 //Access the ditau vertex at gen level (from CV)
 if(tau_gentauh_match){
  double genditauvtx_x = 0;
  double genditauvtx_y = 0;
  double genditauvtx_z = 0;
  Get_genEventVertex(iEvent,genditauvtx_x,genditauvtx_y,genditauvtx_z);
  reco::Vertex::Point genEventVertex(genditauvtx_x,genditauvtx_y,genditauvtx_z);
  tree->genditauvtx_x = genEventVertex.x();
  tree->genditauvtx_y = genEventVertex.y();
  tree->genditauvtx_z = genEventVertex.z();
 } 
 /////
 //   Order Tau by decreasing Pt (Not too elegant sort, improve it)  
 ///// 
 vector<pair<double,int> > TauPtPos;
 for(int taupos=0; taupos<int(PFTaus->size()); taupos++) TauPtPos.push_back(make_pair((*PFTaus)[taupos].pt(),taupos)); 
 sort(TauPtPos.rbegin(),TauPtPos.rend());
 /////
 //   Access tau info
 /////
 for(int taupos=0; taupos<int(PFTaus->size()); taupos++){
  const PFTau& pftau = (*PFTaus)[TauPtPos[taupos].second];
  reco::PFTauRef RefPFTau(PFTaus, taupos);
  //Minimum tau requirements
  float discriminator_value = (*discriminator)[RefPFTau];
  if(!(pftau.pt()>tau_min_pt && fabs(pftau.eta())<tau_max_eta && discriminator_value>0.5)) continue;
  //Matching with Gen Part (posgenmatchedcand = -1 if there is no matching)
  int posgenmatchedcand = -1;
  if(tau_gentauh_match) posgenmatchedcand = MatchRecoTauhGenTauh(pftau,iEvent);
  if(tau_genjet_match)  posgenmatchedcand = MatchRecoTauhGenJet(pftau,iEvent);
  if(posgenmatchedcand==-1) continue;
  //Refit PV removing pftau ch hads
  Vertex unbpv_AVFbs;
  unbpv_AVFbs = Vertex(unbiased_vertex_AVFbs(iEvent,pv,*ttrkbuilder,pftau));
  if(!unbpv_AVFbs.isValid()) unbpv_AVFbs = pv;
  tree->unbpv_AVFbs_nv = unbpv_AVFbs.isValid();
  tree->vtxAVFbs_x = unbpv_AVFbs.position().x();
  tree->vtxAVFbs_y = unbpv_AVFbs.position().y();
  tree->vtxAVFbs_z = unbpv_AVFbs.position().z();
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
  //Direction of the reco tau
  GlobalVector pftaugv(pftau.px(), pftau.py(), pftau.pz());
  tree->pftaugv_x = pftaugv.x();
  tree->pftaugv_y = pftaugv.y();
  tree->pftaugv_z = pftaugv.z();
  //dR between reco tau/jet and gen part 
  math::PtEtaPhiELorentzVector recotaujet_lv(0., 0., 0., 0.);
  math::PtEtaPhiELorentzVector gentaujet_lv(0., 0., 0., 0.);
  if(tau_gentauh_match){
   const GenParticle & genPart = (*genParts)[posgenmatchedcand];
   math::PtEtaPhiELorentzVector tempgen(genPart.pt(), genPart.eta(), genPart.phi(), genPart.energy());
   gentaujet_lv = tempgen;
   math::PtEtaPhiELorentzVector tempreco(pftau.pt(), pftau.eta(), pftau.phi(), pftau.energy());
   recotaujet_lv = tempreco;
  }
  if(tau_genjet_match){
   const GenJet & genJet = (*genjets)[posgenmatchedcand];
   math::PtEtaPhiELorentzVector tempgen(genJet.pt(), genJet.eta(), genJet.phi(), genJet.energy());
   gentaujet_lv = tempgen;
   const PFJetRef& recojettau = pftau.jetRef();
   math::PtEtaPhiELorentzVector tempreco(recojettau->pt(), recojettau->eta(), recojettau->phi(), recojettau->energy());
   recotaujet_lv = tempreco;
  }
  double dR_RecoGen = ROOT::Math::VectorUtil::DeltaR(recotaujet_lv, gentaujet_lv);
  tree->dR_RecoGen  = dR_RecoGen;  
  /////
  //   Reco tau vs its jets (fuctions implemented in interface/TauVSTauJet.h)
  /////
  const PFJetRef& recojettau = pftau.jetRef();
  tree->recojettau_pt  = recojettau->pt();
  tree->recojettau_eta = recojettau->eta();
  tree->recojettau_phi = recojettau->phi();
  tree->recojettau_en  = recojettau->energy();
  tree->tau_pt_DIV_recojettau_pt = pftau.pt()/recojettau->pt();
  tree->tau_en_DIV_recojettau_en = pftau.energy()/recojettau->energy();
  tree->tau_trk0_pt_DIV_recojettau_pt = ratio_tau_trkn_pt_recojettau_pt(pftau,0);
  tree->tau_trk1_pt_DIV_recojettau_pt = ratio_tau_trkn_pt_recojettau_pt(pftau,1);
  tree->tau_trk2_pt_DIV_recojettau_pt = ratio_tau_trkn_pt_recojettau_pt(pftau,2);
  tree->tau_trk0_en_DIV_recojettau_en = ratio_tau_trkn_en_recojettau_en(pftau,0);
  tree->tau_trk1_en_DIV_recojettau_en = ratio_tau_trkn_en_recojettau_en(pftau,1);
  tree->tau_trk2_en_DIV_recojettau_en = ratio_tau_trkn_en_recojettau_en(pftau,2);
  //dR between reco tau and its associated jet
  double dR_recojettau_recotau = deltaR(recojettau->eta(), recojettau->phi(), pftau.eta(), pftau.phi());
  tree->dR_recojettau_recotau = dR_recojettau_recotau;
  /////
  //   Variables related to collimation of tau tracks (fuctions implemented in interface/Collimation.h) 
  /////
  //dR between tau constituent and reco tau
  tree->dR_tautrk_recotau_trk0 = dR_tautrk_recotau(pftau,0);
  tree->dR_tautrk_recotau_trk1 = dR_tautrk_recotau(pftau,1);
  tree->dR_tautrk_recotau_trk2 = dR_tautrk_recotau(pftau,2);
  //dR between tau constituent and reco jet tau
  tree->dR_tautrk_recojettau_trk0 = dR_tautrk_recojettau(pftau,0);
  tree->dR_tautrk_recojettau_trk1 = dR_tautrk_recojettau(pftau,1);
  tree->dR_tautrk_recojettau_trk2 = dR_tautrk_recojettau(pftau,2);
  //Distance between tau constituent and reco tau
  tree->pftauchhads_JTD_val_trk0 = JTD_val_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,0);
  tree->pftauchhads_JTD_val_trk1 = JTD_val_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,1);
  tree->pftauchhads_JTD_val_trk2 = JTD_val_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,2);
  tree->pftauchhads_JTD_sig_trk0 = JTD_sig_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,0);
  tree->pftauchhads_JTD_sig_trk1 = JTD_sig_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,1);
  tree->pftauchhads_JTD_sig_trk2 = JTD_sig_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,2);
  /////
  //   Variables related to tau lifetime (fuctions implemented in interface/LifeTime.h)
  /////
  //IP3D
  tree->pftauchhads_IP3D_val_trk0  = IP3D_val_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,0);
  tree->pftauchhads_IP3D_sig_trk0  = IP3D_sig_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,0);
  tree->pftauchhads_sIP3D_val_trk0 = sIP3D_val_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,0);
  tree->pftauchhads_sIP3D_sig_trk0 = sIP3D_sig_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,0);
  tree->pftauchhads_IP3D_val_trk1  = IP3D_val_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,1);
  tree->pftauchhads_IP3D_sig_trk1  = IP3D_sig_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,1);
  tree->pftauchhads_sIP3D_val_trk1 = sIP3D_val_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,1);
  tree->pftauchhads_sIP3D_sig_trk1 = sIP3D_sig_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,1);
  tree->pftauchhads_IP3D_val_trk2  = IP3D_val_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,2);
  tree->pftauchhads_IP3D_sig_trk2  = IP3D_sig_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,2);
  tree->pftauchhads_sIP3D_val_trk2 = sIP3D_val_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,2);
  tree->pftauchhads_sIP3D_sig_trk2 = sIP3D_sig_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,2);
  //IP2D
  tree->pftauchhads_IP2D_val_trk0  = IP2D_val_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,0);
  tree->pftauchhads_IP2D_sig_trk0  = IP2D_sig_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,0);
  tree->pftauchhads_sIP2D_val_trk0 = sIP2D_val_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,0);
  tree->pftauchhads_sIP2D_sig_trk0 = sIP2D_sig_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,0);
  tree->pftauchhads_IP2D_val_trk1  = IP2D_val_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,1);
  tree->pftauchhads_IP2D_sig_trk1  = IP2D_sig_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,1);
  tree->pftauchhads_sIP2D_val_trk1 = sIP2D_val_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,1);
  tree->pftauchhads_sIP2D_sig_trk1 = sIP2D_sig_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,1);
  tree->pftauchhads_IP2D_val_trk2  = IP2D_val_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,2);
  tree->pftauchhads_IP2D_sig_trk2  = IP2D_sig_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,2);
  tree->pftauchhads_sIP2D_val_trk2 = sIP2D_val_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,2);
  tree->pftauchhads_sIP2D_sig_trk2 = sIP2D_sig_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,2);
  //IP1D (Using the most common definitions of CMSSW (i.e.: with TransverseImpactPointExtrapolator))
  tree->pftauchhads_IP1D_val_trk0  = IP1D_val_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,0);
  tree->pftauchhads_IP1D_sig_trk0  = IP1D_sig_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,0);
  tree->pftauchhads_sIP1D_val_trk0 = sIP1D_val_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,0);
  tree->pftauchhads_sIP1D_sig_trk0 = sIP1D_sig_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,0);
  tree->pftauchhads_IP1D_val_trk1  = IP1D_val_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,1);
  tree->pftauchhads_IP1D_sig_trk1  = IP1D_sig_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,1);
  tree->pftauchhads_sIP1D_val_trk1 = sIP1D_val_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,1);
  tree->pftauchhads_sIP1D_sig_trk1 = sIP1D_sig_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,1);
  tree->pftauchhads_IP1D_val_trk2  = IP1D_val_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,2);
  tree->pftauchhads_IP1D_sig_trk2  = IP1D_sig_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,2);
  tree->pftauchhads_sIP1D_val_trk2 = sIP1D_val_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,2);
  tree->pftauchhads_sIP1D_sig_trk2 = sIP1D_sig_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,2);
  //IP1D (Using AnalyticalImpactPointExtrapolator (below for trk0 only, but it is rather easy to get it for trk1,2))
  tree->pftauchhads_AEIP1D_val_trk0  = AEIP1D_val_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,0);  
  tree->pftauchhads_AEIP1D_sig_trk0  = AEIP1D_sig_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,0);  
  tree->pftauchhads_AEsIP1D_val_trk0 = AEsIP1D_val_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,0);  
  tree->pftauchhads_AEsIP1D_sig_trk0 = AEsIP1D_sig_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,0);   
  //Study on IP1D
  tree->pftauchhads_AEIP1D_x_val_trk0  = AEsIP1D_x_val_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,0);
  tree->pftauchhads_AEIP1D_y_val_trk0  = AEsIP1D_y_val_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,0);
  double IP3DvalAE = 0;
  double IP3DvalTE = 0;
  double IP2DvalAE = 0;
  double IP2DvalTE = 0;
  double IP1DvalAE = 0;
  double IP1DvalTE = 0;
  double tempIPerr = 0; double tempIPsig = 0; double tempsIPval = 0; double tempsIPerr = 0; double tempsIPsig = 0;
  AEIP_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,IP3DvalAE,tempIPerr,tempIPsig,tempsIPval,tempsIPerr,tempsIPsig,3,0);
  TEIP_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,IP3DvalTE,tempIPerr,tempIPsig,tempsIPval,tempsIPerr,tempsIPsig,3,0);
  AEIP_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,IP2DvalAE,tempIPerr,tempIPsig,tempsIPval,tempsIPerr,tempsIPsig,2,0);
  TEIP_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,IP2DvalTE,tempIPerr,tempIPsig,tempsIPval,tempsIPerr,tempsIPsig,2,0);
  AEIP_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,IP1DvalAE,tempIPerr,tempIPsig,tempsIPval,tempsIPerr,tempsIPsig,1,0);
  TEIP_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,IP1DvalTE,tempIPerr,tempIPsig,tempsIPval,tempsIPerr,tempsIPsig,1,0);
  tree->pftauchhads_IP3DvalAE_trk0 = IP3DvalAE;
  tree->pftauchhads_IP3DvalTE_trk0 = IP3DvalTE;
  tree->pftauchhads_IP2DvalAE_trk0 = IP2DvalAE;
  tree->pftauchhads_IP2DvalTE_trk0 = IP2DvalTE;
  tree->pftauchhads_IP1DvalAE_trk0 = IP1DvalAE;
  tree->pftauchhads_IP1DvalTE_trk0 = IP1DvalTE;
  tree->pftauchhads_IP3DvalAE_2DTE_1DAE_trk0 = IP3DvalAE-sqrt(pow(IP2DvalTE,2)+pow(IP1DvalAE,2));
  tree->pftauchhads_IP3DvalAE_2DTE_1DTE_trk0 = IP3DvalAE-sqrt(pow(IP2DvalTE,2)+pow(IP1DvalTE,2));
  tree->pftauchhads_IP3DvalAE_2DAE_1DAE_trk0 = IP3DvalAE-sqrt(pow(IP2DvalAE,2)+pow(IP1DvalAE,2));
  tree->pftauchhads_IP3DvalTE_2DTE_1DTE_trk0 = IP3DvalTE-sqrt(pow(IP2DvalTE,2)+pow(IP1DvalTE,2));
  //DecayLenght3D
  tree->pftauchhads_sDL3D_val_trk0 = sDL3D_val_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,0);
  tree->pftauchhads_sDL3D_val_trk1 = sDL3D_val_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,1);
  tree->pftauchhads_sDL3D_val_trk2 = sDL3D_val_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,2);
  tree->pftauchhads_sDL3D_sig_trk0 = sDL3D_sig_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,0);
  tree->pftauchhads_sDL3D_sig_trk1 = sDL3D_sig_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,1);
  tree->pftauchhads_sDL3D_sig_trk2 = sDL3D_sig_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,2);
  tree->pftauchhads_absDL3D_val_trk0 = absDL3D_val_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,0);
  tree->pftauchhads_absDL3D_val_trk1 = absDL3D_val_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,1);
  tree->pftauchhads_absDL3D_val_trk2 = absDL3D_val_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,2);
  tree->pftauchhads_absDL3D_sig_trk0 = absDL3D_sig_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,0);
  tree->pftauchhads_absDL3D_sig_trk1 = absDL3D_sig_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,1);
  tree->pftauchhads_absDL3D_sig_trk2 = absDL3D_sig_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,2);
  //Tau flight distance
  Vertex taudecvtx = Vertex(get_taudecvtx(*ttrkbuilder,pftau));
  tree->pftau_pvsv_dist3d_val = abs_pvsv_dist3d_val(pftau,*ttrkbuilder,unbpv_AVFbs,taudecvtx,pftaugv);
  tree->pftau_pvsv_dist2d_val = abs_pvsv_dist2d_val(pftau,*ttrkbuilder,unbpv_AVFbs,taudecvtx,pftaugv);
  tree->pftau_pvsv_dist3d_sig = abs_pvsv_dist3d_sig(pftau,*ttrkbuilder,unbpv_AVFbs,taudecvtx,pftaugv);
  tree->pftau_pvsv_dist2d_sig = abs_pvsv_dist2d_sig(pftau,*ttrkbuilder,unbpv_AVFbs,taudecvtx,pftaugv);
  //dR_TauFlightDist_JetAxis
  tree->dR_TauFlightDist_JetAxis = dR_TauFlightDist_JetAxis(pftau,taudecvtx,unbpv_AVFbs,pftaugv);
  //Chi2PV - Chi2PVnonTauTrks
  TransientVertex refittedpv_AVFbs = refitted_vertex_AVFbs(iEvent,pv,*ttrkbuilder);
  TransientVertex unbiasedpv_AVFbs = unbiased_vertex_AVFbs(iEvent,pv,*ttrkbuilder,pftau);
  tree->pv_avfbs_nchi2                        = refittedpv_AVFbs.totalChiSquared()/refittedpv_AVFbs.degreesOfFreedom();
  tree->unbpv_avfbs_nchi2                     = unbiasedpv_AVFbs.totalChiSquared()/unbiasedpv_AVFbs.degreesOfFreedom();
  tree->diff_pv_avfbs_nchi2_unbpv_avfbs_nchi2 = (refittedpv_AVFbs.totalChiSquared()/refittedpv_AVFbs.degreesOfFreedom()) - (unbiasedpv_AVFbs.totalChiSquared()/unbiasedpv_AVFbs.degreesOfFreedom());
  /////
  //   Access tau tracks info
  /////
  const vector<reco::PFCandidatePtr>& sigpfchhadcands = pftau.signalPFChargedHadrCands();
  for(uint p=0; p<sigpfchhadcands.size(); p++){
   PFCandidatePtr cand = sigpfchhadcands[p];
   const Track* candtrk = cand->bestTrack();
   if(!candtrk) continue;
   //if(!is_goodtrk(candtrk,pv)) continue;
   //Kinematic
   tree->pftauchhads_pt[p]  = cand->pt();
   tree->pftauchhads_eta[p] = cand->eta();
   tree->pftauchhads_phi[p] = cand->phi();
   tree->pftauchhads_en[p]  = cand->energy();
   //Charge
   tree->pftauchhads_ch[p]  = cand->charge();
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
/////
//   Checks
/////
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
  //cout<<"Primary  vertex used"<<endl;
  //cout<<"PV   "<<setw(30)<<pv.position().x()<<setw(30)<<pv.position().y()<<setw(30)<<pv.position().z()<<endl;
  //cout<<"unbKVF"<<setw(30)<<unbpv_KVF.position().x()<<setw(30)<<unbpv_KVF.position().y()<<setw(30)<<unbpv_KVF.position().z()<<endl;
  //cout<<"unbKVFbs"<<setw(30)<<unbpv_KVFbs.position().x()<<setw(30)<<unbpv_KVFbs.position().y()<<setw(30)<<unbpv_KVFbs.position().z()<<endl;
  //cout<<"unbAVF"<<setw(30)<<unbpv_AVF.position().x()<<setw(30)<<unbpv_AVF.position().y()<<setw(30)<<unbpv_AVF.position().z()<<endl;
  //cout<<"unbAVFbs"<<setw(30)<<unbpv_AVFbs.position().x()<<setw(30)<<unbpv_AVFbs.position().y()<<setw(30)<<unbpv_AVFbs.position().z()<<endl;
  //cout<<"Their differences"<<endl;
  //cout<<"PV-refPV"<<setw(30)<<fabs(pv.position().x()-refpv.position().x())<<setw(30)<<fabs(pv.position().y()-refpv.position().y())<<setw(30)<<fabs(pv.position().z()-refpv.position().z())<<endl;
  //cout<<"Difference between the different primary vertices"<<endl;
  //cout<<setw(20)<<"Difference"<<setw(20)<<"in x"<<setw(40)<<"in y"<<setw(40)<<"in z"<<endl;
  //cout<<setw(20)<<"PV-unbKVF"<<setw(20)<<fabs(pv.position().x()-unbpv_KVF.position().x())<<setw(40)<<fabs(pv.position().y()-unbpv_KVF.position().y())<<setw(40)<<fabs(pv.position().z()-unbpv_KVF.position().z())<<endl;
  //cout<<setw(20)<<refittedpv_AVFbs.totalChiSquared()<<setw(20)<<tv.totalChiSquared()<<setw(20)<<refittedpv_AVFbs.degreesOfFreedom()<<setw(20)<<tv.degreesOfFreedom()<<endl;
  //cout<<setw(20)<<Vertex(refittedpv_AVFbs).chi2()<<setw(20)<<Vertex(tv).chi2()<<setw(20)<<Vertex(refittedpv_AVFbs).ndof()<<setw(20)<<Vertex(tv).ndof()<<endl;
/////
//   Multiplicity of jet constituent
/////
  //const vector<reco::PFCandidatePtr>& sigpfchhadcands2 = 
  //int recojettau_n = recojettau->getJetConstituents().size();
  //cout<<"Num constituents"<<setw(20)<<sigpfchhadcands.size()<<setw(20)<<recojettau_n<<endl;
/////
//   Study of the definition of IP1D 
/////
  //cout<<"IP "<<setw(20)<<"3D"<<setw(20)<<"1D"<<endl;
  //cout<<"Val"<<setw(20)<<IP3D_val_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,0)<<setw(20)<<IP1D_val_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,0)<<endl;
  //cout<<"Sig"<<setw(20)<<IP3D_sig_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,0)<<setw(20)<<IP1D_sig_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,0)<<endl;
  //cout<<endl;
  /*
  double pftauchhads_AEIP1D_val = 0;
  double pftauchhads_AEIP1D_err = 0;
  double pftauchhads_AEIP1D_sig = 0;
  double pftauchhads_AEsIP1D_val = 0;
  double pftauchhads_AEsIP1D_err = 0;
  double pftauchhads_AEsIP1D_sig = 0; 
  AEIP_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,pftauchhads_AEIP1D_val,pftauchhads_AEIP1D_err,pftauchhads_AEIP1D_sig,pftauchhads_AEsIP1D_val,pftauchhads_AEsIP1D_err,pftauchhads_AEsIP1D_sig,1,0);
  double pftauchhads_TEIP1D_val = 0;
  double pftauchhads_TEIP1D_err = 0;
  double pftauchhads_TEIP1D_sig = 0;
  double pftauchhads_TEsIP1D_val = 0;
  double pftauchhads_TEsIP1D_err = 0;
  double pftauchhads_TEsIP1D_sig = 0; 
  TEIP_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,pftauchhads_TEIP1D_val,pftauchhads_TEIP1D_err,pftauchhads_TEIP1D_sig,pftauchhads_TEsIP1D_val,pftauchhads_TEsIP1D_err,pftauchhads_TEsIP1D_sig,1,0);
  double pftauchhads_zIP1D_val = 0;
  double pftauchhads_zIP1D_err = 0;
  double pftauchhads_zIP1D_sig = 0;
  double pftauchhads_zsIP1D_val = 0;
  double pftauchhads_zsIP1D_err = 0;
  double pftauchhads_zsIP1D_sig = 0; 
  zIP1D(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,pftauchhads_zIP1D_val,pftauchhads_zIP1D_err,pftauchhads_zIP1D_sig,pftauchhads_zsIP1D_val,pftauchhads_zsIP1D_err,pftauchhads_zsIP1D_sig,0);
  cout<<"IP1D"<<setw(20)<<"standard"<<setw(20)<<"TE"<<setw(20)<<"AE"<<setw(20)<<"zIP"<<endl;
  cout<<"Val "<<setw(20)<<sIP1D_val_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,0)<<setw(20)<<pftauchhads_TEsIP1D_val<<setw(20)<<pftauchhads_AEsIP1D_val<<setw(20)<<pftauchhads_zsIP1D_val<<endl;
  cout<<"Sig "<<setw(20)<<sIP1D_sig_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,0)<<setw(20)<<pftauchhads_TEsIP1D_sig<<setw(20)<<pftauchhads_AEsIP1D_sig<<setw(20)<<pftauchhads_zsIP1D_sig<<endl;
  //Different definitions of IP1D used in CMSSW (Be carefull you are using the same vertex)
  //const vector<reco::PFCandidatePtr>& sigpfchhadcands = pftau.signalPFChargedHadrCands();
  //PFCandidatePtr cand = sigpfchhadcands[0];
  //const Track* trk = cand->bestTrack();
  //cout<<setw(20)<<"Perigree"<<setw(20)<<"trk.dz(PV)"<<setw(20)<<"zIP"<<setw(20)<<"trk.vtx.z-vtx.z"<<endl;//trk.vtx.z-vtx.z seems diff
  //cout<<setw(20)<<IP1D_val_trkn(pftau,*ttrkbuilder,pv,0)<<setw(20)<<trk->dz(pv.position())<<setw(20)<<pftauchhads_zIP1D_val<<setw(20)<<abs(trk->vertex().z()-pv.position().z())<<endl;
  double pftauchhads_AEIP2D_val = 0;
  double pftauchhads_AEIP2D_err = 0;
  double pftauchhads_AEIP2D_sig = 0;
  double pftauchhads_AEsIP2D_val = 0;
  double pftauchhads_AEsIP2D_err = 0;
  double pftauchhads_AEsIP2D_sig = 0; 
  AEIP_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,pftauchhads_AEIP2D_val,pftauchhads_AEIP2D_err,pftauchhads_AEIP2D_sig,pftauchhads_AEsIP2D_val,pftauchhads_AEsIP2D_err,pftauchhads_AEsIP2D_sig,2,0);
  double pftauchhads_TEIP2D_val = 0;
  double pftauchhads_TEIP2D_err = 0;
  double pftauchhads_TEIP2D_sig = 0;
  double pftauchhads_TEsIP2D_val = 0;
  double pftauchhads_TEsIP2D_err = 0;
  double pftauchhads_TEsIP2D_sig = 0; 
  TEIP_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,pftauchhads_TEIP2D_val,pftauchhads_TEIP2D_err,pftauchhads_TEIP2D_sig,pftauchhads_TEsIP2D_val,pftauchhads_TEsIP2D_err,pftauchhads_TEsIP2D_sig,2,0);
  cout<<"sIP2D"<<setw(20)<<"standard"<<setw(20)<<"TE"<<setw(20)<<"AE"<<endl;
  cout<<"Val "<<setw(20)<<sIP2D_val_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,0)<<setw(20)<<pftauchhads_TEsIP2D_val<<setw(20)<<pftauchhads_AEsIP2D_val<<endl;
  cout<<"Sig "<<setw(20)<<sIP2D_sig_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,0)<<setw(20)<<pftauchhads_TEsIP2D_sig<<setw(20)<<pftauchhads_AEsIP2D_sig<<endl;
  double pftauchhads_AEIP3D_val = 0;
  double pftauchhads_AEIP3D_err = 0;
  double pftauchhads_AEIP3D_sig = 0;
  double pftauchhads_AEsIP3D_val = 0;
  double pftauchhads_AEsIP3D_err = 0;
  double pftauchhads_AEsIP3D_sig = 0; 
  AEIP_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,pftauchhads_AEIP3D_val,pftauchhads_AEIP3D_err,pftauchhads_AEIP3D_sig,pftauchhads_AEsIP3D_val,pftauchhads_AEsIP3D_err,pftauchhads_AEsIP3D_sig,3,0);
  double pftauchhads_TEIP3D_val = 0;
  double pftauchhads_TEIP3D_err = 0;
  double pftauchhads_TEIP3D_sig = 0;
  double pftauchhads_TEsIP3D_val = 0;
  double pftauchhads_TEsIP3D_err = 0;
  double pftauchhads_TEsIP3D_sig = 0; 
  TEIP_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,pftauchhads_TEIP3D_val,pftauchhads_TEIP3D_err,pftauchhads_TEIP3D_sig,pftauchhads_TEsIP3D_val,pftauchhads_TEsIP3D_err,pftauchhads_TEsIP3D_sig,3,0);
  cout<<"sIP3D"<<setw(20)<<"standard"<<setw(20)<<"TE"<<setw(20)<<"AE"<<endl;
  cout<<"Val "<<setw(20)<<sIP3D_val_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,0)<<setw(20)<<pftauchhads_TEsIP3D_val<<setw(20)<<pftauchhads_AEsIP3D_val<<endl;
  cout<<"Sig "<<setw(20)<<sIP3D_sig_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,0)<<setw(20)<<pftauchhads_TEsIP3D_sig<<setw(20)<<pftauchhads_AEsIP3D_sig<<endl;
  //AEIP1D   
  tree->pftauchhads_AEIP1D_val_trk0  = AEIP1D_val_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,0);
  tree->pftauchhads_AEIP1D_val_trk1  = AEIP1D_val_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,1);
  tree->pftauchhads_AEIP1D_val_trk2  = AEIP1D_val_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,2);
  tree->pftauchhads_AEIP1D_sig_trk0  = AEIP1D_sig_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,0);
  tree->pftauchhads_AEIP1D_sig_trk1  = AEIP1D_sig_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,1);
  tree->pftauchhads_AEIP1D_sig_trk2  = AEIP1D_sig_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,2);
  tree->pftauchhads_AEsIP1D_val_trk0 = AEsIP1D_val_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,0);
  tree->pftauchhads_AEsIP1D_val_trk1 = AEsIP1D_val_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,1);
  tree->pftauchhads_AEsIP1D_val_trk2 = AEsIP1D_val_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,2);
  tree->pftauchhads_AEsIP1D_sig_trk0 = AEsIP1D_sig_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,0);
  tree->pftauchhads_AEsIP1D_sig_trk1 = AEsIP1D_sig_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,1);
  tree->pftauchhads_AEsIP1D_sig_trk2 = AEsIP1D_sig_trkn(pftau,*ttrkbuilder,unbpv_AVFbs,pftaugv,2);
  */
