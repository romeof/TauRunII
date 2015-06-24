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
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
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
 /////
 //   Initial information
 /////
 //Vertex collection
 Handle<vector<Vertex> > vertices;
 iEvent.getByLabel("offlinePrimaryVertices", vertices);
 if(vertices->empty()) return; // skip the event if no PV found
 const reco::Vertex &pv = vertices->front();
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
 //Order in Pt (Not too elegant sort, improve it)  
 vector<pair<double,int> > TauPtPos;
 for(int taupos=0; taupos<int(PFTaus->size()); taupos++) TauPtPos.push_back(make_pair((*PFTaus)[taupos].pt(),taupos)); 
 sort(TauPtPos.rbegin(),TauPtPos.rend());
 tree->loop_initialize();
 //Access the ditau vertex at gen level (from CV)
 double genditauvtx_x = 0;
 double genditauvtx_y = 0;
 double genditauvtx_z = 0;
 Get_genEventVertex(iEvent,genditauvtx_x,genditauvtx_y,genditauvtx_z);
 reco::Vertex::Point genEventVertex(genditauvtx_x,genditauvtx_y,genditauvtx_z);
 tree->genditauvtx_x = genEventVertex.x();
 tree->genditauvtx_y = genEventVertex.y();
 tree->genditauvtx_z = genEventVertex.z();
 /////
 //   Access tau info
 /////
 for(int taupos=0; taupos<int(PFTaus->size()); taupos++){
  const PFTau& pftau = (*PFTaus)[TauPtPos[taupos].second];
  reco::PFTauRef RefPFTau(PFTaus, taupos);
  //Minimum tau requirements
  float discriminator_value = (*discriminator)[RefPFTau];
  if(!(pftau.pt()>tau_min_pt && fabs(pftau.eta())<tau_max_eta && discriminator_value>0.5)) continue;
  //pftau.decayMode()!=-1)) continue;
  //Matching with Gen Part (posgenmatchedcand = -1 if there is no matching)
  int posgenmatchedcand = -1;
  if(tau_gentauh_match) posgenmatchedcand = MatchRecoTauhGenTauh(pftau,iEvent);
  if(tau_genjet_match)  posgenmatchedcand = MatchRecoTauhGenJet(pftau,iEvent);
  if(posgenmatchedcand==-1) continue;
  //Now get the info you need for the study
  tree->numrecovtcs = vertices->size();
  //Access also the refitted primary vertex associated to tau, using already implemented methods
  const reco::PFTauTransverseImpactParameter& recTauLifetimeInfo = *(*recTauLifetimeInfos)[RefPFTau];
  const reco::VertexRef pvtauref = recTauLifetimeInfo.primaryVertex();
  //Vertex resolution
  tree->pvtauvtx_genditauvtx_x = pvtauref->position().x()-genEventVertex.x();
  tree->pvtauvtx_genditauvtx_y = pvtauref->position().y()-genEventVertex.y();
  tree->pvtauvtx_genditauvtx_z = pvtauref->position().z()-genEventVertex.z();
  //Refit PV removing pftau ch hads
  Vertex unbpv_KVF;
  unbpv_KVF = Vertex(unbiased_vertex_KVF(pv,*ttrkbuilder,pftau));
  if(!unbpv_KVF.isValid()) unbpv_KVF = pv;
  tree->unbpv_KVF_nv = unbpv_KVF.isValid();
  Vertex unbpv_KVFbs;
  unbpv_KVFbs = Vertex(unbiased_vertex_KVFbs(iEvent,pv,*ttrkbuilder,pftau));
  if(!unbpv_KVFbs.isValid()) unbpv_KVFbs = pv;
  tree->unbpv_KVFbs_nv = unbpv_KVFbs.isValid();
  Vertex unbpv_AVF;
  unbpv_AVF = Vertex(unbiased_vertex_AVF(pv,*ttrkbuilder,pftau));
  if(!unbpv_AVF.isValid()) unbpv_AVF = pv;
  tree->unbpv_AVF_nv = unbpv_AVF.isValid();
  Vertex unbpv_AVFbs;
  unbpv_AVFbs = Vertex(unbiased_vertex_AVFbs(iEvent,pv,*ttrkbuilder,pftau));
  if(!unbpv_AVFbs.isValid()) unbpv_AVFbs = pv;
  tree->unbpv_AVFbs_nv = unbpv_AVFbs.isValid();
  tree->vtxKVF_x   = unbpv_KVF.position().x();
  tree->vtxKVF_y   = unbpv_KVF.position().y();
  tree->vtxKVF_z   = unbpv_KVF.position().z();
  tree->vtxKVFbs_x = unbpv_KVFbs.position().x();
  tree->vtxKVFbs_y = unbpv_KVFbs.position().y();
  tree->vtxKVFbs_z = unbpv_KVFbs.position().z();
  tree->vtxAVF_x   = unbpv_AVF.position().x();
  tree->vtxAVF_y   = unbpv_AVF.position().y();
  tree->vtxAVF_z   = unbpv_AVF.position().z();
  tree->vtxAVFbs_x = unbpv_AVFbs.position().x();
  tree->vtxAVFbs_y = unbpv_AVFbs.position().y();
  tree->vtxAVFbs_z = unbpv_AVFbs.position().z();
  tree->unbpv_KVF_genditauvtx_x = unbpv_KVF.position().x()-genEventVertex.x();
  tree->unbpv_KVF_genditauvtx_y = unbpv_KVF.position().y()-genEventVertex.y();
  tree->unbpv_KVF_genditauvtx_z = unbpv_KVF.position().z()-genEventVertex.z();
  //Kinematic
  tree->tau_pt  = pftau.pt();
  tree->tau_eta = pftau.eta();
  tree->tau_phi = pftau.phi();
  tree->tau_en  = pftau.energy();
  //Charge
  tree->tau_ch  = pftau.charge();
  /////
  //   Reco tau vs its jets
  /////
  const PFJetRef& recojettau = pftau.jetRef();
  //Multiplicity of jet constituent
  //const vector<reco::PFCandidatePtr>& sigpfchhadcands2 = 
  //int recojettau_n = recojettau->getJetConstituents().size();
  //cout<<"Num constituents"<<setw(20)<<sigpfchhadcands.size()<<setw(20)<<recojettau_n<<endl;
  //Ratio between reco tau tracks and jet tracks
  //dR between reco tau and its associated jet
  double dR_recojettau_recotau = deltaR(recojettau->eta(), recojettau->phi(), pftau.eta(), pftau.phi());
  //cout<<"dR between reco tau and its associated jet"<<setw(20)<<dR_recojettau_recotau<<endl;
  tree->dR_recojettau_recotau = dR_recojettau_recotau;
  //Tau IP
  tree->tau_vtxdz   = TMath::Abs(pftau.vertex().z()-pv.position().z());
  tree->tau_vtxdxy  = sqrt(pow(pftau.vertex().x()-pv.position().x(),2)+pow(pftau.vertex().y()-pv.position().y(),2));
  tree->tau_vtxdxyz = sqrt(pow(pftau.vertex().x()-pv.position().x(),2)+pow(pftau.vertex().y()-pv.position().y(),2)+pow(pftau.vertex().z()-pv.position().z(),2));
  //Info of tau constituents (for tau ch had)
  //GlobalVector pftaugv(0.,0.,0.);
  //To do: Use python booleans to choose the direction for the sign the IPs
  //Direction of the reco tau
  GlobalVector pftaugv(pftau.px(), pftau.py(), pftau.pz());
  //Direction of gen tau associated to the reco tau  
  //if(tau_gentauh_match){
  //  const GenParticle & genPart = (*genParts)[posgenmatchedcand];
  //  GlobalVector temp(genPart.px(), genPart.py(), genPart.pz());
  //  pftaugv = temp;
  //}
  //Direction of gen jet associated to the reco tau
  //if(tau_genjet_match){
  //  const GenJet & genJet = (*genjets)[posgenmatchedcand];
  //  GlobalVector temp(genJet.px(), genJet.py(), genJet.pz());
  //  pftaugv = temp;
  //}
  //Direction of the unbiased vertex (unbpv_KVF) including the tau tracks
  //TransientVertex unbpv_KVF_withtautrks = unbiased_vertex_KVF_withtautrks(pv,*ttrkbuilder,pftau);
  //if(unbpv_KVF_withtautrks.isValid()){
  // GlobalVector temp(unbpv_KVF_withtautrks.position().x()-unbpv_KVF.position().x(),
  //                   unbpv_KVF_withtautrks.position().y()-unbpv_KVF.position().y(),
  //                   unbpv_KVF_withtautrks.position().z()-unbpv_KVF.position().z()
  //                  );
  // pftaugv = temp;
  //}else{
  // GlobalVector temp(pftau.px(), pftau.py(), pftau.pz());
  // pftaugv = temp; 
  //}
  tree->pftaugv_x = pftaugv.x();
  tree->pftaugv_y = pftaugv.y();
  tree->pftaugv_z = pftaugv.z();
  //Access tau tracks info
  const vector<reco::PFCandidatePtr>& sigpfchhadcands = pftau.signalPFChargedHadrCands();
  for(uint p=0; p<sigpfchhadcands.size(); p++){
   PFCandidatePtr cand = sigpfchhadcands[p];
   const Track* candtrk = cand->bestTrack();
   if(!candtrk) continue;
   //To do: check results without is_goodtrk condition
   //if(!is_goodtrk(candtrk,pv)) continue;
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
   //To do: temporary initialisation (pftaugv may not always be from vtxDir)
   //double dR_vtxDirGen = deltaR(pftaugv.eta(), pftaugv.phi(), gentaujet_lv.eta(), gentaujet_lv.phi());
   //tree->dR_vtxDirGen  = dR_vtxDirGen;
   //Kinematic
   tree->pftauchhads_pt[p]  = cand->pt();
   tree->pftauchhads_eta[p] = cand->eta();
   tree->pftauchhads_phi[p] = cand->phi();
   tree->pftauchhads_en[p]  = cand->energy();
   //Charge
   tree->pftauchhads_ch[p]  = cand->charge();
   /////
   //   Variables related to tau lifetime
   /////
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
   IP3D2D( candtrk,*ttrkbuilder,unbpv_AVFbs,pftaugv, pftauchhads_IP3D_val,pftauchhads_IP3D_err,pftauchhads_IP3D_sig,  pftauchhads_IP2D_val,pftauchhads_IP2D_err,pftauchhads_IP2D_sig, pftauchhads_sIP3D_val,pftauchhads_sIP3D_err,pftauchhads_sIP3D_sig,  pftauchhads_sIP2D_val,pftauchhads_sIP2D_err,pftauchhads_sIP2D_sig );
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
   //Def fit
   double pftauchhads_IP3Ddv_val = 0;
   double pftauchhads_IP3Ddv_err = 0;
   double pftauchhads_IP3Ddv_sig = 0;
   double pftauchhads_IP2Ddv_val = 0;
   double pftauchhads_IP2Ddv_err = 0;
   double pftauchhads_IP2Ddv_sig = 0;
   double pftauchhads_sIP3Ddv_val = 0;
   double pftauchhads_sIP3Ddv_err = 0;
   double pftauchhads_sIP3Ddv_sig = 0;
   double pftauchhads_sIP2Ddv_val = 0;
   double pftauchhads_sIP2Ddv_err = 0;
   double pftauchhads_sIP2Ddv_sig = 0;
   IP3D2D( candtrk,*ttrkbuilder,pvtauref,pftaugv, pftauchhads_IP3Ddv_val,pftauchhads_IP3Ddv_err,pftauchhads_IP3Ddv_sig,  pftauchhads_IP2Ddv_val,pftauchhads_IP2Ddv_err,pftauchhads_IP2Ddv_sig, pftauchhads_sIP3Ddv_val,pftauchhads_sIP3Ddv_err,pftauchhads_sIP3Ddv_sig,  pftauchhads_sIP2Ddv_val,pftauchhads_sIP2Ddv_err,pftauchhads_sIP2Ddv_sig );
   tree->pftauchhads_IP3Ddv_val[p] = pftauchhads_IP3Ddv_val;
   tree->pftauchhads_IP3Ddv_err[p] = pftauchhads_IP3Ddv_err;
   tree->pftauchhads_IP3Ddv_sig[p] = pftauchhads_IP3Ddv_sig;
   tree->pftauchhads_IP2Ddv_val[p] = pftauchhads_IP2Ddv_val;
   tree->pftauchhads_IP2Ddv_err[p] = pftauchhads_IP2Ddv_err;
   tree->pftauchhads_IP2Ddv_sig[p] = pftauchhads_IP2Ddv_sig;
   tree->pftauchhads_sIP3Ddv_val[p] = pftauchhads_sIP3Ddv_val;
   tree->pftauchhads_sIP3Ddv_err[p] = pftauchhads_sIP3Ddv_err;
   tree->pftauchhads_sIP3Ddv_sig[p] = pftauchhads_sIP3Ddv_sig;
   tree->pftauchhads_sIP2Ddv_val[p] = pftauchhads_sIP2Ddv_val;
   tree->pftauchhads_sIP2Ddv_err[p] = pftauchhads_sIP2Ddv_err;
   tree->pftauchhads_sIP2Ddv_sig[p] = pftauchhads_sIP2Ddv_sig;
   //IP1D
   double pftauchhads_IP1D_val = 0;
   double pftauchhads_IP1D_err = 0;
   double pftauchhads_IP1D_sig = 0;
   double pftauchhads_sIP1D_val = 0;
   double pftauchhads_sIP1D_err = 0;
   double pftauchhads_sIP1D_sig = 0;
   //IP1D( candtrk,*ttrkbuilder,unbpv_KVF,pftaugv, pftauchhads_IP1D_val,pftauchhads_IP1D_err,pftauchhads_IP1D_sig, pftauchhads_sIP1D_val,pftauchhads_sIP1D_err,pftauchhads_sIP1D_sig );
   IP1D( candtrk,*ttrkbuilder,pvtauref,pftaugv, pftauchhads_IP1D_val,pftauchhads_IP1D_err,pftauchhads_IP1D_sig, pftauchhads_sIP1D_val,pftauchhads_sIP1D_err,pftauchhads_sIP1D_sig );
   tree->pftauchhads_IP1D_val[p]  = pftauchhads_IP1D_val;
   tree->pftauchhads_IP1D_err[p]  = pftauchhads_IP1D_err;
   tree->pftauchhads_IP1D_sig[p]  = pftauchhads_IP1D_sig;
   tree->pftauchhads_sIP1D_val[p] = pftauchhads_sIP1D_val;
   tree->pftauchhads_sIP1D_err[p] = pftauchhads_sIP1D_err;
   tree->pftauchhads_sIP1D_sig[p] = pftauchhads_sIP1D_sig;
   /////
   //   Variables related to collimation of tau tracks 
   /////
   //dR between tau constituent and reco tau
   double dR_tautrk_recotau = deltaR(cand->eta(), cand->phi(), pftau.eta(), pftau.phi());
   //cout<<"dR between tau constituent and reco tau"<<setw(20)<<dR_tautrk_recotau<<endl;
   tree->dR_tautrk_recotau[p]     = dR_tautrk_recotau;
   //Distance between tau constituent and reco tau
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
  //cout<<"Primary  vertex used"<<endl;
  //cout<<"PV   "<<setw(30)<<pv.position().x()<<setw(30)<<pv.position().y()<<setw(30)<<pv.position().z()<<endl;
  //cout<<"unbKVF"<<setw(30)<<unbpv_KVF.position().x()<<setw(30)<<unbpv_KVF.position().y()<<setw(30)<<unbpv_KVF.position().z()<<endl;
  //cout<<"unbKVFbs"<<setw(30)<<unbpv_KVFbs.position().x()<<setw(30)<<unbpv_KVFbs.position().y()<<setw(30)<<unbpv_KVFbs.position().z()<<endl;
  //cout<<"unbAVF"<<setw(30)<<unbpv_AVF.position().x()<<setw(30)<<unbpv_AVF.position().y()<<setw(30)<<unbpv_AVF.position().z()<<endl;
  //cout<<"unbAVFbs"<<setw(30)<<unbpv_AVFbs.position().x()<<setw(30)<<unbpv_AVFbs.position().y()<<setw(30)<<unbpv_AVFbs.position().z()<<endl;
  //cout<<"tauPV"<<setw(30)<<pvtauref->position().x()<<setw(30)<<pvtauref->position().y()<<setw(30)<<pvtauref->position().z()<<endl;
  //cout<<"Their differences"<<endl;
  //cout<<"PV-refPV"<<setw(30)<<fabs(pv.position().x()-refpv.position().x())<<setw(30)<<fabs(pv.position().y()-refpv.position().y())<<setw(30)<<fabs(pv.position().z()-refpv.position().z())<<endl;
  //cout<<"Difference between the different primary vertices"<<endl;
  //cout<<setw(20)<<"Difference"<<setw(20)<<"in x"<<setw(40)<<"in y"<<setw(40)<<"in z"<<endl;
  //cout<<setw(20)<<"PV-unbKVF"<<setw(20)<<fabs(pv.position().x()-unbpv_KVF.position().x())<<setw(40)<<fabs(pv.position().y()-unbpv_KVF.position().y())<<setw(40)<<fabs(pv.position().z()-unbpv_KVF.position().z())<<endl;
  //cout<<setw(20)<<"PV-tauPV"<<setw(20)<<fabs(pv.position().x()-pvtauref->position().x())<<setw(40)<<fabs(pv.position().y()-pvtauref->position().y())<<setw(40)<<fabs(pv.position().z()-pvtauref->position().z())<<endl;
  //cout<<setw(20)<<"unbKVF-tauPV"<<setw(20)<<fabs(unbpv_KVF.position().x()-pvtauref->position().x())<<setw(40)<<fabs(unbpv_KVF.position().y()-pvtauref->position().y())<<setw(40)<<fabs(unbpv_KVF.position().z()-pvtauref->position().z())<<endl;
