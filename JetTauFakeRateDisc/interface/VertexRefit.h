#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexUpdator.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexTrackCompatibilityEstimator.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexSmoother.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalTrajectoryExtrapolatorToLine.h"
KalmanVertexFitter   vtx_kvf(true);
//Namespaces
using namespace reco;
using namespace edm;
using namespace std;
/////
//   Get the refitted vertex using vtx tracks 
/////
//Get the transient tracks collection
std::vector<reco::TransientTrack> get_TransientTracks_Vertex(const Event& iEvent, const reco::Vertex vtx, const TransientTrackBuilder& ttrkbuilder){
 std::vector<reco::TransientTrack> vtxttrk;
 //Get Vertex Tracks 
 reco::TrackCollection vertexTracks;
 for(std::vector<reco::TrackBaseRef>::const_iterator vtxTrkRef=vtx.tracks_begin();vtxTrkRef<vtx.tracks_end();vtxTrkRef++){
  vertexTracks.push_back(**vtxTrkRef);
 }
 //Get Transient Tracks
 for(reco::TrackCollection::iterator iter=vertexTracks.begin(); iter!=vertexTracks.end(); ++iter){
  vtxttrk.push_back(ttrkbuilder.build(*iter));
 }
 return vtxttrk;
}
//Refit the vertex using Kalman vertex fitting (KVF) technique
TransientVertex refitted_vertex_KVF(const reco::Vertex vtx, const TransientTrackBuilder& ttrkbuilder){
 TransientVertex unbvtx;
 vector<TransientTrack> vtxttrk;
 for(vector<reco::TrackBaseRef>::const_iterator pvtrk=vtx.tracks_begin(); pvtrk!=vtx.tracks_end(); pvtrk++){
  if(!(pvtrk->isNonnull())) continue;
  vtxttrk.push_back(ttrkbuilder.build(**pvtrk));
 }
 unbvtx = vtx_kvf.vertex(vtxttrk);
 return unbvtx;
}
//Refit the vertex using Kalman vertex fitting (KVF) technique with beam spot constraint
TransientVertex refitted_vertex_KVFbs(const Event& iEvent, const reco::Vertex vtx, const TransientTrackBuilder& ttrkbuilder){
 TransientVertex unbvtx;
 vector<TransientTrack> vtxttrk;
 for(vector<reco::TrackBaseRef>::const_iterator pvtrk=vtx.tracks_begin(); pvtrk!=vtx.tracks_end(); pvtrk++){
  if(!(pvtrk->isNonnull())) continue;
  vtxttrk.push_back(ttrkbuilder.build(**pvtrk));
 }
 edm::Handle<reco::BeamSpot> beamSpot;
 iEvent.getByLabel("offlineBeamSpot",beamSpot);
 unbvtx = vtx_kvf.vertex(vtxttrk, *beamSpot);
 return unbvtx;
}
//Refit the vertex using Adaptive vertex fitting (AVF) technique
TransientVertex refitted_vertex_AVF(const reco::Vertex vtx, const TransientTrackBuilder& ttrkbuilder){
 TransientVertex unbvtx;
 vector<TransientTrack> vtxttrk;
 for(vector<reco::TrackBaseRef>::const_iterator pvtrk=vtx.tracks_begin(); pvtrk!=vtx.tracks_end(); pvtrk++){
  if(!(pvtrk->isNonnull())) continue;
  vtxttrk.push_back(ttrkbuilder.build(**pvtrk));
 }
 AdaptiveVertexFitter vtx_avf;
 vtx_avf.setWeightThreshold(0.1); //weight per track. allow almost every fit, else --> exception
 if(vtxttrk.size()>=2) unbvtx = vtx_avf.vertex(vtxttrk);
 return unbvtx;
}
//Refit the vertex using Adaptive vertex fitting (AVF) technique with beam spot constraint
TransientVertex refitted_vertex_AVFbs(const Event& iEvent, const reco::Vertex vtx, const TransientTrackBuilder& ttrkbuilder){
 TransientVertex unbvtx;
 vector<TransientTrack> vtxttrk;
 for(vector<reco::TrackBaseRef>::const_iterator pvtrk=vtx.tracks_begin(); pvtrk!=vtx.tracks_end(); pvtrk++){
  if(!(pvtrk->isNonnull())) continue;
  vtxttrk.push_back(ttrkbuilder.build(**pvtrk));
 }
 AdaptiveVertexFitter vtx_avf;
 vtx_avf.setWeightThreshold(0.1); //weight per track. allow almost every fit, else --> exception
 edm::Handle<reco::BeamSpot> beamSpot;
 iEvent.getByLabel("offlineBeamSpot",beamSpot);
 if(vtxttrk.size()>=2) unbvtx = vtx_avf.vertex(vtxttrk, *beamSpot);
 return unbvtx;
}
/////
//   Get the unbiased vertex: refit the vertex using the vtx tracks without the tau tracks
/////
//Get the transient tracks collection removing tau tracks
std::vector<reco::TransientTrack> get_TransientTracks_VertexNoTauTracks(const reco::Vertex vtx, const TransientTrackBuilder& ttrkbuilder, const PFTau& pftau){
 std::vector<reco::TransientTrack> vtxttrk;
 //Get the tracks of the pftau
 std::vector<reco::TrackBaseRef> SignalTracks;
 const std::vector<edm::Ptr<reco::PFCandidate> > cands = pftau.signalPFChargedHadrCands();
 for(std::vector<edm::Ptr<reco::PFCandidate> >::const_iterator iter = cands.begin(); iter!=cands.end(); iter++){
  if(iter->get()->trackRef().isNonnull()) SignalTracks.push_back(reco::TrackBaseRef(iter->get()->trackRef()));
  else if(iter->get()->gsfTrackRef().isNonnull()){SignalTracks.push_back(reco::TrackBaseRef(((iter)->get()->gsfTrackRef())));}
 }
 //Get Vertex Tracks without tau tracks
 reco::TrackCollection nonTauTracks;
 for(std::vector<reco::TrackBaseRef>::const_iterator vtxTrkRef=vtx.tracks_begin();vtxTrkRef<vtx.tracks_end();vtxTrkRef++){
  bool nopftautrk = true;
  for(uint t=0; t<SignalTracks.size(); t++) if((*vtxTrkRef)==SignalTracks[t]) nopftautrk = false;
  if(nopftautrk) nonTauTracks.push_back(**vtxTrkRef);
 }
 //Get Transient Tracks
 for(reco::TrackCollection::iterator iter=nonTauTracks.begin(); iter!=nonTauTracks.end(); ++iter) vtxttrk.push_back(ttrkbuilder.build(*iter));
 return vtxttrk;
}
//Refit the unbiased vertex using Kalman vertex fitting (KVF) technique 
TransientVertex unbiased_vertex_KVF(const reco::Vertex vtx, const TransientTrackBuilder& ttrkbuilder, const PFTau& pftau){
 TransientVertex unbvtx;
 std::vector<reco::TransientTrack> vtxttrk = get_TransientTracks_VertexNoTauTracks(vtx,ttrkbuilder,pftau);
 if(vtxttrk.size()>=2) unbvtx = vtx_kvf.vertex(vtxttrk);
 return unbvtx;
}
//Refit the vertex using Kalman vertex fitting (KVF) technique with beam spot constraint
TransientVertex unbiased_vertex_KVFbs(const Event& iEvent, const reco::Vertex vtx, const TransientTrackBuilder& ttrkbuilder, const PFTau& pftau){
 TransientVertex unbvtx;
 std::vector<reco::TransientTrack> vtxttrk = get_TransientTracks_VertexNoTauTracks(vtx,ttrkbuilder,pftau);
 edm::Handle<reco::BeamSpot> beamSpot;
 iEvent.getByLabel("offlineBeamSpot",beamSpot);
 if(vtxttrk.size()>=2) unbvtx = vtx_kvf.vertex(vtxttrk, *beamSpot);
 return unbvtx;
}
//Refit the vertex using Adaptive vertex fitting (AVF) technique 
TransientVertex unbiased_vertex_AVF(const reco::Vertex vtx, const TransientTrackBuilder& ttrkbuilder, const PFTau& pftau){
 TransientVertex unbvtx;
 std::vector<reco::TransientTrack> vtxttrk = get_TransientTracks_VertexNoTauTracks(vtx,ttrkbuilder,pftau);
 AdaptiveVertexFitter vtx_avf;
 vtx_avf.setWeightThreshold(0.1); //weight per track. allow almost every fit, else --> exception
 if(vtxttrk.size()>=2) unbvtx = vtx_avf.vertex(vtxttrk);
 return unbvtx;
}
//Refit the vertex using Adaptive vertex fitting (AVF) technique with beam spot constraint
TransientVertex unbiased_vertex_AVFbs(const Event& iEvent, const reco::Vertex vtx, const TransientTrackBuilder& ttrkbuilder, const PFTau& pftau){
 TransientVertex unbvtx;
 std::vector<reco::TransientTrack> vtxttrk = get_TransientTracks_VertexNoTauTracks(vtx,ttrkbuilder,pftau);
 AdaptiveVertexFitter vtx_avf;
 vtx_avf.setWeightThreshold(0.1); //weight per track. allow almost every fit, else --> exception
 edm::Handle<reco::BeamSpot> beamSpot;
 iEvent.getByLabel("offlineBeamSpot",beamSpot);
 if(vtxttrk.size()>=2) unbvtx = vtx_avf.vertex(vtxttrk, *beamSpot);
 return unbvtx;
}
/////
//   Other implementations
/////
//Original implementation of the TauAlgorithm with a bug
TransientVertex unbiased_vertex_AVFbs_TauAlgo(const Event& iEvent, const reco::Vertex vtx, const TransientTrackBuilder& ttrkbuilder, const PFTau& pftau){
 TransientVertex unbvtx;
 //Get the tracks of the pftau
 std::vector<reco::TrackBaseRef> SignalTracks;
 const std::vector<edm::Ptr<reco::PFCandidate> > cands = pftau.signalPFChargedHadrCands(); 
 for(std::vector<edm::Ptr<reco::PFCandidate> >::const_iterator iter = cands.begin(); iter!=cands.end(); iter++){
  if(iter->get()->trackRef().isNonnull()) SignalTracks.push_back(reco::TrackBaseRef(iter->get()->trackRef()));
  else if(iter->get()->gsfTrackRef().isNonnull()){SignalTracks.push_back(reco::TrackBaseRef(((iter)->get()->gsfTrackRef())));}
 }
 reco::TrackCollection nonTauTracks;
 int vtxtracks = 0;
 for(std::vector<reco::TrackBaseRef>::const_iterator vtxTrkRef=vtx.tracks_begin();vtxTrkRef<vtx.tracks_end();vtxTrkRef++){
  for(unsigned int sigTrk = 0; sigTrk < SignalTracks.size(); sigTrk++){
   if((*vtxTrkRef)!=SignalTracks[sigTrk]){
    nonTauTracks.push_back(**vtxTrkRef);
   }
  }
  vtxtracks++;
 } 
 std::vector<reco::TransientTrack> vtxttrk;
 for(reco::TrackCollection::iterator iter=nonTauTracks.begin(); iter!=nonTauTracks.end(); ++iter){
  vtxttrk.push_back(ttrkbuilder.build(*iter));
 }
 //cout<<"TauAlgo"<<endl;
 //cout<<setw(30)<<SignalTracks.size()<<setw(30)<<vtxttrk.size()<<setw(30)<<vtxtracks<<endl;
 AdaptiveVertexFitter vtx_avf;
 vtx_avf.setWeightThreshold(0.1); //weight per track. allow almost every fit, else --> exception
 edm::Handle<reco::BeamSpot> beamSpot;
 iEvent.getByLabel("offlineBeamSpot",beamSpot);
 if(vtxttrk.size()>=2) unbvtx = vtx_avf.vertex(vtxttrk, *beamSpot);
 return unbvtx;
}
//Get the unbiased vertex as before, but now we want to add tau trks to the unbiased vertex 
TransientVertex unbiased_vertex_KVF_withtautrks(const reco::Vertex vtx, const TransientTrackBuilder& ttrkbuilder, const PFTau& pftau){
 TransientVertex unbvtx;
 //Get the tracks of the pftau
 vector<double> pftautrks_pt;
 const vector<reco::PFCandidatePtr>& sigpfchhadcands = pftau.signalPFChargedHadrCands();
 for(uint p=0; p<sigpfchhadcands.size(); p++){
  PFCandidatePtr cand = sigpfchhadcands[p];
  const Track* candtrk = cand->bestTrack(); 
  if(!candtrk) continue;
  pftautrks_pt.push_back(candtrk->pt());
 }
 //Get the transient track of the pv excluding the tracks of the pftau
 vector<TransientTrack> vtxttrk;
 for(vector<reco::TrackBaseRef>::const_iterator pvtrk=vtx.tracks_begin(); pvtrk!=vtx.tracks_end(); pvtrk++){
  if(!(pvtrk->isNonnull())) continue;
  bool nopftautrk = true;
  for(uint t=0; t<pftautrks_pt.size(); t++) if(pftautrks_pt[t]==(**pvtrk).pt()) nopftautrk = false;
  if(nopftautrk) vtxttrk.push_back(ttrkbuilder.build(**pvtrk));
 }
 //Add tau trks
 for(uint p=0; p<sigpfchhadcands.size(); p++){
  PFCandidatePtr cand = sigpfchhadcands[p];
  const Track* candtrk = cand->bestTrack();
  if(!candtrk) continue;
  vtxttrk.push_back(ttrkbuilder.build(*candtrk));
 } 
 if(vtxttrk.size()>=2) unbvtx = vtx_kvf.vertex(vtxttrk);
 return unbvtx;
}
//Get the tau secondary vertex, if possible
TransientVertex get_taudecvtx(const TransientTrackBuilder& ttrkbuilder, const PFTau& pftau){
 TransientVertex taudecvtx;
 //Get tau transient tracks
 vector<TransientTrack> ttrks;
 const vector<reco::PFCandidatePtr>& sigpfchhadcands = pftau.signalPFChargedHadrCands();
 for(uint p=0; p<sigpfchhadcands.size(); p++){
  PFCandidatePtr cand = sigpfchhadcands[p];
  const Track* candtrk = cand->bestTrack();
  if(!candtrk) continue;
  TransientTrack candttrk = ttrkbuilder.build(&*candtrk);
  ttrks.push_back(candttrk);
 }
 //Get tau sec vtx
 AdaptiveVertexFitter vtx_avf;
 vtx_avf.setWeightThreshold(0.1); //weight per track. allow almost every fit, else --> exception
 if(ttrks.size()>=2) taudecvtx = vtx_avf.vertex(ttrks);
 return taudecvtx;
}
