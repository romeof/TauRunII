//Track builder infos
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
//Access IPTool variables
//The 3D IP value for the pT-leading track of the reco tauh
double IP3D_val_trkn(const PFTau& pftau, const TransientTrackBuilder& ttrkbuilder, TransientVertex vtx, uint n){
 double IP = 0;
 const vector<reco::PFCandidatePtr>& sigpfchhadcands = pftau.signalPFChargedHadrCands();
 if(sigpfchhadcands.size()<=n) return IP;
 PFCandidatePtr cand = sigpfchhadcands[n];
 const Track* trk = cand->bestTrack();
 if(!trk) return IP;
 TransientTrack ttrk = ttrkbuilder.build(&*trk);
 IP = IPTools::absoluteImpactParameter3D(ttrk,vtx).second.value();
 return IP; 
}

/*
double IP3D_val_trkn(const PFTau& pftau, const TransientTrackBuilder& ttrkbuilder, const reco::VertexRef vtx, uint n){
 double IP = 0;
 const vector<reco::PFCandidatePtr>& sigpfchhadcands = pftau.signalPFChargedHadrCands();
 if(sigpfchhadcands.size()<=n) return IP;
 PFCandidatePtr cand = sigpfchhadcands[n];
 iEvent,vtx,ttrkbuilderconst Track* trk = cand->bestTrack();
 if(!trk) return IP;
 TransientTrack ttrk = ttrkbuilder.build(&*trk);
 IP = IPTools::absoluteImpactParameter3D(ttrk,(*vtx)).second.value();
 return IP; 
}
*/

//The 3D IP significance
double IP3D_sig_trkn(const PFTau& pftau, const TransientTrackBuilder& ttrkbuilder, TransientVertex vtx, uint n){
 double IP = 0;
 const vector<reco::PFCandidatePtr>& sigpfchhadcands = pftau.signalPFChargedHadrCands();
 if(sigpfchhadcands.size()<=n) return IP;
 PFCandidatePtr cand = sigpfchhadcands[n];
 const Track* trk = cand->bestTrack();
 if(!trk) return IP;
 TransientTrack ttrk = ttrkbuilder.build(&*trk);
 IP = IPTools::absoluteImpactParameter3D(ttrk,vtx).second.significance();
 return IP;
}

//The 2D IP value
double IP2D_val_trkn(const PFTau& pftau, const TransientTrackBuilder& ttrkbuilder, TransientVertex vtx, uint n){
 double IP = 0;
 const vector<reco::PFCandidatePtr>& sigpfchhadcands = pftau.signalPFChargedHadrCands();
 if(sigpfchhadcands.size()<=n) return IP;
 PFCandidatePtr cand = sigpfchhadcands[n];
 const Track* trk = cand->bestTrack();
 if(!trk) return IP;
 TransientTrack ttrk = ttrkbuilder.build(&*trk);
 IP = IPTools::absoluteTransverseImpactParameter(ttrk,vtx).second.value(); 
 return IP;
}

//The 2D IP significance
double IP2D_sig_trkn(const PFTau& pftau, const TransientTrackBuilder& ttrkbuilder, TransientVertex vtx, uint n){
 double IP = 0;
 const vector<reco::PFCandidatePtr>& sigpfchhadcands = pftau.signalPFChargedHadrCands();
 if(sigpfchhadcands.size()<=n) return IP;
 PFCandidatePtr cand = sigpfchhadcands[n];
 const Track* trk = cand->bestTrack();
 if(!trk) return IP;
 TransientTrack ttrk = ttrkbuilder.build(&*trk);
 IP = IPTools::absoluteTransverseImpactParameter(ttrk,vtx).second.significance();
 return IP;
}

//The 1D IP value
double IP1D_val_trkn(const PFTau& pftau, const TransientTrackBuilder& ttrkbuilder, TransientVertex vtx, uint n){
 double IP = 0;
 const vector<reco::PFCandidatePtr>& sigpfchhadcands = pftau.signalPFChargedHadrCands();
 if(sigpfchhadcands.size()<=n) return IP;
 PFCandidatePtr cand = sigpfchhadcands[n];
 const Track* trk = cand->bestTrack();
 if(!trk) return IP;
 GlobalPoint vert((vtx).position().x(), (vtx).position().y(), (vtx).position().z());
 TransientTrack ttrk = ttrkbuilder.build(&*trk);
 TrajectoryStateClosestToPoint traj = ttrk.trajectoryStateClosestToPoint(vert);
 IP = fabs(traj.perigeeParameters().longitudinalImpactParameter());
 return IP;
}

//The 3D sIP value
double sIP3D_val_trkn(const PFTau& pftau, const TransientTrackBuilder& ttrkbuilder, TransientVertex vtx, GlobalVector gv, uint n){
 double IP = 0;
 const vector<reco::PFCandidatePtr>& sigpfchhadcands = pftau.signalPFChargedHadrCands();
 if(sigpfchhadcands.size()<=n) return IP;
 PFCandidatePtr cand = sigpfchhadcands[n];
 const Track* trk = cand->bestTrack();
 if(!trk) return IP;
 TransientTrack ttrk = ttrkbuilder.build(&*trk);
 IP = IPTools::signedImpactParameter3D(ttrk,gv,vtx).second.value();
 return IP;
}

//The 3D sIP significance
double sIP3D_sig_trkn(const PFTau& pftau, const TransientTrackBuilder& ttrkbuilder, TransientVertex vtx, GlobalVector gv, uint n){
 double IP = 0;
 const vector<reco::PFCandidatePtr>& sigpfchhadcands = pftau.signalPFChargedHadrCands();
 if(sigpfchhadcands.size()<=n) return IP;
 PFCandidatePtr cand = sigpfchhadcands[n];
 const Track* trk = cand->bestTrack();
 if(!trk) return IP;
 TransientTrack ttrk = ttrkbuilder.build(&*trk);
 IP = IPTools::signedImpactParameter3D(ttrk,gv,vtx).second.significance();
 return IP;
}

//The 2D sIP value
double sIP2D_val_trkn(const PFTau& pftau, const TransientTrackBuilder& ttrkbuilder, TransientVertex vtx, GlobalVector gv, uint n){
 double IP = 0;
 const vector<reco::PFCandidatePtr>& sigpfchhadcands = pftau.signalPFChargedHadrCands();
 if(sigpfchhadcands.size()<=n) return IP;
 PFCandidatePtr cand = sigpfchhadcands[n];
 const Track* trk = cand->bestTrack();
 if(!trk) return IP;
 TransientTrack ttrk = ttrkbuilder.build(&*trk);
 IP = IPTools::signedTransverseImpactParameter(ttrk,gv,vtx).second.value();
 return IP;
} 

//The 2D sIP significance
double sIP2D_sig_trkn(const PFTau& pftau, const TransientTrackBuilder& ttrkbuilder, TransientVertex vtx, GlobalVector gv, uint n){
 double IP = 0;
 const vector<reco::PFCandidatePtr>& sigpfchhadcands = pftau.signalPFChargedHadrCands();
 if(sigpfchhadcands.size()<=n) return IP;
 PFCandidatePtr cand = sigpfchhadcands[n];
 const Track* trk = cand->bestTrack();
 if(!trk) return IP;
 TransientTrack ttrk = ttrkbuilder.build(&*trk);
 IP = IPTools::signedTransverseImpactParameter(ttrk,gv,vtx).second.significance();
 return IP;
}

void IP3D2D(const Track* trk,const TransientTrackBuilder& ttrkbuilder, TransientVertex vtx, GlobalVector gv, double& cand_IP3D_val,double& cand_IP3D_err,double& cand_IP3D_sig,  double& cand_IP2D_val,double& cand_IP2D_err,double& cand_IP2D_sig, double& cand_sIP3D_val,double& cand_sIP3D_err,double& cand_sIP3D_sig,  double& cand_sIP2D_val,double& cand_sIP2D_err,double& cand_sIP2D_sig){
 //Get transient track
 TransientTrack ttrk = ttrkbuilder.build(&*trk);
 //3D
 cand_IP3D_val = IPTools::absoluteImpactParameter3D(ttrk,vtx).second.value();
 cand_IP3D_err = IPTools::absoluteImpactParameter3D(ttrk,vtx).second.error();
 cand_IP3D_sig = IPTools::absoluteImpactParameter3D(ttrk,vtx).second.significance();
 //2D 
 cand_IP2D_val = IPTools::absoluteTransverseImpactParameter(ttrk,vtx).second.value();
 cand_IP2D_err = IPTools::absoluteTransverseImpactParameter(ttrk,vtx).second.error();
 cand_IP2D_sig = IPTools::absoluteTransverseImpactParameter(ttrk,vtx).second.significance();
 //s3D
 cand_sIP3D_val = IPTools::signedImpactParameter3D(ttrk,gv,vtx).second.value();
 cand_sIP3D_err = IPTools::signedImpactParameter3D(ttrk,gv,vtx).second.error();
 cand_sIP3D_sig = IPTools::signedImpactParameter3D(ttrk,gv,vtx).second.significance();
 //s2D 
 cand_sIP2D_val = IPTools::signedTransverseImpactParameter(ttrk,gv,vtx).second.value();
 cand_sIP2D_err = IPTools::signedTransverseImpactParameter(ttrk,gv,vtx).second.error();
 cand_sIP2D_sig = IPTools::signedTransverseImpactParameter(ttrk,gv,vtx).second.significance();
}

void IP3D2D(const Track* trk,const TransientTrackBuilder& ttrkbuilder, Vertex vtx, GlobalVector gv, double& cand_IP3D_val,double& cand_IP3D_err,double& cand_IP3D_sig,  double& cand_IP2D_val,double& cand_IP2D_err,double& cand_IP2D_sig, double& cand_sIP3D_val,double& cand_sIP3D_err,double& cand_sIP3D_sig,  double& cand_sIP2D_val,double& cand_sIP2D_err,double& cand_sIP2D_sig){
 //Get transient track
 TransientTrack ttrk = ttrkbuilder.build(&*trk);
 //3D
 cand_IP3D_val = IPTools::absoluteImpactParameter3D(ttrk,vtx).second.value();
 cand_IP3D_err = IPTools::absoluteImpactParameter3D(ttrk,vtx).second.error();
 cand_IP3D_sig = IPTools::absoluteImpactParameter3D(ttrk,vtx).second.significance();
 //2D 
 cand_IP2D_val = IPTools::absoluteTransverseImpactParameter(ttrk,vtx).second.value();
 cand_IP2D_err = IPTools::absoluteTransverseImpactParameter(ttrk,vtx).second.error();
 cand_IP2D_sig = IPTools::absoluteTransverseImpactParameter(ttrk,vtx).second.significance();
 //s3D
 cand_sIP3D_val = IPTools::signedImpactParameter3D(ttrk,gv,vtx).second.value();
 cand_sIP3D_err = IPTools::signedImpactParameter3D(ttrk,gv,vtx).second.error();
 cand_sIP3D_sig = IPTools::signedImpactParameter3D(ttrk,gv,vtx).second.significance();
 //s2D 
 cand_sIP2D_val = IPTools::signedTransverseImpactParameter(ttrk,gv,vtx).second.value();
 cand_sIP2D_err = IPTools::signedTransverseImpactParameter(ttrk,gv,vtx).second.error();
 cand_sIP2D_sig = IPTools::signedTransverseImpactParameter(ttrk,gv,vtx).second.significance();
}

void IP3D2D(const Track* trk,const TransientTrackBuilder& ttrkbuilder, const reco::VertexRef vtx, GlobalVector gv, double& cand_IP3D_val,double& cand_IP3D_err,double& cand_IP3D_sig,  double& cand_IP2D_val,double& cand_IP2D_err,double& cand_IP2D_sig, double& cand_sIP3D_val,double& cand_sIP3D_err,double& cand_sIP3D_sig,  double& cand_sIP2D_val,double& cand_sIP2D_err,double& cand_sIP2D_sig){
 //Get transient track
 TransientTrack ttrk = ttrkbuilder.build(&*trk);
 //3D
 cand_IP3D_val = IPTools::absoluteImpactParameter3D(ttrk,(*vtx)).second.value();
 cand_IP3D_err = IPTools::absoluteImpactParameter3D(ttrk,(*vtx)).second.error();
 cand_IP3D_sig = IPTools::absoluteImpactParameter3D(ttrk,(*vtx)).second.significance();
 //2D 
 cand_IP2D_val = IPTools::absoluteTransverseImpactParameter(ttrk,(*vtx)).second.value();
 cand_IP2D_err = IPTools::absoluteTransverseImpactParameter(ttrk,(*vtx)).second.error();
 cand_IP2D_sig = IPTools::absoluteTransverseImpactParameter(ttrk,(*vtx)).second.significance();
 //s3D
 cand_sIP3D_val = IPTools::signedImpactParameter3D(ttrk,gv,(*vtx)).second.value();
 cand_sIP3D_err = IPTools::signedImpactParameter3D(ttrk,gv,(*vtx)).second.error();
 cand_sIP3D_sig = IPTools::signedImpactParameter3D(ttrk,gv,(*vtx)).second.significance();
 //s2D 
 cand_sIP2D_val = IPTools::signedTransverseImpactParameter(ttrk,gv,(*vtx)).second.value();
 cand_sIP2D_err = IPTools::signedTransverseImpactParameter(ttrk,gv,(*vtx)).second.error();
 cand_sIP2D_sig = IPTools::signedTransverseImpactParameter(ttrk,gv,(*vtx)).second.significance();
}

void IP1D(const Track* trk,const TransientTrackBuilder& ttrkbuilder, const reco::VertexRef vtx, GlobalVector gv, double& pftauchhads_IP1D_val,double& pftauchhads_IP1D_err,double& pftauchhads_IP1D_sig, double& pftauchhads_sIP1D_val,double& pftauchhads_sIP1D_err,double& pftauchhads_sIP1D_sig){
 //Take value
 GlobalPoint vert((*vtx).position().x(), (*vtx).position().y(), (*vtx).position().z());
 TransientTrack ttrk = ttrkbuilder.build(&*trk);
 TrajectoryStateClosestToPoint traj = ttrk.trajectoryStateClosestToPoint(vert);
 pftauchhads_IP1D_val = fabs(traj.perigeeParameters().longitudinalImpactParameter());
 pftauchhads_IP1D_err = traj.perigeeError().longitudinalImpactParameterError(); 
 pftauchhads_IP1D_sig = pftauchhads_IP1D_val/pftauchhads_IP1D_err; 
 //Take sign
 AnalyticalImpactPointExtrapolator extrapolator(ttrk.field());
 //GlobalPoint vertexPosition = RecoVertex::convertPos((*vtx).position());
 TrajectoryStateOnSurface closestIn3DSpaceState = extrapolator.extrapolate(ttrk.impactPointState(),vert);
 GlobalPoint impactPoint = closestIn3DSpaceState.globalPosition();
 GlobalVector IPVec(0,0,impactPoint.z()-(*vtx).position().z());
 double prod = IPVec.dot(gv);
 double sign = (prod>=0) ? 1. : -1.;
 pftauchhads_sIP1D_val = sign*pftauchhads_IP1D_val; 
 pftauchhads_sIP1D_err = sign*pftauchhads_IP1D_err;
 pftauchhads_sIP1D_sig = sign*pftauchhads_IP1D_sig;
}

void IP1D(const Track* trk,const TransientTrackBuilder& ttrkbuilder, TransientVertex vtx, GlobalVector gv, double& pftauchhads_IP1D_val,double& pftauchhads_IP1D_err,double& pftauchhads_IP1D_sig, double& pftauchhads_sIP1D_val,double& pftauchhads_sIP1D_err,double& pftauchhads_sIP1D_sig){
 //Take value
 GlobalPoint vert(vtx.position().x(), vtx.position().y(), vtx.position().z());
 TransientTrack ttrk = ttrkbuilder.build(&*trk);
 TrajectoryStateClosestToPoint traj = ttrk.trajectoryStateClosestToPoint(vert);
 pftauchhads_IP1D_val = fabs(traj.perigeeParameters().longitudinalImpactParameter());
 pftauchhads_IP1D_err = traj.perigeeError().longitudinalImpactParameterError(); 
 pftauchhads_IP1D_sig = pftauchhads_IP1D_val/pftauchhads_IP1D_err; 
 //Take sign
 AnalyticalImpactPointExtrapolator extrapolator(ttrk.field());
 //GlobalPoint vertexPosition = RecoVertex::convertPos(vtx.position());
 TrajectoryStateOnSurface closestIn3DSpaceState = extrapolator.extrapolate(ttrk.impactPointState(),vert);
 GlobalPoint impactPoint = closestIn3DSpaceState.globalPosition();
 GlobalVector IPVec(0,0,impactPoint.z()-vtx.position().z());
 double prod = IPVec.dot(gv);
 double sign = (prod>=0) ? 1. : -1.;
 pftauchhads_sIP1D_val = sign*pftauchhads_IP1D_val; 
 pftauchhads_sIP1D_err = sign*pftauchhads_IP1D_err;
 pftauchhads_sIP1D_sig = sign*pftauchhads_IP1D_sig;
}

//Get the refitted vertex using vtx tracks 
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
//Get the unbiased vertex: refit the Vertex without the tau tracks
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
 int vtxtracks = 0;
 for(std::vector<reco::TrackBaseRef>::const_iterator vtxTrkRef=vtx.tracks_begin();vtxTrkRef<vtx.tracks_end();vtxTrkRef++){
  bool nopftautrk = true;
  for(uint t=0; t<SignalTracks.size(); t++) if((*vtxTrkRef)==SignalTracks[t]) nopftautrk = false;
  if(nopftautrk) nonTauTracks.push_back(**vtxTrkRef);
  vtxtracks++;
 }
 //Get Transient Tracks
 for(reco::TrackCollection::iterator iter=nonTauTracks.begin(); iter!=nonTauTracks.end(); ++iter) vtxttrk.push_back(ttrkbuilder.build(*iter));
 //cout<<"TauAlgoCorr"<<endl;
 //cout<<setw(30)<<SignalTracks.size()<<setw(30)<<vtxttrk.size()<<setw(30)<<vtxtracks<<endl;
 return vtxttrk;
}
TransientVertex unbiased_vertex_KVF(const reco::Vertex vtx, const TransientTrackBuilder& ttrkbuilder, const PFTau& pftau){
 TransientVertex unbvtx;
 std::vector<reco::TransientTrack> vtxttrk = get_TransientTracks_VertexNoTauTracks(vtx,ttrkbuilder,pftau);
 if(vtxttrk.size()>=2) unbvtx = vtx_kvf.vertex(vtxttrk);
 return unbvtx;
}
TransientVertex unbiased_vertex_KVFbs(const Event& iEvent, const reco::Vertex vtx, const TransientTrackBuilder& ttrkbuilder, const PFTau& pftau){
 TransientVertex unbvtx;
 std::vector<reco::TransientTrack> vtxttrk = get_TransientTracks_VertexNoTauTracks(vtx,ttrkbuilder,pftau);
 edm::Handle<reco::BeamSpot> beamSpot;
 iEvent.getByLabel("offlineBeamSpot",beamSpot);
 if(vtxttrk.size()>=2) unbvtx = vtx_kvf.vertex(vtxttrk, *beamSpot);
 return unbvtx;
}
TransientVertex unbiased_vertex_AVF(const reco::Vertex vtx, const TransientTrackBuilder& ttrkbuilder, const PFTau& pftau){
 TransientVertex unbvtx;
 std::vector<reco::TransientTrack> vtxttrk = get_TransientTracks_VertexNoTauTracks(vtx,ttrkbuilder,pftau);
 AdaptiveVertexFitter vtx_avf;
 vtx_avf.setWeightThreshold(0.1); //weight per track. allow almost every fit, else --> exception
 if(vtxttrk.size()>=2) unbvtx = vtx_avf.vertex(vtxttrk);
 return unbvtx;
}
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
//Check that the track is a good track
bool is_goodtrk(const Track* trk, const reco::Vertex& vtx){
 bool isgoodtrk = false;
 if(trk->pt()>1 &&
   trk->hitPattern().numberOfValidHits()>=8 &&
   trk->hitPattern().numberOfValidPixelHits()>=2 &&
   trk->normalizedChi2()<5 &&
   std::abs(trk->dxy(vtx.position()))<0.2 &&
   std::abs(trk->dz(vtx.position()))<17
   ) isgoodtrk = true;
 return isgoodtrk;
}
