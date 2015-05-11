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
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalTrajectoryExtrapolatorToLine.h"
KalmanVertexFitter vtxFitter(true);
//Namespaces
using namespace reco;
using namespace edm;
using namespace std;
//Access IPTool variables
//void IPToolsValues(const Track* trk,const TransientTrackBuilder& ttrkbuilder,const reco::Vertex vtx, GlobalVector gv, double& cand_IP3D_val,double& cand_IP3D_err,double& cand_IP3D_sig,  double& cand_IP2D_val,double& cand_IP2D_err,double& cand_IP2D_sig, double& cand_sIP3D_val,double& cand_sIP3D_err,double& cand_sIP3D_sig,  double& cand_sIP2D_val,double& cand_sIP2D_err,double& cand_sIP2D_sig){
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

//Get the unbiased vertex: refit the PV without the tau tracks
TransientVertex unbiased_vertex(const reco::Vertex vtx, const TransientTrackBuilder& ttrkbuilder, const PFTau& pftau){
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
 if(vtxttrk.size()>=2) unbvtx = vtxFitter.vertex(vtxttrk);
 return unbvtx;
}
//Get the refitted vertex using vtx tracks (to check the vertex is the same of the original one)
TransientVertex refitted_vertex(const reco::Vertex vtx, const TransientTrackBuilder& ttrkbuilder){
 TransientVertex unbvtx;
 vector<TransientTrack> vtxttrk;
 for(vector<reco::TrackBaseRef>::const_iterator pvtrk=vtx.tracks_begin(); pvtrk!=vtx.tracks_end(); pvtrk++){
  if(!(pvtrk->isNonnull())) continue;
  vtxttrk.push_back(ttrkbuilder.build(**pvtrk));
 }
 unbvtx = vtxFitter.vertex(vtxttrk);
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
