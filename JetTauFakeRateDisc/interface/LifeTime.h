#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "RecoBTag/BTagTools/interface/SignedImpactParameter3D.h"
#include "RecoBTag/BTagTools/interface/SignedTransverseImpactParameter.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalTrajectoryExtrapolatorToLine.h"
#include "RecoBTag/SecondaryVertex/interface/SecondaryVertex.h"
//Namespaces
using namespace reco;
using namespace edm;
using namespace std;
/////
//   Implementation of IP variables considering a function for each variables (as required by Tau POG) 
/////
//The 3D IP value for the pT-leading track of the reco tauh
double IP3D_val_trkn(const PFTau& pftau, const TransientTrackBuilder& ttrkbuilder, Vertex vtx, uint n){
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
//The 3D IP significance
double IP3D_sig_trkn(const PFTau& pftau, const TransientTrackBuilder& ttrkbuilder, Vertex vtx, uint n){
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
double IP2D_val_trkn(const PFTau& pftau, const TransientTrackBuilder& ttrkbuilder, Vertex vtx, uint n){
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
double IP2D_sig_trkn(const PFTau& pftau, const TransientTrackBuilder& ttrkbuilder, Vertex vtx, uint n){
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
double IP1D_val_trkn(const PFTau& pftau, const TransientTrackBuilder& ttrkbuilder, Vertex vtx, uint n){
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
//The 1D IP significance
double IP1D_sig_trkn(const PFTau& pftau, const TransientTrackBuilder& ttrkbuilder, Vertex vtx, uint n){
 double IP = -9999;
 const vector<reco::PFCandidatePtr>& sigpfchhadcands = pftau.signalPFChargedHadrCands();
 if(sigpfchhadcands.size()<=n) return IP;
 PFCandidatePtr cand = sigpfchhadcands[n];
 const Track* trk = cand->bestTrack();
 if(!trk) return IP;
 GlobalPoint vert((vtx).position().x(), (vtx).position().y(), (vtx).position().z());
 TransientTrack ttrk = ttrkbuilder.build(&*trk);
 TrajectoryStateClosestToPoint traj = ttrk.trajectoryStateClosestToPoint(vert);
 double IP_val = fabs(traj.perigeeParameters().longitudinalImpactParameter());
 double IP_err = traj.perigeeError().longitudinalImpactParameterError();
 IP = IP_val/IP_err;
 return IP;
}
//The 3D sIP value
double sIP3D_val_trkn(const PFTau& pftau, const TransientTrackBuilder& ttrkbuilder, Vertex vtx, GlobalVector gv, uint n){
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
double sIP3D_sig_trkn(const PFTau& pftau, const TransientTrackBuilder& ttrkbuilder, Vertex vtx, GlobalVector gv, uint n){
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
double sIP2D_val_trkn(const PFTau& pftau, const TransientTrackBuilder& ttrkbuilder, Vertex vtx, GlobalVector gv, uint n){
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
double sIP2D_sig_trkn(const PFTau& pftau, const TransientTrackBuilder& ttrkbuilder, Vertex vtx, GlobalVector gv, uint n){
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
//The 1D sIP value (Our own implementation, to be checked)
double sIP1D_val_trkn(const PFTau& pftau, const TransientTrackBuilder& ttrkbuilder, Vertex vtx, GlobalVector gv, uint n){
 double IP = -9999;
 const vector<reco::PFCandidatePtr>& sigpfchhadcands = pftau.signalPFChargedHadrCands();
 if(sigpfchhadcands.size()<=n) return IP;
 PFCandidatePtr cand = sigpfchhadcands[n];
 const Track* trk = cand->bestTrack();
 if(!trk) return IP;
 GlobalPoint vert((vtx).position().x(), (vtx).position().y(), (vtx).position().z());
 TransientTrack ttrk = ttrkbuilder.build(&*trk);
 TrajectoryStateClosestToPoint traj = ttrk.trajectoryStateClosestToPoint(vert);
 double IP_abs = fabs(traj.perigeeParameters().longitudinalImpactParameter());
 AnalyticalImpactPointExtrapolator extrapolator(ttrk.field());
 TrajectoryStateOnSurface closestIn3DSpaceState = extrapolator.extrapolate(ttrk.impactPointState(),vert);
 GlobalPoint impactPoint = closestIn3DSpaceState.globalPosition();
 GlobalVector IPVec(0,0,impactPoint.z()-(vtx).position().z());
 double prod = IPVec.dot(gv);
 double sign = (prod>=0) ? 1. : -1.;
 IP = sign*IP_abs;
 return IP;
}
//The 1D sIP significance (Our own implementation, to be checked)
double sIP1D_sig_trkn(const PFTau& pftau, const TransientTrackBuilder& ttrkbuilder, Vertex vtx, GlobalVector gv, uint n){
 double IP = -9999;
 const vector<reco::PFCandidatePtr>& sigpfchhadcands = pftau.signalPFChargedHadrCands();
 if(sigpfchhadcands.size()<=n) return IP;
 PFCandidatePtr cand = sigpfchhadcands[n];
 const Track* trk = cand->bestTrack();
 if(!trk) return IP;
 GlobalPoint vert((vtx).position().x(), (vtx).position().y(), (vtx).position().z());
 TransientTrack ttrk = ttrkbuilder.build(&*trk);
 TrajectoryStateClosestToPoint traj = ttrk.trajectoryStateClosestToPoint(vert);
 double IP_val = fabs(traj.perigeeParameters().longitudinalImpactParameter());
 double IP_err = traj.perigeeError().longitudinalImpactParameterError();
 double IP_abs = IP_val/IP_err;
 AnalyticalImpactPointExtrapolator extrapolator(ttrk.field());
 TrajectoryStateOnSurface closestIn3DSpaceState = extrapolator.extrapolate(ttrk.impactPointState(),vert);
 GlobalPoint impactPoint = closestIn3DSpaceState.globalPosition();
 GlobalVector IPVec(0,0,impactPoint.z()-(vtx).position().z());
 double prod = IPVec.dot(gv);
 double sign = (prod>=0) ? 1. : -1.;
 IP = sign*IP_abs;
 return IP;
}
/////
//   Implementation of other IPTool variables considering a function for each variables (as required by Tau POG) 
/////
//The absDL3D value 
double absDL3D_val_trkn(const PFTau& pftau, const TransientTrackBuilder& ttrkbuilder, Vertex vtx, GlobalVector gv, uint n){
 double DL = -9999;
 const vector<reco::PFCandidatePtr>& sigpfchhadcands = pftau.signalPFChargedHadrCands();
 if(sigpfchhadcands.size()<=n) return DL;
 PFCandidatePtr cand = sigpfchhadcands[n];
 const Track* trk = cand->bestTrack();
 if(!trk) return DL;
 TransientTrack ttrk = ttrkbuilder.build(&*trk);
 DL = IPTools::signedDecayLength3D(ttrk,gv,vtx).second.value();
 return fabs(DL);
}
//The absDL3D significance 
double absDL3D_sig_trkn(const PFTau& pftau, const TransientTrackBuilder& ttrkbuilder, Vertex vtx, GlobalVector gv, uint n){
 double DL = -9999;
 const vector<reco::PFCandidatePtr>& sigpfchhadcands = pftau.signalPFChargedHadrCands();
 if(sigpfchhadcands.size()<=n) return DL;
 PFCandidatePtr cand = sigpfchhadcands[n];
 const Track* trk = cand->bestTrack();
 if(!trk) return DL;
 TransientTrack ttrk = ttrkbuilder.build(&*trk);
 DL = IPTools::signedDecayLength3D(ttrk,gv,vtx).second.significance();
 return fabs(DL);
}
//The sDL3D value
double sDL3D_val_trkn(const PFTau& pftau, const TransientTrackBuilder& ttrkbuilder, Vertex vtx, GlobalVector gv, uint n){
 double DL = -9999;
 const vector<reco::PFCandidatePtr>& sigpfchhadcands = pftau.signalPFChargedHadrCands();
 if(sigpfchhadcands.size()<=n) return DL;
 PFCandidatePtr cand = sigpfchhadcands[n];
 const Track* trk = cand->bestTrack();
 if(!trk) return DL;
 TransientTrack ttrk = ttrkbuilder.build(&*trk);
 DL = IPTools::signedDecayLength3D(ttrk,gv,vtx).second.value();
 return DL;
}
//The sDL3D significance
double sDL3D_sig_trkn(const PFTau& pftau, const TransientTrackBuilder& ttrkbuilder, Vertex vtx, GlobalVector gv, uint n){
 double DL = -9999;
 const vector<reco::PFCandidatePtr>& sigpfchhadcands = pftau.signalPFChargedHadrCands();
 if(sigpfchhadcands.size()<=n) return DL;
 PFCandidatePtr cand = sigpfchhadcands[n];
 const Track* trk = cand->bestTrack();
 if(!trk) return DL;
 TransientTrack ttrk = ttrkbuilder.build(&*trk);
 DL = IPTools::signedDecayLength3D(ttrk,gv,vtx).second.significance();
 return DL;
}
//The sDL2D value (Our own implementation, to be checked)
double sDL2D_val_trkn(const PFTau& pftau, const TransientTrackBuilder& ttrkbuilder, Vertex vtx, GlobalVector gv, uint n){
 double DL = -9999;
 const vector<reco::PFCandidatePtr>& sigpfchhadcands = pftau.signalPFChargedHadrCands();
 if(sigpfchhadcands.size()<=n) return DL;
 PFCandidatePtr cand = sigpfchhadcands[n];
 const Track* trk = cand->bestTrack();
 if(!trk) return DL;
 TransientTrack ttrk = ttrkbuilder.build(&*trk);
 TrajectoryStateOnSurface closestToJetState = IPTools::closestApproachToJet(ttrk.impactPointState(), vtx, gv, ttrk.field());
 if(!closestToJetState.isValid()){
  return DL;
 }
 GlobalVector jetDirection = gv.unit();
 GlobalVector jetDirection2(jetDirection.x(),jetDirection.y(),0);
 GlobalPoint vertexPosition(vtx.position().x(),vtx.position().y(),0);
 GlobalPoint Point1(closestToJetState.globalPosition().x(),closestToJetState.globalPosition().y(),0);
 DL = jetDirection2.dot(Point1-vertexPosition);
 return DL;
}
//The sDL2D significance (Our own implementation, to be checked)
double sDL2D_sig_trkn(const PFTau& pftau, const TransientTrackBuilder& ttrkbuilder, Vertex vtx, GlobalVector gv, uint n){
 double DL = -9999;
 const vector<reco::PFCandidatePtr>& sigpfchhadcands = pftau.signalPFChargedHadrCands();
 if(sigpfchhadcands.size()<=n) return DL;
 PFCandidatePtr cand = sigpfchhadcands[n];
 const Track* trk = cand->bestTrack();
 if(!trk) return DL;
 TransientTrack ttrk = ttrkbuilder.build(&*trk);
 TrajectoryStateOnSurface closestToJetState = IPTools::closestApproachToJet(ttrk.impactPointState(), vtx, gv, ttrk.field());
 if(!closestToJetState.isValid()){
  return DL;
 }
 GlobalVector jetDirection = gv.unit();
 GlobalVector jetDirection2(jetDirection.x(),jetDirection.y(),0);
 GlobalPoint vertexPosition(vtx.position().x(),vtx.position().y(),0);
 GlobalPoint Point1(closestToJetState.globalPosition().x(),closestToJetState.globalPosition().y(),0);
 double DL_val = jetDirection2.dot(Point1-vertexPosition);
 AlgebraicVector3 j;
 j[0] = jetDirection.x();
 j[1] = jetDirection.y();
 j[2] = 0.;
 AlgebraicVector6 jj;
 jj[0] = jetDirection.x();
 jj[1] = jetDirection.y();
 jj[2] = 0.;
 jj[3] = 0.;
 jj[4] = 0.;
 jj[5] = 0.;
 reco::Vertex vtx2=vtx;
 //vtx2.position().z()=0;
 //closestToJetState.globalPosition().z()=0;
 double trackError2 = ROOT::Math::Similarity(jj,closestToJetState.cartesianError().matrix());
 double vertexError2 = ROOT::Math::Similarity(j,vtx2.covariance());
 double DL_err = sqrt(trackError2+vertexError2);
 DL = DL_val/DL_err;
 return DL;
}
/////
//   Original implementation of IP variables
/////
//IP 3D and 2D sig and val, with and w/ sign
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
//IP 1D sig and val, with and w/ sign
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
void IP1D(const Track* trk,const TransientTrackBuilder& ttrkbuilder, Vertex vtx, GlobalVector gv, double& pftauchhads_IP1D_val,double& pftauchhads_IP1D_err,double& pftauchhads_IP1D_sig, double& pftauchhads_sIP1D_val,double& pftauchhads_sIP1D_err,double& pftauchhads_sIP1D_sig){
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
//Tau flight distance
double abs_pvsv_dist3d_val(const PFTau& pftau, const TransientTrackBuilder& ttrkbuilder, Vertex vtx, Vertex taudecvtx, GlobalVector gv){
 double pvsv_dist3d_val = 0;
 if(taudecvtx.isValid()) pvsv_dist3d_val = SecondaryVertex::computeDist3d(vtx,taudecvtx,gv,true).value();
 return fabs(pvsv_dist3d_val);
}
double abs_pvsv_dist2d_val(const PFTau& pftau, const TransientTrackBuilder& ttrkbuilder, Vertex vtx, Vertex taudecvtx, GlobalVector gv){
 double pvsv_dist2d_val = 0;
 if(taudecvtx.isValid()) pvsv_dist2d_val = SecondaryVertex::computeDist2d(vtx,taudecvtx,gv,true).value();
 return fabs(pvsv_dist2d_val);
}
double abs_pvsv_dist3d_sig(const PFTau& pftau, const TransientTrackBuilder& ttrkbuilder, Vertex vtx, Vertex taudecvtx, GlobalVector gv){
 double pvsv_dist3d_sig = 0;
 if(taudecvtx.isValid()) pvsv_dist3d_sig = SecondaryVertex::computeDist3d(vtx,taudecvtx,gv,true).significance();
 return fabs(pvsv_dist3d_sig);
}
double abs_pvsv_dist2d_sig(const PFTau& pftau, const TransientTrackBuilder& ttrkbuilder, Vertex vtx, Vertex taudecvtx, GlobalVector gv){
 double pvsv_dist2d_sig = 0;
 if(taudecvtx.isValid()) pvsv_dist2d_sig = SecondaryVertex::computeDist2d(vtx,taudecvtx,gv,true).significance();
 return fabs(pvsv_dist2d_sig);
}
/*
 cout<<"3DVal"<<endl;
 cout<<IPTools::signedImpactParameter3D(ttrk,gv,vtx).second.value()<<setw(20)<<IPTools::signedImpactParameter3D(ttrk,gv,vtx).second.significance()<<endl;
 SignedImpactParameter3D sip3d;
 cout<<sip3d.apply(ttrk,gv,vtx).second.value()<<setw(20)<<sip3d.apply(ttrk,gv,vtx).second.significance()<<endl;
 cout<<"2DVal"<<endl;
 cout<<IP<<setw(20)<<IPTools::signedTransverseImpactParameter(ttrk,gv,vtx).second.significance()<<endl;
 SignedTransverseImpactParameter stip;
 cout<<stip.apply(ttrk,gv,vtx).second.value()<<setw(20)<<stip.apply(ttrk,gv,vtx).second.significance()<<endl;

 cout<<"1DVal"<<endl;
 cout<<sign*IP_val<<setw(20)<<IP<<endl;
 SignedTransverseImpactParameter stip;
 cout<<stip.zImpactParameter(ttrk,gv,vtx).second.value()<<setw(20)<<stip.zImpactParameter(ttrk,gv,vtx).second.significance()<<endl;

 SignedImpactParameter3D sip3d;
 cout<<"Sum"<<endl;
 double i3 = sip3d.apply(ttrk,gv,vtx).second.value();
 double i2 = stip.apply(ttrk,gv,vtx).second.value();
 double i1 = stip.zImpactParameter(ttrk,gv,vtx).second.value();
 cout<<sqrt((i2*i2)+(i1*i1))<<setw(20)<<i3<<endl;

 cout<<"the 3 ip"<<endl;
 cout<<"3d"<<setw(20)<<sip3d.apply(ttrk,gv,vtx).second.value()<<setw(20)<<sip3d.apply(ttrk,gv,vtx).second.error()<<setw(20)<<sip3d.apply(ttrk,gv,vtx).second.significance()<<endl;
 cout<<"2d"<<setw(20)<<stip.apply(ttrk,gv,vtx).second.value()<<setw(20)<<stip.apply(ttrk,gv,vtx).second.error()<<setw(20)<<stip.apply(ttrk,gv,vtx).second.significance()<<endl;
 cout<<"1d"<<setw(20)<<stip.zImpactParameter(ttrk,gv,vtx).second.value()<<setw(20)<<stip.zImpactParameter(ttrk,gv,vtx).second.error()<<setw(20)<<stip.zImpactParameter(ttrk,gv,vtx).second.significance()<<endl;
*/
