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
//The 3D IP value 
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
 double IP = 0;
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
 double IP = 0;
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
//AEIP1D_val_trkn
double AEIP1D_val_trkn(const PFTau& pftau,const TransientTrackBuilder& ttrkbuilder,Vertex vtx,uint n){
 double IP = 0;
 const vector<reco::PFCandidatePtr>& sigpfchhadcands = pftau.signalPFChargedHadrCands();
 if(sigpfchhadcands.size()<=n) return IP;
 PFCandidatePtr cand = sigpfchhadcands[n];
 const Track* trk = cand->bestTrack();
 if(!trk) return IP;
 //Take the value
 TransientTrack ttrk = ttrkbuilder.build(&*trk);
 AnalyticalImpactPointExtrapolator extrapolator(ttrk.field());
 GlobalPoint vert(vtx.position().x(),vtx.position().y(),vtx.position().z());
 TrajectoryStateOnSurface tsos = extrapolator.extrapolate(ttrk.impactPointState(),vert);
 GlobalPoint refPoint          = tsos.globalPosition();
 GlobalError refPointErr       = tsos.cartesianError().position();
 GlobalPoint vertexPosition    = RecoVertex::convertPos(vtx.position());
 GlobalError vertexPositionErr = RecoVertex::convertError(vtx.error());
 GlobalVector diff = refPoint-vertexPosition;
 AlgebraicVector3 vDiff;
 vDiff[0] = 0;
 vDiff[1] = 0;
 vDiff[2] = diff.z();
 IP = sqrt(pow(vDiff[0],2)+pow(vDiff[1],2)+pow(vDiff[2],2));
 return IP;
}
//AEIP1D_sig_trkn
double AEIP1D_sig_trkn(const PFTau& pftau,const TransientTrackBuilder& ttrkbuilder,Vertex vtx,uint n){
 double IP = 0;
 const vector<reco::PFCandidatePtr>& sigpfchhadcands = pftau.signalPFChargedHadrCands();
 if(sigpfchhadcands.size()<=n) return IP;
 PFCandidatePtr cand = sigpfchhadcands[n];
 const Track* trk = cand->bestTrack();
 if(!trk) return IP;
 //Take the value
 TransientTrack ttrk = ttrkbuilder.build(&*trk);
 AnalyticalImpactPointExtrapolator extrapolator(ttrk.field());
 GlobalPoint vert(vtx.position().x(),vtx.position().y(),vtx.position().z());
 TrajectoryStateOnSurface tsos = extrapolator.extrapolate(ttrk.impactPointState(),vert);
 GlobalPoint refPoint          = tsos.globalPosition();
 GlobalError refPointErr       = tsos.cartesianError().position();
 GlobalPoint vertexPosition    = RecoVertex::convertPos(vtx.position());
 GlobalError vertexPositionErr = RecoVertex::convertError(vtx.error());
 GlobalVector diff = refPoint-vertexPosition;
 AlgebraicVector3 vDiff;
 vDiff[0] = 0;
 vDiff[1] = 0;
 vDiff[2] = diff.z();
 double IP_val = sqrt(pow(vDiff[0],2)+pow(vDiff[1],2)+pow(vDiff[2],2));
 //Error
 AlgebraicSymMatrix33 error = refPointErr.matrix()+vertexPositionErr.matrix();
 double err2    = ROOT::Math::Similarity(error,vDiff);
 double IP_err  = 0;
 if(IP_val!=0)  IP_err = sqrt(err2)/IP_val;
 IP = IP_val/IP_err;
 return IP;
}
//AEsIP1D_val_trkn
double AEsIP1D_val_trkn(const PFTau& pftau,const TransientTrackBuilder& ttrkbuilder,Vertex vtx,GlobalVector gv,uint n){
 double IP = 0;
 const vector<reco::PFCandidatePtr>& sigpfchhadcands = pftau.signalPFChargedHadrCands();
 if(sigpfchhadcands.size()<=n) return IP;
 PFCandidatePtr cand = sigpfchhadcands[n];
 const Track* trk = cand->bestTrack();
 if(!trk) return IP;
 //Take the value
 TransientTrack ttrk = ttrkbuilder.build(&*trk);
 AnalyticalImpactPointExtrapolator extrapolator(ttrk.field());
 GlobalPoint vert(vtx.position().x(),vtx.position().y(),vtx.position().z());
 TrajectoryStateOnSurface tsos = extrapolator.extrapolate(ttrk.impactPointState(),vert);
 GlobalPoint refPoint          = tsos.globalPosition();
 GlobalError refPointErr       = tsos.cartesianError().position();
 GlobalPoint vertexPosition    = RecoVertex::convertPos(vtx.position());
 GlobalError vertexPositionErr = RecoVertex::convertError(vtx.error());
 GlobalVector diff = refPoint-vertexPosition;
 AlgebraicVector3 vDiff;
 vDiff[0] = 0;
 vDiff[1] = 0;
 vDiff[2] = diff.z();
 double IP_val = sqrt(pow(vDiff[0],2)+pow(vDiff[1],2)+pow(vDiff[2],2));
 //Sign
 double IPVec_x = refPoint.x()-vtx.position().x();
 double IPVec_y = refPoint.y()-vtx.position().y();
 double IPVec_z = refPoint.z()-vtx.position().z();
 GlobalVector IPVec(IPVec_x,IPVec_y,IPVec_z);
 double prod  = IPVec.dot(gv);
 double sign  = (prod>=0) ? 1. : -1.;
 IP = sign*IP_val;
 return IP;
}
//AEsIP1D_sig_trkn
double AEsIP1D_sig_trkn(const PFTau& pftau,const TransientTrackBuilder& ttrkbuilder,Vertex vtx,GlobalVector gv,uint n){
 double IP = 0;
 const vector<reco::PFCandidatePtr>& sigpfchhadcands = pftau.signalPFChargedHadrCands();
 if(sigpfchhadcands.size()<=n) return IP;
 PFCandidatePtr cand = sigpfchhadcands[n];
 const Track* trk = cand->bestTrack();
 if(!trk) return IP;
 //Take the value
 TransientTrack ttrk = ttrkbuilder.build(&*trk);
 AnalyticalImpactPointExtrapolator extrapolator(ttrk.field());
 GlobalPoint vert(vtx.position().x(),vtx.position().y(),vtx.position().z());
 TrajectoryStateOnSurface tsos = extrapolator.extrapolate(ttrk.impactPointState(),vert);
 GlobalPoint refPoint          = tsos.globalPosition();
 GlobalError refPointErr       = tsos.cartesianError().position();
 GlobalPoint vertexPosition    = RecoVertex::convertPos(vtx.position());
 GlobalError vertexPositionErr = RecoVertex::convertError(vtx.error());
 GlobalVector diff = refPoint-vertexPosition;
 AlgebraicVector3 vDiff;
 vDiff[0] = 0;
 vDiff[1] = 0;
 vDiff[2] = diff.z();
 double IP_val = sqrt(pow(vDiff[0],2)+pow(vDiff[1],2)+pow(vDiff[2],2));
 //Error
 AlgebraicSymMatrix33 error = refPointErr.matrix()+vertexPositionErr.matrix();
 double err2    = ROOT::Math::Similarity(error,vDiff);
 double IP_err  = 0;
 if(IP_val!=0)  IP_err = sqrt(err2)/IP_val;
 double IP_sig  = IP_val/IP_err;
 //Sign
 double IPVec_x = refPoint.x()-vtx.position().x();
 double IPVec_y = refPoint.y()-vtx.position().y();
 double IPVec_z = refPoint.z()-vtx.position().z();
 GlobalVector IPVec(IPVec_x,IPVec_y,IPVec_z);
 double prod  = IPVec.dot(gv);
 double sign  = (prod>=0) ? 1. : -1.;
 IP = sign*IP_sig;
 return IP;
}
double AEsIP1D_x_val_trkn(const PFTau& pftau,const TransientTrackBuilder& ttrkbuilder,Vertex vtx,GlobalVector gv,uint n){
 double IP = 0;
 const vector<reco::PFCandidatePtr>& sigpfchhadcands = pftau.signalPFChargedHadrCands();
 if(sigpfchhadcands.size()<=n) return IP;
 PFCandidatePtr cand = sigpfchhadcands[n];
 const Track* trk = cand->bestTrack();
 if(!trk) return IP;
 //Take the value
 TransientTrack ttrk = ttrkbuilder.build(&*trk);
 AnalyticalImpactPointExtrapolator extrapolator(ttrk.field());
 GlobalPoint vert(vtx.position().x(),vtx.position().y(),vtx.position().z());
 TrajectoryStateOnSurface tsos = extrapolator.extrapolate(ttrk.impactPointState(),vert);
 GlobalPoint refPoint          = tsos.globalPosition();
 GlobalError refPointErr       = tsos.cartesianError().position();
 GlobalPoint vertexPosition    = RecoVertex::convertPos(vtx.position());
 GlobalError vertexPositionErr = RecoVertex::convertError(vtx.error());
 GlobalVector diff = refPoint-vertexPosition;
 AlgebraicVector3 vDiff;
 vDiff[0] = diff.x();
 vDiff[1] = 0;
 vDiff[2] = 0;
 double IP_val = sqrt(pow(vDiff[0],2)+pow(vDiff[1],2)+pow(vDiff[2],2));
 //Sign
 double IPVec_x = refPoint.x()-vtx.position().x();
 double IPVec_y = 0;
 double IPVec_z = 0;
 GlobalVector IPVec(IPVec_x,IPVec_y,IPVec_z);
 double prod  = IPVec.dot(gv);
 double sign  = (prod>=0) ? 1. : -1.;
 IP = sign*IP_val;
 return IP;
}
double AEsIP1D_y_val_trkn(const PFTau& pftau,const TransientTrackBuilder& ttrkbuilder,Vertex vtx,GlobalVector gv,uint n){
 double IP = 0;
 const vector<reco::PFCandidatePtr>& sigpfchhadcands = pftau.signalPFChargedHadrCands();
 if(sigpfchhadcands.size()<=n) return IP;
 PFCandidatePtr cand = sigpfchhadcands[n];
 const Track* trk = cand->bestTrack();
 if(!trk) return IP;
 //Take the value
 TransientTrack ttrk = ttrkbuilder.build(&*trk);
 AnalyticalImpactPointExtrapolator extrapolator(ttrk.field());
 GlobalPoint vert(vtx.position().x(),vtx.position().y(),vtx.position().z());
 TrajectoryStateOnSurface tsos = extrapolator.extrapolate(ttrk.impactPointState(),vert);
 GlobalPoint refPoint          = tsos.globalPosition();
 GlobalError refPointErr       = tsos.cartesianError().position();
 GlobalPoint vertexPosition    = RecoVertex::convertPos(vtx.position());
 GlobalError vertexPositionErr = RecoVertex::convertError(vtx.error());
 GlobalVector diff = refPoint-vertexPosition;
 AlgebraicVector3 vDiff;
 vDiff[0] = 0;
 vDiff[1] = diff.y();
 vDiff[2] = 0;
 double IP_val = sqrt(pow(vDiff[0],2)+pow(vDiff[1],2)+pow(vDiff[2],2));
 //Sign
 double IPVec_x = 0;
 double IPVec_y = refPoint.y()-vtx.position().y();
 double IPVec_z = 0;
 GlobalVector IPVec(IPVec_x,IPVec_y,IPVec_z);
 double prod  = IPVec.dot(gv);
 double sign  = (prod>=0) ? 1. : -1.;
 IP = sign*IP_val;
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
 TrajectoryStateOnSurface closestIn3DSpaceState = extrapolator.extrapolate(ttrk.impactPointState(),vert);
 GlobalPoint impactPoint = closestIn3DSpaceState.globalPosition();
 GlobalVector IPVec(0,0,impactPoint.z()-vtx.position().z());
 double prod = IPVec.dot(gv);
 double sign = (prod>=0) ? 1. : -1.;
 pftauchhads_sIP1D_val = sign*pftauchhads_IP1D_val; 
 pftauchhads_sIP1D_err = sign*pftauchhads_IP1D_err;
 pftauchhads_sIP1D_sig = sign*pftauchhads_IP1D_sig;
}
void zIP1D(const PFTau& pftau,const TransientTrackBuilder& ttrkbuilder,Vertex vtx,GlobalVector gv,double& pftauchhads_IP1D_val,double& pftauchhads_IP1D_err,double& pftauchhads_IP1D_sig, double& pftauchhads_sIP1D_val,double& pftauchhads_sIP1D_err,double& pftauchhads_sIP1D_sig,uint n){
 const vector<reco::PFCandidatePtr>& sigpfchhadcands = pftau.signalPFChargedHadrCands();
 if(sigpfchhadcands.size()>n){
  PFCandidatePtr cand = sigpfchhadcands[n];
  const Track* trk = cand->bestTrack();  
  if(trk){
   TransientTrack ttrk = ttrkbuilder.build(&*trk);
   SignedTransverseImpactParameter stip;
   pftauchhads_IP1D_val  = fabs(stip.zImpactParameter(ttrk,gv,vtx).second.value());  
   pftauchhads_IP1D_err  = fabs(stip.zImpactParameter(ttrk,gv,vtx).second.error());  
   pftauchhads_IP1D_sig  = fabs(stip.zImpactParameter(ttrk,gv,vtx).second.significance());  
   pftauchhads_sIP1D_val = stip.zImpactParameter(ttrk,gv,vtx).second.value();  
   pftauchhads_sIP1D_err = stip.zImpactParameter(ttrk,gv,vtx).second.error();  
   pftauchhads_sIP1D_sig = stip.zImpactParameter(ttrk,gv,vtx).second.significance();  
  }//if(trk)
 }//if(sigpfchhadcands.size()>=n)
}
void TEIP_trkn(const PFTau& pftau, const TransientTrackBuilder& ttrkbuilder,Vertex vtx,GlobalVector gv,double& pftauchhads_IP_val,double& pftauchhads_IP_err,double& pftauchhads_IP_sig, double& pftauchhads_sIP_val,double& pftauchhads_sIP_err,double& pftauchhads_sIP_sig,int dimension,uint n){
 //Take value
 const vector<reco::PFCandidatePtr>& sigpfchhadcands = pftau.signalPFChargedHadrCands();
 if(sigpfchhadcands.size()>n){
  PFCandidatePtr cand = sigpfchhadcands[n];
  const Track* trk = cand->bestTrack();  
  if(trk){
   TransientTrack ttrk = ttrkbuilder.build(&*trk);
   TransverseImpactPointExtrapolator extrapolator(ttrk.field());
   GlobalPoint vert(vtx.position().x(),vtx.position().y(),vtx.position().z());
   TrajectoryStateOnSurface tsos = extrapolator.extrapolate(ttrk.impactPointState(),vert);
   GlobalPoint refPoint          = tsos.globalPosition();
   GlobalError refPointErr       = tsos.cartesianError().position();
   GlobalPoint vertexPosition    = RecoVertex::convertPos(vtx.position());
   GlobalError vertexPositionErr = RecoVertex::convertError(vtx.error());
   //Error
   AlgebraicSymMatrix33 error = refPointErr.matrix()+vertexPositionErr.matrix();
   GlobalVector diff          = refPoint-vertexPosition;
   AlgebraicVector3 vDiff;
   double IPVec_x = 0;
   double IPVec_y = 0;
   double IPVec_z = 0;
   if(dimension==1){
    vDiff[0] = 0;
    vDiff[1] = 0;
    vDiff[2] = diff.z();
    IPVec_z = refPoint.z()-vtx.position().z();
   }else if(dimension==2){
    vDiff[0] = diff.x();
    vDiff[1] = diff.y();
    vDiff[2] = 0;
    IPVec_x = refPoint.x()-vtx.position().x();
    IPVec_y = refPoint.y()-vtx.position().y();
   }else if(dimension==3){
    vDiff[0] = diff.x();
    vDiff[1] = diff.y();
    vDiff[2] = diff.z();
    IPVec_x = refPoint.x()-vtx.position().x();
    IPVec_y = refPoint.y()-vtx.position().y();
    IPVec_z = refPoint.z()-vtx.position().z();
   }
   pftauchhads_IP_val = sqrt(pow(vDiff[0],2)+pow(vDiff[1],2)+pow(vDiff[2],2));
   double err2  = ROOT::Math::Similarity(error,vDiff);
   if(pftauchhads_IP_val!=0) pftauchhads_IP_err = sqrt(err2)/pftauchhads_IP_val;
   pftauchhads_IP_sig  = pftauchhads_IP_val/pftauchhads_IP_err;
   //Sign
   GlobalVector IPVec(IPVec_x,IPVec_y,IPVec_z);
   double prod = IPVec.dot(gv);
   double sign = (prod>=0) ? 1. : -1.;
   pftauchhads_sIP_val = sign*pftauchhads_IP_val;
   pftauchhads_sIP_err = sign*pftauchhads_IP_err;
   pftauchhads_sIP_sig = sign*pftauchhads_IP_sig;
  }//if(trk)
 }//if(sigpfchhadcands.size()>=n)
}
void AEIP_trkn(const PFTau& pftau, const TransientTrackBuilder& ttrkbuilder,Vertex vtx,GlobalVector gv,double& pftauchhads_IP_val,double& pftauchhads_IP_err,double& pftauchhads_IP_sig, double& pftauchhads_sIP_val,double& pftauchhads_sIP_err,double& pftauchhads_sIP_sig,int dimension,uint n){
 const vector<reco::PFCandidatePtr>& sigpfchhadcands = pftau.signalPFChargedHadrCands();
 if(sigpfchhadcands.size()>n){
  PFCandidatePtr cand = sigpfchhadcands[n];
  const Track* trk = cand->bestTrack();  
  if(trk){
   //Take value
   TransientTrack ttrk = ttrkbuilder.build(&*trk);
   AnalyticalImpactPointExtrapolator extrapolator(ttrk.field());
   GlobalPoint vert(vtx.position().x(),vtx.position().y(),vtx.position().z());
   TrajectoryStateOnSurface tsos = extrapolator.extrapolate(ttrk.impactPointState(),vert);
   GlobalPoint refPoint          = tsos.globalPosition();
   GlobalError refPointErr       = tsos.cartesianError().position();
   GlobalPoint vertexPosition    = RecoVertex::convertPos(vtx.position());
   GlobalError vertexPositionErr = RecoVertex::convertError(vtx.error());
   //Error
   AlgebraicSymMatrix33 error = refPointErr.matrix()+vertexPositionErr.matrix();
   GlobalVector diff          = refPoint-vertexPosition;
   AlgebraicVector3 vDiff;
   double IPVec_x = 0;
   double IPVec_y = 0;
   double IPVec_z = 0;
   if(dimension==1){
    vDiff[0] = 0;
    vDiff[1] = 0;
    vDiff[2] = diff.z();
    IPVec_z = refPoint.z()-vtx.position().z();
   }else if(dimension==2){
    vDiff[0] = diff.x();
    vDiff[1] = diff.y();
    vDiff[2] = 0;
    IPVec_x = refPoint.x()-vtx.position().x();
    IPVec_y = refPoint.y()-vtx.position().y();
   }else if(dimension==3){
    vDiff[0] = diff.x();
    vDiff[1] = diff.y();
    vDiff[2] = diff.z();
    IPVec_x = refPoint.x()-vtx.position().x();
    IPVec_y = refPoint.y()-vtx.position().y();
    IPVec_z = refPoint.z()-vtx.position().z();
   }
   pftauchhads_IP_val = sqrt(pow(vDiff[0],2)+pow(vDiff[1],2)+pow(vDiff[2],2));
   double err2  = ROOT::Math::Similarity(error,vDiff);
   if(pftauchhads_IP_val!=0) pftauchhads_IP_err = sqrt(err2)/pftauchhads_IP_val;
   pftauchhads_IP_sig  = pftauchhads_IP_val/pftauchhads_IP_err;
   //Sign
   GlobalVector IPVec(IPVec_x,IPVec_y,IPVec_z);
   double prod = IPVec.dot(gv);
   double sign = (prod>=0) ? 1. : -1.;
   pftauchhads_sIP_val = sign*pftauchhads_IP_val;
   pftauchhads_sIP_err = sign*pftauchhads_IP_err;
   pftauchhads_sIP_sig = sign*pftauchhads_IP_sig;
  }//if(trk)
 }//if(sigpfchhadcands.size()>=n)
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
