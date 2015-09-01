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
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalTrajectoryExtrapolatorToLine.h"
//Namespaces
using namespace reco;
using namespace edm;
using namespace std;
//dR between tau constituent and reco tau
double dR_tautrk_recotau(const PFTau& pftau, uint n){
 double dR = -9999;
 const vector<reco::PFCandidatePtr>& sigpfchhadcands = pftau.signalPFChargedHadrCands();
 if(sigpfchhadcands.size()<=n) return dR;
 PFCandidatePtr cand = sigpfchhadcands[n];
 dR = deltaR(cand->eta(), cand->phi(), pftau.eta(), pftau.phi());
 return dR;
}
//dR between tau constituent and reco jet tau
double dR_tautrk_recojettau(const PFTau& pftau, uint n){
 double dR = 0;
 const vector<reco::PFCandidatePtr>& sigpfchhadcands = pftau.signalPFChargedHadrCands();
 if(sigpfchhadcands.size()<=n) return dR;
 PFCandidatePtr cand = sigpfchhadcands[n];
 const PFJetRef& recojettau = pftau.jetRef();
 dR = deltaR(cand->eta(), cand->phi(), recojettau->eta(), recojettau->phi());
 return dR;
}
//dR between TauFlight_dist and JetAxis
double dR_TauFlightDist_JetAxis(const PFTau& pftau, Vertex decvtx, Vertex Pvtx, GlobalVector JetAxis){
 double dR = 0;
 const vector<reco::PFCandidatePtr>& sigpfchhadcands = pftau.signalPFChargedHadrCands();
 if(sigpfchhadcands.size()<=1) return dR;
 GlobalVector TauFlightDist(decvtx.position().x()-Pvtx.position().x(), decvtx.position().y()-Pvtx.position().y(), decvtx.position().z()-Pvtx.position().z());
 dR = deltaR(TauFlightDist, JetAxis);
 return dR;
}
//The JTD value
double JTD_val_trkn(const PFTau& pftau, const TransientTrackBuilder& ttrkbuilder, Vertex vtx, GlobalVector gv, uint n){
 double JTD = -9999;
 const vector<reco::PFCandidatePtr>& sigpfchhadcands = pftau.signalPFChargedHadrCands();
 if(sigpfchhadcands.size()<=n) return JTD;
 PFCandidatePtr cand = sigpfchhadcands[n];
 const Track* trk = cand->bestTrack();
 if(!trk) return JTD;
 TransientTrack ttrk = ttrkbuilder.build(&*trk);
 JTD = IPTools::jetTrackDistance(ttrk,gv,vtx).second.value();
 return JTD;
}
//The JTD significance
double JTD_sig_trkn(const PFTau& pftau, const TransientTrackBuilder& ttrkbuilder, Vertex vtx, GlobalVector gv, uint n){
 double JTD = -9999;
 const vector<reco::PFCandidatePtr>& sigpfchhadcands = pftau.signalPFChargedHadrCands();
 if(sigpfchhadcands.size()<=n) return JTD;
 PFCandidatePtr cand = sigpfchhadcands[n];
 const Track* trk = cand->bestTrack();
 if(!trk) return JTD;
 TransientTrack ttrk = ttrkbuilder.build(&*trk);
 JTD = IPTools::jetTrackDistance(ttrk,gv,vtx).second.significance();
 return JTD;
}
