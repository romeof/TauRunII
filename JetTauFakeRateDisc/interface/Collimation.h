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
