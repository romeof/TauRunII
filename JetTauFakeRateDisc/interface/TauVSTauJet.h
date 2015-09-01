//Namespaces
using namespace reco;
using namespace edm;
using namespace std;
//Tau_trkn_pt/recojettau_pt
double ratio_tau_trkn_pt_recojettau_pt(const PFTau& pftau, uint n){
 double R = 0;
 const vector<reco::PFCandidatePtr>& sigpfchhadcands = pftau.signalPFChargedHadrCands();
 if(sigpfchhadcands.size()<=n) return R;
 PFCandidatePtr cand = sigpfchhadcands[n];
 const PFJetRef& recojettau = pftau.jetRef();
 R = cand->pt()/recojettau->pt();
 return R;
}
//Tau_trkn_en/recojettau_en
double ratio_tau_trkn_en_recojettau_en(const PFTau& pftau, uint n){
 double R = 0;
 const vector<reco::PFCandidatePtr>& sigpfchhadcands = pftau.signalPFChargedHadrCands();
 if(sigpfchhadcands.size()<=n) return R;
 PFCandidatePtr cand = sigpfchhadcands[n];
 const PFJetRef& recojettau = pftau.jetRef();
 R = cand->energy()/recojettau->energy();
 return R;
}
