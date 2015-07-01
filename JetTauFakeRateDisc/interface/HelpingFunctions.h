//Math
#include <algorithm>
using namespace reco;
using namespace edm;
using namespace std;
/////
//   Matching between reco candidates and gen counterpart
/////
//RecoTauh - GenTauh matching
int MatchRecoTauhGenTauh(const PFTau& pftau, const Event& iEvent){
 int posgenmatchedcand = -1;
 math::PtEtaPhiELorentzVector recotau(pftau.pt(), pftau.eta(), pftau.phi(), pftau.energy());
 double dRmin = 0.3;
 Handle<vector<GenParticle> > genParts;
 iEvent.getByLabel("genParticles", genParts);
 for(size_t ngenPart=0; ngenPart<genParts->size(); ngenPart++){
  const GenParticle & genPart = (*genParts)[ngenPart];
  if((abs(genPart.pdgId())!=15) || (genPart.status()==3)) continue;
  math::PtEtaPhiELorentzVector gentau(genPart.pt(), genPart.eta(), genPart.phi(), genPart.energy());
  math::PtEtaPhiELorentzVector genvistau(0.,0.,0.,0.);
  int neutrinos = 0;
  for(unsigned int ndaugh=0; ndaugh<genPart.numberOfDaughters(); ndaugh++){
   const Candidate * daughter = genPart.daughter(ndaugh);
   if(abs(daughter->pdgId())==12 || abs(daughter->pdgId())==14 || abs(daughter->pdgId())==16){
    neutrinos++;
    math::PtEtaPhiELorentzVector genneu(daughter->pt(), daughter->eta(), daughter->phi(), daughter->energy());
    if(neutrinos==1) genvistau = gentau - genneu;
    if(neutrinos==2) genvistau = genvistau - genneu;
   }
  }
  double dR = ROOT::Math::VectorUtil::DeltaR(recotau,genvistau);
  if(dR<dRmin && neutrinos==1){
   dRmin = dR;
   posgenmatchedcand = ngenPart;
  }
 }//End Matching with VisGenTau
 return posgenmatchedcand;
}
//RecoTauh - GenJet matching
int MatchRecoTauhGenJet(const PFTau& pftau, const Event& iEvent){
 int posgenmatchedcand = -1;
 math::PtEtaPhiELorentzVector recotau(pftau.pt(),pftau.eta(),pftau.phi(),pftau.energy());
 edm::Handle<vector<GenJet> > genjets;
 iEvent.getByLabel("ak4GenJets", genjets);
 double dRmin = 0.3;
 for(size_t ngenjet=0; ngenjet<genjets->size(); ngenjet++){
  const GenJet & genjet = (*genjets)[ngenjet];
  math::PtEtaPhiELorentzVector genjetlv(genjet.pt(),genjet.eta(),genjet.phi(),genjet.energy());
  double dR = ROOT::Math::VectorUtil::DeltaR(recotau,genjetlv);
  if(dR<dRmin){
   dRmin = dR;
   posgenmatchedcand = ngenjet;
  }
 }
 return posgenmatchedcand;  
}
/////
//   Access the ditau vertex at gen level (from Christian Veelken)
/////
void Get_genEventVertex(const Event& iEvent, double& genditauvtx_x, double& genditauvtx_y, double& genditauvtx_z){
 Handle<vector<GenParticle> > genParts;
 iEvent.getByLabel("genParticles", genParts);
 const reco::GenParticle* genTauPlus  = 0;
 const reco::GenParticle* genTauMinus = 0;
 for ( reco::GenParticleCollection::const_iterator genParticle = genParts->begin();
       genParticle != genParts->end(); ++genParticle ) {
   if ( genParticle->pdgId() == -15 && !genTauPlus  ) genTauPlus  = &(*genParticle);
   if ( genParticle->pdgId() == +15 && !genTauMinus ) genTauMinus = &(*genParticle);
   if ( genTauPlus && genTauMinus ) break;
 }
 if ( !(genTauPlus && genTauMinus) ) {
   edm::LogError ("Tau3ProngTrackRecoAnalyzer::analyze")
     << " Failed to find generator level Tau leptons --> skipping !!";
   return;
 }
 assert(TMath::Abs(genTauPlus->vertex().x() - genTauMinus->vertex().x()) < 1.e-3);
 assert(TMath::Abs(genTauPlus->vertex().y() - genTauMinus->vertex().y()) < 1.e-3);
 assert(TMath::Abs(genTauPlus->vertex().z() - genTauMinus->vertex().z()) < 1.e-2);
 genditauvtx_x = 0.5*(genTauPlus->vertex().x() + genTauMinus->vertex().x());
 genditauvtx_y = 0.5*(genTauPlus->vertex().y() + genTauMinus->vertex().y());
 genditauvtx_z = 0.5*(genTauPlus->vertex().z() + genTauMinus->vertex().z());
}
/////
//   Track quality requirements 
/////
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
