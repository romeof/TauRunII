//Math
#include <algorithm>
using namespace reco;
using namespace edm;
using namespace std;
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