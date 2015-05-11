/**
This Macro
1. Uses a chain to merge the content of the trees that have been output by the analyzer and save it on a new rootpla 

Need to specify
0. See Declare constants
1. Do "voms-proxy-init --voms cms" if you read remote files
TTH. lcg-ls srm://stormfe1.pi.infn.it:8444/srm/managerv2?SFN=/cms/store/user/fromeo/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/crab_TTHbb_
TTH. lcg-ls srm://storm-se-01.ba.infn.it:8444/srm/managerv2?SFN=/cms/store/user/fromeo/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/crab_TTHbb_ 
TTJets. lcg-ls srm://stormfe1.pi.infn.it:8444/srm/managerv2?SFN=/cms/store/user/fromeo/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/crab_TTJetsbb_
TTJets. lcg-ls srm://storm-se-01.ba.infn.it:8444/srm/managerv2?SFN=/cms/store/user/fromeo/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/crab_TTJetsbb_ 
3. Do 0,$s/\/cms/root:\/\/xrootd-cms.infn.it\//g in the file specified in name_file
*/
/////
//   To run: root -l ChainTree.cc+  
/////
/////
//   Prepare Root and Roofit
/////
#include "TTree.h"
#include "TTreePlayer.h"
#include "TFile.h"
#include "TChain.h"
#include "TFileCollection.h"
#include <iostream>
using namespace std;
////
//   Declare constants
/////
const string tree_name    = "demo/tree"; //The name of the tree you have defined in your rootplas
//tth
const string name_file    = "Sig.txt"; //List with the path of the rootfiles
const string name_rootple = "Signal.root"; //The name of the new rootpla
/////
//   Main function
/////
void ChainTree(){
 //Put here the tree
 TTree *tree = new TTree("tree","tree");
 tree->SetMaxTreeSize(99000000000);
 //Create new file
 TFile *newfile = new TFile(name_rootple.c_str(),"recreate");
 newfile->cd();
 //Create chain, merge trees 
 TChain * chain      = new TChain(tree_name.c_str(),"");
 TFileCollection* fc = new TFileCollection("list", "list",name_file.c_str());
 chain->AddFileInfoList((TCollection*)fc->GetList());
 //Save it
 tree    = chain->CopyTree("");
 newfile = tree->GetCurrentFile();
 tree    = NULL;
 newfile->Write();
 newfile->Close();
 delete newfile; 
}
