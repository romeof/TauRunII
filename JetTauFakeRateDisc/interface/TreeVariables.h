//New class
#ifndef TREEVARIABLES
#define TREEVARIABLES
/////
//   Headers
/////
#include <TTree.h>
#include <string>
#include <map>
#include <cmath>
#include <iostream>
/////
//   Constants
/////
#define DEF_SIZE1D 100
#define DEF_VAL_INT -9999
#define DEF_VAL_FLOAT -9999.0f
#define DEF_VAL_DOUBLE -9999.0d
#define FLOAT_EPS 0.0000001f
#define DOUBLE_EPS 0.0000001d
/////
//   Functions
/////
#define INIT_1DARRAY(x,n,y) for(int i=0;i<n;i++) {x[i]=y;}
#define INIT_2DARRAY(x,n,m,y) for(int i=0;i<n;i++) { for(int j=0;j<m;j++) { x[i][j]=y; } }
inline bool is_undef(int x) { return x==DEF_VAL_INT; };
inline bool is_undef(float x) { return fabs(x-DEF_VAL_FLOAT) < FLOAT_EPS; };
inline bool is_undef(double x) { return fabs(x-DEF_VAL_DOUBLE) < DOUBLE_EPS; }
/////
//   Class declaration
/////
class CTree{
public:
 CTree(TTree* _tree) { tree = _tree; };
 TTree* tree;
 /////
 //   Helper functions for accessing branches
 /////
 template <typename T>
 T get_address(const std::string name){
  auto* br = tree->GetBranch(name.c_str());
  if(br==0){
   std::cerr << "ERROR: get_address CTree " << "branch " << name << " does not exist" << std::endl;
   throw std::exception();
  }
  auto* p = br->GetAddress();
  return reinterpret_cast<T>(p);
 }
 /////
 //   Declare variables
 /////
 //Usefull
 double duno;
 //Kinematic
 double tau_pt, tau_eta, tau_phi, tau_en, tau_ch;
 //IP
 double tau_vtxdz, tau_vtxdxy, tau_vtxdxyz;
 //pftau ch hads
 int pftauchhads_numv;
 double pftauchhads_num; 
 double pftauchhads_pt[DEF_SIZE1D], pftauchhads_eta[DEF_SIZE1D], pftauchhads_phi[DEF_SIZE1D], pftauchhads_en[DEF_SIZE1D], pftauchhads_ch[DEF_SIZE1D];
 double pftauchhads_IP3D_val[DEF_SIZE1D], pftauchhads_IP3D_err[DEF_SIZE1D], pftauchhads_IP3D_sig[DEF_SIZE1D], pftauchhads_IP2D_val[DEF_SIZE1D], pftauchhads_IP2D_err[DEF_SIZE1D], pftauchhads_IP2D_sig[DEF_SIZE1D], pftauchhads_sIP3D_val[DEF_SIZE1D], pftauchhads_sIP3D_err[DEF_SIZE1D], pftauchhads_sIP3D_sig[DEF_SIZE1D], pftauchhads_sIP2D_val[DEF_SIZE1D], pftauchhads_sIP2D_err[DEF_SIZE1D], pftauchhads_sIP2D_sig[DEF_SIZE1D];
 double pftauchhads_IP1D_val[DEF_SIZE1D], pftauchhads_IP1D_err[DEF_SIZE1D], pftauchhads_IP1D_sig[DEF_SIZE1D], pftauchhads_sIP1D_val[DEF_SIZE1D], pftauchhads_sIP1D_err[DEF_SIZE1D], pftauchhads_sIP1D_sig[DEF_SIZE1D];
 /////
 //   Initialise
 /////
 void loop_initialize(void){
  //Usefull
  duno        = 1.;
  //Kinematic
  tau_pt  = DEF_VAL_DOUBLE;
  tau_eta = DEF_VAL_DOUBLE;
  tau_phi = DEF_VAL_DOUBLE;
  tau_en  = DEF_VAL_DOUBLE;
  tau_ch  = DEF_VAL_DOUBLE;
  //IP
  tau_vtxdz   = DEF_VAL_DOUBLE;
  tau_vtxdxy  = DEF_VAL_DOUBLE;
  tau_vtxdxyz = DEF_VAL_DOUBLE;
  //pftau ch hads
  pftauchhads_numv = DEF_VAL_INT;
  pftauchhads_num  = DEF_VAL_DOUBLE;
  INIT_1DARRAY(pftauchhads_pt,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(pftauchhads_eta,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(pftauchhads_phi,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(pftauchhads_en,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(pftauchhads_ch,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(pftauchhads_IP3D_val,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(pftauchhads_IP3D_err,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(pftauchhads_IP3D_sig,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(pftauchhads_IP2D_val,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(pftauchhads_IP2D_err,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(pftauchhads_IP2D_sig,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(pftauchhads_sIP3D_val,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(pftauchhads_sIP3D_err,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(pftauchhads_sIP3D_sig,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(pftauchhads_sIP2D_val,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(pftauchhads_sIP2D_err,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(pftauchhads_sIP2D_sig,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(pftauchhads_IP1D_val,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(pftauchhads_IP1D_err,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(pftauchhads_IP1D_sig,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(pftauchhads_sIP1D_val,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(pftauchhads_sIP1D_err,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(pftauchhads_sIP1D_sig,DEF_SIZE1D,DEF_VAL_DOUBLE);
 } 
 /////
 //   Set branches
 /////
 void make_branches(void){
  //Usefull
  tree->Branch("duno", &duno, "duno/D");
  //Kinematic
  tree->Branch("tau_pt", &tau_pt, "tau_pt/D");
  tree->Branch("tau_eta", &tau_eta, "tau_eta/D");
  tree->Branch("tau_phi", &tau_phi, "tau_phi/D");
  tree->Branch("tau_en", &tau_en, "tau_en/D");
  tree->Branch("tau_ch", &tau_ch, "tau_ch/D");
  //IP
  tree->Branch("tau_vtxdz", &tau_vtxdz, "tau_vtxdz/D");
  tree->Branch("tau_vtxdxy", &tau_vtxdxy, "tau_vtxdxy/D");
  tree->Branch("tau_vtxdxyz", &tau_vtxdxyz, "tau_vtxdxyz/D");
  //pftau ch hads
  tree->Branch("pftauchhads_numv", &pftauchhads_numv, "pftauchhads_numv/I");
  tree->Branch("pftauchhads_num", &pftauchhads_num, "pftauchhads_num/D");
  tree->Branch("pftauchhads_pt", &pftauchhads_pt, "pftauchhads_pt[pftauchhads_numv]/D");
  tree->Branch("pftauchhads_eta", &pftauchhads_eta, "pftauchhads_eta[pftauchhads_numv]/D");
  tree->Branch("pftauchhads_phi", &pftauchhads_phi, "pftauchhads_phi[pftauchhads_numv]/D");
  tree->Branch("pftauchhads_en", &pftauchhads_en, "pftauchhads_en[pftauchhads_numv]/D");
  tree->Branch("pftauchhads_ch", &pftauchhads_ch, "pftauchhads_ch[pftauchhads_numv]/D");
  tree->Branch("pftauchhads_IP3D_val", &pftauchhads_IP3D_val, "pftauchhads_IP3D_val[pftauchhads_numv]/D");
  tree->Branch("pftauchhads_IP3D_err", &pftauchhads_IP3D_err, "pftauchhads_IP3D_err[pftauchhads_numv]/D");
  tree->Branch("pftauchhads_IP3D_sig", &pftauchhads_IP3D_sig, "pftauchhads_IP3D_sig[pftauchhads_numv]/D");
  tree->Branch("pftauchhads_IP2D_val", &pftauchhads_IP2D_val, "pftauchhads_IP2D_val[pftauchhads_numv]/D");
  tree->Branch("pftauchhads_IP2D_err", &pftauchhads_IP2D_err, "pftauchhads_IP2D_err[pftauchhads_numv]/D");
  tree->Branch("pftauchhads_IP2D_sig", &pftauchhads_IP2D_sig, "pftauchhads_IP2D_sig[pftauchhads_numv]/D");
  tree->Branch("pftauchhads_sIP3D_val", &pftauchhads_sIP3D_val, "pftauchhads_sIP3D_val[pftauchhads_numv]/D");
  tree->Branch("pftauchhads_sIP3D_err", &pftauchhads_sIP3D_err, "pftauchhads_sIP3D_err[pftauchhads_numv]/D");
  tree->Branch("pftauchhads_sIP3D_sig", &pftauchhads_sIP3D_sig, "pftauchhads_sIP3D_sig[pftauchhads_numv]/D");
  tree->Branch("pftauchhads_sIP2D_val", &pftauchhads_sIP2D_val, "pftauchhads_sIP2D_val[pftauchhads_numv]/D");
  tree->Branch("pftauchhads_sIP2D_err", &pftauchhads_sIP2D_err, "pftauchhads_sIP2D_err[pftauchhads_numv]/D");
  tree->Branch("pftauchhads_sIP2D_sig", &pftauchhads_sIP2D_sig, "pftauchhads_sIP2D_sig[pftauchhads_numv]/D");
  tree->Branch("pftauchhads_IP1D_val", &pftauchhads_IP1D_val, "pftauchhads_IP1D_val[pftauchhads_numv]/D");
  tree->Branch("pftauchhads_IP1D_err", &pftauchhads_IP1D_err, "pftauchhads_IP1D_err[pftauchhads_numv]/D");
  tree->Branch("pftauchhads_IP1D_sig", &pftauchhads_IP1D_sig, "pftauchhads_IP1D_sig[pftauchhads_numv]/D");
  tree->Branch("pftauchhads_sIP1D_val", &pftauchhads_sIP1D_val, "pftauchhads_sIP1D_val[pftauchhads_numv]/D");
  tree->Branch("pftauchhads_sIP1D_err", &pftauchhads_sIP1D_err, "pftauchhads_sIP1D_err[pftauchhads_numv]/D");
  tree->Branch("pftauchhads_sIP1D_sig", &pftauchhads_sIP1D_sig, "pftauchhads_sIP1D_sig[pftauchhads_numv]/D");
 }
 /////
 //   Set branch address
 /////
 //Connects the branches of an existing TTree to variables used when loading the file
 void set_branch_addresses(void){
  //Usefull
  tree->SetBranchAddress("duno", &duno);
  //Kinematic
  tree->SetBranchAddress("tau_pt", &tau_pt);
  tree->SetBranchAddress("tau_eta", &tau_eta);
  tree->SetBranchAddress("tau_phi", &tau_phi);
  tree->SetBranchAddress("tau_en", &tau_en);
  tree->SetBranchAddress("tau_ch", &tau_ch);
  //IP
  tree->SetBranchAddress("tau_vtxdz", &tau_vtxdz);
  tree->SetBranchAddress("tau_vtxdxy", &tau_vtxdxy);
  tree->SetBranchAddress("tau_vtxdxyz", &tau_vtxdxyz);
  //pftau ch hads 
  tree->SetBranchAddress("pftauchhads_numv", &pftauchhads_numv);
  tree->SetBranchAddress("pftauchhads_num", &pftauchhads_num);
  tree->SetBranchAddress("pftauchhads_pt", &pftauchhads_pt);
  tree->SetBranchAddress("pftauchhads_eta", &pftauchhads_eta);
  tree->SetBranchAddress("pftauchhads_phi", &pftauchhads_phi);
  tree->SetBranchAddress("pftauchhads_en", &pftauchhads_en);
  tree->SetBranchAddress("pftauchhads_ch", &pftauchhads_ch);
  tree->SetBranchAddress("pftauchhads_IP3D_val", &pftauchhads_IP3D_val);
  tree->SetBranchAddress("pftauchhads_IP3D_err", &pftauchhads_IP3D_err);
  tree->SetBranchAddress("pftauchhads_IP3D_sig", &pftauchhads_IP3D_sig);
  tree->SetBranchAddress("pftauchhads_IP2D_val", &pftauchhads_IP2D_val);
  tree->SetBranchAddress("pftauchhads_IP2D_err", &pftauchhads_IP2D_err);
  tree->SetBranchAddress("pftauchhads_IP2D_sig", &pftauchhads_IP2D_sig);
  tree->SetBranchAddress("pftauchhads_sIP3D_val", &pftauchhads_sIP3D_val);
  tree->SetBranchAddress("pftauchhads_sIP3D_err", &pftauchhads_sIP3D_err);
  tree->SetBranchAddress("pftauchhads_sIP3D_sig", &pftauchhads_sIP3D_sig);
  tree->SetBranchAddress("pftauchhads_sIP2D_val", &pftauchhads_sIP2D_val);
  tree->SetBranchAddress("pftauchhads_sIP2D_err", &pftauchhads_sIP2D_err);
  tree->SetBranchAddress("pftauchhads_sIP2D_sig", &pftauchhads_sIP2D_sig);
  tree->SetBranchAddress("pftauchhads_IP1D_val", &pftauchhads_IP1D_val);
  tree->SetBranchAddress("pftauchhads_IP1D_err", &pftauchhads_IP1D_err);
  tree->SetBranchAddress("pftauchhads_IP1D_sig", &pftauchhads_IP1D_sig);
  tree->SetBranchAddress("pftauchhads_sIP1D_val", &pftauchhads_sIP1D_val);
  tree->SetBranchAddress("pftauchhads_sIP1D_err", &pftauchhads_sIP1D_err);
  tree->SetBranchAddress("pftauchhads_sIP1D_sig", &pftauchhads_sIP1D_sig);
 }
};
#endif
