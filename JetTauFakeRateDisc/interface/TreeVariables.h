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
 double numrecovtcs;
 //Vertex resolution
 double genditauvtx_x, genditauvtx_y, genditauvtx_z;
 double unbpv_KVF_nv, unbpv_KVFbs_nv, unbpv_AVF_nv, unbpv_AVFbs_nv;
 double vtxKVF_x, vtxKVF_y, vtxKVF_z, vtxKVFbs_x, vtxKVFbs_y, vtxKVFbs_z, vtxAVF_x, vtxAVF_y, vtxAVF_z, vtxAVFbs_x, vtxAVFbs_y, vtxAVFbs_z; 
 double pvtauvtx_genditauvtx_x, pvtauvtx_genditauvtx_y, pvtauvtx_genditauvtx_z, unbpv_KVF_genditauvtx_x, unbpv_KVF_genditauvtx_y, unbpv_KVF_genditauvtx_z;
 //RecoJetTau
 double recojettau_pt, recojettau_eta, recojettau_phi, recojettau_en; 
 double tau_pt_DIV_recojettau_pt, tau_en_DIV_recojettau_en, tau_trk0_pt_DIV_recojettau_pt, tau_trk1_pt_DIV_recojettau_pt, tau_trk2_pt_DIV_recojettau_pt, tau_trk0_en_DIV_recojettau_en, tau_trk1_en_DIV_recojettau_en, tau_trk2_en_DIV_recojettau_en;
 //dR
 double dR_recojettau_recotau;
 double dR_tautrk_recojettau_trk0, dR_tautrk_recojettau_trk1, dR_tautrk_recojettau_trk2;
 double dR_RecoGen;
 double dR_vtxDirGen;
 double dR_tautrk_recotau_trk0, dR_tautrk_recotau_trk1, dR_tautrk_recotau_trk2;
 double dR_TauFlightDist_JetAxis;
 //pftaugv
 double pftaugv_x;
 double pftaugv_y;
 double pftaugv_z;
 //Kinematic
 double tau_pt, tau_eta, tau_phi, tau_en, tau_ch;
 //IP
 double tau_vtxdz, tau_vtxdxy, tau_vtxdxyz;
 double pftauchhads_IP3D_val_trk0, pftauchhads_IP3D_val_trk1, pftauchhads_IP3D_val_trk2;
 double pftauchhads_IP3D_sig_trk0, pftauchhads_IP3D_sig_trk1, pftauchhads_IP3D_sig_trk2;
 double pftauchhads_IP2D_val_trk0, pftauchhads_IP2D_val_trk1, pftauchhads_IP2D_val_trk2;
 double pftauchhads_IP2D_sig_trk0, pftauchhads_IP2D_sig_trk1, pftauchhads_IP2D_sig_trk2;
 double pftauchhads_IP1D_val_trk0, pftauchhads_IP1D_val_trk1, pftauchhads_IP1D_val_trk2;
 double pftauchhads_IP1D_sig_trk0, pftauchhads_IP1D_sig_trk1, pftauchhads_IP1D_sig_trk2;
 double pftauchhads_AEIP1D_val_trk0, pftauchhads_AEIP1D_val_trk1, pftauchhads_AEIP1D_val_trk2;
 double pftauchhads_AEIP1D_sig_trk0, pftauchhads_AEIP1D_sig_trk1, pftauchhads_AEIP1D_sig_trk2;
 double pftauchhads_sIP3D_val_trk0, pftauchhads_sIP3D_val_trk1, pftauchhads_sIP3D_val_trk2;
 double pftauchhads_sIP3D_sig_trk0, pftauchhads_sIP3D_sig_trk1, pftauchhads_sIP3D_sig_trk2;
 double pftauchhads_sIP2D_val_trk0, pftauchhads_sIP2D_val_trk1, pftauchhads_sIP2D_val_trk2;
 double pftauchhads_sIP2D_sig_trk0, pftauchhads_sIP2D_sig_trk1, pftauchhads_sIP2D_sig_trk2;
 double pftauchhads_sIP1D_val_trk0, pftauchhads_sIP1D_val_trk1, pftauchhads_sIP1D_val_trk2;
 double pftauchhads_sIP1D_sig_trk0, pftauchhads_sIP1D_sig_trk1, pftauchhads_sIP1D_sig_trk2;
 double pftauchhads_AEsIP1D_val_trk0, pftauchhads_AEsIP1D_val_trk1, pftauchhads_AEsIP1D_val_trk2;
 double pftauchhads_AEsIP1D_sig_trk0, pftauchhads_AEsIP1D_sig_trk1, pftauchhads_AEsIP1D_sig_trk2;
 double pftauchhads_sDL3D_val_trk0, pftauchhads_sDL3D_val_trk1, pftauchhads_sDL3D_val_trk2;
 double pftauchhads_sDL3D_sig_trk0, pftauchhads_sDL3D_sig_trk1, pftauchhads_sDL3D_sig_trk2;
 double pftauchhads_absDL3D_val_trk0, pftauchhads_absDL3D_val_trk1, pftauchhads_absDL3D_val_trk2;
 double pftauchhads_absDL3D_sig_trk0, pftauchhads_absDL3D_sig_trk1, pftauchhads_absDL3D_sig_trk2;
 //Compare Analytical Extrapolator (AE) vs Transverse Extrapolator (TE)
 double pftauchhads_AEIP1D_x_val_trk0, pftauchhads_AEIP1D_y_val_trk0, pftauchhads_IP3DvalAE_trk0, pftauchhads_IP3DvalTE_trk0, pftauchhads_IP2DvalAE_trk0, pftauchhads_IP2DvalTE_trk0, pftauchhads_IP1DvalAE_trk0, pftauchhads_IP1DvalTE_trk0;
 double pftauchhads_IP3DvalAE_2DTE_1DAE_trk0, pftauchhads_IP3DvalAE_2DTE_1DTE_trk0, pftauchhads_IP3DvalAE_2DAE_1DAE_trk0, pftauchhads_IP3DvalTE_2DTE_1DTE_trk0;
 //Collimation
 double pftauchhads_JTD_val_trk0, pftauchhads_JTD_val_trk1, pftauchhads_JTD_val_trk2;
 double pftauchhads_JTD_sig_trk0, pftauchhads_JTD_sig_trk1, pftauchhads_JTD_sig_trk2;
 double pftau_pvsv_dist3d_val, pftau_pvsv_dist2d_val, pftau_pvsv_dist3d_sig, pftau_pvsv_dist2d_sig; 
 double pv_avfbs_nchi2, unbpv_avfbs_nchi2, diff_pv_avfbs_nchi2_unbpv_avfbs_nchi2;
 //pftau ch hads
 int pftauchhads_numv;
 double pftauchhads_num; 
 double pftauchhads_pt[DEF_SIZE1D], pftauchhads_eta[DEF_SIZE1D], pftauchhads_phi[DEF_SIZE1D], pftauchhads_en[DEF_SIZE1D], pftauchhads_ch[DEF_SIZE1D];
 double pftauchhads_IP3D_val[DEF_SIZE1D], pftauchhads_IP3D_err[DEF_SIZE1D], pftauchhads_IP3D_sig[DEF_SIZE1D], pftauchhads_IP2D_val[DEF_SIZE1D], pftauchhads_IP2D_err[DEF_SIZE1D], pftauchhads_IP2D_sig[DEF_SIZE1D], pftauchhads_sIP3D_val[DEF_SIZE1D], pftauchhads_sIP3D_err[DEF_SIZE1D], pftauchhads_sIP3D_sig[DEF_SIZE1D], pftauchhads_sIP2D_val[DEF_SIZE1D], pftauchhads_sIP2D_err[DEF_SIZE1D], pftauchhads_sIP2D_sig[DEF_SIZE1D];
 double pftauchhads_IP3Ddv_val[DEF_SIZE1D], pftauchhads_IP3Ddv_err[DEF_SIZE1D], pftauchhads_IP3Ddv_sig[DEF_SIZE1D], pftauchhads_IP2Ddv_val[DEF_SIZE1D], pftauchhads_IP2Ddv_err[DEF_SIZE1D], pftauchhads_IP2Ddv_sig[DEF_SIZE1D], pftauchhads_sIP3Ddv_val[DEF_SIZE1D], pftauchhads_sIP3Ddv_err[DEF_SIZE1D], pftauchhads_sIP3Ddv_sig[DEF_SIZE1D], pftauchhads_sIP2Ddv_val[DEF_SIZE1D], pftauchhads_sIP2Ddv_err[DEF_SIZE1D], pftauchhads_sIP2Ddv_sig[DEF_SIZE1D];
 double pftauchhads_IP1D_val[DEF_SIZE1D], pftauchhads_IP1D_err[DEF_SIZE1D], pftauchhads_IP1D_sig[DEF_SIZE1D], pftauchhads_sIP1D_val[DEF_SIZE1D], pftauchhads_sIP1D_err[DEF_SIZE1D], pftauchhads_sIP1D_sig[DEF_SIZE1D];
 //Collimation of tau tracks
 double dR_tautrk_recotau[DEF_SIZE1D];
 /////
 //   Initialise
 /////
 void loop_initialize(void){
  //Usefull
  duno        = 1.;
  numrecovtcs = DEF_VAL_DOUBLE;
  //Vertex resolution
  genditauvtx_x = DEF_VAL_DOUBLE;
  genditauvtx_y = DEF_VAL_DOUBLE;
  genditauvtx_z = DEF_VAL_DOUBLE;  
  unbpv_KVF_nv = DEF_VAL_DOUBLE;
  unbpv_KVFbs_nv = DEF_VAL_DOUBLE;
  unbpv_AVF_nv = DEF_VAL_DOUBLE;
  unbpv_AVFbs_nv = DEF_VAL_DOUBLE;
  vtxKVF_x = DEF_VAL_DOUBLE;
  vtxKVF_y = DEF_VAL_DOUBLE;
  vtxKVF_z = DEF_VAL_DOUBLE;
  vtxKVFbs_x = DEF_VAL_DOUBLE;
  vtxKVFbs_y = DEF_VAL_DOUBLE;
  vtxKVFbs_z = DEF_VAL_DOUBLE;
  vtxAVF_x = DEF_VAL_DOUBLE;
  vtxAVF_y = DEF_VAL_DOUBLE;
  vtxAVF_z = DEF_VAL_DOUBLE;
  vtxAVFbs_x = DEF_VAL_DOUBLE;
  vtxAVFbs_y = DEF_VAL_DOUBLE;
  vtxAVFbs_z = DEF_VAL_DOUBLE;
  pvtauvtx_genditauvtx_x = DEF_VAL_DOUBLE;
  pvtauvtx_genditauvtx_y = DEF_VAL_DOUBLE;
  pvtauvtx_genditauvtx_z = DEF_VAL_DOUBLE;
  unbpv_KVF_genditauvtx_x = DEF_VAL_DOUBLE;
  unbpv_KVF_genditauvtx_y = DEF_VAL_DOUBLE;
  unbpv_KVF_genditauvtx_z = DEF_VAL_DOUBLE;  
  //RecoJetTau
  recojettau_pt  = DEF_VAL_DOUBLE;
  recojettau_eta = DEF_VAL_DOUBLE;
  recojettau_phi = DEF_VAL_DOUBLE;
  recojettau_en  = DEF_VAL_DOUBLE;
  tau_pt_DIV_recojettau_pt = DEF_VAL_DOUBLE;
  tau_en_DIV_recojettau_en = DEF_VAL_DOUBLE;
  tau_trk0_pt_DIV_recojettau_pt = DEF_VAL_DOUBLE;
  tau_trk1_pt_DIV_recojettau_pt = DEF_VAL_DOUBLE;
  tau_trk2_pt_DIV_recojettau_pt = DEF_VAL_DOUBLE;
  tau_trk0_en_DIV_recojettau_en = DEF_VAL_DOUBLE;
  tau_trk1_en_DIV_recojettau_en = DEF_VAL_DOUBLE;
  tau_trk2_en_DIV_recojettau_en = DEF_VAL_DOUBLE;
  //dR
  dR_recojettau_recotau = DEF_VAL_DOUBLE;
  dR_tautrk_recojettau_trk0 = DEF_VAL_DOUBLE;
  dR_tautrk_recojettau_trk1 = DEF_VAL_DOUBLE;
  dR_tautrk_recojettau_trk2 = DEF_VAL_DOUBLE;
  dR_RecoGen = DEF_VAL_DOUBLE;
  dR_vtxDirGen = DEF_VAL_DOUBLE;
  dR_tautrk_recotau_trk0 = DEF_VAL_DOUBLE;
  dR_tautrk_recotau_trk1 = DEF_VAL_DOUBLE;
  dR_tautrk_recotau_trk2 = DEF_VAL_DOUBLE;
  dR_TauFlightDist_JetAxis = DEF_VAL_DOUBLE;
  //pftaugv
  pftaugv_x = DEF_VAL_DOUBLE;
  pftaugv_y = DEF_VAL_DOUBLE;
  pftaugv_z = DEF_VAL_DOUBLE;
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
  pftauchhads_IP3D_val_trk0 = DEF_VAL_DOUBLE;
  pftauchhads_IP3D_val_trk1 = DEF_VAL_DOUBLE;
  pftauchhads_IP3D_val_trk2 = DEF_VAL_DOUBLE;
  pftauchhads_IP3D_sig_trk0 = DEF_VAL_DOUBLE;
  pftauchhads_IP3D_sig_trk1 = DEF_VAL_DOUBLE;
  pftauchhads_IP3D_sig_trk2 = DEF_VAL_DOUBLE;
  pftauchhads_IP2D_val_trk0 = DEF_VAL_DOUBLE;
  pftauchhads_IP2D_val_trk1 = DEF_VAL_DOUBLE;
  pftauchhads_IP2D_val_trk2 = DEF_VAL_DOUBLE;
  pftauchhads_IP2D_sig_trk0 = DEF_VAL_DOUBLE;
  pftauchhads_IP2D_sig_trk1 = DEF_VAL_DOUBLE;
  pftauchhads_IP2D_sig_trk2 = DEF_VAL_DOUBLE;
  pftauchhads_IP1D_val_trk0 = DEF_VAL_DOUBLE;
  pftauchhads_IP1D_val_trk1 = DEF_VAL_DOUBLE;
  pftauchhads_IP1D_val_trk2 = DEF_VAL_DOUBLE;
  pftauchhads_IP1D_sig_trk0 = DEF_VAL_DOUBLE;
  pftauchhads_IP1D_sig_trk1 = DEF_VAL_DOUBLE;
  pftauchhads_IP1D_sig_trk2 = DEF_VAL_DOUBLE;
  pftauchhads_AEIP1D_val_trk0 = DEF_VAL_DOUBLE;
  pftauchhads_AEIP1D_val_trk1 = DEF_VAL_DOUBLE;
  pftauchhads_AEIP1D_val_trk2 = DEF_VAL_DOUBLE;
  pftauchhads_AEIP1D_sig_trk0 = DEF_VAL_DOUBLE;
  pftauchhads_AEIP1D_sig_trk1 = DEF_VAL_DOUBLE;
  pftauchhads_AEIP1D_sig_trk2 = DEF_VAL_DOUBLE;
  pftauchhads_sIP3D_val_trk0 = DEF_VAL_DOUBLE;
  pftauchhads_sIP3D_val_trk1 = DEF_VAL_DOUBLE;
  pftauchhads_sIP3D_val_trk2 = DEF_VAL_DOUBLE;
  pftauchhads_sIP3D_sig_trk0 = DEF_VAL_DOUBLE;
  pftauchhads_sIP3D_sig_trk1 = DEF_VAL_DOUBLE;
  pftauchhads_sIP3D_sig_trk2 = DEF_VAL_DOUBLE;
  pftauchhads_sIP2D_val_trk0 = DEF_VAL_DOUBLE;
  pftauchhads_sIP2D_val_trk1 = DEF_VAL_DOUBLE;
  pftauchhads_sIP2D_val_trk2 = DEF_VAL_DOUBLE;
  pftauchhads_sIP2D_sig_trk0 = DEF_VAL_DOUBLE;
  pftauchhads_sIP2D_sig_trk1 = DEF_VAL_DOUBLE;
  pftauchhads_sIP2D_sig_trk2 = DEF_VAL_DOUBLE;
  pftauchhads_sIP1D_val_trk0 = DEF_VAL_DOUBLE;
  pftauchhads_sIP1D_val_trk1 = DEF_VAL_DOUBLE;
  pftauchhads_sIP1D_val_trk2 = DEF_VAL_DOUBLE;
  pftauchhads_sIP1D_sig_trk0 = DEF_VAL_DOUBLE;
  pftauchhads_sIP1D_sig_trk1 = DEF_VAL_DOUBLE;
  pftauchhads_sIP1D_sig_trk2 = DEF_VAL_DOUBLE;
  pftauchhads_AEsIP1D_val_trk0 = DEF_VAL_DOUBLE;
  pftauchhads_AEsIP1D_val_trk1 = DEF_VAL_DOUBLE;
  pftauchhads_AEsIP1D_val_trk2 = DEF_VAL_DOUBLE;
  pftauchhads_AEsIP1D_sig_trk0 = DEF_VAL_DOUBLE;
  pftauchhads_AEsIP1D_sig_trk1 = DEF_VAL_DOUBLE;
  pftauchhads_AEsIP1D_sig_trk2 = DEF_VAL_DOUBLE;
  pftauchhads_sDL3D_val_trk0 = DEF_VAL_DOUBLE;
  pftauchhads_sDL3D_val_trk1 = DEF_VAL_DOUBLE;
  pftauchhads_sDL3D_val_trk2 = DEF_VAL_DOUBLE;
  pftauchhads_sDL3D_sig_trk0 = DEF_VAL_DOUBLE;
  pftauchhads_sDL3D_sig_trk1 = DEF_VAL_DOUBLE;
  pftauchhads_sDL3D_sig_trk2 = DEF_VAL_DOUBLE;
  pftauchhads_absDL3D_val_trk0 = DEF_VAL_DOUBLE;
  pftauchhads_absDL3D_val_trk1 = DEF_VAL_DOUBLE;
  pftauchhads_absDL3D_val_trk2 = DEF_VAL_DOUBLE;
  pftauchhads_absDL3D_sig_trk0 = DEF_VAL_DOUBLE;
  pftauchhads_absDL3D_sig_trk1 = DEF_VAL_DOUBLE;
  pftauchhads_absDL3D_sig_trk2 = DEF_VAL_DOUBLE;
  //Compare Analytical Extrapolator (AE) vs Transverse Extrapolator (TE)
  pftauchhads_AEIP1D_x_val_trk0 = DEF_VAL_DOUBLE;
  pftauchhads_AEIP1D_y_val_trk0 = DEF_VAL_DOUBLE;
  pftauchhads_IP3DvalAE_trk0 = DEF_VAL_DOUBLE;
  pftauchhads_IP3DvalTE_trk0 = DEF_VAL_DOUBLE;
  pftauchhads_IP2DvalAE_trk0 = DEF_VAL_DOUBLE;
  pftauchhads_IP2DvalTE_trk0 = DEF_VAL_DOUBLE;
  pftauchhads_IP1DvalAE_trk0 = DEF_VAL_DOUBLE;
  pftauchhads_IP1DvalTE_trk0 = DEF_VAL_DOUBLE;
  pftauchhads_IP3DvalAE_2DTE_1DAE_trk0 = DEF_VAL_DOUBLE;
  pftauchhads_IP3DvalAE_2DTE_1DTE_trk0 = DEF_VAL_DOUBLE;
  pftauchhads_IP3DvalAE_2DAE_1DAE_trk0 = DEF_VAL_DOUBLE;
  pftauchhads_IP3DvalTE_2DTE_1DTE_trk0 = DEF_VAL_DOUBLE;
  //Collimation
  pftauchhads_JTD_val_trk0 = DEF_VAL_DOUBLE;
  pftauchhads_JTD_val_trk1 = DEF_VAL_DOUBLE;
  pftauchhads_JTD_val_trk2 = DEF_VAL_DOUBLE;
  pftauchhads_JTD_sig_trk0 = DEF_VAL_DOUBLE;
  pftauchhads_JTD_sig_trk1 = DEF_VAL_DOUBLE;
  pftauchhads_JTD_sig_trk2 = DEF_VAL_DOUBLE;
  pftau_pvsv_dist3d_val = DEF_VAL_DOUBLE;
  pftau_pvsv_dist2d_val = DEF_VAL_DOUBLE;
  pftau_pvsv_dist3d_sig = DEF_VAL_DOUBLE;
  pftau_pvsv_dist2d_sig = DEF_VAL_DOUBLE;
  pv_avfbs_nchi2                        = DEF_VAL_DOUBLE;
  unbpv_avfbs_nchi2                     = DEF_VAL_DOUBLE;
  diff_pv_avfbs_nchi2_unbpv_avfbs_nchi2 = DEF_VAL_DOUBLE;
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
  INIT_1DARRAY(pftauchhads_IP3Ddv_val,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(pftauchhads_IP3Ddv_err,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(pftauchhads_IP3Ddv_sig,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(pftauchhads_IP2Ddv_val,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(pftauchhads_IP2Ddv_err,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(pftauchhads_IP2Ddv_sig,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(pftauchhads_sIP3Ddv_val,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(pftauchhads_sIP3Ddv_err,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(pftauchhads_sIP3Ddv_sig,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(pftauchhads_sIP2Ddv_val,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(pftauchhads_sIP2Ddv_err,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(pftauchhads_sIP2Ddv_sig,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(pftauchhads_IP1D_val,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(pftauchhads_IP1D_err,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(pftauchhads_IP1D_sig,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(pftauchhads_sIP1D_val,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(pftauchhads_sIP1D_err,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(pftauchhads_sIP1D_sig,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(dR_tautrk_recotau,DEF_SIZE1D,DEF_VAL_DOUBLE);
 } 
 /////
 //   Set branches
 /////
 void make_branches(void){
  //Usefull
  tree->Branch("duno", &duno, "duno/D");
  tree->Branch("numrecovtcs", &numrecovtcs, "numrecovtcs/D");
  //Vertex resolution
  tree->Branch("genditauvtx_x", &genditauvtx_x, "genditauvtx_x/D");
  tree->Branch("genditauvtx_y", &genditauvtx_y, "genditauvtx_y/D");
  tree->Branch("genditauvtx_z", &genditauvtx_z, "genditauvtx_z/D");
  tree->Branch("unbpv_KVF_nv", &unbpv_KVF_nv, "unbpv_KVF_nv/D");
  tree->Branch("unbpv_KVFbs_nv", &unbpv_KVFbs_nv, "unbpv_KVFbs_nv/D");
  tree->Branch("unbpv_AVF_nv", &unbpv_AVF_nv, "unbpv_AVF_nv/D");
  tree->Branch("unbpv_AVFbs_nv", &unbpv_AVFbs_nv, "unbpv_AVFbs_nv/D");
  tree->Branch("vtxKVF_x", &vtxKVF_x, "vtxKVF_x/D");
  tree->Branch("vtxKVF_y", &vtxKVF_y, "vtxKVF_y/D");
  tree->Branch("vtxKVF_z", &vtxKVF_z, "vtxKVF_z/D");
  tree->Branch("vtxKVFbs_x", &vtxKVFbs_x, "vtxKVFbs_x/D");
  tree->Branch("vtxKVFbs_y", &vtxKVFbs_y, "vtxKVFbs_y/D");
  tree->Branch("vtxKVFbs_z", &vtxKVFbs_z, "vtxKVFbs_z/D");
  tree->Branch("vtxAVF_x", &vtxAVF_x, "vtxAVF_x/D");
  tree->Branch("vtxAVF_y", &vtxAVF_y, "vtxAVF_y/D");
  tree->Branch("vtxAVF_z", &vtxAVF_z, "vtxAVF_z/D");
  tree->Branch("vtxAVFbs_x", &vtxAVFbs_x, "vtxAVFbs_x/D");
  tree->Branch("vtxAVFbs_y", &vtxAVFbs_y, "vtxAVFbs_y/D");
  tree->Branch("vtxAVFbs_z", &vtxAVFbs_z, "vtxAVFbs_z/D");
  tree->Branch("pvtauvtx_genditauvtx_x", &pvtauvtx_genditauvtx_x, "pvtauvtx_genditauvtx_x/D");
  tree->Branch("pvtauvtx_genditauvtx_y", &pvtauvtx_genditauvtx_y, "pvtauvtx_genditauvtx_y/D");
  tree->Branch("pvtauvtx_genditauvtx_z", &pvtauvtx_genditauvtx_z, "pvtauvtx_genditauvtx_z/D");
  tree->Branch("unbpv_KVF_genditauvtx_x", &unbpv_KVF_genditauvtx_x, "unbpv_KVF_genditauvtx_x/D");
  tree->Branch("unbpv_KVF_genditauvtx_y", &unbpv_KVF_genditauvtx_y, "unbpv_KVF_genditauvtx_y/D");
  tree->Branch("unbpv_KVF_genditauvtx_z", &unbpv_KVF_genditauvtx_z, "unbpv_KVF_genditauvtx_z/D");
  //RecoJetTau
  tree->Branch("recojettau_pt", &recojettau_pt, "recojettau_pt/D");
  tree->Branch("recojettau_eta", &recojettau_eta, "recojettau_eta/D");
  tree->Branch("recojettau_phi", &recojettau_phi, "recojettau_phi/D");
  tree->Branch("recojettau_en", &recojettau_en, "recojettau_en/D");
  tree->Branch("tau_pt_DIV_recojettau_pt", &tau_pt_DIV_recojettau_pt, "tau_pt_DIV_recojettau_pt/D");
  tree->Branch("tau_en_DIV_recojettau_en", &tau_en_DIV_recojettau_en, "tau_en_DIV_recojettau_en/D");
  tree->Branch("tau_trk0_pt_DIV_recojettau_pt", &tau_trk0_pt_DIV_recojettau_pt, "tau_trk0_pt_DIV_recojettau_pt/D");
  tree->Branch("tau_trk1_pt_DIV_recojettau_pt", &tau_trk1_pt_DIV_recojettau_pt, "tau_trk1_pt_DIV_recojettau_pt/D");
  tree->Branch("tau_trk2_pt_DIV_recojettau_pt", &tau_trk2_pt_DIV_recojettau_pt, "tau_trk2_pt_DIV_recojettau_pt/D");
  tree->Branch("tau_trk0_en_DIV_recojettau_en", &tau_trk0_en_DIV_recojettau_en, "tau_trk0_en_DIV_recojettau_en/D");
  tree->Branch("tau_trk1_en_DIV_recojettau_en", &tau_trk1_en_DIV_recojettau_en, "tau_trk1_en_DIV_recojettau_en/D");
  tree->Branch("tau_trk2_en_DIV_recojettau_en", &tau_trk2_en_DIV_recojettau_en, "tau_trk2_en_DIV_recojettau_en/D");
  //dR
  tree->Branch("dR_recojettau_recotau", &dR_recojettau_recotau, "dR_recojettau_recotau/D");
  tree->Branch("dR_tautrk_recojettau_trk0", &dR_tautrk_recojettau_trk0, "dR_tautrk_recojettau_trk0/D");
  tree->Branch("dR_tautrk_recojettau_trk1", &dR_tautrk_recojettau_trk1, "dR_tautrk_recojettau_trk1/D");
  tree->Branch("dR_tautrk_recojettau_trk2", &dR_tautrk_recojettau_trk2, "dR_tautrk_recojettau_trk2/D");
  tree->Branch("dR_RecoGen", &dR_RecoGen, "dR_RecoGen/D");
  tree->Branch("dR_vtxDirGen", &dR_vtxDirGen, "dR_vtxDirGen/D");
  tree->Branch("dR_tautrk_recotau_trk0", &dR_tautrk_recotau_trk0, "dR_tautrk_recotau_trk0/D");
  tree->Branch("dR_tautrk_recotau_trk1", &dR_tautrk_recotau_trk1, "dR_tautrk_recotau_trk1/D");
  tree->Branch("dR_tautrk_recotau_trk2", &dR_tautrk_recotau_trk2, "dR_tautrk_recotau_trk2/D");
  tree->Branch("dR_TauFlightDist_JetAxis", &dR_TauFlightDist_JetAxis, "dR_TauFlightDist_JetAxis/D");
  //pftaugv
  tree->Branch("pftaugv_x", &pftaugv_x, "pftaugv_x/D");
  tree->Branch("pftaugv_y", &pftaugv_y, "pftaugv_y/D");
  tree->Branch("pftaugv_z", &pftaugv_z, "pftaugv_z/D");
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
  tree->Branch("pftauchhads_IP3D_val_trk0", &pftauchhads_IP3D_val_trk0, "pftauchhads_IP3D_val_trk0/D");
  tree->Branch("pftauchhads_IP3D_val_trk1", &pftauchhads_IP3D_val_trk1, "pftauchhads_IP3D_val_trk1/D");
  tree->Branch("pftauchhads_IP3D_val_trk2", &pftauchhads_IP3D_val_trk2, "pftauchhads_IP3D_val_trk2/D");
  tree->Branch("pftauchhads_IP3D_sig_trk0", &pftauchhads_IP3D_sig_trk0, "pftauchhads_IP3D_sig_trk0/D");
  tree->Branch("pftauchhads_IP3D_sig_trk1", &pftauchhads_IP3D_sig_trk1, "pftauchhads_IP3D_sig_trk1/D");
  tree->Branch("pftauchhads_IP3D_sig_trk2", &pftauchhads_IP3D_sig_trk2, "pftauchhads_IP3D_sig_trk2/D");
  tree->Branch("pftauchhads_IP2D_val_trk0", &pftauchhads_IP2D_val_trk0, "pftauchhads_IP2D_val_trk0/D");
  tree->Branch("pftauchhads_IP2D_val_trk1", &pftauchhads_IP2D_val_trk1, "pftauchhads_IP2D_val_trk1/D");
  tree->Branch("pftauchhads_IP2D_val_trk2", &pftauchhads_IP2D_val_trk2, "pftauchhads_IP2D_val_trk2/D");
  tree->Branch("pftauchhads_IP2D_sig_trk0", &pftauchhads_IP2D_sig_trk0, "pftauchhads_IP2D_sig_trk0/D");
  tree->Branch("pftauchhads_IP2D_sig_trk1", &pftauchhads_IP2D_sig_trk1, "pftauchhads_IP2D_sig_trk1/D");
  tree->Branch("pftauchhads_IP2D_sig_trk2", &pftauchhads_IP2D_sig_trk2, "pftauchhads_IP2D_sig_trk2/D");
  tree->Branch("pftauchhads_IP1D_val_trk0", &pftauchhads_IP1D_val_trk0, "pftauchhads_IP1D_val_trk0/D");
  tree->Branch("pftauchhads_IP1D_val_trk1", &pftauchhads_IP1D_val_trk1, "pftauchhads_IP1D_val_trk1/D");
  tree->Branch("pftauchhads_IP1D_val_trk2", &pftauchhads_IP1D_val_trk2, "pftauchhads_IP1D_val_trk2/D");
  tree->Branch("pftauchhads_IP1D_sig_trk0", &pftauchhads_IP1D_sig_trk0, "pftauchhads_IP1D_sig_trk0/D");
  tree->Branch("pftauchhads_IP1D_sig_trk1", &pftauchhads_IP1D_sig_trk1, "pftauchhads_IP1D_sig_trk1/D");
  tree->Branch("pftauchhads_IP1D_sig_trk2", &pftauchhads_IP1D_sig_trk2, "pftauchhads_IP1D_sig_trk2/D");
  tree->Branch("pftauchhads_AEIP1D_val_trk0", &pftauchhads_AEIP1D_val_trk0, "pftauchhads_AEIP1D_val_trk0/D");
  tree->Branch("pftauchhads_AEIP1D_val_trk1", &pftauchhads_AEIP1D_val_trk1, "pftauchhads_AEIP1D_val_trk1/D");
  tree->Branch("pftauchhads_AEIP1D_val_trk2", &pftauchhads_AEIP1D_val_trk2, "pftauchhads_AEIP1D_val_trk2/D");
  tree->Branch("pftauchhads_AEIP1D_sig_trk0", &pftauchhads_AEIP1D_sig_trk0, "pftauchhads_AEIP1D_sig_trk0/D");
  tree->Branch("pftauchhads_AEIP1D_sig_trk1", &pftauchhads_AEIP1D_sig_trk1, "pftauchhads_AEIP1D_sig_trk1/D");
  tree->Branch("pftauchhads_AEIP1D_sig_trk2", &pftauchhads_AEIP1D_sig_trk2, "pftauchhads_AEIP1D_sig_trk2/D");
  tree->Branch("pftauchhads_sIP3D_val_trk0", &pftauchhads_sIP3D_val_trk0, "pftauchhads_sIP3D_val_trk0/D");
  tree->Branch("pftauchhads_sIP3D_val_trk1", &pftauchhads_sIP3D_val_trk1, "pftauchhads_sIP3D_val_trk1/D");
  tree->Branch("pftauchhads_sIP3D_val_trk2", &pftauchhads_sIP3D_val_trk2, "pftauchhads_sIP3D_val_trk2/D");
  tree->Branch("pftauchhads_sIP3D_sig_trk0", &pftauchhads_sIP3D_sig_trk0, "pftauchhads_sIP3D_sig_trk0/D");
  tree->Branch("pftauchhads_sIP3D_sig_trk1", &pftauchhads_sIP3D_sig_trk1, "pftauchhads_sIP3D_sig_trk1/D");
  tree->Branch("pftauchhads_sIP3D_sig_trk2", &pftauchhads_sIP3D_sig_trk2, "pftauchhads_sIP3D_sig_trk2/D");
  tree->Branch("pftauchhads_sIP2D_val_trk0", &pftauchhads_sIP2D_val_trk0, "pftauchhads_sIP2D_val_trk0/D");
  tree->Branch("pftauchhads_sIP2D_val_trk1", &pftauchhads_sIP2D_val_trk1, "pftauchhads_sIP2D_val_trk1/D");
  tree->Branch("pftauchhads_sIP2D_val_trk2", &pftauchhads_sIP2D_val_trk2, "pftauchhads_sIP2D_val_trk2/D");
  tree->Branch("pftauchhads_sIP2D_sig_trk0", &pftauchhads_sIP2D_sig_trk0, "pftauchhads_sIP2D_sig_trk0/D");
  tree->Branch("pftauchhads_sIP2D_sig_trk1", &pftauchhads_sIP2D_sig_trk1, "pftauchhads_sIP2D_sig_trk1/D");
  tree->Branch("pftauchhads_sIP2D_sig_trk2", &pftauchhads_sIP2D_sig_trk2, "pftauchhads_sIP2D_sig_trk2/D");
  tree->Branch("pftauchhads_sIP1D_val_trk0", &pftauchhads_sIP1D_val_trk0, "pftauchhads_sIP1D_val_trk0/D");
  tree->Branch("pftauchhads_sIP1D_val_trk1", &pftauchhads_sIP1D_val_trk1, "pftauchhads_sIP1D_val_trk1/D");
  tree->Branch("pftauchhads_sIP1D_val_trk2", &pftauchhads_sIP1D_val_trk2, "pftauchhads_sIP1D_val_trk2/D");
  tree->Branch("pftauchhads_sIP1D_sig_trk0", &pftauchhads_sIP1D_sig_trk0, "pftauchhads_sIP1D_sig_trk0/D");
  tree->Branch("pftauchhads_sIP1D_sig_trk1", &pftauchhads_sIP1D_sig_trk1, "pftauchhads_sIP1D_sig_trk1/D");
  tree->Branch("pftauchhads_sIP1D_sig_trk2", &pftauchhads_sIP1D_sig_trk2, "pftauchhads_sIP1D_sig_trk2/D");
  tree->Branch("pftauchhads_AEsIP1D_val_trk0", &pftauchhads_AEsIP1D_val_trk0, "pftauchhads_AEsIP1D_val_trk0/D");
  tree->Branch("pftauchhads_AEsIP1D_val_trk1", &pftauchhads_AEsIP1D_val_trk1, "pftauchhads_AEsIP1D_val_trk1/D");
  tree->Branch("pftauchhads_AEsIP1D_val_trk2", &pftauchhads_AEsIP1D_val_trk2, "pftauchhads_AEsIP1D_val_trk2/D");
  tree->Branch("pftauchhads_AEsIP1D_sig_trk0", &pftauchhads_AEsIP1D_sig_trk0, "pftauchhads_AEsIP1D_sig_trk0/D");
  tree->Branch("pftauchhads_AEsIP1D_sig_trk1", &pftauchhads_AEsIP1D_sig_trk1, "pftauchhads_AEsIP1D_sig_trk1/D");
  tree->Branch("pftauchhads_AEsIP1D_sig_trk2", &pftauchhads_AEsIP1D_sig_trk2, "pftauchhads_AEsIP1D_sig_trk2/D");
  tree->Branch("pftauchhads_sDL3D_val_trk0", &pftauchhads_sDL3D_val_trk0, "pftauchhads_sDL3D_val_trk0/D");
  tree->Branch("pftauchhads_sDL3D_val_trk1", &pftauchhads_sDL3D_val_trk1, "pftauchhads_sDL3D_val_trk1/D");
  tree->Branch("pftauchhads_sDL3D_val_trk2", &pftauchhads_sDL3D_val_trk2, "pftauchhads_sDL3D_val_trk2/D");
  tree->Branch("pftauchhads_sDL3D_sig_trk0", &pftauchhads_sDL3D_sig_trk0, "pftauchhads_sDL3D_sig_trk0/D");
  tree->Branch("pftauchhads_sDL3D_sig_trk1", &pftauchhads_sDL3D_sig_trk1, "pftauchhads_sDL3D_sig_trk1/D");
  tree->Branch("pftauchhads_sDL3D_sig_trk2", &pftauchhads_sDL3D_sig_trk2, "pftauchhads_sDL3D_sig_trk2/D");
  tree->Branch("pftauchhads_absDL3D_val_trk0", &pftauchhads_absDL3D_val_trk0, "pftauchhads_absDL3D_val_trk0/D");
  tree->Branch("pftauchhads_absDL3D_val_trk1", &pftauchhads_absDL3D_val_trk1, "pftauchhads_absDL3D_val_trk1/D");
  tree->Branch("pftauchhads_absDL3D_val_trk2", &pftauchhads_absDL3D_val_trk2, "pftauchhads_absDL3D_val_trk2/D");
  tree->Branch("pftauchhads_absDL3D_sig_trk0", &pftauchhads_absDL3D_sig_trk0, "pftauchhads_absDL3D_sig_trk0/D");
  tree->Branch("pftauchhads_absDL3D_sig_trk1", &pftauchhads_absDL3D_sig_trk1, "pftauchhads_absDL3D_sig_trk1/D");
  tree->Branch("pftauchhads_absDL3D_sig_trk2", &pftauchhads_absDL3D_sig_trk2, "pftauchhads_absDL3D_sig_trk2/D");
  //Compare Analytical Extrapolator (AE) vs Transverse Extrapolator (TE)
  tree->Branch("pftauchhads_AEIP1D_x_val_trk0", &pftauchhads_AEIP1D_x_val_trk0, "pftauchhads_AEIP1D_x_val_trk0/D");
  tree->Branch("pftauchhads_AEIP1D_y_val_trk0", &pftauchhads_AEIP1D_y_val_trk0, "pftauchhads_AEIP1D_y_val_trk0/D");
  tree->Branch("pftauchhads_IP3DvalAE_trk0", &pftauchhads_IP3DvalAE_trk0, "pftauchhads_IP3DvalAE_trk0/D");
  tree->Branch("pftauchhads_IP3DvalTE_trk0", &pftauchhads_IP3DvalTE_trk0, "pftauchhads_IP3DvalTE_trk0/D");
  tree->Branch("pftauchhads_IP2DvalAE_trk0", &pftauchhads_IP2DvalAE_trk0, "pftauchhads_IP2DvalAE_trk0/D");
  tree->Branch("pftauchhads_IP2DvalTE_trk0", &pftauchhads_IP2DvalTE_trk0, "pftauchhads_IP2DvalTE_trk0/D");
  tree->Branch("pftauchhads_IP1DvalAE_trk0", &pftauchhads_IP1DvalAE_trk0, "pftauchhads_IP1DvalAE_trk0/D");
  tree->Branch("pftauchhads_IP1DvalTE_trk0", &pftauchhads_IP1DvalTE_trk0, "pftauchhads_IP1DvalTE_trk0/D");
  tree->Branch("pftauchhads_IP3DvalAE_2DTE_1DAE_trk0", &pftauchhads_IP3DvalAE_2DTE_1DAE_trk0, "pftauchhads_IP3DvalAE_2DTE_1DAE_trk0/D");
  tree->Branch("pftauchhads_IP3DvalAE_2DTE_1DTE_trk0", &pftauchhads_IP3DvalAE_2DTE_1DTE_trk0, "pftauchhads_IP3DvalAE_2DTE_1DTE_trk0/D");
  tree->Branch("pftauchhads_IP3DvalAE_2DAE_1DAE_trk0", &pftauchhads_IP3DvalAE_2DAE_1DAE_trk0, "pftauchhads_IP3DvalAE_2DAE_1DAE_trk0/D");
  tree->Branch("pftauchhads_IP3DvalTE_2DTE_1DTE_trk0", &pftauchhads_IP3DvalTE_2DTE_1DTE_trk0, "pftauchhads_IP3DvalTE_2DTE_1DTE_trk0/D");
  //Collimation
  tree->Branch("pftauchhads_JTD_val_trk0", &pftauchhads_JTD_val_trk0, "pftauchhads_JTD_val_trk0/D");
  tree->Branch("pftauchhads_JTD_val_trk1", &pftauchhads_JTD_val_trk1, "pftauchhads_JTD_val_trk1/D");
  tree->Branch("pftauchhads_JTD_val_trk2", &pftauchhads_JTD_val_trk2, "pftauchhads_JTD_val_trk2/D");
  tree->Branch("pftauchhads_JTD_sig_trk0", &pftauchhads_JTD_sig_trk0, "pftauchhads_JTD_sig_trk0/D");
  tree->Branch("pftauchhads_JTD_sig_trk1", &pftauchhads_JTD_sig_trk1, "pftauchhads_JTD_sig_trk1/D");
  tree->Branch("pftauchhads_JTD_sig_trk2", &pftauchhads_JTD_sig_trk2, "pftauchhads_JTD_sig_trk2/D");
  tree->Branch("pftau_pvsv_dist3d_val", &pftau_pvsv_dist3d_val, "pftau_pvsv_dist3d_val/D");
  tree->Branch("pftau_pvsv_dist2d_val", &pftau_pvsv_dist2d_val, "pftau_pvsv_dist2d_val/D");
  tree->Branch("pftau_pvsv_dist3d_sig", &pftau_pvsv_dist3d_sig, "pftau_pvsv_dist3d_sig/D");
  tree->Branch("pftau_pvsv_dist2d_sig", &pftau_pvsv_dist2d_sig, "pftau_pvsv_dist2d_sig/D");
  tree->Branch("pv_avfbs_nchi2", &pv_avfbs_nchi2, "pv_avfbs_nchi2/D");
  tree->Branch("unbpv_avfbs_nchi2", &unbpv_avfbs_nchi2, "unbpv_avfbs_nchi2/D");
  tree->Branch("diff_pv_avfbs_nchi2_unbpv_avfbs_nchi2", &diff_pv_avfbs_nchi2_unbpv_avfbs_nchi2, "diff_pv_avfbs_nchi2_unbpv_avfbs_nchi2/D");
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
  tree->Branch("pftauchhads_IP3Ddv_val", &pftauchhads_IP3Ddv_val, "pftauchhads_IP3Ddv_val[pftauchhads_numv]/D");
  tree->Branch("pftauchhads_IP3Ddv_err", &pftauchhads_IP3Ddv_err, "pftauchhads_IP3Ddv_err[pftauchhads_numv]/D");
  tree->Branch("pftauchhads_IP3Ddv_sig", &pftauchhads_IP3Ddv_sig, "pftauchhads_IP3Ddv_sig[pftauchhads_numv]/D");
  tree->Branch("pftauchhads_IP2Ddv_val", &pftauchhads_IP2Ddv_val, "pftauchhads_IP2Ddv_val[pftauchhads_numv]/D");
  tree->Branch("pftauchhads_IP2Ddv_err", &pftauchhads_IP2Ddv_err, "pftauchhads_IP2Ddv_err[pftauchhads_numv]/D");
  tree->Branch("pftauchhads_IP2Ddv_sig", &pftauchhads_IP2Ddv_sig, "pftauchhads_IP2Ddv_sig[pftauchhads_numv]/D");
  tree->Branch("pftauchhads_sIP3Ddv_val", &pftauchhads_sIP3Ddv_val, "pftauchhads_sIP3Ddv_val[pftauchhads_numv]/D");
  tree->Branch("pftauchhads_sIP3Ddv_err", &pftauchhads_sIP3Ddv_err, "pftauchhads_sIP3Ddv_err[pftauchhads_numv]/D");
  tree->Branch("pftauchhads_sIP3Ddv_sig", &pftauchhads_sIP3Ddv_sig, "pftauchhads_sIP3Ddv_sig[pftauchhads_numv]/D");
  tree->Branch("pftauchhads_sIP2Ddv_val", &pftauchhads_sIP2Ddv_val, "pftauchhads_sIP2Ddv_val[pftauchhads_numv]/D");
  tree->Branch("pftauchhads_sIP2Ddv_err", &pftauchhads_sIP2Ddv_err, "pftauchhads_sIP2Ddv_err[pftauchhads_numv]/D");
  tree->Branch("pftauchhads_sIP2Ddv_sig", &pftauchhads_sIP2Ddv_sig, "pftauchhads_sIP2Ddv_sig[pftauchhads_numv]/D");
  tree->Branch("pftauchhads_IP1D_val", &pftauchhads_IP1D_val, "pftauchhads_IP1D_val[pftauchhads_numv]/D");
  tree->Branch("pftauchhads_IP1D_err", &pftauchhads_IP1D_err, "pftauchhads_IP1D_err[pftauchhads_numv]/D");
  tree->Branch("pftauchhads_IP1D_sig", &pftauchhads_IP1D_sig, "pftauchhads_IP1D_sig[pftauchhads_numv]/D");
  tree->Branch("pftauchhads_sIP1D_val", &pftauchhads_sIP1D_val, "pftauchhads_sIP1D_val[pftauchhads_numv]/D");
  tree->Branch("pftauchhads_sIP1D_err", &pftauchhads_sIP1D_err, "pftauchhads_sIP1D_err[pftauchhads_numv]/D");
  tree->Branch("pftauchhads_sIP1D_sig", &pftauchhads_sIP1D_sig, "pftauchhads_sIP1D_sig[pftauchhads_numv]/D");
  tree->Branch("dR_tautrk_recotau", &dR_tautrk_recotau, "dR_tautrk_recotau[pftauchhads_numv]/D");
 }
 /////
 //   Set branch address
 /////
 //Connects the branches of an existing TTree to variables used when loading the file
 void set_branch_addresses(void){
  //Usefull
  tree->SetBranchAddress("duno", &duno);
  tree->SetBranchAddress("numrecovtcs", &numrecovtcs);
  //Vertex resolution
  tree->SetBranchAddress("genditauvtx_x", &genditauvtx_x);
  tree->SetBranchAddress("genditauvtx_y", &genditauvtx_y);
  tree->SetBranchAddress("genditauvtx_z", &genditauvtx_z);
  tree->SetBranchAddress("unbpv_KVF_nv", &unbpv_KVF_nv);
  tree->SetBranchAddress("unbpv_KVFbs_nv", &unbpv_KVFbs_nv);
  tree->SetBranchAddress("unbpv_AVF_nv", &unbpv_AVF_nv);
  tree->SetBranchAddress("unbpv_AVFbs_nv", &unbpv_AVFbs_nv);
  tree->SetBranchAddress("vtxKVF_x", &vtxKVF_x);
  tree->SetBranchAddress("vtxKVF_y", &vtxKVF_y);
  tree->SetBranchAddress("vtxKVF_z", &vtxKVF_z);
  tree->SetBranchAddress("vtxKVFbs_x", &vtxKVFbs_x);
  tree->SetBranchAddress("vtxKVFbs_y", &vtxKVFbs_y);
  tree->SetBranchAddress("vtxKVFbs_z", &vtxKVFbs_z);
  tree->SetBranchAddress("vtxAVF_x", &vtxAVF_x);
  tree->SetBranchAddress("vtxAVF_y", &vtxAVF_y);
  tree->SetBranchAddress("vtxAVF_z", &vtxAVF_z);
  tree->SetBranchAddress("vtxAVFbs_x", &vtxAVFbs_x);
  tree->SetBranchAddress("vtxAVFbs_y", &vtxAVFbs_y);
  tree->SetBranchAddress("vtxAVFbs_z", &vtxAVFbs_z);
  tree->SetBranchAddress("pvtauvtx_genditauvtx_x", &pvtauvtx_genditauvtx_x);
  tree->SetBranchAddress("pvtauvtx_genditauvtx_y", &pvtauvtx_genditauvtx_y);
  tree->SetBranchAddress("pvtauvtx_genditauvtx_z", &pvtauvtx_genditauvtx_z);
  tree->SetBranchAddress("unbpv_KVF_genditauvtx_x", &unbpv_KVF_genditauvtx_x);
  tree->SetBranchAddress("unbpv_KVF_genditauvtx_y", &unbpv_KVF_genditauvtx_y);
  tree->SetBranchAddress("unbpv_KVF_genditauvtx_z", &unbpv_KVF_genditauvtx_z);
  //RecoJetTau
  tree->SetBranchAddress("recojettau_pt", &recojettau_pt);
  tree->SetBranchAddress("recojettau_eta", &recojettau_eta);
  tree->SetBranchAddress("recojettau_phi", &recojettau_phi);
  tree->SetBranchAddress("recojettau_en", &recojettau_en);
  tree->SetBranchAddress("tau_pt_DIV_recojettau_pt", &tau_pt_DIV_recojettau_pt);
  tree->SetBranchAddress("tau_en_DIV_recojettau_en", &tau_en_DIV_recojettau_en);
  tree->SetBranchAddress("tau_trk0_pt_DIV_recojettau_pt", &tau_trk0_pt_DIV_recojettau_pt);
  tree->SetBranchAddress("tau_trk1_pt_DIV_recojettau_pt", &tau_trk1_pt_DIV_recojettau_pt);
  tree->SetBranchAddress("tau_trk2_pt_DIV_recojettau_pt", &tau_trk2_pt_DIV_recojettau_pt);
  tree->SetBranchAddress("tau_trk0_en_DIV_recojettau_en", &tau_trk0_en_DIV_recojettau_en);
  tree->SetBranchAddress("tau_trk1_en_DIV_recojettau_en", &tau_trk1_en_DIV_recojettau_en);
  tree->SetBranchAddress("tau_trk2_en_DIV_recojettau_en", &tau_trk2_en_DIV_recojettau_en);
  //dR
  tree->SetBranchAddress("dR_recojettau_recotau", &dR_recojettau_recotau);
  tree->SetBranchAddress("dR_tautrk_recojettau_trk0", &dR_tautrk_recojettau_trk0);
  tree->SetBranchAddress("dR_tautrk_recojettau_trk1", &dR_tautrk_recojettau_trk1);
  tree->SetBranchAddress("dR_tautrk_recojettau_trk2", &dR_tautrk_recojettau_trk2);
  tree->SetBranchAddress("dR_RecoGen", &dR_RecoGen);
  tree->SetBranchAddress("dR_vtxDirGen", &dR_vtxDirGen);
  tree->SetBranchAddress("dR_tautrk_recotau_trk0", &dR_tautrk_recotau_trk0);
  tree->SetBranchAddress("dR_tautrk_recotau_trk1", &dR_tautrk_recotau_trk1);
  tree->SetBranchAddress("dR_tautrk_recotau_trk2", &dR_tautrk_recotau_trk2);
  tree->SetBranchAddress("dR_TauFlightDist_JetAxis", &dR_TauFlightDist_JetAxis);
  //pftaugv
  tree->SetBranchAddress("pftaugv_x", &pftaugv_x);
  tree->SetBranchAddress("pftaugv_y", &pftaugv_y);
  tree->SetBranchAddress("pftaugv_z", &pftaugv_z);
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
  tree->SetBranchAddress("pftauchhads_IP3D_val_trk0", &pftauchhads_IP3D_val_trk0);
  tree->SetBranchAddress("pftauchhads_IP3D_val_trk1", &pftauchhads_IP3D_val_trk1);
  tree->SetBranchAddress("pftauchhads_IP3D_val_trk2", &pftauchhads_IP3D_val_trk2);
  tree->SetBranchAddress("pftauchhads_IP3D_sig_trk0", &pftauchhads_IP3D_sig_trk0);
  tree->SetBranchAddress("pftauchhads_IP3D_sig_trk1", &pftauchhads_IP3D_sig_trk1);
  tree->SetBranchAddress("pftauchhads_IP3D_sig_trk2", &pftauchhads_IP3D_sig_trk2);
  tree->SetBranchAddress("pftauchhads_IP2D_val_trk0", &pftauchhads_IP2D_val_trk0);
  tree->SetBranchAddress("pftauchhads_IP2D_val_trk1", &pftauchhads_IP2D_val_trk1);
  tree->SetBranchAddress("pftauchhads_IP2D_val_trk2", &pftauchhads_IP2D_val_trk2);
  tree->SetBranchAddress("pftauchhads_IP2D_sig_trk0", &pftauchhads_IP2D_sig_trk0);
  tree->SetBranchAddress("pftauchhads_IP2D_sig_trk1", &pftauchhads_IP2D_sig_trk1);
  tree->SetBranchAddress("pftauchhads_IP2D_sig_trk2", &pftauchhads_IP2D_sig_trk2);
  tree->SetBranchAddress("pftauchhads_IP1D_val_trk0", &pftauchhads_IP1D_val_trk0);
  tree->SetBranchAddress("pftauchhads_IP1D_val_trk1", &pftauchhads_IP1D_val_trk1);
  tree->SetBranchAddress("pftauchhads_IP1D_val_trk2", &pftauchhads_IP1D_val_trk2);
  tree->SetBranchAddress("pftauchhads_IP1D_sig_trk0", &pftauchhads_IP1D_sig_trk0);
  tree->SetBranchAddress("pftauchhads_IP1D_sig_trk1", &pftauchhads_IP1D_sig_trk1);
  tree->SetBranchAddress("pftauchhads_IP1D_sig_trk2", &pftauchhads_IP1D_sig_trk2);
  tree->SetBranchAddress("pftauchhads_AEIP1D_val_trk0", &pftauchhads_AEIP1D_val_trk0);
  tree->SetBranchAddress("pftauchhads_AEIP1D_val_trk1", &pftauchhads_AEIP1D_val_trk1);
  tree->SetBranchAddress("pftauchhads_AEIP1D_val_trk2", &pftauchhads_AEIP1D_val_trk2);
  tree->SetBranchAddress("pftauchhads_AEIP1D_sig_trk0", &pftauchhads_AEIP1D_sig_trk0);
  tree->SetBranchAddress("pftauchhads_AEIP1D_sig_trk1", &pftauchhads_AEIP1D_sig_trk1);
  tree->SetBranchAddress("pftauchhads_sIP3D_val_trk0", &pftauchhads_sIP3D_val_trk0);
  tree->SetBranchAddress("pftauchhads_sIP3D_val_trk1", &pftauchhads_sIP3D_val_trk1);
  tree->SetBranchAddress("pftauchhads_sIP3D_val_trk2", &pftauchhads_sIP3D_val_trk2);
  tree->SetBranchAddress("pftauchhads_sIP3D_sig_trk0", &pftauchhads_sIP3D_sig_trk0);
  tree->SetBranchAddress("pftauchhads_sIP3D_sig_trk1", &pftauchhads_sIP3D_sig_trk1);
  tree->SetBranchAddress("pftauchhads_sIP3D_sig_trk2", &pftauchhads_sIP3D_sig_trk2);
  tree->SetBranchAddress("pftauchhads_sIP2D_val_trk0", &pftauchhads_sIP2D_val_trk0);
  tree->SetBranchAddress("pftauchhads_sIP2D_val_trk1", &pftauchhads_sIP2D_val_trk1);
  tree->SetBranchAddress("pftauchhads_sIP2D_val_trk2", &pftauchhads_sIP2D_val_trk2);
  tree->SetBranchAddress("pftauchhads_sIP2D_sig_trk0", &pftauchhads_sIP2D_sig_trk0);
  tree->SetBranchAddress("pftauchhads_sIP2D_sig_trk1", &pftauchhads_sIP2D_sig_trk1);
  tree->SetBranchAddress("pftauchhads_sIP2D_sig_trk2", &pftauchhads_sIP2D_sig_trk2);
  tree->SetBranchAddress("pftauchhads_sIP1D_val_trk0", &pftauchhads_sIP1D_val_trk0);
  tree->SetBranchAddress("pftauchhads_sIP1D_val_trk1", &pftauchhads_sIP1D_val_trk1);
  tree->SetBranchAddress("pftauchhads_sIP1D_val_trk2", &pftauchhads_sIP1D_val_trk2);
  tree->SetBranchAddress("pftauchhads_sIP1D_sig_trk0", &pftauchhads_sIP1D_sig_trk0);
  tree->SetBranchAddress("pftauchhads_sIP1D_sig_trk1", &pftauchhads_sIP1D_sig_trk1);
  tree->SetBranchAddress("pftauchhads_sIP1D_sig_trk2", &pftauchhads_sIP1D_sig_trk2);
  tree->SetBranchAddress("pftauchhads_AEIP1D_sig_trk2", &pftauchhads_AEIP1D_sig_trk2);
  tree->SetBranchAddress("pftauchhads_AEsIP1D_val_trk0", &pftauchhads_AEsIP1D_val_trk0);
  tree->SetBranchAddress("pftauchhads_AEsIP1D_val_trk1", &pftauchhads_AEsIP1D_val_trk1);
  tree->SetBranchAddress("pftauchhads_AEsIP1D_val_trk2", &pftauchhads_AEsIP1D_val_trk2);
  tree->SetBranchAddress("pftauchhads_AEsIP1D_sig_trk0", &pftauchhads_AEsIP1D_sig_trk0);
  tree->SetBranchAddress("pftauchhads_AEsIP1D_sig_trk1", &pftauchhads_AEsIP1D_sig_trk1);
  tree->SetBranchAddress("pftauchhads_AEsIP1D_sig_trk2", &pftauchhads_AEsIP1D_sig_trk2);
  tree->SetBranchAddress("pftauchhads_sDL3D_val_trk0", &pftauchhads_sDL3D_val_trk0);
  tree->SetBranchAddress("pftauchhads_sDL3D_val_trk1", &pftauchhads_sDL3D_val_trk1);
  tree->SetBranchAddress("pftauchhads_sDL3D_val_trk2", &pftauchhads_sDL3D_val_trk2);
  tree->SetBranchAddress("pftauchhads_sDL3D_sig_trk0", &pftauchhads_sDL3D_sig_trk0);
  tree->SetBranchAddress("pftauchhads_sDL3D_sig_trk1", &pftauchhads_sDL3D_sig_trk1);
  tree->SetBranchAddress("pftauchhads_sDL3D_sig_trk2", &pftauchhads_sDL3D_sig_trk2);
  tree->SetBranchAddress("pftauchhads_absDL3D_val_trk0", &pftauchhads_absDL3D_val_trk0);
  tree->SetBranchAddress("pftauchhads_absDL3D_val_trk1", &pftauchhads_absDL3D_val_trk1);
  tree->SetBranchAddress("pftauchhads_absDL3D_val_trk2", &pftauchhads_absDL3D_val_trk2);
  tree->SetBranchAddress("pftauchhads_absDL3D_sig_trk0", &pftauchhads_absDL3D_sig_trk0);
  tree->SetBranchAddress("pftauchhads_absDL3D_sig_trk1", &pftauchhads_absDL3D_sig_trk1);
  tree->SetBranchAddress("pftauchhads_absDL3D_sig_trk2", &pftauchhads_absDL3D_sig_trk2);
  //Compare Analytical Extrapolator (AE) vs Transverse Extrapolator (TE)
  tree->SetBranchAddress("pftauchhads_AEIP1D_x_val_trk0", &pftauchhads_AEIP1D_x_val_trk0);
  tree->SetBranchAddress("pftauchhads_AEIP1D_y_val_trk0", &pftauchhads_AEIP1D_y_val_trk0);
  tree->SetBranchAddress("pftauchhads_IP3DvalAE_trk0", &pftauchhads_IP3DvalAE_trk0);
  tree->SetBranchAddress("pftauchhads_IP3DvalTE_trk0", &pftauchhads_IP3DvalTE_trk0);
  tree->SetBranchAddress("pftauchhads_IP2DvalAE_trk0", &pftauchhads_IP2DvalAE_trk0);
  tree->SetBranchAddress("pftauchhads_IP2DvalTE_trk0", &pftauchhads_IP2DvalTE_trk0);
  tree->SetBranchAddress("pftauchhads_IP1DvalAE_trk0", &pftauchhads_IP1DvalAE_trk0);
  tree->SetBranchAddress("pftauchhads_IP1DvalTE_trk0", &pftauchhads_IP1DvalTE_trk0);
  tree->SetBranchAddress("pftauchhads_IP3DvalAE_2DTE_1DAE_trk0", &pftauchhads_IP3DvalAE_2DTE_1DAE_trk0);
  tree->SetBranchAddress("pftauchhads_IP3DvalAE_2DTE_1DTE_trk0", &pftauchhads_IP3DvalAE_2DTE_1DTE_trk0);
  tree->SetBranchAddress("pftauchhads_IP3DvalAE_2DAE_1DAE_trk0", &pftauchhads_IP3DvalAE_2DAE_1DAE_trk0);
  tree->SetBranchAddress("pftauchhads_IP3DvalTE_2DTE_1DTE_trk0", &pftauchhads_IP3DvalTE_2DTE_1DTE_trk0);
  //Collimation
  tree->SetBranchAddress("pftauchhads_JTD_val_trk0", &pftauchhads_JTD_val_trk0);
  tree->SetBranchAddress("pftauchhads_JTD_val_trk1", &pftauchhads_JTD_val_trk1);
  tree->SetBranchAddress("pftauchhads_JTD_val_trk2", &pftauchhads_JTD_val_trk2);
  tree->SetBranchAddress("pftauchhads_JTD_sig_trk0", &pftauchhads_JTD_sig_trk0);
  tree->SetBranchAddress("pftauchhads_JTD_sig_trk1", &pftauchhads_JTD_sig_trk1);
  tree->SetBranchAddress("pftauchhads_JTD_sig_trk2", &pftauchhads_JTD_sig_trk2);
  tree->SetBranchAddress("pftau_pvsv_dist3d_val", &pftau_pvsv_dist3d_val);
  tree->SetBranchAddress("pftau_pvsv_dist2d_val", &pftau_pvsv_dist2d_val);
  tree->SetBranchAddress("pftau_pvsv_dist3d_sig", &pftau_pvsv_dist3d_sig);
  tree->SetBranchAddress("pftau_pvsv_dist2d_sig", &pftau_pvsv_dist2d_sig);
  tree->SetBranchAddress("pv_avfbs_nchi2", &pv_avfbs_nchi2);
  tree->SetBranchAddress("unbpv_avfbs_nchi2", &unbpv_avfbs_nchi2);
  tree->SetBranchAddress("diff_pv_avfbs_nchi2_unbpv_avfbs_nchi2", &diff_pv_avfbs_nchi2_unbpv_avfbs_nchi2);
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
  tree->SetBranchAddress("pftauchhads_IP3Ddv_val", &pftauchhads_IP3Ddv_val);
  tree->SetBranchAddress("pftauchhads_IP3Ddv_err", &pftauchhads_IP3Ddv_err);
  tree->SetBranchAddress("pftauchhads_IP3Ddv_sig", &pftauchhads_IP3Ddv_sig);
  tree->SetBranchAddress("pftauchhads_IP2Ddv_val", &pftauchhads_IP2Ddv_val);
  tree->SetBranchAddress("pftauchhads_IP2Ddv_err", &pftauchhads_IP2Ddv_err);
  tree->SetBranchAddress("pftauchhads_IP2Ddv_sig", &pftauchhads_IP2Ddv_sig);
  tree->SetBranchAddress("pftauchhads_sIP3Ddv_val", &pftauchhads_sIP3Ddv_val);
  tree->SetBranchAddress("pftauchhads_sIP3Ddv_err", &pftauchhads_sIP3Ddv_err);
  tree->SetBranchAddress("pftauchhads_sIP3Ddv_sig", &pftauchhads_sIP3Ddv_sig);
  tree->SetBranchAddress("pftauchhads_sIP2Ddv_val", &pftauchhads_sIP2Ddv_val);
  tree->SetBranchAddress("pftauchhads_sIP2Ddv_err", &pftauchhads_sIP2Ddv_err);
  tree->SetBranchAddress("pftauchhads_sIP2Ddv_sig", &pftauchhads_sIP2Ddv_sig);
  tree->SetBranchAddress("pftauchhads_IP1D_val", &pftauchhads_IP1D_val);
  tree->SetBranchAddress("pftauchhads_IP1D_err", &pftauchhads_IP1D_err);
  tree->SetBranchAddress("pftauchhads_IP1D_sig", &pftauchhads_IP1D_sig);
  tree->SetBranchAddress("pftauchhads_sIP1D_val", &pftauchhads_sIP1D_val);
  tree->SetBranchAddress("pftauchhads_sIP1D_err", &pftauchhads_sIP1D_err);
  tree->SetBranchAddress("pftauchhads_sIP1D_sig", &pftauchhads_sIP1D_sig);
  tree->SetBranchAddress("dR_tautrk_recotau", &dR_tautrk_recotau);
 }
};
#endif
