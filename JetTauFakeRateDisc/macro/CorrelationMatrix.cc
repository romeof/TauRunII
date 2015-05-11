/**
This Macro   
1. Studies correlation among some variables
    Plots scatter plots and correlation matrix  

Need to specify
1. See Declare Constants
*/
/////
//   To run: root -l CorrelationMatrix.cc+
/////
/////
//   Prepare Root and Roofit
/////
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF2.h"
#include "THStack.h"
#include "TLegend.h"
//#include "TPaletteAxis.h"
#include "TGaxis.h"
#include "TLine.h"
#include "TColor.h"
#include "TTree.h"
#include "TTreePlayer.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TEfficiency.h"
#include "TEllipse.h"
#include "TArc.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <iostream>
using namespace std;
/////
//   Declare constants
/////
//Path - channel - samples - selection  
const string path       = "../Rootplas/";
const char *samples[]   = {"bkg", "GGHLL_b"};
const string selection  = "";
const double Luminosity = 19703.225;//19779.362;//19600;; //pb^-1
const bool noLumiNorm   = true; //true means NO luminosity normalization done
const bool noPUcorr     = true; //true means NO PU corr done done
const bool save_plots   = true;
const int  numVar = 100;

const int  vini   = 0;
const int  vfin   = 7;
const int  vvini  = 0;
const int  vvfin  = 7;

const string corrmat = "CorrelationMat_";
const char *variablesNum[]    = {"pftauchhads_IP2D_val", "pftauchhads_IP2D_sig", "pftauchhads_IP3D_val", "pftauchhads_IP3D_sig", "pftauchhads_sIP2D_val", "pftauchhads_sIP2D_sig", "pftauchhads_sIP3D_val", "pftauchhads_sIP3D_sig"};//, "pftauchhads_sIP2D_val", "pftauchhads_sIP2D_sig", "pftauchhads_sIP3D_val", "pftauchhads_IP3D_sig"};
//Meaning of vtNum, vtDen
//int (i): 0, double (d): 1, int arr (iarr): arr pos = abs(- val) -2, double arr (darr): arr pos = val +2
const int  vtNum[numVar]      = {2, 2, 2, 2, 2, 2, 2, 2};
const char *variablesDen[]    = {"duno", "duno", "duno", "duno", "duno", "duno", "duno", "duno"};
const int  vtDen[numVar]      = {1, 1, 1, 1, 1, 1, 1, 1};
const char *titles[]          = {"IP2D_val", "IP2D_sig", "IP3D_val", "IP3D_sig", "sIP2D_val", "sIP2D_sig", "sIP3D_val", "sIP3D_sig"};
const char *titlespdf[]       = {"IP2D_val", "IP2D_sig", "IP3D_val", "IP3D_sig", "sIP2D_val", "sIP2D_sig", "sIP3D_val", "sIP3D_sig"};
const double inRange[numVar]  = {0, 0, 0, 0, -1.5, -100, -30, -200};
const double endRange[numVar] = {3.5, 100, 35, 200, 1.5, 100, 30, 200};
const int    bin[numVar]      = {100, 100, 100, 100, 100, 100, 100, 100};
/////
//   Declare functions 
/////
TFile* Call_TFile(string rootpla);
//XNum_XDen_YNum_YDen
//iarr_d
TH2F*  iarr_d_darr_darr(string sample_vars, TTree* tree, int r,int v,int vv);
TH2F*  iarr_d_darr_iarr(string sample_vars, TTree* tree, int r,int v,int vv);
TH2F*  iarr_d_iarr_d(string sample_vars, TTree* tree, int r,int v,int vv);
//darr_d
TH2F*  darr_d_iarr_d(string sample_vars, TTree* tree, int r,int v,int vv);
TH2F*  darr_d_darr_d(string sample_vars, TTree* tree, int r,int v,int vv);
TH2F*  darr_d_darr_iarr(string sample_vars, TTree* tree, int r,int v,int vv);
TH2F*  darr_d_darr_darr(string sample_vars, TTree* tree, int r,int v,int vv);
//darr_darr
TH2F*  darr_darr_darr_darr(string sample_vars, TTree* tree, int r,int v,int vv);
TH2F*  darr_darr_darr_iarr(string sample_vars, TTree* tree, int r,int v,int vv);
//darr_iarr
TH2F*  darr_iarr_darr_iarr(string sample_vars, TTree* tree, int r,int v,int vv);
TH2F*  darr_iarr_darr_darr(string sample_vars, TTree* tree, int r,int v,int vv);

void   drawscattplot(TH2F* histo);
void   drawcorrmat(TH2F* histo);
TH2F*  emptyhisto(string sample_vars,int r,int v,int vv);
void setTDRStyle();
/////
//   Main function
/////
void CorrelationMatrix(){
 setTDRStyle();
 vector<string> varNum(variablesNum, variablesNum + sizeof(variablesNum)/sizeof(variablesNum[0]));
 vector<string> varDen(variablesDen, variablesDen + sizeof(variablesDen)/sizeof(variablesDen[0]));
 vector<string> varTitles(titles, titles + sizeof(titles)/sizeof(titles[0]));
 vector<string> varTitlespdf(titlespdf, titlespdf + sizeof(titlespdf)/sizeof(titlespdf[0]));
 //Loop over samples
 vector<string> rootplas(samples, samples + sizeof(samples)/sizeof(samples[0]));
 for(uint r=0; r<rootplas.size(); r++){
  cout<<rootplas[r]<<endl;
  string s_corrmat = "Correlation matrix";
  if(rootplas[r]=="bkg_tot") s_corrmat = s_corrmat+" Jet";
  if(rootplas[r]=="GGHLL")    s_corrmat = s_corrmat+" Tau";
  TH2F* h_corrmat = new TH2F(s_corrmat.c_str(), s_corrmat.c_str(), (vfin-vini)+1,0,(vfin-vini)+1, (vvfin-vvini)+1,0,(vvfin-vvini)+1);
  //Loop on pairs of variables
  for(int v=vini,posv=1; v<=vfin; v++,posv++){
   for(int vv=vvini,posvv=1; vv<=vvfin; vv++,posvv++){  
    if(vv==v){
     h_corrmat->SetBinContent(posv,posvv,1.0);
     h_corrmat->GetXaxis()->SetBinLabel(posv,varTitles[v].c_str());
     h_corrmat->GetYaxis()->SetBinLabel(posvv,varTitles[vv].c_str());
    } 
    if(vv<=v) continue;
    //Call tree  
    TFile* f = Call_TFile(rootplas[r]); TTree* tree; f->GetObject("tree",tree);
    cout<<setiosflags(ios::fixed)<<setprecision(10);
    //Canvas
    const string sample_vars = "Correlation_"+rootplas[r]+"_"+varTitles[v]+"_"+varTitles[vv];
    TCanvas* c1 = new TCanvas(sample_vars.c_str(),sample_vars.c_str(),100,100,700,700);
    //Tak correlation plot for the histo
    TH2F *h_scattplot;
    //iarr_d
    if(vtNum[v]<-1 && vtDen[v]==1 && vtNum[vv]>1 && vtDen[vv]>1)   h_scattplot = iarr_d_darr_darr(sample_vars,tree,r,v,vv);
    if(vtNum[v]<-1 && vtDen[v]==1 && vtNum[vv]>1 && vtDen[vv]<-1)  h_scattplot = iarr_d_darr_iarr(sample_vars,tree,r,v,vv);
    if(vtNum[v]<-1 && vtDen[v]==1 && vtNum[vv]<-1 && vtDen[vv]==1) h_scattplot = iarr_d_iarr_d(sample_vars,tree,r,v,vv);
    //darr_d
    if(vtNum[v]>1 && vtDen[v]==1 && vtNum[vv]<-1 && vtDen[vv]==1) h_scattplot = darr_d_iarr_d(sample_vars,tree,r,v,vv);
    if(vtNum[v]>1 && vtDen[v]==1 && vtNum[vv]>1 && vtDen[vv]==1)  h_scattplot = darr_d_darr_d(sample_vars,tree,r,v,vv);
    if(vtNum[v]>1 && vtDen[v]==1 && vtNum[vv]>1 && vtDen[vv]<-1)  h_scattplot = darr_d_darr_iarr(sample_vars,tree,r,v,vv);
    if(vtNum[v]>1 && vtDen[v]==1 && vtNum[vv]>1 && vtDen[vv]>1)   h_scattplot = darr_d_darr_darr(sample_vars,tree,r,v,vv);
    //darr_darr
    if(vtNum[v]>1 && vtDen[v]>1 && vtNum[vv]>1 && vtDen[vv]>1)    h_scattplot = darr_darr_darr_darr(sample_vars,tree,r,v,vv);
    if(vtNum[v]>1 && vtDen[v]>1  && vtNum[vv]>1 && vtDen[vv]<-1)  h_scattplot = darr_darr_darr_iarr(sample_vars,tree,r,v,vv);
    //darr_iarr
    if(vtNum[v]>1 && vtDen[v]<-1 && vtNum[vv]>1 && vtDen[vv]<-1)  h_scattplot = darr_iarr_darr_iarr(sample_vars,tree,r,v,vv);
    if(vtNum[v]>1 && vtDen[v]<-1 && vtNum[vv]>1 && vtDen[vv]>1)   h_scattplot = darr_iarr_darr_darr(sample_vars,tree,r,v,vv);
    drawscattplot(h_scattplot);
    //Take correlation values
    double d_corrval = h_scattplot->GetCorrelationFactor(); 
    cout<<varTitles[v]<<setw(50)<<varTitles[vv]<<setw(50)<<d_corrval<<endl;
    h_corrmat->SetBinContent(posv,posvv,d_corrval); 
    h_corrmat->GetXaxis()->SetBinLabel(posv,varTitles[v].c_str());
    h_corrmat->GetYaxis()->SetBinLabel(posvv,varTitles[vv].c_str());
    const string sample_vars_pdf = "Correlation_"+rootplas[r]+"___"+varTitlespdf[v]+"___"+varTitlespdf[vv];
    if(save_plots) c1->SaveAs((sample_vars_pdf+".pdf").c_str());
    delete tree;
   } 
  }//Combinatoric of variables 
  string s_corr    = corrmat+rootplas[r];
  TCanvas* c0 = new TCanvas(s_corr.c_str(),s_corr.c_str(),100,100,700,600);
  c0->SaveAs((s_corr+".pdf").c_str());
  drawcorrmat(h_corrmat);
 }//Loop on samples 
}
/////
//   Call TFile to be read
/////
TFile* Call_TFile(string rootpla){
 string dotroot   = ".root";
 string file_name = path+rootpla+selection+dotroot;
 TFile* f = new TFile(file_name.c_str(),"update");
 return f;
}
/////
//   Return the scatter plot of variables in position v and vv
/////
//iarr_d
TH2F* iarr_d_darr_darr(string sample_vars, TTree* tree, int r,int v,int vv){
 TH2F *h_scattplot = emptyhisto(sample_vars,r,v,vv);
 vector<string> varNum(variablesNum, variablesNum + sizeof(variablesNum)/sizeof(variablesNum[0]));
 vector<string> varDen(variablesDen, variablesDen + sizeof(variablesDen)/sizeof(variablesDen[0]));
 int    varXNum[100];
 double varXDen;
 double varYNum[100];
 double varYDen[100];
 tree->SetBranchAddress(varNum[v].c_str(),&varXNum);
 tree->SetBranchAddress(varDen[v].c_str(),&varXDen);
 tree->SetBranchAddress(varNum[vv].c_str(),&varYNum);
 tree->SetBranchAddress(varDen[vv].c_str(),&varYDen);
 double wgt_lumi = 0; double wgt_pu = 0;
 //Fill histo
 double weight = 0.;
 for(int en=0; en<tree->GetEntries(); en++)
 //for(int en=0; en<10; en++)
 {
  tree->GetEntry(en);
  weight = wgt_lumi;
  weight = weight*Luminosity/100.;
  if(noLumiNorm) weight = 1.;
  if(noPUcorr)   wgt_pu = 1.;
  h_scattplot->Fill(varXNum[abs(vtNum[v])-2]/double(varXDen),varYNum[abs(vtNum[vv])-2]/double(varYDen[abs(vtDen[vv])-2]),weight*wgt_pu);
 }
 return h_scattplot;
}
TH2F*  iarr_d_darr_iarr(string sample_vars, TTree* tree, int r,int v,int vv){
 TH2F *h_scattplot = emptyhisto(sample_vars,r,v,vv);
 vector<string> varNum(variablesNum, variablesNum + sizeof(variablesNum)/sizeof(variablesNum[0]));
 vector<string> varDen(variablesDen, variablesDen + sizeof(variablesDen)/sizeof(variablesDen[0]));
 int    varXNum[100];  ///int///
 double varXDen;
 double varYNum[100];
 int    varYDen[100];  ///int///
 tree->SetBranchAddress(varNum[v].c_str(),&varXNum);
 tree->SetBranchAddress(varDen[v].c_str(),&varXDen);
 tree->SetBranchAddress(varNum[vv].c_str(),&varYNum);
 tree->SetBranchAddress(varDen[vv].c_str(),&varYDen);
 double wgt_lumi = 0; double wgt_pu = 0;
 //Fill histo
 double weight = 0.;
 for(int en=0; en<tree->GetEntries(); en++)
 {
  tree->GetEntry(en);
  weight = wgt_lumi;
  weight = weight*Luminosity/100.;
  if(noLumiNorm) weight = 1.;
  if(noPUcorr)   wgt_pu = 1.;
  h_scattplot->Fill(varXNum[abs(vtNum[v])-2]/double(varXDen),varYNum[abs(vtNum[vv])-2]/double(varYDen[abs(vtDen[vv])-2]),weight*wgt_pu);
 }
 return h_scattplot;
}
TH2F*  iarr_d_iarr_d(string sample_vars, TTree* tree, int r,int v,int vv){
 TH2F *h_scattplot = emptyhisto(sample_vars,r,v,vv);
 vector<string> varNum(variablesNum, variablesNum + sizeof(variablesNum)/sizeof(variablesNum[0]));
 vector<string> varDen(variablesDen, variablesDen + sizeof(variablesDen)/sizeof(variablesDen[0]));
 int    varXNum[100];   ///int///
 double varXDen;
 tree->SetBranchAddress(varNum[v].c_str(),&varXNum);
 tree->SetBranchAddress(varDen[v].c_str(),&varXDen);
 double wgt_lumi = 0; double wgt_pu = 0;
 //Fill histo
 double weight = 0.;
 for(int en=0; en<tree->GetEntries(); en++)
 {
  tree->GetEntry(en);
  weight = wgt_lumi;
  weight = weight*Luminosity/100.;
  if(noLumiNorm) weight = 1.;
  if(noPUcorr)   wgt_pu = 1.;
  h_scattplot->Fill(varXNum[abs(vtNum[v])-2]/double(varXDen),varXNum[abs(vtNum[vv])-2]/double(varXDen),weight*wgt_pu);
 }
 return h_scattplot;
}
//darr_d
TH2F*  darr_d_iarr_d(string sample_vars, TTree* tree, int r,int v,int vv){
 TH2F *h_scattplot = emptyhisto(sample_vars,r,v,vv);
 vector<string> varNum(variablesNum, variablesNum + sizeof(variablesNum)/sizeof(variablesNum[0]));
 vector<string> varDen(variablesDen, variablesDen + sizeof(variablesDen)/sizeof(variablesDen[0]));
 double varXNum[100];
 double varXDen;
 int    varYNum[100];  
 tree->SetBranchAddress(varNum[v].c_str(),&varXNum);
 tree->SetBranchAddress(varDen[v].c_str(),&varXDen);
 tree->SetBranchAddress(varNum[vv].c_str(),&varYNum);
 double wgt_lumi = 0; double wgt_pu = 0;
 //Fill histo
 double weight = 0.;
 for(int en=0; en<tree->GetEntries(); en++)
 //for(int en=0; en<10; en++)
 {
  tree->GetEntry(en);
  weight = wgt_lumi;
  weight = weight*Luminosity/100.;
  if(noLumiNorm) weight = 1.;
  if(noPUcorr)   wgt_pu = 1.;
  if(varDen[v]==varDen[vv]){
   h_scattplot->Fill(varXNum[abs(vtNum[v])-2]/double(varXDen),varYNum[abs(vtNum[vv])-2]/double(varXDen),weight*wgt_pu);
  }else{
   double varYDen;
   tree->SetBranchAddress(varDen[vv].c_str(),&varYDen);
   h_scattplot->Fill(varXNum[abs(vtNum[v])-2]/double(varXDen),varYNum[abs(vtNum[vv])-2]/double(varYDen),weight*wgt_pu);
  }
 }
 return h_scattplot;
}
TH2F*  darr_d_darr_d(string sample_vars, TTree* tree, int r,int v,int vv){
 TH2F *h_scattplot = emptyhisto(sample_vars,r,v,vv);
 vector<string> varNum(variablesNum, variablesNum + sizeof(variablesNum)/sizeof(variablesNum[0]));
 vector<string> varDen(variablesDen, variablesDen + sizeof(variablesDen)/sizeof(variablesDen[0]));
 double varXNum[100];
 tree->SetBranchAddress(varNum[v].c_str(),&varXNum);
 double varYNum[100];   ///int///
 tree->SetBranchAddress(varNum[vv].c_str(),&varYNum);

 double wgt_lumi = 0; double wgt_pu = 0;
 //Fill histo
 double weight = 0.;
 for(int en=0; en<
 //100; en++)
 tree->GetEntries(); en++)
 {
  tree->GetEntry(en);
  weight = wgt_lumi;
  weight = weight*Luminosity/100.;
  if(noLumiNorm) weight = 1.;
  if(noPUcorr)   wgt_pu = 1.;
  h_scattplot->Fill(varXNum[abs(vtNum[v])-2],varYNum[abs(vtNum[vv])-2],weight*wgt_pu);
 }
 return h_scattplot;
}
TH2F*  darr_d_darr_iarr(string sample_vars, TTree* tree, int r,int v,int vv){
 TH2F *h_scattplot = emptyhisto(sample_vars,r,v,vv);
 vector<string> varNum(variablesNum, variablesNum + sizeof(variablesNum)/sizeof(variablesNum[0]));
 vector<string> varDen(variablesDen, variablesDen + sizeof(variablesDen)/sizeof(variablesDen[0]));
 double varXNum[100];
 double varXDen;
 double varYNum[100];
 int    varYDen[100];  ///int///
 tree->SetBranchAddress(varNum[v].c_str(),&varXNum);
 tree->SetBranchAddress(varDen[v].c_str(),&varXDen);
 tree->SetBranchAddress(varNum[vv].c_str(),&varYNum);
 tree->SetBranchAddress(varDen[vv].c_str(),&varYDen);
 double wgt_lumi = 0; double wgt_pu = 0;
 //Fill histo
 double weight = 0.;
 for(int en=0; en<tree->GetEntries(); en++)
 {
  tree->GetEntry(en);
  weight = wgt_lumi;
  weight = weight*Luminosity/100.;
  if(noLumiNorm) weight = 1.;
  if(noPUcorr)   wgt_pu = 1.;
  h_scattplot->Fill(varXNum[abs(vtNum[v])-2]/double(varXDen),varYNum[abs(vtNum[vv])-2]/double(varYDen[abs(vtDen[vv])-2]),weight*wgt_pu);
 }
 return h_scattplot;
}
TH2F*  darr_d_darr_darr(string sample_vars, TTree* tree, int r,int v,int vv){
 TH2F *h_scattplot = emptyhisto(sample_vars,r,v,vv);
 vector<string> varNum(variablesNum, variablesNum + sizeof(variablesNum)/sizeof(variablesNum[0]));
 vector<string> varDen(variablesDen, variablesDen + sizeof(variablesDen)/sizeof(variablesDen[0]));
 double varXNum[100];
 double varXDen;
 double varYNum[100];
 double varYDen[100];
 tree->SetBranchAddress(varNum[v].c_str(),&varXNum);
 tree->SetBranchAddress(varDen[v].c_str(),&varXDen);
 tree->SetBranchAddress(varNum[vv].c_str(),&varYNum);
 tree->SetBranchAddress(varDen[vv].c_str(),&varYDen);
 double wgt_lumi = 0; double wgt_pu = 0;
 //Fill histo
 double weight = 0.;
 for(int en=0; en<tree->GetEntries(); en++)
 {
  tree->GetEntry(en);
  weight = wgt_lumi;
  weight = weight*Luminosity/100.;
  if(noLumiNorm) weight = 1.;
  if(noPUcorr)   wgt_pu = 1.;
  h_scattplot->Fill(varXNum[abs(vtNum[v])-2]/double(varXDen),varYNum[abs(vtNum[vv])-2]/double(varYDen[abs(vtDen[vv])-2]),weight*wgt_pu);
 }
 return h_scattplot;
}
//darr_darr
TH2F*  darr_darr_darr_darr(string sample_vars, TTree* tree, int r,int v,int vv){
 TH2F *h_scattplot = emptyhisto(sample_vars,r,v,vv);
 vector<string> varNum(variablesNum, variablesNum + sizeof(variablesNum)/sizeof(variablesNum[0]));
 vector<string> varDen(variablesDen, variablesDen + sizeof(variablesDen)/sizeof(variablesDen[0]));
 double varXNum[100];
 double varXDen[100];
 tree->SetBranchAddress(varNum[v].c_str(),&varXNum);
 tree->SetBranchAddress(varDen[v].c_str(),&varXDen);
 double wgt_lumi = 0; double wgt_pu = 0;
 //Fill histo
 double weight = 0.;
 for(int en=0; en<tree->GetEntries(); en++)
 {
  tree->GetEntry(en);
  weight = wgt_lumi;
  weight = weight*Luminosity/100.;
  if(noLumiNorm) weight = 1.;
  if(noPUcorr)   wgt_pu = 1.;
  if(varNum[v]==varNum[vv] && varDen[v]==varDen[vv]){ 
   h_scattplot->Fill(varXNum[abs(vtNum[v])-2]/double(varXDen[abs(vtDen[v])-2]),varXNum[abs(vtNum[vv])-2]/double(varXDen[abs(vtDen[vv])-2]),weight*wgt_pu);
  }else if(varNum[v]!=varNum[vv] && varDen[v]==varDen[vv]){
   double varYNum[100];
   tree->SetBranchAddress(varNum[vv].c_str(),&varYNum);
   h_scattplot->Fill(varXNum[abs(vtNum[v])-2]/double(varXDen[abs(vtDen[v])-2]),varYNum[abs(vtNum[vv])-2]/double(varXDen[abs(vtDen[vv])-2]),weight*wgt_pu);
  }else{
   double varYNum[100];
   double varYDen[100];
   tree->SetBranchAddress(varDen[vv].c_str(),&varYDen);
   tree->SetBranchAddress(varNum[vv].c_str(),&varYNum);
   h_scattplot->Fill(varXNum[abs(vtNum[v])-2]/double(varXDen[abs(vtDen[v])-2]),varYNum[abs(vtNum[vv])-2]/double(varYDen[abs(vtDen[vv])-2]),weight*wgt_pu);
  }
 }
 return h_scattplot;
}
TH2F*  darr_darr_darr_iarr(string sample_vars, TTree* tree, int r,int v,int vv){
 TH2F *h_scattplot = emptyhisto(sample_vars,r,v,vv);
 vector<string> varNum(variablesNum, variablesNum + sizeof(variablesNum)/sizeof(variablesNum[0]));
 vector<string> varDen(variablesDen, variablesDen + sizeof(variablesDen)/sizeof(variablesDen[0]));
 double varXNum[100];
 double varXDen[100];
 double varYNum[100];
 int    varYDen[100];   ///int///
 tree->SetBranchAddress(varNum[v].c_str(),&varXNum);
 tree->SetBranchAddress(varDen[v].c_str(),&varXDen);
 tree->SetBranchAddress(varNum[vv].c_str(),&varYNum);
 tree->SetBranchAddress(varDen[vv].c_str(),&varYDen);
 double wgt_lumi = 0; double wgt_pu = 0;
 //Fill histo
 double weight = 0.;
 for(int en=0; en<tree->GetEntries(); en++)
 {
  tree->GetEntry(en);
  weight = wgt_lumi;
  weight = weight*Luminosity/100.;
  if(noLumiNorm) weight = 1.;
  if(noPUcorr)   wgt_pu = 1.;
  h_scattplot->Fill(varXNum[abs(vtNum[v])-2]/double(varXDen[abs(vtDen[v])-2]),varYNum[abs(vtNum[vv])-2]/double(varYDen[abs(vtDen[vv])-2]),weight*wgt_pu);
 }
 return h_scattplot;
}
//darr_iarr
TH2F*  darr_iarr_darr_iarr(string sample_vars, TTree* tree, int r,int v,int vv){
 TH2F *h_scattplot = emptyhisto(sample_vars,r,v,vv);
 vector<string> varNum(variablesNum, variablesNum + sizeof(variablesNum)/sizeof(variablesNum[0]));
 vector<string> varDen(variablesDen, variablesDen + sizeof(variablesDen)/sizeof(variablesDen[0]));
 double varXNum[100];
 int    varXDen[100];   ///int///
 double varYNum[100];
 tree->SetBranchAddress(varNum[v].c_str(),&varXNum);
 tree->SetBranchAddress(varDen[v].c_str(),&varXDen);
 tree->SetBranchAddress(varNum[vv].c_str(),&varYNum);
 double wgt_lumi = 0; double wgt_pu = 0;
 //Fill histo
 double weight = 0.;
 for(int en=0; en<tree->GetEntries(); en++)
 {
  tree->GetEntry(en);
  weight = wgt_lumi;
  weight = weight*Luminosity/100.;
  if(noLumiNorm) weight = 1.;
  if(noPUcorr)   wgt_pu = 1.;
  if(varDen[v]==varDen[vv]){
   h_scattplot->Fill(varXNum[abs(vtNum[v])-2]/double(varXDen[abs(vtDen[v])-2]),varYNum[abs(vtNum[vv])-2]/double(varXDen[abs(vtDen[vv])-2]),weight*wgt_pu);
  }else{
   int    varYDen[100];   ///int///
   tree->SetBranchAddress(varDen[vv].c_str(),&varYDen);
   h_scattplot->Fill(varXNum[abs(vtNum[v])-2]/double(varXDen[abs(vtDen[v])-2]),varYNum[abs(vtNum[vv])-2]/double(varYDen[abs(vtDen[vv])-2]),weight*wgt_pu);
  }
 }
 return h_scattplot;
}
TH2F*  darr_iarr_darr_darr(string sample_vars, TTree* tree, int r,int v,int vv){
 TH2F *h_scattplot = emptyhisto(sample_vars,r,v,vv);
 vector<string> varNum(variablesNum, variablesNum + sizeof(variablesNum)/sizeof(variablesNum[0]));
 vector<string> varDen(variablesDen, variablesDen + sizeof(variablesDen)/sizeof(variablesDen[0]));
 double varXNum[100];
 int    varXDen[100];   ///int///
 double varYDen[100];
 tree->SetBranchAddress(varNum[v].c_str(),&varXNum);
 tree->SetBranchAddress(varDen[v].c_str(),&varXDen);
 tree->SetBranchAddress(varDen[vv].c_str(),&varYDen);
 double wgt_lumi = 0; double wgt_pu = 0;
 //Fill histo
 double weight = 0.;
 for(int en=0; en<tree->GetEntries(); en++)
 {
  tree->GetEntry(en);
  weight = wgt_lumi;
  weight = weight*Luminosity/100.;
  if(noLumiNorm) weight = 1.;
  if(noPUcorr)   wgt_pu = 1.;
  if(varNum[v]!=varNum[vv]){
   double varYNum[100];
   tree->SetBranchAddress(varNum[vv].c_str(),&varYNum);
   h_scattplot->Fill(varXNum[abs(vtNum[v])-2]/double(varXDen[abs(vtDen[v])-2]),varYNum[abs(vtNum[vv])-2]/double(varYDen[abs(vtDen[vv])-2]),weight*wgt_pu);
  }else{ 
   h_scattplot->Fill(varXNum[abs(vtNum[v])-2]/double(varXDen[abs(vtDen[v])-2]),varXNum[abs(vtNum[vv])-2]/double(varYDen[abs(vtDen[vv])-2]),weight*wgt_pu);
  }
 }
 return h_scattplot;
}
/////
//   Draw scatter plot histo 
/////
void drawscattplot(TH2F* h){
 const Int_t NRGBs = 7;
 const Int_t NCont = 255;
 Double_t stops[NRGBs] = { 0.00, 0.15, 0.30, 0.45, 0.65, 0.85, 1.00 };
 Double_t red[NRGBs]   = { 0.60, 0.30, 0.00, 0.00, 0.60, 0.40, 0.00 };
 Double_t green[NRGBs] = { 1.00, 0.90, 0.80, 0.75, 0.20, 0.00, 0.00 };
 Double_t blue[NRGBs]  = { 1.00, 1.00, 1.00, 0.30, 0.00, 0.00, 0.00 };
 TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
 gStyle->SetNumberContours(NCont);
 TGaxis::SetMaxDigits(3);
 h->Draw("Scat");
 h->Draw("COLZ");
}
/////
//   Draw correlation matrix 
/////
void drawcorrmat(TH2F* h){
 const Int_t NRGBs = 7;
 const Int_t NCont = 255;
 Double_t stops[NRGBs] = { 0.00, 0.15, 0.30, 0.45, 0.65, 0.85, 1.00 };
 Double_t red[NRGBs]   = { 0.60, 0.30, 0.00, 0.00, 0.60, 0.40, 0.00 };
 Double_t green[NRGBs] = { 1.00, 0.90, 0.80, 0.75, 0.20, 0.00, 0.00 };
 Double_t blue[NRGBs]  = { 1.00, 1.00, 1.00, 0.30, 0.00, 0.00, 0.00 };
 TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
 gStyle->SetNumberContours(NCont);
 h->Draw("Scat");
 h->Draw("COLZ");
 h->Draw("textsame");
}
/////
//   Return an empty histos with all sets on axis
/////
TH2F* emptyhisto(string sample_vars,int r,int v,int vv){
 vector<string> varTitles(titles, titles + sizeof(titles)/sizeof(titles[0]));
 TH2F* h_scattplot = new TH2F(sample_vars.c_str(),sample_vars.c_str(),bin[v],inRange[v],endRange[v],bin[vv],inRange[vv],endRange[vv]);
 string titleY;
 if(r==0) titleY = "TTJets";
 if(r==1) titleY = "TTH";
 h_scattplot->SetTitle(titleY.c_str());
 h_scattplot->SetMarkerStyle(10);
 h_scattplot->SetMarkerColor(1);
 char bin_size_c[1000]; float bin_size_f = ((endRange[v]-inRange[v])/bin[v]); sprintf(bin_size_c,"%.2f",bin_size_f);
 string titlexaxis = varTitles[v];//+"/"+(string) bin_size_c;
 h_scattplot->GetXaxis()->SetTitle(titlexaxis.c_str());
 bin_size_f = ((endRange[vv]-inRange[vv])/bin[vv]); sprintf(bin_size_c,"%.2f",bin_size_f);
 string titleYaxis = varTitles[vv];//+"/"+(string) bin_size_c;
 h_scattplot->GetYaxis()->SetTitle(titleYaxis.c_str()); 
 return h_scattplot;
}
/////
//   Set setTDRStyle_modified (from link https://twiki.cern.ch/twiki/pub/CMS/TRK10001/setTDRStyle_modified.C)
/////
void setTDRStyle(){
  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

  // For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600); //Height of canvas
  tdrStyle->SetCanvasDefW(600); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);
  // For the Pad:
  tdrStyle->SetPadBorderMode(0);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);
  // For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);

  // For the histo:
  tdrStyle->SetHistFillColor(0);
  // tdrStyle->SetHistFillStyle(0);
  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);
  tdrStyle->SetErrorX(0.);

  tdrStyle->SetMarkerStyle(20);

  //For the fit/function:
  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

  //For the date:
  tdrStyle->SetOptDate(0);

  // For the statistics box:
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat(""); // To display the mean and RMS:   SetOptStat("mr");
  tdrStyle->SetStatColor(kGray);
  tdrStyle->SetStatFont(42);

  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(0);
  tdrStyle->SetStatX(0.9); //Starting position on X axis
  tdrStyle->SetStatY(1.); //Starting position on Y axis
  tdrStyle->SetStatFontSize(0.025); //Vertical Size
  tdrStyle->SetStatW(0.15); //Horizontal size 

  // Margins:
  tdrStyle->SetPadTopMargin(0.1);
  tdrStyle->SetPadBottomMargin(0.2);
  tdrStyle->SetPadLeftMargin(0.16);
  tdrStyle->SetPadRightMargin(0.14);

  // For the Global title:

  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);
  tdrStyle->SetTitleBorderSize(0);

  // For the axis titles:
  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.06, "XYZ");
  tdrStyle->SetTitleXOffset(1.4);
  tdrStyle->SetTitleYOffset(1.2);

  // For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.05, "XYZ");

  // For the axis:

  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

  // Change for log plots:
  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);
  tdrStyle->SetPaintTextFormat("4.2f");

  tdrStyle->cd();
}

