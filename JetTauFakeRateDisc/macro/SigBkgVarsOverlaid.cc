/**
This Macro   
1. Plots variables of signal and background overlaid

Need to specify
1. See Declare Constants
2. There is a difference between int first[arraycomp]; int second[arraycomp]; and double first[arraycomp]; double second[arraycomp];
   There is also a difference to tak into account in the code in case you plot array or simple double/int
*/
/////
//   To run: root -l SigBkgVarsOverlaid.cc+
/////
/////
//   Prepare Root and Roofit
/////
#include "TFile.h"
#include "TH1F.h"
#include "THStack.h"
#include "TLegend.h"
#include "TTree.h"
#include "TTreePlayer.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TEfficiency.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <iostream>
using namespace std;
/////
//   Declare constants
/////
const string path      = "/home/francescoromeovb/Shared_win_ubu_/Work/JetTauFakeRateDisc/AEIP1D_3/";
//const char *samples[]  = {"Qcd", "Tau"};
const char *samples[]  = {"Tau"};
const string selection = "";
const int numvar       = 100;

/*
const char *varfirst[]        = {"pftauchhads_AEsIP1D_val_trk0", "pftauchhads_AEsIP1D_sig_trk0"};
const char *varsecond[]       = {"1", "1"};
const char *vartitle[]        = {"signed IP1D val (z-axis)", "signed IP1D sig (x-axis)"};
const double inRange[numvar]  = {-0.05,-50};
const double endRange[numvar] = {0.05,100};
const int    bin[numvar]      = {100,150};
bool ylogscale    = true;
double setminimum = 0.1;
double setmaximum = 100;
*/

/*
const char *varfirst[]        = {"pftauchhads_AEsIP1D_val_trk0", "pftauchhads_AEIP1D_x_val_trk0", "pftauchhads_AEIP1D_y_val_trk0"};
const char *varsecond[]       = {"1", "1", "1"};
const char *vartitle[]        = {"signed IP1D (z-axis)", "signed IP1D (x-axis)", "signed IP1D (y-axis)"};
const double inRange[numvar]  = {-0.05,-0.05,-0.05};
const double endRange[numvar] = {0.05,0.05,0.05};
const int    bin[numvar]      = {100,100,100};
bool ylogscale    = true;
double setminimum = 0.1;
double setmaximum = 100;
*/

//Difference of IP3D minus IP2D+IP1D
const char *varfirst[]        = {"pftauchhads_IP3DvalAE_2DAE_1DAE_trk0", "pftauchhads_IP3DvalAE_2DTE_1DAE_trk0", "pftauchhads_IP3DvalAE_2DTE_1DTE_trk0"};
const char *varsecond[]       = {"pftauchhads_IP3DvalAE_trk0", "pftauchhads_IP3DvalAE_trk0", "pftauchhads_IP3DvalAE_trk0"};
const char *vartitle[]        = {"#frac{AE_IP3D - #sqrt{AE_IP2D^{2}+AE_IP1D^{2}}}{AE_IP3D}", "#frac{AE_IP3D - #sqrt{TE_IP2D^{2}+AE_IP1D^{2}}}{AE_IP3D}", "#frac{AE_IP3D - #sqrt{TE_IP2D^{2}+TE_IP1D^{2}}}{AE_IP3D}"};
const double inRange[numvar]  = {-0.1,-0.1,-7};
const double endRange[numvar] = {0.1,  1, 3};
const int    bin[numvar]      = {50,  50,50};
bool ylogscale    = true;
double setminimum = 0.1;
double setmaximum = 100;
/*
const char *varfirst[]        = {"pftauchhads_IP3DvalAE_2DAE_1DAE_trk0", "pftauchhads_IP3DvalAE_2DTE_1DAE_trk0", "pftauchhads_IP3DvalAE_2DTE_1DTE_trk0"};
const char *varsecond[]       = {"1", "1", "1"};
const char *vartitle[]        = {"AE_IP3D - #sqrt{AE_IP2D^{2}+AE_IP1D^{2}}", "AE_IP3D - #sqrt{TE_IP2D^{2}+AE_IP1D^{2}}", "AE_IP3D - #sqrt{TE_IP2D^{2}+TE_IP1D^{2}}"};
const double inRange[numvar]  = {-5,-5,-5};
const double endRange[numvar] = {5,5,5};
const int    bin[numvar]      = {50,50,50};
bool ylogscale    = true;
double setminimum = 0.1;
double setmaximum = 100;
*/

//TauVSTauJet
/*
const char *varfirst[]        = {"tau_pt_DIV_recojettau_pt"};
const char *varsecond[]       = {"1"};
const char *vartitle[]        = {"pT_{#tau}/pT_{jet#tau}"};
const double inRange[numvar]  = {0};
const double endRange[numvar] = {1};
const int    bin[numvar]      = {50};
bool ylogscale    = true;
double setminimum = 0.1;
double setmaximum = 100;
*/
//Impact parameters
/*
const char *varfirst[]        = {"pftauchhads_IP3D_val_trk0", "pftauchhads_AEIP1D_val_trk0", "pftauchhads_IP3D_sig_trk0", "pftauchhads_AEIP1D_sig_trk0"};
const char *varsecond[]       = {"1", "1", "1", "1", "1", "1"};
const char *vartitle[]        = {"sIP3D val trk0", "AEsIP1D val trk0", "sIP3D sig trk0", "AEsIP1D sig trk0"};
const double inRange[numvar]  = {0,0,0,0,0,0};
//const double inRange[numvar]  = {-0.05,-0.05,-10,-10};
const double endRange[numvar] = {0.05,0.05,10,10};
const int    bin[numvar]      = {100,100,100,100};
bool ylogscale    = true;
double setminimum = 0.1;
double setmaximum = 100;
*/

/*
const char *varfirst[]        = {"pftauchhads_sIP3D_val_trk0", "pftauchhads_sIP2D_val_trk0", "pftauchhads_sIP1D_val_trk0", "pftauchhads_sIP3D_sig_trk0", "pftauchhads_sIP2D_sig_trk0", "pftauchhads_sIP1D_sig_trk0"};
const char *varsecond[]       = {"1", "1", "1", "1", "1", "1"};
const char *vartitle[]        = {"sIP3D val trk0", "sIP2D val trk0", "sIP1D val trk0", "sIP3D sig trk0", "sIP2D sig trk0", "sIP1D sig trk0"};
//const double inRange[numvar]  = {0,0,0,0,0,0};
const double inRange[numvar]  = {-0.05,-0.05,-0.05,-10,-10,-10};
const double endRange[numvar] = {0.05,0.05,0.05,10,10,10};
const int    bin[numvar]      = {100,100,100,100,100,100};
bool ylogscale    = true;
double setminimum = 0.1;
double setmaximum = 100;
*/

/*
//Decay lenght
const char *varfirst[]        = {"pftauchhads_absDL3D_val_trk0", "pftauchhads_absDL3D_val_trk1", "pftauchhads_absDL3D_val_trk2", "pftauchhads_absDL3D_sig_trk0", "pftauchhads_absDL3D_sig_trk1", "pftauchhads_absDL3D_sig_trk2"};
const char *varsecond[]       = {"1", "1", "1", "1", "1", "1"};
const char *vartitle[]        = {"absDL3D val trk0", "absDL3D val trk1", "absDL3D val trk2", "absDL3D sig trk0", "absDL3D sig trk1", "absDL3D sig trk2"};
const double inRange[numvar]  = {0,0,0,0,0,0};
const double endRange[numvar] = {7,7,7,7000,7000,7000};
const int    bin[numvar]      = {100,100,100,100,100,100};
bool ylogscale    = true;
double setminimum = 0.1;
double setmaximum = 100;
*/

/*
//Flight distance
const char *varfirst[]        = {"pftau_pvsv_dist3d_val", "pftau_pvsv_dist3d_sig", "pftau_pvsv_dist2d_val", "pftau_pvsv_dist2d_sig"};
const char *varsecond[]       = {"1", "1", "1", "1", "1", "1"};
const char *vartitle[]        = {"Flight Dist 3D val", "Flight Dist 3D sig", "Flight Dist 2D val", "Flight Dist 2D sig"};
const double inRange[numvar]  = {0,0,0,0,0,0};
const double endRange[numvar] = {5,15,5,15};
const int    bin[numvar]      = {100,100,100,100};
bool ylogscale    = true;
double setminimum = 0.1;
double setmaximum = 100;
*/

//Tau Collimation
/*
const char *varfirst[]        = {"dR_tautrk_recotau_trk0", "dR_tautrk_recotau_trk1", "dR_tautrk_recotau_trk2"};//, "pftauchhads_JTD_val_trk0", "pftauchhads_JTD_val_trk1", "pftauchhads_JTD_val_trk2"};
const char *varsecond[]       = {"1", "1", "1", "1", "1", "1"};
const char *vartitle[]        = {"dR(Trk0,RecoTau)", "dR(Trk1,RecoTau)", "dR(Trk2,RecoTau)", "dist(Trk0,RecoTau)", "dist(Trk1,RecoTau)", "dist(Trk2,RecoTau)",};
const double inRange[numvar]  = {0,0,0,0,0,0};
const double endRange[numvar] = {0.15,0.15,0.15,0.05,0.05,0.05};
const int    bin[numvar]      = {100,100,100,100,100,100,100};
bool ylogscale    = true;
double setminimum = 0.1;
double setmaximum = 100;
*/
//Vertex compatibility
/*
const char *varfirst[]        = {"diff_pv_avfbs_nchi2_unbpv_avfbs_nchi2"};
const char *varsecond[]       = {"1", "1", "1", "1", "1", "1"};
const char *vartitle[]        = {"#chi^{2}_{PV}-#chi^{2}_{PV-TauTrks}"};
const double inRange[numvar]  = {-0.1,0,0,0,0,0};
const double endRange[numvar] = {0.15,0.15,0.15,0.05,0.05,0.05};
const int    bin[numvar]      = {100,100,100,100,100,100,100};
bool ylogscale    = true;
double setminimum = 0.1;
double setmaximum = 100;
*/

/*
const char *varfirst[]        = {"pftau_pvsv_dist3d_val", "pftau_pvsv_dist3d_sig"};
const char *varsecond[]       = {"1", "1"};
const char *vartitle[]        = {"pftau_pvsv_dist3d_val", "pftau_pvsv_dist3d_sig"};
const double inrange[numvar]  = {0,0};
const double endrange[numvar] = {5,10};
const int    bin[numvar]      = {50,100};
bool ylogscale    = true;
double setminimum = 0.1;
double setmaximum = 100;
*/

/*
const char *varfirst[]       = {"pftauchhads_pt", "pftauchhads_pt", "pftauchhads_pt"};//, "pftauchhads_IP2D_sig"};
const char *varsecond[]      = {"1","1","1", "1","1","1"};
const char *vartitle[]       = {"tau pt 1", "tau pt 2", "tau pt 3"};
const double inRange[numvar]  = {0, 0, 0};
const double endRange[numvar] = {500, 500, 500};
const int    bin[numvar]      = {100, 100, 100};
bool ylogscale    = true;
double setminimum = 0.1;
double setmaximum = 100;
*/

/*
const char *varfirst[]       = {"pftauchhads_IP3D_sig", "pftauchhads_IP2D_sig"};
const char *varsecond[]      = {"1","1","1", "1","1","1"};
const char *vartitle[]       = {"IP3D_sig trk1", "IP2D_sig trk1"};
const double inRange[numvar]  = {0, 0, 0};
const double endRange[numvar] = {9, 8, 10};
const int    bin[numvar]      = {100, 100, 100};
bool ylogscale    = true;
double setminimum = 0.1;
double setmaximum = 100;
*/

/*
const char *varfirst[]        = {"dR_tautrk_recotau"};
const char *varsecond[]       = {"1"};
const char *vartitle[]        = {"dR_tautrk_recotau"};
const double inRange[numvar]  = {0};
const double endRange[numvar] = {0.05};
const int    bin[numvar]      = {150};
bool ylogscale    = true;
double setminimum = 0.1;
double setmaximum = 100;
*/

/*
const char *varfirst[]        = {"pvtauvtx_genditauvtx_x", "pvtauvtx_genditauvtx_y", "pvtauvtx_genditauvtx_z"};//, "unbpv_KVF_genditauvtx_x", "unbpv_KVF_genditauvtx_y", "unbpv_KVF_genditauvtx_z"};
const char *varsecond[]       = {"1", "1", "1", "1", "1", "1"};
const char *vartitle[]        = {"(recoVtx(AVF)-genVtx).x", "(recoVtx(AVF)-genVtx).y", "(recoVtx(AVF)-genVtx).z", "(recoVtx(KVF)-genVtx).x", "(recoVtx(KVF)-genVtx).y", "(recoVtx(KVF)-genVtx).z"};
const double inRange[numvar]  = {-0.1,-0.1,-20,-0.1,-0.1,-20};
const double endRange[numvar] = {0.1,0.1,20,0.1,0.1,20};
const int    bin[numvar]      = {200,200,400,200,200,400};
bool ylogscale    = true;
double setminimum = 0.01;
double setmaximum = 100;
*/

//Other options
bool saveplots = true;

/////
//   Declare functions 
/////
TFile* Call_TFile(string rootpla);
string convertToString(double number);
void setTDRStyle();
/////
//   Main function
/////
void SigBkgVarsOverlaid(){
 //Preliminarly
 setTDRStyle();
 //Loop over variables
 vector<string> varfirsts(varfirst, varfirst + sizeof(varfirst)/sizeof(varfirst[0])); 
 vector<string> varseconds(varsecond, varsecond + sizeof(varsecond)/sizeof(varsecond[0]));
 vector<string> vartitles(vartitle, vartitle + sizeof(vartitle)/sizeof(vartitle[0])); 
 const uint variables_size = varfirsts.size();
 for(uint vars=0; vars<variables_size; vars++)
 {
  TCanvas* c1 = new TCanvas(vartitles[vars].c_str(),vartitles[vars].c_str(),200,200,800,800);
  if(ylogscale) c1->SetLogy();
  TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);
  //Do plots
  //Loop over samples
  vector<string> rootplas(samples, samples + sizeof(samples)/sizeof(samples[0]));
  const uint rootplas_size = rootplas.size();
  for(uint smp=0; smp<rootplas_size; smp++){
   cout<<"Rootpla "<<rootplas[smp]<<endl;
   //Declare histo
   TH1F *hist = new TH1F(rootplas[smp].c_str(),rootplas[smp].c_str(),bin[vars],inRange[vars],endRange[vars]); 
   hist->SetTitle("");
   hist->SetMarkerStyle(1);
   hist->GetXaxis()->SetTitle(vartitles[vars].c_str());
   hist->GetYaxis()->SetTitle("Percentage");
   if(smp==0){
    hist->SetMarkerColor(kRed);
    hist->SetLineColor(kRed);
   }else if(smp==1){
    hist->SetMarkerColor(kBlue);
    hist->SetLineColor(kBlue);
   }
   //Make plot
   TFile* f = Call_TFile((rootplas[smp]).c_str()); TTree* tree; f->GetObject("tree",tree);
   double first;
   double second;
   double pftauchhads_num;
   //const int arraycomp = 10;
   //double first[arraycomp];
   //double second[arraycomp];
   tree->SetBranchAddress(varfirsts[vars].c_str(),&first);
   tree->SetBranchAddress(varseconds[vars].c_str(),&second);
   tree->SetBranchAddress("pftauchhads_num", &pftauchhads_num);
   //for(int en=0; en<tree->GetEntries(); en++)
   for(int en=0; en<100000; en++)
   {
    tree->GetEntry(en);
    //if(first==-999 || first==0) continue;
    if(first==-999) continue;
    if(pftauchhads_num==0) continue;       
    if((varfirsts[vars]=="dR_tautrk_recotau_trk0" && pftauchhads_num!=1) ||
       (varfirsts[vars]=="dR_tautrk_recotau_trk1" && pftauchhads_num!=2) ||
       (varfirsts[vars]=="dR_tautrk_recotau_trk2" && pftauchhads_num!=3)
      ) continue;
    if((varfirsts[vars]=="pftauchhads_JTD_val_trk0" || varfirsts[vars]=="pftauchhads_JTD_val_trk1" ||varfirsts[vars]=="pftauchhads_JTD_val_trk2")) first = -first;
    if(varseconds[vars]=="1") hist->Fill(first);         
    if(varseconds[vars]!="1") hist->Fill(double(first)/double(second));               
    //cout<<"Val "<<en<<" "<<first[vars]<<" "<<second[vars]<<" "<<double(second[vars])<<endl;
    //if(first[vars]==-999 || first[vars]==0) continue;
    //if(pftauchhads_num==0) continue;
    //if(varseconds[vars]=="1") hist->Fill(first[vars]);  
    //if(varseconds[vars]!="1") hist->Fill(double(first[vars])/double(second[vars]));  
   }
   cout<<"Entries of "<<rootplas[smp]<<" is "<<hist->Integral()<<endl;
   double scale = hist->Integral();
   hist->Scale(1/scale*100);
   //Draw plot
   gStyle->SetOptStat("mr");     
   gStyle->SetStatColor(kWhite);
   gStyle->SetStatX(0.9); //Starting position on X axis
   gStyle->SetStatW(0.2); //Horizontal size 
   //jet_csv
   //gStyle->SetStatX(0.7); //Starting position on X axis
   //gStyle->SetStatW(0.2); //Horizontal size 
   if(smp==0){
    gStyle->SetStatTextColor(kRed);
    gStyle->SetStatY(0.7); //Starting position on Y axis
    gStyle->SetStatFontSize(0.1); //Vertical Size
   }else if(smp==1){
    gStyle->SetStatTextColor(kBlue);
    gStyle->SetStatY(0.5); //Starting position on Y axis
    gStyle->SetStatFontSize(0.1); //Vertical Size
   }
   if(smp==0){
    leg->AddEntry(hist,"Tau","L"); 
    hist->SetMinimum(setminimum);
    hist->SetMaximum(setmaximum);
    hist->Draw(""); 
   }else{
    leg->AddEntry(hist,"Tau","L"); 
    hist->SetMinimum(setminimum);
    hist->SetMaximum(setmaximum);
    hist->Draw("sames");
   }
  } 
  leg->Draw();
  string num = convertToString(double(vars));
  string namefile = "SigBkgVarsOverlaid_"+varfirsts[vars]+".pdf";
  if(saveplots) c1->SaveAs(namefile.c_str());
 }
}
/////
//   Call TFile to be read
/////
TFile* Call_TFile(string rootpla){
 string dotroot   = ".root";
 string file_name = path+rootpla+dotroot;
 cout<<"file_name "<<file_name<<endl;
 TFile* f = new TFile(file_name.c_str(),"update");
 return f;
}
/////
//   Convert number to string 
/////
string convertToString(double number){
 char buffer[5];
 sprintf(buffer,"%g",number);
 return string(buffer);
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
  // tdrStyle->SetPadBorderSize(Width_t size = 1);
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
  // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
  // tdrStyle->SetNumberContours(Int_t number = 20);
//  tdrStyle->SetEndErrorSize(0);
  tdrStyle->SetErrorX(0.);
//  tdrStyle->SetErrorMarker(20);

  tdrStyle->SetMarkerStyle(20);

  //For the fit/function:
  tdrStyle->SetOptFit(0);
  tdrStyle->SetFitFormat("5.4g");
  //tdrStyle->SetFuncColor(1);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

  //For the date:
  tdrStyle->SetOptDate(0);
  // tdrStyle->SetDateX(Float_t x = 0.01);
  // tdrStyle->SetDateY(Float_t y = 0.01);

  // For the statistics box:
  tdrStyle->SetOptFile(0);
  //tdrStyle->SetOptStat("mr"); // To display the mean and RMS:   SetOptStat("mr");
  //tdrStyle->SetStatColor(kWhite);
  //tdrStyle->SetStatColor(kGray);
  //tdrStyle->SetStatFont(42);

  //tdrStyle->SetTextSize(11);
  tdrStyle->SetTextAlign(11);

  //tdrStyle->SetStatTextColor(1);
  //tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(0);
  //tdrStyle->SetStatX(1.); //Starting position on X axis
  //tdrStyle->SetStatY(1.); //Starting position on Y axis
  //tdrStyle->SetStatFontSize(0.025); //Vertical Size
  //tdrStyle->SetStatW(0.25); //Horizontal size 
  // tdrStyle->SetStatStyle(Style_t style = 1001);

  // Margins:
  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.165);
  tdrStyle->SetPadLeftMargin(0.1);
  tdrStyle->SetPadRightMargin(0.05);

  // For the Global title:
  //  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.5);
  tdrStyle->SetTitleH(0.05); // Set the height of the title box
  //tdrStyle->SetTitleW(0); // Set the width of the title box
  tdrStyle->SetTitleX(0.15); // Set the position of the title box
  tdrStyle->SetTitleY(1.0); // Set the position of the title box
  // tdrStyle->SetTitleStyle(Style_t style = 1001);
  tdrStyle->SetTitleBorderSize(0);

  // For the axis titles:
  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.045, "X");
  tdrStyle->SetTitleSize(0.06, "YZ");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  tdrStyle->SetTitleXOffset(1.5);
  //tdrStyle->SetTitleYOffset(1.0);
  tdrStyle->SetTitleOffset(0.75, "Y"); // Another way to set the Offset

  // For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.04, "XYZ");

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
  // Postscript options:
  // tdrStyle->SetPaperSize(15.,15.);
  // tdrStyle->SetLineScalePS(Float_t scale = 3);
  // tdrStyle->SetLineStyleString(Int_t i, const char* text);
  // tdrStyle->SetHeaderPS(const char* header);
  // tdrStyle->SetTitlePS(const char* pstitle);
  // tdrStyle->SetBarOffset(Float_t baroff = 0.5);
  // tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
  // tdrStyle->SetPaintTextFormat(const char* format = "g");
  // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // tdrStyle->SetTimeOffset(Double_t toffset);
  // tdrStyle->SetHistMinimumZero(kTRUE);

  tdrStyle->cd();
}

