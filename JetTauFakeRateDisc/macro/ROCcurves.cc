/**
This Macro   
1. Plots ROC curves 

Need to specify
1. See Declare Constants
*/
/////
//   To run: root -l ROCcurves.cc+
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
#include "TMath.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <iostream>
using namespace std;
/////
//   Declare constants
/////
const string path     = "../Rootplas/";
const char *samples[] = {"bkg", "GGHLL_b"};
const string selection = "";
const int numvar       = 100;

const char *varfirst[]        = {"pftauchhads_IP2D_val[0]", "pftauchhads_IP3D_val[0]"};
const char *varsecond[]       = {"1", "1"}; 
//,"jettrksallch_num2vno2v[0]","jettrksallch_num2vno2v[1]","jettrksallch_num2vno2v[2]","jettrksallch_num2vno2v[3]","jettrksallch_num2vno2v[4]"};
const char *vartitle[]        = {"IP2D_val", "IP3D_val"};
const double cutIni[numvar]  = {0,0,};
const double cutFin[numvar]  = {3.5,35};
const int    npoints[numvar] = {100,100};
const string namefile        = "ROC_IP_val.pdf";


/////
//   Declare functions 
/////
TFile* Call_TFile(string rootpla);
void setTDRStyle();
/////
//   Main function
/////
void ROCcurves(){
 //Preliminarly
 setTDRStyle();
 vector<string> varfirsts(varfirst, varfirst + sizeof(varfirst)/sizeof(varfirst[0])); 
 vector<string> varseconds(varsecond, varsecond + sizeof(varsecond)/sizeof(varsecond[0]));
 vector<string> vartitles(vartitle, vartitle + sizeof(vartitle)/sizeof(vartitle[0])); 
 const uint variables_size = varfirsts.size();
 TCanvas* c1 = new TCanvas("ROC","ROC",200,200,1100,700);
 TLegend *leg = new TLegend(0.65, 0.2, 0.9, 0.7);
 leg->SetHeader("Variables");
 leg->SetBorderSize(0);
 leg->SetTextSize(0.035);
 //Loop over variables
 for(uint vars=0; vars<variables_size; vars++) 
 //for(uint vars=0; vars<1; vars++)
 {
  cout<<"Var "<<varfirsts[vars]<<endl;
  double cutStep = (cutFin[vars]-cutIni[vars])/npoints[vars];
  vector<string> rootplas(samples, samples + sizeof(samples)/sizeof(samples[0]));
  const uint rootplas_size = rootplas.size();
  const int npointsint = npoints[vars];   
  double effvals[rootplas_size][npointsint];
  double efferrs[rootplas_size][npointsint]; 
  //Loop over samples
  for(uint smp=0; smp<rootplas_size; smp++){
   cout<<rootplas[smp]<<endl;
   TFile* f = Call_TFile((rootplas[smp]).c_str()); TTree* tree; f->GetObject("tree",tree);
   double evt_bef = tree->GetEntries();
   //double evt_bef = tree->GetEntries(Form("%s>0",varfirsts[vars].c_str()));
   int cutpos = 0;
   for(double cut=cutIni[vars]; cut<=cutFin[vars]; cut+=cutStep){
    double evt_aft = 0;
    if(varseconds[vars]=="1") evt_aft = tree->GetEntries(Form("%s>=%f",varfirsts[vars].c_str(),cut));
    if(varseconds[vars]!="1") evt_aft = tree->GetEntries(Form("%s/%s>=%f",varfirsts[vars].c_str(),varseconds[vars].c_str(),cut));
    //if(varseconds[vars]=="1" && varfirsts[vars].find("chipval")>varfirsts[vars].length())
    // evt_aft = tree->GetEntries(Form("%s>=%f",varfirsts[vars].c_str(),cut));
    //if(varseconds[vars]=="1" && varfirsts[vars].find("chipval")<varfirsts[vars].length())
    // evt_aft = tree->GetEntries(Form("(-log(%s)>=%f && %s>0) || (%s<=0)",varfirsts[vars].c_str(),cut,varfirsts[vars].c_str(),varfirsts[vars].c_str()));
     //evt_aft = tree->GetEntries(Form("(-log(%s)>=%f && %s>0)",varfirsts[vars].c_str(),cut,varfirsts[vars].c_str()));
    double eff = evt_aft/evt_bef;
    effvals[smp][cutpos] = eff; 
    efferrs[smp][cutpos] = 0;  
    //cout<<effvals[smp][cutpos]<<" ";
    cutpos++;
   }
   //cout<<endl;
  }//Loop over samples
  TGraphErrors *rocurves = new TGraphErrors(npointsint,&effvals[0][0],&effvals[1][0],&efferrs[0][0],&efferrs[1][0]);  
  rocurves->SetTitle(0);
  rocurves->GetXaxis()->SetTitle("Jet efficiency");
  rocurves->GetYaxis()->SetTitle("Tau efficiency");
  rocurves->SetMarkerSize(0.3);
  rocurves->SetMarkerColor(1+vars);
  rocurves->SetLineColor(1+vars);
  if(vars==0){
   rocurves->Draw("APL");
  }else{
   rocurves->Draw("same");
  }
  leg->AddEntry(rocurves,vartitles[vars].c_str(),"L"); 
 }//Loop over variables
 leg->Draw();
 gPad->RedrawAxis();
 c1->SaveAs(namefile.c_str());
}
/////
//   Call TFile to be read
/////
TFile* Call_TFile(string rootpla){
 string dotroot   = ".root";
 string file_name = path+rootpla+dotroot;
 //cout<<"file_name "<<file_name<<endl;
 TFile* f = new TFile(file_name.c_str(),"update");
 return f;
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
  tdrStyle->SetPadTopMargin(0.005);
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
  //tdrStyle->SetTitleSize(0.045, "X");
  tdrStyle->SetTitleSize(0.06, "XYZ");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  tdrStyle->SetTitleXOffset(1.);
  //tdrStyle->SetTitleYOffset(1.0);
  tdrStyle->SetTitleOffset(0.75, "Y"); // Another way to set the Offset

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
