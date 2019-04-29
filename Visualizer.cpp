#include "Visualizer.hh"

// C/C++ includes
#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;

// ROOT includes
#include <TFrame.h>
#include <TKey.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TTimer.h>
#include <Getline.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TMath.h>
#include <TF1.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <Math/DistFunc.h>
#include <TPRegexp.h>
#include <THashTable.h>
#include <TGraph.h>

Visualizer::Visualizer(string *_files,int _n_files,int _bin,int _flags) : files(_files),n_files(_n_files),bin(_bin),flags(_flags),canvas(NULL) {
  io=new IO*[n_files];
  for(int i=0;i<n_files;i++) io[i]=new IO(files[i]);
}

Visualizer::~Visualizer() {
  
}

void Visualizer::prompt() {
  TTimer  *timer = new TTimer("gSystem->ProcessEvents();",50,kFALSE);
  TString input = "";
  while (input != "exit" && input != "quit") {
    string chrom = "",tmp="";
    int start,end;
    istringstream sin(input.Data());
    sin>>tmp;
    if(tmp!="") {
      cout << 1 << endl;
      int pos=tmp.find(":");
      chrom=tmp.substr(0,pos);
      tmp.erase(0,pos+1);
      pos=tmp.find("-");
      start=stoi(tmp.substr(0,pos));
      tmp.erase(0,pos+1);
      end=stoi(tmp);
      generateViewBAF(chrom,start,end);
    }
    timer->TurnOn();
    timer->Reset();
    input = Getline(">");
    input.ReplaceAll("\n","\0");
    input = input.Remove(TString::kBoth,' ');
    timer->TurnOff();
  }
  delete timer;
}

void Visualizer::generateViewBAF(string chrom,int start,int end) {
  int win = 1*(end - start + 1);
  if (win < 10000) win = 10000;
  
  if (canvas) delete canvas;
  TStyle *st = new TStyle("st","st");
  st->SetOptStat(false);  // No box with statistics
  st->SetOptTitle(false); // No box with title
  gROOT->SetStyle("st");
  canvas = new TCanvas("canv","canv",900,900);
  canvas->SetFillColor(kWhite);
  canvas->SetBorderMode(0); // No borders
  canvas->Divide(1,3);
  
  TString title = chrom; title += ":";
  title += start; title += "-";
  title += end;
  canvas->SetTitle(title);
  
  TH1 *raw=io[0]->getSignal(chrom,bin,"RD raw",0);
  TH1 *his=io[0]->getSignal(chrom,bin,"RD",0);
  TH1 *hisc=io[0]->getSignal(chrom,bin,"RD",flags);
  TH1 *hisp=io[0]->getSignal(chrom,bin,"RD partition",flags);
  TH1 *hism=io[0]->getSignal(chrom,bin,"RD call",flags);
  cout << his << endl;
  TVirtualPad *pad = canvas->cd(1);
  pad->SetFillColor(kWhite);
  pad->SetLineColor(kWhite);
  pad->SetFrameLineColor(kWhite);
  pad->SetFrameBorderMode(0);
  drawHistograms(chrom,start,end,win,"",pad,raw,his,hisc,hisp,hism);

//  if (hisc) {
//    double mean,sigma;
//    TH1 *rd_his = getHistogram(getDistrName("chr1",bin_size,
//                                            useATcorr,useGCcorr));
//    getMeanSigma(rd_his,mean,sigma);
//    hisc->GetYaxis()->SetRangeUser(0,2*mean);
//    hisc->SetLineWidth(1);
//  }

  TH1 *hisl=io[0]->getSignal(chrom,bin,"SNP likelihood",flags);
  pad = canvas->cd(2);
  pad->SetFillColor(kWhite);
  pad->SetLineColor(kWhite);
  pad->SetFrameLineColor(kWhite);
  pad->SetFrameBorderMode(0);
  drawHistogram2D(chrom,start,end,win,"",pad,hisl);
  
  TH1 *hisbaf=io[0]->getSignal(chrom,bin,"SNP i1",flags);
  TH1 *hisbafc=io[0]->getSignal(chrom,bin,"SNP i1 partition",flags);
  pad = canvas->cd(3);
  pad->SetFillColor(kWhite);
  pad->SetLineColor(kWhite);
  pad->SetFrameLineColor(kWhite);
  pad->SetFrameBorderMode(0);
  drawHistogramsBAF(chrom,start,end,win,"",pad,hisbaf,hisbafc);

  canvas->cd();
  canvas->Update();
}

void Visualizer::drawHistograms(TString chrom,int start,int end,
            int win,TString title,
            TVirtualPad *pad,
            TH1* raw,TH1 *his,TH1 *hisc,TH1 *hisp,TH1 *hism)
{
  TH1 *main = hisc;
  if (!main) main = his;
  if (!main) main = raw;
  if (!main) return;

  int s  = start - win, e  = end + win;
  int n_bins = main->GetNbinsX();
  int bs = s/bin - 1; if (bs < 0) bs = 1;
  int be = e/bin + 1; if (be > n_bins) be = n_bins;
  double max = 0;
  for (int i = bs;i <= be;i++) {
    if (raw  && raw->GetBinContent(i)  > max) max = raw->GetBinContent(i);
    if (his  && his->GetBinContent(i)  > max) max = his->GetBinContent(i);
    if (hisc && hisc->GetBinContent(i) > max) max = hisc->GetBinContent(i);
  }
  max *= 1.05;

  main->Draw();
  main->GetXaxis()->SetRangeUser(s,e);
  main->GetXaxis()->SetTitle(chrom);
  main->GetYaxis()->SetRangeUser(0,max);
  main->GetYaxis()->SetTitle(title);
  main->SetLineWidth(3);
  if (raw) {
    raw->Draw("same hist");
    raw->SetLineColor(kYellow);
  }
  if (his) {
    his->Draw("same hist");
    his->SetLineColor(kGray);
  }
  if (hisc) {
    hisc->Draw("same hist");
    hisc->SetLineColor(kBlack);
  }
  if (hisp) {
    hisp->Draw("same");
    hisp->SetLineColor(kRed);
    hisp->SetLineWidth(3);
  }
  if (hism) {
    hism->Draw("same hist");
    hism->SetLineColor(kGreen);
    hism->SetLineWidth(3);
  }
  TLine *line1 = new TLine(0,0,0,0),*line2 = new TLine(0,0,0,0);
  line1->SetX1(start); line1->SetX2(start);
  line1->SetY1(0);     line1->SetY2(max);
  line2->SetX1(end);   line2->SetX2(end);
  line2->SetY1(0);     line2->SetY2(max);
  line1->SetLineColor(kCyan);
  line2->SetLineColor(kCyan);
  line1->SetLineWidth(3);
  line2->SetLineWidth(3);
  line1->Draw(); line2->Draw();
}

void Visualizer::drawHistogram2D(TString chrom,int start,int end,
            int win,TString title,
            TVirtualPad *pad,
            TH1* his2d)
{
  TH1 *main = his2d;

  int s  = start - win, e  = end + win;
  int n_bins = main->GetNbinsX();
  int bs = s/bin - 1; if (bs < 0) bs = 1;
  int be = e/bin + 1; if (be > n_bins) be = n_bins;
  double max = 1.05;

  main->Draw("CONT4");
  main->GetXaxis()->SetRangeUser(s,e);
  main->GetXaxis()->SetTitle(chrom);
  main->GetYaxis()->SetRangeUser(0,max);
  main->GetYaxis()->SetTitle(title);
  main->SetLineWidth(3);

  TLine *line1 = new TLine(0,0,0,0),*line2 = new TLine(0,0,0,0);
  line1->SetX1(start); line1->SetX2(start);
  line1->SetY1(0);     line1->SetY2(max);
  line2->SetX1(end);   line2->SetX2(end);
  line2->SetY1(0);     line2->SetY2(max);
  line1->SetLineColor(kCyan);
  line2->SetLineColor(kCyan);
  line1->SetLineWidth(3);
  line2->SetLineWidth(3);
  line1->Draw(); line2->Draw();
}

void Visualizer::drawHistogramsBAF(TString chrom,int start,int end,
            int win,TString title,
            TVirtualPad *pad,
            TH1* hiss,TH1 *hisc)
{
  TH1 *main = hiss;
  if (!main) main = hisc;

  int s  = start - win, e  = end + win;
  int n_bins = main->GetNbinsX();
  int bs = s/bin - 1; if (bs < 0) bs = 1;
  int be = e/bin + 1; if (be > n_bins) be = n_bins;
  double max = 0.5;

  main->Draw();
  main->GetXaxis()->SetRangeUser(s,e);
  main->GetXaxis()->SetTitle(chrom);
  main->GetYaxis()->SetRangeUser(0,max);
  main->GetYaxis()->SetTitle(title);
  main->SetLineWidth(3);
  if (hisc) {
    hisc->Draw("same hist");
    hisc->SetLineColor(kGreen);
    hisc->SetLineWidth(3);
  }
  TLine *line1 = new TLine(0,0,0,0),*line2 = new TLine(0,0,0,0);
  line1->SetX1(start); line1->SetX2(start);
  line1->SetY1(0);     line1->SetY2(max);
  line2->SetX1(end);   line2->SetX2(end);
  line2->SetY1(0);     line2->SetY2(max);
  line1->SetLineColor(kCyan);
  line2->SetLineColor(kCyan);
  line1->SetLineWidth(3);
  line2->SetLineWidth(3);
  line1->Draw(); line2->Draw();
}
