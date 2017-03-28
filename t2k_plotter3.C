
{

  gROOT->SetStyle("Plain");
  gStyle->SetTitleBorderSize(0);
  gStyle->SetOptStat("");

  gStyle->SetLabelFont(102,"");
  gStyle->SetLabelSize(0.06,"");
  gStyle->SetLabelFont(102,"xyz");
  gStyle->SetLabelSize(0.034,"xyz");
  gStyle->SetLabelOffset(0.001,"x");
  gStyle->SetLabelOffset(0.01,"y");

  gStyle->SetTitleFont(42,"xyz");
  gStyle->SetTitleFontSize(0.06);
  gStyle->SetTitleFont(42);
  gStyle->SetTitleSize(0.06,"xyz");
  gStyle->SetTitleOffset(1.05,"x");
  gStyle->SetTitleOffset(0.8,"y");

  gStyle->SetStripDecimals(kFALSE);
  
  gStyle->SetPadLeftMargin(0.20);
  gStyle->SetPadBottomMargin(0.15);

  gStyle->SetStatW(0.35);
  gStyle->SetStatH(0.25);

  gStyle->SetPadTickX(kTRUE);
  gStyle->SetPadTickY(kTRUE);

  gStyle->SetPalette(1);
  gStyle->SetNumberContours(99);

  gStyle->SetHistLineWidth(2);
  gStyle->SetFrameLineWidth(2);
  gStyle->SetFuncWidth(2);

  gStyle->SetStatFont(42);
  gStyle->SetOptFit(1111);
  gStyle->SetOptStat(0);
gStyle->SetOptStat(000000);


TFile f1("t2k_scan_1.root");
TH2D* contour1 = (TH2D*)f1.Get("contourplot");

TFile f2("t2k_scan_2.root");
TH2D* contour2 = (TH2D*)f2.Get("contourplot");


TCanvas *c7 = new TCanvas("c7","Correlation_a",820,700);

contour2 ->Draw("hist");
contour1 ->Draw("same");
contour2->SetLineColor(1);
contour1->SetLineColor(2);
contour2->SetLineWidth(3);
contour1->SetLineWidth(3);

TLegend *leg = new TLegend(0.6,.75,.9,.9);
leg->SetFillColor(0);
leg->AddEntry(contour2, "1 and 2 #sigma with SK ", " l");
leg->AddEntry(contour1, "1 and 2 #sigma without SK ", "l ");
leg->Draw();



}




