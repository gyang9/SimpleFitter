{

  gROOT->SetStyle("Plain");
  gStyle->SetTitleBorderSize(0);
  gStyle->SetOptStat("");

  gStyle->SetLabelFont(102,"");
  gStyle->SetLabelSize(0.06,"");
  gStyle->SetLabelFont(102,"xyz");
  gStyle->SetLabelSize(0.06,"xyz");
  gStyle->SetLabelOffset(0.001,"x");
  gStyle->SetLabelOffset(0.01,"y");

  gStyle->SetTitleFont(42,"xyz");
  gStyle->SetTitleFontSize(0.06);
  gStyle->SetTitleFont(42);
  gStyle->SetTitleSize(0.06,"xyz");
  gStyle->SetTitleOffset(1.05,"x");
  gStyle->SetTitleOffset(1.35,"y");

  gStyle->SetStripDecimals(kFALSE);
  
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);

  gStyle->SetStatW(0.35);
  gStyle->SetStatH(0.25);

  gStyle->SetPadTickX(kTRUE);
  gStyle->SetPadTickY(kTRUE);

  gStyle->SetPalette(1);
  gStyle->SetNumberContours(99);

  gStyle->SetHistLineWidth(3);
  gStyle->SetFrameLineWidth(3);
  gStyle->SetFuncWidth(3);

  gStyle->SetStatFont(42);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);


  Int_t nBins = 100;

  ifstream in;
  in.open(Form("bin%d.txt",nBins));

TH1D* MCDraw =  new TH1D("","",nBins,0,2);
TH1D* dataDraw = new TH1D("","",nBins,0,2);
int i=0;
double s,a;
while(1){in>>s>>a;
if(!in.good())
{break;}
MCDraw->SetBinContent(i+1,s);
dataDraw->SetBinContent(i+1,a);
i++;
}
in.close();

MCDraw->Rebin(1);
dataDraw->Rebin(1);

TCanvas *c1 = new TCanvas();
c1->SetFillColor(0);
c1->SetGrid();
MCDraw->SetLineColor(2);
dataDraw->SetLineColor(4);
MCDraw->SetLineWidth(2.5);
dataDraw->SetLineWidth(2.5);
MCDraw->Draw("e");
dataDraw->Draw("same e");
TLegend *leg = new TLegend(0.68,.68,.9,.9);
leg->SetFillColor(0);
leg->AddEntry(dataDraw, "Data", " l");
leg->AddEntry(MCDraw, "MC", "l ");
leg->Draw();






}