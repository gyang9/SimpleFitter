
{

  gROOT->SetStyle("Plain");
  gStyle->SetTitleBorderSize(0);
  gStyle->SetOptStat("");

  gStyle->SetLabelFont(102,"");
  gStyle->SetLabelSize(0.06,"");
  gStyle->SetLabelFont(102,"xyz");
  gStyle->SetLabelSize(0.04,"xyz");
  gStyle->SetLabelOffset(0.001,"x");
  gStyle->SetLabelOffset(0.01,"y");

  gStyle->SetTitleFont(42,"xyz");
  gStyle->SetTitleFontSize(0.06);
  gStyle->SetTitleFont(42);
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetTitleOffset(1.05,"x");
  gStyle->SetTitleOffset(1.05,"y");

  gStyle->SetStripDecimals(kFALSE);
  
  gStyle->SetPadLeftMargin(0.15);
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


double s=0;
double a=0;
double b=0; 
double c=0;
double d=0;
double e=0;
double f=0;
int ccon=0;
double cchh[200]={};
for(int i=0;i<200;i++){cchh[i]=1000.;}
double icounter=0,icounter2=0;

double bin_x[101]={};
double bin_y[101]={};
double bin_v[101]={};
for(int i=100;i>=0;i--){bin_x[100-i]=TMath::Power(10.,(-3.*i/100.));}
for(int i=0;i<=100;i++){bin_y[i]=TMath::Power(10.,(-4 + 5.*i/100.));}
for(int i=0;i<=100;i++){bin_v[i]=TMath::Power(10.,(-4 + 4.*i/100.));}

double h0[10000][10000]={};
double h1[10000][10000]={};
double h0_CL[10000][10000]={};
double referee[10000]={};
double res[10000]={};
double crit1[10000]={};

ifstream in;
ifstream in2;
ifstream in3;
double min=800;
in.open("./test308/all.txt");

TH2D* p1 = new TH2D("","",50,2,102,100,0.05,5.05);
TH2D* p2 = new TH2D("","",50,2,102,100,0.05,5.05);
TH2D* pp1 = new TH2D("pp1","pp1",100,bin_x,100,bin_y);
TH2D* p1_CL = new TH2D("p1_CL","p1_CL",100,bin_x,100,bin_y);
double p_n[100]={};


TH1D* cp1 = new TH1D("cp1","cp1",100,0,1);
TH1D* cp2 = new TH1D("cp2","cp2",100,0,1);
double re_s[20000];
double re_a[20000];
double re_b[20000];
double re_c[20000];
double ree_s[51];
double ree_a[51];
double ree_b[51];
double ree_c[51];
for(int i=0;i<51;i++){
ree_s[i]=100;
ree_b[i]=100;
}
for(int i=0;i<20000;i++){
re_s[i]=1000;
re_a[i]=1000;
re_b[i]=1000;
re_c[i]=1000;
}
double min_s,min_a,min_b,min_c,minn;
double pre=0;

double recording[4]={};
double line1[100]={};
double line2[100]={};
double line3[100]={};
double line4[100]={};

double lineS1[100]={};
double lineS2[100]={};
double lineS3[100]={};
double lineS4[100]={};

int counter=0;
while(1){in>>s>>a>>b>>c>>d>>e>>f;
//cout<<a<<"  "<<b<<"   "<<c<<"   "<<f<<endl;
if(!in.good())
//if(b==3)
{break;}

//if(counter > 10) break;
//counter++;
int sss = (s+0.000000001)/10.;
int aaa = (a+0.000000001)/10.;
int bbb = (b+0.000000001)/10.;
int ccc = (c+0.000000001)/1.;
double ddd = f;

if(aaa==0 && bbb==0 && ccc==0) {line1[sss]=ddd; }
if(sss==0 && bbb==0 && ccc==0) {line2[aaa]=ddd; }
if(sss==0 && aaa==0 && ccc==0) {line3[bbb]=ddd; }
if(sss==0 && aaa==0 && bbb==0) {line4[ccc]=ddd; }
//cout<<"here2"<<endl;
if(ddd > lineS1[sss] && ddd<200 ) { lineS1[sss]=ddd; }
if(ddd > lineS2[aaa] && ddd<200) { lineS2[aaa]=ddd; }
if(ddd > lineS3[bbb] && ddd<200) { lineS3[bbb]=ddd; }
if(ddd > lineS4[ccc] && ddd<200) { lineS4[ccc]=ddd; }
//cout<<"here"<<endl;
}
in.close();

cout<<"let's see.."<<endl;
for(Int_t i=0;i<10;i++){cout<<line1[i]<<" "<<lineS1[i]<<" "<<line2[i]<<" "<<lineS2[i]<<endl;}

cout<<"done in1"<<endl;

double bin_x[10],bin_xx[100];
for(Int_t i=0;i<10;i++){bin_x[i]=i*0.1;}
for(Int_t i=0;i<100;i++){bin_xx[i]=i*0.01;}

TCanvas *c1 = new TCanvas("c1","lines");
c1->SetFillColor(0);
c1->SetGrid();
c1->GetFrame()->SetFillColor(21);
c1->GetFrame()->SetBorderSize(12);
TGraphErrors *gr20 = new TGraphErrors(10,bin_x,line1,0,0);
TGraphErrors *gr21 = new TGraphErrors(10,bin_x,line2,0,0);
TGraphErrors *gr22 = new TGraphErrors(10,bin_x,line3,0,0);
TGraphErrors *gr23 = new TGraphErrors(100,bin_xx,line4,0,0);

gr20->GetXaxis()->SetTitle("correlation");
gr20->GetYaxis()->SetTitle("#Delta#chi^{2}");
gr20->SetTitle("correlation test");
gr20->SetMarkerColor(1);
gr21->SetMarkerColor(2);
gr22->SetMarkerColor(3);
gr23->SetMarkerColor(4);
gr20->SetLineColor(1);
gr21->SetLineColor(2);
gr22->SetLineColor(3);
gr23->SetLineColor(4);
gr20->SetLineWidth(4);
gr21->SetLineWidth(4);
gr22->SetLineWidth(4);
gr23->SetLineWidth(4);
gr20->SetMarkerStyle(20);
gr21->SetMarkerStyle(20);
gr22->SetMarkerStyle(20);
gr23->SetMarkerStyle(20);
gr20->Draw("AL");
gr21->Draw("same");
gr22->Draw("same");
gr23->Draw("same");

TLegend *leg = new TLegend(0.6,.7,.9,.9);
leg->SetFillColor(0);
leg->AddEntry(gr20, "FD bin-bin", " l");
leg->AddEntry(gr21, "FD/ND channel-channel", "l ");
leg->AddEntry(gr22, "ND bin-bin", " l");
leg->AddEntry(gr23, "FD-ND", "l ");
leg->Draw();


TCanvas *c2 = new TCanvas("c2","lines2");
c2->SetFillColor(0);
c2->SetGrid();
c2->GetFrame()->SetFillColor(21);
c2->GetFrame()->SetBorderSize(12);
TGraphErrors *gr20 = new TGraphErrors(10,bin_x,lineS1,0,0);
TGraphErrors *gr21 = new TGraphErrors(10,bin_x,lineS2,0,0);
TGraphErrors *gr22 = new TGraphErrors(10,bin_x,lineS3,0,0);
TGraphErrors *gr23 = new TGraphErrors(100,bin_xx,lineS4,0,0);

gr20->GetXaxis()->SetTitle("correlation");
gr20->GetYaxis()->SetTitle("#Delta#chi^{2}");
gr20->SetTitle("correlation variation");
gr20->SetMarkerColor(1);
gr21->SetMarkerColor(2);
gr22->SetMarkerColor(3);
gr23->SetMarkerColor(4);
gr20->SetMarkerStyle(20);
gr20->Draw("APL");
gr21->Draw("same");
gr22->Draw("same");
gr23->Draw("same");






}
