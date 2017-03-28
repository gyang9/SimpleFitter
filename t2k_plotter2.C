
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


int s=0;
int a=0;
int b=0; 
double c=0;
double d=0;
double e=0;
int ccon=0;
double cchh[200]={};
for(int i=0;i<200;i++){cchh[i]=1000.;}
double icounter=0,icounter2=0;

double bin_x[50]={};
double bin_y[50]={};
double bin_v[101]={};
//for(int i=100;i>=0;i--){bin_x[100-i]=TMath::Power(10.,(-3.*i/100.));}
//for(int i=0;i<=100;i++){bin_y[i]=TMath::Power(10.,(-4 + 5.*i/100.));}
//for(int i=0;i<=100;i++){bin_v[i]=TMath::Power(10.,(-4 + 4.*i/100.));}
for(int i=49;i>=0;i--){bin_x[49-i]=TMath::Power(10.,(-3.*i*2/100.));}
for(int i=0;i<50;i++){bin_y[i]=TMath::Power(10,(-4 + 5.*i*2/100.));}
for(int i=0;i<50;i++){cout<<bin_x[i]<<endl; if(i>0&&bin_x[i]<bin_x[i-1]){cout<<"hey"<<endl;} }
for(int i=0;i<50;i++){cout<<bin_y[i]<<endl; if(i>0&&bin_y[i]<bin_y[i-1]){cout<<"hey"<<endl;}}

double h0[10000][10000]={};
double h1[10000][10000]={};
double h0_CL[10000][10000]={};
double referee[10000]={};
double res[10000]={};
double crit1[10000]={};

ifstream in;
ifstream in2;
ifstream in3;
ifstream in4;
ifstream in5;
ifstream in6;
ifstream in7;
ifstream in8;
ifstream in9;
double min=800;
in.open("./t2k1D_delta.txt");
in2.open("./t2kSk1D_delta.txt");
in4.open("./t2k2D_delta.txt");
in5.open("./t2kSk2D_delta.txt");

double p_n[100]={};

cout<<"done in1"<<endl;
TH1D* cp1 = new TH1D("cp1","cp1",100,-10,10);
TH1D* cp2 = new TH1D("cp2","cp2",100,-10,10);
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
ree_a[i]=100;
ree_b[i]=100;
ree_c[i]=100;
}
for(int i=0;i<20000;i++){
re_s[i]=1000;
re_a[i]=1000;
re_b[i]=1000;
re_c[i]=1000;
}
double min_s,min_a,min_b,min_c,minn;

double mind1[101]={};
double mind2[101]={};
double mind3[101]={};
double mind4[101]={};
for(int mm=0;mm<51;mm++){
mind1[mm]=1000;
mind2[mm]=1000;
mind3[mm]=1000;
mind3[mm]=1000;
}

double min1[101][101]={};
double min2[101][101]={};
double min3[101][101]={};
double min4[101][101]={};
for(int mm=0;mm<21;mm++){
for(int nn=0;nn<21;nn++){
min1[mm][nn]=1000;
min2[mm][nn]=1000;
min3[mm][nn]=1000;
min4[mm][nn]=1000;
}
}
double marker=1000;
double minn1=1000;
double minn2=1000;
double minn3=1000;
double minn4=1000;

TH2D* p1_CL = new TH2D("p1_CL","",50,0.0,0.3,   300, -3,3);
TH2D* p1_CL2 = new TH2D("p1_CL2","",50,0.0,0.3,   300, -3,3);

while(1){in>>s>>d>>e;
if(!in.good())
{break;}
s = (s+0.000001);
re_s[s] = e; 
}
in.close();

while(1){in2>>s>>d>>e;
if(!in2.good())
{break;}
s = (s+0.000001);
re_a[s] = e; 
}
in2.close();



while(1){in4>>s>>c>>d>>e;
if(!in4.good())
{break;}
if (e<minn1 ){minn1=e;}
int shere=((c+0.00000001) - 0 )/(0.3/50.);
int ahere=((d+0.00000001) + 3)/(6./300.);
min1[shere][ahere]=e;
}
in4.close();


while(1){in5>>s>>c>>d>>e;
if(!in5.good())
{break;}
if (e<minn2 ){minn2=e;}
int shere=((c+0.00000001) - 0 )/(0.3/50.);
int ahere=((d+0.00000001) + 3)/(6./300.);
min2[shere][ahere]=e;
}
in5.close();




double xx[50]={};
double xx1[50]={};
for(int i=0;i<50;i++){xx[i] = -3 + i*(6/50.); xx1[i] = -3 + i*(6/50.); }


TCanvas *c1 = new TCanvas("c1","sensitivity ");
c1->SetFillColor(0);
c1->SetGrid();
c1->GetFrame()->SetFillColor(21);
c1->GetFrame()->SetBorderSize(12);
 TGraphErrors *gr22 = new TGraphErrors(50,xx,re_s,0,0);
 TGraphErrors *gr23 = new TGraphErrors(50,xx,re_a,0,0);

gr22->SetLineColor(1);
gr22->SetLineStyle(1);
gr22->SetLineWidth(4);
gr23->SetLineColor(2);
gr23->SetLineStyle(1);
gr23->SetLineWidth(4);

TMultiGraph *mg = new TMultiGraph();
    mg->Add(gr22);
     mg->Add(gr23);
 // mg->SetMinimum(0.02);
 // mg->SetMaximum(0.04);
mg->Draw("ALP");
mg->GetXaxis()->SetTitle("#delta_{CP}");
mg->GetYaxis()->SetTitle("#Delta #chi^{2}");
//c9->BuildLegend();

TLegend *leg = new TLegend(0.5,0.7,.9,.9);
  leg->SetBorderSize(2);
  leg->SetTextFont(62);
  leg->SetLineColor(1);
  leg->SetShadowColor(0);
  leg->SetLineStyle(1);
  leg->SetLineWidth(2);
  leg->SetFillColor(0);
  leg->SetFillStyle(1001);

leg->AddEntry(gr23, "with SK", "l");
leg->AddEntry(gr22, "without SK", "l");

leg->Draw();

/////////////////////////////////////////////////////////////////////////////

//b,c
for(int mm=0;mm<50;mm++){
for(int nn=0;nn<300;nn++){
p1_CL->SetBinContent(mm,nn,min1[mm][nn]); //cout<<minbc[mm][nn]<<endl;
p1_CL2->SetBinContent(mm,nn,min2[mm][nn]);
}
}


cout<<"done in1 "<<minn1<<" "<<minn2<<endl;

   double conto3[1]={minn1+20};
   double conto4[1]={minn2+20};
   double conto5[1]={minn1+9.};

   p1_CL2->SetContour(1,conto4);
   p1_CL2->SetLineWidth(5);
   p1_CL2->SetLineColor(1);

   p1_CL->SetContour(1,conto3);
   p1_CL->SetLineWidth(5);
   p1_CL->SetLineColor(2);



TCanvas *c7 = new TCanvas("c7","Correlation_a",820,700);

p1_CL2->GetXaxis()->SetTitle("#theta_{13}");// sin^{2}(2 #theta_{13})");
p1_CL2->GetYaxis()->SetTitle("#delta_{CP} ");
p1_CL2->GetYaxis()->SetTitleOffset(1.4);
p1_CL2->Draw("cont3");
p1_CL->Draw("cont3 same");

TLegend *leg = new TLegend(0.6,.75,.9,.9);
leg->SetFillColor(0);
leg->AddEntry(p1_CL2, "90/% with SK ", " l");
leg->AddEntry(p1_CL, "90/% without SK ", "l ");
leg->Draw();

}




