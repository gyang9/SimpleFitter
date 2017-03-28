
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
  gStyle->SetTitleOffset(0.5,"y");

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
int c=0;
int d=0;
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
for(int i=49;i>=0;i--){bin_x[49-i]=TMath::Power(10.,(-5.*i*2/100.));}
for(int i=0;i<50;i++){bin_y[i]=TMath::Power(10,(-2 + 5.*i*2/100.));}
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
double minn=800;
in.open("./all.txt");

//TH2D* p1_CL2 = new TH2D("","",50,0.05,0.15,50,0.0015,0.0035);
//TH2D* p1_CL3 = new TH2D("","",50,0.05,0.15,50,0.0005,0.0045);
//TH2D* p1_CL3 = new TH2D("","",50,0.05,0.15,50,0.0005,0.0045);
//TH2D* p1_CL = new TH2D("","",31,(-0.7+1)*0.97+0.604+0.0701,(0.5+1)*0.97+0.604+0.0701,31,(-0.5+1)*0.95+4.334+1.55,(2.1+1)*0.95+4.334+1.55);

TH2D* p1_CL = new TH2D("p1_CL","",49,bin_x,49,bin_y);
TH2D* p1_CL2 = new TH2D("p1_CL2","sterile sensitivity",49,bin_x,49,bin_y);
TH2D* p1_CL3 = new TH2D("p1_CL3","",49,bin_x,49,bin_y);
//TH1D* cp1 = new TH1D("cp1","cp1",200,-100,100);
//TH1D* cp2 = new TH1D("cp2","cp2",200,-100,100);
//TH2D* p_n = new TH2D("p_n","p_n",100,bin_x,100,bin_v);
//TH2D* p_n2 = new TH2D("p_n2","p_n2",100,bin_y,100,bin_v);
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
ree_b[i]=100;
}
for(int i=0;i<20000;i++){
re_s[i]=1000;
re_a[i]=1000;
re_b[i]=1000;
re_c[i]=1000;
}
double min_s,min_a,min_b,min_c,minn;

double minab[51][51]={};
double minbc[51][51]={};
double minac[51][51]={};
for(int mm=0;mm<51;mm++){
for(int nn=0;nn<51;nn++){
minbc[mm][nn]=1000;
minac[mm][nn]=1000;
minab[mm][nn]=1000;
}
}
double marker=1000;

while(1){in>>c>>a>>b>>s>>e;
//cout<<a<<"  "<<b<<"   "<<c<<"   "<<d<<endl;
if(!in.good())
//if(b==3)
{break;}
if(e<minn){minn=e;}
if (e<re_s[s]){re_s[s] = e; }
if (e<re_b[b]){re_b[b] = e; }

int bhere=(b+0.00000001)/1.;
int chere=(s+0.00000001)/1.;
minbc[bhere][chere]=e;
//if(e<minbc[bhere][chere]){minbc[bhere][chere]=e;}

}
in.close();



//b,c
for(int mm=0;mm<50;mm++){
for(int nn=0;nn<50;nn++){
p1_CL->SetBinContent(mm,nn,minbc[mm][nn]);
p1_CL2->SetBinContent(mm,nn,minbc[mm][nn]);
p1_CL3->SetBinContent(mm,nn,minbc[mm][nn]);}
}

cout<<"done in1 "<<minn<<endl;
/*
double xx1[51],xx2[51],xx3[51],xx4[51];
for(int i=0;i<31;i++){xx1[i]=(i*0.004+0.03)/1.; }
for(int i=0;i<31;i++){xx2[i]=((i*4)/100.- 0.7 +1)*0.97+0.0701+0.604;}
for(int i=0;i<31;i++){xx3[i]=((i*8)/100.-0.5 +1 )*0.95+4.334+1.55;}


for(int i=0;i<31;i++){
ree_b[i]= re_b[i*8];
}
for(int i=0;i<31;i++){
ree_a[i]= re_a[i*4];
}
for(int i=0;i<31;i++){
ree_s[i]= re_s[i*4+30];
}


for(int i=1;i<50; i++){if(ree_s[i]<=ree_s[i+1] && ree_s[i]<=ree_s[i-1]){min_s = i; minn = ree_s[i];} }
for(int i=1;i<50; i++){if(ree_a[i]<ree_a[i+1]&& ree_a[i]<ree_a[i-1]){min_a = i;} 
if(ree_b[i]<ree_b[i+1]&& ree_b[i]<ree_b[i-1]){min_b = i;}  }
*/


//   double conto3[3] = {0.1+1.,0.1+4.,0.1+9.};
   double conto3[1]={minn+2.3};
   double conto4[1]={minn+4.6};
   double conto5[1]={minn+9.2};

   p1_CL2->SetContour(1,conto3);
   p1_CL2->SetLineWidth(4);
   p1_CL2->SetLineColor(1);

   p1_CL->SetContour(1,conto4);
   p1_CL->SetLineWidth(4);
   p1_CL->SetLineColor(2);

   p1_CL3->SetContour(1,conto5);
   p1_CL3->SetLineWidth(4);
   p1_CL3->SetLineColor(4);

TCanvas *c7 = new TCanvas("c7","s2t14 vs.dms14",700,700);

c7->SetLogx();
c7->SetLogy();
p1_CL2->GetXaxis()->SetTitle("sin^{2}2#theta_{14}");
p1_CL2->GetYaxis()->SetTitle("#Delta m^{2}_{41}");
p1_CL2->GetYaxis()->SetTitleOffset(1.4);
//p1_CL3->Draw("col ");
p1_CL2->Draw("cont3");
p1_CL->Draw("cont3 same");
p1_CL3->Draw("cont3 same");
//p3_CL->Draw("cont3 same");

TLegend *leg = new TLegend(0.7,.7,.9,.9);
leg->SetFillColor(0);
leg->AddEntry(p1_CL2, "1 #sigma ", " l");
leg->AddEntry(p1_CL, "2 #sigma ", "l ");
leg->AddEntry(p1_CL3, "3 #sigma ", "l ");
leg->Draw();
c7->SaveAs("contour3.png");

}




