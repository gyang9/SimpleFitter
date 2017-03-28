
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
in.open("./test311/all.txt");

TH1D* curve1[5];
for(Int_t i=0;i<5;i++){ curve1[i] = new TH1D("","",101, 0.04, 4.04);}
TH1D* curve2[5];
for(Int_t i=0;i<5;i++){ curve2[i] = new TH1D("","",101, 2, 202);}

TH2D* p1 = new TH2D("","",201,1,201,201,0.02,2.02);
TH2D* pp1 = new TH2D("pp1","pp1",100,bin_x,100,bin_y);
TH2D* p1_CL = new TH2D("p1_CL","p1_CL",100,bin_x,100,bin_y);
double p_n[100]={};


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
double pre=0;

while(1){in>>s>>a>>b>>c>>d>>e>>f;
//cout<<a<<"  "<<b<<"   "<<c<<"   "<<d<<endl;
if(!in.good()){break;}
//if(b==3){break;}
int sss = (s+0.000000001-0.)/200.;
int aaa = (c+0.000000001-0.)/4.;
//if(e<0.0001){ e=0; }
p1->SetBinContent(sss,aaa,f);

if(sss == 5 )curve1[0]->SetBinContent(aaa,f);
if(sss == 10)curve1[1]->SetBinContent(aaa,f);
if(sss == 15)curve1[2]->SetBinContent(aaa,f);
if(sss == 20)curve1[3]->SetBinContent(aaa,f);
if(sss == 25)curve1[4]->SetBinContent(aaa,f);

if(aaa == 10)curve2[0]->SetBinContent(sss,f);
if(aaa == 20)curve2[1]->SetBinContent(sss,f);
if(aaa == 30)curve2[2]->SetBinContent(sss,f);
if(aaa == 40)curve2[3]->SetBinContent(sss,f);
if(aaa == 50)curve2[4]->SetBinContent(sss,f);

}
in.close();
/*
for(Int_t i=2;i<p1->GetNbinsX();i++) { 
  for (Int_t j=2;j<p1->GetNbinsY();j++){
       if(TMath::Abs(p1->GetBinContent(i,j)-p1->GetBinContent(i+1,j))>2 && TMath::Abs(p1->GetBinContent(i,j)-p1->GetBinContent(i-1,j)) >2 ) p1->SetBinContent(i,j,(p1->GetBinContent(i+1,j)+p1->GetBinContent(i-1,j))/2.);
     }
  }
for(Int_t i=2;i<p1->GetNbinsX();i++) { 
  for (Int_t j=2;j<p1->GetNbinsY();j++){
       if(TMath::Abs(p1->GetBinContent(i,j)-p1->GetBinContent(i,j+1))>2 && TMath::Abs(p1->GetBinContent(i,j)-p1->GetBinContent(i,j-1)) >2 ) p1->SetBinContent(i,j,(p1->GetBinContent(i,j+1)+p1->GetBinContent(i,j-1))/2.);
     }
  }

for(Int_t i=4;i<p1->GetNbinsX()-3;i++) { 
  for (Int_t j=4;j<p1->GetNbinsY()-3;j++){
        if(p1->GetBinContent(i,j)>20)  p1->SetBinContent(i,j, (p1->GetBinContent(i+3,j+3)+p1->GetBinContent(i-3,j-3))/2.);
     }
  }
for(Int_t i=4;i<p1->GetNbinsX()-3;i++) { 
  for (Int_t j=4;j<p1->GetNbinsY()-3;j++){
        if(p1->GetBinContent(i,j)>20)  p1->SetBinContent(i,j, (p1->GetBinContent(i+6,j+6)+p1->GetBinContent(i-6,j-6))/2.);
     }
  }
for(Int_t i=4;i<p1->GetNbinsX()-3;i++) { 
  for (Int_t j=4;j<p1->GetNbinsY()-3;j++){
        if(p1->GetBinContent(i,j)>20)  p1->SetBinContent(i,j, (p1->GetBinContent(i+7,j+7)+p1->GetBinContent(i-7,j-7))/2.);
     }
  }

for(Int_t i=3;i<p1->GetNbinsX()-1;i++) { 
  for (Int_t j=3;j<p1->GetNbinsY()-1;j++){
       if(TMath::Abs(p1->GetBinContent(i,j)-p1->GetBinContent(i+2,j))>2 && TMath::Abs(p1->GetBinContent(i,j)-p1->GetBinContent(i-2,j)) >2 ) p1->SetBinContent(i,j,(p1->GetBinContent(i+2,j)+p1->GetBinContent(i-2,j))/2.);
     }
  }
for(Int_t i=3;i<p1->GetNbinsX()-1;i++) { 
  for (Int_t j=3;j<p1->GetNbinsY()-1;j++){
       if(TMath::Abs(p1->GetBinContent(i,j)-p1->GetBinContent(i,j+2))>2 && TMath::Abs(p1->GetBinContent(i,j)-p1->GetBinContent(i,j-2)) >2 ) p1->SetBinContent(i,j,(p1->GetBinContent(i,j+2)+p1->GetBinContent(i,j-2))/2.);
     }
  }

for(Int_t i=40;i<p1->GetNbinsX()-20;i++) { 
  for (Int_t j=40;j<p1->GetNbinsY()-20;j++){
       if(TMath::Abs(p1->GetBinContent(i,j)-p1->GetBinContent(i+6,j))>2 && TMath::Abs(p1->GetBinContent(i,j)-p1->GetBinContent(i-6,j)) >2 ) p1->SetBinContent(i,j,(p1->GetBinContent(i+6,j)+p1->GetBinContent(i-6,j))/2.);
     }
  }
for(Int_t i=40;i<p1->GetNbinsX()-20;i++) { 
  for (Int_t j=40;j<p1->GetNbinsY()-20;j++){
       if(TMath::Abs(p1->GetBinContent(i,j)-p1->GetBinContent(i,j+6))>2 && TMath::Abs(p1->GetBinContent(i,j)-p1->GetBinContent(i,j-6)) >2 ) p1->SetBinContent(i,j,(p1->GetBinContent(i,j+6)+p1->GetBinContent(i,j-6))/2.);
     }
  }
*/
cout<<"done in1"<<endl;

TCanvas *c1 = new TCanvas("c1","",750,700);
c1->SetFillColor(0);
c1->SetGrid();
p1->GetXaxis()->SetTitle("Number of bins");
p1->GetYaxis()->SetTitle("Number of years");
p1->GetXaxis()->CenterTitle(); 
p1->GetYaxis()->CenterTitle(); 
p1->GetYaxis()->SetTitleOffset(1);
p1->Draw("colz");


TCanvas *c2 = new TCanvas("c2","",750,700);
c2->SetFillColor(0);
c2->SetGrid();
curve1[1]->GetXaxis()->SetTitle("Statistics");
curve1[1]->GetYaxis()->SetTitle("#Delta#chi^{2}");
curve1[1]->GetXaxis()->CenterTitle(); 
curve1[1]->GetYaxis()->CenterTitle(); 
curve1[1]->GetYaxis()->SetTitleOffset(1);
curve1[1]->Draw("hist");
for(Int_t i=2;i<5;i++){
curve1[i]->Draw("same");
}

for(Int_t i=0;i<5;i++){curve1[i]->SetLineColor(i+1);}

  TLegend* leg = new TLegend(0.6,0.68,0.9,0.9,NULL,"brNDC");

  //leg->AddEntry(curve1[0] ,"10 bins","l");
  leg->AddEntry(curve1[1] ,"20 bins","l");
  leg->AddEntry(curve1[2] ,"30 bins","l");
  leg->AddEntry(curve1[3] ,"40 bins","l");
  leg->AddEntry(curve1[4] ,"50 bins","l");

  leg->SetBorderSize(2);
  leg->SetTextFont(62);
  leg->SetFillColor(1);
  leg->SetLineColor(1);
  leg->SetShadowColor(0);
  leg->SetLineStyle(1);
  leg->SetLineWidth(2);
  leg->SetFillColor(0);
  leg->SetFillStyle(1001);
  leg->Draw();



TCanvas *c3 = new TCanvas("c3","",750,700);
c3->SetFillColor(0);
c3->SetGrid();
curve2[0]->GetXaxis()->SetTitle("Number of bins");
curve2[0]->GetYaxis()->SetTitle("Number of years");
curve2[0]->GetXaxis()->CenterTitle(); 
curve2[0]->GetYaxis()->CenterTitle(); 
curve2[0]->GetYaxis()->SetTitleOffset(1);
curve2[0]->Draw("hist");
for(Int_t i=1;i<5;i++){
curve2[i]->Draw("same");
}

for(Int_t i=0;i<5;i++){curve2[i]->SetLineColor(i+1);}

  TLegend* leg = new TLegend(0.6,0.68,0.9,0.9,NULL,"brNDC");

  leg->AddEntry(curve2[0] ,"0.8 times statistics","l");
  leg->AddEntry(curve2[1] ,"1.2 times statistics","l");
  leg->AddEntry(curve2[2] ,"1.6 times statistics","l");
  leg->AddEntry(curve2[3] ,"2.0 times statistics","l");
  leg->AddEntry(curve2[4] ,"2.4 times statistics","l");

  leg->SetBorderSize(2);
  leg->SetTextFont(62);
  leg->SetFillColor(1);
  leg->SetLineColor(1);
  leg->SetShadowColor(0);
  leg->SetLineStyle(1);
  leg->SetLineWidth(2);
  leg->SetFillColor(0);
  leg->SetFillStyle(1001);
  leg->Draw();




}
