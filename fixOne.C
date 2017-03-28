{
   gStyle->SetOptStat(0);

   int NBins1 = 200;  double bin_min1 = 0. ; double bin_max1 = 850;  
   int NBins2 = 200;  double bin_min2 = 0.e3; double bin_max2 = 3e3; 


   TH1F* h1 = new TH1F("fixE","fixE",NBins1, bin_min1, bin_max1); 
   TH1F* h2 = new TH1F("fixL","fixL",NBins2, bin_min2, bin_max2); 
    
   //global setup for sin2t
   Double_t sin2t = 0.1;

   // i is baseline here,  energy in GeV
   double dm2 = 20.; double Energy = 700;
   for(int i=0;i<NBins1;i++){ 
   Double_t value =  sin2t * TMath::Power(TMath::Sin(dm2 * (bin_min1+(double)(i/(NBins1+0.00001))*(bin_max1-bin_min1)) /(4*Energy)),2);
   h1->SetBinContent(i,value);}

   // i is energy here , baseline in m
   dm2 = 20.; double baseline = 600;
   for(int i=0;i<NBins2;i++){ 
   Double_t value =  sin2t * TMath::Power(TMath::Sin(dm2 * baseline /(4* (bin_min2+(double)(i/(NBins2+0.000001))*(bin_max2-bin_min2)) )),2);
   h2->SetBinContent(i,value);}


   //Draw a straight line(nonoscillated)..
   TF1* non1 = new TF1("non1","[0]+[1]*x",bin_min1,bin_max1);
   non1 ->SetParameter(0,1);
   non1 ->SetParameter(1,0);

   TF1* non2 = new TF1("non2","[0]+[1]*x",bin_min2,bin_max2);
   non2 ->SetParameter(0,1);
   non2 ->SetParameter(1,0);

/////////////////////////////////////////////////////Draw a plot

TCanvas *c1 = new TCanvas();
c1->SetFillColor(0);
c1->SetGrid();
c1->GetFrame()->SetFillColor(21);
c1->GetFrame()->SetBorderSize(12);
h1->GetXaxis()->SetTitle("baseline [m]");
h1->GetYaxis()->SetTitle(" Probability");
h1->SetLineColor(4);
non1->SetLineColor(1);

h1->SetMarkerSize(0.6);
h1->SetMarkerStyle(20);
h1->Draw("hist");
non1->Draw("same");

TLegend *leg = new TLegend(0.2,.4,.6,.6);
leg->AddEntry(h1, "Oscillated ", "l");
leg->AddEntry(non1, "Non oscillated ", "l");
//leg->Draw();

///////////////////////////////////////////////////////Draw another plot

TCanvas *c2 = new TCanvas();
c2->SetFillColor(0);
c2->SetGrid();
c2->GetFrame()->SetFillColor(21);
c2->GetFrame()->SetBorderSize(12);
h2->GetXaxis()->SetTitle("Energy [GeV]");
h2->GetYaxis()->SetTitle(" Probability");
h2->SetLineColor(4);
non2->SetLineColor(1);

h2->SetMarkerSize(0.6);
h2->SetMarkerStyle(20);
h2->Draw("hist");
non2->Draw("same");

TLegend *leg = new TLegend(0.2,.4,.6,.6);
leg->AddEntry(h2, "Oscillated ", "l");
leg->AddEntry(non2, "Non oscillated ", "l");
//leg->Draw();












}
