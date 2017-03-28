{
//this macro reads asci data file and plots data

#include<iostream.h>
#include<fstream.h>
#include<iomanip.h>
#include<string>
#include<vector>
#include<strstream>
#include<ctime>
#include<math.h>


gROOT->Reset();
gStyle->SetOptFit(1111);
gStyle->SetOptStat(1111);


 double s22theta = 0.1; //this is sin^2 (2*theta), where theta is the mixing angle 
 double dm2=0.4; //this is delta(m^2) in [eV] 



//muon to electron neutrino oscillation probability, P_me= s22theta x sin^2(1.27 dm2 L/Enu_MeV)
//we may calculate oscillation probability in several different cases

//simple case 1: fix the energy of neutrino to value L_f, then calculate oscillation probability as a function of distance
//L[n_data_1]

 int n_data_1=1000; //number of points
 double E_f=700.0; //energy in [MeV]
 double L[n_data_1] = {}; //distance will be calculated in [m]
 double P_me_1[n_data_1] = {}; //oscillation probability will be calculated

 //simple case 2: fox the distance L between neutrino source and the detector, then calculate the oscillation probability
 //as a function of neutrino energy E[n_data_2] 
 
 int n_data_2 = 1000;//number of points
 double L_f=600.0; //distance in [m]
 double E[n_data_2] = {}; //distance will be calculated in [m]
 double P_me_2[n_data_2] = {}; //oscillation probability will be calculated



//calculate probability as a function of distance for fixed energy
for(int i=0; i<n_data_1; i++) {
  L[i] = i*1; //each point is 1 m apart
  if(E_f!=0) P_me_1[i] = s22theta*( sin(1.27*dm2*L[i]/E_f) * sin(1.27*dm2*L[i]/E_f) ); 
  //cout<< "L = " << L[i] << "; " << "P_me_1 = " << P_me_1[i] << " ; " << endl; 
  
 }


//calculate probability at fixed disyatne as a function of energy
for(int i=0; i<n_data_2; i++) {
  E[i] = i*1; //each point is 1 MeV apart
  if(E[i]!=0) P_me_2[i] = s22theta*( sin(1.27*dm2*L_f/E[i]) * sin(1.27*dm2*L_f/E[i]) ); 
  //cout<< "E = " << E[i] << "; " << "P_me_2 = " << P_me_2[i] << " ; " << endl; 
 
 }

/////////////////////////////
//make graphs of calculated points

TCanvas *c1 = new TCanvas("c1", "Survival Probability",200, 30, 700, 500);
//c1->SetFillColor(20);
c1->SetGrid();
c1->GetFrame()->SetFillColor(21);
c1->GetFrame()->SetBorderSize(12);

TGraph *gr_P_me_1 = new TGraph(n_data_1, L, P_me_1); 
gr_P_me_1->SetMarkerColor(2);
gr_P_me_1->SetLineColor(2);
gr_P_me_1->SetMarkerStyle(20);
gr_P_me_1->Draw("AC");
gr_P_me_1->GetXaxis() ->SetTitle("Travel Distance [m]");
gr_P_me_1->GetYaxis() ->SetTitle("Oscillation Probability");

/*
TLatex *latex1 = new TLatex();
latex1->SetTextColor(1);
latex1->SetTextSize(0.033);
latex1->DrawLatex(6. , 0.9 ,"P(E_{#nu})= 1 - sin^{2}2#theta_{12}sin^{2}(#Delta m_{12}^{2} L/E_{#nu})");
latex1->DrawLatex(6. , 0.8 ,"sin^{2}2#theta_{12} = 1.0");
latex1->DrawLatex(6. , 0.7 ,"#Delta m_{12}^{2} = 6.9x10^{-5} eV^{2}");
latex1->DrawLatex(6. , 0.6 ,"L = 160459.9 m");
*/
c1->Update();



TCanvas *c2 = new TCanvas("c2", "Survival Probability 2",200, 30, 700, 500);
//c2->SetFillColor(20);
c2->SetGrid();
c2->GetFrame()->SetFillColor(21);
c2->GetFrame()->SetBorderSize(12);

TGraph *gr_P_me_2 = new TGraph(n_data_2, E, P_me_2); 
gr_P_me_2->SetMarkerColor(4);
gr_P_me_2->SetLineColor(4);
gr_P_me_2->SetMarkerStyle(20);
gr_P_me_2->Draw("AC");
gr_P_me_2->GetXaxis() ->SetTitle("Energy [MeV]");
gr_P_me_2->GetYaxis() ->SetTitle("Oscillation Probability");
c2->Update();

}
