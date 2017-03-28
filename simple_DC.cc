/*
 *  t2k simple combined fit
 *
 *  Author: Guang Yang
 */
#include "simple_t2k.hh"
#include "TMath.h"

#include "RooArgList.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include <TVectorD.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TFrame.h>

using namespace std;

  Sterile ::Sterile (const char* name) 
  : RooAbsReal(name,name)
{

// there will be: pull 0-2 s12, s23, s13, pull 3, CP delta, pull 4-6 dm2(21,32,31), pull 7-10 numu, nue Xsec and numu, nue selections

  _pulls     = new RooListProxy("_pulls","_pulls",this);
  RooRealVar* Par1 = new RooRealVar("s12","par1",TMath::ASin(TMath::Sqrt(0.85))/2.,0,100);
  RooRealVar* Par2 = new RooRealVar("s23","par2",TMath::ASin(TMath::Sqrt(0.95))/2.,0,100);
  RooRealVar* Par3 = new RooRealVar("s13","par3",0.1,0,100);
  RooRealVar* Par4 = new RooRealVar("delta","par4",-1.5,-10,10);
  RooRealVar* Par5 = new RooRealVar("dm21","par5",0.000075,-10,10);
  RooRealVar* Par6 = new RooRealVar("dm32","par6",0.00238,-10,10);
  RooRealVar* Par7 = new RooRealVar("dm31","par7",0.00244,-10,10);
  RooRealVar* Par8 = new RooRealVar("numuX","par8",1,0.,100);
  RooRealVar* Par9 = new RooRealVar("nueX","par9",1,0.,100);
  RooRealVar* Par10 = new RooRealVar("numuSel","par10",1,0.,100);
  RooRealVar* Par11 = new RooRealVar("nueSel","par11",1,0.,100);


Par1->setConstant(false);
Par2->setConstant(false);
Par3->setConstant(false);
Par4->setConstant(false);
Par5->setConstant(false);
Par6->setConstant(false);
Par7->setConstant(false);
Par8->setConstant(false);
Par9->setConstant(false);
Par10->setConstant(false);
Par11->setConstant(false);


_parlist.add(*Par1);
_parlist.add(*Par2);
_parlist.add(*Par3);
_parlist.add(*Par4);
_parlist.add(*Par5);
_parlist.add(*Par6);
_parlist.add(*Par7);
_parlist.add(*Par8);
_parlist.add(*Par9);
_parlist.add(*Par10);
_parlist.add(*Par11);
_pulls->add(_parlist);

  this->addServerList(*_pulls);

   };



  Sterile ::~Sterile ()
  {;}

TMatrixD* Sterile::prepareMatrix(TH2D* conv) const
{
Int_t nbin = conv->GetNbinsX();
TMatrixD* convMat = new TMatrixD(nbin,nbin);
Double_t sum[nbin];
    for(Int_t i=0;i<nbin;i++){
     for(Int_t j=0;j<nbin;j++){
       sum[i] += conv->GetBinContent(i,j);
                               }
                             }

    for(Int_t i=0;i<nbin;i++){
      for(Int_t j=0;j<nbin; j++){
        (*convMat)(i,j) = conv->GetBinContent(i,j)/sum[i];
                                }
                             }  

    return convMat;
}


TMatrixD* Sterile::prepareT2kCovMatrix(TMatrixD* covM_t2k, TVectorD* fVec_t2k) const
{
Int_t nbin = 30;
TMatrixD* covmatT2k = new TMatrixD(nbin,nbin);
    for(Int_t i=0;i<nbin;i++){
      for(Int_t j=0;j<nbin; j++){
        (*covmatT2k )(i,j) = (*covM_t2k)(i,j) * (*fVec_t2k)[i] * (*fVec_t2k)[j];
                                }
                             }  
    for(Int_t i=0;i<nbin;i++){
          (*covmatT2k )(i,i) += (*fVec_t2k)[i];                           
                             }
std::cout<<"t2k matrix sum "<<covmatT2k->Sum()<<std::endl;
return covmatT2k ;
}


TMatrixD* Sterile::prepareSkCovMatrix(TMatrixD* covM_sk, TVectorD* fVec_sk) const
{
Int_t nbin = 10;
TMatrixD* covmatSk = new TMatrixD(nbin,nbin);
    for(Int_t i=0;i<nbin;i++){
      for(Int_t j=0;j<nbin; j++){
        (*covmatSk )(i,j) = (*covM_sk)(i,j) * (*fVec_sk)[i] * (*fVec_sk)[j];
                                }
                             }  
return covmatSk ;
}

TMatrixD* Sterile::prepareJunoCovMatrix(Int_t nBins, TMatrixD* covM_JUNO, TVectorD* fVec_JUNO) const
{
Int_t nbin = nBins;
TMatrixD* covmatJUNO = new TMatrixD(nbin,nbin);
    for(Int_t i=0;i<nbin;i++){
      for(Int_t j=0;j<nbin; j++){
        (*covmatJUNO )(i,j) = (*covM_JUNO)(i,j) * (*fVec_JUNO)[i] * (*fVec_JUNO)[j];
                                }
                             }  
    for(Int_t i=0;i<nbin;i++){
          (*covmatJUNO )(i,i) += (37./36.)*(*fVec_JUNO)[i];                           
                             }
//std::cout<<"JUNO matrix sum "<<covmatJUNO ->Sum()<<std::endl;
return covmatJUNO ;
}

Double_t Sterile::surv_JUNO(Double_t E, RooListProxy* _pulls, Double_t L) const
{
  Double_t delta_21 = ((RooAbsReal*)_pulls->at(4))->getVal()* (L/E)*1.27 ;
  Double_t delta_32 = ((RooAbsReal*)_pulls->at(5))->getVal()* (L/E)*1.27 ;
  Double_t delta_31 = ((RooAbsReal*)_pulls->at(6))->getVal()* (L/E)*1.27 ;

  Double_t prob = 1 -  ((RooAbsReal*)_pulls->at(2))->getVal() *  TMath::Power(TMath::Sin(delta_31),2);

    return prob;       
}

Double_t Sterile::surv_JUNO(Double_t E, TVectorD* parVec, Double_t L) const
{
  Double_t delta_21 = (*parVec)[4]* (L/E)*1.27 ;
  Double_t delta_32 = (*parVec)[5]* (L/E)*1.27 ;
  Double_t delta_31 = (*parVec)[6]* (L/E)*1.27 ;

  Double_t prob = 1 - (*parVec)[2] *  TMath::Power(TMath::Sin(delta_31),2);

    return prob;       
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



Double_t Sterile ::FillEv( RooListProxy* _pulls ) const
{

   Int_t nBins = _Bins;
   Double_t junoSpec[nBins];
   Double_t junoSpecN[nBins];
  
   Double_t intercept = 3.;
   Double_t endPoint = 8.;

      //  [0 -> 750] ~ 6 years
      TF1* fnFunc = new TF1("fnFunc","[0]*exp(-0.5 * TMath::Power((x-[1])/[2],2))",1,intercept );

      fnFunc->SetParameter(0,20000);
      fnFunc->SetParameter(1,intercept );
      fnFunc->SetParameter(2,1);

      TF1* fnFunc2 = new TF1("fnFunc2","[0]*exp(-0.5 * TMath::Power((x-[1])/[2],2))",intercept ,endPoint );

      fnFunc2->SetParameter(0,20000);
      fnFunc2->SetParameter(1,intercept );
      fnFunc2->SetParameter(2,2);


//      TF1* sumFunc = new TF1("sumFunc","fnFunc + fnFunc2");

      TF1* sumFunc = new TF1("sumFunc","[0]*exp(-0.5 * TMath::Power((x-[1])/[2],2)) + [3]*exp(-0.5 * TMath::Power((x-[4])/[5],2))",1.,endPoint );

      sumFunc->SetParameter(0,4.*6000 * _time);
      sumFunc->SetParameter(1,intercept );
      sumFunc->SetParameter(2,0.65);
      sumFunc->SetParameter(3,4.*8000 * _time);
      sumFunc->SetParameter(4,intercept );
      sumFunc->SetParameter(5,2);


         for(Int_t iBin=0; iBin<nBins ; ++iBin) {
	Double_t content = sumFunc ->Integral(1+iBin*7./nBins, 1+(iBin+1)*7./nBins);
	junoSpec[iBin]=content;
	junoSpecN[iBin]=content*6.;
         }

/*
         for(Int_t iBin=0; iBin<nBins *(intercept /endPoint ); ++iBin) {
	Double_t content = fnFunc->Integral(iBin*intercept /(nBins *(intercept /endPoint )), (iBin+1)*intercept /(nBins *(intercept /endPoint )));
	junoSpec[iBin]=content;
      }
         for(Int_t iBin=0; iBin<nBins *((endPoint - intercept)/endPoint ); ++iBin) {
        Int_t inBinNum = iBin+nBins *(intercept /endPoint );
	Double_t content = fnFunc->Integral(inBinNum *(endPoint - intercept)/(nBins *((endPoint - intercept)/endPoint )), (inBinNum +1)*(endPoint - intercept)/(nBins *((endPoint - intercept)/endPoint )));
	junoSpec[inBinNum]=content;
   //     std::cout<<inBinNum<<" "<<junoSpec[inBinNum]<<std::endl;
      }


   TH1D* tranBin = new TH1D("","",nBins,1.,endPoint );
   for(Int_t i=0;i<700;i++){tranBin->SetBinContent(i+1, junoSpec[i]);}
   tranBin->Rebin(700/nBins);
   for(Int_t i=0;i<nBins;i++){junoSpec[i]=tranBin->GetBinContent(i+1);} 
*/
   TMatrixD* covM_JUNO = new TMatrixD(nBins,nBins);
   for(Int_t i=0;i<nBins;i++){ (*covM_JUNO)(i,i) = 0.007* 0.017 * 0.017; }

   for(Int_t i=0;i<nBins;i++){ (*covM_JUNO)(i,i) += 0.001 * 0.001; }

   for(Int_t i=0;i<nBins;i++){
       for(Int_t j=0;j<nBins;j++){
         (*covM_JUNO)(i,j) += 0.869* 0.007 * 0.007; 
             }
       }

   Int_t nSpec = 1;
   TH1D* allVec[nSpec];
   allVec[0] = new TH1D("","",nBins,1.,endPoint );
   TVectorD* fVec_JUNO = new TVectorD(nBins);
   TVectorD* fVecShadow_JUNO = new TVectorD(nBins);
   TVectorD* fVecShadow2_JUNO = new TVectorD(nBins);

   for(Int_t i=0;i<nBins;i++){   
         Double_t thisE = (((RooAbsReal*)_pulls->at(8))->getVal()-1.) + ((RooAbsReal*)_pulls->at(7))->getVal() * allVec[0]->GetBinCenter(i+1);
         (*fVec_JUNO)[i] = (1/6.) * junoSpec[i] * this->surv_JUNO( (thisE+0.8) , _pulls, 1050) /  this->surv_JUNO( (thisE+0.8), _pulls, 400) ;
         (*fVecShadow2_JUNO)[i] = junoSpec[i];
                        }


   for(Int_t i=0;i<nBins;i++){   
         Double_t thisE = (((RooAbsReal*)_pulls->at(8))->getVal()-1.) + ((RooAbsReal*)_pulls->at(7))->getVal() * allVec[0]->GetBinCenter(i+1);
         (*fVecShadow_JUNO)[i] = (1/6.) * junoSpecN[i] * this->surv_JUNO( (thisE+0.8) , _pulls, 400) ;
                        }
   
                           
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  TVectorD* dataIn = new TVectorD(11);
  (*dataIn )[0] = TMath::ASin(TMath::Sqrt(0.85))/2.;
  (*dataIn )[1] = TMath::ASin(TMath::Sqrt(0.95))/2.;
  (*dataIn )[2] = 0.1;
  (*dataIn )[3] = -1.5;
  (*dataIn )[4] = 0.000075;
  (*dataIn )[5] = 0.00238;
  (*dataIn )[6] = 0.00244;
  (*dataIn )[7] = 1;
  (*dataIn )[8] = 1;
  (*dataIn )[9] = 1;
  (*dataIn )[10] = 1;

 //  std::cout<<"setting up data "<<std::endl;
   TVectorD* fData_t2k = new TVectorD(30);
   TVectorD* fData_sk = new TVectorD(10);
   TVectorD* fData_JUNO = new TVectorD(nBins);

   for(Int_t i=0;i<nBins;i++){ 
         Double_t thisE2 =  allVec[0]->GetBinCenter(i+1);
        (*fData_JUNO)[i] = (1/6.) * junoSpec[i] * this->surv_JUNO( (thisE2+0.8), dataIn, 1050 ) / this->surv_JUNO( (thisE2+0.8), dataIn, 400 )  ;  
 //     std::cout<<i<<" "<<(*fData_JUNO)[i]<<" "<<(*fVec_JUNO)[i]<<std::endl;
              }


TH1D* hShadow_JUNO =  new TH1D("","",nBins,1,8);
TH1D* hData_JUNO = new TH1D("","",nBins,1,8);
for(Int_t i=0;i<nBins; i++){
hShadow_JUNO->SetBinContent(i+1,(*fVecShadow2_JUNO)(i));
hData_JUNO->SetBinContent(i+1,(*fVec_JUNO)(i));
}

  ofstream out;
  out.open(Form("bin%d.txt",nBins),ios::out);
     for(Int_t i=0;i<nBins;i++){ out<<hShadow_JUNO->GetBinContent(i+1)<<" "<<6*hData_JUNO->GetBinContent(i+1)<<endl; }


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


   TMatrixD* covMatJUNO = this->prepareJunoCovMatrix(nBins ,covM_JUNO, fVec_JUNO);

 //  fVec_JUNO->Print();
 //  std::cout<<"Ok"<<std::endl;
 //  fData_JUNO->Print();
   
   for(Int_t i=0;i<nBins; i++){
       (*fVec_JUNO)[i] -= (*fData_JUNO)[i];
         if( (*covMatJUNO)(i,i) ==0 ){(*covMatJUNO)(i,i) = 1000000000;}
     }

   covMatJUNO->Invert();

   TVectorD mulVec(*fVec_JUNO);
   mulVec *= (*covMatJUNO);

   Double_t JunoResult = TMath::Abs(mulVec*(*fVec_JUNO));
   std::cout<<"chi2 sans pull "<<JunoResult<<std::endl;

   return (Double_t) JunoResult ; 
}


Double_t Sterile ::evaluate() const
{ 

Double_t matPart = this->FillEv(_pulls);

Double_t extraPull = this -> ExtraPull (_pulls);
Double_t tot = matPart + extraPull; //If needed, add pull terms here.

return tot;

}

Double_t Sterile ::ExtraPull (RooListProxy* _pulls) const
{
Double_t pullAdd = 0;
for(Int_t i=0;i<11;i++){
 pullAdd += TMath::Power(( ((RooAbsReal*)_pulls->at(i))->getVal() - (*pullCV)[i] ),2) / TMath::Power( (*pullUnc)[i],2) ;
    }
 std::cout<<"extra pull penalty: "<<pullAdd<<std::endl;
 return pullAdd;
}


int main(int argc, char**argv){


RooFitResult* res;
Sterile * rep = new Sterile ("_rep");
  char formula[10];


rep->setAtmBaseline(0);
rep->setDensity(0);
TH1D* vecInput1 = new TH1D("","",11,0,11);
TH1D* vecInput2 = new TH1D("","",11,0,11);
//vecInput1->Print();

vecInput1 ->SetBinContent(1, TMath::ASin(TMath::Sqrt(0.85))/2.);
vecInput1 ->SetBinContent(2,TMath::ASin(TMath::Sqrt(0.95))/2.);
vecInput1 ->SetBinContent(3,0.1);
vecInput1 ->SetBinContent(4,-1.5);
vecInput1 ->SetBinContent(5,0.000075);
vecInput1 ->SetBinContent(6,0.00238);
vecInput1 ->SetBinContent(7,0.00244);
vecInput1 ->SetBinContent(8,1);
vecInput1 ->SetBinContent(9,1);
vecInput1 ->SetBinContent(10,1);
vecInput1 ->SetBinContent(11,1);

vecInput2 ->SetBinContent(1,TMath::ASin(TMath::Sqrt(0.85+0.021))/2.-TMath::ASin(TMath::Sqrt(0.85))/2.);
vecInput2 ->SetBinContent(2,2);
vecInput2 ->SetBinContent(3,100);
vecInput2 ->SetBinContent(4,2);
vecInput2 ->SetBinContent(5,0.0000018);
vecInput2 ->SetBinContent(6,0.0001);
vecInput2 ->SetBinContent(7,0.00009);
vecInput2 ->SetBinContent(8,0.009);
vecInput2 ->SetBinContent(9,0.03);
vecInput2 ->SetBinContent(10,0.04); 
vecInput2 ->SetBinContent(11,0.027);

rep->setPull(vecInput1); 
rep->setPullUnc(vecInput2);
rep->addSK(true);

rep->setNBins(200);
rep->setTime(6);
rep->FillEv(rep->getPullList());


    const int npts = 100;
    double bin_x[npts]={};
    double low_x = 1, high_x = 10;
    double step_x = (high_x-low_x)/(double)npts;  
    for(int i=0;i<npts ;i++){bin_x[i]= low_x + i * step_x;}

  int aaa = (atof(argv[5])+0.000001)*100;
  ofstream out2;
  out2.open(Form("s%d_a%d_JUNO2D.txt",atoi(argv[2]),aaa ),ios::app);
  Int_t binSetup = atoi(argv[2]);

RooArgList list("list");
list.add(*rep);
  sprintf(formula,"%s","@0");
RooFormulaVar* fcn = new RooFormulaVar("fit","fit",formula,list);

//rep->getParVar(3)->setVal(2.5);
//rep->getParVar(3)->setConstant(true);
//rep->getParVar(0)->setConstant(true);
//rep->getParVar(1)->setConstant(true);
//rep->getParVar(2)->setConstant(true);
//rep->getParVar(3)->setConstant(true);
//rep->getParVar(7)->setConstant(true);
rep->getParVar(8)->setConstant(true);
//rep->getParVar(9)->setConstant(true);

//    for(Int_t nn= 0; nn< npts ; ++nn) {

rep->setNBins(binSetup );
rep->setTime(atof(argv[5]));

RooMinuit m(*fcn);
m.setStrategy(2);
 Double_t callsEDM[2] = {10500., 1.e-6};
 Int_t irf = 0;

 gMinuit->mnexcm("MIGRAD",callsEDM,2,irf);
m.migrad();
//m.hesse();
m.minos(); 
res = m.save();
double bestFit = res->minNll(); 

std::cout<<"result list: "<<std::endl;
std::cout<<"chi2: "<<bestFit <<std::endl;

  double bb =  rep->getParVar(2)->getAsymErrorLo();
  double dd =  rep->getParVar(2)->getAsymErrorHi();

cout<<"errors are "<<bb<<" "<<dd<<endl;

for(Int_t i=0;i<11;i++){std::cout<<" "<<rep->getPar(i)<<std::endl;}

  int aa = (atof(argv[5])+0.000001)*100;

  out2<<atoi(argv[2])<<" "<<aa<<" "<<dd<<endl; 
  cout<<atoi(argv[2])<<" "<<aa<<" "<<dd<<endl; 

//  }
}



Double_t Sterile ::getPar(int i) {
(((RooAbsReal*)_pulls->at(i))->getVal());
}

RooRealVar* Sterile ::getParVar(int i) {
return ((RooRealVar*)_pulls->at(i));
}

TVectorD* Sterile::getT2kVec() {
return fVec_t2k;
}

TVectorD* Sterile::getJunoVec() {
return fVecShadow_JUNO;
}

TVectorD* Sterile::getJunoData() {
return fData_JUNO;
}

TH1D* Sterile::getJunoHist() {
return hShadow_JUNO;
}

TH1D* Sterile::getJunoDataHist() {
return hData_JUNO;
}

TVectorD* Sterile::getSkVec() {
return fVec_sk;
}

void Sterile :: setSyst(Double_t syst){
_syst = syst;
}

void Sterile :: setdm2CV(Double_t dm2CV){
_dm2CV = dm2CV;
}

void Sterile :: setdm2Unc(Double_t dm2Unc){
_dm2Unc = dm2Unc;
}

void Sterile :: addSK(Bool_t wSK){
withSK = wSK;
}

void Sterile :: setAtmBaseline(Double_t AtmBaseline){
_AtmBaseline = AtmBaseline;
}

void Sterile :: setDensity(Double_t Density){
_Density = Density;
}

void Sterile :: setNBins(Double_t Bins){
_Bins= Bins;
}

void Sterile :: setTime(Double_t time){
_time= time;
}

void Sterile :: setPull(TH1D* pullvecCV){
pullCV = new TVectorD(11);
for(Int_t i=0;i<11;i++){
(*pullCV)[i] =  pullvecCV->GetBinContent(i+1);
    }
}

void Sterile :: setPullUnc(TH1D* pullvecUnc){
pullUnc = new TVectorD(11);
for(Int_t i=0;i<11;i++){
(*pullUnc)[i] = pullvecUnc->GetBinContent(i+1);
    }
}

Double_t Sterile::getPullUnc(Int_t pN){
return (*pullUnc)[pN];
}

void Sterile::DataSwitch(Bool_t dataSwitch) const
{
Bool_t _dataSwitch = dataSwitch;
}

Bool_t Sterile::getDataSwitch() const
{
return _dataSwitch;
}

RooListProxy* Sterile::getPullList() const
{
return _pulls;
}

