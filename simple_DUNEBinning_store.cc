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
  RooRealVar* Par2 = new RooRealVar("s23","par2",TMath::ASin(TMath::Sqrt(1.))/2.,0,100);
  RooRealVar* Par3 = new RooRealVar("s13","par3",TMath::ASin(TMath::Sqrt(0.1))/2.,0,100);
  RooRealVar* Par4 = new RooRealVar("delta","par4",-1.5,-10,10);
  RooRealVar* Par5 = new RooRealVar("dm21","par5",0.000075,-10,10);
  RooRealVar* Par6 = new RooRealVar("dm32","par6",0.00244,-10,10);
  RooRealVar* Par7 = new RooRealVar("dm31","par7",0.00252,-10,10);
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


TMatrixD* Sterile::prepareT2kCovMatrix(TMatrixD* covM_t2k, TVectorD* fVec_t2k, Int_t nBins) const
{
Int_t nbin = nBins*1.5;
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


Double_t Sterile::surv_t2k(Double_t L, Double_t E, RooListProxy* _pulls) const
{
  Double_t prob = 1 - 4*( TMath::Power(TMath::Sin(((RooAbsReal*)_pulls->at(0))->getVal()),2) * TMath::Power(TMath::Cos(((RooAbsReal*)_pulls->at(1))->getVal()),2)  +  TMath::Power(TMath::Sin(((RooAbsReal*)_pulls->at(2))->getVal()),2) * TMath::Power(TMath::Sin(((RooAbsReal*)_pulls->at(1))->getVal()),2) * TMath::Power(TMath::Cos(((RooAbsReal*)_pulls->at(0))->getVal()),2)   + 2 * TMath::Sin(((RooAbsReal*)_pulls->at(0))->getVal()) * TMath::Sin(((RooAbsReal*)_pulls->at(2))->getVal()) * TMath::Sin(((RooAbsReal*)_pulls->at(1))->getVal()) * TMath::Cos(((RooAbsReal*)_pulls->at(0))->getVal()) * TMath::Cos(((RooAbsReal*)_pulls->at(1))->getVal()) * TMath::Cos(((RooAbsReal*)_pulls->at(3))->getVal())  )
                  * TMath::Power(TMath::Sin(((RooAbsReal*)_pulls->at(1))->getVal()),2) * TMath::Power(TMath::Cos(((RooAbsReal*)_pulls->at(2))->getVal()),2) * TMath::Power(TMath::Sin(((RooAbsReal*)_pulls->at(6))->getVal()* L*1.27 /E ),2)
   
                  - 4*( TMath::Power(TMath::Cos(((RooAbsReal*)_pulls->at(0))->getVal()),2) * TMath::Power(TMath::Cos(((RooAbsReal*)_pulls->at(1))->getVal()),2)  +  TMath::Power(TMath::Sin(((RooAbsReal*)_pulls->at(2))->getVal()),2) * TMath::Power(TMath::Sin(((RooAbsReal*)_pulls->at(1))->getVal()),2) * TMath::Power(TMath::Sin(((RooAbsReal*)_pulls->at(0))->getVal()),2)   - 2 * TMath::Sin(((RooAbsReal*)_pulls->at(0))->getVal()) * TMath::Sin(((RooAbsReal*)_pulls->at(2))->getVal()) * TMath::Sin(((RooAbsReal*)_pulls->at(1))->getVal()) * TMath::Cos(((RooAbsReal*)_pulls->at(0))->getVal()) * TMath::Cos(((RooAbsReal*)_pulls->at(1))->getVal()) * TMath::Cos(((RooAbsReal*)_pulls->at(3))->getVal())  )
                  * TMath::Power(TMath::Sin(((RooAbsReal*)_pulls->at(1))->getVal()),2) * TMath::Power(TMath::Cos(((RooAbsReal*)_pulls->at(2))->getVal()),2) * TMath::Power(TMath::Sin(((RooAbsReal*)_pulls->at(5))->getVal()* L*1.27 /E ),2)

                  - 4*( TMath::Power(TMath::Sin(((RooAbsReal*)_pulls->at(0))->getVal()),2) * TMath::Power(TMath::Cos(((RooAbsReal*)_pulls->at(1))->getVal()),2)  +  TMath::Power(TMath::Sin(((RooAbsReal*)_pulls->at(2))->getVal()),2) * TMath::Power(TMath::Sin(((RooAbsReal*)_pulls->at(1))->getVal()),2) * TMath::Power(TMath::Cos(((RooAbsReal*)_pulls->at(0))->getVal()),2)   + 2 * TMath::Sin(((RooAbsReal*)_pulls->at(0))->getVal()) * TMath::Sin(((RooAbsReal*)_pulls->at(2))->getVal()) * TMath::Sin(((RooAbsReal*)_pulls->at(1))->getVal()) * TMath::Cos(((RooAbsReal*)_pulls->at(0))->getVal()) * TMath::Cos(((RooAbsReal*)_pulls->at(1))->getVal()) * TMath::Cos(((RooAbsReal*)_pulls->at(3))->getVal())  )
                  * ( TMath::Power(TMath::Cos(((RooAbsReal*)_pulls->at(0))->getVal()),2) * TMath::Power(TMath::Cos(((RooAbsReal*)_pulls->at(1))->getVal()),2)  +  TMath::Power(TMath::Sin(((RooAbsReal*)_pulls->at(2))->getVal()),2) * TMath::Power(TMath::Sin(((RooAbsReal*)_pulls->at(1))->getVal()),2) * TMath::Power(TMath::Sin(((RooAbsReal*)_pulls->at(0))->getVal()),2)   - 2 * TMath::Sin(((RooAbsReal*)_pulls->at(0))->getVal()) * TMath::Sin(((RooAbsReal*)_pulls->at(2))->getVal()) * TMath::Sin(((RooAbsReal*)_pulls->at(1))->getVal()) * TMath::Cos(((RooAbsReal*)_pulls->at(0))->getVal()) * TMath::Cos(((RooAbsReal*)_pulls->at(1))->getVal()) * TMath::Cos(((RooAbsReal*)_pulls->at(3))->getVal())  )
                  * TMath::Power(TMath::Sin(((RooAbsReal*)_pulls->at(4))->getVal()* L*1.27 /E ),2) ; 
 //   std::cout<<prob<<std::endl;
    return prob;       
}



Double_t Sterile::app_t2k(Double_t L, Double_t E,  Double_t density, RooListProxy* _pulls) const
{
  Double_t prob = 4 * TMath::Power(TMath::Cos(((RooAbsReal*)_pulls->at(2))->getVal()),2) * TMath::Power(TMath::Sin(((RooAbsReal*)_pulls->at(2))->getVal()),2) * TMath::Power(TMath::Sin(((RooAbsReal*)_pulls->at(1))->getVal()),2) * TMath::Power(TMath::Sin(((RooAbsReal*)_pulls->at(6))->getVal()* L*1.27 /E ),2) * (1+ (2*7.56e-5 * density * E / ((RooAbsReal*)_pulls->at(6))->getVal()) * (1 - 2* TMath::Power(TMath::Sin(((RooAbsReal*)_pulls->at(2))->getVal()),2)))
                + 8 *  TMath::Power(TMath::Cos(((RooAbsReal*)_pulls->at(2))->getVal()),2) * TMath::Sin(((RooAbsReal*)_pulls->at(0))->getVal()) * TMath::Sin(((RooAbsReal*)_pulls->at(2))->getVal()) * TMath::Sin(((RooAbsReal*)_pulls->at(1))->getVal())
                * (TMath::Cos(((RooAbsReal*)_pulls->at(0))->getVal()) * TMath::Cos(((RooAbsReal*)_pulls->at(1))->getVal()) * TMath::Cos(((RooAbsReal*)_pulls->at(3))->getVal()) - TMath::Sin(((RooAbsReal*)_pulls->at(0))->getVal()) * TMath::Sin(((RooAbsReal*)_pulls->at(2))->getVal()) * TMath::Sin(((RooAbsReal*)_pulls->at(1))->getVal()) )
                * TMath::Cos(((RooAbsReal*)_pulls->at(5))->getVal()* L*1.27 /E ) * TMath::Sin(((RooAbsReal*)_pulls->at(6))->getVal()* L*1.27 /E ) * TMath::Sin(((RooAbsReal*)_pulls->at(4))->getVal()* L*1.27 /E )

                - 8 * TMath::Power(TMath::Cos(((RooAbsReal*)_pulls->at(2))->getVal()),2) * TMath::Cos(((RooAbsReal*)_pulls->at(0))->getVal()) * TMath::Cos(((RooAbsReal*)_pulls->at(1))->getVal()) * TMath::Sin(((RooAbsReal*)_pulls->at(0))->getVal()) * TMath::Sin(((RooAbsReal*)_pulls->at(2))->getVal()) * TMath::Sin(((RooAbsReal*)_pulls->at(1))->getVal()) * TMath::Sin(((RooAbsReal*)_pulls->at(3))->getVal())* TMath::Sin(((RooAbsReal*)_pulls->at(5))->getVal()* L*1.27 /E ) * TMath::Sin(((RooAbsReal*)_pulls->at(6))->getVal()* L*1.27 /E ) * TMath::Sin(((RooAbsReal*)_pulls->at(4))->getVal()* L*1.27 /E )
                
                + 4* TMath::Power(TMath::Sin(((RooAbsReal*)_pulls->at(0))->getVal()),2) * TMath::Power(TMath::Cos(((RooAbsReal*)_pulls->at(2))->getVal()),2)
                * (  TMath::Power(TMath::Cos(((RooAbsReal*)_pulls->at(0))->getVal()),2) * TMath::Power(TMath::Cos(((RooAbsReal*)_pulls->at(1))->getVal()),2) + TMath::Power(TMath::Sin(((RooAbsReal*)_pulls->at(0))->getVal()),2) * TMath::Power(TMath::Sin(((RooAbsReal*)_pulls->at(1))->getVal()),2) * TMath::Power(TMath::Sin(((RooAbsReal*)_pulls->at(2))->getVal()),2) - 2* TMath::Cos(((RooAbsReal*)_pulls->at(0))->getVal()) * TMath::Cos(((RooAbsReal*)_pulls->at(1))->getVal()) * TMath::Sin(((RooAbsReal*)_pulls->at(0))->getVal()) 
                * TMath::Sin(((RooAbsReal*)_pulls->at(1))->getVal()) * TMath::Sin(((RooAbsReal*)_pulls->at(2))->getVal()) * TMath::Cos(((RooAbsReal*)_pulls->at(3))->getVal()) ) * TMath::Power(TMath::Sin(((RooAbsReal*)_pulls->at(4))->getVal()* L*1.27 /E ),2)
        
                - 8 * TMath::Power(TMath::Cos(((RooAbsReal*)_pulls->at(2))->getVal()),2) * TMath::Power(TMath::Sin(((RooAbsReal*)_pulls->at(2))->getVal()),2) * TMath::Power(TMath::Sin(((RooAbsReal*)_pulls->at(1))->getVal()),2) * (1 - 2 * TMath::Power(TMath::Sin(((RooAbsReal*)_pulls->at(2))->getVal()),2))
                * 7.56e-5 * density * E  * L / (4. * E) * TMath::Cos(((RooAbsReal*)_pulls->at(5))->getVal()* L*1.27 /E ) * TMath::Sin(((RooAbsReal*)_pulls->at(6))->getVal()* L*1.27 /E );
       
    return prob;
}

Double_t Sterile::survNue_sk(Double_t L, Double_t E, Double_t density, RooListProxy* _pulls) const
{
  Double_t prob = 1 + (TMath::Power(TMath::Sin(2 * ((RooAbsReal*)_pulls->at(0))->getVal()),2) * TMath::Sin(((RooAbsReal*)_pulls->at(4))->getVal()* L*1.27 /E )) * (2 * TMath::Power(TMath::Cos(((RooAbsReal*)_pulls->at(1))->getVal()),2) - 1) ;
    return prob;       
}



Double_t Sterile::survNumu_sk(Double_t L, Double_t E,   RooListProxy* _pulls) const
{
  Double_t phi = (((RooAbsReal*)_pulls->at(6))->getVal() + TMath::Power(TMath::Sin(((RooAbsReal*)_pulls->at(0))->getVal()),2) * ((RooAbsReal*)_pulls->at(4))->getVal()) * L*2.54 /  E;

  Double_t P_ex = (TMath::Power(TMath::Sin(2 * ((RooAbsReal*)_pulls->at(0))->getVal()),2) * TMath::Sin(((RooAbsReal*)_pulls->at(4))->getVal()* L*1.27 /E ));

  Double_t prob = 1 - (TMath::Power(TMath::Sin(2 * ((RooAbsReal*)_pulls->at(0))->getVal()),2) * TMath::Sin(((RooAbsReal*)_pulls->at(4))->getVal()* L*1.27 /E )) * (0.5 * TMath::Power(TMath::Cos(((RooAbsReal*)_pulls->at(1))->getVal()),2) ) * (2 * TMath::Power(TMath::Cos(((RooAbsReal*)_pulls->at(1))->getVal()),2) - 1) 
                    - 0.5 * TMath::Power(TMath::Sin(2 * ((RooAbsReal*)_pulls->at(1))->getVal()),2) * (1- TMath::Sqrt(1 - P_ex) * TMath::Cos(phi));
      
  return prob;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Double_t Sterile::surv_t2k(Double_t L, Double_t E, TVectorD* parVec) const
{
  Double_t prob = 1 - 4*( TMath::Power(TMath::Sin((*parVec)[0]),2) * TMath::Power(TMath::Cos((*parVec)[1]),2)  +  TMath::Power(TMath::Sin((*parVec)[2]),2) * TMath::Power(TMath::Sin((*parVec)[1]),2) * TMath::Power(TMath::Cos((*parVec)[0]),2)   + 2 * TMath::Sin((*parVec)[0]) * TMath::Sin((*parVec)[2]) * TMath::Sin((*parVec)[1]) * TMath::Cos((*parVec)[0]) * TMath::Cos((*parVec)[1]) * TMath::Cos((*parVec)[3])  )
                  * TMath::Power(TMath::Sin((*parVec)[1]),2) * TMath::Power(TMath::Cos((*parVec)[2]),2) * TMath::Power(TMath::Sin((*parVec)[6]* L*1.27 /E ),2)
   
                  - 4*( TMath::Power(TMath::Cos((*parVec)[0]),2) * TMath::Power(TMath::Cos((*parVec)[1]),2)  +  TMath::Power(TMath::Sin((*parVec)[2]),2) * TMath::Power(TMath::Sin((*parVec)[1]),2) * TMath::Power(TMath::Sin((*parVec)[0]),2)   - 2 * TMath::Sin((*parVec)[0]) * TMath::Sin((*parVec)[2]) * TMath::Sin((*parVec)[1]) * TMath::Cos((*parVec)[0]) * TMath::Cos((*parVec)[1]) * TMath::Cos((*parVec)[3])  )
                  * TMath::Power(TMath::Sin((*parVec)[1]),2) * TMath::Power(TMath::Cos((*parVec)[2]),2) * TMath::Power(TMath::Sin((*parVec)[5]* L*1.27  /E ),2)

                  - 4*( TMath::Power(TMath::Sin((*parVec)[0]),2) * TMath::Power(TMath::Cos((*parVec)[1]),2)  +  TMath::Power(TMath::Sin((*parVec)[2]),2) * TMath::Power(TMath::Sin((*parVec)[1]),2) * TMath::Power(TMath::Cos((*parVec)[0]),2)   + 2 * TMath::Sin((*parVec)[0]) * TMath::Sin((*parVec)[2]) * TMath::Sin((*parVec)[1]) * TMath::Cos((*parVec)[0]) * TMath::Cos((*parVec)[1]) * TMath::Cos((*parVec)[3])  )
                  * ( TMath::Power(TMath::Cos((*parVec)[0]),2) * TMath::Power(TMath::Cos((*parVec)[1]),2)  +  TMath::Power(TMath::Sin((*parVec)[2]),2) * TMath::Power(TMath::Sin((*parVec)[1]),2) * TMath::Power(TMath::Sin((*parVec)[0]),2)   - 2 * TMath::Sin((*parVec)[0]) * TMath::Sin((*parVec)[2]) * TMath::Sin((*parVec)[1]) * TMath::Cos((*parVec)[0]) * TMath::Cos((*parVec)[1]) * TMath::Cos((*parVec)[3])  )
                  * TMath::Power(TMath::Sin((*parVec)[4]* L*1.27  /E ),2) ; 
 //   std::cout<<prob<<std::endl;
    return prob;       
}



Double_t Sterile::app_t2k(Double_t L, Double_t E,  Double_t density, TVectorD* parVec) const
{
  Double_t prob = 4 * TMath::Power(TMath::Cos((*parVec)[2]),2) * TMath::Power(TMath::Sin((*parVec)[2]),2) * TMath::Power(TMath::Sin((*parVec)[1]),2) * TMath::Power(TMath::Sin((*parVec)[6]* L*1.27  /E ),2) * (1+ (2*7.56e-5 * density * E / (*parVec)[6]) * (1 - 2* TMath::Power(TMath::Sin((*parVec)[2]),2)))
                + 8 *  TMath::Power(TMath::Cos((*parVec)[2]),2) * TMath::Sin((*parVec)[0]) * TMath::Sin((*parVec)[2]) * TMath::Sin((*parVec)[1])
                * (TMath::Cos((*parVec)[0]) * TMath::Cos((*parVec)[1]) * TMath::Cos((*parVec)[3]) - TMath::Sin((*parVec)[0]) * TMath::Sin((*parVec)[2]) * TMath::Sin((*parVec)[1]) )
                * TMath::Cos((*parVec)[5]* L*1.27  /E ) * TMath::Sin((*parVec)[6]* L*1.27  /E ) * TMath::Sin((*parVec)[4]* L*1.27  /E )

                - 8 * TMath::Power(TMath::Cos((*parVec)[2]),2) * TMath::Cos((*parVec)[0]) * TMath::Cos((*parVec)[1]) * TMath::Sin((*parVec)[0]) * TMath::Sin((*parVec)[2]) * TMath::Sin((*parVec)[1]) * TMath::Sin((*parVec)[3])* TMath::Sin((*parVec)[5]* L*1.27  /E ) * TMath::Sin((*parVec)[6]* L*1.27  /E ) * TMath::Sin((*parVec)[4]* L*1.27  /E )
                
                + 4* TMath::Power(TMath::Sin((*parVec)[0]),2) * TMath::Power(TMath::Cos((*parVec)[2]),2)
                * (  TMath::Power(TMath::Cos((*parVec)[0]),2) * TMath::Power(TMath::Cos((*parVec)[1]),2) + TMath::Power(TMath::Sin((*parVec)[0]),2) * TMath::Power(TMath::Sin((*parVec)[1]),2) * TMath::Power(TMath::Sin((*parVec)[2]),2) - 2* TMath::Cos((*parVec)[0]) * TMath::Cos((*parVec)[1]) * TMath::Sin((*parVec)[0]) 
                * TMath::Sin((*parVec)[1]) * TMath::Sin((*parVec)[2]) * TMath::Cos((*parVec)[3]) ) * TMath::Power(TMath::Sin((*parVec)[4]* L*1.27  /E ),2)
        
                - 8 * TMath::Power(TMath::Cos((*parVec)[2]),2) * TMath::Power(TMath::Sin((*parVec)[2]),2) * TMath::Power(TMath::Sin((*parVec)[1]),2) * (1 - 2 * TMath::Power(TMath::Sin((*parVec)[2]),2))
                * 7.56e-5 * density * E  * L / (4. * E) * TMath::Cos((*parVec)[5]* L*1.27  /E ) * TMath::Sin((*parVec)[6]* L*1.27  /E );
       
    return prob;
}

Double_t Sterile::survNue_sk(Double_t L, Double_t E, Double_t density, TVectorD* parVec) const
{
  Double_t prob = 1 + (TMath::Power(TMath::Sin(2 * (*parVec)[0]),2) * TMath::Sin((*parVec)[4]* L*1.27  /E )) * (2 * TMath::Power(TMath::Cos((*parVec)[1]),2) - 1) ;
    return prob;       
}



Double_t Sterile::survNumu_sk(Double_t L, Double_t E,   TVectorD* parVec) const
{
  Double_t phi = ((*parVec)[6] + TMath::Power(TMath::Sin((*parVec)[0]),2) * (*parVec)[4]) * L*2.54 /  E;

  Double_t P_ex = (TMath::Power(TMath::Sin(2 * (*parVec)[0]),2) * TMath::Sin((*parVec)[4]* L*1.27  /E ));

  Double_t prob = 1 - (TMath::Power(TMath::Sin(2 * (*parVec)[0]),2) * TMath::Sin((*parVec)[4]* L*1.27  /E )) * (0.5 * TMath::Power(TMath::Cos((*parVec)[1]),2) ) * (2 * TMath::Power(TMath::Cos((*parVec)[1]),2) - 1) 
                    - 0.5 * TMath::Power(TMath::Sin(2 * (*parVec)[1]),2) * (1- TMath::Sqrt(1 - P_ex) * TMath::Cos(phi));
      
  return prob;
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



Double_t Sterile ::FillEv( RooListProxy* _pulls ) const
{

   Int_t nBins = _Bins;
   Double_t nCorr1 = _Corr1;
   Double_t nCorr2 = _Corr2;
   Double_t nCorr3 = _Corr3;
   Double_t nCorr4 = _Corr4;
  
   Double_t intPoint = 0.;
   Double_t intercept = 3.;
   Double_t endPoint = 16.;


   Double_t numuVec_t2k[nBins];
   Double_t nueVec_t2k[nBins];

      TF1* sumFunc = new TF1("sumFunc","[0]*exp(-0.5 * TMath::Power((x-[1])/[2],2)) + [3]*exp(-0.5 * TMath::Power((x-[4])/[5],2)) ",0.,endPoint );
 //+ [6] * TMath::Erf((x-0.8)*2) ",0.,endPoint );

      sumFunc->SetParameter(0,1600 * _time);
      sumFunc->SetParameter(1,intercept );
      sumFunc->SetParameter(2,1.);
      sumFunc->SetParameter(3,2100 * _time);
      sumFunc->SetParameter(4,intercept );
      sumFunc->SetParameter(5,1.5);
      sumFunc->SetParameter(6,20);


         for(Int_t iBin=0; iBin<nBins ; ++iBin) {
	Double_t content = sumFunc ->Integral(iBin*(endPoint - intPoint)/nBins, (iBin+1)*(endPoint - intPoint)/nBins);
	numuVec_t2k[iBin]=content;
	nueVec_t2k[iBin]=0;

	numuVec_t2k_nd[iBin]=content;
	nueVec_t2k_nd[iBin]=0;
         }




   TMatrixD* covM_t2k = new TMatrixD(nBins+nBins,nBins+nBins);

   for(Int_t i=0;i<nBins;i++){ (*covM_t2k)(i,i) = 0.077*0.077 ; }
   for(Int_t i=nBins;i<nBins*2;i++){ (*covM_t2k)(i,i) = 0.068*0.068 ; }    

   // nCorr1 and nCorr2 were set to 0.5 and 0.12 respectively.
   //cout<<nCorr1<<" "<<nCorr2<<endl;

   for(Int_t i=0;i<nBins*2;i++){
      for(Int_t j=0;j<nBins*2;j++){
          if(i!=j) (*covM_t2k)(i,j) = nCorr1  * TMath::Sqrt((*covM_t2k)(i,i)*(*covM_t2k)(j,j)) ;
       }
    }


   for(Int_t i=0;i<nBins ;i++){
      for(Int_t j=nBins ;j<nBins *2;j++){
        (*covM_t2k)(i,j) = nCorr2  * TMath::Sqrt((*covM_t2k)(i,j-nBins )*(*covM_t2k)(i+nBins ,j)) ;
        (*covM_t2k)(j,i) = (*covM_t2k)(i,j) ;
       }
    }
   
   covM_t2k->ResizeTo(nBins *1.5,nBins *1.5);


   TMatrixD* covM_t2k_nd = new TMatrixD(nBins+nBins,nBins+nBins);

   for(Int_t i=0;i<nBins;i++){ (*covM_t2k_nd)(i,i) = 0.077*0.077 ; }
   for(Int_t i=nBins;i<nBins*2;i++){ (*covM_t2k_nd)(i,i) = 0.068*0.068 ; }    

   for(Int_t i=0;i<nBins*2;i++){
      for(Int_t j=0;j<nBins*2;j++){
          if(i!=j) (*covM_t2k_nd)(i,j) = nCorr3  * TMath::Sqrt((*covM_t2k_nd)(i,i)*(*covM_t2k_nd)(j,j)) ;
       }
    }


   for(Int_t i=0;i<nBins ;i++){
      for(Int_t j=nBins ;j<nBins *2;j++){
        (*covM_t2k_nd)(i,j) = nCorr2  * TMath::Sqrt((*covM_t2k_nd)(i,j-nBins )*(*covM_t2k_nd)(i+nBins ,j)) ;
        (*covM_t2k_nd)(j,i) = (*covM_t2k_nd)(i,j) ;
       }
    }
   
   covM_t2k_nd->ResizeTo(nBins *1.5,nBins *1.5);



   TMatrixD* covM_t2k_tot = new TMatrixD(nBins+nBins+nBins,nBins+nBins+nBins);
      for(Int_t i=0;i<nBins*1.5 ;i++){
      for(Int_t j=0 ;j<nBins *1.5;j++){
        (*covM_t2k_tot)(i,j) = (*covM_t2k)(i,j) ;
       }
    }

      for(Int_t i=nBins*1.5;i<nBins*3 ;i++){
      for(Int_t j=nBins*1.5 ;j<nBins *3;j++){
        (*covM_t2k_tot)(i,j) = (*covM_t2k_nd)(i,j) ;
       }
    }

   for(Int_t i=0;i<nBins*1.5 ;i++){
      for(Int_t j=nBins*1.5 ;j<nBins *3;j++){
        (*covM_t2k_tot)(i,j) = nCorr4  * TMath::Sqrt((*covM_t2k_tot)(i,j-nBins )*(*covM_t2k_tot)(i+nBins ,j)) ;
        (*covM_t2k_tot)(j,i) = (*covM_t2k_tot)(i,j) ;
       }
    }




   Int_t nSpec = 2;
   TH1D* allVec[nSpec];
   allVec[0] = new TH1D("","",nBins ,0.,endPoint );
   allVec[1] = new TH1D("","",nBins /2,0.,endPoint/2);

   TVectorD* fVec_t2k = new TVectorD(nBins *1.5);
   TVectorD* fVec_t2kShadow = new TVectorD(nBins);

   for(Int_t i=0;i<nBins ;i++){ (*fVec_t2k)[i] = numuVec_t2k[i] * 
        this->surv_t2k(1300, allVec[0]->GetBinCenter(i+1) , _pulls) * ((RooAbsReal*)_pulls->at(7))->getVal() * ((RooAbsReal*)_pulls->at(9))->getVal() ; 
                      }
   for(Int_t i=nBins ;i<nBins *1.5;i++){ (*fVec_t2k)[i] = numuVec_t2k[i-nBins ] * this->app_t2k(1300, allVec[1]->GetBinCenter(i+1-nBins ) , _Density,  _pulls) * ((RooAbsReal*)_pulls->at(8))->getVal() * ((RooAbsReal*)_pulls->at(10))->getVal() ; 
        (*fVec_t2kShadow)[i-nBins] = (*fVec_t2k)[i];
                      }

                           
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  TVectorD* dataIn = new TVectorD(11);
  (*dataIn )[0] = TMath::ASin(TMath::Sqrt(0.85))/2.;
  (*dataIn )[1] = TMath::ASin(TMath::Sqrt(1.))/2.;
  (*dataIn )[2] = TMath::ASin(TMath::Sqrt(0.1))/2.;
  (*dataIn )[3] = 0;
  (*dataIn )[4] = 0.000075;
  (*dataIn )[5] = 0.00244;
  (*dataIn )[6] = 0.00252;
  (*dataIn )[7] = 1;
  (*dataIn )[8] = 1;
  (*dataIn )[9] = 1;
  (*dataIn )[10] = 1;

 //  std::cout<<"setting up data "<<std::endl;
   TVectorD* fData_t2k = new TVectorD(nBins *1.5);

   for(Int_t i=0;i<nBins ;i++){ (*fData_t2k)[i] = numuVec_t2k[i] * this->surv_t2k(1300, allVec[0]->GetBinCenter(i+1) , dataIn )  ;  }
   for(Int_t i=nBins ;i<nBins *1.5;i++){ (*fData_t2k)[i] = numuVec_t2k[i-nBins ] * this->app_t2k(1300, allVec[1]->GetBinCenter(i+1-nBins ), _Density , dataIn ) ; }

 //   for(Int_t i=nBins ;i<nBins *1.5;i++){  cout<<(*dataIn)[3]<< " "<< ((RooAbsReal*)_pulls->at(3))->getVal() <<" "<< this->app_t2k(1300, allVec[1]->GetBinCenter(i+1-nBins ), _Density , dataIn )<<" "<<this->app_t2k(1300, allVec[1]->GetBinCenter(i+1-nBins ) , _Density,  _pulls)<<endl; }
 //   cout<<endl;


TH1D* hShadow_t2k  =  new TH1D("","",nBins,0,endPoint );
TH1D* hData_t2k  = new TH1D("","",nBins,0,endPoint );
for(Int_t i=nBins;i<1.5*nBins; i++){
hShadow_t2k  ->SetBinContent(i-nBins+1,(*fVec_t2kShadow  )(i-nBins));
hData_t2k  ->SetBinContent(i-nBins+1,(*fData_t2k  )(i));
}

  ofstream out;
  out.open(Form("bin%d.txt",nBins),ios::out);
     for(Int_t i=nBins;i<nBins*1.5;i++){ out<<hShadow_t2k  ->GetBinContent(i-nBins+1)<<" "<<hData_t2k  ->GetBinContent(i-nBins+1)<<endl; }


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


   TMatrixD* covMatT2k = this->prepareT2kCovMatrix(covM_t2k, fVec_t2k, nBins);

 //  fVec_JUNO->Print();
 //  std::cout<<"Ok"<<std::endl;
 //  fData_JUNO->Print();
   
   for(Int_t i=0;i<nBins *1.5; i++){
       (*fVec_t2k)[i] -= (*fData_t2k)[i];
         if( (*covMatT2k)(i,i) ==0 ){(*covMatT2k)(i,i) = 1000000000;}
     }
 //   fVec_t2k->Print();

   covMatT2k->Invert();

   TVectorD mulVec(*fVec_t2k);
   mulVec *= (*covMatT2k);

   Double_t t2kResult = TMath::Abs(mulVec*(*fVec_t2k));
  std::cout<<"sans pull t2k: "<<t2kResult<<std::endl;

   return (Double_t) t2kResult;  
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


rep->setAtmBaseline(10000);
rep->setDensity(3.);
TH1D* vecInput1 = new TH1D("","",11,0,11);
TH1D* vecInput2 = new TH1D("","",11,0,11);
//vecInput1->Print();

vecInput1 ->SetBinContent(1, TMath::ASin(TMath::Sqrt(0.85))/2.);
vecInput1 ->SetBinContent(2,TMath::ASin(TMath::Sqrt(1.))/2.);
vecInput1 ->SetBinContent(3,TMath::ASin(TMath::Sqrt(0.1))/2.);
vecInput1 ->SetBinContent(4,-1.5);
vecInput1 ->SetBinContent(5,0.000075);
vecInput1 ->SetBinContent(6,0.00244);
vecInput1 ->SetBinContent(7,0.00252);
vecInput1 ->SetBinContent(8,1);
vecInput1 ->SetBinContent(9,1);
vecInput1 ->SetBinContent(10,1);
vecInput1 ->SetBinContent(11,1);

vecInput2 ->SetBinContent(1,TMath::ASin(TMath::Sqrt(0.85+0.021))/2.-TMath::ASin(TMath::Sqrt(0.85))/2.);
vecInput2 ->SetBinContent(2,2);
vecInput2 ->SetBinContent(3,TMath::ASin(TMath::Sqrt(0.1+0.005))/2.-TMath::ASin(TMath::Sqrt(0.1))/2.);
vecInput2 ->SetBinContent(4,10);
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

  Int_t binSetup = atoi(argv[2]);
  Double_t corr1Setup = atof(argv[3]);
  Double_t corr2Setup = atof(argv[4]);

rep->setNBins(binSetup);
rep->setCorr1(corr1Setup);
rep->setCorr2(corr2Setup);

rep->setTime(1);
rep->FillEv(rep->getPullList());


    const int npts = 100;
    double bin_x[npts]={};
    double low_x = 1, high_x = 10;
    double step_x = (high_x-low_x)/(double)npts;  
    for(int i=0;i<npts ;i++){bin_x[i]= low_x + i * step_x;}

  int aaa = (atof(argv[5])+0.000001)*100;
  ofstream out2;
  out2.open(Form("s%d_a%d_JUNO2D.txt",atoi(argv[2]),aaa ),ios::app);

RooArgList list("list");
list.add(*rep);
  sprintf(formula,"%s","@0");
RooFormulaVar* fcn = new RooFormulaVar("fit","fit",formula,list);

//rep->getParVar(3)->setVal(2.5);
//rep->getParVar(3)->setConstant(true);
//rep->getParVar(0)->setConstant(true);
//rep->getParVar(1)->setConstant(true);
//rep->getParVar(2)->setConstant(true);
rep->getParVar(3)->setConstant(true);
rep->getParVar(7)->setConstant(true);
rep->getParVar(8)->setConstant(true);
rep->getParVar(9)->setConstant(true);
rep->getParVar(10)->setConstant(true);

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
//m.minos(); 
res = m.save();
double bestFit = res->minNll(); 

std::cout<<"result list: "<<std::endl;
std::cout<<"chi2: "<<bestFit <<std::endl;

  double bb =  rep->getParVar(1)->getError();
  double dd =  rep->getParVar(3)->getError();

cout<<"errors are "<<bb<<" "<<dd<<endl;

for(Int_t i=0;i<11;i++){std::cout<<" "<<rep->getPar(i)<<std::endl;}

  int aa = (atof(argv[5])+0.000001)*100;

  out2<<atoi(argv[2])<<" "<<aa<<" "<<bb<<" "<<bestFit<<endl; 
  cout<<atoi(argv[2])<<" "<<aa<<" "<<bb<<" "<<bestFit <<endl; 

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

void Sterile :: setCorr1(Double_t Corr1){
_Corr1= Corr1;
}

void Sterile :: setCorr2(Double_t Corr2){
_Corr2= Corr2;
}

void Sterile :: setCorr3(Double_t Corr3){
_Corr3= Corr3;
}

void Sterile :: setCorr4(Double_t Corr4){
_Corr4= Corr4;
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

