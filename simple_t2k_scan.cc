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


  Sterile ::Sterile (const char* name) 
  : RooAbsReal(name,name)
{

// there will be: pull 0-2 s12, s23, s13, pull 3, CP delta, pull 4-6 dm2(21,32,31), pull 7-10 numu, nue Xsec and numu, nue selections

  _pulls     = new RooListProxy("_pulls","_pulls",this);
  RooRealVar* Par1 = new RooRealVar("s12","par1",TMath::ASin(TMath::Sqrt(0.85))/2.,0,100);
  RooRealVar* Par2 = new RooRealVar("s23","par2",TMath::ASin(TMath::Sqrt(0.95))/2.,0,100);
  RooRealVar* Par3 = new RooRealVar("s13","par3",TMath::ASin(TMath::Sqrt(0.1))/2.,0,100);
  RooRealVar* Par4 = new RooRealVar("delta","par4",-1.5,-10,10);
  RooRealVar* Par5 = new RooRealVar("dm21","par5",0.0000753,-10,10);
  RooRealVar* Par6 = new RooRealVar("dm32","par6",0.00244,-10,10);
  RooRealVar* Par7 = new RooRealVar("dm31","par7",0.0025,-10,10);
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
//std::cout<<"t2k matrix sum "<<covmatT2k->Sum()<<std::endl;
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

   Double_t numuVec_t2k[20] = {0,0,6,36,66,85,84,66,41,28,  18,12,10,8,7,5,4,3,2,2};
   Double_t nueVec_t2k[10] = {0,0,0,0,0,0,0,0,0,0};

   TMatrixD* covM_t2k = new TMatrixD(40,40);

   if(withSK){
   for(Int_t i=0;i<20;i++){ (*covM_t2k)(i,i) = 0.077*0.077 - 0.05*0.05 - 0.04*0.04; }
   for(Int_t i=20;i<40;i++){ (*covM_t2k)(i,i) = 0.068*0.068 - 0.047*.0047 - 0.027*0.027; }
      }

     else{
   for(Int_t i=0;i<20;i++){ (*covM_t2k)(i,i) = 0.077*0.077 ; }
   for(Int_t i=20;i<40;i++){ (*covM_t2k)(i,i) = 0.068*0.068 ; }    
          }

   for(Int_t i=0;i<40;i++){
      for(Int_t j=0;j<40;j++){
          if(i!=j) (*covM_t2k)(i,j) = 0.5 * TMath::Sqrt((*covM_t2k)(i,i)*(*covM_t2k)(j,j)) ;
       }
    }


   for(Int_t i=0;i<20;i++){
      for(Int_t j=20;j<40;j++){
        (*covM_t2k)(i,j) = 0.12 * TMath::Sqrt((*covM_t2k)(i,j-20)*(*covM_t2k)(i+20,j)) ;
        (*covM_t2k)(j,i) = (*covM_t2k)(i,j) ;
       }
    }
   
   covM_t2k->ResizeTo(30,30);

   Double_t numuVec_sk[5]={1.58e-2, 1.77e-2, 1.86e-2,1.68e-2, 1.38e-2};
   Double_t nueVec_sk[5]={1.21e-2,1.46e-2,1.50e-2,1.37e-2,1.16e-2};

   TMatrixD* covM_sk = new TMatrixD(10,10);
   (*covM_sk)(0,0) = 0.21*0.21;
   (*covM_sk)(1,1) = 0.16*0.16;
   (*covM_sk)(2,2) = 0.15*0.15;
   (*covM_sk)(3,3) = 0.16*0.16;
   (*covM_sk)(4,4) = 0.18*0.18;
   (*covM_sk)(5,5) = 0.18*0.18;
   (*covM_sk)(6,6) = 0.17*0.17;
   (*covM_sk)(7,7) = 0.16*0.16;
   (*covM_sk)(8,8) = 0.15*0.15;
   (*covM_sk)(9,9) = 0.17*0.17;

   (*covM_sk)(0,1) = 0.84;
   (*covM_sk)(0,2) = 0.699;
   (*covM_sk)(0,3) = 0.574;
   (*covM_sk)(0,4) = 0.596;
   (*covM_sk)(1,2) = 0.934;
   (*covM_sk)(1,3) = 0.872;
   (*covM_sk)(1,4) = 0.816;
   (*covM_sk)(2,3) = 0.942;
   (*covM_sk)(2,4) = 0.889;
   (*covM_sk)(3,4) = 0.912;

   (*covM_sk)(1,0) = 0.84;
   (*covM_sk)(2,0) = 0.699;
   (*covM_sk)(3,0) = 0.574;
   (*covM_sk)(4,0) = 0.596;
   (*covM_sk)(2,1) = 0.934;
   (*covM_sk)(3,1) = 0.872;
   (*covM_sk)(4,1) = 0.816;
   (*covM_sk)(3,2) = 0.942;
   (*covM_sk)(4,2) = 0.889;
   (*covM_sk)(4,3) = 0.912;



   (*covM_sk)(5,6) = 0.948;
   (*covM_sk)(5,7) = 0.720;
   (*covM_sk)(5,8) = 0.512;
   (*covM_sk)(5,9) = 0.396;
   (*covM_sk)(6,7) = 0.88;
   (*covM_sk)(6,8) = 0.678;
   (*covM_sk)(6,9) = 0.58;
   (*covM_sk)(7,8) = 0.914;
   (*covM_sk)(7,9) = 0.836;
   (*covM_sk)(8,9) = 0.931;

   (*covM_sk)(6,5) = 0.948;
   (*covM_sk)(7,5) = 0.720;
   (*covM_sk)(8,5) = 0.512;
   (*covM_sk)(9,5) = 0.396;
   (*covM_sk)(7,6) = 0.88;
   (*covM_sk)(8,6) = 0.678;
   (*covM_sk)(9,6) = 0.58;
   (*covM_sk)(8,7) = 0.914;
   (*covM_sk)(9,7) = 0.836;
   (*covM_sk)(9,8) = 0.931;





   (*covM_sk)(0,5) = 0.843;
   (*covM_sk)(0,6) = 0.888;
   (*covM_sk)(0,7) = 0.809;
   (*covM_sk)(0,8) = 0.682;
   (*covM_sk)(0,9) = 0.537;
   (*covM_sk)(1,5) = 0.692;
   (*covM_sk)(1,6) = 0.832;
   (*covM_sk)(1,7) = 0.940;
   (*covM_sk)(1,8) = 0.899;
   (*covM_sk)(1,9) = 0.842;
   (*covM_sk)(2,5) = 0.520;
   (*covM_sk)(2,6) = 0.686;
   (*covM_sk)(2,7) = 0.896;
   (*covM_sk)(2,8) = 0.944;
   (*covM_sk)(2,9) = 0.900;
   (*covM_sk)(3,5) = 0.423;
   (*covM_sk)(3,6) = 0.597;
   (*covM_sk)(3,7) = 0.839;
   (*covM_sk)(3,8) = 0.915;
   (*covM_sk)(3,9) = 0.941;
   (*covM_sk)(4,5) = 0.400;
   (*covM_sk)(4,6) = 0.544;
   (*covM_sk)(4,7) = 0.771;
   (*covM_sk)(4,8) = 0.889;
   (*covM_sk)(4,9) = 0.874;

   (*covM_sk)(5,0) = 0.843;
   (*covM_sk)(6,0) = 0.888;
   (*covM_sk)(7,0) = 0.809;
   (*covM_sk)(8,0) = 0.682;
   (*covM_sk)(9,0) = 0.537;
   (*covM_sk)(5,1) = 0.692;
   (*covM_sk)(6,1) = 0.832;
   (*covM_sk)(7,1) = 0.940;
   (*covM_sk)(8,1) = 0.899;
   (*covM_sk)(9,1) = 0.842;
   (*covM_sk)(5,2) = 0.520;
   (*covM_sk)(6,2) = 0.686;
   (*covM_sk)(7,2) = 0.896;
   (*covM_sk)(8,2) = 0.944;
   (*covM_sk)(9,2) = 0.900;
   (*covM_sk)(5,3) = 0.423;
   (*covM_sk)(6,3) = 0.597;
   (*covM_sk)(7,3) = 0.839;
   (*covM_sk)(8,3) = 0.915;
   (*covM_sk)(9,3) = 0.941;
   (*covM_sk)(5,4) = 0.400;
   (*covM_sk)(6,4) = 0.544;
   (*covM_sk)(7,4) = 0.771;
   (*covM_sk)(8,4) = 0.889;
   (*covM_sk)(9,4) = 0.874;

for(Int_t i=0;i<5;i++){(*covM_sk)(i,i) -= (0.04*0.04+0.05*0.05); (*covM_sk)(i+5,i+5) -= (0.027*0.027+0.047*0.047); }

for(Int_t i=0;i<5;i++){
    for(Int_t j=0;j<5;j++){
          if (i != j) { (*covM_sk)(i,j) *= TMath::Sqrt( (*covM_sk)(i,i) * (*covM_sk)(j,j) ); (*covM_sk)(i,j) -= 0.05*0.05; }
       }
  }

for(Int_t i=5;i<10;i++){
    for(Int_t j=5;j<10;j++){
          if (i != j) { (*covM_sk)(i,j) *= TMath::Sqrt( (*covM_sk)(i,i) * (*covM_sk)(j,j) ); (*covM_sk)(i,j) -= 0.05*0.05; }
       }
  }

for(Int_t i=0;i<5;i++){
    for(Int_t j=5;j<10;j++){
          if (i != j) { (*covM_sk)(i,j) *= TMath::Sqrt( (*covM_sk)(i,j-5) * (*covM_sk)(i+5,j) ); (*covM_sk)(i,j) -= 0.05*0.05; (*covM_sk)(j,i) = (*covM_sk)(i,j); }
       }
  }
   



   Int_t nSpec = 4;
   TH1D* allVec[nSpec];
   allVec[0] = new TH1D("","",20,0.,2.);
   allVec[1] = new TH1D("","",10,0.,1.);
   Double_t bin1[6]={TMath::Power(10,-0.8),TMath::Power(10,-0.6),TMath::Power(10,-0.4),TMath::Power(10,-0.2),TMath::Power(10,0),TMath::Power(10,0.2)};
   Double_t bin2[6]={TMath::Power(10,-0.6),TMath::Power(10,-0.4),TMath::Power(10,-0.2),TMath::Power(10,0),TMath::Power(10,0.2),TMath::Power(10,0.4)};
   allVec[2] = new TH1D("","",5,bin1);
   allVec[3] = new TH1D("","",5,bin2);

   TVectorD* fVec_t2k = new TVectorD(30);
   TVectorD* fVec_sk = new TVectorD(10);

   for(Int_t i=0;i<20;i++){ (*fVec_t2k)[i] = numuVec_t2k[i] * this->surv_t2k(295, allVec[0]->GetBinCenter(i+1) , _pulls) * ((RooAbsReal*)_pulls->at(7))->getVal() * ((RooAbsReal*)_pulls->at(9))->getVal() ; }
   for(Int_t i=20;i<30;i++){ (*fVec_t2k)[i] = numuVec_t2k[i-20] * this->app_t2k(295, allVec[1]->GetBinCenter(i+1-20) , _Density,  _pulls) * ((RooAbsReal*)_pulls->at(8))->getVal() * ((RooAbsReal*)_pulls->at(10))->getVal() ; }

   for(Int_t i=0;i<5;i++){ (*fVec_sk)[i] = numuVec_sk[i] * this->survNumu_sk(_AtmBaseline, allVec[2]->GetBinCenter(i+1) , _pulls) * ((RooAbsReal*)_pulls->at(7))->getVal() * ((RooAbsReal*)_pulls->at(9))->getVal() ; }
   for(Int_t i=5;i<10;i++){ (*fVec_sk)[i] = nueVec_sk[i-5] * this->survNue_sk(_AtmBaseline, allVec[3]->GetBinCenter(i+1-5) , _Density,  _pulls) * ((RooAbsReal*)_pulls->at(8))->getVal() * ((RooAbsReal*)_pulls->at(10))->getVal() ; }

 // std::cout<<this->getDataSwitch()<<std::endl;
  Int_t counter;
 
  TVectorD* dataIn = new TVectorD(11);
  (*dataIn )[0] = TMath::ASin(TMath::Sqrt(0.85))/2.;
  (*dataIn )[1] = TMath::ASin(TMath::Sqrt(0.95))/2.;
  (*dataIn )[2] = TMath::ASin(TMath::Sqrt(0.1))/2.;
  (*dataIn )[3] = -1.5;
  (*dataIn )[4] = 0.0000753;
  (*dataIn )[5] = 0.00244;
  (*dataIn )[6] = 0.0025;
  (*dataIn )[7] = 1;
  (*dataIn )[8] = 1;
  (*dataIn )[9] = 1;
  (*dataIn )[10] = 1;

 //  std::cout<<"setting up data "<<std::endl;
   TVectorD* fData_t2k = new TVectorD(30);
   TVectorD* fData_sk = new TVectorD(10);

   for(Int_t i=0;i<20;i++){ (*fData_t2k)[i] = numuVec_t2k[i] * this->surv_t2k(295, allVec[0]->GetBinCenter(i+1) , dataIn )  ;  }
   for(Int_t i=20;i<30;i++){ (*fData_t2k)[i] = numuVec_t2k[i-20] * this->app_t2k(295, allVec[1]->GetBinCenter(i+1-20), _Density , dataIn ) ; }

   for(Int_t i=0;i<5;i++){ (*fData_sk)[i] = numuVec_sk[i] * this->survNumu_sk(_AtmBaseline, allVec[2]->GetBinCenter(i+1) , dataIn ) ; }
   for(Int_t i=5;i<10;i++){ (*fData_sk)[i] = nueVec_sk[i-5] * this->survNue_sk(_AtmBaseline, allVec[3]->GetBinCenter(i+1-5), _Density , dataIn ) ; }


   counter++;

   TMatrixD* covMatT2k = this->prepareT2kCovMatrix(covM_t2k, fVec_t2k);
   TMatrixD* covMatSk = this->prepareSkCovMatrix(covM_sk, fVec_sk);

 //  fVec_t2k->Print();
 //  std::cout<<"Ok"<<std::endl;
 //  fData_t2k->Print();
   
   for(Int_t i=0;i<30; i++){
       (*fVec_t2k)[i] -= (*fData_t2k)[i];
         if( (*covMatT2k)(i,i) ==0 ){(*covMatT2k)(i,i) = 1000000000;}
     }

   covMatT2k->Invert();

   TVectorD mulVec(*fVec_t2k);
   mulVec *= (*covMatT2k);

   Double_t t2kResult = TMath::Abs(mulVec*(*fVec_t2k));


//////////////////////////////////////////////////////////////////////////

   for(Int_t i=0;i<10; i++){
       (*fVec_sk)[i] -= (*fData_sk)[i];
         if( (*covMatSk)(i,i) ==0 ){(*covMatSk)(i,i) = 1000000000;}
     }

   covMatSk->Invert();

   TVectorD mulVec2(*fVec_sk);
   mulVec2 *= (*covMatSk);

   Double_t skResult = TMath::Abs(mulVec2*(*fVec_sk));

//  std::cout<<"sans pull t2k: "<<t2kResult<<"   sk: "<<skResult<<std::endl;

   if(withSK) return (Double_t) (skResult + t2kResult); 
   else  return (Double_t) t2kResult; 
}


Double_t Sterile ::evaluate() const
{ 

Double_t matPart = this->FillEv( _pulls);

Double_t extraPull = this -> ExtraPull (_pulls);
Double_t tot = matPart + extraPull; //If needed, add pull terms here.

return tot;

}

Double_t Sterile ::ExtraPull (RooListProxy* _pulls) const
{
Double_t pullAdd = 0;
if(withSK){
for(Int_t i=0;i<11;i++){
 pullAdd += TMath::Power(( ((RooAbsReal*)_pulls->at(i))->getVal() - (*pullCV)[i] ),2) / TMath::Power( (*pullUnc)[i],2) ;
    }
  }
else{
for(Int_t i=0;i<7;i++){
 pullAdd += TMath::Power(( ((RooAbsReal*)_pulls->at(i))->getVal() - (*pullCV)[i] ),2) / TMath::Power( (*pullUnc)[i],2) ;
    }
  }
// std::cout<<"extra pull penalty: "<<pullAdd<<std::endl;
 return pullAdd;
}

Double_t Sterile ::getPar(int i) {
(((RooAbsReal*)_pulls->at(i))->getVal());
}

RooRealVar* Sterile ::getParVar(int i) {
return ((RooRealVar*)_pulls->at(i));
}

TVectorD* Sterile ::getT2kVec() {
return fVec_t2k;
}

TVectorD* Sterile ::getSkVec() {
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

int main(int argc, char**argv){


RooFitResult* res;
Sterile * rep = new Sterile ("_rep");
  char formula[10];


rep->setAtmBaseline(atof(argv[1]));
rep->setDensity(atof(argv[2]));
TH1D* vecInput1 = new TH1D("","",11,0,11);
TH1D* vecInput2 = new TH1D("","",11,0,11);
//vecInput1->Print();

vecInput1 ->SetBinContent(1, TMath::ASin(TMath::Sqrt(0.85))/2.);
vecInput1 ->SetBinContent(2,TMath::ASin(TMath::Sqrt(0.95))/2.);
vecInput1 ->SetBinContent(3,TMath::ASin(TMath::Sqrt(0.1))/2.);
vecInput1 ->SetBinContent(4,-1.5);
vecInput1 ->SetBinContent(5,0.0000753);
vecInput1 ->SetBinContent(6,0.00244);
vecInput1 ->SetBinContent(7,0.0025);
vecInput1 ->SetBinContent(8,1);
vecInput1 ->SetBinContent(9,1);
vecInput1 ->SetBinContent(10,1);
vecInput1 ->SetBinContent(11,1);

vecInput2 ->SetBinContent(1,TMath::ASin(TMath::Sqrt(0.85+0.021))/2.-TMath::ASin(TMath::Sqrt(0.85))/2.);
vecInput2 ->SetBinContent(2,10);
vecInput2 ->SetBinContent(3,10);
vecInput2 ->SetBinContent(4,10);
vecInput2 ->SetBinContent(5,0.0000018);
vecInput2 ->SetBinContent(6,0.00006);
vecInput2 ->SetBinContent(7,0.00006);
vecInput2 ->SetBinContent(8,0.05);
vecInput2 ->SetBinContent(9,0.047);
vecInput2 ->SetBinContent(10,0.04);
vecInput2 ->SetBinContent(11,0.027);

rep->setPull(vecInput1);
rep->setPullUnc(vecInput2);
rep->addSK(true);

rep->FillEv(rep->getPullList());

RooArgList list("list");
list.add(*rep);
  sprintf(formula,"%s","@0");
RooFormulaVar* fcn = new RooFormulaVar("fit","fit",formula,list);

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

for(Int_t i=0;i<11;i++){std::cout<<" "<<rep->getPar(i)<<std::endl;}

 RooPlot *frame;

 rep ->getParVar(0)->setConstant(false);
 rep ->getParVar(1)->setConstant(false);
 frame = m.contour(*rep->getParVar(2),*rep->getParVar(3),1,2); 
 std::string contourpath = "./";
 contourpath += "t2k_scan_2.root";
 frame->SetLineColor(1);
 frame->SetLineWidth(1);
 frame->SetYTitle("#delta_{CP}");
 frame->SetXTitle("#theta_{13}");
 //frame->Draw();

 TFile frameout(contourpath.c_str(),"recreate");
 frame->Write("contourplot");
 frameout.Close();


/*
    ofstream out1;
    out1.open("t2kSk1D_s23.txt");

      const int npts = 50;
      double scanx[npts], scany[npts], deltascany[npts];
      const double xlo = 0., xhi = 1.57,
	gap = (xhi-xlo)/npts;
     
      RooFitResult *scanres;
      for(Int_t i = 0; i < npts; ++i) {
	scanx[i] = xlo + i*gap;
	rep ->getParVar(1)->setVal(xlo + i*gap);
	rep ->getParVar(1)->setConstant(true);
       
	RooMinuit m(*fcn);
	m.setStrategy(2);
       
	Double_t callsEDM[2] = {10500., 1.e-6};
	Int_t irf = 0;
	gMinuit->mnexcm("MIGRAD",callsEDM,2,irf);  
	scanres = m.save();
       
	scany[i] = scanres->minNll();
	deltascany[i] = scany[i] - bestFit ;

      out1<<i<<" "<<scanx[i]<<" "<<deltascany[i]<<std::endl;
      }
     
*/


}
