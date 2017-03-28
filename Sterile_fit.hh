/*
 *  Sterile for NoVA header file.
 *
 *  Author: Guang Yang
 */

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <iomanip>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <sstream>
#include <TList.h>

#include <TROOT.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TH1.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TFile.h>
#include <TRint.h>
#include <TH2.h>
#include <TFormula.h>
#include <TF1.h>
#include <TF2.h>
#include <TMath.h>
#include <Math/DistFunc.h>
#include <TLine.h>
#include <TTree.h>
#include <TGraph.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TGraphErrors.h>
#include <TVirtualFFT.h>
#include <TFoamIntegrand.h>
#include <TMatrixD.h>
#include <TVectorT.h>
#include <TDecompChol.h>

#include <RooFit.h>
#include "RooGlobalFunc.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooDataSet.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooWorkspace.h"
#include "RooMinuit.h"
#include "RooNLLVar.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooDataSet.h"
#include "RooExtendPdf.h"
#include "RooChi2Var.h"
#include "RooMinuit.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooRandom.h"
#include <RooMsgService.h>
#include <RooHist.h>
#include <RooTrace.h>
#include <RooCategory.h>
#include "RooConstVar.h"
#include "RooBinning.h"

#include "TStopwatch.h"
#include "TFile.h"
#include "TMinuit.h"

#include "RooFit.h"
#include "RooMinuit.h"
#include "RooFitResult.h"
#include "TMinuit.h"
#include <RooRealVar.h>

using namespace RooFit;


  class Sterile : public RooAbsReal {

  public:

    Sterile (const char* name);

    Sterile (const Sterile & other, const char* name = 0): RooAbsReal(other,name) {};
    virtual TObject* clone(const char* newname) const {return new Sterile (*this, newname);};
    virtual ~Sterile () ;

    Sterile (const Sterile & Sterile );

    void randomGo(int newIsotope, RooListProxy* _pulls);

    Sterile & operator=(const Sterile & rhs);

    RooFormulaVar* Chi2() ;

    Double_t FillEv(RooListProxy* _pulls) const;

    TMatrixD* prepareMatrix(TH2D* conv) const;
 
    void setSyst(Double_t syst) ;
    void setdm2CV(Double_t dm2CV) ;
    void setdm2Unc(Double_t dm2Unc) ;

    TMatrixD* prepareCovMatrix(TH2D* conv, TVectorD* vec, Double_t syst) const;

    RooRealVar Par1 ;
    RooRealVar Par2 ;
    RooRealVar Par3 ;
    RooRealVar Par4 ;

    Double_t _par1;
    Double_t _par2;

    Double_t getPar(int i) ;

    RooRealVar* getParVar(int i) ;
    RooListProxy* getParVar() ;

    RooArgList _parlist;
    RooListProxy* _pulls;

   TH2D* conv;
   TH1D* nueAfter;
   TH1D* nueBefore;
   TH1D* numuTrue;

   Double_t _syst;
   Double_t _dm2CV;
   Double_t _dm2Unc;

    //  private:
   
  virtual  Double_t evaluate() const ;

  protected:


  };


