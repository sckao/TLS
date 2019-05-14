#ifndef Module_H
#define Module_H

#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <TMath.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <THStack.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TRandom3.h>
#include <TLeaf.h>
#include <TString.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TFrame.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TMatrixD.h>
#include <TMatrixDEigen.h>

#include <algorithm>

#include "AnaInput.h"
#include "MathTools.h"
#include "Match.h"

struct swcs {
   int i ;
   int j ;
   double x ;
   double y ;
   bool b ;
};

class Module {

public:

   Module();     
   ~Module();
  
   void Analysis() ; 
   void AnalyzeStepFAST() ; 
   void CompareXY() ;
   void GetPitch( vector<vec>& output ) ;
   //void GetXY( vector<vec>& output ) ;
   void GetXY( vector<vec>& input, vector<swcs>& output ) ;
   void AxisCorrection( vector<vec>& dataV, double theta, double phi ) ;
   void CenterCorrection( vector<vec>& wA, double& dX0, double& dY0 ) ;
   void CompareFFLY_PSET( vector<vec>& fflyV, vector<vec>& psetV0, double sx = 0, double sy = 0) ;

   pair<double,double> gsl_vector_mean( vector<double>& v ) ;
   pair<double,double> gsl_vector_linear_fit( vector<double>& vx, vector<double>& vy ) ;

private:

   AnaInput*    Input ;
   Match*       match ;

   string hfolder  ;
   string hFileName ;
   string plotType ;
   string plotName0 ;
   string outFileName ;
   string logFileName ;
   string logFileName1 ;
   string logFileName2 ;
   string logFileName3 ;
   string logFileName4 ;

   int debug ;
   int site ;
   int die ;
   int tob ;

   //ClassDef(Module, 1);
};

//#if !defined(__CINT__)
//    ClassImp(Module);
#endif

