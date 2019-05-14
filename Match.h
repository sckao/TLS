#ifndef Match_H
#define Match_H

#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <math.h>

/*
#include <TMath.h>
#include <TH2.h>
#include <THStack.h>
#include <TTree.h>
#include <TBranch.h>
#include <TRandom3.h>
#include <TLeaf.h>
#include <TString.h>
#include <TSystem.h>
#include <TLegend.h>
#include <TMatrixD.h>
#include <TMatrixDEigen.h>
*/

#include <TFrame.h>
#include <TStyle.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TCanvas.h>

#include <gsl/gsl_fit.h>
#include <gsl/gsl_statistics.h>

#include "AnaInput.h"
#include "MathTools.h"

class Match {

public:

   Match();     
   ~Match();   
   void RepeatRun() ; 
   void RepeatAna() ; 
   void OutputFakePSET() ;
   void ReworkPSET() ;
   void Matching() ;
   void Minimum_GdR( FILE* logfile, vector<vec>& vNom, vector<vec>& vMes ) ; 
   void Rotation( double& x, double& y, double theta ) ;
   void Orth( double& x, double& y, double theta ) ;
   double ReadWaferAlignment( vector<vec>& data, int nSkip = 1  ) ;
   void ReadXCD( vector<vec>& data, int nSkip = 1  ) ;
   void ReadXCDGroup( vector<vec>& data, int nSkip = 1  ) ;
   void ReadDPM( vector<vec>& data, int nSkip = 1  ) ;
   void ReadPset( vector<vec>& data  ) ;
   void ReadPsetAlign( vector<vec>& data  ) ;
   void ReadPsetCFG( vector<vec>& data, int nSkip = 3  ) ;
   void SelectPsetCFG( vector<vec>& data  ) ;
   void ReadPickPlace( vector<vec>& data  ) ;
   void ReadnSort( vector<vec>& data, int axis1, int axis2, int nSkip = 1, bool rise1 = true, bool rise2 = true) ;
   void ReadMask( vector<vec>& data  ) ;
   void ReadStepFAST( vector<vec>& data ) ;  // 0:St,1:Die,2:Tob,3:Valid, 6:X,7:Y,8:dX,9:dY
   void AlignJetStep() ;
   void AlignPickPlace() ;
   void AlignMask() ;
   void AlignDECA() ;
   vector<vec> Sort2D(vector<vec>& data, int axis1, int axis2, bool rise1 = true, bool rise2 = true ) ;


private:

   AnaInput*     Input;

   string cfolder  ;
   string cFileName ;
   string pFileName ;
   string outFileName ;
   string path  ;
   string logFileName ;
   string rootFileName ;

   int debug ;
   double mindX ;

   int die ;
   int site ;
   int tob ;

   //ClassDef(Match, 1);
};

//#if !defined(__CINT__)
//    ClassImp(Match);
#endif

