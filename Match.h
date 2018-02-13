#ifndef Match_H
#define Match_H

#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

/*
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
*/

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
   void OutputFakePSET( int mode = 0 ) ;
   void Matching() ;
   void Minimum_GdR( FILE* logfile, vector<vec>& vNom, vector<vec>& vMes ) ; 
   void Rotation( double& x, double& y, double theta ) ;
   void ReadXCD( vector<vec>& data, int nSkip = 1  ) ;
   void ReadDPM( vector<vec>& data, int nSkip = 23  ) ;
   void ReadPset( vector<vec>& data  ) ;

private:

   AnaInput*     Input;

   string cfolder  ;
   string cFileName ;
   string pFileName ;
   string outFileName ;
   string path  ;
   string logFileName ;

   int debug ;
   double mindX ;

   int die ;
   int site ;

   //ClassDef(Match, 1);
};

//#if !defined(__CINT__)
//    ClassImp(Match);
#endif

