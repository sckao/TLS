#ifndef Reader_H
#define Reader_H

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
#include <math.h>
#include <algorithm>

#include "AnaInput.h"
#include "MathTools.h"

struct shotMes {
 
  shotMes() : all(true) {} 

  bool all ; 
  double x11, x12, x21, x22, x31, x32, x41, x42 ;
  double y11, y12, y21, y22, y31, y32, y41, y42 ;
  double dy21, dy22, dx21, dx22 ;
  double dy41, dy42, dx41, dx42 ;

  void mes() {
       dy21 = y21 - y11 ;
       dy22 = y22 - y12 ;
       dx21 = x21 - x11 - 5.1 - 8.8 ;
       dx22 = x22 - y12 - 5.1 - 8.8 ;
       dy41 = y11 - y41 - 9.8 - 4.2 ;
       dy42 = y12 - y42 - 9.8 - 4.2 ;
       dx41 = x41 - x11 - 5.1 - 8.8 ;
       dx42 = x42 - y12 - 5.1 - 8.8 ;
  }
      
} ;

class Reader {

public:

   Reader( );     
   ~Reader();   
  
   void GetDataXCD() ; 
   void SkimXCD( FILE* logfile, string csvfile ) ; 
   void GetDataFromSINF() ; 
   void GetDataFromMapCorrection() ; 
   void GetDataFromASML() ; 
   void GetDataFromUPROI() ;
   void GetDataFromUPROI( FILE* logfile, int site ) ; 
   void GetDataFromTLSDPM() ;
   void GetDataFromTLSDPM( FILE* logfile, int site ) ;
   void GetDataFromSiteCorrectionMap( vector<vec>& data ) ; 
   void MatchDPM_SiteCorrMap( vector<vec>& data ) ; 
   void MatchDPM_MapCorr( vector<vec>& data ) ; 
   void GetPsetFromPTI() ;
   void LableDie( vector<vec>& data ) ; 
   void MakePlots( vector<vec> data, string plotname ) ;
   void DieRotation(vector<vec>& data ) ;
  
   double Pull(double* xa, double* ya, int size ) ;
   pair<double,double> DieTranslation( vector<vec>& v1, vector<vec>& v2, int ix, int iy ) ;

   void GetDieStatistic() ;

private:

   AnaInput*     Input;

   string cfolder  ;
   string cFileName ;
   string mapFileName ;
   string mapCorrFileName ;
   string hfolder  ;
   string hFileName ;
   string outFileName ;
   string plotType ;
   string plotName0 ;
   string plotName1 ;

   int debug ;
   int nbin ;

   bool Y2X ;

   //ClassDef(Reader, 1);
};


//#if !defined(__CINT__)
//    ClassImp(Reader);
#endif

