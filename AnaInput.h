#ifndef AnaInput_H
#define AnaInput_H

#include "TObject.h"
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <math.h>

#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <THStack.h>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TString.h>
#include <TSystem.h>
#include <TLorentzVector.h>

using namespace std;

struct mes{

   double time ;
   double p1 ;
   double p2 ;
   double p3 ;
   double c  ; 
   double np ; // normalized p2
   int    idx  ;
} ;

struct df {

   float x ;
   float y ;
   float s ;
};

typedef vector<double> vec ;
typedef pair<int,int> addr;
//class AnaInput : public TObject {
class AnaInput {

public:

   ~AnaInput();     
   static AnaInput* Instance() ;
  
   void SetDatacard( string datacardInput ) ;

   void GetParameters( string paraName, int* thePara, string cfgFile = "" );
   void GetParameters( string paraName, double* thePara, string cfgFile ="" );
   void GetParameters( string paraName, string* thePara, string cfgFile ="" );
   void GetParameters( string paraName, vector<double>* thePara, string cfgFile = "" );
   void GetParameters( string paraName, vector<float>* thePara, string cfgFile = "" );
   void GetParameters( string paraName, vector<string>* thePara, string cfgFile = "" );
   void GetParameters( string paraName, vector<int>* thePara, string cfgFile = "" );
   
   void ReadFreeCSV( string fileName, vector<vec>& dataV, int nSkipLine = 0 ) ;
   int ReadDPM( string fileName, vector<vec>& dataV ) ;
   void ReadParameter( string paraName, double* thePara, string mapFile );
   void ReadMapCorrection( string mapFile, vector<vec>& data  );
   int ReadMap( string fileName, vector<vec>& dataV, int nSkipLine =0  ) ;

   void ReadSINF( string fileName, vector<vec>& dataV  ) ;
   int ReadKLARF( string fileName, vector<df>& dataV, double lowLimit = 0 ) ;
   //bool Exclusion( df& defect )  ;
   //void SetExclusion( float xL = -1 , float xR = 999999, float yU= 999999, float yD = -1, float threshold = 50 ) ;
   //bool SetExclusion( df& defect, float xL = -1 , float xR = 999999, float yU= 999999, float yD = -1, float threshold = 50 ) ;
   //bool SetExclusion( df& defect, float x_center= 0, float y_center = 0, float radius = 150000, float threshold = 50 ) ;
   bool SetExclusion( df& defect, vector<float>* limit, float threshold = 50 ) ;
   vector<vec> Sort2D( vector<vec>& data, int axis1, int axis2 ) ;


private:

   //AnaInput( string datacardInput = "DataCard.txt" );     
   AnaInput();     

   string datacardfile ;

   static AnaInput* m_Instance;

   float xLimit_l, xLimit_r, yLimit_u, yLimit_d, threshold_k ;


   //ClassDef(AnaInput, 1);

};

//#if !defined(__CINT__)
//    ClassImp(AnaInput);

#endif

