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

class Module {

public:

   Module();     
   ~Module();
  
   void Analysis() ; 

private:

   AnaInput*     Input;

   string hfolder  ;
   string hFileName ;
   string plotType ;
   string plotName0 ;
   string outFileName ;
   string logFileName ;

   int debug ;

   //ClassDef(Module, 1);
};

//#if !defined(__CINT__)
//    ClassImp(Module);
#endif

