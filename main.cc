#include <iostream> 
#include <vector>
#include <fstream>
#include <string>
#include <algorithm>
#include <stdio.h>
#include <TString.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TMinuit.h>

//#include <pthread.h>
//#include <unistd.h>
//#include "TThread.h"

#include "AnaInput.h"
#include "Reader.h"
#include "Study.h"
#include "Match.h"
#include "Module.h"

using namespace std; 


int main( int argc, const char* argv[] ) { 

  string datacardfile = ( argc > 1 ) ? argv[1] : "DataCard.txt";
  //AnaInput  *Input = new AnaInput( datacardfile );
  AnaInput  *Input = AnaInput::Instance() ;
  Input->SetDatacard( datacardfile ) ;
  

  int module = -1 ;
  Input->GetParameters( "Module", & module ) ;

  if ( module == 0 ) {
     Reader  *reader  = new Reader( ) ;
     //reader->GetDataFromMap();
     reader->GetPsetFromPTI();
     delete reader ;
  }

  if ( module == 2 ) {
     vector<vec> data ;
     Reader  *reader  = new Reader( ) ;
     reader->GetDataFromSiteCorrectionMap( data );
     delete reader ;
  }

  if ( module == 3 ) {
     Reader  *reader  = new Reader( ) ;
     reader->GetDataFromMapCorrection();
     delete reader ;
  }

  if ( module == 4 ) {
     Study  *study  = new Study( ) ;
     study->FromDPM_Matching( );
     delete study ;
  }

  if ( module == 5 ) {
     Module  *mod  = new Module( ) ;
     mod->Analysis( );
     delete mod ;
  }

  if ( module == 6 ) {
     Reader  *reader  = new Reader( ) ;
     reader->GetDataFromSINF( );
     delete reader ;
  }

  if ( module == 7 ) {
     Match  *matcher  = new Match( ) ;
     matcher->AlignJetStep( );
     delete matcher ;
  }

  if ( module == 8 ) {
     Reader  *reader  = new Reader( ) ;
     reader->GetDataFromTLSDPM( );
     delete reader ;
  }

  if ( module == 9 ) {
     Match  *matcher  = new Match( ) ;
     //matcher->AlignPickPlace( );
     //matcher->AlignMask( );
     matcher->AlignDECA( );
     delete matcher ;
  }

  if ( module == 10 ) {
     Match  *matcher  = new Match( ) ;
     //matcher->OutputFakePSET( );
     matcher->ReworkPSET( );
     delete matcher ;
  }

  if ( module == 11 ) {
     Match  *matcher  = new Match( ) ;
     vector<vec> data; 
     matcher->ReadXCD( data, 1 );
     delete matcher ;
  }

  if ( module == 12 ) {
     Match  *matcher  = new Match( ) ;
     matcher->OutputFakePSET( );
     delete matcher ;
  }

  if ( module == 13 ) {
     Reader  *reader  = new Reader( ) ;
     reader->GetDieStatistic();
     delete reader ;
  }

  if (module == 14 ) {
     Module  *mod  = new Module( ) ;
     mod->AnalyzeStepFAST( );
     delete mod ;
  }

  if (module == 15 ) {
     Module  *mod  = new Module( ) ;
     mod->CompareXY( );
     delete mod ;
  }


  delete Input ;
  cout<<" Finished !!!"<<endl ;

  return 0;

}

