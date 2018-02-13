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

  if ( module == 1 ) {
     Reader  *reader  = new Reader( ) ;
     reader->GetDataFromNSX();
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
     //vector<vec> data ;
     //reader->MatchDPM_SiteCorrMap( data );
     //reader->MatchDPM_MapCorr( data );
     reader->GetDataFromMapCorrection();
     delete reader ;
  }

  if ( module == 4 ) {
     Study  *study  = new Study( ) ;
     study->FromDPM_Matching( );
     delete study ;
  }

  if ( module == 5 ) {
     Study  *study  = new Study( ) ;
     study->FromCDM_Matching( );
     delete study ;
  }

  if ( module == 6 ) {
     Reader  *reader  = new Reader( ) ;
     reader->GetDataFromSINF( );
     delete reader ;
  }

  if ( module == 7 ) {
     Reader  *reader  = new Reader( ) ;
     reader->GetDataFromASML( );
     //reader->GetDataFromUPROI( );
     delete reader ;
  }

  if ( module == 8 ) {
     Reader  *reader  = new Reader( ) ;
     reader->GetDataFromTLSDPM( );
     delete reader ;
  }

  if ( module == 9 ) {
     Reader  *reader  = new Reader( ) ;
     reader->GetDataFromUPROI( );
     delete reader ;
  }

  if ( module == 10 ) {
     int mode = 0 ;
     Input->GetParameters("Mode", &mode ) ;
     Match  *matcher  = new Match( ) ;
     //matcher->Matching( );
     matcher->OutputFakePSET( mode );
     //matcher->RepeatRun( );
     delete matcher ;
  }

  if ( module == 11 ) {
     Reader  *reader  = new Reader( ) ;
     reader->GetDataXCD( );
     delete reader ;
  }

  if ( module == 12 ) {
     Match  *matcher  = new Match( ) ;
     matcher->RepeatAna( );
     delete matcher ;
  }

  if ( module == 13 ) {
     Reader  *reader  = new Reader( ) ;
     reader->GetDieStatistic();
     delete reader ;
  }

  delete Input ;
  cout<<" Finished !!!"<<endl ;

  return 0;

}

