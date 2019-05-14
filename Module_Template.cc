#include "Module.h"
#include "Match.h"

Module::Module( ) {

  Input = AnaInput::Instance() ;

  Input->GetParameters("PlotType",      &plotType ) ;
  Input->GetParameters("Path",          &hfolder ) ;

  Input->GetParameters("HFileName",     &hFileName );
  Input->GetParameters("debug",         &debug );

  Input->GetParameters("PlotName0",     &plotName0 );

  Input->GetParameters("OutputFile",    &outFileName );
  logFileName = hfolder +  outFileName ;

}

Module::~Module(){

  //delete Input ;
    cout<<" done ! "<<endl ;
}


void Module::Analysis() {

   ///// ============ DPM Output ============================================================ 
   //       0    1   2     3     4     5   6    7     8   9   10  11   12     13     14   
   //double x0, y0, xoff, yoff, toff, dx1, dy1, dx2, dy2, xc, yc, dA, diedx, diedy, diedR  ;

   // Logfile 
   FILE* logfile = fopen( logFileName.c_str() ,"w");
   // Read data
   vector<vec> data ;
   Match  *match  = new Match( ) ;
   match->ReadnSort( data, 4, 5, 1 );

   for (size_t i=0; i < data.size(); i++) {
       fprintf( logfile, "%1.0f, %3.0f, %3.0f, %4.6f, %4.6f, %4.6f, %4.6f \n",
                           data[i][1], data[i][4], data[i][5], data[i][11], data[i][12], data[i][13], data[i][14] ) ;

   } 

   fclose( logfile ) ;

  // Book histograms
  //Open Root file to store result
  /*
  string hfile = hfolder + hFileName + ".root" ;
  TFile* rootFile = new TFile( hfile.c_str() , "RECREATE" );
  rootFile->cd() ;

  // Define histograms
  TH1D* h_Yoff = new TH1D("offsetY" , "Y_offset",   20, -500, 500 )  ;
  TH2D* h_r_Toff = new TH2D("r_Toff" , "T_offset",  50, 0, 150,  20, -20, 20 )    ; 
  TGraph2D* g_Xoff = new TGraph2D() ;
  */


 
  // Writing histograms and root files
  /*
  h_Yoff->Write() ;
  h_r_Toff->Write() ;
  g_Xoff->Write() ;
  rootFile->Close() ;
  */
   delete match ;
}

