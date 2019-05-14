#include "Module.h"
#include "Match.h"
#include <gsl/gsl_fit.h>
#include <gsl/gsl_statistics.h>


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

   vector<vec> m1;
   vector<vec> m2;
   for (size_t i=0; i < data.size(); i++) {
       if ( data[i][1] == -1 ) continue ;

       if ( data[i][1] == 0 ) m1.push_back( data[i] ) ;
       if ( data[i][1] == 3 ) m2.push_back( data[i] ) ;

       //fprintf( logfile, "%1.0f, %3.0f, %3.0f, %4.6f, %4.6f, %4.6f, %4.6f \n",
       //                    data[i][1], data[i][4], data[i][5], data[i][11], data[i][12], data[i][13], data[i][14] ) ;
   } 

   // Calculate angle, pitch
   double dx1,dx2,dy1,dy2, ag1, ag2, dA ;
   double px(0), py(0) ;
   bool newColumn =  true ;
   vector<double> pitchX ;
   vector<double> pitchY ;
   for ( size_t i=0; i < m1.size(); i++ ) {

       // calculate angle of die 
       dx1 = m1[i][13] - m2[i][13] ;
       dy1 = m1[i][14] - m2[i][14] ;
       ag1 = atan2(dy1,dx1)*180/3.1415926 ;
       dx2 = m1[i][13] - m2[i][13] + m1[i][11] - m2[i][11] ;
       dy2 = m1[i][14] - m2[i][14] + m1[i][12] - m2[i][12] ;
       ag2 = atan2(dy2,dx2)*180/3.1415926 ;
       dA = ag1 - ag2 ;

       if ( (i+14) < m1.size() ) {
          px = m1[i+14][13] + m1[i+14][11] - m1[i][13] - m1[i][11] ;
          pitchX.push_back( px ) ;
       } else {
          px = 0 ;
       }
       if ( !newColumn ) {

            if ( m1[i][4] != m1[i-1][4] ) { 
               newColumn = true ;
            } else  {
               py = m1[i][14] + m1[i][12] - m1[i-1][14] - m1[i-1][12] ;
               pitchY.push_back( py ) ;
            }
       } else {
            py = 0 ;
            newColumn = false ;
       }

       fprintf( logfile, "%3.0f, %3.0f, %4.6f, %4.6f, %4.6f, %4.6f, %4.6f, %4.6f, %4.6f, %4.6f, %2.6f, %2.6f, %2.6f\n",
                           m1[i][4], m1[i][5], m1[i][11], m1[i][12], m1[i][13], m1[i][14], 
                                               m2[i][11], m2[i][12], m2[i][13], m2[i][14], dA, py, px             ) ;
   }

   const int szx = (int)pitchX.size() ;
   const int szy = (int)pitchY.size() ;
   double pxa[szx], pya[szy] ; 
   cout<<" xsz = "<< szx <<endl ;
   cout<<" ysz = "<< szy <<endl ;
   for ( int i=0; i< szx; i++) {
       pxa[i] = pitchX[i] ;
   }
   for ( int i=0; i< szy; i++) {
       pya[i] = pitchY[i] ;
   }
   double mean_pX = gsl_stats_mean( pxa, 1, szx );
   double std_pX  = gsl_stats_sd(   pxa, 1, szx );
   double mean_pY = gsl_stats_mean( pya, 1, szy );
   double std_pY  = gsl_stats_sd(   pya, 1, szy );

   fprintf( logfile, " , , , , , , , , , , , %2.6f, %2.6f\n", mean_pY, mean_pX ) ; 
   fprintf( logfile, " , , , , , , , , , , , %2.6f, %2.6f\n", std_pY, std_pX ) ; 

   fclose( logfile ) ;
   delete match ;
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
}

