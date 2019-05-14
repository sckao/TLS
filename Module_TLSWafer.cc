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

   // Get Data organized in different site , label top and bottom marks 
   vector<vec> s1t;
   vector<vec> s2t;
   vector<vec> s3t;
   vector<vec> s4t;
   vector<vec> s1b;
   vector<vec> s2b;
   vector<vec> s3b;
   vector<vec> s4b;
   int st(0) ;
   for (size_t i=0; i < data.size(); i++) {
       if ( data[i][1] == -1 ) continue ;
       // Get Site Id
       st = SiteID( data[i] ) ;

       if ( data[i][1] == 0 && st == 1 ) s1b.push_back( data[i] ) ;
       if ( data[i][1] == 1 && st == 1 ) s1t.push_back( data[i] ) ;
       if ( data[i][1] == 0 && st == 2 ) s2b.push_back( data[i] ) ;
       if ( data[i][1] == 1 && st == 2 ) s2t.push_back( data[i] ) ;
       if ( data[i][1] == 0 && st == 3 ) s3b.push_back( data[i] ) ;
       if ( data[i][1] == 1 && st == 3 ) s3t.push_back( data[i] ) ;
       if ( data[i][1] == 0 && st == 4 ) s4b.push_back( data[i] ) ;
       if ( data[i][1] == 1 && st == 4 ) s4t.push_back( data[i] ) ;

       //fprintf( logfile, "%1.0f, %3.0f, %3.0f, %4.6f, %4.6f, %4.6f, %4.6f \n",
       //                    data[i][1], data[i][4], data[i][5], data[i][11], data[i][12], data[i][13], data[i][14] ) ;
   } 

   // Calculate angle, pitch
   double dA1(0), dA3(0) ;
   double p12x(0), p12y(0), p14x(0) ;
   vector<double> vA1 ;
   vector<double> vA3 ;
   vector<double> vP12x ;
   vector<double> vP12y ;
   vector<double> vP14x ;
   fprintf( logfile, "x1, y1, x2, y2, x3, y3, x4, y4, dA1, dA3, dx12, dy12, dx14\n") ;
   for ( size_t i=0; i < s1t.size(); i++ ) {

       // calculate angle of die 
       dA1 = dAngle( s1t[i], s1b[i] )  ;
       dA3 = dAngle( s3t[i], s3b[i] )  ;
       p12x = s1t[i][13] + s1t[i][11] - s2t[i][13] - s2t[i][11] ; 
       p12y = s1t[i][14] + s1t[i][12] - s2t[i][14] - s2t[i][12] ; 
       p14x = s1t[i][13] + s1t[i][11] - s4t[i][13] - s4t[i][11] ;

       // Record to vectors
       vA1.push_back(dA1) ;
       vA3.push_back(dA3) ;
       vP12x.push_back(p12x) ;
       vP12y.push_back(p12y) ;
       vP14x.push_back(p14x) ;

       // print out to file 
       fprintf( logfile, "%4.5f, %4.5f, %4.5f, %4.5f, %4.5f, %4.5f, %4.5f, %4.5f, %2.5f, %2.5f, %2.5f, %2.5f, %2.5f\n",
                             s1t[i][11]+ s1t[i][13], s1t[i][12]+ s1t[i][14], 
                             s2t[i][11]+ s2t[i][13], s2t[i][12]+ s2t[i][14], 
                             s3t[i][11]+ s3t[i][13], s3t[i][12]+ s3t[i][14], 
                             s4t[i][11]+ s4t[i][13], s4t[i][12]+ s4t[i][14], 
                             dA1, dA3, p12x, p12y, p14x  ) ;
      
   }

   const int sz = (int)vA1.size() ;
   double AA1[sz], AA3[sz], AP12x[sz], AP12y[sz], AP14x[sz] ;
   for ( int i=0; i< sz; i++) {
       AA1[i] = vA1[i] ;
       AA3[i] = vA3[i] ;
       AP12x[i] = vP12x[i] ;
       AP12y[i] = vP12y[i] ;
       AP14x[i] = vP14x[i] ;
   }
   double mean_a1 = gsl_stats_mean( AA1, 1, sz );
   double std_a1  = gsl_stats_sd(   AA1, 1, sz );
   double mean_a3 = gsl_stats_mean( AA3, 1, sz );
   double std_a3  = gsl_stats_sd(   AA3, 1, sz );
   double mean_p12x = gsl_stats_mean( AP12x, 1, sz );
   double std_p12x  = gsl_stats_sd(   AP12x, 1, sz );
   double mean_p12y = gsl_stats_mean( AP12y, 1, sz );
   double std_p12y  = gsl_stats_sd(   AP12y, 1, sz );
   double mean_p14x = gsl_stats_mean( AP14x, 1, sz );
   double std_p14x  = gsl_stats_sd(   AP14x, 1, sz );

   fprintf( logfile, " , , , , , , , Mean, %2.6f, %2.6f, %2.6f, %2.6f, %2.6f\n", mean_a1, mean_a3, mean_p12x+14, mean_p12y, mean_p14x+14 ) ; 
   fprintf( logfile, " , , , , , , , Stdev, %2.6f, %2.6f, %2.6f, %2.6f, %2.6f\n", std_a1, std_a3, std_p12x, std_p12y, std_p14x ) ; 
   /*  
   fprintf( logfile, " , , , , , , , , , , , %2.6f, %2.6f\n", mean_a1, std_a1 ) ; 
   fprintf( logfile, " , , , , , , , , , , , %2.6f, %2.6f\n", mean_a3, std_a3 ) ; 
   fprintf( logfile, " , , , , , , , , , , , %2.6f, %2.6f\n", mean_p12x+14, std_p12x ) ; 
   fprintf( logfile, " , , , , , , , , , , , %2.6f, %2.6f\n", mean_p14x+14, std_p14x ) ; 
   fprintf( logfile, " , , , , , , , , , , , %2.6f, %2.6f\n", mean_p12y, std_p12y ) ; 
   */

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

int Module::SiteID( vec data ){

     int st = 0 ;
     if ( (int)(data[5])%2 == 0 && (int)(data[4])%2 == 1 ) st  = 1 ;
     if ( (int)(data[5])%2 == 0 && (int)(data[4])%2 == 0 ) st  = 2 ;
     if ( (int)(data[5])%2 == 1 && (int)(data[4])%2 == 1 ) st  = 3 ;
     if ( (int)(data[5])%2 == 1 && (int)(data[4])%2 == 0 ) st  = 4 ;

     return st ;

}
 
//Return angle 
double Module::dAngle( vec pt1, vec pt2, bool degree ) {


       // calculate angle of die 
       double dx1 = pt1[13] - pt2[13] ;
       double dy1 = pt1[14] - pt2[14] ;
       double ag1 = atan2(dy1,dx1) ;
       double dx2 = pt1[13] - pt2[13] + pt1[11] - pt2[11] ;
       double dy2 = pt1[14] - pt2[14] + pt1[12] - pt2[12] ;
       double ag2 = atan2(dy2,dx2) ;
       double dA = (ag1 - ag2)*1000 ;

       if ( degree ) dA = dA*180/3.1415926 ;

       return dA ;
}

