#include "Match.h"
#include <gsl/gsl_fit.h>
#include <gsl/gsl_statistics.h>

Match::Match( ) {

  Input = AnaInput::Instance() ;

  Input->GetParameters("Path",          &path ) ;
  Input->GetParameters("CSVDIR",        &cfolder ) ;
  Input->GetParameters("CSVFile",       &cFileName );
  Input->GetParameters("PSETFile",      &pFileName );
  Input->GetParameters("OutputFile",    &outFileName );
  Input->GetParameters("HFileName",     &rootFileName );
  Input->GetParameters("Die",           &die ) ;
  Input->GetParameters("Site",          &site ) ;
  Input->GetParameters("Tob",           &tob ) ;

  Input->GetParameters("debug",         &debug );

  logFileName = path +  outFileName ;
  mindX = 1 ;
}

Match::~Match(){

  //delete Input ;
    cout<<" done ! "<<endl ;
}

extern int is ;
extern bool ScanSort( vec s1, vec s2) { return (s1[is] < s2[is]) ; }
extern bool ScanSortN( vec s1, vec s2) { return (s1[is] > s2[is]) ; }

void Match::RepeatRun() {

   for (int r=1; r <11; r++ ) {
       for (int i=1; i < 5; i++ ) {
           for (int j=1; j<3; j++  ) {
               site = i ;
               die = j ;
 
               char buff[32] ;
               sprintf(buff, "FF_PSET_S%dD%d_L2_R%d.txt", site, die, r ) ;
	       logFileName = path +  buff ;

               char puff[32] ;
               sprintf(puff, "PSET-S%d-D%d.csv", site, die ) ;
               pFileName = puff ;

               char cuff[32] ;
               sprintf(cuff, "TLS_XCD_L2_%d.csv", r ) ;
               cFileName = cuff ;

               printf(" PSET:%s ,  XCD:%s , Output: %s\n", pFileName.c_str(), cFileName.c_str(), logFileName.c_str() ) ;

               OutputFakePSET() ;

           }
       }
   }
}

// Repeatability study
void Match::RepeatAna() {

   // Record the result
   FILE* logfile = fopen( logFileName.c_str() ,"w");

   int ax1 = 11 ;
   int ax2 = 12 ;
   int nSkip = 1 ;

   // Number of run
   int n_run ; 
   Input->GetParameters("NRun",           &n_run ) ;
   const int rn = n_run ;
   // Number of measurement
   const int sz = 4624 ;
   double x[sz][rn], y[sz][rn] ;

   // Read data from each run for specific parameter
   for (int r=0; r < rn; r++ ) {
 
        // Setup data file name
	char cuff[32] ;
	sprintf(cuff, "PTI_Pretest_Run%d.csv", r+1 ) ;
	cFileName = cuff ;

        // Read Data for each run
        vector<vec> data ; 
	string csvfile = cfolder + cFileName ;
	printf(" Data file = %s \n", csvfile.c_str() ) ;
	Input->ReadFreeCSV( csvfile, data, nSkip ) ;

        //fprintf(logfile, " == Run %d == \n", r+1 ) ;

        for (size_t i=0; i< sz; i++) { 
 
            x[i][r] = data[i][ax1] ;            
            y[i][r] = data[i][ax2] ; 
           
            //fprintf(logfile, " %04.6f  % 04.6f  %04.6f  %04.6f  %04.6f  %04.6f\n", 
            //        psetV[i][ax1], psetV[i][ax2], xcdV[i][ax1], xcdV[i][ax2], dpmV[i][ax1], dpmV[i][ax2] ) ;
            //fprintf(logfile, " %04.6f  % 04.6f \n", data[i][ax1], data[i][ax2] ) ;

        }
   }

   fprintf(logfile, "\n\n === statistic ==== \n") ;
   
   for ( int i=0; i< sz; i++) {

       double xa[rn],ya[rn] ; 
       for ( int j=0; j< rn; j++) {
           xa[j] = x[i][j] ;
           ya[j] = y[i][j] ;
       }

       double mn_x  = gsl_stats_mean( xa, 1, rn );
       double st_x  = gsl_stats_sd(xa, 1, rn );
       double mn_y  = gsl_stats_mean( ya, 1, rn );
       double st_y  = gsl_stats_sd(ya, 1, rn );

       fprintf(logfile, " %04.6f, %04.6f, %04.6f, %04.6f \n",
                           mn_x,  st_x,  mn_y,  st_y  ) ;
   }
   
   fclose( logfile );

}


// Create FakePSET using FFLY RRF output
void Match::OutputFakePSET() {

    int nAM = 2 ;
    Input->GetParameters("NAlignMark",  &nAM ) ;

    vector<vec> ffV ;
    FILE* logfile = fopen( logFileName.c_str() ,"w");
    ReadnSort( ffV, 4, 5, 1 );


    double x(0), y(0), epx(0), epy(0) ;
    double dx = 0 ;
    double dy = 0 ;
    fprintf(logfile, "INDEX, X, Y, AlignType, EXPOSURE_X, EXPOSURE_Y, SiteRelativeX, SiteRelativeY\n") ;    
    int k = 0 ;
    for (int i=0; i< (int)ffV.size(); i++ ) {
          
	  if ( ffV[i][1] < 0 ) continue ;
	  if ( tob > -1 &&  ffV[i][1] != tob ) continue ;
	  int x_i = (int) ffV[i][14] ;
	  x = ( x_i > 0 ) ? (double) x_i : (double) (x_i -1) ;
          int y_i = (int) (-1*ffV[i][13]) ; 
	  y = ( y_i > 0 ) ? (double) y_i : (double) (y_i -1) ;
          
          epx =  x - dx; 
          epy =  y - dy; 
        // PSET CFG Format
        fprintf(logfile, "%d, % 04.6f, % 04.6f, M, % 04.6f, % 04.6f, % 04.6f, % 04.6f\n", 
                     k+nAM+1,       x,       y,        epx,     epy,      dx,     dy ) ;
        // Debug Format 
        //fprintf(logfile, " (%.0f, %.0f) %04.6f  % 04.6f, (%.0f, %.0f)  %04.6f  %04.6f,  %04.6f  %04.6f\n", 
        //                 psetV[i][5], psetV[i][6], x, y, ffV[i][1], ffV[i][2],ffV[i][ax1], ffV[i][ax2], dx, dy ) ;
	k++ ;
    } 
    fclose( logfile );

}

void Match::ReworkPSET() {

    string csvfile = cfolder + pFileName ;
    printf(" PSET file = %s \n", csvfile.c_str() ) ;

    vector<vec> data ;
    //Input->ReadFreeCSV( csvfile, data, 1 ) ;
    //
    // 0: cl, 1: dX, 2:dY, 3:X, 4:Y, 5:st, 6:de, 7:X_shot, 8:Y_shot 
    ReadPsetCFG( data, 3 ) ;

    // Record the result
    FILE* logfile = fopen( logFileName.c_str() ,"w");

    fprintf(logfile, "INDEX, X, Y, AlignType, EXPOSURE_X, EXPOSURE_Y, SiteRelativeX, SiteRelativeY\n") ;  

    for (int i=0; i< (int)data.size(); i++ ) {

	/*     
        if (i<3) {
           fprintf(logfile, "%d, % 04.6f, % 04.6f, A, % 04.6f, % 04.6f, % 04.6f, % 04.6f\n", 
                             i+1, data[i][1], data[i][2], data[i][4],  data[i][5], data[i][6], data[i][7] ) ;
	} else {
           fprintf(logfile, "%d, % 04.6f, % 04.6f, M, % 04.6f, % 04.6f, % 04.6f, % 04.6f\n", 
                             i+1, data[i][1], data[i][2], data[i][4],  data[i][5], data[i][6], data[i][7] ) ;
        }
	*/
        fprintf(logfile, "%d, % 04.6f, % 04.6f, M, % 04.6f, % 04.6f, % 04.6f, % 04.6f\n", 
                           i+6, data[i][3], data[i][4], data[i][7],  data[i][8], data[i][1], data[i][2] ) ;

    } 
    fclose( logfile );

}

void Match::Matching() {

    vector<vec> xcdV ;
    ReadDPM( xcdV ) ;

    vector<vec> psetV ;
    ReadPsetCFG( psetV ) ;

    // Record the result
    FILE* logfile = fopen( logFileName.c_str() ,"w");

    Minimum_GdR( logfile, psetV, xcdV ) ;
    fclose( logfile );

}


void Match::AlignJetStep() {

     // log file
     FILE* logfile = fopen( logFileName.c_str() ,"w");

     // Get FireFly Wafer alignment
     vector<vec> ffV ;
     double angle = ReadWaferAlignment( ffV, 1 );
     double dxc = (ffV[0][1] - ffV[1][1])    ;  
     double dyc = (ffV[0][2] - ffV[1][2])    ;  
     //double angle = atan(dyc/dxc) ; 
     double x1 = ffV[0][1] ;
     double y1 = ffV[0][2] ;
     double x2 = ffV[1][1] ;
     double y2 = ffV[1][2] ;
     Rotation( x1, y1, -1*angle );
     Rotation( x2, y2, -1*angle );
     double xc = (x1 + x2) / 2.   ;  
     double yc = (y1 + y2) / 2.   ;  
     fprintf( logfile, " FireFly \n" ) ;     
     fprintf( logfile, "Point1: %04.6f, %04.6f ->  %04.6f, %04.6f \n", ffV[0][1], ffV[0][2], x1, y1  ) ;     
     fprintf( logfile, "Point2: %04.6f, %04.6f ->  %04.6f, %04.6f \n", ffV[1][1], ffV[1][2], x2, y2  ) ;     
     fprintf( logfile, "Center :  %04.6f, %04.6f \n\n",          xc,        yc  ) ;     
     fprintf( logfile, "dX = %04.6f , dY = %04.6f \n",       dxc,       dyc  ) ;     
     fprintf( logfile, " %f uRad, %f deg\n",  atan(dyc/dxc)*1000000, atan(dyc/dxc)*180/3.1415926 ) ;    

     // Get JetStep Wafer alignment
     vector<vec> jsA ;
     ReadPsetAlign(jsA ) ;
     double jxc = (jsA[0][1] + jsA[1][1]) / 2.   ;  
     double jyc = (jsA[0][2] + jsA[1][2]) / 2.   ;  
     fprintf( logfile, "\n JetStep Wafer Alignment\n" ) ;     
     fprintf( logfile, "Point1: %04.6f %04.6f \n", jsA[0][1], jsA[0][2]  ) ;     
     fprintf( logfile, "Point2: %04.6f %04.6f \n", jsA[1][1], jsA[1][2]  ) ;     
     fprintf( logfile, "Center: %04.6f  %04.6f \n\n",       jxc,        jyc  ) ;
     fprintf( logfile, "dX = %04.6f , dY = %04.6f \n",      jsA[0][1] - jsA[1][1],   jsA[0][2] - jsA[1][2]  ) ;     
  
     // Get JetStep Die positions - Reading pset.CFG file  
     vector<vec> jsV ;
     ReadPsetCFG(jsV ) ;
     SelectPsetCFG( jsV ) ;
     // Get JetStep Corrections - Reading pset file
     // Only if JetStep provide it
     /*
     vector<vec> jsP ;
     ReadPset(jsP ) ;
     */

     // JetStep - FFLY shift -> Using the center of two alignment marks
     double dx_jf = jxc -xc ;
     double dy_jf = jyc -yc ;
     fprintf( logfile, "Center Shift\n" ) ;     
     fprintf( logfile, " %04.6f  %04.6f \n\n",       dx_jf,        dy_jf ) ;     

     // Get FFLY Die Alignment
     vector<vec> ffD0 ;
     ReadDPM( ffD0  ) ;

     // Booking root file and histograms 
     //Open Root file to store result
     string hfile = path + rootFileName + ".root" ;
     TFile* rootFile = new TFile( hfile.c_str() , "RECREATE" );
     rootFile->cd() ;
     // Book histogram
     TGraph2D* g_dX = new TGraph2D()  ;
     TGraph2D* g_dY = new TGraph2D()  ;
     g_dX->SetName("dX") ;
     g_dY->SetName("dY") ;

     double xr,yr ;
     const int sz = (int)jsV.size(); 
     double xa[sz], ya[sz] ;
     fprintf(logfile, "       CFG_X,      CFG_Y,          FF_X,       FF_Y,      CorrFF_X,   CorrFF_Y,        dX,        dY\n") ;
     for (size_t i=0 ; i < ffD0.size(); i++ ) {
         xr = ffD0[i][3] ;
         yr = ffD0[i][4] ;
         Rotation( xr, yr, -1*angle );
         xr = xr + dx_jf ;
         yr = yr + dy_jf ;
         fprintf(logfile, " %4.6f, %4.6f,   %4.6f, %4.6f,   %4.6f, %4.6f,  %4.6f, %4.6f\n", 
                        jsV[i][3], jsV[i][4], ffD0[i][3], ffD0[i][4], xr, yr, jsV[i][3]-xr, jsV[i][4]-yr ) ;
         //fprintf(logfile, " %4.6f, %4.6f,   %4.6f, %4.6f,   %4.6f, %4.6f,  %4.6f, %4.6f,  %4.6f, %4.6f,  %4.6f, %4.6f\n", 
         //               jsV[i][3], jsV[i][4], ffD0[i][3], ffD0[i][4], xr, yr, jsV[i][3]-xr, jsV[i][4]-yr, 
         //                                                         jsP[i][0], jsP[i][1], jsP[i][2], jsP[i][3] ) ;

         xa[i] = jsV[i][3]-xr ; 
         ya[i] = jsV[i][4]-yr ; 

         g_dX->SetPoint( i , jsV[i][3], jsV[i][4], jsV[i][3]-xr );
         g_dY->SetPoint( i , jsV[i][3], jsV[i][4], jsV[i][4]-yr );

     }

     double mean_dX = gsl_stats_mean( xa, 1, sz );
     double std_dX  = gsl_stats_sd(   xa, 1, sz );
     double mean_dY = gsl_stats_mean( ya, 1, sz );
     double std_dY  = gsl_stats_sd(   ya, 1, sz );
     fprintf(logfile, " %d, %4.6f,  %4.6f \n ", sz, xa[sz-1], ya[sz-1] ) ;
     fprintf(logfile, " %4.6f +/- %4.6f \n", mean_dX, std_dX ) ;
     fprintf(logfile, " %4.6f +/- %4.6f \n", mean_dY, std_dY ) ;

     fclose( logfile );


     // Plot dX &  dY and save the root file
     /*
     TCanvas* c0 = new TCanvas("c0","", 0,0, 800, 600);
     c0->SetFillColor(10);
     c0->SetFillColor(10);
     //gPad->SetGridx();
     c0->Clear() ;

     gStyle->SetPalette(1);
     gStyle->SetNumberContours(10);
     g_dX->Draw("surf1z");
     c0->Update();
     TString gPlotname1 = path + "dX.png"  ;
     c0->Print( gPlotname1 ) ;

     c0->Clear() ;
     g_dY->Draw("surf1z");
     c0->Update();
     gPlotname1 = path + "dY.png"  ;
     c0->Print( gPlotname1 ) ;

     delete c0; 
     */

     g_dX->Write() ;
     g_dY->Write() ;
     rootFile->Close() ;

}

void Match::AlignPickPlace() {

     // log file
     FILE* logfile = fopen( logFileName.c_str() ,"w");

     // Get FireFly Wafer alignment
     
     vector<vec> ffV ;
     double angle = ReadWaferAlignment( ffV, 1 );
     double xc = (ffV[0][1] + ffV[1][1]) / 2.   ;  
     double yc = (ffV[0][2] + ffV[1][2]) / 2.   ;  
     double dxc = (ffV[0][1] - ffV[1][1])    ;  
     double dyc = (ffV[0][2] - ffV[1][2])    ;
     //double angle = atan(dyc/dxc) ; 
     double x1 = ffV[0][1] ;
     double y1 = ffV[0][2] ;
     double x2 = ffV[1][1] ;
     double y2 = ffV[1][2] ;
     Rotation( x1, y1, -1*angle );
     Rotation( x2, y2, -1*angle );
 
     fprintf( logfile, " FireFly \n" ) ;     
     fprintf( logfile, " %04.6f %04.6f \n", ffV[0][1], ffV[0][2]  ) ;     
     fprintf( logfile, " %04.6f %04.6f \n", ffV[1][1], ffV[1][2]  ) ;     
     fprintf( logfile, " %04.6f %04.6f \n",        xc,        yc  ) ;     
     fprintf( logfile, " %04.6f %04.6f  %f uRad, %f deg\n",     dxc,     dyc, atan(dyc/dxc)*1000000, atan(dyc/dxc)*180/3.1415926 ) ;    
     fprintf( logfile, " %04.6f %04.6f \n",          x1,        y1  ) ;     
     fprintf( logfile, " %04.6f %04.6f \n\n",        x2,        y2  ) ;     
      

     // Get FFLY Die Alignment
     vector<vec> ffD0 ;
     ReadDPM( ffD0  ) ;
     vector<vec> ffD = Sort2D( ffD0, 3, 4 ) ;

     // Get Pick & Place Design Position
     vector<vec> ppD ;
     ReadPickPlace( ppD ) ;

     const int sz = (int)ppD.size(); 
     double xa[sz], ya[sz] ;
     double xr,yr ;
     int j = 0 ;
     for (int i= 0; i< sz; i++ ) { 
         if ( i%18 ==0 && i!= 0 ) j+=18 ;
         //fprintf( logfile, "(%d) %4.6f, %4.6f, (%d) %4.6f, %4.6f\n", i, ppD[i][4], ppD[i][5], j, ffD[j][3],   ffD[j][4] ) ;
         //fprintf( logfile, "(%d) %4.6f, %4.6f, (%d) %4.6f, %4.6f\n", i, ppD[i][4], ppD[i][5], j+18, ffD[j+18][3], ffD[j+18][4] ) ;
         
         // For mark 1
         //fprintf( logfile, "%4.6f, %4.6f, %4.6f, %4.6f, %4.6f, %4.6f\n", ppD[i][4], ppD[i][5], ffD[j][3], ffD[j][4], ffD[j][5], ffD[j][6] ) ;
         // For mark 2
         //fprintf( logfile, "%4.6f, %4.6f, %4.6f, %4.6f, %4.6f, %4.6f, %.0f\n", ppD[i][4], ppD[i][5], ffD[j+18][3], ffD[j+18][4], ffD[j+18][5], ffD[j+18][6], ffD[j+18][0] ) ;
         //xa[i] = ffD[j+18][5] ;
         //ya[i] = ffD[j+18][6] ;
         //xa[i] = ffD[j+18][3] - ppD[i][4] ;
         //ya[i] = ffD[j+18][4] - ppD[i][5] ;
         xr = ffD[j+18][3] ;
         yr = ffD[j+18][4] ;
         Rotation( xr, yr, -1*angle ) ;
         xa[i] = xr - ppD[i][4] ;
         ya[i] = yr - ppD[i][5] ;
  
         fprintf( logfile, "%4.6f, %4.6f, %4.6f, %4.6f, %4.6f, %4.6f, %4.6f, %4.6f, %4.6f, %4.6f, %.0f\n", 
                          ppD[i][4], ppD[i][5], ffD[j+18][3], ffD[j+18][4], ffD[j+18][5], ffD[j+18][6], xr, yr, xa[i], ya[i], ffD[j+18][0] ) ;
         j++ ;
     }
     double mean_dX = gsl_stats_mean( xa, 1, sz );
     double std_dX  = gsl_stats_sd(   xa, 1, sz );
     double mean_dY = gsl_stats_mean( ya, 1, sz );
     double std_dY  = gsl_stats_sd(   ya, 1, sz );
     fprintf(logfile, " %4.6f +/- %4.6f \n", mean_dX, std_dX ) ;
     fprintf(logfile, " %4.6f +/- %4.6f \n", mean_dY, std_dY ) ;
  

     fclose( logfile );
}

// In this case, FFLY X,Y is the same as Mask Coordiante X, Y
void Match::AlignMask() {

     // log file
     FILE* logfile = fopen( logFileName.c_str() ,"w");

     // Get FireFly Wafer alignment
     vector<vec> ffV ;
     double angle = ReadWaferAlignment( ffV, 1 );
     double xc = (ffV[0][2] + ffV[1][2]) / 2.   ;  
     double yc = (ffV[0][1] + ffV[1][1]) / 2.   ;  
     double dxc = (ffV[0][2] - ffV[1][2])    ;  
     double dyc = (ffV[0][1] - ffV[1][1])    ;
     //double angle = atan(dyc/dxc) ; 
     double x1 = ffV[0][2] ;
     double y1 = ffV[0][1] ;
     double x2 = ffV[1][2] ;
     double y2 = ffV[1][1] ;
     Rotation( x1, y1, -1*angle );
     Rotation( x2, y2, -1*angle );
 
     fprintf( logfile, " FireFly \n" ) ;     
     fprintf( logfile, " %04.6f %04.6f \n", ffV[0][2], ffV[0][1]  ) ;     
     fprintf( logfile, " %04.6f %04.6f \n", ffV[1][2], ffV[1][1]  ) ;     
     fprintf( logfile, " %04.6f %04.6f \n",        xc,        yc  ) ;     
     fprintf( logfile, " %04.6f %04.6f  %f uRad, %f deg\n",     dxc,     dyc, atan(dyc/dxc)*1000000, atan(dyc/dxc)*180/3.1415926 ) ;    
     fprintf( logfile, " %04.6f %04.6f \n",          x1,        y1  ) ;     
     fprintf( logfile, " %04.6f %04.6f \n\n",        x2,        y2  ) ;     

     // Mask Wafer Alignment
     double m1[2] = { -3.65, 45 } ;       
     double m2[2] = { 81.15, 45 } ;       
     double mcx = (m1[0] + m2[0]) /2 ;
     double mcy = (m1[1] + m2[1]) /2 ;
     // X, Y Shift
     double dx_mf = mcx - ((x1 + x2)/2) ;
     double dy_mf = mcy - ((y1 + y2)/2) ;
     fprintf( logfile, " center shift \n") ; 
     fprintf( logfile, " %04.6f %04.6f \n\n",       dx_mf,        dy_mf  ) ;     


     // Get FFLY Die Alignment
     vector<vec> ffD0 ;
     ReadDPM( ffD0  ) ;
     vector<vec> ffD = Sort2D( ffD0, 4, 3 ) ;

     // Get Pick & Place Design Position
     vector<vec> ppD ;
     ReadMask( ppD ) ;

     const int sz = (int)ppD.size(); 
     double xa[sz], ya[sz] ;
     double xr,yr ;
     for (int i= 0; i< sz; i++ ) { 
         xr = ffD[i][4] ;
         yr = ffD[i][3] ;
         Rotation( xr, yr, -1*angle ) ;
         xr = xr + dx_mf ;
         yr = yr + dy_mf ;
         xa[i] = xr - ppD[i][3] ;
         ya[i] = yr - ppD[i][4] ;
  
         fprintf( logfile, "%4.6f, %4.6f, %4.6f, %4.6f, %4.6f, %4.6f, %4.6f, %4.6f, %4.6f, %4.6f \n", 
                          ppD[i][3], ppD[i][4], ffD[i][4], ffD[i][3], xr, yr, xa[i], ya[i], ffD[i][4]-ppD[i][3], ffD[i][3]-ppD[i][4] ) ;
         //fprintf( logfile, "%4.6f, %4.6f, %4.6f, %4.6f, %4.6f, %4.6f, %4.6f, %4.6f \n", 
         //                 ppD[i][3], ppD[i][4], ffD[i][4], ffD[i][3], xr, yr, xa[i], ya[i]  ) ;
         // MM Pset                
         //fprintf( logfile, "%4.6f, %4.6f, %4.6f, %4.6f\n", 
         //                 ppD[i][3], ppD[i][4], ffD[i][4] - ppD[i][3] , ffD[i][3] - ppD[i][4] ) ;
     }
     //double mean_dX = gsl_stats_mean( xa, 1, sz );
     //double std_dX  = gsl_stats_sd(   xa, 1, sz );
     //double mean_dY = gsl_stats_mean( ya, 1, sz );
     //double std_dY  = gsl_stats_sd(   ya, 1, sz );
     //fprintf(logfile, " %4.6f +/- %4.6f \n", mean_dX, std_dX ) ;
     //fprintf(logfile, " %4.6f +/- %4.6f \n", mean_dY, std_dY ) ;
  

     fclose( logfile );
}

void Match::AlignDECA() {

     // log file
     FILE* logfile = fopen( logFileName.c_str() ,"w");

     // Read Data 
     vector<vec> ppD ;
     ReadnSort( ppD, 2,3,1 ) ;

     // Calculate Rotation
     //double x1 = ppD[3][4] + ppD[3][6] ;
     //double y1 = ppD[3][5] + ppD[3][7];
     //double x2 = ppD[3494][4] + ppD[3494][6] ;
     //double y2 = ppD[3494][5] + ppD[3494][7];
     double x1 = ppD[15][4] + ppD[15][6] ;
     double y1 = ppD[15][5] + ppD[15][7];
     double x2 = ppD[3482][4] + ppD[3482][6] ;
     double y2 = ppD[3482][5] + ppD[3482][7];
     double dxc = ( x1 - x2 )    ;  
     double dyc = ( y1 - y2 )    ;  
     double angle = atan(dyc/dxc) ; 
     Rotation( x1, y1, -1*angle );
     Rotation( x2, y2, -1*angle );
     fprintf( logfile, " Rotation from two alignment marks\n" ) ;     
     fprintf( logfile, "  %f uRad, %f deg\n", atan(dyc/dxc)*1000000, atan(dyc/dxc)*180/3.1415926 ) ;    
     fprintf( logfile, " %04.6f %04.6f => ", ppD[15][4] + ppD[15][6], ppD[15][5] + ppD[15][7]  ) ;     
     fprintf( logfile, " %04.6f %04.6f \n",          x1,        y1  ) ;     
     fprintf( logfile, " %04.6f %04.6f => ", ppD[3482][4] + ppD[3482][6], ppD[3482][5] + ppD[3482][7]  ) ;     
     fprintf( logfile, " %04.6f %04.6f \n",          x2,        y2  ) ;     
     // Calculate Shift
     // Center shift frim two wafer alignment marks
     double xc = (x1 + x2) / 2. ;
     double yc = (y1 + y2) / 2. ;
     double xc0 = (ppD[15][0] + ppD[3482][0])  / 2000 ; 
     double yc0 = (ppD[15][1] + ppD[3482][1])  / 2000 ; 
     double dx = xc - xc0 ;
     double dy = yc - yc0 ;
     fprintf( logfile, "\nCenter Shift from two alignment marks\n" ) ;     
     fprintf( logfile, "Measurement Center:  %04.6f %04.6f \n", xc, yc) ;
     fprintf( logfile, "Nominal Center:      %04.6f %04.6f \n", xc0, yc0) ;
     fprintf( logfile, "Differece:           %04.6f %04.6f \n", dx, dy ) ;

     const int sz = (int) ppD.size() ;     
     double xa[sz], ya[sz], xm1[sz], ym1[sz], theta[sz] ;
     //double xm2[sz], ym2[sz] ;
     // Calculate Angle
     double vx1,vy1, vx2, vy2 ;
     MathTools *math = new MathTools() ;
     for (int i= 0; i< sz; i++ ) { 
         // record original measurement data
         xm1[i] = ppD[i][4] + ppD[i][6] ;
         ym1[i] = ppD[i][5] + ppD[i][7] ;
         //xm2[i] = ppD[i][8] + ppD[i][10] ;
         //ym2[i] = ppD[i][9] + ppD[i][11] ;
         // Nominal data 
         xa[i] = ppD[i][0]/1000. ;
         ya[i] = ppD[i][1]/1000. ;

         // ffly expected - vector1 from two ffly expected postions
         vx1 = ppD[i][6] - ppD[i][10] ;
         vy1 = ppD[i][7] - ppD[i][11] ;
         // ffly measured - vector2 from two ffly measured positions
         vx2 = ppD[i][4] + ppD[i][6] - ppD[i][8] - ppD[i][10] ; 
         vy2 = ppD[i][5] + ppD[i][7] - ppD[i][9] - ppD[i][11] ;     
         // calculated theta 
         theta[i] = math->angleAB( vx1, vy1, vx2, vy2 ) ;
     }

     double mean_dX0 = gsl_stats_mean( xa, 1, sz );
     double std_dX0  = gsl_stats_sd(   xa, 1, sz );
     double mean_dY0 = gsl_stats_mean( ya, 1, sz );
     double std_dY0  = gsl_stats_sd(   ya, 1, sz );
     double mean_mdx = gsl_stats_mean( xm1, 1, sz );
     double mean_mdy = gsl_stats_mean( ym1, 1, sz );
     //dx = mean_mdx - mean_dX0 ;
     //dy = mean_mdy - mean_dY0 ;
     fprintf(logfile, "\nAverage XY from all nominal positions \n" ) ;
     fprintf(logfile, " %4.6f +/- %4.6f \n", mean_dX0, std_dX0 ) ;
     fprintf(logfile, " %4.6f +/- %4.6f \n", mean_dY0, std_dY0 ) ;
     fprintf(logfile, "\nAverage XY from all measurement positions \n" ) ;
     fprintf(logfile, " %4.6f  \n", mean_mdx ) ;
     fprintf(logfile, " %4.6f  \n", mean_mdy ) ;
     fprintf( logfile, "Differece:   %04.6f %04.6f \n", dx, dy ) ;

    
     fprintf( logfile, "\nDie Measurement List\n" ) ;     
     fprintf( logfile, "\nNominalX, NominalY, DieX, DieY, X1, Y1, X2, Y2, dX1, dY1, Theta\n" ) ;     
     double xr1[sz], yr1[sz], xr2[sz], yr2[sz];
     //double  xb[sz], yb[sz] ;
     double sx2(0), sy2(0) ;
     double xf1,yf1, xf2, yf2 ; 
     double xd[sz], yd[sz] ;
     double xd2(0), yd2(0) ;


     for (int i= 0; i< sz; i++ ) { 
         
         // for calculating R-squre 
         sx2 += (xm1[i] - mean_mdx)*(xm1[i] - mean_mdx) ;
         sy2 += (ym1[i] - mean_mdy)*(ym1[i] - mean_mdy) ;
 
         // Rotating measurement data      
         x1 = ppD[i][4] + ppD[i][6] ;
         y1 = ppD[i][5] + ppD[i][7] ;
         x2 = ppD[i][8] + ppD[i][10] ;
         y2 = ppD[i][9] + ppD[i][11] ;

	 Rotation( x1, y1, -1*angle );
	 Rotation( x2, y2, -1*angle );

         xr1[i] = x1 ;
         yr1[i] = y1 ;
         xr2[i] = x2 ;
         yr2[i] = y2 ;
         //xb[i] = xm1[i] - xa[i] ; 
         //yb[i] = ym1[i] - ya[i] ; 

         // get the final coordinates with rotation and shift 
         xf1 = xr1[i] - dx ;
         yf1 = yr1[i] - dy ;
         xf2 = xr2[i] - dx ;
         yf2 = yr2[i] - dy ;

         xd[i] = xf1 - ppD[i][0]/1000. ;
         yd[i] = yf1 - ppD[i][1]/1000. ;

         // for calculating R-squre
         xd2 += ((xm1[i] - xf1)*(xm1[i] - xf1))  ;
         yd2 += ((ym1[i] - yf1)*(ym1[i] - yf1)) ;

         fprintf( logfile, "%4.4f, %4.4f, %3.0f, %3.0f, %4.6f, %4.6f, %4.6f, %4.6f, %4.6f, %4.6f, %1.6f\n", 
                          ppD[i][0]/1000., ppD[i][1]/1000., ppD[i][2], ppD[i][3], xf1, yf1,  xf2, yf2, xd[i], yd[i], theta[i]*180/3.1415926 ) ;
     } 
     //double mean_dX = gsl_stats_mean( xr1, 1, sz );
     //double std_dX  = gsl_stats_sd(   xr1, 1, sz );
     //double mean_dY = gsl_stats_mean( yr1, 1, sz );
     //double std_dY  = gsl_stats_sd(   yr1, 1, sz );
     //fprintf(logfile, "\nAverage XY from measurement \n" ) ;
     //fprintf(logfile, "meanX: %4.6f +/- %4.6f \n", mean_dX, std_dX ) ;
     //fprintf(logfile, "meanY: %4.6f +/- %4.6f \n", mean_dY, std_dY ) ;

     double RsqureX = 1 - (xd2/sx2) ;
     double RsqureY = 1 - (yd2/sy2) ;
     double mean_dX = gsl_stats_mean( xd, 1, sz );
     double std_dX  = gsl_stats_sd(   xd, 1, sz );
     double mean_dY = gsl_stats_mean( yd, 1, sz );
     double std_dY  = gsl_stats_sd(   yd, 1, sz );
     fprintf(logfile, "\n Average XY from measurement \n" ) ;
     fprintf(logfile, "meanX: %4.6f +/- %4.6f \n", mean_dX, std_dX ) ;
     fprintf(logfile, "meanY: %4.6f +/- %4.6f \n", mean_dY, std_dY ) ;
     fprintf(logfile, "R-squared X: %4.6f, Y: %4.6f \n", RsqureX, RsqureY ) ;

     
     fclose( logfile );
}

// Output id, x, y, dx, dy
double Match::ReadWaferAlignment( vector<vec>& data1, int nSkip  ) {

    data1.clear() ;

    vector<vec> data ;
    data.clear() ;

    string aFileName ;
    Input->GetParameters("WaferAlign",       &aFileName );
    string csvfile = cfolder + aFileName ;
    printf(" WaferAlignment file = %s \n", csvfile.c_str() ) ;
    Input->ReadFreeCSV( csvfile, data, nSkip ) ;

    // Remove alignment mode -1 and miss-measured sites
    // printf(" size = %d \n", (int)data.size() ) ;
    // printf(" row size 1   = %d \n", (int)data[0].size() ) ;
    // printf(" row size end = %d \n", (int)data[ (int)data.size() -1 ].size() ) ;
    for ( std::vector<vec>::iterator it = data.begin() ; it != data.end();  ) {
        //if ( it->size() < 10 || (*it)[3] < 0 ) { 
        if ( (*it)[1] < 0 ) { 
           data.erase( it ) ; 
           //continue ;
        } else {
           it++ ;
        }
    }
    
    string debugFile = path + "debug_wafer.csv" ;
    FILE* logfile = fopen( debugFile.c_str() ,"w");

    is = 6 ;
    if ( data.size() > 1 ) {
        sort( data.begin(), data.end(), ScanSort );
    }

    vec mdata ;
    for (size_t j=0; j < data.size(); j++ ) {

        mdata.push_back( double(j+1) ) ;
	mdata.push_back( data[j][14] + data[j][12]  ) ;
	mdata.push_back( -1*data[j][13] - data[j][11]  ) ;
	mdata.push_back( data[j][12] ) ;  // dx - new TLS version
	mdata.push_back( -1*data[j][11] ) ;  // dy - new TLS version

	data1.push_back( mdata ) ;
        mdata.clear() ;

	fprintf(logfile, "%d, % 04.6f, % 04.6f, % 04.6f, % 04.6f \n", 
                        (int)j+1, data[j][14], -1*data[j][13], data[j][12], -1*data[j][11]   ) ;

    } 

    double dxc = (data1[0][1] - data1[1][1])    ;  
    double dyc = (data1[0][2] - data1[1][2])    ;  
    double angle = atan(dyc/dxc) ; 

    double x0 = data1[0][1] ;
    double y0 = data1[0][2] ;
    double x1 = data1[1][1] ;
    double y1 = data1[1][2] ;
    printf(" Before rotation (%.6f,%.6f) - (%.6f,%.6f) \n", x0, y0, x1, y1 ) ;
    Rotation( x0, y0, -1*angle ) ;
    Rotation( x1, y1, -1*angle ) ;
    printf(" After rotation  (%.6f,%.6f) - (%.6f,%.6f) \n", x0, y0, x1, y1 ) ;

    printf(" Data Size = %d entries - %d items\n", (int)data1.size(), (int)data1[0].size() ) ;
    printf(" Angle = %.6f  Rad \n", angle ) ;
    fclose( logfile ) ;

    return angle ;
}

// Output id, site, die, x, y, dx, dy
void Match::ReadDPM( vector<vec>& data1, int nSkip  ) {

    data1.clear() ;

    vector<vec> data ;
    data.clear() ;

    string csvfile = cfolder + cFileName ;
    printf(" DPM file = %s \n", csvfile.c_str() ) ;
    Input->ReadFreeCSV( csvfile, data, nSkip ) ;

    // Remove alignment mode -1 and miss-measured sites
    // printf(" size = %d \n", (int)data.size() ) ;
    // printf(" row size 1   = %d \n", (int)data[0].size() ) ;
    // printf(" row size end = %d \n", (int)data[ (int)data.size() -1 ].size() ) ;
    for ( std::vector<vec>::iterator it = data.begin() ; it != data.end();  ) {
        //if ( it->size() < 10 || (*it)[3] < 0 ) { 
        if ( (*it)[1] < 0 ) { 
           data.erase( it ) ; 
           //continue ;
        } else {
           it++ ;
        }
    }
    
    string debugFile = path + "debug.csv" ;
    FILE* logfile = fopen( debugFile.c_str() ,"w");

    int axis1 = 6 ;
    int axis2 = ( axis1 == 7 ) ? 6 : 7  ;
    if ( data.size() > 1 ) {
        is = axis1 ;
        sort( data.begin(), data.end(), ScanSort );
    }

    vector<vec> col ;
    vec mdata ;
    double val = data[0][axis1] ;
    bool flush = false ;
    is = axis2 ;
    int cl(0), st(0), de(0) ;

    for (size_t i=0; i < data.size(); i++ ) {

        if ( fabs( val - data[i][axis1] ) < mindX ) {
           col.push_back( data[i] ) ;
           if ( i == data.size() -1 ) flush = true ;
        } else {
           flush = true ;
        }

        if ( flush ) {
           sort( col.begin(), col.end(), ScanSort ) ;
           for ( size_t j=0; j< col.size(); j++ ) {

               // Giving die and site number - Only work if data taking contains all 4 sites
               /* Old data type
               if ( cl%2 == 0 ) de = 1 ;
               if ( cl%2 == 1 ) de = 2 ;

               if ( j%4  < 2 && cl%4  < 2 ) st = 3 ;
               if ( j%4 >= 2 && cl%4  < 2 ) st = 1 ;
               if ( j%4  < 2 && cl%4 >= 2 ) st = 4 ;
               if ( j%4 >= 2 && cl%4 >= 2 ) st = 2 ;
               */
               // new RFF data type 
               // 34x68 die arry
               if ( (int)(col[j][4])%2 == 0 && (int)(col[j][5])%4  < 2 ) st  = 1 ; 
               if ( (int)(col[j][4])%2 == 0 && (int)(col[j][5])%4 >= 2 ) st  = 2 ; 
               if ( (int)(col[j][4])%2 == 1 && (int)(col[j][5])%4  < 2 ) st  = 3 ; 
               if ( (int)(col[j][4])%2 == 1 && (int)(col[j][5])%4 >= 2 ) st  = 4 ; 
               if ( (int)(col[j][5])%2 == 0 ) de = 1 ;
               if ( (int)(col[j][5])%2 == 1 ) de = 2 ;   
               

               // Only turn on this if measuring all the sites and dies
               //if ( de != die || site != st ) continue; 
               //printf(" %d-%d-%d \n", cl, st, de ) ;

               //printf("[%d, %d] %.5f  %.5f \n", (int)j, cl, col[j][axis1], col[j][axis2] ) ;
               // For JetStep
               
               mdata.push_back( (double)cl ) ;
               mdata.push_back( (double)st ) ;
               mdata.push_back( (double)de ) ; 
               mdata.push_back( col[j][14] + col[j][12]  ) ;
               mdata.push_back( -1*col[j][13] - col[j][11]  ) ;
               mdata.push_back( col[j][12] ) ;  // dx - new TLS version
               mdata.push_back( -1*col[j][11] ) ;  // dy - new TLS version
               
               // For pick and placed
               /*
               mdata.push_back( col[j][0] ) ;
               mdata.push_back( col[j][4] ) ;
               mdata.push_back( col[j][5] ) ; 
               mdata.push_back( -1*col[j][14] - col[j][12]  ) ;
               mdata.push_back(    col[j][13] + col[j][11]  ) ;
               mdata.push_back( -1*col[j][12] ) ;  // dx - pick&place version
               mdata.push_back(    col[j][11] ) ;  // dy - pick&place version
               */

               // Store data
               data1.push_back( mdata ) ;
               mdata.clear() ;

               fprintf(logfile, "%d, %d, %d,  % 04.6f, % 04.6f \n", 
                                 cl, st, de, col[j][14], -1*col[j][13]   ) ;
           } 

           flush = false ;
	   val = data[i][axis1];
	   col.clear() ;
	   col.push_back(data[i] ) ;
	   cl++ ;
        }
    } 

    printf(" DPM Size = %d - %d \n", (int)data1.size(), (int)data1[0].size() ) ;
    fclose( logfile ) ;
}

// Read the target XCD items
void Match::ReadXCD( vector<vec>& data1, int nSkip  ) {

    int roi ;
    Input->GetParameters("ROI",           &roi ) ;

    vector<vec> data0 ;
    data0.clear() ;

    string csvfile = cfolder + cFileName ;
    printf(" CSV file = %s \n", csvfile.c_str() ) ;
    Input->ReadFreeCSV( csvfile, data0, nSkip ) ;

    vector<vec> data ;
    data.clear() ;
    data = Sort2D( data0, 4, 5 ) ;

    data1.clear() ;
    string debugFile = path + "debug_xcd.csv" ;
    FILE* logfile = fopen( debugFile.c_str() ,"w");
    //int icfg = 7 ;
    for (size_t i=0; i< data.size(); i++ ) {

        if ( roi  >=0 && (int)(data[i][1]) !=  roi ) continue ; 
        //if ( icfg >=0 && (int)(data[i][3]) !=  icfg ) continue ; 
        if ( data[i][10] >=0 ) continue ;  // is defect ?

        data1.push_back( data[i] ) ;
        fprintf( logfile, "%.0f %.0f %.0f  %.0f  %.0f %4.6f \n", data[i][0], data[i][1], data[i][3], data[i][4], data[i][5],  data[i][11] ) ;

    }

    const int sz = data1.size() ;
    double xcd[sz] ;
    for (size_t i=0 ; i< data1.size(); i++) {
        xcd[i] = data1[i][11] ;
    }
    double mean = gsl_stats_mean( xcd, 1, sz );
    double std  = gsl_stats_sd( xcd, 1, sz );
    fprintf(logfile, " mean = %4.6f  +/-  %4.6f \n", mean, std ) ;

    fclose( logfile ) ;

    printf(" XCD Size = %d - %d \n", (int)data1.size(), (int)data1[0].size() ) ;

}

// Read the target Group XCD items
void Match::ReadXCDGroup( vector<vec>& data1, int nSkip  ) {

    vector<int> roi ;
    Input->GetParameters("ROI",           &roi ) ;

    vector<vec> data0 ;
    data0.clear() ;

    string csvfile = cfolder + cFileName ;
    printf(" CSV file = %s \n", csvfile.c_str() ) ;
    Input->ReadFreeCSV( csvfile, data0, nSkip ) ;

    // Sort data based on die index
    vector<vec> data ;
    data.clear() ;
    data = Sort2D( data0, 4, 5 ) ;

    data1.clear() ;
    string debugFile = path + "debug_xcdgroup.csv" ;
    FILE* logfile = fopen( debugFile.c_str() ,"w");
    /*
    for (size_t i=0; i< data.size(); i++ ) {
        if ( roi  >=0 && (int)(data[i][1]) !=  roi ) continue ; 
        if ( (int)data[i][15] ==0 && (int)data[i][16] ==0  ) continue ;  // invalid measurement, max = min 
        data1.push_back( data[i] ) ;
        fprintf( logfile, "%.0f %.0f %.0f  %.0f  %.0f %4.6f \n", data[i][0], data[i][1], data[i][8], data[i][4], data[i][5],  data[i][12] ) ;
    }
    */

    //Open Root file to store result
    string hfile = path +"h_XCDGroup.root" ;
    TFile* rootFile = new TFile( hfile.c_str() , "RECREATE" );
    rootFile->cd() ;
    // Book histogram
    TH1D* h_dX = new TH1D("h_dX" , "dX",   100, -12.5, 12.5 )  ;
    TH1D* h_dY = new TH1D("h_dY" , "dY",   100, -12.5, 12.5 )  ;
    TH1D* h_dR = new TH1D("h_dR" , "dR",   100, -12.5, 12.5 )  ;

    data1.clear() ;
    vector<double> dx ;
    vector<double> dy ;
    vector<double> dr ;
    typedef pair<int,int> dieXY ;
    vector<dieXY> dd ;
   
    for (size_t i=0; i< data.size(); i++ ) {
        bool select = false ;
        for ( size_t j=0 ; j < roi.size() ; j++ ) {  
            if ( roi[j] == -1 ) select = true ;     
            if ( roi[j] == data[i][1] ) select = true ; 
        } 
        if ( !select ) continue ; 

	// roi sequence is dx, dy, dr
	if (data[i][1] == roi[0] ) dx.push_back( data[i][12]);
	if (data[i][1] == roi[1] ) dy.push_back( data[i][12]);
	//if (data[i][1] == roi[2] ) dr.push_back( data[i][12]);
	if (data[i][1] == roi[0] ) dd.push_back( make_pair( (int)data[i][4], (int)data[i][5]) );

       
	//if (data[i][1] == roi[0] ) h_dX->Fill( data[i][12] ) ;
	//if (data[i][1] == roi[1] ) h_dY->Fill( data[i][12] ) ;
	//if (data[i][1] == roi[2] ) h_dR->Fill( data[i][12] ) ;

        data1.push_back( data[i] ) ;
        //fprintf( logfile, "%.0f %.0f  %.0f  %.0f %4.6f \n", data[i][0], data[i][1], data[i][4], data[i][5],  data[i][12] ) ;
    }

    fprintf(logfile, "\n summary \n\n" ) ;
    
    cout<<" dr size = "<< dr.size() <<", "<< dx.size() <<", "<< dy.size() <<endl; 
    cout<<" dd size = "<< dd.size() <<" data size = "<< data1.size() << endl ; 

    vector<vec> fdata ; 
    vector<double> mdata ;
    double dr_ ;
    for (size_t i=0 ; i< dx.size(); i++) {
        if ( dx[i] == 0 || dy[i] == 0 ) continue ;

        dr_ = sqrt( (dx[i]*dx[i]) + (dy[i]*dy[i]) ) ;

        h_dX->Fill( dx[i] ) ;
        h_dY->Fill( dy[i] ) ;
        h_dR->Fill( dr_ ) ;
        mdata.push_back( (double) dd[i].first ) ;
        mdata.push_back( (double) dd[i].second ) ;
        mdata.push_back( dx[i] ) ;
        mdata.push_back( dy[i] ) ;
        mdata.push_back( dr_ ) ;
 
        fdata.push_back( mdata ) ;
        mdata.clear() ;

        fprintf( logfile, " %d, %d,  %4.6f, %4.6f, %4.6f\n", dd[i].first, dd[i].second, dx[i], dy[i], dr_ ) ;
    }

    //Calculate final mean and standard deviation
    const int sz = (int)fdata.size() ;
    double dxA[sz], dyA[sz], drA[sz] ;
    for (size_t i=0; i< fdata.size() ; i++ ) {
        dxA[i] = fdata[i][2] ;
        dyA[i] = fdata[i][3] ;
        drA[i] = fdata[i][4] ;
    }
    double mean_x = gsl_stats_mean( dxA, 1, sz );
    double std_x  = gsl_stats_sd( dxA, 1, sz );
    double mean_y = gsl_stats_mean( dyA, 1, sz );
    double std_y  = gsl_stats_sd( dyA, 1, sz );
    double mean_r = gsl_stats_mean( drA, 1, sz );
    double std_r  = gsl_stats_sd( drA, 1, sz );
    fprintf(logfile, " mean X = %4.6f  +/-  %4.6f \n", mean_x, std_x ) ;
    fprintf(logfile, " mean Y = %4.6f  +/-  %4.6f \n", mean_y, std_y ) ;
    fprintf(logfile, " mean R = %4.6f  +/-  %4.6f \n", mean_r, std_r ) ;

    fclose( logfile ) ;

    printf(" XCD Size = %d - %d \n", (int)data1.size(), (int)data1[0].size() ) ;

 
   TCanvas* c0 = new TCanvas("c_0","", 0,0, 800, 600);
   c0->SetFillColor(10);
   c0->SetFillColor(10);
   gPad->SetGridx();

   h_dX->Draw();
   c0->Update();
   TString gPlotname = path + "h_dX.png"  ;
   c0->Print( gPlotname ) ;

   h_dY->Draw();
   c0->Update();
   gPlotname = path + "h_dY.png"  ;
   c0->Print( gPlotname ) ;
 
   h_dR->Draw();
   c0->Update();
   gPlotname = path + "h_dR.png"  ;
   c0->Print( gPlotname ) ;

   h_dX->Write() ;
   h_dY->Write() ;
   h_dR->Write() ;
   rootFile->Close() ;
}


// Output id, x, y
// Just Read JetStep Wafer alignment info - first two entries
void Match::ReadPsetAlign( vector<vec>& data1 ) {


    vector<vec> data0 ;
    data0.clear() ;

    string csvfile = cfolder + pFileName ;
    printf(" PAlign file = %s \n", csvfile.c_str() ) ;
    Input->ReadFreeCSV( csvfile, data0, 1 ) ;

    vec mdata ;

    for (size_t i=0; i < 2; i++ ) {

        mdata.push_back( (double)(i+1) ) ;
	mdata.push_back( data0[i][1] ) ;
	mdata.push_back( data0[i][2] ) ;

	data1.push_back( mdata ) ;
	mdata.clear() ;
    } 
 
    printf(" PAlign Size = %d - %d\n", (int)data1.size(), (int)data1[0].size() ) ;
}

// Output id, dx, dy, x, y, st, de
void Match::ReadPsetCFG( vector<vec>& data1, int nSkip ) {

    data1.clear() ;

    vector<vec> data ;
    data.clear() ;

    // Input CFG file  
    string csvfile = cfolder + pFileName ;
    printf(" PCAD file = %s \n", csvfile.c_str() ) ;
    Input->ReadFreeCSV( csvfile, data, nSkip ) ;

    string debugFile = path + "debug_pcad.csv" ;
    FILE* logfile = fopen( debugFile.c_str() ,"w");

    int axis1 = 2 ;
    int axis2 = 1 ;

    if ( data.size() > 1 ) {
        is = axis1 ;
        sort( data.begin(), data.end(), ScanSortN );
    }

    vector<vec> col ;
    vec mdata ;
    double val = data[0][axis1] ;
    bool flush = false ;
    is = axis2 ;
    int cl(0) ;
    int st(0), de(0) ;

    for (size_t i=0; i < data.size(); i++ ) {

        if ( fabs( val - data[i][axis1] ) < mindX ) {
           col.push_back( data[i] ) ;
           if ( i == data.size() -1 ) flush = true ;
        } else {
           flush = true ;
        }

        if ( flush ) {
           sort( col.begin(), col.end(), ScanSort ) ;
           for ( size_t j=0; j< col.size(); j++ ) {

               // Giving die and site number
               //printf("[%d, %d] %.5f  %.5f\n", (int)col[j][2], (int)col[j][3], col[j][axis1], col[j][axis2] ) ;
               de = (j%2)+1 ;
               if ( cl%4 < 2 ) {
                  st =  (j%4 < 2 ) ? 1 : 2 ;
               } else {
                  st =  (j%4 < 2 ) ? 3 : 4 ;
               }

	       int tb = cl%2 ;
	       bool selectSite = ( st == site || site == -1 ) ? true : false ;
	       bool selectDie  = ( de == die  || die  == -1 ) ? true : false ;
               bool selectTob  = ( tb == tob  || tob  == -1 ) ? true : false ;

	       // 0: cl, 1: dX, 2:dY, 3:X, 4:Y, 5:st, 6:de, 7:X_shot, 8:Y_shot 
               if ( selectSite && selectDie && selectTob ) { 
                  mdata.push_back( (double)cl ) ;
		  mdata.push_back( col[j][6] ) ;  // dX
		  mdata.push_back( col[j][7] ) ;  // dY
		  mdata.push_back( col[j][axis2] ) ; // X
		  mdata.push_back( col[j][axis1] ) ; // Y
		  mdata.push_back( st ) ;
		  mdata.push_back( de ) ;
		  mdata.push_back( col[j][4] ) ;  // X_shot
		  mdata.push_back( col[j][5] ) ;  // Y_shot

		  data1.push_back( mdata ) ;
                  mdata.clear() ;

                  fprintf(logfile, "%d, %d, %d, % 04.6f, % 04.6f \n", 
                                    cl, st, de, col[j][1], col[j][2]   ) ;
	       }
           } 

           flush = false ;
	   val = data[i][axis1];
	   col.clear() ;
	   col.push_back(data[i] ) ;
	   cl++ ;
        }
    } 

    fclose( logfile ) ;
    printf(" PCAD Size = %d - %d\n", (int)data1.size(), (int)data1[0].size() ) ;
}

void Match::SelectPsetCFG( vector<vec>& data ) {

    //string debugFile = path + "debug_pcad1.csv" ;
    string debugFile = path + "PTIPANEL.PR2.CFG" ;
    FILE* logfile = fopen( debugFile.c_str() ,"w");

    int idx = 3 ;
    for ( std::vector<vec>::iterator it = data.begin() ; it != data.end();  ) {
        //if ( it->size() < 10 || (*it)[3] < 0 ) { 
        int col = (int) (*it)[0] ;
        if ( col%2 != 0 ) { 
           data.erase( it ) ; 
           //continue ;
        } else {
           //fprintf(logfile, "%.0f, % 04.6f, % 04.6f \n", 
           //                     (*it)[0],  (*it)[3], (*it)[4]   ) ;
           fprintf(logfile, "%d, % 04.6f,% 04.6f, M, % 04.6f,% 04.6f, % 04.6f,% 04.6f \n", 
                             idx,  (*it)[3], (*it)[4],  (*it)[7], (*it)[8], (*it)[1], (*it)[2]   ) ;
           it++ ;
           idx++ ;
        }
    }

    fclose( logfile ) ;
    printf(" PCAD1 Size = %d - %d\n", (int)data.size(), (int)data[0].size() ) ;
}

// Output id, x, y, dx, dy
// This is the correction file for jetstep
void Match::ReadPset( vector<vec>& data1 ) {

    vector<vec> data0 ;
    data0.clear() ;

    string PsetName = "pset1.csv" ;

    string csvfile = cfolder + PsetName ;
    printf(" PSET file = %s \n", csvfile.c_str() ) ;
    Input->ReadFreeCSV( csvfile, data0, 2 ) ;

    data1.clear() ;
    data1 = Sort2D( data0, 1, 0, false, true ) ;

    string debugFile = path + "debug_pset0.csv" ;
    FILE* logfile = fopen( debugFile.c_str() ,"w");
  
    for (size_t i=0; i< data1.size(); i++ ) {
        fprintf(logfile, "%d, %04.6f, %04.6f, %04.6f, %04.6f \n", 
                         (int)i, data1[i][0],  data1[i][1],  data1[i][2],  data1[i][3]  ) ;
    }


    fclose( logfile ) ;
    printf(" PSET Size = %d - %d\n", (int)data1.size(), (int)data1[0].size() ) ;
}

void Match::Minimum_GdR( FILE* logfile, vector<vec>& vNom, vector<vec>& vMes ) {

    int ax1 = 3 ;
    int ax2 = ax1 + 1 ;

    // Determine the range for minimum
    double c0,c1,cov00, cov01, cov11, sumsq ;
    double x_max(0), y_max(0), a_max(0) ;
    double val = vMes[0][ax1] ;
    bool flush = false ;
    vector<vec> colm ;
    vector<vec> coln ;
    double pos1(0) ; 
    for ( size_t i=0 ; i< vMes.size() ; i++ ) {

        if ( fabs( val - vMes[i][ax1] ) < mindX ) {
           colm.push_back( vMes[i] ) ;
           coln.push_back( vNom[i] ) ;
           pos1 = vNom[i][ax1] ; 
           if ( i == vMes.size() -1 ) flush = true ;
        } else {
           flush = true ;
        }

        if ( flush ) {
           const int csz = (int) colm.size() ;
	   double xa[csz], ya[csz], dxa[csz] ;
	   for ( int j=0; j < csz; j++ ) {
               // For each vec, their x values are similar but y are different
               // Therefor, flipping x and y for fitting 
               //printf("%d, %d, %d - %.5f, %d - %.5f \n",(int)i , csz, (int)j, col[j][ax2], (int)i-csz+j, vNom[i-csz+j][ax2] ) ;
               //printf("(%d) %.5f , %.5f |  %.5f , %.5f \n",(int)i , colm[j][ax1], colm[j][ax2], coln[j][ax1], coln[j][ax2] ) ;
               //dxa[j] = col[j][ax2] - vNom[i-csz+j][ax2] ;
               dxa[j] = colm[j][ax2] - coln[j][ax2] ;
               xa[j]  = colm[j][ax2] ;
               ya[j]  = colm[j][ax1] ;
           }
           gsl_fit_linear( xa, 1, ya, 1, csz, &c0, &c1, &cov00, &cov01, &cov11, &sumsq ) ;
	   double dx_ = c0 - pos1 ;
           dx_ = ( dx_ > 0 ) ? dx_ + sqrt(cov00) : dx_ - sqrt(cov00) ;
	   x_max = ( dx_ > x_max ) ? dx_ : x_max ;

           double mean_dY = gsl_stats_mean( dxa, 1, csz );
           double std_dY  = gsl_stats_sd( dxa, 1, csz );
           double dy_ = ( mean_dY > 0 ) ? mean_dY + std_dY : mean_dY - std_dY ;
	   y_max = ( fabs(dy_) > fabs(y_max) ) ? dy_ : y_max ; 

           printf(" theta = %f +/- %f \n", c1, sqrt(cov11) ) ;
           double da_ = ( c1 > 0 ) ? c1 + sqrt(cov11) : c1 - sqrt(cov11) ;
           a_max = ( fabs(da_) > fabs(a_max) ) ?  da_ : a_max ;

           flush = false ;
	   val = vMes[i][ax1];
	   coln.clear() ;
	   colm.clear() ;
	   coln.push_back( vNom[i] ) ;
	   colm.push_back( vMes[i] ) ;
        }
    }

    printf(" xmax: %.5f, ymax: %.5f , amax: %.9f \n", x_max, y_max, a_max ) ;
    //Global dR minimization
    int nStep = 100 ;
    double gdR = 9999999. ;
    double x_(0), y_(0) ;
    double dx(0), dy(0), da(0) ; 
    double xs(0), ys(0), as(0) ;

    // Tripple loop for x, y , angle shift
    for (int i=0; i < nStep ; i++ ) {
        dx +=  ( x_max/ (double)nStep )  ;
        dy = 0 ;
        for (int j=0; j < nStep ; j++ ) {
            dy +=  ( y_max/ (double)nStep )  ;
            da = -1.*a_max ;
            for (int k=0; k < nStep ; k++ ) {
                da += ( 2*a_max/ (double)nStep )  ;

                double gdR_ = 0 ;
                for ( size_t m = 0; m< vNom.size(); m++ ) {
		    x_ = vMes[m][ax1] - dx ;
		    y_ = vMes[m][ax2] - dy ;
		    Rotation( x_, y_, da ) ;

                    double da1 = x_ - vNom[m][ax1] ;
                    double da2 = y_ - vNom[m][ax2] ;
                    double dR2 = (da1*da1) + (da2*da2) ;

                    gdR_ += dR2 ; 
                }
                //printf(" GdR = %f, shift = %f %f %f \n", gdR, dx, dy, da ) ;
                if ( gdR_ < gdR ) {
                   gdR = gdR_ ;
                   xs = dx ;
                   ys = dy ;
                   as = da ;
                }


            }
        }
    } 
    printf(" Matching Done \n") ;
    printf(" Final GdR = %f, shift = %f %f %f \n", gdR, xs, ys,as ) ;

    x_ = 0 ;
    y_ = 0 ;
    gdR = 0 ;
    fprintf(logfile, "%f, %f, %f \n", xs, ys,as ) ;
    for ( size_t m = 0; m< vMes.size(); m++ ) {
        x_ = vMes[m][ax1] - xs ;
	y_ = vMes[m][ax2] - ys ;
	Rotation( x_, y_, as ) ;

	double da1 = x_ - vNom[m][ax1] ;
	double da2 = y_ - vNom[m][ax2] ;
	double dR2 = (da1*da1) + (da2*da2) ;
        fprintf(logfile, " %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f \n", 
                          vNom[m][ax1], vNom[m][ax2], vMes[m][ax1], vMes[m][ax2], x_, y_, da1, da2, vNom[m][1], vNom[m][2]  ) ;

	gdR += dR2 ; 
    }
}

// Rotate (x,y) CounterClockwise 
//  |x'|    [ cos   -sin ] | x | 
//  |y'| =  [ sin    cos ] | y |
void Match::Rotation( double& x, double& y, double theta ) {

     double x_ = x*cos(theta) - y*sin(theta) ;
     double y_ = x*sin(theta) + y*cos(theta) ;

     x = x_ ;
     y = y_ ;
}

// theta is the angle between y and y' 
//  |x'|    [ 1   Sin ] | x | 
//  |y'| =  [ 0   Cos ] | y |
void Match::Orth( double& x, double& y, double theta ) {

     double x_ = x + y*sin(theta) ;
     double y_ =     y*cos(theta) ;
     x = x_ ;
     y = y_ ;
}


// Sort data on axis1 , within same axis1, sort by axis2
vector<vec> Match::Sort2D( vector<vec>& data, int axis1, int axis2, bool rise1, bool rise2 ) {

    vector<vec> data1 ;

    // Sort the data - Example
    is = axis1 ;
    if ( data.size() > 1 ) {
       if ( rise1) sort( data.begin(), data.end(), ScanSort );
       if (!rise1) sort( data.begin(), data.end(), ScanSortN );
    }

    vector<vec> col ;
    col.clear() ;
    bool nextCol = false ;
    double x_ = data[0][axis1] ;
    is = axis2 ;

    for (size_t i=0; i< data.size(); i++ ) {
        // Accumulate the same column data
        if ( fabs(data[i][axis1] - x_ ) < mindX ) {
           col.push_back( data[i] ) ;
           if ( i == data.size() -1 ) nextCol = true ;
        } else {
           nextCol = true ;
        }

        // Sort the column
        if ( nextCol ) {
           //cout<<" col size = "<< col.size() <<endl ;
           if ( rise2)  sort( col.begin(), col.end(), ScanSort ) ;
           if (!rise2)  sort( col.begin(), col.end(), ScanSortN ) ;

           // Fill up the new sorted data
           for ( size_t j=0; j < col.size() ; j++ ) {
               data1.push_back( col[j] ) ;
           }
           col.clear() ;
           col.push_back( data[i] ) ;
           nextCol = false ;
        }
        x_ = data[i][axis1] ;
    }

    return data1 ;

}

void Match::ReadPickPlace( vector<vec>& data1 ) {

    data1.clear() ;

    vector<vec> data ;
    data.clear() ;

    // Read the design file name
    string csvfile = cfolder + pFileName ;
    printf(" PickPlace file = %s \n", csvfile.c_str() ) ;
    Input->ReadFreeCSV( csvfile, data, 1 ) ;

    data1 = Sort2D( data, 4, 5) ;

}

// Read and sort data 
void Match::ReadnSort( vector<vec>& data1, int axis1, int axis2, int nSkip, bool rise1, bool rise2  ) {

    data1.clear() ;

    vector<vec> data ;
    data.clear() ;

    // Read the design file name
    string csvfile = cfolder + cFileName ;
    printf(" Input CSV file = %s \n", csvfile.c_str() ) ;
    Input->ReadFreeCSV( csvfile, data, nSkip ) ;

    data1 = Sort2D( data, axis1, axis2, rise1, rise2) ;

}

void Match::ReadMask( vector<vec>& data1 ) {

    data1.clear() ;

    vector<vec> data ;
    data.clear() ;

    // Read the design file name
    string csvfile = cfolder + pFileName ;
    printf(" Mask Design file = %s \n", csvfile.c_str() ) ;
    Input->ReadFreeCSV( csvfile, data, 1 ) ;

    data1 = Sort2D( data, 3, 4) ;

}

// Read StepFAST panel after StepFAST correction
//  0     1    2     3      4         5      6  7   8   9
// Site, Die, Tob, Valid, Column, NoMeaning, X, Y, dX, dY
void Match::ReadStepFAST( vector<vec>& output ) {

     // log file
     FILE* logfile = fopen( logFileName.c_str() ,"w");

     // Read Data 
     vector<vec> sfV ;
     ReadnSort( sfV, 11, 10, 1, false ) ;


     // Output 
     int cl(0) ;
     vector<vec> col ;
     vector<vec> dataV ;
     bool flush = false ;
     double val = sfV[0][11] ;
     size_t data_sz = sfV[0].size() ;

     for ( int i=0; i< (int)sfV.size() ; i++ ) {
         // debug for the entry - some entries have missing data
	 while ( sfV[i].size() < data_sz ) {
               sfV[i].push_back(0.) ;
         }
         if ( fabs( val - sfV[i][11] ) < mindX ) {
 
            col.push_back( sfV[i] ) ;
            if ( i == (int)sfV.size() -1 ) flush = true ;
         } else {
            flush = true ;
    	    val = sfV[i][11] ;
         }

        if ( flush ) {
           for ( int j=0 ; j < (int)col.size(); j++ ) {

               if ( j%4 == 0 && cl%4 < 2 )  { 
		  col[j].push_back(1) ; // st  1
		  col[j].push_back(1) ; // die 1
		  col[j].push_back(cl%2) ; // top or buttom 
		  col[j].push_back(cl) ; // column number
	       }
               if ( j%4 == 1 && cl%4 < 2 ) {
		  col[j].push_back(1) ; // st  1
		  col[j].push_back(2) ; // die 2
		  col[j].push_back(cl%2) ; // top or buttom 
		  col[j].push_back(cl) ; // column number
	       }
               if ( j%4 == 2 && cl%4 < 2 )  { 
		  col[j].push_back(2) ; // st  2
		  col[j].push_back(1) ; // die 1
		  col[j].push_back(cl%2) ; // top or buttom 
		  col[j].push_back(cl) ; // column number
	       }
               if ( j%4 == 3 && cl%4 < 2 ) {
		  col[j].push_back(2) ; // st  2
		  col[j].push_back(2) ; // die 2
		  col[j].push_back(cl%2) ; // top or buttom 
		  col[j].push_back(cl) ; // column number
	       }

               if ( j%4 == 0 && cl%4 >= 2 )  { 
		  col[j].push_back(3) ; // st  3
		  col[j].push_back(1) ; // die 1
		  col[j].push_back(cl%2) ; // top or buttom 
		  col[j].push_back(cl) ; // column number
	       }
               if ( j%4 == 1 && cl%4 >= 2 ) {
		  col[j].push_back(3) ; // st  3
		  col[j].push_back(2) ; // die 2
		  col[j].push_back(cl%2) ; // top or buttom 
		  col[j].push_back(cl) ; // column number
	       }
               if ( j%4 == 2 && cl%4 >= 2 )  { 
		  col[j].push_back(4) ; // st  4
		  col[j].push_back(1) ; // die 1
		  col[j].push_back(cl%2) ; // top or buttom 
		  col[j].push_back(cl) ; // column number
	       }
               if ( j%4 == 3 && cl%4 >= 2 ) {
		  col[j].push_back(4) ; // st  4
		  col[j].push_back(2) ; // die 2
		  col[j].push_back(cl%2) ; // top or buttom 
		  col[j].push_back(cl) ; // column number
               }

	       // ----  Debug ----
	       if (col[j].size() != 20 ) {
		  printf(" %d, %d, ", j, cl ) ;
                  for (size_t k=0; k< col[j].size(); k++) { 
	  	      printf(" %.4f,", col[j][k] ) ;
	          }
	          printf("= > %d \n", (int)col[j].size() ) ;
	       }
               // ----------------

	       dataV.push_back( col[j] ) ;
	   }
	   cl++ ;
	   col.clear() ;
           col.push_back( sfV[i] ) ;
	   flush = false ;
        } 
     }
 
     vec mdata ;
     for ( int i=0; i< (int)dataV.size() ; i++ ) {
	    double valid = ( dataV[i][12] == 0.0 && dataV[i][13] == 0.0 ) ? 0. : 1. ; 
            mdata.clear() ;
	    mdata.push_back( dataV[i][16] ) ;     // site
	    mdata.push_back( dataV[i][17] ) ;     // die
	    mdata.push_back( dataV[i][18] ) ;     // top or buttom
	    mdata.push_back( valid  ) ;           // valid or not
	    mdata.push_back( dataV[i][19] ) ;     // column number 
	    mdata.push_back( dataV[i][5] ) ;      // No meaning
	    mdata.push_back( -1*dataV[i][11] ) ;  // X
	    mdata.push_back( dataV[i][10] ) ;     // Y
	    mdata.push_back( -1*dataV[i][13] ) ;  // dX
	    mdata.push_back( dataV[i][12] ) ;     // dY
            output.push_back( mdata ) ;
	    if ( dataV[i][16] == site && dataV[i][17] == die && dataV[i][18] == tob ) {
               fprintf( logfile, "%.0f,%.0f,%.0f, %.0f, %.0f, %.4f, %.4f, %.6f, %.6f\n", 
	             dataV[i][16], dataV[i][17], dataV[i][18], valid, dataV[i][19], -1.*dataV[i][11], dataV[i][10], -1.*dataV[i][13], dataV[i][12] ) ; 
	    }
     }
     fclose( logfile ) ;
}
