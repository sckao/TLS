#include "Module.h"
#include "Match.h"

Module::Module( ) {

  Input = AnaInput::Instance() ;

  Input->GetParameters("PlotType",      &plotType ) ;
  Input->GetParameters("PlotName0",     &plotName0 );
  Input->GetParameters("HFileName",     &hFileName );

  Input->GetParameters("Path",          &hfolder ) ;
  Input->GetParameters("debug",         &debug );
  Input->GetParameters("Site",          &site );
  Input->GetParameters("Die",           &die );
  Input->GetParameters("Tob",           &tob );

  Input->GetParameters("OutputFile",    &outFileName );
  logFileName = hfolder +  outFileName ;
  logFileName1 = hfolder + "pch_" + outFileName ;
  logFileName2 = hfolder + "xy_" + outFileName ;

}

Module::~Module(){

  //delete Input ;
    cout<<" done ! "<<endl ;
}


void Module::Analysis() {

   ///// ============ DPM Output ============================================================ 
   //   0    1      2         3         4        5       6       7        8         9    10  11   12     13     14   
   //  ID,ROI_ID,ROI_Name,Inspex_ID,XDie_Idx,YDie_Idx,XCenter,YCenter,IsValid?,IsDefective?,Defect ID,OffsetX_MM,OffsetY_MM,XExpectedLocationSWCS_MM,YExpectedLocationSWCS_MM,AngleDegrees,WasDieFound,IsInterpolated,XInspectionSiteIndex,YInspectionSiteIndex

   // Logfile 
   FILE* logfile = fopen( logFileName.c_str() ,"w");
   // Read data
   vector<vec> data ;
   Match  *match  = new Match( ) ;
   // Sort Data by die address X and Y
   match->ReadnSort( data, 4, 5, 1 );
   printf(" Read Site %d Die %d Tob %d \n", site, die, tob) ;

   // Figure out Number of columns and rows 
   int nCol = ((int)data[ data.size()-1 ][4] +1) /2 ;
   int nRow = ((int)data[ data.size()-1 ][5] +1) /4 ;

   // ud indicates top or bottom fiducial  
   int st(0),de(0),tb(9), vd(1) ;
   vector<vec> output ;
   vec mdata ;
   for (size_t i=0; i < data.size(); i++) {

       // label the data with site,die and top/bottom fiducial
       if ( (int)data[i][5]%2 == 0 ) de = 1 ; 
       if ( (int)data[i][5]%2 == 1 ) de = 2 ; 
       if ( (int)data[i][4]%2 == 0 && (int)data[i][5]%4  < 2 ) st  = 1 ; 
       if ( (int)data[i][4]%2 == 0 && (int)data[i][5]%4 >= 2 ) st  = 2 ; 
       if ( (int)data[i][4]%2 == 1 && (int)data[i][5]%4  < 2 ) st  = 3 ; 
       if ( (int)data[i][4]%2 == 1 && (int)data[i][5]%4 >= 2 ) st  = 4 ; 
       if ( i%2 == 0 &&  (int)data[i][4] == (int)data[i+1][4] && (int)data[i][5] == (int)data[i+1][5] ) {
          if ( data[i][6] < data[i+1][6] ) tb = 0 ; 
	  else tb = 1 ;
       }
       if ( i%2 == 1 &&  (int)data[i][4] == (int)data[i-1][4] && (int)data[i][5] == (int)data[i-1][5] ) {
          if ( data[i][6] > data[i-1][6] ) tb = 1 ; 
	  else tb = 0 ;
       }
       // Is data Valid ? 
       if ( data[i][11] ==  0.0 && data[i][12] == 0.0 ) {
	   vd = 0 ;
       } else {
	   vd = 1  ;
       }

       // Store useful information
       if ( st == site && de == die && tb == tob ) {
          mdata.clear() ;
	  mdata.push_back( st ) ;
	  mdata.push_back( de ) ;
	  mdata.push_back( tb ) ;
	  mdata.push_back( vd ) ;
	  mdata.push_back( data[i][4] ) ;
	  mdata.push_back( data[i][5] ) ;
	  mdata.push_back( data[i][11] ) ;
	  mdata.push_back( data[i][12] ) ;
	  mdata.push_back( data[i][13] ) ;
	  mdata.push_back( data[i][14] ) ;
          output.push_back( mdata ) ; 

          fprintf( logfile, "%1d, %1d, %1d, %1d, %3.0f, %3.0f, %4.6f, %4.6f, %4.6f, %4.6f \n",
                          st,de,tb, vd, data[i][4], data[i][5], data[i][11], data[i][12], data[i][13], data[i][14] ) ;

       }
   } 

   printf(" Col (%d) x Row (%d) \n", nCol, nRow ) ;
   fclose( logfile ) ;

   GetPitch( output ) ;

   GetXY( output ) ;

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

void Module::GetPitch( vector<vec>& data ) {

     FILE* logfile = fopen( logFileName1.c_str() ,"w");
     const int xsz = 17 ;
     const int ysz = 17 ;

     // Fill X Y position in two arrays
     double xA[xsz][ysz] ;	
     double yA[xsz][ysz] ;
     bool   bA[xsz+1][ysz+1] ;
     int ix = -1 ;
     int iy = 0 ;     
     for ( size_t i=0 ; i < data.size() ; i++ ) {
	 iy = i% ysz ;
	 if (iy == 0 ) ix++ ; 
         xA[ix][iy] = data[i][6] + data[i][8] ; 
         yA[ix][iy] = data[i][7] + data[i][9] ;

	 if ( data[i][3] == 0 ) bA[ix][iy] = false ;
	 else                   bA[ix][iy] = true  ;
	 //printf(" (%d,%d) = [ %4.5f, %4.5f ]\n", ix, iy, xA[ix][iy], yA[ix][iy] ) ;
     }
     for ( int i=0; i< xsz+1; i++) { bA[i][ysz] = false  ; }
     for ( int i=0; i< ysz+1; i++) { bA[xsz][i] = false  ; }

     // Check dY w.r.t. different X
     double dy_  = 0 ;
     double dx_  = 0 ;
     double dxA[xsz][ysz], dyA[xsz][ysz] ;

     for ( int j=0 ; j < ysz ; j++ ) {
         for ( int i= 0; i < xsz; i++ ) {
             if (bA[i][j+1] && bA[i][j]) {
                dy_ = yA[i][j+1] - yA[i][j] ;
		dyA[i][j] = dy_ ;
	     } else  {
		dyA[i][j] = 0 ;
	     }

             if (bA[i+1][j] && bA[i][j] ) {
                dx_ = xA[i+1][j] - xA[i][j] ;
		dxA[i][j] = dx_ ;
	     }  else {
		dxA[i][j] = 0 ;
	     }
         }
     }

     // Print out dX 
     fprintf(logfile," dx \n" ) ;
     vector<double> dxV ;
     for ( int j=ysz-1 ; j >=0  ; j-- ) {
         for ( int i= 0; i < xsz; i++ ) {
             fprintf( logfile, " %.4f,", dxA[i][j] ) ;
	     if (dxA[i][j] != 0 ) dxV.push_back( dxA[i][j] ) ;
      	 }
	 pair<double,double> stats = gsl_vector_mean(dxV );
         fprintf(logfile, " %.5f, %.5f \n", stats.first, stats.second ) ;
         dxV.clear() ;
     }
     // Calculate mean and stdev w.r.t. each X
     vector<double> stdV ;
     for ( int i=0 ; i < xsz ; i++ ) {
         for ( int j= 0; j < ysz; j++ ) {
	     if (dxA[i][j] != 0 ) dxV.push_back( dxA[i][j] ) ;
      	 }
	 pair<double,double> stats = gsl_vector_mean(dxV );
	 // output mean for each column(x)
         fprintf(logfile, " %.5f,", stats.first ) ;
	 stdV.push_back( stats.second ) ;
         dxV.clear() ;
     }
     // output stdev for each column(x)
     fprintf(logfile,"\n") ;
     for ( size_t i =0; i< stdV.size(); i++ ) {
	 fprintf(logfile, " %.5f,", stdV[i] ) ;
     }
     stdV.clear() ;

     // Print out dY
     fprintf(logfile,"\n dy \n" ) ;
     vector<double> dyV ;
     for ( int j=ysz-1 ; j >=0  ; j-- ) {
         for ( int i= 0; i < xsz; i++ ) {
             fprintf( logfile, " %.4f,", dyA[i][j] ) ;
	     if (dyA[i][j] != 0 ) dyV.push_back( dyA[i][j] ) ;
      	 }
	 pair<double,double> stats = gsl_vector_mean(dyV );
         fprintf(logfile, " %.5f, %.5f \n", stats.first, stats.second ) ;
         dyV.clear() ;
     }
     for ( int i=0 ; i < xsz ; i++ ) {
         for ( int j= 0; j < ysz; j++ ) {
	     if (dyA[i][j] != 0 ) dyV.push_back( dyA[i][j] ) ;
      	 }
	 pair<double,double> stats = gsl_vector_mean(dyV );
	 // output mean for each column(x)
         fprintf(logfile, " %.5f,", stats.first ) ;
	 stdV.push_back( stats.second ) ;
         dyV.clear() ;
     }
     // output stdev for each column(x)
     fprintf(logfile,"\n") ;
     for ( size_t i =0; i< stdV.size(); i++ ) {
	 fprintf(logfile, " %.5f,", stdV[i] ) ;
     }
     stdV.clear() ;

     fclose( logfile ) ;
}

void Module::GetXY( vector<vec>& data ) {

     FILE* logfile = fopen( logFileName2.c_str() ,"w");
     const int xsz = 17 ;
     const int ysz = 17 ;

     // Fill X Y position in two arrays
     double xA[xsz][ysz] ;	
     double yA[xsz][ysz] ;
     bool   bA[xsz][ysz] ;
     int ix = -1 ;
     int iy = 0 ;     
     for ( size_t i=0 ; i < data.size() ; i++ ) {
	 iy = i% ysz ;
	 if (iy == 0 ) ix++ ; 
         xA[ix][iy] = data[i][6] + data[i][8] ; 
         yA[ix][iy] = data[i][7] + data[i][9] ;

	 if ( data[i][3] == 0 ) bA[ix][iy] = false ;
	 else                   bA[ix][iy] = true  ;
	 //printf(" (%d,%d) = [ %4.5f, %4.5f ]\n", ix, iy, xA[ix][iy], yA[ix][iy] ) ;
     }

     // Print out X 
     fprintf(logfile," X \n" ) ;
     for ( int j=ysz-1 ; j >=0  ; j-- ) {
         for ( int i= 0; i < xsz; i++ ) {
             if (bA[i][j]) fprintf( logfile, " %.4f,", xA[i][j] ) ;
	     else          fprintf( logfile, " %.4f,", 0.0 ) ;
      	 }
	 fprintf( logfile, "\n" ) ;
     }

     // Fit each column for X position X = c0 + c1*Y
     vector<double> xV ;
     vector<double> yV ;
     vector<double> c1V ;
     for ( int i=0 ; i < xsz ; i++ ) {
         for ( int j= 0; j < ysz; j++ ) {
	     if ( bA[i][j] ) {
		xV.push_back( xA[i][j] ) ;
		yV.push_back( yA[i][j] ) ;
	     }
      	 }
	 pair<double,double> stats = gsl_vector_linear_fit( yV, xV );
	 // output mean for each column(x)
         fprintf(logfile, " %.4f,", stats.first ) ;
	 c1V.push_back( stats.second ) ;
         xV.clear() ;
         yV.clear() ;
     }
     // output stdev for each column(x)
     fprintf(logfile,"\n") ;
     for ( size_t i =0; i< c1V.size(); i++ ) {
	 fprintf(logfile, " %.3f,", c1V[i]*1000000.0 ) ;
     }
     c1V.clear() ;
 
     // Print out Y
     // Fit each row for Y position Y = c0 + c1*x
     fprintf(logfile,"\n Y \n" ) ;
     for ( int j=ysz-1 ; j >=0 ; j-- ) {
         for ( int i= 0; i < xsz; i++ ) {
             if ( bA[i][j] ) { 
		fprintf( logfile, " %.4f,", yA[i][j] ) ;
	        yV.push_back( yA[i][j] ) ;
	        xV.push_back( xA[i][j] ) ;
	     } else {
		fprintf( logfile, " %.4f,", 0.0 ) ;
	     }

      	 }
	 pair<double,double> stats = gsl_vector_linear_fit( xV, yV );
         fprintf(logfile, " %.4f, %.3f \n", stats.first, stats.second*1000000.0 ) ;
         yV.clear() ;
         xV.clear() ;
     }

     fclose( logfile ) ;
}

// Return mean and stdev pair 
pair<double,double> Module::gsl_vector_mean( vector<double>& v ) {
      

      const int sz = v.size() ;
      double vA[sz] ;
      for (int i = 0; i < sz ; i++ ) {
          vA[i] = v[i] ;
	  //printf("%.4f, ", vA[i]);
      }    
      double mean_v =  gsl_stats_mean( vA, 1, sz );
      double std_v  =  gsl_stats_sd( vA, 1, sz );
      pair<double,double> ms = make_pair( mean_v , std_v ) ;

      return ms ;
}


pair<double,double> Module::gsl_vector_linear_fit( vector<double>& vx, vector<double>& vy ) {
      

      const int sz = vx.size() ;
      double ax[sz], ay[sz] ;
      double c0,c1,cov00, cov01, cov11, sumsq ;
      for (int i = 0; i < sz ; i++ ) {
          ax[i] = vx[i] ;
          ay[i] = vy[i] ;
	  //printf("%.4f, ", vA[i]);
      }   
      gsl_fit_linear( ax, 1, ay, 1, sz, &c0, &c1, &cov00, &cov01, &cov11, &sumsq ) ;

      pair<double,double> ms = make_pair( c0 , c1 ) ;

      return ms ;
}
