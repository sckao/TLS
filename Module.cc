#include "Module.h"

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
  logFileName3 = hfolder + "dXY_" + outFileName ;
  logFileName4 = hfolder + "diff_" + outFileName ;

  match  = new Match( ) ;
}

Module::~Module(){

  //delete Input ;
  delete match ;
  cout<<" done ! "<<endl ;

}


void Module::Analysis() {

   ///// ============ DPM Output ==================================== 
   //   0    1    2     3     4      5    6   7   8   9    
   // site, die, tob, valid, DieX, DieY,  X,  Y, dX, dY

   // Logfile 
   FILE* logfile = fopen( logFileName.c_str() ,"w");
   // Read data
   vector<vec> data ;
   //Match  *match  = new Match( ) ;
   // Sort Data by die address X and Y
   match->ReadnSort( data, 4, 5, 1 );
   printf(" Read Site %d Die %d Tob %d Data Size = %d \n", site, die, tob, (int)data.size() ) ;
  
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

       bool select_site = ( st == site || site == -1 ) ;
       bool select_die  = ( de == die  || die  == -1 ) ;
       bool select_tob  = ( tb == tob  || tob  == -1 ) ;
       // Store useful information
       if ( select_site && select_die && select_tob ) {
          mdata.clear() ;
	  mdata.push_back( st ) ;
	  mdata.push_back( de ) ;
	  mdata.push_back( tb ) ;
	  mdata.push_back( vd ) ;
	  mdata.push_back( data[i][4] ) ;
	  mdata.push_back( data[i][5] ) ;
	  mdata.push_back( data[i][13] ) ;  // X
	  mdata.push_back( data[i][14] ) ;  // Y
	  mdata.push_back( data[i][11] ) ;  // dX
	  mdata.push_back( data[i][12] ) ;  // dY
          output.push_back( mdata ) ; 

          fprintf( logfile, "%1d, %1d, %1d, %1d, %3.0f, %3.0f, %4.6f, %4.6f, %4.6f, %4.6f \n",
                          st,de,tb, vd, data[i][4], data[i][5], data[i][11], data[i][12], data[i][13], data[i][14] ) ;

       }
   } 

   fclose( logfile ) ;
   printf(" Col (%d) x Row (%d), output size = %d \n\n", nCol, nRow,  (int)output.size() ) ;

   // Get Pitch  
   GetPitch( output ) ;

   // Get Raw XY  
   vector<swcs> xyV ;
   GetXY( output, xyV ) ;

   // Get Wafer Alignment data 
   // This function use JetStep Coordinate system
   vector<vec> wA ;
   double theta = match->ReadWaferAlignment( wA, 1 );

   double dX0(0), dY0(0) ;
   CenterCorrection( wA, dX0, dY0 ) ;

   // Rotation angle theta and orth angle phi
   AxisCorrection( output, (theta*1.) , -0.0000 ) ;

   // Read PSET file
   vector<vec> psetV ;
   match->ReadPsetCFG( psetV ) ;

   // Compare PSET and FFLY_Corrected
   CompareFFLY_PSET( output, psetV, -1*dY0, dX0 ) ;

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
   //delete match ;
}

// Compare Corrected (rotation, orth, shift) FFLY data with PSET
void Module::CompareFFLY_PSET( vector<vec>& fflyV, vector<vec>& psetV0, double sx, double sy ){

  // Rotate and sort psetV coordinate system 
  // Y_pset = - X_ffly  
  // X_pset =   Y_ffly 
  vector<vec> psetV ;
  psetV = match->Sort2D( psetV0, 4, 3, false, true ) ; 
  
  FILE* logfile = fopen( logFileName4.c_str() ,"w");
  double fflyX(0), fflyY(0), dx(0), dy(0) ;
  for (size_t i=0; i< psetV.size(); i++ ) {
      
      fflyX = fflyV[i][6] + fflyV[i][8] + sx;
      fflyY = fflyV[i][7] + fflyV[i][9] + sy;
      dx = (-1*psetV[i][4]) - fflyX ;
      dy =    psetV[i][3] - fflyY ;

      fprintf( logfile, " %.4f, %.4f, %.6f, %.6f, %.6f, %.6f \n", 
                  -1*psetV[i][4], psetV[i][3], fflyX, fflyY,  dx, dy ) ;
  }
 
  fclose( logfile ) ;

}

// Analyze data after StepFAST correction
void Module::AnalyzeStepFAST() {

   // Read data
   vector<vec> output ;
   //Match  *match  = new Match( ) ;
   // Sort Data by die address X and Y
   match->ReadStepFAST( output );
   printf(" Read Site %d Die %d Tob %d \n", site, die, tob) ;

   vector<vec> sample ;

   bool select_site = ( site == -1 ) ? true : false ;
   bool select_die  = ( die  == -1 ) ? true : false ;
   bool select_tob  = ( tob  == -1 ) ? true : false ;
   
   for ( size_t i =0 ; i < output.size() ; i++ ) {
       select_site = ( output[i][0] == site || site == -1 ) ;
       select_die  = ( output[i][1] == die  || die  == -1 ) ;
       select_tob  = ( output[i][2] == tob  || tob  == -1 ) ;
       if ( select_site && select_die && select_tob ) {
          sample.push_back( output[i] ) ;
       }
   }


   GetPitch( sample ) ;

   vector<swcs> xyV ;
   GetXY( sample, xyV ) ;

}

void Module::CompareXY( ) {

   int rSite(1), rDie(1), rTob(0) ;
   Input->GetParameters("RefSite",          &rSite );
   Input->GetParameters("RefDie",           &rDie );
   Input->GetParameters("RefTob",           &rTob );
   int szx = 17 ;
   int szy = 17 ;


   // Read data
   vector<vec> output ;
   //Match  *match  = new Match( ) ;
   // Sort Data by die address X and Y
   match->ReadStepFAST( output );
   printf(" Read Site %d Die %d Tob %d \n", rSite, rDie, rTob) ;

   vector<vec> sample ;
   vector<vec> ref ;

   for ( size_t i =0 ; i < output.size() ; i++ ) {
       if ( output[i][0] == rSite && output[i][1] == rDie && output[i][2] == rTob ) {
          ref.push_back( output[i] ) ;
       }
       if ( output[i][0] == site && output[i][1] == die && output[i][2] == tob ) {
          sample.push_back( output[i] ) ;
       }
   }

   printf(" ref size  %d , sample size %d \n", (int)ref.size() , (int)sample.size() ) ;
   vector<swcs> s1d1V ;
   GetXY( ref, s1d1V ) ;

   vector<swcs> sNdMV ;
   GetXY( sample, sNdMV ) ;

   double dxA[szx][szy], dyA[szx][szy] ;
   bool bA[szx][szy] ;
   int idx(0), idy(0) ;
   for ( size_t i=0 ; i< s1d1V.size() ; i++ ) {
       idx = s1d1V[i].i ;
       idy = s1d1V[i].j ;
       dxA[idx][idy] = sNdMV[i].x - s1d1V[i].x ;
       dyA[idx][idy] = sNdMV[i].y - s1d1V[i].y ;
       bA[idx][idy] = ( sNdMV[i].b == 0 || s1d1V[i].b == 0) ? false : true ;
       printf(" %.4f, %.4f, %.5f, %.4f, %.4f, %.5f \n", 
       	 	sNdMV[i].x, s1d1V[i].x, dxA[idx][idy], sNdMV[i].y, s1d1V[i].y, dyA[idx][idy] ) ;
   }

   // Print out difference
   FILE* logfile = fopen( logFileName3.c_str() ,"w");

   fprintf(logfile," dX \n" ) ;
   vector<double> dxV ;
   for ( int j=szy-1 ; j >=0  ; j-- ) {
       for ( int i= 0; i < szx; i++ ) {
	   if ( bA[i][j] ) {
	      fprintf( logfile, " %.5f,", dxA[i][j] ) ;
	      dxV.push_back( dxA[i][j] ) ;
	   } else {
	      fprintf( logfile, " %.5f,", 0.0 ) ;
	   }
       }
       pair<double,double> stats = gsl_vector_mean(dxV );
       fprintf(logfile, " %.5f, %.5f \n", stats.first, stats.second ) ;
       dxV.clear() ;
   }
   // Calculate mean and stdev w.r.t. each X
   vector<double> stdV ;
   for ( int i=0 ; i < szx ; i++ ) {
       for ( int j= 0; j < szy; j++ ) {
           if ( bA[i][j] ) dxV.push_back( dxA[i][j] ) ;
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

   fprintf(logfile,"\n dY \n" ) ;
   vector<double> dyV ;
   for ( int j=szy-1 ; j >=0  ; j-- ) {
       for ( int i= 0; i < szx; i++ ) {
	   if ( bA[i][j] ) {
	      fprintf( logfile, " %.5f,", dyA[i][j] ) ;
	      dyV.push_back( dyA[i][j] ) ;
	   } else {
	      fprintf( logfile, " %.5f,", 0.0 ) ;
	   }
       }
       pair<double,double> stats = gsl_vector_mean(dyV );
       fprintf(logfile, " %.5f, %.5f \n", stats.first, stats.second ) ;
   }
   // Calculate mean and stdev w.r.t. each X
   for ( int i=0 ; i < szx ; i++ ) {
       for ( int j= 0; j < szy; j++ ) {
           if ( bA[i][j] ) dyV.push_back( dyA[i][j] ) ;
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

   fclose(logfile) ;
}

// Only can look at one die and one site at a time 
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
                dy_ = sqrt( (yA[i][j+1] - yA[i][j])*(yA[i][j+1] - yA[i][j]) +  (xA[i][j+1] - xA[i][j])*(xA[i][j+1] - xA[i][j]) )  ;
                //dy_ = yA[i][j+1] - yA[i][j] ;
		dyA[i][j] = dy_ ;
	     } else  {
		dyA[i][j] = 0 ;
	     }

             if (bA[i+1][j] && bA[i][j] ) {
                dx_ = sqrt( (xA[i+1][j] - xA[i][j])*(xA[i+1][j] - xA[i][j]) + (yA[i+1][j] - yA[i][j])*(yA[i+1][j] - yA[i][j]) ) ;
                //dx_ = xA[i+1][j] - xA[i][j] ;
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

void Module::GetXY( vector<vec>& data, vector<swcs>& output ) {

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

	 swcs XYc ;
	 if ( data[i][3] == 0 ) {
            bA[ix][iy] = false ;
	    XYc.i = ix ;
	    XYc.j = iy ;
	    XYc.x = 0 ;
	    XYc.y = 0 ;
	    XYc.b = false ;
	 } else { 
            bA[ix][iy] = true  ;
	    XYc.i = ix ;
	    XYc.j = iy ;
	    XYc.x = xA[ix][iy] ;
	    XYc.y = yA[ix][iy] ;
	    XYc.b = true ;
	 }
	 output.push_back( XYc ) ;
	 //printf(" (%d,%d) = [ %4.5f, %4.5f ]\n", ix, iy, xA[ix][iy], yA[ix][iy] ) ;
     }

     // Print out X 
     fprintf(logfile," X \n" ) ;
     for ( int j=ysz-1 ; j >=0  ; j-- ) {
         for ( int i= 0; i < xsz; i++ ) {
             if (bA[i][j])   fprintf( logfile, " %.4f,", xA[i][j] ) ;
	     else            fprintf( logfile, " %.4f,", 0.0 ) ;
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

// Ratation angle : theta , orthogonality angle : phi
void Module::AxisCorrection( vector<vec>& dataV, double theta, double phi ) {

     double x(0), y(0), dx(0), dy(0) ;
     for (size_t i=0 ; i < dataV.size() ; i++ ) {
	
	  x = dataV[i][6] ;
         dx = dataV[i][8] ;
	  y = dataV[i][7] ;
	 dy = dataV[i][9] ;
         match->Rotation(  x,  y, -1*(theta) );
         match->Rotation( dx, dy, -1*(theta) );
         match->Orth(  x,  y, phi ) ;
         match->Orth( dx, dy, phi ) ;
	 dataV[i][6] =  x ;
	 dataV[i][7] =  y ;
	 dataV[i][8] = dx ;
	 dataV[i][9] = dy ;
     }

}

void Module::CenterCorrection( vector<vec>& wA, double& dX0, double& dY0 ) {

   // PCAD Positions
   double pCad0[2] = { -233.67, 11.9 } ;
   double pCad1[2] = {  214.33, 11.9 } ;
   double centerX = (pCad0[0] + pCad1[0]) / 2. ;
   double centerY = (pCad0[1] + pCad1[1]) / 2. ;
   double ff_X0 = (wA[0][1] + wA[0][3] + wA[1][1] + wA[1][3]) / 2. ;
   double ff_Y0 = (wA[0][2] + wA[0][4] + wA[1][2] + wA[1][4]) / 2. ;
   dX0 = centerX - ff_X0 ;
   dY0 = centerY - ff_Y0 ;

   printf(" PCAD0 : %.6f , %.6f \n ", pCad0[0], pCad0[1] ) ;
   printf(" PCAD1 : %.6f , %.6f \n ", pCad1[0], pCad1[1] ) ;
   printf(" FFLY0 : %.6f , %.6f \n ", wA[0][1] + wA[0][3], wA[0][2] + wA[0][4] ) ;
   printf(" FFLY1 : %.6f , %.6f \n ", wA[1][1] + wA[1][3], wA[1][2] + wA[1][4] ) ;
   printf(" dX = %.6f , dY = %.6f \n", dX0, dY0 ) ;
}

