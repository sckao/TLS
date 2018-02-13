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
  Input->GetParameters("Die",           &die ) ;
  Input->GetParameters("Site",          &site ) ;

  Input->GetParameters("debug",         &debug );

  logFileName = path +  outFileName ;
  mindX = 0.1 ;
}

Match::~Match(){

  //delete Input ;
    cout<<" done ! "<<endl ;
}

extern int is ;
extern bool ScanSort( vec s1, vec s2) { return (s1[is] < s2[is]) ; }

void Match::RepeatRun() {

   for (int r=1; r <11; r++ ) {
       for (int i=1; i < 5; i++ ) {
           for (int j=1; j<3; j++  ) {
               site = i ;
               die = j ;
 
               char buff[22] ;
               sprintf(buff, "FF_PSET_S%dD%d_L2_R%d.txt", site, die, r ) ;
	       logFileName = path +  buff ;

               char puff[14] ;
               sprintf(puff, "PSET-S%d-D%d.csv", site, die ) ;
               pFileName = puff ;

               char cuff[14] ;
               sprintf(cuff, "TLS_XCD_L2_%d.csv", r ) ;
               cFileName = cuff ;

               printf(" PSET:%s ,  XCD:%s , Output: %s\n", pFileName.c_str(), cFileName.c_str(), logFileName.c_str() ) ;

               OutputFakePSET() ;

           }
       }
   }
}

void Match::RepeatAna() {


   // Record the result
   FILE* logfile = fopen( logFileName.c_str() ,"w");

   vector<vec> psetV ;
   ReadPset( psetV ) ;
   int ax1 = 3 ;
   int ax2 = ax1 + 1 ;

   const int sz = (int) psetV.size() ;
   double xc[sz][10], yc[sz][10], xd[sz][10], yd[sz][10] ;

   for (int r=0; r <1; r++ ) {
 
	char cuff[14] ;
	sprintf(cuff, "Report_%d.csv", r ) ;
	cFileName = cuff ;

	vector<vec> dpmV ;
	//ReadDPM( dpmV ) ;
	cFileName = "Report_150118_054626.csv" ;
	ReadXCD( dpmV, 182 ) ;

	vector<vec> xcdV ;
	cFileName = "Report_150118_055105.csv" ;
	ReadXCD( xcdV, 75 ) ;

        
	//fprintf(logfile, " Report:%s  - pset size = %d, xcd size = %d, dpm size= %d\n", 
        //        cFileName.c_str(), (int)psetV.size(), (int)xcdV.size(), (int)dpmV.size() ) ;
        fprintf(logfile, "   CAD_X       CAD_Y      NAXCD_X     NAXCD_Y        XCD_X       XCD_Y\n" ) ;

        for (size_t i=0; i<psetV.size(); i++) { 

            xc[i][r] = xcdV[i][ax1] ;            
            yc[i][r] = xcdV[i][ax2] ;            
            xd[i][r] = dpmV[i][ax1] ;            
            yd[i][r] = dpmV[i][ax2] ; 
           
            //fprintf(logfile, " %04.6f  % 04.6f  %04.6f  %04.6f  %04.6f  %04.6f\n", 
            //        psetV[i][ax1], psetV[i][ax2], xcdV[i][ax1], xcdV[i][ax2], dpmV[i][ax1], dpmV[i][ax2] ) ;
            fprintf(logfile, " %04.6f  % 04.6f  %04.6f  %04.6f  %04.6f  %04.6f  %04.6f  %04.6f  %04.6f  %04.6f  %04.6f  %04.6f\n", 
                    psetV[i][ax1], psetV[i][ax2], xcdV[i][ax1], xcdV[i][ax2], dpmV[i][ax1], dpmV[i][ax2], xc[i][r] - xd[i][r], yc[i][r] - yd[i][r], xc[i][r]-psetV[i][ax1], yc[i][r]-psetV[i][ax2], xd[i][r]-psetV[i][ax1], yd[i][r]-psetV[i][ax2] ) ;

        }
   }

   fprintf(logfile, "\n\n === statistic ==== \n") ;
   /*
   for ( int i=0; i< sz; i++) {
       double xa[10],ya[10] ; 
       double xb[10],yb[10] ; 

       for ( int j=0; j< 10; j++) {
           xa[j] = xc[i][j] ;
           ya[j] = yc[i][j] ;
           xb[j] = xd[i][j] ;
           yb[j] = yd[i][j] ;
       }

       double mn_x_cd  = gsl_stats_mean( xa, 1, 10 );
       double st_x_cd   = gsl_stats_sd(xa, 1, 10 );
       double mn_y_cd  = gsl_stats_mean( ya, 1, 10 );
       double st_y_cd   = gsl_stats_sd(ya, 1, 10 );
   
       double mn_x_dp  = gsl_stats_mean( xb, 1, 10 );
       double st_x_dp   = gsl_stats_sd(xb, 1, 10 );
       double mn_y_dp  = gsl_stats_mean( yb, 1, 10 );
       double st_y_dp   = gsl_stats_sd(yb, 1, 10 );

       fprintf(logfile, " %04.6f  %04.6f  %04.6f  %04.6f,  %04.6f  %04.6f  %04.6f  %04.6f \n",
                         mn_x_cd, st_x_cd, mn_y_cd, st_y_cd, mn_x_dp, st_x_dp, mn_y_dp, st_y_dp ) ;
   }
   */
   fclose( logfile );

}


// Compare Stepper Actual to FF Actual
// Using DPM or XCD
void Match::OutputFakePSET( int mode ) {

    int xcdPos = 0 ;
    Input->GetParameters("XCDPos",         &xcdPos );

    vector<vec> psetV ;
    ReadPset( psetV ) ;

    vector<vec> ffV ;
    if ( mode == 0 ) ReadDPM( ffV ) ;
    if ( mode == 1 ) ReadXCD( ffV, xcdPos ) ;

    // Record the result
    FILE* logfile = fopen( logFileName.c_str() ,"w");

    int ax1 = 3 ;
    int ax2 = ax1 + 1 ;

    double x(0),y(0), dx(0),dy(0);
    fprintf(logfile, "Align Job, TLS-FFLY\nExpose Job, TLS-FFLY\nLot ID, GG978\nSubstrate ID, GG978-00\nDate, 2018-01-10 13:38:12\n") ;    
    fprintf(logfile, "Points,     % 03d\n", (int)psetV.size()+2 );
    for (int i=0; i< (int)psetV.size(); i++ ) {

         x = psetV[i][ax1] ; 
         y = psetV[i][ax2] ; 
         if ( mode == 0 ) {  // DPM
            //dx = ffV[i][5] ;
            //dy = ffV[i][6] ;
            dx = ffV[i][ax1] - x ; 
            dy = ffV[i][ax2] - y ; 
         }
         if ( mode == 1 ) {  // XCD
            dx = ffV[i][ax1] - x ; 
            dy = ffV[i][ax2] - y ; 
         }
        // PSET Format
        //fprintf(logfile, "%.6f	%.6f	%.6f	%.6f\n", x, y, dx, dy ) ;
        //fprintf(logfile, "%.6f	%.6f	%.6f %.6f	%.6f	%.6f\n", x, y, ffV[i][ax1], ffV[i][ax2], dx, dy ) ;
        // MAP Format 
        fprintf(logfile, "%3d,  % 04.6f, % 03.6f, % 04.3f, % 04.3f, M,    0.000,    0.000\n", 
                          i+2,        x,       y, dx*1000,  dy*1000   ) ;
        // Debug Format 
        //fprintf(logfile, " %04.6f  % 04.6f  %04.6f  %04.6f  %04.6f  %04.6f\n", x, y, ffV[i][ax1], ffV[i][ax2], dx, dy ) ;

    } 
    fclose( logfile );

}

void Match::Matching() {

    vector<vec> xcdV ;
    ReadXCD( xcdV ) ;

    vector<vec> psetV ;
    ReadPset( psetV ) ;

    // Record the result
    FILE* logfile = fopen( logFileName.c_str() ,"w");

    Minimum_GdR( logfile, psetV, xcdV ) ;
    fclose( logfile );

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
        if ( it->size() < 10 || (*it)[3] < 0 ) { 
           data.erase( it ) ; 
           //continue ;
        } else {
           it++ ;
        }
    }
    /*  
    printf(" size1 = %d \n", (int)data.size() ) ;
    for ( std::vector<vec>::iterator it = data.begin() ; it != data.end(); it++  ) {
        printf(" %.0f %.0f %.0f %.0f %.4f %.4f %.4f %.4f %.4f %.4f \n",
              (*it)[0], (*it)[1], (*it)[2], (*it)[3], (*it)[4], (*it)[5], (*it)[6], (*it)[7], (*it)[8],(*it)[9]  ) ;
    }
    */

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
               if ( cl%2 == 0 ) de = 1 ;
               if ( cl%2 == 1 ) de = 2 ;

               if ( j%4  < 2 && cl%4  < 2 ) st = 3 ;
               if ( j%4 >= 2 && cl%4  < 2 ) st = 1 ;
               if ( j%4  < 2 && cl%4 >= 2 ) st = 4 ;
               if ( j%4 >= 2 && cl%4 >= 2 ) st = 2 ;

               // Only turn on this if measuring all the sites and dies
               //if ( de != die || site != st ) continue; 
               //printf(" %d-%d-%d \n", cl, st, de ) ;

               //printf("[%d, %d] %.5f  %.5f \n", (int)j, cl, col[j][axis1], col[j][axis2] ) ;
               mdata.push_back( (double)cl ) ;
               mdata.push_back( (double)st ) ;
               mdata.push_back( (double)de ) ;
               mdata.push_back( col[j][axis1] ) ;
               mdata.push_back( col[j][axis2] ) ;
               mdata.push_back( col[j][8] ) ;  // dx
               mdata.push_back( col[j][9] ) ;  // dy

               data1.push_back( mdata ) ;
               mdata.clear() ;
           } 

           flush = false ;
	   val = data[i][axis1];
	   col.clear() ;
	   col.push_back(data[i] ) ;
	   cl++ ;
        }
    } 

    printf(" DPM Size = %d - %d \n", (int)data1.size(), (int)data1[0].size() ) ;

}

// Output id, site, die, x, y
void Match::ReadXCD( vector<vec>& data1, int nSkip  ) {

    data1.clear() ;

    vector<vec> data ;
    data.clear() ;

    string csvfile = cfolder + cFileName ;
    printf(" CSV file = %s \n", csvfile.c_str() ) ;
    Input->ReadFreeCSV( csvfile, data, nSkip ) ;

    for ( std::vector<vec>::iterator it = data.begin() ; it != data.end();  ) {
        if ( it->size() != 4  ) { 
           data.erase( it ) ; 
           //continue ;
        } else {
           it++ ;
        }
    }

    int axis1 = 1 ;
    int axis2 = ( axis1 == 2 ) ? 1 : 2  ;

    if ( data.size() > 1 ) {
        is = axis1 ;
        sort( data.begin(), data.end(), ScanSort );
	//for (size_t k =0 ; k< data.size(); k++ ) {
        //    printf(" %.5f  %.5f \n", data[k][axis1], data[k][axis2] ) ;
	//}
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

               // Giving die and site number
               if ( cl%2 == 0 ) de = 1 ;
               if ( cl%2 == 1 ) de = 2 ;

               if (  j%4  < 2 && cl%4  < 2 ) st = 3 ;
               if (  j%4 >= 2 && cl%4  < 2 ) st = 1 ;
               if (  j%4  < 2 && cl%4 >= 2 ) st = 4 ;
               if (  j%4 >= 2 && cl%4 >= 2 ) st = 2 ;

               //printf(" %d-%d-%d \n", cl, st, de ) ;
               //if ( de != die || site != st ) continue; 

               //printf("[%d, %d] %.5f  %.5f \n", (int)j, cl, col[j][axis1], col[j][axis2] ) ;
               mdata.push_back( (double)cl ) ;
               mdata.push_back( (double)st ) ;
               mdata.push_back( (double)de ) ;
               mdata.push_back( col[j][axis1] ) ;
               mdata.push_back( col[j][axis2] ) ;

               data1.push_back( mdata ) ;
               mdata.clear() ;
           } 

           flush = false ;
	   val = data[i][axis1];
	   col.clear() ;
	   col.push_back(data[i] ) ;
	   cl++ ;
        }
    } 

    printf(" XCD Size = %d - %d \n", (int)data1.size(), (int)data1[0].size() ) ;

}

// Output id, dx, dy, x, y
void Match::ReadPset( vector<vec>& data1 ) {

    data1.clear() ;

    vector<vec> data ;
    data.clear() ;

    string csvfile = cfolder + pFileName ;
    printf(" PSET file = %s \n", csvfile.c_str() ) ;
    Input->ReadFreeCSV( csvfile, data, 1 ) ;

    int axis1 = 0 ;
    int axis2 = ( axis1 == 1 ) ? 0 : 1  ;

    if ( data.size() > 1 ) {
        is = axis1 ;
        sort( data.begin(), data.end(), ScanSort );
    }

    vector<vec> col ;
    vec mdata ;
    double val = data[0][axis1] ;
    bool flush = false ;
    is = axis2 ;
    int cl(0) ;

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
               mdata.push_back( (double)cl ) ;
               mdata.push_back( col[j][2] ) ;
               mdata.push_back( col[j][3] ) ;
               mdata.push_back( col[j][axis1] ) ;
               mdata.push_back( col[j][axis2] ) ;

               data1.push_back( mdata ) ;
               mdata.clear() ;
           } 

           flush = false ;
	   val = data[i][axis1];
	   col.clear() ;
	   col.push_back(data[i] ) ;
	   cl++ ;
        }
    } 

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

void Match::Rotation( double& x, double& y, double theta ) {

     double x_ = x*cos(theta) - y*sin(theta) ;
     double y_ = x*sin(theta) + y*cos(theta) ;

     x = x_ ;
     y = y_ ;
}


