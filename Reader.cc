#include "Reader.h"
#include <gsl/gsl_fit.h>
#include <gsl/gsl_statistics.h>

//static bool dIncreasing( df s1, df s2) { return ( s1.s < s2.s ) ; }
static bool pIncreasing( vec s1, vec s2) { 
 
  if ( s1[2] < s2[2] ) {
     return true ;
  } 
  else if ( s1[2] == s2[2] ) {
     return ( s1[1] < s2[1] ) ; 
  } else {
    return false ;
  } 
}

static bool ScanSortXP( vec s1, vec s2) { return ( s1[4] < s2[4] ) ; }
static bool ScanSortXN( vec s1, vec s2) { return ( s1[4] > s2[4] ) ; }

static bool ScanSortY( vec s1, vec s2) { 
   if ( ( s1[5] - s2[5]) > 2. ) {
      return true ;
   } else {
      return false ;
   }
}

// sort x from small to large
static bool ScanSortYP( vec s1, vec s2) { return ( s1[5] < s2[5] ) ; }
static bool ScanSortYN( vec s1, vec s2) { return ( s1[5] > s2[5] ) ; }
static bool ScanSortX( vec s1, vec s2) { 
   if (  s1[4] < s2[4] ) {
      return true ;
   } else {
      return false ;
   }
}

int is = 3 ;
static bool ScanSort( vec s1, vec s2) { return (s1[is] < s2[is]) ; }   
 

Reader::Reader( ) {

  Input = AnaInput::Instance() ;

  Input->GetParameters("PlotType",      &plotType ) ; 
  Input->GetParameters("Path",          &hfolder ) ; 
  Input->GetParameters("CSVDIR",        &cfolder ) ; 
  Input->GetParameters("CSVFile",       &cFileName );

  Input->GetParameters("MapFile",       &mapFileName );
  Input->GetParameters("MapCorrFile",   &mapCorrFileName );

  Input->GetParameters("OutputFile",    &outFileName );
  Input->GetParameters("HFileName",     &hFileName );
  Input->GetParameters("debug",         &debug );

  Input->GetParameters("PlotName0",     &plotName0 );
  Input->GetParameters("PlotName1",     &plotName1 );

}

Reader::~Reader(){

  //delete Input ;
  cout<<" done ! "<<endl ;
}

void Reader::GetDataFromMap( ) {

    vector<vec> data ;
    data.clear() ;

    // Read data from Map file
    string kfile = cfolder + mapCorrFileName ;
    printf(" file = %s \n", kfile.c_str() ) ;
    //Input->ReadMap( kfile, data ) ;
    //Input->ReadNSX( kfile, data ) ;
    //Input->ReadDPM( kfile, data ) ;
    Input->ReadMapCorrection( kfile ,data ) ;
    cout<<" Read MapCorr File !! [ "<< data.size()<<" ]" <<endl ;

    // Record the result
    string logFileName = hfolder +  outFileName ;
    FILE* logfile = fopen( logFileName.c_str() ,"w");

    // Sort the data - Y first then X , small to large
    //if ( data.size() > 1 ) { 
    //   sort( data.begin(), data.end(), pIncreasing );
       //printf(" small d = %.1f ~ %.1f \n", data[0], data[ data.size() - 1 ] ) ;
    //}
    for ( size_t j=0 ; j< data.size(); j++ ) {
        fprintf(logfile, "%.0f,  %.4f,  %.4f,  %.4f,  %.4f,  %.9f\n", 
                           data[j][0], data[j][1], data[j][2], data[j][3], data[j][4], data[j][5] ) ;


    }
    fclose( logfile ) ;
 
}

// Sort Y then X 
void Reader::GetDataFromNSX( ) {

    Y2X = false ;
    vector<vec> data ;
    data.clear() ;

    // Read data from Map file
    string kfile = cfolder + cFileName ;
    printf(" file = %s \n", kfile.c_str() ) ;
    Input->ReadNSX( kfile, data ) ;

    // Record the result
    string logFileName = hfolder +  outFileName ;
    FILE* logfile = fopen( logFileName.c_str() ,"w");
    vector<vec> data1 ;

    if ( Y2X ) {
       // Sort the data - Y first then X , small to large
       if ( data.size() > 1 ) {
	  cout<<" data length = "<< data[ data.size() -1 ].size() <<endl ;
	  sort( data.begin(), data.end(), ScanSortY );
       }
       /*
       for (size_t i=0 ; i < data.size() ; i++ ) {
           printf(" == %.4f  %.4f \n", data[i][4], data[i][5] ) ;
        }
        cout<<" =========== End of sorting Y =========== "<<endl ;
        */

       vector<vec> vecY ;

       double nomY = data[0][5] ;
       bool scanP = true ;
       for (size_t i=0 ; i < data.size() ; i++ ) {
           //printf(" == %.4f  %.4f \n", data[i][4], data[i][5] ) ;
           if ( nomY - data[i][5] < 1.5 && i < data.size()-2 ) { 
              vecY.push_back( data[i] ) ; 
           } else {

              if ( i == data.size()- 1 ) vecY.push_back( data[i] ) ;
              if ( scanP ) { 
                 sort( vecY.begin(), vecY.end(), ScanSortXP );
              } else {
                 sort( vecY.begin(), vecY.end(), ScanSortXN );
              }
              for ( size_t j=0 ; j < vecY.size() ; j++ ) { 
                  data1.push_back( vecY[j] ) ;
		  //printf(" == %.4f  %.4f \n", vecY[j][4], vecY[j][5] ) ;
		  fprintf(logfile, " %.4f  %.4f \n",  vecY[j][4], vecY[j][5] ) ;
              }
              vecY.clear() ;
              scanP = !scanP ;
              nomY = data[i][5] ;
              vecY.push_back( data[i] ) ; 
           }
       }
    } else {

       // Sort the data - X first then Y , small to large
       if ( data.size() > 1 ) {
	  cout<<" data length = "<< data[ data.size() -1 ].size() <<endl ;
	  sort( data.begin(), data.end(), ScanSortX );
       }
       /*for (size_t i=0 ; i < data.size() ; i++ ) {
           printf(" == %.4f  %.4f \n", data[i][4], data[i][5] ) ;
       }
        cout<<" =========== End of sorting X =========== "<<endl ;
       */
 
       vector<vec> vecX ;

       double nomX = data[0][4] ;
       bool scanP = false ;
       for (size_t i=0 ; i < data.size() ; i++ ) {
           //printf(" == %.4f  %.4f \n", data[i][4], data[i][5] ) ;
           if ( fabs(nomX - data[i][4]) < 1.5 && i < (data.size()-2) ) { 
              vecX.push_back( data[i] ) ;
              //printf(" --> %.4f  %.4f %.4f \n", nomX, data[i][4], nomX - data[i][4]  ) ;
             
           } else {

              //cout<<" vecX size = "<< vecX.size() <<endl ; 
              if ( i == data.size()- 1 ) vecX.push_back( data[i] ) ;

              if ( scanP ) { 
                 sort( vecX.begin(), vecX.end(), ScanSortYP );
              } else {
                 sort( vecX.begin(), vecX.end(), ScanSortYN );
              }
              for ( size_t j=0 ; j < vecX.size() ; j++ ) { 
                  data1.push_back( vecX[j] ) ;
		  //printf(" == %.4f  %.4f \n", vecX[j][4], vecX[j][5] ) ;
		  fprintf(logfile, " %.4f  %.4f \n",  vecX[j][4], vecX[j][5] ) ;
              }
              vecX.clear() ;
              scanP = !scanP ;
              nomX = data[i][4] ;
              vecX.push_back( data[i] ) ; 
           }
       }
    }

    fclose( logfile ) ;
    // Check the sortitng 
    if (debug ) {
       for ( size_t i=0 ; i < data.size() ; i++ ) {
           int match = 0 ;
           for ( size_t j=0 ; j < data1.size() ; j++ ) {
               if ( data1[j] == data[i] ) match += 1 ; 
           }
           if ( match != 1 ) {
               printf(" >>>  %.4f  %.4f match %d times !! \n", data[i][4], data[i][5], match ) ;
           }
       }
    }

}

// Read MapCorrection file and match with NSX CDM
void Reader::GetDataFromSiteCorrectionMap( vector<vec>& outputV ) {

    /// 1 Read data from Map Correction file
    string kfile = cfolder + mapCorrFileName ;
    printf(" Map Corr file = %s \n", kfile.c_str() ) ;

    // 1.1 Collect shot data
    vector<vec> data ;
    data.clear() ;
    Input->ReadMapCorrection( kfile ,data ) ;

    /// 2 Read data from Map file
    string mfile = cfolder + mapFileName ;
    printf(" Map file = %s \n", kfile.c_str() ) ;
    // 2.1 Collect data
    vector<vec> mdata ;
    mdata.clear() ;
    Input->ReadMap( mfile, mdata ) ;

    /// 3 Read data from NSX/FF csv file
    string xfile = cfolder + cFileName ;
    printf(" NSX ile = %s \n", xfile.c_str() ) ;
    // 3.1 Collect data
    vector<vec> xdata ;
    xdata.clear() ;
    Input->ReadNSX( xfile, xdata ) ;
    

    // Record the result
    string logFileName = hfolder +  outFileName ;
    FILE* logfile = fopen( logFileName.c_str() ,"w");
    vec output_i ;

    double mX1(0), mY1(0), mX2(0), mY2(0) ;
    fprintf(logfile, "Site Xcoord  Ycoord    offsetX    offsetY     offsetT      dx1     dy1     dx2     dy2        ds1     dx1     dy1       ds2      dx2     dy2\n" ) ; 
    for ( size_t j= 0 ; j< data.size(); j++ ) {

        //fprintf(logfile, " %.0f  %.1f  %.1f  %.9f  %.9f  %.9f \n", 
        //             data[j][0], data[j][1], data[j][2], data[j][3], data[j][4], data[j][5] ) ;

        // Find the two measurement points
        mX1 =  data[j][1] - 5.  ;
        mY1 =  data[j][2] - 5.  ;
        mX2 =  data[j][1] + 5.  ;
        mY2 =  data[j][2] + 5.  ;

        // Match two meas. points of the shot
        int k1(-1), k2(-1) ;
        for ( size_t k=0; k< mdata.size(); k++ ) {
            if ( mdata[k][1] == mX1 && mdata[k][2] == mY1 )  k1 = (int)k ;
            if ( mdata[k][1] == mX2 && mdata[k][2] == mY2 )  k2 = (int)k ;
            if ( k1 > 0 && k2 > 0 ) break ;
        }

        // Match two meas. points of the shot
        int f1(-1), f2(-1) ;
        for ( size_t k=0; k< xdata.size(); k++ ) {
            if ( fabs(xdata[k][4] - mX1) < 2.  && fabs(xdata[k][5] - mY1) < 2. )  f1 = (int)k ;
            if ( fabs(xdata[k][4] - mX2) < 2.  && fabs(xdata[k][5] - mY2) < 2. )  f2 = (int)k ;
            if ( f1 > 0 && f2 > 0 ) break ;
        }

        // Output data
        fprintf(logfile, "%.0f  %.1f  %.1f  %.9f  %.9f  %.9f  %.3f  %.3f  %.3f  %.3f  %.4f %.4f %.4f %.4f %.4f %.4f\n", 
                           data[j][0], data[j][1], data[j][2], data[j][3], data[j][4], data[j][5], 
                           mdata[k1][3], mdata[k1][4], mdata[k2][3], mdata[k2][4],
                           xdata[f1][3], xdata[f1][4], xdata[f1][5], xdata[f2][3], xdata[f2][4], xdata[f2][5] ) ;
        
        output_i.clear() ;
        output_i.push_back( data[j][1] ) ;
        output_i.push_back( data[j][2] ) ;
        output_i.push_back( data[j][3] ) ;
        output_i.push_back( data[j][4] ) ;
        output_i.push_back( data[j][5] ) ;
        output_i.push_back( mdata[k1][3] ) ;
        output_i.push_back( mdata[k1][4] ) ;
        output_i.push_back( mdata[k2][3] ) ;
        output_i.push_back( mdata[k2][4] ) ;
        output_i.push_back( xdata[f1][3] ) ;
        output_i.push_back( xdata[f1][4] ) ;
        output_i.push_back( xdata[f1][5] ) ;
        output_i.push_back( xdata[f2][3] ) ;
        output_i.push_back( xdata[f2][4] ) ;
        output_i.push_back( xdata[f2][5] ) ;

        if ( output_i.size() != 15 ) {
            printf(" !!! scan fail (%d) \n", (int)output_i.size() ) ;
            break;
        }
        outputV.push_back( output_i ) ;
    }
    fclose( logfile ) ;
    

} 


void Reader::MatchDPM_SiteCorrMap( vector<vec>& outputV ) {

    /// 1 Read data from Map Correction file
    string kfile = cfolder + mapCorrFileName ;
    printf(" Map Corr file = %s \n", kfile.c_str() ) ;
    // 1.1 Collect shot data
    vector<vec> data ;
    data.clear() ;
    Input->ReadMapCorrection( kfile ,data ) ;
    cout<<" Read MapCorr File !! [ "<< data.size()<<" ]" <<endl ;

    /// 2 Read data from Map file
    string mfile = cfolder + mapFileName ;
    printf(" Map file = %s \n", mfile.c_str() ) ;
    // 2.1 Collect data
    vector<vec> mdata ;
    mdata.clear() ;
    //Input->ReadMap( mfile, mdata ) ;
    Input->ReadFreeCSV( mfile, mdata, 3 ) ;
    cout<<" Read PSET File !! [ "<< mdata.size() <<" ]" <<endl ;

    /// 3 Read data from NSX/FF DPM file
    string xfile = cfolder + cFileName ;
    printf(" DPM file = %s \n", xfile.c_str() ) ;
    // 3.1 Collect data
    vector<vec> xdata ;
    xdata.clear() ;
    //Input->ReadDPM( xfile, xdata ) ;
    Input->ReadFreeCSV( xfile, xdata, 1 ) ;
    cout<<" Read DPM File !! [ "<< xdata.size() <<" ]" <<endl ;
    
    //4  Record the result
    string logFileName = hfolder +  outFileName ;
    FILE* logfile = fopen( logFileName.c_str() ,"w");

    vec output_i ;

    //double mX1(0), mY1(0), mX2(0), mY2(0) ;
    fprintf(logfile, "Site Xcoord  Ycoord    offsetX    offsetY     offsetT      dx1     dy1     dx2     dy2       xc        yc     dA   die_dX   die_dY   die_dR\n" ) ; 
    for ( size_t j= 0 ; j< data.size(); j++ ) {

        //fprintf(logfile, " %.0f  %.1f  %.1f  %.9f  %.9f  %.9f \n", 
        //             data[j][0], data[j][1], data[j][2], data[j][3], data[j][4], data[j][5] ) ;

        // Find the two measurement points
        //mX1 =  data[j][1] - 5.  ;
        //mY1 =  data[j][2] - 5.  ;
        //mX2 =  data[j][1] + 5.  ;
        //mY2 =  data[j][2] + 5.  ;
        /*
        mX1 =  data[j][1] - 2.775  ;
        mY1 =  data[j][2] - 2.775  ;
        mX2 =  data[j][1] + 2.775  ;
        mY2 =  data[j][2] + 2.775  ;
        */

        // Match two meas. points of the shot
        int k1(-1), k2(-1) ;
        for ( size_t k=0; k< mdata.size(); k++ ) {
            //printf("%d) %.4f %.4f  - %.4f %.4f = %.4f %.4f \n", 
            //       (int)k, data[j][1], data[j][2], mdata[k][1], mdata[k][2], data[j][1] - mdata[k][1], data[j][2]-mdata[k][2] ) ;
            //if ( mdata[k][1] == mX1 && mdata[k][2] == mY1 )  k1 = (int)k ;
            //if ( mdata[k][1] == mX2 && mdata[k][2] == mY2 )  k2 = (int)k ;
            double dx_ = data[j][1] - mdata[k][1] ;
            double dy_ = data[j][2] - mdata[k][2] ;
            if ( fabs(dx_) < 2.776 && fabs(dy_) < 2.776 && dx_ > 0 && dy_ > 0 )  k1 = (int)k ;
            if ( fabs(dx_) < 2.776 && fabs(dy_) < 2.776 && dx_ < 0 && dy_ < 0 )  k2 = (int)k ;
            if ( k1 > 0 && k2 > 0 ) break ;
        }
        printf(" k = %d %d \n", k1, k2 ) ;

        // Match DPM point
        int f1(-1) ;
        for ( size_t k=0; k< xdata.size(); k++ ) {
            printf("%d) %.4f %.4f  - %.4f %.4f = %.4f %.4f \n", 
                   (int)k, data[j][1], data[j][2], xdata[k][2]/1000, xdata[k][3]/1000, data[j][1] - (xdata[k][2]/1000), data[j][2] - (xdata[k][3]/1000) ) ;
            //if ( fabs(xdata[k][7] - data[j][1]) < 2.  && fabs(xdata[k][8] - data[j][2]) < 2. )  f1 = (int)k ;
            if ( fabs(xdata[k][2]/1000 - data[j][1]) < 3.  && fabs(xdata[k][3]/1000 - data[j][2]) < 3. )  f1 = (int)k ;
            if ( f1 > 0 ) break ;
        }
        printf(" f1 = %d \n", f1 ) ;

        // Output data
        output_i.clear() ;
        fprintf(logfile, "%.0f  %.1f  %.1f  %.9f  %.9f  %.9f  %.3f  %.3f  %.3f  %.3f  %.4f %.4f %.4f %.4f %.4f %.4f\n", 
                           data[j][0], data[j][1], data[j][2], data[j][3], data[j][4], data[j][5], 
                           mdata[k1][3], mdata[k1][4], mdata[k2][3], mdata[k2][4],
                           //xdata[f1][7], xdata[f1][8], xdata[f1][6], xdata[f1][3], xdata[f1][4], xdata[f1][5] ) ;
                           xdata[f1][2]/1000, xdata[f1][3]/1000, xdata[f1][6], xdata[f1][4], xdata[f1][5], xdata[f1][7] ) ;

        output_i.push_back( data[j][1] ) ;
        output_i.push_back( data[j][2] ) ;
        output_i.push_back( data[j][3] ) ;
        output_i.push_back( data[j][4] ) ;
        output_i.push_back( data[j][5] ) ;
        output_i.push_back( mdata[k1][3] ) ;
        output_i.push_back( mdata[k1][4] ) ;
        output_i.push_back( mdata[k2][3] ) ;
        output_i.push_back( mdata[k2][4] ) ;
        output_i.push_back( xdata[f1][2] ) ;
        output_i.push_back( xdata[f1][3] ) ;
        output_i.push_back( xdata[f1][6] ) ;
        output_i.push_back( xdata[f1][4] ) ;
        output_i.push_back( xdata[f1][5] ) ;
        output_i.push_back( xdata[f1][7] ) ;
        //output_i.push_back( xdata[f1][7] ) ;
        //output_i.push_back( xdata[f1][8] ) ;
        //output_i.push_back( xdata[f1][6] ) ;
        //output_i.push_back( xdata[f1][3] ) ;
        //output_i.push_back( xdata[f1][4] ) ;
        //output_i.push_back( xdata[f1][5] ) ;
       
        if ( output_i.size() != 15 ) {
            printf(" !!! scan fail (%d) \n", (int)output_i.size() ) ;
            break;
        }
        outputV.push_back( output_i ) ;
    }

    fclose( logfile ) ;
    
} 

void Reader::MatchDPM_MapCorr( vector<vec>& outputV ) {

    /// 1 Read data from Map Correction file
    string kfile = cfolder + mapCorrFileName ;
    printf(" Map Corr file = %s \n", kfile.c_str() ) ;
    // 1.1 Collect shot data
    vector<vec> data ;
    data.clear() ;
    Input->ReadMapCorrection( kfile ,data ) ;
    cout<<" Read MapCorr File !! [ "<< data.size()<<" ]" <<endl ;

    /// 2 Read data from NSX/FF DPM file
    string xfile = cfolder + cFileName ;
    printf(" DPM file = %s \n", xfile.c_str() ) ;
    // 2.1 Collect data
    vector<vec> xdata ;
    xdata.clear() ;
    Input->ReadFreeCSV( xfile, xdata, 1 ) ;
    cout<<" Read DPM File !! [ "<< xdata.size() <<" ]" <<endl ;
    
    //3  Record the result
    string logFileName = hfolder +  outFileName ;
    FILE* logfile = fopen( logFileName.c_str() ,"w");

    vec output_i ;

    fprintf(logfile, "Site Xcoord  Ycoord    offsetX    offsetY     offsetT      dx1     dy1     dx2     dy2       xc        yc     dA   die_dX   die_dY   die_dR\n" ) ; 
    for ( size_t j= 0 ; j< data.size(); j++ ) {

        //fprintf(logfile, " %.0f  %.1f  %.1f  %.9f  %.9f  %.9f \n", 
        //             data[j][0], data[j][1], data[j][2], data[j][3], data[j][4], data[j][5] ) ;

        // Match DPM point
        int f1(-1) ;
        for ( size_t k=0; k< xdata.size(); k++ ) {
            printf("%d) %.4f %.4f  - %.4f %.4f = %.4f %.4f \n", 
                   (int)k, data[j][1], data[j][2], xdata[k][2]/1000, xdata[k][3]/1000, data[j][1] - (xdata[k][2]/1000), data[j][2] - (xdata[k][3]/1000) ) ;
            if ( fabs(xdata[k][2]/1000 - data[j][1]) < 5.  && fabs(xdata[k][3]/1000 - data[j][2]) < 5. )  f1 = (int)k ;
            if ( f1 > 0 ) break ;
        }
        printf(" f1 = %d \n", f1 ) ;

        // Output data
        output_i.clear() ;
        fprintf(logfile, "%.0f  %.1f  %.1f  %.9f  %.9f  %.9f  %.4f %.4f %.4f %.4f %.4f %.4f\n", 
                           data[j][0], data[j][1], data[j][2], data[j][3], data[j][4], data[j][5], 
                           xdata[f1][2]/1000, xdata[f1][3]/1000, xdata[f1][6], xdata[f1][4], xdata[f1][5], xdata[f1][7] ) ;

        output_i.push_back( data[j][1] ) ;
        output_i.push_back( data[j][2] ) ;
        output_i.push_back( data[j][3] ) ;
        output_i.push_back( data[j][4] ) ;
        output_i.push_back( data[j][5] ) ;
        output_i.push_back( xdata[f1][2] ) ;
        output_i.push_back( xdata[f1][3] ) ;
        output_i.push_back( xdata[f1][6] ) ;
        output_i.push_back( xdata[f1][4] ) ;
        output_i.push_back( xdata[f1][5] ) ;
        output_i.push_back( xdata[f1][7] ) ;
       
        if ( output_i.size() != 11 ) {
            printf(" !!! scan fail (%d) \n", (int)output_i.size() ) ;
            break;
        }
        outputV.push_back( output_i ) ;
    }

    fclose( logfile ) ;
    
} 


void Reader::GetDataFromMapCorrection( ) {

    /// 1 Read data from Map Correction file
    string kfile = cfolder + mapCorrFileName ;
    printf(" Map Corr file = %s \n", kfile.c_str() ) ;
    // 1 Collect shot data
    vector<vec> data ;
    data.clear() ;
    Input->ReadMapCorrection( kfile ,data ) ;
    cout<<" Read MapCorr File !! [ "<< data.size()<<" ]" <<endl ;

    //2  Record the result
    string logFileName = hfolder +  outFileName ;
    FILE* logfile = fopen( logFileName.c_str() ,"w");

    vec output_i ;

    fprintf(logfile, "Site  X      Y      offsetX       offsetY      offset_Theta\n" ) ; 
    for ( size_t j= 0 ; j< data.size(); j++ ) {

        fprintf(logfile, "%.0f, % 04.2f, % 04.2f, %.9f, %.9f, %.9f\n", 
                     data[j][0], data[j][1], data[j][2], data[j][3], data[j][4], data[j][5] ) ;


    }

    fclose( logfile ) ;
    
} 



void Reader::GetDataFromSINF( ) {

    string sfolder, sinfFile1, sinfFile2 ;
    double gd1, gd2 ;
    vector<int> ad1, ad2, ad2_xf ;
    vector<double> dfcodes ;
    Input->GetParameters("SINFDIR",        &sfolder ) ; 
    Input->GetParameters("SINFFile1",      &sinfFile1 );
    Input->GetParameters("SINFFile2",      &sinfFile2 );
    Input->GetParameters("GoodDie1",       &gd1 );
    Input->GetParameters("GoodDie2",       &gd2 );
    Input->GetParameters("AlignDie1",      &ad1 );
    Input->GetParameters("AlignDie2",      &ad2 );
    Input->GetParameters("AlignDie2_XF",   &ad2_xf );
    Input->GetParameters("DFCodes",        &dfcodes );
    int di = ad1[1] - ad2[1] ;
    int dj = ad1[0] - ad2[0] ;
    int di_xf = ad2[1] - ad2_xf[1] ;
    int dj_xf = ad2[0] + ad2_xf[0] ;
 
    vector<vec> data1 ;
    data1.clear() ;
    vector<vec> data2 ;
    data2.clear() ;

    // Read data from SINF files
    string sfile1 = sfolder + sinfFile1 ;
    string sfile2 = sfolder + sinfFile2 ;
    printf(" sinf file1 = %s \n", sfile1.c_str() ) ;
    printf(" sinf file2 = %s \n", sfile2.c_str() ) ;

    Input->ReadSINF( sfile1, data1 );
    printf(" N1 of row = %d \n", (int)data1.size() ) ;
    printf(" N1 of col = %d \n", (int) data1[0].size() ) ;

    Input->ReadSINF( sfile2, data2 );
    printf(" N2 of row = %d \n", (int)data2.size() ) ;
    printf(" N2 of col = %d \n", (int) data2[0].size() ) ;

    // Check Reading
    bool ReadError = false ;
    if ( data1.size() > 0 && data2.size() > 0 ) {
       if ( data1.size() != data2.size() || data1[0].size() != data2[0].size() ) {
          printf(" Caution : Two SINF files have different numbers of row and column !!!\n" ) ;
          //ReadError = false ;
       }
    } else {
        printf(" Error : SINF files Reading Error !!!\n" ) ;
        ReadError = true ;
    }

    if ( ReadError ) return ;

    // Record the result
    string logFileName = hfolder +  outFileName ;
    FILE* logfile = fopen( logFileName.c_str() ,"w");

    // Defect Codes 
    double gdf  = dfcodes[0] ;
    double bdf  = dfcodes[1] ;
    double g1b2 = dfcodes[2] ;
    double g2b1 = dfcodes[3] ;
    double mm   = dfcodes[4] ;

    fprintf(logfile, "File 1 = %s \n", sfile1.c_str() );
    fprintf(logfile, "File 2 = %s \n", sfile2.c_str() );
    fprintf(logfile,  "Good Die Code : %02.0f , Bad Die Code : %02.0f, Mis-matched Code : %02.0f\n", gdf, bdf, mm ) ;
    fprintf(logfile,  "Good1 Bad2 : %02.0f , Good2 Bad1 : %02.0f \n\n", g1b2, g2b1 ) ;

    // Compare two files
    vector<vec> dataC ; 
    vec ctmp ;
    int nGD = 0 ;
    int nBD = 0 ;
    int nB1 = 0 ;
    int nB2 = 0 ;
    int nMM = 0 ;
    vector<addr> llist ;
    // Loop of row (Y-axis)
    for ( int j=0 ; j< (int)data1.size(); j++ ) {

        // Just print the die x-index
        if ( j- dj == -1 ) {
           fprintf(logfile, "    ") ; 
           for ( int k=0; k< (int)data1[j].size(); k++ ) {
               fprintf(logfile, " %2d", k-di+1 ) ; 
           }
           fprintf(logfile, "\n") ; 
        }

        if ( j - dj > (int)data2.size() - 1 || j  < dj  ) continue ;
        //printf("[%d] ", j ) ;
        // Print die y-index
        fprintf(logfile, "[%2d] ", 89-j ) ; 

        // Loop of column (X-axis)
        ctmp.clear() ;
        for ( int i=0; i< (int)data1[j].size(); i++ ) {

            if ( i - di > (int)data2[j-dj].size() - 1 || i < di ) continue ;
            //printf(" data 1-2 = %02.0f  - %02.0f \n", data1[j][i] ,  data2[j-dj][i-di]  ) ;
            if ( data2[j-dj][i-di] < 0 && data1[j][i] < 0 ) {
               ctmp.push_back( -1 ) ; 
               printf("-1 " ) ;
               fprintf(logfile, "-1 " ) ; 
            } else {
               if ( data1[j][i] == gd1 && data2[j-dj][i-di] == gd2 ) {
                  ctmp.push_back( gdf ) ;
                  printf("%02.0f ", gdf ) ;
                  fprintf(logfile, "%02.0f ", gdf ) ; 
                  nGD++ ;
               } else if ( (data1[j][i] != gd1) && (data2[j-dj][i-di] != gd2) && data1[j][i] > 0 && data2[j-dj][i-di] > 0 ) {
                  ctmp.push_back( bdf ) ;
                  printf("%02.0f ", bdf ) ;
                  fprintf(logfile, "%02.0f ", bdf ) ; 
                  nBD++ ;
               } else if ( (data1[j][i] != gd1) && (data2[j-dj][i-di] == gd2) && data1[j][i] > 0 && data2[j-dj][i-di] > 0 ) {
                  ctmp.push_back( g2b1 ) ;
                  printf("%02.0f ", g2b1 ) ;
                  fprintf(logfile, "%02.0f ", g2b1 ) ;
                  llist.push_back( make_pair(i - di - di_xf , dj_xf - j + dj ) )  ;
                  nB1++ ;
               } else if ( (data1[j][i] == gd1) && (data2[j-dj][i-di] != gd2) && data1[j][i] > 0 && data2[j-dj][i-di] > 0 ) {
                  ctmp.push_back( g1b2 ) ;
                  printf("%02.0f ", g1b2 ) ;
                  fprintf(logfile, "%02.0f ", g1b2 ) ; 
                  nB2++ ;
               } else {
                  ctmp.push_back( mm ) ;
                  printf("%02.0f ", mm ) ;
                  fprintf(logfile, "%02.0f ", mm ) ; 
                  nMM++ ;
               }

            }            
        }
        printf("\n") ;
        fprintf(logfile, "\n" ) ;
        dataC.push_back( ctmp ) ;
    }

    printf(" Nu of row = %d \n", (int)dataC.size() ) ;
    printf(" Nu of col = %d \n", (int)dataC[0].size() ) ;
    printf(" Nu of g2b1 = %d \n", (int)llist.size() ) ;
    fprintf(logfile, "\n\n N of Good Die = %d \n", nGD ) ;
    fprintf(logfile, " N of Bad Die = %d \n", nBD ) ;
    fprintf(logfile, " N of Good1 Bad2 = %d \n", nB2 ) ;
    fprintf(logfile, " N of Good2 Bad1 = %d \n", nB1 ) ;
    fprintf(logfile, " N of MisMatched Die = %d \n", nMM ) ;

    fprintf(logfile, "\n List of Good2 Bad1 dies address (x,y) \n") ;
    for (size_t i=0; i < llist.size(); i++ ) {
        printf(" [%d] ( %d, %d ) \n", (int)i+1, llist[i].first, llist[i].second ) ;
        fprintf(logfile, " [%d] ( %d, %d ) \n", (int)i+1, llist[i].first, llist[i].second ) ;
    }

    fclose( logfile ) ;

}

void Reader::GetDataFromCSV( ) {

    vector<vec> data ;
    data.clear() ;

    // Read data from CSV file
    string kfile = cfolder + cFileName ;
    Input->ReadCSV( kfile, data ) ;
    printf(" file = %s \n", kfile.c_str() ) ;

    // Record the result
    string logFileName = hfolder +  outFileName ;
    FILE* logfile = fopen( logFileName.c_str() ,"w");

    // Sort the data - Y first then X , small to large
    if ( data.size() > 1 ) { 
       sort( data.begin(), data.end(), pIncreasing );
       //printf(" small d = %.1f ~ %.1f \n", data[0], data[ data.size() - 1 ] ) ;
    }
 
    string plotname0_ = hfolder + plotName0 + "." + plotType ;
    MakePlots( data, plotname0_ ) ;

    LableDie( data ) ;
    DieRotation( data ) ;
    //sort( data.begin(), data.end(), pIncreasing );

    int row = 0 ; 
    for ( size_t j= 0 ; j< data.size(); j++ ) {
        //if (data[j].size() != 8 ) fprintf(logfile, " wrong entry !! \n") ;

        if ( j%132 == 0 ) row++ ;
        if ( j > 1 && data[j][1] == data[j-1][1] && data[j][2] == data[j-1][2] ) {
           if (data[j].size() != 8 ) fprintf(logfile, " duplicate entry !! \n") ;
        }
        int dieNu = (int) data[j][6] ;
        if ( row%2 == 1 ) continue ;
        if ( dieNu % 2 == 0 ) continue ;
        fprintf(logfile, " %.0f  %.6f  %.6f   %.6f  %.6f  %.1f  %.0f %.0f %.6f \n", 
                     data[j][0], data[j][1], data[j][2], data[j][3], data[j][4], data[j][5], data[j][6], data[j][7], data[j][8] ) ;
    }


    fclose( logfile ) ;

    string plotname1_ = hfolder + plotName1 + "." + plotType ;
    MakePlots( data, plotname1_ ) ;

}


void Reader::LableDie( vector<vec>& data ) {
    
    double goodM = 1 ; 
    double iDie = 0 ;
    //double theX = data[0][1] ;
    //double theY = data[0][2] ;
    int k = 0 ;

    for ( size_t j= 0 ; j< data.size() ; j++ ) {
        k++ ; 
        if ( k > 23000 ) { 
           cout<<" endless problem !! sz = "<< data.size() <<endl ;
           break ;
        }

        if ( j == 0 ) {
           if ( data[j][0] == 0 ) {
              iDie += 1 ;
              data[j].push_back( iDie ) ;  
              data[j].push_back( goodM ) ; 
           } else {
              iDie += 1 ;
              data[j].push_back( iDie ) ;  
              data[j].push_back( goodM ) ; 
              vec miss0 = data[j] ; // copy values from mark1 
              miss0[0] = 0 ;        // change mark to 1
	      miss0[1] = data[j][1] - 5.55 ;          // change x position 
	      miss0[7] = 0 ;        // mark it as a vitual mark   
	      data.insert( data.begin() ,miss0) ;    

           }
           continue ;
        }
        if ( j == data.size()-1 ) {
           if ( data[j][0] == 1 ) {
              data[j].push_back( iDie ) ;  
              data[j].push_back( goodM ) ; 
           } else {
              iDie += 1 ;
              data[j].push_back( iDie ) ;  
              data[j].push_back( goodM ) ; 
              vec miss1 = data[j] ; // copy values from mark1 
              miss1[0] = 0 ;        // change mark to 1
	      miss1[1] = data[j][1] + 5.55 ;          // change x position 
	      miss1[7] = 0 ;        // mark it as a vitual mark   
	      data.insert( data.end() ,miss1) ;    
           }
           continue ;
        }

        // Possible X nominal values : 5.5 or 6.9 or 1.35 or 22 or 16.45
        double dX_b =  data[j][1]   - data[j-1][1] ;          
	double dX_a =  data[j+1][1] - data[j][1] ;    
        // Possible Y nominal values : 6.9 or 18.2
        double dY_b =  data[j][2]   - data[j-1][2] ;
        double dY_a =  data[j+1][2] - data[j][2] ;
       
        if ( data[j][0] == 0 ) {
           if ( data[j].size() == 6 ) {
              iDie += 1 ;
	      data[j].push_back( iDie ) ;  
	      data[j].push_back( goodM ) ;
           } 
           // missing next mark1 
           if ( data[j+1][0] != 1 || fabs(dX_a - 5.55) > 0.2 || fabs(dY_a) > 0.5  ) {
              if ( data[j][1] > -243 && data[j][1] < 240  ) {
                 vec miss1 = data[j] ; // copy values from mark0 
		 miss1[0] = 1 ;        // change mark to 1
		 miss1[1] = data[j][1] + 5.55 ;          // change x position 
		 miss1[7] = 0 ;         // mark it as a vitual mark   
		 data.insert( data.begin()+j+1,miss1) ; 
              }
           }
           continue ; 
        }
     
        if ( data[j][0] == 1 ) { 
           if ( fabs(dX_b -5.55 ) < 0.2 ) {
              if ( data[j].size() == 6 ) {
                 data[j].push_back( iDie ) ;  
                 data[j].push_back( goodM ) ;  
              }
              // Edge of the shot
              if ( fabs( data[j][1] + 91.675 ) < 0.5  ||  fabs( data[j][1] - 75.225 ) < 0.5 )  {
                 if ( data[j+1][0] != 0 || fabs(dX_a - 16.45) > 0.2 ) {
                    iDie += 1 ;
		    vec miss0 = data[j] ; // copy values from mark1 
		    miss0[0] = 0 ;        // change mark to 0
		    miss0[1] = data[j][1] + 16.45 ;          // change x position 
		    miss0[6] = iDie ;        // next die, die id +1   
		    miss0[7] = 0 ;        // mark it as a vitual mark   
		    data.insert( data.begin()+j+1 ,miss0) ;  
                    j-- ;
                 }
              // missing next mark0 
              } else if ( data[j+1][0] != 0 || fabs(dX_a - 1.35) > 0.2 || fabs(dY_b) > 0.5  ) {
              //} else if ( data[j+1][0] != 0 || fabs(dX_a - 1.35) > 0.2 ) {
                 if ( data[j][1] > -243 && data[j][1] < 240  ) {
                    iDie +=1 ; 
		    vec miss0 = data[j] ; // copy values from mark1 
		    miss0[0] = 0 ;        // change mark to 1
		    miss0[1] = data[j][1] + 1.35 ;          // change x position 
		    miss0[6] = iDie ;        // next die, die id +1   
		    miss0[7] = 0 ;        // mark it as a vitual mark   
		    data.insert( data.begin()+j+1 ,miss0) ;    
		    j-- ; 
                 }  
              } 
              continue ;
           } else {
              // missing previous mark0
              if ( data[j].size() == 6 ) {
                 iDie += 1 ; 
		 data[j].push_back( iDie ) ;  
		 data[j].push_back( goodM ) ; 
		 data[j][6] = iDie ; 
              }
	      vec miss0 = data[j] ; // copy values from mark1 
	      miss0[0] = 0 ;        // change mark to 0
	      miss0[1] = data[j][1] - 5.55 ;          // change x position 
	      miss0[6] = iDie ;        // next die, die id +1   
	      miss0[7] = 0 ;        // mark it as a vitual mark   
	      data.insert( data.begin()+j ,miss0) ;
	      continue ;
           }

        }

    }
}


void Reader::DieRotation( vector<vec>& data ) {

    sort( data.begin(), data.end(), pIncreasing );

    double dx = 0 ;
    double dy = 0 ;
    double theta = -999999. ;
    for ( size_t j= 0 ; j< data.size(); j++ ) {

        if ( j == data.size()-1 ) {
           data[j].push_back( theta ) ;
           continue ;
        }

        if ( (int)data[j][0] == 0 && (int)data[j+1][0] == 1 && data[j][6] == data[j+1][6] ) {
           dx = data[j+1][1] - data[j][1] + data[j+1][3] - data[j][3] ;      
           dy = data[j+1][2] - data[j][2] + data[j+1][4] - data[j][4] ;   
           if ( dx != 0.0 ) theta = atan( dy/dx ) ;   
           data[j].push_back( theta ) ;
        } else {
           data[j].push_back( theta ) ;
        }
    }

}

pair<double,double> Reader::DieTranslation( vector<vec>& v1, vector<vec>& v2, int ix, int iy ) {
 
     double x1(0), y1(0) ; 
     for (size_t i=0 ; i< v1.size(); i++ ) {
         x1 += v1[i][ix] ;
         y1 += v1[i][iy] ;
     }
     x1 = x1 / (double)v1.size() ;
     y1 = y1 / (double)v1.size() ;

     double x2(0), y2(0) ; 
     for (size_t i=0 ; i< v2.size(); i++ ) {
         x2 += v2[i][ix] ;
         y2 += v2[i][iy] ;
     }
     x2 = x2 / (double)v2.size() ;
     y2 = y2 / (double)v2.size() ;

     double dx = x1 - x2 ;
     double dy = y1 - y2 ;

     pair<double,double> dxy = make_pair( dx,dy ) ;
     return dxy ;

}

void Reader::MakePlots( vector<vec> data, string plotname ) {

 
     // Open a root file to store graph
     string hfile = hfolder + hFileName + ".root" ;
     TFile* rootFile = new TFile( hfile.c_str() , "RECREATE" );
     rootFile->cd() ;

     int sz0_(0),sz1_(0), sz2_(0),sz3_(0);
     for ( size_t i =0; i < data.size() ; i++ ) { 
         if ( data[i].size() > 6 ) { 
            if ( data[i][0] == 0 && data[i][7] == 1 ) sz0_++ ;
	    if ( data[i][0] == 1 && data[i][7] == 1 ) sz1_++ ;
	    if ( data[i][0] == 0 && data[i][7] == 0 ) sz2_++ ;
	    if ( data[i][0] == 1 && data[i][7] == 0 ) sz3_++ ;
         } else {
            data[i].push_back(0) ;
            data[i].push_back(1) ;
            if ( data[i][0] == 0  ) sz0_++ ;
	    if ( data[i][0] == 1  ) sz1_++ ;
	    if ( data[i][0] == 0  ) sz2_++ ;
	    if ( data[i][0] == 1  ) sz3_++ ;
         }
     }

     const int sz0 = sz0_ ;
     const int sz1 = sz1_ ;
     const int sz2 = sz2_ ;
     const int sz3 = sz3_ ;
     double x0[sz0], y0[sz0], x1[sz1], y1[sz1]  ; 
     double x2[sz2], y2[sz2], x3[sz3], y3[sz3]  ; 
     cout<<" size0 = "<< sz0 <<" size1 = "<< sz1 <<endl ;
     int k(0),j(0),m(0),n(0) ;
     for ( size_t i =0; i < data.size() ; i++ ) {
          if ( data[i][0] == 0 && data[i][7] == 1 ) {
             x0[j] = data[i][1] ; 
             y0[j] = data[i][2] ;
             j++ ;
          } 
          if ( data[i][0] == 1 && data[i][7] == 1 ) {
             x1[k] = data[i][1] ; 
             y1[k] = data[i][2] ;
             k++ ;
          }
          if ( data[i][0] == 0 && data[i][7] == 0 ) {
             x2[m] = data[i][1] ; 
             y2[m] = data[i][2] ;
             m++ ;
          } 
          if ( data[i][0] == 1 && data[i][7] == 0 ) {
             x3[n] = data[i][1] ; 
             y3[n] = data[i][2] ;
             n++ ;
          }
     }

     TCanvas* c0 = new TCanvas("c_0","", 1200, 1200);
     c0->SetFillColor(10);
     c0->SetFillColor(10);
     //gPad->SetGridx();
     //gPad->SetGridy();
 
     TGraph* g0 = new TGraph( sz0, x0, y0 ) ; 
     g0->GetXaxis()->SetLimits(-260, 260);
     g0->SetName(" Mark0 ") ;
     g0->SetTitle("Mark0 Distribution") ;
     g0->SetMarkerStyle(20) ;
     g0->SetMarkerSize(0.3) ;
     g0->SetMarkerColor(1) ;
     g0->Draw("AP") ;
     c0->Update() ;
     g0->Write() ;

     TGraph* g1 = new TGraph( sz1, x1, y1 ) ; 
     g1->GetXaxis()->SetLimits(-260, 260);
     g1->SetName("Mark 1") ;
     g1->SetTitle("Mark1 Distribution") ;
     g1->SetMarkerStyle(21) ;
     g1->SetMarkerSize(0.3) ;
     g1->SetMarkerColor(4) ;
     g1->Draw("SAMEP") ;
     c0->Update() ;
     g1->Write() ;
    
     TGraph* g2 = new TGraph( sz2, x2, y2 ) ; 
     g2->SetMarkerStyle(20) ;
     g2->SetMarkerSize(0.3) ;
     g2->SetMarkerColor(2) ;
     g2->Draw("SAMEP") ;
     c0->Update() ;
     g2->Write() ;

     TGraph* g3 = new TGraph( sz3, x3, y3 ) ; 
     g3->SetMarkerStyle(21) ;
     g3->SetMarkerSize(0.3) ;
     g3->SetMarkerColor(6) ;
     g3->Draw("SAMEP") ;
     c0->Update() ;
     g3->Write() ;

     //plotname = hfolder + "mark1" + ".png" ;
     c0->Print( plotname.c_str() ) ;

     rootFile->Close() ;
 
     delete c0 ;
}


void Reader::GetDataFromASML( ) {

    vector<vec> data ;
    data.clear() ;

    // Read data from the CSV file
    string csvfile = cfolder + cFileName ;
    printf(" CSV file = %s \n", csvfile.c_str() ) ;
    // output a double vector of 5 elements
    // ID, XSize YSize, CenterX, CenterY
    Input->ReadFreeCSV( csvfile, data, 1 ) ;  

    int axis = 4 ;
    Input->GetParameters("Axis",         &axis );

    // Record the result
    string logFileName = hfolder +  outFileName ;
    FILE* logfile = fopen( logFileName.c_str() ,"w");

    // Sorting sequence X first (axis = 3) or Y first (axis = 4)
    // Only Axis1 need to be changed, Datacard output filename also need to be changed.
    int axis1 = axis ;
    int axis2 = ( axis1 == 4 ) ? 3 : 4  ;

    // sorting data on Axis1
    if ( data.size() > 1 ) { 
       is = axis1 ;
       sort( data.begin(), data.end(), ScanSort );
       //printf(" small d = %.1f ~ %.1f \n", data[0], data[ data.size() - 1 ] ) ;
    }
    printf(" data size = %d \n", (int)data.size() ) ;

    double val = data[0][axis1] ;
    vector<vec> col ;
    is = axis2  ;
    vector<vec> lineSum ;
    bool flush = false ; 
    for (size_t i=0; i < data.size(); i++ ) {
        if ( data[i][1] < 250. || data[i][2] < 250. ) continue ;


        if ( fabs(val - data[i][axis1] ) < 0.05 ) {
            col.push_back( data[i] ) ;
            if ( i == data.size() -1 ) flush = true ;
        } else {
            flush = true ;
        }

        if ( flush ) {
            // sorting data on X for each Y values
            sort( col.begin(), col.end(), ScanSort ) ;
            // For Fitting
            const int csz = (int) col.size() ;
            double xa[csz], ya[csz] ;
            double c0,c1,cov00, cov01, cov11, sumsq ;
            for ( size_t j=0; j < col.size(); j++ ) {
                xa[j] = col[j][axis2] ;
                ya[j] = col[j][axis1] ;

                if (j ==0 ) fprintf(logfile, "%d, %.4f, %.4f, 0 \n",    (int)j, col[j][3], col[j][4] ) ;
                else        fprintf(logfile, "%d, %.4f, %.4f, %.4f\n", (int)j, col[j][3], col[j][4], col[j][axis2] - col[j-1][axis2] ) ;
            }
            // Do some GSL functions
            gsl_fit_linear( xa, 1, ya, 1, csz, &c0, &c1, &cov00, &cov01, &cov11, &sumsq ) ;
            double mean  = gsl_stats_mean( ya, 1, csz );
            double std = gsl_stats_sd(ya, 1, csz );
            vec ls ;
            ls.push_back( mean ) ;
            ls.push_back( std ) ;
            ls.push_back( c0 ) ;
            ls.push_back( c1 ) ;
            ls.push_back( (double)csz ) ;
            //fprintf(logfile, "[%d]  Y= %.5f +/- %.5f %f  %f \n", 
            //                  csz,     mean,    std, c0, c1  ) ;
            fprintf(logfile, " --- \n" ) ;
 
            // Get the real dx and dy
            /*
            double ye = 0 ;
            double dyfit   = 0 ; 
            double fitstd = 0 ; 
            for (size_t j=0; j< col.size(); j++ ) {
                ye = (c1*xa[j]) + c0 ;
                dyfit = ya[j] - ye ;
                fitstd += (dyfit*dyfit);
                //fprintf(logfile," %.5f %.5f \n", ye, dyfit ) ;
            }
            double dof = (double)col.size() - 1 ;
            fitstd = sqrt( fitstd/dof  ) ;
            //fprintf(logfile, "\n %.5f \n", fitstd ) ;
            ls.push_back( fitstd );
            */
            double pull = Pull( xa, ya, csz ) ;
            ls.push_back( pull );
  
            lineSum.push_back( ls ) ;

            col.clear() ;
            col.push_back( data[i] ) ;
            val = data[i][axis1] ;
            flush = false ;
        }
    }
    if ( axis1 == 4 ) fprintf(logfile, " == Y Summary == \n" ) ;
    if ( axis1 == 3 ) fprintf(logfile, " == X Summary == \n" ) ;
    fprintf(logfile, "#,     Ave,     stdev,    Intercept,   Angle,  DoF,  res\n" ) ;
    for ( int i=0; i< (int)lineSum.size(); i++ ) {
        fprintf(logfile,"%d, %.5f, %.5f,  %f, %f, %.0f, %.5f\n", i, lineSum[i][0], lineSum[i][1], lineSum[i][2], atan2(lineSum[i][3],1), lineSum[i][4], lineSum[i][5] ) ;

    }

    fclose(logfile);
}

double Reader::Pull( double* xa, double* ya, int csz ) {

     const int sz = csz -1 ;
     double x[sz], y[sz] ; 
     double c0,c1,cov00, cov01, cov11, sumsq ;
     int k = 0 ;
     double ye(0), dy(0), stdy(0) ;
     for ( int i=0; i< csz ; i++) {
         k = 0 ;
         for ( int j= 0; j< sz; j++ ) {
             if ( j == i ) k = j+1  ;
             x[j] = xa[k] ;
             y[j] = ya[k] ;
             k++ ;
         }
         gsl_fit_linear( x, 1, y, 1, sz, &c0, &c1, &cov00, &cov01, &cov11, &sumsq ) ;
         ye = (c1*xa[i]) + c0 ;
         dy = (ya[i] - ye)*(ya[i] - ye) ;
         stdy += dy ;
     }
     double dof = (double)(csz -1) ;
     stdy = sqrt( stdy/ dof  ) ;
     return stdy ;
}

// Read data from Unpatterned ROI
void Reader::GetDataFromUPROI( ) {

    int site = 0 ;
    Input->GetParameters("Site",      &site ) ; 
     
    // Record the result
    string logFileName = hfolder +  outFileName ;
    FILE* logfile = fopen( logFileName.c_str() ,"w");

    if ( site ==0 ) { 
       for ( int i=1 ; i<5; i++ ) {
           GetDataFromUPROI( logfile, i ) ;
       }
    } else {
          GetDataFromUPROI( logfile, site )  ;
    }

    fclose(logfile) ;
}

void Reader::GetDataFromUPROI( FILE* logfile, int site_ ) {

    vector<vec> data ;
    data.clear() ;

    // Read data from the CSV file
    string csvfile = cfolder + cFileName ;
    printf(" CSV file = %s \n", csvfile.c_str() ) ;
    // output a double vector of 5 elements
    // ID, DieX DieY, CenterX, CenterY
    Input->ReadFreeCSV( csvfile, data, 1 ) ;  

    // Die 1 or Die 2     
    int die = 1 ;
    Input->GetParameters("Die",      &die ) ; 

    // Sorting sequence DieX first (axis = 5) or DieY first (axis = 6)
    int axis1 = 5 ;
    int axis2 = ( axis1 == 6 ) ? 5 : 6  ;

    // sorting data on Y
    if ( data.size() > 1 ) { 
       is = axis1 ;
       sort( data.begin(), data.end(), ScanSort );
       //printf(" small d = %.1f ~ %.1f \n", data[0], data[ data.size() - 1 ] ) ;
    }
   
    printf(" data size = %d \n", (int)data.size() ) ;

    // Giving Shot ID, lable site and die numbers
    vector<vec> col ;
    vector<vec> sumV ;
    vector<vec> s1V ;
    vec mdata ;
    double val = data[0][axis1] ;
    bool flush = false ;
    is = axis2 ;
    int cl(0), st(0), de(0) ;
    double dy = 2.2 ;
    double dx(0) ;
    for (size_t i=0; i < data.size(); i++ ) { 

        if ( fabs( val - data[i][axis1] ) < 0.1 ) {
           col.push_back( data[i] ) ;
           if ( i == data.size() -1 ) flush = true ;
        } else {
           flush = true ;
        }

        if ( flush ) { 
           sort( col.begin(), col.end(), ScanSort ) ;
           for ( int j=0; j< (int)col.size(); j++ ) {
               if ( j == 0 ) { 
                   dx = 0 ; 
                   dy = 2.2 ;
               } else {
                 dx = col[j][5] - col[j-1][5] ;
                 dy = col[j][6] - col[j-1][6] ;
               }
               if ( cl%2 == 0 ) de = 1 ;
               if ( cl%2 == 1 ) de = 2 ;

               if (  j%4  < 2 && cl%4  < 2 ) st = 3 ;
               if (  j%4 >= 2 && cl%4  < 2 ) st = 1 ;
               if (  j%4  < 2 && cl%4 >= 2 ) st = 4 ;
               if (  j%4 >= 2 && cl%4 >= 2 ) st = 2 ;

               //fprintf(logfile, "[%d,%d] %.0f %.4f %.4f  %.5f %.5f %.5f\n", st,de, (col[j][1]*10)+col[j][2], col[j][5], col[j][6], dx, dy, atan2(dx,dy)*1000 ) ;
               if ( de != die ) continue  ;
               // Only pick the top mark
               if ( dy > 8  ) {
                  mdata.clear() ;
		  mdata.push_back( st*1.0 );  
		  mdata.push_back( de*1.0 );  
		  mdata.push_back( (col[j][1]*10)+col[j][2] );  
		  mdata.push_back( col[j][5] );  
		  mdata.push_back( col[j][6] );  
		  mdata.push_back( dx );  
		  mdata.push_back( dy ); 
		  if ( st == site_ ) {
                     sumV.push_back( mdata ) ; 
                     //fprintf(logfile, "[%d,%d] %.0f %.4f %.4f  %.5f %.5f %.5f\n", st,de, (col[j][1]*10)+col[j][2], col[j][5], col[j][6], dx, dy, atan2(dx,dy)*1000 ) ;
                  }
		  if ( st == 1 )    s1V.push_back( mdata ) ; 
               }
           }
           flush = false ;
           val = data[i][axis1]; 
           col.clear() ;
           col.push_back(data[i] ) ;
           cl++ ;
        }
    }

    printf( " data1 size = %d \n", (int)sumV.size() ) ;
    printf( " data1 size = %d \n", (int)s1V.size() ) ;
    const int sz = 7 ;
    double xa[sz], ya[sz], dxa[sz], dya[sz], dA[sz] ;
    val = sumV[0][3] ; 
    int j = 0 ;

    // For Fitting
    double c0,c1,cov00, cov01, cov11, sumsq ;
    vector<vec> colsum ;

    fprintf(logfile, "Shot, S,D,   X,        Y,         dX,      dY,      dA\n");
    for ( size_t i=0; i< sumV.size() ; i++ ) {
        if ( fabs( sumV[i][3] - val ) < 0.5 ) { 
           xa[j]  = sumV[i][3] ; 
           ya[j]  = sumV[i][4] ;
           dA[j]  = atan2( sumV[i][5], sumV[i][6])*1000  ;
           if ( sumV[i][0] == 1 ) {
              dxa[j] = sumV[i][5] ; 
              dya[j] = sumV[i][6] ; 
           }
           if ( sumV[i][0] == 2 ) {
              dxa[j] = sumV[i][3] - s1V[i][3] - 14 ; 
              dya[j] = sumV[i][4] - s1V[i][4] ; 
           }
           if ( sumV[i][0] == 4 ) {
              dxa[j] = sumV[i][3] - s1V[i][3] - 14 ; 
              dya[j] = sumV[i][4] - s1V[i][4] + 14 ; 
           }
           if ( sumV[i][0] == 3 ) {
              dxa[j] = sumV[i][3] - s1V[i][3] ; 
              dya[j] = sumV[i][4] - s1V[i][4] + 14 ; 
           }
           fprintf(logfile, "%.0f,  %.0f,%.0f,  %.4f, %.4f,  %.5f, %.5f, %.5f\n", 
                     sumV[i][2], sumV[i][0], sumV[i][1], sumV[i][3], sumV[i][4], dxa[j], dya[j], dya[j] ) ;
           j++ ;
           if ( i == sumV.size()-1) flush = true ;
        } else {
            flush = true ;
        }
 
        val = sumV[i][3] ;
        if ( flush) {
            // Do some GSL functions
            gsl_fit_linear( ya, 1, xa, 1, sz, &c0, &c1, &cov00, &cov01, &cov11, &sumsq ) ;
            double Xmean  = gsl_stats_mean(  xa, 1, sz );
            double dXmean = gsl_stats_mean( dxa, 1, sz );
            double dYmean = gsl_stats_mean( dya, 1, sz );
            double dAmean  = gsl_stats_mean(  dA, 1, sz );
            double X_std  = gsl_stats_sd(  xa, 1, sz );
            double dX_std = gsl_stats_sd( dxa, 1, sz );
            double dY_std = gsl_stats_sd( dya, 1, sz );
            double dA_std = gsl_stats_sd(  dA, 1, sz );
            vec ls ;
            ls.push_back( Xmean );
            ls.push_back( X_std );
            ls.push_back( dXmean );
            ls.push_back( dX_std );
            ls.push_back( dYmean );
            ls.push_back( dY_std );
            ls.push_back( dAmean );
            ls.push_back( dA_std );
            ls.push_back( c0 ) ;
            ls.push_back( c1 ) ;
            colsum.push_back( ls );
            flush = false ;
            j= 0 ;
            //fprintf(logfile, "> %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f \n", xa[0], xa[1], xa[2], xa[3], xa[4], xa[5], xa[6]) ;
            //fprintf(logfile, " [%d]  %.5f, %.5f,  %.5f, %.5f,  %.5f, %.5f,  %.5f, %.5f,  %f, %f\n", 
            //   j, Xmean, X_std, dXmean, dX_std, dYmean, dY_std, dAmean, dA_std, c0, c1 );

	    xa[j]  = sumV[i][3] ; 
	    ya[j]  = sumV[i][4] ;
	    dA[j]  = atan2( sumV[i][5], sumV[i][6])*1000  ;
            if ( sumV[i][0] == 1 ) {
               dxa[j] = sumV[i][5] ; 
               dya[j] = sumV[i][6] ; 
            }
	    if ( sumV[i][0] == 2 ) {
		    dxa[j] = sumV[i][3] - s1V[i][3] -14. ; 
		    dya[j] = sumV[i][4] - s1V[i][4] ; 
	    }
	    if ( sumV[i][0] == 4 ) {
		    dxa[j] = sumV[i][3] - s1V[i][3] -14 ; 
		    dya[j] = sumV[i][4] - s1V[i][4] +14 ; 
	    }
	    if ( sumV[i][0] == 3 ) {
		    dxa[j] = sumV[i][3] - s1V[i][3] ; 
		    dya[j] = sumV[i][4] - s1V[i][4] +14 ; 
	    }
            fprintf(logfile, "\n" ) ;
            if ( i < sumV.size() -1 ) fprintf(logfile, "%.0f,  %.0f,%.0f,  %.4f, %.4f,  %.5f, %.5f, %.5f\n", 
                      sumV[i][2], sumV[i][0], sumV[i][1], sumV[i][3], sumV[i][4], dxa[j], dya[j], dya[j] ) ;
            j++ ;
        }
    }

    fprintf(logfile, "\n == Summary of Site %d == \n", site_ ) ;
    fprintf(logfile, "col,     X,      std_X,      dX,    std_dX,     dY,    std_dY,      dA,    std_dA,    Intercept,   Slope\n" ) ;
    for ( int i=0; i< (int)colsum.size(); i++ ) {
        fprintf(logfile, " %d,  %.5f, %.5f,  %.5f, %.5f,  %.5f, %.5f,  %.5f, %.5f,  %f, %f\n", 
        i, colsum[i][0], colsum[i][1], colsum[i][2], colsum[i][3], colsum[i][4], colsum[i][5], colsum[i][6], colsum[i][7], colsum[i][8], colsum[i][9] );
    }
    fprintf(logfile, "\n" ) ;

}


void Reader::GetDataFromTLSDPM( ) {
       
    int site = 0 ;
    Input->GetParameters("Site",      &site ) ; 
     
    // Record the result
    string logFileName = hfolder +  outFileName ;
    FILE* logfile = fopen( logFileName.c_str() ,"w");

    if ( site == 0 ) {
       for (int i=1 ; i< 5; i++ ) {
           GetDataFromTLSDPM( logfile, i ) ;
       }
    } else {
           GetDataFromTLSDPM( logfile, site ) ;
    }

    fclose(logfile) ;
}

void Reader::GetDataFromTLSDPM( FILE* logfile, int site ) {

    vector<vec> data ;
    data.clear() ;

    // Read data from the CSV file
    string csvfile = cfolder + cFileName ;
    printf(" CSV file = %s \n", csvfile.c_str() ) ;
    // output a double vector of 5 elements
    // ID, DieX, DieY, CenterX, CenterY, XError, YError, ThetaError
    Input->ReadFreeCSV( csvfile, data, 9 ) ;  

    // Sorting sequence X first (axis = 3) or Y first (axis = 4)
    int axis1 = 3 ;
    int axis2 = ( axis1 == 4 ) ? 3 : 4  ;

    // sorting data on Y
    if ( data.size() > 1 ) { 
       is = axis1 ;
       sort( data.begin(), data.end(), ScanSort );
       //printf(" small d = %.1f ~ %.1f \n", data[0], data[ data.size() - 1 ] ) ;
    }
    printf(" data size = %d \n", (int)data.size() ) ;


    // classify 4 different sites
    int dieX(0), dieY(0) ;
    vector<vec> slist ;
    vector<vec> s1list ;
    slist.clear() ;
    s1list.clear() ;
    for (size_t i=0; i < data.size(); i++ ) { 
        dieX = (int)(data[i][1]) ;
        dieY = (int)(data[i][2]) ;
        if ( dieX % 2 == 1 && dieY%2 == 0              ) s1list.push_back( data[i]) ;  
        if ( dieX % 2 == 1 && dieY%2 == 0 && site == 1 )  slist.push_back( data[i]) ;  
        if ( dieX % 2 == 0 && dieY%2 == 0 && site == 2 )  slist.push_back( data[i]) ;  
        if ( dieX % 2 == 1 && dieY%2 == 1 && site == 3 )  slist.push_back( data[i]) ;  
        if ( dieX % 2 == 0 && dieY%2 == 1 && site == 4 )  slist.push_back( data[i]) ;  
    }
    printf(" site list size = %d \n", (int)slist.size() ) ;
    

    fprintf(logfile, " Site %d \n", site ) ;
    fprintf(logfile, "col,DieX,DieY   X        Y         dX          dY         dA\n" ) ; 
    double val = slist[0][axis1] ;
    vector<vec> col ;
    vector<vec> col1 ;
    col.clear() ;
    col1.clear() ;
    is = axis2  ;
    vector<vec> lineSum ;
    lineSum.clear() ;
    bool flush = false ; 
    for (size_t i=0; i < slist.size(); i++ ) {
        if ( fabs(val - slist[i][axis1] ) < 0.1 ) {
            col.push_back( slist[i] ) ;
            col1.push_back( s1list[i] ) ;
            if ( i == slist.size() -1 ) flush = true ;
        } else {
            flush = true ;
        }

        if ( flush ) {
            // sorting data on X for `each Y values
            sort( col.begin(), col.end(), ScanSort ) ;
            sort( col1.begin(), col1.end(), ScanSort ) ;
            // For Fitting
            const int csz = (int) col.size() ;
            double xa[csz], ya[csz], ta[csz], tx[csz], ty[csz], ref[csz], refx[csz], refy[csz] ;
            double c0,c1,cov00, cov01, cov11, sumsq ;
            for ( size_t j=0; j < col.size(); j++ ) {

                xa[j]   = col[j][axis2] ;
                ya[j]   = col[j][axis1] ;

                ta[j]   = col[j][7]*3141.5926/180.0 ;
                ty[j]   = col[j][6] ;
                tx[j]   = col[j][5] ;
                ref[j]  = col1[j][7]*3141.5926/180.0 ;
                refy[j] = col1[j][6] ;
                refx[j] = col1[j][5] ;

                fprintf(logfile, " %d, %02.0f,%02.0f, %.4f, %.4f,  %.6f,  %.6f,  %.6f\n", 
                       (int)j, col[j][1], col[j][2], col[j][3], col[j][4], tx[j], ty[j], ta[j] ) ;
                //fprintf(logfile, "1 %d, %02.0f,%02.0f, %.4f, %.4f,  %.6f,  %.6f,  %.6f\n", 
                //       (int)j, col1[j][1], col1[j][2], col1[j][3], col1[j][4], refx[j], refy[j], ref[j] ) ;
                tx[j] = tx[j] - refx[j] ;
                ty[j] = ty[j] - refy[j] ;
                ta[j] = ta[j] - ref[j] ;
            }
            // Do some GSL functions
            gsl_fit_linear( xa, 1, ya, 1, csz, &c0, &c1, &cov00, &cov01, &cov11, &sumsq ) ;
            double meanX  = gsl_stats_mean( tx, 1, csz );
            double meanY  = gsl_stats_mean( ty, 1, csz );
            double meanA  = gsl_stats_mean( ta, 1, csz );
            double stdX   = gsl_stats_sd(tx, 1, csz );
            double stdY   = gsl_stats_sd(ty, 1, csz );
            double stdA   = gsl_stats_sd(ta, 1, csz );
            vec ls ;
            ls.push_back( meanX ) ;
            ls.push_back( stdX ) ;
            ls.push_back( meanY ) ;
            ls.push_back( stdY ) ;
            ls.push_back( meanA ) ;
            ls.push_back( stdA ) ;
            ls.push_back( c0 ) ;
            ls.push_back( c1 ) ;
            ls.push_back( (double)csz ) ;
            lineSum.push_back( ls ) ;
            fprintf(logfile, " End of Column %d \n", (int)lineSum.size() ) ;

            col.clear() ;
            col.push_back( slist[i] ) ;
            col1.clear() ;
            col1.push_back( s1list[i] ) ;
            val = slist[i][axis1] ;
            flush = false ;
        }
    }
    fprintf(logfile, "\n == Summary of Site %d == \n", site ) ;
    fprintf(logfile, "col,    X,       dX,        Y,      dY,         A,       dA,      Intercept,    Slope,  (N_entries)\n" ) ;
    for ( int i=0; i< (int)lineSum.size(); i++ ) {
        fprintf(logfile, " %d,  %.5f, %.5f,  %.5f, %.5f,  %.5f, %.5f,   %f,  %f,  %.0f\n", 
             i, lineSum[i][0], lineSum[i][1], lineSum[i][2], lineSum[i][3], lineSum[i][4], lineSum[i][5], lineSum[i][6], lineSum[i][7], lineSum[i][8]);

    }
    fprintf(logfile, "\n" ) ;

}

// PSET File
void Reader::GetPsetFromPTI() {

    vector<vec> data ;
    data.clear() ;
    vector<vec> data1 ;
    data1.clear() ;

    // Read data from the CSV file
    string csvfile = cfolder + mapFileName ;
    printf(" CSV file = %s \n", csvfile.c_str() ) ;
    // output a double vector of 5 elements
    // ID, XSize YSize, CenterX, CenterY
    Input->ReadFreeCSV( csvfile, data, 3 ) ;  

    int axis = 1 ;
    Input->GetParameters("Axis",         &axis );

    // Sorting sequence X first (axis = 1) or Y first (axis = 2)
    // Only Axis1 need to be changed, Datacard output filename also need to be changed.
    int axis1 = axis ;
    int axis2 = ( axis1 == 2 ) ? 1 : 2 ;

    // Record the result
    string logFileName = hfolder +  outFileName ;
    FILE* logfile = fopen( logFileName.c_str() ,"w");

    // sorting data on Axis1
    if ( data.size() > 1 ) { 
       is = axis1 ;
       sort( data.begin(), data.end(), ScanSort );
       //printf(" small d = %.1f ~ %.1f \n", data[0], data[ data.size() - 1 ] ) ;
    }
    printf(" data size = %d \n", (int)data.size() ) ;

    // Lable Die
    double x0 = -242.125 ;
    double y0 = -239.575 ;
    double y_pitch = 6.9 ;
    double x_ = x0 ;
    double dieX(0), dieY(0) ;
   
    bool nextCol = false ;
    vector<vec> col ;
    for ( size_t i =0; i < data.size(); i++ ) {

        double dx_ = fabs(  data[i][axis1] - x_ ) ;
        
        if ( dx_ < 0.2 ) {
           col.push_back( data[i] ) ;
           if ( i == data.size() -1 ) nextCol = true ;
        } else {
          nextCol = true ;
        }

        if ( nextCol ) {
           is = axis2 ;
           sort( col.begin(), col.end(), ScanSort );
           
	   for( size_t j=0; j< col.size(); j++ ) {
              if ( col[j][axis2] > 0 ) y0 = 6.325 ;
              else                     y0 = -239.575 ;

              double frac = modf( fabs( col[j][axis2] - y0 )/y_pitch , &dieY ) ;
              if ( frac > 0.9 ) dieY  = dieY + 1 ;
              if ( col[j][axis2] > 0 ) dieY = dieY + 34 ;
                
              col[j].push_back(dieX);
              col[j].push_back(dieY);
              col[j][0] = (dieX*100)+dieY ;
              //fprintf(logfile, " %.0f, %.4f, %.4f, %.4f, %.4f, %.0f, %.0f, %.0f \n", 
              //                   col[j][0], col[j][1], col[j][2], col[j][3], col[j][4], col[j][6], col[j][7], col[j][8] ) ; 
              data1.push_back( col[j] ) ;
	   }

           if ( dx_ < 2. || dx_ > 8. ) dieX += 1 ; 
           col.clear() ;
           nextCol = false ;
        } 
        x_ = data[i][axis1] ;      
  
    } 
    //fprintf(logfile, "===================\n") ;
    if ( data.size() > 1 ) { 
       MathTools *math = new MathTools() ;
       is = 0 ;
       sort( data1.begin(), data1.end(), ScanSort );
       for (size_t j=0; j< data1.size(); j++) {

            double theta = 0 ;
            if ( j > 0 && data1[j][0] == data1[j-1][0] ) {
               double v1x = data1[j-1][1] + data1[j-1][3] ;
	       double v1y = data1[j-1][2] + data1[j-1][4] ;
	       double v2x = data1[j][1] + data1[j][3] ;
	       double v2y = data1[j][2] + data1[j][4] ;
	       theta = math->angleAB( v1x, v1y, v2x, v2y ) ;
            }
          
            fprintf(logfile, " %.0f, %.4f, %.4f, %.4f, %.4f, %.0f, %.0f, %.5f \n", 
                                 data1[j][0], data1[j][1], data1[j][2], data1[j][3], data1[j][4], data1[j][6], data1[j][7], theta ) ; 
       }
       //printf(" small d = %.1f ~ %.1f \n", data[0], data[ data.size() - 1 ] ) ;
    }

    fclose( logfile ) ;
}

void Reader::SkimXCD( FILE* logfile, string csvfile ) {

    // vector to store data
    vector<vec> data ;
    data.clear() ;

    // output a double vector of 5 elements
    Input->ReadFreeCSV( csvfile, data, 1 ) ;  

    fprintf(logfile, "ID,dieX,dieY,XSize,YSize,X,Y \n") ;
    for (size_t i=0; i < data.size(); i++ ) {

        //printf(" -> [%.0f] size = %d \n", data[i][8], (int)data[i].size() ) ;
        if ( data[i][8] != 231 ) continue ; 
        fprintf( logfile, "%.0f,%.0f,%.0f,%f,%f,%f,%f,%.0f\n", 
                       data[i][0], data[i][1], data[i][2], data[i][12],data[i][13],data[i][14],data[i][15],data[i][17] ) ;
 
    }
    //fclose( logfile ) ;

}

void Reader::GetDataXCD() {

  int nLoad = 1 ;
  Input->GetParameters("N_Load",        &nLoad );
  int nRun = 1 ;
  Input->GetParameters("N_Run",         &nRun );

  if ( nLoad ==1 && nRun == 1 ) {
     // Record the result
     string logFileName = hfolder +  outFileName ;
     FILE* logfile = fopen( logFileName.c_str() ,"w");

     // Read data from the CSV file
     string csvfile = cfolder + cFileName ;
     printf(" CSV file = %s \n", csvfile.c_str() ) ;

     SkimXCD( logfile, csvfile ) ;
     fclose( logfile ) ; 
  } else {

     for ( int i=1; i< nLoad+1; i++ ) {
         for ( int j=1; j< nRun+1; j++ ) {

             char cfname[20] ;
	     sprintf(cfname, "TLS_XCD_LUL%d_%d.csv", i, j ) ;
	     string csvFileName = cfolder +  cfname ;

	     char logname[20] ;
	     sprintf(logname, "TLS_XCD_L%d_%d.csv", i, j ) ;
	     string logFileName = hfolder +  logname ;
             FILE* logfile = fopen( logFileName.c_str() ,"w");

             printf(" CSV file = %s  => %s \n", csvFileName.c_str(), logFileName.c_str() ) ;
             SkimXCD( logfile, csvFileName ) ;
             fclose( logfile ) ;
         }
     }

  } 
}

void Reader::GetDieStatistic( ) {

    vector<vec> data ;
    data.clear() ;

    // Read data from the CSV file
    string csvfile = cfolder + cFileName ;
    printf(" CSV file = %s \n", csvfile.c_str() ) ;
    // output a double vector of 5 elements
    // ID, DieX, DieY, Size, SizeX, SizeY, CenterX, CenterY, 
    Input->ReadFreeCSV( csvfile, data, 1 ) ;  

    // Sorting sequence X first (axis = 1) or Y first (axis = 2)
    int axis1 = 1 ;
    int axis2 = ( axis1 == 2 ) ? 1 : 2  ;

    // sorting data on Y
    if ( data.size() > 1 ) { 
       is = axis1 ;
       sort( data.begin(), data.end(), ScanSort );
       //printf(" small d = %.1f ~ %.1f \n", data[0], data[ data.size() - 1 ] ) ;
    }
    printf(" data size = %d \n", (int)data.size() ) ;

    //Open Root file to store result
    string hfile = hfolder + hFileName + ".root" ;
    TFile* rootFile = new TFile( hfile.c_str() , "RECREATE" );
    rootFile->cd() ;

    // Define histograms
    TH1D* hDcnt = new TH1D("hDcnt" , "CF defect count per die",   30, 0, 30 )  ;

    int ndef = 0 ;
    double dieX(0), dieY(0) ;
    for ( size_t i=0; i< data.size() ; i++ ) {

        
        if (  dieX == data[i][1] && dieY == data[i][2] ) {
           ndef++ ;
        } else {
           hDcnt->Fill( ndef ) ;
           //printf(" DCount = %d \n", ndef ) ;
           ndef = 1 ;
           dieX = data[i][1] ;
           dieY = data[i][2] ;
        }
        //printf(" [%.0f][%.0f] \n", data[i][1], data[i][2] ) ;
    }
    printf(" ======== Defect Count per die =========== \n") ;
    for ( int i=0; i<32; i++) {
        double bc = hDcnt->GetBinContent(i) ; 
        //printf("[%d] = [%.0f] \n", i, bc ) ;
        printf(" %d, %.0f \n", i, bc ) ;
    }


    hDcnt->Write() ;
    rootFile->Close() ;
 
}
