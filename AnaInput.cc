#include "AnaInput.h"

AnaInput* AnaInput::m_Instance = NULL ;

//AnaInput::AnaInput( string datacardInput ) {
AnaInput::AnaInput() {

  datacardfile = "DataCard.txt" ;
  xLimit_l = -1 ; 
  xLimit_r = 999999 ; 
  yLimit_d = -1 ; 
  yLimit_u = 999999 ; 

}

AnaInput::~AnaInput(){

   cout<<" close input "<<endl ;

}

AnaInput* AnaInput::Instance() {

    if (! m_Instance ) {
       m_Instance = new AnaInput( ) ;
    } 
       
    return m_Instance ;

}

int iis;
bool SortI( vec s1, vec s2) { return (s1[iis] < s2[iis]) ; }

void AnaInput::SetDatacard( string datacardInput ) {
 
     datacardfile = datacardInput ;
}


// Methods to read DataCard.txt
void AnaInput::GetParameters(string paraName, int* thePara, string cfgFile ){

     if ( cfgFile == "" ) cfgFile = datacardfile ;

     fstream paraFile( cfgFile.c_str() );
     if ( !paraFile.is_open() )  cout<<" file opened error => check file path and the folder "<<endl;
     string  line;
     string  getName;
     string  getValue;
     size_t  pos ;
     size_t  vpos ;

     bool gotIt = false ;
     while ( getline(paraFile, line) ) {
           if ( line[0] == '#' ) continue ;

           pos = line.find( paraName );
           vpos = pos + paraName.size() + 2;
           if ( pos < line.npos ) {
              string str_end = line.substr(vpos-1, 1) ;
              if ( str_end == ' ' || str_end == '=') {
                 getName  = line.substr( pos, paraName.size() );
                 getValue = line.substr( vpos );
                 *thePara = atoi( getValue.c_str() );
                 //cout<< paraName <<" = "<< *thePara << endl;
                 gotIt = true;
              }
           }
           if ( gotIt ) break ;
     }
     paraFile.close();
}

void AnaInput::GetParameters(string paraName, double* thePara, string cfgFile ){

     if ( cfgFile == "" ) cfgFile = datacardfile ;

     fstream paraFile( cfgFile.c_str() );
     if ( !paraFile.is_open() )  cout<<" file opened error => check file path and the folder "<<endl;
     string  line;
     string  getName;
     string  getValue;
     size_t  pos ;
     size_t  vpos ;

     bool gotIt = false ;
     while ( getline(paraFile, line) ) {
           if ( line[0] == '#' ) continue ;

           pos = line.find( paraName );
           vpos = pos + paraName.size() + 2;
           // exclude the case when this paraName is the sub-string of other paraName 
           if ( pos > 0 && pos < 999 && line[pos-1] != ' ' ) continue ;

           if ( pos < line.npos ) {
              string str_end = line.substr(vpos-1, 1) ;
              if ( str_end == ' ' || str_end == '=') {
                 getName  = line.substr( pos, paraName.size() );
                 getValue = line.substr( vpos );
                 *thePara = atof( getValue.c_str() );
                 //cout<< paraName <<" = "<< *thePara << endl;
                 gotIt = true ;
              }
           }
           if ( gotIt ) break ;
     }
     paraFile.close();
}

void AnaInput::GetParameters(string paraName, string* thePara, string cfgFile ){

     if ( cfgFile == "" ) cfgFile = datacardfile ;

     fstream paraFile( cfgFile.c_str() );
     if ( !paraFile.is_open() )  cout<<" file opened error => check file path and the folder "<<endl;
     string  line;
     string  getName;
     size_t  pos ;
     size_t  vpos ;

     bool gotIt = false ;
     while ( getline(paraFile, line) ) {
           if ( line[0] == '#' ) continue ;

           pos = line.find( paraName );
           vpos = pos + paraName.size() + 2;
           // exclude the case when this paraName is the sub-string of other paraName 
           if ( pos > 0 && pos < 999 && line[pos-1] != ' ' ) continue ;

           if ( pos < line.npos ) {
              string str_end = line.substr(vpos-1, 1) ;
              if ( str_end == ' ' || str_end == '=') {
                 //cout<<" pos = "<< pos <<endl;
                 getName  = line.substr( pos, paraName.size() );
                 //*thePara = line.substr( vpos );
                 //cout<< paraName <<" = "<< *thePara << endl;
                 string strTmp = line.substr( vpos );
                 for (string::iterator it = strTmp.begin(); it< strTmp.end(); it++) {
                     if ( (*it) != ',' && (*it) != ' ' && (*it) != '(' && (*it) != ')' && (*it) != '=') thePara->push_back( *it );
                 }
                 gotIt = true ;
              }
           }
           if ( gotIt ) break;
     }
     paraFile.close();
}

void AnaInput::GetParameters(string paraName, vector<double>* thePara, string cfgFile ){

     if ( cfgFile == "" ) cfgFile = datacardfile ;

     fstream paraFile( cfgFile.c_str() );
     if ( !paraFile.is_open() )  cout<<" file opened error => check file path and the folder "<<endl;
     string  line;
     string  getName;
     string  getValue;
     size_t  pos(0) ;
     size_t  vpos(0) ;
     vector<double>  vvec;

     while ( getline(paraFile, line) ) {
           if ( line[0] == '#' ) continue ;

           pos = line.find( paraName );
           vpos = pos + paraName.size() + 1;
           // exclude the case when this paraName is the sub-string of other paraName 
           if ( pos > 0 && pos < 999 && line[pos-1] != ' ' ) continue ;

           if ( pos < line.npos ) {
              getName  = line.substr( pos, paraName.size() );
              string arrVal = line.substr( vpos );
              if ( arrVal[0] != '=' && arrVal[0] != ' ' ) continue;
	      int vidx = 0;
	      string vtemp ;
	      //cout<< paraName <<" = ( " ;
              for (string::iterator it = arrVal.begin(); it< arrVal.end(); it++) {
                  if ( (*it) != ',' && (*it) != ' ' && (*it) != '(' && (*it) != ')' && (*it) != '=') vtemp.push_back( *it );
                  if ( (*it) == ',' || (*it) == ')' ) { 
                     if ( vtemp.size() > 0 ) vvec.push_back( atof( vtemp.c_str() ) ) ;
		     //cout<< vtemp << *it;
		     vidx++ ;
		     vtemp.clear() ;
                  }
              }
              *thePara = vvec ;
           }
     }
     //cout<<""<<endl ;
     paraFile.close();

} 

void AnaInput::GetParameters(string paraName, vector<float>* thePara, string cfgFile ){

     if ( cfgFile == "" ) cfgFile = datacardfile ;

     fstream paraFile( cfgFile.c_str() );
     if ( !paraFile.is_open() )  cout<<" file opened error => check file path and the folder "<<endl;
     string  line;
     string  getName;
     string  getValue;
     size_t  pos(0) ;
     size_t  vpos(0) ;
     vector<float>  vvec;

     while ( getline(paraFile, line) ) {
           if ( line[0] == '#' ) continue ;

           pos = line.find( paraName );
           vpos = pos + paraName.size() + 1;
           // exclude the case when this paraName is the sub-string of other paraName 
           if ( pos > 0 && pos < 999 && line[pos-1] != ' ' ) continue ;

           if ( pos < line.npos ) {
              getName  = line.substr( pos, paraName.size() );
              string arrVal = line.substr( vpos );
              if ( arrVal[0] != '=' && arrVal[0] != ' ' ) continue;
	      int vidx = 0;
	      string vtemp ;
	      //cout<< paraName <<" = ( " ;
              for (string::iterator it = arrVal.begin(); it< arrVal.end(); it++) {
                  if ( (*it) != ',' && (*it) != ' ' && (*it) != '(' && (*it) != ')' && (*it) != '=') vtemp.push_back( *it );
                  if ( (*it) == ',' || (*it) == ')' ) { 
                     if ( vtemp.size() > 0 ) vvec.push_back( atof( vtemp.c_str() ) ) ;
		     //cout<< vtemp << *it;
		     vidx++ ;
		     vtemp.clear() ;
                  }
              }
              *thePara = vvec ;
           }
     }
     //cout<<""<<endl ;
     paraFile.close();

} 

void AnaInput::GetParameters(string paraName, vector<string>* thePara, string cfgFile ){

     if ( cfgFile == "" ) cfgFile = datacardfile ;

     fstream paraFile( cfgFile.c_str() );

     if ( !paraFile.is_open() )  cout<<" file opened error => check file path and the folder "<<endl;
     string  line;
     string  getName;
     string  getValue;
     size_t  pos ;
     size_t  vpos ;
     vector<string>  vvec;

     while ( getline(paraFile, line) ) {
           if ( line[0] == '#' ) continue ;

           pos = line.find( paraName );
           vpos = pos + paraName.size() ;
           if ( pos < line.npos ) {
              getName  = line.substr( pos, paraName.size() );
              string arrVal = line.substr( vpos );
              if ( arrVal[0] != '=' && arrVal[0] != ' ' ) continue;
	      int vidx = 0;
	      string vtemp ;
	      //cout<< paraName <<" = ( " ;
              for (string::iterator it = arrVal.begin(); it< arrVal.end(); it++) {
                  if ( (*it) != ',' && (*it) != ' ' && (*it) != '(' && (*it) != ')' && (*it) != '=') vtemp.push_back( *it );
                  if ( (*it) == ',' || (*it) == ')' ) { 
                     if ( vtemp.size() > 0 ) vvec.push_back( vtemp ) ;
		     //cout<< vtemp << *it;
		     vidx++ ;
		     vtemp.clear() ;
                  }
              }
              *thePara = vvec ;
           }
     }
     paraFile.close();

}
 
void AnaInput::GetParameters(string paraName, vector<int>* thePara, string cfgFile ){

     if ( cfgFile == "" ) cfgFile = datacardfile ;

     fstream paraFile( cfgFile.c_str() );
     if ( !paraFile.is_open() )  cout<<" file opened error => check file path and the folder "<<endl;
     string  line;
     string  getName;
     string  getValue;
     size_t  pos ;
     size_t  vpos ;
     vector<int>  vvec;

     while ( getline(paraFile, line) ) {
           if ( line[0] == '#' ) continue ;

           pos = line.find( paraName );
           vpos = pos + paraName.size() ;
           if ( pos < line.npos ) {
              getName  = line.substr( pos, paraName.size() );
              string arrVal = line.substr( vpos );
              if ( arrVal[0] != '=' && arrVal[0] != ' ' ) continue;
	      //int vidx = 0;
	      string vtemp ;
	      //cout<< paraName <<" = ( " ;
              for (string::iterator it = arrVal.begin(); it< arrVal.end(); it++) {
                  if ( (*it) != ',' && (*it) != ' ' && (*it) != '(' && (*it) != ')' && (*it) != '=') vtemp.push_back( *it );
                  if ( (*it) == ',' || (*it) == ')' ) { 
                     if ( vtemp.size() > 0 ) vvec.push_back( atoi( vtemp.c_str() ) ) ;
		     //cout<< vtemp << *it;
		     //vidx++ ;
		     vtemp.clear() ;
                  }
              }
              *thePara = vvec ;
           }
     }
     paraFile.close();

}

// Read KLARF files - Only record defect Size
int AnaInput::ReadKLARF( string fileName, vector<df>& data, double lowLimit ) {

  // Open data file
  //printf(" File:%s \n", fileName.c_str() ) ;
  fstream logfile( fileName.c_str() );

  // Open KLARF file to read
  if ( !logfile.is_open() )  { 
      printf(" file open error \n") ;
      return 1 ;
  }
 
  // Variables to hold data 
  int id, xi,yi,xs,ys,da, clN  ;
  float x,y,dsize ;
  int fo = -1 ;
  int debug = 0 ;
  GetParameters( "debug", &debug ) ;
  // data collection
  bool get = false ;
  df entry ;

  string  line;
  while ( getline( logfile, line) ) {
        string dflist  = line.substr( 0, 10 );
        //cout<< dflist.c_str() << endl ;
        if ( strcmp( dflist.c_str() , "DefectList" ) == 0 ) { 
           get = true ;
           continue ;
        } 
        if ( !get ) continue ;
 
        // ID X Y XINDEX YINDEX XSIZE YSIZE DEFECTAREA DSIZE CLASSNUMBER TEST 
        fo = sscanf( line.c_str() , "%d %f %f %d %d %d %d %d %f %d" , &id, &x, &y, &xi, &yi, &xs, &ys, &da, &dsize, &clN );
        if (debug) printf(  " >> %d %.0f %.0f %d %d %d %d %d (%.1f) %d \n" , id, x, y, xi, yi, xs, ys, da, dsize, clN );
 
        if ( da < 1 ) {
            entry.x = x ;
            entry.y = y ;
            entry.s = (dsize > 4096 ) ? 4096 : dsize ;
            //if ( Exclusion( entry ) )  data.push_back( entry ) ;
            if ( entry.s > lowLimit ) data.push_back( entry ) ;
        }

        if ( fo != 10 ) {
          // printf(" !!! scan done (%d) \n", fo ) ;
           break;
        }
  }
  logfile.close() ;

  //printf(" data size == %d , %d measurement \n", (int)data.size(), id+1 ) ;
  if ( data.size() < 1 )  {
     return 1 ;
  } else {
     return 0 ;
  }

}

// Read the map file from Stepper 
int AnaInput::ReadMap( string fileName, vector<vec>& data, int nSkipLine ) {

  // Open data file
  //printf(" File:%s \n", fileName.c_str() ) ;
  fstream logfile( fileName.c_str() );

  // Open CSV file to read
  if ( !logfile.is_open() )  { 
      printf(" file open error \n") ;
      return 1 ;
  }
 
  // Variables to hold data 
  double id, x,y,dx,dy  ;
  char mk ;
  int fo = -1 ;
  int debug = 0 ;
  GetParameters( "debug", &debug ) ;
  // data collection
  vec entry ;

  string  line;
  int ii = 0 ;
  while ( getline( logfile, line) ) {
        if ( ii < nSkipLine ) continue ;
        ii++ ;
 
        if ( line.size() < 72 ) continue ;
        // read data : ID X Y dX dY type
        entry.clear() ;
        fo = sscanf( line.c_str() , " %lf, %lf, %lf, %lf, %lf, %c, " , &id, &x, &y, &dx, &dy, &mk );
        if ( mk != 'M' ) continue ;
        if (debug) printf(  " (%d) %.0f %.3f %.3f %.5f %.5f %c \n" ,fo,  id, x, y, dx, dy, mk );
        entry.push_back(id) ;
        entry.push_back(x) ;
        entry.push_back(y) ;
        entry.push_back(dx) ;
        entry.push_back(dy) ;

        if ( fo != 6 || entry.size() != 5 ) {
           printf(" !!! scan fail (%d) \n", fo ) ;
           break;
        } 
        data.push_back( entry ) ;
  }
  logfile.close() ;

  printf(" data size == %d \n", (int)data.size() ) ;
  if ( data.size() < 1 )  {
     return 1 ;
  } else {
     return 0 ;
  }

}

// Read site correction map 
void AnaInput::ReadParameter(string paraName, double* thePara, string cfgFile ){

     fstream paraFile( cfgFile.c_str() );
     if ( !paraFile.is_open() )  cout<<" file opened error => check file path and the folder "<<endl;
     string  line;
     string  getName;
     string  getValue;
     size_t  pos ;
     size_t  vpos ;

     while ( getline(paraFile, line) ) {
           if ( line[0] != '%' ) continue ;

           pos = line.find( paraName );
           vpos = pos + paraName.size() + 1 ;
//           cout<< paraName <<" pos = "<< pos <<" vpos "<< vpos <<  endl;
           
           // exclude the case when this paraName is the sub-string of other paraName 
           //if ( pos > 0 && pos < 999 && line[pos-1] != ' ' ) continue ;

           if ( pos < line.npos ) {
              string str_end = line.substr(vpos-1, 1) ;
              if ( str_end == ' ' || str_end == ':' || str_end == '	' ) {
                 getName  = line.substr( pos, paraName.size() );
                 getValue = line.substr( vpos );
                 *thePara = atof( getValue.c_str() );
                 cout<< paraName <<" = "<< *thePara << endl;
              }
           }
     }
     paraFile.close();
}

// Read Data from Map Correction Files
void AnaInput::ReadMapCorrection( string cfgFile, vector<vec>& data ){

     fstream paraFile( cfgFile.c_str() );
     if ( !paraFile.is_open() )  cout<<" file opened error => check file path and the folder "<<endl;
     string  line;
     string  getName;
     string  getValue;
     size_t  pos ;
     size_t  vpos ;

     string paraList[6] = { "SS        Site", "Site Xcoord(horizontal)", "Site Ycoord(vertical)", 
                                             "Site offsetX(horizontal)", "Site offsetY(vertical)", "Site offsetT(theta)" } ;

     double val = 0.0 ;
     vec entry ;
     entry.clear() ;

     string paraName ;
     while ( getline(paraFile, line) ) {
           if ( line[0] != '%' ) continue ;

           for ( int i=0; i< 6; i++) { 
               paraName = paraList[i] ;
    
               pos = line.find( paraName );
               vpos = pos + paraName.size() + 1;
               //cout<< paraName <<" pos = "<< pos <<" vpos "<< vpos <<  endl;

               if ( pos < line.npos ) {
                  string str_end = line.substr(vpos-1, 1) ;
                  if ( str_end == ' ' || str_end == ':' ||  str_end == '	' ) {
                     getName  = line.substr( pos, paraName.size() );
		     getValue = line.substr( vpos );
		     val = atof( getValue.c_str() );
		     //cout<< paraName <<" = "<< val << endl;
                     entry.push_back( val ) ;
                  }
               }
               if ( entry.size() == 6 && i == 5  ) {
                  //printf(" [%d]-> (%d)!!!\n", i , (int)entry.size() ) ;
                  if ( entry[3] == 0.0 && entry[4] == 0.0 && entry[5] == 0. ) {
                     entry.clear() ; 
                  } else {
                     data.push_back( entry ) ;
                     entry.clear() ;
                  }
               }
           }
     }
     paraFile.close();
}


int AnaInput::ReadDPM( string fileName, vector<vec>& data ) {

  // Open data file
  //printf(" File:%s \n", fileName.c_str() ) ;
  fstream logfile( fileName.c_str() );

  // Open CSV file to read
  if ( !logfile.is_open() )  { 
      printf(" file open error \n") ;
      return 1 ;
  }
 
  // Variables to hold data 
  double id, dieX, dieY, diedX, diedY, diedR, diedA, x,y, sth  ;
  int fo(-1) ;
  int debug = 0 ;
  GetParameters( "debug", &debug ) ;
  // data collection
  vec entry ;

  bool startRead = false ;
  string  line ;
  while ( getline( logfile, line) ) {
        if ( !startRead ) {
            if ( line[0] == 'I' && line[1] == 'D' ) {
              startRead = true ;
            }
            continue ;
        } 

        // read data : ID X Y dX dY type
        entry.clear() ;

        // sscanf need fixed length for the string reading so the whole line need to be broken into two 
        // in order to read the 2nd part of information
        fo = sscanf( line.c_str() , "%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, ,%lf", 
                                    &id, &dieX, &dieY, &diedX, &diedY, &diedR, &diedA, &x, &y, &sth );
        if (debug) printf(  " (%d) %.0f %.0f %.0f %.4f %.4f %.4f %.4f %.4f %.4f %.6f \n" ,
                               fo,  id, dieX, dieY, diedX, diedY, diedR, diedA, x, y, sth );

        entry.push_back(id) ;
        entry.push_back(dieX) ;
        entry.push_back(dieY) ;
        entry.push_back(diedX) ;
        entry.push_back(diedY) ;
        entry.push_back(diedR) ;
        entry.push_back(diedA) ;
        entry.push_back(x/1000.) ;
        entry.push_back(y/1000.) ;
        entry.push_back(sth) ;
        
        if ( fo != 10 || entry.size() != 10 ) {
           printf(" !!! scan fail (%d) \n", fo ) ;
           break;
        } 
        
        data.push_back( entry ) ;
  }
  logfile.close() ;

  printf(" data size == %d \n", (int)data.size() ) ;
  if ( data.size() < 1 )  {
     return 1 ;
  } else {
     return 0 ;
  }

}


// unit for x,y coordinate -> um
//bool AnaInput::SetExclusion( df& defect, float xl, float xr, float yu, float yd, float threshold ) {
bool AnaInput::SetExclusion( df& defect, vector<float>* limit, float threshold ) {
      
     bool InArea = true ;
     if ( (*limit).size() == 4 ) {

        float xl = (*limit)[0]  ;
	float xr = (*limit)[1]  ;
	float yu = (*limit)[2]  ;
	float yd = (*limit)[3]  ;

	bool inX =  ( defect.x > xl && defect.x < xr ) ;
	bool inY =  ( defect.y > yd && defect.y < yu ) ;
        if ( !inX || !inY ) InArea = false ;

     } else if ( (*limit).size() == 3 ) {

        float x_center = (*limit)[0] ;
        float y_center = (*limit)[1] ;
        float radius   = (*limit)[2] ;

        float dx = defect.x - x_center ;
	float dy = defect.y - y_center ;
	float dr = sqrt( (dx*dx) + (dy*dy) ) ;
	InArea =  ( dr < radius ) ? true : false ;
        
        //printf("limit x,y : (%.f , %.f) , r =  %.f - %.f \n", dx, dy , dr, radius ) ;
     } 


     bool inSize = ( defect.s > threshold_k ) ? true : false ;

     //printf("limit x : %.f - %.f , y: %.f - %.f \n", xLimit_l, xLimit_r , yLimit_d, yLimit_u  ) ;
     //printf("defect x : %.f , %.f \n", defect.x, defect.y  ) ;

     bool pass = ( InArea && inSize ) ? true  : false ; 
     //if ( pass ) cout<<"pass " <<endl ;
     return pass ;
}

// Read the output SINF file
void AnaInput::ReadSINF( string fileName, vector<vec>& data ) {

  // Open data file
  //printf(" File:%s \n", fileName.c_str() ) ;
  fstream logfile( fileName.c_str() );

  // Open CSV file to read
  if ( !logfile.is_open() )  { 
      printf(" file open error \n") ;
      //return 1 ;
  }
 
  // Variables to hold data 
  int nRow(0), nCol(0)  ;
  int debug = 0 ;

  GetParameters( "debug", &debug ) ;
  // data collection
  vec rowV ;
  rowV.clear() ;

  bool startRead = false ;
  string  line, line1, str_ ;
  int ii =0 ;
  while ( getline( logfile, line) ) {

        ii++ ;
        str_ = line.substr(0,6) ; 

        if ( str_  == "ROWCT:" ) {
           string row_str = line.substr(6) ;
	   //nRow =  stoi( row_str ) ;
	   nRow =  atoi( row_str.c_str()  ) ;
           printf(" row = %d \n", nRow ) ;
           continue ;
        }
	if ( str_  == "COLCT:" ) {
           string col_str = line.substr(6) ;
	   //nCol =  stoi( col_str ) ;
	   nCol =  atoi( col_str.c_str() ) ;
           printf(" col = %d \n", nCol ) ;
           continue ;
	}
        if ( str_  != "RowDat" ) continue ;

        if ( str_  == "RowDat" ) {
           startRead = true ;
           line1 = line.substr(8) ;   
        }

        // read data - Die info
        if ( startRead ) {
	   int vidx = 0 ;
	   double dfID = 0 ;
	   string vtemp ;
	   for (string::iterator it = line1.begin(); it< line1.end(); it++) {
               if ( (*it) != ',' && (*it) != ' ' && (*it) != ':') {
                  vtemp.push_back( *it );
               }
               if ( it == line1.end() -1 ) {
                  vtemp.push_back( *it );
                  //printf(" (%s) \n", vtemp.c_str() ) ;
               }
               if ( (*it) == ' ' || it == line1.end()-1 ) {
                  //printf("%s ", vtemp.c_str() ) ;
                  if ( vtemp[0] == '_' ) {
                     dfID = -9 ;
                  } else if ( vtemp[0] == '@' ) { 
                     dfID = -1 ;
                  } else {
                     int hx_ = strtol( vtemp.c_str(), NULL, 16 ) ;
                     dfID = (float)hx_ ;
                     //dfID = atof(vtemp.c_str() ) ;
                  }
                  printf("%02.0f ", dfID ) ;
                  rowV.push_back( dfID ) ;
                  vtemp.clear() ;
                  vidx++ ;
               }
 
           }
	   printf("\n" ) ;
           data.push_back( rowV ) ;
           rowV.clear() ;
        }
  }
  logfile.close() ;

}

void AnaInput::ReadFreeCSV( string fileName, vector<vec>& data, int nSkipLine ) {

  // Open data file
  //printf(" File:%s \n", fileName.c_str() ) ;
  fstream logfile( fileName.c_str() );

  // Open CSV file to read
  if ( !logfile.is_open() )  { 
      printf(" file open error \n") ;
      //return 1 ;
  }
 
  // Variables to hold data 
  int debug = 0 ;

  GetParameters( "debug", &debug ) ;
  // data collection
  vec rowV ;
  rowV.clear() ;

  string  line ;
  int ii =0 ;
  double val_ ;
  while ( getline( logfile, line) ) {

        ii++ ;
        if ( ii <= nSkipLine ) continue ;

        // read data - Die info
	string vtemp ;
	for (string::iterator it = line.begin(); it< line.end(); it++) {
            //if ( (*it) != ',' && (*it) != ' ' ) {
            if ( (*it) != ','  ) {
               vtemp.push_back( *it );
            }
            if ( (*it) == ',' || it == line.end()-1 ) {
               if ( vtemp.size() < 1 ) continue ;
               //printf("%s ", vtemp.c_str() ) ;
               val_ = atof(vtemp.c_str() ) ;
	       if (debug == 1)  printf("%f ", val_ ) ;
	       rowV.push_back( val_ ) ;
	       vtemp.clear() ;
            }
        }
	if (debug == 1) printf("\n" ) ;
        //printf(" row size = %d \n", (int)rowV.size() ) ;
	data.push_back( rowV ) ;
	rowV.clear() ;
  }
  logfile.close() ;

}

vector<vec> AnaInput::Sort2D( vector<vec>& data, int axis1, int axis2 ) {

    vector<vec> data1 ;

    // Sort the data - Example
    iis = axis1 ;
    if ( data.size() > 1 ) {
       sort( data.begin(), data.end(), SortI );
    }

    vector<vec> col ;
    col.clear() ;
    bool nextCol = false ;
    double x_ = data[0][axis1] ;
    iis = axis2 ;

    for (size_t i=0; i< data.size(); i++ ) {
        // Accumulate the same column data
        if ( fabs(data[i][axis1] - x_ ) < 0.1 ) {
           col.push_back( data[i] ) ;
           if ( i == data.size() -1 ) nextCol = true ;
        } else {
           nextCol = true ;
        }

        // Sort the column
        if ( nextCol ) {
           cout<<" col size = "<< col.size() <<endl ;
           sort( col.begin(), col.end(), SortI ) ;

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

