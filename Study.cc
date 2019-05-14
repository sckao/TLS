#include "Study.h"
#include "Reader.h"

Study::Study( ) {

  Input = AnaInput::Instance() ;

  Input->GetParameters("PlotType",      &plotType ) ;
  Input->GetParameters("Path",          &hfolder ) ;
  Input->GetParameters("CSVDIR",        &cfolder ) ;
  Input->GetParameters("CSVFile",       &cFileName );

  Input->GetParameters("MapFile",       &mapFileName );
  Input->GetParameters("MapCorrFile",   &mapCorrFileName );

  Input->GetParameters("HFileName",     &hFileName );
  Input->GetParameters("debug",         &debug );

  Input->GetParameters("PlotName0",     &plotName0 );
  Input->GetParameters("PlotName1",     &plotName1 );

  Input->GetParameters("OutputFile",    &outFileName );
  logFileName = hfolder +  outFileName ;

}

Study::~Study(){

  //delete Input ;
    cout<<" done ! "<<endl ;
}

void Study::VerificationStats() {

   vector<double> targetXY ;
   Input->GetParameters("TargetXY",     &targetXY );


   FILE* logfile = fopen( logFileName.c_str() ,"w");
   printf(" Log File = %s \n", logFileName.c_str() ) ;

   vector<vec> data ;
   for (int i=1; i<8; i++ ) {

       char cuff[32] ;
       sprintf(cuff, "AxisCalibrationResults_Y_Camera%d.csv", i ) ;
       cFileName = cuff ;
       string cvsfile = cfolder + cFileName ;

       ReadCalib( cvsfile, data );
       fprintf(logfile, " file %d , size = %d \n", i, (int)data.size() ) ;
       for (size_t j=0 ; j < data.size(); j++ ) {
           //if ( fabs(data[j][2] - targetXY[0]) < 0.0001 && fabs(data[j][3] - targetXY[1]) < 0.0001 ) {
           if ( fabs(data[j][2] - targetXY[0]) < 2 ) {
              fprintf( logfile, "%f, %f, %f, %f, %f \n",  data[j][2], data[j][3], data[j][19], data[j][20], data[j][26] ) ;
              //break ;
           }
       }
       fprintf(logfile, "========= \n" ) ;
   }


}



// Read the verification of calibration , trim down those invalid measurements
void Study::ReadCalib( string fileName, vector<vec>& data1  ) {

    data1.clear() ;
    vector<vec> data ;
    data.clear() ;

    printf(" CSV file = %s \n", fileName.c_str() ) ;
    Input->ReadFreeCSV( fileName, data, 1 ) ;

    // skim those invalid entries
    for ( size_t i=0; i < data.size(); i++ ) {
        if ( data[i][9] == -1  ) continue ;
        if ( data[i][17] == 0 || data[i][18] == 0 || data[i][19] == 0 || data[i][20] == 0 ) continue ;
        data1.push_back( data[i] ) ;   
    }

}


void Study::FromDPM_Matching() {

   ///// ============ DPM Output ============================================================ 
   //       0    1   2     3     4     5   6    7     8   9   10  11   12     13     14   
   //double x0, y0, xoff, yoff, toff, dx1, dy1, dx2, dy2, xc, yc, dA, diedx, diedy, diedR  ;
   vector<vec> data ;
   Reader  *reader  = new Reader( ) ;
   reader->MatchDPM_SiteCorrMap( data );
   delete reader ;
 
  // Book histograms
  //Open Root file to store result
  string hfile = hfolder + hFileName + ".root" ;
  TFile* rootFile = new TFile( hfile.c_str() , "RECREATE" );
  rootFile->cd() ;

  // Define histograms
  TH1D* h_Xoff = new TH1D("offsetX" , "X_offset",   20, -500, 500 )  ;
  TH1D* h_Yoff = new TH1D("offsetY" , "Y_offset",   20, -500, 500 )  ;
  TH1D* h_Toff = new TH1D("offsetT" , "T_offset",   50, -25, 25 )    ; 
  
  TH1D* h_dxc = new TH1D("h_dxc" , "dx from FF DPM",  190, -0.8, -0.65 )    ; 
  TH1D* h_dyc = new TH1D("h_dyc" , "dy from FF DPM",  100, 0.3, 0.4 )    ; 

  TH2D* h_r_Toff = new TH2D("r_Toff" , "T_offset",  50, 0, 150,  20, -20, 20 )    ; 
   
  TGraph2D* g_Xoff = new TGraph2D() ;

  TGraph* g_xoff1 = new TGraph();
  TGraph* g_xoff2 = new TGraph();
  TGraph* g_xoff3 = new TGraph();
  TGraph* g_xoff4 = new TGraph();
  TGraph* g_xoff5 = new TGraph();

  double r= 0 ;
  int sz = (int) data.size() ;
  int k[5] = { 0 } ;
  double theX(0), theY(0) ;
  for (int i=0 ; i< sz; i++ ) {

      r = sqrt( (data[i][0]*data[i][0]) + (data[i][1]*data[i][1]) ) ;

      h_Xoff->Fill(data[i][2]*1000000 ) ;
      h_Yoff->Fill(data[i][3]*1000000 ) ;
      h_Toff->Fill(data[i][4]*1000000 ) ;
      h_dxc->Fill( data[i][0] - data[i][9] ) ;
      h_dyc->Fill( data[i][1] - data[i][10] ) ;
      h_r_Toff->Fill( r, data[i][4]*1000000 ) ;
      //printf(" %.6f , %.6f , %.6f \n", data[i][3], data[i][4], data[i][5] ) ;

      g_Xoff->SetPoint( i , data[i][0], data[i][1], data[i][2]*1000000 );

      theX = data[i][0] ;
      theY = data[i][1] - data[i][10]  ;

      if ( data[i][1] == -120 ) { 
         g_xoff1->SetPoint(k[0], theX , theY ) ;
         k[0]++ ;
      }
      if ( data[i][1] == -100 )  { 
         g_xoff2->SetPoint(k[1], theX , theY ) ;
         k[1]++ ;
      }
      if ( data[i][1] ==  -40 ) {
         g_xoff3->SetPoint(k[2], theX , theY ) ;
         k[2]++ ;
      } 
      if ( data[i][1] ==  -80 ) { 
         g_xoff4->SetPoint(k[3], theX , theY ) ; 
         k[3]++ ;
      }
      if ( data[i][1] ==  -60 ) { 
         g_xoff5->SetPoint(k[4], theX , theY ) ;
         k[4]++ ; 
      }

  }  


  TCanvas* c0 = new TCanvas("c_0","", 0,0, 800, 600);
  c0->SetFillColor(10);
  c0->SetFillColor(10);
  gPad->SetGridx();

  gStyle->SetPalette(1);
  gStyle->SetNumberContours(10);
  g_Xoff->Draw("surf1z");
  c0->Update();
  TString gPlotname1 = hfolder + "Xoff.png"  ;
  c0->Print( gPlotname1 ) ;

  TString gPlotname2 = hfolder + "dx1_1D.png"  ;
  c0->Clear() ;
  //g_xoff3->GetXaxis()->SetLimits( 30., 150.);
  g_xoff3->SetMaximum(0.4) ;
  g_xoff3->SetMinimum(0.3) ;
  g_xoff3->SetMarkerStyle(20) ;
  g_xoff3->SetMarkerSize(1) ;
  g_xoff3->SetMarkerColor(1) ;
  g_xoff3->SetLineColor(1) ;
  g_xoff3->SetLineWidth(2) ;
  g_xoff3->Draw("ALP") ;
  c0->Update() ;
  g_xoff3->Write() ;
 
  g_xoff1->SetMarkerStyle(21) ;
  g_xoff1->SetMarkerSize(1) ;
  g_xoff1->SetMarkerColor(2) ;
  g_xoff1->SetLineColor(2) ;
  g_xoff1->SetLineWidth(2) ;
  g_xoff1->Draw("SAMELP") ;
  c0->Update() ;
  g_xoff1->Write() ;

  g_xoff2->SetMarkerStyle(21) ;
  g_xoff2->SetMarkerSize(1) ;
  g_xoff2->SetMarkerColor(kOrange) ;
  g_xoff2->SetLineColor(kOrange) ;
  g_xoff2->SetLineWidth(2) ;
  g_xoff2->Draw("SAMELP") ;
  c0->Update() ;
  g_xoff2->Write() ;

  g_xoff4->SetMarkerStyle(22) ;
  g_xoff4->SetMarkerSize(1) ;
  g_xoff4->SetMarkerColor(8) ;
  g_xoff4->SetLineColor(8) ;
  g_xoff4->SetLineWidth(2) ;
  g_xoff4->Draw("SAMELP") ;
  c0->Update() ;
  g_xoff4->Write() ;

  g_xoff5->SetMarkerStyle(21) ;
  g_xoff5->SetMarkerSize(1) ;
  g_xoff5->SetMarkerColor(4) ;
  g_xoff5->SetLineColor(4) ;
  g_xoff5->SetLineWidth(2) ;
  g_xoff5->Draw("SAMELP") ;
  c0->Update() ;
  g_xoff5->Write() ;
  c0->Print( gPlotname2 ) ;



  delete c0 ;
 
  h_Xoff->Write() ;
  h_Yoff->Write() ;
  h_Toff->Write() ;
  h_dxc->Write() ;
  h_dyc->Write() ;
  h_r_Toff->Write() ;
  
  g_Xoff->Write() ;
  rootFile->Close() ;
}

void Study::FromCDM_Matching( ) {

   ///// ============ CDM Output ============================================================ 
   //       0    1   2     3     4     5   6    7     8    9   10   11   12   13   14   
   //double x0, y0, xoff, yoff, toff, dx1, dy1, dx2, dy2, ds1, x1,  y1, ds2,  x2,  y2;
   vector<vec> data ;
   Reader  *reader  = new Reader( ) ;
   reader->GetDataFromSiteCorrectionMap( data );
   delete reader ;
 
  // Book histograms
  //Open Root file to store result
  string hfile = hfolder + hFileName + ".root" ;
  TFile* rootFile = new TFile( hfile.c_str() , "RECREATE" );
  rootFile->cd() ;

  // Define histograms
  TH1D* h_Xoff = new TH1D("offsetX" , "X_offset",   20, -500, 500 )  ;
  TH1D* h_Yoff = new TH1D("offsetY" , "Y_offset",   20, -500, 500 )  ;
  TH1D* h_Toff = new TH1D("offsetT" , "T_offset",   50, -25, 25 )    ; 

  TH2D* h_r_Toff = new TH2D("r_Toff" , "T_offset",  50, 0, 150,  20, -20, 20 )    ; 
   
  TGraph2D* g_Xoff = new TGraph2D() ;

  TGraph* g_xoff1 = new TGraph();
  TGraph* g_xoff2 = new TGraph();
  TGraph* g_xoff3 = new TGraph();
  TGraph* g_xoff4 = new TGraph();
  TGraph* g_xoff5 = new TGraph();

  double r= 0 ;
  int sz = (int) data.size() ;
  int k[5] = { 0 } ;
  double theX(0), theY(0) ;
  for (int i=0 ; i< sz; i++ ) {

      r = sqrt( (data[i][0]*data[i][0]) + (data[i][1]*data[i][1]) ) ;

      h_Xoff->Fill(data[i][2]*1000000 ) ;
      h_Yoff->Fill(data[i][3]*1000000 ) ;
      h_Toff->Fill(data[i][4]*1000000 ) ;
      h_r_Toff->Fill( r, data[i][4]*1000000 ) ;
      //printf(" %.6f , %.6f , %.6f \n", data[i][3], data[i][4], data[i][5] ) ;

      g_Xoff->SetPoint( i , data[i][0], data[i][1], data[i][2]*1000000 );

      theX = data[i][0] ;
      theY = data[i][0] - 5.0 - data[i][10]  ;

      if ( data[i][1] == -120 ) { 
         g_xoff1->SetPoint(k[0], theX , theY ) ;
         k[0]++ ;
      }
      if ( data[i][1] == -100 )  { 
         g_xoff2->SetPoint(k[1], theX , theY ) ;
         k[1]++ ;
      }
      if ( data[i][1] ==  -40 ) {
         g_xoff3->SetPoint(k[2], theX , theY ) ;
         k[2]++ ;
      } 
      if ( data[i][1] ==  -80 ) { 
         g_xoff4->SetPoint(k[3], theX , theY ) ; 
         k[3]++ ;
      }
      if ( data[i][1] ==  -60 ) { 
         g_xoff5->SetPoint(k[4], theX , theY ) ;
         k[4]++ ; 
      }

  }  


  TCanvas* c0 = new TCanvas("c_0","", 0,0, 800, 600);
  c0->SetFillColor(10);
  c0->SetFillColor(10);
  gPad->SetGridx();

  gStyle->SetPalette(1);
  gStyle->SetNumberContours(10);
  g_Xoff->Draw("surf1z");
  c0->Update();
  TString gPlotname1 = hfolder + "Xoff.png"  ;
  c0->Print( gPlotname1 ) ;

  TString gPlotname2 = hfolder + "CDM_dx1_1D.png"  ;
  c0->Clear() ;
  //g_xoff3->GetXaxis()->SetLimits( 30., 150.);
  g_xoff3->SetMaximum(-0.45) ;
  g_xoff3->SetMinimum(-0.65) ;
  g_xoff3->SetMarkerStyle(20) ;
  g_xoff3->SetMarkerSize(1) ;
  g_xoff3->SetMarkerColor(1) ;
  g_xoff3->SetLineColor(1) ;
  g_xoff3->SetLineWidth(2) ;
  g_xoff3->Draw("ALP") ;
  c0->Update() ;
  g_xoff3->Write() ;
 
  g_xoff1->SetMarkerStyle(21) ;
  g_xoff1->SetMarkerSize(1) ;
  g_xoff1->SetMarkerColor(2) ;
  g_xoff1->SetLineColor(2) ;
  g_xoff1->SetLineWidth(2) ;
  g_xoff1->Draw("SAMELP") ;
  c0->Update() ;
  g_xoff1->Write() ;

  g_xoff2->SetMarkerStyle(21) ;
  g_xoff2->SetMarkerSize(1) ;
  g_xoff2->SetMarkerColor(kOrange) ;
  g_xoff2->SetLineColor(kOrange) ;
  g_xoff2->SetLineWidth(2) ;
  g_xoff2->Draw("SAMELP") ;
  c0->Update() ;
  g_xoff2->Write() ;

  g_xoff4->SetMarkerStyle(22) ;
  g_xoff4->SetMarkerSize(1) ;
  g_xoff4->SetMarkerColor(8) ;
  g_xoff4->SetLineColor(8) ;
  g_xoff4->SetLineWidth(2) ;
  g_xoff4->Draw("SAMELP") ;
  c0->Update() ;
  g_xoff4->Write() ;

  g_xoff5->SetMarkerStyle(21) ;
  g_xoff5->SetMarkerSize(1) ;
  g_xoff5->SetMarkerColor(4) ;
  g_xoff5->SetLineColor(4) ;
  g_xoff5->SetLineWidth(2) ;
  g_xoff5->Draw("SAMELP") ;
  c0->Update() ;
  g_xoff5->Write() ;
  c0->Print( gPlotname2 ) ;

  delete c0 ;
 
  h_Xoff->Write() ;
  h_Yoff->Write() ;
  h_Toff->Write() ;
  h_r_Toff->Write() ;
  
  g_Xoff->Write() ;
  rootFile->Close() ;
}

/*
void Study::MindR() {

   ///// ============ CDM Output ============================================================ 
   //       0    1   2     3     4     5   6    7     8    9   10   11   12   13   14   
   //double x0, y0, xoff, yoff, toff, dx1, dy1, dx2, dy2, ds1, x1,  y1, ds2,  x2,  y2;
   vector<vec> data ;
   Reader  *reader  = new Reader( ) ;
   reader->GetDataFromSiteCorrectionMap( data );
   delete reader ;



}
*/


void Study::MatrixTest() {

    TMatrixD a(2,2) ;
    a[0][0] = 25 ;
    a[0][1] = -15 ;
    a[1][0] = -15 ;
    a[1][1] = 25 ;

    TMatrixDEigen e(a) ;

    TMatrixD eM = e.GetEigenValues() ;
  
    printf(" eigen value = %f , %f  -  %f %f  \n", eM[0][0], eM[1][1], eM[0][1], eM[1][0] ) ;
}


