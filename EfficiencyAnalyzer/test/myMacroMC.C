#include "tdrstyle.C"
#include "CMS_lumi.C"
#include "TH1.h"
#include "TH1F.h"
 
TCanvas* muHadMass_MC_plot( int iPeriod, int iPos, const char* MTBin );
TCanvas* HLT_eff_MC_plot( int iPeriod, int iPos, const char* prodMode );
TCanvas* tau_eff_MC_plot( int iPeriod, int iPos, const char* var );
 
void myMacroMC(const char* MTBin, const char* prodMode)
{
  // gROOT->LoadMacro("tdrstyle.C");
  setTDRStyle();

  // gROOT->LoadMacro("CMS_lumi.C");

  writeExtraText = true;       // if extra text
  extraText  = "Simulation";  // default extra text is "Preliminary"
  // lumi_8TeV  = "19.1 fb^{-1}"; // default is "19.7 fb^{-1}"
  // lumi_7TeV  = "4.9 fb^{-1}";  // default is "5.1 fb^{-1}"
  // lumi_sqrtS = "13 TeV";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

  int iPeriod = 12;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 0=free form (uses lumi_sqrtS)

  // second parameter in muHadMass_plot is iPos, which drives the position of the CMS logo in the plot
  // iPos=11 : top-left, left-aligned
  // iPos=33 : top-right, right-aligned
  // iPos=22 : center, centered
  // mode generally : 
  //   iPos = 10*(alignement 1/2/3) + position (1/2/3 = left/center/right)

  // muHadMass_MC_plot( iPeriod, 0, MTBin );   // out of frame (in exceptional cases)
  // HLT_eff_MC_plot( iPeriod, 0, prodMode );   // out of frame (in exceptional cases)
  tau_eff_MC_plot( iPeriod, 0, "PT" );   // out of frame (in exceptional cases)
  tau_eff_MC_plot( iPeriod, 0, "Eta" );   // out of frame (in exceptional cases)
  // muHadMass_MC_plot( iPeriod, 11, MTBin );  // left-aligned
  //  muHadMass_plot( iPeriod, 33 );  // right-aligned

  //  writeExtraText = false;       // remove Preliminary
  
  //  muHadMass_plot( iPeriod, 0 );   // out of frame (in exceptional cases)

  //  muHadMass_plot( iPeriod, 11 );  // default: left-aligned
  //  muHadMass_plot( iPeriod, 22 );  // centered
  //  muHadMass_plot( iPeriod, 33 );  // right-aligned  
}

TCanvas* muHadMass_MC_plot( int iPeriod, int iPos, const char* MTBin )
{ 
  //  if( iPos==0 ) relPosX = 0.12;

  int W = 800;
  int H = 600;

  // 
  // Simple example of macro: plot with CMS name and lumi text
  //  (this script does not pretend to work in all configurations)
  // iPeriod = 1*(0/1 7 TeV) + 2*(0/1 8 TeV)  + 4*(0/1 13 TeV) 
  // For instance: 
  //               iPeriod = 3 means: 7 TeV + 8 TeV
  //               iPeriod = 7 means: 7 TeV + 8 TeV + 13 TeV 
  // Initiated by: Gautier Hamel de Monchenault (Saclay)
  // Updated by:   Dinko Ferencek (Rutgers)
  //
  int H_ref = 600; 
  int W_ref = 800; 

  // references for T, B, L, R
  float T = 0.08*H_ref;
  float B = 0.12*H_ref; 
  float L = 0.12*W_ref;
  float R = 0.04*W_ref;

  TString canvName = "sigVsBkg_muHadMass_";
  // canvName += W;
  // canvName += "-";
  // canvName += H;
  // canvName += "_";  
  // canvName += iPeriod;
  // if( writeExtraText ) canvName += "-prelim";
  // if( iPos%10==0 ) canvName += "-out";
  // else if( iPos%10==1 ) canvName += "-left";
  // else if( iPos%10==2 )  canvName += "-center";
  // else if( iPos%10==3 )  canvName += "-right";
  canvName += MTBin;
  canvName += "_v87";

  TCanvas* canv = new TCanvas(canvName,canvName,50,50,W,H);
  canv->SetFillColor(0);
  canv->SetBorderMode(0);
  canv->SetFrameFillStyle(0);
  canv->SetFrameBorderMode(0);
  canv->SetLeftMargin( L/W );
  canv->SetRightMargin( R/W );
  canv->SetTopMargin( T/H );
  canv->SetBottomMargin( B/H );
  canv->SetTickx(0);
  canv->SetTicky(0);
  canv->SetLogy();

  TH1* h = new TH1F("h","h",40,0.0, 11.0);
  h->GetXaxis()->SetNdivisions(6,5,0);
  h->GetXaxis()->SetTitle("m_{#mu+X} (GeV)");  
  h->GetYaxis()->SetNdivisions(6,5,0);
  h->GetYaxis()->SetTitleOffset(1);
  h->GetYaxis()->SetTitle("Events / bin");  

  h->SetMaximum( 260 );
  if( iPos==1 ) h->SetMaximum( 300 );
  h->Draw();
  h->GetYaxis()->SetRangeUser(0.01, 1000000.0);

  int histLineColor = kBlue - 10;
  int histFillColor = kBlue - 10;
  float markerSize  = 1.0;

  {
    TLatex latex;
				
    int nLeft_ = /*2*/4;
    int nRight_ = /*2*/8;

    float x1_l = /*0.92*/0.95;
    float y1_l = /*0.60*//*0.9*//*0.85*/0.9175;

    float dx_l = 0.30;
    float dy_l = /*0.18*//*0.54*//*0.3*/0.45;
    float x0_l = x1_l-dx_l;
    float x0_l_legendLeft = x1_l-2.0*dx_l;
    float y0_l = y1_l-dy_l;

    float ar_l = dy_l/dx_l;
		
    float x_l[1];
    float ex_l[1];
    float y_l_left[1];
    float ey_l_left[1];
    float y_l_right[1];
    float ey_l[1];
		
    //    float gap_ = 0.09/ar_l;
    float gapLeft_ = 1./(nLeft_+1);
    float gapRight_ = 1./(nRight_+1);
		
    float bwx_ = 0.12;
    float bwyLeft_ = gapLeft_/1.5;
    float bwyRight_ = gapRight_/1.5;
		
    x_l[0] = 1.2*bwx_;
    //    y_l_left[0] = 1-(1-0.10)/ar_l;
    y_l_left[0] = 1-gapLeft_;
    y_l_right[0] = 1-gapRight_;
    ex_l[0] = 0;
    ey_l[0] = 0.04/ar_l;

    float xx_ = x_l[0];
    float yyLeft_ = y_l_left[0];

    TPad* legendLeft = new TPad("legendLeft_0","legendLeft_0",
				/*x0_l*/x0_l_legendLeft,y0_l,/*x1_l*/x0_l+0.0*dx_l, y1_l );
    // legendLeft->SetFillColor( kGray );
    legendLeft->Draw();
    legendLeft->cd();

    // TGraph* gr_l = new TGraphErrors(1, x_l, y_l_left, ex_l, ey_l );
		
    // gStyle->SetEndErrorSize(0);
    // gr_l->SetMarkerSize(0.9);
    // gr_l->Draw("0P");
		
    latex.SetTextFont(42);
    latex.SetTextAngle(0);
    latex.SetTextColor(kBlack);    
    latex.SetTextSize(/*0.25*/0.09375);    
    latex.SetTextAlign(12); 
		
    TLine line_;
    TBox  box_;
    // latex.DrawLatex(xx_+1.*bwx_,yyLeft_,"Data");

    line_.SetLineStyle( kDashed );
    line_.SetLineWidth( 4 );
    line_.SetLineColor( kGreen - 2 );
    line_.DrawLine( xx_-bwx_/2, yyLeft_, xx_+bwx_/2, yyLeft_ );
    latex.DrawLatex(xx_+1.*bwx_,yyLeft_,"ggH m_{a} = 9 GeV");

    yyLeft_ -= gapLeft_;
    line_.SetLineStyle( kDashed );
    line_.SetLineWidth( 4 );
    line_.SetLineColor( kBlue );
    line_.DrawLine( xx_-bwx_/2, yyLeft_, xx_+bwx_/2, yyLeft_ );
    latex.DrawLatex(xx_+1.*bwx_,yyLeft_,"WH m_{a} = 9 GeV");

    yyLeft_ -= gapLeft_;
    line_.SetLineStyle( kDashed );
    line_.SetLineWidth( 4 );
    line_.SetLineColor( kMagenta + 2 );
    line_.DrawLine( xx_-bwx_/2, yyLeft_, xx_+bwx_/2, yyLeft_ );
    latex.DrawLatex(xx_+1.*bwx_,yyLeft_,"VBF m_{a} = 9 GeV");

    yyLeft_ -= gapLeft_;
    line_.SetLineStyle( kDashed );
    line_.SetLineWidth( 4 );
    line_.SetLineColor( kOrange + 7 );
    line_.DrawLine( xx_-bwx_/2, yyLeft_, xx_+bwx_/2, yyLeft_ );
    latex.DrawLatex(xx_+1.*bwx_,yyLeft_,"ZH m_{a} = 9 GeV");

    canv->cd();
	
    TPad* legendRight = new TPad("legendRight_0","legendRight_0",x0_l,y0_l,x1_l, y1_l );
    // legendRight->SetFillColor( kGray );
    legendRight->Draw();
    legendRight->cd();

    float yyRight_ = y_l_right[0];
    box_.SetLineStyle( kSolid );
    box_.SetLineWidth( 0 );
    box_.SetFillColor( kMagenta + 2 );
    box_.DrawBox( xx_-bwx_/2, yyRight_-bwyRight_/2, xx_+bwx_/2, yyRight_+bwyRight_/2 );
    box_.SetFillStyle(1001);
    box_.DrawBox( xx_-bwx_/2, yyRight_-bwyRight_/2, xx_+bwx_/2, yyRight_+bwyRight_/2 );
    latex.DrawLatex(xx_+1.*bwx_,yyRight_,"WW");

    yyRight_ -= gapRight_;
    box_.SetLineStyle( kSolid );
    box_.SetLineWidth( 0 );
    box_.SetFillColor( kCyan + 2 );
    box_.DrawBox( xx_-bwx_/2, yyRight_-bwyRight_/2, xx_+bwx_/2, yyRight_+bwyRight_/2 );
    box_.SetFillStyle(1001);
    box_.DrawBox( xx_-bwx_/2, yyRight_-bwyRight_/2, xx_+bwx_/2, yyRight_+bwyRight_/2 );
    latex.DrawLatex(xx_+1.*bwx_,yyRight_,"ZZ");

    yyRight_ -= gapRight_;
    box_.SetLineStyle( kSolid );
    box_.SetLineWidth( 0 );
    box_.SetFillColor( kRed + 2 );
    box_.DrawBox( xx_-bwx_/2, yyRight_-bwyRight_/2, xx_+bwx_/2, yyRight_+bwyRight_/2 );
    box_.SetFillStyle(1001);
    box_.DrawBox( xx_-bwx_/2, yyRight_-bwyRight_/2, xx_+bwx_/2, yyRight_+bwyRight_/2 );
    latex.DrawLatex(xx_+1.*bwx_,yyRight_,"WZ");

    yyRight_ -= gapRight_;
    box_.SetLineStyle( kSolid );
    box_.SetLineWidth( 0 );
    box_.SetFillColor( kYellow );
    box_.DrawBox( xx_-bwx_/2, yyRight_-bwyRight_/2, xx_+bwx_/2, yyRight_+bwyRight_/2 );
    box_.SetFillStyle(1001);
    box_.DrawBox( xx_-bwx_/2, yyRight_-bwyRight_/2, xx_+bwx_/2, yyRight_+bwyRight_/2 );
    latex.DrawLatex(xx_+1.*bwx_,yyRight_,"W + #geq1 jet");

    yyRight_ -= gapRight_;
    box_.SetLineStyle( kSolid );
    box_.SetLineWidth( 0 );
    box_.SetFillColor( kViolet - 7 );
    box_.DrawBox( xx_-bwx_/2, yyRight_-bwyRight_/2, xx_+bwx_/2, yyRight_+bwyRight_/2 );
    box_.SetFillStyle(1001);
    box_.DrawBox( xx_-bwx_/2, yyRight_-bwyRight_/2, xx_+bwx_/2, yyRight_+bwyRight_/2 );
    latex.DrawLatex(xx_+1.*bwx_,yyRight_,"Single top");

    yyRight_ -= gapRight_;
    box_.SetLineStyle( kSolid );
    box_.SetLineWidth( 0 );
    box_.SetFillColor( kSpring + 4 );
    box_.DrawBox( xx_-bwx_/2, yyRight_-bwyRight_/2, xx_+bwx_/2, yyRight_+bwyRight_/2 );
    box_.SetFillStyle(1001);
    box_.DrawBox( xx_-bwx_/2, yyRight_-bwyRight_/2, xx_+bwx_/2, yyRight_+bwyRight_/2 );
    latex.DrawLatex(xx_+1.*bwx_,yyRight_,"t#bar{t} + jets");

    yyRight_ -= gapRight_;
    box_.SetLineStyle( kSolid );
    box_.SetLineWidth( 0 );
    box_.SetFillColor( kBlue + 1 );
    box_.DrawBox( xx_-bwx_/2, yyRight_-bwyRight_/2, xx_+bwx_/2, yyRight_+bwyRight_/2 );
    box_.SetFillStyle(1001);
    box_.DrawBox( xx_-bwx_/2, yyRight_-bwyRight_/2, xx_+bwx_/2, yyRight_+bwyRight_/2 );
    latex.DrawLatex(xx_+1.*bwx_,yyRight_,"Drell-Yan + jets");

    yyRight_ -= gapRight_;
    box_.SetLineStyle( kSolid );
    box_.SetLineWidth( 0 );
    box_.SetFillColor( kGray + 2 );
    box_.DrawBox( xx_-bwx_/2, yyRight_-bwyRight_/2, xx_+bwx_/2, yyRight_+bwyRight_/2 );
    box_.SetFillStyle(1001);
    box_.DrawBox( xx_-bwx_/2, yyRight_-bwyRight_/2, xx_+bwx_/2, yyRight_+bwyRight_/2 );
    latex.DrawLatex(xx_+1.*bwx_,yyRight_,"QCD (from data)");

    canv->cd();

    TLatex latexTitle;
    latexTitle.SetTextFont(42);
    latexTitle.SetTextAngle(0);
    latexTitle.SetTextColor(kBlack);    
    latexTitle.SetTextSize(0.0524476);    
    latexTitle.SetTextAlign(12);     
    Double_t userX = 0.814011;
    Double_t userY = 9158;
    string label;
    if (string(MTBin) == "lowMT") label = "Low M_{T}";
    if (string(MTBin) == "highMT") label = "High M_{T}";
    latexTitle.DrawLatex(userX,userY,label.c_str());
  }

  {
    // Observed data
    string fileName("../sigVsBkgQCDFromData_muHadIsoAnalysis_");
    fileName = fileName + string(MTBin);
    fileName = fileName + string("_a9_19p7fb-1_v214.root");
    TFile file_(fileName.c_str(),"READ");

    TCanvas *canvas   = static_cast<TCanvas*>(file_.Get("muHadMassCanvas")->Clone());

    canvas->Draw("goff");
    THStack *stack   =
      static_cast<THStack*>(canvas->cd(1)->GetPrimitive("muHadMassStack")->Clone());
    stack->SetMinimum(0.01);
    stack->SetMaximum(1000000.0);
    canv->cd();
    stack->Draw("same");
    
    TList* sigs = canvas->GetListOfPrimitives();
    for (Int_t i = 2; i < 6; ++i) {
      TH1F *sig   = static_cast<TH1F*>(sigs->At(i)->Clone());
      sig->SetDirectory(0);
      int sigCol = kBlack;
      switch (i) {
      case 2:
      	sigCol = kBlue;
	break;
      case 3:
      	sigCol = kGreen - 2;
	break;
      case 4:
      	sigCol = kOrange + 7;
	break;
      case 5:
      	sigCol = kMagenta + 2;
	break;
      default:
      	break;
      }
      canv->cd();
      sig->SetLineColor(sigCol);
      sig->Draw("histsame");
    }

    file_.Close();
  }

  // writing the lumi information and the CMS "logo"
  CMS_lumi( canv, iPeriod, iPos );

  canv->Update();
  canv->RedrawAxis();
  canv->GetFrame()->Draw();

  canv->Print(canvName+".pdf",".pdf");
  canv->Print(canvName+".png",".png");

  return canv;
}

TCanvas* HLT_eff_MC_plot( int iPeriod, int iPos, const char* prodMode )
{ 
  //  if( iPos==0 ) relPosX = 0.12;

  int W = /*(1.0+0.115)**/800;
  int H = 600;

  // 
  // Simple example of macro: plot with CMS name and lumi text
  //  (this script does not pretend to work in all configurations)
  // iPeriod = 1*(0/1 7 TeV) + 2*(0/1 8 TeV)  + 4*(0/1 13 TeV) 
  // For instance: 
  //               iPeriod = 3 means: 7 TeV + 8 TeV
  //               iPeriod = 7 means: 7 TeV + 8 TeV + 13 TeV 
  // Initiated by: Gautier Hamel de Monchenault (Saclay)
  // Updated by:   Dinko Ferencek (Rutgers)
  //
  int H_ref = 600; 
  int W_ref = 800; 

  // references for T, B, L, R
  float T = 0.08*H_ref;
  float B = 0.12*H_ref; 
  float L = 0.12*W_ref;
  float R = 0.04/*0.115*/*W_ref;

  TString canvName = "HLT_eff_Wmu_";
  // canvName += W;
  // canvName += "-";
  // canvName += H;
  // canvName += "_";  
  // canvName += iPeriod;
  // if( writeExtraText ) canvName += "-prelim";
  // if( iPos%10==0 ) canvName += "-out";
  // else if( iPos%10==1 ) canvName += "-left";
  // else if( iPos%10==2 )  canvName += "-center";
  // else if( iPos%10==3 )  canvName += "-right";
  canvName += prodMode;

  TCanvas* canv = new TCanvas(canvName,canvName,50,50,W,H);
  canv->SetFillColor(0);
  canv->SetBorderMode(0);
  canv->SetFrameFillStyle(0);
  canv->SetFrameBorderMode(0);
  canv->SetLeftMargin( L/W );
  canv->SetRightMargin( R/W );
  canv->SetTopMargin( T/H );
  canv->SetBottomMargin( B/H );
  canv->SetTickx(0);
  canv->SetTicky(0);
  
  // TPad* padEff = new TPad("padEff", "padEff", 0, 0, 1, 1);
  // padEff->SetFillStyle(4000);
  // padEff->SetLeftMargin( L/W );
  // padEff->SetRightMargin( R/W );
  // padEff->SetTopMargin( T/H );
  // padEff->SetBottomMargin( B/H );

  // canv->cd();
  // padEff->Draw();
  // padEff->cd();
  // canv->cd();
  // padEff->cd();
  // padEff->SetTicks(0, 0);

  TH1* hEff = new TH1F("hEff","hEff",40,0.0, 100.0);
  hEff->GetXaxis()->SetNdivisions(6,5,0);
  hEff->GetXaxis()->SetTitle("p_{T} (GeV)");  
  hEff->GetYaxis()->SetNdivisions(6,5,0);
  hEff->GetYaxis()->SetTitleOffset(1/*0.9*/);
  hEff->GetYaxis()->SetTitle("Simulated efficiency");  

  hEff->SetMaximum( 260 );
  if( iPos==1 ) hEff->SetMaximum( 300 );
  hEff->GetYaxis()->SetRangeUser(0.0, 1.0);
  hEff->Draw();

  // padEff->SetTicks(0, 0);
  // padEff->Update();
  // canv->cd();

  // TPad* padDist = new TPad("padDist", "padDist", 0, 0, 1, 1);
  // padDist->SetFillStyle(4000);
  // padDist->SetLeftMargin( L/W );
  // padDist->SetRightMargin( R/W );
  // padDist->SetTopMargin( T/H );
  // padDist->SetBottomMargin( B/H );
  
  // padDist->Draw();
  // padDist->cd();
  // padDist->SetTicks(0, 0);

  // TH1* hDist = new TH1F("hDist","hDist",40,0.0, 100.0);
  // hDist->GetXaxis()->SetNdivisions(6,5,0);
  // hDist->GetXaxis()->SetTitle("p_{T} (GeV)");  
  // hDist->GetYaxis()->SetNdivisions(6,5,0);
  // hDist->GetYaxis()->SetTitleOffset(/*1*/0.9);
  // hDist->GetYaxis()->SetTitle("Events / 5 GeV");  

  // Double_t yAxisMax = 0.0;
  // if (string(prodMode) == "WH") yAxisMax = 120.0;
  // if (string(prodMode) == "ggH") yAxisMax = 1300.0;
  
  // hDist->SetMaximum( 260 );
  // if( iPos==1 ) hDist->SetMaximum( 300 );
  // hDist->Draw("Y+");
  // hDist->GetYaxis()->SetRangeUser(0.0, yAxisMax);

  // padDist->SetTicks(0, 0);

  int histLineColor = kBlue - 10;
  int histFillColor = kBlue - 10;
  float markerSize  = 1.0;

  {
    TLatex latex;
				
    int nLeft_ = 2;
    int nRight_ = /*2*/8;

    float x1_l = 0.65;
    float y1_l = 0.4;

    float dx_l = 0.30;
    float dy_l = 0.18;
    float x0_l = x1_l-dx_l;
    float x0_l_legendLeft = x1_l-2.0*dx_l;
    float y0_l = y1_l-dy_l;

    float ar_l = dy_l/dx_l;
		
    float x_l[1];
    float ex_l[1];
    float y_l_left[1];
    float ey_l_left[1];
    float y_l_right[1];
    float ey_l[1];
		
    //    float gap_ = 0.09/ar_l;
    float gapLeft_ = 1./(nLeft_+1);
    float gapRight_ = 1./(nRight_+1);
		
    float bwx_ = 0.12;
    float bwyLeft_ = gapLeft_/1.5;
    float bwyRight_ = gapRight_/1.5;
		
    x_l[0] = 1.2*bwx_;
    //    y_l_left[0] = 1-(1-0.10)/ar_l;
    y_l_left[0] = 1-gapLeft_;
    y_l_right[0] = 1-gapRight_;
    ex_l[0] = 0;
    ey_l[0] = 0.04/ar_l;

    float xx_ = x_l[0];
    float yyLeft_ = y_l_left[0];

    TPad* legendLeft = new TPad("legendLeft_0","legendLeft_0",
  				x0_l/*x0_l_legendLeft*/,y0_l,x1_l/*x0_l+0.0*dx_l*/, y1_l );
    // legendLeft->SetFillColor( kGray );
    // legendLeft->Draw();
    // legendLeft->cd();

    // TGraph* gr_l = new TGraphErrors(1, x_l, y_l_left, ex_l, ey_l );
		
    // gStyle->SetEndErrorSize(0);
    // gr_l->SetMarkerSize(0.9);
    // gr_l->Draw("0P");
		
    // latex.SetTextFont(42);
    // latex.SetTextAngle(0);
    // latex.SetTextColor(kBlack);    
    // latex.SetTextSize(0.25/*0.09375*/);    
    // latex.SetTextAlign(12); 
		
    // TLine line_;
    // TBox  box_;
    // latex.DrawLatex(xx_+1.*bwx_,yyLeft_,"Efficiency");

    // yyLeft_ -= gapLeft_;
    // line_.SetLineStyle( kSolid );
    // line_.SetLineWidth( 2 );
    // line_.SetLineColor( kBlack );
    // line_.DrawLine( xx_-bwx_/2, yyLeft_, xx_+bwx_/2, yyLeft_ );
    // latex.DrawLatex(xx_+1.*bwx_,yyLeft_,"p_{T} distribution");

    canv->cd();
	
    // TPad* legendRight = new TPad("legendRight_0","legendRight_0",x0_l,y0_l,x1_l, y1_l );
    // // legendRight->SetFillColor( kGray );
    // legendRight->Draw();
    // legendRight->cd();

    // float yyRight_ = y_l_right[0];
    // box_.SetLineStyle( kSolid );
    // box_.SetLineWidth( 0 );
    // box_.SetFillColor( kMagenta + 2 );
    // box_.DrawBox( xx_-bwx_/2, yyRight_-bwyRight_/2, xx_+bwx_/2, yyRight_+bwyRight_/2 );
    // box_.SetFillStyle(1001);
    // box_.DrawBox( xx_-bwx_/2, yyRight_-bwyRight_/2, xx_+bwx_/2, yyRight_+bwyRight_/2 );
    // latex.DrawLatex(xx_+1.*bwx_,yyRight_,"WW");

    // yyRight_ -= gapRight_;
    // box_.SetLineStyle( kSolid );
    // box_.SetLineWidth( 0 );
    // box_.SetFillColor( kCyan + 2 );
    // box_.DrawBox( xx_-bwx_/2, yyRight_-bwyRight_/2, xx_+bwx_/2, yyRight_+bwyRight_/2 );
    // box_.SetFillStyle(1001);
    // box_.DrawBox( xx_-bwx_/2, yyRight_-bwyRight_/2, xx_+bwx_/2, yyRight_+bwyRight_/2 );
    // latex.DrawLatex(xx_+1.*bwx_,yyRight_,"ZZ");

    // yyRight_ -= gapRight_;
    // box_.SetLineStyle( kSolid );
    // box_.SetLineWidth( 0 );
    // box_.SetFillColor( kRed + 2 );
    // box_.DrawBox( xx_-bwx_/2, yyRight_-bwyRight_/2, xx_+bwx_/2, yyRight_+bwyRight_/2 );
    // box_.SetFillStyle(1001);
    // box_.DrawBox( xx_-bwx_/2, yyRight_-bwyRight_/2, xx_+bwx_/2, yyRight_+bwyRight_/2 );
    // latex.DrawLatex(xx_+1.*bwx_,yyRight_,"WZ");

    // yyRight_ -= gapRight_;
    // box_.SetLineStyle( kSolid );
    // box_.SetLineWidth( 0 );
    // box_.SetFillColor( kYellow );
    // box_.DrawBox( xx_-bwx_/2, yyRight_-bwyRight_/2, xx_+bwx_/2, yyRight_+bwyRight_/2 );
    // box_.SetFillStyle(1001);
    // box_.DrawBox( xx_-bwx_/2, yyRight_-bwyRight_/2, xx_+bwx_/2, yyRight_+bwyRight_/2 );
    // latex.DrawLatex(xx_+1.*bwx_,yyRight_,"W + #geq1 jet");

    // yyRight_ -= gapRight_;
    // box_.SetLineStyle( kSolid );
    // box_.SetLineWidth( 0 );
    // box_.SetFillColor( kViolet - 7 );
    // box_.DrawBox( xx_-bwx_/2, yyRight_-bwyRight_/2, xx_+bwx_/2, yyRight_+bwyRight_/2 );
    // box_.SetFillStyle(1001);
    // box_.DrawBox( xx_-bwx_/2, yyRight_-bwyRight_/2, xx_+bwx_/2, yyRight_+bwyRight_/2 );
    // latex.DrawLatex(xx_+1.*bwx_,yyRight_,"Single top");

    // yyRight_ -= gapRight_;
    // box_.SetLineStyle( kSolid );
    // box_.SetLineWidth( 0 );
    // box_.SetFillColor( kSpring + 4 );
    // box_.DrawBox( xx_-bwx_/2, yyRight_-bwyRight_/2, xx_+bwx_/2, yyRight_+bwyRight_/2 );
    // box_.SetFillStyle(1001);
    // box_.DrawBox( xx_-bwx_/2, yyRight_-bwyRight_/2, xx_+bwx_/2, yyRight_+bwyRight_/2 );
    // latex.DrawLatex(xx_+1.*bwx_,yyRight_,"t#bar{t} + jets");

    // yyRight_ -= gapRight_;
    // box_.SetLineStyle( kSolid );
    // box_.SetLineWidth( 0 );
    // box_.SetFillColor( kBlue + 1 );
    // box_.DrawBox( xx_-bwx_/2, yyRight_-bwyRight_/2, xx_+bwx_/2, yyRight_+bwyRight_/2 );
    // box_.SetFillStyle(1001);
    // box_.DrawBox( xx_-bwx_/2, yyRight_-bwyRight_/2, xx_+bwx_/2, yyRight_+bwyRight_/2 );
    // latex.DrawLatex(xx_+1.*bwx_,yyRight_,"Drell-Yan + jets");

    // yyRight_ -= gapRight_;
    // box_.SetLineStyle( kSolid );
    // box_.SetLineWidth( 0 );
    // box_.SetFillColor( kGray + 2 );
    // box_.DrawBox( xx_-bwx_/2, yyRight_-bwyRight_/2, xx_+bwx_/2, yyRight_+bwyRight_/2 );
    // box_.SetFillStyle(1001);
    // box_.DrawBox( xx_-bwx_/2, yyRight_-bwyRight_/2, xx_+bwx_/2, yyRight_+bwyRight_/2 );
    // latex.DrawLatex(xx_+1.*bwx_,yyRight_,"QCD (from data)");

    // canv->cd();

    TLatex latexTitle;
    latexTitle.SetTextFont(42);
    latexTitle.SetTextAngle(0);
    latexTitle.SetTextColor(kBlack);    
    latexTitle.SetTextSize(0.0524476);    
    latexTitle.SetTextAlign(12);
    latexTitle.SetNDC();
    Double_t userX = 0.135;
    Double_t userY = 0.875;
    string label;
    if (string(prodMode) == "ggH") label = "gg#rightarrowH#rightarrowaa";
    if (string(prodMode) == "WH") label = "pp#rightarrowWH#rightarrow#mu#nuaa";
    latexTitle.DrawLatex(userX,userY,label.c_str());
  }

  {
    // Observed data
    string fileName("../");
    fileName = fileName + string(prodMode);
    fileName = fileName + string("WMuHLTEff.root");
    TFile file_(fileName.c_str(),"READ");

    TH1F *denominator   = static_cast<TH1F*>(file_.Get("denominatorPT")->Clone());
    denominator->SetDirectory(0);
    denominator->SetLineColor(kBlack);
    denominator->SetFillStyle(0);

    TH1F *numerator   = static_cast<TH1F*>(file_.Get("numeratorPT")->Clone());
    numerator->SetDirectory(0);
    numerator->SetLineColor(kBlack);

    TGraphAsymmErrors *eff = new TGraphAsymmErrors(numerator, denominator);
    eff->SetMarkerStyle(20);
    eff->SetMarkerSize(0.9);

    for (Int_t i = 0; i < eff->GetN(); ++i) {
      eff->SetPointEXhigh(i, 0.0);
      eff->SetPointEXlow(i, 0.0);
    }
    
    // padEff->cd();
    // padEff->SetTicks(0, 0);
    eff->Draw("P");
    // padEff->SetTicks(0, 0);
    // padEff->Update();

    // canv->cd();
    // padDist->cd();
    // padDist->SetTicks(0, 0);
    // Double_t weight = 0.0;
    // if (string(prodMode) == "WH") weight = 0.01507435332;
    // if (string(prodMode) == "ggH") weight = 0.14914894205;
    // denominator->Scale(weight);
    // denominator->Draw("histsame");
    // padDist->SetTicks(0, 0);
    
    file_.Close();
  }

  // writing the lumi information and the CMS "logo"
  CMS_lumi( canv, iPeriod, iPos );

  canv->Update();
  canv->RedrawAxis();
  canv->GetFrame()->Draw();

  canv->Print(canvName+".pdf",".pdf");
  canv->Print(canvName+".png",".png");

  return canv;
}

TCanvas* tau_eff_MC_plot( int iPeriod, int iPos, const char* var )
{ 
  int W = 800;
  int H = 600;

  int H_ref = 600; 
  int W_ref = 800; 

  // references for T, B, L, R
  float T = 0.08*H_ref;
  float B = 0.12*H_ref; 
  float L = 0.12*W_ref;
  float R = 0.04*W_ref;

  TString canvName = "tau_eff_";
  canvName += var;

  TCanvas* canv = new TCanvas(canvName,canvName,50,50,W,H);
  canv->SetFillColor(0);
  canv->SetBorderMode(0);
  canv->SetFrameFillStyle(0);
  canv->SetFrameBorderMode(0);
  canv->SetLeftMargin( L/W );
  canv->SetRightMargin( R/W );
  canv->SetTopMargin( T/H );
  canv->SetBottomMargin( B/H );
  canv->SetTickx(0);
  canv->SetTicky(0);
  
  TH1* hEff = new TH1F("hEff","hEff",40,0.0, 100.0);
  hEff->GetXaxis()->SetNdivisions(6,5,0);
  if (string(var) == "PT") hEff->GetXaxis()->SetTitle("#tau_{X} p_{T} (GeV)");
  if (string(var) == "Eta") hEff->GetXaxis()->SetTitle("#tau_{X} #eta");
  hEff->GetYaxis()->SetNdivisions(6,5,0);
  hEff->GetYaxis()->SetTitleOffset(1);
  hEff->GetYaxis()->SetTitle("#tau_{X} efficiency");  

  hEff->SetMaximum( 260 );
  if( iPos==1 ) hEff->SetMaximum( 300 );
  hEff->GetYaxis()->SetRangeUser(0.0, 1.0);
  hEff->Draw();

  int histLineColor = kBlue - 10;
  int histFillColor = kBlue - 10;
  float markerSize  = 1.0;

  {
    TLatex latex;
				
    int nLeft_ = 2;
    int nRight_ = 8;

    float x1_l = 0.65;
    float y1_l = 0.4;

    float dx_l = 0.30;
    float dy_l = 0.18;
    float x0_l = x1_l-dx_l;
    float x0_l_legendLeft = x1_l-2.0*dx_l;
    float y0_l = y1_l-dy_l;

    float ar_l = dy_l/dx_l;
		
    float x_l[1];
    float ex_l[1];
    float y_l_left[1];
    float ey_l_left[1];
    float y_l_right[1];
    float ey_l[1];
		
    float gapLeft_ = 1./(nLeft_+1);
    float gapRight_ = 1./(nRight_+1);
		
    float bwx_ = 0.12;
    float bwyLeft_ = gapLeft_/1.5;
    float bwyRight_ = gapRight_/1.5;
		
    x_l[0] = 1.2*bwx_;
    y_l_left[0] = 1-gapLeft_;
    y_l_right[0] = 1-gapRight_;
    ex_l[0] = 0;
    ey_l[0] = 0.04/ar_l;

    float xx_ = x_l[0];
    float yyLeft_ = y_l_left[0];

    TPad* legendRight = new TPad("legendRight_0","legendRight_0",x0_l,y0_l,x1_l, y1_l );
    legendRight->Draw();
    legendRight->cd();

    TGraph* gr_l = new TGraphErrors(1, x_l, y_l_left, ex_l, ey_l );
		
    gStyle->SetEndErrorSize(0);
    gr_l->SetMarkerSize(0.9);
    gr_l->SetMarkerColor(kRed);
    gr_l->SetMarkerStyle(23);
    gr_l->Draw("0P");
		
    latex.SetTextFont(42);
    latex.SetTextAngle(0);
    latex.SetTextColor(kBlack);    
    latex.SetTextSize(0.09375);    
    latex.SetTextAlign(12); 
		
    latex.DrawLatex(xx_+1.*bwx_,yyLeft_,"Wh, h#rightarrowaa#rightarrow4#tau");

    canv->cd();
  }

  {
    string fileName("WH_vs_ZTauTau_tauX_eff.root");
    TFile file_(fileName.c_str(),"READ");

    TCanvas* canvas = NULL;
    string canvasName("eff_visible");
    canvasName = canvasName + var;
    canvasName = canvasName + string("numerator_over_visible");
    canvasName = canvasName + var;
    canvasName = canvasName + string("denominator");
    file_.GetObject(canvasName.c_str(), canvas);

    if (canvas != NULL) {

      TList* primitives = canvas->GetListOfPrimitives();

      for (unsigned int i = 0; i < 2; ++i) {

	TGraphAsymmErrors *eff = (TGraphAsymmErrors*)primitives->At(i);
    
	// for (Int_t i = 0; i < eff->GetN(); ++i) {
	//   eff->SetPointEXhigh(i, 0.0);
	//   eff->SetPointEXlow(i, 0.0);
	// }
    
	if (i = 0) eff->Draw("P");
	else eff->Draw("PSAME");

      }
      
    }
    
    file_.Close();
  }

  // writing the lumi information and the CMS "logo"
  CMS_lumi( canv, iPeriod, iPos );

  canv->Update();
  canv->RedrawAxis();
  canv->GetFrame()->Draw();

  canv->Print(canvName+".pdf",".pdf");
  canv->Print(canvName+".png",".png");

  return canv;
}
