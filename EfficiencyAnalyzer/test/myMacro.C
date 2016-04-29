#include "tdrstyle.C"
#include "CMS_lumi.C"
#include "TH1.h"
#include "TH1F.h"
 
TCanvas* muHadMass_plot( int iPeriod, int iPos, const char* MTBin );
TCanvas* limit_plot( int iPeriod, int iPos, const char* MTBin );
 
void myMacro(const char* MTBin)
{
  gROOT->LoadMacro("tdrstyle.C");
  setTDRStyle();

  gROOT->LoadMacro("CMS_lumi.C");

  writeExtraText = true;       // if extra text
  extraText  = "Preliminary";  // default extra text is "Preliminary"
  // lumi_8TeV  = "19.1 fb^{-1}"; // default is "19.7 fb^{-1}"
  // lumi_7TeV  = "4.9 fb^{-1}";  // default is "5.1 fb^{-1}"
  // lumi_sqrtS = "13 TeV";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

  int iPeriod = 2;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 0=free form (uses lumi_sqrtS)

  // second parameter in muHadMass_plot is iPos, which drives the position of the CMS logo in the plot
  // iPos=11 : top-left, left-aligned
  // iPos=33 : top-right, right-aligned
  // iPos=22 : center, centered
  // mode generally : 
  //   iPos = 10*(alignement 1/2/3) + position (1/2/3 = left/center/right)

  // muHadMass_plot( iPeriod, 0 );   // out of frame (in exceptional cases)
  muHadMass_plot( iPeriod, 11, MTBin );  // left-aligned
  // limit_plot( iPeriod, 11, MTBin );
  //  muHadMass_plot( iPeriod, 33 );  // right-aligned

  //  writeExtraText = false;       // remove Preliminary
  
  //  muHadMass_plot( iPeriod, 0 );   // out of frame (in exceptional cases)

  //  muHadMass_plot( iPeriod, 11 );  // default: left-aligned
  //  muHadMass_plot( iPeriod, 22 );  // centered
  //  muHadMass_plot( iPeriod, 33 );  // right-aligned  
}

TCanvas* muHadMass_plot( int iPeriod, int iPos, const char* MTBin )
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

  TString canvName = "muHadMassCanvas_final_a9_";
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
  h->GetYaxis()->SetTitle("Events");  

  h->SetMaximum( 260 );
  if( iPos==1 ) h->SetMaximum( 300 );
  h->Draw();
  h->GetYaxis()->SetRangeUser(0.1, 10000.0);

  int histLineColor = kBlue - 10;
  int histFillColor = kBlue - 10;
  float markerSize  = 1.0;

  {
    TLatex latex;
				
    int nLeft_ = /*2*/3;
    int nRight_ = /*2*/4;

    float x1_l = /*0.92*//*0.95*/1.0;
    float y1_l = /*0.60*//*0.9*/0.95/*0.85*/;

    float dx_l = /*0.30*/0.4;
    float dy_l = /*0.18*//*0.54*//*0.3*/0.35;
    float x0_l = x1_l-dx_l;
    float x0_l_legendLeft = x1_l-1.7*dx_l;
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
    ex_l[0] = 0.04/ar_l;
    ey_l[0] = 0.04/ar_l;

    float xx_ = x_l[0];
    float yyLeft_ = y_l_left[0];

    TPad* legendLeft =
      new TPad("legendLeft_0","legendLeft_0",/*x0_l*/x0_l_legendLeft,y0_l,/*x1_l*/x0_l+0.3*dx_l, y1_l );
    // legendLeft->SetFillColor( kGray );
    legendLeft->Draw();
    legendLeft->cd();

    TGraph* gr_l = new TGraphErrors(1, x_l, y_l_left, ex_l, ey_l );
		
    gStyle->SetEndErrorSize(0);
    gr_l->SetMarkerSize(0.9);
    gr_l->Draw("0P");
		
    latex.SetTextFont(42);
    latex.SetTextAngle(0);
    latex.SetTextColor(kBlack);    
    latex.SetTextSize(/*0.25*/0.09375);    
    latex.SetTextAlign(12); 
		
    TLine line_;
    TBox  box_;
    latex.DrawLatex(xx_+1.*bwx_,yyLeft_,"Data");
		
    yyLeft_ -= gapLeft_;
    box_.SetLineStyle( kSolid );
    box_.SetLineWidth( 1 );
    //		box_.SetLineColor( kBlack );
    box_.SetLineColor( histLineColor );
    box_.SetFillColor( histFillColor );
    box_.DrawBox( xx_-bwx_/2, yyLeft_-bwyLeft_/2, xx_+bwx_/2, yyLeft_+bwyLeft_/2 );
    box_.SetFillStyle(0);
    box_.DrawBox( xx_-bwx_/2, yyLeft_-bwyLeft_/2, xx_+bwx_/2, yyLeft_+bwyLeft_/2 );
    latex.DrawLatex(xx_+1.*bwx_,yyLeft_,"Jet fake bkg.");

    yyLeft_ -= gapLeft_;
    box_.SetLineStyle( kSolid );
    box_.SetLineWidth( 1 );
    //		box_.SetLineColor( kBlack );
    box_.SetLineColor( kRed );
    box_.SetFillColor( kRed );
    box_.DrawBox( xx_-bwx_/2, yyLeft_-bwyLeft_/2, xx_+bwx_/2, yyLeft_+bwyLeft_/2 );
    box_.SetFillStyle(3005);
    box_.DrawBox( xx_-bwx_/2, yyLeft_-bwyLeft_/2, xx_+bwx_/2, yyLeft_+bwyLeft_/2 );
    latex.DrawLatex(xx_+1.*bwx_,yyLeft_,"Bkg. syst. error");

    canv->cd();
	
    TPad* legendRight = new TPad("legendRight_0","legendRight_0",x0_l,y0_l,x1_l, y1_l );
    // legendRight->SetFillColor( kGray );
    legendRight->Draw();
    legendRight->cd();

    // xx_+=dx_l;
    float yyRight_ = y_l_right[0];
    // yyLeft_ -= gapLeft_;
    line_.SetLineStyle( kDashed );
    line_.SetLineWidth( 4 );
    //		line_.SetLineColor( kBlack );
    line_.SetLineColor( kGreen - 2 );
    line_.DrawLine( xx_-bwx_/2, yyRight_, xx_+bwx_/2, yyRight_ );
    latex.DrawLatex(xx_+1.*bwx_,yyRight_,"ggH m_{a} = 9 GeV");

    yyRight_ -= gapRight_;
    line_.SetLineStyle( kDashed );
    line_.SetLineWidth( 4 );
    //		line_.SetLineColor( kBlack );
    line_.SetLineColor( kBlue );
    line_.DrawLine( xx_-bwx_/2, yyRight_, xx_+bwx_/2, yyRight_ );
    latex.DrawLatex(xx_+1.*bwx_,yyRight_,"WH m_{a} = 9 GeV");

    yyRight_ -= gapRight_;
    line_.SetLineStyle( kDashed );
    line_.SetLineWidth( 4 );
    //		line_.SetLineColor( kBlack );
    line_.SetLineColor( kMagenta + 2 );
    line_.DrawLine( xx_-bwx_/2, yyRight_, xx_+bwx_/2, yyRight_ );
    latex.DrawLatex(xx_+1.*bwx_,yyRight_,"VBF m_{a} = 9 GeV");

    yyRight_ -= gapRight_;
    line_.SetLineStyle( kDashed );
    line_.SetLineWidth( 4 );
    //		line_.SetLineColor( kBlack );
    line_.SetLineColor( kOrange + 7 );
    line_.DrawLine( xx_-bwx_/2, yyRight_, xx_+bwx_/2, yyRight_ );
    latex.DrawLatex(xx_+1.*bwx_,yyRight_,"ZH m_{a} = 9 GeV");

    canv->cd();

    TLatex latexTitle;
    latexTitle.SetTextFont(42);
    latexTitle.SetTextAngle(0);
    latexTitle.SetTextColor(kBlack);    
    latexTitle.SetTextSize(0.0524476);    
    latexTitle.SetTextAlign(12);     
    Double_t userX = 0.5;
    Double_t userY = 497.869;
    string label;
    if (string(MTBin) == "lowMT") label = "Low M_{T}";
    if (string(MTBin) == "highMT") label = "High M_{T}";
    latexTitle.DrawLatex(userX,userY,label.c_str());

  }

  {
    // Observed data
    string fileName("../final_");
    fileName = fileName + string(MTBin);
    fileName = fileName + string("_a9_v214.root");
    TFile file_(fileName.c_str(),"READ");

    TH1F *data = static_cast<TH1F*>(file_.Get("muHadMassDataSearch")->Clone());
    data->SetDirectory(0);
    data->SetMarkerStyle(20);
    data->SetMarkerSize(markerSize);

    TH1F *bkg   = static_cast<TH1F*>(file_.Get("muHadMassAvg")->Clone());
    bkg->SetDirectory(0);
    bkg->SetLineColor(histLineColor);
    bkg->SetFillColor(histFillColor);
    bkg->SetFillStyle(1001);

    TH1F *systErrBand   = static_cast<TH1F*>(file_.Get("muHadMassSystErr")->Clone());
    systErrBand->SetDirectory(0);
    systErrBand->SetLineColor(kRed);
    systErrBand->SetFillColor(kRed);
    systErrBand->SetFillStyle(3005);

    TH1F *sigWH   = static_cast<TH1F*>(file_.Get("muHadMassWHSearch")->Clone());
    sigWH->SetDirectory(0);
    sigWH->SetLineColor(kBlue);
    sigWH->SetFillColor(kBlue);

    TH1F *sigGGH   = static_cast<TH1F*>(file_.Get("muHadMassggHSearch")->Clone());
    sigGGH->SetDirectory(0);
    sigGGH->SetLineColor(kGreen - 2);
    sigGGH->SetFillColor(kGreen - 2);

    TH1F *sigZH   = static_cast<TH1F*>(file_.Get("muHadMassZHSearch")->Clone());
    sigZH->SetDirectory(0);
    sigZH->SetLineColor(kOrange + 7);
    sigZH->SetFillColor(kOrange + 7);

    TH1F *sigVBF   = static_cast<TH1F*>(file_.Get("muHadMassVBFSearch")->Clone());
    sigVBF->SetDirectory(0);
    sigVBF->SetLineColor(kMagenta + 2);
    sigVBF->SetFillColor(kMagenta + 2);

    bkg->Draw("histsame");
    systErrBand->Draw("e2same");
    sigWH->Draw("histsame");
    sigGGH->Draw("histsame");
    sigZH->Draw("histsame");
    sigVBF->Draw("histsame");
    data->Draw("esame");

    TLine SRBound;
    SRBound.SetLineStyle( kSolid );
    SRBound.SetLineWidth( 2 );
    SRBound.SetLineColor( kBlack );
    SRBound.DrawLine( 4.0, 60.0, 4.0, 200.0 );

    TArrow SRArrow(4.0, 100.0, 5.0, 100.0, 0.025, "|>");
    SRArrow.SetLineStyle( kSolid );
    SRArrow.SetLineWidth( 2 );
    SRArrow.SetLineColor( kBlack );
    SRArrow.Draw();

    TLatex SRLabel;
    SRLabel.SetTextFont(42);
    SRLabel.SetTextAngle(0);
    SRLabel.SetTextColor(kBlack);    
    SRLabel.SetTextSize(/*0.025*/0.03);    
    SRLabel.SetTextAlign(12);     
    Double_t userX = 5.1;
    Double_t userY = 90.0;
    SRLabel.DrawLatex(userX,userY,"Counting experiment signal window m_{#mu+X} #geq 4 GeV");

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

TCanvas* limit_plot( int iPeriod, int iPos, const char* MTBin )
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

  TString canvName = "expLimits_Br_";
  string MTBinSuffix;
  if (string(MTBin) == "lowMT") MTBinSuffix = "ggHVBF";
  if (string(MTBin) == "highMT") MTBinSuffix = "ggHWH";
  if (string(MTBin) == "combined") canvName = "expLimits_Br_20GeV_4sig";
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
  else {
    canvName += MTBin;
    canvName += "_20GeV_" + MTBinSuffix;
  }

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

  TH1* h = new TH1F("h","h",40,4.0, 16.0);
  h->GetXaxis()->SetNdivisions(6,5,0);
  h->GetXaxis()->SetTitle("m_{a} (GeV)");  
  h->GetYaxis()->SetNdivisions(6,5,0);
  h->GetYaxis()->SetTitleOffset(1);
  h->GetYaxis()->SetTitle("95% CL limit on BR(H#rightarrowaa#rightarrow4#tau)");  

  h->SetMaximum( 260 );
  if( iPos==1 ) h->SetMaximum( 300 );
  h->Draw();
  h->GetYaxis()->SetRangeUser(0.05, 10000.0);

  int histLineColor = kBlue - 10;
  int histFillColor = kBlue - 10;
  float markerSize  = 1.0;

  {
    TLatex latex;
				
    int nLeft_ = /*2*/3;
    int nRight_ = /*2*/4;

    float x1_l = /*0.92*/0.95;
    float y1_l = /*0.60*//*0.9*//*0.85*/0.9175;

    float dx_l = 0.30;
    float dy_l = /*0.18*/0.54/*0.3*/;
    float x0_l = x1_l-dx_l;
    float x0_l_legendLeft = x1_l-1.7*dx_l;
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

    // TPad* legendLeft =
    //   new TPad("legendLeft_0","legendLeft_0",/*x0_l*/x0_l_legendLeft,y0_l,/*x1_l*/x0_l+0.3*dx_l, y1_l );
    // // legendLeft->SetFillColor( kGray );
    // legendLeft->Draw();
    // legendLeft->cd();

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
		
    // yyLeft_ -= gapLeft_;
    // box_.SetLineStyle( kSolid );
    // box_.SetLineWidth( 1 );
    // //		box_.SetLineColor( kBlack );
    // box_.SetLineColor( histLineColor );
    // box_.SetFillColor( histFillColor );
    // box_.DrawBox( xx_-bwx_/2, yyLeft_-bwyLeft_/2, xx_+bwx_/2, yyLeft_+bwyLeft_/2 );
    // box_.SetFillStyle(0);
    // box_.DrawBox( xx_-bwx_/2, yyLeft_-bwyLeft_/2, xx_+bwx_/2, yyLeft_+bwyLeft_/2 );
    // latex.DrawLatex(xx_+1.*bwx_,yyLeft_,"Jet fake bkg.");

    // yyLeft_ -= gapLeft_;
    // box_.SetLineStyle( kSolid );
    // box_.SetLineWidth( 1 );
    // //		box_.SetLineColor( kBlack );
    // box_.SetLineColor( kRed );
    // box_.SetFillColor( kRed );
    // box_.DrawBox( xx_-bwx_/2, yyLeft_-bwyLeft_/2, xx_+bwx_/2, yyLeft_+bwyLeft_/2 );
    // box_.SetFillStyle(3005);
    // box_.DrawBox( xx_-bwx_/2, yyLeft_-bwyLeft_/2, xx_+bwx_/2, yyLeft_+bwyLeft_/2 );
    // latex.DrawLatex(xx_+1.*bwx_,yyLeft_,"Bkg. syst. error");

    // canv->cd();
	
    TPad* legendRight = new TPad("legendRight_0","legendRight_0",x0_l,y0_l,x1_l, y1_l );
    // legendRight->SetFillColor( kGray );
    legendRight->Draw();
    legendRight->cd();

    float yyRight_ = y_l_right[0];
    line_.SetLineStyle( kSolid );
    line_.SetLineWidth( 2 );
    line_.SetLineColor( kBlack );
    line_.DrawLine( xx_-bwx_/2, yyRight_, xx_+bwx_/2, yyRight_ );
    latex.DrawLatex(xx_+1.*bwx_,yyRight_,"Observed");

    yyRight_ -= gapRight_;
    line_.SetLineStyle( kDashed );
    line_.SetLineWidth( 2 );
    line_.SetLineColor( kBlack );
    line_.DrawLine( xx_-bwx_/2, yyRight_, xx_+bwx_/2, yyRight_ );
    latex.DrawLatex(xx_+1.*bwx_,yyRight_,"Expected");

    yyRight_ -= gapRight_;
    box_.SetLineStyle( kSolid );
    box_.SetLineWidth( 0 );
    box_.SetFillColor( 3 );
    box_.DrawBox( xx_-bwx_/2, yyRight_-bwyRight_/2, xx_+bwx_/2, yyRight_+bwyRight_/2 );
    box_.SetFillStyle(1001);
    box_.DrawBox( xx_-bwx_/2, yyRight_-bwyRight_/2, xx_+bwx_/2, yyRight_+bwyRight_/2 );
    latex.DrawLatex(xx_+1.*bwx_,yyRight_,"#pm1 #sigma Exp.");

    yyRight_ -= gapRight_;
    box_.SetLineStyle( kSolid );
    box_.SetLineWidth( 0 );
    box_.SetFillColor( 5 );
    box_.DrawBox( xx_-bwx_/2, yyRight_-bwyRight_/2, xx_+bwx_/2, yyRight_+bwyRight_/2 );
    box_.SetFillStyle(1001);
    box_.DrawBox( xx_-bwx_/2, yyRight_-bwyRight_/2, xx_+bwx_/2, yyRight_+bwyRight_/2 );
    latex.DrawLatex(xx_+1.*bwx_,yyRight_,"#pm2 #sigma Exp.");

    canv->cd();
    
    TLatex latexTitle;
    latexTitle.SetTextFont(42);
    latexTitle.SetTextAngle(0);
    latexTitle.SetTextColor(kBlack);    
    latexTitle.SetTextSize(0.0524476);    
    latexTitle.SetTextAlign(12);
    Double_t userX = 0.0;
    Double_t userY = 0.0;
    if (string(MTBin) == "combined") {
      userX = 5.0;
      userY = 80.3009;
    }
    else {
      userX = 9.26705;
      userY = 80.3009;
    }
    string label;
    if (string(MTBin) == "lowMT") label = "Low M_{T}";
    if (string(MTBin) == "highMT") label = "High M_{T}";
    if (string(MTBin) == "combined") label = "Low and high M_{T} combination";
    latexTitle.DrawLatex(userX,userY,label.c_str());
  }

  {
    // Observed data
    string fileName("comb_");
    fileName = fileName + string(MTBin);
    fileName = fileName + string("_limit.root");
    TFile file_(fileName.c_str(),"READ");

    TGraphAsymmErrors *exp2SigmaBand   =
      static_cast<TGraphAsymmErrors*>(file_.Get("band_exp2")->Clone());
    // exp2SigmaBand->SetDirectory(0);
    exp2SigmaBand->SetLineColor(0);
    exp2SigmaBand->SetFillColor(5);
    exp2SigmaBand->SetFillStyle(1001);

    TGraphAsymmErrors *exp1SigmaBand   =
      static_cast<TGraphAsymmErrors*>(file_.Get("band_exp1")->Clone());
    // exp1SigmaBand->SetDirectory(0);
    exp1SigmaBand->SetLineColor(0);
    exp1SigmaBand->SetFillColor(3);
    exp1SigmaBand->SetFillStyle(1001);

    TGraph *obsLimit = static_cast<TGraph*>(file_.Get("limit_obs")->Clone());
    // obsLimit->SetDirectory(0);
    obsLimit->SetLineWidth(2);

    TGraph *expLimit = static_cast<TGraph*>(file_.Get("limit_exp")->Clone());
    // expLimit->SetDirectory(0);
    expLimit->SetLineStyle(2);
    
    // expLimit->Draw("AL");
    exp2SigmaBand->Draw("3same");
    exp1SigmaBand->Draw("3same");
    obsLimit->Draw("Lsame");
    expLimit->Draw("Lsame");
    
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
