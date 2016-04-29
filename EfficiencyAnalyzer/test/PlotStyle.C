#include <string>
#include "TH1F.h"
#include "TH2F.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLegend.h"

//default drawing options
const Double_t defaultXAxisLabelSize = 0.05;
const Int_t defaultCanvasWidth = 600;
const Int_t defaultCanvasHeight = 600;

//set pad drawing options
void setCanvasOptions(TVirtualPad& canvas, const Int_t grid, const Int_t logY, const Int_t logZ)
{
  canvas.SetFillStyle(0);
  canvas.SetFillColor(0);
  canvas.SetGrid(grid, grid);
  canvas.SetLogy(logY);
  canvas.SetLogz(logZ);
}

//set canvas drawing options
void setCanvasOptions(TCanvas& canvas, const Int_t grid, const Int_t logY, const Int_t logZ)
{
  canvas.SetFillStyle(0);
  canvas.SetFillColor(0);
  canvas.SetGrid(grid, grid);
  canvas.SetLogy(logY);
  canvas.SetLogz(logZ);
}

//set pad margins
void setCanvasMargins(TVirtualPad& canvas, const float left, const float top, const float right, 
		      const float bottom)
{
  canvas.cd()->SetLeftMargin(left);
  canvas.cd()->SetTopMargin(top);
  canvas.cd()->SetRightMargin(right);
  canvas.cd()->SetBottomMargin(bottom);
}

//set canvas margins
void setCanvasMargins(TCanvas& canvas, const float left, const float top, const float right, 
		      const float bottom)
{
  canvas.cd()->SetLeftMargin(left);
  canvas.cd()->SetTopMargin(top);
  canvas.cd()->SetRightMargin(right);
  canvas.cd()->SetBottomMargin(bottom);
}

//set axis drawing options
void setAxisOptions(TAxis* axis, const Float_t labelSize, const Float_t titleOffset, 
		    const char* title)
{
  axis->SetLabelFont(42);
  axis->SetLabelOffset(0.007);
  axis->SetLabelSize(labelSize);
  axis->SetTitleFont(42);
  axis->SetTitleSize(0.06);
  axis->SetTitleOffset(titleOffset);
  axis->SetTitle(title);
}

//set axis drawing options incl. title
void setAxisOptions(TAxis* axis, const Float_t labelSize, const Float_t titleSize, 
		    const Float_t titleOffset, const char* title)
{
  axis->SetLabelFont(42);
  axis->SetLabelOffset(0.007);
  axis->SetLabelSize(labelSize);
  axis->SetTitleFont(42);
  axis->SetTitleSize(titleSize);
  axis->SetTitleOffset(titleOffset);
  axis->SetTitle(title);
}

//set axis drawing options excl. title
void setAxisOptions(TAxis* axis, const Float_t labelSize, const Float_t titleSize, 
		    const Float_t titleOffset)
{
  axis->SetLabelFont(42);
  axis->SetLabelOffset(0.007);
  axis->SetLabelSize(labelSize);
  axis->SetTitleFont(42);
  axis->SetTitleSize(titleSize);
  axis->SetTitleOffset(titleOffset);
}

//set axis labels
void setAxisLabels(TAxis* axis, const vector<string>& binLabels)
{
  for (Int_t iBin = 1; iBin <= axis->GetNbins(); ++iBin) {
    if (iBin <= (int)binLabels.size()) axis->SetBinLabel(iBin, binLabels[iBin - 1].c_str());
  }
}

//set graph drawing options
void setGraphOptions(TGraphAsymmErrors& graph, const Color_t color, const Size_t size, 
		     const Style_t style, const char* xAxisTitle, const char* yAxisTitle)
{
  graph.SetMarkerColor(color);
  graph.SetMarkerSize(size);
  graph.SetMarkerStyle(style);
  graph.SetLineColor(color);
  graph.SetTitle("");
  setAxisOptions(graph.GetXaxis(), defaultXAxisLabelSize, 0.9, xAxisTitle);
  setAxisOptions(graph.GetYaxis(), defaultXAxisLabelSize, 1.05, yAxisTitle);
}

//set 1D histogram drawing options
void setHistogramOptions(TH1F* histogram, const Color_t color, const Size_t size, 
			 const Style_t style, const Double_t scale, const char* xAxisTitle, 
			 const char* yAxisTitle, 
			 const Double_t xAxisLabelSize = defaultXAxisLabelSize)
{
  histogram->SetMarkerColor(color);
  histogram->SetMarkerSize(size);
  histogram->SetMarkerStyle(style);
  histogram->SetLineColor(color);
  histogram->SetLineWidth(2);
  histogram->SetFillStyle(0);
  setAxisOptions(histogram->GetXaxis(), xAxisLabelSize, 0.9, xAxisTitle);
  setAxisOptions(histogram->GetYaxis(), defaultXAxisLabelSize, 1.05, yAxisTitle);
  histogram->Scale(scale);
}

//set 2D histogram options incl. y-axis title
void setHistogramOptions(TH2F* histogram, const Color_t color, const Size_t size, 
			 const Style_t style, const Float_t yAxisTitleOffset, 
			 const Double_t scale, const char* xAxisTitle, const char* yAxisTitle)
{
  histogram->SetMarkerColor(color);
  histogram->SetMarkerSize(size);
  histogram->SetMarkerStyle(style);
  histogram->SetLineColor(color);
  histogram->SetLineWidth(1);
  histogram->SetFillStyle(0);
  setAxisOptions(histogram->GetXaxis(), 0.04, 0.04, 1.1, xAxisTitle);
  setAxisOptions(histogram->GetYaxis(), 0.04, 0.04, yAxisTitleOffset, yAxisTitle);
  setAxisOptions(histogram->GetZaxis(), 0.04, 0.04, histogram->GetZaxis()->GetTitleOffset(), "");
  histogram->Scale(scale);
}

//set 2D histogram options excl. y-axis title
void setHistogramOptions(TH2F* histogram, const Color_t color, const Size_t size, 
			 const Style_t style, const Float_t yAxisTitleOffset, const Double_t scale)
{
  histogram->SetMarkerColor(color);
  histogram->SetMarkerSize(size);
  histogram->SetMarkerStyle(style);
  histogram->SetLineColor(color);
  histogram->SetLineWidth(1);
  histogram->SetFillStyle(0);
  setAxisOptions(histogram->GetXaxis(), 0.04, 0.04, 1.1);
  setAxisOptions(histogram->GetYaxis(), 0.04, 0.04, yAxisTitleOffset);
  setAxisOptions(histogram->GetZaxis(), 0.04, 0.04, histogram->GetZaxis()->GetTitleOffset(), "");
  histogram->Scale(scale);
}

//set legend drawing options
void setLegendOptions(TLegend& legend, const char* header)
{
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.SetTextFont(42);
  legend.SetHeader(header);
}
