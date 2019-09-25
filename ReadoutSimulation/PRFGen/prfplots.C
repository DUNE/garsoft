#include<TH2F.h>
#include<TCanvas.h>
#include<TMath.h>
#include<TAxis.h>
#include<TStyle.h>
#include<TROOT.h>
#include<TString.h>
#include<TFile.h>

void plotrespfunc();
void plotrocpad(TString rocname="iroc", float padwidth=0.4, float padlength=7.5, float xpwid=1.0, float ypwid=1.0,
		int nbx=50, int nby=50, float wireoffset = 0);

float rfunc(float x, float y);
float rfuncroc(float x, float y, float padwidth, float padlength, float wireoffset);
float rocint(float x, float y, float padwidth, float padlength);

void norm2f(TH2F *h);

void prfplots(TString rootfile="mpdtpcprf_v1.root")
{
  TFile outrootfile(rootfile,"RECREATE");
  gStyle->SetOptStat(0);
  plotrespfunc();
  plotrocpad("IROC", 0.4,0.75,1.0, 1.0, 50,50, 0);
  plotrocpad("IOROC",0.6,1.0, 1.0, 1.0, 50,50, 0.5);
  plotrocpad("OOROC",0.6,1.5, 1.0, 1.0, 50,50, 0.5);
  outrootfile.Close();
}

// normalize to unit area a la fig 7.7 in the ALICE TDR
void norm2f(TH2F *h)
{
  int nbinsx = h->GetNbinsX();
  int nbinsy = h->GetNbinsY();
  float xmin = h->GetXaxis()->GetBinLowEdge(1);
  float xmax = h->GetXaxis()->GetBinUpEdge(nbinsx);
  float ymin = h->GetYaxis()->GetBinLowEdge(1);
  float ymax = h->GetYaxis()->GetBinUpEdge(nbinsy);
  float integral = h->Integral();
  float fact = 1.0/((ymax-ymin)*(xmax-xmin)*integral/(nbinsx*nbinsy));
  fact /= 4;  // not sure -- maybe four pads total
  h->Scale(fact);
}

void plotrespfunc()
{
  TCanvas *mycanvas = (TCanvas*) new TCanvas("rsc","Response Function",700,460);
  TH2F *respfunc = (TH2F*) new TH2F("crf",";pad direction (cm);pad row direction (cm)",30,-1,1,30,-1,1);
  for (int ibinx=1; ibinx <= respfunc->GetNbinsX(); ++ibinx)
    {
      float x = respfunc->GetXaxis()->GetBinCenter(ibinx);
      for (int ibiny=1; ibiny <= respfunc->GetNbinsY(); ++ibiny)
        {
	  float y = respfunc->GetYaxis()->GetBinCenter(ibiny);
	  respfunc->SetBinContent(ibinx,ibiny,rfunc(x,y));
	}
    }
  int ndiv = -404;
  respfunc->GetXaxis()->SetNdivisions(ndiv);
  respfunc->GetYaxis()->SetNdivisions(ndiv);
  respfunc->SetLineColor(1);
  respfunc->SetLineWidth(2);
  respfunc->Draw("lego");
  respfunc->SetDirectory(0);
  respfunc->Write();
  mycanvas->Print("response_function.pdf");
}

float rfunc(float x, float y)
{
  float resp = 0.07*TMath::Exp( -(x*x + y*y)/0.03 )
     + 0.01*TMath::Exp( -(x*x + y*y)/0.08 );
  return resp;
}


void plotrocpad(TString rocname, float padwidth, float padlength, float xpwid, float ypwid, int nbx, int nby,
                float wireoffset)
{
  TString cname="canvas";
  cname += rocname;
  TString hname="resp";
  hname += rocname;
  TString ctitle=rocname;
  ctitle += " Response Function";
  TString htitle=rocname;
  htitle += " Pad Response;pad direction (cm);pad row direction (cm)";
  TCanvas *mycanvas = (TCanvas*) new TCanvas(cname,ctitle,700,460);
  TH2F *roc = (TH2F*) new TH2F(hname,htitle,nbx,-xpwid,xpwid,nby,-ypwid,ypwid);
  for (int ibinx=1; ibinx <= roc->GetNbinsX(); ++ibinx)
    {
      float x = roc->GetXaxis()->GetBinCenter(ibinx);
      for (int ibiny=1; ibiny <= roc->GetNbinsY(); ++ibiny)
        {
	  float y = roc->GetYaxis()->GetBinCenter(ibiny);
	  roc->SetBinContent(ibinx,ibiny,rfuncroc(x,y,padwidth,padlength,wireoffset));
	}
    }
  int ndiv = -404;
  roc->GetXaxis()->SetNdivisions(ndiv);
  roc->GetYaxis()->SetNdivisions(ndiv);
  roc->SetLineColor(1);
  roc->SetLineWidth(2);
  norm2f(roc);
  roc->Draw("lego");  
  roc->SetDirectory(0);
  roc->Write();

  TString outfilename = rocname;
  outfilename += "_pad_response.pdf";
  mycanvas->Print(outfilename);

}

float rfuncroc(float x, float y, float padwidth, float padlength, float wireoffset)
{
  float wpitch = 0.25;  // wire pitch in cm
  int iwire = nearbyint( (y/wpitch) - wireoffset );
  float tmpy = (iwire+wireoffset)*wpitch;
  return rocint(x, tmpy, padwidth, padlength);
}


// numerically integrate the pad response over the pad
float rocint(float x, float y, float padwidth, float padlength)
{
  float xlim = padwidth/2.0;
  float ylim  = padlength/2.0;
  float dx = 0.01;
  float dy = 0.01;

  double sum = 0;
  for (float xprime = -xlim; xprime<=xlim; xprime += dx)
    {
      for (float yprime = -ylim; yprime < ylim; yprime += dy)
	{
	  sum += rfunc(x-xprime,y-yprime);
	}
    }
  return sum*dx*dy;
}

