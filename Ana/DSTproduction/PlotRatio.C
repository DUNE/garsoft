// Data/MC ratio plotter.  Histograms must be of class TH1D; the histogram
// RAT will exist after execution of doRat or doEff.  These methods do the
// ratio with uncertainties either by sum-in-quadrature based on event weights
// or by computing the relative uncertainty from the number of effective events
// with the binomial formula and scaling that to the ratio from sum-of-weights.
// Entries with zero or negative uncertainties are incorporated anyway in the
// fit, starting in v5 of ROOT, so this macro will will not give the correct
// fit in that situation, and there isn't a simple fix.  Because ROOT is so
// great.



#include <iostream>
#include <fstream>
#include <sstream>
#include <deque>

#include "Rtypes.h"
#include "TMath.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TStyle.h"

using namespace std;


//---   class definition here   ----------------------------------------
class PlotRatio {
public:



	PlotRatio() {
		ratioCanvas = new TCanvas("ratioCanvas","ratioCanvas",600,900);
		minYaxis = maxYaxis = 0;
		Slope = false;
	};
	PlotRatio(const char* canvasname, int Xwide=600,int Yhigh=900, int Xpos=50,int Ypos=50) {
		ratioCanvas = new TCanvas(canvasname,canvasname,Xpos,Ypos,Xwide,Yhigh);
		minYaxis = maxYaxis = 0;
		Slope = false;
	};
	~PlotRatio() {
		delete ratioCanvas;
	}

	int doRat(char* histNum, char* histDen, bool Log=true, bool Scale=false) {
		TH1D* hNum = (TH1D*)gDirectory->Get(histNum);
		TH1D* hDen = (TH1D*)gDirectory->Get(histDen);
		return doFullRatio(hNum,hDen,false,Log,Scale);
	}
	int doEff(char* histNum, char* histDen, bool Log=true, bool Scale=false) {
		TH1D* hNum = (TH1D*)gDirectory->Get(histNum);
		TH1D* hDen = (TH1D*)gDirectory->Get(histDen);
		return doFullRatio(hNum,hDen,true,Log,Scale);
	}

	int doRat(TH1D* hNum, TH1D* hDen, bool Log=true, bool Scale=false) {
		return doFullRatio(hNum,hDen,false,Log,Scale);
	}
	int doEff(TH1D* hNum, TH1D* hDen, bool Log=true, bool Scale=false) {
		return doFullRatio(hNum,hDen,true,Log,Scale);
	}

	double minYaxis,maxYaxis;
	bool Slope;


	
private:
	TCanvas* ratioCanvas;



	int doFullRatio(TH1D* hNum, TH1D* hDen,
					bool binomial, bool Log, bool Scale) {
		Int_t	  Nbin[2];
		Double_t  lowX[2],highX[2];

		// Check the binning
		Nbin[0] = hNum->GetNbinsX();
		Nbin[1] = hDen->GetNbinsX();
		if (Nbin[0] != Nbin[1]) {
			cout << "PlotRatio: Binning mismatch (Nbin)" << endl;
			return 1;
		}
		lowX[0]  = hNum->GetBinLowEdge(1);
		lowX[1]  = hDen->GetBinLowEdge(1);
		if ( lowX[0] != lowX[1] ) {
			cout << "PlotRatio: Binning mismatch (low edge)" << endl;
			return 1;
		}

		highX[0] = hNum->GetBinLowEdge(Nbin[0]);
		highX[1] = hDen->GetBinLowEdge(Nbin[0]);
		if ( highX[0] != highX[1] ) {
			cout << "PlotRatio: Binning mismatch (high edge)" << endl;
			return 1;
		}

		// I think one needs to use Sumw2() in all scenarios
		if (hNum->GetSumw2N()==0) {
			cout << "Numerator never heard you say Sumw2()\n";
			return 1;
		}
		if (hDen->GetSumw2N()==0) {
			cout << "Denominator never heard you say Sumw2()\n";
			return 1;
		}

		// Create RAT for ratio
		TH1D* RAT = (TH1D*)hNum->Clone();
		RAT->SetName("RAT");
		RAT->Divide(hDen);

		// Create SCL = scaled denominator for plotting, Chi2 only
		TH1D* SCL = (TH1D*)hDen->Clone();		SCL->Sumw2();
		Double_t ScaleFactor = 1.0;
		if (Scale) {
			ScaleFactor = hNum->Integral(1,Nbin[0])
						 /hDen->Integral(1,Nbin[0]);
			// cout << "Scaling MC by " << ScaleFactor;
		}
		SCL->Scale(ScaleFactor);


		// Find the correct uncertainties.  Default is sum-in-quad of relative
		// uncertainties; optionally, do a binomial ratio
		if (binomial) {
			// Not right for weighted histograms
			Double_t num, den, err;
			int badN=0, badD=0;
			for (Int_t i=0; i<=Nbin[0]+1; ++i) {
				num   = hNum->GetBinContent(i);
				den   = hDen->GetBinContent(i);

				if ( num > den ) {
					RAT->SetBinContent(i, 0.0);
					RAT->SetBinError  (i, 0.0);
					if ( badN==0 ) {
						cout << "N(numerator) > N(denominator) ";
						cout << "starting at bin " << i << endl;
						badN = i;
					}
				}
				if ( den <= 0.0 ) {
					RAT->SetBinContent(i, 0.0);
					RAT->SetBinError(i, 0.0);
					if ( i!=0 && i!=Nbin[0] && badD==0 ) {
						cout << "N(denominator) <= 0 ";
						cout << "starting at bin " << i << endl;
						badD = i;
					}
				}
				if ( num<=den && den>0 ) {
					err = num*(1 - num/den);
					err = sqrt(err) / den;
					RAT->SetBinError(i,err);
				}
			}
		}


		// Plot the dang thing.
		ratioCanvas->Divide(1,2);
		ratioCanvas->cd(1);				// Plot the MC over the data on top
		Double_t macks,mackSave;

		if ( hNum->GetMaximum() > SCL->GetMaximum() ){
			macks = hNum->GetMaximum();
		} else {
			macks = SCL->GetMaximum();
		}
		if (maxYaxis!=0) macks = maxYaxis;
		if (minYaxis!=0) cout << "Yah, gotta implement minYaxis someday" << endl;
		if (macks > 0.0 && Log) {
			gPad->SetLogy(1);
		} else {
			gPad->SetLogy(0);
		}
		mackSave = hNum->GetMaximum();
		hNum->SetMaximum(1.2*macks);
		hNum->SetMarkerStyle(20);
		hNum->SetMarkerSize(0.8);
		hNum->SetMarkerColor(kBlack);
		hNum->SetLineColor(kBlack);
		hNum->DrawCopy("PE");
		//gStyle->SetOptStat("ourMi");
		//gPad->Update();
		SCL->SetLineWidth(2);
		SCL->SetLineColor(kRed);
		SCL->DrawCopy("SAME HIST");

		ratioCanvas->cd(2);				// Fit the ratio on the bottom
		hNum->SetMaximum(mackSave);
		gPad->SetLogy(0);
		RAT->SetMarkerStyle(20);
		RAT->SetMarkerSize(0.8);
		RAT->SetLineWidth(1);
		RAT->SetMarkerColor(kBlack);
		// TH1.Divide(TH1) screws up max & min; this fixes it.
		// Because ROOT is so great.
		RAT->SetMaximum();	RAT->SetMinimum();
		RAT->GetMaximum();	RAT->GetMinimum();
		RAT->SetMaximum();	RAT->SetMinimum();
		RAT->GetMaximum();	RAT->GetMinimum();
		if (Slope) {
			RAT->Fit("pol1","QEM","PE9");
		} else {
			RAT->Fit("pol0","QEM","PE9");
		}
		gStyle->SetOptStat("");
		gStyle->SetOptFit(1011);
		gPad->Update();


		Double_t Chi2,pull,err2;	// Print quality of fit message
		Chi2 = 0.0;
		Int_t NDoF = -2;			// 'cuz 2 DoF in line Fit
		if (Scale) ++NDoF;
		Int_t OVR  = Nbin[0]+1;
		for (Int_t i=1; i<=Nbin[0]; ++i) {
			if (hNum->GetBinContent(i)!=0.0 || SCL->GetBinContent(i)!=0.0) {
				pull =  hNum->GetBinContent(i) -SCL->GetBinContent(i);
				err2 =  hNum->GetBinError(i)*hNum->GetBinError(i)
					   +SCL->GetBinError(i)*SCL->GetBinError(i);
				//cout << "Bin " << i << " Pull " << pull/sqrt(err2) << endl;
				if (err2<=0.0) {
					cout << "Zero error, non-zero content at bin "
						 << i << endl;
					err2 = hNum->GetBinContent(i) +SCL->GetBinContent(i);			
				}
				Chi2 += pull*pull/err2;
				++NDoF;
			}
		}
		cout << "\nPlotRatio: " << hNum->GetName() << " / " << hDen->GetName() << "\n";
		cout << "Chi2/NDoF: " << Chi2 << " / " << NDoF << "   C.L.: " 
			 << TMath::Prob(Chi2,NDoF) << endl;
		if (hNum->GetBinContent(0)!=0.0 || SCL->GetBinContent(0)!=0.0) {
			pull =  hNum->GetBinContent(0) -SCL->GetBinContent(0);
			err2 =  hNum->GetBinError(0)*hNum->GetBinError(0)
				   +SCL->GetBinError(0)*SCL->GetBinError(0);
			//cout << "Bin " << i << " Pull " << pull/sqrt(err2) << endl;
			if (err2<=0.0) {
				cout << "Zero error, non-zero content at UNDERFLOW\n";
				err2 = hNum->GetBinContent(0) +SCL->GetBinContent(0);			
			}
			Chi2 += pull*pull/err2;
			++NDoF;
		}
		if (hNum->GetBinContent(OVR)!=0.0 || SCL->GetBinContent(OVR)!=0.0) {
			pull =  hNum->GetBinContent(OVR) -SCL->GetBinContent(OVR);
			err2 =  hNum->GetBinError(OVR)*hNum->GetBinError(OVR)
				   +SCL->GetBinError(OVR)*SCL->GetBinError(OVR);
			//cout << "Bin " << i << " Pull " << pull/sqrt(err2) << endl;
			if (err2<=0.0) {
				cout << "Zero error, non-zero content at OVERFLOW\n";
				err2 = hNum->GetBinContent(OVR) +SCL->GetBinContent(OVR);
			}			
			Chi2 += pull*pull/err2;
			++NDoF;
		}
		cout << "Chi2/NDoF: " << Chi2 << " / " << NDoF << "   C.L.: " 
			 << TMath::Prob(Chi2,NDoF);
		cout << "  with over/underflows\n";
		return 0;
	}



};
