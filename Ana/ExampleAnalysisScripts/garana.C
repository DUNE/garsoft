#define garana_cxx
#include "garana.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TString.h>
#include <TVector3.h>
#include <TH2F.h>
#include <TMath.h>
#include <TROOT.h>
#include <TPad.h>

void garana::Loop()
{
//   In a ROOT session, you can do:
//      root> .L garana.C
//      root> garana t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   TVector3 detcent(0,218.196,1124.02);
   TVector3 xhat(1,0,0);
   float detR = 250;

   TH1F *allmuonr = new TH1F("allMuonR","Radius of Start of Muon from Center of Detector",30,0,300);
   TH1F *recomuonr = new TH1F("recoMuonR","Radius of Start of Muon from Center of Detector",30,0,300);

   TH1F *allmuonpt = new TH1F("allMuonPt",";Pt of muon (GeV)",30,0,2);
   TH1F *allmuonpx = new TH1F("allMuonPx",";Px of muon (GeV)",30,-2,2);
   TH1F *allmuonsinx = new TH1F("allMuonSinx","Px/P all muons",30,-1,1); 
   TH1F *recomuonpt = new TH1F("recoMuonPt",";Pt of muon (GeV)",30,0,2);
   TH1F *recomuonpx = new TH1F("recoMuonPx",";Px of muon (GeV)",30,-2,2);
   TH1F *recomuonsinx = new TH1F("recoMuonSinx","Px/P reco muons",30,-1,1); 
   TH1F *nchprimh = new TH1F("nchprim","Number of P.V. charged particles pt gt 10 MeV",30,-0.5,29.5);
   TH1F *othermcpt = new TH1F("otherMCpt",";Pt of non-muon charged particles (GeV)",30,0,2);
   TH1F *othermcpx = new TH1F("otherMCpx",";Px of non-muon charged particles (GeV)",30,-2,2);

   TH1F *alltrackradstart = new TH1F("allTrackR","Radius of Track Start/End from Center of Detector",30,0,300);
   TH1F *alltrackpt = new TH1F("alltrackpt",";Track Pt (GeV)",30,0,3);
   TH1F *alltrackpx = new TH1F("alltarckpx",";Track Px (GeV)",30,-3,3);

   TH2F *primyz = new TH2F("primyz","Primary vertex;Z(cm);Y(cm)",50,detcent.Z()-detR,detcent.Z()+detR,50,detcent.Y()-detR,detcent.Y()+detR);
   TH2F *tendpoint = new TH2F("tendpoint","Track Endpoints;Z(cm);Y(cm)",50,detcent.Z()-detR,detcent.Z()+detR,50,detcent.Y()-detR,detcent.Y()+detR);
   TH2F *recovtxyz = new TH2F("recovtxyz","Reco Vertices;Z(cm);Y(cm)",50,detcent.Z()-detR,detcent.Z()+detR,50,detcent.Y()-detR,detcent.Y()+detR);

   TH1F *minryzvertex = new TH1F("minryzvertex","min YZ distance of reco vertex to primary;cm",30,-0.11,3); 

   TH1F *ntracksreco = new TH1F("ntracksreco",";Reconstructed Tracks",30,-0.5,29.5);
   TH1F *nvtxreco = new TH1F("nvtxreco",";Reconstructed Vertices",10,-0.5,9.5);

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if ( CCNC->at(0) !=0 ) continue;  // only do charged current for now.
      if ( std::abs( PDG->at(0) ) != 13 ) continue;  // only do muons

      TVector3 primvtx(MCPStartX->at(0),MCPStartY->at(0),MCPStartZ->at(0));
      TVector3 muonp(MCPPX->at(0),MCPPY->at(0),MCPPZ->at(0));
      allmuonr->Fill( (primvtx-detcent).Mag() );
      //cout << muonp.Perp() << endl;
      allmuonpt->Fill( (muonp.Cross(xhat)).Mag() );
      allmuonpx->Fill( muonp.X() );
      primyz->Fill(primvtx.Z(),primvtx.Y());
      float pxoverp = muonp.X()/muonp.Mag();
      allmuonsinx->Fill(pxoverp);

      // fill MC histograms

      int nchprim = 0;
      for (size_t i=0; i< MCPStartX->size(); ++i)
	{
	  TVector3 mcstart(MCPStartX->at(i),MCPStartY->at(i),MCPStartZ->at(i));
	  TVector3 mcp(MCPPX->at(i),MCPPY->at(i),MCPPZ->at(i));
	  float mcpt = (mcp.Cross(xhat)).Mag();
	  float mcpx = mcp.X();
	  int ipdg = TMath::Abs(PDG->at(i));
	  if ( (mcstart-primvtx).Mag() < 0.1 && mcpt > 0.01)
	    {
	      if (ipdg == 11 || ipdg == 13 || ipdg == 211 || ipdg == 2212 || ipdg == 321)
		{
		  ++nchprim;
		  if (ipdg != 13)
		    {
		      othermcpt->Fill(mcpt);
		      othermcpx->Fill(mcpx);
		    }
		}
	    }
	}
      nchprimh->Fill(nchprim);

      int ntreco = 0;
      for (size_t i=0; i< TrackStartX->size(); ++i)
	{
	  ++ntreco;
	  TVector3 trackstart(TrackStartX->at(i),TrackStartY->at(i),TrackStartZ->at(i));
	  TVector3 trackend(TrackEndX->at(i),TrackEndY->at(i),TrackEndZ->at(i));
	  TVector3 trackp(TrackStartPX->at(i),TrackStartPY->at(i),TrackStartPZ->at(i));
	  TVector3 trackpend(TrackEndPX->at(i),TrackEndPY->at(i),TrackEndPZ->at(i));
	  std::cout << "Track start " << trackstart.X() << " " << trackstart.Y() << " " << trackstart.Z() << 
	    " P: " << trackp.X() << " " << trackp.Y() << " " << trackp.Z() << std::endl;
	  std::cout << "Track end   " << trackend.X() << " " << trackend.Y() << " " << trackend.Z() << 
	    " P: " << trackpend.X() << " " << trackpend.Y() << " " << trackpend.Z() << std::endl;
	  std::cout << "Primvtx " << primvtx.X() << " " << primvtx.Y() << " " << primvtx.Z() << 
	    " Pmuon: " << muonp.X() << " " << muonp.Y() << " " << muonp.Z() << std::endl;

	  // cosine of the angle of the track and the muon

	  float ctmb = trackp.Dot(muonp)/(trackp.Mag()*muonp.Mag());
	  float ctme = trackpend.Dot(muonp)/(trackpend.Mag()*muonp.Mag());
	  float missb = (trackp.Cross(trackstart-primvtx)).Mag()/(trackp.Mag());
	  float misse = (trackpend.Cross(trackend-primvtx)).Mag()/(trackpend.Mag());
	  std::cout << "ctmb: " << ctmb << " ctme " << ctme << " missb " << missb << " misse: " << misse << std::endl << std::endl;
 
	  if ( (missb < 10 && ctmb > 0.95) || (misse < 10 && ctme > 0.95) )
	    {
	      recomuonr->Fill( (primvtx - detcent).Mag() );
	      recomuonpt->Fill( (muonp.Cross(xhat)).Mag() );
              float recopxoverp = muonp.X()/muonp.Mag();
              recomuonsinx->Fill(pxoverp);

	      break;
	    }
	  alltrackradstart->Fill( (trackstart - detcent).Mag() );
	  alltrackradstart->Fill( (trackend - detcent).Mag() );
	  alltrackpt->Fill( (trackp.Cross(xhat)).Mag() );
	  alltrackpx->Fill( trackp.X() );
	  alltrackpt->Fill( (trackpend.Cross(xhat)).Mag() );
	  alltrackpx->Fill( trackpend.X() );
	  tendpoint->Fill(trackstart.Z(),trackstart.Y());
	  tendpoint->Fill(trackend.Z(),trackstart.Y());
	}
      ntracksreco->Fill(ntreco);

      int nvreco = 0;
      float minryz=-0.1;
      for (size_t i=0; i< VertX->size(); ++i)
	{
	  nvreco++;
	  recovtxyz->Fill(VertZ->at(i),VertY->at(i));
	  TVector3 rvpos(VertX->at(i),VertY->at(i),VertZ->at(i));
	  float ryz = ((rvpos-primvtx).Cross(xhat)).Mag();
	  //cout << ryz << endl;
	  if (i==0) 
	    {
	      minryz = ryz;
	    }
	  else
	    {
	      minryz = TMath::Min(ryz,minryz);
	    }
	}
      minryzvertex->Fill(minryz);

   } // end loop over tree entries

   TCanvas *mccanvas1 = new TCanvas("mccanvas1","",1000,800);
   mccanvas1->Divide(3,3);
   mccanvas1->cd(1);
   allmuonr->Draw("hist");
   recomuonr->SetLineColor(2);
   recomuonr->Draw("hist,same");
   mccanvas1->cd(2);
   allmuonpt->Draw("hist");
   recomuonpt->SetLineColor(2);
   recomuonpt->Draw("hist,same");
   mccanvas1->cd(3);
   nchprimh->Draw("hist");
   mccanvas1->cd(4);
   allmuonpx->Draw("hist");
   mccanvas1->cd(5);
   othermcpt->Draw("hist");
   mccanvas1->cd(6);
   othermcpx->Draw("hist");
   mccanvas1->cd(7);
   allmuonsinx->Draw("hist");
   recomuonsinx->SetLineColor(2);
   recomuonsinx->Draw("hist,same");
   mccanvas1->Print("mccanvas1.png");

   TCanvas *trackcanvas1 = new TCanvas("trackcanvas1","",1000,800);
   trackcanvas1->Divide(2,2);
   trackcanvas1->cd(1);
   ntracksreco->Draw("hist");
   trackcanvas1->cd(2);
   alltrackradstart->Draw("hist");
   trackcanvas1->cd(3);
   alltrackpt->Draw("hist");
   trackcanvas1->cd(4);
   alltrackpx->Draw("hist");
   trackcanvas1->Print("trackcanvas1.png");

   TCanvas *vertexcanvas = new TCanvas("Vertices","Vertices",1000,800);
   vertexcanvas->Divide(2,2);
   vertexcanvas->cd(1);
   primyz->Draw("box");
   vertexcanvas->cd(2);
   recovtxyz->Draw("box");
   vertexcanvas->cd(3);
   tendpoint->Draw("box");
   vertexcanvas->cd(4);
   gPad->SetLogy();
   minryzvertex->Draw("hist");

   vertexcanvas->Print("vertexcanvas.png");

}
