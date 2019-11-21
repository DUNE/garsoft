void checkmuons()
{
    // CaliceStyle();

    //Compare between ptrue <-> preco and angle <-> anglereco
    TChain *chain = new TChain("caf");

    // for(int i = 5000; i < 5050; i++){
    //     TString filename = TString::Format("../Cafs/caf_%i.root", i);
    //     chain->Add(filename);
    // }

    chain->Add("../caf.root");

    std::vector<float> *truep = 0;
    std::vector<float> *preco = 0;
    std::vector<int> *truepdg = 0;
    std::vector<int> *recopid = 0;
    std::vector<float> *angle = 0;
    std::vector<float> *anglereco = 0;
    std::vector<float> *erecon = 0;

    chain->SetBranchAddress("truep", &truep);
    chain->SetBranchAddress("preco", &preco);
    chain->SetBranchAddress("truepdg", &truepdg);
    chain->SetBranchAddress("recopid", &recopid);
    chain->SetBranchAddress("angle", &angle);
    chain->SetBranchAddress("anglereco", &anglereco);
    chain->SetBranchAddress("erecon", &erecon);

    TH1F* hPullMomentum = new TH1F("hPullMomentum", "hPullMomentum", 300, -1, 1);
    hPullMomentum->SetLineWidth(2);
    hPullMomentum->SetLineColor(kBlack);

    TH2F* hMomemtumTruevsReco = new TH2F("hMomemtumTruevsReco", "hMomemtumTruevsReco", 100, 0, 5, 100, 0, 5);

    TH1F* hPullAngle = new TH1F("hPullAngle", "hPullAngle", 100, -1, 1);
    hPullAngle->SetLineWidth(2);
    hPullAngle->SetLineColor(kBlack);

    TH2F* hAngleTruevsReco = new TH2F("hAngleTruevsReco", "hAngleTruevsReco", 100, -1, 1, 100, -1, 1);

    TH1F* hPullNeutrons = new TH1F("hPullNeutrons", "hPullNeutrons", 100, -1, 1);
    hPullNeutrons->SetLineWidth(2);
    hPullNeutrons->SetLineColor(kBlack);

    for(int itree = 0; itree < chain->GetEntries(); itree++)
    {
        chain->GetEntry(itree);

        for(int i = 0; i < truep->size(); i++)
        {
            //Muons
            if( std::abs(truepdg->at(i)) == 13 && std::abs(recopid->at(i)) == 13){
                hPullMomentum->Fill( (truep->at(i) - preco->at(i)) / truep->at(i) );
                hMomemtumTruevsReco->Fill(truep->at(i), preco->at(i));
                hPullAngle->Fill( (angle->at(i) - anglereco->at(i)) / angle->at(i) );
                hAngleTruevsReco->Fill(angle->at(i), anglereco->at(i));
            }
        }
    }

    gStyle->SetOptStat(1110);
    gStyle->SetOptFit(1);
    hPullMomentum->Scale(1./hPullMomentum->Integral());

    hPullMomentum->SetStats(kTRUE);
    hPullMomentum->Fit("gaus", "EMR+", "", -0.8, 0.8);
    TCanvas *c1 = new TCanvas("c1", "Pull p", 800, 600);
    hPullMomentum->GetXaxis()->SetTitle( "(p_{true} - p_{reco}) / p_{true}" );
    hPullMomentum->GetYaxis()->SetTitle( "Normalized Events" );
    hPullMomentum->Draw("hist");
    hPullMomentum->GetFunction("gaus")->SetLineColor(kRed);
    hPullMomentum->GetFunction("gaus")->SetLineStyle(2);
    hPullMomentum->GetFunction("gaus")->Draw("same");
    // c1->SaveAs("PullMomentumMuons.pdf");

    TCanvas *c2 = new TCanvas("c2", "Mom truth vs reco", 800, 600);
    gPad->SetRightMargin(0.12);
    hMomemtumTruevsReco->GetXaxis()->SetTitle( "p_{true} [GeV]" );
    hMomemtumTruevsReco->GetYaxis()->SetTitle( "p_{reco} [GeV]" );
    hMomemtumTruevsReco->Draw("COLZ");
    // c2->SaveAs("MomentumMuons_TruevsReco.pdf");

    hPullAngle->Scale(1./hPullAngle->Integral());
    hPullAngle->SetStats(kTRUE);

    TCanvas *c3 = new TCanvas("c3", "Pull angle", 800, 600);
    hPullAngle->GetXaxis()->SetTitle( "(angle_{true} - angle_{reco}) / angle_{true}" );
    hPullAngle->GetYaxis()->SetTitle( "Normalized Events" );
    hPullAngle->Draw("hist");
    // c3->SaveAs("PullAngleMuons.pdf");

    TCanvas *c4 = new TCanvas("c4", "Angle truth vs reco", 800, 600);
    gPad->SetRightMargin(0.12);
    hAngleTruevsReco->GetXaxis()->SetTitle( "angle_{true} [GeV]" );
    hAngleTruevsReco->GetYaxis()->SetTitle( "angle_{reco} [GeV]" );
    hAngleTruevsReco->Draw("COLZ");
    // c4->SaveAs("AngleMuons_TruevsReco.pdf");

    // TCanvas *c5 = new TCanvas("c5", "Pull Neutrons", 800, 600);
    // hPullNeutrons->GetXaxis()->SetTitle( "(E_{true} - E_{reco}) / E_{true}" );
    // hPullNeutrons->GetYaxis()->SetTitle( "Normalized Events" );
    // hPullNeutrons->Draw("hist");
    // c5->SaveAs("PullEnergyNeutrons.pdf");

}
