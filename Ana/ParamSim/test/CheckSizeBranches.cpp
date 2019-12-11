void CheckSizeBranches()
{
    // CaliceStyle();

    //Compare between ptrue <-> preco and angle <-> anglereco
    TChain *chain = new TChain("caf");

    // for(int i = 5000; i < 5050; i++){
    //     TString filename = TString::Format("../Cafs/caf_%i.root", i);
    //     chain->Add(filename);
    // }

    chain->Add("caf.root");

    int evt;
    std::vector<double>* mctime = 0;
    std::vector<int> *_nFSP = 0;
    std::vector<int> *detected = 0;
    std::vector<int> *pdgmother = 0;
    std::vector<int> *truepdg = 0;
    std::vector<int> * _MCPStartX = 0;
    std::vector<int> *_MCPStartY = 0;
    std::vector<int> *_MCPStartZ = 0;
    std::vector<int> *_MCPEndX = 0;
    std::vector<int> *_MCPEndY = 0;
    std::vector<int> *_MCPEndZ = 0;
    std::vector<std::string> *_MCProc = 0;
    std::vector<std::string> *_MCEndProc = 0;
    std::vector<float> *trkLen = 0;
    std::vector<float> *trkLenPerp = 0;
    std::vector<float> *truep = 0;
    std::vector<float> *truepx = 0;
    std::vector<float> *truepy = 0;
    std::vector<float> *truepz = 0;
    std::vector<float> *_angle = 0;
    //Reco values
    std::vector<int> *recopid = 0;
    std::vector<float> *prob_arr = 0;
    std::vector<float> *partereco = 0;
    std::vector<float> *anglereco = 0;
    std::vector<float> *_preco = 0;
    std::vector<float> *erecon = 0;

    chain->SetBranchAddress("Event", &evt);
    chain->SetBranchAddress("nFSP", &_nFSP);
    chain->SetBranchAddress("detected", &detected);
    chain->SetBranchAddress("pdgmother", &pdgmother);
    chain->SetBranchAddress("MCPTime", &mctime);
    chain->SetBranchAddress("MCPStartX", &_MCPStartX);
    chain->SetBranchAddress("MCPStartY", &_MCPStartY);
    chain->SetBranchAddress("MCPStartZ", &_MCPStartZ);
    chain->SetBranchAddress("truepx", &truepx);
    chain->SetBranchAddress("truepy", &truepy);
    chain->SetBranchAddress("truepz", &truepz);
    chain->SetBranchAddress("MCPEndX", &_MCPEndX);
    chain->SetBranchAddress("MCPEndY", &_MCPEndY);
    chain->SetBranchAddress("MCPEndZ", &_MCPEndZ);
    chain->SetBranchAddress("MCProc", &_MCProc);
    chain->SetBranchAddress("MCEndProc", &_MCEndProc);
    chain->SetBranchAddress("angle", &_angle);
    chain->SetBranchAddress("truep", &truep);
    chain->SetBranchAddress("truepdg", &truepdg);
    //Reco info
    chain->SetBranchAddress("recopid", &recopid);
    chain->SetBranchAddress("trkLen", &trkLen);
    chain->SetBranchAddress("trkLenPerp", &trkLenPerp);
    chain->SetBranchAddress("preco", &_preco);
    chain->SetBranchAddress("anglereco", &anglereco);
    chain->SetBranchAddress("erecon", &erecon);
    chain->SetBranchAddress("prob_arr", &prob_arr);

    for(int itree = 0; itree < 5; itree++)
    {
        chain->GetEntry(itree);

        std::cout << "Event " << evt << std::endl;

        std::cout << "Size pdgmother: " << pdgmother->size() << std::endl;
        std::cout << "Size MCPTime: " << mctime->size() << std::endl;
        std::cout << "Size MCPStartX: " << _MCPStartX->size() << std::endl;
        std::cout << "Size MCPStartY: " << _MCPStartY->size() << std::endl;
        std::cout << "Size MCPStartZ: " << _MCPStartZ->size() << std::endl;
        std::cout << "Size truepx: " << truepx->size() << std::endl;
        std::cout << "Size truepy: " << truepy->size() << std::endl;
        std::cout << "Size truepz: " << truepz->size() << std::endl;
        std::cout << "Size MCPEndX: " << _MCPEndX->size() << std::endl;
        std::cout << "Size MCPEndY: " << _MCPEndY->size() << std::endl;
        std::cout << "Size MCPEndZ: " << _MCPEndZ->size() << std::endl;
        std::cout << "Size _MCProc: " << _MCProc->size() << std::endl;
        std::cout << "Size _MCEndProc: " << _MCEndProc->size() << std::endl;
        std::cout << "Size _angle: " << _angle->size() << std::endl;
        std::cout << "Size truep: " << truep->size() << std::endl;
        std::cout << "Size truepdg: " << truepdg->size() << std::endl;
        std::cout << "Size recopid: " << recopid->size() << std::endl;
        std::cout << "Size trkLen: " << trkLen->size() << std::endl;
        std::cout << "Size trkLenPerp: " << trkLenPerp->size() << std::endl;
        std::cout << "Size _preco: " << _preco->size() << std::endl;
        std::cout << "Size anglereco: " << anglereco->size() << std::endl;
        std::cout << "Size erecon: " << erecon->size() << std::endl;
        std::cout << "Size prob_arr: " << prob_arr->size() << std::endl;
        std::cout << "Size _nFSP: " << _nFSP->size() << " value " << _nFSP->at(0) << std::endl;
        std::cout << "Size detected: " << detected->size() << std::endl;
    }

}
