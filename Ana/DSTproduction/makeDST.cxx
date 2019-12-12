#include <list>

#include "globbing.h"
#include "treeReader.h"
#include "includeROOT.h"

#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "ReconstructionDataProducts/Track.h"
#include "ReconstructionDataProducts/IDNumberGen.h"

void TransMogrifyTree(int fileNumber, bool firstCall);



// In a fiducial corresponding to 1 metric ton of Ar nuclei?  (Corresponds to
// about 35cm from edge of the TPC).
bool inFiducial(double x, double y, double z) {
	double const Rcut = 445.0 /2.0;		double const Xcut = 430.0 /2.0;
	double r = sqrt( z*z + y*y );
	if ( r > Rcut ) return false;
	if (fabs(x) > Xcut ) return false;
	return true;
}



//==============================================================================
//==============================================================================
#define writeECAL
bool  viaXrootd = true;
float  PmagCut = 0.150;		// In GeV.



TFile* outFile;            TFile* inFile;
TTree* outTree;            TTree* inTree;



enum CounterTag   {nEvent,      nDatainFid,
                   nCCQEinFid,  nNCQEinFid,
				   nRESinFid,   nCCDISinFid,
				   nNCDISinFid, nCCCOHinFid,
				   nNCCOHinFid, nNuEinFid,
nCounterTag};  // Must be last enum
int Counter[nCounterTag];



//==============================================================================
//==============================================================================
int main(int argc , const char* argv[]){
    int maxNfiles = -1;		string outName, SearchDir;	bool gotOutDir = false;
    int iArg = 1;
    while ((iArg < argc) && (argv[iArg][0] == '-')){
        switch (argv[iArg][1]) {

        case 'f':
            maxNfiles = atoi(argv[iArg+1]);
			cout << "Only processing " << maxNfiles << " files." << endl;
            if (maxNfiles < 1) {
                cout << "Wait! " << maxNfiles << "files!  That's not right!" << endl;
                exit(1);
            }
            ++iArg;
            break;

		case 'c':
			PmagCut = std::stof(argv[iArg+1]);
			cout << "Removing tracks and vertices with tracks with P below " << 
				PmagCut << " MeV." << endl;
			PmagCut /= 1000.0;
			++iArg;
			break;

		case 'L':
	    	SearchDir =  argv[iArg+1];
			viaXrootd = false;
			++iArg;
			break;

		case 'o':
			// Cleverly, argv comes to us with environmental variables expanded
			// and we don't need wordexp.h.
			outName = argv[iArg+1];
			gotOutDir = true;			
			++iArg;
			break;

        case 'h':
        default:
            cout << "[-h]       This message" << endl;
            cout << "[-o name]  Name for output files; TTree will be in <name>.root," << endl;
			cout << "               statistics will be in <name>.txt.  Environmental" << endl;
			cout << "               variables in <name> will be expanded." << endl;
            cout << "[-f n]     Analyze n background files, max" << endl;
            cout << "[-c m]     Set Ptrk cut to Ptrk > m MeV for output." << endl;
            cout << "               The default cut is 150 MeV." << endl;
            cout << "[-L name]  Input files from directory <name> without using" << endl;
			cout << "               xrootd rather than from files in $PWD, with" << endl;
			cout << "               xrootd, which is the default." << endl;
            exit(0);
        }
        ++iArg;
    }
	if (!gotOutDir) {
		cout << "Did not see -o <dir> so no output directory selected." << endl;
		exit(1);
	}



    // Get list of input files
    string LookFor;
	vector<string> fileList;
    int nFiles;

	if (viaXrootd){ 
		SearchDir = "";
		LookFor   = "anatree_*.root";
		fileList  = globbing(SearchDir +LookFor);
		nFiles    = (int)fileList.size();
		string newAccess = "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr";
		string PWD = getenv("PWD");		PWD.erase(0,5);		PWD += "/";
		newAccess += PWD;

		for (int iFile=0; iFile<nFiles; ++iFile) {
			fileList[iFile].insert(0,newAccess);
		}
	} else {
		LookFor   = "anatree_*.root";
 		fileList  = globbing(SearchDir +LookFor);
		nFiles = (int)fileList.size();
   }
    cout << ( nFiles = (int)fileList.size() ) << " total files found.\n";
    if (nFiles==0) exit(1);

    if ( maxNfiles>0 && maxNfiles<nFiles) {
        nFiles = maxNfiles;
    }
    cout << endl;



    // Output file made and opened here
	string RootOutName  = outName + ".root";
    if (!viaXrootd && !filehamna(RootOutName)) {
        cout << "Already existing " << RootOutName << endl;
        exit(2);
    }
    outFile = TFile::Open(RootOutName.c_str(),"NEW");
    if (outFile->IsZombie()) {cout << "Eat you brains!" << endl; return 1;};
    outFile->mkdir("anatree");
    outFile->cd("anatree");
    outTree = new TTree("GArAnaTree","GArAnaTree");



    for (int iFile=0; iFile<nFiles; ++iFile) {
        cout << "Reading file " << iFile+1 << ": " << fileList[iFile].c_str() << " ";
        inFile = TFile::Open(fileList[iFile].c_str(),"READ");
        if (!inFile || !inFile->IsOpen()) {
            cout << "...but that file can not be opened!" << endl;
            exit(3);
        }
		
		size_t strtLoc = fileList[iFile].find("_") +1;
		size_t lenLoc  = fileList[iFile].find(".root") -strtLoc;
		int fileNumber;
		string numberAsString = fileList[iFile].substr(strtLoc,lenLoc);
		sscanf(numberAsString.c_str(), "%d", &fileNumber);

        TDirectory* dir = (TDirectory*)(inFile)->Get("anatree");
        dir->GetObject("GArAnaTree",inTree);
        TransMogrifyTree(fileNumber, iFile==0);
        inFile->Close();
    }
	outFile->Write();
	outFile->Close();



    string TextOutName = outName + ".txt";
    #define OUTDEV TextFileOut
    std::ofstream TextFileOut(TextOutName.c_str());
	// #define OUTDEV cout
	cout << endl << endl;

    OUTDEV << Counter[nEvent]      << "\tevents" << endl;

	OUTDEV << Counter[nDatainFid]  << "\tevents with MC vertex in fiducial"
		<< endl;

	OUTDEV << Counter[nCCQEinFid]  << "\tCCQE events with MC vertex in fiducial"
		<< endl;

	OUTDEV << Counter[nNCQEinFid]  << "\tNCQE events with MC vertex in fiducial"
		<< endl;

	OUTDEV << Counter[nRESinFid]   << "\tRES events with MC vertex in fiducial"
		<< endl;

	OUTDEV << Counter[nCCDISinFid] << "\tCC DIS events with MC vertex in fiducial"
		<< endl;

	OUTDEV << Counter[nNCDISinFid] << "\tNC DIS events with MC vertex in fiducial"
		<< endl;

	OUTDEV << Counter[nCCCOHinFid] << "\tCC coherent pi events with MC vertex in fiducial"
		<< endl;

	OUTDEV << Counter[nNCCOHinFid] << "\tNC coherent pi events with MC vertex in fiducial"
		<< endl;

	OUTDEV << Counter[nNuEinFid]   << "\tElastic nu-e events with MC vertex in fiducial"
		<< endl;

    TextFileOut.close();

    exit(0);
}





//==============================================================================
//==============================================================================
void TransMogrifyTree(int fileNumber, bool firstCall) {



    Long64_t iEntry;

    Int_t                      fFileNumber;
	if (firstCall) outTree->Branch("FileNumber",&fFileNumber,"FileNumber/I");

    scalarFromTree<Int_t>       Event (inTree,"Event",&iEntry);
    Int_t                      fEvent;
	if (firstCall) outTree->Branch("Event",&fEvent,"Event/I");



    vectorFromTree<Int_t>       NType(inTree,"NType",&iEntry);
    vector<Int_t>              fNType;
	if (firstCall) outTree->Branch("NType",&fNType);

    vectorFromTree<Int_t>       Mode(inTree,"Mode",&iEntry);
    vector<Int_t>              fMode;
	if (firstCall) outTree->Branch("Mode",&fMode);

    vectorFromTree<int>         InterT(inTree,"InterT",&iEntry);
    vector<int>                 fInterT;
	if (firstCall) outTree->Branch("InterT",&fInterT);

    vectorFromTree<Float_t>     MCVertX(inTree,"MCVertX",&iEntry);
    vector<Float_t>            fMCVertX;
	if (firstCall) outTree->Branch("MCVertX",&fMCVertX);

    vectorFromTree<Float_t>     MCVertY(inTree,"MCVertY",&iEntry);
    vector<Float_t>            fMCVertY;
	if (firstCall) outTree->Branch("MCVertY",&fMCVertY);

    vectorFromTree<Float_t>     MCVertZ(inTree,"MCVertZ",&iEntry);
    vector<Float_t>            fMCVertZ;
	if (firstCall) outTree->Branch("MCVertZ",&fMCVertZ);

    vectorFromTree<Float_t>     MCNuPx(inTree,"MCNuPx",&iEntry);
    vector<Float_t>            fMCNuPx;
	if (firstCall) outTree->Branch("MCNuPx",&fMCNuPx);

    vectorFromTree<Float_t>     MCNuPy(inTree,"MCNuPy",&iEntry);
    vector<Float_t>            fMCNuPy;
	if (firstCall) outTree->Branch("MCNuPy",&fMCNuPy);

    vectorFromTree<Float_t>     MCNuPz(inTree,"MCNuPz",&iEntry);
    vector<Float_t>            fMCNuPz;
	if (firstCall) outTree->Branch("MCNuPz",&fMCNuPz);

    vectorFromTree<Float_t>     MC_T(inTree,"MC_T",&iEntry);
    vector<Float_t>            fMC_T;
	if (firstCall) outTree->Branch("MC_T",&fMC_T);



    vectorFromTree<Int_t>       PDG(inTree,"PDG",&iEntry);
    vector<Int_t>              fPDG;
	if (firstCall) outTree->Branch("PDG",&fPDG);

    vectorFromTree<Int_t>       PDGMother(inTree,"PDGMother",&iEntry);
    vector<Int_t>              fPDGMother;
	if (firstCall) outTree->Branch("PDGMother",&fPDGMother);

    vectorFromTree<Float_t>     MCPStartX(inTree,"MCPStartX",&iEntry);
    vector<Float_t>            fMCPStartX;
	if (firstCall) outTree->Branch("MCPStartX",&fMCPStartX);

    vectorFromTree<Float_t>     MCPStartY(inTree,"MCPStartY",&iEntry);
    vector<Float_t>            fMCPStartY;
	if (firstCall) outTree->Branch("MCPStartY",&fMCPStartY);

    vectorFromTree<Float_t>     MCPStartZ(inTree,"MCPStartZ",&iEntry);
    vector<Float_t>            fMCPStartZ;
	if (firstCall) outTree->Branch("MCPStartZ",&fMCPStartZ);

    vectorFromTree<Float_t>     MCPStartPX(inTree,"MCPStartPX",&iEntry);
    vector<Float_t>            fMCPStartPX;
	if (firstCall) outTree->Branch("MCPStartPX",&fMCPStartPX);

    vectorFromTree<Float_t>     MCPStartPY(inTree,"MCPStartPY",&iEntry);
    vector<Float_t>            fMCPStartPY;
	if (firstCall) outTree->Branch("MCPStartPY",&fMCPStartPY);

    vectorFromTree<Float_t>     MCPStartPZ(inTree,"MCPStartPZ",&iEntry);
    vector<Float_t>            fMCPStartPZ;
	if (firstCall) outTree->Branch("MCPStartPZ",&fMCPStartPZ);



    vectorFromTree<gar::rec::IDNumber>
	                            TrackIDNumber(inTree,"TrackIDNumber",&iEntry);
    vector<gar::rec::IDNumber> fTrackIDNumber;
	if (firstCall) outTree->Branch("TrackIDNumber",&fTrackIDNumber);

    vectorFromTree<Float_t>     TrackStartX(inTree,"TrackStartX",&iEntry);
    vector<Float_t>            fTrackStartX;
	if (firstCall) outTree->Branch("TrackStartX",&fTrackStartX);

    vectorFromTree<Float_t>     TrackStartY(inTree,"TrackStartY",&iEntry);
    vector<Float_t>            fTrackStartY;
	if (firstCall) outTree->Branch("TrackStartY",&fTrackStartY);

    vectorFromTree<Float_t>     TrackStartZ(inTree,"TrackStartZ",&iEntry);
    vector<Float_t>            fTrackStartZ;
	if (firstCall) outTree->Branch("TrackStartZ",&fTrackStartZ);

    vectorFromTree<Float_t>     TrackStartPX(inTree,"TrackStartPX",&iEntry);
    vector<Float_t>            fTrackStartPX;
	if (firstCall) outTree->Branch("TrackStartPX",&fTrackStartPX);

    vectorFromTree<Float_t>     TrackStartPY(inTree,"TrackStartPY",&iEntry);
    vector<Float_t>            fTrackStartPY;
	if (firstCall) outTree->Branch("TrackStartPY",&fTrackStartPY);

    vectorFromTree<Float_t>     TrackStartPZ(inTree,"TrackStartPZ",&iEntry);
    vector<Float_t>            fTrackStartPZ;
	if (firstCall) outTree->Branch("TrackStartPZ",&fTrackStartPZ);

    vectorFromTree<Int_t>       TrackStartQ(inTree,"TrackStartQ",&iEntry);
    vector<Int_t>              fTrackStartQ;
	if (firstCall) outTree->Branch("TrackStartQ",&fTrackStartQ);

    vectorFromTree<Float_t>     TrackEndX(inTree,"TrackEndX",&iEntry);
    vector<Float_t>            fTrackEndX;
	if (firstCall) outTree->Branch("TrackEndX",&fTrackEndX);

    vectorFromTree<Float_t>     TrackEndY(inTree,"TrackEndY",&iEntry);
    vector<Float_t>            fTrackEndY;
	if (firstCall) outTree->Branch("TrackEndY",&fTrackEndY);

    vectorFromTree<Float_t>     TrackEndZ(inTree,"TrackEndZ",&iEntry);
    vector<Float_t>            fTrackEndZ;
	if (firstCall) outTree->Branch("TrackEndZ",&fTrackEndZ);

    vectorFromTree<Float_t>     TrackEndPX(inTree,"TrackEndPX",&iEntry);
    vector<Float_t>            fTrackEndPX;
	if (firstCall) outTree->Branch("TrackEndPX",&fTrackEndPX);

    vectorFromTree<Float_t>     TrackEndPY(inTree,"TrackEndPY",&iEntry);
    vector<Float_t>            fTrackEndPY;
	if (firstCall) outTree->Branch("TrackEndPY",&fTrackEndPY);

    vectorFromTree<Float_t>     TrackEndPZ(inTree,"TrackEndPZ",&iEntry);
    vector<Float_t>            fTrackEndPZ;
	if (firstCall) outTree->Branch("TrackEndPZ",&fTrackEndPZ);

    vectorFromTree<Int_t>       TrackEndQ(inTree,"TrackEndQ",&iEntry);
    vector<Int_t>              fTrackEndQ;
	if (firstCall) outTree->Branch("TrackEndQ",&fTrackEndQ);



    vectorFromTree<Float_t>     TrackLenF(inTree,"TrackLenF",&iEntry);
    vector<Float_t>            fTrackLenF;
	if (firstCall) outTree->Branch("TrackLenF",&fTrackLenF);

    vectorFromTree<Float_t>     TrackLenB(inTree,"TrackLenB",&iEntry);
    vector<Float_t>            fTrackLenB;
	if (firstCall) outTree->Branch("TrackLenB",&fTrackLenB);

    vectorFromTree<Int_t>       NTPCClustersOnTrack(inTree,"NTPCClustersOnTrack",&iEntry);
    vector<Int_t>              fNTPCClustersOnTrack;
	if (firstCall) outTree->Branch("NTPCClustersOnTrack",&fNTPCClustersOnTrack);

    vectorFromTree<Float_t>     TrackAvgIonF(inTree,"TrackAvgIonF",&iEntry);
    vector<Float_t>            fTrackAvgIonF;
	if (firstCall) outTree->Branch("TrackAvgIonF",&fTrackAvgIonF);

    vectorFromTree<Float_t>     TrackAvgIonB(inTree,"TrackAvgIonB",&iEntry);
    vector<Float_t>            fTrackAvgIonB;
	if (firstCall) outTree->Branch("TrackAvgIonB",&fTrackAvgIonB);



    vectorFromTree<gar::rec::IDNumber>
	                            VertIDNumber(inTree,"VertIDNumber",&iEntry);
    vector<gar::rec::IDNumber> fVertIDNumber;
	if (firstCall) outTree->Branch("VertIDNumber",&fVertIDNumber);

    vectorFromTree<Float_t>     VertX(inTree,"VertX",&iEntry);
    vector<Float_t>            fVertX;
	if (firstCall) outTree->Branch("VertX",&fVertX);

    vectorFromTree<Float_t>     VertY(inTree,"VertY",&iEntry);
    vector<Float_t>            fVertY;
	if (firstCall) outTree->Branch("VertY",&fVertY);

    vectorFromTree<Float_t>     VertZ(inTree,"VertZ",&iEntry);
    vector<Float_t>            fVertZ;
	if (firstCall) outTree->Branch("VertZ",&fVertZ);

    vectorFromTree<Int_t>       VertN(inTree,"VertN",&iEntry);
    vector<Int_t>              fVertN;
	if (firstCall) outTree->Branch("VertN",&fVertN);

    vectorFromTree<Int_t>       VertQ(inTree,"VertQ",&iEntry);
    vector<Int_t>              fVertQ;
	if (firstCall) outTree->Branch("VertQ",&fVertQ);



    vectorFromTree<gar::rec::IDNumber>
	                            VT_VertIDNumber(inTree,"VT_VertIDNumber",&iEntry);
    vector<gar::rec::IDNumber> fVT_VertIDNumber;
	if (firstCall) outTree->Branch("VT_VertIDNumber",&fVT_VertIDNumber);

    vectorFromTree<gar::rec::IDNumber>
	                            VT_TrackIDNumber(inTree,"VT_TrackIDNumber",&iEntry);
    vector<gar::rec::IDNumber> fVT_TrackIDNumber;
	if (firstCall) outTree->Branch("VT_TrackIDNumber",&fVT_TrackIDNumber);

    vectorFromTree<gar::rec::TrackEnd>
	                            VT_TrackEnd(inTree,"VT_TrackEnd",&iEntry);
    vector<gar::rec::TrackEnd> fVT_TrackEnd;
	if (firstCall) outTree->Branch("VT_TrackEnd",&fVT_TrackEnd);



	#ifdef writeECAL
	// Change the name to get it right
    vectorFromTree<gar::rec::IDNumber>
	                            ClusterIDNumber(inTree,"ClusterIDNumber",&iEntry);
    vector<gar::rec::IDNumber> fClusterIDNumber;
	if (firstCall) outTree->Branch("ClusterIDNumber",&fClusterIDNumber);

    vectorFromTree<UInt_t>      ClusterNhits(inTree,"ClusterNhits",&iEntry);
    vector<UInt_t>             fClusterNhits;
	if (firstCall) outTree->Branch("ClusterNhits",&fClusterNhits);

    vectorFromTree<Float_t>     ClusterEnergy(inTree,"ClusterEnergy",&iEntry);
    vector<Float_t>            fClusterEnergy;
	if (firstCall) outTree->Branch("ClusterEnergy",&fClusterEnergy);



    vectorFromTree<gar::rec::IDNumber>
								ECALAssn_ClusIDNumber(inTree,"ECALAssn_ClusIDNumber",&iEntry);
    vector<gar::rec::IDNumber> fECALAssn_ClusIDNumber;
	if (firstCall) outTree->Branch("ECALAssn_ClusIDNumber",&fECALAssn_ClusIDNumber);

    vectorFromTree<gar::rec::IDNumber>
								ECALAssn_TrackIDNumber(inTree,"ECALAssn_TrackIDNumber",&iEntry);
    vector<gar::rec::IDNumber> fECALAssn_TrackIDNumber;
	if (firstCall) outTree->Branch("ECALAssn_TrackIDNumber",&fECALAssn_TrackIDNumber);

    vectorFromTree<gar::rec::TrackEnd>
								ECALAssn_TrackEnd(inTree,"ECALAssn_TrackEnd",&iEntry);
    vector<gar::rec::TrackEnd> fECALAssn_TrackEnd;
	if (firstCall) outTree->Branch("ECALAssn_TrackEnd",&fECALAssn_TrackEnd);
	#endif




    Long64_t nEntry = inTree->GetEntries();
    cout << "nEntry:\t" << nEntry << "\tRun Number: " << fileNumber << endl;
    for (iEntry=0; iEntry<nEntry; ++iEntry) {

        ++Counter[nEvent];



        inTree->GetEntry(iEntry);

		fFileNumber             = fileNumber;
		fEvent                  = Event.getData();

		fNType                  = NType.getData();
		fMode                   = Mode.getData();
		fInterT                 = InterT.getData();
		fMCVertX                = MCVertX.getData();
		fMCVertY                = MCVertY.getData();
		fMCVertZ                = MCVertZ.getData();
		fMCNuPx                 = MCNuPx.getData();
		fMCNuPy                 = MCNuPy.getData();
		fMCNuPz                 = MCNuPz.getData();
		fMC_T                   = MC_T.getData();

		fPDG                    = PDG.getData();
		fPDGMother              = PDGMother.getData();
		fMCPStartX              = MCPStartX.getData();
		fMCPStartY              = MCPStartY.getData();
		fMCPStartZ              = MCPStartZ.getData();
		fMCPStartPX             = MCPStartPX.getData();
		fMCPStartPY             = MCPStartPY.getData();
		fMCPStartPZ             = MCPStartPZ.getData();

		fTrackIDNumber          = TrackIDNumber.getData();
		fTrackStartX		    = TrackStartX.getData();
		fTrackStartY		    = TrackStartY.getData();
		fTrackStartZ		    = TrackStartZ.getData();
		fTrackStartPX           = TrackStartPX.getData();
		fTrackStartPY           = TrackStartPY.getData();
		fTrackStartPZ           = TrackStartPZ.getData();
		fTrackStartQ            = TrackStartQ.getData();
		fTrackEndX  		    = TrackEndX.getData();
		fTrackEndY  		    = TrackEndY.getData();
		fTrackEndZ  		    = TrackEndZ.getData();
		fTrackEndPX             = TrackEndPX.getData();
		fTrackEndPY             = TrackEndPY.getData();
		fTrackEndPZ             = TrackEndPZ.getData();
		fTrackEndQ              = TrackEndQ.getData();

		fTrackLenF              = TrackLenF.getData();
		fTrackLenB              = TrackLenB.getData();
		fNTPCClustersOnTrack    = NTPCClustersOnTrack.getData();
		fTrackAvgIonF           = TrackAvgIonF.getData();
		fTrackAvgIonB           = TrackAvgIonB.getData();

		fVertIDNumber           = VertIDNumber.getData();
		fVertX                  = VertX.getData();
		fVertY                  = VertY.getData();
		fVertZ                  = VertZ.getData();
		fVertN                  = VertN.getData();
		fVertQ                  = VertQ.getData();

		fVT_VertIDNumber        = VT_VertIDNumber.getData();
		fVT_TrackIDNumber       = VT_TrackIDNumber.getData();
		fVT_TrackEnd            = VT_TrackEnd.getData();

		#ifdef writeECAL
		fClusterIDNumber        = ClusterIDNumber.getData();
		fClusterNhits           = ClusterNhits.getData();
		fClusterEnergy          = ClusterEnergy.getData();

		fECALAssn_ClusIDNumber  = ECALAssn_ClusIDNumber.getData();
		fECALAssn_TrackIDNumber = ECALAssn_TrackIDNumber.getData();
		fECALAssn_TrackEnd      = ECALAssn_TrackEnd.getData();
		#endif




		Float_t XvertMC, YvertMC, ZvertMC;
		int nMCVerts = MCVertX.size();
		for (int iMCVert=0; iMCVert<nMCVerts; ++iMCVert) {
			XvertMC = MCVertX.getData(iMCVert);
			YvertMC = MCVertY.getData(iMCVert);
			ZvertMC = MCVertZ.getData(iMCVert);
			if ( inFiducial(XvertMC,YvertMC,ZvertMC) ) {
				++Counter[nDatainFid];
				int interType = InterT.getData(iMCVert);
				if (interType==simb::kCCQE)						++Counter[nCCQEinFid];
				if (interType==simb::kNCQE)						++Counter[nNCQEinFid];
				if (interType==simb::kNuanceOffset)				++Counter[nRESinFid];
				if (simb::kResCCNuProtonPiPlus<=interType &&
					interType<=simb::kResCCNuBarProtonPi0Pi0)	++Counter[nRESinFid];
				if (interType==simb::kCCDIS)					++Counter[nCCDISinFid];
				if (interType==simb::kNCDIS)					++Counter[nNCDISinFid];
				if (interType==simb::kCCCOH)					++Counter[nCCCOHinFid];
				if (interType==simb::kNCCOH)					++Counter[nNCCOHinFid];
				if (interType==simb::kNuElectronElastic)		++Counter[nNuEinFid];				
			}
		}





		if (fVertIDNumber.size()==0) continue;

		// Setup iterators for the vector.erase function on the tracks
		vector<gar::rec::IDNumber>::iterator itTrackIDNumber       = fTrackIDNumber.begin();
		vector<Float_t>::iterator            itTrackStartX  	   = fTrackStartX.begin();
		vector<Float_t>::iterator            itTrackStartY  	   = fTrackStartY.begin();
		vector<Float_t>::iterator            itTrackStartZ  	   = fTrackStartZ.begin();
		vector<Float_t>::iterator            itTrackStartPX        = fTrackStartPX.begin();
		vector<Float_t>::iterator            itTrackStartPY        = fTrackStartPY.begin();
		vector<Float_t>::iterator            itTrackStartPZ        = fTrackStartPZ.begin();
		vector<Int_t>::iterator              itTrackStartQ         = fTrackStartQ.begin();
		vector<Float_t>::iterator            itTrackEndX		   = fTrackEndX.begin();
		vector<Float_t>::iterator            itTrackEndY		   = fTrackEndY.begin();
		vector<Float_t>::iterator            itTrackEndZ		   = fTrackEndZ.begin();
		vector<Float_t>::iterator            itTrackEndPX          = fTrackEndPX.begin();
		vector<Float_t>::iterator            itTrackEndPY          = fTrackEndPY.begin();
		vector<Float_t>::iterator            itTrackEndPZ          = fTrackEndPZ.begin();
		vector<Int_t>::iterator              itTrackEndQ           = fTrackEndQ.begin();
		vector<Float_t>::iterator            itTrackLenF           = fTrackLenF.begin();
		vector<Float_t>::iterator            itTrackLenB           = fTrackLenB.begin();
		vector<Int_t>::iterator              itNTPCClustersOnTrack = fNTPCClustersOnTrack.begin();
		vector<Float_t>::iterator            itTrackAvgIonF        = fTrackAvgIonF.begin();
		vector<Float_t>::iterator            itTrackAvgIonB        = fTrackAvgIonB.begin();

		// Step through the vectors, deleting as you go (this is slow, ok) be sure to keep
		// track of the reduced lengths
		int nTracks = fTrackIDNumber.size();
		for (int iTrack=0; iTrack<nTracks; ++iTrack) {
			Float_t Pstart = Qadd(fTrackStartPX[iTrack],fTrackStartPY[iTrack],fTrackStartPZ[iTrack]);
			Float_t Pend   = Qadd(fTrackEndPX[iTrack],  fTrackEndPY[iTrack],  fTrackEndPZ[iTrack]);
			if ( Pstart<PmagCut && Pend<PmagCut ) {
				itTrackIDNumber       = fTrackIDNumber.erase(itTrackIDNumber);
				itTrackStartX		  = fTrackStartX.erase(itTrackStartX);
				itTrackStartY		  = fTrackStartY.erase(itTrackStartY);
				itTrackStartZ		  = fTrackStartZ.erase(itTrackStartZ);
				itTrackStartPX        = fTrackStartPX.erase(itTrackStartPX);
				itTrackStartPY        = fTrackStartPY.erase(itTrackStartPY);
				itTrackStartPZ        = fTrackStartPZ.erase(itTrackStartPZ);
				itTrackStartQ         = fTrackStartQ.erase(itTrackStartQ);
				itTrackEndX 		  = fTrackEndX.erase(itTrackEndX);
				itTrackEndY 		  = fTrackEndY.erase(itTrackEndY);
				itTrackEndZ 		  = fTrackEndZ.erase(itTrackEndZ);
				itTrackEndPX          = fTrackEndPX.erase(itTrackEndPX);
				itTrackEndPY          = fTrackEndPY.erase(itTrackEndPY);
				itTrackEndPZ          = fTrackEndPZ.erase(itTrackEndPZ);
				itTrackEndQ           = fTrackEndQ.erase(itTrackEndQ);
				itTrackLenF           = fTrackLenF.erase(itTrackLenF);
				itTrackLenB           = fTrackLenB.erase(itTrackLenB);
				itNTPCClustersOnTrack = fNTPCClustersOnTrack.erase(itNTPCClustersOnTrack);
				itTrackAvgIonF        = fTrackAvgIonF.erase(itTrackAvgIonF);
				itTrackAvgIonB        = fTrackAvgIonB.erase(itTrackAvgIonB);
				--nTracks;	--iTrack;	// for (... iTrack<nTrack...) will break when iTrack==-1, nTrack==0
			} else {
				++itTrackIDNumber;
				++itTrackStartX;	++itTrackStartY;	++itTrackStartZ;
				++itTrackStartPX;	++itTrackStartPY;	++itTrackStartPZ;	++itTrackStartQ;
				++itTrackEndX;		++itTrackEndY;		++itTrackEndZ;
				++itTrackEndPX;		++itTrackEndPY;		++itTrackEndPZ;		++itTrackEndQ;
				++itTrackLenF;		++itTrackLenB;		++itNTPCClustersOnTrack;
				++itTrackAvgIonF;	++itTrackAvgIonB;
			}
		}
		if (fTrackIDNumber.size()==0) continue; 



		// Find the vertices with vert-track associations that are missing one or more tracks
		std::list<gar::rec::IDNumber> deadVertex;
		int nAssns = fVT_VertIDNumber.size();
		for (int iAssn=0; iAssn<nAssns; ++iAssn) {
			gar::rec::IDNumber thisIDNumber = fVT_TrackIDNumber[iAssn];
			vector<gar::rec::IDNumber>::iterator found;
			found = std::find(fTrackIDNumber.begin(),fTrackIDNumber.end(), thisIDNumber);
			if ( found == fTrackIDNumber.end() ) {
				// This vertex-track association is dead.  Save the dead vertex number
				deadVertex.push_back(fVT_VertIDNumber[iAssn]);
			}
		}
		if (deadVertex.size()>0) {
			// Get the list cleaned up
			std::list<gar::rec::IDNumber>::iterator deadStart = deadVertex.begin();
			std::list<gar::rec::IDNumber>::iterator deadEnd;
			deadVertex.sort();
			deadEnd = std::unique(deadVertex.begin(),deadVertex.end());



			// Delete the dead vertices
			vector<gar::rec::IDNumber>::iterator itVertIDNumber          = fVertIDNumber.begin();
			vector<Float_t>::iterator            itVertX                 = fVertX.begin();
			vector<Float_t>::iterator            itVertY                 = fVertY.begin();
			vector<Float_t>::iterator            itVertZ                 = fVertZ.begin();
			vector<Int_t>::iterator              itVertN                 = fVertN.begin();
			vector<Int_t>::iterator              itVertQ                 = fVertQ.begin();
			int nVerts = fVertIDNumber.size();
			for (int iVert=0; iVert<nVerts; ++iVert) {
				gar::rec::IDNumber thisIDNumber = fVertIDNumber[iVert];
				std::list<gar::rec::IDNumber>::iterator found;
				found = std::find(deadStart,deadEnd, thisIDNumber);
				if ( found != deadEnd ) {
					// Indeed, this is a dead vertex.  Delete it				
					itVertIDNumber        = fVertIDNumber.erase(itVertIDNumber);
					itVertX               = fVertX.erase(itVertX);
					itVertY               = fVertY.erase(itVertY);
					itVertZ               = fVertZ.erase(itVertZ);
					itVertN               = fVertN.erase(itVertN);
					itVertQ               = fVertQ.erase(itVertQ);
					--nVerts;	--iVert;
				} else {
					++itVertIDNumber;
					++itVertX;		++itVertY;		++itVertZ;		++itVertN;		++itVertQ;
				}
			}
			if (fVertIDNumber.size()==0) continue; 



			// Delete the dead vert-track associations
			vector<gar::rec::IDNumber>::iterator itVT_VertIDNumber       = fVT_VertIDNumber.begin();
			vector<gar::rec::IDNumber>::iterator itVT_TrackIDNumber      = fVT_TrackIDNumber.begin();
			vector<gar::rec::TrackEnd>::iterator itVT_TrackEnd           = fVT_TrackEnd.begin();
			for (int iAssn=0; iAssn<nAssns; ++iAssn) {
				gar::rec::IDNumber thisIDNumber = fVT_VertIDNumber[iAssn];
				std::list<gar::rec::IDNumber>::iterator found;
				found = std::find(deadStart,deadEnd, thisIDNumber);
				if ( found != deadEnd ) {
					// Indeed, this is a dead track-vertex association			
					itVT_VertIDNumber     = fVT_VertIDNumber.erase(itVT_VertIDNumber);
					itVT_TrackIDNumber    = fVT_TrackIDNumber.erase(itVT_TrackIDNumber);
					itVT_TrackEnd         = fVT_TrackEnd.erase(itVT_TrackEnd);
					--nAssns;	--iAssn;
				} else {
					++itVT_VertIDNumber;	++itVT_TrackIDNumber;		++itVT_TrackEnd;
				}
			}



			#ifdef writeECAL
			// Delete the dead ECAL-track associations, but not the ECAL cluster info
			vector<gar::rec::IDNumber>::iterator itECALAssn_ClusIDNumber  = fECALAssn_ClusIDNumber.begin();
			vector<gar::rec::IDNumber>::iterator itECALAssn_TrackIDNumber = fECALAssn_TrackIDNumber.begin();
			vector<gar::rec::TrackEnd>::iterator itECALAssn_TrackEnd      = fECALAssn_TrackEnd.begin();
			nAssns = fECALAssn_ClusIDNumber.size();
			for (int iAssn=0; iAssn<nAssns; ++iAssn) {
				gar::rec::IDNumber thisIDNumber = fECALAssn_TrackIDNumber[iAssn];
				vector<gar::rec::IDNumber>::iterator found;
				found = std::find(fTrackIDNumber.begin(),fTrackIDNumber.end(), thisIDNumber);
				if ( found == fTrackIDNumber.end() ) {
					// Indeed, this is a dead track-vertex association			
					itECALAssn_ClusIDNumber  = fECALAssn_ClusIDNumber.erase(itECALAssn_ClusIDNumber);
					itECALAssn_TrackIDNumber = fECALAssn_TrackIDNumber.erase(itECALAssn_TrackIDNumber);
					itECALAssn_TrackEnd      = fECALAssn_TrackEnd.erase(itECALAssn_TrackEnd);
					--nAssns;	--iAssn;
				} else {
					++itECALAssn_ClusIDNumber;	++itECALAssn_TrackIDNumber;		++itECALAssn_TrackEnd;
				}
			}
			#endif
		}




        outTree->Fill();
    }
    return;
}
