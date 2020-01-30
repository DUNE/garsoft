// Templates to simplify reading TTrees.  Constructor-getter pattern, modest
// efficiency considerations considered.
//
// The use of try-catch structures is for decorative purposes only; root v6
// does not hurl, it segfaults.  See https://xkcd.com/371/
// The signature for TTree::SetBranchAddress in the vectorFromTree construtor is
// template <class T> Int_t SetBranchAddress(const char*, T**, TBranch** = 0)
// which means that you could call with a std::deque or perhaps some other
// container in the 2nd argument rather than a std::vector.  And it would still
// compile.  But then getData() actually will getGarbage().  Because root is so
// great.
//
// Honestly, I wouldn't try resizing the vector returned from getData() directly.
// maybe try a deep copy first or something.  I suspect this only works because
// vectors are guaranteed to be contiguous storage, unlike other containers.
//
// Leo Bellantoni, 2019 copyright FRA.
//
// Use like so:
//
//		Long64_t iEntry = 0;
//		scalarFromTree<int> readRun(Oak,"Run",&iEntry);
//		vectorFromTree<int> readVertN(Oak,"VertN",&iEntry);
//		vectorFromTree<int> readVertQ(Oak,"VertQ",&iEntry);
//		Int_t Run    = readRun();	// cuz' you initialized iEntry to 0
//		Int_t nEntry = Oak->GetEntries();
//		for (iEntry=0; iEntry<nEntry; ++iEntry) {
//			vector<int> VertN = readVertN();
//			size_t nVertInEvent = VertN.size();
//			vector<int> TwoTrackVerts = readVertN.findIndices(2);
//			for (auto aTwoTrackVertIndex : TwoTrackVerts) {
//				int Qvert = readVertQ.getData(aTwoTrackVertIndex);
//				... whatever
//			}
//		}



#ifndef treeReader_h
#define treeReader_h



#include <string>
using std::string;
#include <vector>
using std::vector;
//#include <algorithm>

#include "Rtypes.h"
#include "TTree.h"
#include "TLeaf.h"





// Read a scalar ===============================================================
// =============================================================================
template<typename T> class scalarFromTree {
public:



	scalarFromTree(TTree* whichTree, string varname, Long64_t* iEntry_in) :
		thisTree(whichTree),thisVarname(varname),iEntry_local(iEntry_in) {
		lastEntry = -1;
	}

	// Compiler won't make a default constructor.  Gotta do it here.
	scalarFromTree() {
		thisTree = NULL;		iEntry_local = NULL;	thisVarname = "";
		lastEntry = -2;
	}



	// Get that scalar!  Get it!   Here is the preferred signature
	T operator()() {return getData();};

	// This signature deprecated but is nicer than (*scalarFromTree)(i)
	T getData() {
		if (*iEntry_local != lastEntry) {
			lastEntry = *iEntry_local;
			try {
				thisTree->GetEntry(*iEntry_local);
				sData = thisTree->GetLeaf(thisVarname.c_str())->GetValue(0);
			} catch (...) {
				throw std::runtime_error("Can not GetLeaf in scalarFromTree.getData()");
			}
		}
		return sData;
	}



private:
	TTree* thisTree;
	string thisVarname;
	Long64_t* iEntry_local;
	T sData;
	Long64_t lastEntry;
};





// Read a std::vector ==========================================================
// =============================================================================
template<typename T> class vectorFromTree {
public:



	// Constructor only calls SetBranchAddress =================================
	vectorFromTree(TTree* whichTree, string varname, Long64_t* iEntry_in) :
		thisTree(whichTree),thisVarname(varname),iEntry_local(iEntry_in),
		vData(NULL),bData(NULL) {
		// EQUIVALENCE was a bad idea in FORTRAN IV, never mind C++14.
		try {
			thisTree->SetBranchAddress(thisVarname.c_str(),&vData,&bData);		
		} catch (...) {
			throw std::runtime_error("SetBranchAddress fails in vectorFromTree constructor");
		}
		lastEntry = -1;
	}

	// Compiler won't make a default constructor.  Gotta do it here. ===========
	vectorFromTree() {
		thisTree = NULL;	iEntry_local = NULL;
		vData = NULL;		bData = NULL;
		thisVarname = "";
		lastEntry = -2;
	}



	// Get that vector!  Get it! ===============================================
	vector<T>& getDataVector() {
		if (*iEntry_local != lastEntry) {
			lastEntry = *iEntry_local;
			try {
				Long64_t tEntry = thisTree->LoadTree(*iEntry_local);
				bData->GetEntry(tEntry);
			} catch (...) {
				throw std::runtime_error("Can not GetEntry in vectorFromTree.getData()");
			}
		}
		return (*vData);
	}



	// Get something from that vector!  Get that! ==============================
	// This is the preferred signature
	T operator()(int i) {
		return getDataVector().at(i);
	}
	// This signature deprecated but is nicer than (*vectorFromTree)(i)
	T getData(int i) {
		return getDataVector().at(i);
	}



	// Get that vector's size!  Get him! =======================================
	int size() {
		return getDataVector().size();
	}



	// Construct a vector with all indices where data matches searchval ========
	vector<int> findIndices(T searchval) {
		vector<int> retval;
		int nData = getDataVector().size();
		for (int iDatum=0; iDatum<nData; ++iDatum) {
			if ( getDataVector().at(iDatum)==searchval ) retval.push_back(iDatum);
		}
		return retval;
	}



private:
	TTree* thisTree;
	string thisVarname;
	Long64_t* iEntry_local;
	vector<T>* vData;
	TBranch*   bData;
	Long64_t lastEntry;
};



#endif
