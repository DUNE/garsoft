// A selection of includes & inlines for analysis



#ifndef includeROOT_h
#define includeROOT_h



#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <sys/stat.h>
#include <string>
using std::cout;		using std::endl;
using std::string;		using std::vector;



// Includes for common ROOT things
#include "Rtypes.h"         // Look in $ROOT_DIR/source/root-<version>/core/base/inc
#include "TMath.h"
#include "TVector3.h"
#include "TCut.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TFile.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TMinuit.h"



// Math made E-Z.  Too E-Z
#define BIG     1.79768e+308
#define BIGINT  2147483647
#define BIGSIZE 18446744073709551615

// Why isn't this in cmath?
#define M_PI 3.14159265358979323846

#define MIN(X,Y) ( (X)<(Y) ? (X) : (Y) )
#define MAX(X,Y) ( (X)>(Y) ? (X) : (Y) )
#define NINT(a)	 ( int(floor((a)+0.5)) )
#define SGN(a)   ( (a)<0 ? -1 : ( (a)>0 ? +1 : 0 ) )

inline int			SQR(int b)			{return b*b;};
inline float		SQR(float b)		{return b*b;}
inline double		SQR(double b)		{return b*b;}
inline long double	SQR(long double b)	{return b*b;}

inline int			CUBE(int b)			{return b*b*b;}
inline float		CUBE(float b)		{return b*b*b;}
inline double		CUBE(double b)		{return b*b*b;}
inline long double	CUBE(long double b)	{return b*b*b;}

inline int			QUAD(int b)			{int		 p = b*b; return p*p;}
inline float		QUAD(float b)		{float		 p = b*b; return p*p;}
inline double		QUAD(double b)		{double		 p = b*b; return p*p;}
inline long double	QUAD(long double b)	{long double p = b*b; return p*p;}



// std::hypot appeared in C++17.
inline double Qadd(double const x, double const y){
	return std::hypot(x,y);
}
inline double Qadd(double const x, double const y, double const z){
    double absX,absY,absZ;
    absX = fabs(x);     absY = fabs(y);     absZ = fabs(z);
    if ( absX>absY && absX>absZ ) return std::hypot(x,std::hypot(y,z));
    if ( absY>absX && absY>absZ ) return std::hypot(y,std::hypot(x,z));
    return Qadd(z,Qadd(x,y));
}
inline double Qsub(double const x, double const y) {
	if (std::fabs(x) < std::fabs(y)) {
		std::cout << "Fatality in Qsub(double,double): Invalid arguments";
		std::cout << std::endl;
		std::exit(1);
	}
	long double x2 = x*x;		long double y2 = y*y;
	return double(std::sqrt( x2 - y2 )); 
}



// Primitive earth algorithm for uncertainty in binomial-distributed variables
inline double binerr(double const passed, double const tested) {
    // Garbage input
    if ((passed<0)   || (tested<=0) || (passed>tested)) {
		std::cout << "Fatality in binerr(double,double): Invalid arguments ";
		std::cout << passed << "/" << tested << std::endl;
		std::exit(1);
	}
	double err;
	// Fairly ad-hoc statistical uncertainties for low stats
	if ( tested==1.0 ) {
		err = 0.393;
	} else if ( passed==0.0 || passed==tested ) {
		err = 0.421 / pow(tested,0.94);
	} else {
		err = std::sqrt( passed*(tested-passed)/tested );
	}
	err /= tested;
    return err;
}



// A string for std::cout that is a ratio and its binomial uncertainty
char* stringEffErr(int num, int den) {
	// Yea, its a memory leak.  Fight me.
	char* temp = new char[32];
	double rat = 1000.0*double(num)/double(den);
	double err = 1000.0*binerr(double(num),double(den));
// Assuming UTF-8 encoding on your terminal or whatever
	sprintf(temp,"%5.2f %c%c %5.2f m%c%c", rat, 0xC2,0xB1, err, 0xC2,0xB0);
	return temp;
};



// Yea does that file even exist?  Probably doesn't work for xrood access.
inline bool filehamna(const std::string& filename) {
    struct stat buf;
    int retval = stat(filename.c_str(), &buf);
    if (retval==0 && !S_ISREG(buf.st_mode)) {
        cout << "filehamna(string) : Not a regular file" << endl;
    }
    // -1 is the error condition, meaning that the file isn't there
    return (retval == -1);
}



int const nothingPDG   =    0;
int const protonPDG    = 2212;
int const neutronPDG   = 2112;
int const electronPDG  =   11;	// Means an electron not a positron
int const muonPDG      =   13;
int const nuePDG       =   12;
int const numuPDG      =   14;
int const pichPDG      =  211;	// Means a pi+; pi- is -211
int const pi0PDG       =  111;
int const KchPDG       =  321;	// Means a K+; K- is -321
int const KlongPDG     =  130;
int const KshortPDG    =  310;
int const KzeroPDG     =  311;
int const gammaPDG     =   22;

int const ArgonPDG     = 1000180400;
int const CarbonPDG    = 1000060120;
int const DeuteronPDG  = 1000010020;
int const HeliumPDG    = 1000020040;

int const lambdaPDG    = 3122;
int const sigmaM_PDG   = 3112;
int const sigmaP_PDG   = 3222;
int const cascadeM_PDG = 3312;
int const cascade0_PDG = 3322;



double const Melectron = 0.00051100;	// in GeV
double const Mmuon     = 0.10565837;
double const Mpi0      = 0.13497770;
double const MpiCH     = 0.13957061;
double const MchK      = 0.493677;
double const MneuK     = 0.497611;
double const Mprot     = 0.93827208;
double const Mneut     = 0.93956541;
double const Mlambda   = 1.115683;
double const MsigmaNeg = 1.197449;
double const MsigmaPos = 1.18937;
double const Mcasc0    = 1.31486;
double const McascNeg  = 1.32171;



// Mass of this particle, in GeV, from PDG code.  Negative value
// means unstable particle!
double massPDG(int PDGcode) {

	if (PDGcode > 1000000000) {
		int A = (PDGcode -1000000000)/10;
		A = A%1000;
		return A * (39.948 / 40.0);
	}

	switch (PDGcode){
		case +electronPDG:	case -electronPDG:
			return Melectron;
		case +muonPDG:		case -muonPDG:
			return Mmuon;
		case  gammaPDG:
		case +nuePDG:		case -nuePDG:
		case +numuPDG:		case -numuPDG:
			return 0.000000;
		case +pichPDG:		case -pichPDG:
			return Mpi0;
		case +KchPDG:		case -KchPDG:
			return MpiCH;
		case +protonPDG:	case -protonPDG:
			return Mprot;
		case  pi0PDG:
			return Mpi0;
		case  KlongPDG:
			return MneuK;
		case +neutronPDG:	case -neutronPDG:
			return Mneut;

		case  KshortPDG:
		case +KzeroPDG:		case -KzeroPDG:
			return -MneuK;
		case +lambdaPDG:	case -lambdaPDG:
			return -Mlambda;
		case +sigmaM_PDG:	case -sigmaM_PDG:
			return -MsigmaNeg;
		case +sigmaP_PDG:	case -sigmaP_PDG:
			return -MsigmaPos;
		case +cascade0_PDG:
			return -Mcasc0;
		case +cascadeM_PDG:	case -cascadeM_PDG:
			return -McascNeg;

		default:
			// Should be some wide, fast-decaying resonance.
			return -2468;
	}
};



// Is this particle charged, from the PDG code
bool chargedPDG(int PDGcode) {

	if (PDGcode > 1000000000) return true;

	switch (PDGcode){
		case +electronPDG:	case -electronPDG:
		case +muonPDG:		case -muonPDG:
		case +pichPDG:		case -pichPDG:
		case +KchPDG:		case -KchPDG:
		case +protonPDG:	case -protonPDG:
		case +sigmaM_PDG:	case -sigmaM_PDG:
		case +sigmaP_PDG:	case -sigmaP_PDG:
		case +cascadeM_PDG:	case -cascadeM_PDG:
			return true;

		case  gammaPDG:
		case +nuePDG:		case -nuePDG:
		case +numuPDG:		case -numuPDG:
		case  pi0PDG:
		case  KlongPDG:		case  KshortPDG:
		case +KzeroPDG:		case -KzeroPDG:
		case +neutronPDG:	case -neutronPDG:
		case +lambdaPDG:	case -lambdaPDG:
		case +cascade0_PDG:
			return false;

		default:
			// Should be some wide, fast-decaying resonance.
			return false;
	}
};



#endif
