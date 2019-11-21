#ifndef CAF_h
#define CAF_h

#include "TFile.h"
#include "TTree.h"

#include "Utils.h"

class CAF {

    typedef std::pair<float, std::string> P;

public:

    /* Default Constructor */
    CAF();

    /* Constructor */
    CAF( std::string infile, std::string filename );

    /* Default Copy Constructor */
    CAF(const CAF &) = default;

    /* Destructor */
    ~CAF();

    /* Method to book the output TFile */
    bool BookTFile();

    /* Method to fill the output TTree */
    void FillTTree();

    /* Method to write the output TTree */
    void WriteTTree();

    /* Method to close the output TTree */
    void CloseTFile();

    /* Method to clear the event variables */
    void ClearVectors();

    /* Method to loop */
    void loop();

private:
    std::string _outputFile; ///< The output TFile name */
    TFile * cafFile; ///< The output TFile pointer */
    TTree * cafMVA; ///< The output TTree pointer */

    std::string _inputfile;
    TFile* _intfile;
    TTree* _inttree;

    Utils *_util;

    //Event-wise values
    int _Run, _Event, _SubRun;
    //Generator values
    std::vector<int> mode, ccnc, ntype, gint, weight, tgtpdg, gt_t, intert;
    std::vector<double> q2, w, y, x, theta, t, mctime, mcnupx, mcnupy, mcnupz, vertx, verty, vertz;
    //MC Particle Values
    std::vector<int> _nFSP;
    std::vector<int> mctrkid, mother, pdgmother, truepdg, _MCPStartX, _MCPStartY, _MCPStartZ, _MCPEndX, _MCPEndY, _MCPEndZ;
    std::vector<std::string> _MCProc, _MCEndProc;
    std::vector<float> trkLen, trkLenPerp, truep, truepx, truepy, truepz, _angle;
    //Reco values
    std::vector<int> recopid;
    std::vector<float> prob_arr, partereco, anglereco, _preco, erecon;
};

#endif
