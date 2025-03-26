//
//    AnaUtils.h
//
//    Created by C. Hilgenberg

//#ifndef GAR_ANAUTILS_cxx
//#define GAR_ANAUTILS_cxx

#include "AnaUtils.h"
#include <iostream>

using std::pair;
using std::vector;

namespace gar {

  garana::GTruth MakeAnaGTruth(const simb::GTruth& gt, const int& vtxregion)
  {

    garana::GTruth outtruth;
    FillGTruth(gt, outtruth);

    outtruth.fVertexRegion = vtxregion;

    return outtruth;
  }

  garana::FSParticle MakeFSParticle(const simb::MCParticle& mcp)
  {
    return garana::FSParticle(mcp.TrackId(), mcp.PdgCode(), mcp.Position(), mcp.Momentum());
  }

  garana::G4Particle MakeG4Particle(const simb::MCParticle& mcp,
                                    int parentPdg,
                                    int progenitorPdg,
                                    int progenitorTrackId,
                                    const vector<pair<TLorentzVector, TLorentzVector>>& positions,
                                    const vector<pair<TLorentzVector, TLorentzVector>>& momenta,
                                    const vector<int>& regions,
                                    const vector<size_t>& nptsPerRegion)
  {
    return garana::G4Particle(mcp.NumberTrajectoryPoints(),
                              mcp.PdgCode(),
                              parentPdg,
                              progenitorPdg,
                              mcp.TrackId(),
                              mcp.Mother(),
                              progenitorTrackId,
                              gar::ProcessNameToCode(mcp.Process()),
                              gar::ProcessNameToCode(mcp.EndProcess()),
                              positions,
                              momenta,
                              regions,
                              nptsPerRegion);
  }

  //////////////////////////////////////////////
  garana::CaloCluster MakeAnaCalCluster(const rec::Cluster& clust,
                                        const int& region,
                                        const std::vector<std::pair<int, float>>& edeps)
  {

    TLorentzVector pos(clust.Position()[0], clust.Position()[1], clust.Position()[2], clust.Time());

    //trade in our arrays for vectors/TVector3s
    const float* eigs = clust.EigenVectors();
    std::vector<TVector3> vecs;
    for (int i = 0; i < 3; i++) {
      TVector3 eig(eigs[i], eigs[i + 1], eigs[i + 2]);
      vecs.push_back(eig);
    }

    garana::CaloCluster outclust(pos,
                                 region,
                                 clust.Energy(),
                                 clust.EnergyError(),
                                 clust.TimeDiffFirstLast(),
                                 clust.Shape(),
                                 clust.ITheta(),
                                 clust.IPhi(),
                                 vecs,
                                 edeps);

    return outclust;
  } //

  ///////////////////////////////////////////////
  garana::Vee MakeAnaVee(const rec::Vee& vee)
  {

    TLorentzVector vtx(vee.Position()[0], vee.Position()[1], vee.Position()[2], vee.Time());
    std::vector<TLorentzVector> vecs;
    for (int i = 0; i < 3; i++) {
      vecs.push_back(vee.FourMomentum(i));
    }

    garana::Vee outvee(vtx, vecs, vee.Chisq(), vee.VertexCov());

    return outvee;
  } //
  ///////////////////////////////////////////////
  garana::Vertex MakeAnaVtx(const rec::Vertex& vtx)
  {

    TLorentzVector pos(vtx.Position()[0], vtx.Position()[1], vtx.Position()[2], vtx.Time());
    garana::Vertex outvtx(pos, vtx.CovMat());

    return outvtx;
  } //

  ////////////////////////////////////////////////
  std::unordered_map<std::string, int> processMap{{"unknown", 0},
                                                  {"primary", 1},
                                                  {"compt", 2},
                                                  {"phot", 3},
                                                  {"annihil", 4},
                                                  {"eIoni", 5},
                                                  {"eBrem", 6},
                                                  {"conv", 7},
                                                  {"muIoni", 8},
                                                  {"muMinusCaptureAtRest", 9},
                                                  {"neutronInelastic", 10},
                                                  {"nCapture", 11},
                                                  {"hadElastic", 12},
                                                  {"Decay", 13},
                                                  {"CoulombScat", 14},
                                                  {"muPairProd", 15},
                                                  {"muBrems", 16},
                                                  {"muPairProd", 17},
                                                  {"PhotonInelastic", 18},
                                                  {"hIoni", 19},
                                                  {"protonInelastic", 20},
                                                  {"pi+Inelastic", 21},
                                                  {"CHIPSNuclearCaptureAtRest", 22},
                                                  {"pi-Inelastic", 23},
                                                  {"hBertiniCaptureAtRest", 24},
                                                  {"photonNuclear", 25},
                                                  {"muonNuclear", 26},
                                                  {"kaon-Inelastic", 27},
                                                  {"kaon+Inelastic", 28},
                                                  {"Transportation", 29},
                                                  {"nKiller", 30},
                                                  {"ionIoni", 31},
                                                  {"msc", 32},
                                                  {"dInelastic", 33},
                                                  {"kaon0LInelastic", 34},
                                                  {"electronNuclear", 35},
                                                  {"tInelastic", 36},
                                                  {"kaon0SInelastic", 37},
                                                  {"sigma+Inelastic", 38},
                                                  {"lambdaInelastic", 39},
                                                  {"anti_protonInelastic", 40},
                                                  {"alphaInelastic", 41}};

  //convert process name from G4 to numeric value for
  // easy storage in TTrees
  int ProcessNameToCode(std::string const& p)
  {
    if (processMap.find(p) == processMap.end()) {
      std::cout << "WARNING(AnaUtils::ProcessNameToCode): "
                << "unknown process name, '" << p << "'" << std::endl;
      return -9999;
    }

    else
      return processMap[p];
  }

  // make a copy of the ghep record and write to empty garana product
  void FillGTruth(const simb::GTruth& gt, garana::GTruth& outtruth)
  {

    outtruth.fVertex = gt.fVertex;
    outtruth.fVertex.SetVect(
      100. * gt.fVertex.Vect()); // GENIE uses meters - convert to cm (same as MCParticle)
    outtruth.fVertex.SetT(
      1.e9 * gt.fVertex.T()); // GENIE uses seconds - convert to ns (same as MCParticle)
    outtruth.fweight = gt.fweight;
    outtruth.fprobability = gt.fprobability;
    outtruth.fXsec = gt.fXsec;
    outtruth.fDiffXsec = gt.fDiffXsec;
    outtruth.fGPhaseSpace = gt.fGPhaseSpace;

    outtruth.fProbePDG = gt.fProbePDG;
    outtruth.fProbeP4 = gt.fProbeP4;
    outtruth.fTgtP4 = gt.fTgtP4;

    outtruth.ftgtZ = gt.ftgtZ;
    outtruth.ftgtA = gt.ftgtA;
    outtruth.ftgtPDG = gt.ftgtPDG;
    outtruth.fHitNucPDG = gt.fHitNucPDG;
    outtruth.fHitQrkPDG = gt.fHitQrkPDG;
    outtruth.fIsSeaQuark = gt.fIsSeaQuark;
    outtruth.fHitNucP4 = gt.fHitNucP4;
    outtruth.fHitNucPos = gt.fHitNucPos;

    outtruth.fGscatter = gt.fGscatter;
    outtruth.fGint = gt.fGint;

    outtruth.fgQ2 = gt.fgQ2;
    outtruth.fgq2 = gt.fgq2;
    outtruth.fgW = gt.fgW;
    outtruth.fgT = gt.fgT;
    outtruth.fgX = gt.fgX;
    outtruth.fgY = gt.fgY;
    outtruth.fFSleptonP4 = gt.fFSleptonP4;
    outtruth.fFShadSystP4 = gt.fFShadSystP4;

    outtruth.fIsCharm = gt.fIsCharm;
    outtruth.fCharmHadronPdg = gt.fCharmHadronPdg;
    outtruth.fIsStrange = gt.fIsStrange;
    outtruth.fStrangeHadronPdg = gt.fStrangeHadronPdg;
    outtruth.fNumProton = gt.fNumProton;
    outtruth.fNumNeutron = gt.fNumNeutron;
    outtruth.fNumPi0 = gt.fNumPi0;
    outtruth.fNumPiPlus = gt.fNumPiPlus;
    outtruth.fNumPiMinus = gt.fNumPiMinus;
    outtruth.fResNum = gt.fResNum;
    outtruth.fDecayMode = gt.fDecayMode;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  const garana::Track MakeAnaTrack(const rec::Track& trk,
                                   const vector<pair<int, float>>& pidf,
                                   const vector<pair<int, float>>& pidb,
                                   float ionf,
                                   float ionb,
                                   const vector<pair<UInt_t, TLorentzVector>>& posBeg,
                                   const vector<pair<UInt_t, TLorentzVector>>& posEnd,
                                   const vector<pair<UInt_t, TLorentzVector>>& momBeg,
                                   const vector<pair<UInt_t, TLorentzVector>>& momEnd,
                                   const vector<pair<int, float>>& edeps)
  {

    const float* tmp = trk.Vertex();
    TLorentzVector vtx(tmp[0], tmp[1], tmp[2], trk.Time());

    tmp = trk.End();
    TLorentzVector end(
      tmp[0], tmp[1], tmp[2], trk.Time()); // eventually track end times should be different

    tmp = trk.VtxDir();
    TVector3 vtxDir(tmp[0], tmp[1], tmp[2]);

    tmp = trk.EndDir();
    TVector3 endDir(tmp[0], tmp[1], tmp[2]);

    garana::Track anatrk(trk.LengthForward(),
                         trk.LengthBackward(),
                         trk.Momentum_beg(),
                         trk.Momentum_end(),
                         vtx,
                         end,
                         vtxDir,
                         endDir,
                         trk.ChisqForward(),
                         trk.ChisqBackward(),
                         trk.NHits(),
                         trk.TrackParBeg(),
                         trk.TrackParEnd(),
                         trk.CovMatBegPacked(),
                         trk.CovMatEndPacked(),
                         trk.ChargeBeg(),
                         trk.ChargeEnd(),
                         pidf,
                         pidb,
                         ionf,
                         ionb,
                         posBeg,
                         posEnd,
                         momBeg,
                         momEnd,
                         edeps);

    return anatrk;
  }

  /*void ChangeToTpcCoords(garana::GTruth& gt,         const TLorentzVector& origin) {
        gt
    }

    void ChangeToTpcCoords(garana::Track& track,       const TLorentzVector& origin) {

    }

    void ChangeToTpcCoords(garana::CaloCluster& clust, const TLorentzVector& origin) {

    }

    void ChangeToTpcCoords(garana::Vertex& vtx,        const TLorentzVector& origin) {

    }

    void ChangeToTpcCoords(garana::Vee& vee,           const TLorentzVector& origin) {

    }

    void ChangeToTpcCoords(garana::G4Particle& g4p,    const TLorentzVector& origin) {

    }

    void ChangeToTpcCoords(garana::FSParticle& g4p,    const TLorentzVector& origin) {

    }*/

} //namespace

//#endif
