#ifndef GAR_ANAUTILS_H
#define GAR_ANAUTILS_H

#include <string>
#include <unordered_map>
#include <vector>

#include <TLorentzVector.h>

#include "garana/DataProducts/GTruth.h"
#include "nusimdata/SimulationBase/GTruth.h"

#include "garana/DataProducts/FSParticle.h"
#include "garana/DataProducts/G4Particle.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include "ReconstructionDataProducts/Track.h"
#include "garana/DataProducts/Track.h"

#include "ReconstructionDataProducts/Cluster.h"
#include "garana/DataProducts/CaloCluster.h"

#include "ReconstructionDataProducts/Vertex.h"
#include "garana/DataProducts/Vertex.h"

#include "ReconstructionDataProducts/Vee.h"
#include "garana/DataProducts/Vee.h"

using std::pair;
using std::vector;

namespace gar {

  int ProcessNameToCode(std::string const& p);
  void FillGTruth(const simb::GTruth& gt, garana::GTruth& outtruth);

  garana::GTruth MakeAnaGTruth(const simb::GTruth& gt, const int& vtxregion);
  garana::FSParticle MakeFSParticle(const simb::MCParticle& mcp);
  garana::G4Particle MakeG4Particle(const simb::MCParticle& mcp,
                                    int parentPdg,
                                    int progenitorPdg,
                                    int progenitorTrackId,
                                    const vector<pair<TLorentzVector, TLorentzVector>>& positions,
                                    const vector<pair<TLorentzVector, TLorentzVector>>& momenta,
                                    const vector<int>& regions,
                                    const vector<size_t>& nptsPerRegion);

  const garana::Track MakeAnaTrack(const rec::Track& trk,
                                   const vector<pair<int, float>>& pidf,
                                   const vector<pair<int, float>>& pidb,
                                   float ionf,
                                   float ionb,
                                   const vector<pair<UInt_t, TLorentzVector>>& posBeg,
                                   const vector<pair<UInt_t, TLorentzVector>>& posEnd,
                                   const vector<pair<UInt_t, TLorentzVector>>& momBeg,
                                   const vector<pair<UInt_t, TLorentzVector>>& momEnd,
                                   const vector<pair<int, float>>& edeps);

  garana::CaloCluster MakeAnaCalCluster(const rec::Cluster& clust,
                                        const int& region,
                                        const vector<pair<int, float>>& edeps);
  garana::Vee MakeAnaVee(const rec::Vee& vee);
  garana::Vertex MakeAnaVtx(const rec::Vertex& vtx);

  /*void ChangeToTpcCoords(garana::GTruth& gt,         const TLorentzVector& origin);
    void ChangeToTpcCoords(garana::Track& track,       const TLorentzVector& origin);
    void ChangeToTpcCoords(garana::CaloCluster& clust, const TLorentzVector& origin);
    void ChangeToTpcCoords(garana::Vertex& vtx,        const TLorentzVector& origin);
    void ChangeToTpcCoords(garana::Vee& vee,           const TLorentzVector& origin);
    void ChangeToTpcCoords(garana::G4Particle& g4p,    const TLorentzVector& origin);
    void ChangeToTpcCoords(garana::FSParticle& g4p,    const TLorentzVector& origin);*/
}

#endif
