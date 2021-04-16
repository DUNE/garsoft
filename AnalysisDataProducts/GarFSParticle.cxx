#include "garsoft/AnalysisDataProducts/GarFSParticle.h"

using namespace gar;

GarFSParticle::GarFSParticle(simb::MCParticle mcp) {
   fFSP = new garana::FSParticle(mcp.TrackId(), mcp.PdgCode(),mcp.Position(),mcp.Momentum()); 
}
