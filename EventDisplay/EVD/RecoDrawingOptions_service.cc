////////////////////////////////////////////////////////////////////////
/// \file RecoDrawingOption_plugin.cc
///
/// \version $Id: RecoDrawingOptions_plugin.cc,v 1.1 2010/11/11 18:11:22 p-novaart Exp $
/// \author  brebel@fnal.gov

// Framework includes

/// GArSoft includes
#include "EventDisplay/EVD/RecoDrawingOptions.h"

#include <iostream>

namespace gar {
namespace evd {

//......................................................................
RecoDrawingOptions::RecoDrawingOptions(fhicl::ParameterSet const& pset,
                                       art::ActivityRegistry& /* reg */)
  : evdb::Reconfigurable{pset}
{
    this->reconfigure(pset);
}

//......................................................................
RecoDrawingOptions::~RecoDrawingOptions()
{
}

//......................................................................
void RecoDrawingOptions::reconfigure(fhicl::ParameterSet const& pset)
{
    fDrawHits                  = pset.get< int                      >("DrawHits"                 );
    fDrawTPCClusters           = pset.get< int                      >("DrawTPCClusters"          );
    fTPCClusterMarker          = pset.get< int                      >("TPCClusterMarker"         );
    fTPCClusterMarkerSize      = pset.get< int                      >("TPCClusterMarkerSize"     );
    fDrawCaloClusters          = pset.get< int                      >("DrawCaloClusters"         );
    fDrawCaloHits              = pset.get< int                      >("DrawCaloHits"             );
    fDrawTracks      	       = pset.get< int                      >("DrawTracks"     	         );
    fTrackWidth      	       = pset.get< int                      >("TrackWidth"     	         );
    fDrawTrackTrajectoryPoints = pset.get< int                      >("DrawTrackTrajectoryPoints");
    fDrawShowers     	       = pset.get< int                      >("DrawShowers"    	         );
    fDrawVecHits    	       = pset.get< int                      >("DrawVecHits"   	     	 );
    fDrawVertices    	       = pset.get< int                      >("DrawVertices"   	     	 );
    fDrawEvents      	       = pset.get< int                      >("DrawEvents"     	     	 );
    fHitLabels                 = pset.get< std::vector<std::string> >("HitModuleLabels"          );
    fTPCClusterLabels          = pset.get< std::vector<std::string> >("TPCClusterModuleLabels"   );
    fCaloClusterLabels         = pset.get< std::vector<std::string> >("CaloClusterModuleLabels"  );
    fCaloHitLabels             = pset.get< std::vector<std::string> >("CaloHitModuleLabels"      );
    fTrackLabels      	       = pset.get< std::vector<std::string> >("TrackModuleLabels"     	 );
    fShowerLabels     	       = pset.get< std::vector<std::string> >("ShowerModuleLabels"    	 );
    fVecHitLabels     	       = pset.get< std::vector<std::string> >("VecHitModuleLabels"    	 );
    fVertexLabels     	       = pset.get< std::vector<std::string> >("VertexModuleLabels"    	 );
    fEventLabels      	       = pset.get< std::vector<std::string> >("EventModuleLabels"     	 );
    fCaloClusterScale          = pset.get< float                    >("CaloClusterScale"         );
    fCaloHitScale              = pset.get< float                    >("CaloHitScale"             );
  }

}

namespace evd {

  DEFINE_ART_SERVICE(RecoDrawingOptions)

} // namespace evd
}
////////////////////////////////////////////////////////////////////////
