#ifndef GARSOFT_ANAUTIL_H
#define GARSOFT_ANAUTIL_H


////////////////////////////////////////////////////////////////////////
/// \file  anautil.h
/// \brief Code to compute quantities that will go into the anatree
///
/// \version $Id:  $
/// \author  bellanto@fnal.gov
////////////////////////////////////////////////////////////////////////


#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "ReconstructionDataProducts/Hit.h"
#include "ReconstructionDataProducts/Track.h"
#include "ReconstructionDataProducts/TrackIoniz.h"
#include "ReconstructionDataProducts/Vertex.h"



namespace gar {

    // Converte signal values from ADCs into calibrated dE/dX values
    void processIonizationInfo( rec::TrackIoniz& ion );

    // Average dE/dX for a track
    float AverageIonization( rec::TrackIoniz ion );

    // Compute T for coherent pion analysis
    float computeT( simb::MCTruth theMCTruth );

}



#endif
