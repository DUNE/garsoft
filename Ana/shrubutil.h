#ifndef GARSOFT_SHRUBUTIL_H
#define GARSOFT_SHRUBUTIL_H


////////////////////////////////////////////////////////////////////////
/// \file  shrubutil.h
/// \brief Code to compute quantities that will go into the anatree, shrub version
///
//
// The functions are:
//      float gar::processIonizationInfo( rec::TrackIoniz& ion, float ionizeTruncate )
//      float gar::computeT( simb::MCTruth theMCTruth )
//
/// \version $Id:  $
/// \author  bellanto@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "ReconstructionDataProducts/Hit.h"
#include "ReconstructionDataProducts/Track.h"
#include "ReconstructionDataProducts/TrackIoniz.h"
#include "ReconstructionDataProducts/Vertex.h"



namespace gar {



    // Convert signal values from ADCs into calibrated dE/dX values
    void processIonizationInfo_shrub( rec::TrackIoniz& ion, float ionizeTruncate,
                                float& forwardIonVal, float& backwardIonVal );
    float processOneDirection_shrub(std::vector<std::pair<float,float>> SigData, float ionizeTruncate);
    // Helper method for processOneDirection
    bool lessThan_byE_shrub(std::pair<float,float> a, std::pair<float,float> b)
        {return a.first < b.first;}



    // Compute T for coherent pion analysis
    float computeT_shrub( simb::MCTruth theMCTruth );



}



#endif
