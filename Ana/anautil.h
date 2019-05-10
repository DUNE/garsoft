#ifndef GARSOFT_ANAUTIL_H
#define GARSOFT_ANAUTIL_H


////////////////////////////////////////////////////////////////////////
/// \file  anautil.h
/// \brief Code to compute quantities that will go into the anatree
///
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
  float processIonizationInfo( rec::TrackIoniz& ion, float ionizeTruncate );

  // Compute T for coherent pion analysis
  float computeT( simb::MCTruth theMCTruth );
}



#endif
