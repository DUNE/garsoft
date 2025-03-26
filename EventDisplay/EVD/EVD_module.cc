////////////////////////////////////////////////////////////////////////
/// \file EVD_plugin.cc
///
/// \version $Id: EVD_plugin.cc,v 1.1 2010/11/10 22:49:25 p-novaart Exp $
/// \author  jpaley@anl.gov

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"

#ifndef EVD_EVD_H
#define EVD_EVD_H

// Framework Includes
#include "art/Framework/Core/EDAnalyzer.h"

#include "TH1D.h"
#include <string>

extern "C" {
#include <sys/stat.h>
#include <sys/types.h>
}

//GArSoft includes
#include "EventDisplay/EVD/CalorView.h"
#include "EventDisplay/EVD/Display3DView.h"
#include "nuevdb/EventDisplayBase/DisplayWindow.h"

// Framework includes
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

/// The Event Display
// Builder for the Display3D view
static evdb::Canvas* mk_display3d_canvas(TGMainFrame* mf)
{
  return new gar::evd::Display3DView(mf);
}

// Builder for the Calor view
static evdb::Canvas* mk_calor_canvas(TGMainFrame* mf)
{
  return new gar::evd::CalorView(mf);
}

// // Builder for the MCTruth view
// static evdb::ObjListCanvas* mk_mctrue_canvas(TGMainFrame* mf)
// {
//   return new evd::MCTrueView(mf);
// }

namespace gar {
  namespace evd {

    /// a class for transporting photons in a roughly realistic way
    class EVD : public art::EDAnalyzer {
    public:
      explicit EVD(fhicl::ParameterSet const& pset);
      virtual ~EVD();

      void analyze(art::Event const& evt);
      void beginJob();

    private:
      //unused bool fWindowsDrawn; ///< flag for whether windows are already drawn
    };

    //----------------------------------------------------
    EVD::EVD(fhicl::ParameterSet const& pset) : EDAnalyzer(pset)
    //, fWindowsDrawn(false)
    {}

    //----------------------------------------------------
    EVD::~EVD() {}

    //----------------------------------------------------
    void EVD::beginJob()
    {
      evdb::DisplayWindow::Register("Display3D", "Display3D", 700, 700, mk_display3d_canvas);

      evdb::DisplayWindow::Register("Calorimetry", "Calorimetry", 700, 700, mk_calor_canvas);

      //     evdb::ListWindow::Register("MC Particle List",
      // 			       "MC Particle List",
      // 			       400,
      // 			       800,
      // 			       mk_mctrue_canvas);

      // Open up the main display window and run
      evdb::DisplayWindow::OpenWindow(0);
    }

    //----------------------------------------------------
    void EVD::analyze(const art::Event& /*evt*/) {}

    DEFINE_ART_MODULE(EVD)

  } // namespace evd
} // namespace gar

#endif // EVD
////////////////////////////////////////////////////////////////////////
