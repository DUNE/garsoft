////////////////////////////////////////////////////////////////////////
// $Id: InfoTransfer.h,v 1.15 2010/08/30 21:33:24 andrzejs Exp $
//
// Transfer hitlist and run info into a producer module.
// Do not copy this code without contacting Andrzej Szelc and Brian Rebel first.
//
// \author andrzej.szelc@yale.edu
// \author ellen.klein@yale.edu
////////////////////////////////////////////////////////////////////////
#ifndef INFOTRANSFER_H
#define INFOTRANSFER_H
#ifndef __CINT__
#include <iostream>
#include <string>
#include <vector>

#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

#include "ReconstructionDataProducts/Cluster.h"
#include "ReconstructionDataProducts/Hit.h"

namespace util {
  class PxPoint;
  class PxLine;
}

namespace evd {
  class InfoTransfer {
  public:
    InfoTransfer(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
    ~InfoTransfer();

    // The Rebuild function rebuilds the various maps we need to pickup hits.
    // It is called automatically before each event is processed. For jobs involving
    // Monte Carlo generation, this is too soon. So, we'll call rebuild after those data
    // products are put into the event in LArG4.  This is the least bad way of ensuring the
    // InfoTransfer works in jobs that combine MC production and reconstruction analysis based
    // on MC truth.  Don't copy this design pattern without talking to brebel@fnal.gov first
    void Rebuild(const art::Event& evt);

    void reconfigure(fhicl::ParameterSet const& pset);
    void SetTestFlag(int value) { testflag = value; };
    int GetTestFlag() { return testflag; };
    void SetRunNumber(int value) { fRun = value; };
    int GetRunNumber() { return fRun; };
    void SetSubRunNumber(int value) { fSubRun = value; };
    int GetSubRunNumber() { return fSubRun; };
    void SetEvtNumber(int value) { fEvt = value; };
    int GetEvtNumber() { return fEvt; };

    void SetHitList(unsigned int p, std::vector<art::Ptr<recob::Hit>> hits_to_save)
    {
      fSelectedHitlist[p].clear();
      fSelectedHitlist[p] = hits_to_save;
    };

    std::vector<art::Ptr<recob::Hit>> GetHitList(unsigned int plane)
    {
      return fRefinedHitlist[plane];
    };

    std::vector<art::Ptr<recob::Hit>> GetSelectedHitList(unsigned int plane)
    {
      return fSelectedHitlist[plane];
    };

    void ClearSelectedHitList(int plane)
    {
      if (fSelectedHitlist.size() == 0) {
        return;
        std::cout << "no size" << std::endl;
      }
      fSelectedHitlist[plane].clear();
      for (unsigned int i = 0; i < fRefStartHit.size(); i++) {
        fRefStartHit[i] = NULL;
        fRefEndHit[i] = NULL;
      }
      return;
    }

    void SetStartHit(unsigned int p, recob::Hit* starthit) { fStartHit[p] = starthit; };

    recob::Hit* GetStartHit(unsigned int plane) { return fRefStartHit[plane]; };

    void SetEndHit(unsigned int p, recob::Hit* endhit) { fEndHit[p] = endhit; };

    recob::Hit* GetEndHit(unsigned int plane) { return fRefEndHit[plane]; };

    std::vector<double> GetStartHitCoords(unsigned int plane) { return refstarthitout[plane]; };
    std::vector<double> GetEndHitCoords(unsigned int plane) { return refendhitout[plane]; };

    void SetStartHitCoords(unsigned int plane, std::vector<double> starthitin)
    {
      starthitout[plane].clear();
      starthitout[plane].resize(2);
      starthitout[plane] = starthitin;
    };
    void SetEndHitCoords(unsigned int plane, std::vector<double> endhitin)
    {
      endhitout[plane].clear();
      endhitout[plane].resize(2);
      endhitout[plane] = endhitin;
    };

    void SetSeedList(std::vector<util::PxLine> seedlines);

    std::vector<util::PxLine> GetSeedList();

  private:
    void FillStartEndHitCoords(unsigned int plane);

    int testflag;
    int fEvt;
    int fRun;
    int fSubRun;
    std::vector<std::vector<art::Ptr<recob::Hit>>>
      fSelectedHitlist; ///< the list selected by the GUI (one for each plane)
    std::vector<std::vector<art::Ptr<recob::Hit>>>
      fRefinedHitlist; ///< the refined hitlist after rebuild (one for each plane)
    std::vector<art::Ptr<recob::Hit>> fFullHitlist; ///< the full Hit list from the Hitfinder.
    std::string fHitModuleLabel;                    ///< label for geant4 module

    std::vector<recob::Hit*> fStartHit;    ///< The Starthit
    std::vector<recob::Hit*> fRefStartHit; ///< The Refined Starthit

    std::vector<recob::Hit*> fEndHit;    ///< The Starthit
    std::vector<recob::Hit*> fRefEndHit; ///< The Refined Starthit

    std::vector<util::PxLine> fSeedList;

    std::vector<std::vector<double>> starthitout;
    std::vector<std::vector<double>> endhitout;

    std::vector<std::vector<double>> refstarthitout;
    std::vector<std::vector<double>> refendhitout;
  };
} //namespace
#endif // __CINT__
DECLARE_ART_SERVICE(evd::InfoTransfer, LEGACY)
#endif
