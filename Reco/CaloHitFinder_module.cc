////////////////////////////////////////////////////////////////////////
// Class:       CaloHitFinder
// Plugin Type: producer (art v2_11_02)
// File:        CompressedHitFinder_module.cc
//
// Generated at Wed Jun 13 15:27:13 2018 by Thomas Junk using cetskelgen
// from cetlib version v3_03_01.
// Report the compressed raw digit blocks as hits -- just a threshold
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "RawDataProducts/CaloRawDigit.h"
#include "RawDataProducts/raw.h"
#include "ReconstructionDataProducts/CaloHit.h"
#include "Geometry/Geometry.h"
#include "DetectorInfo/DetectorClocksService.h"
#include "DetectorInfo/DetectorPropertiesService.h"

#include "ReadoutSimulation/ECALUtils.h"

#include <memory>

namespace gar {
  namespace rec {

    class CaloHitFinder : public art::EDProducer {
    public:
      explicit CaloHitFinder(fhicl::ParameterSet const & p);
      // The compiler-generated destructor is fine for non-base
      // classes without bare pointers or other resource use.

      // Plugins should not be copied or assigned.
      CaloHitFinder(CaloHitFinder const &) = delete;
      CaloHitFinder(CaloHitFinder &&) = delete;
      CaloHitFinder & operator = (CaloHitFinder const &) = delete;
      CaloHitFinder & operator = (CaloHitFinder &&) = delete;

      // Required functions.
      void produce(art::Event & e) override;
      float CalibrateToMIP(unsigned int ADC);
      float CalibratetoMeV(double MIP);

    private:

      // Declare member data here.

      float fMIPThreshold;   ///< zero-suppression threshold (in case the raw digits need to be zero-suppressed)
      int   fClusterHits;    ///< hit clustering algorithm number

      std::string fRawDigitLabel;  ///< label to find the right raw digits
      const detinfo::DetectorProperties*  fDetProp;      ///< detector properties
      const geo::GeometryCore*            fGeo;          ///< pointer to the geometry
      std::unique_ptr<rosim::ECALUtils>          fECALUtils;    ///< pointer to the Util fcn for the ECAL
    };


    CaloHitFinder::CaloHitFinder(fhicl::ParameterSet const & p)
    // :
    {
      fMIPThreshold = p.get<float>("MIPThreshold", 0.25);
      fRawDigitLabel = p.get<std::string>("RawDigitLabel", "daqecal");
      fClusterHits  = p.get<int>("ClusterHits", 0);

      fGeo     = gar::providerFrom<geo::Geometry>();
      fDetProp = gar::providerFrom<detinfo::DetectorPropertiesService>();
      fECALUtils = std::make_unique<rosim::ECALUtils>(fDetProp->EffectivePixel(), 0.95);

      produces< std::vector<rec::CaloHit> >();
    }

    void CaloHitFinder::produce(art::Event & e)
    {
      // create an emtpy output hit collection -- to add to.
      std::unique_ptr<std::vector<CaloHit> > hitCol (new std::vector<CaloHit> );
      std::unique_ptr<std::vector<CaloHit> > hitClusterCol (new std::vector<CaloHit> );

      // the input raw digits
      auto rdCol = e.getValidHandle< std::vector<raw::CaloRawDigit> >(fRawDigitLabel);

      for (size_t ird = 0; ird < rdCol->size(); ++ ird)
      {
        auto const& rd = (*rdCol)[ird];
        unsigned int hitADC = rd.ADC();
        float hitTime = rd.Time();
        float x = rd.X();
        float y = rd.Y();
        float z = rd.Z();
        unsigned int id = rd.CaloID();

        //Do Calibration of the hit in MIPs
        float hitMIP = this->CalibrateToMIP(hitADC);

        if(hitMIP < fMIPThreshold)
        {
          LOG_DEBUG("CaloHitFinder") << "Signal under the " << fMIPThreshold << " MIP threshold" << std::endl;
          continue;
        }

        double sat_px = hitMIP * fDetProp->LightYield();
        //DeSaturate
        double unsat_px = fECALUtils->DeSaturate(sat_px);

        //Calibrate to the MeV scale
        double unsat_energy = unsat_px / fDetProp->LightYield();
        float energy = this->CalibratetoMeV(unsat_energy);
        float pos[3] = {x, y, z};

        hitCol->emplace_back(energy, hitTime, pos, id);
      }

      // cluster hits if requested

      if (fClusterHits == 0)
      {
        e.put(std::move(hitCol));
        return;
      }
      else{
        LOG_DEBUG("CaloHitFinder") << "Clustering of calo hit not done yet!" << std::endl;
        return;
      }
    }

    float CaloHitFinder::CalibrateToMIP(unsigned int ADC)
    {
      if(ADC <= 0) return ADC;

      return ADC / (fDetProp->LightYield() * fDetProp->SiPMGain());
    }

    float CaloHitFinder::CalibratetoMeV(double MIP)
    {
      if(MIP <= 0) return MIP;

      return MIP * 0.814;//MeV
    }

    DEFINE_ART_MODULE(CaloHitFinder)

  } // namespace rec
} // namespace gar
