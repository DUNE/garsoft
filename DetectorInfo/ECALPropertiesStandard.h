/**
 * @file    ECALPropertiesStandard.h
 * @brief   Service provider with utility ECAL functions
 * @see     ECALPropertiesStandard.cxx LArPropertiesStandardTestHelpers.h
 *
 * The provider detinfo::LArProperiesStandard supports simple setup for testing
 * environment, by including in your test:
 *
 *     #include "DetectorInfo/LArPropertiesStandardTestHelpers.h"
 *
 */

#ifndef DETECTORINFO_ECALPROPERTIESSTANDARD_H
#define DETECTORINFO_ECALPROPERTIESSTANDARD_H


// GArSoft libraries
#include "DetectorInfo/ECALProperties.h"

// FHiCL libraries
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"

// C/C++ standard libraries
#include <string>
#include <vector>
#include <map>
#include <set>


namespace gar {
  namespace detinfo {
    /**
     * @brief Properties related to liquid argon environment in the detector
     *
     * This class can access databases via DatabaseUtil service.
     *
     * @note Some of the database connection properties are established before
     * the beginning of the job and if they change this service will not be
     * aware of it. These properties petrain, so far, only the connection mode
     * and not any content of the databases themselves.
     * @note 2: the database connection features for this base class have been removed
     */
    class ECALPropertiesStandard : public ECALProperties {
    public:
      ECALPropertiesStandard();
      explicit ECALPropertiesStandard(fhicl::ParameterSet const& pset,
                                     std::set<std::string>      ignore_params = {});
      ECALPropertiesStandard(ECALPropertiesStandard const&) = delete;
      virtual ~ECALPropertiesStandard() = default;

      /**
       * @brief Configures the provider
       * @param p configuration parameter set
       * @param ignore_params parameters to be ignored (optional)
       *
       * This method will validate the parameter set (except for the parameters
       * it's explicitly told to ignore) and extract the useful information out
       * of it.
       */
      bool   Configure(fhicl::ParameterSet const& pset,
                       std::set<std::string>      ignore_params = {});
      bool   Update(uint64_t ts=0);

      /// SiPM Number of effective pixels (px)
      virtual double EffectivePixel()       const override { return fNeffPx; } ///< g/cm^2

      /// Light yield of the tile (px/MIP)
      virtual double LightYield()          const override { return fLY; }

      /// SiPM Gain (ADC/px)
      virtual double SiPMGain()            const override { return fGain; }

      /// Birks constant (mm/MeV)
      virtual double ScintBirksConstant()            const override { return fBirks; }

      void SetEffectivePixel      (double effpx) { fNeffPx  = effpx; }
      void SetLightYield         (double ly ) { fLY = ly;                 }
      void SetSiPMGain           (double gain ) { fGain = gain;                 }
      void SetScintBirksConstant (double birks) { fBirks = birks; }

    private:
    protected:

        /// structure with all configuration parameters
      struct Configuration_t {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;

        fhicl::Atom<double> EffectivePixel     { Name("EffectivePixel" ),     Comment("SiPM Number of effective pixels (px)")     };
        fhicl::Atom<double> LightYield         { Name("LightYield"    ),      Comment("Light yield of the tile (px/MIP)")         };
        fhicl::Atom<double> SiPMGain           { Name("SiPMGain"      ),      Comment("SiPM Gain (ADC/px)")                       };
        fhicl::Atom<double> ScintBirksConstant { Name("ScintBirksConstant"),  Comment("ScintBirksConstant (mm/MeV)")              };

      }; // Configuration_t


      bool   fIsConfigured;

      double fNeffPx;  ///< px
      double fLY;                ///< px/MIP
      double fGain;                ///< ADC/px
      double fBirks;  ///< mm/MeV

    public:
      // expose the configuration object for framework service
      using ConfigurationParameters_t = Configuration_t;

    }; // class ECALPropertiesStandard
  } //namespace detinfo
} // gar
#endif // DETECTORINFO_ECALPROPERTIESSTANDARD_H
