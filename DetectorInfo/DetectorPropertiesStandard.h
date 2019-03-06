////////////////////////////////////////////////////////////////////////
// \file DetectorProperties.h
//
// \brief service to contain information about detector electronics, etc
//
// \author brebel@fnal.gov
//
// Separation of service from Detector info class:
// jpaley@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef DETINFO_DETECTORPROPERTIESSTD_H
#define DETINFO_DETECTORPROPERTIESSTD_H

// GArSoft libraries
#include "Geometry/GeometryCore.h"
#include "CoreUtils/ProviderPack.h"
#include "DetectorInfo/GArProperties.h"
#include "DetectorInfo/ECALProperties.h"
#include "DetectorInfo/DetectorClocks.h"
#include "DetectorInfo/DetectorProperties.h"

// framework libraries
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"

// C/C++ standard libraries
#include <set>

///General GArSoft Utilities
namespace gar {
  namespace detinfo{

    class DetectorPropertiesStandard : public DetectorProperties {
    public:
        /// List of service providers we depend on
      using providers_type = gar::ProviderPack<geo::GeometryCore,
                                               detinfo::GArProperties,
                                               detinfo::ECALProperties,
                                               detinfo::DetectorClocks>;

        /// Structure for configuration parameters
      struct Configuration_t {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;

        fhicl::Sequence<double> Efield{
          Name   ("Efield"),
          Comment("electric field in front of each wire plane (the last one is the big one!) [kV/cm]")
        };

        fhicl::Atom<double> Electronlifetime{
          Name   ("Electronlifetime"),
          Comment("electron lifetime in gaseous argon [us]")
        };
        fhicl::Atom<double> Temperature{
          Name   ("Temperature"),
          Comment("argon temperature [K]")
        };
        fhicl::Atom<double> ElectronsToADC{
          Name   ("ElectronsToADC"),
          Comment("conversion factor: (ADC counts)/(ionization electrons)")
        };
        fhicl::Atom<unsigned int> NumberTimeSamples{
          Name   ("NumberTimeSamples"),
          Comment("number of TPC readout TDC clock ticks per event (= readout window)")
        };

        fhicl::Atom<double> SternheimerA{
          Name   ("SternheimerA"),
          Comment("parameter a of Sternheimer correction delta = 2log(10) x - cbar + { a (x1-x)^k } theta(x1-x), x = log10(p/m)")
        };
        fhicl::Atom<double> SternheimerK{
          Name   ("SternheimerK"),
          Comment("parameter k of Sternheimer correction delta = 2log(10) x - cbar + { a (x_1-x)^k } theta(x1-x), x = log10(p/m)")
        };
        fhicl::Atom<double> SternheimerX0{
          Name   ("SternheimerX0"),
          Comment("minimum x = log10(p/m) for the application of Sternheimer correction")
        };
        fhicl::Atom<double> SternheimerX1{
          Name   ("SternheimerX1"),
          Comment("parameter x_1 of Sternheimer correction delta = 2log(10) x - cbar + { a (x_1-x)^k } theta(x1-x), x = log10(p/m)")
        };
        fhicl::Atom<double> SternheimerCbar{
          Name   ("SternheimerCbar"),
          Comment("parameter cbar of Sternheimer correction delta = 2log(10) x - cbar + { a (x_1-x)^k } theta(x1-x), x = log10(p/m)")
        };

      }; // Configuration_t

      DetectorPropertiesStandard();
      DetectorPropertiesStandard(fhicl::ParameterSet    const&  pset,
                                 const geo::GeometryCore     *  geo,
                                 const detinfo::GArProperties*  gp,
                                 const detinfo::ECALProperties* ecalp,
                                 const detinfo::DetectorClocks* c,
                                 std::set<std::string>   const& ignore_params = {});

      /**
       * @brief Constructs the provider and sets up the dependencies
       * @param pset FHiCL parameter set for provider configuration
       * @param providers pack of providers DetectorPropertiesStandard depends on
       * @see Setup()
       */
      DetectorPropertiesStandard(fhicl::ParameterSet const& pset,
                                 providers_type providers,
                                 std::set<std::string> const& ignore_params = {});
      DetectorPropertiesStandard(DetectorPropertiesStandard const&) = delete;
      virtual ~DetectorPropertiesStandard() = default;

      /**
       * @brief Configures the provider, first validating the configuration
       * @param p configuration parameter set
       * @param ignore_params parameters to be ignored (optional)
       *
       * This method will validate the parameter set (except for the parameters
       * it's explicitly told to ignore) and extract the useful information out
       * of it.
       */
      void ValidateAndConfigure(fhicl::ParameterSet   const& p,
                                std::set<std::string> const& ignore_params = {});


        /// Extracts the relevant configuration from the specified object
      void Configure(Configuration_t const& config);

      /**
       * @brief Validates the specified configuration
       * @param p configuration parameter set
       * @param ignore_params parameters to be ignored (optional)
       * @return a parsed configuration object
       * @see ValidateAndConfigure(), Configure()
       *
       * This method will validate the parameter set (except for the parameters
       * it's explicitly told to ignore) and it returns an object ready to
       * be used with Configure().
       */
      Configuration_t ValidateConfiguration(fhicl::ParameterSet   const& p,
                                            std::set<std::string> const& ignore_params = {});

      bool Update(uint64_t ts);
      bool UpdateClocks(const detinfo::DetectorClocks* clks);

      /**
       * @brief Sets all the providers at once
       * @param providers the pack of service providers we depend on
       *
       * Example:
       *
       *     gar::DetectorPropertiesStandard::providers_type providers;
       *     providers.set(gar::providerFrom<geo::Geometry>());
       *     providers.set(gar::providerFrom<detinfo::GArPropertiesService>());
       *     providers.set(gar::providerFrom<detinfo::DetectorClocksService>());
       *     detprop->Setup(providers);
       *
       */
      void Setup(providers_type providers);

      void SetGeometry      (const geo::GeometryCore* g)          { fGeo    = g;    }
      void SetGArProperties (const detinfo::GArProperties* gp)    { fGP     = gp;   }
      void SetECALProperties (const detinfo::ECALProperties* ecalp)    { fECALP     = ecalp;   }
      void SetDetectorClocks(const detinfo::DetectorClocks* clks) { fClocks = clks; }

      void SetNumberTimeSamples(unsigned int nsamp) { fNumberTimeSamples=nsamp;}
        // Accessors.

      virtual double Efield(unsigned int planegap=0) const override; ///< kV/cm

      virtual double DriftVelocity(double efield=0.,
                                   double temperature=0.,
                                   bool   cmPerns=true) const override;  ///< cm/ns if true, otherwise cm/us

      /// dQ/dX in electrons/cm, returns dE/dX in MeV/cm.
      virtual double BirksCorrection(double dQdX) const override;
      virtual double ModBoxCorrection(double dQdX) const override;

      virtual double ElectronLifetime()      const override { return fElectronlifetime;     }   //< microseconds

      /**
       * @brief Returns argon density at a given temperature
       * @param temperature the temperature in kelvin
       * @return argon density in g/cm^3
       *
       * Density is nearly a linear function of temperature.
       * See the NIST tables for details
       * Slope is between -6.2 and -6.1, intercept is 1928 kg/m^3.
       * This parameterization will be good to better than 0.5%.
       */
      virtual double Density(double temperature) const override;                          ///< g/cm^3

        // need to provide a definition, since the override above hides the inherited one
      virtual double Density() const override { return Density(Temperature()); }

        /// In kelvin.
      virtual double Temperature()                   const override { return fTemperature; }

      /**
       * @brief Restricted mean energy loss (dE/dx)
       * @param mom  momentum of incident particle [GeV/c]
       * @param mass mass of incident particle [GeV/c^2]
       * @param tcut maximum kinetic energy of delta rays [MeV]; 0 for unlimited
       * @return the restricted mean energy loss (dE/dx) in units of MeV/cm
       *
       * Returned value is always positive.
       * For unrestricted mean energy loss, set tcut = 0 (special case),
       * or tcut large.
       *
       * Based on Bethe-Bloch formula as contained in particle data book.
       * Material parameters are from the configuration.
       */
      virtual double Eloss(double mom, double mass, double tcut) const override;

      /**
       * @brief Energy loss fluctuation (@f$ \sigma_{E}^2 / x @f$)
       * @param mom  momentum of incident particle in [GeV/c]
       * @return energy loss fluctuation in MeV^2/cm
       *
       * Based on Bichsel formula referred to but not given in pdg.
       */
      virtual double ElossVar(double mom, double mass) const override;

      virtual double       SamplingRate()      const override { return fTPCClock.TickPeriod() * 1.e3; }
      virtual double       ElectronsToADC()    const override { return fElectronsToADC; }
      virtual unsigned int NumberTimeSamples() const override { return fNumberTimeSamples; }
      virtual int          TriggerOffset()     const override;

      virtual double       ConvertXToTicks(double X)     const override;
      virtual double       ConvertTicksToX(double ticks) const override;

      // The following methods convert between TDC counts (SimChannel time) and
      // ticks (RawDigit/Wire time).
      virtual double       ConvertTDCToTicks(double tdc) const override;
      virtual double       ConvertTicksToTDC(double ticks) const override;

      //ECAL Properties
      virtual double        EffectivePixel() const override { return fECALP->EffectivePixel(); }
      virtual double        LightYield() const override { return fECALP->LightYield(); }
      virtual double        SiPMGain() const override { return fECALP->SiPMGain(); }
      virtual double        IntercalibrationFactor() const override { return fECALP->IntercalibrationFactor(); }
      virtual double        ADCSaturation() const override { return fECALP->ADCSaturation(); }
      virtual double        TimeResolution() const override { return fECALP->TimeResolution(); }
      virtual double        MeVtoMIP() const override { return fECALP->MeVtoMIP(); }

      /// Verifies that the provider is in a fully configured status
      /// @throw cet::exception (category DetectorPropertiesStandard) if not ok
      void CheckIfConfigured() const;

    protected:


        /// Parameters for Sternheimer density effect corrections
      struct SternheimerParameters_t {
        double a;               ///< parameter a
        double k;               ///< parameter k
        double x0;              ///< parameter x0
        double x1;              ///< parameter x1
        double cbar;            ///< parameter Cbar
      }; //  SternheimerParameters_t

      void         CalculateXTicksParams();

      // service providers we depend on;
      // in principle could be replaced by a single providerpacl_type.
      const detinfo::GArProperties*  fGP;
      const detinfo::ECALProperties*  fECALP;
      const detinfo::DetectorClocks* fClocks;
      const geo::GeometryCore*       fGeo;

      std::vector< double >          fEfield;                ///< kV/cm (per inter-plane volume)
      double                         fElectronlifetime;      ///< microseconds
      double                         fTemperature;           ///< kelvin
      double                         fSamplingRate;          ///< in ns
      double 	                       fElectronsToADC;        ///< conversion factor for # of ionization electrons to 1 ADC count
      unsigned int                   fNumberTimeSamples;     ///< number of clock ticks per event (= readout window)

      SternheimerParameters_t        fSternheimerParameters; ///< Sternheimer parameters

      double                         fXTicksCoefficient;     ///< Parameters for x<-->ticks

      detinfo::ElecClock             fTPCClock;              ///< TPC electronics clock
    }; // class DetectorPropertiesStandard
  } //namespace detinfo
} // gar

#endif // DETINFO_DETECTOR_PROPERTIES_H
