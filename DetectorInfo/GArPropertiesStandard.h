/**
 * @file    LArPropertiesStandard.h
 * @brief   Service provider with utility LAr functions
 * @see     LArPropertiesStandard.cxx LArPropertiesStandardTestHelpers.h
 *
 * The provider detinfo::LArProperiesStandard supports simple setup for testing
 * environment, by including in your test:
 *
 *     #include "DetectorInfo/LArPropertiesStandardTestHelpers.h"
 *
 */

#ifndef DETECTORINFO_GARPROPERTIESSTANDARD_H
#define DETECTORINFO_GARPROPERTIESSTANDARD_H

// GArSoft libraries
#include "DetectorInfo/GArProperties.h"

// FHiCL libraries
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"

// C/C++ standard libraries
#include <map>
#include <set>
#include <string>
#include <vector>

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
    class GArPropertiesStandard : public GArProperties {
    public:
      GArPropertiesStandard();
      explicit GArPropertiesStandard(fhicl::ParameterSet const& pset,
                                     std::set<std::string> ignore_params = {});
      GArPropertiesStandard(GArPropertiesStandard const&) = delete;
      virtual ~GArPropertiesStandard() = default;

      /**
       * @brief Configures the provider
       * @param p configuration parameter set
       * @param ignore_params parameters to be ignored (optional)
       *
       * This method will validate the parameter set (except for the parameters
       * it's explicitly told to ignore) and extract the useful information out
       * of it.
       */
      bool Configure(fhicl::ParameterSet const& pset, std::set<std::string> ignore_params = {});
      bool Update(uint64_t ts = 0);

      virtual double RadiationLength() const override { return fRadiationLength; } ///< g/cm^2

      /// Ar atomic number
      virtual double AtomicNumber() const override { return fZ; }

      /// Ar atomic mass (g/mol)
      virtual double AtomicMass() const override { return fA; }

      /// Ar mean excitation energy (eV)
      virtual double ExcitationEnergy() const override { return fI; }

      /// Fano factor
      virtual double FanoFactor() const override { return fF; }

      /// Diffusion constants
      virtual double LongitudinalDiffusion() const override { return fLongDiff; }
      virtual double TransverseDiffusion() const override { return fTransDiff; }

      void SetRadiationLength(double rl) { fRadiationLength = rl; }
      void SetAtomicNumber(double z) { fZ = z; }
      void SetAtomicMass(double a) { fA = a; }
      void SetMeanExcitationEnergy(double e) { fI = e; }
      void SetFanoFactor(double f) { fF = f; }
      void SetLongitudinalDiffusion(double ld) { fLongDiff = ld; }
      void SetTransverseDiffusion(double td) { fTransDiff = td; }

    private:
    protected:
      /// structure with all configuration parameters
      struct Configuration_t {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;

        fhicl::Atom<double> RadiationLength{Name("RadiationLength"),
                                            Comment("radiation length [g/cm^2]")};
        fhicl::Atom<double> AtomicNumber{Name("AtomicNumber"),
                                         Comment("atomic number (yes, yes, it's 18...)")};
        fhicl::Atom<double> AtomicMass{Name("AtomicMass"), Comment("atomic mass [g/mol]")};
        fhicl::Atom<double> MeanExcitationEnergy{Name("ExcitationEnergy"),
                                                 Comment("mean excitation energy [eV]")};
        fhicl::Atom<double> FanoFactor{Name("FanoFactor"),
                                       Comment("controls fluctuations in ionization")};
        fhicl::Atom<double> LongitudinalDiffusion{Name("LongitudinalDiffusion"),
                                                  Comment("controls fluctuations in ionization")};
        fhicl::Atom<double> TransverseDiffusion{Name("TransverseDiffusion"),
                                                Comment("controls fluctuations in ionization")};

      }; // Configuration_t

      bool fIsConfigured;

      double fRadiationLength; ///< g/cm^2
      double fZ;               ///< Ar atomic number
      double fA;               ///< Ar atomic mass (g/mol)
      double fI;               ///< Ar mean excitation energy (eV)
      double fF;               ///< Ar Fano factor
      double fLongDiff;        ///< longitudinal diffusion in cm^2/ns
      double fTransDiff;       ///< transverse diffusion in cm^2/ns

    public:
      // expose the configuration object for framework service
      using ConfigurationParameters_t = Configuration_t;

    }; // class GArPropertiesStandard
  }    //namespace detinfo
} // gar
#endif // GARPROPERTIES_H
