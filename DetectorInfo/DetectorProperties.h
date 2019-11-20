////////////////////////////////////////////////////////////////////////
// \file DetectorProperties.h
//
// \brief pure virtual base interface for detector properties
//
// \author jpaley@fnal.gov
//
////////////////////////////////////////////////////////////////////////
#ifndef DETINFO_IDETECTORPROPERTIES_H
#define DETINFO_IDETECTORPROPERTIES_H

#include "fhiclcpp/ParameterSet.h"
#include "Geometry/Geometry.h"

///General GArSoft Utilities
namespace gar {
  namespace detinfo{

    // Conversion for energy deposited in GeV to number of ionization electrons produced
    constexpr double kGeVToElectrons = 3.788e7;  ///< 26.4 eV per ion pair, 1e9 eV/GeV

    class DetectorProperties {
    public:

      DetectorProperties(const DetectorProperties &) = delete;
      DetectorProperties(DetectorProperties &&) = delete;
      DetectorProperties& operator = (const DetectorProperties &) = delete;
      DetectorProperties& operator = (DetectorProperties &&) = delete;
      virtual ~DetectorProperties() = default;

      /**
       * @brief Returns the nominal electric field in the specified volume
       * @param planegap volume specification (default: 0, the big drift volume)
       * @return electric field in the volume, in kV/cm
       *
       * The electric field is "nominal", i.e., a completely uniform field is
       * assumed.
       *
       * The planegap argument identifies which volume to return the field value
       * for. The relation between planegap and readout plane is not perfectly
       * formalized yet. In general, a good rule is that planegap N describes
       * the volume on the cathode side of wire plane N. This rule is formally
       * valid also for ArgoNeuT/LArIAT, where three wire planes are present.
       * But only two of them are instrumented and read, that are called
       * "readout plane 0" and "readout plane 1", but effectively correspond to
       * planegap 1 and 2.
       *
       * Note that all TPCs are assumed to have the same electric field values.
       */
      virtual double Efield(unsigned int planegap=0) const = 0;

      virtual double DriftVelocity(double efield=0.,
                                   double temperature=0.,
                                   bool   cmPerns=true) const = 0;

      virtual double ElectronLifetime() const = 0;

      /**
       * @brief Returns argon density at a given temperature
       * @param temperature the temperature in kelvin
       * @return argon density in g/cm^3
       */
      virtual double Density(double temperature) const = 0;
      virtual double Temperature() const = 0;

      /**
       * @brief Restricted mean energy loss (@f$ dE/dx @f$)
       * @param mom  momentum of incident particle [GeV/c]
       * @param mass mass of incident particle [GeV/c^2]
       * @param tcut maximum kinetic energy of delta rays [MeV]; 0 for unlimited
       * @return the restricted mean energy loss (dE/dx) in units of MeV/cm
       *
       * Returned value is always positive.
       * For unrestricted mean energy loss, set tcut = 0 (special case),
       * or tcut large.
       */
      virtual double Eloss(double mom,
                           double mass,
                           double tcut)                   const = 0;

      /**
       * @brief Energy loss fluctuation (@f$ \sigma_{E}^2 / x @f$)
       * @param mom  momentum of incident particle in [GeV/c]
       * @return energy loss fluctuation in MeV^2/cm
       */
      virtual double ElossVar(double mom,
                              double mass)                 const = 0;

        /// Returns argon density at the temperature from Temperature()
      virtual double Density() const { return Density(Temperature()); }

      virtual double       SamplingRate()                  const = 0;
      virtual double       ElectronsToADC()                const = 0;
      virtual unsigned int NumberTimeSamples()             const = 0;
      virtual int          TriggerOffset()                 const = 0;

      virtual double       ConvertXToTicks(double X)       const = 0;
      virtual double       ConvertTicksToX(double ticks)   const = 0;

      // The following methods convert between TDC counts (SimChannel time) and
      // ticks (RawDigit time).
      virtual double       ConvertTDCToTicks(double tdc)   const = 0;
      virtual double       ConvertTicksToTDC(double ticks) const = 0;

      //ECAL Properties
      virtual double        EffectivePixel() const = 0;
      virtual double        LightYield() const = 0;
      virtual double        SiPMGain() const = 0;
      virtual double        IntercalibrationFactor() const = 0;
      virtual double        ADCSaturation() const = 0;
      virtual double        TimeResolution() const = 0;
      virtual double        MeVtoMIP() const = 0;
      virtual double        NoisePx() const = 0;

    protected:
      DetectorProperties() = default;

    }; // class DetectorProperties
  } //namespace detinfo
}
#endif // DETINFO_IDETECTORPROPERTIES_H
