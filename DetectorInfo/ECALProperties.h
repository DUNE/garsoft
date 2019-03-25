////////////////////////////////////////////////////////////////////////
// \file ECALProperties.h
//
// \brief pure virtual base interface for ECAL properties
//
// \author eldwan.brianne@desy.de
//
////////////////////////////////////////////////////////////////////////
#ifndef DETINFO_ECALProperties_H
#define DETINFO_ECALProperties_H

// C/C++ standard libraries
#include <map>
#include <string>


///General GArSoft Utilities
namespace gar {
  namespace detinfo{

    class ECALProperties {
    public:

      ECALProperties(const ECALProperties &)              = delete;
      ECALProperties(ECALProperties &&)                   = delete;
      ECALProperties& operator = (const ECALProperties &) = delete;
      ECALProperties& operator = (ECALProperties &&)      = delete;
      virtual ~ECALProperties()                          = default;

      //Number of effective pixels on the SiPM
      virtual double EffectivePixel()        const = 0;

      //Light Yield of the ECAL tile (px/MIP)
      virtual double LightYield()        const = 0;

      //Gain of the ECAL SiPM (ADC/px)
      virtual double SiPMGain()        const = 0;

      //Birks coeficient for Scintillator (mm/MeV)
      virtual double ScintBirksConstant()        const = 0;

      //Intercalibration High/Low Gain factor
      virtual double IntercalibrationFactor()        const = 0;

      //ADC saturation value
      virtual double ADCSaturation()        const = 0;

      //Time resolution of the ECAL in ns
      virtual double TimeResolution()        const = 0;

      //MeV to MIP factor for 5 mm scintillator
      virtual double MeVtoMIP()        const = 0;

      //Noise in px
      virtual double NoisePx()        const = 0;

    protected:
      ECALProperties() = default;

    }; // class ECALProperties
  } //namespace detinfo
} // gar

#endif // DETINFO_ECALProperties_H
