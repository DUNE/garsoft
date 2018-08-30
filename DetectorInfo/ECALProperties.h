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

    protected:
      ECALProperties() = default;

    }; // class ECALProperties
  } //namespace detinfo
} // gar

#endif // DETINFO_ECALProperties_H
