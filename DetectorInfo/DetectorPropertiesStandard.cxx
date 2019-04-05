////////////////////////////////////////////////////////////////////////
//
//  \file DetectorProperties.cxx
//
// Separation of service from Detector info class:
// jpaley@fnal.gov
//
////////////////////////////////////////////////////////////////////////
// Framework includes

#include <cassert>

// GArSoft includes
#include "DetectorInfo/DetectorPropertiesStandard.h"
#include "CoreUtils/ProviderUtil.h" // gar::IgnorableProviderConfigKeys()
#include "CoreUtils/DebugUtils.h"   // gar::IgnorableProviderConfigKeys()
#include "Geometry/Geometry.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

// Art includes
#include "fhiclcpp/make_ParameterSet.h"

namespace {

  template <typename T>
  inline T sqr(T v) { return v*v; }

} // local namespace

namespace gar {
  namespace detinfo{

    //--------------------------------------------------------------------
    DetectorPropertiesStandard::DetectorPropertiesStandard()
    : fGP(0)
    , fECALP(0)
    , fClocks(0)
    , fGeo(0)
    {

    }

    //--------------------------------------------------------------------
    DetectorPropertiesStandard::DetectorPropertiesStandard(fhicl::ParameterSet      const& pset,
                                                           const geo::GeometryCore*        geo,
                                                           const detinfo::GArProperties*   gp,
                                                           const detinfo::ECALProperties*  ecalp,
                                                           const detinfo::DetectorClocks*  c,
                                                           std::set<std::string> const&    ignore_params /* = {} */)
    : fGP(gp)
    , fECALP(ecalp)
    , fClocks(c)
    , fGeo(geo)
    {
      {
        mf::LogInfo debug("setupProvider<DetectorPropertiesStandard>");

        debug << "Asked to ignore " << ignore_params.size() << " keys:";
        for (auto const& key: ignore_params) debug << " '" << key << "'";
      }

      ValidateAndConfigure(pset, ignore_params);

      fTPCClock = fClocks->TPCClock();

      LOG_WARNING("DetectorPropertiesStandard")
      << "DetectorPropertiesStandard Warning!"
      << "\t The following functions need to be verified for gaseous argon:"
      << "\t\t BirksCorrection"
      << "\t\t Density"
      << "\t\t DriftVelocity"
      << "\t\t EField";


    }

    //--------------------------------------------------------------------
    DetectorPropertiesStandard::DetectorPropertiesStandard(fhicl::ParameterSet const& pset,
                                                           providers_type             providers,
                                                           std::set<std::string> const& ignore_params /* = {} */)
    :DetectorPropertiesStandard(pset,
                                providers.get<geo::GeometryCore>(),
                                providers.get<detinfo::GArProperties>(),
                                providers.get<detinfo::ECALProperties>(),
                                providers.get<detinfo::DetectorClocks>(),
                                ignore_params)
    {}

    //--------------------------------------------------------------------
    bool DetectorPropertiesStandard::Update(uint64_t t)
    {

      CalculateXTicksParams();
      return true;
    }

    //--------------------------------------------------------------------
    bool DetectorPropertiesStandard::UpdateClocks(const detinfo::DetectorClocks* clks)
    {
      fClocks = clks;

      fTPCClock = fClocks->TPCClock();
      CalculateXTicksParams();
      return true;
    }

    //------------------------------------------------------------
    double DetectorPropertiesStandard::ConvertTDCToTicks(double tdc) const
    {
      return fClocks->TPCTDC2Tick(tdc);
    }

    //--------------------------------------------------------------
    double DetectorPropertiesStandard::ConvertTicksToTDC(double ticks) const
    {
      return fClocks->TPCTick2TDC(ticks);
    }


    //--------------------------------------------------------------------
    void DetectorPropertiesStandard::Configure(Configuration_t const& config) {

      fEfield                     = config.Efield();
      fElectronlifetime           = config.Electronlifetime();
      fTemperature                = config.Temperature();
      fDriftVelocity              = config.DriftVelocity();
      fElectronsToADC             = config.ElectronsToADC();
      fNumberTimeSamples          = config.NumberTimeSamples();

      fSternheimerParameters.a    = config.SternheimerA();
      fSternheimerParameters.k    = config.SternheimerK();
      fSternheimerParameters.x0   = config.SternheimerX0();
      fSternheimerParameters.x1   = config.SternheimerX1();
      fSternheimerParameters.cbar = config.SternheimerCbar();

      CalculateXTicksParams();

    } // DetectorPropertiesStandard::Configure()

    //--------------------------------------------------------------------
    DetectorPropertiesStandard::Configuration_t
    DetectorPropertiesStandard::ValidateConfiguration(fhicl::ParameterSet   const& p,
                                                      std::set<std::string> const& ignore_params /* = {} */)
    {
      std::set<std::string> ignorable_keys = gar::IgnorableProviderConfigKeys();
      ignorable_keys.insert(ignore_params.begin(), ignore_params.end());

      gar::debug::printBacktrace(mf::LogDebug("DetectorPropertiesStandard"), 15);

        // parses and validates the parameter set:
      fhicl::Table<Configuration_t> config_table { p, ignorable_keys };

      return std::move(config_table());

    } // DetectorPropertiesStandard::ValidateConfiguration()

    //--------------------------------------------------------------------
    void DetectorPropertiesStandard::ValidateAndConfigure(fhicl::ParameterSet    const& p,
                                                          std::set<std::string> const& ignore_params /* = {} */)
    {
      Configure(ValidateConfiguration(p, ignore_params));
    } // ValidateAndConfigure()


    //------------------------------------------------------------------------------------//
    void DetectorPropertiesStandard::Setup(providers_type providers) {

      SetGeometry(providers.get<geo::GeometryCore>());
      SetGArProperties(providers.get<detinfo::GArProperties>());
      SetECALProperties(providers.get<detinfo::ECALProperties>());
      SetDetectorClocks(providers.get<detinfo::DetectorClocks>());

    } // DetectorPropertiesStandard::Setup()


    //------------------------------------------------------------------------------------//
    double DetectorPropertiesStandard::Efield(unsigned int planegap) const
    {
      if(planegap >= fEfield.size())
        throw cet::exception("DetectorPropertiesStandard")
        << "requesting Electric field in a plane gap that is not defined\n";

      return fEfield[planegap];
    }


      //------------------------------------------------
    double DetectorPropertiesStandard::Density(double temperature) const
    {
      // Default temperature use internal value.
      if(temperature == 0.)
        temperature = Temperature();

      double density = -0.00615*temperature + 1.928;

      return density;
    } // DetectorPropertiesStandard::Density()


    //----------------------------------------------------------------------------------
    // Restricted mean energy loss (dE/dx) in units of MeV/cm.
    //
    // For unrestricted mean energy loss, set tcut = 0, or tcut large.
    //
    // Arguments:
    //
    // mom  - Momentum of incident particle in GeV/c.
    // mass - Mass of incident particle in GeV/c^2.
    // tcut - Maximum kinetic energy of delta rays (MeV).
    //
    // Returned value is positive.
    //
    // Based on Bethe-Bloch formula as contained in particle data book.
    // Material parameters (stored in GArProperties.fcl) are taken from
    // pdg web site http://pdg.lbl.gov/AtomicNuclearProperties/
    //
    double DetectorPropertiesStandard::Eloss(double mom,
                                             double mass,
                                             double tcut) const
    {
      // Some constants.

      double K = 0.307075;     // 4 pi N_A r_e^2 m_e c^2 (MeV cm^2/mol).
      double me = 0.510998918; // Electron mass (MeV/c^2).

      // Calculate kinematic quantities.

      double bg    = mom / mass;          // beta*gamma.
      double gamma = sqrt(1. + bg*bg);    // gamma.
      double beta  = bg / gamma;          // beta (velocity).
      double mer   = 0.001 * me / mass;   // electron mass / mass of incident particle.
      double tmax  = 2.*me* bg*bg / (1. + 2.*gamma*mer + mer*mer);  // Maximum delta ray energy (MeV).

      // Make sure tcut does not exceed tmax.

      if(tcut == 0. || tcut > tmax)
        tcut = tmax;

      // Calculate density effect correction (delta).

      double x = std::log10(bg);
      double delta = 0.;
      if(x >= fSternheimerParameters.x0) {
        delta = 2. * std::log(10.) * x - fSternheimerParameters.cbar;
        if(x < fSternheimerParameters.x1)
          delta += fSternheimerParameters.a * std::pow(fSternheimerParameters.x1 - x, fSternheimerParameters.k);
      }

      // Calculate stopping number.

      double B = 0.5 * std::log(2.*me*bg*bg*tcut / (1.e-12 * sqr(fGP->ExcitationEnergy())))
      - 0.5 * beta*beta * (1. + tcut / tmax) - 0.5 * delta;

      // Don't let the stopping number become negative.

      if(B < 1.)
        B = 1.;

        // Calculate dE/dx.

      double dedx = Density() * K*fGP->AtomicNumber()*B / (fGP->AtomicMass() * beta*beta);

      // Done.

      return dedx;
    } // DetectorPropertiesStandard::Eloss()

    //----------------------------------------------------------------------------------
    double DetectorPropertiesStandard::ElossVar(double mom,
                                                double mass) const
    {
      // Some constants.

      double K  = 0.307075;     // 4 pi N_A r_e^2 m_e c^2 (MeV cm^2/mol).
      double me = 0.510998918; // Electron mass (MeV/c^2).

      // Calculate kinematic quantities.

      double bg     = mom / mass;     // beta*gamma.
      double gamma2 = 1. + bg*bg;     // gamma^2.
      double beta2  = bg*bg / gamma2; // beta^2.

      // Calculate final result.

      double result = gamma2 * (1. - 0.5 * beta2) * me * (fGP->AtomicNumber() / fGP->AtomicMass()) * K * Density();
      return result;
    } // DetectorPropertiesStandard::ElossVar()


    //------------------------------------------------------------------------------------//
    double DetectorPropertiesStandard::DriftVelocity(double efield,
                                                     double temperature,
                                                     bool   cmPerns) const
    {

      // Efield should have units of kV/cm
      // Temperature should have units of Kelvin

      // Default Efield, use internal value.
      if(efield == 0.)
        efield = Efield();

      //
      //if(efield > 4.0)
      //  LOG_WARNING("DetectorPropertiesStandard")
      //  << "DriftVelocity Warning! : E-field value of "
      //  << efield
      //  << " kV/cm is outside of range covered by drift"
      // << " velocity parameterization. Returned value"
      //<< " may not be correct";

      // Default temperature use internal value.
      if(temperature == 0.)
        temperature = Temperature();

      // read in from fcl parameter

      double vd = fDriftVelocity; // cm/us.  For now just take it out of a fcl parameter.  Calcualted with magboltz and it's a strong function of gas composition

      if(cmPerns) return vd * 1.e-3; // cm/ns

      return vd; // in cm/us
    }

    //----------------------------------------------------------------------------------
    // The below function assumes that the user has applied the lifetime correction and
    // effective pitch between the wires (usually after 3D reconstruction). Using with
    // mean wire pitch will not give correct results.
    // parameters:
    //  dQdX in electrons/cm, charge (amplitude or integral obtained) divided by
    //         effective pitch for a given 3D track.
    // returns dEdX in MeV/cm
    double DetectorPropertiesStandard::BirksCorrection(double dQdx) const
    {
      // Correction for charge quenching using parameterization from
      // S.Amoruso et al., NIM A 523 (2004) 275
      double  A3t    = detinfo::kRecombA;
      double  K3t    = detinfo::kRecombk;                  // in KV/cm*(g/cm^2)/MeV
      double  rho    = Density();                          // GAr density in g/cm^3
      double Wion    = 1000./detinfo::kGeVToElectrons;     // 23.6 eV = 1e, Wion in MeV/e
      double E_field = Efield();                           // Electric Field in the drift region in KV/cm
      K3t           /= rho;                                // KV/MeV
      double dEdx    = dQdx/(A3t/Wion-K3t/E_field*dQdx);   //MeV/cm

      return dEdx;
    }

    //----------------------------------------------------------------------------------
    // Modified Box model correction
    double DetectorPropertiesStandard::ModBoxCorrection(double dQdx) const
    {
      // Modified Box model correction has better behavior than the Birks
      // correction at high values of dQ/dx.
      double rho     = Density();                          // GAr density in g/cm^3
      double Wion    = 1000./detinfo::kGeVToElectrons;     // 23.6 eV = 1e, Wion in MeV/e
      double E_field = Efield();                           // Electric Field in the drift region in KV/cm
      double Beta    = detinfo::kModBoxB / (rho * E_field);
      double Alpha   = detinfo::kModBoxA;
      double dEdx    = (std::exp(Beta * Wion * dQdx ) - Alpha) / Beta;

      return dEdx;

    }

    //------------------------------------------------------------------------------------//
    int  DetectorPropertiesStandard::TriggerOffset()     const
    {
      return fTPCClock.Ticks(fClocks->TriggerOffsetTPC() * -1.);
    }


    //--------------------------------------------------------------------
    // Take an X coordinate, and convert to a number of ticks, the
    // charge deposit occured at t=0
    double DetectorPropertiesStandard::ConvertXToTicks(double X) const
    {
      return (X / fXTicksCoefficient);
    }

    //-------------------------------------------------------------------
    // Take a cooridnate in ticks, and convert to an x position
    // assuming event deposit occured at t=0
    double  DetectorPropertiesStandard::ConvertTicksToX(double ticks) const
    {
      return ticks * fXTicksCoefficient;
    }


    //--------------------------------------------------------------------
    void DetectorPropertiesStandard::CheckIfConfigured() const
    {
      if (!fGeo)    throw cet::exception(__FUNCTION__) << "Geometry is uninitialized!";
      if (!fGP)     throw cet::exception(__FUNCTION__) << "GArPropertiesStandard is uninitialized!";
      if (!fECALP)  throw cet::exception(__FUNCTION__) << "ECALPropertiesStandard is uninitialized!";
      if (!fClocks) throw cet::exception(__FUNCTION__) << "DetectorClocks is uninitialized!";
    }

    //--------------------------------------------------------------------
    // Recalculte x<-->ticks conversion parameters from detector constants
    void DetectorPropertiesStandard::CalculateXTicksParams()
    {
      CheckIfConfigured();

      double samplingRate   = SamplingRate();
      double efield         = Efield();
      double temperature    = Temperature();
      double driftVelocity  = DriftVelocity(efield, temperature);

      fXTicksCoefficient    = 0.001 * driftVelocity * samplingRate;

      return;
    }

  } // namespace
} // gar
