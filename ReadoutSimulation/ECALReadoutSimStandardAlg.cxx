//
//  ECALReadoutSimStandardAlg.cxx
//
//  Created by Eldwan Brianne on 8/27/18
//

#include "fhiclcpp/ParameterSet.h"

#include "ReadoutSimulation/ECALReadoutSimStandardAlg.h"
#include "ReadoutSimulation/IonizationAndScintillation.h"
#include "DetectorInfo/DetectorPropertiesService.h"
#include "Geometry/Geometry.h"
#include "CoreUtils/ServiceUtil.h"
#include "RawDataProducts/raw.h"
#include "SimulationDataProducts/CaloDeposit.h"

#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandBinomial.h"

namespace gar {
  namespace rosim{

    //----------------------------------------------------------------------------
    ECALReadoutSimStandardAlg::ECALReadoutSimStandardAlg(CLHEP::HepRandomEngine      & engine,
      fhicl::ParameterSet    const& pset)
      : ECALReadoutSimAlg(engine, pset)
      {
        fGeo = gar::providerFrom<geo::Geometry>();
        fDetProp = gar::providerFrom<detinfo::DetectorPropertiesService>();
        this->reconfigure(pset);

        return;
      }

      //----------------------------------------------------------------------------
      ECALReadoutSimStandardAlg::~ECALReadoutSimStandardAlg()
      {
        return;
      }

      //----------------------------------------------------------------------------
      void ECALReadoutSimStandardAlg::reconfigure(fhicl::ParameterSet const& pset)
      {
        fAddNoise = pset.get<bool>("AddNoise", false);
        fCellSize = pset.get<double>("CellSize", 2.); // CellSize in cm
        fSaturation = pset.get<bool>("Saturation", false);
        fTimeSmearing = pset.get<bool>("TimeSmearing", false);
        fTimeResolution = pset.get<double>("TimeResolution", 1); // default time resolution in ns
        fADCSaturation = pset.get<int>("ADCSaturation", 4095); // default to 12-bit ADC's
        fMeVtoMIP = pset.get<double>("MeVtoMIP", 0.814); // default MeVtoMIP for 5 mm Scint

        fECALUtils = std::make_unique<ECALUtils>(fDetProp->EffectivePixel(), 0.95);

        return;
      }

      //----------------------------------------------------------------------------
      void ECALReadoutSimStandardAlg::CreateCaloRawDigits(std::vector<sdp::CaloDeposit> SimCaloVec, std::vector<raw::CaloRawDigit> &digCol)
      {
        LOG_DEBUG("ECALReadoutSimStandardAlg") << "CreateCaloRawDigits()";

        std::map<unsigned long long int, std::vector<sdp::CaloDeposit> > m_SimCaloHits;

        for(auto const &SimCaloHit : SimCaloVec)
        {
          //Create the raw digit hit
          sdp::CaloDeposit SegSimCaloHit = this->Segmentation(SimCaloHit);
          //push the hit into a map with the cellID as key
          m_SimCaloHits[SegSimCaloHit.CellID()].push_back(SegSimCaloHit);
        }

        std::vector<sdp::CaloDeposit> SimCaloHitVec;
        this->AddSubHits(m_SimCaloHits, SimCaloHitVec);

        for(auto const &SimCaloHit : SimCaloHitVec)
        {
          unsigned long long int cellID = SimCaloHit.CellID();
          double energy = SimCaloHit.Energy();
          double time = SimCaloHit.Time();
          double x = SimCaloHit.X();
          double y = SimCaloHit.Y();
          double z = SimCaloHit.Z();
          double id = SimCaloHit.CaloID();
          unsigned int layer = SimCaloHit.Layer();

          this->DoPhotonStatistics(energy);

          if(fTimeSmearing)
          this->DoTimeSmearing(time);

          raw::CaloRawDigit digithit = raw::CaloRawDigit(cellID, energy, time, x, y, z, layer, id);
          digCol.emplace_back(digithit);
        }
      }

      //----------------------------------------------------------------------------
      sdp::CaloDeposit ECALReadoutSimStandardAlg::Segmentation(sdp::CaloDeposit const &SimCaloHit)
      {
        LOG_DEBUG("ECALReadoutSimStandardAlg") << "Segmentation()";

        unsigned long long int cellID = 0;
        unsigned int layer = 0;
        double x = SimCaloHit.X();
        double y = SimCaloHit.Y();
        double z = SimCaloHit.Z();
        int id = SimCaloHit.CaloID();

        //Do the Segmentation per type of hits based on the id (endcap and barrel separately)

        //Axes definition, z is in the plane
        //        y
        //        ^
        //        |
        // x <-----

        //TO DO (getting these values from geo)
        double layer_thickness = 0.8; //in cm

        //Case 1 Endcap (endcap is similar in y/z)
        if(id == 3)
        {
          double offset = 0.5 * fCellSize; //in cm
          double xEndcapStart = 260.0; //in cm

          //The Endcap plane is in yz.. x defines the depth of the layer
          layer = std::floor( (std::abs(x) - xEndcapStart) / layer_thickness );//starts at 0

          //Get the bin position of the hit based on the cellsize
          int zbin = fECALUtils->PositionToBin(z, fCellSize, offset);
          int ybin = fECALUtils->PositionToBin(y, fCellSize, offset);

          //Create a unique cellID based on the bins in y, z and the layer
          cellID = fECALUtils->MakeCellID(id, zbin, ybin, layer);

          //Get the new position based on the bin number
          if(x > 0)
          x = xEndcapStart + layer * layer_thickness + 0.45;
          else
          x = - (xEndcapStart + layer * layer_thickness + 0.45);

          y = fECALUtils->BinToPosition(ybin, fCellSize, offset);
          z = fECALUtils->BinToPosition(zbin, fCellSize, offset);
        }

        //Case 2 Inner Barrel
        if(id == 1)
        {
          double Rinner = 274; //in cm
          double offset = 0.5 * fCellSize; //in cm
          double dPhi = std::atan(fCellSize / Rinner); //in rad

          double R = std::sqrt(y*y + z*z);
          layer = std::floor( (std::abs(R) - Rinner) / layer_thickness );//starts at 0

          double phi_angle = std::acos(z/R); //in rad

          int xbin = fECALUtils->PositionToBin(x, fCellSize, offset);
          int phibin = fECALUtils->PositionToBin(phi_angle, dPhi, 0.5 * dPhi);

          //Create a unique cellID based on the bins in x, phi and the layer
          cellID = fECALUtils->MakeCellID(id, xbin, phibin, layer);

          //Get the new position based on the bin number
          x = fECALUtils->BinToPosition(xbin, fCellSize, offset);
          y = R * std::sin( fECALUtils->BinToPosition(phibin, dPhi, 0.5 * dPhi) );
          z = R * std::cos( fECALUtils->BinToPosition(phibin, dPhi, 0.5 * dPhi) );
        }

        //Case 3 Outer Barrel
        if(id == 2)
        {
          double Router = 300; //in cm
          double offset = 0.5 * fCellSize; //in cm
          double dPhi = std::atan(fCellSize / Router); //in rad

          double R = std::sqrt(y*y + z*z);
          layer = std::floor( (std::abs(R) - Router) / layer_thickness );//starts at 0

          double phi_angle = std::acos(z/R); //in rad

          int xbin = fECALUtils->PositionToBin(x, fCellSize, offset);
          int phibin = fECALUtils->PositionToBin(phi_angle, dPhi, 0.5 * dPhi);

          //Create a unique cellID based on the bins in x, phi and the layer
          cellID = fECALUtils->MakeCellID(id, xbin, phibin, layer);

          //Get the new position based on the bin number
          x = fECALUtils->BinToPosition(xbin, fCellSize, offset);
          y = R * std::sin( fECALUtils->BinToPosition(phibin, dPhi, 0.5 * dPhi) );
          z = R * std::cos( fECALUtils->BinToPosition(phibin, dPhi, 0.5 * dPhi) );
        }

        return sdp::CaloDeposit(SimCaloHit.TrackID(), SimCaloHit.Time(), SimCaloHit.Energy(), x, y, z, id, cellID, layer);
      }

      //----------------------------------------------------------------------------
      void ECALReadoutSimStandardAlg::AddSubHits(std::map<unsigned long long int, std::vector<sdp::CaloDeposit> > m_SimCaloHits, std::vector<sdp::CaloDeposit> &SimCaloHitVec)
      {
        LOG_DEBUG("ECALReadoutSimStandardAlg") << "AddSubHits()";

        for(std::map<unsigned long long int, std::vector<sdp::CaloDeposit> >::iterator it = m_SimCaloHits.begin(); it != m_SimCaloHits.end(); ++it)
        {
          double EsubhitSum = 0.;
          double TsubhitSum = 0.;
          unsigned int nsubhits = 0;
          double x = 0.;
          double y = 0.;
          double z = 0.;
          unsigned int layer = 0;
          int id = 0;
          unsigned long long int cellID = it->first;

          for(std::vector<sdp::CaloDeposit>::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2)
          {
            sdp::CaloDeposit SimSubHit = *it2;
            double energy = SimSubHit.Energy();
            double time = SimSubHit.Time();

            x = SimSubHit.X();
            y = SimSubHit.Y();
            z = SimSubHit.Z();
            id = SimSubHit.CaloID();
            layer = SimSubHit.Layer();

            EsubhitSum += energy;
            TsubhitSum += time;
            nsubhits++;
          }

          sdp::CaloDeposit SimHit = sdp::CaloDeposit(-1, TsubhitSum/nsubhits, EsubhitSum/nsubhits, x, y, z, id, cellID, layer);
          SimCaloHitVec.push_back(SimHit);
        }

        return;
      }

      //----------------------------------------------------------------------------
      void ECALReadoutSimStandardAlg::DoPhotonStatistics(double &energy)
      {
        LOG_DEBUG("ECALReadoutSimStandardAlg") << "DoPhotonStatistics()";
        CLHEP::RandBinomial BinomialRand(fEngine);

        //convertion from GeV to MIP
        energy /= fMeVtoMIP * CLHEP::MeV / CLHEP::GeV;
        //conversion to px
        energy *= fDetProp->LightYield();
        //Saturation
        if(fSaturation)
        energy = fECALUtils->Saturate(energy);

        //Binomial Smearing
        double prob = energy / fDetProp->EffectivePixel();
        energy += BinomialRand.shoot(fDetProp->EffectivePixel(), prob);

        //Convertion to ADC
        if(energy > 0)
        energy *= fDetProp->SiPMGain();
        else
        energy = 0.;

        if(energy > fADCSaturation)
        energy = fADCSaturation;

        return;
      }

      //----------------------------------------------------------------------------
      void ECALReadoutSimStandardAlg::DoTimeSmearing(double &time)
      {
        LOG_DEBUG("ECALReadoutSimStandardAlg") << "DoTimeSmearing()";
        CLHEP::RandGauss GausRand(fEngine);
        time = GausRand.fire(time, fTimeResolution);
        return;
      }

      //----------------------------------------------------------------------------
      void ECALReadoutSimStandardAlg::CreateNoiseDigits(std::vector<raw::CaloRawDigit> & digits)
      {
        LOG_DEBUG("ECALReadoutSimStandardAlg") << "CreateNoiseDigits()";
        if(!fAddNoise) return;

        //TO DO

        return;
      }

    } // rosim
  } // gar
