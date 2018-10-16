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
        fNoiseMeV = pset.get<float>("NoiseMeV", 0.1);
        fCellSize = pset.get<double>("DefaultCellSize", 2.); // CellSize in cm
        fSaturation = pset.get<bool>("Saturation", false);
        fTimeSmearing = pset.get<bool>("TimeSmearing", false);
        fTimeResolution = pset.get<double>("TimeResolution", 1.); // default time resolution in ns
        fADCSaturation = pset.get<int>("ADCSaturation", 4096); // default to 12-bit ADC's
        fMeVtoMIP = pset.get<double>("MeVtoMIP", 0.814); // default MeVtoMIP for 5 mm Scint

        //MultiSegmentation
        fMultiSeg = pset.get<bool>("MultiSegmentation", false);
        fMultiSegLayer = pset.get<unsigned int>("MultiSegmentation_FromLayer");
        fMultiSegCellSize = pset.get<float>("MultiSegmentation_CellSize");

        if(fMultiSeg)
        {
          std::cout << "ECALReadoutSimStandardAlg::reconfigure() - MultiSegmentation put to " << fMultiSeg << std::endl;
          std::cout << "ECALReadoutSimStandardAlg::reconfigure() - From layer " << fMultiSegLayer
          << " change of the granularity from " << fCellSize
          << " cm to " << fMultiSegCellSize << " cm" << std::endl;
        }
        else{
          std::cout << "ECALReadoutSimStandardAlg::reconfigure() - Normal segmentation of " << fCellSize << " cm" << std::endl;
        }

        fECALUtils = std::make_unique<ECALUtils>(fDetProp->EffectivePixel(), 0.95);

        fECALLayerThickness = fGeo->GetECALLayerThickness(); // in cm
        fEndcapStartXPosition = fGeo->GetECALEndcapStartPosition(); // in cm
        fECALRinner = fGeo->GetECALInnerBarrelRadius(); // in cm
        fECALRouter = fGeo->GetECALOuterBarrelRadius(); // in cm
        fECALPVThickness = fGeo->GetPVThickness(); // in cm

        fAbsorberThickness = fGeo->GetECALAbsorberThickness();
        fActiveMatThickness = fGeo->GetECALActiveMatThickness();

        fTPCOriginX = fGeo->TPCXCent();
        fTPCOriginY = fGeo->TPCYCent();
        fTPCOriginZ = fGeo->TPCZCent();

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
          float energy = SimCaloHit.Energy();
          float time = SimCaloHit.Time();
          float x = SimCaloHit.X();
          float y = SimCaloHit.Y();
          float z = SimCaloHit.Z();
          unsigned int id = SimCaloHit.CaloID();
          unsigned int layer = SimCaloHit.Layer();

          if(fAddNoise)
          this->AddElectronicNoise(energy);

          this->DoPhotonStatistics(energy);

          if(fTimeSmearing)
          this->DoTimeSmearing(time);

          raw::CaloRawDigit digithit = raw::CaloRawDigit(static_cast<unsigned int>(energy), time, x, y, z, id, cellID, layer);
          digCol.emplace_back(digithit);
        }
      }

      //----------------------------------------------------------------------------
      sdp::CaloDeposit ECALReadoutSimStandardAlg::Segmentation(sdp::CaloDeposit const &SimCaloHit)
      {
        LOG_DEBUG("ECALReadoutSimStandardAlg") << "Segmentation()";

        unsigned long long int cellID = 0;
        unsigned int layer = 0;
        //Need to change the origin from world to local with the center of the TPC as origin
        float x = SimCaloHit.X() - fTPCOriginX;
        float y = SimCaloHit.Y() - fTPCOriginY;
        float z = SimCaloHit.Z() - fTPCOriginZ;
        unsigned int id = SimCaloHit.CaloID();

        //Do the Segmentation per type of hits based on the id (endcap and barrel separately)

        //Axes definition, z is in the plane
        //        y
        //        ^
        //        |
        // x <-----

        //Case 1 Endcap (endcap is similar in y/z)
        if(id == 3)
        {
          //The Endcap plane is in yz.. x defines the depth of the layer
          layer = std::floor( (std::abs(x) - fEndcapStartXPosition) / fECALLayerThickness );//starts at 0

          if(fMultiSeg && layer >= fMultiSegLayer) fCellSize = fMultiSegCellSize;

          double offset = 0.5 * fCellSize; //in cm

          //Get the bin position of the hit based on the cellsize
          int zbin = fECALUtils->PositionToBin(z, fCellSize, offset);
          int ybin = fECALUtils->PositionToBin(y, fCellSize, offset);

          //Create a unique cellID based on the bins in y, z and the layer
          cellID = fECALUtils->MakeCellID(id, zbin, ybin, layer);

          //Get the new position based on the bin number
          if(x > 0)
          x = fEndcapStartXPosition + layer * fECALLayerThickness + (fAbsorberThickness + fActiveMatThickness / 2.0);
          else
          x = - (fEndcapStartXPosition + layer * fECALLayerThickness + (fAbsorberThickness + fActiveMatThickness / 2.0));

          y = fECALUtils->BinToPosition(ybin, fCellSize, offset);
          z = fECALUtils->BinToPosition(zbin, fCellSize, offset);
        }

        //Case 2 Inner Barrel
        if(id == 1)
        {
          float R = std::sqrt(y*y + z*z);
          layer = std::floor( (std::abs(R) - fECALRinner) / fECALLayerThickness );//starts at 0

          if(fMultiSeg && layer >= fMultiSegLayer) fCellSize = fMultiSegCellSize;

          double offset = 0.5 * fCellSize; //in cm
          float dPhi = std::atan(fCellSize / fECALRinner); //in rad

          float phi_angle = std::acos(z/R); //in rad

          int xbin = fECALUtils->PositionToBin(x, fCellSize, offset);
          int phibin = fECALUtils->PositionToBin(phi_angle, dPhi, 0.5 * dPhi);

          //Create a unique cellID based on the bins in x, phi and the layer
          cellID = fECALUtils->MakeCellID(id, xbin, phibin, layer);

          //Get the new position based on the bin number
          x = fECALUtils->BinToPosition(xbin, fCellSize, offset);
          if(y < 0)
          y = - R * std::sin( fECALUtils->BinToPosition(phibin, dPhi, 0.5 * dPhi) );
          else
          y = R * std::sin( fECALUtils->BinToPosition(phibin, dPhi, 0.5 * dPhi) );
          z = R * std::cos( fECALUtils->BinToPosition(phibin, dPhi, 0.5 * dPhi) );
        }

        //Case 3 Outer Barrel
        if(id == 2)
        {
          float R = std::sqrt(y*y + z*z);
          layer = std::floor( (std::abs(R) - (fECALRinner + fECALPVThickness)) / fECALLayerThickness );//starts at last InnerBarrelLayer + 1

          if(fMultiSeg && layer >= fMultiSegLayer) fCellSize = fMultiSegCellSize;

          double offset = 0.5 * fCellSize; //in cm
          float dPhi = std::atan(fCellSize / fECALRouter); //in rad

          float phi_angle = std::acos(z/R); //in rad

          int xbin = fECALUtils->PositionToBin(x, fCellSize, offset);
          int phibin = fECALUtils->PositionToBin(phi_angle, dPhi, 0.5 * dPhi);

          //Create a unique cellID based on the bins in x, phi and the layer
          cellID = fECALUtils->MakeCellID(id, xbin, phibin, layer);

          //Get the new position based on the bin number
          x = fECALUtils->BinToPosition(xbin, fCellSize, offset);
          if(y < 0)
          y = - R * std::sin( fECALUtils->BinToPosition(phibin, dPhi, 0.5 * dPhi) );
          else
          y = R * std::sin( fECALUtils->BinToPosition(phibin, dPhi, 0.5 * dPhi) );
          z = R * std::cos( fECALUtils->BinToPosition(phibin, dPhi, 0.5 * dPhi) );
        }

        return sdp::CaloDeposit(SimCaloHit.TrackID(), SimCaloHit.Time(), SimCaloHit.Energy(), x + fTPCOriginX, y + fTPCOriginY, z + fTPCOriginZ, id, cellID, layer);
      }

      //----------------------------------------------------------------------------
      void ECALReadoutSimStandardAlg::AddSubHits(std::map<unsigned long long int, std::vector<sdp::CaloDeposit> > m_SimCaloHits, std::vector<sdp::CaloDeposit> &SimCaloHitVec)
      {
        LOG_DEBUG("ECALReadoutSimStandardAlg") << "AddSubHits()";

        for(std::map<unsigned long long int, std::vector<sdp::CaloDeposit> >::iterator it = m_SimCaloHits.begin(); it != m_SimCaloHits.end(); ++it)
        {
          float EsubhitSum = 0.;
          float TsubhitSum = 100000.;
          float x = 0.;
          float y = 0.;
          float z = 0.;
          unsigned int layer = 0;
          unsigned int id = 0;
          unsigned long long int cellID = it->first;

          for(std::vector<sdp::CaloDeposit>::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2)
          {
            sdp::CaloDeposit SimSubHit = *it2;
            float energy = SimSubHit.Energy();
            float time = SimSubHit.Time();

            x = SimSubHit.X();
            y = SimSubHit.Y();
            z = SimSubHit.Z();
            id = SimSubHit.CaloID();
            layer = SimSubHit.Layer();

            EsubhitSum += energy;

            if(time < TsubhitSum)
            TsubhitSum = time;
          }

          sdp::CaloDeposit SimHit = sdp::CaloDeposit(-1, TsubhitSum, EsubhitSum, x, y, z, id, cellID, layer);
          SimCaloHitVec.push_back(SimHit);
        }

        return;
      }

      //----------------------------------------------------------------------------
      void ECALReadoutSimStandardAlg::DoPhotonStatistics(float &energy)
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
      void ECALReadoutSimStandardAlg::DoTimeSmearing(float &time)
      {
        LOG_DEBUG("ECALReadoutSimStandardAlg") << "DoTimeSmearing()";
        CLHEP::RandGauss GausRand(fEngine);
        time = GausRand.fire(time, fTimeResolution);
        return;
      }

      //----------------------------------------------------------------------------
      void ECALReadoutSimStandardAlg::AddElectronicNoise(float &energy)
      {
        LOG_DEBUG("ECALReadoutSimStandardAlg") << "AddElectronicNoise()";

        //Random noise shooting (Gaussian electronic noise)
        CLHEP::RandGauss GaussRand(fEngine);
        energy += GaussRand.shoot(0., fNoiseMeV * CLHEP::MeV / CLHEP::GeV);

        return;
      }

    } // rosim
  } // gar
