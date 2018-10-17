////////////////////////////////////////////////////////////////////////
/// \file MaterialPropertyLoader.cxx
//
/// \author  bjpjones@mit.edu
////////////////////////////////////////////////////////////////////////
// Class to set material properties for different materials in
// the detector. Currently mainly used to set optical properties
// for LAr and other optical components
//



#include "GArG4/MaterialPropertyLoader.h"
#include "CoreUtils/ServiceUtil.h"

#include "Geant4/G4Material.hh"
#include "Geant4/G4MaterialPropertiesTable.hh"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "DetectorInfo/ECALProperties.h"
#include "DetectorInfo/ECALPropertiesService.h"
#include "CoreUtils/ServiceUtil.h"

namespace gar {
  namespace garg4 {

    //----------------------------------------------
    void MaterialPropertyLoader::SetMaterialProperty(std::string              const& Material,
                                                     std::string              const& Property,
                                                     std::map<double, double> const& PropertyVector,
                                                     double                          Unit)
    {
      std::map<double,double> PropVectorWithUnit;
      for(auto const& itr : PropertyVector){
        PropVectorWithUnit[itr.first * CLHEP::eV] = itr.second * Unit;
      }

      fPropertyList[Material][Property] = PropVectorWithUnit;

      LOG_INFO("MaterialPropertyLoader")
      <<"Added property "
      << Material<< "  "
		  << Property;
    }

    //----------------------------------------------
    void MaterialPropertyLoader::SetMaterialConstProperty(std::string const& Material,
                                                          std::string const& Property,
                                                          double             PropertyValue,
                                                          double             Unit)
    {
      fConstPropertyList[Material][Property]=PropertyValue*Unit;

      LOG_INFO("MaterialPropertyLoader")
      << "Added const property "
      << Material
      << "  "
      << Property
      << " = "
      << PropertyValue;
    }

    //----------------------------------------------
    void MaterialPropertyLoader::SetBirksConstant(std::string const& Material,
                                                  double             PropertyValue,
                                                  double             Unit)
    {
      fBirksConstants[Material] = PropertyValue * Unit;
      LOG_INFO("MaterialPropertyLoader")
      << "Set Birks constant "
      << Material;
    }

    //----------------------------------------------
    void MaterialPropertyLoader::UpdateGeometry(G4LogicalVolumeStore * lvs)
    {
      std::map<std::string,G4MaterialPropertiesTable*> MaterialTables;
      std::map<std::string,bool> MaterialsSet;

      LOG_INFO("MaterialPropertyLoader")
      << "UPDATING GEOMETRY";

      // Loop over each material with a property vector and create a new material table for it
      for(auto const& itr : fPropertyList) {
        std::string Material     = itr.first;
        MaterialsSet[Material]   = true;
        MaterialTables[Material] = new G4MaterialPropertiesTable;
      }

      // Loop over each material with a const property,
      // if material table does not exist, create one
      for(auto const& itr : fConstPropertyList){
        std::string Material = itr.first;
        if(!MaterialsSet[Material]){
          MaterialsSet[Material]   = true;
          MaterialTables[Material] = new G4MaterialPropertiesTable;
        }
      }

      // For each property vector, convert to an array of g4doubles and
      // feed to materials table Lots of firsts and seconds!  See annotation
      // in MaterialPropertyLoader.h to follow what each element is

      for(auto const& itr : fPropertyList){
        std::string Material = itr.first;
        for(auto const& jitr : itr.second){
          std::string Property = jitr.first;
          std::vector<G4double> g4MomentumVector;
          std::vector<G4double> g4PropertyVector;

          for(auto const& kitr : jitr.second){
            g4MomentumVector.push_back(kitr.first);
            g4PropertyVector.push_back(kitr.second);
          }
          int NoOfElements = g4MomentumVector.size();
          MaterialTables[Material]->AddProperty(Property.c_str(),
                                                &g4MomentumVector[0],
                                                &g4PropertyVector[0],
                                                NoOfElements);
          LOG_INFO("MaterialPropertyLoader")
          << "Added property "
          << Property
          << " to material table "
          << Material;
        }
      }

      //Add each const property element
      for(auto const& itr : fConstPropertyList){
        std::string Material = itr.first;
        for(auto const& jitr : itr.second){
          std::string Property   = jitr.first;
          G4double PropertyValue = jitr.second;
          MaterialTables[Material]->AddConstProperty(Property.c_str(), PropertyValue);

          LOG_INFO("MaterialPropertyLoader")
          << "Added const property "
          << Property
          << " to material table "
          << Material;
        }
      }

      //Loop through geometry elements and apply relevant material table where materials match
      for ( G4LogicalVolumeStore::iterator i = lvs->begin(); i != lvs->end(); ++i ){
        G4LogicalVolume* volume = (*i);
        G4Material* TheMaterial = volume->GetMaterial();
        std::string Material = TheMaterial->GetName();
        for(auto const& jitr : MaterialTables){
          if(Material == jitr.first){
            TheMaterial->SetMaterialPropertiesTable(jitr.second);
            //Birks Constant, for some reason, must be set separately
            if(fBirksConstants[Material]!=0)
              TheMaterial->GetIonisation()->SetBirksConstant(fBirksConstants[Material]);
            volume->SetMaterial(TheMaterial);
          }
        }
      }
    }

    //--------------------------------------------------------------------------
    void MaterialPropertyLoader::SetReflectances(std::string                                     const& /*Material*/,
                                                 std::map<std::string,std::map<double, double> > const& Reflectances,
                                                 std::map<std::string,std::map<double, double> > const& DiffuseFractions)
    {
      std::map<double, double> ReflectanceToStore;
      std::map<double, double> DiffuseToStore;

      for(auto const& itMat : Reflectances){
        std::string ReflectancePropName = std::string("REFLECTANCE_") + itMat.first;
        ReflectanceToStore.clear();
        for(auto const& itEn : itMat.second){
          ReflectanceToStore[itEn.first] = itEn.second;
        }
        SetMaterialProperty("GAr", ReflectancePropName, ReflectanceToStore,1);
      }

      for(auto const& itMat : DiffuseFractions){
        std::string DiffusePropName = std::string("DIFFUSE_REFLECTANCE_FRACTION_") + itMat.first;
        DiffuseToStore.clear();
        for(auto const& itEn : itMat.second){
          DiffuseToStore[itEn.first] = itEn.second;
        }
        SetMaterialProperty("GAr",
                            DiffusePropName,
                            DiffuseToStore,
                            1);
      }

    }

    //--------------------------------------------------------------------------
    void MaterialPropertyLoader::GetPropertiesFromServices()
    {
//      const detinfo::GArProperties*      GarProp = gar::providerFrom<detinfo::GArPropertiesService>();
//      const detinfo::DetectorProperties* DetProp = gar::providerFrom<detinfo::DetectorPropertiesService>();

        const detinfo::ECALProperties*     EcalProp = gar::providerFrom<detinfo::ECALPropertiesService>();

        SetMaterialConstProperty("Scintillator", "DENSITY", 1.06, CLHEP::g/CLHEP::cm3);
        SetBirksConstant("Scintillator", EcalProp->ScintBirksConstant(), CLHEP::mm/CLHEP::MeV);
    }

  } // garg4
} // gar
