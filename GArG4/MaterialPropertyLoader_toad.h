////////////////////////////////////////////////////////////////////////
/// \file MaterialPropertyLoader.h
//
/// \author  bjpjones@mit.edu
////////////////////////////////////////////////////////////////////////
// Class to set material properties for different materials in
// the detector. Currently used to supply optical properties to LAr
// and other optical elements of the detector.
//

#ifndef GARG4MaterialPropertyLoaderTOAD_h
#define GARG4MaterialPropertyLoaderTOAD_h

#include "Geant4/G4LogicalVolumeStore.hh"
#include <map>

namespace gar {
  namespace garg4 {
    
    class MaterialPropertyLoaderToad
    {
    public:
      
      MaterialPropertyLoaderToad() {}
      ~MaterialPropertyLoaderToad() {}
      
      
      
    public:
      
      //Accessors
      std::map<double,double> GetMaterialProperty(std::string const& Material,
                                                  std::string const& Property)
      {return fPropertyList[Material][Property];}
      
      double GetMaterialConstProperty(std::string const& Material,
                                      std::string const& Property)
      {return fConstPropertyList[Material][Property];}
      
      std::map<std::string,double> GetMaterialConstProperties(std::string const& Material)
      {return fConstPropertyList[Material];}
      
      std::map<std::string,std::map<double,double> >  GetMaterialProperties(std::string const& Material)
      {return fPropertyList[Material];}
      
      
      // Methods to set material properties
      void SetMaterialProperty(std::string             const& Material,
                               std::string             const& Property,
                               std::map<double,double> const& Values,
                               double                         Unit);
      void SetMaterialConstProperty(std::string const& Material,
                                    std::string const& Property,
                                    double             Value,
                                    double             Unit);
      
      // Method to set GArG4 Birks constant
      void SetBirksConstant(std::string const&,
                            double,
                            double);
      
      void SetReflectances(std::string                                      const&,
                           std::map<std::string, std::map<double, double> > const&,
                           std::map<std::string, std::map<double, double> > const&);
      
      // Set material properties supplied from services
      void GetPropertiesFromServices();
      
      // Geometry updater.  Generate Geant4 material property tables and attach to detector geometry
      void UpdateGeometry( G4LogicalVolumeStore* );
      
      
      
    private:
      
        //         materials                properties            values
      std::map < std::string , std::map < std::string,double> > fConstPropertyList;
      
        //         materials                properties               energies  values
      std::map < std::string , std::map < std::string , std::map < double ,  double > > > fPropertyList;
      
      std::map<std::string, double> fBirksConstants;
      
      
      
      
    };  
    
  }
} // gar

#endif
