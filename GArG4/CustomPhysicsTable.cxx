////////////////////////////////////////////////////////////////////////
/// \file CustomPhysicsTable.cxx
//
////////////////////////////////////////////////////////////////////////
// See header file for full description

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "GArG4/CustomPhysicsTable.hh"
#include "GArG4/CustomPhysicsFactory.hh"

#include "Geant4/G4VPhysicsConstructor.hh"

namespace gar {
  namespace garg4 {
    
    CustomPhysicsTable* TheCustomPhysicsTable{nullptr};
    
    //-----------------------------------------------------------------
    CustomPhysicsTable::CustomPhysicsTable(CustomPhysicsFactoryBase * Factory)
    {
      if(!TheCustomPhysicsTable){
        TheCustomPhysicsTable = new CustomPhysicsTable;
        TheCustomPhysicsTable->AddPhysics(Factory);
      }
      else{
        TheCustomPhysicsTable->AddPhysics(Factory);
      }
    }
    
    //-----------------------------------------------------------------
    std::vector<std::string> CustomPhysicsTable::GetAvailablePhysicsList()
    {
      std::vector<std::string> ReturnVector;
      for(std::map<std::string,CustomPhysicsFactoryBase* >::const_iterator i = TheCustomPhysicsTable->theTable.begin(); i!=TheCustomPhysicsTable->theTable.end(); i++){
        ReturnVector.push_back((*i).first);
      }
      return ReturnVector;
    }
    
      //-----------------------------------------------------------------
    bool CustomPhysicsTable::IsPhysicsAvailable(std::string PhysicsName)
    {
      if(!TheCustomPhysicsTable->theTable[PhysicsName])
        return false;
      else
        return true;
    }
    
      //-----------------------------------------------------------------
    G4VPhysicsConstructor * CustomPhysicsTable::GetPhysicsConstructor(std::string PhysicsName)
    {
      if(IsPhysicsAvailable(PhysicsName)){
        G4VPhysicsConstructor * G4VPC =  TheCustomPhysicsTable->theTable[PhysicsName]->Build();
        return G4VPC;
      }
      else{
        G4VPhysicsConstructor * G4VPC=0;
        return G4VPC;
      }
      
    }
    
    //-----------------------------------------------------------------
    void CustomPhysicsTable::AddPhysics(CustomPhysicsFactoryBase * Factory)
    {
      
      if(IsPhysicsAvailable(Factory->GetName()))
        mf::LogWarning("CustomPhysicsTable") << "Physics constructor being overwritten"
					   << " in CustomPhysicsTable";
      TheCustomPhysicsTable->theTable[Factory->GetName()]=Factory;
      MF_LOG_DEBUG("CustomPhysicsTable")<<"CustomPhysicsTable : Physics Table registering new physics "
      << Factory->GetName();
    }
    
  }
  
} // gar

//Ben Jones, MIT, Sept 2009
