/**
 * @file   garsoft/Utilities/GeneratedParticleInfo.h
 * @brief  Contains data associated to particles from detector simulation.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   January 10, 2018
 * 
 */

#ifndef GARSOFT_SIMULATION_GENERATEDPARTICLEINFO_H
#define GARSOFT_SIMULATION_GENERATEDPARTICLEINFO_H

#include "nusimdata/SimulationBase/simb.h" // sim::GeneratedParticleIndex_t, ...


namespace sim {
  
  /**
   * @brief Contains information about a generated particle.
   * 
   * This information is designed to be added as metadata to an association
   * between a `simb::MCParticle` produced by the simulation of particle
   * propagation through the detector (e.g. `GArG4`) with the generator truth
   * record that particle derives from:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * art::Assns<simb::MCParticle, simb::MCTruth, sim::GeneratedParticleInfo>
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * 
   * The associated `simb::MCParticle` is expected to be the generated particle
   * of the generator truth record, as processed by the simulation.
   * The provided information is currently:
   * * index of the generated particle within the generator record; the particle
   *   itself can be obtained as:
   *   @code{.cpp}
   *   mctruthPtr->GetParticle(genParticleInfo.generatedParticleIndex());
   *   @endcode
   * 
   */
  class GeneratedParticleInfo {
    
      public:
    
    //--- BEGIN Generated particle index ---------------------------------------
    /// @{
    /// @name Generated particle index
    
    /// Type of particle index in the generator truth record (`simb::MCTruth`).
    using GeneratedParticleIndex_t = simb::GeneratedParticleIndex_t;
    
    /// Constant representing the absence of generator truth information.
    static constexpr GeneratedParticleIndex_t NoGeneratedParticleIndex
      = simb::NoGeneratedParticleIndex;
    
    /**
     * @brief  Returns the generated particle index.
     * @return generated particle index (`NoGeneratedParticleIndex` if none)
     */
    GeneratedParticleIndex_t generatedParticleIndex() const
      { return fGeneratedParticleIndex; }
    
    /// Returns whether the specified one is an acceptable generator index.
    bool hasGeneratedParticleIndex() const
      { return simb::isGeneratedParticleIndex(fGeneratedParticleIndex); }
    
    /// @}
    //--- END Generated particle index -----------------------------------------
    
    
    //--- BEGIN Constructors ---------------------------------------------------
    /// Default constructor (for ROOT I/O only).
    GeneratedParticleInfo() = default;
    
    /**
     * @brief Constructor: assigns all the values.
     * @param index generated particle index (`NoGeneratedParticleIndex` if none)
     */
    GeneratedParticleInfo(GeneratedParticleIndex_t index)
      : fGeneratedParticleIndex(index)
      {}
    //--- END Constructors ---------------------------------------------------
    
    
      private:
    
    /// Index of the generated particle in its `simb::MCTruth` record.
    GeneratedParticleIndex_t fGeneratedParticleIndex = NoGeneratedParticleIndex;
    
  }; // class GeneratedParticleInfo
  
  
  
} // namespace sim

#endif // GARSOFT_SIMULATION_GENERATEDPARTICLEINFO_H
