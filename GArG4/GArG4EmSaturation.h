////////////////////////////////////////////////////////////////////////////
// Adapted from GEANT4 to accommodate Birks scintillator parameterization
// Search for GAr to find modified lines
////////////////////////////////////////////////////////////////////////////
//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id$
//
//
#ifndef GArG4EmSaturation_h
#define GArG4EmSaturation_h
// -------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4EmSaturation
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 18.02.2008
//
// Modifications:
//
//
// Class Description:
//   Compution on saturation effect, which reduce visible energy
//   deposition at the step. Default implementation takes into
//   account Birks effect. Birks coefficients for some materials
//   from G4 database on materials are provided
//
// -------------------------------------------------------------
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Geant4/globals.hh"
#include "Geant4/G4Step.hh"
#include "Geant4/G4ParticleDefinition.hh"
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4LossTableManager;
class G4NistManager;
class G4MaterialCutsCouple;
class G4Material;

namespace gar {
  namespace garg4 {
    class GArG4EmSaturation
    {
    public:

      GArG4EmSaturation();

      virtual ~GArG4EmSaturation();

      G4double VisibleEnergyDeposition(const G4Track *track,
        const G4ParticleDefinition*,
        const G4MaterialCutsCouple*,
        G4double length,
        G4double edepTotal,
        G4double edepNIEL = 0.0);

        inline G4double VisibleEnergyDeposition(const G4Step* step);

        // find and Birks coefficient
        G4double FindG4BirksCoefficient(const G4Material*);

        // dump coeffitients used in run time
        void DumpBirksCoefficients();

        // dump G4 list
        void DumpG4BirksCoefficients();

      private:

        // hide assignment operator
        GArG4EmSaturation & operator=(const GArG4EmSaturation &right);

        GArG4EmSaturation(const GArG4EmSaturation&);

        G4double FindBirksCoefficient(const G4Material*);

        void Initialise();

        const G4ParticleDefinition* electron;
        const G4ParticleDefinition* proton;
        G4LossTableManager*         manager;
        G4NistManager*              nist;

        // cash
        const G4Material*           curMaterial;
        G4double                    curBirks;
        G4double                    curChou;  // Nova specific!!!
        G4double                    curRatio;
        G4double                    curChargeSq;

        G4int    nMaterials;
        G4int    nG4Birks;

        // list of materials used in run time
        std::vector<const G4Material*>    matPointers;
        std::vector<G4String>             matNames;
        std::vector<G4double>             massFactors;
        std::vector<G4double>             effCharges;

        // list of G4 materials
        std::vector<G4double>             g4MatData;
        std::vector<G4String>             g4MatNames;
      };

    } // garg4
  } // gar

  inline G4double gar::garg4::GArG4EmSaturation::VisibleEnergyDeposition(const G4Step* step)
  {
    const G4Track* track = step->GetTrack();
    return VisibleEnergyDeposition(track,
      track->GetParticleDefinition(),
      track->GetMaterialCutsCouple(),
      step->GetStepLength(),
      step->GetTotalEnergyDeposit(),
      step->GetNonIonizingEnergyDeposit());
    }

    #endif
