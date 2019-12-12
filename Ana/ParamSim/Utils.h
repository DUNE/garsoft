#ifndef UTILS_H
#define UTILS_H 1

#include "TVector3.h"
#include "TRandom3.h"

class Utils
{
public:

    /* Default constructor */
    Utils();

    /* Copy constructor */
    Utils(const Utils &) = default;

    /* Destructor */
    ~Utils();

    /* Set the seed */
    void SetSeed(int seed);

    /* Check if MCP started in tracker */
    bool hasOriginInTracker(TVector3 spoint);

    /* Check if MCP decayed in calo */
    bool hasDecayedInCalo(TVector3 epoint);

    /* Check if the mcp is a backscatter */
    bool isBackscatter(TVector3 spoint, TVector3 epoint);

    /* Check if it is a Bremsstrahlung photon */
    bool isBremsstrahlung(TVector3 spoint, int pdg, int motherpdg);

    float GetRamdomNumber() { return _rando->Rndm(); }

    float GaussianSmearing(float mean, float sigma) { return _rando->Gaus(mean, sigma); }

private:
    unsigned long int _seed;         ///< seed
    TRandom3 *_rando;                ///< random generator
};

#endif
