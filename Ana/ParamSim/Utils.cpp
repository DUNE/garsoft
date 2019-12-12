#include "Utils.h"

Utils::Utils()
: _seed(0), _rando(new TRandom3(0))
{
}

Utils::~Utils()
{
    delete _rando;
}

void Utils::SetSeed(int seed)
{
    _seed = seed;
    _rando->SetSeed(_seed);
}

bool Utils::hasOriginInTracker(TVector3 spoint)
{
    //TPC Volume radius 2600 mm
    //TPC full length 5 m
    bool hasOriginInTracker = true;

    float r_point = std::sqrt(spoint.Y()*spoint.Y() + spoint.Z()*spoint.Z());
    //in the Barrel
    if( r_point > 260 ) hasOriginInTracker = false;
    //in the Endcap
    if( r_point < 260 && std::abs(spoint.X()) > 250 ) hasOriginInTracker = false;

    return hasOriginInTracker;
}

bool Utils::hasDecayedInCalo(TVector3 epoint)
{
    return !hasOriginInTracker(epoint);
}

bool Utils::isBackscatter(TVector3 spoint, TVector3 epoint)
{
    bool isBackscatter = false;

    //check if started in the calo but made it to the tracker
    float r_spoint = std::sqrt(spoint.Y()*spoint.Y() + spoint.Z()*spoint.Z());
    float r_epoint = std::sqrt(epoint.Y()*epoint.Y() + epoint.Z()*epoint.Z());

    //in the Barrel
    if( r_spoint > 260 && r_epoint < 260 ) isBackscatter = true;
    //in the Endcap
    if( (r_spoint < 260 && r_epoint < 260) && ( std::abs(spoint.X()) > 250 && std::abs(epoint.X()) < 250 ) ) isBackscatter = true;

    return isBackscatter;
}

bool Utils::isBremsstrahlung(TVector3 spoint, int pdg, int motherpdg)
{
    bool isBremsstrahlung = false;

    //Check if it has origin in the tracker and that the pdg is photon and mother is electron/positron (most probable)
    if(hasOriginInTracker(spoint) && pdg == 22 && std::abs(motherpdg) == 11) isBremsstrahlung = true;

    return isBremsstrahlung;
}
