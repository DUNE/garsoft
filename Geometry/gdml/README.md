Description of the gdml files

MPD Volumes by descending order are:

- volWorld
---- rockBox_lv (if in Hall)
---- volDetEnclosure
------- volMPD
--------- volNDHPgTPC
------------ volGArTPC
--------------- volTPCChamber
------------------ volTPCGas
------------ volPVBarrel
------------ volPVEndcap
------------ volBarrelECal
--------------- BarrelECal_stave%i_module%i_vol (staves 1 to 8 and modules 1 to 5)
------------ volEndcapECal
--------------- EndcapECal_stave%i_module%i_vol (staves 1 to 4 and modules 0 and 6)
--------- volMagnet
------------ volMagnet_Coil%i (if Coils)
--------- volYokeBarrel (if Yoke)
------------ YokeBarrel_stave%i_module%i (staves 1 to 8 and modules 1 and 2)
--------- volYokeEndcap (if Yoke)

=============================================================

MPD only in world volume

- MPD_ECalOctagon_60l_UniformMagnet.gdml

Baseline MPD with HPgTPC and PV, an ECAL in Octagonal shape (8 sides) with 60 layers in Barrel and Endcap of 2 mm Cu, 5 mm Sc and a uniform cylinder magnet of Aluminium of around 100t.

- MPD_ECalDodecagon_80l_UniformMagnet.gdml

MPD with HPgTPC and PV, an ECAL in Dodecagonal shape (12 sides) with 80 layers in Barrel and 60 layers in the Endcap of 2 mm Cu, 5 mm Sc and a uniform cylinder magnet of Aluminium of around 100t.

=============================================================

MPD with LAr in ND Hall

- MPD_LAr_Hall_2CoilsMagnet_NoYoke.gdml

LArTPC and Baseline MPD with HPgTPC and PV, an ECAL in Octagonal shape (8 sides) with 60 layers in Barrel and Endcap of 2 mm Cu, 5 mm Sc and a 2 Helmholtz coil configuration w/o PRY

- MPD_LAr_Hall_4CoilsMagnet.gdml

LArTPC and Baseline MPD with HPgTPC and PV, an ECAL in Octagonal shape (8 sides) with 60 layers in Barrel and Endcap of 2 mm Cu, 5 mm Sc and a 4 Helmholtz coil configuration (removed central coil)

- MPD_LAr_Hall_5CoilsMagnet.gdml

LArTPC and Baseline MPD with HPgTPC and PV, an ECAL in Octagonal shape (8 sides) with 60 layers in Barrel and Endcap of 2 mm Cu, 5 mm Sc and a 5 Helmholtz coil configuration

- MPD_LAr_Hall_SPYMagnet.gdml

LArTPC and Baseline MPD with HPgTPC and PV, an ECAL in Octagonal shape (8 sides) with 60 layers in Barrel and Endcap of 2 mm Cu, 5 mm Sc and the SPY Magnet configuration and a PRY according to the ND WS at DESY (https://indico.fnal.gov/event/21340/session/5/contribution/40/material/slides/0.pdf)

- MPD_LAr_Hall_UniformMagnet.gdml

LArTPC and Baseline MPD with HPgTPC and PV, an ECAL in Octagonal shape (8 sides) with 60 layers in Barrel and Endcap of 2 mm Cu, 5 mm Sc and a uniform cylindrical magnet of Aluminium of around 100t

=============================================================

=============================================================

Full ND Complex

- ND_Baseline_Hall.gdml

Full ND Complex baseline in the hall

- ND_onlyMPD_Hall.gdml

Only MPD baseline in the hall

=============================================================
