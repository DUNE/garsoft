# garsoft: Software for ND-GAr

## Quick-Start: Using pre-installed GArSoft in CVMFS

If your computer has CVMFS installed and running, you can type

```
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup garsoft v02_20_00 -q e26:prof
```

and then proceed to the section below on running a sample of GArSoft events.  
The version in the setup command above will go out of date. To see what versions are available, type

```
ups list -aK+ garsoft
```

after sourcing the setup script. This list is not expected to be sorted.  Versions starting with v02_20_00 use the e26 (gcc v12.1.0) and the c14 (clang 14.0.6) compilers.  Most people on DUNE run production with gcc, though performance is similar.


## Building GArSoft from Source -- all other dependencies in CVMFS

To set up a new test release for doing development in your own area, execute the following commands (or put them in a shell script) to set up a new development area.  These instructions are tested to work on the current version in the head of develop. If you are working with a private test release made before September 27, 2021, you may need to use the old mrb version.  After sourcing the DUNE setup script, unsetup mrb and setup mrb -o.  Or make a new test release with the latest software.  GArSoft v02_13_00 works with mrb 5.  v02_12_00 is buggy; don't use it.

```
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
mkdir <new empty directory>
cd <new emtpy directory>
export MRB_PROJECT=garsoft
export COMPILER=e26
export BUILDTYPE=prof
mrb newDev -v develop -q ${COMPILER}:${BUILDTYPE}
source localProducts*/setup
mkdir work
cd srcs
mrb g garsoft
cd $MRB_BUILDDIR
mrbsetenv
mrb i -j4
mrbslp
```

The COMPILER variable above is a UPS qualifier indicating which version of which compiler to use. See https://cdcvs.fnal.gov/redmine/projects/cet-is-public/wiki/AboutQualifiers for the meanings of these qualifiers.

The BUILDTYPE variable above is either "prof" or "debug". Both kinds of builds include debug symbols, but "prof" turns on optimization, which can make using a debugger more challenging, but will make the code run faster.

The -j4 argument on the mrb i command refers to the number of concurrent build processes (like g++ compiling a single source file).  A good rule of thum is to match it to the number of CPU threads available on the machine you are working.  DUNE gpvms have four cores and around 10 GB of memory, so use -j4 on those.  DUNE build nodes have 16 cores each and enough memory to support running the compiler on each one.

Each time you log in and want to work with your garsoft test release, execute (or write a script), do the following to set up a build environment:

```
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
export MRB_PROJECT=garsoft
source <directory you created from the creation step above>/localProducts*/setup
cd ${MRB_BUILDDIR}
mrbsetenv
```

Note: mrbsetenv sets up a build environment, and environment variables will point to the build directory. mrbslp will set up installed local products for running. Sometimes this makes a difference -- if a build fails, the build directory can have an incomplete set of things in it. It may also be missing some items (pandora XML files is a common one in far detector sim/reco that only get picked up with mrbslp).

See the DUNE computing [training](https://dune.github.io/computing-basics/) for a tutorial on how to develop software.  Some older examples from dunetpc on how to work with git, check out code, develop on branches, and pushing code: [dunetpc git/mrb tutorial](https://cdcvs.fnal.gov/redmine/projects/dunetpc/wiki/_Tutorial_)

## Run a 1000-event MC sample 

After building or setting up code:

```
art -n 1000 -c prodgenie.fcl
art -c readoutsimjob.fcl genie_gen.root
art -c recojob.fcl readoutsim.root         (sim_readout.root for older checkouts of the fcl files)
art -c evd.fcl reco.root
```

or, insetead of the last line, use

```
art -c evd3D.fcl reco.root
```

for a nice eve-based event display Eldwan provided. Note that the prodgenie.fcl job requires flux files. Currently the location of those flux files is set to a pnfs directory at FNAL. (If you are running somewhere other than the FNAL cluster, you might want to try prodsingle.fcl instead, or you will need to edit the search path and copy method in the prodgenie.fcl file. The fcl files are in $FHICL_FILE_PATH.)  In order for art to find the flux files, you need a grid proxy or an unexpired Kerberos ticket from which ifdh can make a grid proxy.  setup_fnal_security (in duneutil) also gives you a grid proxy based on your Kerberos ticket.

You can now make an analysis root tree with

```
art -c anajob.fcl reco.root
```

There is an example root macro that reads the anatree.root file made by anajob.fcl. It's in garsoft/Ana/ExampleAnalysisScripts/garana.C and garana.h. There's no CMakeLists.txt file in there so the script doesn't get installed anywhere by default, so it just lives in the git repository. Copy it to your working directory along with the header file and run it to make some plots.

## Run TOAD Simulation and Reconstruction
These commands will generate a muon track, simulate detector readout, reconstruct and produce an event display. The toadfgen.fcl points to a event generation file that can be adapted as needed.

```
# event generation
art -n 2 -c toadtfgen.fcl
# readout simulation
art -c readoutsimjob_toad.fcl text_gen.root
# reconstruction
art -c recojob_toad.fcl readoutsim.root
# event display
art -c evd3D_toad.fcl readoutsim.root
```

DUNE-DAQ hdf5 files can be read and reconstructed as follows:

```
#convert hdf5 to root
art -c runTOADRawDecoder.fcl <hdf5 file > -o <output decoded root file name>
# reconstruction
art -c recojob_toad.fcl <output decoded root file name> -o <reco root file name>
#event display
art -c evd3D_toad.fcl <reco root file name>
```

## Tips, tricks, and gotchas

During development, testing and debugging, it is good to have three login sessions up on your screen simultaneously.

* A session used for editing code. No special environment is needed for this
* A session for compiling code. Use the setup instructions above and use mrbsetenv to set up the build environment.
* A session used for running code. Use the setup instructions above and use mrbslp to set up the run environment.

## Event Display tips
The event display starts a ROOT TApplication, and reads rootlogon.C. Users may have a rootlogon.C which is incompatible with the event display -- the symptom is a segfault starting up the event display. It works fine with a blank rootlogon.C, which you can put in the directory you're running the event display from so as not to clobber your regular rootlogon.C.

## Debugger tips
The forge_tools debugger from arm (setup forge_tools to get it) comes with customized versions of gdb and several system libraries, like libGL.so. The system library replacements, probably libGL.so, interfere with EVE, so you have to unsetup forge_tools in order to get evd3D.fcl to work (or any other EVE-based event display for that matter). I like to have a separate login session set up with the same running environment as the one the debugger is running in so I can search for source code that ddt cannot find automatically. You may have to use UPS-defined environment variables to find source code in CVMFS.

## Computing tips
Lots of computing how-to information is here: https://wiki.dunescience.org/wiki/DUNE_Computing/Computing_How-To_Documentation

## Copyright and Licensing
Copyright © 2023 FERMI NATIONAL ACCELERATOR LABORATORY for the benefit of the DUNE Collaboration.

This repository, and all software contained within, except where noted within the individual source files, is licensed under
the Apache License, Version 2.0 (the "License"); you may not use this
file except in compliance with the License. You may obtain a copy of
the License at

    http://www.apache.org/licenses/LICENSE-2.0

Copyright is granted to FERMI NATIONAL ACCELERATOR LABORATORY on behalf
of the Deep Underground Neutrino Experiment (DUNE). Unless required by
applicable law or agreed to in writing, software distributed under the
License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for
the specific language governing permissions and limitations under the
License.
