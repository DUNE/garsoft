// ROOT libraries

// cheats to get ROOT include headers to work with c14

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wsizeof-pointer-div"

#include "TROOT.h"
#include "TApplication.h"
#include "TGClient.h"
#include "TSystem.h"
#include "TVirtualX.h"
#include "TGX11.h" // this header currently triggers a compiler error,
// that we had to disable via compiler flag -Wno-variadic-macros

#pragma GCC diagnostic pop

// C/C++ standard libraries
#include <iostream>
#include <string>
#include <stdexcept> // std::runtime_error
#include <cstdlib> // getenv()

namespace { // local namespace
/**
* This class implements a work-around to initialize ROOT graphics system
* before some other code messes with it (that is, we are the first ones
* to mess with it, and of course we do it right(TM)).
* To ensure this code is executed as soon as possible, it's bound to
* the construction of a static variable. While it is not predictable
* when this code will be executed, it is expected to be executed before
* main() function of the executable if linked statically, or as soon
* as a library statically linked with this one is dynamically loaded.
*
* The work around consists of making sure there is an active TApplication
* (but not the ones that pull in an interactive prompt, e.g. TRint),
* pulling ROOT out of batch mode, and creating ROOT's X11 graphics client
* (that, for example, TEve needs).
*/
struct RootGraphicsEnablerClass
{
    /// Default constructor: quiet
    RootGraphicsEnablerClass() { EnableRootGraphics(); }

    /// Constructor: enable messages on specified stream
    RootGraphicsEnablerClass(std::ostream& out) { EnableRootGraphics(&out); }

    /// Enacts the tricks to enable the graphics
/// @param out pointer to output stream (default: none, be quiet)
static void EnableRootGraphics(std::ostream* out = nullptr);
}; // RootGraphicEnablerClass

// static instance that is here only to ensure the initialization function is called;
// this is currently verbose
RootGraphicsEnablerClass RootGraphicsEnabler{
    std::cout
};

void RootGraphicsEnablerClass::EnableRootGraphics(std::ostream* out /* = nullptr */)
{
    if (out) (*out) << "RootGraphicsEnablerClass hacking its way forth." << std::endl;

    //======================================================================
    // Setup the root environment for a program started with no arguments
    //======================================================================
    if (out) (*out) << "  ==> get the current TApplication (and make sure gROOT is valid)" << std::endl;
    TApplication* app = ROOT::GetROOT()->GetApplication();

    // ROOT::GetROOT() should initialize gROOT.
    if (!gROOT)
    throw std::runtime_error("RootGraphicsEnabler: no ROOT global pointer");

    if (!app) {
        if (out) (*out) << "  ==> create a TApplication" << std::endl;
        int    argc = 0;
        char** argv = nullptr;
        new TApplication("TApplicationFromRootGraphicsEnabler", &argc, argv);
    } // if no application

    if (out) (*out) << "  ==> set batch mode off (now it is " << (gROOT->IsBatch()? "on)": "already off)") << std::endl;

    gROOT->SetBatch(kFALSE);

    if (!gClient)
    {
        if (out) (*out) << "  ==> creating a TGClient" << std::endl;
        if (out) (*out) << "      ==> loading graphics library (X11)" << std::endl;
        int res = gSystem->Load("libGX11.so");
        if (out) {
            switch (res) {
            case  0: break; // successfully loaded
            case  1: (*out) << "          (it was already)" << std::endl; break;
            case -1: (*out) << "          ERROR: not found!" << std::endl; break;
            case -2: (*out) << "          ERROR: version mismatch!" << std::endl; break;
            default: (*out) << "          ERROR: undocumented (code=" << res << ")" << std::endl; break;
            } // switch
        }

        if (out) (*out) << "      ==> creating TVirtualX" << std::endl;

        gVirtualX = new TGX11("X11", "X11 session");
        std::string const DISPLAY = getenv("DISPLAY");

        if (out) (*out) << "      ==> creating the TGClient (DISPLAY='" << DISPLAY << "')" << std::endl;

        new TGClient(DISPLAY.c_str());
    } // if no graphics client

    if (out) (*out) << "RootGraphicsEnablerClass hacking compleled." << std::endl;

} // RootGraphicsEnabler::EnableRootGraphics()

} // local namespace
