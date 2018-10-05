///
/// \file    SimulationDrawer.cxx
/// \brief   Render the objects from the Simulation package
/// \author  messier@indiana.edu
/// \version $Id: SimulationDrawer.cxx,v 1.2 2010/11/11 18:11:22 p-novaart Exp $
///

#include <iomanip>
#include <algorithm>

#include "TParticle.h"
#include "TLatex.h"
#include "TPolyMarker3D.h"
#include "TPolyMarker.h"
#include "TPolyLine.h"
#include "TPolyLine3D.h"
#include "TDatabasePDG.h"

#include "EventDisplay/SimulationDrawer.h"
#include "nutools/EventDisplayBase/View2D.h"
#include "nutools/EventDisplayBase/View3D.h"
#include "CoreUtils/ServiceUtil.h"
#include "Geometry/Geometry.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "EventDisplay/Style.h"
#include "EventDisplay/SimulationDrawingOptions.h"
#include "EventDisplay/RawDrawingOptions.h"
#include "DetectorInfo/DetectorPropertiesService.h"
#include "DetectorInfo/DetectorProperties.h"
#include "DetectorInfo/DetectorClocksService.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/View.h"
#include "art/Framework/Principal/Event.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace {
  // Utility function to make uniform error messages.
  void writeErrMsg(const char* fcn,
		   cet::exception const& e)
  {
    mf::LogWarning("SimulationDrawer") << "SimulationDrawer::" << fcn
				       << " failed with message:\n"
				       << e;
  }
}

namespace gar {
namespace evd{

  SimulationDrawer::SimulationDrawer() 
  {
    art::ServiceHandle<geo::Geometry> geom;
    minx =  1e9;
    maxx = -1e9;
    miny =  1e9;
    maxy = -1e9;
    minz =  1e9;
    maxz = -1e9;
    double world[3] = {geom->TPCXCent(),geom->TPCYCent(),geom->TPCZCent()};

    if (minx>world[0]-geom->DetHalfWidth())
      minx = world[0]-geom->DetHalfWidth();
    if (maxx<world[0]+geom->DetHalfWidth())
      maxx = world[0]+geom->DetHalfWidth();
    if (miny>world[1]-geom->DetHalfHeight())
      miny = world[1]-geom->DetHalfHeight();
    if (maxy<world[1]+geom->DetHalfHeight())
      maxy = world[1]+geom->DetHalfHeight();
    if (minz>world[2]-geom->DetLength()/2.)
      minz = world[2]-geom->DetLength()/2.;
    if (maxz<world[2]+geom->DetLength()/2.)
      maxz = world[2]+geom->DetLength()/2.;
  }

  //......................................................................

  SimulationDrawer::~SimulationDrawer()
  {
  }

  //......................................................................
  void SimulationDrawer::MCTruthShortText(const art::Event& evt,
                                          evdb::View2D*     view)
  {

    if( evt.isRealData() ) return;

    art::ServiceHandle<evd::SimulationDrawingOptions> drawopt;
    // Skip drawing if option is turned off
    if (!drawopt->fShowMCTruthText) return;

    std::vector<const simb::MCTruth*> mctruth;
    this->GetMCTruth(evt, mctruth);
    
    for (unsigned int i = 0; i<mctruth.size(); ++i) {
        std::string mctext;
        bool firstin  = true;
        bool firstout = true;
        std::string origin;
        std::string incoming;
        std::string outgoing;
        // Label cosmic rays -- others are pretty obvious
        if (mctruth[i]->Origin()==simb::kCosmicRay)  origin = "c-ray: ";
        int jmax = TMath::Min(20,mctruth[i]->NParticles());
        for (int j=0; j<jmax; ++j) {
            const simb::MCParticle& p = mctruth[i]->GetParticle(j);
            char buff[1024];
            if (p.P()>0.05) {
                sprintf(buff,"#color[%d]{%s #scale[0.75]{[%.1f GeV/c]}}",
                        Style::ColorFromPDG(p.PdgCode()),
                        Style::LatexName(p.PdgCode()),
                        p.P());
            }
            else {
                sprintf(buff,"#color[%d]{%s}",
                        Style::ColorFromPDG(p.PdgCode()),
                        Style::LatexName(p.PdgCode()));
            }
            if (p.StatusCode()==0) {
                if (firstin==false) incoming += " + ";
                incoming += buff;
                firstin = false;
            }
            if (p.StatusCode()==1) {
                if (firstout==false) outgoing += " + ";
                outgoing += buff;
                firstout = false;
            }
        } // loop on j particles
        if (origin=="" && incoming=="") {
            mctext = outgoing;
        }
        else {
            mctext = origin+incoming+" #rightarrow "+outgoing;
        }
        TLatex& latex = view->AddLatex(0.03, 0.2, mctext.c_str());
        latex.SetTextSize(0.6);

    } // loop on i mctruth objects
  }

  //......................................................................
  void SimulationDrawer::MCTruthLongText(const art::Event& evt,
                                         evdb::View2D* /*view*/)
  {
    if( evt.isRealData() ) return;
    
    art::ServiceHandle<evd::SimulationDrawingOptions> drawopt;
    // Skip drawing if option is turned off
    if (!drawopt->fShowMCTruthText) return;
    
    std::vector<const simb::MCTruth*> mctruth;
    this->GetMCTruth(evt, mctruth);
    std::cout<<"\nMCTruth Ptcl trackID            PDG      P      T   Moth  Process\n";
    for (unsigned int i=0; i<mctruth.size(); ++i) {
      for (int j=0; j<mctruth[i]->NParticles(); ++j) {
        const simb::MCParticle& p = mctruth[i]->GetParticle(j);
        if(p.StatusCode() == 0 || p.StatusCode() == 1) {
          int KE = 1000 * (p.E() - p.Mass());
          std::cout<<std::right<<std::setw(7)<<i<<std::setw(5)<<j
          <<std::setw(8)<<p.TrackId()
          <<" "<<std::setw(14)<<Style::LatexName(p.PdgCode())
          <<std::setw(7)<<int(1000 * p.P())
          <<std::setw(7)<<KE<<std::setw(7)<<p.Mother()
          <<" "<<p.Process()
          <<"\n";
        }
      } // loop on j particles in list
    }
    std::cout<<"Note: Momentum, P, and kinetic energy, T, in MeV/c\n";
  } // MCTruthLongText


  //......................................................................
  //this method draws the true particle trajectories in 3D
  void SimulationDrawer::MCTruth3D(const art::Event& evt,
                                   evdb::View3D*     view)
  {
    if( evt.isRealData() ) return;

    art::ServiceHandle<evd::SimulationDrawingOptions> drawopt;
    // If the option is turned off, there's nothing to do
    if (!drawopt->fShowMCTruthTrajectories) return;

    art::ServiceHandle<geo::Geometry> geo;
    double xcent = geo->TPCXCent();
    double ycent = geo->TPCYCent();
    double zcent = geo->TPCZCent();
    xcent = 0;
    ycent = 0;
    zcent = 0;

    // get the particles from the Geant4 step
    std::vector<const simb::MCParticle*> plist;
    this->GetParticle(evt, plist);
      
    // Useful variables
    double xMinimum(-1.*(maxx-minx));
    double xMaximum( 2.*(maxx-minx));
      
    // Define a couple of colors for neutrals and if we gray it out...
    int neutralColor(12);
    int grayedColor(15);
    int neutrinoColor(38);
    
    // Using the voxel information can be slow (see previous implementation of this code).
    // In order to speed things up we have modified the strategy:
    // 1) Make one pass through the list of voxels
    // 2) For each voxel, keep track of the MCParticle contributing energy to it and it's position
    //    which is done by keeping a map between the MCParticle and a vector of positions
    // 3) Then loop through the map to draw the particle trajectories.
    // One caveat is the need for MCParticles... and the voxels contain the track ids. So we'll need one
    // more loop to make a map of track id's and MCParticles.

    // First up is to build the map between track id's and associated MCParticles so we can recover when looping over voxels
    std::map<int, const simb::MCParticle*> trackToMcParticleMap;
      
    // Should we display the trajectories too?
    double minPartEnergy(0.01);
      
    for(size_t p = 0; p < plist.size(); ++p){
        trackToMcParticleMap[plist[p]->TrackId()] = plist[p];
        
        // Quick loop through to drawn trajectories...
        if (drawopt->fShowMCTruthTrajectories){
          // Is there an associated McTrajectory?
          const simb::MCParticle*   mcPart = plist[p];
          const simb::MCTrajectory& mcTraj = mcPart->Trajectory();
          
          int           pdgCode(mcPart->PdgCode());
          int           colorIdx(evd::Style::ColorFromPDG(mcPart->PdgCode()));
          TParticlePDG* partPDG(TDatabasePDG::Instance()->GetParticle(pdgCode));
          double        partCharge = partPDG ? partPDG->Charge() : 0.;
          double        partEnergy = mcPart->E();
          
          if (!drawopt->fShowMCTruthColors) colorIdx = grayedColor;
          
          if (!mcTraj.empty() && partEnergy > minPartEnergy && mcPart->TrackId() < 100000000){
              // The following is meant to get the correct offset for drawing the particle trajectory
              // In particular, the cosmic rays will not be correctly placed without this
            
              // collect the points from this particle
            int numTrajPoints = mcTraj.size();
            
            std::unique_ptr<double[]> hitPositions(new double[3*numTrajPoints]);
            int                       hitCount(0);
            
            for(int hitIdx = 0; hitIdx < numTrajPoints; hitIdx++){
              double xPos = mcTraj.X(hitIdx) + xcent;
              double yPos = mcTraj.Y(hitIdx) + ycent;
              double zPos = mcTraj.Z(hitIdx) + zcent;
              
                // If the original simulated hit did not occur in the TPC volume then don't draw it
              if (xPos < minx || xPos > maxx || yPos < miny || yPos > maxy|| zPos < minz || zPos > maxz) continue;
              
                // Check fiducial limits
              if (xPos > xMinimum && xPos < xMaximum){
                hitPositions[3*hitCount    ] = xPos;
                hitPositions[3*hitCount + 1] = yPos;
                hitPositions[3*hitCount + 2] = zPos;
                hitCount++;
              }
            }
            
            TPolyLine3D& pl(view->AddPolyLine3D(1, colorIdx, 1, 1));
            
              // Draw neutrals as a gray dotted line to help fade into background a bit...
            if (partCharge == 0.){
              pl.SetLineColor(neutralColor);
              pl.SetLineStyle(3);
              pl.SetLineWidth(1);
            }
            pl.SetPolyLine(hitCount, hitPositions.get(), "");
          }
        }
    }
    
    // Now we set up and build the map between MCParticles and a vector of positions obtained from the voxels
    std::map<const simb::MCParticle*, std::vector<std::vector<double> > > partToPosMap;
      
    // Finally ready for the main event! Simply loop through the map between MCParticle and positions to
    // draw the trajectories
    std::map<const simb::MCParticle*, std::vector<std::vector<double> > >::iterator partToPosMapItr;
      
    for(partToPosMapItr = partToPosMap.begin(); partToPosMapItr != partToPosMap.end(); partToPosMapItr++)
    {
        // Recover the McParticle, we'll need to access several data members so may as well dereference it
        const simb::MCParticle* mcPart = partToPosMapItr->first;
        
        // Apparently, it can happen that we get a null pointer here or maybe no points to plot
        if (!mcPart || partToPosMapItr->second.empty()) continue;
        
        // The following is meant to get the correct offset for drawing the particle trajectory
        // In particular, the cosmic rays will not be correctly placed without this
      
        int colorIdx(evd::Style::ColorFromPDG(mcPart->PdgCode()));
        int markerIdx(kFullDotSmall);
        int markerSize(2);
        
        if (!drawopt->fShowMCTruthFullSize){
            colorIdx   = grayedColor;
            markerIdx  = kDot;
            markerSize = 1;
        }
        
        TPolyMarker3D& pm = view->AddPolyMarker3D(partToPosMapItr->second.size(), colorIdx, markerIdx, markerSize); 
        
        // Now loop over points and add to trajectory
        for(size_t posIdx = 0; posIdx < partToPosMapItr->second.size(); posIdx++){
            const std::vector<double>& posVec = partToPosMapItr->second[posIdx];
            
            double xCoord = posVec[0];
            
            if (xCoord > xMinimum && xCoord < xMaximum)
                pm.SetPoint(posIdx, xCoord, posVec[1], posVec[2]);
        }
    }
      
    // Finally, let's see if we can draw the incoming particle from the MCTruth information
    std::vector<const simb::MCTruth*> mctruth;
    this->GetMCTruth(evt, mctruth);

    // Loop through the MCTruth vector
    for (unsigned int idx = 0; idx < mctruth.size(); idx++){
        // Go through each MCTruth object in the list
      for (int particleIdx = 0; particleIdx < mctruth[idx]->NParticles(); particleIdx++){
        const simb::MCParticle& mcPart = mctruth[idx]->GetParticle(particleIdx);
        
          // A negative mother id indicates the "primary" particle
        if(mcPart.Mother() == -1 && mcPart.StatusCode() == 0){
          mf::LogDebug("SimulationDrawer") << mcPart << std::endl;
          
            // Get position vector
          TVector3 particlePosition(mcPart.Vx(),mcPart.Vy(),mcPart.Vz());
          
            // Get direction vector (in opposite direction)
          TVector3 oppPartDir(-mcPart.Px(),-mcPart.Py(),-mcPart.Pz());
          
          if (oppPartDir.Mag2() > 0.) oppPartDir.SetMag(1.);
          
          double arcLenToDraw = -particlePosition.Z() / oppPartDir.CosTheta();
          
            // No point in drawing if arc length is zero (e.g. Ar nucleus)
          if (arcLenToDraw > 0.){
              // Draw the line, use an off color to be unique
            TPolyLine3D& pl(view->AddPolyLine3D(2, neutrinoColor, 1, 2));
            
            pl.SetPoint(0,particlePosition.X()+xcent,particlePosition.Y()+ycent,particlePosition.Z()+zcent);
            
            particlePosition += std::min(arcLenToDraw + 10.,1000.) * oppPartDir;
            
            pl.SetPoint(1,particlePosition.X()+xcent,particlePosition.Y()+ycent,particlePosition.Z()+zcent);
          }
        }
          // The particles we want to draw will be early in the list so break out if we didn't find them
        else break;
      } // loop on particles in list
    }
    
    return;
  }

  //......................................................................
  int SimulationDrawer::GetParticle(const art::Event&                     evt,
                                    std::vector<const simb::MCParticle*>& plist)
  {
    plist.clear();

    if( evt.isRealData() ) return 0;

    art::ServiceHandle<evd::SimulationDrawingOptions> drawopt;

    std::vector<const simb::MCParticle*> temp;

    art::View<simb::MCParticle> plcol;
    // use get by Type because there should only be one collection of these in the event
    try{
      evt.getView(drawopt->fG4ModuleLabel, plcol);
      for(unsigned int i = 0; i < plcol.vals().size(); ++i){
        temp.push_back(plcol.vals().at(i));
      }
      temp.swap(plist);
    }
    catch(cet::exception& e){
      writeErrMsg("GetRawDigits", e);
    }
  
    return plist.size();

  }

  //......................................................................
  int SimulationDrawer::GetMCTruth(const art::Event& evt,
                                   std::vector<const simb::MCTruth*>& mcvec)
  {
    mcvec.clear();

    if( evt.isRealData() ) return 0;

    std::vector<const simb::MCTruth*> temp;

    std::vector< art::Handle< std::vector<simb::MCTruth> > > mctcol;
    // use get by Type because there should only be one collection of these in the event
    try{
      evt.getManyByType(mctcol);
      for(size_t mctc = 0; mctc < mctcol.size(); ++mctc){
        art::Handle< std::vector<simb::MCTruth> > mclistHandle = mctcol[mctc];
        
        for(size_t i = 0; i < mclistHandle->size(); ++i){
          temp.push_back(&(mclistHandle->at(i)));
        }
      }
      temp.swap(mcvec);
    }
    catch(cet::exception& e){
      writeErrMsg("GetMCTruth", e);
    }
  
    return mcvec.size();
  }

  //......................................................................
  void SimulationDrawer::HiLite(int trkId, bool dohilite)
  {
    fHighlite[trkId] = dohilite;
  }
}
}// namespace
////////////////////////////////////////////////////////////////////////
