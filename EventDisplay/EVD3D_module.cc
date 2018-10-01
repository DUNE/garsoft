//
// - ROOT-based 3D event display for ART toy experiment.  Requires
//   EvtDisplayUtils, NavState, and EvtDisplayService.
//

// art includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Utilities/make_tool.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/View.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"

// ROOT includes
// ... libCore
#include <TApplication.h>
#include <TString.h>
#include <TSystem.h>
#include <TROOT.h>
// ... libRIO
#include <TFile.h>
// ... libGui
#include <TGString.h>
#include <TGLabel.h>
#include <TGButton.h>
#include <TGButtonGroup.h>
#include <TGTextEntry.h>
#include <TGTextView.h>
#include <TGLayout.h>
#include <TGTab.h>
#include <TG3DLine.h>
// ... libGeom
#include <TGeoManager.h>
#include <TGeoTube.h>
#include <TGeoCompositeShape.h>
#include <TGeoBoolNode.h>
// ... libEG
#include <TParticle.h>
// ... libRGL
#include <TGLViewer.h>
// ... libEve
#include <TEveManager.h>
#include <TEveEventManager.h>
#include <TEveBrowser.h>
#include <TEveGeoNode.h>
#include <TEveViewer.h>
#include <TEveScene.h>
#include <TEveProjectionManager.h>
#include <TEveProjectionAxes.h>
#include <TEvePointSet.h>
#include <TEveTrack.h>
#include <TEveTrackPropagator.h>
#include <TEveLine.h>

#include "nutools/EventDisplayBase/NavState.h"
#include "Geometry/Geometry.h"
#include "EventDisplay3DUtils.h"
#include "EventDisplay/Style.h"
#include "CoreUtils/ServiceUtil.h"
#include "DetectorInfo/DetectorPropertiesService.h"
#include "DetectorInfo/DetectorClocksService.h"
#include "DetectorInfo/GArPropertiesService.h"
#include "DetectorInfo/ECALPropertiesService.h"

#include "RawDataProducts/CaloRawDigit.h"
#include "RawDataProducts/RawDigit.h"
#include "RawDataProducts/raw.h"

#include "ReconstructionDataProducts/Hit.h"
#include "ReconstructionDataProducts/CaloHit.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include "SimulationDataProducts/CaloDeposit.h"
#include "SimulationDataProducts/EnergyDeposit.h"

#include "cetlib_except/demangle.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Utilities/Exception.h"

#include "TEveManager.h"
#include "TEveEventManager.h"
#include "TEveScene.h"
#include "TEveViewer.h"
#include <TEveBrowser.h>
#include "TEveRGBAPaletteOverlay.h"
#include "TEveFrameBox.h"
#include "TEveQuadSet.h"
#include "TEveTrans.h"
#include "TEveProjectionAxes.h"
#include "TEveProjectionManager.h"
#include "TEveWindow.h"
#include "TGLViewer.h"
#include "TGLCameraOverlay.h"
#include "TGTab.h"
#include "TRandom.h"

#include <sstream>
#include <iostream>
#include <iomanip>

#include "TDatabasePDG.h"
#include "TParticle.h"

namespace {
  // Utility function to make uniform error messages.
  void writeErrMsg(const char* fcn,
    cet::exception const& e)
    {
      mf::LogWarning("EventDisplay3D") << "EventDisplay3D::" << fcn
      << " failed with message:\n"
      << e;
    }
  }

  namespace gar{
    namespace evd
    {
      class EventDisplay3D : public art::EDAnalyzer
      {
      public:

        explicit EventDisplay3D(fhicl::ParameterSet const& pset);
        virtual ~EventDisplay3D();

        void beginJob() override;
        void endJob()   override;
        void beginRun(const art::Run& run)   override;
        void analyze(const art::Event& event) override;

      private:

        // Keep track of fhicl parameter set
        const fhicl::ParameterSet fParamSet;

        // Set by parameter set variables.
        bool fDrawSimHits;
        bool fDrawRawHits;
        bool fDrawRecoHits;
        bool fDrawMCTruth;

        const gar::geo::GeometryCore* fGeometry = gar::providerFrom<geo::Geometry>();
        const gar::detinfo::DetectorClocks* fTime = gar::providerFrom<detinfo::DetectorClocksService>();
        const gar::detinfo::DetectorProperties* fDetProp = gar::providerFrom<detinfo::DetectorPropertiesService>();

        std::unique_ptr<evd::EventDisplay3DUtils> fEvtDisplayUtil;

        TEveManager*           fEve;
        TEveGeoTopNode*        fEveGeoTopNode;

        TGTextEntry      *fTeRun,*fTeEvt;
        TGLabel          *fTlRun,*fTlEvt;

        TEveElementList* fCaloSimHitList;
        TEveElementList* fTPCSimHitList;
        TEveElementList* fCaloRawHitList;
        TEveElementList* fTPCRawHitList;
        TEveElementList* fCaloRecoHitList;
        TEveElementList* fTPCRecoHitList;
        TEveElementList* fTrajectoryList;

        void makeNavPanel();
        void getSimHits(std::vector<const sdp::EnergyDeposit*> &simTPC, std::vector<const sdp::CaloDeposit*> &simCalo, const art::Event& event);
        void getRawHits(std::vector<const raw::RawDigit*> &digitTPC, std::vector<const raw::CaloRawDigit*> &digitCalo, const art::Event& event);
        void getRecoHits(std::vector<const rec::Hit*> &recoTPC, std::vector<const rec::CaloHit*> &recoCalo, const art::Event& event);
        void getMCTruth(std::vector<const simb::MCTruth*>& mcvec, const art::Event& event);
        void getParticle(std::vector<const simb::MCParticle*>& plist, const art::Event& event);
        void cleanEvt();
      };

      EventDisplay3D::EventDisplay3D(fhicl::ParameterSet const& pset):
      art::EDAnalyzer(pset),
      fDrawSimHits ( pset.get<bool> ("drawSimHits", true) ),
      fDrawRawHits ( pset.get<bool> ("drawRawHits", false) ),
      fDrawRecoHits ( pset.get<bool> ("drawRecoHits", false) ),
      fDrawMCTruth ( pset.get<bool> ("drawMCTruth", true) ),
      fEve(0),
      fEveGeoTopNode(0),
      fTeRun(0),
      fTeEvt(0),
      fTlRun(0),
      fTlEvt(0),
      fCaloSimHitList(0),
      fTPCSimHitList(0),
      fCaloRawHitList(0),
      fTPCRawHitList(0),
      fCaloRecoHitList(0),
      fTPCRecoHitList(0),
      fTrajectoryList(0)
      {
        fEvtDisplayUtil = std::make_unique<evd::EventDisplay3DUtils>();
      }

      //----------------------------------------------------
      EventDisplay3D::~EventDisplay3D()
      {
      }

      void EventDisplay3D::makeNavPanel()
      {
        // Create control panel for event navigation
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        TEveBrowser* browser = fEve->GetBrowser();
        browser->StartEmbedding(TRootBrowser::kLeft); // insert nav frame as new tab in left pane

        TGMainFrame* frmMain = new TGMainFrame(gClient->GetRoot(), 1000, 600);
        frmMain->SetWindowName("EVT NAV");
        frmMain->SetCleanup(kDeepCleanup);

        TGHorizontalFrame* navFrame = new TGHorizontalFrame(frmMain);
        TGVerticalFrame* evtidFrame = new TGVerticalFrame(frmMain);

        TString icondir(TString::Format("%s/icons/", gSystem->Getenv("ROOTSYS")) );
        TGPictureButton* b = 0;

        // ... Create back button and connect to "PrevEvent" rcvr in visutils
        b = new TGPictureButton(navFrame, gClient->GetPicture(icondir + "GoBack.gif"));
        navFrame->AddFrame(b);
        b->Connect("Clicked()", "gar::evd::EventDisplay3DUtils", fEvtDisplayUtil.get(), "PrevEvent()");

        // ... Create forward button and connect to "NextEvent" rcvr in visutils
        b = new TGPictureButton(navFrame, gClient->GetPicture(icondir + "GoForward.gif"));
        navFrame->AddFrame(b);
        b->Connect("Clicked()", "gar::evd::EventDisplay3DUtils", fEvtDisplayUtil.get(), "NextEvent()");

        // ... Create run num text entry widget and connect to "GotoEvent" rcvr in visutils
        TGHorizontalFrame* runoFrame = new TGHorizontalFrame(evtidFrame);
        fTlRun = new TGLabel(runoFrame,"Run Number");
        fTlRun->SetTextJustify(kTextLeft);
        fTlRun->SetMargins(5,5,5,0);
        runoFrame->AddFrame(fTlRun);

        fTeRun = new TGTextEntry(runoFrame, fEvtDisplayUtil->fTbRun = new TGTextBuffer(5), 1);
        fEvtDisplayUtil->fTbRun->AddText(0, "1");
        fTeRun->Connect("ReturnPressed()","gar::evd::EventDisplay3DUtils", fEvtDisplayUtil.get(),"GotoEvent()");
        runoFrame->AddFrame(fTeRun,new TGLayoutHints(kLHintsExpandX));

        // ... Create evt num text entry widget and connect to "GotoEvent" rcvr in visutils
        TGHorizontalFrame* evnoFrame = new TGHorizontalFrame(evtidFrame);
        fTlEvt = new TGLabel(evnoFrame,"Evt Number");
        fTlEvt->SetTextJustify(kTextLeft);
        fTlEvt->SetMargins(5,5,5,0);
        evnoFrame->AddFrame(fTlEvt);

        fTeEvt = new TGTextEntry(evnoFrame, fEvtDisplayUtil->fTbEvt = new TGTextBuffer(5), 1);
        fEvtDisplayUtil->fTbEvt->AddText(0, "1");
        fTeEvt->Connect("ReturnPressed()","gar::evd::EventDisplay3DUtils", fEvtDisplayUtil.get(),"GotoEvent()");
        evnoFrame->AddFrame(fTeEvt,new TGLayoutHints(kLHintsExpandX));

        // ... Add horizontal run & event number subframes to vertical evtidFrame
        evtidFrame->AddFrame(runoFrame,new TGLayoutHints(kLHintsExpandX));
        evtidFrame->AddFrame(evnoFrame,new TGLayoutHints(kLHintsExpandX));

        // ... Add navFrame and evtidFrame to MainFrame
        frmMain->AddFrame(navFrame);
        TGHorizontal3DLine *separator = new TGHorizontal3DLine(frmMain);
        frmMain->AddFrame(separator, new TGLayoutHints(kLHintsExpandX));
        frmMain->AddFrame(evtidFrame);

        frmMain->MapSubwindows();
        frmMain->Resize();
        frmMain->MapWindow();

        browser->StopEmbedding();
        browser->SetTabTitle("Event Nav", 0);
      }


      void EventDisplay3D::beginJob()
      {
        // Initialize global Eve application manager (return gEve)
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        gROOT->SetBatch(kFALSE);
        fEve = TEveManager::Create();

        makeNavPanel();

        TGLViewer* glViewer = gEve->GetDefaultGLViewer();
        glViewer->SetCurrentCamera(TGLViewer::kCameraPerspXOZ);
        glViewer->SetGuideState(TGLUtil::kAxesEdge,kTRUE,kFALSE,0);
        glViewer->SetDrawCameraCenter(kTRUE);

        // Recover the root geometry
        TGeoManager* rootGeoManager = fGeometry->ROOTGeoManager();
        rootGeoManager->DefaultColors();
        auto top = rootGeoManager->GetTopNode();
        auto det = new TEveGeoTopNode(rootGeoManager, top);

        det->SetVisLevel(3);
        fEve->AddGlobalElement(det);
        fEve->Redraw3D(kTRUE);

        std::cout << "Event viewer successfully drawn" << std::endl;

        return;
      }

      void EventDisplay3D::beginRun(const art::Run& run)
      {
        return;
      }

      void EventDisplay3D::analyze(const art::Event& event)
      {
        this->cleanEvt();

        // ... Update the run and event numbers in the TGTextEntry widgets in the Navigation panel
        std::ostringstream sstr;
        sstr << event.id().run();

        fEvtDisplayUtil->fTbRun->Clear();
        fEvtDisplayUtil->fTbRun->AddText(0, sstr.str().c_str());
        gClient->NeedRedraw(fTeRun);

        sstr.str("");
        sstr << event.id().event();
        fEvtDisplayUtil->fTbEvt->Clear();
        fEvtDisplayUtil->fTbEvt->AddText(0, sstr.str().c_str());
        gClient->NeedRedraw(fTeEvt);

        if(fDrawSimHits)
        {
          //Get Sim hits
          std::vector<const sdp::EnergyDeposit*> simTPC;
          std::vector<const sdp::CaloDeposit*> simECAL;
          this->getSimHits(simTPC, simECAL, event);

          //TPC
          TEveElementList* TPChitList = new TEveElementList("TPC", "Energy Deposit for TPC");
          fTPCSimHitList = new TEveElementList("SimTPCHit", "Sim TPC hits");
          fTPCSimHitList->SetMainColor(kRed);
          fTPCSimHitList->SetMainAlpha(1.0);

          for(auto hit: simTPC)
          {
            const double energy = hit->Energy();
            const double x = hit->X();
            const double y = hit->Y();
            const double z = hit->Z();

            TEveLine* eveHit = new TEveLine(2);
            eveHit->SetName("TPCSimHit");
            std::ostringstream title;
            title << "TPC Hit ";
            title << std::fixed << std::setprecision(2)
            << " " << energy << " GeV";
            title << " at (" << x << " cm"
            << "," <<  y << " cm"
            << "," <<  z << " cm"
            << ")";

            eveHit->SetTitle(title.str().c_str());
            eveHit->SetLineWidth(5);
            eveHit->SetLineColor(fEvtDisplayUtil->LogColor(energy, 0, 10, 3));
            eveHit->SetPoint(0, x, y, z);//cm
            eveHit->SetPoint(1, x+1, y+1, z+1);//cm
            TPChitList->AddElement(eveHit);
          }

          fTPCSimHitList->AddElement(TPChitList);
          fEve->AddElement(fTPCSimHitList);

          //ECAL
          TEveElementList* hitList = new TEveElementList("ECAL", "Energy Deposit for ECAL");
          fCaloSimHitList = new TEveElementList("SimCaloHit", "Sim ECAL hits");
          fCaloSimHitList->SetMainColor(kRed);
          fCaloSimHitList->SetMainAlpha(1.0);

          for(auto hit: simECAL)
          {
            const double energy = hit->Energy();
            const double x = hit->X();
            const double y = hit->Y();
            const double z = hit->Z();

            TEveLine* eveHit = new TEveLine(2);
            eveHit->SetName("CaloSimHit");
            std::ostringstream title;
            title << "Hit ";
            title << std::fixed << std::setprecision(2)
            << " " << energy << " GeV";
            title << " at (" << x << " cm"
            << "," <<  y << " cm"
            << "," <<  z << " cm"
            << ")";

            eveHit->SetTitle(title.str().c_str());
            eveHit->SetLineWidth(5);
            eveHit->SetLineColor(fEvtDisplayUtil->LogColor(energy, 0, 10, 3));
            eveHit->SetPoint(0, x, y, z);//cm
            eveHit->SetPoint(1, x+1, y+1, z+1);//cm
            hitList->AddElement(eveHit);
          }

          fCaloSimHitList->AddElement(hitList);
          fEve->AddElement(fCaloSimHitList);
        }

        if(fDrawRawHits)
        {
          //Get Raw hits
          std::vector<const raw::RawDigit*> digitTPC;
          std::vector<const raw::CaloRawDigit*> digitECAL;
          this->getRawHits(digitTPC, digitECAL, event);

          //TPC
          TEveElementList* TPChitList = new TEveElementList("TPC", "Energy Deposit for TPC");
          fTPCRawHitList = new TEveElementList("RawTPCHit", "Raw TPC hits");
          fTPCRawHitList->SetMainColor(kRed);
          fTPCRawHitList->SetMainAlpha(1.0);

          double driftVelocity = fDetProp->DriftVelocity(fDetProp->Efield(), fDetProp->Temperature());
          gar::raw::ADCvector_t uncompressed;

          for(auto hit: digitTPC)
          {
            float xyz[3];
            fGeometry->ChannelToPosition((unsigned int)hit->Channel(), xyz);
            double chanposx = xyz[0];
            uncompressed.resize(hit->Samples());
            short ped = hit->Pedestal();
            gar::raw::Uncompress(hit->ADCs(), uncompressed, ped, hit->Compression());

            for(size_t t = 0; t < hit->Samples(); ++t)
            {
              if(uncompressed[t] < 5) continue;

              double driftdistance = fTime->TPCTick2Time(t) * driftVelocity;

              if (chanposx < 0)
              xyz[0] = chanposx + driftdistance;
              else
              xyz[0] = chanposx - driftdistance;

              TEveLine* eveHit = new TEveLine(2);
              eveHit->SetName("TPCRawHit");
              std::ostringstream title;
              title << "TPC Hit ";
              title << std::fixed << std::setprecision(2)
              << " " << uncompressed[t] << " ADC";
              title << " at (" << xyz[0] << " cm"
              << "," <<  xyz[1] << " cm"
              << "," <<  xyz[2] << " cm"
              << ")";

              eveHit->SetTitle(title.str().c_str());
              eveHit->SetLineWidth(5);
              eveHit->SetLineColor(fEvtDisplayUtil->LogColor(uncompressed[t], 0, 4096, 3));
              eveHit->SetPoint(0, xyz[0], xyz[1], xyz[2]);//cm
              eveHit->SetPoint(1, xyz[0]+1, xyz[1]+1, xyz[2]+1);//cm
              TPChitList->AddElement(eveHit);
            }
          }

          fTPCRawHitList->AddElement(TPChitList);
          fEve->AddElement(fTPCRawHitList);

          //ECAL
          TEveElementList* hitList = new TEveElementList("ECAL", "Energy Deposit for ECAL");
          fCaloRawHitList = new TEveElementList("RawCaloHit", "Raw ECAL hits");
          fCaloRawHitList->SetMainColor(kRed);
          fCaloRawHitList->SetMainAlpha(1.0);

          for(auto hit: digitECAL)
          {
            const double energy = hit->ADC();
            const double x = hit->X();
            const double y = hit->Y();
            const double z = hit->Z();

            TEveLine* eveHit = new TEveLine(2);
            eveHit->SetName("CaloRawHit");
            std::ostringstream title;
            title << "Hit ";
            title << std::fixed << std::setprecision(2)
            << " " << energy << " ADC";
            title << " at (" << x << " mm"
            << "," <<  y << " mm"
            << "," <<  z << " mm"
            << ")";

            eveHit->SetTitle(title.str().c_str());
            eveHit->SetLineWidth(5);
            eveHit->SetLineColor(fEvtDisplayUtil->LogColor(energy, 0, 4096, 3));
            eveHit->SetPoint(0, x, y, z);//cm
            eveHit->SetPoint(1, x+1, y+1, z+1);//cm
            hitList->AddElement(eveHit);
          }

          fCaloRawHitList->AddElement(hitList);
          fEve->AddElement(fCaloRawHitList);
        }

        if(fDrawRecoHits)
        {
          //Get Reco hits
          std::vector<const rec::Hit*> recoTPC;
          std::vector<const rec::CaloHit*> recoECAL;
          this->getRecoHits(recoTPC, recoECAL, event);

          //TPC
          TEveElementList* TPChitList = new TEveElementList("TPC", "Energy Deposit for TPC");
          fTPCRecoHitList = new TEveElementList("RecoTPCHit", "Reco TPC hits");
          fTPCRecoHitList->SetMainColor(kRed);
          fTPCRecoHitList->SetMainAlpha(1.0);

          for(auto itr : recoTPC)
          {
            auto const signal = itr->Signal();
            auto const* pos = itr->Position();

            TEveLine* eveHit = new TEveLine(2);
            eveHit->SetName("TPCRecoHit");
            std::ostringstream title;
            title << "Hit ";
            title << std::fixed << std::setprecision(2)
            << " " << signal << " ADC";
            title << " at (" << pos[0] << " cm"
            << "," <<  pos[1] << " cm"
            << "," <<  pos[2] << " cm"
            << ")";

            eveHit->SetTitle(title.str().c_str());
            eveHit->SetLineWidth(5);
            eveHit->SetLineColor(fEvtDisplayUtil->LogColor(signal, 0, 4096, 3));
            eveHit->SetPoint(0, pos[0], pos[1], pos[2]);//cm
            eveHit->SetPoint(1, pos[0]+1, pos[1]+1, pos[2]+1);//cm
            TPChitList->AddElement(eveHit);
          }

          fTPCRecoHitList->AddElement(TPChitList);
          fEve->AddElement(fTPCRecoHitList);

          //ECAL
          TEveElementList* CalohitList = new TEveElementList("ECAL", "Energy Deposit for ECAL");
          fCaloRecoHitList = new TEveElementList("RecoECALHit", "Reco ECAL hits");
          fCaloRecoHitList->SetMainColor(kRed);
          fCaloRecoHitList->SetMainAlpha(1.0);

          for(auto itr : recoECAL)
          {
            auto const energy = itr->Energy();
            auto const* pos = itr->Position();

            TEveLine* eveHit = new TEveLine(2);
            eveHit->SetName("CaloRecoHit");
            std::ostringstream title;
            title << "Hit ";
            title << std::fixed << std::setprecision(2)
            << " " << energy << " MeV";
            title << " at (" << pos[0] << " mm"
            << "," <<  pos[1] << " mm"
            << "," <<  pos[2] << " mm"
            << ")";

            eveHit->SetTitle(title.str().c_str());
            eveHit->SetLineWidth(5);
            eveHit->SetLineColor(fEvtDisplayUtil->LogColor(energy, 0, 10, 3));
            eveHit->SetPoint(0, pos[0], pos[1], pos[2]);//cm
            eveHit->SetPoint(1, pos[0]+1, pos[1]+1, pos[2]+1);//cm
            CalohitList->AddElement(eveHit);
          }

          fCaloRecoHitList->AddElement(CalohitList);
          fEve->AddElement(fCaloRecoHitList);
        }

        if(fDrawMCTruth)
        {
          std::vector<const simb::MCParticle*> plist;
          this->getParticle(plist, event);

          fTrajectoryList = new TEveElementList("g4Trajectories", "Geant4 Trajectories");
          fTrajectoryList->SetMainColor(kRed);
          fTrajectoryList->SetMainAlpha(1.0);

          // First up is to build the map between track id's and associated MCParticles so we can recover when looping over voxels
          std::map<int, const simb::MCParticle*> trackToMcParticleMap;

          for(size_t p = 0; p < plist.size(); ++p)
          {
            trackToMcParticleMap[plist[p]->TrackId()] = plist[p];

            // Is there an associated McTrajectory?
            const simb::MCParticle*   mcPart = plist[p];
            const simb::MCTrajectory& mcTraj = mcPart->Trajectory();

            int           pdgCode(mcPart->PdgCode());
            int           colorIdx(evd::Style::ColorFromPDG(mcPart->PdgCode()));
            double        partEnergy = mcPart->E();

            if (!mcTraj.empty() && partEnergy > 0.01 && mcPart->TrackId() < 100000000)
            {
              // collect the points from this particle
              int numTrajPoints = mcTraj.size();

              std::ostringstream label;
              label << "Particle pdg " << pdgCode;
              label << " (Energy " << partEnergy << " MeV)";

              TEveLine *track = new TEveLine();
              track->SetName("trajectory");
              track->SetTitle(label.str().c_str());
              track->SetLineColor(colorIdx);
              track->SetLineStyle(1);
              track->SetLineWidth(3);

              for(int hitIdx = 0; hitIdx < numTrajPoints; hitIdx++)
              {
                double xPos = mcTraj.X(hitIdx);
                double yPos = mcTraj.Y(hitIdx);
                double zPos = mcTraj.Z(hitIdx);

                track->SetPoint(hitIdx, xPos, yPos, zPos);
              }

              fTrajectoryList->AddElement(track);
            }
          }

          fEve->AddElement(fTrajectoryList);
        }

        return;
      } // end EventDisplay3D::EventDisplay3D::analyze

      void EventDisplay3D::getSimHits(std::vector<const sdp::EnergyDeposit*> &simTPC, std::vector<const sdp::CaloDeposit*> &simCalo, const art::Event& event)
      {
        simTPC.clear();
        auto TPChandle = event.getValidHandle< std::vector<sdp::EnergyDeposit> >("geant");
        std::vector<sdp::EnergyDeposit> const& simTPCtemp(*TPChandle);

        for (sdp::EnergyDeposit const& prod: simTPCtemp) {
          // do something with prod
          simTPC.push_back(&prod);
        }

        simCalo.clear();
        auto Calohandle = event.getValidHandle< std::vector<sdp::CaloDeposit> >("geant");
        std::vector<sdp::CaloDeposit> const& simCalotemp(*Calohandle);

        for (sdp::CaloDeposit const& prod: simCalotemp) {
          // do something with prod
          simCalo.push_back(&prod);
        }

      }

      void EventDisplay3D::getRawHits(std::vector<const raw::RawDigit*> &digitTPC, std::vector<const raw::CaloRawDigit*> &digitCalo, const art::Event& event)
      {
        digitTPC.clear();
        std::vector<const raw::RawDigit*> tempTPC;

        try
        {
          event.getView("daq", tempTPC);
          for(size_t t = 0; t < tempTPC.size(); ++t)
          digitTPC.push_back(tempTPC[t]);
        }
        catch(cet::exception& e){
          writeErrMsg("GetDigits TPC", e);
        }

        digitCalo.clear();
        std::vector<const raw::CaloRawDigit*> tempCalo;

        try
        {
          event.getView("daqecal", tempCalo);
          for(size_t t = 0; t < tempCalo.size(); ++t)
          digitCalo.push_back(tempCalo[t]);
        }
        catch(cet::exception& e){
          writeErrMsg("GetDigits Calo", e);
        }
      }

      void EventDisplay3D::getRecoHits(std::vector<const rec::Hit*> &recoTPC, std::vector<const rec::CaloHit*> &recoCalo, const art::Event& event)
      {
        recoTPC.clear();
        std::vector<const rec::Hit*> tempTPC;

        try
        {
          event.getView("hit", tempTPC);
          for(size_t t = 0; t < tempTPC.size(); ++t)
          recoTPC.push_back(tempTPC[t]);
        }
        catch(cet::exception& e){
          writeErrMsg("GetHits TPC", e);
        }

        recoCalo.clear();
        std::vector<const rec::CaloHit*> tempCalo;

        try
        {
          event.getView("calohit", tempCalo);
          for(size_t t = 0; t < tempCalo.size(); ++t)
          recoCalo.push_back(tempCalo[t]);
        }
        catch(cet::exception& e){
          writeErrMsg("GetHits Calo", e);
        }
      }

      void EventDisplay3D::getMCTruth(std::vector<const simb::MCTruth*>& mcvec, const art::Event& event)
      {
        mcvec.clear();

        if( event.isRealData() ) return;

        std::vector<const simb::MCTruth*> temp;
        std::vector< art::Handle< std::vector<simb::MCTruth> > > mctcol;
        // use get by Type because there should only be one collection of these in the event
        try
        {
          event.getManyByType(mctcol);

          for(size_t mctc = 0; mctc < mctcol.size(); ++mctc)
          {
            art::Handle< std::vector<simb::MCTruth> > mclistHandle = mctcol[mctc];
            for(size_t i = 0; i < mclistHandle->size(); ++i)
            temp.push_back(&(mclistHandle->at(i)));
          }

          temp.swap(mcvec);
        }
        catch(cet::exception& e){
          writeErrMsg("GetMCTruth", e);
        }

        return;
      }

      void EventDisplay3D::getParticle(std::vector<const simb::MCParticle*>& plist, const art::Event& event)
      {
        plist.clear();

        if( event.isRealData() ) return;

        std::vector<const simb::MCParticle*> temp;
        art::View<simb::MCParticle> plcol;
        // use get by Type because there should only be one collection of these in the event
        try
        {
          event.getView("geant", plcol);
          for(unsigned int i = 0; i < plcol.vals().size(); ++i)
          temp.push_back(plcol.vals().at(i));

          temp.swap(plist);
        }
        catch(cet::exception& e){
          writeErrMsg("GetParticle", e);
        }

        return;
      }

      void EventDisplay3D::cleanEvt()
      {
        //check if any element is present in the event and remove it
        if(fEve->GetCurrentEvent() != nullptr)
        fEve->GetCurrentEvent()->DestroyElements();
      }

      void EventDisplay3D::endJob(){

      }

      DEFINE_ART_MODULE(EventDisplay3D)

    }//end namespace evd
  }//end namespace gar
