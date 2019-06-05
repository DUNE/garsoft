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
#include <TGeoNode.h>
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
#include <TEveTrackPropagator.h>
#include <TEveLine.h>
#include "TEveTrack.h"
#include "TEveTrackPropagator.h"
#include "TEvePathMark.h"
#include "TTimeStamp.h"

#include "nutools/EventDisplayBase/NavState.h"
#include "Geometry/Geometry.h"

#include "EventDisplay3DUtils.h"
#include "EventDisplay/Style.h"
#include "EventDisplay/eventdisplay.h"

#include "CoreUtils/ServiceUtil.h"
#include "DetectorInfo/DetectorPropertiesService.h"
#include "DetectorInfo/DetectorClocksService.h"
#include "DetectorInfo/GArPropertiesService.h"
#include "DetectorInfo/ECALPropertiesService.h"

#include "Utilities/TrackPropagator.h"

#include "RawDataProducts/CaloRawDigit.h"
#include "RawDataProducts/RawDigit.h"
#include "RawDataProducts/raw.h"

#include "ReconstructionDataProducts/Hit.h"
#include "ReconstructionDataProducts/CaloHit.h"

#include "ReconstructionDataProducts/Cluster.h"
#include "ReconstructionDataProducts/Track.h"
#include "ReconstructionDataProducts/Vertex.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include "cetlib_except/demangle.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Utilities/Exception.h"

#include "TEveManager.h"
#include "TEveEventManager.h"
#include "TEveScene.h"
#include "TEveViewer.h"
#include "TEveBrowser.h"
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
#include "TEveArrow.h"
#include "TGLAnnotation.h"
#include "TGLFontManager.h"

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

            void reconfigure(fhicl::ParameterSet const& pset);

        private:

            // Keep track of fhicl parameter set
            const fhicl::ParameterSet fParamSet;

            const gar::geo::GeometryCore* fGeometry = gar::providerFrom<geo::Geometry>();
            const gar::detinfo::DetectorClocks* fTime = gar::providerFrom<detinfo::DetectorClocksService>();
            const gar::detinfo::DetectorProperties* fDetProp = gar::providerFrom<detinfo::DetectorPropertiesService>();

            //Event display 3D pointer
            std::unique_ptr<evd::EventDisplay3DUtils> fEvtDisplayUtil;

            // Set by parameter set variables.
            int fDrawECALRawHits;
            int fDrawECALRecoHits;
            int fDrawECALClusters;
            int fDrawTracks;
            int fDrawVertices;
            int fDrawMCTruth;

            std::string fG4Label;                                       ///< module label that produced G4 hits
            std::vector<std::string> fSimHitLabels;     		        ///< module labels that produced sim hits
            std::vector<std::string> fRawHitLabels;     		        ///< module labels that produced raw hits
            std::vector<std::string> fRecoHitLabels;     		        ///< module labels that produced reco hits
            std::vector<std::string> fCaloClusterLabels;    	        ///< module labels that produced Calorimeter Clusters
            std::vector<std::string> fTrackLabels;   		            ///< module labels that produced tracks
            std::vector<std::string> fVertexLabels;  		            ///< module labels that produced vertices

            std::vector<std::string> fVolumesToShow;                    ///< list of volumes to show in the Eve display

            TGeoManager*  fGeoManager;
            TEveManager*  fEve;
            TGLViewer* glViewer;
            TGLAnnotation* ann;

            TGTextEntry      *fTeRun, *fTeEvt;
            TGLabel          *fTlRun, *fTlEvt;

            //Calo specific
            TEveElementList* fCaloRawHitList;
            TEveElementList* fCaloRecoHitList;

            //High level reco (clusters and tracks)
            TEveElementList* fCaloClusterList;
            TEveElementList* fTrackList;
            TEveElementList* fVertexList;

            //MC Track and incoming particle
            TEveElementList* fMCTrajectoryList;

            //scaling factor for cluster axis
            float fScalingfactor;
            //draw track calo intersections
            int fDrawIntersection;

            //Create the navigation panel
            void makeNavPanel();

            void PickVolumes(TEveElementList* &list);

            void UpdateHeader(const art::Event& event);

            void DrawMCTruth(const art::Event& event);

            void DrawRawHits(const art::Event& event);

            void DrawRecoHits(const art::Event& event);

            void DrawHighLevelReco(const art::Event& event);

            void DrawTrack(const gar::rec::Track* trk, TEveElementList* &eve_list, int counter, bool forward);

            void DrawIntersections(const gar::rec::Track* trk, TEveElementList* &fTrackList, int counter);

            void DrawHelix3D(const float *trackpar, const float xpar, const float xother, TEveLine* &eve_track, int color);

            void DrawArrow3D(const float *fVertex, const float *fDir, int color, TEveArrow* &arrow);

            //get raw hits from the event handle
            void GetRawCaloHits(std::vector<const raw::CaloRawDigit*> &digitCalo, const art::Event& event);

            //get reco hits from the event handle
            void GetRecoCaloHits(std::vector<const rec::CaloHit*> &recoCalo, const art::Event& event);

            //get calo clusters
            void GetCaloClusters(std::vector<const rec::Cluster*> &caloCluster, const art::Event& event);

            //get tracks clusters
            void GetTracks(std::vector<const rec::Track*> &track, const art::Event& event);

            //get reco vertexes
            void GetVertices(std::vector<const rec::Vertex*> &vertex, const art::Event& event);

            //get mc truth from the event handle
            void GetMCTruth(std::vector<const simb::MCTruth*>& mcvec, const art::Event& event);

            //get particle list from the event handle
            void GetParticle(std::vector<const simb::MCParticle*>& plist, const art::Event& event);

            //clean the event before going to next
            void cleanEvt();
        };

        //......................................................................
        EventDisplay3D::EventDisplay3D(fhicl::ParameterSet const& pset):
            art::EDAnalyzer(pset),
            fEve(0),
            fTeRun(0),
            fTeEvt(0),
            fTlRun(0),
            fTlEvt(0),
            fCaloRawHitList(0),
            fCaloRecoHitList(0),
            fCaloClusterList(0),
            fTrackList(0),
            fVertexList(0),
            fMCTrajectoryList(0),
            fScalingfactor(1.0)
            {
                this->reconfigure(pset);

                fEvtDisplayUtil = std::make_unique<evd::EventDisplay3DUtils>();

                fGeoManager = fGeometry->ROOTGeoManager();
                fGeoManager->DefaultColors();
            }

            //......................................................................
            void EventDisplay3D::reconfigure(fhicl::ParameterSet const& pset)
            {
                fDrawECALRawHits           = pset.get<int                       > ("drawECALRawHits"      , 0);
                fDrawECALRecoHits          = pset.get<int                       > ("drawECALRecoHits"     , 0);
                fDrawECALClusters          = pset.get<int                       > ("drawECALClusters"     , 1);
                fDrawTracks                = pset.get<int                       > ("drawTracks"           , 0);
                fDrawVertices              = pset.get<int                       > ("drawVertices"         , 1);
                fDrawMCTruth               = pset.get<int                       > ("drawMCTruth"          , 1);
                fVolumesToShow             = pset.get< std::vector<std::string> > ("VolumesToShow"           );

                fG4Label                   = pset.get< std::string              > ("G4ModuleLabel"           );
                fRawHitLabels              = pset.get< std::vector<std::string> > ("RawHitModuleLabels"      );
                fRecoHitLabels             = pset.get< std::vector<std::string> > ("RecoHitModuleLabels"     );
                fCaloClusterLabels         = pset.get< std::vector<std::string> > ("CaloClusterModuleLabels" );
                fTrackLabels      	       = pset.get< std::vector<std::string> > ("TrackModuleLabels"     	 );
                fVertexLabels     	       = pset.get< std::vector<std::string> > ("VertexModuleLabels"    	 );

                fScalingfactor      	   = pset.get<float                     > ("Scalingfactor",       1.0);
                fDrawIntersection          = pset.get<int                       > ("drawIntersection"     , 0);
            }

            //----------------------------------------------------
            EventDisplay3D::~EventDisplay3D()
            {
            }

            //----------------------------------------------------
            void EventDisplay3D::makeNavPanel()
            {
                // Create control panel for event navigation
                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                TEveBrowser* browser = fEve->GetBrowser();
                browser->StartEmbedding(TRootBrowser::kLeft); // insert nav frame as new tab in left pane

                TGMainFrame* frmMain = new TGMainFrame(gClient->GetRoot(), 800, 600);
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

            //----------------------------------------------------
            void EventDisplay3D::beginJob()
            {
                // Initialize global Eve application manager (return gEve)
                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                gROOT->SetBatch(kFALSE);
                fEve = TEveManager::Create();

                TEveElementList* simple = new TEveElementList("simplifiedGeometry");
                this->PickVolumes(simple);

                //Add the simplifiedGeometry
                fEve->AddGlobalElement(simple);

                //Create the navigation panel
                makeNavPanel();

                //Sets the GL viewer parameters
                glViewer = fEve->GetDefaultGLViewer();
                double ref[3] = {fGeometry->TPCXCent(), fGeometry->TPCYCent(), fGeometry->TPCZCent()};

                glViewer->ColorSet().Background().SetColor(kWhite);
                glViewer->SetCurrentCamera(TGLViewer::kCameraPerspXOZ);
                glViewer->SetGuideState(TGLUtil::kAxesEdge, kTRUE, kTRUE, ref);
                glViewer->SetDrawCameraCenter(kTRUE);

                std::cout << "Drawing at reference (" << ref[0] << ", " << ref[1] << ", " << ref[2] << ")" << std::endl;

                ann = new TGLAnnotation(glViewer, "", 0.75, 0.85);
                ann->SetRole(TGLOverlayElement::kViewer);
                ann->SetUseColorSet(true);
                ann->SetTextSize(0.04);// % of window diagonal
                ann->SetTextAlign(TGLFont::kLeft);
                ann->SetAllowClose(false);
                ann->SetTextColor(kBlack);

                //Draw the display
                fEve->Redraw3D(kTRUE);

                return;
            }

            //----------------------------------------------------
            void EventDisplay3D::beginRun(const art::Run& run)
            {
                return;
            }

            //----------------------------------------------------
            void EventDisplay3D::analyze(const art::Event& event)
            {
                this->cleanEvt();

                this->UpdateHeader(event);

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

                //Draw MCTruth
                if(fDrawMCTruth)
                this->DrawMCTruth(event);

                //Draw raw hits
                if(fDrawECALRawHits)
                this->DrawRawHits(event);

                //Draw reco hits
                if(fDrawECALRecoHits)
                this->DrawRecoHits(event);

                //Draw high level reco
                if(fDrawECALClusters || fDrawTracks || fDrawVertices)
                this->DrawHighLevelReco(event);

                return;
            } // end EventDisplay3D::EventDisplay3D::analyze

            //----------------------------------------------------
            void EventDisplay3D::cleanEvt()
            {
                //check if any element is present in the event and remove it
                if(fEve->GetCurrentEvent() != nullptr)
                fEve->GetCurrentEvent()->DestroyElements();
            }

            //----------------------------------------------------
            void EventDisplay3D::endJob(){

            }

            //----------------------------------------------------
            void EventDisplay3D::PickVolumes(TEveElementList* &list)
            {
                //Check if the top volume is Enclosure or ND
                bool hasEnclosure = false;
                std::vector<const TGeoNode*> topnodes = fGeometry->FindVolumePath("volDetEnclosure_0");
                TGeoScale nullmatgm;
                TGeoMatrix* topEnclosure = &nullmatgm;
                if(topnodes.size() != 0) hasEnclosure = true;

                std::cout << "The top Volume is Enclosure? " << hasEnclosure << std::endl;

                //Matrix Enclosure
                if(hasEnclosure)
                {
                    for(unsigned int i = 0; i < topnodes.size(); i++)
                    {
                        std::string nodename(topnodes.at(i)->GetName());
                        if(nodename.find("volDetEnclosure_0") == std::string::npos) continue;

                        const TGeoNode *top = topnodes.at(i);
                        topEnclosure = top->GetMatrix();
                    }
                }

                //Matrix ND
                TGeoMatrix* topND = &nullmatgm;
                //Get the matrix of the top volume (rotation of the full ND)
                topnodes = fGeometry->FindVolumePath("volNDHPgTPC_0");
                for(unsigned int i = 0; i < topnodes.size(); i++)
                {
                    std::string nodename(topnodes.at(i)->GetName());
                    if(nodename.find("volNDHPgTPC_0") == std::string::npos) continue;

                    const TGeoNode *top = topnodes.at(i);
                    topND = top->GetMatrix();
                }

                for(auto volname : fVolumesToShow)
                {
                    std::vector<const TGeoNode*> nodevec = fGeometry->FindVolumePath(volname);

                    for(unsigned int i = 0; i < nodevec.size(); i++)
                    {
                        std::string nodename(nodevec.at(i)->GetName());
                        if(nodename.find(volname) == std::string::npos) continue;

                        const TGeoNode *node = nodevec.at(i);

                        //If Barrel or Endcap -> show daughters
                        if(nodename.find("BarrelECal") != std::string::npos || nodename.find("EndcapECal") != std::string::npos)
                        {
                            TGeoScale nullmatgm;
                            TGeoMatrix* topECal = &nullmatgm;
                            topECal = node->GetMatrix();

                            for(int idaugh = 0; idaugh < node->GetNdaughters(); idaugh++)
                            {
                                const TGeoNode *daugh_node = node->GetDaughter(idaugh);
                                TGeoShape *daugh_shape = daugh_node->GetVolume()->GetShape();
                                TGeoShape* daugh_clonedShape = dynamic_cast<TGeoShape*> (daugh_shape->Clone("fakeShape"));
                                TEveGeoShape *daugh_fakeShape = new TEveGeoShape(daugh_shape->GetName());
                                daugh_fakeShape->SetShape(daugh_clonedShape);

                                if(nodename.find("BarrelECal") != std::string::npos)
                                {
                                    daugh_fakeShape->SetMainColor(kOrange+10);
                                    daugh_fakeShape->SetMainTransparency(50);
                                }
                                else if(nodename.find("EndcapECal") != std::string::npos)
                                {
                                    daugh_fakeShape->SetMainColor(kOrange+10);
                                    daugh_fakeShape->SetMainTransparency(70);
                                }

                                TGeoMatrix* currMat = daugh_node->GetMatrix();
                                TGeoMatrix* mat = currMat->MakeClone();
                                TGeoHMatrix *m = new TGeoHMatrix(*mat);
                                m->MultiplyLeft(topECal);
                                m->MultiplyLeft(topND);
                                if(hasEnclosure) m->MultiplyLeft(topEnclosure);
                                daugh_fakeShape->SetTransMatrix(*m);

                                list->AddElement(daugh_fakeShape);
                            }
                        }
                        else{
                            TGeoShape *shape = node->GetVolume()->GetShape();
                            TGeoShape* clonedShape = dynamic_cast<TGeoShape*> (shape->Clone("fakeShape"));

                            TEveGeoShape *fakeShape = new TEveGeoShape(shape->GetName());
                            fakeShape->SetShape(clonedShape);

                            if(nodename.find("TPC") != std::string::npos)
                            {
                                fakeShape->SetMainColor(kBlue+8);
                                fakeShape->SetMainTransparency(50);
                            }
                            else if(nodename.find("PV") != std::string::npos)
                            {
                                fakeShape->SetMainColor(kGray+1);
                                fakeShape->SetMainTransparency(70);
                            }
                            else{
                                fakeShape->SetMainColor(node->GetVolume()->GetLineColor());
                                fakeShape->SetMainTransparency(70);
                            }

                            TGeoMatrix* currMat = node->GetMatrix();
                            TGeoMatrix* mat = currMat->MakeClone();
                            TGeoHMatrix *m = new TGeoHMatrix(*mat);
                            m->MultiplyLeft(topND);
                            if(hasEnclosure) m->MultiplyLeft(topEnclosure);
                            fakeShape->SetTransMatrix(*m);
                            // fakeShape->SetTransMatrix(*mat);

                            list->AddElement(fakeShape);
                        }
                    }
                    nodevec.clear();
                }

                return;
            }

            //----------------------------------------------------
            void EventDisplay3D::DrawMCTruth(const art::Event& event)
            {
                // Define a couple of colors for neutrals and if we gray it out...
                int neutralColor(12);
                int neutrinoColor(38);

                std::vector<const simb::MCParticle*> plist;
                this->GetParticle(plist, event);

                fMCTrajectoryList = new TEveElementList("g4Trajectories", "Geant4 Trajectories");
                fMCTrajectoryList->SetMainColor(kRed);
                fMCTrajectoryList->SetMainAlpha(1.0);

                double minPartEnergy(0.01);

                for(unsigned int p = 0; p < plist.size(); ++p)
                {
                    // Is there an associated McTrajectory?
                    const simb::MCParticle*   mcPart = plist[p];
                    const simb::MCTrajectory& mcTraj = mcPart->Trajectory();

                    int           pdgCode(mcPart->PdgCode());
                    int           colorIdx(evd::Style::ColorFromPDG(mcPart->PdgCode()));
                    TParticlePDG* partPDG(TDatabasePDG::Instance()->GetParticle(pdgCode));
                    double        partCharge = partPDG ? partPDG->Charge() : 0.;
                    double        partEnergy = mcPart->E();

                    if(nullptr == partPDG) continue;

                    if (!mcTraj.empty() && partEnergy > minPartEnergy && mcPart->TrackId() < 100000000)
                    {
                        // collect the points from this particle
                        int numTrajPoints = mcTraj.size();

                        std::ostringstream label;
                        label << "Particle pdg " << pdgCode << "\n";
                        label << "Charge: " << partCharge << "\n";
                        label << "Energy: " << partEnergy << " GeV\n";
                        label << "Momentum (" << mcPart->Px() << ", " << mcPart->Py() << ", " << mcPart->Pz() << " )";

                        TEveLine *MCtrack = new TEveLine(numTrajPoints);
                        MCtrack->SetName("MC trajectory");
                        MCtrack->SetTitle(label.str().c_str());
                        if (partCharge == 0.){
                            MCtrack->SetLineColor(neutralColor);
                            MCtrack->SetLineWidth(2);
                            MCtrack->SetLineStyle(3);
                        }
                        else{
                            MCtrack->SetLineColor(colorIdx);
                            MCtrack->SetLineWidth(2);
                            MCtrack->SetLineStyle(1);
                        }

                        for(int hitIdx = 0; hitIdx < numTrajPoints; hitIdx++)
                        {
                            double xPos = mcTraj.X(hitIdx);
                            double yPos = mcTraj.Y(hitIdx);
                            double zPos = mcTraj.Z(hitIdx);

                            // If the original simulated hit did not occur in the enclosure volume then don't draw it
                            TVector3 point(xPos, yPos, zPos);
                            if (!fGeometry->PointInDetEnclosure(point)) continue;

                            MCtrack->SetPoint(hitIdx, xPos, yPos, zPos);
                        }

                        fMCTrajectoryList->AddElement(MCtrack);
                    }
                }

                // Finally, let's see if we can draw the incoming particle from the MCTruth information
                std::vector<const simb::MCTruth*> mctruth;
                this->GetMCTruth(mctruth, event);

                // Loop through the MCTruth vector
                for (unsigned int idx = 0; idx < mctruth.size(); idx++)
                {
                    // Go through each MCTruth object in the list
                    for (int particleIdx = 0; particleIdx < mctruth[idx]->NParticles(); particleIdx++)
                    {
                        const simb::MCParticle& mcPart = mctruth[idx]->GetParticle(particleIdx);

                        // A negative mother id indicates the "primary" particle
                        if(mcPart.Mother() == -1 && mcPart.StatusCode() == 0)
                        {
                            mf::LogDebug("EVD3D") << mcPart << std::endl;

                            // Get position vector
                            TVector3 particlePosition(mcPart.Vx(), mcPart.Vy(), mcPart.Vz());

                            // Get direction vector (in opposite direction)
                            TVector3 oppPartDir(-mcPart.Px(), -mcPart.Py(), -mcPart.Pz());

                            if (oppPartDir.Mag2() > 0.) oppPartDir.SetMag(1.);

                            double arcLenToDraw = 50.0;  // always draw 50 cm of neutrinos

                            // Draw the line, use an off color to be unique
                            TEveLine *motherPart = new TEveLine(2);
                            motherPart->SetLineColor(neutrinoColor);
                            motherPart->SetLineWidth(2);
                            motherPart->SetLineStyle(2);

                            motherPart->SetPoint(0, particlePosition.X(), particlePosition.Y(), particlePosition.Z());
                            particlePosition += std::min(arcLenToDraw + 10., 1000.) * oppPartDir;
                            motherPart->SetPoint(1,particlePosition.X(), particlePosition.Y(), particlePosition.Z());

                            fMCTrajectoryList->AddElement(motherPart);
                        }
                        // The particles we want to draw will be early in the list so break out if we didn't find them
                        else break;
                    } // loop on particles in list
                }

                fEve->AddElement(fMCTrajectoryList);
            }

            //----------------------------------------------------
            void EventDisplay3D::DrawRawHits(const art::Event& event)
            {
                std::vector<const raw::CaloRawDigit*> rawlist;
                this->GetRawCaloHits(rawlist, event);

                fCaloRawHitList = new TEveElementList("ECAL Raw Calo Hits", "Digitized ECAL hits");
                fCaloRawHitList->SetMainColor(kRed);
                fCaloRawHitList->SetMainAlpha(1.0);

                for(unsigned int p = 0; p < rawlist.size(); ++p)
                {
                    const raw::CaloRawDigit* rawHit = rawlist[p];

                    std::ostringstream label;
                    label << "Digi Hit " << p << "\n";
                    label << "Energy: " << rawHit->ADC() << " ADC\n";
                    label << "Position (" << rawHit->X() << ", " << rawHit->X() << ", " << rawHit->Z() << " ) cm\n";
                    label << "CellID: " << rawHit->CellID();

                    TEvePointSet *evehit = new TEvePointSet(1);
                    evehit->SetName(TString::Format("ECAL digi hit %i", p).Data());
                    evehit->SetTitle(label.str().c_str());
                    evehit->SetMarkerSize(0.5);
                    evehit->SetMarkerStyle(20);
                    evehit->SetMarkerColor(fEvtDisplayUtil->LogColor(rawHit->ADC(), 0, 4096, 5));
                    evehit->SetPoint(0, rawHit->X(), rawHit->Y(), rawHit->Z());//cm
                    fCaloRawHitList->AddElement(evehit);
                }

                fEve->AddElement(fCaloRawHitList);
            }

            //----------------------------------------------------
            void EventDisplay3D::DrawRecoHits(const art::Event& event)
            {
                std::vector<const rec::CaloHit*> recolist;
                this->GetRecoCaloHits(recolist, event);

                fCaloRecoHitList = new TEveElementList("ECAL Reco Calo Hits", "Reconstructed ECAL hits");
                fCaloRecoHitList->SetMainColor(kRed);
                fCaloRecoHitList->SetMainAlpha(1.0);

                for(unsigned int p = 0; p < recolist.size(); ++p)
                {
                    const rec::CaloHit* recoHit = recolist[p];

                    std::ostringstream label;
                    label << "Reco Hit " << p << "\n";
                    label << "Energy: " << recoHit->Energy() * 1000 << " MeV\n";
                    label << "Position (" << recoHit->Position()[0] << ", " << recoHit->Position()[1] << ", " << recoHit->Position()[2] << " ) cm\n";
                    label << "CellID: " << recoHit->CellID();

                    TEvePointSet *evehit = new TEvePointSet(1);
                    evehit->SetName(TString::Format("ECAL reco hit %i", p));
                    evehit->SetTitle(label.str().c_str());
                    evehit->SetMarkerSize(0.5);
                    evehit->SetMarkerStyle(20);
                    evehit->SetMarkerColor(fEvtDisplayUtil->LogColor(recoHit->Energy() * 1000, 0, 1000, 5));
                    evehit->SetPoint(0, recoHit->Position()[0], recoHit->Position()[1], recoHit->Position()[2]);//cm

                    fCaloRecoHitList->AddElement(evehit);
                }

                fEve->AddElement(fCaloRecoHitList);
            }

            //----------------------------------------------------
            void EventDisplay3D::DrawHighLevelReco(const art::Event& event)
            {
                if(fDrawECALClusters)
                {
                    std::vector<const gar::rec::Cluster*> clusters;
                    this->GetCaloClusters(clusters, event);

                    fCaloClusterList = new TEveElementList("ECAL Clusters", "Clustered ECAL hits");
                    fCaloClusterList->SetMainColor(kRed);
                    fCaloClusterList->SetMainAlpha(1.0);

                    for(unsigned int p = 0; p < clusters.size(); ++p)
                    {
                        const rec::Cluster* clus = clusters[p];

                        auto const* pos = clus->Position();
                        auto const* shape = clus->Shape();
                        auto const* dir = clus->EigenVectors();
                        // float length = clus->Energy();   // todo -- make the scale factor a fcl parameter

                        //position of the 6 hit of that define the outline of the cluster in x, y, z
                        float poslist[6][3];

                        for (int i = 0; i < 3; ++i)
                        {
                            poslist[0][i] = pos[i] + (shape[0] * dir[i]) * fScalingfactor;
                            poslist[1][i] = pos[i] - (shape[1] * dir[i]) * fScalingfactor;

                            poslist[2][i] = pos[i] + (shape[2] * dir[3+i]) * fScalingfactor;
                            poslist[3][i] = pos[i] - (shape[2] * dir[3+i]) * fScalingfactor;

                            poslist[4][i] = pos[i] + (shape[3] * dir[6+i]) * fScalingfactor;
                            poslist[5][i] = pos[i] - (shape[3] * dir[6+i]) * fScalingfactor;
                        }

                        std::ostringstream label;
                        label << "Cluster " << p << "\n";
                        label << "Energy: " << clus->Energy() * 1000 << " MeV\n";
                        label << "Position (" << clus->Position()[0] << ", " << clus->Position()[1] << ", " << clus->Position()[2] << " ) cm\n";
                        label << "Shape: r_forw " << shape[0] << ", r_bck " << shape[1] << ", r2 " << shape[2] << ", r3 " << shape[3] << ", vol " << shape[4] << ", width " << shape[5];

                        TEveLine *eve_r1cluster = new TEveLine(6);
                        eve_r1cluster->SetName(TString::Format("ECAL cluster %i", p));
                        eve_r1cluster->SetTitle(label.str().c_str());
                        eve_r1cluster->SetLineWidth(2);
                        eve_r1cluster->SetLineStyle(1);
                        eve_r1cluster->SetLineColor(p%evd::kNCOLS);
                        eve_r1cluster->SetPoint(0, poslist[0][0], poslist[0][1], poslist[0][2]);
                        eve_r1cluster->SetPoint(1, poslist[1][0], poslist[1][1], poslist[1][2]);
                        eve_r1cluster->SetPoint(2, poslist[2][0], poslist[2][1], poslist[2][2]);
                        eve_r1cluster->SetPoint(3, poslist[3][0], poslist[3][1], poslist[3][2]);
                        eve_r1cluster->SetPoint(4, poslist[4][0], poslist[4][1], poslist[4][2]);
                        eve_r1cluster->SetPoint(5, poslist[5][0], poslist[5][1], poslist[5][2]);

                        fCaloClusterList->AddElement(eve_r1cluster);
                    }

                    fEve->AddElement(fCaloClusterList);
                }

                if(fDrawTracks)
                {
                    std::vector<const gar::rec::Track*> tracks;
                    this->GetTracks(tracks, event);

                    fTrackList = new TEveElementList("Fitted Tracks", "Tracks in the TPC");
                    fTrackList->SetMainColor(kRed);
                    fTrackList->SetMainAlpha(1.0);

                    unsigned int icounter = 0;
                    for (auto tv = tracks.begin(); tv != tracks.end(); ++tv)
                    {
                        //Forward
                        this->DrawTrack((*tv), fTrackList, icounter, true);
                        //Backward
                        this->DrawTrack((*tv), fTrackList, icounter, false);

                        if(fDrawIntersection)
                        {
                            this->DrawIntersections((*tv), fTrackList, icounter);
                        }

                        ++icounter;
                    }

                    fEve->AddElement(fTrackList);
                }

                if(fDrawVertices)
                {
                    std::vector<const gar::rec::Vertex*> vertexs;
                    this->GetVertices(vertexs, event);

                    fVertexList = new TEveElementList("Reconstructed vertex", "Reconstructed vertex");
                    fVertexList->SetMainColor(kRed);
                    fVertexList->SetMainAlpha(1.0);

                    for(unsigned int p = 0; p < vertexs.size(); ++p)
                    {
                        const rec::Vertex* vx = vertexs[p];

                        auto const* pos = vx->Position();

                        TEvePointSet *evevx = new TEvePointSet(1);
                        evevx->SetMarkerSize(1);
                        evevx->SetMarkerStyle(20);
                        evevx->SetMarkerColor(5);
                        evevx->SetPoint(0, pos[0], pos[1], pos[2]);

                        fVertexList->AddElement(evevx);
                    }

                    fEve->AddElement(fVertexList);
                }
            }

            //----------------------------------------------------
            void EventDisplay3D::DrawTrack(const gar::rec::Track* trk, TEveElementList* &eve_list, int counter, bool forward)
            {
                int color  = evd::kColor[counter%evd::kNCOLS];

                const float *fTrackpar = (forward == true) ? trk->TrackParBeg() : trk->TrackParEnd();
                const float fStartX = (forward == true) ? trk->Vertex()[0] : trk->End()[0];
                const float fEndX = (forward == true) ? trk->End()[0] : trk->Vertex()[0];
                const float fMomentum = (forward == true) ? trk->Momentum_beg() : trk->Momentum_end();
                TString trackname = (forward == true) ? TString::Format("Track Fit forward %i", counter) : TString::Format("Track Fit backward %i", counter);

                TEveLine *eve_track = new TEveLine();
                this->DrawHelix3D(fTrackpar, fStartX, fEndX, eve_track, color);

                eve_track->SetName(trackname.Data());
                std::ostringstream label;
                label << "Track " << counter << "\n";
                label << "Vertex (" << trk->Vertex()[0] << ", " << trk->Vertex()[1] << ", " << trk->Vertex()[2] << ") cm\n";
                label << "End (" << trk->End()[0] << ", " << trk->End()[1] << ", " << trk->End()[2] << ") cm\n";
                label << "Omega " << fTrackpar[2] << "\n";
                label << "phi " << fTrackpar[3] << "\n";
                label << "TanLambda " << std::tan(fTrackpar[4]) << "\n";
                label << "Momentum " << fMomentum << " GeV\n";
                label << "Charge " << fTrackpar[2] / std::fabs(fTrackpar[2]);
                eve_track->SetTitle(label.str().c_str());

                const float *fVertex = (forward == true) ? trk->Vertex() : trk->End();
                const float *fDir = (forward == true) ? trk->VtxDir() : trk->EndDir();

                TEveArrow *arrow = new TEveArrow();
                this->DrawArrow3D(fVertex, fDir, color, arrow);

                eve_list->AddElement(eve_track);
                eve_list->AddElement(arrow);

                return;
            }

            //----------------------------------------------------
            void EventDisplay3D::DrawIntersections(const gar::rec::Track* trk, TEveElementList* &fTrackList, int counter)
            {
                int color  = evd::kColor[counter%evd::kNCOLS];
                //Check the intersection with the ECAL Barrel
                float xyz_intersection[3] = {0., 0., 0.};

                //Forward case
                int result = util::TrackPropagator::PropagateToCylinder(trk->TrackParEnd(), fGeometry->GetECALEndcapStartX(), fGeometry->GetECALInnerBarrelRadius(), fGeometry->TPCYCent(), fGeometry->TPCZCent(), trk->End()[0], xyz_intersection);
                std::ostringstream label;

                if(result == 0)
                {
                    TEvePointSet *intersection1 = new TEvePointSet(1);
                    intersection1->SetName(TString::Format("Track forward Calo Intersection Barrel %i", counter).Data());
                    label << "Barrel Intersection forward - Track " << counter << "\n";
                    label << "Position (" << xyz_intersection[0] << ", " << xyz_intersection[1] << ", " << xyz_intersection[2] << ") cm";
                    intersection1->SetTitle(label.str().c_str());
                    intersection1->SetMarkerSize(1.3);
                    intersection1->SetMarkerStyle(22);
                    intersection1->SetMarkerColor(color);
                    intersection1->SetPoint(0, xyz_intersection[0], xyz_intersection[1], xyz_intersection[2]);//cm
                    fTrackList->AddElement(intersection1);
                }
                else{
                    //Try Endcap
                    result = util::TrackPropagator::PropagateToXPlane(trk->TrackParEnd(), (trk->End()[0] > 0) ? fGeometry->GetECALEndcapStartX() : -fGeometry->GetECALEndcapStartX(), fGeometry->GetECALOuterBarrelRadius(), trk->End()[0], xyz_intersection);

                    if(result == 0)
                    {
                        TEvePointSet *intersection3 = new TEvePointSet(1);
                        intersection3->SetName(TString::Format("Track forward Calo Intersection Endcap %i", counter).Data());
                        label.str("");
                        label.clear();
                        label << "Endcap Intersection forward - Track " << counter << "\n";
                        label << "Position (" << xyz_intersection[0] << ", " << xyz_intersection[1] << ", " << xyz_intersection[2] << ") cm";
                        intersection3->SetTitle(label.str().c_str());
                        intersection3->SetMarkerSize(1.3);
                        intersection3->SetMarkerStyle(29);
                        intersection3->SetMarkerColor(color);
                        intersection3->SetPoint(0, xyz_intersection[0], xyz_intersection[1], xyz_intersection[2]);//cm
                        fTrackList->AddElement(intersection3);
                    }
                }

                //Backward case
                result = util::TrackPropagator::PropagateToCylinder(trk->TrackParBeg(), fGeometry->GetECALEndcapStartX(), fGeometry->GetECALInnerBarrelRadius(), fGeometry->TPCYCent(), fGeometry->TPCZCent(), trk->Vertex()[0], xyz_intersection);

                if(result == 0)
                {
                    TEvePointSet *intersection2 = new TEvePointSet(1);
                    intersection2->SetName(TString::Format("Track backward Calo Intersection Barrel %i", counter).Data());
                    label.str("");
                    label.clear();
                    label << "Barrel Intersection backward - Track " << counter << "\n";
                    label << "Position (" << xyz_intersection[0] << ", " << xyz_intersection[1] << ", " << xyz_intersection[2] << ") cm";
                    intersection2->SetTitle(label.str().c_str());
                    intersection2->SetMarkerSize(1.3);
                    intersection2->SetMarkerStyle(23);
                    intersection2->SetMarkerColor(color);
                    intersection2->SetPoint(0, xyz_intersection[0], xyz_intersection[1], xyz_intersection[2]);//cm
                    fTrackList->AddElement(intersection2);
                }
                else{
                    //Try Endcap
                    result = util::TrackPropagator::PropagateToXPlane(trk->TrackParBeg(), (trk->Vertex()[0] > 0) ? fGeometry->GetECALEndcapStartX() : -fGeometry->GetECALEndcapStartX(), fGeometry->GetECALOuterBarrelRadius(), trk->Vertex()[0], xyz_intersection);

                    if(result == 0)
                    {
                        TEvePointSet *intersection4 = new TEvePointSet(1);
                        intersection4->SetName(TString::Format("Track backward Calo Intersection Endcap %i", counter).Data());
                        label.str("");
                        label.clear();
                        label << "Endcap Intersection backward - Track " << counter << "\n";
                        label << "Position (" << xyz_intersection[0] << ", " << xyz_intersection[1] << ", " << xyz_intersection[2] << ") cm";
                        intersection4->SetTitle(label.str().c_str());
                        intersection4->SetMarkerSize(1.3);
                        intersection4->SetMarkerStyle(29);
                        intersection4->SetMarkerColor(color);
                        intersection4->SetPoint(0, xyz_intersection[0], xyz_intersection[1], xyz_intersection[2]);//cm
                        fTrackList->AddElement(intersection4);
                    }
                }

                return;
            }

            //----------------------------------------------------
            void EventDisplay3D::DrawHelix3D(const float *trackpar, const float xpar, const float xother, TEveLine* &eve_track, int color)
            {
                float r = 0;
                if (trackpar[2] != 0) r=1.0/trackpar[2];

                float si = TMath::Tan(trackpar[4]);
                float ycc = trackpar[0] + r*TMath::Cos(trackpar[3]);
                float zcc = trackpar[1] - r*TMath::Sin(trackpar[3]);
                float dphimax = TMath::Pi();
                float phi1 = 0;
                float phi2 = (xother-xpar)*trackpar[2]/si;

                if (phi2-phi1>dphimax) phi2 = phi1+dphimax;
                if (phi2-phi1<-dphimax) phi2 = phi1-dphimax;
                if (phi2-phi1==0) phi2=phi1+0.01;

                int nptshelix=100;
                float dphihelix=(phi2-phi1)/((float) nptshelix);

                for (int ipoint=0;ipoint<nptshelix;++ipoint)
                {
                    float philoc = phi1 + ipoint*dphihelix;
                    float xl = xpar + r*si*(philoc);
                    float yl = ycc - r*TMath::Cos(philoc + trackpar[3]);
                    float zl = zcc + r*TMath::Sin(philoc + trackpar[3]);

                    TVector3 point(xl, yl, zl);
                    if(!fGeometry->PointInGArTPC(point)) continue;

                    eve_track->SetPoint(ipoint, xl, yl, zl);
                }

                eve_track->SetLineWidth(2);
                eve_track->SetLineStyle(1);
                eve_track->SetLineColor(color);
            }

            //----------------------------------------------------
            void EventDisplay3D::DrawArrow3D(const float *fVertex, const float *fDir, int color, TEveArrow* &arrow)
            {
                //Draw arraw
                TVector3 cent(fGeometry->TPCXCent(), fGeometry->TPCYCent(), fGeometry->TPCZCent());
                TVector3 spv(fVertex);
                TVector3 spcv = spv + cent;
                TVector3 av(fDir);

                arrow->SetOrigin(spcv.X(), spcv.Y(), spcv.Z());
                arrow->SetVector(av.X()*10, av.Y()*10, av.Z()*10);
                arrow->SetMainColor(color);
                arrow->SetTubeR(0.1);
                arrow->SetConeL(0.5);
                arrow->SetConeR(0.2);

                return;
            }

            //----------------------------------------------------
            void EventDisplay3D::GetRawCaloHits(std::vector<const raw::CaloRawDigit*> &digitCalo, const art::Event& event)
            {
                digitCalo.clear();
                std::vector<const raw::CaloRawDigit*> tempCalo;

                try
                {
                    event.getView(fRawHitLabels.at(0), tempCalo);
                    for(size_t t = 0; t < tempCalo.size(); ++t)
                    digitCalo.push_back(tempCalo[t]);
                }
                catch(cet::exception& e){
                    writeErrMsg("GetDigits Calo", e);
                }
            }

            //----------------------------------------------------
            void EventDisplay3D::GetRecoCaloHits(std::vector<const rec::CaloHit*> &recoCalo, const art::Event& event)
            {
                recoCalo.clear();
                std::vector<const rec::CaloHit*> tempCalo;

                try
                {
                    event.getView(fRecoHitLabels.at(0), tempCalo);
                    for(size_t t = 0; t < tempCalo.size(); ++t)
                    recoCalo.push_back(tempCalo[t]);
                }
                catch(cet::exception& e){
                    writeErrMsg("GetHits Calo", e);
                }
            }

            //----------------------------------------------------
            void EventDisplay3D::GetCaloClusters(std::vector<const rec::Cluster*> &caloCluster, const art::Event& event)
            {
                caloCluster.clear();

                std::vector<const gar::rec::Cluster*> temp;

                try{
                    event.getView(fCaloClusterLabels.at(0), temp);
                    for(size_t t = 0; t < temp.size(); ++t){
                        caloCluster.push_back(temp[t]);
                    }
                }
                catch(cet::exception& e){
                    writeErrMsg("GetCaloClusters", e);
                }

                return;
            }

            //----------------------------------------------------
            void EventDisplay3D::GetTracks(std::vector<const rec::Track*> &track, const art::Event& event)
            {
                track.clear();

                std::vector<const gar::rec::Track*> temp;

                try{
                    event.getView(fTrackLabels.at(0), temp);
                    for(size_t t = 0; t < temp.size(); ++t){
                        track.push_back(temp[t]);
                    }
                }
                catch(cet::exception& e){
                    writeErrMsg("GetTracks", e);
                }

                return;
            }

            //----------------------------------------------------
            void EventDisplay3D::GetVertices(std::vector<const rec::Vertex*> &vertex, const art::Event& event)
            {
                vertex.clear();

                std::vector<const gar::rec::Vertex*> temp;

                try{
                    event.getView(fVertexLabels.at(0), temp);
                    for(size_t t = 0; t < temp.size(); ++t){
                        vertex.push_back(temp[t]);
                    }
                }
                catch(cet::exception& e){
                    writeErrMsg("GetVertices", e);
                }

                return;
            }

            //----------------------------------------------------
            void EventDisplay3D::GetMCTruth(std::vector<const simb::MCTruth*>& mcvec, const art::Event& event)
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

            //----------------------------------------------------
            void EventDisplay3D::GetParticle(std::vector<const simb::MCParticle*>& plist, const art::Event& event)
            {
                plist.clear();

                if( event.isRealData() ) return;

                std::vector<const simb::MCParticle*> temp;
                art::View<simb::MCParticle> plcol;
                // use get by Type because there should only be one collection of these in the event
                try
                {
                    event.getView(fG4Label, plcol);
                    for(unsigned int i = 0; i < plcol.vals().size(); ++i)
                    temp.push_back(plcol.vals().at(i));

                    temp.swap(plist);
                }
                catch(cet::exception& e){
                    writeErrMsg("GetParticle", e);
                }

                return;
            }

            //----------------------------------------------------
            void EventDisplay3D::UpdateHeader(const art::Event& event)
            {
                unsigned long long int tsval = event.time().value();
                int run   = event.run();
                int srun  = event.subRun();
                int evt = event.id().event();

                unsigned int year, month, day, dayofweek;
                unsigned int hour, minute, second;
                int          nano;

                const unsigned long int mask32 = 0xFFFFFFFFUL;
                unsigned long int lup = ( tsval >> 32 ) & mask32;
                unsigned long int llo = tsval & mask32;
                TTimeStamp ts(lup, (int)llo);

                ts.GetDate(kTRUE,0,&year,&month,&day);
                ts.GetTime(kTRUE,0,&hour,&minute,&second);
                nano = ts.GetNanoSec();
                dayofweek = ts.GetDayOfWeek();
                char eventbuff[256];
                char runbuff[256];
                char datebuff[256];
                char timebuff[256];

                // Skip first one since ROOT returns these numbers starting from 1 not 0
                static const char* days[] = {"","Mon","Tue","Wed","Thu","Fri","Sat","Sun"};
                static const char* months[] = {"","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"};

                sprintf(runbuff,  "Run:   %d/%d\n", run, srun);
                sprintf(eventbuff,"Event: %d\n", evt);
                sprintf(datebuff, "UTC %s %s %d, %d\n", days[dayofweek], months[month], day, year);
                sprintf(timebuff, "%.2d:%.2d:%2.9f\n", hour, minute, (float)second+(float)nano/1.0E9);

                //Final char
                char annotation[256];
                sprintf(annotation, "DUNE ND HPgTPC\n%s%s%s%s", runbuff, eventbuff, datebuff, timebuff);
                ann->SetText(annotation);
                ann->SetRole(TGLOverlayElement::kViewer);
                ann->SetUseColorSet(true);
                ann->SetTextSize(0.04);// % of window diagonal
                ann->SetTextAlign(TGLFont::kLeft);
                ann->SetAllowClose(false);
                ann->SetTextColor(kBlack);
            }

            DEFINE_ART_MODULE(EventDisplay3D)

        }//end namespace evd
    }//end namespace gar
