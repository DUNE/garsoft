////////////////////////////////////////////////////////////////////////
// Class:       no class at all; trailer trash.
/// \file       anautil.cxx
/// \brief Code to compute quantities that will go into the anatree
///
/// \version $Id:  $
/// \author  bellanto@fnal.gov
//
// This is a separate file of functions used by anatree_module.cc
// It exists to keep that file of manageable size.
//
// The functions are:
//     void  gar::anatree::ClearVectors()
//     void  gar::anatree::FillVectors(art::Event const & e)
//     float gar::anatree::processIonizationInfo(rec::TrackIoniz& ion, float ionizeTruncate,
//                                           float& forwardIonVal, float& backwardIonVal)
//     float gar::anatree::computeT(simb::MCTruth theMCTruth)
//
////////////////////////////////////////////////////////////////////////



//#include "Ana/anautil.h"
#include "anatree_module.cc"



namespace gar {
//    namespace anatree {



        //==========================================================================
        void anatree::ClearVectors() {
            // clear out all our vectors
            if (fWriteMCinfo) {
                fNeutrinoType.clear();
                fCCNC.clear();
                fMode.clear();
                fInteractionType.clear();
                fQ2.clear();
                fW.clear();
                fX.clear();
                fY.clear();
                fTheta.clear();
                if (fWriteCohInfo) fT.clear();
                fMCVertexX.clear();
                fMCVertexY.clear();
                fMCVertexZ.clear();
                fMCnuPx.clear();
                fMCnuPy.clear();
                fMCnuPz.clear();
        
                fGint.clear();
                fTgtPDG.clear();
                fWeight.clear();
                fgT.clear();
        
                fMCPTrkID.clear();
                fMCPDG.clear();
                fMCMother.clear();
                fMCPDGMother.clear();
                fMCPStartX.clear();
                fMCPStartY.clear();
                fMCPStartZ.clear();
                fMCPTime.clear();
                fMCPStartPX.clear();
                fMCPStartPY.clear();
                fMCPStartPZ.clear();
                fMCPEndX.clear();
                fMCPEndY.clear();
                fMCPEndZ.clear();
                fMCPEndPX.clear();
                fMCPEndPY.clear();
                fMCPEndPZ.clear();
                fMCPProc.clear();
                fMCPEndProc.clear();
        
                if (fWriteMCPTrajectory) {
                    fTrajMCPX.clear();
                    fTrajMCPY.clear();
                    fTrajMCPZ.clear();
                    fTrajMCPT.clear();
                    fTrajMCPE.clear();
                    fTrajMCPTrajIndex.clear();
                }
        
                fSimnHits = 0;
                fSimHitX.clear();
                fSimHitY.clear();
                fSimHitZ.clear();
                fSimHitTime.clear();
                fSimHitEnergy.clear();
                fSimHitTrackID.clear();
                fSimHitCellID.clear();
                fSimEnergySum = 0.;
            }
        
            if (fWriteHits) {
                fHitX.clear();
                fHitY.clear();
                fHitZ.clear();
                fHitSig.clear();
                fHitRMS.clear();
            }
        
            if (fWriteTPCClusters) {
                fTPCClusterX.clear();
                fTPCClusterY.clear();
                fTPCClusterZ.clear();
                fTPCClusterSig.clear();
                fTPCClusterRMS.clear();
            }
        
            if (fWriteAllTracks) {
                fTrackStartX.clear();
                fTrackStartY.clear();
                fTrackStartZ.clear();
                fTrackStartPX.clear();
                fTrackStartPY.clear();
                fTrackStartPZ.clear();
                fTrackEndX.clear();
                fTrackEndY.clear();
                fTrackEndZ.clear();
                fTrackEndPX.clear();
                fTrackEndPY.clear();
                fTrackEndPZ.clear();
                fTrackLenF.clear();
                fTrackLenB.clear();
                fTrackAvgIonF.clear();
                fTrackAvgIonB.clear();
                fNTPCClustersOnTrack.clear();
            }
        
            if (fWriteTPCClustersInTracks) {
                fTrkTPCClusterX.clear();
                fTrkTPCClusterY.clear();
                fTrkTPCClusterZ.clear();
                fTrkTPCClusterSig.clear();
                fTrkTPCClusterRMS.clear();
                fTrkTPCClusterTrkIndex.clear();
            }
        
            if (fWriteVertsNtracks) {
                fVertexX.clear();
                fVertexY.clear();
                fVertexZ.clear();
                fVertexT.clear();
                fVertexN.clear();
                fVertexQ.clear();
        
                fVTAssn_Vertex.clear();
                fVTAssn_TrkEndPx.clear();
                fVTAssn_TrkEndPy.clear();
                fVTAssn_TrkEndPz.clear();
                fVTAssn_TrkEndLen.clear();
                fVTAssn_TrkEndQ.clear();
                fVTAssn_TrkEndChiSq.clear();
                fVTAssn_TrkEndIon.clear();
                fVTAssn_TrkEndNclus.clear();
            }
        
            if (fWriteCaloDigits) {
                fDiginHits = 0;
                fDigiHitX.clear();
                fDigiHitY.clear();
                fDigiHitZ.clear();
                fDigiHitTime.clear();
                fDigiHitADC.clear();
                fDigiHitCellID.clear();
            }
        
            if (fWriteCaloInfo) {
                fReconHits = 0;
                fRecoHitX.clear();
                fRecoHitY.clear();
                fRecoHitZ.clear();
                fRecoHitTime.clear();
                fRecoHitEnergy.clear();
                fRecoHitCellID.clear();
                fRecoEnergySum = 0.;
        
                fnCluster = 0;
                fClusterNhits.clear();
                fClusterEnergy.clear();
                fClusterX.clear();
                fClusterY.clear();
                fClusterZ.clear();
                fClusterTheta.clear();
                fClusterPhi.clear();
                fClusterPID.clear();
                // fClusterShape.clear();
                fClusterMainAxisX.clear();
                fClusterMainAxisY.clear();
                fClusterMainAxisZ.clear();
        
                if (fWriteAllTracks) {
                    fECALAssn_Cluster.clear();
                    fECALAssn_TrkWhich.clear();
                    fECALAssn_TrkStartX.clear();
                    fECALAssn_TrkStartY.clear();
                    fECALAssn_TrkStartZ.clear();
                    fECALAssn_TrkStartPX.clear();
                    fECALAssn_TrkStartPY.clear();
                    fECALAssn_TrkStartPZ.clear();
                    fECALAssn_TrkEndX.clear();
                    fECALAssn_TrkEndY.clear();
                    fECALAssn_TrkEndZ.clear();
                    fECALAssn_TrkEndPX.clear();
                    fECALAssn_TrkEndPY.clear();
                    fECALAssn_TrkEndPZ.clear();
                    fECALAssn_TrkLenF.clear();
                    fECALAssn_TrkLenB.clear();
                    fECALAssn_TrkAvgIonF.clear();
                    fECALAssn_TrkAvgIonB.clear();
                    fECALAssn_TrkNTPCClustersOnTrack.clear();
                }
        
            }
        } // end gar::anatree::ClearVectors



        //==========================================================================
        void FillVectors(art::Event const & e) {
        
            fRun    = e.run();
            fSubRun = e.subRun();
            fEvent  = e.id().event();
        
            // Get a good grip on that data!
        
            std::vector< art::Handle< std::vector<simb::MCTruth> > > mcthandlelist;
            std::vector< art::Handle< std::vector<simb::GTruth> > > gthandlelist;
        
            art::Handle< std::vector<simb::MCParticle> > MCPHandle;
            art::Handle< std::vector<gar::sdp::CaloDeposit> > SimHitHandle;
        
            if (fWriteMCinfo) {
                // Handles for MC info might be there, might not
                if (fGeneratorLabels.size()<1)
                {
                    // get them all (even if there are none)
                    e.getManyByType(mcthandlelist);
                }
                else
                {
                    mcthandlelist.resize(fGeneratorLabels.size());
                    for (size_t i=0; i< fGeneratorLabels.size(); ++i)
                    {
                        // complain if we wanted a specific one but didn't find it
                        if (!e.getByLabel(fGeneratorLabels.at(i),mcthandlelist.at(i)))
                        throw cet::exception("anatree") << " No simb::MCTruth branch - " << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
                    }
                }
        
                if (fGENIEGeneratorLabels.size()<1)
                {
                    e.getManyByType(gthandlelist);  // get them all (even if there are none)
                }
                else
                {
                    gthandlelist.resize(fGENIEGeneratorLabels.size());
                    for (size_t i=0; i< fGENIEGeneratorLabels.size(); ++i)
                    {
                        // complain if we wanted a specific one but didn't find it
                        if (!e.getByLabel(fGENIEGeneratorLabels.at(i),gthandlelist.at(i)))
                        throw cet::exception("anatree")
                        << " No simb::GTruth branch - "
                        << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        
                    }
                }
        
                if (!e.getByLabel(fGeantLabel, MCPHandle)) {
                    throw cet::exception("anatree")
                    << " No simb::MCParticle branch - "
                    << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
                }
        
                if (!e.getByLabel(fGeantLabel, SimHitHandle)) {
                    throw cet::exception("anatree")
                    << " No gar::sdp::CaloDeposit branch - "
                    << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
                }
            }
        
            art::Handle< std::vector<gar::rec::Hit> > HitHandle;
            art::Handle< std::vector<gar::rec::TPCCluster> > TPCClusterHandle;
            art::Handle< std::vector<gar::rec::Track> > TrackHandle;
            art::Handle< std::vector<gar::rec::TrackIoniz> > TrackIonHandle;
            art::Handle< std::vector<gar::rec::Vertex> > VertexHandle;
            if (fWriteHits) {
                if (!e.getByLabel(fHitLabel, HitHandle)) {
                    throw cet::exception("anatree")
                    << " No gar::rec::Hit branch - "
                    << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
                }
            }
            
            if (fWriteTPCClusters) {
                if (!e.getByLabel(fTPCClusterLabel, TPCClusterHandle)) {
                    throw cet::exception("anatree")
                    << " No gar::rec::TPCCluster branch - "
                    << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
                }
            }
            
            if(fWriteAllTracks)
            {
                if (!e.getByLabel(fTrackLabel, TrackHandle)) {
                    throw cet::exception("anatree")
                    << " No gar::rec::Track branch - "
                    << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
                }
            
                if (!e.getByLabel(fTrackLabel, TrackIonHandle)) {
                    throw cet::exception("anatree")
                    << " No gar::rec::TrackIoniz branch - "
                    << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
                }
            }
            
            if(fWriteVertsNtracks)
            {
                if (!e.getByLabel(fVertexLabel, VertexHandle)) {
                    throw cet::exception("anatree")
                    << " No gar::rec::Vertex branch - "
                    << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
                }
            }
            
            art::Handle< std::vector<gar::raw::CaloRawDigit> > RawHitHandle;
            art::Handle< std::vector<gar::rec::CaloHit> >      RecoHitHandle;
            art::Handle< std::vector<gar::rec::Cluster> >      RecoClusterHandle;
            if (fWriteCaloInfo) {
                if (!e.getByLabel(fRawCaloHitLabel, RawHitHandle)) {
                    throw cet::exception("anatree")
                    << " No gar::raw::CaloRawDigit branch - "
                    << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
                }
            
                if (!e.getByLabel(fCaloHitLabel, RecoHitHandle)) {
                    throw cet::exception("anatree")
                    << " No gar::rec::CaloHit branch - "
                    << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
                }
            
                if (!e.getByLabel(fClusterLabel, RecoClusterHandle)) {
                    throw cet::exception("anatree")
                    << " No gar::rec::Cluster branch - "
                    << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
                }
            }
            
            
            
            // Pull on the handles now!
            if (fWriteMCinfo) {
                // save MCTruth info
            
                for (size_t imchl = 0; imchl < mcthandlelist.size(); ++imchl)
                {
                    for ( auto const& mct : (*mcthandlelist.at(imchl)) ) {
                        if (mct.NeutrinoSet())
                        {
                            simb::MCNeutrino nuw = mct.GetNeutrino();
                            fNeutrinoType.push_back(nuw.Nu().PdgCode());
                            fCCNC.push_back(nuw.CCNC());
                            fMode.push_back(nuw.Mode());
                            fInteractionType.push_back(nuw.InteractionType());
                            fQ2.push_back(nuw.QSqr());
                            fW.push_back(nuw.W());
                            fX.push_back(nuw.X());
                            fY.push_back(nuw.Y());
                            fTheta.push_back(nuw.Theta());
                            if (fWriteCohInfo) {
                                double getT = gar::anatree::computeT(mct);
                                fT.push_back( static_cast<Float_t>(getT) );
                            }
                            fMCVertexX.push_back(nuw.Nu().EndX());
                            fMCVertexY.push_back(nuw.Nu().EndY());
                            fMCVertexZ.push_back(nuw.Nu().EndZ());
                            fMCnuPx.push_back(nuw.Nu().Px());
                            fMCnuPy.push_back(nuw.Nu().Py());
                            fMCnuPz.push_back(nuw.Nu().Pz());
                        }  // end MC info from MCTruth
                    }
                }
        
                // save GTruth info
                for (size_t igthl = 0; igthl < gthandlelist.size(); ++igthl)
                {
                    for ( auto const& gt : (*gthandlelist.at(igthl)) ) {
                        fGint.push_back(gt.fGint);
                        fTgtPDG.push_back(gt.ftgtPDG);
                        fWeight.push_back(gt.fweight);
                        fgT.push_back(gt.fgT);
                    }
                }
        
                // save MCParticle info (post-GEANT; for pre-GEANT, maybe try MCTruth?)
                Int_t mcpIndex = 0;
                for ( auto const& mcp : (*MCPHandle) ) {
        
                    fMCPTrkID.push_back(mcp.TrackId());
                    fMCPDG.push_back(mcp.PdgCode());
                    // if mcp.Mother() == 0, particle is from initial vertex; MCParticles in
                    // the event per se are indexed from 1 via TrackID.
                    // (*MCPHandle) is a vector<simb::MCParticle> and so is indexed from 0.
                    int momNumber = mcp.Mother();    int momPDG = 0;
                    // Short loop to find the mother.  MCParticle.fmother appears to be the TrackID of
                    // of the mother particle, minus 1.  However, not all TrackID values appear
                    // sequentially in the event record; if they did you could look at the MCParticle
                    // (*MCPHandle)[momNumber-1].  Too bad!  BTW, StatusCode==1 at this point, no point
                    // checking that.
                    if (momNumber>0) {
                        for ( auto const& mcp_short : (*MCPHandle) ) {
                            if (mcp_short.TrackId() == momNumber) {
                                momPDG = mcp_short.PdgCode();
                                break;
                            }
                        }
                    }
                    fMCMother.push_back(momNumber);
                    fMCPDGMother.push_back(momPDG);
        
                    const TLorentzVector& position = mcp.Position(0);
                    const TLorentzVector& momentum = mcp.Momentum(0);
                    fMCPStartX.push_back(position.X());
                    fMCPStartY.push_back(position.Y());
                    fMCPStartZ.push_back(position.Z());
                    fMCPTime.push_back(mcp.T());
                    fMCPStartPX.push_back(momentum.Px());
                    fMCPStartPY.push_back(momentum.Py());
                    fMCPStartPZ.push_back(momentum.Pz());
        
                    const TLorentzVector& positionEnd = mcp.EndPosition();
                    const TLorentzVector& momentumEnd = mcp.EndMomentum();
                    fMCPEndX.push_back(positionEnd.X());
                    fMCPEndY.push_back(positionEnd.Y());
                    fMCPEndZ.push_back(positionEnd.Z());
                    fMCPEndPX.push_back(momentumEnd.Px());
                    fMCPEndPY.push_back(momentumEnd.Py());
                    fMCPEndPZ.push_back(momentumEnd.Pz());
                    fMCPProc.push_back(mcp.Process());
                    fMCPEndProc.push_back(mcp.EndProcess());
            
                    if(fWriteMCPTrajectory) {
                        const TDatabasePDG* databasePDG = TDatabasePDG::Instance();
                        const TParticlePDG* definition = databasePDG->GetParticle( mcp.PdgCode() );
                        //No charge don't store the trajectory
                        if(nullptr == definition || definition->Charge() == 0) continue;
            
                        for(uint iTraj=0; iTraj < mcp.Trajectory().size(); iTraj++) {
                            fTrajMCPX.push_back(mcp.Trajectory().X(iTraj));
                            fTrajMCPY.push_back(mcp.Trajectory().Y(iTraj));
                            fTrajMCPZ.push_back(mcp.Trajectory().Z(iTraj));
                            fTrajMCPT.push_back(mcp.Trajectory().T(iTraj));
                            fTrajMCPE.push_back(mcp.Trajectory().E(iTraj));
                            fTrajMCPTrajIndex.push_back(mcpIndex);
                        }
                        mcpIndex++;
                    }
                }  // end MC info from MCParticle
            
                // Save simulation hit info
                for ( auto const& SimHit : (*SimHitHandle) ) {
                    fSimnHits++;
                    fSimHitX.push_back(SimHit.X());
                    fSimHitY.push_back(SimHit.Y());
                    fSimHitZ.push_back(SimHit.Z());
                    fSimHitTime.push_back(SimHit.Time());
                    fSimHitEnergy.push_back(SimHit.Energy());
                    fSimHitTrackID.push_back(SimHit.TrackID());
                    fSimHitCellID.push_back(SimHit.CellID());
                    fSimEnergySum += SimHit.Energy();
                }
            } // end if MCWriteInfo
        
            // save TPC Hit info
            if (fWriteHits) {
                for ( auto const& Hit : (*HitHandle) ) {
                    fHitX.push_back(Hit.Position()[0]);
                    fHitY.push_back(Hit.Position()[1]);
                    fHitZ.push_back(Hit.Position()[2]);
                    fHitSig.push_back(Hit.Signal());
                    fHitRMS.push_back(Hit.RMS());
                }
            }
        
            // save TPCCluster info
            if (fWriteTPCClusters) {
                for ( auto const& TPCCluster : (*TPCClusterHandle) ) {
                    fTPCClusterX.push_back(TPCCluster.Position()[0]);
                    fTPCClusterY.push_back(TPCCluster.Position()[1]);
                    fTPCClusterZ.push_back(TPCCluster.Position()[2]);
                    fTPCClusterSig.push_back(TPCCluster.Signal());
                    fTPCClusterRMS.push_back(TPCCluster.RMS());
                }
            }
        
            // save Track info.  Need track-TPCCluster-Ionization Assns outside if () for later use
            if (fWriteAllTracks) {
                const art::FindManyP<gar::rec::TPCCluster> findManyTPCClusters(TrackHandle,e,fTrackLabel);
                const art::FindOneP<gar::rec::TrackIoniz>  findIonization(TrackHandle,e,fTrackLabel);
        
                size_t iTrack = 0;
                for ( auto const& track : (*TrackHandle) ) {
        
                    // track is a gar::rec::Track, not a gar::rec::TrackPar
                    fTrackStartX.push_back(track.Vertex()[0]);
                    fTrackStartY.push_back(track.Vertex()[1]);
                    fTrackStartZ.push_back(track.Vertex()[2]);
                    fTrackStartPX.push_back(track.Momentum_beg()*track.VtxDir()[0]);
                    fTrackStartPY.push_back(track.Momentum_beg()*track.VtxDir()[1]);
                    fTrackStartPZ.push_back(track.Momentum_beg()*track.VtxDir()[2]);
        
                    fTrackEndX.push_back(track.End()[0]);
                    fTrackEndY.push_back(track.End()[1]);
                    fTrackEndZ.push_back(track.End()[2]);
                    fTrackEndPX.push_back(track.Momentum_end()*track.EndDir()[0]);
                    fTrackEndPY.push_back(track.Momentum_end()*track.EndDir()[1]);
                    fTrackEndPZ.push_back(track.Momentum_end()*track.EndDir()[2]);
        
                    fTrackLenF.push_back(track.LengthForward());
                    fTrackLenB.push_back(track.LengthBackward());
        
                    if (findIonization.isValid()) {
                        // No calibration for now.  Someday this should all be in reco
                        rec::TrackIoniz ionization = *(findIonization.at(iTrack));
                        float avgIonF, avgIonB;
                        gar::anatree::processIonizationInfo(ionization, fIonizTruncate, avgIonF, avgIonB);
                        fTrackAvgIonF.push_back( avgIonF );
                        fTrackAvgIonB.push_back( avgIonB );
                    } else {
                        // must push_back something so that fTrackAvgIonF,B are of correct size.
                        fTrackAvgIonF.push_back( 0.0 );
                        fTrackAvgIonB.push_back( 0.0 );
                    }
        
                    fNTPCClustersOnTrack.push_back(track.NHits());
        
                    if (fWriteTPCClustersInTracks) {
                        size_t iTrackedTPCCluster=0;
                        for (; iTrackedTPCCluster<track.NHits(); iTrackedTPCCluster++) {
                            auto const& TPCCluster = *(findManyTPCClusters.at(iTrack).at(iTrackedTPCCluster));
                            fTrkTPCClusterX.push_back(TPCCluster.Position()[0]);
                            fTrkTPCClusterY.push_back(TPCCluster.Position()[1]);
                            fTrkTPCClusterZ.push_back(TPCCluster.Position()[2]);
                            fTrkTPCClusterSig.push_back(TPCCluster.Signal());
                            fTrkTPCClusterRMS.push_back(TPCCluster.RMS());
                            fTrkTPCClusterTrkIndex.push_back((Int_t)iTrack);
                        }
                    }
                    iTrack++;
                } // end loop over TrackHandle
            }
        
            // save Vertex and Track-Vertex association info
            if (fWriteVertsNtracks) {
                const art::FindManyP<gar::rec::TPCCluster> findManyTPCClusters(TrackHandle,e,fTrackLabel);
                const art::FindOneP<gar::rec::TrackIoniz>  findIonization(TrackHandle,e,fTrackLabel);
                const art::FindMany<gar::rec::Track, gar::rec::TrackEnd> findManyTrack(VertexHandle,e,fVertexLabel);
                size_t iVertex = 0;
                for ( auto const& vertex : (*VertexHandle) ) {
                    fVertexX.push_back(vertex.Position()[0]);
                    fVertexY.push_back(vertex.Position()[1]);
                    fVertexZ.push_back(vertex.Position()[2]);
                    fVertexT.push_back(vertex.Time());
        
                    int nVertexedTracks = 0;
                    if ( findManyTrack.isValid() ) {
                        // tracks is a vector of gar::rec::Track
                        auto const& tracks = findManyTrack.at(iVertex);
                        nVertexedTracks = tracks.size();
                    }
                    fVertexN.push_back(nVertexedTracks);
        
                    int vertexCharge = 0;
                    for (int iVertexedTrack=0; iVertexedTrack<nVertexedTracks; ++iVertexedTrack) {
                        fVTAssn_Vertex.push_back(iVertex);  // 1 push of iVertex for each track in vertex
        
                        // Get this vertexed track.  findManyTrack.at(iVertex) is a
                        // vector<const gar::rec::Track*, std::allocator<const gar::rec::Track*> >
                        gar::rec::Track track = *(findManyTrack.at(iVertex).at(iVertexedTrack));
        
                        // Get the end of the track in the vertex.  It isn't that odd for the end
                        // of the track not used in the vertex to be closer to the vertex than the
                        // one actually used; you might have a very short stub track in a 3 track
                        // vertex with small opening angles and the other 2 tracks might pull the
                        // vertex some distance towards the far end of the stub track
                        gar::rec::TrackEnd fee = *(findManyTrack.data(iVertex).at(iVertexedTrack));
        
                        if (fee==gar::rec::TrackEndBeg) {
                            fVTAssn_TrkEndPx.push_back( track.Momentum_beg()*track.VtxDir()[0] );
                            fVTAssn_TrkEndPy.push_back( track.Momentum_beg()*track.VtxDir()[1] );
                            fVTAssn_TrkEndPz.push_back( track.Momentum_beg()*track.VtxDir()[2] );
                            fVTAssn_TrkEndLen.push_back( track.LengthForward() );
                            fVTAssn_TrkEndQ.push_back( track.ChargeBeg() );
                            fVTAssn_TrkEndChiSq.push_back( track.ChisqForward() );
                            vertexCharge += track.ChargeBeg();
                        } else {
                            fVTAssn_TrkEndPx.push_back( track.Momentum_end()*track.EndDir()[0] );
                            fVTAssn_TrkEndPy.push_back( track.Momentum_end()*track.EndDir()[1] );
                            fVTAssn_TrkEndPz.push_back( track.Momentum_end()*track.EndDir()[2] );
                            fVTAssn_TrkEndLen.push_back( track.LengthBackward() );
                            fVTAssn_TrkEndQ.push_back( track.ChargeEnd() );
                            fVTAssn_TrkEndChiSq.push_back( track.ChisqBackward() );
                            vertexCharge += track.ChargeEnd();
                        }
        
                        // We have assns from Vert to Track and Track to Ionize (and TPCClusters).  At this
                        // point in the code we have a Track from a Vert of interest and want the Ionize for
                        // it.  But FindOne, FindMany, FindOneP and FindManyP can't be indexed by the values
                        // they store; the at() and get() and for that matter, the data() methods only have
                        // size_t based indexing.
                        int iTrack = 0;   bool found = false;
                        for ( auto const& soughtTrack : (*TrackHandle) ) {
                            if ( track==soughtTrack ) {
                                found = true;
                                break;
                            }
                            ++iTrack;
                        }
                        if (!found) {
                            throw cet::exception("anatree")
                            << " Fail to find vertexed track in TrackHandle "
                            << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
                        }
        
                        if (found && findIonization.isValid()) {
                            // No calibration for now
                            rec::TrackIoniz ionization = *(findIonization.at(iTrack));
                            float avgIonF, avgIonB;
                            gar::anatree::processIonizationInfo(ionization, fIonizTruncate, avgIonF, avgIonB);
                            fTrackAvgIonF.push_back( avgIonF );
                            fTrackAvgIonB.push_back( avgIonB );
                            if (fee==gar::rec::TrackEndBeg) {
                                fVTAssn_TrkEndIon.push_back( avgIonF );
                            } else {
                                fVTAssn_TrkEndIon.push_back( avgIonB );
                            }
                        } else {
                            fVTAssn_TrkEndIon.push_back( 0.0 );
                        }
        
                        if (found && findManyTPCClusters.isValid()) {
                            // TPCClusters is a vector of gar::rec::TPCCluster
                            auto const& TPCClusters = findManyTPCClusters.at(iTrack);
                            fVTAssn_TrkEndNclus.push_back( TPCClusters.size() );
                        } else {
                            fVTAssn_TrkEndNclus.push_back( 0 );
                        }
                    }
        
                    fVertexQ.push_back(vertexCharge);
                    ++iVertex;
                } // end loop over VertexHandle
            }
        
            // save calorimetry info
            if (fWriteCaloDigits) {
                //Save Digit Hit info
                for ( auto const& DigiHit : (*RawHitHandle) ) {
                    fDiginHits++;
                    fDigiHitX.push_back(DigiHit.X());
                    fDigiHitY.push_back(DigiHit.Y());
                    fDigiHitZ.push_back(DigiHit.Z());
                    fDigiHitTime.push_back( (DigiHit.Time().first + DigiHit.Time().second) / 2.0 );
                    fDigiHitADC.push_back(DigiHit.ADC());
                    fDigiHitCellID.push_back(DigiHit.CellID());
                }
            }
        
            if (fWriteCaloInfo) {
                //Save Reco Hit info
                for ( auto const& Hit : (*RecoHitHandle) ) {
                    fReconHits++;
                    fRecoHitX.push_back(Hit.Position()[0]);
                    fRecoHitY.push_back(Hit.Position()[1]);
                    fRecoHitZ.push_back(Hit.Position()[2]);
                    fRecoHitTime.push_back(Hit.Time());
                    fRecoHitEnergy.push_back(Hit.Energy());
                    fRecoHitCellID.push_back(Hit.CellID());
                    fRecoEnergySum += Hit.Energy();
                }
        
                // save Cluster info
                for ( auto const& Cluster : (*RecoClusterHandle) ) {
                    fnCluster++;
                    fClusterNhits.push_back(Cluster.CalorimeterHits().size());
                    fClusterEnergy.push_back(Cluster.Energy());
                    fClusterX.push_back(Cluster.Position()[0]);
                    fClusterY.push_back(Cluster.Position()[1]);
                    fClusterZ.push_back(Cluster.Position()[2]);
                    fClusterTheta.push_back(Cluster.ITheta());
                    fClusterPhi.push_back(Cluster.IPhi());
                    fClusterPID.push_back(Cluster.ParticleID());
                    // fClusterShape.push_back(Cluster.Shape());
                    fClusterMainAxisX.push_back(Cluster.EigenVectors()[0]);
                    fClusterMainAxisY.push_back(Cluster.EigenVectors()[1]);
                    fClusterMainAxisZ.push_back(Cluster.EigenVectors()[2]);
        
                    // Write matching track info
                    if (fWriteAllTracks) {
                        const art::FindManyP<gar::rec::Track, gar::rec::TrackEnd> findManyECALTrk(RecoClusterHandle,e,fECALAssnLabel);
                        int nECALedTracks = 0;
                        if ( findManyECALTrk.isValid() ) {
                            nECALedTracks = findManyECALTrk.at(fnCluster-1).size();
                        }
                        if (nECALedTracks==0) continue;    // If ECAL in generation vol, most clusters have no tracks
        
                        const art::FindOneP<gar::rec::TrackIoniz>  findIonization(TrackHandle,e,fTrackLabel);
                        const art::FindManyP<gar::rec::TPCCluster> findManyTPCClusters(TrackHandle,e,fTrackLabel);
                        for (int iECALedTrack=0; iECALedTrack<nECALedTracks; ++iECALedTrack) {
                            fECALAssn_Cluster.push_back(fnCluster-1);    // Cluster which this Assn belongs to
                            gar::rec::TrackEnd fee = *(findManyECALTrk.data(fnCluster-1).at(iECALedTrack));
                            fECALAssn_TrkWhich.push_back(fee);    // The gar::rec::TrackEnd (see Track.h) that extrapolated to cluster
                            gar::rec::Track track = *(findManyECALTrk.at(fnCluster-1).at(iECALedTrack));
                            fECALAssn_TrkStartX.push_back( track.Vertex()[0] );
                            fECALAssn_TrkStartY.push_back( track.Vertex()[1] );
                            fECALAssn_TrkStartZ.push_back( track.Vertex()[2] );
                            fECALAssn_TrkStartPX.push_back(track.Momentum_beg()*track.VtxDir()[0]);
                            fECALAssn_TrkStartPY.push_back(track.Momentum_beg()*track.VtxDir()[1]);
                            fECALAssn_TrkStartPZ.push_back(track.Momentum_beg()*track.VtxDir()[2]);
                            fECALAssn_TrkEndX.push_back(track.End()[0]);
                            fECALAssn_TrkEndY.push_back(track.End()[1]);
                            fECALAssn_TrkEndZ.push_back(track.End()[2]);
                            fECALAssn_TrkEndPX.push_back(track.Momentum_end()*track.EndDir()[0]);
                            fECALAssn_TrkEndPY.push_back(track.Momentum_end()*track.EndDir()[1]);
                            fECALAssn_TrkEndPZ.push_back(track.Momentum_end()*track.EndDir()[2]);
                            fECALAssn_TrkLenF.push_back(track.LengthForward());
                            fECALAssn_TrkLenB.push_back(track.LengthBackward());
        
                            // Loop again to find index of this track
                            int iTrack = 0;   bool found = false;
                            for ( auto const& soughtTrack : (*TrackHandle) ) {
                                if ( track==soughtTrack ) {
                                    found = true;
                                    break;
                                }
                                ++iTrack;
                            }
                            if (!found) {
                                throw cet::exception("anatree")
                                << " Fail to find ECAL-matched track in TrackHandle "
                                << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
                            }
        
                            if (found && findIonization.isValid()) {
                                rec::TrackIoniz ionization = *(findIonization.at(iTrack));
                                float avgIonF, avgIonB;
                                gar::anatree::processIonizationInfo(ionization, fIonizTruncate, avgIonF, avgIonB);
                                 fECALAssn_TrkAvgIonF.push_back( avgIonF );
                                 fECALAssn_TrkAvgIonB.push_back( avgIonB );
                            } else {
                                // must push_back something so that fTrackAvgIon is of correct size.
                                fECALAssn_TrkAvgIonF.push_back( 0.0 );
                                fECALAssn_TrkAvgIonB.push_back( 0.0 );
                            }
        
                            int nTrackedTPCClusters = 0;
                            if (found && findManyTPCClusters.isValid()) {
                                // TPCClusters is a vector of gar::rec::TPCCluster
                                auto const& TPCClusters = findManyTPCClusters.at(iTrack);
                                nTrackedTPCClusters = TPCClusters.size();
                            }
                            fECALAssn_TrkNTPCClustersOnTrack.push_back(nTrackedTPCClusters);
                        }
                    }
                }
        
        
            } // end branch on fWriteCaloInfo
        } // end gar::anatree::FillVectors



        //==========================================================================
        // Process ionization.  Eventually this moves into the reco code.
        void processIonizationInfo(rec::TrackIoniz& ion, float ionizeTruncate,
                                   float& forwardIonVal, float& backwardIonVal) {
            // Get the ADC data
            std::vector<std::pair<float,float>> SigData = ion.getFWD_dSigdXs();

            // NO CALIBRATION SERVICE FOR NOW

            forwardIonVal = processOneDirection(SigData, ionizeTruncate);

            SigData = ion.getBAK_dSigdXs();
            backwardIonVal = processOneDirection(SigData, ionizeTruncate);

            return;
        }



        float processOneDirection(std::vector<std::pair<float,float>> SigData, float ionizeTruncate) {

            std::vector<std::pair<float,float>> dEvsX;    // Will be the ionization vs distance along track

            // The first hit on the track never had its ionization info stored.  Not a problem
            // really.  Each pair is a hit and the step along the track that ends at the hit
            // For the last hit, just take the step from the n-1 hit; don't guess some distance to
            // (nonexistant!) n+1 hit.  Using pointer arithmetic because you are a real K&R C nerd!
            float distAlongTrack = 0;
            std::vector<std::pair<float,float>>::iterator littlebit = SigData.begin();
            for (; littlebit<(SigData.end()-1); ++littlebit) {
                float dE =   std::get<0>(*littlebit);
               // tpctrackfit2_module.cc fills the TrackIoniz data product so that 
               // this quantity is really dL > 0 not dX, a coordinate on the drift axis
                float dX  = std::get<1>(*littlebit);
                distAlongTrack += dX;    // But count full step to get hit position on track
                // Take dX to be 1/2 the previous + last segment
                dX += std::get<1>(*(littlebit+1));
                float dEdX = dE/(0.5*dX);

                std::pair pushme = std::make_pair(dEdX,distAlongTrack);
                dEvsX.push_back( pushme );
            }

            // Get the truncated mean; first sort then take mean
            std::sort(dEvsX.begin(),dEvsX.end(), lessThan_byE);

            // Get the dEdX vs length data, truncated.
            int goUpTo = ionizeTruncate * dEvsX.size() +0.5;
            if (goUpTo > (int)dEvsX.size()) goUpTo = dEvsX.size();
            int i = 1;        float returnvalue = 0;
            littlebit = dEvsX.begin();
            for (; littlebit<dEvsX.end(); ++littlebit) {
              returnvalue += std::get<0>(*littlebit);
              ++i;
              if (i>goUpTo) break;
            }
            returnvalue /= goUpTo;
            return returnvalue;
        }



        //==========================================================================
        // Coherent pion analysis specific code
        float computeT( simb::MCTruth theMCTruth ) {
            // Warning.  You probably want the absolute value of t, not t.
            int nPart = theMCTruth.NParticles();
            enum { nu, mu, pi};
            float E[3], Px[3], Py[3], Pz[3];
            E[nu] = E[mu] = E[pi] = -1e42;

            for (int i=0; i<3;++i) {
                Px[i] = 0; 
                Py[i] = 0;
                Pz[i] = 0;
                E[i]  = 0;
            }
            // Find t from the MCParticles via the
            for (int iPart=0; iPart<nPart; iPart++) {
                simb::MCParticle Part = theMCTruth.GetParticle(iPart);
                int code = Part.PdgCode();
                int mom  = Part.Mother();

                // get the neutrino
                if ( abs(code) == 12 || abs(code) == 14 || abs(code) == 16 ) {
                    if (mom == -1) {
                        E[nu] = Part.E();   Px[nu] = Part.Px();   Py[nu] = Part.Py();   Pz[nu] = Part.Pz();
                    }
                }

                // get the lepton
                if ( abs(code) == 11 || abs(code) == 13 || abs(code) == 15 ) {
                    if (mom == 0) {
                        E[mu] = Part.E();   Px[mu] = Part.Px();   Py[mu] = Part.Py();   Pz[mu] = Part.Pz();
                    }
                }

                // get the pion
                if ( code==111 || abs(code)==211 ) {
                    if (mom == 1) {
                        E[pi] = Part.E();   Px[pi] = Part.Px();   Py[pi] = Part.Py();   Pz[pi] = Part.Pz();
                    }
                }

                // get outa here
                if ( E[nu]!=0 && E[mu]!=0 && E[pi]!=0) break;

            }

            // Compute t; reuse nu 4-vector to get first q, then t.
            E[nu] -= E[mu];   Px[nu] -= Px[mu];   Py[nu] -= Py[mu];   Pz[nu] -= Pz[mu];
            E[nu] -= E[pi];   Px[nu] -= Px[pi];   Py[nu] -= Py[pi];   Pz[nu] -= Pz[pi];
            float t = E[nu]*E[nu] -Px[nu]*Px[nu] -Py[nu]*Py[nu] -Pz[nu]*Pz[nu];
            return t;
        }



//    } // end namespace ana
} // end namespace gar
