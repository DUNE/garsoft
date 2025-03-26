#ifndef TRACK_CREATOR_H
#define TRACK_CREATOR_H 1

#include "Api/PandoraApi.h"
#include "ReconstructionDataProducts/Track.h"
#include "RotationTransformation.h"
#include "art/Framework/Principal/Event.h"

namespace gar {
  namespace gar_pandora {

    typedef std::vector<art::Ptr<gar::rec::Track>> TrackVector;
    typedef std::vector<gar::rec::Track> RawTrackVector;
    typedef std::set<const gar::rec::Track*> TrackList;

    class TrackCreator {
    public:
      class Settings {
      public:
        Settings();

        std::string m_trackCollection; ///< The reconstructed track collection
        std::string m_V0Collection;    ///< The vees collection

        unsigned int m_minTrackHits; ///< Track quality cut: the minimum number of track hits
        unsigned int m_maxTrackHits; ///< Track quality cut: the maximum number of track hits

        float
          m_d0TrackCut; ///< Track d0 cut used to determine whether track can be used to form pfo
        float
          m_z0TrackCut; ///< Track z0 cut used to determine whether track can be used to form pfo
        float m_unmatchedVertexTrackMaxEnergy;
        float m_d0UnmatchedVertexTrackCut;
        float m_z0UnmatchedVertexTrackCut;

        float
          m_minTrackECalDistanceFromIp; ///< Sanity check on separation between ip and track projected ecal position
        float m_maxTrackSigmaPOverP; ///< Track fraction momentum error cut

        float m_bField;                ///< The bfield
        int m_eCalBarrelInnerSymmetry; ///< ECal barrel inner symmetry order
        float m_eCalBarrelInnerPhi0;   ///< ECal barrel inner phi 0
        float m_eCalBarrelInnerR;      ///< ECal barrel inner radius
        float m_eCalEndCapInnerZ;      ///< ECal endcap inner z
        float m_GArCenterY;
        float m_GArCenterZ;
      };

      TrackCreator(const Settings& settings,
                   const pandora::Pandora* const pPandora,
                   const RotationTransformation* const pRotation);

      virtual ~TrackCreator();

      pandora::StatusCode CreateTracks(const art::Event& pEvent);
      pandora::StatusCode CreateTrackAssociations(const art::Event& pEvent);

      const TrackVector& GetTrackVector() const;

      void Reset();

    private:
      pandora::StatusCode CreateTracks();
      void GetTrackStates(const gar::rec::Track* const pTrack,
                          PandoraApi::Track::Parameters& trackParameters) const;
      bool PassesQualityCuts(const gar::rec::Track* const pTrack,
                             const PandoraApi::Track::Parameters& trackParameters) const;
      void TrackReachesECAL(const gar::rec::Track* const pTrack,
                            PandoraApi::Track::Parameters& trackParameters) const;
      void DefineTrackPfoUsage(const gar::rec::Track* const pTrack,
                               PandoraApi::Track::Parameters& trackParameters) const;
      void CalculateTrackStateAtCalo(const gar::rec::Track* const pTrack,
                                     pandora::CartesianVector& posAtCalo) const;
      void CalculateTimeAtCalo(const gar::rec::Track* const pTrack, float& timeAtCalo) const;

      pandora::StatusCode CollectTracks(const art::Event& pEvent,
                                        const std::string& label,
                                        TrackVector& trkVector);

      pandora::StatusCode ExtractKinks(const art::Event& pEvent);
      pandora::StatusCode ExtractProngsAndSplits(const art::Event& pEvent);
      pandora::StatusCode ExtractV0s(const art::Event& pEvent);

      bool IsV0(const gar::rec::Track* const pTrack) const;
      bool IsParent(const gar::rec::Track* const pTrack) const;
      bool IsDaughter(const gar::rec::Track* const pTrack) const;

      const Settings m_settings;
      const pandora::Pandora& m_pandora;
      const RotationTransformation& m_rotation;

      TrackVector artTrkVector;
      TrackVector m_trackVector;
      TrackList m_v0TrackList;       ///< The list of v0 tracks
      TrackList m_parentTrackList;   ///< The list of parent tracks
      TrackList m_daughterTrackList; ///< The list of daughter tracks
    };

    //------------------------------------------------------------------------------------------------------------------------------------------

    inline const TrackVector& TrackCreator::GetTrackVector() const
    {
      return m_trackVector;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    inline void TrackCreator::Reset()
    {
      artTrkVector.clear();
      m_trackVector.clear();
      m_v0TrackList.clear();
      m_parentTrackList.clear();
      m_daughterTrackList.clear();
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    inline bool TrackCreator::IsV0(const gar::rec::Track* const pTrack) const
    {
      return (m_v0TrackList.end() != m_v0TrackList.find(pTrack));
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    inline bool TrackCreator::IsParent(const gar::rec::Track* const pTrack) const
    {
      return (m_parentTrackList.end() != m_parentTrackList.find(pTrack));
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    inline bool TrackCreator::IsDaughter(const gar::rec::Track* const pTrack) const
    {
      return (m_daughterTrackList.end() != m_daughterTrackList.find(pTrack));
    }
  }
}

#endif // #ifndef TRACK_CREATOR_H
