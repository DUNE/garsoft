////////////////////////////////////////////////////////////////////////
/// \file  IonizationReadout_module.cc
/// \brief Create the readout from the ionization signal
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

#ifndef GAR_READOUTSIMULATION_IONIZATIONREADOUT
#define GAR_READOUTSIMULATION_IONIZATIONREADOUT

// C++ Includes
#include <memory>
#include <vector> // std::ostringstream
#include <iostream>

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/Assns.h"
#include "cetlib_except/exception.h"
#include "cetlib/search_path.h"

// nutools extensions
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nutools/ParticleNavigation/ParticleList.h"
#include "nutools/RandomUtils/NuRandomService.h"

// GArSoft Includes
#include "DetectorInfo/DetectorClocksService.h"
#include "ReadoutSimulation/IonizationAndScintillation.h"
#include "ReadoutSimulation/ElectronDriftStandardAlg.h"
#include "ReadoutSimulation/TPCReadoutSimStandardAlg.h"
#include "Utilities/AssociationUtil.h"
#include "SimulationDataProducts/EnergyDeposit.h"
#include "RawDataProducts/RawDigit.h"
#include "Geometry/Geometry.h"
#include "CoreUtils/ServiceUtil.h"

// Forward declarations

///Geant4 interface
namespace gar {
  namespace rosim {

    /**
     * @brief Runs Readout simulation including propagation of electrons and photons to readout
     *
     *
     * Randomness
     * -----------
     *
     * The random number generators used by this process are:
     * - 'GEANT' instance: used by Geant4
     * - 'propagation' instance: used in electron propagation
     *
     */
    class IonizationReadout : public ::art::EDProducer{
    public:

      /// Standard constructor and destructor for an FMWK module.
      explicit IonizationReadout(fhicl::ParameterSet const& pset);
      virtual ~IonizationReadout();

      void produce (::art::Event& evt);
      void beginJob();
      void beginRun(::art::Run& run);
      void reconfigure(fhicl::ParameterSet const& pset);

    private:

      void DriftElectronsToReadout(std::vector<sdp::EnergyDeposit> const& edepCol,
                                   std::vector<edepIDE>                 & edepIDEs);
      void CombineIDEs(std::vector<edepIDE>                 & edepIDEs,
                       std::vector<sdp::EnergyDeposit> const& edepCol);
      void CreateSignalDigit(unsigned int                                     const& channel,
                             std::vector<float>                                    & electrons,
                             std::set<size_t>                                      & eDepLocs,
                             std::vector<raw::RawDigit>                            & digCol,
                             ::art::ValidHandle< std::vector<sdp::EnergyDeposit> > & eDepCol,
                             ::art::Assns<sdp::EnergyDeposit, raw::RawDigit>       & erassn,
                             ::art::Event                                          & evt);
      void CheckChannelToEnergyDepositMapping(unsigned int       const& channel,
                                              sdp::EnergyDeposit const& edep,
                                              std::string        const& id);

      std::string                         fG4Label;    ///< label of G4 module
      std::unique_ptr<ElectronDriftAlg>   fDriftAlg;   ///< algorithm to drift ionization electrons
      const gar::detinfo::DetectorClocks* fTime;       ///< electronics clock
      std::unique_ptr<TPCReadoutSimAlg>   fROSimAlg;   ///< algorithm to simulate the electronics
      fhicl::ParameterSet                 fISCalcPars; ///< parameter set for the IS calculator
      size_t                              fNumTicks;   ///< number of TDC samples
      const gar::geo::GeometryCore*       fGeo;        ///< geometry information
      bool                                fCheckChan;  ///< flag to check mapping of energy deposits to channels
    };

  } // namespace rosim

  namespace rosim {

    //----------------------------------------------------------------------
    // Constructor
    IonizationReadout::IonizationReadout(fhicl::ParameterSet const& pset)
    {
      fTime  = gar::providerFrom<detinfo::DetectorClocksService>();

      fNumTicks = gar::providerFrom<detinfo::DetectorPropertiesService>()->NumberTimeSamples();

      fGeo = gar::providerFrom<geo::Geometry>();

      // setup the random number service
      // obtain the random seed from NuRandomService
      ::art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, "HepJamesRandom", "ionization", pset, "IonizationSeed");

      this->reconfigure(pset);

      produces< std::vector<raw::RawDigit>                      >();
      produces< ::art::Assns<sdp::EnergyDeposit, raw::RawDigit> >();

      return;
    }

    //----------------------------------------------------------------------
    // Destructor
    IonizationReadout::~IonizationReadout()
    {
    }

    //----------------------------------------------------------------------
    void IonizationReadout::reconfigure(fhicl::ParameterSet const& pset)
    {
      LOG_DEBUG("IonizationReadout") << "Debug: IonizationReadout()";
      ::art::ServiceHandle<::art::RandomNumberGenerator> rng;

      fISCalcPars = pset.get<fhicl::ParameterSet>("ISCalcPars"                 );
      fG4Label    = pset.get<std::string        >("G4ModuleLabel",      "geant");
      fCheckChan  = pset.get<bool               >("CheckChannelMapping", false );

      auto driftAlgPars = pset.get<fhicl::ParameterSet>("ElectronDriftAlgPars");
      auto driftAlgName = driftAlgPars.get<std::string>("DriftAlgType");

      if(driftAlgName.compare("Standard") == 0)
        fDriftAlg = std::make_unique<gar::rosim::ElectronDriftStandardAlg>(rng->getEngine("ionization"),
                                                                           driftAlgPars);
      else
        throw cet::exception("IonizationReadout")
	  << "Unable to determine which electron drift algorithm to use, bail";

      auto tpcROAlgPars = pset.get<fhicl::ParameterSet>("TPCReadoutSimAlgPars");
      auto tpcROAlgName = tpcROAlgPars.get<std::string>("TPCReadoutSimType");

      if(tpcROAlgName.compare("Standard") == 0)
        fROSimAlg = std::make_unique<gar::rosim::TPCReadoutSimStandardAlg>(rng->getEngine("ionization"),
                                                                           tpcROAlgPars);
      else
        throw cet::exception("IonizationReadout")
	  << "Unable to determine which TPC readout simulation algorithm to use, bail";

      return;
    }

    //----------------------------------------------------------------------
    void IonizationReadout::beginJob()
    {
      auto* rng = &*(::art::ServiceHandle<::art::RandomNumberGenerator>());

      // create the ionization and scintillation calculator;
      // this is a singleton (!) so we just need to make the instance in one
      // location
      IonizationAndScintillation::CreateInstance(rng->getEngine("ionization"),
                                                 fISCalcPars);

      return;
    }

    //--------------------------------------------------------------------------
    void IonizationReadout::beginRun(::art::Run& run)
    {
      return;
    }

    //--------------------------------------------------------------------------
    void IonizationReadout::produce(::art::Event& evt)
    {
      LOG_DEBUG("IonizationReadout") << "produce()";

      // loop over the lists and put the particles and voxels into the event as collections
      std::unique_ptr< std::vector<raw::RawDigit>                      > rdCol (new std::vector<raw::RawDigit>                     );
      std::unique_ptr< ::art::Assns<sdp::EnergyDeposit, raw::RawDigit> > erassn(new ::art::Assns<sdp::EnergyDeposit, raw::RawDigit>);

      // first get the energy deposits from the event record
      auto eDepCol = evt.getValidHandle< std::vector<sdp::EnergyDeposit> >(fG4Label);

      if(eDepCol->size() > 0){

        std::vector<edepIDE> eDepIDEs;

	// drift the ionization electrons to the readout and create edepIDE objects
	this->DriftElectronsToReadout(*eDepCol, eDepIDEs);

	//if (eDepIDEs.size()>0)
	{

	  // The IDEs should have been combined already so that there are no repeat
	  // TDC values for any channel
	  unsigned int       prevChan = eDepIDEs.front().Channel;
	  std::set<size_t>   digitEDepLocs;
	  std::vector<float> electrons(fNumTicks, 0.);

	  // make the signal raw digits and set their associations to the energy deposits
	  for(auto edide : eDepIDEs){

	    LOG_DEBUG("IonizationReadout")
	      << "Current eDepIDE channel is "
	      << edide.Channel
	      << " previous channel is "
	      << prevChan;

	    if(edide.Channel != prevChan){
	      LOG_DEBUG("IonizationReadout")
		<< "There are  "
		<< digitEDepLocs.size()
		<< " locations for "
		<< edide.Channel
		<< " rdCol size is currently "
		<< rdCol->size();

	      // this method clears the electrons and digitEDepLocs collections
	      // after creating the RawDigit
	      this->CreateSignalDigit(prevChan,
				      electrons,
				      digitEDepLocs,
				      *rdCol,
				      eDepCol,
				      *erassn,
				      evt);

	      // reset the previous channel info
	      prevChan = edide.Channel;

	    }

	    // put overflow times in the last bin.  Is this okay?  TODO

	    size_t esize = electrons.size();
	    if (esize>0)
	      {
		if (edide.TDC >= esize)
		  {
		    electrons[esize - 1] = edide.NumElect;
		  }
		else
		  {
		    electrons[edide.TDC] = edide.NumElect;
		  }
	      }

	    for(auto loc : edide.edepLocs) digitEDepLocs.insert(loc);

	  } // end loop to fill signal raw digit vector and make EnergyDeposit associations

	  // still one more digit to make because we ran out of channels to compare against
	  this->CreateSignalDigit(eDepIDEs.back().Channel,
				  electrons,
				  digitEDepLocs,
				  *rdCol,
				  eDepCol,
				  *erassn,
				  evt);

	  LOG_DEBUG("IonizationReadout")
	    << "Created "
	    << rdCol->size()
	    << " raw digits from signal";

	  // now make the noise digits
	  // to do -- only make noise digits on channels we haven't
	  // yet considered for noise digits, but which may have
	  // been entirely zero-suppressed -- may need to keep a
	  // list of channels and pass it in
	  fROSimAlg->CreateNoiseDigits(*rdCol);

	} // end if the EdepIDEs have any size

      } // end if there were energy deposits to use

      evt.put(std::move(rdCol));
      evt.put(std::move(erassn));

      return;
    } // IonizationReadout::produce()

    //--------------------------------------------------------------------------
    void IonizationReadout::DriftElectronsToReadout(std::vector<sdp::EnergyDeposit> const& edepCol,
                                                    std::vector<edepIDE>                 & edepIDEs)
    {
      auto geo = gar::providerFrom<geo::Geometry>();

      float        xyz[3] = {0.};
      unsigned int chan   = 0;

      // now instantiate an ElectronDriftInfo object to keep track of the
      // drifted locations of each electron cluster from each energy deposit
      rosim::ElectronDriftInfo driftInfo;

      // loop over the energy deposits
      for(size_t e = 0; e < edepCol.size(); ++e){

        // get the positions, arrival times, and electron cluster sizes
        // for this energy deposition
        fDriftAlg->DriftElectronsToReadout(edepCol[e], driftInfo);

        auto clusterXPos = driftInfo.ClusterXPos();
        auto clusterYPos = driftInfo.ClusterYPos();
        auto clusterZPos = driftInfo.ClusterZPos();
        auto clusterTime = driftInfo.ClusterTime();
        auto clusterSize = driftInfo.ClusterSize();

        // the vectors should all have the same size by the time we get them
        // here (verified by the ElectronDriftInfo object when they are filled)
        for(size_t c = 0; c < clusterXPos.size(); ++c){

          xyz[0] = clusterXPos[c];
          xyz[1] = clusterYPos[c];
          xyz[2] = clusterZPos[c];

          // see if this cluster of electrons can be mapped to a channel.
          // if not, move on to the next one
          try{
            chan = geo->NearestChannel(xyz);
          }
          catch(cet::exception &except){
            continue;
          }

          edepIDEs.emplace_back(clusterSize[c],
                                chan,
                                fTime->TPCG4Time2TDC(clusterTime[c]),
                                e);

          LOG_DEBUG("IonizationReadout")
	    << "cluster time: "
	    << clusterTime[c]
	    << " TDC "
	    << fTime->TPCG4Time2TDC(clusterTime[c])
	    << " "
	    << fTime->G4ToElecTime(clusterTime[c])
	    << " "
	    << fTime->TPCClock().TickPeriod();

          this->CheckChannelToEnergyDepositMapping(edepIDEs.back().Channel,
                                                   edepCol[e],
                                                   "DriftElectronsToReadout");

        }

      } // end loop over deposit collections

      this->CombineIDEs(edepIDEs, edepCol);

      return;
    }

    //--------------------------------------------------------------------------
    void IonizationReadout::CombineIDEs(std::vector<edepIDE>                 & edepIDEs,
                                        std::vector<sdp::EnergyDeposit> const& edepCol)
    {
      LOG_DEBUG("IonizationReadout")
	<< "starting with "
	<< edepIDEs.size()
	<< " energy deposits";

      std::vector<edepIDE> temp;

      // sort the edepIDE objects.  This is the sorting by channel and TDC
      // because that is all the < operator of edepIDEs does
      std::sort(edepIDEs.begin(), edepIDEs.end());

      for(auto itr : edepIDEs){
        for(auto edloc : itr.edepLocs)
          this->CheckChannelToEnergyDepositMapping(itr.Channel,
                                                   edepCol[edloc],
                                                   "CombineIDEsAfterSort");
      }

      edepIDE prev = edepIDEs.front();
      edepIDE sum(edepIDEs.front());
      edepIDE cur(edepIDEs.front());

      // now loop over the sorted vector and combine any edepIDEs
      // with the same channel and tdc values for the IDEs
      // start with entry 1 as sum is already holding the information
      // from entry 0
      for(size_t e = 1; e < edepIDEs.size(); ++e){
        cur = edepIDEs[e];

        LOG_DEBUG("IonizationReadout")
	  << "current edepIDE: "
	  << cur.NumElect
	  << " "
	  << cur.Channel
	  << " "
	  << cur.TDC
	  << " "
	  << cur.edepLocs.size();

        if(cur != prev){
          LOG_DEBUG("IonizationReadout")
	    << "storing edepIDE sum: "
	    << sum.NumElect
	    << " "
	    << sum.Channel
	    << " "
	    << sum.TDC
	    << " "
	    << sum.edepLocs.size();

          if(fCheckChan)
            for(auto edloc : sum.edepLocs)
              this->CheckChannelToEnergyDepositMapping(sum.Channel,
                                                       edepCol[edloc],
                                                       "CombineIDEsStore");

          // put the summed edepIDE into the temp vector
          temp.push_back(sum);

          // start over with a fresh sum
          sum  = cur;
          prev = cur;
        }
        else{
          LOG_DEBUG("IonizationReadout")
	    << "summing current edepIDE";

          sum  += cur;
          prev  = cur;
        }

      } // end loop to sum edepIDEs

      // now swap the input vector with the temp vector
      temp.swap(edepIDEs);

      LOG_DEBUG("IonizationReadout")
	<< "ending with "
	<< edepIDEs.size()
	<< " energy deposits";

      return;
    }

    //--------------------------------------------------------------------------
    void IonizationReadout::CreateSignalDigit(unsigned int                                     const& channel,
                                              std::vector<float>                                    & electrons,
                                              std::set<size_t>                                      & eDepLocs,
                                              std::vector<raw::RawDigit>                            & digCol,
                                              ::art::ValidHandle< std::vector<sdp::EnergyDeposit> > & eDepCol,
                                              ::art::Assns<sdp::EnergyDeposit, raw::RawDigit>       & erassn,
                                              ::art::Event                                          & evt)
    {

      // could be that all the adc's fell below threshold, so test if we got a raw digit at all.

      bool todrop=false;
      raw::RawDigit tmpdigit = fROSimAlg->CreateRawDigit(channel, electrons, todrop);
      if (!todrop)
	{
	  digCol.emplace_back(tmpdigit);

	  LOG_DEBUG("IonizationReadout")
	    << "Associating "
	    << eDepLocs.size()
	    << " energy deposits to digit for channel "
	    << channel;

	  // loop over the locations in the eDepCol to make the associations
	  for(auto ed : eDepLocs){
	    auto const ptr = art::Ptr<sdp::EnergyDeposit>(eDepCol, ed);

	    this->CheckChannelToEnergyDepositMapping(channel,
						     *ptr,
						     "CreateSignalDigit");

	    util::CreateAssn(*this, evt, digCol, ptr, erassn);
	  }
	}

      eDepLocs.clear();
      electrons.clear();
      electrons.resize(fNumTicks, 0);

      return;
    }

    //--------------------------------------------------------------------------
    void IonizationReadout::CheckChannelToEnergyDepositMapping(unsigned int       const& channel,
                                                               sdp::EnergyDeposit const& edep,
                                                               std::string        const& id)
    {
      if(!fCheckChan) return;

      // check that the channel for this cluster is close to what we expect
      // for the energy deposit

      float xyz[3] = {0.};
      fGeo->ChannelToPosition(channel, xyz);

      if(std::abs(edep.Y() - xyz[1]) > 1 ||
         std::abs(edep.Z() - xyz[2]) > 1){
        LOG_VERBATIM("IonizationReadout")
	  << "In function "
	  << id
	  << ": Channel "
	  << channel
	  << " is off from the energy deposit: ("
	  << xyz[1]
	  << ", "
	  << xyz[2]
	  << ") vs ("
	  << edep.Y()
	  << ", "
	  << edep.Z()
	  << ") ";
      }

      return;
    }

  } // namespace rosim

  namespace rosim {

    DEFINE_ART_MODULE(IonizationReadout)

  } // namespace rosim
} // gar
#endif // GAR_READOUTSIMULATION_IONIZATIONREADOUT
