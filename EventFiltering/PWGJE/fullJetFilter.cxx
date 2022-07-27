// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
// O2 includes

#include "ReconstructionDataFormats/Track.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
//FK #include "Common/Core/PID/PIDResponse.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGJE/Core/JetFinder.h"

// #include "Commmon/CCDB/TriggerAliases.h"
#include "/Users/gijsvanweelden/alice/O2/DataFormats/Detectors/EMCAL/include/DataFormatsEMCAL/AnalysisCluster.h"
#include "/Users/gijsvanweelden/alice/O2/DataFormats/Detectors/EMCAL/include/DataFormatsEMCAL/Cluster.h"

#include "../filterTables.h"

#include "Framework/HistogramRegistry.h"

#include <cmath>
#include <string>
#include <TMath.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

static const std::vector<std::string> matchInfoNames {"MinimumBias", "Patch", "MatchedJet"};

struct fullJetFilter {
  enum { kMatchedJet = 0,
         kCategories };

  //event selection cuts
  Configurable<float> selectionHighPtTrack{"selectionHighPtTrack", 10., "Minimum track pT trigger threshold"};       //we want to keep all events having a track with pT above this
  Configurable<float> selectionJetPt{"selectionJetPt", 25., "Minimum charged jet pT trigger threshold"}; //we want to keep all events having a charged jet with pT above this
  Configurable<int> fMatchDist{"fMatchDist", 5, "Matching distance in number of towers"}; // Maximum clusters/patch distance
  // Acceptance cuts
  Configurable<float> cfgVertexCut{"cfgVertexCut", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgTrackEtaCut{"cfgTrackEtaCut", 0.9f, "Eta range for tracks"};
  Configurable<float> cfgTrackLowPtCut{"cfgTrackLowPtCut", 0.1f, "Minimum constituent pT"};
  // TODO: Complete jet should be inside the EMCal. Use fiducial cuts.
  // EMCAL: phi (1.40, 3.26), |eta| < 0.7
  // DCAL: [phi (320, 327) deg, |eta| < 0.7] && [phi (4.54, 5.70), 0.23 < |eta| < 0.7]
  Configurable<float> cfgJetEtaCut{"cfgJetEtaCut", 0.7f, "Maximum (absolute) jet eta, dimensions of EMCal"};
  Configurable<float> cfgJetPhiMin{"cfgJetPhiMin", 1.396f, "Minimum jet phi, dimensions of EMCal"}; // 80 degrees
  Configurable<float> cfgJetPhiMax{"cfgJetPhiMax", 3.264f, "Maximum jet phi, dimensions of EMCal"}; // 187 degrees
  // Soft drop settings
  Configurable<float> cfg_zCut{"cfg_zCut", 0.1, "soft drop z cut"};
  Configurable<float> cfg_beta{"cfg_beta", 0.0, "soft drop beta"};
  Configurable<float> cfg_jetR{"cfg_jetR", 0.4, "jet resolution parameter"}; //possible to get configurable from another task? jetR
  Configurable<bool> cfg_DoConstSub{"cfg_DoConstSub", false, "do constituent subtraction"};

  Produces<aod::FJetFilters> tags;

  std::vector<fastjet::PseudoJet> jetConstituents;
  std::vector<fastjet::PseudoJet> jetReclustered;
  JetFinder jetReclusterer;
  HistogramRegistry spectra{"spectra", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  void init(o2::framework::InitContext&)
  {
    // Bins here
    Int_t fgkNPtBins = 200; // Should go from -50 to 150 for corr pt jet
    Float_t kMinPt   = 0.;
    Float_t kMaxPt   = 200.;

    Int_t fgkNPhiBins = 18*8;
    Float_t kMinPhi   = 0.;
    Float_t kMaxPhi   = 2.*TMath::Pi();

    Int_t fgkNEtaBins = 100;
    Float_t kMinEta = -1.;
    Float_t kMaxEta =  1.;

    // Add histograms
    spectra.add("fCollZpos", "collision z position", HistType::kTH1F, {{200, -20., +20., "#it{z}_{vtx} position (cm)"}});
    spectra.add("fTrackPtSelected", "pT of selected high pT tracks", HistType::kTH1F, {{150, 0., +150., "track #it{p}_{T} (GeV/#it{c})"}});
    spectra.add("fJetChPtSelected", "pT of selected high pT charged jets", HistType::kTH1F, {{150, 0., +150., "charged jet #it{p}_{T} (GeV/#it{c})"}});
    //------------------------
    // Minimum Bias histograms
    //------------------------
        spectra.add("fMBTrack", "MBTracks", HistType::kTH3F, {
        {150, 0., +150., "track #it{p}_{T} (GeV/#it{c})"},
        {fgkNPhiBins, kMinPhi, kMaxPhi, "#phi"},
        {fgkNEtaBins, kMinEta, kMaxEta, "#eta"}
        });
    spectra.add("fMBJet", "MB Jet properties", HistType::kTH3F, {
        {fgkNPtBins, kMinPt, kMaxPt, "#it{E} (GeV)"},
        {fgkNPhiBins, kMinPhi, kMaxPhi, "#phi"},
        {fgkNEtaBins, kMinEta, kMaxEta, "#eta"}
        });
    // Jet substructure
    spectra.add("fMBJetZg", "MB Jet zg", HistType::kTH2F, {
        {fgkNPtBins, kMinPt, kMaxPt, "#it{p}_{T}^{jet} (GeV)"},
        {20, 0., 1., "#it{z}_{g}"}
        });
    spectra.add("fMBJetRg", "MB Jet Rg", HistType::kTH2F, {
        {fgkNPtBins, kMinPt, kMaxPt, "#it{p}_{T}^{jet} (GeV)"},
        {20, 0., 1., "#it{R}_{g}"}
        });
    spectra.add("fMBJetnsd", "MB Jet nsd", HistType::kTH2F, {
        {fgkNPtBins, kMinPt, kMaxPt, "#it{p}_{T}^{jet} (GeV)"},
        {50, 0, 50, "#it{n}_{SD}"}
        });
    //------------------------
    // Selection histograms
    //------------------------
    spectra.add("fTrack", "Tracks", HistType::kTH3F, {
        {150, 0., +150., "track #it{p}_{T} (GeV/#it{c})"},
        {fgkNPhiBins, kMinPhi, kMaxPhi, "#phi"},
        {fgkNEtaBins, kMinEta, kMaxEta, "#eta"}
        });
    spectra.add("fMaxPatch", "Max Patch", HistType::kTH3F, {
        {fgkNPtBins, kMinPt, kMaxPt, "#it{E} (GeV)"},
        {fgkNPhiBins, kMinPhi, kMaxPhi, "#phi"},
        {fgkNEtaBins, kMinEta, kMaxEta, "#eta"}
        }); // Max Patch info
    spectra.add("fCluster", "Cluster", HistType::kTH3F, {
        {fgkNPtBins, kMinPt, kMaxPt, "#it{E} (GeV)"},
        {fgkNPhiBins, kMinPhi, kMaxPhi, "#phi"},
        {fgkNEtaBins, kMinEta, kMaxEta, "#eta"}
        }); // Cluster info
    spectra.add("fJet", "Jet properties", HistType::kTH3F, {
        {fgkNPtBins, kMinPt, kMaxPt, "#it{E} (GeV)"},
        {fgkNPhiBins, kMinPhi, kMaxPhi, "#phi"},
        {fgkNEtaBins, kMinEta, kMaxEta, "#eta"}
        }); // Jet pt, eta, phi before matching
    spectra.add("fMBJet", "MB Jet properties", HistType::kTH3F, {
        {fgkNPtBins, kMinPt, kMaxPt, "#it{E} (GeV)"},
        {fgkNPhiBins, kMinPhi, kMaxPhi, "#phi"},
        {fgkNEtaBins, kMinEta, kMaxEta, "#eta"}
        }); // Jet pt, eta, phi before matching
    spectra.add("fJetTrigger", "Match info", HistType::kTH3F, {
        {fgkNPtBins, kMinPt, kMaxPt, "#it{E}^{p} (GeV)"},
        {fgkNPtBins, kMinPt, kMaxPt, "#it{E}^{cl} (GeV)"},
        {fgkNPtBins, kMinPt, kMaxPt, "#it{p}_{T}^{jet} (GeV)"}
        }); // Match info
    spectra.add("fJetZg", "Jet zg", HistType::kTH2F, {
        {fgkNPtBins, kMinPt, kMaxPt, "#it{p}_{T}^{jet} (GeV)"},
        {20, 0., 1., "#it{z}_{g}"}
        });
    spectra.add("fJetRg", "Jet Rg", HistType::kTH2F, {
        {fgkNPtBins, kMinPt, kMaxPt, "#it{p}_{T}^{jet} (GeV)"},
        {20, 0., 1., "#it{R}_{g}"}
        });
    spectra.add("fJetnsd", "Jet nsd", HistType::kTH2F, {
        {fgkNPtBins, kMinPt, kMaxPt, "#it{p}_{T}^{jet} (GeV)"},
        {50, 0, 50, "#it{n}_{SD}"}
        });
    // Tracks #MB events, #Patches, #Matched in one histogram
    auto scalers{std::get<std::shared_ptr<TH1>>(spectra.add("fProcessedEvents", ";;Number of filtered events", HistType::kTH1F, {{kCategories, -0.5, kCategories + 0.5}}))};
    for (uint32_t iS{1}; iS <= matchInfoNames.size(); ++iS) {
      scalers->GetXaxis()->SetBinLabel(iS, matchInfoNames[iS - 1].data());
    }

    jetReclusterer.isReclustering = true;
    jetReclusterer.algorithm = fastjet::JetAlgorithm::cambridge_algorithm;
    jetReclusterer.jetR = cfg_jetR * 2.5; // Why is this *2.5?
  } // init()

  // Declare filters
  Filter collisionFilter = nabs(aod::collision::posZ) < cfgVertexCut;
  // TODO: Should be fiducial jet cut. I.e.: eta_jet < eta_max - R
  Filter jetEtaFilter = (nabs(aod::jet::eta) < cfgJetEtaCut);
  Filter jetPhiFilter = ((aod::jet::phi < cfgJetPhiMax) && (aod::jet::phi > cfgJetPhiMin));
  Filter trackFilter = (nabs(aod::track::eta) < cfgTrackEtaCut) && (aod::track::isGlobalTrack == (uint8_t) true) && (aod::track::pt > cfgTrackLowPtCut);

  //using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection>>;

  void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision,
               //TrackCandidates const& tracks,
               aod::Jets const& jets
              //  aod::Tracks const& tracks,
              //  aod::JetTrackConstituents const& constituents,
              //  aod::JetConstituentsSub const& constituentsSub
              )
  {
    // collision process loop
    bool keepEvent[kCategories]{false};

    // TODO: Check that there is an emcal trigger for this event
    // Can we implement this as a filter on the event table?

    // selectedClusters const& clusters
    // for (const auto& cluster : clusters) {
    for (auto& jet : jets)
    {
        if (jet.pt() > selectionJetPt){
          spectra.fill(HIST("fJet"), jet.pt());
          keepEvent[kMatchedJet] = true;
          break;
        }
    }

    //count events which passed the selections
    for (int iDecision{0}; iDecision < kCategories; ++iDecision) {
      if (keepEvent[iDecision]) {
        spectra.fill(HIST("fProcessedEvents"), iDecision);
      }
    }

    tags(keepEvent[0]);
  } // process()
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{
    adaptAnalysisTask<fullJetFilter>(cfg)};
}
