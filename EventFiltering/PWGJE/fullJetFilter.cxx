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
using selectedClusters = o2::soa::Filtered<o2::aod::EMCALClusters>;

static const std::vector<std::string> matchInfoNames {"Patch", "MatchedJet"}; //{"MB", "Patch", "MatchedJet"};

struct fullJetFilter {
  enum { //kMinimumBias = 0,
    kPatch = 0,
    kMatchedJet,
    kCategories };

  //event selection cuts
  Configurable<float> selectionHighPtTrack{"selectionHighPtTrack", 10., "Minimum track pT trigger threshold"};       //we want to keep all events having a track with pT above this
  Configurable<float> selectionJetPt{"selectionJetPt", 25., "Minimum full jet pT trigger threshold"}; //we want to keep all events having a charged jet with pT above this
  Configurable<float> selectionHighECluster{"selectionHighPtCluster", 5., "Minimum cluster E trigger threshold"}; //we want to keep all events having a charged jet with pT above this
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

  Configurable<std::string> mClusterDefinition{"clusterDefinition", "kV3Default", "cluster definition to be selected, e.g. V3Default"};

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
    spectra.add("fJetPtSelected", "pT of selected high pT jets", HistType::kTH1F, {{150, 0., +150., "jet #it{p}_{T} (GeV/#it{c})"}});
    spectra.add("fClusterESelected", "E of selected high E clusters", HistType::kTH1F, {{150, 0., +150., "cluster #it{E} (GeV/#it{c})"}});
    spectra.add("fClusterDistSelected", "Distance of selected  cluster to jet", HistType::kTH2F, {
      {150, 0., +150., "cluster #it{E} (GeV/#it{c})"},
      {100, 0., +10., "dist"}
      });

    spectra.add("fJetMB", "pT of MB jets", HistType::kTH1F, {{150, 0., +150., "jet #it{p}_{T} (GeV/#it{c})"}});
    spectra.add("fClusterMB", "E of MB clusters", HistType::kTH1F, {{150, 0., +150., "cluster #it{E} (GeV/#it{c})"}});

    spectra.add("fCluster", "Cluster", HistType::kTH3F, {
        {fgkNPtBins, kMinPt, kMaxPt, "#it{E} (GeV)"},
        {fgkNPhiBins, kMinPhi, kMaxPhi, "#phi"},
        {fgkNEtaBins, kMinEta, kMaxEta, "#eta"}
        }); // Cluster info
    spectra.add("fJet", "Jet properties", HistType::kTH3F, {
        {fgkNPtBins, kMinPt, kMaxPt, "#it{p}_T (GeV)"},
        {fgkNPhiBins, kMinPhi, kMaxPhi, "#phi"},
        {fgkNEtaBins, kMinEta, kMaxEta, "#eta"}
        }); // Jet pt, eta, phi before matching

    // Tracks #MB events, #Patches, #Matched in one histogram
    auto scalers{std::get<std::shared_ptr<TH1>>(spectra.add("fProcessedEvents", ";;Number of filtered events", HistType::kTH1F, {{kCategories, -0.5, kCategories - 0.5}}))};
    for (uint32_t iS{1}; iS <= matchInfoNames.size(); ++iS) {
      scalers->GetXaxis()->SetBinLabel(iS, matchInfoNames[iS - 1].data());
    }

    jetReclusterer.isReclustering = true;
    jetReclusterer.algorithm = fastjet::JetAlgorithm::cambridge_algorithm;
    jetReclusterer.jetR = cfg_jetR * 2.5; // Why is this *2.5?
  } // init()

  // Declare filters
  o2::aod::EMCALClusterDefinition clusDef = o2::aod::emcalcluster::getClusterDefinitionFromString(mClusterDefinition.value);
  Filter clusterDefinitionSelection = o2::aod::emcalcluster::definition == static_cast<int>(clusDef);

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgVertexCut;
  // TODO: Should be fiducial jet cut. I.e.: eta_jet < eta_max - R
  // Filter jetEtaFilter = (nabs(aod::jet::eta) < cfgJetEtaCut);
  // Filter jetPhiFilter = ((aod::jet::phi < cfgJetPhiMax) && (aod::jet::phi > cfgJetPhiMin));
  // Filter trackFilter = (nabs(aod::track::eta) < cfgTrackEtaCut) && (aod::track::isGlobalTrack == (uint8_t) true) && (aod::track::pt > cfgTrackLowPtCut);

  //using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection>>;

  void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision,
               aod::Jets const& jets,
               aod::JetClusterConstituents const& clusterConstituents,
               selectedClusters const& clusters
              )
  {
    // TODO:
    // Check if EMCAL trigger is fired. Skip otherwise.
    // -> Filter on collision table?
    // Construct patch proxy from cells
    // -> Maybe do this in a different task?
    // Match cluster with cells
    // Set flags
    // collision process loop
    bool keepEvent[kCategories]{false};
    int nMC = 0; // number of matched clusters
    int nKC = 0; // number of kept clusters
    // keepEvent[kMinimumBias] = true;
    for (const auto& jet : jets){
      spectra.fill(HIST("fJet"), 0., jet.phi(), jet.eta());
      spectra.fill(HIST("fJetMB"), jet.pt());
    }
    for (const auto& cluster : clusters){
      spectra.fill(HIST("fCluster"), cluster.energy(), cluster.phi(), cluster.eta());
      spectra.fill(HIST("fClusterMB"), cluster.energy());
    }

    for (auto clusterConstituent : clusterConstituents){
      auto cluster = clusterConstituent.cluster();
      // if (cluster.energy() > selectionHighECluster){
        keepEvent[kPatch] = true;
        spectra.fill(HIST("fClusterESelected"), cluster.energy());
        break;
      // }
    }

    for (const auto& cluster : clusters){
      if (cluster.energy() > selectionHighECluster){
        keepEvent[kPatch] = true;
        spectra.fill(HIST("fClusterESelected"), cluster.energy());
        nKC++;
        break;
      }
    }

    for (const auto& jet : jets){
      for (const auto& cluster : clusters){
        double dist = TMath::Sqrt(TMath::Power(cluster.phi() - jet.phi(), 2) + TMath::Power(cluster.eta() - jet.eta(), 2));
        if (dist > jet.r()) continue;
        if (cluster.energy() > selectionHighECluster){
          keepEvent[kPatch] = true;
          spectra.fill(HIST("fClusterDistSelected"), cluster.energy(), dist);
          nKC++;
          break;
        }
      }
    }

    for (const auto& jet : jets){
      if (jet.pt() < selectionJetPt) continue;
      if (TMath::Abs(jet.eta() > cfgJetEtaCut || jet.phi() < cfgJetPhiMin || jet.phi() > cfgJetPhiMax)) continue;
      keepEvent[kMatchedJet] = true;
      spectra.fill(HIST("fJetPtSelected"), jet.pt());
      nMC++;
    }

    for (int iDecision{0}; iDecision < kCategories; iDecision++){
      if (keepEvent[iDecision]) {
        spectra.fill(HIST("fProcessedEvents"), iDecision);
      }
    }
    // if (nKC - nMC != 0) std::cout << "nKC (" << nKC << ") - nMC (" << nMC << ") = " << nKC - nMC << std::endl;

    tags(keepEvent[0], keepEvent[1]);//, keepEvent[2]);
  } // process()
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{
    adaptAnalysisTask<fullJetFilter>(cfg)};
}
