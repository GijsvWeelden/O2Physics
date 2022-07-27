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

// jet trigger QA task (subscribing to jet finder task)
//
// Author: Gijs van Weelden
//

#ifndef O2_ANALYSIS_JETTRIGGER_QA_H
#define O2_ANALYSIS_JETTRIGGER_QA_H

#include "TH1F.h"
#include "TTree.h"

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/RunningWorkflowInfo.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/Core/JetFinder.h"
#include "EventFiltering/filterTables.h"
// #include "EventFiltering/PWGJE/jetFilter.cxx"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/Core/JetFinder.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
namespace jettriggerqa
{
DECLARE_SOA_COLUMN(Zg, zg, float);
DECLARE_SOA_COLUMN(Rg, rg, float);
DECLARE_SOA_COLUMN(Nsd, nsd, float);
} // namespace jettriggerqa
DECLARE_SOA_TABLE(JetTriggerQA, "AOD", "JETTRIGGERQA", jettriggerqa::Zg, jettriggerqa::Rg, jettriggerqa::Nsd);
} // namespace o2::aod

struct JetTriggerQA {
  Produces<aod::JetTriggerQA> jetTriggerQA;
  HistogramRegistry spectra{"spectra", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  Configurable<float> f_jetPtMin{"f_jetPtMin", 0.0, "minimum jet pT cut"};
  Configurable<float> f_SD_zCut{"f_SD_zCut", 0.1, "soft drop z cut"};
  Configurable<float> f_SD_beta{"f_SD_beta", 0.0, "soft drop beta"};
  Configurable<float> f_jetR{"f_jetR", 0.4, "jet resolution parameter"}; //possible to get configurable from another task? jetR. Or use jet.r()?

  Configurable<std::vector<float>> f_ang_kappa{"f_ang_kappa", {1.0, 1.0, 2.0}, "angularity momentum exponent"}; // This is usually 1.0
  Configurable<std::vector<float>> f_ang_alpha{"f_ang_alpha", {1.0, 2.0, 1.0}, "angularity angle exponent"};

  Configurable<bool> b_DoConstSub{"b_DoConstSub", false, "do constituent subtraction"};

  std::vector<fastjet::PseudoJet> jetConstituents;
  std::vector<fastjet::PseudoJet> jetClusterConstituents;
  std::vector<fastjet::PseudoJet> jetReclustered;
  JetFinder jetReclusterer;

  void init(InitContext& initContext)
  {
    jetReclusterer.isReclustering = true;
    jetReclusterer.algorithm = fastjet::JetAlgorithm::cambridge_algorithm;
    jetReclusterer.jetR = f_jetR * 2.5; // Use larger R for CA reclustering to prevent losing particles

    auto& workflows = initContext.services().get<RunningWorkflowInfo const>();
    for (DeviceSpec device : workflows.devices) {
      // std::cout << "Device: " << device.name << std::endl; // Which devices do we have?
      if (device.name == "jet-finder-data") {
        for (auto option : device.options) {
          std::cout << "Option: " << option.name << std::endl; // Which options are available?
          if (option.name == "jetPtMin") {
            // LOGF(info, "Retrieved jet finder vertex z cut: %s.%s", device.name, option.name);
            std::cout << "Jet Pt Min: " << option.defaultValue.get<float>() << std::endl;
          }
          if (option.name == "jetR") {
            std::cout << "Default jet R: " << option.defaultValue << std::endl;
            // f_jetR = option.defaultValue.get<float*>();
            std::cout << "New jet R: " << *(option.defaultValue.get<double*>()) << std::endl;
          }
        }
      }
    }

    Int_t fgkNPtBins = 200;
    Float_t kMinPt   = 0.;
    Float_t kMaxPt   = 200.;

    // Histograms
    spectra.add("hJetZg", "Jet zg", HistType::kTH2F, {
      {fgkNPtBins, kMinPt, kMaxPt, "#it{p}_{T}^{jet} (GeV)"},
      {10, 0., .5, "#it{z}_{g}"}
      });
    spectra.add("hJetRg", "Jet Rg", HistType::kTH2F, {
      {fgkNPtBins, kMinPt, kMaxPt, "#it{p}_{T}^{jet} (GeV)"},
      {20, 0., 1., "#it{R}_{g}"}
      });
    spectra.add("hJetnsd", "Jet nsd", HistType::kTH2F, {
      {fgkNPtBins, kMinPt, kMaxPt, "#it{p}_{T}^{jet} (GeV)"},
      {10, 0., 10., "#it{n}_{SD}"}
      });
    spectra.add("hJetMass", "Jet mass", HistType::kTH2F, {
      {fgkNPtBins, kMinPt, kMaxPt, "#it{p}_{T}^{jet} (GeV)"},
      {50, 0., 50., "#it{m} (GeV)"}
      });
    spectra.add("hJetNEF", "Jet NEF", HistType::kTH2F, {
      {fgkNPtBins, kMinPt, kMaxPt, "#it{p}_{T}^{jet} (GeV)"},
      {50, 0., 1., "NEF"}
      });
    spectra.add("hJetnconst", "Jet nconst", HistType::kTH2F, {
      {fgkNPtBins, kMinPt, kMaxPt, "#it{p}_{T}^{jet} (GeV)"},
      {50, 0., 50., "#it{n}_{const}"}
      });
    spectra.add("hJetPtD", "Jet ptD", HistType::kTH2F, {
      {fgkNPtBins, kMinPt, kMaxPt, "#it{p}_{T}^{jet} (GeV)"},
      {50, 0., 1., "#it{p}_{T}^{D}"}
      });
    spectra.add("hJetScaledMass", "Jet scaled mass", HistType::kTH2F, {
      {fgkNPtBins, kMinPt, kMaxPt, "#it{p}_{T}^{jet} (GeV)"},
      {200, 0., 2., "#it{m}/#it{p}_{T}"}
      });
    spectra.add("hJetChargeFrag", "Jet Charge Fragmentation", HistType::kTH2F, {
      {fgkNPtBins, kMinPt, kMaxPt, "#it{p}_{T}^{jet} (GeV)"},
      {100, 0., 1., "#it{F}^{jet}"}
      });
    spectra.add("hJetAngularities_11", "Jet Angularities 1,1: girth", HistType::kTH2F, {
      {fgkNPtBins, kMinPt, kMaxPt, "#it{p}_{T}^{jet} (GeV)"},
      {20, 0., 1., "#it{z}^{1}#theta^{1}"}
      });
    spectra.add("hJetAngularities_21", "Jet Angularities 2,1", HistType::kTH2F, {
      {fgkNPtBins, kMinPt, kMaxPt, "#it{p}_{T}^{jet} (GeV)"},
      {20, 0., 1., "#it{z}^{2}#theta{1}"}
      });
    spectra.add("hJetAngularities_12", "Jet Angularities 1,2", HistType::kTH2F, {
      {fgkNPtBins, kMinPt, kMaxPt, "#it{p}_{T}^{jet} (GeV)"},
      {20, 0., 1., "#it{z}^{1}#theta{2}"}
      });
    spectra.add("hClusterE", "Cluster E", HistType::kTH2F, {
      {fgkNPtBins, kMinPt, kMaxPt, "#it{p}_{T}^{jet} (GeV)"},
      {fgkNPtBins, kMinPt, kMaxPt, "#it{E}"}
      });
  }

  void SoftDrop(aod::Jet const& jet,
                std::vector<fastjet::PseudoJet> const& jetReclustered,
                float zg = 0,
                float rg = 0,
                int nsd = 0);
  void NEF(aod::Jet const& jet, std::vector<fastjet::PseudoJet> const& jetClusterConstituents);
  void ChargeFragmentation(aod::Jet const& jet,
                           std::vector<fastjet::PseudoJet> const& jetConstituents);
  float Angularity(aod::Jet const& jet,
                   std::vector<fastjet::PseudoJet> const& jetConstituents,
                   std::vector<fastjet::PseudoJet> const& jetClusterConstituents,
                   float alpha = 1.0,
                   float kappa = 1.0);
  void ptD(aod::Jet const& jet,
           std::vector<fastjet::PseudoJet> const& jetConstituents,
           std::vector<fastjet::PseudoJet> const& jetClusterConstituents);

  void process(aod::Jet const& jet,
               aod::Tracks const& tracks,
               aod::JetTrackConstituents const& constituents,
               aod::JetClusterConstituents const& clusterConstituents,
               aod::EMCALClusters const& emcalClusters,
               aod::JetConstituentsSub const& constituentsSub)//,
              //  aod::JetFilters const& jetFilter)
  {
    spectra.fill(HIST("hJetMass"), jet.pt(), jet.mass());
    spectra.fill(HIST("hJetScaledMass"), jet.pt(), jet.mass() / jet.pt());

    jetConstituents.clear();
    jetClusterConstituents.clear();
    jetReclustered.clear();
    if (b_DoConstSub) {
      for (const auto& constituent : constituentsSub) {
        fillConstituents(constituent, jetConstituents);
      }
    } else {
      for (const auto& constituentIndex : constituents) {
        auto constituent = constituentIndex.track();
        fillConstituents(constituent, jetConstituents);
      }
      int count = 0;
      for (const auto& clusterConstituent : clusterConstituents) {
        // This is based on https://github.com/AliceO2Group/O2Physics/blob/90b01104988b5697ac108e51ea4d60429d2e9e40/PWGJE/TableProducer/jetfinder.cxx#L231
        count++;
        auto cluster = clusterConstituent.cluster();
        double pt = cluster.energy() / std::cosh(cluster.eta());
        jetClusterConstituents.emplace_back(
          pt * std::cos(cluster.phi()),
          pt * std::sin(cluster.phi()),
          pt * std::sinh(cluster.eta()),
          cluster.energy()
        );
      }
    }

    /*
    int nJ = 0;
    int nT = 0;

    for (const auto& filter : jetFilters) {
      if (filter.hasHighPtTrack()) {
        // std::cout << "High pt track ";
        ++nT;
        if (filter.hasJetChHighPt()) {
          // std::cout << " has high pt jet." << std::endl;
          ++nJ;
        }
        // else std::cout << " has no high pt jet." << std::endl;
      }
    }
    std::cout << "nJ/nT = " << nJ << "/" << nT << std::endl;
    */

    spectra.fill(HIST("hJetnconst"), jet.pt(), jetConstituents.size() + jetClusterConstituents.size());
    ChargeFragmentation(jet, jetConstituents);
    NEF(jet, jetClusterConstituents);
    ptD(jet, jetConstituents, jetClusterConstituents);

    if (f_ang_alpha->size() != f_ang_kappa->size()){
      std::cout << "Warning: different amount of alpha (" << f_ang_alpha->size() << "), kappa (" << f_ang_kappa->size() << ") values. Truncating to shortest array." << std::endl;
    }
    float ang = -999.;
    ang = Angularity(jet, jetConstituents,
                      jetClusterConstituents,
                      f_ang_alpha->at(0), f_ang_kappa->at(0));
    spectra.fill(HIST("hJetAngularities_11"), jet.pt(), ang);
    ang = Angularity(jet, jetConstituents,
                      jetClusterConstituents,
                      f_ang_alpha->at(1), f_ang_kappa->at(1));
    spectra.fill(HIST("hJetAngularities_12"), jet.pt(), ang);
    ang = Angularity(jet, jetConstituents,
                      jetClusterConstituents,
                      f_ang_alpha->at(2), f_ang_kappa->at(2));
    spectra.fill(HIST("hJetAngularities_21"), jet.pt(), ang);

    std::vector<fastjet::PseudoJet> allJetConstituents;
    allJetConstituents.reserve(jetConstituents.size() + jetClusterConstituents.size());
    allJetConstituents.insert(allJetConstituents.end(), jetConstituents.begin(), jetConstituents.end());
    allJetConstituents.insert(allJetConstituents.end(), jetClusterConstituents.begin(), jetClusterConstituents.end());
    fastjet::ClusterSequenceArea clusterSeq(jetReclusterer.findJets(allJetConstituents, jetReclustered)); // TODO: Include cluster constituents here
    jetReclustered = sorted_by_pt(jetReclustered);
    // bool softDropped = false;
    int nsd = 0.0;
    float zg = -1.0;
    float rg = -1.0;
    SoftDrop(jet, jetReclustered, zg, rg, nsd);

    jetTriggerQA(zg, rg, nsd);
  }
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<JetTriggerQA>(cfgc, TaskName{"jet-trigger-qa"})};
}

#endif
