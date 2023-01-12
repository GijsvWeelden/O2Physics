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

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/CCDB/TriggerAliases.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGJE/Core/JetFinder.h"

// FIXME: Make these paths non-global!
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

  OutputObj<TH1I> hProcessedEvents{"hProcessedEvents"};
  OutputObj<TH2F> hClusterComparison{"hClusterComparison"};

  OutputObj<TH1F> hCellAmp{"hCellAmp"};
  OutputObj<TH1F> hCellMaxAmp{"hCellMaxAmp"};

  // MB histograms
  OutputObj<TH1I> hNMBEvents{"hNMBEvents"};
  OutputObj<TH1I> hNMBJets{"hNMBJets"};
  OutputObj<TH1I> hNMBClusters{"hNMBClusters"};
  OutputObj<TH1F> hMBJetPt{"hMBJetPt"};
  OutputObj<TH2F> hMBJetEtaPhi{"hMBJetEtaPhi"};
  OutputObj<TH1F> hMBClusterPt{"hMBClusterPt"};
  OutputObj<TH2F> hMBClusterEtaPhi{"hMBClusterEtaPhi"};
  OutputObj<TH2F> hMBMaxSpectrum{"hMBMaxSpectrum"};

  // EMCAL histograms
  OutputObj<TH1I> hNEmcEvents{"hNEmcEvents"};
  OutputObj<TH1I> hNEmcJets{"hNEmcJets"};
  OutputObj<TH1I> hNEmcClusters{"hNEmcClusters"};
  OutputObj<TH1F> hEmcJetPt{"hEmcJetPt"};
  OutputObj<TH2F> hEmcJetEtaPhi{"hEmcJetEtaPhi"};
  OutputObj<TH1F> hEmcClusterPt{"hEmcClusterPt"};
  OutputObj<TH2F> hEmcClusterEtaPhi{"hEmcClusterEtaPhi"};
  OutputObj<TH2F> hEmcMaxSpectrum{"hEmcMaxSpectrum"};

  // Selected histograms
  OutputObj<TH1I> hNSelectedEvents{"hNSelectedEvents"};
  OutputObj<TH1I> hNSelectedJets{"hNSelectedJets"};
  OutputObj<TH1I> hNSelectedClusters{"hNSelectedClusters"};
  OutputObj<TH1F> hSelectedJetPt{"hSelectedJetPt"};
  OutputObj<TH2F> hSelectedJetEtaPhi{"hSelectedJetEtaPhi"};
  OutputObj<TH1F> hSelectedClusterPt{"hSelectedClusterPt"};
  OutputObj<TH2F> hSelectedClusterEtaPhi{"hSelectedClusterEtaPhi"};
  OutputObj<TH2F> hSelectedMaxSpectrum{"hSelectedMaxSpectrum"};

  // Configurables
  Configurable<float> f_jetPtMin{"f_jetPtMin", 0.0, "minimum jet pT cut"};
  Configurable<float> cfgVertexCut{"cfgVertexCut", 10.0f, "Accepted z-vertex range"};
  Configurable<double> minCellAmplitude{"minCellAmplitude", 0.1, "Minimum cell amplitude for histograms."};
  Configurable<std::string> mClusterDefinition{"clusterDefinition", "kV3Default", "cluster definition to be selected, e.g. V3Default"};
  // Configurable<std::string> geometry{"geometry", "EMCAL", "geometry to be selected on, e.g. EMCAL"};

  void init(o2::framework::InitContext&)
  {
    // Bins here
    Int_t nPtBins = 200;
    Float_t kMinPt   = 0.;
    Float_t kMaxPt   = 200.;

    Int_t nPhiBins = 18*8;
    Float_t kMinPhi   = 0.;
    Float_t kMaxPhi   = 2.*TMath::Pi();

    Int_t nEtaBins = 100;
    Float_t kMinEta = -1.;
    Float_t kMaxEta =  1.;

    // Misc. histograms
    hCellAmp.setObject(
      new TH1F("hCellAmp",
               "Emc cell amp;amp",
               10*nPtBins, kMinPt, kMaxPt
               ));
    hCellMaxAmp.setObject(
      new TH1F("hCellMaxAmp",
               "Highest amp emc cell;amp",
               10*nPtBins, kMinPt, kMaxPt
               ));
    hClusterComparison.setObject(
      new TH2F("hClusterComparison",
               "MB jet cluster/event cluster spectra;#it{p}_{T}^{jet cl} (GeV);#it{p}_{T}^{evt cl} (GeV)",
               nPtBins, kMinPt, kMaxPt/2, nPtBins, kMinPt, kMaxPt/2
               ));

    // MB histograms
    hNMBEvents.setObject( new TH1I("hNMBEvents", "Number of events", 1, 0, 1));
    hNMBJets.setObject( new TH1I("hNMBJets", "Number of jets per event", 100, 0, 100));
    hNMBClusters.setObject( new TH1I("hNMBClusters", "Number of clusters per event", 100, 0, 100));
    hMBJetPt.setObject(
      new TH1F("hMBJetPt",
               "All jet #it{p}_{T};#it{p}_{T}",
               nPtBins, kMinPt, kMaxPt
               ));
    hMBJetEtaPhi.setObject(
      new TH2F("hMBJetEtaPhi",
               "MB jet eta phi;#eta;#phi",
               nEtaBins, kMinEta, kMaxEta, nPhiBins, kMinPhi, kMaxPhi
               ));
    hMBClusterPt.setObject(
      new TH1F("hMBClusterPt",
               "All cluster #it{E}_{T};#it{E}_{T}",
               nPtBins, kMinPt, kMaxPt/2
               ));
    hMBClusterEtaPhi.setObject(
      new TH2F("hMBClusterEtaPhi",
               "MB cluster eta phi;#eta;#phi",
               nEtaBins, kMinEta, kMaxEta, nPhiBins, kMinPhi, kMaxPhi
               ));
    hMBMaxSpectrum.setObject(
      new TH2F("hMBMaxSpectrum",
               "MB jet/cluster spectra;#it{p}_{T}^{jet} (GeV);#it{p}_{T}^{cl} (GeV)",
               nPtBins, kMinPt, kMaxPt, nPtBins, kMinPt, kMaxPt/2
               ));

    // EMCAL histograms
    hNEmcEvents.setObject( new TH1I("hNEmcEvents", "Number of events where EMCAL is live", 1, 0, 1));
    hNEmcJets.setObject( new TH1I("hNEmcJets", "Number of jets per emc event", 100, 0, 100));
    hNEmcClusters.setObject( new TH1I("hNEmcClusters", "Number of clusters per emc event", 100, 0, 100));
    hEmcJetPt.setObject(
      new TH1F("hEmcJetPt",
               "Emc jet #it{p}_{T};#it{p}_{T}",
               nPtBins, kMinPt, kMaxPt
               ));
    hEmcJetEtaPhi.setObject(
      new TH2F("hEmcJetEtaPhi",
               "Emc jet eta phi;#eta;#phi",
               nEtaBins, kMinEta, kMaxEta, nPhiBins, kMinPhi, kMaxPhi
               ));
    hEmcClusterPt.setObject(
      new TH1F("hEmcClusterPt",
               "AlEmcl cluster #it{E}_{T};#it{E}_{T}",
               nPtBins, kMinPt, kMaxPt/2
               ));
    hEmcClusterEtaPhi.setObject(
      new TH2F("hEmcClusterEtaPhi",
               "Emc cluster eta phi;#eta;#phi",
               nEtaBins, kMinEta, kMaxEta, nPhiBins, kMinPhi, kMaxPhi
               ));
    hEmcMaxSpectrum.setObject(
      new TH2F("hEmcMaxSpectrum",
               "Emc jet/cluster spectra;#it{p}_{T}^{jet} (GeV);#it{p}_{T}^{cl} (GeV)",
               nPtBins, kMinPt, kMaxPt, nPtBins, kMinPt, kMaxPt/2
               ));

    // Selected histograms
    hNSelectedEvents.setObject( new TH1I("hNSelectedEvents", "Number of selected events", 1, 0, 1));
    hSelectedJetPt.setObject(
      new TH1F("hSelectedJetPt",
               "Selected jet #it{p}_{T};#it{p}_{T}",
               nPtBins, kMinPt, kMaxPt
               ));
    hSelectedJetEtaPhi.setObject(
      new TH2F("hSelectedJetEtaPhi",
               "Selected jet eta phi;#eta;#phi",
               nEtaBins, kMinEta, kMaxEta, nPhiBins, kMinPhi, kMaxPhi
               ));
    hSelectedClusterPt.setObject(
      new TH1F("hSelectedClusterPt",
               "Selected cluster #it{E}_{T};#it{E}_{T}",
               nPtBins, kMinPt, kMaxPt/2
               ));
    hSelectedClusterEtaPhi.setObject(
      new TH2F("hSelectedClusterEtaPhi",
               "MB cluster eta phi;#eta;#phi",
               nEtaBins, kMinEta, kMaxEta, nPhiBins, kMinPhi, kMaxPhi
               ));
    hSelectedMaxSpectrum.setObject(
      new TH2F("hSelectedMaxSpectrum",
               "Selected jet/cluster spectra;#it{p}_{T}^{jet} (GeV);#it{p}_{T}^{cl} (GeV)",
               nPtBins, kMinPt, kMaxPt, nPtBins, kMinPt, kMaxPt/2
               ));
  } // init()

  // Declare filters
  o2::aod::EMCALClusterDefinition clusDef = o2::aod::emcalcluster::getClusterDefinitionFromString(mClusterDefinition.value);
  Filter clusterDefinitionSelection = o2::aod::emcalcluster::definition == static_cast<int>(clusDef);

  Bool_t isJetInEmcal(aod::Jet const& jet){
    double emcalEtaMin = -0.7, emcalEtaMax = 0.7, emcalPhiMin = 1.40, emcalPhiMax = 3.26; // Phi: 80 - 187 deg
    double R = jet.r() * 1e-2; // Jet R is saved as round(100*R)
    if ( (jet.eta() >= emcalEtaMin + R) && (jet.eta() <= emcalEtaMax - R)
          && (jet.phi() >= emcalPhiMin + R) && (jet.phi() <= emcalPhiMax - R) ){
      return true;
    }
    return false;
  }

  // Filter collisionFilter = nabs(aod::collision::posZ) < cfgVertexCut;

  // In the backend, the process function matches the tables via their shared identifiers. Here: the bcID
  // This means we only get events that also have cells
  // We can circumvent this by using two process functions: one that knows about the cells and one that doesn't
  void process(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision,
               aod::Jets const& jets,
               selectedClusters const& clusters
               )
  {
    bool keepEvent[kCategories]{false};

    double L0emcal = 2.5, L1emcal = 19; // EMCAL hardware JE trigger thresholds
    double maxJetPt = -1., maxClusterPt = -1., maxSelectedJetPt = -1.;
    // maxJetClusterPt = -1., maxCellAmp = -1.;
    int nJets = jets.size(), nClusters = clusters.size();

    hNMBEvents->Fill(0.5); // All events, not just ones with jets and clusters
    hNMBJets->Fill(nJets);
    hNMBClusters->Fill(nClusters);

    // /*
    for (const auto& jet : jets){
      hMBJetPt->Fill(jet.pt());
      hMBJetEtaPhi->Fill(jet.eta(), jet.phi());
      if (jet.pt() > maxJetPt){
        maxJetPt = jet.pt();
      }
      if (nClusters > 0){
        hEmcJetPt->Fill(jet.pt());
        hEmcJetEtaPhi->Fill(jet.eta(), jet.phi());
        if (isJetInEmcal(jet)){
          if (jet.pt() > maxSelectedJetPt) maxSelectedJetPt = jet.pt();
        }
      }
    }
    // */
    // /*
    for (const auto& cluster : clusters){
      double clusterPt = cluster.energy() / std::cosh(cluster.eta());
      hMBClusterPt->Fill(clusterPt);
      hMBClusterEtaPhi->Fill(cluster.eta(), cluster.phi());
      if (clusterPt > maxClusterPt){
        maxClusterPt = clusterPt;
      }
      hEmcClusterPt->Fill(clusterPt);
      hEmcClusterEtaPhi->Fill(cluster.eta(), cluster.phi());
    }
    // */

    // /*
    hMBMaxSpectrum->Fill(maxJetPt, maxClusterPt);
    // double emcalEtaMin = -0.7, emcalEtaMax = 0.7, emcalPhiMin = 1.40, emcalPhiMax = 3.26; // Phi: 80 - 187 deg
    if (nClusters > 0){
      hNEmcEvents->Fill(0.5);
      hEmcMaxSpectrum->Fill(maxJetPt, maxClusterPt);
      if (maxSelectedJetPt >= 0 && maxClusterPt >= 0){
        hNSelectedEvents->Fill(0.5);
        hSelectedMaxSpectrum->Fill(maxSelectedJetPt, maxClusterPt);
        for (const auto& jet : jets){
          if (isJetInEmcal(jet)){
            hSelectedJetPt->Fill(jet.pt());
            hSelectedJetEtaPhi->Fill(jet.eta(), jet.phi());
          }
        }
        for (const auto& cluster : clusters){
          double clusterPt = cluster.energy() / std::cosh(cluster.eta());
          hSelectedClusterPt->Fill(clusterPt);
          hSelectedClusterEtaPhi->Fill(cluster.eta(), cluster.phi());
        }
      }
    }
    // */

    /*
    // For now, use events with a cluster as nEmcEvents
    if (nClusters > 0){
      hNEmcEvents->Fill(0.5);
      hEmcMaxSpectrum->Fill(maxJetPt, maxClusterPt);
      if (maxSelectedJetPt >= 0 && maxSelectedClusterPt >= 0){
        hNSelectedEvents->Fill(0.5);
        hSelectedMaxSpectrum->Fill(maxSelectedJetPt, maxSelectedClusterPt);
      }
      // TODO: enforce EMCAL geometry on jets here?
      for (const auto& jet : jets){
        hEmcJetPt->Fill(jet.pt());
        hEmcJetEtaPhi->Fill(jet.eta(), jet.phi());
        if (maxJetPt >= 0 && maxClusterPt >= 0){
          if (isJetInEmcal(jet)){
            hSelectedJetPt->Fill(jet.pt());
            hSelectedJetEtaPhi->Fill(jet.eta(), jet.phi());
          }
        }
      }
      for (const auto& cluster : clusters){
        double clusterPt = cluster.energy() / std::cosh(cluster.eta());
        hEmcClusterPt->Fill(clusterPt);
        hEmcClusterEtaPhi->Fill(cluster.eta(), cluster.phi());
        if (maxJetPt >= 0 && maxClusterPt >= 0){
          hSelectedClusterPt->Fill(clusterPt);
          hSelectedClusterEtaPhi->Fill(cluster.eta(), cluster.phi());
        }
      }
    }
    // */

    // /*
    for (int iDecision{0}; iDecision < kCategories; iDecision++){
      if (keepEvent[iDecision]) {
        hProcessedEvents->Fill(iDecision);
      }
    }
    // */
    // tags(keepEvent[0], keepEvent[1]);
  } // process()
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{
    adaptAnalysisTask<fullJetFilter>(cfg)};
}
