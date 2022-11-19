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

///
/// \file   qaTrackSelection.cxx
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \since  04/11/2022
/// \brief  Task to check how many tracks pass the cuts
///

// O2 includes
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/TrackSelectionTables.h"

using namespace o2::framework;

struct QaTrackCuts {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
#define fillHistogram(name, ...) histos.fill(HIST(name), __VA_ARGS__)

  static constexpr int nhist = 10;
  static constexpr std::string_view hselection[nhist] = {"NoEvSel/alltracks", "NoEvSel/hastof", "NoEvSel/hastpc", "NoEvSel/hasits", "NoEvSel/hastrd",
                                                         "sel8/alltracks", "sel8/hastof", "sel8/hastpc", "sel8/hasits", "sel8/hastrd"};
  static constexpr std::string_view hisQualityTrack[2] = {"NoEvSel/isQualityTrack",
                                                          "sel8/isQualityTrack"};
  static constexpr std::string_view htrackTypes[2] = {"NoEvSel/tracktypes",
                                                      "sel8/tracktypes"};

  void init(InitContext&)
  {
    const AxisSpec axisSelections{30, 0.5, 30.5f, "Selection"};
    // histos.add("events", "events", kTH1F, {axisSelections});
    for (int i = 0; i < nhist; i++) {
      auto h = histos.add<TH1>(hselection[i].data(), "", kTH1F, {axisSelections});
      h->SetTitle(hselection[i].data());
      h->GetXaxis()->SetBinLabel(1, "Tracks read");
      h->GetXaxis()->SetBinLabel(2, "isGlobalTrackSDD");
      h->GetXaxis()->SetBinLabel(3, "passedTrackType");
      h->GetXaxis()->SetBinLabel(4, "passedPtRange");
      h->GetXaxis()->SetBinLabel(5, "passedEtaRange");
      h->GetXaxis()->SetBinLabel(6, "passedTPCNCls");
      h->GetXaxis()->SetBinLabel(7, "passedTPCCrossedRows");
      h->GetXaxis()->SetBinLabel(8, "passedTPCCrossedRowsOverNCls");
      h->GetXaxis()->SetBinLabel(9, "passedTPCChi2NDF");
      h->GetXaxis()->SetBinLabel(10, "passedTPCRefit");
      h->GetXaxis()->SetBinLabel(11, "passedITSNCls");
      h->GetXaxis()->SetBinLabel(12, "passedITSChi2NDF");
      h->GetXaxis()->SetBinLabel(13, "passedITSRefit");
      h->GetXaxis()->SetBinLabel(14, "passedITSHits");
      h->GetXaxis()->SetBinLabel(15, "passedGoldenChi2");
      h->GetXaxis()->SetBinLabel(16, "passedDCAxy");
      h->GetXaxis()->SetBinLabel(17, "passedDCAz");
      h->GetXaxis()->SetBinLabel(18, "isQualityTrack");
      h->GetXaxis()->SetBinLabel(19, "isPrimaryTrack");
      h->GetXaxis()->SetBinLabel(20, "isInAcceptanceTrack");
      h->GetXaxis()->SetBinLabel(21, "isGlobalTrack");
      h->GetXaxis()->SetBinLabel(22, "isGlobalTrackWoPtEta");
      h->GetXaxis()->SetBinLabel(23, "isGlobalTrackWoDCA");
    }

    for (int i = 0; i < 2; i++) {
      auto h = histos.add<TH1>(hisQualityTrack[i].data(), "Tracks selection for isQualityTrack", kTH1F, {axisSelections});
      h->GetXaxis()->SetBinLabel(1, "Tracks read");
      h->GetXaxis()->SetBinLabel(2, "passedTrackType");
      h->GetXaxis()->SetBinLabel(3, "passedTPCNCls");
      h->GetXaxis()->SetBinLabel(4, "passedTPCCrossedRows");
      h->GetXaxis()->SetBinLabel(5, "passedTPCCrossedRowsOverNCls");
      h->GetXaxis()->SetBinLabel(6, "passedTPCChi2NDF");
      h->GetXaxis()->SetBinLabel(7, "passedTPCRefit");
      h->GetXaxis()->SetBinLabel(8, "passedITSNCls");
      h->GetXaxis()->SetBinLabel(9, "passedITSChi2NDF");
      h->GetXaxis()->SetBinLabel(10, "passedITSRefit");
      h->GetXaxis()->SetBinLabel(11, "passedITSHits");
    }

    for (int i = 0; i < 2; i++) {
      auto h = histos.add<TH1>(htrackTypes[i].data(), "Tracks types seen", kTH1F, {axisSelections});
      h->GetXaxis()->SetBinLabel(1, "TrackIU");
      h->GetXaxis()->SetBinLabel(2, "Track");
      h->GetXaxis()->SetBinLabel(3, "Run2Track");
      h->GetXaxis()->SetBinLabel(4, "Run2Tracklet");
      h->GetXaxis()->SetBinLabel(5, "Undefined");
    }
  }

  template <int index, typename TrackType>
  void fillHistoCuts(const TrackType& track)
  {
    fillHistogram(hselection[index].data(), 1.f);
    if (track.isGlobalTrackSDD()) {
      fillHistogram(hselection[index].data(), 2.f);
    }
    if (track.passedTrackType()) {
      fillHistogram(hselection[index].data(), 3.f);
    }
    if (track.passedPtRange()) {
      fillHistogram(hselection[index].data(), 4.f);
    }
    if (track.passedEtaRange()) {
      fillHistogram(hselection[index].data(), 5.f);
    }
    if (track.passedTPCNCls()) {
      fillHistogram(hselection[index].data(), 6.f);
    }
    if (track.passedTPCCrossedRows()) {
      fillHistogram(hselection[index].data(), 7.f);
    }
    if (track.passedTPCCrossedRowsOverNCls()) {
      fillHistogram(hselection[index].data(), 8.f);
    }
    if (track.passedTPCChi2NDF()) {
      fillHistogram(hselection[index].data(), 9.f);
    }
    if (track.passedTPCRefit()) {
      fillHistogram(hselection[index].data(), 10.f);
    }
    if (track.passedITSNCls()) {
      fillHistogram(hselection[index].data(), 11.f);
    }
    if (track.passedITSChi2NDF()) {
      fillHistogram(hselection[index].data(), 12.f);
    }
    if (track.passedITSRefit()) {
      fillHistogram(hselection[index].data(), 13.f);
    }
    if (track.passedITSHits()) {
      fillHistogram(hselection[index].data(), 14.f);
    }
    if (track.passedGoldenChi2()) {
      fillHistogram(hselection[index].data(), 15.f);
    }
    if (track.passedDCAxy()) {
      fillHistogram(hselection[index].data(), 16.f);
    }
    if (track.passedDCAz()) {
      fillHistogram(hselection[index].data(), 17.f);
    }
    if (track.isQualityTrack()) {
      fillHistogram(hselection[index].data(), 18.f);
    }
    if (track.isPrimaryTrack()) {
      fillHistogram(hselection[index].data(), 19.f);
    }
    if (track.isInAcceptanceTrack()) {
      fillHistogram(hselection[index].data(), 20.f);
    }
    if (track.isGlobalTrack()) {
      fillHistogram(hselection[index].data(), 21.f);
    }
    if (track.isGlobalTrackWoPtEta()) {
      fillHistogram(hselection[index].data(), 22.f);
    }
    if (track.isGlobalTrackWoDCA()) {
      fillHistogram(hselection[index].data(), 23.f);
    }
  }

  template <int index, typename TrackType>
  void fillIsQualityTrack(const TrackType& track)
  {
    fillHistogram(hisQualityTrack[index], 1.f);
    if (track.passedTrackType()) {
      fillHistogram(hisQualityTrack[index], 2.f);
      if (track.passedTPCNCls()) {
        fillHistogram(hisQualityTrack[index], 3.f);
        if (track.passedTPCCrossedRows()) {
          fillHistogram(hisQualityTrack[index], 4.f);
          if (track.passedTPCCrossedRowsOverNCls()) {
            fillHistogram(hisQualityTrack[index], 5.f);
            if (track.passedTPCChi2NDF()) {
              fillHistogram(hisQualityTrack[index], 6.f);
              if (track.passedTPCRefit()) {
                fillHistogram(hisQualityTrack[index], 7.f);
                if (track.passedITSNCls()) {
                  fillHistogram(hisQualityTrack[index], 8.f);
                  if (track.passedITSChi2NDF()) {
                    fillHistogram(hisQualityTrack[index], 9.f);
                    if (track.passedITSRefit()) {
                      fillHistogram(hisQualityTrack[index], 10.f);
                      if (track.passedITSHits()) {
                        fillHistogram(hisQualityTrack[index], 11.f);
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  template <int index, typename TrackType>
  void fillTrackType(const TrackType& track)
  {
    switch (track.trackType()) {
      case o2::aod::track::TrackTypeEnum::TrackIU:
        fillHistogram(htrackTypes[index].data(), 1.f);
        break;
      case o2::aod::track::TrackTypeEnum::Track:
        fillHistogram(htrackTypes[index].data(), 2.f);
        break;
      case o2::aod::track::TrackTypeEnum::Run2Track:
        fillHistogram(htrackTypes[index].data(), 3.f);
        break;
      case o2::aod::track::TrackTypeEnum::Run2Tracklet:
        fillHistogram(htrackTypes[index].data(), 4.f);
        break;
      default:
        fillHistogram(htrackTypes[index].data(), 5.f);
        break;
    }
  }

  using AnalysisTracks = o2::soa::Join<o2::aod::Tracks, o2::aod::TracksExtra, o2::aod::TrackSelection>;
  using AnalysisColls = o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels>;
  void process(const AnalysisTracks& tracks,
               const AnalysisColls&)
  {
    for (const auto& track : tracks) {
      fillHistoCuts<0>(track);
      if (track.hasTOF()) {
        fillHistoCuts<1>(track);
      }
      if (track.hasTPC()) {
        fillHistoCuts<2>(track);
      }
      if (track.hasITS()) {
        fillHistoCuts<3>(track);
      }
      if (track.hasTRD()) {
        fillHistoCuts<4>(track);
      }

      fillIsQualityTrack<0>(track);
      fillTrackType<0>(track);

      if (!track.has_collision()) {
        continue;
      }
      if (!track.collision_as<AnalysisColls>().sel8()) {
        continue;
      }

      fillHistoCuts<5>(track);
      if (track.hasTOF()) {
        fillHistoCuts<6>(track);
      }
      if (track.hasTPC()) {
        fillHistoCuts<7>(track);
      }
      if (track.hasITS()) {
        fillHistoCuts<8>(track);
      }
      if (track.hasTRD()) {
        fillHistoCuts<9>(track);
      }

      fillIsQualityTrack<1>(track);
      fillTrackType<1>(track);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<QaTrackCuts>(cfgc)};
}
