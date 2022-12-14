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
//
//  *+-+*+-+*+-+*+-+*+-+*+-+*
//  Strangeness builder task
//  *+-+*+-+*+-+*+-+*+-+*+-+*
//
//  This task loops over a set of V0 and cascade indices and
//  creates the corresponding analysis tables that contain
//  the typical information required for analysis.
//
//  PERFORMANCE WARNING: this task includes several track
//  propagation calls that are intrinsically heavy. Please
//  also be cautious when adjusting selections: these can
//  increase / decrease CPU consumption quite significantly.
//
//  IDEAL USAGE: if you are interested in taking V0s and
//  cascades and propagating TrackParCovs based on these,
//  please do not re-propagate the daughters. Instead,
//  the tables generated by this builder task can be used
//  to instantiate a TrackPar object (default operation)
//  or even a TrackParCov object (for which you will
//  need to enable the option of producing the V0Cov and
//  CascCov tables too).
//
//    Comments, questions, complaints, suggestions?
//    Please write to:
//    david.dobrigkeit.chinellato@cern.ch
//

#include <cmath>
#include <array>
#include <cstdlib>
#include <map>
#include <iterator>
#include <utility>

#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "DetectorsVertexing/DCAFitterN.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"

#include "TFile.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TLorentzVector.h"
#include "Math/Vector4D.h"
#include "TPDGCode.h"
#include "TDatabasePDG.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

// use parameters + cov mat non-propagated, aux info + (extension propagated)
using FullTracksExt = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TracksDCA>;
using FullTracksExtIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA>;
using LabeledTracks = soa::Join<aod::Tracks, aod::McTrackLabels>;
using V0full = soa::Join<aod::V0Datas, aod::V0Covs>;

//*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
// Builder task: rebuilds multi-strange candidates
struct multistrangeBuilder {
  Produces<aod::CascData> cascdata;
  Produces<aod::CascCovs> casccovs; // if requested by someone
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  Configurable<bool> d_UseAutodetectMode{"d_UseAutodetectMode", true, "Autodetect requested topo sels"};

  // Configurables related to table creation
  Configurable<int> createCascCovMats{"createCascCovMats", -1, {"Produces V0 cov matrices. -1: auto, 0: don't, 1: yes. Default: auto (-1)"}};

  // Topological selection criteria
  Configurable<int> mincrossedrows{"mincrossedrows", 70, "min crossed rows"};
  Configurable<int> tpcrefit{"tpcrefit", 0, "demand TPC refit"};
  Configurable<float> dcabachtopv{"dcabachtopv", .05, "DCA Bach To PV"};
  Configurable<float> cascradius{"cascradius", 0.9, "cascradius"};
  Configurable<float> casccospa{"casccospa", 0.95, "casccospa"};
  Configurable<float> dcacascdau{"dcacascdau", 1.0, "DCA cascade Daughters"};
  Configurable<float> lambdaMassWindow{"lambdaMassWindow", .01, "Distance from Lambda mass"};

  // Operation and minimisation criteria
  Configurable<double> d_bz_input{"d_bz", -999, "bz field, -999 is automatic"};
  Configurable<bool> d_UseAbsDCA{"d_UseAbsDCA", true, "Use Abs DCAs"};
  Configurable<bool> d_UseWeightedPCA{"d_UseWeightedPCA", false, "Vertices use cov matrices"};
  Configurable<int> useMatCorrType{"useMatCorrType", 0, "0: none, 1: TGeo, 2: LUT"};
  Configurable<int> rejDiffCollTracks{"rejDiffCollTracks", 0, "rejDiffCollTracks"};

  // CCDB options
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};

  int mRunNumber;
  float d_bz;
  float maxSnp;  // max sine phi for propagation
  float maxStep; // max step size (cm) for propagation
  o2::base::MatLayerCylSet* lut = nullptr;

  // Define o2 fitter, 2-prong, active memory (no need to redefine per event)
  o2::vertexing::DCAFitterN<2> fitter;
  enum cascstep { kCascAll = 0,
                  kCascLambdaMass,
                  kBachTPCrefit,
                  kBachCrossedRows,
                  kBachDCAxy,
                  kCascDCADau,
                  kCascCosPA,
                  kCascRadius,
                  kNCascSteps };

  // Helper struct to pass cascade information
  struct {
    int v0Id;
    int bachelorId;
    int collisionId;
    int charge;
    std::array<float, 3> pos;
    std::array<float, 3> bachP;
    float dcacascdau;
    float bachDCAxy;
    float cosPA;
    float cascradius;
  } cascadecandidate;

  o2::track::TrackParCov lBachelorTrack;
  o2::track::TrackParCov lV0Track;

  // Helper struct to do bookkeeping of building parameters
  struct {
    std::array<long, kNCascSteps> cascstats;
    long exceptions;
    long eventCounter;
  } statisticsRegistry;

  HistogramRegistry registry{
    "registry",
    {{"hEventCounter", "hEventCounter", {HistType::kTH1F, {{1, 0.0f, 1.0f}}}},
     {"hCaughtExceptions", "hCaughtExceptions", {HistType::kTH1F, {{1, 0.0f, 1.0f}}}},
     {"hCascadeCriteria", "hCascadeCriteria", {HistType::kTH1F, {{10, -0.5f, 9.5f}}}}}};

  void resetHistos()
  {
    statisticsRegistry.exceptions = 0;
    statisticsRegistry.eventCounter = 0;
    for (Int_t ii = 0; ii < kNCascSteps; ii++)
      statisticsRegistry.cascstats[ii] = 0;
  }

  void fillHistos()
  {
    registry.fill(HIST("hEventCounter"), 0.0, statisticsRegistry.eventCounter);
    registry.fill(HIST("hCaughtExceptions"), 0.0, statisticsRegistry.exceptions);
    for (Int_t ii = 0; ii < kNCascSteps; ii++)
      registry.fill(HIST("hCascadeCriteria"), ii, statisticsRegistry.cascstats[ii]);
  }

  void init(InitContext& context)
  {
    resetHistos();

    mRunNumber = 0;
    d_bz = 0;
    maxSnp = 0.85f;  // could be changed later
    maxStep = 2.00f; // could be changed later

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(lutPath));
    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      ccdb->get<TGeoManager>(geoPath);
    }

    if (doprocessRun2 == false && doprocessRun3 == false) {
      LOGF(fatal, "Neither processRun2 nor processRun3 enabled. Please choose one.");
    }
    if (doprocessRun2 == true && doprocessRun3 == true) {
      LOGF(fatal, "Cannot enable processRun2 and processRun3 at the same time. Please choose one.");
    }

    if (d_UseAutodetectMode) {
      // Checking for subscriptions to:
      double loosest_casccospa = 100;
      float loosest_dcacascdau = -100;
      float loosest_dcabachtopv = 100;
      float loosest_dcav0topv = 100;
      float loosest_radius = 100;
      float loosest_v0masswindow = -100;

      double detected_casccospa = 100;
      float detected_dcacascdau = -100;
      float detected_dcabachtopv = 100;
      float detected_dcav0topv = 100;
      float detected_radius = 100;
      float detected_v0masswindow = -100;

      LOGF(info, "*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*");
      LOGF(info, " Multi-strange builder self-configuration");
      LOGF(info, "*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*");
      auto& workflows = context.services().get<RunningWorkflowInfo const>();
      for (DeviceSpec const& device : workflows.devices) {
        // Step 1: check if this device subscribed to the V0data table
        for (auto const& input : device.inputs) {
          if (device.name.compare("cascade-initializer") == 0)
            continue; // don't listen to the initializer, it's just to extend stuff
          const std::string CascDataName = "CascData";
          const std::string CascDataExtName = "CascDataExt";
          if (input.matcher.binding == CascDataName || input.matcher.binding == CascDataExtName) {
            LOGF(info, "Device named %s has subscribed to CascData table! Will now scan for desired settings...", device.name);
            for (auto const& option : device.options) {

              // 5 V0 topological selections + 1 mass
              if (option.name.compare("cascadesetting_cospa") == 0) {
                detected_casccospa = option.defaultValue.get<double>();
                LOGF(info, "%s requested cascade cospa = %f", device.name, detected_casccospa);
                if (detected_casccospa < loosest_casccospa)
                  loosest_casccospa = detected_casccospa;
              }
              if (option.name.compare("cascadesetting_dcacascdau") == 0) {
                detected_dcacascdau = option.defaultValue.get<float>();
                LOGF(info, "%s requested DCA cascade daughters = %f", device.name, detected_dcacascdau);
                if (detected_dcacascdau > loosest_dcacascdau)
                  loosest_dcacascdau = detected_dcacascdau;
              }
              if (option.name.compare("cascadesetting_dcabachtopv") == 0) {
                detected_dcabachtopv = option.defaultValue.get<float>();
                LOGF(info, "%s requested DCA bachelor daughter = %f", device.name, detected_dcabachtopv);
                if (detected_dcabachtopv < loosest_dcabachtopv)
                  loosest_dcabachtopv = detected_dcabachtopv;
              }
              if (option.name.compare("cascadesetting_cascradius") == 0) {
                detected_radius = option.defaultValue.get<float>();
                LOGF(info, "%s requested  to PV = %f", device.name, detected_radius);
                if (detected_radius < loosest_radius)
                  loosest_radius = detected_radius;
              }
              if (option.name.compare("cascadesetting_mindcav0topv") == 0) {
                detected_dcav0topv = option.defaultValue.get<float>();
                LOGF(info, "%s requested minimum V0 DCA to PV = %f", device.name, detected_dcav0topv);
                if (detected_dcav0topv < loosest_dcav0topv)
                  loosest_dcav0topv = detected_dcav0topv;
              }
              if (option.name.compare("cascadesetting_v0masswindow") == 0) {
                detected_v0masswindow = option.defaultValue.get<float>();
                LOGF(info, "%s requested minimum V0 mass window (GeV/c^2) = %f", device.name, detected_v0masswindow);
                if (detected_v0masswindow > loosest_v0masswindow)
                  loosest_v0masswindow = detected_v0masswindow;
              }
            }
          }
          const std::string CascCovsName = "CascCovs";
          if (input.matcher.binding == CascCovsName) {
            LOGF(info, "Device named %s has subscribed to CascCovs table! Enabling.", device.name);
            createCascCovMats.value = 1;
          }
        }
      }

      LOGF(info, "Self-configuration finished! Decided on selections:");
      LOGF(info, " -+*> Cascade cospa ............: %.6f", loosest_casccospa);
      LOGF(info, " -+*> DCA cascade daughters ....: %.6f", loosest_dcacascdau);
      LOGF(info, " -+*> DCA bachelor daughter ....: %.6f", loosest_dcabachtopv);
      LOGF(info, " -+*> Min DCA V0 to PV .........: %.6f", loosest_dcav0topv);
      LOGF(info, " -+*> Min cascade decay radius .: %.6f", loosest_radius);
      LOGF(info, " -+*> V0 mass window ...........: %.6f", loosest_v0masswindow);

      casccospa.value = loosest_casccospa;
      dcacascdau.value = loosest_dcacascdau;
      dcabachtopv.value = loosest_dcabachtopv;
      // dcav0dau.value = loosest_dcav0topv;
      cascradius.value = loosest_radius;
      lambdaMassWindow.value = loosest_v0masswindow;
    }

    //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
    LOGF(info, "Strangeness builder configuration:");
    if (doprocessRun2 == true) {
      LOGF(info, "Run 2 processing enabled. Will subscribe to Tracks table.");
    };
    if (doprocessRun3 == true) {
      LOGF(info, "Run 3 processing enabled. Will subscribe to TracksIU table.");
    };
    if (createCascCovMats > 0) {
      LOGF(info, "-> Will produce cascade cov mat table");
    };
    //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*

    // initialize O2 2-prong fitter (only once)
    fitter.setPropagateToPCA(true);
    fitter.setMaxR(200.);
    fitter.setMinParamChange(1e-3);
    fitter.setMinRelChi2Change(0.9);
    fitter.setMaxDZIni(1e9);
    fitter.setMaxChi2(1e9);
    fitter.setUseAbsDCA(d_UseAbsDCA);
    fitter.setWeightedFinalPCA(d_UseWeightedPCA);

    // Material correction in the DCA fitter
    o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
    if (useMatCorrType == 1)
      matCorr = o2::base::Propagator::MatCorrType::USEMatCorrTGeo;
    if (useMatCorrType == 2)
      matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
    fitter.setMatCorrType(matCorr);
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    auto run3grp_timestamp = bc.timestamp();

    o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grp_timestamp);
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (grpo) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      if (d_bz_input < -990) {
        // Fetch magnetic field from ccdb for current collision
        d_bz = grpo->getNominalL3Field();
        LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
      } else {
        d_bz = d_bz_input;
      }
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3grp_timestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << run3grp_timestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
      if (d_bz_input < -990) {
        // Fetch magnetic field from ccdb for current collision
        d_bz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
        LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
      } else {
        d_bz = d_bz_input;
      }
    }
    o2::base::Propagator::Instance()->setMatLUT(lut);
    mRunNumber = bc.runNumber();
    // Set magnetic field value once known
    fitter.setBz(d_bz);
  }

  template <class TTracksTo, typename TV0Object>
  bool buildCascadeCandidate(aod::Collision const& collision, TTracksTo const& bachTrack, TV0Object const& v0)
  {
    // value 0.5: any considered cascade
    statisticsRegistry.cascstats[kCascAll]++;

    // Overall cascade charge
    cascadecandidate.charge = bachTrack.signed1Pt() > 0 ? +1 : -1;
    cascadecandidate.bachDCAxy = bachTrack.dcaXY();

    // check also against charge
    if (cascadecandidate.charge < 0 && TMath::Abs(v0.mLambda() - 1.116) > lambdaMassWindow)
      return false;
    if (cascadecandidate.charge > 0 && TMath::Abs(v0.mAntiLambda() - 1.116) > lambdaMassWindow)
      return false;
    statisticsRegistry.cascstats[kCascLambdaMass]++;

    if (tpcrefit) {
      if (!(bachTrack.trackType() & o2::aod::track::TPCrefit)) {
        return false;
      }
    }
    statisticsRegistry.cascstats[kBachTPCrefit]++;
    if (bachTrack.tpcNClsCrossedRows() < mincrossedrows) {
      return false;
    }
    statisticsRegistry.cascstats[kBachCrossedRows]++;

    // bachelor DCA track to PV
    if (TMath::Abs(bachTrack.dcaXY()) < dcabachtopv)
      return false;
    statisticsRegistry.cascstats[kBachDCAxy]++;

    // Do actual minimization
    lBachelorTrack = getTrackParCov(bachTrack);

    // Set up covariance matrices (should in fact be optional)
    std::array<float, 21> covV = {0.};
    constexpr int MomInd[6] = {9, 13, 14, 18, 19, 20}; // cov matrix elements for momentum component
    for (int i = 0; i < 6; i++) {
      covV[MomInd[i]] = v0.momentumCovMat()[i];
      covV[i] = v0.positionCovMat()[i];
    }
    lV0Track = o2::track::TrackParCov(
      {v0.x(), v0.y(), v0.z()},
      {v0.pxpos() + v0.pxneg(), v0.pypos() + v0.pyneg(), v0.pzpos() + v0.pzneg()},
      covV, cascadecandidate.charge, true);

    //---/---/---/
    // Move close to minima
    int nCand = 0;
    try {
      nCand = fitter.process(lV0Track, lBachelorTrack);
    } catch (...) {
      registry.fill(HIST("hCaughtExceptions"), 0.5f);
      LOG(error) << "Exception caught in DCA fitter process call!";
      return false;
    }
    if (nCand == 0)
      return false;

    lV0Track = fitter.getTrack(0);
    lBachelorTrack = fitter.getTrack(1);

    // DCA between cascade daughters
    cascadecandidate.dcacascdau = TMath::Sqrt(fitter.getChi2AtPCACandidate());
    if (cascadecandidate.dcacascdau > dcacascdau)
      return false;
    statisticsRegistry.cascstats[kCascDCADau]++;

    fitter.getTrack(1).getPxPyPzGlo(cascadecandidate.bachP);
    // get decay vertex coordinates
    const auto& vtx = fitter.getPCACandidate();
    for (int i = 0; i < 3; i++) {
      cascadecandidate.pos[i] = vtx[i];
    }

    cascadecandidate.cosPA = RecoDecay::cpa(
      array{collision.posX(), collision.posY(), collision.posZ()},
      array{cascadecandidate.pos[0], cascadecandidate.pos[1], cascadecandidate.pos[2]},
      array{v0.pxpos() + v0.pxneg() + cascadecandidate.bachP[0], v0.pypos() + v0.pyneg() + cascadecandidate.bachP[1], v0.pzpos() + v0.pzneg() + cascadecandidate.bachP[2]});
    if (cascadecandidate.cosPA < casccospa) {
      return false;
    }
    statisticsRegistry.cascstats[kCascCosPA]++;

    // Cascade radius
    cascadecandidate.cascradius = RecoDecay::sqrtSumOfSquares(cascadecandidate.pos[0], cascadecandidate.pos[1]);
    if (cascadecandidate.cascradius < cascradius)
      return false;
    statisticsRegistry.cascstats[kCascRadius]++;

    return true;
  }

  template <class TTracksTo>
  void buildStrangenessTables(aod::Collision const& collision, aod::Cascades const& cascades, TTracksTo const& tracks)
  {
    statisticsRegistry.eventCounter++;

    for (auto& cascade : cascades) {
      // Track casting
      auto bachTrackCast = cascade.template bachelor_as<TTracksTo>();
      auto v0index = cascade.v0_as<o2::aod::V0sLinked>();
      if (!(v0index.has_v0Data())) {
        // cascdataLink(-1);
        continue; // skip those cascades for which V0 doesn't exist
      }
      auto v0 = v0index.v0Data_as<V0full>(); // de-reference index to correct v0data in case it exists
      //
      bool validCascadeCandidate = buildCascadeCandidate(collision, bachTrackCast, v0);
      if (!validCascadeCandidate)
        continue; // doesn't pass cascade selections

      cascdata(v0index.globalIndex(),
               bachTrackCast.globalIndex(),
               cascade.collisionId(),
               cascadecandidate.charge,
               cascadecandidate.pos[0], cascadecandidate.pos[1], cascadecandidate.pos[2],
               v0.x(), v0.y(), v0.z(),
               v0.pxpos(), v0.pypos(), v0.pzpos(),
               v0.pxneg(), v0.pyneg(), v0.pzneg(),
               cascadecandidate.bachP[0], cascadecandidate.bachP[1], cascadecandidate.bachP[2],
               v0.dcaV0daughters(), cascadecandidate.dcacascdau,
               v0.dcapostopv(), v0.dcanegtopv(),
               cascadecandidate.bachDCAxy);
    }
    // En masse filling at end of process call
    fillHistos();
    resetHistos();
  }

  void processRun2(aod::Collision const& collision, aod::V0sLinked const& V0s, V0full const&, aod::Cascades const& cascades, FullTracksExt const& tracks, aod::BCsWithTimestamps const&)
  {
    /* check the previous run number */
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);

    // do v0s, typecase correctly into tracks (Run 2 use case)
    buildStrangenessTables<FullTracksExt>(collision, cascades, tracks);
  }
  PROCESS_SWITCH(multistrangeBuilder, processRun2, "Produce Run 2 V0 tables", true);

  void processRun3(aod::Collision const& collision, aod::V0sLinked const& V0s, V0full const&, aod::Cascades const& cascades, FullTracksExtIU const& tracks, aod::BCsWithTimestamps const&)
  {
    /* check the previous run number */
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);

    // do v0s, typecase correctly into tracksIU (Run 3 use case)
    buildStrangenessTables<FullTracksExtIU>(collision, cascades, tracks);
  }
  PROCESS_SWITCH(multistrangeBuilder, processRun3, "Produce Run 3 V0 tables", false);
};

//*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
struct multistrangeLabelBuilder {
  Produces<aod::McCascLabels> casclabels; // MC labels for cascades
  // for bookkeeping purposes: how many V0s come from same mother etc

  void init(InitContext const&) {}

  void processDoNotBuildLabels(aod::Collisions::iterator const& collision)
  {
    // dummy process function - should not be required in the future
  }
  PROCESS_SWITCH(multistrangeLabelBuilder, processDoNotBuildLabels, "Do not produce MC label tables", true);

  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  // build cascade labels if requested to do so
  void processBuildCascadeLabels(aod::Collision const& collision, aod::CascDataExt const& casctable, aod::V0sLinked const&, aod::V0Datas const& v0table, LabeledTracks const&, aod::McParticles const&)
  {
    for (auto& casc : casctable) {
      // Loop over those that actually have the corresponding V0 associated to them
      auto v0 = casc.v0_as<o2::aod::V0sLinked>();
      if (!(v0.has_v0Data())) {
        casclabels(-1);
        continue; // skip those cascades for which V0 doesn't exist
      }
      auto v0data = v0.v0Data(); // de-reference index to correct v0data in case it exists
      int lLabel = -1;

      // Acquire all three daughter tracks, please
      auto lBachTrack = casc.bachelor_as<LabeledTracks>();
      auto lNegTrack = v0data.negTrack_as<LabeledTracks>();
      auto lPosTrack = v0data.posTrack_as<LabeledTracks>();

      // Association check
      // There might be smarter ways of doing this in the future
      if (lNegTrack.has_mcParticle() && lPosTrack.has_mcParticle() && lBachTrack.has_mcParticle()) {
        auto lMCBachTrack = lBachTrack.mcParticle_as<aod::McParticles>();
        auto lMCNegTrack = lNegTrack.mcParticle_as<aod::McParticles>();
        auto lMCPosTrack = lPosTrack.mcParticle_as<aod::McParticles>();

        // Step 1: check if the mother is the same, go up a level
        if (lMCNegTrack.has_mothers() && lMCPosTrack.has_mothers()) {
          for (auto& lNegMother : lMCNegTrack.mothers_as<aod::McParticles>()) {
            for (auto& lPosMother : lMCPosTrack.mothers_as<aod::McParticles>()) {
              if (lNegMother == lPosMother) {
                // if we got to this level, it means the mother particle exists and is the same
                // now we have to go one level up and compare to the bachelor mother too
                for (auto& lV0Mother : lNegMother.mothers_as<aod::McParticles>()) {
                  for (auto& lBachMother : lMCBachTrack.mothers_as<aod::McParticles>()) {
                    if (lV0Mother == lBachMother) {
                      lLabel = lV0Mother.globalIndex();
                    }
                  }
                } // end conditional V0-bach pair
              }   // end neg = pos mother conditional
            }
          } // end loop neg/pos mothers
        }   // end conditional of mothers existing
      }     // end association check
      // Construct label table (note: this will be joinable with CascDatas)
      casclabels(
        lLabel);
    } // end casctable loop
  }
  PROCESS_SWITCH(multistrangeLabelBuilder, processBuildCascadeLabels, "Produce cascade MC label tables", false);
  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
};

/// Extends the cascdata table with expression columns
struct cascadeInitializer {
  Spawns<aod::CascDataExt> cascdataext;
  void init(InitContext const&) {}
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<multistrangeBuilder>(cfgc),
    adaptAnalysisTask<multistrangeLabelBuilder>(cfgc),
    adaptAnalysisTask<cascadeInitializer>(cfgc)};
}
