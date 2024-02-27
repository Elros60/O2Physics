#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include <cmath>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>
#include <ctime>
#include <iostream>
#include <gsl/span>

#include <Math/Vector4D.h>
#include <TCanvas.h>
#include <TDatabasePDG.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TMatrixD.h>
#include <TParameter.h>
#include <TSystem.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TChain.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLine.h>
#include <TSystem.h>
#include <TStopwatch.h>

#include "CommonConstants/LHCConstants.h"
#include "CommonUtils/NameConf.h"
#include "DataFormatsMCH/Cluster.h"
#include "DataFormatsMCH/Digit.h"
#include "DataFormatsMCH/ROFRecord.h"
#include "DataFormatsMCH/TrackMCH.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/GRPGeomHelper.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/Logger.h"
#include "Framework/CallbackService.h"
#include "MCHGeometryTransformer/Transformations.h"
#include "MCHMappingInterface/Segmentation.h"
#include "MCHTracking/TrackExtrap.h"
#include "MCHTracking/TrackParam.h"
#include "ReconstructionDataFormats/TrackMCHMID.h"
#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"
#include "DPLUtils/RootTreeWriter.h"

#include "SimulationDataFormat/MCCompLabel.h"
#include "DataFormatsMCH/Cluster.h"
#include "MCHTracking/Track.h"
#include "MCHTracking/TrackExtrap.h"
#include "MCHTracking/TrackParam.h"
#include "MCHTracking/TrackFitter.h"
#include "MCHAlign/Aligner.h"

#include "DetectorsCommonDataFormats/AlignParam.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "DetectorsCommonDataFormats/DetectorNameConf.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace std;

const int fgNCh = 10;
const int fgNDetElemCh[fgNCh] = {4, 4, 4, 4, 18, 18, 26, 26, 26, 26};
const int fgSNDetElemCh[fgNCh + 1] = {0, 4, 8, 12, 16, 34, 52, 78, 104, 130, 156};

using Container = std::vector<o2::mch::Track>;

struct mchClustersAOD {

  ccdb::CcdbApi ccdbApi;
  Service<ccdb::BasicCCDBManager> ccdb;
  parameters::GRPMagField* grpmag;
  TGeoManager* geo;
  mch::TrackFitter trackFitter;
  mch::geo::TransformationCreator transformation;
  mch::Aligner mAlign{};
  TStopwatch sw;
  Double_t weightRecord{1.0};

  double Reso_X = 0.0;
  double Reso_Y = 0.0;
  double ImproveCut = 6.0;

  map<int, math_utils::Transform3D> transformOld;
  map<int, math_utils::Transform3D> transformNew;

  string inputConfig = fmt::format("rofs:MCH/CLUSTERROFS;clusters:MCH/CLUSTERS");
  map<string, string> metadataRCT, headers;
  uint64_t ts{};

  Configurable<string> fConfigCcdbUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<string> fConfigColType{"collision-type", "pp", "Resolution specification for trackfitter"};
  Configurable<string> fFixChamber{"fix-chamber", "", "Fixing chamber"};
  Configurable<int> fRunNumber{"run-number", 539483, "Run number"};
  Configurable<bool> fDoNewGeo{"do-realign", false, "Transform to a given new geometry"};
  Configurable<bool> fDoEvaluation{"do-evaluation", false, "Enable storage of residuals"};
  Configurable<string> fConfigNewGeoFile{"new-geo", "o2sim_geometry-aligned.root", "New geometry for transformation"};
  Configurable<string> fOutputRecFileName{"outfile-data", "recDataFile.root", "Name of output data record file"};
  Configurable<string> fOutputConsFileName{"outfile-constraint", "recConsFile.root", "Name of output constraint record file"};

  void init(InitContext& ic)
  {

    // Load field and geometry informations here
    ccdbApi.init(fConfigCcdbUrl.value);
    ccdb->setURL(fConfigCcdbUrl.value);
    ccdb->setCaching(true);
    headers = ccdbApi.retrieveHeaders(Form("RCT/Info/RunInformation/%i", fRunNumber.value), metadataRCT, -1);
    ts = atol(headers["SOR"].c_str());
    grpmag = ccdb->getForTimeStamp<parameters::GRPMagField>("GLO/Config/GRPMagField", ts);
    base::Propagator::initFieldFromGRP(grpmag);
    mch::TrackExtrap::setField();
    mch::TrackExtrap::useExtrapV2();
    mAlign.SetBFieldOn(mch::TrackExtrap::isFieldON());
    trackFitter.smoothTracks(true);

    if (fConfigColType.value == "pp") {
      Reso_X = 0.4;
      Reso_Y = 0.4;
      ImproveCut = 6.0;
      LOG(info) << "Using pp parameter set for TrackFitter";
    } else {
      Reso_X = 0.2;
      Reso_Y = 0.2;
      ImproveCut = 4.0;
      LOG(info) << "Using PbPb parameter set for TrackFitter";
    }

    trackFitter.setChamberResolution(Reso_X, Reso_Y);
    trackFitter.useChamberResolution();

    mAlign.SetDoEvaluation(fDoEvaluation.value);
    mAlign.SetAllowedVariation(0, 2.0);
    mAlign.SetAllowedVariation(1, 0.3);
    mAlign.SetAllowedVariation(2, 0.002);
    mAlign.SetAllowedVariation(3, 2.0);

    // Fix chambers
    auto chambers = fFixChamber.value;
    for (int i = 0; i < chambers.length(); ++i) {
      if (chambers[i] == ',')
        continue;
      int chamber = chambers[i] - '0';
      LOG(info) << Form("%s%d", "Fixing chamber: ", chamber);
      mAlign.FixChamber(chamber);
    }

    mAlign.init(fOutputRecFileName.value, fOutputConsFileName.value);

    // Load reference geometry
    geo = ccdb->getForTimeStamp<TGeoManager>("GLO/Config/GeometryAligned", ts);
    transformation = mch::geo::transformationFromTGeoManager(*geo);
    for (int i = 0; i < 156; i++) {
      int iDEN = GetDetElemId(i);
      transformOld[iDEN] = transformation(iDEN);
    }

    if (fDoNewGeo.value) {
      // Load new geometry with which we want to check
      base::GeometryManager::loadGeometry(fConfigNewGeoFile.value);
      transformation = mch::geo::transformationFromTGeoManager(*gGeoManager);
      for (int i = 0; i < 156; i++) {
        int iDEN = GetDetElemId(i);
        transformNew[iDEN] = transformation(iDEN);
      }
    }

    sw.Stop();
    sw.Start(false);

    ic.services().get<CallbackService>().set<CallbackService::Id::Stop>([this]() {
      LOG(info) << "Saving records into ROOT file";

      if (fDoEvaluation.value) {
        drawHisto(*(mAlign.GetResTree()));
      }

      mAlign.terminate();
      sw.Stop();
      LOG(info) << "CPU time: " << sw.CpuTime();
    });
  }

  //_________________________________________________________________________________________________
  Int_t GetDetElemNumber(Int_t iDetElemId)
  {
    // get det element number from ID
    // get chamber and element number in chamber
    const Int_t iCh = iDetElemId / 100;
    const Int_t iDet = iDetElemId % 100;

    // make sure detector index is valid
    if (!(iCh > 0 && iCh <= fgNCh && iDet < fgNDetElemCh[iCh - 1])) {
      LOG(fatal) << "Invalid detector element id: " << iDetElemId;
    }

    // add number of detectors up to this chamber
    return iDet + fgSNDetElemCh[iCh - 1];
  }

  //_________________________________________________________________________________________________
  Int_t GetDetElemId(Int_t iDetElemNumber)
  {
    // make sure detector number is valid
    if (!(iDetElemNumber >= fgSNDetElemCh[0] &&
          iDetElemNumber < fgSNDetElemCh[fgNCh])) {
      LOG(fatal) << "Invalid detector element number: " << iDetElemNumber;
    }
    /// get det element number from ID
    // get chamber and element number in chamber
    int iCh = 0;
    int iDet = 0;
    for (int i = 1; i <= fgNCh; i++) {
      if (iDetElemNumber < fgSNDetElemCh[i]) {
        iCh = i;
        iDet = iDetElemNumber - fgSNDetElemCh[i - 1];
        break;
      }
    }

    // make sure detector index is valid
    if (!(iCh > 0 && iCh <= fgNCh && iDet < fgNDetElemCh[iCh - 1])) {
      LOG(fatal) << "Invalid detector element id: " << 100 * iCh + iDet;
    }

    // add number of detectors up to this chamber
    return 100 * iCh + iDet;
  }

  //_________________________________________________________________________________________________
  bool RemoveTrack(mch::Track& track, double ImproveCut)
  {
    double maxChi2Cluster = 2 * ImproveCut * ImproveCut;
    bool removeTrack = false;

    try {
      trackFitter.fit(track, false);
    } catch (exception const& e) {
      removeTrack = true;
      return removeTrack;
    }

    auto itStartingParam = std::prev(track.rend());

    while (true) {
      try {
        trackFitter.fit(track, true, false, (itStartingParam == track.rbegin()) ? nullptr : &itStartingParam);
      } catch (exception const&) {
        removeTrack = true;
        break;
      }

      double worstLocalChi2 = -1.0;

      track.tagRemovableClusters(0x1F, false);
      auto itWorstParam = track.end();

      for (auto itParam = track.begin(); itParam != track.end(); ++itParam) {
        if (itParam->getLocalChi2() > worstLocalChi2) {
          worstLocalChi2 = itParam->getLocalChi2();
          itWorstParam = itParam;
        }
      }

      if (worstLocalChi2 < maxChi2Cluster)
        break;

      if (!itWorstParam->isRemovable()) {
        removeTrack = true;
        track.removable();
        break;
      }

      auto itNextParam = track.removeParamAtCluster(itWorstParam);
      auto itNextToNextParam = (itNextParam == track.end()) ? itNextParam : std::next(itNextParam);
      itStartingParam = track.rbegin();

      if (track.getNClusters() < 10) {
        removeTrack = true;
        break;
      } else {
        while (itNextToNextParam != track.end()) {
          if (itNextToNextParam->getClusterPtr()->getChamberId() != itNextParam->getClusterPtr()->getChamberId()) {
            itStartingParam = std::make_reverse_iterator(++itNextParam);
            break;
          }
          ++itNextToNextParam;
        }
      }
    }

    if (!removeTrack) {
      for (auto& param : track) {
        param.setParameters(param.getSmoothParameters());
        param.setCovariances(param.getSmoothCovariances());
      }
    }

    return removeTrack;
  }

  //_________________________________________________________________________________________________
  void drawHisto(TTree& Res_Tree, string outFileName = "Residual_QA")
  {

    double deNumber[156];

    // Saving plots
    string PlotFiles_name = Form("%s%s", outFileName.c_str(), "_results.root");
    TFile* PlotFiles = TFile::Open(PlotFiles_name.c_str(), "RECREATE");

    int RefClDetElem;
    int RefClDetElemNumber;
    float RefClusterX;
    float RefClusterY;
    float RefTrackX;
    float RefTrackY;
    float RefTrackSlopeX;
    float RefTrackSlopeY;
    Res_Tree.SetBranchAddress("fClusterX", &RefClusterX);
    Res_Tree.SetBranchAddress("fClusterY", &RefClusterY);
    Res_Tree.SetBranchAddress("fTrackX", &RefTrackX);
    Res_Tree.SetBranchAddress("fTrackY", &RefTrackY);
    Res_Tree.SetBranchAddress("fTrackSlopeX", &RefTrackSlopeX);
    Res_Tree.SetBranchAddress("fTrackSlopeY", &RefTrackSlopeY);
    Res_Tree.SetBranchAddress("fClDetElem", &RefClDetElem);
    Res_Tree.SetBranchAddress("fClDetElemNumber", &RefClDetElemNumber);

    TH1F* Histos_Res[2][11];
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 11; j++) {
        if (i == 0) {
          Histos_Res[i][j] = new TH1F(Form("%s%d", "Residual_X_Ch", j), Form("%s%d", "Residual_x_Ch", j), 200, -5, 5);
        }

        if (i == 1) {
          Histos_Res[i][j] = new TH1F(Form("%s%d", "Residual_Y_Ch", j), Form("%s%d", "Residual_y_Ch", j), 200, -5, 5);
        }
      }
    }

    TH1F* Histos_DETRes[2][156];
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 156; j++) {
        if (i == 0)
          Histos_DETRes[i][j] = new TH1F(Form("%s%d", "Hist_x_DET", j + 1), Form("%s%d", "Hist_x_DET", j + 1), 200, -5, 5);
        if (i == 1)
          Histos_DETRes[i][j] = new TH1F(Form("%s%d", "Hist_y_DET", j + 1), Form("%s%d", "Hist_y_DET", j + 1), 200, -5, 5);
      }
    }

    int Ref_NbEntries = Res_Tree.GetEntries();
    for (int i = 0; i < Ref_NbEntries; i++) {

      Res_Tree.GetEntry(i);
      double Res_X = RefClusterX - RefTrackX;
      double Res_Y = RefClusterY - RefTrackY;

      for (int iCh = 0; iCh < 11; iCh++) {
        if (iCh == 0) {
          Histos_Res[0][iCh]->Fill(Res_X);
          Histos_Res[1][iCh]->Fill(Res_Y);
        } else {
          if (iCh == int(RefClDetElem / 100)) {
            Histos_Res[0][iCh]->Fill(Res_X);
            Histos_Res[1][iCh]->Fill(Res_Y);
          }
        }
      }

      for (int iDEN = 0; iDEN < 156; iDEN++) {
        if (RefClDetElemNumber == iDEN) {
          Histos_DETRes[0][iDEN]->Fill(Res_X);
          Histos_DETRes[1][iDEN]->Fill(Res_Y);
        }
      }
    }

    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 11; j++) {
        if (i == 0)
          PlotFiles->WriteObjectAny(Histos_Res[i][j], "TH1F", Form("%s%d", "Residual_X_Ch", j));
        if (i == 1)
          PlotFiles->WriteObjectAny(Histos_Res[i][j], "TH1F", Form("%s%d", "Residual_Y_Ch", j));
      }
    }

    double ResX_mean[156] = {};
    double ResX_err[156] = {};
    double ResY_mean[156] = {};
    double ResY_err[156] = {};
    for (int iDEN = 0; iDEN < 156; iDEN++) {
      deNumber[iDEN] = iDEN + 0.5;

      ResX_mean[iDEN] = Histos_DETRes[0][iDEN]->GetMean();
      ResX_err[iDEN] = Histos_DETRes[0][iDEN]->GetRMS();

      ResY_mean[iDEN] = Histos_DETRes[1][iDEN]->GetMean();
      ResY_err[iDEN] = Histos_DETRes[1][iDEN]->GetRMS();
    }

    TGraphErrors* graphResX = new TGraphErrors(156, deNumber, ResX_mean, nullptr, ResX_err);
    TGraphErrors* graphResY = new TGraphErrors(156, deNumber, ResY_mean, nullptr, ResY_err);

    TCanvas* graphRes = new TCanvas("graphRes", "graphRes", 1200, 1600);
    TLine limLine(4, -5, 4, 5);
    TH1F* gHisto = new TH1F("gHisto", "Residuals", 161, 0, 160);
    gHisto->SetXTitle("Det. Elem. Number");

    graphRes->Divide(1, 2);

    graphRes->cd(1);
    gHisto->SetYTitle("TrackX - ClusterX (cm)");
    gHisto->GetYaxis()->SetRangeUser(-5.0, 5.0);
    gHisto->DrawCopy("goff");
    graphResX->SetMarkerStyle(8);
    graphResX->SetMarkerSize(0.7);
    graphResX->SetLineColor(kBlue);
    graphResX->Draw("PZsame goff");
    limLine.DrawLine(4, -5, 4, 5);
    limLine.DrawLine(8, -5, 8, 5);
    limLine.DrawLine(12, -5, 12, 5);
    limLine.DrawLine(16, -5, 16, 5);
    limLine.DrawLine(16 + 18, -5, 16 + 18, 5);
    limLine.DrawLine(16 + 2 * 18, -5, 16 + 2 * 18, 5);
    limLine.DrawLine(16 + 2 * 18 + 26, -5, 16 + 2 * 18 + 26, 5);
    limLine.DrawLine(16 + 2 * 18 + 2 * 26, -5, 16 + 2 * 18 + 2 * 26, 5);
    limLine.DrawLine(16 + 2 * 18 + 3 * 26, -5, 16 + 2 * 18 + 3 * 26, 5);

    graphRes->cd(2);
    gHisto->SetYTitle("TrackY - ClusterY (cm)");
    gHisto->GetYaxis()->SetRangeUser(-5.0, 5.0);
    gHisto->DrawCopy("goff");
    graphResY->SetMarkerStyle(8);
    graphResY->SetMarkerSize(0.7);
    graphResY->SetLineColor(kBlue);
    graphResY->Draw("PZsame goff");
    limLine.DrawLine(4, -5, 4, 5);
    limLine.DrawLine(8, -5, 8, 5);
    limLine.DrawLine(12, -5, 12, 5);
    limLine.DrawLine(16, -5, 16, 5);
    limLine.DrawLine(16 + 18, -5, 16 + 18, 5);
    limLine.DrawLine(16 + 2 * 18, -5, 16 + 2 * 18, 5);
    limLine.DrawLine(16 + 2 * 18 + 26, -5, 16 + 2 * 18 + 26, 5);
    limLine.DrawLine(16 + 2 * 18 + 2 * 26, -5, 16 + 2 * 18 + 2 * 26, 5);
    limLine.DrawLine(16 + 2 * 18 + 3 * 26, -5, 16 + 2 * 18 + 3 * 26, 5);

    PlotFiles->WriteObjectAny(graphRes, "TCanvas", "ResGraph");
    PlotFiles->Close();
  }

  template <typename TTrack, typename TClusters>
  void runProcessTracks(TTrack const& track, TClusters const& clusters)
  {

    int clIndex = -1;
    mch::Track convertedTrack;

    for (auto& cluster : clusters) {

      clIndex += 1;

      mch::Cluster* mch_cluster = new mch::Cluster();
      mch_cluster->x = cluster.x();
      mch_cluster->y = cluster.y();
      mch_cluster->z = cluster.z();

      if (fDoNewGeo.value) {
        math_utils::Point3D<double> local;
        math_utils::Point3D<double> master;

        master.SetXYZ(cluster.x(), cluster.y(), cluster.z());

        transformOld[cluster.deId()].MasterToLocal(master, local);
        transformNew[cluster.deId()].LocalToMaster(local, master);

        mch_cluster->x = master.x();
        mch_cluster->y = master.y();
        mch_cluster->z = master.z();
      }

      uint32_t ClUId = mch::Cluster::buildUniqueId(int(cluster.deId() / 100) - 1, cluster.deId(), clIndex);
      mch_cluster->uid = ClUId;

      mch_cluster->ex = cluster.isGoodX() ? 0.2 : 10.0;
      mch_cluster->ey = cluster.isGoodY() ? 0.2 : 10.0;

      convertedTrack.createParamAtCluster(*mch_cluster);
    }

    if (convertedTrack.getNClusters() >= 9) {
      // Erase removable track
      if (!RemoveTrack(convertedTrack, ImproveCut)) {
        mAlign.ProcessTrack(convertedTrack, transformation, true, weightRecord);
      }
    }
  }

  void processTracks(aod::FwdTracks::iterator const& track, aod::FwdTrkCls const& clusters)
  {
    runProcessTracks(track, clusters);
  }

  PROCESS_SWITCH(mchClustersAOD, processTracks, "Process tracks", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<mchClustersAOD>(cfgc)};
}
