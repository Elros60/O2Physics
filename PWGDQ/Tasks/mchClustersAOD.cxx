#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"

#include <cmath>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>
#include <ctime>
#include <iostream>

#include <gsl/span>

// #include <TSystem.h>
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

#include "CommonConstants/LHCConstants.h"
#include "CommonUtils/NameConf.h"
#include "DataFormatsMCH/Cluster.h"
#include "DataFormatsMCH/Digit.h"
#include "DataFormatsMCH/ROFRecord.h"
#include "DataFormatsMCH/TrackMCH.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/Logger.h"
#include "MCHGeometryTransformer/Transformations.h"
#include "MCHMappingInterface/Segmentation.h"
#include "MCHTracking/TrackExtrap.h"
#include "MCHTracking/TrackParam.h"
#include "ReconstructionDataFormats/TrackMCHMID.h"

#include "SimulationDataFormat/MCCompLabel.h"

//  Test with alignment codes
#include "DataFormatsMCH/Cluster.h"
#include "MCHTracking/Track.h"
#include "MCHTracking/TrackExtrap.h"
#include "MCHTracking/TrackParam.h"
#include "MCHTracking/TrackFitter.h"

#include "DetectorsCommonDataFormats/AlignParam.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "DetectorsCommonDataFormats/DetectorNameConf.h"


using namespace o2;
using namespace o2::framework;

struct mchClustersAOD
{
	void process(aod::FwdTracks const& tracks, aod::FwdTrkCls const& clusters){
		for(auto& cluster : clusters){
			std::cout << std::endl;
			std::cout << std::endl;
			//LOG(info) << Form("%s: %f","Track index",cluster.fwdtrack());
			LOG(info) << Form("%s: %f","pos x",cluster.x());
			LOG(info) << Form("%s: %f","pos y",cluster.y());
			LOG(info) << Form("%s: %f","pos z",cluster.z());
			LOG(info) << Form("%s: %d","DEId",cluster.deId());

		}

	}
	
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<mchClustersAOD>(cfgc)};
}