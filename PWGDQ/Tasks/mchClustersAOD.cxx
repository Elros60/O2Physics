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
#include "DetectorsBase/GRPGeomHelper.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/Logger.h"
#include "MCHGeometryTransformer/Transformations.h"
#include "MCHMappingInterface/Segmentation.h"
#include "MCHTracking/TrackExtrap.h"
#include "MCHTracking/TrackParam.h"
#include "ReconstructionDataFormats/TrackMCHMID.h"
#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"

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
using namespace std;

const int fgNCh = 10;
const int fgNDetElemCh[fgNCh] = {4, 4, 4, 4, 18, 18, 26, 26, 26, 26};
const int fgSNDetElemCh[fgNCh + 1] = {0,  4,  8,   12,  16, 34,
                                      52, 78, 104, 130, 156}; 



Int_t GetDetElemNumber(Int_t iDetElemId);
Int_t GetDetElemId(Int_t iDetElemNumber);


struct mchClustersAOD
{
	vector<mch::Track*> mch_tracks;
	mch::Track* convertedTrack;
	ccdb::CcdbApi ccdbApi;
	Service<ccdb::BasicCCDBManager> ccdb;
	mch::geo::TransformationCreator transformation;
	parameters::GRPMagField* grpmag;
	TGeoManager* geo;
	mch::TrackFitter *trackFitter;
	string inputConfig = fmt::format("rofs:MCH/CLUSTERROFS;clusters:MCH/CLUSTERS");
	map<string, string> metadataRCT, headers;
	uint64_t ts{};

	map<int, math_utils::Transform3D> transformOld;
	map<int, math_utils::Transform3D> transformNew;

	double Reso_X = 0.4;
    double Reso_Y = 0.4;
    int runNumber = 539483;
	
	void init(InitContext&){

		//Load field and geometry informations here
		ccdbApi.init("http://alice-ccdb.cern.ch");
    	ccdb->setURL("http://alice-ccdb.cern.ch");
    	ccdb->setCaching(true);

    	auto inputs = select(inputConfig.c_str());
		auto ccdbRequest = make_shared<base::GRPGeomRequest>(		false,                             // orbitResetTime
                                                              		false,                             // GRPECS=true
                                                                	false,                             // GRPLHCIF
                                                                	false,                             // GRPMagField
                                                                	false,                             // askMatLUT
                                                                	base::GRPGeomRequest::Aligned, // geometry
                                                                	inputs);

		headers = ccdbApi.retrieveHeaders(Form("RCT/Info/RunInformation/%i", runNumber), metadataRCT, -1);
      	ts = atol(headers["SOR"].c_str());

    	grpmag = ccdb->getForTimeStamp<parameters::GRPMagField>("GLO/Config/GRPMagField", ts);
		base::Propagator::initFieldFromGRP(grpmag);
		mch::TrackExtrap::setField();
  		mch::TrackExtrap::useExtrapV2();
 
  		trackFitter = new mch::TrackFitter();
  		trackFitter->initField(grpmag->getL3Current(), grpmag->getDipoleCurrent());
  		trackFitter->smoothTracks(true);
  		trackFitter->setChamberResolution(Reso_X, Reso_Y);
  		trackFitter->useChamberResolution();

  		//Load reference geometry
  		geo = ccdb->getForTimeStamp<TGeoManager>("GLO/Config/GeometryAligned", ts);
  		//base::GRPGeomHelper::instance().setRequest(ccdbRequest);
  		transformation = mch::geo::transformationFromTGeoManager(*geo);
  		for (int i = 0; i < 156; i++) { 
    		int iDEN = GetDetElemId(i);
    	transformOld[iDEN] = transformation(iDEN);
  		}

  		//Load new geometry with which we want to check
  		base::GeometryManager::loadGeometry("o2sim_geometry-aligned.root");
  		transformation = mch::geo::transformationFromTGeoManager(*gGeoManager);
    	for (int i = 0; i < 156; i++) {
      		int iDEN = GetDetElemId(i);
      		transformNew[iDEN] = transformation(iDEN);
    	}

    	//Print loaded geometries to check
    	/*
    	for (int i = 0; i < 156; i++) { 
			int iDEN = GetDetElemId(i);
			auto transform3D_Old = transformOld[iDEN];
			auto transform3D_New = transformNew[iDEN];


			TMatrixD MTransOld(3,4);
			TMatrixD MTransNew(3,4);

			transform3D_Old.GetTransformMatrix(MTransOld);
			transform3D_New.GetTransformMatrix(MTransNew);
			cout << "===================================================" <<endl;
			cout << "DET ID: " << iDEN <<endl;
			MTransOld.Print();
			MTransNew.Print();
			cout << "===================================================" <<endl;
			cout << endl;

	    }
	    */

	};
	

	void process(aod::FwdTracks::iterator const& track, aod::FwdTrkCls const& clusters){

		convertedTrack = new mch::Track();
		int clIndex = -1;

		for(auto& cluster : clusters){

			LOG(info) << Form("%s%lld","Processing track: ",cluster.fwdtrack().globalIndex());
			clIndex += 1;
			LOG(info) << Form("%s%d","DEId: ",cluster.deId());
			LOG(info) << Form("%s%f","pos x: ",cluster.x());
			LOG(info) << Form("%s%f","pos y: ",cluster.y());
			LOG(info) << Form("%s%f","pos z: ",cluster.z());

			math_utils::Point3D<double> local;
			math_utils::Point3D<double> master;

			master.SetXYZ(cluster.x(), cluster.y(), cluster.z());

			transformOld[cluster.deId()].MasterToLocal(master, local);
			transformNew[cluster.deId()].LocalToMaster(local, master);

			//Refit
			mch::Cluster* mch_cluster = new mch::Cluster();
			mch_cluster->x = master.x();
			mch_cluster->y = master.y();
			mch_cluster->z = master.z();
			uint32_t ClUId = mch::Cluster::buildUniqueId(int(cluster.deId()/100),cluster.deId()%100,clIndex);
			mch_cluster->uid = ClUId;
			//std::cout << ClUId << std::endl;
			std::cout << "Cluster index: " << mch::Cluster::getClusterIndex(ClUId) << std::endl;
			//std::cout << "Chamber ID: " << mch::Cluster::getChamberId(ClUId) << std::endl;
			//std::cout << "DEId: " << mch::Cluster::getDEId(ClUId) << std::endl;
			mch_cluster->ex = cluster.isGoodX() ? 0.2 : 10.0;
			mch_cluster->ey = cluster.isGoodY() ? 0.2 : 10.0;
			convertedTrack->createParamAtCluster(*mch_cluster);
	
		}


		if(convertedTrack && convertedTrack->getNClusters()!=0){
			LOG(info) << "Track before refit: ";
			convertedTrack->print();
			mch_tracks.push_back(convertedTrack);
			try{
				LOG(info) << "Trying to refit current track...";
    			trackFitter->fit(*convertedTrack, false);
  			}catch(exception const& e){
    			LOG(info) << "Track to be removed!";
  			}
  			LOG(info) << "Track after refit: ";
  			convertedTrack->print();

		}


	}
	
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<mchClustersAOD>(cfgc)};
}



Int_t GetDetElemNumber(Int_t iDetElemId) {
  /// get det element number from ID
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


Int_t GetDetElemId(Int_t iDetElemNumber) {
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
