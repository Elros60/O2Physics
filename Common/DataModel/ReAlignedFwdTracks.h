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
/// \file   ReAlignedFwdTracks.h
/// \author Chi Zhang <chi.zhang@cern.ch>
///
/// \brief  Declaration of the table for the re-aligned forward muon tracks.
///

#ifndef COMMON_DATAMODEL_REALIGNEDFWDTRACKS_H_
#define COMMON_DATAMODEL_REALIGNEDFWDTRACKS_H_

#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace realignedfwdtracks
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);                                              //!
DECLARE_SOA_COLUMN(TrackType, trackType, uint8_t);                                           //! Type of track. See enum ForwardTrackTypeEnum
DECLARE_SOA_COLUMN(X, x, float);                                                             //! TrackParFwd parameter x
DECLARE_SOA_COLUMN(Y, y, float);                                                             //! TrackParFwd parameter y
DECLARE_SOA_COLUMN(Z, z, float);                                                             //! TrackParFwd propagation parameter z
DECLARE_SOA_COLUMN(Phi, phi, float);                                                         //! TrackParFwd parameter phi; (i.e. pt pointing direction)
DECLARE_SOA_COLUMN(Tgl, tgl, float);                                                         //! TrackParFwd parameter tan(\lamba); (\lambda = 90 - \theta_{polar})
DECLARE_SOA_COLUMN(Signed1Pt, signed1Pt, float);                                             //! TrackParFwd parameter: charged inverse transverse momentum; (q/pt)
DECLARE_SOA_COLUMN(NClusters, nClusters, int8_t);                                            //! Number of clusters
DECLARE_SOA_COLUMN(MFTClusterSizesAndTrackFlags, mftClusterSizesAndTrackFlags, uint64_t);    //! Cluster sizes per track, stored per layer (each 6 bits). Remaining 4 bits for MFT flags
DECLARE_SOA_COLUMN(Chi2, chi2, float);                                                       //! Track chi^2
DECLARE_SOA_COLUMN(PDca, pDca, float);                                                       //! PDca for MUONStandalone
DECLARE_SOA_COLUMN(RAtAbsorberEnd, rAtAbsorberEnd, float);                                   //! RAtAbsorberEnd for MUONStandalone tracks and GlobalMuonTrackstracks
DECLARE_SOA_COLUMN(Chi2MatchMCHMID, chi2MatchMCHMID, float);                                 //! MCH-MID Match Chi2 for MUONStandalone tracks
DECLARE_SOA_COLUMN(Chi2MatchMCHMFT, chi2MatchMCHMFT, float);                                 //! MCH-MFT Match Chi2 for GlobalMuonTracks
DECLARE_SOA_COLUMN(MatchScoreMCHMFT, matchScoreMCHMFT, float);                               //! MCH-MFT Machine Learning Matching Score for GlobalMuonTracks
DECLARE_SOA_SELF_INDEX_COLUMN_FULL(MCHTrack, matchMCHTrack, int, "FwdTracks_MatchMCHTrack"); //! Index of matching MCH track for GlobalMuonTracks and GlobalForwardTracks
DECLARE_SOA_COLUMN(MCHBitMap, mchBitMap, uint16_t);                                          //! Fired muon trackig chambers bitmap
DECLARE_SOA_COLUMN(MIDBitMap, midBitMap, uint8_t);                                           //! MID bitmap: non-bending plane (4bit), bending plane (4bit)
DECLARE_SOA_COLUMN(MIDBoards, midBoards, uint32_t);                                          //! Local boards on each MID plane (8 bits per plane)
DECLARE_SOA_COLUMN(TrackTime, trackTime, float);                                             //! Estimated time of the track in ns wrt collision().bc() or ambiguoustrack.bcSlice()[0]
DECLARE_SOA_COLUMN(TrackTimeRes, trackTimeRes, float);                                       //! Resolution of the track time in ns
DECLARE_SOA_DYNAMIC_COLUMN(Sign, sign,                                                       //! Sign of the track eletric charge
                           [](float signed1Pt) -> short { return (signed1Pt > 0) ? 1 : -1; });
DECLARE_SOA_DYNAMIC_COLUMN(IsCA, isCA, //! Returns true if used track-finding algorithm was Cellular Automaton
                           [](uint64_t mftClusterSizesAndTrackFlags) -> bool { return mftClusterSizesAndTrackFlags & (0x1ULL << 60); });
DECLARE_SOA_EXPRESSION_COLUMN(Eta, eta, float, //!
                              -1.f * nlog(ntan(o2::constants::math::PIQuarter - 0.5f * natan(aod::fwdtrack::tgl))));
DECLARE_SOA_EXPRESSION_COLUMN(Pt, pt, float, //!
                              ifnode(nabs(aod::fwdtrack::signed1Pt) < o2::constants::math::Almost0, o2::constants::math::VeryBig, nabs(1.f / aod::fwdtrack::signed1Pt)));
DECLARE_SOA_EXPRESSION_COLUMN(P, p, float, //!
                              ifnode((nabs(aod::fwdtrack::signed1Pt) < o2::constants::math::Almost0) || (nabs(o2::constants::math::PIQuarter - 0.5f * natan(aod::fwdtrack::tgl)) < o2::constants::math::Almost0), o2::constants::math::VeryBig, 0.5f * (ntan(o2::constants::math::PIQuarter - 0.5f * natan(aod::fwdtrack::tgl)) + 1.f / ntan(o2::constants::math::PIQuarter - 0.5f * natan(aod::fwdtrack::tgl))) / nabs(aod::fwdtrack::signed1Pt)));
DECLARE_SOA_DYNAMIC_COLUMN(Px, px, //!
                           [](float pt, float phi) -> float {
                             return pt * std::cos(phi);
                           });
DECLARE_SOA_DYNAMIC_COLUMN(Py, py, //!
                           [](float pt, float phi) -> float {
                             return pt * std::sin(phi);
                           });
DECLARE_SOA_DYNAMIC_COLUMN(Pz, pz, //!
                           [](float pt, float tgl) -> float {
                             return pt * tgl;
                           });
DECLARE_SOA_DYNAMIC_COLUMN(MIDBoardCh1, midBoardCh1, //!
                           [](uint32_t midBoards) -> int {
                             return static_cast<int>(midBoards & 0xFF);
                           });
DECLARE_SOA_DYNAMIC_COLUMN(MIDBoardCh2, midBoardCh2, //!
                           [](uint32_t midBoards) -> int {
                             return static_cast<int>((midBoards >> 8) & 0xFF);
                           });
DECLARE_SOA_DYNAMIC_COLUMN(MIDBoardCh3, midBoardCh3, //!
                           [](uint32_t midBoards) -> int {
                             return static_cast<int>((midBoards >> 16) & 0xFF);
                           });
DECLARE_SOA_DYNAMIC_COLUMN(MIDBoardCh4, midBoardCh4, //!
                           [](uint32_t midBoards) -> int {
                             return static_cast<int>((midBoards >> 24) & 0xFF);
                           });

DECLARE_SOA_COLUMN(SigmaX, sigmaX, float);        //! Covariance matrix
DECLARE_SOA_COLUMN(SigmaY, sigmaY, float);        //! Covariance matrix
DECLARE_SOA_COLUMN(SigmaPhi, sigmaPhi, float);    //! Covariance matrix
DECLARE_SOA_COLUMN(SigmaTgl, sigmaTgl, float);    //! Covariance matrix
DECLARE_SOA_COLUMN(Sigma1Pt, sigma1Pt, float);    //! Covariance matrix
DECLARE_SOA_COLUMN(RhoXY, rhoXY, int8_t);         //! Covariance matrix in compressed form
DECLARE_SOA_COLUMN(RhoPhiX, rhoPhiX, int8_t);     //! Covariance matrix in compressed form
DECLARE_SOA_COLUMN(RhoPhiY, rhoPhiY, int8_t);     //! Covariance matrix in compressed form
DECLARE_SOA_COLUMN(RhoTglX, rhoTglX, int8_t);     //! Covariance matrix in compressed form
DECLARE_SOA_COLUMN(RhoTglY, rhoTglY, int8_t);     //! Covariance matrix in compressed form
DECLARE_SOA_COLUMN(RhoTglPhi, rhoTglPhi, int8_t); //! Covariance matrix in compressed form
DECLARE_SOA_COLUMN(Rho1PtX, rho1PtX, int8_t);     //! Covariance matrix in compressed form
DECLARE_SOA_COLUMN(Rho1PtY, rho1PtY, int8_t);     //! Covariance matrix in compressed form
DECLARE_SOA_COLUMN(Rho1PtPhi, rho1PtPhi, int8_t); //! Covariance matrix in compressed form
DECLARE_SOA_COLUMN(Rho1PtTgl, rho1PtTgl, int8_t); //! Covariance matrix in compressed form

DECLARE_SOA_EXPRESSION_COLUMN(CXX, cXX, float, //!
                              aod::fwdtrack::sigmaX* aod::fwdtrack::sigmaX);
DECLARE_SOA_EXPRESSION_COLUMN(CXY, cXY, float, //!
                              (aod::fwdtrack::rhoXY / 128.f) * (aod::fwdtrack::sigmaX * aod::fwdtrack::sigmaY));
DECLARE_SOA_EXPRESSION_COLUMN(CYY, cYY, float, //!
                              aod::fwdtrack::sigmaY* aod::fwdtrack::sigmaY);
DECLARE_SOA_EXPRESSION_COLUMN(CPhiX, cPhiX, float, //!
                              (aod::fwdtrack::rhoPhiX / 128.f) * (aod::fwdtrack::sigmaPhi * aod::fwdtrack::sigmaX));
DECLARE_SOA_EXPRESSION_COLUMN(CPhiY, cPhiY, float, //!
                              (aod::fwdtrack::rhoPhiY / 128.f) * (aod::fwdtrack::sigmaPhi * aod::fwdtrack::sigmaY));
DECLARE_SOA_EXPRESSION_COLUMN(CPhiPhi, cPhiPhi, float, //!
                              aod::fwdtrack::sigmaPhi* aod::fwdtrack::sigmaPhi);
DECLARE_SOA_EXPRESSION_COLUMN(CTglX, cTglX, float, //!
                              (aod::fwdtrack::rhoTglX / 128.f) * (aod::fwdtrack::sigmaTgl * aod::fwdtrack::sigmaX));
DECLARE_SOA_EXPRESSION_COLUMN(CTglY, cTglY, float, //!
                              (aod::fwdtrack::rhoTglY / 128.f) * (aod::fwdtrack::sigmaTgl * aod::fwdtrack::sigmaY));
DECLARE_SOA_EXPRESSION_COLUMN(CTglPhi, cTglPhi, float, //!
                              (aod::fwdtrack::rhoTglPhi / 128.f) * (aod::fwdtrack::sigmaTgl * aod::fwdtrack::sigmaPhi));
DECLARE_SOA_EXPRESSION_COLUMN(CTglTgl, cTglTgl, float, //!
                              aod::fwdtrack::sigmaTgl* aod::fwdtrack::sigmaTgl);
DECLARE_SOA_EXPRESSION_COLUMN(C1PtY, c1PtY, float, //!
                              (aod::fwdtrack::rho1PtY / 128.f) * (aod::fwdtrack::sigma1Pt * aod::fwdtrack::sigmaY));
DECLARE_SOA_EXPRESSION_COLUMN(C1PtX, c1PtX, float, //!
                              (aod::fwdtrack::rho1PtX / 128.f) * (aod::fwdtrack::sigma1Pt * aod::fwdtrack::sigmaX));
DECLARE_SOA_EXPRESSION_COLUMN(C1PtPhi, c1PtPhi, float, //!
                              (aod::fwdtrack::rho1PtPhi / 128.f) * (aod::fwdtrack::sigma1Pt * aod::fwdtrack::sigmaPhi));
DECLARE_SOA_EXPRESSION_COLUMN(C1PtTgl, c1PtTgl, float, //!
                              (aod::fwdtrack::rho1PtTgl / 128.f) * (aod::fwdtrack::sigma1Pt * aod::fwdtrack::sigmaTgl));
DECLARE_SOA_EXPRESSION_COLUMN(C1Pt21Pt2, c1Pt21Pt2, float, //!
                              aod::fwdtrack::sigma1Pt* aod::fwdtrack::sigma1Pt);

} // namespace realignedfwdtracks

} // namespace o2::aod

#endif // COMMON_DATAMODEL_REALIGNEDFWDTRACKS_H_