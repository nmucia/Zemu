// -*- C++ -*-
////////////////////////////////////////////////////////////////
// Package:    Spots
// Class:      Spots
//
// Description:  This is an analysis for looking at beam spots^
//
// Implementation:
//     [Notes on implementation]
//
//
// Original Author:  Dale Stentz
//         Created:  Tue Aug 14 17:14:59 CDT 2012
//
// ^Disclaimer:  Description may be misleading or false
//
// $Id$
////////////////////////////////////////////////////////////////

//c/c++
// system include files
#include <memory>
#include <vector>
#include <string>
#include <map>
#include <utility>
#include <math.h>
#include <cmath>
#include <iostream>

//CMSSW
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"        //djs
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/Common/interface/View.h"             //djs
#include "DataFormats/Common/interface/RefToBaseVector.h"  //djs

#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryVector/interface/LocalVector.h"
#include "DataFormats/Math/interface/Vector.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "DataFormats/Math/interface/deltaPhi.h"    //djs -removed

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "DataFormats/Candidate/interface/Candidate.h"          //djs
#include "DataFormats/Candidate/interface/CandidateFwd.h"

//djs - Event Shapes
#include <PhysicsTools/CandUtils/interface/EventShapeVariables.h>
#include "PhysicsTools/CandUtils/interface/EventShapeVariables.h"
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TVectorD.h>
#include <Math/VectorUtil.h>


//trigger info
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "DataFormats/TrackReco/interface/TrackBase.h"   //djs
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"  //djs
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"

#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "SimDataFormats/JetMatching/interface/JetFlavour.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"
#include "SimDataFormats/JetMatching/interface/MatchedPartons.h"
#include "SimDataFormats/JetMatching/interface/JetMatchedPartons.h"

#include "DataFormats/JetReco/interface/GenJet.h"             //djs
#include "DataFormats/JetReco/interface/GenJetCollection.h"   //djs
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/JetID.h"
#include "DataFormats/JetReco/interface/JPTJetCollection.h"
#include "DataFormats/JetReco/interface/JPTJet.h"

#include "DataFormats/JetReco/interface/JetExtendedAssociation.h"  //djs

#include "DataFormats/VertexReco/interface/Vertex.h"     //djs
#include "DataFormats/VertexReco/interface/VertexFwd.h"  //djs
#include "DataFormats/BeamSpot/interface/BeamSpot.h"     //djs

//djs - jec
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

//djs - MET
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"

#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
// electron stuff
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "RecoParticleFlow/PFProducer/interface/Utils.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "EgammaAnalysis/ElectronTools/interface/EGammaMvaEleEstimator.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "EgammaAnalysis/ElectronTools/interface/PFIsolationEstimator.h"
#include "EgammaAnalysis/ElectronTools/interface/ElectronEffectiveArea.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "CommonTools/ParticleFlow/test/PFIsoReaderDemo.h"

//andy vertex info

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TrackingTools/TransientTrack/interface/TrackTransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"


#include <algorithm>

#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TString.h"     //djs
#include "TStopwatch.h"  //djs
#include "TProfile.h"
#include "TRandom.h"
#include "TTree.h"
#include "TNtuple.h"
#include "Math/LorentzVector.h"
#include "TVector2.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TVector2.h"
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include <TH1F.h>
#include <TH2F.h>
////////////////////////////////////////////////////////////////

using namespace edm;
using namespace std;
using namespace reco;
//math?

////////////////////////////////////////////////////////////////

const char
  PROCESS_no = '0',
  PROCESS_di = '2',
  PROCESS_lj = '1',
  DILEP_no = 'X',
  DILEP_ee = 'e',
  DILEP_mm = 'm',
  DILEP_em = 'b',
  LEPJET_no = 'X',
  LEPJET_ej = 'e',
  LEPJET_mj = 'm',
  SIGN_os = 'o',
  SIGN_ss = 's',
  SIGN_na = 'n';

////////////////////////////////////////////////////////////////

const bool
  IS_DATA = false,
  SKIP_SS_CUT = false,
  SKIP_Z_MASS_CUT = true,
  SKIP_MET_CUT = false,
  SKIP_HT_CUT = true;

////////////////////////////////////////////////////////////////

const int
  N_DEBUG = 10,
  N_HIST = 329,  //Updated for Nick Code
  N_JETS = 20, 
  N_TAGS = 10,
  CUT_TRACK_VHITS = 10,
  ID_0 = 0,
  ID_U = 1,
  ID_D = 2,
  ID_S = 3,
  ID_C = 4,
  ID_B = 5,
  ID_T = 6,
  ID_G = 21,
  ID_H = 25,
  ID_W = 24,
  MET_NBINS = 60,
  PHI_NBINS = 32,
  DPHI_NBINS = 20,
  LPT_NBINS = 60,
  JPT_NBINS = 50,
  BPT_NBINS = 75,
  HT_NBINS = 50,
  MT_NBINS = 50,
  M3J_NBINS = 60,
  ETA_NBINS = 25,
  DETA_NBINS = 25,
  DELTAR_NBINS = 40,
  ISO_NBINS = 25,
  VZ_NBINS = 40,
  DZ_NBINS = 50,
  DISCRIM_NBINS = 50,
  MASS_NBINS = 40;

////////////////////////////////////////////////////////////////

//check all these!
const double   
  R40 = 0.40,
  JUNK = 999999.9,
  DEFAULT_WEIGHT = 1.,
  CUT_LEPTON_TRIG_PT = 25.,     
  CUT_LEPTON_LOOSE_PT = 25.,
  CUT_ELEC_ETA = 2.4,
  CUT_MUON_ETA = 2.4,
  //CUT_CENTRAL_ETA = 1.4,   
  CUT_ELEC_ISO = 0.17,
  CUT_MUON_ISO = 0.15,
  CUT_MUON_CHI2 = 10.,
  CUT_ELEC_DXY = 0.04,
  CUT_MUON_DXY = 0.02,
  CUT_DELTAZ = 0.05,
  CUT_JET_PT = 20.,
  CUT_JET_ETA = 2.4,
  CUT_NEUTR_EM_EN_FRAC = 0.99,
  CUT_CHARG_EM_EN_FRAC = 0.99,
  CUT_NEUTR_HAD_EN_FRAC = 0.99,
  CUT_CHARG_HAD_EN_FRAC = 0.,
  CUT_1ST_BJET_PT = 50.,
//  CUT_2ND_BJET_PT = CUT_JET_PT,
//  CUT_BTAG_PT = CUT_JET_PT,
  CUT_BTAG_ETA = CUT_JET_ETA,
  CUT_BTAG_DISCRIM = 3.3+0.7,
  CUT_HT_LEPJET = 200.,
  CUT_HT_DILEP = 200.,
  CUT_MET_LEPJET = 20.,
  CUT_MET_DILEP = 15.,
  CUT_DELTA_R_LEPJET = 0.30,
  MATCH_DELTA_R_LEPJET = 0.30,
  MATCH_DELTA_R_BJET = 0.10,
  MATCH_DELTA_R_FLVR = 0.475,
  MATCH_VERTEX_DZ = 0.05,
  MATCH_DELTA_R_E_MU = 0.10,
  HT_MIN = 100.,
  HT_MAX = 1600.,  //1100.
  MT_MIN = 0.,
  MT_MAX = 250.,  //200.
  M3J_MIN = 50.,
  M3J_MAX = 1550.,   //960. & 60.
  MET_MAX = 360.,
  PHI_MAX = M_PI,
  PHI_MIN = -1.*PHI_MAX,
  DPHI_MAX = M_PI,
  LPT_MAX = 360.,  //250.
  JPT_MAX = 250.,
  BPT_MAX = 525.,
  ETA_MAX = 2.5,
  ETA_MIN = -1.*ETA_MAX,
  DETA_MAX = 2.*ETA_MAX,
  DELTAR_MAX = 6.,
  ISO_MAX = 0.25,
  VZ_MAX = 20.,
  VZ_MIN = -1.*VZ_MAX,
  DZ_MAX = 0.25,
  DISCRIM_MAX = 50.,
  DISCRIM_MIN = -25.,
  MASS_MIN = 10.,
  MASS_MAX = 510.;
string TriggerDoubleMu= "HLT_Mu17_TkMu8_v";
string TriggerDoubleEl= "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v";
string TriggerMuonElectron= "HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v";
string TriggerMuonElectron2= "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v";
string TriggerMuoniso="HLT_IsoMu24_eta2p1_v";

string TriggerSingleMuon="HLT_Mu40_eta2p1_v";

string TriggerMuoniso2="HLT_IsoMu24_v";
string TriggerSingleElectron="HLT_Ele27_WP80_v";
// these will change depending on the MC used

// for 2011 madgraph ZMuMu

     const double weights[]={
	2.18044,
	0.75958,
	0.932347,
	1.23022,
	1.01107,
	0.583834,
	1.09618,
	0.611194,
	0.791631,
0.933951,
0.710509,
0.670534,
0.74651,
0.833824,
1.05458,
2.9503,
2.55892,
2.29829,
1.5941,
	1.51827,
	1.41812,
	1.22111,
	1.03856,
	0.996799,
	0.926675,
	0.868843,
	0.838478,
	0.848388,
	0.851921,
0.826438,
0.816059,
0.835708,
0.813565,
0.845544,
0.910387,
1.02427,
1.10238,
1.09722,
1.01079,
	0.948719,
	1.07007,
	1.10456,
	1.00248,
	1.13232,
	1.17993,
	1.30207,
	1.43889,
	1.24978,
	1.18103,
1.60285,
1.35181,
1.26355,
1.30403,
1.77235,
1.32171,
2.06278,
1.48103,
1.04302,
1.50724,
	1.87661,
	10.0205,
	1.53306,
	1.28852,
	1.71814,
	0.783741,
	0.750644,
	1.18673,
	1.50129,
	1.40299,
0.869793,
1.69193,
2.31746,
4.28939,
1.99576,
3.72939,
0.,
3.20513,
2.76427,
1.01277,
	0.,
	0.,
	0.,
	0.786388,
	2.43066,
	0.474216,
	0.619579,
	1.97789,
	0.893623,
	0.,
1.44171};


/*
const double weights [50] = {
1.07652,
0.0430093,
0.0737551,
0.834337,
3.93636,
0.339041,
0.0329401,
0.0112451,
0.0104611,
0.0557284,
0.347516,
1.06351,
1.63179,
1.75244,
1.70012,
1.64782,
1.68977,
1.72174,
1.70717,
1.60452,
1.41793,
1.20383,
1.02947,
0.882494,
0.762143,
0.649683,
0.53972,
0.438161,
0.349207,
0.275196,
0.213219,
0.164608,
0.125573,
0.0965567,
0.0723858,
0.0538551,
0.0394291,
0.0286018,
0.0205512,
0.0147488,
0.0107557,
0.00767967,
0.00551422,
0.00382477,
0.002607,
0.00187566,
0.00124657,
0.000838946,
0.000566705,
0.00034506
};
*/
// this should be the good one
/*
const double weights [50] = {
  0.210543,
  0.260915,
  0.410351,
  0.293983,
  0.348181,
  0.548456,
  0.429515,
  0.435769,
  0.615016,
  0.911685,
  1.29756,
  1.65689,
  1.71087,
  1.51894,
  1.30673,
  1.13692,
  1.06478,
  1.03777,
  1.0563,
  1.10139,
  1.14746,
  1.17283,
  1.19067,
  1.18594,
  1.16977,
  1.1265,
  1.05372,
  0.961634,
  0.854111,
  0.736257,
  0.608833,
  0.489204,
  0.3792,
  0.289076,
  0.209037,
  0.145429,
  0.09628,
  0.0610844,
  0.0372639,
  0.0221993,
  0.01328,
  0.00779486,
  0.00469022,
  0.00283262,
  0.00177902,
  0.00126465,
  0.000894619,
  0.000687725,
  0.000564471,
  0.000439988

};
*/
//WJets
/*
const double weights [70] ={
0,
0,
0,
0,
0.108669,
0.439995,
0.513769,
0.476766,
0.615551,
1.01029,
1.31398,
1.65008,
1.66161,
1.62107,
1.26026,
1.17285,
1.12439,
1.07138,
1.05086,
1.11121,
1.16064,
1.09344,
1.22514,
1.21456,
1.13627,
1.12642,
1.04075,
0.913217,
0.90782,
0.728936,
0.636317,
0.482431,
0.395955,
0.263238,
0.20603,
0.14702,
0.0768348,
0.0595063,
0.0369977,
0.0203208,
0.0127694,
0.00680935,
0.00422929,
0.00243946,
0.0018024,
0.001134,
0.000946013,
0.000735173,
0.000597588,
0.000580238,
0.000358262,
0.000279389,
0.000175174,
0.000219788,
0.000274887,
8.54002e-05,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0

};
*/
//PH mm 2012
/*
const double weights[70]={
  0.232996,
0.188308,
0.28514,
0.335045,
0.312058,
0.525378,
0.409806,
0.422886,
0.590811,
0.882024,
1.2593,
1.60519,
1.66225,
1.49387,
1.2798,
1.1223,
1.04738,
1.02414,
1.04358,
1.09524,
1.13834,
1.17427,
1.19407,
1.18958,
1.17696,
1.13587,
1.06656,
0.974953,
0.868722,
0.751653,
0.625525,
0.508479,
0.393395,
0.301064,
0.216132,
0.151867,
0.100262,
0.0646608,
0.0389621,
0.0238347,
0.0143091,
0.00826623,
0.00495892,
0.00303034,
0.00197947,
0.00137222,
0.00100112,
0.000768244,
0.00062902,
0.000507785,
0.000431537,
0.000346621,
0.000339226,
0.000310266,
0.000237263,
0.00023888,
0.000214654,
0.000193391,
0.000181547,
0.000286683,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0
};
*/
//new DY file
/*
const double weights [70] ={
  0.253016,
  0.223964,
  0.345833,
  0.331658,
  0.282931,
  0.559827,
  0.418742,
  0.438717,
  0.60052,
  0.90624,
  1.30459,
  1.652,
  1.72323,
  1.52817,
  1.30963,
  1.13671,
  1.06103,
  1.03612,
  1.07005,
  1.10396,
  1.14369,
  1.17449,
  1.19561,
  1.19068,
  1.16651,
  1.13267,
  1.05971,
  0.964067,
  0.850157,
  0.739615,
  0.607219,
  0.488804,
  0.376253,
  0.290718,
  0.211292,
  0.146061,
  0.0951318,
  0.0607154,
  0.0371998,
  0.021883,
  0.0129877,
  0.00775003,
  0.00475692,
  0.00297684,
  0.00183914,
  0.00125111,
  0.000861426,
  0.000754962,
  0.000594288,
  0.000454728,
  0.000363417,
  0.000340092,
  0.000256547,
  0.000257508,
  0.000299594,
  0.000222348,
  0.000352027,
  0.000272738,
  0.000150207,
  0.000152482,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0
};
*/
//for ttbar
/*
const double weights [70] ={

0.0269952,
0.0358432,
0,
0.371551,
0.271683,
0.347377,
0.409936,
0.408317,
0.571691,
0.869095,
1.22879,
1.50614,
1.56116,
1.46834,
1.32692,
1.10557,
1.01358,
1.00059,
1.00235,
1.12025,
1.13414,
1.18642,
1.14954,
1.168,
1.21584,
1.13614,
1.06994,
0.979156,
0.895114,
0.76494,
0.642517,
0.530911,
0.403718,
0.314578,
0.233107,
0.158755,
0.107273,
0.0762124,
0.0447324,
0.0248056,
0.0139894,
0.00997956,
0.00550258,
0.00328401,
0.0019117,
0.00167198,
0.00122333,
0.000808718,
0.000664009,
0.000805914,
0.000447843,
0.000349249,
0.000159254,
0.000274745,
0.00022908,
0.000213508,
0.000131457,
0.000160047,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0
};
*/
//for WW
/*
const double weights[70] = {
0,
0.381452,
0.539932,
0.329511,
0.301178,
0.501717,
0.443817,
0.435142,
0.571672,
0.898414,
1.29985,
1.60998,
1.67085,
1.48369,
1.29085,
1.11663,
1.02958,
1.01307,
1.04105,
1.08806,
1.13889,
1.17112,
1.18441,
1.18678,
1.20948,
1.12153,
1.06566,
0.981956,
0.890417,
0.753286,
0.623932,
0.511752,
0.391036,
0.306247,
0.214427,
0.156375,
0.0997561,
0.0627456,
0.0404164,
0.0235465,
0.0144495,
0.00833542,
0.0049903,
0.00309074,
0.00188589,
0.00143973,
0.000934532,
0.000694077,
0.000648967,
0.000448781,
0.000372347,
0.000362613,
0.000245303,
0.000224915,
0.00020316,
0.000568049,
0.000349747,
0.000141938,
0.000255831,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0


};
*/
//for WZ
/*
const double weights[70] = {
  0.130163,
0,
0.978513,
0.447876,
0.327493,
0.623999,
0.459899,
0.455912,
0.607679,
0.953522,
1.28324,
1.61292,
1.72884,
1.53767,
1.30753,
1.13918,
1.04842,
1.04482,
1.04262,
1.11005,
1.1629,
1.16546,
1.18935,
1.19237,
1.17411,
1.1367,
1.05228,
0.969686,
0.869565,
0.733506,
0.623114,
0.499794,
0.390109,
0.27882,
0.20676,
0.140371,
0.0959233,
0.0597092,
0.0351317,
0.0222691,
0.0126518,
0.00775241,
0.0044295,
0.00275076,
0.00175828,
0.00128331,
0.000909883,
0.000744159,
0.000523907,
0.000431764,
0.0003374,
0.000364102,
0.000312839,
0.000331184,
0.000207104,
0.000171578,
0.000158461,
0.000154339,
0.00011591,
0.000274551,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0

};
*/
//for ZZ
/*
const double weights[70] = {
0.260077,
0.345321,
0.279308,
0.550707,
0.408976,
0.508698,
0.476565,
0.446181,
0.622815,
0.939315,
1.313,
1.69165,
1.74857,
1.5471,
1.33365,
1.14501,
1.07306,
1.04933,
1.06078,
1.10655,
1.15263,
1.16331,
1.17791,
1.17792,
1.16014,
1.10863,
1.05264,
0.961601,
0.855172,
0.723787,
0.606992,
0.491839,
0.377372,
0.281063,
0.202851,
0.142971,
0.0972845,
0.0594489,
0.0361231,
0.0225534,
0.0131023,
0.0076494,
0.00439161,
0.00285959,
0.00196059,
0.00121513,
0.000823589,
0.000647122,
0.000530643,
0.000415946,
0.000337079,
0.000309401,
0.000263706,
0.000264694,
0.000254654,
0.000228553,
0.00021108,
0.00038548,
0.000154399,
0.000548579,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0

};
*/
//for WJets, multiply below by ttbar
/*
const double WJetsCorrection [50] = {1.24207,
				     0.84968,
				     1.01623,
				     0.985822,
				     0.992387,
				     0.981417,
				     1.04359,
				     1.01589,
				     0.948213,
				     1.01962,
				     0.970325,
				     0.963572,
				     0.957001,
				     1.00634,
				     1.0132,
				     0.999944,
				     0.874279,
				     1.08339,
				     1.06702,
				     0.96785,
				     0.916341,
				     1.14305,
				     0.981721,
				     0.967427,
				     0.911003,
				     0.963762,
				     0.919244,
				     1.17332,
				     0.884276,
				     1.22269,
				     0.906702,
				     1.24275,
				     0.436754,
				     1.89711,
				     1.10911,
				     0,
				     1.0363,
				     0,
				     0,
				     0.641599,
				     0,
				     0,
				     0,
				     0,
				     0,
				     0,
				     0,
				     0,
				     0,
0};
*/
//for 2011 powheg TupleMaker
/*
const double weights [50] = {1.25501,
    1.07428,
    1.25551,
    1.15666,
    1.22246,
    1.15657,
    1.24121,
    1.28427,
    1.23065,
    1.34526,
    1.38593,
    1.50142,
    1.26444,
    1.08908,
    0.939453,
    0.839133,
    0.875633,
    0.568408,
    0.523315,
    0.395852,
    0.299744,
    0.220639,
    0.174965,
    0.1181,
    0.0872553,
    0.0686035,
    0.0534855,
    0.0667177,
    0.032515,
    0.018864,
    0.0125302,
    0.0139801,
    0.011766,
    0.0297955,
    0.0049456,
    0.00724493,
    0.00173486,
    0.0382712,
    0,
    0.000172399,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0};
*/
//CUT_BTAG_DISCRIM options:
//TCHE: 3.30 ? [djs fudge factor 0.7]
//JP:   0.70 ?

////////////////////////////////////////////////////////////////

//jet collections:  ak5CaloJets & ak5PFJets
//old btag collections:  trackCountingHighEffBJetTags, jetProbabilityBJetTags
//new btag collections:  newJetProbabilityBJetTags, newTrackCountingHighEffBJetTags
//flavor collections:  JetbyValAlgo, GenJetbyValAlgo, GenJetbyValPhys

const char* MET_COLLECTION = "pfMet";
//const char* BEAMSPOT_COLLECTION = "offlineBeamSpot";
const char* VERTEX_COLLECTION = "offlinePrimaryVertices";
const char* ELEC_COLLECTION = "gsfElectrons";
const char* MUON_COLLECTION = "muons";
//const char* FLAVOUR_COLLECTION = "GenJetbyValAlgo";
const char* JET_COLLECTION = "ak5PFJetsL1FastL2L3";
//const char* JET_COLLECTION = "ak5PFJets";
//const char* BTAG_COLLECTION = "newTrackCountingHighEffBJetTags";
//const char* CORR_COLLECTION = "AK5PF";
//const char* JET_CORR_L1 = "ak5PFL1Fastjet";
//const char* JET_CORR_L2 = "ak5PFL2Relative";
//const char* JET_CORR_L3 = "ak5PFL3Absolute";
//const char* JET_CORR_MC = "ak5PFResidual";
const reco::PFCandidateCollection * pfCandidates;  

const char* 
  LINE = 
  "-----------------------------------------------------------------";

const math::XYZTLorentzVector ZERO_P4(0., 0., 0., 0.);

////////////////////////////////////////////////////////////////

int get_flvr_index(int);

char get_flavor(int);

TString get_true_false(bool);

double DeltaPhiX(double, double);
double DeltaPhiX(math::XYZVector, math::XYZVector);
double DeltaRX(double, double, double, double);
double DeltaRX(double, double);
double DeltaRX(math::XYZVector, math::XYZVector);
double get_mass(math::XYZTLorentzVector);
double get_rapid(TLorentzVector);
double acolinearity(math::XYZTLorentzVector, math::XYZTLorentzVector);
double get_pt(math::XYZVector);
double get_abs_p(math::XYZTLorentzVector);
double get_MT(double, double, double);
double get_MT(math::XYZVector, math::XYZVector);

////////////////////////////////////////////////////////////////

class TupleMaker : public EDAnalyzer {

public:
  explicit TupleMaker(const ParameterSet&);
  ~TupleMaker();
  static void fillDescriptions(ConfigurationDescriptions& descriptions);
  
  Handle<BeamSpot> beam_spots;
  Handle<PFMETCollection> MET;
  //Handle<GenMETCollection> metHandle;
  Handle<vector<GsfElectron> > electrons;
  Handle<vector<Muon> > muons;
  Handle<PFJetCollection> jets;
  Handle<BeamSpot> bsHandle;
  //  Handle<ConversionCollection> hConversions;

  


Handle<reco::VertexCollection> vtx_h;
Handle<double> rhoIso_h;
  // Handle< GenEventInfoProduct > GenInfoHandle;
  // Handle< HepMCProduct > EvtHandle ;

  //comment out for data
  //  Handle<GenParticleCollection> genParticles;
  //  Handle<JetTagCollection> btags;
   Handle<VertexCollection> vertices;
  //  Handle<JetFlavourMatchingCollection> flavours;
  
  //  ESHandle<JetCorrectorParametersCollection> jet_corrections;

private:
  //cms
  virtual void beginJob();
  virtual void analyze(const Event&, const EventSetup&);
  virtual void endJob();
  
  virtual void beginRun(Run const&, EventSetup const&);
  virtual void endRun(Run const&, EventSetup const&);
  virtual void beginLuminosityBlock(LuminosityBlock const&, EventSetup const&);
  virtual void endLuminosityBlock(LuminosityBlock const&, EventSetup const&);
  bool trainTrigPresel(const reco::GsfElectron& ele);
  //  virtual void evaluate_mvas(const Event&, const  EventSetup&);
  //djs
  edm::InputTag inputTagPhotons_;
  edm::InputTag inputTagPFCandidateMap_;
  edm::InputTag triggerResultsTag_;
  edm::InputTag triggerEventTag_; 
 std::vector<edm::InputTag> inputTagIsoDepElectrons_;
  std::vector<edm::InputTag> inputTagIsoDepPhotons_;
  //  std::vector<edm::InputTag> inputTagIsoValElectronsNoPFId_;
  std::vector<edm::InputTag> inputTagIsoValElectronsPFId_;   
  std::vector<edm::InputTag> inputTagIsoValPhotonsPFId_;   
  //MVA
  ParameterSet conf_;

  EGammaMvaEleEstimator* myMVATrigV0;
  EGammaMvaEleEstimator* myMVATrigNoIPV0;
  EGammaMvaEleEstimator* myMVANonTrigV0;
  
  TMVA::Reader             *myTMVAReader;
  Float_t                   myMVAVar_fbrem;
  Float_t                   myMVAVar_kfchi2;
  Float_t                   myMVAVar_kfhits;
  Float_t                   myMVAVar_gsfchi2;
  
  Float_t                   myMVAVar_deta;
  Float_t                   myMVAVar_dphi;
  Float_t                   myMVAVar_detacalo;
  Float_t                   myMVAVar_dphicalo;

  Float_t                   myMVAVar_see;
  Float_t                   myMVAVar_spp;
  Float_t                   myMVAVar_etawidth;
  Float_t                   myMVAVar_phiwidth;
  Float_t                   myMVAVar_e1x5e5x5;
  Float_t                   myMVAVar_R9;
  Float_t                   myMVAVar_nbrems;

  Float_t                   myMVAVar_HoE;
  Float_t                   myMVAVar_EoP;
  Float_t                   myMVAVar_IoEmIoP;
  Float_t                   myMVAVar_eleEoPout;
  Float_t                   myMVAVar_PreShowerOverRaw;
  Float_t                   myMVAVar_EoPout;

  Float_t                   myMVAVar_d0;
  Float_t                   myMVAVar_ip3d;

  Float_t                   myMVAVar_eta;
  Float_t                   myMVAVar_pt;

  double _Rho;
  unsigned int ev;




  void setup_stopwatch();
  void setup_names();
  void setup_counters();
  void setup_histograms(const ParameterSet&);
  void setup_flags();

  
  void time_stamp();
  
  bool cut_electron(vector<GsfElectron>::const_iterator);
  bool cut_tight_electron(int);
  bool cut_muon(vector<Muon>::const_iterator);
  bool cut_tight_muon(vector<Muon>::const_iterator);
  bool cut_mu_electron(double, double);
  bool cut_btag(JetTagCollection::const_iterator);
  bool cut_bjets();
  bool cut_jet(PFJet);
  bool cut_jet_vtx(PFJet);
  bool cut_Z_window();
  bool cut_2or3_jets();
  bool cut_3_bjets();
  bool cut_same_sign();
  bool cut_MET_dilep();
  bool cut_MET_lepjet();
  bool cut_HT_dilep();
  bool cut_HT_lepjet();
  bool lepton_jet(PFJetCollection::const_iterator);
  bool noniso_jet(PFJet);
  double deltaR_ljet(PFJet);
  double deltaR_ejet(PFJet, int);
  double deltaR_mjet(PFJet, int);
  double leptonVertex(reco::TrackRef, reco::TrackRef, const EventSetup&);  
  double leptonVertex(reco::GsfTrackRef, reco::GsfTrackRef, const EventSetup&);  
  double leptonVertex(reco::TrackRef, reco::GsfTrackRef, const EventSetup&);  
 
 void veto_event();
  
  void fill(TH1D*, int);
  void fill(TH1D*, double);
  void set_h_bin(TH1D*, int, int);
  void set_h_bin(TH1D*, int, double);

  // Nick UPDATE
  void FinalCuts(double&, double&, double&); 
  void NickNPlots(double&, double&, double&);
  // end
  // void MC_info(double&, double&, double&, const Event&);
  void MC_tau(const Event&, double, double);
  void order_elecs(int, vector<GsfElectron>::const_iterator);
  void order_muons(int, vector<Muon>::const_iterator);
  void order_jets(int, PFJet);
  //void order_btags(int, JetTagCollection::const_iterator);
  //void jet_addons(int);
  
  void update_MET();
  void set_sign();
  
  PFJet correct_jet(PFJetCollection::const_iterator, const Event&, const EventSetup&);

  void match_to_jet(JetTagCollection::const_iterator);  
  //void match_to_bjet(JetTagCollection::const_iterator);  //old method
  bool match_to_vtx(PFJet);
  
  void load_prep();
  void load_beam(const Event&);
  void load_rMET(const Event&);
  void load_leps(const Event&);
  int  load_elec(const Event&, const EventSetup&, int num_m);
  int  load_muon(const Event&);
  void load_vtxs(const Event&);
  void load_corr(const Event&, const EventSetup&);
  void load_jets(const Event&, const EventSetup&);
  void load_btag(const Event&);
  void load_flvr(const Event&);
  void load_triggers(const Event&);
  virtual bool triggerDecision(edm::Handle<edm::TriggerResults>& hltR, int iTrigger);
      void load_PhiStar(const Event&);
  void load_Wacceptance(const Event&, const EventSetup&);

  void load_acceptance(const Event&);
  void load_dileptons();
  void load_leptonjet();
  void save_counters();

  //comment out for data
           void load_pileup(const Event&);

virtual TLorentzVector GetReducedMET(TLorentzVector sumJet, TLorentzVector lep1, TLorentzVector lep2, TLorentzVector metP4, int version);
  void load_NickAnalysis();
  void lepton_search(const Event&, const EventSetup&);

  void fill_pre_ll();
  void fill_pre_lj();
  void fill_pre_jets();
  void fill_pre_bjet();   //void fill_pre_bjet(int j, JetTagCollection::const_iterator i_tag);
  void fill_1b_ll();
  void fill_2b_ll();
  void fill_1b_lj();
  void fill_2b_lj();
  void fill_1b_ll_jets();
  void fill_2b_ll_jets();
  void fill_1b_lj_jets();
  void fill_2b_lj_jets();
  
  double get_MET();
  double get_METphi();
  double get_METsig();
  double get_dphi_ll();
  double get_dphi_lv();
  double get_MT_lv();
  double get_HT();
  double get_H();
  double get_centrality();
  double get_sphericity();
  double get_aplanarity();  
  double get_dphi_vj();
  double get_dphi_lv_jj();
  double get_dphi_lvb_bjj();
  double get_dphi_lv_b();
  double get_e_iso(vector<GsfElectron>::const_iterator);
  double get_m_iso(vector<Muon>::const_iterator);
  double get_sum_track_pt(PFJet);
  double get_lj_M_bjj();
  double get_MT_lvb();
  double get_dR_bj();
  double get_lep_vz();
  double get_weight();
 
  vector<math::XYZVector> get_evs_vector();
  
  math::XYZVector get_lep_vector(int which = 1);

  void get_delta_bb(double&, double&, double&);
  void get_delta_jj(double&, double&, double&);
  void get_delta_b_bjj(double&, double&, double&);
  void get_delta_l_jj(double&, double&, double&);
  void get_delta_lb_bjj(double&, double&, double&);  
  
  void get_S_and_A(double&, double&);
  
  //void print_debug(const char*, int n_check=10, int var = (int)(JUNK));
  void print_debug(const char*, int n_check=N_DEBUG, double var = JUNK);
  
  //varibles
  //const Event event_;
  //const EventSetup setup_;
  int nPUVertices;
  int nPUVerticesTrue;
  //int num_btags;
  int num_bjets;
  int num_jets;
  
  int muCounter;

  double isquark, motherId1, motherId2;

  int num_events;  
  int num_trilep;
  int num_no_lep;
  
  int num_2jet_veto;
  int num_bjet_veto;
  int num_isoj_veto;
  int num_lldz_veto;
  int num_Z_ll_veto;
  int num_3jet_veto;
  int num_mET1_veto;
  int num_mET2_veto;
  int num_HT1_veto;
  int num_HT2_veto;
  
  int num_pass_ee;
  int num_pass_em;
  int num_pass_mm;
  int num_pass_ej;
  int num_pass_mj;
  
  int num_good_ee;
  int num_good_em;
  int num_good_mm;
  int num_good_ej;
  int num_good_mj;
  
  int num_ss_ee;
  int num_ss_mm;
  int num_ss_em;
  int num_3b_xx;
  
  int num_2b_ee;
  int num_2b_mm;
  int num_2b_em;
  int num_2b_ej;
  int num_2b_mj;  
  
  int num_miss_b;
  int num_unkn_b;
  int num_true_b;
  int num_falseb;
  int num_unkn_j;
  
  int num_uds_g;
  int num_unfound;
  int num_0flavor;

  double MCrapid;
  int numberjets;
  int num_m, num_e;
  int nummuon;
  int eventid, runid;
  double ZpT;
  double InvMass;
  double weight;
  double mother, mother2, upz, downpz;
  int extraMuons;
  int extraElectrons;
  int vtxCount;
  double lep1pt, lep2pt;
  double_t meta1, meta2;
  double normChi2, d_xy, m_iso;
  int numPH, numVH;
  int jflavor[N_JETS];
  double lep1dxy, lep2dxy;
  double METSignificance;
  double lep1Q, lep2Q;
  double elec1Q, elec2Q;
  double liso1, liso2, eiso1, eiso2;
  int num_mm, num_ee, num_me;
  double MCpTdiff;
  double lep1ptE, lep2ptE;
  double mujetd1, mujetd2;

  // these are variables for the W study with isolation 8/2/14
  double Npt, Neta, Nphi, Nisolation, Mpt, Mphi, Meta, Njetd, Wpass, Mq, WqT;

  bool is_bjet[N_JETS];
  bool DoubleMuonTrigger;
  bool DoubleElectronTrigger;
  bool MuonElectronTrigger;
  bool SingleMuonTrigger40;
  bool SingleMuonTriggeriso;
  bool SingleMuonTriggeriso2;
  bool SingleElectronTrigger;
  bool lep1ID;
  bool lep2ID;
  bool elec1ID;
  bool elec2ID;
  int istautau;
  bool passacc;
  int eff1;
  int eff2;
  int eff3;
  int eff4;
  int eleceff3;
  int eleceff4;
  int muoneff3;
  int muoneff4;
  int aftertrigger;
   vector<string> hlNames;
  char process_;
  char dilep_;
  char lepjet_;
  char sign_;
  bool acceptancemm;
  bool acceptanceee;
  bool acceptanceem;
  double normChi21, normChi22;
  int numPH21,numPH22, hitPattern1, hitPattern2, numberStations1, numberStations2, hits1, hits2;

  TString Name_h[N_HIST];
  TString Title_h[N_HIST];
  string hlTriggerResults_; 
  double mET_x;
  double mET_y;
  double mET;
  double mETphi;
  double elec1dxy, elec2dxy;
  //double bdiscrim_1;
  //double bdiscrim_2;
  PFIsolationEstimator eleIsolator;
  double dReEm1[3];
  double dReEm2[3];

  float vertexProb;
  reco::GsfTrackRef elec1track, elec2track;
  double acol;
  double muon1iso;
  double muon2iso;
  double elecReg1, elecReg2;
  double lepton1eta, lepton2eta;
  double lepton1phi, lepton2phi;
  double elec1eta, elec1phi, elec2phi, elec2eta;
  double redMETtotal;
  double rMETparallel;
  double rMETperp;
  double jet_pt[N_JETS];
  double Emuon_pt[3];
  double Emuon_eta[3];
  double Eelectron_pt[3];
  double Eelectron_eta[3];

  double MCinvmass, MC1phi, MC1eta, MC1pt, MC2phi, MC2eta, MC2pt, MC1Q, MC2Q, MCN1Q, MCN2Q, MCB1Q, MCB2Q, MCD1Q, MCD2Q;
  int MC1a, MC2a;
  double MC1truth, MC2truth;
  int DiLeptonType;

  double jdiscrim[N_JETS];
  double jet_vz[N_JETS];

  double elec1pt, elec2pt;
  math::XYZVector vec_mET;
  math::XYZPoint beam_xyz;
  math::XYZTLorentzVector jflvr_p[N_JETS];
  TLorentzVector rMET;
  TLorentzVector sumJet;
  TLorentzVector lep1;
  TLorentzVector lep2;
  TLorentzVector* metP4;
  TLorentzVector* l14vector;
  TLorentzVector* l24vector;
  TLorentzVector* MC14vector;
  TLorentzVector* MC24vector;
  TLorentzVector* MCN14vector;
  TLorentzVector* MCN24vector;
  TLorentzVector* MCB14vector;
  TLorentzVector* MCB24vector;
  TLorentzVector* MCD14vector;
  TLorentzVector* MCD24vector;


  math::XYZTLorentzVector elec1P4;
  math::XYZTLorentzVector elec2P4;
  TLorentzVector redMET;
  const JetCorrector* corr_L1;
  const JetCorrector* corr_L2;     
  const JetCorrector* corr_L3;
  const JetCorrector* corr_MC; 
  VertexCollection::const_iterator lep_vtx;
  VertexCollection::const_iterator jet_vtx;
  //VertexCollection::const_iterator jvtx_[N_JETS];
  vector<GsfElectron>::const_iterator elec_1;
  vector<GsfElectron>::const_iterator elec_2;
  vector<Muon>::const_iterator muon_1;
  vector<Muon>::const_iterator muon_2;
  // Handle<GenParticleCollection> genParticles;
  //  const GenParticle & MCmuon1;
  // const GenParticle & MCmuon2;
  PFJet bjet_1;         //PFJetCollection::const_iterator
  PFJet bjet_2;         //PFJetCollection::const_iterator
  PFJet jet_[N_JETS];   //PFJetCollection::const_iterator jet_[N_JETS];
  //JetTagCollection::const_iterator btag_[N_TAGS];
  
  //RefToBase<reco::Jet> jref_[N_JETS];
  
  TStopwatch* watch;
  TTree* tree;
  TNtuple* ntuple;
  TH1D* h_base;
  TH1D* h_passacc;
  TH1D* h_passrecoeff;
  TH1D* h_passseleceff;
  TH1D* h_num_events;
  TH1D* h_num_leps;
  TH1D* h_num_jets;
  TH1D* h_num_bets;
  TH1D* h_num_pass_ee;
  TH1D* h_num_pass_mm;
  TH1D* h_num_pass_em;
  TH1D* h_num_pass_ej;  
  TH1D* h_num_pass_mj;
  TH1D* h_num_good_ee;
  TH1D* h_num_good_mm;
  TH1D* h_num_good_em;
  TH1D* h_num_good_ej;  
  TH1D* h_num_good_mj;
  
  TH1D* h_MCM1pt;
  TH1D* h_MCM2pt;
  TH1D* h_MCM1eta;
  TH1D* h_MCM2eta;
  TH1D* h_MCMIM;
  TH1D* h_MCMrapid;
  TH1D* h_MCE1pt;
  TH1D* h_MCE2pt;
  TH1D* h_MCE1eta;
  TH1D* h_MCE2eta;
  TH1D* h_MCEIM;
  TH1D* h_MCErapid;

  TH1D* h_MCS1pt;
  TH1D* h_MCS2pt;
  TH1D* h_MCS1eta;
  TH1D* h_MCS2eta;
  TH1D* h_MCSIM;
  TH1D* h_MCSrapid;

  TH1D* h_MCmu1dxy;
  TH1D* h_MCmu1dz;
  TH1D* h_MCmu1Global;
  TH1D* h_MCmu1PF;
  TH1D* h_MCmu1Chi2;
  TH1D* h_MCmu1hits;
  TH1D* h_MCmu1Stations;
  TH1D* h_MCmu1Phits;
  TH1D* h_MCmu1Thits;
  TH1D* h_MCmu1iso;

  TH1D* h_MCmu2dxy;
  TH1D* h_MCmu2dz;
  TH1D* h_MCmu2Global;
  TH1D* h_MCmu2PF;
  TH1D* h_MCmu2Chi2;
  TH1D* h_MCmu2hits;
  TH1D* h_MCmu2Stations;
  TH1D* h_MCmu2Phits;
  TH1D* h_MCmu2Thits;
  TH1D* h_MCmu2iso;

  TH1D* h_uppz;
  TH1D* h_uppx;
  TH1D* h_uppy;

  TH1D* h_downpz;
  TH1D* h_downpx;
  TH1D* h_downpy;


  TH1D* h_num_uds_g;
  TH1D* h_num_unfound;
  TH1D* h_num_0flavor;
  TH1D* h_num_HT1_veto;
  TH1D* h_num_HT2_veto;
  
  TH1D* h_num_2jet_veto;
  TH1D* h_num_bjet_veto;
  TH1D* h_num_isoj_veto;
  TH1D* h_num_lldz_veto;
  TH1D* h_num_Z_ll_veto;
  TH1D* h_num_3jet_veto;
  //
  TH1D* h_num_mET1_veto;
  TH1D* h_num_mET2_veto;
  TH1D* h_num_ss_ee;
  TH1D* h_num_ss_mm;
  TH1D* h_num_ss_em;
  TH1D* h_num_3b_xx;
  TH1D* h_num_2b_ee;
  TH1D* h_num_2b_mm;
  TH1D* h_num_2b_em;
  TH1D* h_num_2b_ej;
  TH1D* h_num_2b_mj;
  //
  TH1D* h_num_miss_b;
  TH1D* h_num_unkn_b;
  TH1D* h_num_true_b;
  TH1D* h_num_falseb;
  TH1D* h_num_unkn_j;

  TH1D* h_pre_b_jpt;
  TH1D* h_pre_b_fpt;
  TH1D* h_pre_missb_jpt;
  TH1D* h_pre_missb_fpt;
  TH1D* h_pre_b_dptfj;
  TH1D* h_pre_b_rptfj;
  TH1D* h_pre_b_fmass;

  TH1D* jet_cut_pt;

  TH1D* h_jet_pt;
  TH1D* h_jet_pt_50;
  TH1D* h_jet_pt_100;
  TH1D* h_jet_pt_200;
  TH1D* h_jet_pt_300;
  TH1D* h_jet_pt_500;
  TH1D* h_jet_pt_700;
  
  TH1D* h_jet_lead_pt;
  TH1D* h_jet_lead_pt_50;
  TH1D* h_jet_lead_pt_100;
  TH1D* h_jet_lead_pt_200;
  TH1D* h_jet_lead_pt_300;
  TH1D* h_jet_lead_pt_500;
  TH1D* h_jet_lead_pt_700;

  TH1D* h_ZpT_log;
  TH1D* h_ZpT;
  TH1D* h_ZpT_50;
  TH1D* h_ZpT_100;
  TH1D* h_ZpT_200;
  TH1D* h_ZpT_300;
  TH1D* h_ZpT_500;
  TH1D* h_ZpT_700;

  TH1D* h_InvMass;

  TH1D* h_numJets_50;
  TH1D* h_numJets_100;
  TH1D* h_numJets_200;
  TH1D* h_numJets_300;
  TH1D* h_numJets_500;
  TH1D* h_numJets_700;
  TH1D* h_numJets;

  TH1D* h_mET;
  TH1D* h_mET_50;
  TH1D* h_mET_100;
  TH1D* h_mET_200;
  TH1D* h_mET_300;
  TH1D* h_mET_500;
  TH1D* h_mET_700;

  TH1D* h_extra_electrons;
  TH1D* h_extra_muons;
  TH1D* h_extra_electrons50;
  TH1D* h_extra_muons50;
  TH1D* h_extra_electrons100;
  TH1D* h_extra_muons100;
  TH1D* h_extra_electrons200;
  TH1D* h_extra_muons200;
  TH1D* h_extra_electrons300;
  TH1D* h_extra_muons300;
  TH1D* h_extra_electrons500;
  TH1D* h_extra_muons500;
  TH1D* h_extra_electrons700;
  TH1D* h_extra_muons700;

  TH1D* h_nPUVertices;
  TH1D* h_nPUVerticesTrue;
  TH1D* h_vtxCount;

  TH1F* h_mva_nonTrig;
};

////////////////////////////////////////////////////////////////
