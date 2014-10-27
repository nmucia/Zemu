//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Aug 29 11:48:53 2014 by ROOT version 5.32/00
// from TTree T/T
// found on file: /uscms/home/nmucia2/nobackup/825/SMA.root
//////////////////////////////////////////////////////////

#ifndef AnalyzerZemu_h
#define AnalyzerZemu_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include "TMath.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH1.h"
#include "THStack.h"
#include "TChain.h"
#include "TMath.h"
#include <TBranch.h>
#include <TROOT.h>
#include <TRint.h>
#include <vector>
#include <math.h>
#include <cmath>
#include <string>
#include "TNtuple.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include <rochcor2012jan22.C>

#include <TRandom3.h>
#include <fstream>
#include <set>
#include <algorithm>
#include <sstream>
#include <utility>
// Header file for the classes stored in the TTree if any.
#include <TLorentzVector.h>

const double
  PTCUT1 = 20,
  PTCUT2 = 30,
  ETACUTMUON = 2.4,
  ETACUTELECTRON = 2.4;
      float qter = 1.0;  

//MG Ful 2012

   const double jet_weight[]={1.06,0.935,0.888,0.874,0.82,0.91,0.90,0.98,1.,1.,1.};


    const double weights[]={
      0.291545,
      0.227707,
      0.335468,
      0.302465,
      0.297379,
      0.522843,
      0.415166,
      0.420852,
      0.594346,
      0.885749,
      1.26379,
      1.60555,
      1.66451,
      1.48817,
      1.28037,
      1.12133,
      1.04736,
      1.02398,
      1.04326,
      1.09364,
      1.13989,
      1.17122,
      1.18878,
      1.19121,
      1.17746,
      1.13489,
      1.06959,
      0.979838,
      0.872769,
      0.751919,
      0.626006,
      0.505877,
      0.393746,
      0.300135,
      0.217398,
      0.150956,
      0.100377,
      0.064226,
      0.0392567,
      0.0238096,
      0.014123,
      0.00832263,
      0.0049879,
      0.00305456,
      0.00197542,
      0.0013562,
      0.000976246,
      0.000750864,
      0.000606018,
      0.000491019,
      0.00043391,
      0.000375931,
      0.0003277,
      0.000306953,
      0.000284736,
      0.000246177,
      0.000211898,
      0.000204152,
      0.000158951,
      0.000341641,
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


// PH mm Full 2012
/*
    const double weights[]={
      0.203918,
0.280091,
0.38974,
0.331493,
0.298871,
0.520424,
0.411246,
0.419909,
0.591984,
0.89478,
1.27703,
1.61219,
1.66609,
1.50169,
1.28773,
1.1272,
1.04556,
1.02812,
1.04819,
1.09716,
1.14222,
1.17487,
1.1923,
1.19424,
1.17847,
1.1297,
1.06873,
0.977302,
0.868024,
0.746,
0.622732,
0.504105,
0.392784,
0.298446,
0.213667,
0.149793,
0.0996344,
0.0627497,
0.0389887,
0.0232458,
0.0138045,
0.00809082,
0.00485627,
0.00298571,
0.0019049,
0.00131002,
0.000948511,
0.000722666,
0.000561008,
0.000491535,
0.000403371,
0.000347891,
0.000298485,
0.000287582,
0.000240339,
0.000219928,
0.000222314,
0.000164859,
0.000175731,
0.000268827,
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
//PHtt Full 2012
/*
    const double weights[]={
      0,
      0,
      0,
      0.0950039,
      0.138936,
      1.12509,
      0.394119,
      0.405007,
      0.52372,
      0.829041,
      1.18514,
      1.60989,
      1.59931,
      1.32782,
      1.14246,
      1.1097,
      0.965382,
      0.986245,
      1.00461,
      1.06121,
      1.10876,
      1.10133,
      1.18908,
      1.22836,
      1.22658,
      1.13908,
      1.14918,
      1.09031,
      0.885819,
      0.866404,
      0.674173,
      0.550461,
      0.419268,
      0.31741,
      0.289055,
      0.174235,
      0.113551,
      0.0735609,
      0.05367,
      0.02649,
      0.0163259,
      0.00954841,
      0.00574519,
      0.00304637,
      0.00238975,
      0.00175508,
      0.00151187,
      0.00103393,
      0.000679137,
      0.000741847,
      0.000572557,
      0.000476274,
      0.000447927,
      0.000562008,
      0,
      0.000218372,
      6.72257e-05,
      0,
      4.91738e-05,
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
//ttbar
/*
const double weights[]={
0.137415,
0.121636,
0.258258,
0.630441,
0.329276,
0.479959,
0.445799,
0.409461,
0.58484,
0.917204,
1.24911,
1.61568,
1.62394,
1.52692,
1.30098,
1.12894,
1.04155,
1.00775,
1.05177,
1.08248,
1.13499,
1.16518,
1.21073,
1.19456,
1.17526,
1.12911,
1.06196,
0.977012,
0.875399,
0.7461,
0.628265,
0.510214,
0.40089,
0.295261,
0.224149,
0.147102,
0.0973747,
0.0664454,
0.0401442,
0.0247427,
0.013659,
0.00832307,
0.00491053,
0.00294335,
0.00185624,
0.00141246,
0.000940566,
0.000762343,
0.000573969,
0.000576897,
0.000356199,
0.000302604,
0.000228648,
0.000174818,
0.00023322,
0.000181138,
0.000148702,
9.05215e-05,
0.000163157,
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
//WW
/*
    const double weights[]={
      0,
      0.199772,
      0.377027,
      0.31859,
      0.302844,
      0.510915,
      0.429537,
      0.404236,
      0.570785,
      0.884785,
      1.28372,
      1.57436,
      1.66374,
      1.46586,
      1.28369,
      1.11761,
      1.02756,
      1.01015,
      1.04675,
      1.08372,
      1.13925,
      1.1707,
      1.18621,
      1.19017,
      1.19504,
      1.12823,
      1.0563,
      0.976784,
      0.893089,
      0.761551,
      0.638589,
      0.516794,
      0.397877,
      0.312445,
      0.214385,
      0.157585,
      0.101942,
      0.0655436,
      0.0413829,
      0.0236911,
      0.0148912,
      0.00834674,
      0.00519986,
      0.00303758,
      0.00185546,
      0.00150179,
      0.000964531,
      0.000713192,
      0.000686755,
      0.000481259,
      0.000396198,
      0.000458008,
      0.000250349,
      0.000226857,
      0.000255355,
      0.000297496,
      0.000293068,
      0.000178404,
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

//WZ
/*
	const double weights[]={
	  0.222879,
0,
0.418879,
0.383451,
0.287574,
0.518976,
0.420828,
0.43071,
0.606598,
0.917725,
1.25825,
1.60116,
1.68683,
1.49599,
1.28237,
1.12865,
1.02961,
1.02369,
1.03271,
1.10677,
1.12371,
1.17705,
1.20257,
1.20393,
1.1676,
1.14101,
1.06661,
0.973378,
0.871803,
0.759375,
0.627964,
0.506604,
0.394131,
0.291779,
0.213238,
0.147696,
0.101908,
0.0646001,
0.0378155,
0.0237596,
0.01391,
0.00839865,
0.0048184,
0.00311466,
0.00191139,
0.00136989,
0.000973103,
0.000806397,
0.00057372,
0.000499036,
0.000366089,
0.000397721,
0.000336355,
0.000362937,
0.000246697,
0.000220346,
0.000197334,
0.000220231,
0.000158779,
0.000470117,
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
//ZZ
/*
	const double weights[]={
	  0.201523,
	  0.356767,
	  0.27545,
	  0.34671,
	  0.38267,
	  0.513241,
	  0.4329,
	  0.434912,
	  0.590702,
	  0.903646,
	  1.2679,
	  1.61852,
	  1.65907,
	  1.48784,
	  1.29171,
	  1.12613,
	  1.05215,
	  1.02827,
	  1.05083,
	  1.09413,
	  1.15072,
	  1.16145,
	  1.17603,
	  1.1804,
	  1.16907,
	  1.13198,
	  1.07001,
	  0.977125,
	  0.863408,
	  0.74564,
	  0.626435,
	  0.504273,
	  0.39217,
	  0.298674,
	  0.215021,
	  0.150707,
	  0.102891,
	  0.0643313,
	  0.0381549,
	  0.0243249,
	  0.0142938,
	  0.00806898,
	  0.00487656,
	  0.00306445,
	  0.0020701,
	  0.00135217,
	  0.000932531,
	  0.000749032,
	  0.000573792,
	  0.000458868,
	  0.000407709,
	  0.000359613,
	  0.000315119,
	  0.000321727,
	  0.000213765,
	  0.000335551,
	  0.000206598,
	  0.000298693,
	  0.000205093,
	  0.000425072,
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
/*
const double weights[]={
	0.749946,
	0.604585,
	0.644087,
	0.336204,
	0.222737,
	0.307487,
	0.294162,
	0.423543,
	0.704628,
	1.05018,
	1.4612,
	1.81815,
	1.83947,
	1.59679,
	1.34102,
	1.16494,
	1.09254,
	1.07301,
	1.09734,
	1.14941,
	1.1865,
	1.19629,
	1.18603,
	1.15932,
	1.1147,
	1.03965,
	0.943131,
	0.829977,
	0.711761,
	0.592483,
	0.477177,
	0.371968,
	0.277421,
	0.200643,
	0.136119,
	0.0870989,
	0.0523611,
	0.0296754,
	0.0157536,
	0.00816188,
	0.00408293,
	0.0020068,
	0.000988836,
	0.0004862,
	0.00024331,
	0.000123126,
	6.19031e-05,
	3.15919e-05,
	1.62003e-05,
	8.06197e-06,
	4.26594e-06,
	2.17158e-06,
	1.09619e-06,
	5.87595e-07,
	3.08713e-07,
	1.49736e-07,
	7.16594e-08,
	3.80626e-08,
	1.62112e-08,
	1.89263e-08,
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

/*
const double weights[]={
	0.749946,
	0.604585,
	0.644087,
	0.336204,
	0.222737,
	0.307487,
	0.294162,
	0.423543,
	0.704628,
	1.05018,
	1.4612,
	1.81815,
	1.83947,
	1.59679,
	1.34102,
	1.16494,
	1.09254,
	1.07301,
	1.09734,
	1.14941,
	1.1865,
	1.19629,
	1.18603,
	1.15932,
	1.1147,
	1.03965,
	0.943131,
	0.829977,
	0.711761,
	0.592483,
	0.477177,
	0.371968,
	0.277421,
	0.200643,
	0.136119,
	0.0870989,
	0.0523611,
	0.0296754,
	0.0157536,
	0.00816188,
	0.00408293,
	0.0020068,
	0.000988836,
	0.0004862,
	0.00024331,
	0.000123126,
	6.19031e-05,
	3.15919e-05,
	1.62003e-05,
	8.06197e-06,
	4.26594e-06,
	2.17158e-06,
	1.09619e-06,
	5.87595e-07,
	3.08713e-07,
	1.49736e-07,
	7.16594e-08,
	3.80626e-08,
	1.62112e-08,
	1.89263e-08,
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
// Fixed size dimensions of array or collections stored in the TTree if any.

class AnalyzerZemu : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   bool IS_DATA;
   double lep1eta, lep2eta, lep1phi, lep2phi;
   // double MC1eta, MC1phi, MC1pt, MC2phi, MC2eta, MC2pt, MC1Q, MC2Q;
   double leadeta;
   int numa, numb;
   int filetype;
   int produceForest;
   bool notdoublecounted;
 vector<pair<long, long> > oldevents;
   std::pair <long, long> foo;
 ostringstream ss;
 ostringstream forest;
  std::pair <long, long> filler;
  long eventida, runida;
  ofstream myfile;
  ofstream forestfile;
  int vtxCount;
  double InvMass3;
  double InvMass2;
  double delphi;
  double weighta, weightb, weightc;
   Long64_t larks;
   // Declaration of leaf types
   Double_t        lep1pt;
   Double_t        lep2pt;
   Int_t           numberjets;
   Double_t        leadjetpt;
   //   Double_t        2ndjetpt;
   // Double_t        3rdjetpt;
   // Double_t        4thjetpt;
   // Double_t        5thjetpt;
   Double_t        met;
   Double_t        ZpT;
   Double_t        InvMass;
   Int_t           Vertices;
   Int_t           extraElectrons;
   Int_t           extraMuons;
   Double_t        weight;
   Int_t           nPUVerticesTrue;
   //  Double_t        1stEmuonpt;
   //  Double_t        2ndEmuonpt;
   // Double_t        3rdEmuonpt;
   // Double_t        1stEmuoneta;
   //Double_t        2ndEmuoneta;
   //Double_t        3rdEmuoneta;
   //Double_t        1stEelectronpt;
   //Double_t        2ndEelectronpt;
   //Double_t        3rdEelectronpt;
   //Double_t        1stEelectroneta;
   // Double_t        2ndEelectroneta;
   // Double_t        3rdEelectroneta;
   Double_t        rMETparallel;
   Double_t        rMETperp;
   Double_t        redMETtotal;
   Double_t        lepton1phi;
   Double_t        lepton2phi;
   Double_t        lepton1eta;
   Double_t        lepton2eta;
   Int_t           DiLeptonType;
   Double_t        METSignificance;
   Double_t        acol;
   Double_t        lep1dxy;
   Double_t        lep2dxy;
   Double_t        lep1Q;
   Double_t        lep2Q;
   Float_t         vertexProb;
   Double_t        MCinvmass;
   Double_t        MC1pt;
   Double_t        MC1eta;
   Double_t        MC1phi;
   Double_t        MC2pt;
   Double_t        MC2eta;
   Double_t        MC2phi;
   Double_t        MC1Q;
   Double_t        MC2Q;
   Double_t        MCpTdiff;
   bool ZpTcut, METcut, Sumptcut, NumJetscut, acolcut, vtxprobcut, invmasscutabove, invmasscutbelow;
   Bool_t          DoubleMuonTrigger;
   Bool_t          DoubleElectronTrigger;
   Bool_t          MuonElectronTrigger;
   Bool_t          SingleMuonTrigger40;
   Bool_t          SingleMuonTriggeriso;
   Bool_t          SingleMuonTriggeriso2;
   Bool_t          SingleElectronTrigger;
   Int_t           aftertrigger;
   Double_t        liso1;
   Double_t        liso2;
   TLorentzVector  *l14vector;
   TLorentzVector  *l24vector;
   Double_t        lep1ptE;
   Double_t        lep2ptE;
   Double_t        elecReg1;
   Double_t        elecReg2;
   Int_t           eventid;
   Int_t           runid;
   Bool_t          istautau;
   TLorentzVector  *MC14vector;
   TLorentzVector  *MC24vector;
   Double_t        mujetd1;
   Double_t        mujetd2;
   Double_t        isquark;
   Double_t        motherId1;
   Double_t        motherId2;
   Double_t        Npt;
   Double_t        Neta;
   Double_t        Nphi;
   Double_t        Mpt;
   Double_t        Meta;
   Double_t        Mphi;
   Double_t        Nisolation;
   Double_t        Njetd;
   Double_t        Mq;
   Double_t        WqT;
   Int_t           eff1;



  rochcor2012 *muCorrector;
  TGraphErrors *_isorewe[4];
  TGraphErrors *_idrewe[4];
  TH2F *_elrewe;

   //list of histograms

  TH1D* h_jet_pt;

  
  TH1D* h_jet_lead_pt;


  TH1D* h_ZpT_log;
  TH1D* h_ZpT;
  TH1D* h_ZpTee;
  TH1D* h_ZpTmm;
  TH1D* h_ZpTem;
  TH1D* h_ZpTCutee;
  TH1D* h_ZpTCutmm;
  TH1D* h_ZpTCutee2;
  TH1D* h_ZpTCutmm2;
  TH1D* h_ZpTCutem;
  TH1D* h_CutQTee;
  TH1D* h_CutQTmm;
  TH1D* h_CutQTem;

  TH1D* h_IMemnocuts;
  TH1D* h_IMemJetres;
  TH1D* h_IMemVtexres;
  TH1D* h_IMemacolres;
  TH1D* h_IMemSumptres;
  TH1D* h_NJemJetres;
  TH1D* h_Vtexemres;
  TH1D* h_acolemres;
  TH1D* h_sumptemres;

 TH1D* h_IMeenocuts;
  TH1D* h_IMeeJetres;
  TH1D* h_IMeeVtexres;
  TH1D* h_IMeeacolres;
  TH1D* h_IMeeSumptres;
  TH1D* h_NJeeJetres;
  TH1D* h_Vtexeeres;
  TH1D* h_acoleeres;
  TH1D* h_sumpteeres;


 TH1D* h_IMmmnocuts;
  TH1D* h_IMmmJetres;
  TH1D* h_IMmmVtexres;
  TH1D* h_IMmmacolres;
  TH1D* h_IMmmSumptres;
  TH1D* h_NJmmJetres;
  TH1D* h_Vtexmmres;
  TH1D* h_acolmmres;
  TH1D* h_sumptmmres;


  TH1D* h_numJetsem;
  TH1D* h_numJetseminner;
  TH1D* h_numJetseminner2;

  TH1D* h_numJetsemouter;
  TH1D* h_numJetseeinner;
  TH1D* h_numJetseeouter;

  TH1D* h_numJetsmminner;
  TH1D* h_numJetsmmouter;

  TH1D* h_numJetsmm;
  TH1D* h_numJetsee;
  TH1D* h_numJetsmm2;
  TH1D* h_numJetsee2;

  TH1D* h_InvMass;
  TH1D* h_InvMassbase;
  TH1D* h_InvMassbase2;
  TH1D* h_InvMassbase3;
  TH1D* h_InvMassbase4;
  TH1D* h_InvMassnone;
  TH1D* h_MCinvmass;



  TH1D* h_numJets;

  TH1D* h_mET;
  TH1D* h_mETmm;
  TH1D* h_mETee;
  TH1D* h_mETem;
  TH1D* h_mETem2;
  TH1D* h_mETeminner;
  TH1D* h_mETeminner2;
  TH1D* h_mETemouter;
  TH1D* h_mETemouter2;

  TH1D* h_redMETtotalmm;
  TH1D* h_redMETtotalee;
  TH1D* h_redMETtotalem;


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
  TH1D* h_vtxCountmm;
  TH1D* h_vtxCountee;


  TH2D* h_METcomponents;
  TH2D* h_METcomparison;  
  TH1D* h_InvMassee; 
  TH1D* h_InvMassmm; 
  TH1D* h_InvMassem; 
  TH1D* h_InvMasseeinner; 
  TH1D* h_InvMassmminner; 
  TH1D* h_InvMasseminner; 
  TH1D* h_InvMasseminner2; 

  TH1D* h_InvMasseeouter; 
  TH1D* h_InvMassmmouter; 
  TH1D* h_InvMassemouter; 
  TH1D* h_InvMassbelow; 
  TH1D* h_InvMassabove; 
  TH1D* h_InvMassbelowmm; 
  TH1D* h_InvMassabovemm; 
  TH1D* h_InvMassbelowee; 
  TH1D* h_InvMassaboveee; 
  TH1D* h_InvMassnoIso; 
  TH2D* h_InvmassvsMET;
  TH2D* h_InvmassvsMsig;
  TH2D* h_qTvsM;
  TH1D* h_X;

  TH1D* h_METresem;
  TH1D* h_NJresem;
  TH1D* h_SPresem;
  TH1D* h_VPresem;

  TH1D* h_METresee;
  TH1D* h_NJresee; 
  TH1D* h_SPresee;
  TH1D* h_VPresee;

  TH1D* h_METresmm;
  TH1D* h_NJresmm;
  TH1D* h_SPresmm;
  TH1D* h_VPresmm;

  TH1D* h_Zrapid;

  TH1D* h_METSignificance;
  TH1D* h_METSignificancemm;
  TH1D* h_METSignificanceee;
  TH1D* h_METSignificancemm2;
  TH1D* h_METSignificanceee2;
  TH1D* h_METSignificanceem;
  TH1D* h_METSignificance2;
  TH1D* h_Cutmm;
  TH1D* h_Cutee;
  TH1D* h_Cutmm2;

  TH1D* h_Cutmm3;
  TH1D* h_Cutee3;
  TH1D* h_Cutee2;

  TH1D* h_Cutem;
  TH1D* h_Cutem2;
  TH1D* h_acol;
  TH1D* h_acolee;
  TH1D* h_acolmm;
  TH1D* h_acolee2;
  TH1D* h_acolmm2;
  TH1D* h_acolem;
  TH1D* h_Sumpt;
  TH1D* h_Sumptee;
  TH1D* h_Sumptmm;
  TH1D* h_Sumptee2;
  TH1D* h_Sumptmm2;
  TH1D* h_Sumptem;
  TH1D* h_delphi;
  TH1D* h_delphiee;
  TH1D* h_delphimm;
  TH1D* h_delphiem;
  TH1D* h_lep1dxy;
  TH1D* h_lep2dxy;
  TH1D* h_lep1dxyee;
  TH1D* h_lep2dxyee;
  TH1D* h_lep1dxymm;
  TH1D* h_lep2dxymm;
  TH1D* h_lep1dxyem;
  TH1D* h_lep2dxyem;
  TH1D* h_Michael1;
  TH1D* h_Michael2;
  TH1D* h_lepton1pt;
  TH1D* h_lepton2pt;
  TH1D* h_leadptem;
  TH1D* h_subptem;
  TH1D* h_leadptee;
  TH1D* h_subptee;
  TH1D* h_leadptmm;
  TH1D* h_subptmm;
  TH1D* h_leptonbothpt;
  TH1D* h_leptonsmalleta;
  TH1D* h_lepton1smalleta;
  TH1D* h_lepton2smalleta;
  TH1D* h_vertexProb;
  TH1D* h_vertexProbmm;
 TH1D* h_vertexProbmm2;
 TH1D* h_vertexProbee2;
  TH1D* h_vertexProbem;
  //  TH1D* h_liso1;
  TH1D* h_liso1mm;
  TH1D* h_liso2mm;
  TH1D* h_leta1mm;
  TH1D* h_leta2mm;
  TH1D* h_lphi1mm;
  TH1D* h_lphi2mm;
  TH1D* h_leta1ee;
  TH1D* h_leta2ee;
  TH1D* h_lphi1ee;
  TH1D* h_lphi2ee;
  TH1D* h_liso1ee;
  TH1D* h_liso2ee;
  TH1D* h_InvMassee2;

   // List of branches
   TBranch        *b_lep1pt;   //!
   TBranch        *b_lep2pt;   //!
   TBranch        *b_numberjets;   //!
   TBranch        *b_leadjetpt;   //!
   // TBranch        *b_2ndjetpt;   //!
   // TBranch        *b_3rdjetpt;   //!
   // TBranch        *b_4thjetpt;   //!
   // TBranch        *b_5thjetpt;   //!
   TBranch        *b_met;   //!
   TBranch        *b_ZpT;   //!
   TBranch        *b_InvMass;   //!
   TBranch        *b_Vertices;   //!
   TBranch        *b_extraElectrons;   //!
   TBranch        *b_extraMuons;   //!
   TBranch        *b_weight;   //!
   TBranch        *b_nPUVerticesTrue;   //!
   //   TBranch        *b_1stEmuonpt;   //!
   // TBranch        *b_2ndEmuonpt;   //!
   // TBranch        *b_3rdEmuonpt;   //!
   // TBranch        *b_1stEmuoneta;   //!
   // TBranch        *b_2ndEmuoneta;   //!
   // TBranch        *b_3rdEmuoneta;   //!
   // TBranch        *b_1stEelectronpt;   //!
   // TBranch        *b_2ndEelectronpt;   //!
   // TBranch        *b_3rdEelectronpt;   //!
   // TBranch        *b_1stEelectroneta;   //!
   // TBranch        *b_2ndEelectroneta;   //!
   // TBranch        *b_3rdEelectroneta;   //!
   TBranch        *b_rMETparallel;   //!
   TBranch        *b_rMERperp;   //!
   TBranch        *b_redMETtotal;   //!
   TBranch        *b_lepton1phi;   //!
   TBranch        *b_lepton2phi;   //!
   TBranch        *b_lepton1eta;   //!
   TBranch        *b_lepton2eta;   //!
   TBranch        *b_DiLeptonType;   //!
   TBranch        *b_METSignificance;   //!
   TBranch        *b_acol;   //!
   TBranch        *b_lep1dxy;   //!
   TBranch        *b_lep2dxy;   //!
   TBranch        *b_lep1Q;   //!
   TBranch        *b_lep2Q;   //!
   TBranch        *b_vertexProb;   //!
   TBranch        *b_MCinvmass;   //!
   TBranch        *b_MC1pt;   //!
   TBranch        *b_MC1eta;   //!
   TBranch        *b_MC1phi;   //!
   TBranch        *b_MC2pt;   //!
   TBranch        *b_MC2eta;   //!
   TBranch        *b_MC2phi;   //!
   TBranch        *b_MC1Q;   //!
   TBranch        *b_MC2Q;   //!
   TBranch        *b_MCpTdiff;   //!
   TBranch        *b_DoubleMuonTrigger;   //!
   TBranch        *b_DoubleElectronTrigger;   //!
   TBranch        *b_MuonElectronTrigger;   //!
   TBranch        *b_SingleMuonTrigger40;   //!
   TBranch        *b_SingleMuonTriggeriso;   //!
   TBranch        *b_SingleMuonTriggeriso2;   //!
   TBranch        *b_SingleElectronTrigger;   //!
   TBranch        *b_aftertrigger;   //!
   TBranch        *b_liso1;   //!
   TBranch        *b_liso2;   //!
   TBranch        *b_l14vector;   //!
   TBranch        *b_l24vector;   //!
   TBranch        *b_lep1ptE;   //!
   TBranch        *b_lep2ptE;   //!
   TBranch        *b_elecReg1;   //!
   TBranch        *b_elecReg2;   //!
   TBranch        *b_eventid;   //!
   TBranch        *b_runid;   //!
   TBranch        *b_istautau;   //!
   TBranch        *b_MC14vector;   //!
   TBranch        *b_MC24vector;   //!
   TBranch        *b_mujetd1;   //!
   TBranch        *b_mujetd2;   //!
   TBranch        *b_isquark;   //!
   TBranch        *b_motherId1;   //!
   TBranch        *b_motherId2;   //!
   TBranch        *b_Npt;   //!
   TBranch        *b_Neta;   //!
   TBranch        *b_Nphi;   //!
   TBranch        *b_Mpt;   //!
   TBranch        *b_Meta;   //!
   TBranch        *b_Mphi;   //!
   TBranch        *b_Nisolation;   //!
   TBranch        *b_Njetd;   //!
   TBranch        *b_Mq;   //!
   TBranch        *b_WqT;   //!
   TBranch        *b_eff1;   //!
   TFile *f;
   TFile *f_muisorewe;
   TFile *f_muidrewe;
   TFile *f_elrewe;

   AnalyzerZemu(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~AnalyzerZemu() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();
    void fill(TH1D*, int);
    void fill(TH1D*, double);
    //double get_bin(TH1D*, int ibin=1);
    double get_mass(TLorentzVector);
    //    double get_eta(TLorentzVector);
    double get_rapid(TLorentzVector);
    double get_pT(TLorentzVector);
     double DeltaPhiX(double, double);
   // double spreader(double, double);
    float GetMuEff(double, double) const;
    float GetElEff(double, double);

    float TriggerMuEff(double, double);
    float TriggerElEff(double, double);

   ClassDef(AnalyzerZemu,0);
 private:
   // TFile* f;
};

#endif

#ifdef AnalyzerZemu_cxx
void AnalyzerZemu::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   l14vector = 0;
   l24vector = 0;
   MC14vector = 0;
   MC24vector = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("lep1pt", &lep1pt, &b_lep1pt);
   fChain->SetBranchAddress("lep2pt", &lep2pt, &b_lep2pt);
   fChain->SetBranchAddress("numberjets", &numberjets, &b_numberjets);
   fChain->SetBranchAddress("leadjetpt", &leadjetpt, &b_leadjetpt);
   //   fChain->SetBranchAddress("2ndjetpt", &2ndjetpt, &b_2ndjetpt);
   // fChain->SetBranchAddress("3rdjetpt", &3rdjetpt, &b_3rdjetpt);
   // fChain->SetBranchAddress("4thjetpt", &4thjetpt, &b_4thjetpt);
   // fChain->SetBranchAddress("5thjetpt", &5thjetpt, &b_5thjetpt);
   fChain->SetBranchAddress("met", &met, &b_met);
   fChain->SetBranchAddress("ZpT", &ZpT, &b_ZpT);
   fChain->SetBranchAddress("InvMass", &InvMass, &b_InvMass);
   fChain->SetBranchAddress("Vertices", &Vertices, &b_Vertices);
   fChain->SetBranchAddress("extraElectrons", &extraElectrons, &b_extraElectrons);
   fChain->SetBranchAddress("extraMuons", &extraMuons, &b_extraMuons);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("nPUVerticesTrue", &nPUVerticesTrue, &b_nPUVerticesTrue);
   //  fChain->SetBranchAddress("1stEmuonpt", &1stEmuonpt, &b_1stEmuonpt);
   // fChain->SetBranchAddress("2ndEmuonpt", &2ndEmuonpt, &b_2ndEmuonpt);
   // fChain->SetBranchAddress("3rdEmuonpt", &3rdEmuonpt, &b_3rdEmuonpt);
   // fChain->SetBranchAddress("1stEmuoneta", &1stEmuoneta, &b_1stEmuoneta);
   // fChain->SetBranchAddress("2ndEmuoneta", &2ndEmuoneta, &b_2ndEmuoneta);
   // fChain->SetBranchAddress("3rdEmuoneta", &3rdEmuoneta, &b_3rdEmuoneta);
   // fChain->SetBranchAddress("1stEelectronpt", &1stEelectronpt, &b_1stEelectronpt);
   // fChain->SetBranchAddress("2ndEelectronpt", &2ndEelectronpt, &b_2ndEelectronpt);
   // fChain->SetBranchAddress("3rdEelectronpt", &3rdEelectronpt, &b_3rdEelectronpt);
   // fChain->SetBranchAddress("1stEelectroneta", &1stEelectroneta, &b_1stEelectroneta);
   //fChain->SetBranchAddress("2ndEelectroneta", &2ndEelectroneta, &b_2ndEelectroneta);
   //fChain->SetBranchAddress("3rdEelectroneta", &3rdEelectroneta, &b_3rdEelectroneta);
   fChain->SetBranchAddress("rMETparallel", &rMETparallel, &b_rMETparallel);
   fChain->SetBranchAddress("rMETperp", &rMETperp, &b_rMERperp);
   fChain->SetBranchAddress("redMETtotal", &redMETtotal, &b_redMETtotal);
   fChain->SetBranchAddress("lepton1phi", &lepton1phi, &b_lepton1phi);
   fChain->SetBranchAddress("lepton2phi", &lepton2phi, &b_lepton2phi);
   fChain->SetBranchAddress("lepton1eta", &lepton1eta, &b_lepton1eta);
   fChain->SetBranchAddress("lepton2eta", &lepton2eta, &b_lepton2eta);
   fChain->SetBranchAddress("DiLeptonType", &DiLeptonType, &b_DiLeptonType);
   fChain->SetBranchAddress("METSignificance", &METSignificance, &b_METSignificance);
   fChain->SetBranchAddress("acol", &acol, &b_acol);
   fChain->SetBranchAddress("lep1dxy", &lep1dxy, &b_lep1dxy);
   fChain->SetBranchAddress("lep2dxy", &lep2dxy, &b_lep2dxy);
   fChain->SetBranchAddress("lep1Q", &lep1Q, &b_lep1Q);
   fChain->SetBranchAddress("lep2Q", &lep2Q, &b_lep2Q);
   fChain->SetBranchAddress("vertexProb", &vertexProb, &b_vertexProb);
   fChain->SetBranchAddress("MCinvmass", &MCinvmass, &b_MCinvmass);
   fChain->SetBranchAddress("MC1pt", &MC1pt, &b_MC1pt);
   fChain->SetBranchAddress("MC1eta", &MC1eta, &b_MC1eta);
   fChain->SetBranchAddress("MC1phi", &MC1phi, &b_MC1phi);
   fChain->SetBranchAddress("MC2pt", &MC2pt, &b_MC2pt);
   fChain->SetBranchAddress("MC2eta", &MC2eta, &b_MC2eta);
   fChain->SetBranchAddress("MC2phi", &MC2phi, &b_MC2phi);
   fChain->SetBranchAddress("MC1Q", &MC1Q, &b_MC1Q);
   fChain->SetBranchAddress("MC2Q", &MC2Q, &b_MC2Q);
   fChain->SetBranchAddress("MCpTdiff", &MCpTdiff, &b_MCpTdiff);
   fChain->SetBranchAddress("DoubleMuonTrigger", &DoubleMuonTrigger, &b_DoubleMuonTrigger);
   fChain->SetBranchAddress("DoubleElectronTrigger", &DoubleElectronTrigger, &b_DoubleElectronTrigger);
   fChain->SetBranchAddress("MuonElectronTrigger", &MuonElectronTrigger, &b_MuonElectronTrigger);
   fChain->SetBranchAddress("SingleMuonTrigger40", &SingleMuonTrigger40, &b_SingleMuonTrigger40);
   fChain->SetBranchAddress("SingleMuonTriggeriso", &SingleMuonTriggeriso, &b_SingleMuonTriggeriso);
   fChain->SetBranchAddress("SingleMuonTriggeriso2", &SingleMuonTriggeriso2, &b_SingleMuonTriggeriso2);
   fChain->SetBranchAddress("SingleElectronTrigger", &SingleElectronTrigger, &b_SingleElectronTrigger);
   fChain->SetBranchAddress("aftertrigger", &aftertrigger, &b_aftertrigger);
   fChain->SetBranchAddress("liso1", &liso1, &b_liso1);
   fChain->SetBranchAddress("liso2", &liso2, &b_liso2);
   fChain->SetBranchAddress("l14vector", &l14vector, &b_l14vector);
   fChain->SetBranchAddress("l24vector", &l24vector, &b_l24vector);
   fChain->SetBranchAddress("lep1ptE", &lep1ptE, &b_lep1ptE);
   fChain->SetBranchAddress("lep2ptE", &lep2ptE, &b_lep2ptE);
   fChain->SetBranchAddress("elecReg1", &elecReg1, &b_elecReg1);
   fChain->SetBranchAddress("elecReg2", &elecReg2, &b_elecReg2);
   fChain->SetBranchAddress("eventid", &eventid, &b_eventid);
   fChain->SetBranchAddress("runid", &runid, &b_runid);
   fChain->SetBranchAddress("istautau", &istautau, &b_istautau);
   fChain->SetBranchAddress("MC14vector", &MC14vector, &b_MC14vector);
   fChain->SetBranchAddress("MC24vector", &MC24vector, &b_MC24vector);
   fChain->SetBranchAddress("mujetd1", &mujetd1, &b_mujetd1);
   fChain->SetBranchAddress("mujetd2", &mujetd2, &b_mujetd2);
   fChain->SetBranchAddress("isquark", &isquark, &b_isquark);
   fChain->SetBranchAddress("motherId1", &motherId1, &b_motherId1);
   fChain->SetBranchAddress("motherId2", &motherId2, &b_motherId2);
   fChain->SetBranchAddress("Npt", &Npt, &b_Npt);
   fChain->SetBranchAddress("Neta", &Neta, &b_Neta);
   fChain->SetBranchAddress("Nphi", &Nphi, &b_Nphi);
   fChain->SetBranchAddress("Mpt", &Mpt, &b_Mpt);
   fChain->SetBranchAddress("Meta", &Meta, &b_Meta);
   fChain->SetBranchAddress("Mphi", &Mphi, &b_Mphi);
   fChain->SetBranchAddress("Nisolation", &Nisolation, &b_Nisolation);
   fChain->SetBranchAddress("Njetd", &Njetd, &b_Njetd);
   fChain->SetBranchAddress("Mq", &Mq, &b_Mq);
   fChain->SetBranchAddress("WqT", &WqT, &b_WqT);\
   fChain->SetBranchAddress("eff1", &eff1, &b_eff1);
}

Bool_t AnalyzerZemu::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef AnalyzerZemu_cxx
