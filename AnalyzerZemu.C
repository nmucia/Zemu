#define AnalyzerZemu_cxx
// The class definition in AnalyzerZemu.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("AnalyzerZemu.C")
// Root > T->Process("AnalyzerZemu.C","some options")
// Root > T->Process("AnalyzerZemu.C+")
//

#include "AnalyzerZemu.h"
#include <TH2.h>
#include <TStyle.h>


void AnalyzerZemu::Begin(TTree * /*tree*/)
{
  filetype=0;
  IS_DATA=false;
  larks=0;
  produceForest=0;
  muCorrector = new rochcor2012();
   numa=numb=0;
   

 if(filetype==2)
    {
      long a,b;
      cout<<"ddddddddd"<<endl;
      std::ifstream infile("/uscms_data/d3/nmucia2/dblcountfiles/SMA.txt");
      if(infile.is_open())
	{
	  while (infile>>a>>b)  
	    {
	      // 	      cout<<a<<" "<<b<<endl;
	      foo=make_pair(a,b);
	      oldevents.push_back(foo);
	    }
	}
      sort(oldevents.begin(), oldevents.end());
    }
  f_muisorewe = new TFile("MuonEfficiencies_ISO_Run_2012ReReco_53X.root", "OPEN"); 
  _isorewe[0] = (TGraphErrors*)f_muisorewe->Get("DATA_over_MC_combRelIsoPF04dBeta<012_Tight_pt_abseta<0.9");
  _isorewe[1] = (TGraphErrors*)f_muisorewe->Get("DATA_over_MC_combRelIsoPF04dBeta<012_Tight_pt_abseta0.9-1.2");
  _isorewe[2] = (TGraphErrors*)f_muisorewe->Get("DATA_over_MC_combRelIsoPF04dBeta<012_Tight_pt_abseta1.2-2.1");
  _isorewe[3] = (TGraphErrors*)f_muisorewe->Get("DATA_over_MC_combRelIsoPF04dBeta<012_Tight_pt_abseta2.1-2.4");
  f_muisorewe->Close();

  f_muidrewe = new TFile("MuonEfficiencies_Run2012ReReco_53X.root", "OPEN"); 
  _idrewe[0] = (TGraphErrors*)f_muidrewe->Get("DATA_over_MC_Tight_pt_abseta<0.9");
  _idrewe[1] = (TGraphErrors*)f_muidrewe->Get("DATA_over_MC_Tight_pt_abseta0.9-1.2");
  _idrewe[2] = (TGraphErrors*)f_muidrewe->Get("DATA_over_MC_Tight_pt_abseta1.2-2.1");
      _idrewe[3] = (TGraphErrors*)f_muidrewe->Get("DATA_over_MC_Tight_pt_abseta2.1-2.4");  
  f_muidrewe->Close();       
  f_elrewe = new TFile("electrons_scale_factors.root", "OPEN"); 
  _elrewe = (TH2F*)f_elrewe->Get("electronsDATAMCratio_FO_ID_ISO");

 //Madgraph forr 2012 DoubleMu C
  

  f = new TFile("rootfiles/outfileTW1.root","RECREATE");           






  //  float xlow[8] = {.0035,.0035,.030,.025,.025,.025,.025,.025,.025,.025,.025,.025,.025,.025,.025,.025,.025,.025,.025,.025,.025,.025,.01,.01,.01,.01,.01,.05,.1,.5,1};

  // Float_t xlow[35] = {0.00,0.004,0.008,0.012,0.016,0.020,0.024,0.029,0.034,0.039,0.045,0.051,0.057,0.064,0.072,0.081,0.091,0.102,0.114,0.128,0.145,0.165,0.189,0.219,0.258,0.312,.391,.524,.695,.918,1.153,1.496,1.947,2.522,3.277};
  // Float_t invmassmm[35] ={52,56,60,64,68,72,76,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,104,108,112,116,120,124,128};
  h_jet_pt = new TH1D("h_jet_pt","jet pt", 100, 0, 1000.);
  // h_jet_pt_50 = new TH1D("h_jet_pt_50","jet pt, qT>50", 100, 0, 1000.);
  // h_jet_pt_100 = new TH1D("h_jet_pt_100","jet pt, qT>100", 100, 0, 1000.);
  // h_jet_pt_200 = new TH1D("h_jet_pt_200","jet pt, qT>200", 100, 0, 1000.);
  // h_jet_pt_300 = new TH1D("h_jet_pt_300","jet pt, qT>300", 100, 0, 1000.);
  // h_jet_pt_500 = new TH1D("h_jet_pt_500","jet pt, qT>500", 100, 0, 1000.);
  //  h_jet_pt_700 = new TH1D("h_jet_pt_700","jet pt, qT>700", 100, 0, 1000.);
  h_jet_lead_pt = new TH1D("h_jet_lead_pt","lead jet pt", 100, 0, 1000.);
  // h_jet_lead_pt_50 = new TH1D("h_jet_lead_pt_50","lead jet pt, qT>50", 100, 0, 1000.);
  // h_jet_lead_pt_100 = new TH1D("h_jet_lead_pt_100","lead jet pt, qT>100", 100, 0, 1000.);
  // h_jet_lead_pt_200 = new TH1D("h_jet_lead_pt_200","lead jet pt, qT>200", 100, 0, 1000.);
  // h_jet_lead_pt_300 = new TH1D("h_jet_lead_pt_300","lead jet pt, qT>300", 100, 0, 1000.);
  // h_jet_lead_pt_500 = new TH1D("h_jet_lead_pt_500","lead jet pt, qT>500", 100, 0, 1000.);
  //  h_jet_lead_pt_700 = new TH1D("h_jet_lead_pt_700","lead jet pt, qT>700", 100, 0, 1000.);

   h_mET = new TH1D("h_mET","mET", 50, 0, 300.);
   h_mETmm = new TH1D("h_mETmm","mET", 50, 0, 300.);
   h_mETee = new TH1D("h_mETee","mET", 50, 0, 300.);
   h_mETem = new TH1D("h_mETem","mET", 50, 0, 300.);
   h_mETemouter = new TH1D("h_mETemouter","mET", 50, 0, 300.);
   h_mETemouter2 = new TH1D("h_mETemouter2","mET", 50, 0, 300.);

   h_mETem2 = new TH1D("h_mETem2","mET", 50, 0, 300.);

   h_mETeminner = new TH1D("h_mETeminner","mET", 50, 0, 300.);
   h_mETeminner2 = new TH1D("h_mETeminner2","mET", 50, 0, 300.);


   h_redMETtotalmm = new TH1D("h_redMETtotalmm", "reduced met", 100, 0, 100);
   h_redMETtotalee = new TH1D("h_redMETtotalee", "reduced met", 100, 0, 100);
   h_redMETtotalem = new TH1D("h_redMETtotalem", "reduced met", 100, 0, 100);

  h_ZpT_log = new TH1D("h_ZpT_log","log plot of qT", 140, 0, 700);
  h_InvMass = new TH1D("h_InvMass","Invariant Mass of the Z", 100, 50, 150);
  h_MCinvmass = new TH1D("h_MCinvmass","Invariant Mass of the Z", 100, 50, 150);

  h_ZpT = new TH1D("h_ZpT","ZpT", 50, 0., 100);
  h_numJets = new TH1D("h_numJets","number of jets", 11, -.5, 10.5);


  h_extra_electrons = new TH1D("h_extra_electrons", "# of extra electrons", 6, -.5, 5.5);
  h_extra_muons = new TH1D("h_extra_muons", "# of extra muons", 6, -.5, 5.5);
  
  //  h_electron_pt = new TH1D("h_electron_pt", "electron pt", 200,0,200);
  // h_muon_pt = new TH1D("h_muon_pt", "muon pt", 200,0,200);



  
  h_nPUVertices = new TH1D("PUVertices", "PUVertices", 70, 0, 70);
  h_nPUVerticesTrue = new TH1D("PUVerticesTrue", "PUVerticesTrue", 70, 0, 70);
 h_vtxCount = new TH1D("vertex count", "vertex count", 70, 0, 70); 
 h_vtxCountmm= new TH1D("vertex countmm", "vertex count mm", 70, 0, 70); 
 h_vtxCountee= new TH1D("vertex countee", "vertex count mm", 70, 0, 70); 

 h_InvMassnone = new TH1D("h_InvMassnone", "InvMass loose cuts", 100, 50, 250);

 // h_dReEm1 = new TH1D("h_dReEm1", "dR between elec and first muon", 100, 0, 6);
 //h_dReEm2 = new TH1D("h_dReEm2", "dR between elec and second muon", 100, 0, 6);
 
 // h_rMET = new TH1D("h_rMET", "reduced MET", 50,0,300);
 //h_rMETx = new TH1D("h_rMETx", "reduced MET parallel", 100,-100,300);
 //h_rMETy = new TH1D("h_rMETy", "reduced MET perpendiculat", 50,-100,100);
 // h_changeMET = new TH1D("h_changeMET", "change between MET and rMET", 100,-50,50);
 h_METcomponents = new TH2D("h_METcomponents","MET perp vs MET parallel", 200,-100,100,200,-100,100);
 h_METcomparison = new TH2D("h_METcomparison","rMET vs MET ", 200,0,100,200,0,100);
 h_InvMassnoIso = new TH1D("h_InvMassnoIso", "InvMass with no iso", 100,50,250);
 h_InvmassvsMET = new TH2D("h_InvmassvsMET", "MET vs Inv Mass",30,60,120,40,0,80);
 h_InvmassvsMsig = new TH2D("h_InvmassvsMsig", "MET vs Inv Mass",60,60,120,10,0,10);
 h_qTvsM = new TH2D("h_qTvsM", "qT vs M", 200, 0, 1000, 200, 0, 1000);
 h_X = new TH1D("h_X","sqrt(M*M + qT*qT)", 500, 0, 1000);
 // h_XoverM = new TH1D("XoverM"," X over M", 100, 0, 10);


 h_Zrapid = new TH1D("h_Zrapid","Z rapidity",60, -3, 3);
 //h_BarreltestNick = new TH1D("h_BarreltestNick","Invariant Mass Barrel mm",34, invmassmm);
 // h_BarreltestNickBB = new TH1D("h_BarreltestNickBB","Invariant Mass Barrel mm",34, invmassmm);
 // h_BarreltestNickBE = new TH1D("h_BarreltestNickBE","Invariant Mass Barrel mm",34, invmassmm);
 // h_BarreltestNickEE = new TH1D("h_BarreltestNickEE","Invariant Mass Barrel mm",34, invmassmm);

  h_InvMassee = new TH1D("h_InvMassee","Invariant Mass of the Z ee", 70, 60, 130);
 h_InvMasseeinner = new TH1D("h_InvMasseeinner","Invariant Mass of the Z ee", 70, 60, 130);
 h_InvMasseeouter = new TH1D("h_InvMasseeouter","Invariant Mass of the Z ee", 70, 60, 130);
  h_InvMassee2 = new TH1D("h_InvMassee2","Invariant Mass of the Z ee", 70, 60, 130);
  h_InvMassbase = new TH1D("h_InvMassbase","Invariant Mass of the Z mm", 70, 60, 130);
 h_InvMassbase2 = new TH1D("h_InvMassbase2","Invariant Mass of the Z mm", 70, 60, 130);
 h_InvMassbase3 = new TH1D("h_InvMassbase3","Invariant Mass of the Z mm", 70, 60, 130);
 h_InvMassbase4 = new TH1D("h_InvMassbase4","Invariant Mass of the Z mm", 70, 60, 130);
  h_InvMassmm = new TH1D("h_InvMassmm","Invariant Mass of the Z mm", 70, 60, 130);
  h_InvMassmminner = new TH1D("h_InvMassmminner","Invariant Mass of the Z mm", 70, 60, 130);
  h_InvMassmmouter = new TH1D("h_InvMassmmouter","Invariant Mass of the Z mm", 70, 60, 130);

  h_InvMassem = new TH1D("h_InvMassem","Invariant Mass of the Z em", 70, 60, 130);
  h_InvMasseminner = new TH1D("h_InvMasseminner","Invariant Mass of the Z em", 70, 60, 130); 
  h_InvMasseminner2 = new TH1D("h_InvMasseminner2","Invariant Mass of the Z em", 70, 60, 130); 

  h_InvMassemouter = new TH1D("h_InvMassemouter","Invariant Mass of the Z em", 70, 60, 130);
  h_InvMassbelow = new TH1D("h_InvMassbelow","Invariant Mass of the Z em", 70, 60, 130);
  h_InvMassabove = new TH1D("h_InvMassabove","Invariant Mass of the Z em", 70, 60, 130);
  h_InvMassbelowee = new TH1D("h_InvMassbelowee","Invariant Mass of the Z em", 70, 60, 130);
  h_InvMassaboveee= new TH1D("h_InvMassaboveee","Invariant Mass of the Z em", 70, 60, 130);  
  h_InvMassbelowmm = new TH1D("h_InvMassbelowmm","Invariant Mass of the Z em", 70, 60, 130);
  h_InvMassabovemm = new TH1D("h_InvMassabovemm","Invariant Mass of the Z em", 70, 60, 130);

  h_IMmmnocuts = new TH1D("h_IMmmnocuts","Invariant Mass of the Z em", 70, 60, 130);
  h_IMmmJetres = new TH1D("h_IMmmJetres","Invariant Mass of the Z em", 70, 60, 130);
  h_IMmmVtexres = new TH1D("h_IMmmVtexres","Invariant Mass of the Z em", 70, 60, 130);
  h_IMmmacolres = new TH1D("h_IMmmacolres","Invariant Mass of the Z em", 70, 60, 130);
  h_IMmmSumptres = new TH1D("h_IMmmSumptres","Invariant Mass of the Z em", 70, 60, 130);
  h_NJmmJetres = new TH1D("h_NJmmJerres","number of jets", 11, -.5, 10.5);
  h_Vtexmmres = new TH1D("h_Vtexmmres","vertex probability", 80, 0, 4);
  h_acolmmres = new TH1D("h_acolmmres","pi - acolineariry", 32, 0, 3.2);
   h_sumptmmres = new TH1D("h_Sumptmmres","scalar sum of both lepton pt", 100, 0, 500.);

   h_IMeenocuts = new TH1D("h_IMeenocuts","Invariant Mass of the Z em", 70, 60, 130);
  h_IMeeJetres = new TH1D("h_IMeeJetres","Invariant Mass of the Z em", 70, 60, 130);
  h_IMeeVtexres = new TH1D("h_IMeeVtexres","Invariant Mass of the Z em", 70, 60, 130);
  h_IMeeacolres = new TH1D("h_IMeeacolres","Invariant Mass of the Z em", 70, 60, 130);
  h_IMeeSumptres = new TH1D("h_IMeeSumptres","Invariant Mass of the Z em", 70, 60, 130);
  h_NJeeJetres = new TH1D("h_NJeeJerres","number of jets", 11, -.5, 10.5);
  h_Vtexeeres = new TH1D("h_Vtexeeres","vertex probability", 80, 0, 4);
  h_acoleeres = new TH1D("h_acoleeres","pi - acolineariry", 32, 0, 3.2);
   h_sumpteeres = new TH1D("h_Sumpteeres","scalar sum of both lepton pt", 100, 0, 500.);

 h_IMemnocuts = new TH1D("h_IMemnocuts","Invariant Mass of the Z em", 70, 60, 130);
 h_IMemJetres = new TH1D("h_IMemJetres","Invariant Mass of the Z em", 70, 60, 130);
 h_IMemVtexres = new TH1D("h_IMemVtexres","Invariant Mass of the Z em", 70, 60, 130);
  h_IMemacolres = new TH1D("h_IMemacolres","Invariant Mass of the Z em", 70, 60, 130);
  h_IMemSumptres = new TH1D("h_IMemSumptres","Invariant Mass of the Z em", 70, 60, 130);
  h_NJemJetres = new TH1D("h_NJemJerres","number of jets", 11, -.5, 10.5);
  h_Vtexemres = new TH1D("h_Vtexemres","vertex probability", 80, 0, 4);
  h_acolemres = new TH1D("h_acolemres","pi - acolineariry", 32, 0, 3.2);
   h_sumptemres = new TH1D("h_Sumptemres","scalar sum of both lepton pt", 100, 0, 500.);


   h_ZpTee = new TH1D("h_ZpTee","ZpT", 50, 0., 100);
   h_ZpTmm = new TH1D("h_ZpTmm","ZpT", 50, 0., 100);
   h_ZpTem = new TH1D("h_ZpTem","ZpT", 50, 0., 100);
   h_ZpTCutee = new TH1D("h_ZpTCutee","ZpT", 50, 0., 100);
   h_ZpTCutmm = new TH1D("h_ZpTCutmm","ZpT", 50, 0., 100);
   h_ZpTCutee2 = new TH1D("h_ZpTCutee2","ZpT", 50, 0., 100);
   h_ZpTCutmm2 = new TH1D("h_ZpTCutmm2","ZpT", 50, 0., 100);
   h_ZpTCutem = new TH1D("h_ZpTCutem","ZpT", 50, 0., 100);
  h_numJetsee = new TH1D("h_numJetsee","number of jets", 11, -.5, 10.5);
  h_numJetsmm = new TH1D("h_numJetsmm","number of jets", 11, -.5, 10.5);
  h_numJetsee2 = new TH1D("h_numJetsee2","number of jets", 11, -.5, 10.5);
  h_numJetsmm2 = new TH1D("h_numJetsmm2","number of jets", 11, -.5, 10.5);
  h_numJetsem = new TH1D("h_numJetsem","number of jets", 11, -.5, 10.5);
  h_numJetseminner = new TH1D("h_numJetseminner","number of jets", 11, -.5, 10.5);
  h_numJetseminner2 = new TH1D("h_numJetseminner2","number of jets", 11, -.5, 10.5);

  h_numJetsemouter = new TH1D("h_numJetsemouter","number of jets", 11, -.5, 10.5);
 h_numJetseeinner = new TH1D("h_numJetseeinner","number of jets", 11, -.5, 10.5);
  h_numJetseeouter = new TH1D("h_numJetseeouter","number of jets", 11, -.5, 10.5);
 h_numJetsmminner = new TH1D("h_numJetsmminner","number of jets", 11, -.5, 10.5);
  h_numJetsmmouter = new TH1D("h_numJetsmmouter","number of jets", 11, -.5, 10.5);
   h_METSignificance = new TH1D("h_METSignificance","MET significance", 51, -.5, 50.5);
   h_METSignificanceee = new TH1D("h_METSignificanceee","MET significance", 51, -.5, 50.5);
   h_METSignificancemm = new TH1D("h_METSignificancemm","MET significance", 51, -.5, 50.5);
   h_METSignificanceee2 = new TH1D("h_METSignificanceee2","MET significance", 51, -.5, 50.5);
   h_METSignificancemm2 = new TH1D("h_METSignificancemm2","MET significance", 51, -.5, 50.5);
   h_METSignificanceem = new TH1D("h_METSignificanceem","MET significance", 51, -.5, 50.5);
   h_METSignificance2 = new TH1D("h_METSignificance2","MET significance2", 51, -.5, 50.5);
 h_Cutmm = new TH1D("h_Cutmm","Invariant Mass with cuts MM", 60, 60, 120);
 h_Cutee = new TH1D("h_Cutee","Invariant Mass with cuts ee", 60, 60, 120);
 h_Cutem = new TH1D("h_Cutem","Invariant Mass with cuts eM", 60, 60, 120);

 h_Cutmm2 = new TH1D("h_Cutmm2","Invariant Mass with cuts MM", 60, 60, 120);
 h_Cutee2 = new TH1D("h_Cutee2","Invariant Mass with cuts ee", 60, 60, 120);
 h_Cutem2 = new TH1D("h_Cutem2","Invariant Mass with cuts eM", 60, 60, 120);
 h_CutQTmm = new TH1D("h_CutQTmm","Invariant Mass with cuts MM", 60, 60, 120);
 h_CutQTee = new TH1D("h_CutQTee","Invariant Mass with cuts ee", 60, 60, 120);
 h_CutQTem = new TH1D("h_CutQTem","Invariant Mass with cuts eM", 60, 60, 120);
 h_Cutmm3 = new TH1D("h_Cutmm3","Invariant Mass with cuts eM", 70, 60, 120);
 h_Cutee3 = new TH1D("h_Cutee3","Invariant Mass with cuts eM", 70, 60, 120);

h_acol = new TH1D("h_acol","pi - acolineariry", 32, 0, 3.2);
h_acolee = new TH1D("h_acolee","pi - acolineariry", 32, 0, 3.2);
h_acolmm = new TH1D("h_acolmm","pi - acolineariry", 32, 0, 3.2);
h_acolee2 = new TH1D("h_acolee2","pi - acolineariry", 32, 0, 3.2);
h_acolmm2 = new TH1D("h_acolmm2","pi - acolineariry", 32, 0, 3.2);
h_acolem = new TH1D("h_acolem","pi - acolineariry", 32, 0, 3.2);
h_delphi = new TH1D("h_delphi","phi angle between leps", 50, 0, 3.2);
h_delphiee = new TH1D("h_delphiee","phi angle between leps", 50, 0, 3.2);
h_delphimm = new TH1D("h_delphimm","phi angle between leps", 50, 0, 3.2);
h_delphiem = new TH1D("h_delphiem","phi angle between leps", 50, 0, 3.2);
   h_Sumpt = new TH1D("h_Sumpt","scalar sum of both lepton pt", 100, 0, 1000.);
   h_Sumptee = new TH1D("h_Sumptee","scalar sum of both lepton pt", 100, 0, 1000.);
   h_Sumptmm = new TH1D("h_Sumptmm","scalar sum of both lepton pt", 100, 0, 1000.);
   h_Sumptem = new TH1D("h_Sumptem","scalar sum of both lepton pt", 60, 20, 140.);
   h_lep1dxy = new TH1D("h_lep1dxy","beam spot information", 100, 0, 0.5);
   h_lep2dxy = new TH1D("h_lep2dxy","beam spot information", 100, 0, 0.5);
   h_lep1dxyee = new TH1D("h_lep1dxyee","beam spot information", 100, 0, 0.5);
   h_lep2dxyee = new TH1D("h_lep2dxyee","beam spot information", 100, 0, 0.5);
   h_lep1dxymm = new TH1D("h_lep1dxymm","beam spot information", 100, 0, 0.05);
   h_lep2dxymm = new TH1D("h_lep2dxymm","beam spot information", 100, 0, 0.05);
   h_lep1dxyem = new TH1D("h_lep1dxyem","beam spot information", 100, 0, 0.05);
   h_lep2dxyem = new TH1D("h_lep2dxyem","beam spot information", 100, 0, 0.5); 
   h_Michael1 = new TH1D("h_Michael1","Invariant Mass of the mumu", 70, 60, 130); 
   h_Michael2 = new TH1D("h_Michael2","Invariant Mass of the emu", 70, 60, 130);
   h_lepton1pt = new TH1D("h_lepton1pt","lepton pt", 60, 0, 60.);
   h_lepton2pt = new TH1D("h_lepton2pt","lepton pt", 60, 0, 60.);
   h_leadptem = new TH1D("h_leadptem","lepton pt", 100, 0, 100.);
   h_subptem = new TH1D("h_subptem","lepton pt", 100, 0, 100.);
   h_leadptmm = new TH1D("h_leadptmm","lepton pt", 100, 0, 100.);
   h_subptmm = new TH1D("h_subptmm","lepton pt", 100, 0, 100.); 
   h_leadptee = new TH1D("h_leadptee","lepton pt", 100, 0, 100.);
   h_subptee = new TH1D("h_subptee","lepton pt", 100, 0, 100.);
   h_leptonbothpt = new TH1D("h_leptonbothpt","lepton pt", 60, 0, 60.);
   h_leptonsmalleta = new TH1D("h_leptonsmalleta","lepton pt", 60, 0, 60.);
   h_lepton1smalleta = new TH1D("h_lepton1smalleta","lepton pt", 60, 0, 60.);
   h_lepton2smalleta = new TH1D("h_lepton2smalleta","lepton pt", 60, 0, 60.);
   h_vertexProb = new TH1D("h_vertexProb","vertex probability", 80, 0, 4);
   h_vertexProbmm = new TH1D("h_vertexProbmm","vertex probability", 80, 0, 4);
   h_vertexProbmm2 = new TH1D("h_vertexProbmm2","vertex probability", 80, 0, 4);
   h_vertexProbee2 = new TH1D("h_vertexProbee2","vertex probability", 80, 0, 4);

   h_vertexProbem = new TH1D("h_vertexProbem","vertex probability", 80, 0, 4);
   h_liso1mm = new TH1D("h_liso1mm","liso1 mm", 100, 0.0, 1.0);
   h_liso2mm = new TH1D("h_liso2mm","liso2 mm", 100, 0.0, 1.0);
   h_leta1mm = new TH1D("h_leta1mm","leta1", 100, -2.5, 2.5);
   h_leta2mm = new TH1D("h_leta2mm","leta2", 100, -2.5, 2.5);
   h_lphi1mm = new TH1D("h_lphi1mm","lphi1", 100, 0.0, 3.2);
   h_lphi2mm = new TH1D("h_lphi2mm","lphi2", 100, 0.0, 3.2);
   h_leta1ee = new TH1D("h_leta1ee","leta1", 100, -2.5, 2.5);
   h_leta2ee = new TH1D("h_leta2ee","leta2", 100, -2.5, 2.5);
   h_lphi1ee = new TH1D("h_lphi1ee","lphi1", 100, 0.0, 3.2);
   h_lphi2ee = new TH1D("h_lphi2ee","lphi2", 100, 0.0, 3.2);
   h_liso1ee = new TH1D("h_liso1ee","liso1 ee", 100, 0.0, 0.15);
   h_liso2ee = new TH1D("h_liso2ee","liso2 ee", 100, 0.0, 0.15);

  h_METresee = new TH1D("h_METresee","MET significance", 51, -.5, 50.5);
  h_VPresee = new TH1D("h_VPresee","vertex probability", 80, 0, 4);
 h_NJresee = new TH1D("h_NJresee","number of jets", 11, -.5, 10.5);
 h_SPresee = new TH1D("h_SPresee","scalar sum of both lepton pt", 100, 0, 500.);

 h_METresmm = new TH1D("h_METresmm","MET significance", 51, -.5, 50.5);
  h_VPresmm = new TH1D("h_VPresmm","vertex probability", 80, 0, 4);
 h_NJresmm = new TH1D("h_NJresmm","number of jets", 11, -.5, 10.5);
 h_SPresmm = new TH1D("h_SPresmm","scalar sum of both lepton pt", 100, 0, 500.);

 h_METresem = new TH1D("h_METresem","MET significance", 51, -.5, 50.5);
  h_VPresem = new TH1D("h_VPresem","vertex probability", 80, 0, 4);
 h_NJresem = new TH1D("h_NJresem","number of jets", 11, -.5, 10.5);
 h_SPresem = new TH1D("h_SPresem","scalar sum of both lepton pt", 100, 0, 500.);

  cout<<"Ok so far ..."<<endl;
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

void AnalyzerZemu::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t AnalyzerZemu::Process(Long64_t entry)
{
  GetEntry(entry);
	 notdoublecounted=true;
  ZpTcut=METcut=invmasscutabove=invmasscutbelow=Sumptcut=NumJetscut=acolcut=vtxprobcut= false;
  weighta=weightb=weightc=1.0;
  lep1pt=l14vector->Pt();
  lep2pt=l24vector->Pt();
  //SingleElectronTrigger=0;
  larks++;
      if((larks%10000)==0)
	 cout<<larks<<endl;
  if(-log10(vertexProb)>4)
    vertexProb=.00009;
    if(!(IS_DATA))
      {    
		       		weight=weights[nPUVerticesTrue];
     
    /////// jet-reweight
	//    weight=weight*(jet_weight[numberjets]/1.03);
    // cout<<numberjets<<" , "<<jet_weight[numberjets]<<endl;
      }
 h_nPUVerticesTrue->Fill(nPUVerticesTrue, 1);  	
 
   if(acol<5)
     acolcut=true;
   if((lep1pt+lep2pt)<95)
     Sumptcut=true;
   if(met<30)
     METcut=true;
   if(numberjets==0)
     NumJetscut=true;
   if((-log10(vertexProb))<3.0)
     vtxprobcut=true;
   if(InvMass>88)
     invmasscutbelow=true;
   if(InvMass<96)
     InvMass2=true;
   if(ZpT<30)
     ZpTcut=true;
  //  h_vtxCount->Fill(vtxCount, weight);
   // fill(h_InvMassmm, InvMass);
   //   delphi=3.14159265359-DeltaPhiX(lep1phi,lep2phi);

   //This is for the acceptance study, not needed normally
    
      InvMass2=get_mass(*l14vector+*l24vector);
      //              cout<<MC14vector->M()<<" , "<<MC14vector->Pt()<<endl;
      /*      InvMass3=get_mass(*MC14vector+*MC24vector);
       cout<<eff1<<endl;
       if(eff1==1)
	 {
	   fill(h_InvMassbase, InvMass3);
	   if(MC14vector->Pt()>30 && MC24vector->Pt()>20 && fabs(MC14vector->Eta())<2.1 && fabs(MC24vector->Eta()) <2.1 &&InvMass3<130 && InvMass3>60)
	     {
	       fill(h_InvMassbase2, InvMass3);
	       if(DiLeptonType==5)
		   {
		     fill(h_InvMassbase3, InvMass3);
		     if((SingleMuonTriggeriso==1||SingleElectronTrigger==1))
		     {
		       fill(h_InvMassbase4, InvMass3);

		     }
		   }
	     }
	 }
      */
      
      if(DiLeptonType==5)
	{
	  if(liso1<0.12 && liso2<0.15 && l14vector->Pt()>20 && l24vector->Pt()>20&&(l14vector->Pt()>30 ||l24vector->Pt()>30 ) )
	     {
	       if(l14vector->Pt()>l24vector->Pt())
		 {
		   leadeta=l14vector->Eta(); 
		 }
	       else
		 {
		   leadeta=l24vector->Eta(); 
		 }
	       if(fabs(leadeta)<2.1 && fabs(l14vector->Eta())<2.4 && fabs(l14vector->Eta())<2.4 )
		 {
		   if(fabs(lep2eta)<1.44 || fabs(lep2eta)>1.57)
		     {
		       InvMass2=get_mass(*l14vector+*l24vector);
		       numa++;
		       if(SingleMuonTriggeriso==1||SingleElectronTrigger==1)
			 {
			   numb++;
			 }
		     }
		 }
	     }
	  
	}
      
      //    if( (SingleMuonTriggeriso==1||SingleElectronTrigger==1))
      if(lep1Q==-lep2Q && (SingleMuonTriggeriso==1||SingleElectronTrigger==1))
	 {
  if(filetype==1)
	{
	  //  cout<<runid<<" "<<eventid<<endl;
	 	 ss << runid << " " << eventid << "\n";
		 //    string mySting = ss.str();
			//string hellos;
			//hellos <<runid;
			//cout<<hellos<<endl;
		 //			       		 cout<<mySting<<endl;
	 }
       //checks if event was used before
   

          if(filetype==2)
	 {
	       runida=runid;
	     eventida=eventid;

	     	   foo=make_pair(runida,eventida);
		   //cout<<foo.first<<" "<<foo.second<<endl;
		   //for use in unsorted array
		   //		   if(find(oldevents.begin(), oldevents.end(), foo) != oldevents.end())
		   //for use if oldevents is sorted
		   if(binary_search(oldevents.begin(), oldevents.end(), foo))
	      {notdoublecounted=false;
		cout<<larks<<" duplicate event found"<<endl;
	      }


	 }

	        if(notdoublecounted)
		  // if(true)
	 {
	 	   /////////////////////////////////////muons/////////////////////////////////////////////////

   if(DiLeptonType==1) 
     {
       if(IS_DATA)
	 {
	    muCorrector->momcor_data(*l14vector,lep1Q,0,qter); 
	    muCorrector->momcor_data(*l24vector,lep2Q,0,qter);  
	   lep1pt=l14vector->Pt();
	   lep2pt=l24vector->Pt();
	   lep1eta=l14vector->Eta();
	   lep2eta=l24vector->Eta();
	 }
       else
	 {
	   //	  cout<<endl<<endl<<weight<<endl;
	   
	   //	    cout<<weighta<<" , "<<l14vector->Pt()<<" , "<<l14vector->Eta()<<" , "<<weight<<endl;
	   
	   
	    muCorrector->momcor_mc(*l14vector,lep1Q,0,qter); 
	   muCorrector->momcor_mc(*l24vector,lep2Q,0,qter); 
	   //  cout<<l14vector->Pt();
	   lep1pt=l14vector->Pt();
	   lep2pt=l24vector->Pt();
	   lep1eta=l14vector->Eta();
	   lep2eta=l24vector->Eta();	    
	   
	   	    weighta=GetMuEff(lep1pt, lep1eta);
	        	weightb=GetMuEff(lep2pt, lep2eta);
			//  weightc=triggerMuEff(lep1pt, lep1eta);
	   //	      cout<<weight<<" , ";
			//apply single muon or single electron or both trigger weights to the MC files
			      if(SingleMuonTriggeriso==1 && !(SingleElectronTrigger==1))
	{
	  weightc=TriggerMuEff(lep1pt, lep1eta);
	}
      else if(!(SingleMuonTriggeriso==1) && (SingleElectronTrigger==1))
	{
	  	  weightc=TriggerElEff(lep1pt, lep1eta);
	    }
      else if(SingleMuonTriggeriso==1 && (SingleElectronTrigger==1))
	{
	  	  weightc=TriggerMuEff(lep1pt, lep1eta)*TriggerElEff(lep1pt, lep1eta);
	}


	         weight=weight*weighta*weightb*weightc;
		 // if(weight<0.000000001)
		   //cout<<lep1eta<<" , "<<lep2eta<<" , "<<endl;
	   // weight=1.0;
	   //cout<<weight<<endl;

	 }
            InvMass2=get_mass(*l14vector+*l24vector);
       if( fabs(lep1eta) <2.1 && fabs(lep2eta)<2.4 && lep1pt>30 && lep2pt>20 && InvMass2>60&& InvMass2<120 && extraElectrons<1 && extraMuons<1 && vertexProb>-.001)
	 {
	   if(extraElectrons>1)
	     {
	       //cout<<extraElectrons<<endl;
	     }
	       fill(h_liso1mm, liso1);
	       fill(h_liso2mm, liso2);
	      

	   if(liso1<0.12 &&liso2<0.12)
	     {
	       fill(h_leta1mm, lep1eta);
	       fill(h_leta2mm, lep2eta);
	       fill(h_lphi1mm, l14vector->Phi());
	       fill(h_lphi2mm, l24vector->Phi());
	       
	       fill(h_vtxCountmm, Vertices);
	       fill(h_leadptmm, lep1pt);
	       fill(h_subptmm, lep2pt);
	      //   fill(h_lepton1pt, lep1pt);
	      //    fill(h_lepton2pt, lep2pt);
	      //    fill(h_leptonbothpt, lep1pt);
	      //   fill(h_leptonbothpt, lep2pt);
	      // if(fabs(lep1eta)<1)
	      //	{	
	      //  fill(h_leptonsmalleta, lep1pt);
	      //	  fill(h_lepton1smalleta, lep1pt);
	      //	}
	      //   if(fabs(lep2eta)<1)
	      //	{
	      //  fill(h_leptonsmalleta, lep2pt);
	      //  fill(h_lepton2smalleta, lep2pt);
	      //	} 
	      fill(h_InvMassmm, InvMass2);

	     	      fill(h_numJetsmm, numberjets);
		      //double hello= sqrt((l14vector->Px()+l24vector->Px())**2+(l14vector->Py()+l24vector->Py())**2);
	      
		        fill(h_ZpTmm,ZpT);
			fill(h_mETmm, met);
			fill(h_redMETtotalmm, redMETtotal);
	       fill(h_METSignificancemm, met);
	       fill(h_delphimm, delphi); 
	        fill(h_Sumptmm, (lep1pt+lep2pt));
	      //  fill(h_lep1dxymm, lep1dxy);
	           fill(h_vertexProbmm, -log10(vertexProb));
	      //fill(h_lep2dxymm, lep2dxy);
	       fill(h_acolmm, acol);
	       if(Sumptcut&&vtxprobcut&&NumJetscut&&(InvMass2<74 || InvMass2>106))
		 {
		   fill(h_METresmm, METSignificance);
		 }
	       if((met>50 ||( InvMass2<74 || InvMass2>106)))
		 {
		   fill(h_IMmmnocuts, InvMass2);
		   if(Sumptcut&&vtxprobcut)
		     {
		       fill(h_NJresmm, numberjets);
		     }
		   if(Sumptcut&&NumJetscut)
		     {
		       fill(h_SPresmm, lep1pt+lep2pt);
		     }
		   if(NumJetscut&&vtxprobcut)
		     {
		       fill(h_VPresmm, -log10(vertexProb));
		     }
		   if(Sumptcut)
		     { 
		       if(acolcut&&ZpTcut)
			 {
			   if(vtxprobcut)
			     {
			       fill(h_numJetsmmouter, numberjets);
			       if(NumJetscut)
				 {
			       fill(h_InvMassmmouter, InvMass2);
				 }
			       else
				 {
				   fill(h_IMmmJetres, InvMass2);
				   fill(h_NJmmJetres, numberjets);
				 }
			     }
			   else
			     {
			       fill(h_IMmmVtexres, InvMass2);
			       fill(h_Vtexmmres, -log10(vertexProb));
			     }
			   
			 }
		       else
			 {
			   fill(h_IMmmacolres, InvMass2);
			   fill(h_acolmmres, acol);
			 }
		     }
		   else
		     {
		       fill(h_IMmmSumptres, InvMass2);
		       fill(h_sumptmmres, lep1pt+lep2pt);
		     }

		 }


	       if(((met>30 && met<50) || (InvMass2<88 ||InvMass2>96)) && InvMass2>74 && InvMass2<106)
			 {  
			   if(Sumptcut && acolcut && vtxprobcut && ZpTcut)
			     {
			   fill(h_numJetsmminner, numberjets);
			   fill(h_InvMassmminner, InvMass2);
			     }
			 }
	     			   
 if(METcut && Sumptcut && acolcut && vtxprobcut && InvMass2>88 && InvMass2<96 && ZpTcut)
	{
	  fill(h_numJetsmm2, numberjets);
	}
	       
     if(METcut && Sumptcut&& NumJetscut && vtxprobcut && InvMass2>88 && InvMass2<96 && ZpTcut)
       {
	 fill(h_acolmm2, acol);
       }
     
     if(METcut && Sumptcut && NumJetscut && acolcut && InvMass2>88 && InvMass2<96 && ZpTcut)
       {
      	fill(h_vertexProbmm2, -log10(vertexProb));
       }
     
      // fill(h_delphiem, delphi); 
     //   fill(h_lep1dxyem, lep1dxy);
     // fill(h_lep2dxyem, lep2dxy);
     
     //   if(numberjets==0 && InvMass>60 && InvMass<130 && METSignificance<4 && acol<2.5 && (lep1pt+lep2pt)>50
     if(METcut && Sumptcut && NumJetscut && acolcut && vtxprobcut && InvMass2>88 && ZpTcut)
       {
	 fill(h_InvMassabovemm, InvMass2);
       }
     
     if(METcut && Sumptcut && NumJetscut && acolcut && vtxprobcut && InvMass2<96 && ZpTcut)
       {
	 fill(h_InvMassbelowmm, InvMass2);
	     }
     
     if( Sumptcut && acolcut && vtxprobcut && InvMass2>88 && InvMass2<96 && ZpTcut)
       {
	 fill(h_METSignificancemm2, met);
	 
       }
     if(METcut && Sumptcut && NumJetscut && acolcut && vtxprobcut && ZpTcut)
       {
	 fill(h_Cutmm3, InvMass2);
       }
     
     if(METcut && Sumptcut && NumJetscut && acolcut && vtxprobcut && InvMass2>88 && InvMass2<96)
       
       {
	 fill(h_ZpTCutmm2, ZpT);
       }
		   if(METcut && Sumptcut && NumJetscut && acolcut && vtxprobcut && InvMass2>88 && InvMass2<96 && ZpTcut)
			{
			  fill(h_Cutmm, InvMass2);
			  fill(h_ZpTCutmm, ZpT);
			  if(ZpT<ZpTcut)
			    {
			      fill(h_CutQTmm, InvMass);
			      }
			}
	     }
	 }
     }
	 

   ///////////////////////////////////////////////////////electrons//////////////////////////////////////////////////////
       //  if(DiLeptonType==3  && aftertrigger==3)
             
       if(DiLeptonType==3)
	 {
	   
	   lep1eta=l14vector->Eta();
	   lep2eta=l24vector->Eta();
	   lep1phi=l14vector->Phi();
	   lep2phi=l24vector->Phi();
	   lep1pt=l14vector->Pt();
	   lep2pt=l24vector->Pt();

	   if(!(IS_DATA))
	     {
	       weighta=GetElEff(lep1pt, lep1eta);
	       weightb=GetElEff(lep2pt, lep2eta);
			//apply single muon or single electron or both trigger weights to the MC files

      if(SingleMuonTriggeriso==1 && !(SingleElectronTrigger==1))
	{
	  weightc=TriggerMuEff(lep1pt, lep1eta);
	}
      else if(!(SingleMuonTriggeriso==1) && (SingleElectronTrigger==1))
	{
	  weightc=TriggerElEff(lep1pt, lep1eta);
	    }
      else if(SingleMuonTriggeriso==1 && (SingleElectronTrigger==1))
	{
	  weightc=TriggerMuEff(lep1pt, lep1eta)*TriggerElEff(lep1pt, lep1eta);
	}

	       weight=weight*weighta*weightb*weightc;
	     }	   

	   // cout<<lep1eta<<endl;
	   if(liso1<0.15 &&liso2<0.15 && fabs(lep1eta)<2.1 && fabs(lep2eta)<2.4 )
	     {
	       	   
	       
	       if(fabs(lep1eta)<1.44 || fabs(lep1eta)>1.57)
		 {
		   if(fabs(lep2eta)<1.44 || fabs(lep2eta)>1.57)
		     // && leadpt>20 &&subpt>10)
		     
		     {
		       //cout<<lep1eta<<" , "<<lep2eta<<endl;
		       //	   cout<<elecReg1<<endl;
		       // cout<<l14vector->Pt()<<endl;
		       // cout<<l24vector->Pt()<<endl;
		       
		       //       if(!(IS_DATA))
		       //	 {
			   //	   l14vector->SetPtEtaPhiE(l14vector->Pt()*elecReg1/l14vector->E(),lep1eta,lep1phi,elecReg1);
			   //   l24vector->SetPtEtaPhiE(l24vector->Pt()*elecReg2/l24vector->E(),lep2eta,lep2phi,elecReg2);
		       //	 }
		   	       
		       if(l14vector->Pt()>20 && l24vector->Pt()>20&&(l14vector->Pt()>30 ||l24vector->Pt()>30 ) && InvMass2<120 && InvMass2>60 && extraElectrons <1 && extraMuons<1)
			 {
			   if(Sumptcut&&vtxprobcut&&NumJetscut&&(InvMass2<74 || InvMass2>106))
		 {
		   fill(h_METresee, met);
		 }
			   if((met>50) ||(InvMass2<74 || InvMass2>106))
		 {
		   fill(h_IMeenocuts, InvMass2);
		   if(Sumptcut&&vtxprobcut)
		     {
		       fill(h_NJresee, numberjets);
		     }
		   if(Sumptcut&&NumJetscut)
		     {
		       fill(h_SPresee, lep1pt+lep2pt);
		     }
		   if(NumJetscut&&vtxprobcut)
		     {
		       fill(h_VPresee, -log10(vertexProb));
		     }
		   if(Sumptcut)
		     { 
		       if(acolcut&&ZpTcut)
			 {
			   if(vtxprobcut)
			     {
			       fill(h_numJetseeouter, numberjets);
			       if(NumJetscut)
				 {
			       fill(h_InvMasseeouter, InvMass2);
				 }
			       else
				 {
				   fill(h_IMeeJetres, InvMass2);
				   fill(h_NJeeJetres, numberjets);
				 }
			     }
			   else
			     {
			       fill(h_IMeeVtexres, InvMass2);
			       fill(h_Vtexeeres, -log10(vertexProb));
			     }
			   
			 }
		       else
			 {
			   fill(h_IMeeacolres, InvMass2);
			   fill(h_acoleeres, acol);
			 }
		     }
		   else
		     {
		       fill(h_IMeeSumptres, InvMass2);
		       fill(h_sumpteeres, lep1pt+lep2pt);
		     }

		 }
			   if(Sumptcut && acolcut && vtxprobcut&&ZpTcut)
			     {
			       
			     
			       if(((met>30 && met<50) || (InvMass2<88 ||InvMass2>96)) && InvMass2>74 && InvMass2<106)

				 {     
				   fill(h_InvMasseeinner, InvMass2);
				   fill(h_numJetseeinner, numberjets);
				 }
			     }
			 
			 
			   // cout<<" , "<<l14vector->Pt()<<endl;
			   InvMass2=get_mass(*l14vector+*l24vector);
			   fill(h_InvMassee2, InvMass2);
			   	       fill(h_vtxCountee, Vertices);

			      fill(h_liso1ee, liso1);
			      fill(h_liso2ee, liso2);
			    fill(h_leta1ee, lep1eta);
			    fill(h_leta2ee, lep2eta);
			    fill(h_lphi1ee, lep1phi);
			  
			   
			   
			    fill(h_lphi2ee, lep2phi);
			   fill(h_leadptee, l14vector->Pt());
			   fill(h_subptee, l24vector->Pt());
			   // fill(h_vertexProb, -log10(vertexProb));
			   fill(h_InvMassee, InvMass2);
			   fill(h_ZpTee, ZpT);
			   fill(h_numJetsee, numberjets);
			   fill(h_METSignificanceee, met);
			   fill(h_mETee, met);
			fill(h_redMETtotalee, redMETtotal);
			   //	   fill(h_delphiee, delphi); 
			   	   fill(h_Sumptee, (lep1pt+lep2pt));
			   //	   fill(h_lep1dxyee, lep1dxy);
			   //	   fill(h_lep2dxyee, lep2dxy);
			   fill(h_acolee, acol);
			   
 if(METcut && Sumptcut && acolcut && vtxprobcut && InvMass2>88 && InvMass2<96 && ZpTcut)
	{
	  fill(h_numJetsee2, numberjets);
	}
	       
     if(METcut && Sumptcut&& NumJetscut && vtxprobcut && InvMass2>88 && InvMass2<96 && ZpTcut)
       {
	 fill(h_acolee2, acol);
       }
     
     if(METcut && Sumptcut && NumJetscut && acolcut && InvMass2>88 && InvMass2<96 && ZpTcut)
       {
      	fill(h_vertexProbee2, -log10(vertexProb));
       }
     
      // fill(h_delphiem, delphi); 
     //   fill(h_lep1dxyem, lep1dxy);
     // fill(h_lep2dxyem, lep2dxy);
     
     //   if(numberjets==0 && InvMass>60 && InvMass<130 && METSignificance<4 && acol<2.5 && (lep1pt+lep2pt)>50
     if(METcut && Sumptcut && NumJetscut && acolcut && vtxprobcut && InvMass2>88 && ZpTcut)
       {
	 fill(h_InvMassaboveee, InvMass2);
       }
     
     if(METcut && Sumptcut && NumJetscut && acolcut && vtxprobcut && InvMass2<96 && ZpTcut)
       {
	 fill(h_InvMassbelowee, InvMass2);
	     }
     
     if( Sumptcut && acolcut && vtxprobcut && InvMass2>88 && InvMass2<96 && ZpTcut)
       {
	 fill(h_METSignificanceee2, met);
	 
       }
     if(METcut && Sumptcut && NumJetscut && acolcut && vtxprobcut && ZpTcut)
       {
	 fill(h_Cutee3, InvMass2);
       }
     
     if(METcut && Sumptcut && NumJetscut && acolcut && vtxprobcut && InvMass2>88 && InvMass2<96)
       
       {
	 fill(h_ZpTCutee2, ZpT);
       }
			   if(METcut && Sumptcut && NumJetscut && acolcut && vtxprobcut && InvMass2>88 && InvMass2<96 && ZpTcut)
			     {
			       fill(h_Cutee, InvMass2);
			       fill(h_ZpTCutee, ZpT);
			       if(ZpT<ZpTcut)
				 {
				   fill(h_CutQTee, InvMass);
				 }
			       
			     }
			 }
		     }  
		 }
	     }
	 }
	 
  //  if(DiLeptonType==5 && aftertrigger==5)
       if(DiLeptonType==5)
	 {
	   weighta=weightb=weightc=1.0;
      if(IS_DATA)
	{
	  muCorrector->momcor_data(*l14vector,lep1Q,0,qter); 
	   lep1pt=l14vector->Pt();

	}
      else
	{
	  
	  muCorrector->momcor_mc(*l14vector,lep1Q,0,qter);
	  lep1pt=l14vector->Pt();
	   lep1eta=l14vector->Eta();

	  weighta=GetMuEff(lep1pt, lep1eta);
	  weightb=GetElEff(lep2pt, lep2eta);
	  //  weightc=triggerEff(lep1pt, lep1eta, lep2pt, lep2eta);
	  //	      cout<<weight<<" , ";
			//apply single muon or single electron or both trigger weights to the MC files

      if(SingleMuonTriggeriso==1 && !(SingleElectronTrigger==1))
	{
	  weightc=TriggerMuEff(lep1pt, lep1eta);
	}
      else if(!(SingleMuonTriggeriso==1) && (SingleElectronTrigger==1))
	{
	  weightc=TriggerElEff(lep1pt, lep1eta);
	    }
      else if(SingleMuonTriggeriso==1 && (SingleElectronTrigger==1))
	{
	  weightc=TriggerMuEff(lep1pt, lep1eta)*TriggerElEff(lep1pt, lep1eta);
	}

      weight=weight*weighta*weightb*weightc; 
	}
	   lep1eta=l14vector->Eta();
	   lep2eta=l24vector->Eta();
	   lep1phi=l14vector->Phi();
	   lep2phi=l24vector->Phi();
	   lep2pt=l24vector->Pt();
	   if(liso1<0.12 && liso2<0.5 && l14vector->Pt()>20 && l24vector->Pt()>20&&(l14vector->Pt()>30 ||l24vector->Pt()>30 ) )
	     {
	       if(l14vector->Pt()>l24vector->Pt())
		 {
		   leadeta=l14vector->Eta(); 
		 }
	       else
		 {
		   leadeta=l24vector->Eta(); 
		 }
	       if(fabs(leadeta)<2.1 && fabs(l14vector->Eta())<2.4 && fabs(l14vector->Eta())<2.4 )
		 {
	   if(fabs(lep2eta)<1.44 || fabs(lep2eta)>1.57)
	     {
             InvMass2=get_mass(*l14vector+*l24vector);

	     if(produceForest==1)
	       
	       {
		 if(liso1<0.12 && liso2<0.5)
		   {
		     forest<<"1"<<" "<<l24vector->Pt() << " " <<l14vector->Pt()<< " "<<ZpT<< " "<<met<<" "<<DeltaPhiX(lep1phi,lep2phi)<<" "<<"0"<<" "<<numberjets<<" "<<"0"<<" "<<"0"<<" "<<"0"<<" "<<"0"<<" "<<"0"<<" "<<"0"<<" "<<"0"<<" "<<"1"<<" "<<InvMass2<<" "<<(extraElectrons+extraMuons)<<" "<<"0"<<" "<<"0"<<" "<<met<<" "<<"0"<<" "<<lep1pt+lep2pt<<" "<<METSignificance<<"\n";
		     
		   
		   }
	       }

	     if( InvMass2>60&& InvMass2<120 && extraElectrons<1 && extraMuons<1)
	       {
		 

		   fill(h_mETem, met);

		     fill(h_InvMassem, InvMass2);

		 if(Sumptcut && acolcut && vtxprobcut)
		   {
		     //fill(h_InvMassem, InvMass2);
		   }
		 fill(h_ZpTem, ZpT);

	   fill(h_leadptem, l14vector->Pt());
	   fill(h_subptem, l24vector->Pt());


		       
     //  cout<<"hello"<<endl;
	   if(Sumptcut&&vtxprobcut&&NumJetscut&&ZpTcut&&(InvMass2<74 || InvMass2>106))
		 {
		   fill(h_METresem, met);
		    fill(h_mETem2, met);
		   fill(h_redMETtotalem, redMETtotal);
		 }
  
	       
		   //   	   if(METcut && NumJetscut && acolcut && vtxprobcut && InvMass2>88 && InvMass2<96)
		   //    {
		   //     fill(h_Sumptem, (lep1pt+lep2pt));
		   //  }
	       
		   if((met>50) ||(InvMass2<74 || InvMass2>106))
		 {
		   fill(h_IMemnocuts, InvMass2);
		   fill(h_mETemouter, met);
		   if(Sumptcut&&vtxprobcut)
		     {
		       fill(h_NJresem, numberjets);
		     }
		   if(Sumptcut&&NumJetscut)
		     {
		       fill(h_SPresem, lep1pt+lep2pt);
		     }
		   if(NumJetscut&&vtxprobcut)
		     {
		       fill(h_VPresem, -log10(vertexProb));
		     }
		   if(Sumptcut)
		     { 
		       if(acolcut&&ZpTcut)
			 {
			   if(vtxprobcut)
			     {
			       fill(h_numJetsemouter, numberjets);
			       if(NumJetscut)
				 {
				   fill(h_mETemouter2, met);
				   
				   fill(h_InvMassemouter, InvMass2);
				 }
			       else
				 {
				   
				      fill(h_IMemJetres, InvMass2);
				   fill(h_NJemJetres, numberjets);
				 }
			     }
			   else
			     {
			       fill(h_IMemVtexres, InvMass2);
			       fill(h_Vtexemres, -log10(vertexProb));
			     }
			   
			 }
		       else
			 {
			   fill(h_IMemacolres, InvMass2);
			   fill(h_acolemres, acol);
			 }
		     }
		   else
		     {
		       fill(h_IMemSumptres, InvMass2);
		       fill(h_sumptemres, lep1pt+lep2pt);
		     }

		 }		       
		       if(((met>30 && met<50) || (InvMass2<88 ||InvMass2>96)) && InvMass2>74 && InvMass2<106)
			 
			 { 
			   fill(h_mETeminner, met);
			   fill(h_InvMasseminner, InvMass2);
			   fill(h_numJetseminner, numberjets);
			   if(Sumptcut && acolcut && vtxprobcut && ZpTcut)
			     {
			       fill(h_numJetseminner2, numberjets);  
			       if(NumJetscut)
				 {
				   fill(h_mETeminner2, met);
				   fill(h_InvMasseminner2, InvMass2);
				 }
			     }
			 }			     
      if(METcut && Sumptcut && acolcut && vtxprobcut && InvMass2>88 && InvMass2<96 && ZpTcut)
	{

	  fill(h_numJetsem, numberjets);
	}
	       
     if(METcut && Sumptcut&& NumJetscut && vtxprobcut && InvMass2>88 && InvMass2<96 && ZpTcut)
       {

	 fill(h_acolem, acol);
       }
     
     if(METcut && Sumptcut && NumJetscut && acolcut && InvMass2>88 && InvMass2<96 && ZpTcut)
       {

      	fill(h_vertexProbem, -log10(vertexProb));
       }
     
      // fill(h_delphiem, delphi); 
     //   fill(h_lep1dxyem, lep1dxy);
     // fill(h_lep2dxyem, lep2dxy);
     
     //   if(numberjets==0 && InvMass>60 && InvMass<130 && METSignificance<4 && acol<2.5 && (lep1pt+lep2pt)>50
     if(METcut && Sumptcut && NumJetscut && acolcut && vtxprobcut && InvMass2>88 && ZpTcut)
       {

	 fill(h_InvMassabove, InvMass2);
       }
     
     if(METcut && Sumptcut && NumJetscut && acolcut && vtxprobcut && InvMass2<96 && ZpTcut)
       {

	 fill(h_InvMassbelow, InvMass2);
	     }
     
     if( Sumptcut && acolcut && vtxprobcut && InvMass2>88 && InvMass2<96 && ZpTcut && NumJetscut)
       {

	 fill(h_METSignificanceem, met);
	 
       }
     if(METcut && Sumptcut && NumJetscut && acolcut && vtxprobcut && ZpTcut)
       {

	 fill(h_Cutem2, InvMass2);
       }
     
     if(METcut && Sumptcut && NumJetscut && acolcut && vtxprobcut && InvMass2>88 && InvMass2<96)       
       {

	 fill(h_ZpTCutem, ZpT);
       }		  

     if(Sumptcut && NumJetscut && acolcut && vtxprobcut && ZpTcut)
       {

	    h_InvmassvsMET->Fill(InvMass, met, weight);
	 // h_InvmassvsMsig->Fill(InvMass, met, weight);
	 if(METcut  && InvMass2>88 && InvMass2<96)
	   {
	 fill(h_Cutem, InvMass2);
	 
	 fill(h_METSignificance2, met);
	   }

	 if(ZpT<ZpTcut)
	   {

	     fill(h_CutQTem, InvMass);
	   }
	 
       }
	       }
	     }     
	   
		 }
	     }
	 }
       //end emu
	 }
	 }
	 
  // cout<<eventid<<endl;
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either AnalyzerZemu::GetEntry() or TBranch::GetEntry()
   // to read either all or the required parts of the data. When processing
   // keyed objects with PROOF, the object is already loaded and is available
   // via the fObject pointer.
   //
   // This function should contain the "body" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.
	 
   return kTRUE;
}

void AnalyzerZemu::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void AnalyzerZemu::Terminate()
{

   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
  double ratio=numa/numb;
  cout<<numa<<" , "<<numb<<" "<<ratio<<endl;
  if(produceForest==1)
    {
      forestfile.open("/uscms_data/d3/nmucia2/dblcountfiles/ForestEMU2.txt");
      string Mystring =forest.str();
      forestfile<< Mystring;
      
      forestfile.close();
    }

  if(filetype==1)
	 {
	   myfile.open("/uscms_data/d3/nmucia2/dblcountfiles/SMA.txt");
      string mySting = ss.str();

       myfile<< mySting;

      myfile.close();
	 }
  f_elrewe->Close();

   f->Write();
   f->Close();


}
//===================================================================//

void AnalyzerZemu::fill(TH1D* h, int value)
{
  fill(h, (double)(value));
}

//===================================================================//

void AnalyzerZemu::fill(TH1D* h, double value)
{
  //  cout<<weight<<endl;

  h->Fill(value, weight);
}

//===================================================================//

double AnalyzerZemu::get_pT(TLorentzVector v)
{
  return sqrt(v.Px()*v.Px() + v.Py()*v.Py());
}
////
double AnalyzerZemu::get_mass(TLorentzVector pll)
{
  return pll.M();
}
double AnalyzerZemu::get_rapid(TLorentzVector pll)
{
  return pll.Rapidity();
}
float AnalyzerZemu::GetMuEff(double pt, double eta) const
{ // TGraphErrors *_isorewe[4];
  //  TFile* f_muisorewe = new TFile("MuonEfficiencies_ISO_Run_2012ReReco_53X.root", "OPEN"); 
    //    _isorewe[0] = (TGraphErrors*)f_muisorewe->Get("DATA_over_MC_combRelIsoPF04dBeta<012_Tight_pt_abseta<0.9");
     //    _isorewe[1] = (TGraphErrors*)f_muisorewe->Get("DATA_over_MC_combRelIsoPF04dBeta<012_Tight_pt_abseta0.9-1.2");
      //      _isorewe[2] = (TGraphErrors*)f_muisorewe->Get("DATA_over_MC_combRelIsoPF04dBeta<012_Tight_pt_abseta1.2-2.1");
  //   _isorewe[3] = (TGraphErrors*)f_muisorewe->Get("DATA_over_MC_combRelIsoPF04dBeta<012_Tight_pt_abseta2.1-2.4");
    int etaBin = 0;
    float binningEta[] = {0., 0.9, 1.2, 2.1, 2.4};
    float weighta = 1.;
 
    for (int i = 0; i < 4; ++i) {
        if (fabs(eta) > binningEta[i] && fabs(eta) <= binningEta[i+1]) {
            etaBin = i;
            break;
        }
    }

    if (pt < 500.)
      {
	//	cout<<weighta<<endl;
	weighta =_isorewe[etaBin]->Eval(pt);
	//      cout<<weighta<<endl;
	weighta=weighta*_idrewe[etaBin]->Eval(pt);
	//    cout << etaBin << "\t" << lep.Pt() << "\t" << weighta << endl;
      }
    else
      weighta = 1;
    
    return weighta;
}

float AnalyzerZemu::GetElEff(double pt, double eta) 
{

  float weighta=1.0;
  int etaBin=0;
  int ptBin=0;
  float binningEta[] = {0.0,0.8,1.4442,2.0,2.5};
  float binningPt[] = {10.,15.,20.,30.,40. ,50.,9999999999999.,};
  if(pt<15)
    { return 1.0;
    }

  for(int i=0; i<4; i++)
    {
      if(fabs(eta)> binningEta[i] && fabs(eta)<= binningEta[i+1])
	{
	etaBin= i;
	break;
	}
    }
  for(int i=0; i<6; i++)
    {
      if(pt>binningPt[i] && pt<binningPt[i+1])
	{
	  ptBin = i;
	  break;
	}

    }
  weighta = _elrewe->GetBinContent(etaBin+1, ptBin+1);
 


  //  cout<<endl<<endl; 
 //_elrewe->Print("elerewe.eps");

  

  return weighta;

}
//==================================================================//
 float AnalyzerZemu::TriggerMuEff(double pt, double eta)
 {

  if(pt<25)
    return 1.;
 if(fabs(eta)<0.9)
   return 0.9837;
 if(fabs(eta)>0.9 &&fabs(eta)<1.2)
   return 0.9656;
 if(fabs(eta)>1.2)
   return 0.9962;
 else
   return 1.;

   

 }
//===================================================================//
 float AnalyzerZemu::TriggerElEff(double pt, double eta)
 {

  if(pt<30)
    return 1.;
 if(fabs(eta)<0.8)
   {
     if(pt>30 && pt<40)
       return 0.987;
  if(pt>40 && pt<50)
       return 0.997;
  if(pt>50)
       return 0.998;
   }

 if(fabs(eta)>0.8 &&fabs(eta)<1.2)
   {
  if(pt>30 && pt<40)
       return 0.964;
  if(pt>40 && pt<50)
       return 0.980;
  if(pt>50)
       return 0.988;
   }
 if(fabs(eta)>1.2)   
   {
     if(pt>30 && pt<40)
       return 1.004;
     if(pt>40 && pt<50)
       return 1.033;
     if(pt>50)
       return 0.976;
   }
 else
   return 1.;

   
 
 }


double AnalyzerZemu::DeltaPhiX(double phi1, double phi2)
{
  double dphi = (phi1 - phi2);

  while(dphi > 3.14159265359)      dphi -= 2.*3.14159265359;
  while(dphi <= -3.14159265359)    dphi += 2.*3.14159265359;
  
  return fabs(dphi);
}

//====================================================================//
