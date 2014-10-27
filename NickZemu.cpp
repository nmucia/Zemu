#include "TMath.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1.h"
#include "THStack.h"
#include "TChain.h"
#include "TMath.h"
#include <TBranch.h>
#include <TROOT.h>
#include <TRint.h>
#include <vector>
#include <math.h>
#include <TPad.h>
#include <string>
//This is trial code to produce n-1 plots Last Updated 1-30-2013
using namespace std;
const int N_MC =9;

const int N_HIST = 475;

// changed to three files to compare powheg to madgraph
//const char* FILE_NAME[N_MC] = {"HistsDYme.root", "HistsWWme.root", "Histsttbarme.root", "HistsZZme.root","HistsWZme.root", "HistsWWme.root"};
//const char* FILE_NAME[N_MC] = {"rootfiles/met912/MG.root","rootfiles/met912/signal.root","rootfiles/met912/ttbar.root","rootfiles/met912/ZZ.root","rootfiles/met912/WZ.root","rootfiles/met912/WW.root","rootfiles/met912/WW.root","rootfiles/met912/all.root", "rootfiles/913PHtt.root"};

//const char* FILE_NAME[N_MC] = {"w4v/electrons/DY.root","w4v/electrons/DY.root","w4v/electrons/ttbar.root","w4v/electrons/ZZ.root","w4v/electrons/WZ.root","w4v-/electrons/WW.root","w4v/electrons/WJets.root","w4v/electrons/combo.root"};

const char* FILE_NAME[N_MC] = {"rootfiles/924test/MG.root","rootfiles/923m30q30/signal.root","rootfiles/924test/ttbar.root","rootfiles/924test/ZZ.root","rootfiles/924test/WZ.root","rootfiles/924test/WW.root","rootfiles/924test/WW.root","rootfiles/923m30q30/total.root", "rootfiles/924JRWm30q30/WW.root"};

//const char* FILE_NAME[N_MC] = {"rootfiles/923m30q20/MG.root","rootfiles/923m30q20/signal.root","rootfiles/923m30q20/ttbar.root","rootfiles/923m30q20/ZZ.root","rootfiles/923m30q20/WZ.root","rootfiles/923m30q20/WW.root","rootfiles/923m30q20/WW.root","rootfiles/923/total.root", "rootfiles/923m30q20/WW.root"};
//const char* FILE_NAME[N_MC] = {"w4v/electrons/DY.root","w4v/electrons/DY.root","w4v/electrons/ttbar.root","w4v/electrons/ZZ.root","w4v/electrons/WZ.root","w4v-/electrons/WW.root","w4v/electrons/WJets.root","w4v/electrons/combo.root"};
//const char* FILE_NAME[N_MC] = {"w4v/523/PHmm9.root","w4v/523/PHtt.root","w4v/523/ttbar.root","w4v/523/ZZ.root","w4v/523/WZ.root","w4v/523/WW.root","w4v/523/PHtt.root","w4v/523/D.root"};
//const char* FILE_NAME[N_MC] = {"Reg/DYee.root","Reg/DYee.root","Reg/ttbar.root","Reg/ZZ.root","Reg/WZ.root","Reg/WW.root","Reg/PHtt.root","Reg/C.root"};

//const char* FILE_NAME[N_MC] = {"Histsupdate/HistsDYupdate.root", "Histsupdate/HistsDYupdate.root", "Histsupdate/Histttbarupdate5.root","Histsupdate/HistZZupdate5.root","Histsupdate/WZupdate5.root","Histsupdate/WWupdate5.root","Histsupdate/HistWupdate5.root", "Histsupdate/HistsCEupdate5.root"};
//const char* FILE_NAME[1] = {"tprimeNick.root", "tprimeNick.root"};

//
//
//

class NickZemu{

public:
  
  strats();
   ~strats();
  
  void open_files();
  void load_histograms();
  void make_plots();
  void print_summary();
  void n_one_plots();
  void n_one_plots2();
  void yields();
  void plot(int);
  void dataComparisons(int, TString, TString);
  void n95(int, bool, TString);
  void jet_reweighting();
  void PhiStar();
  TH1D* h[N_MC][N_HIST];
  TH1D* g[N_HIST];
  THStack* hs[N_HIST];
  TH2D* MET;
  TH2D* MET1;
  TH1D* has;
  TH1D* has2;
  TH1D* has3;
  //TH2D* h[N_MC];
  TCanvas* c[N_HIST];
  TPad* cc_1;
  TPad* cc_2;
  TFile* f11;

  //TCanvas* c1;
  //TCanvas* c2;
  //TCanvas* c3;
    
private:
  TFile* file[N_MC];
  void Ratio(int, int, int, TString, TString);
  void scale_histogram(int, int);
  double get_bin(TH1D*, int ibin=1);
};
//
//methods
//

NickZemu::NickZemu()
{
  cout<<endl<<"Starting.."<<endl<<endl;
  
}


void NickZemu::open_files()
{
  for(int k=0; k< N_MC; k++)
    {
      cout << "Opening"<<FILE_NAME[k]<< endl;
      
      file[k] = new TFile(FILE_NAME[k], "READ");
      f11 = new TFile("Zemu.root",  "RECREATE");
    }
  
}

//================================================================//

void NickZemu::load_histograms()
{
  
  cout<<"Loading " <<N_HIST << " histograms" << endl;
  
  TString name_h;

  for(int k=0; k < N_MC; k++)
    {
      for(int i=0; i < N_HIST; i++)
	{
	  name_h ="demo/h_";
	  
	  if(i<10) name_h += "00";
	  else if(i<100) name_h +="0";
	  
	  
	  name_h +=i;
	  
	  //	  h[k][i] = (TH1D*)(file[k]->Get(name_h));
	  
	  //	cout<<"hello "<<k<<" "<<i<<endl;
	  
	  //  h[k][i]->Sumw2();
	  
	  
	  
	//	double l=0;
	//	l= h[k][i]->Integral();
	//	cout<<l<<endl;
      }
         
      h[k][100] = (TH1D*)(file[k]->Get("h_METSignificanceee"));
      h[k][101] = (TH1D*)(file[k]->Get("h_METresmm"));
      h[k][102] = (TH1D*)(file[k]->Get("h_NJresmm"));
      h[k][103] = (TH1D*)(file[k]->Get("h_SPresmm"));
      h[k][104] = (TH1D*)(file[k]->Get("h_VPresmm"));
      h[k][105] = (TH1D*)(file[k]->Get("h_METresee"));
      h[k][106] = (TH1D*)(file[k]->Get("h_NJresee"));
      h[k][108] = (TH1D*)(file[k]->Get("h_SPresee"));
      h[k][109] = (TH1D*)(file[k]->Get("h_VPresee"));
      h[k][110] = (TH1D*)(file[k]->Get("h_METresem"));
      h[k][111] = (TH1D*)(file[k]->Get("h_NJresem"));
      h[k][112] = (TH1D*)(file[k]->Get("h_SPresem"));
      h[k][113] = (TH1D*)(file[k]->Get("h_VPresem"));

      //      h[k][101] = (TH1D*)(file[k]->Get("h_jet_pt_50"));
      //  h[k][102] = (TH1D*)(file[k]->Get("h_jet_pt_100"));
      // h[k][103] = (TH1D*)(file[k]->Get("h_jet_pt_200"));
      // h[k][104] = (TH1D*)(file[k]->Get("h_jet_pt_300"));
      //h[k][105] = (TH1D*)(file[k]->Get("h_jet_pt_500"));
      //h[k][106] = (TH1D*)(file[k]->Get("h_jet_pt_700"));
      h[k][107] = (TH1D*)(file[k]->Get("h_jet_lead_pt"));
      // h[k][108] = (TH1D*)(file[k]->Get("h_jet_lead_pt_50"));
      // h[k][109] = (TH1D*)(file[k]->Get("h_jet_lead_pt_100"));
      // h[k][110] = (TH1D*)(file[k]->Get("h_jet_lead_pt_200"));
      //h[k][111] = (TH1D*)(file[k]->Get("h_jet_lead_pt_300"));
      //h[k][112] = (TH1D*)(file[k]->Get("h_jet_lead_pt_500"));
      //h[k][113] = (TH1D*)(file[k]->Get("h_jet_lead_pt_700"));
      
      
      h[k][121] = (TH1D*)(file[k]->Get("h_ZpT_log"));
      h[k][114] = (TH1D*)(file[k]->Get("h_IMemnocuts"));
      h[k][115] = (TH1D*)(file[k]->Get("h_mETemouter"));
      h[k][116] = (TH1D*)(file[k]->Get("h_mETemouter2"));
      h[k][117] = (TH1D*)(file[k]->Get("h_mETeminner"));
      h[k][118] = (TH1D*)(file[k]->Get("h_mETeminner2"));
      h[k][119] = (TH1D*)(file[k]->Get("h_ZpT"));
      h[k][120] = (TH1D*)(file[k]->Get("h_ZpT"));
     
      h[k][150] = (TH1D*)(file[k]->Get("h_InvMass"));
      
      h[k][122] = (TH1D*)(file[k]->Get("h_mETmm"));
      h[k][123] = (TH1D*)(file[k]->Get("h_mETee"));
      h[k][124] = (TH1D*)(file[k]->Get("h_mETem"));
      h[k][125] = (TH1D*)(file[k]->Get("h_numJets"));
      h[k][126] = (TH1D*)(file[k]->Get("h_numJets"));
      h[k][127] = (TH1D*)(file[k]->Get("h_numJets"));
      h[k][128] = (TH1D*)(file[k]->Get("h_numJets"));
	
      h[k][129] = (TH1D*)(file[k]->Get("h_mET"));
      h[k][130] = (TH1D*)(file[k]->Get("h_mET"));
      h[k][131] = (TH1D*)(file[k]->Get("h_mET"));
      h[k][132] = (TH1D*)(file[k]->Get("h_mET"));
      h[k][133] = (TH1D*)(file[k]->Get("h_mET"));
      h[k][134] = (TH1D*)(file[k]->Get("h_mET"));
      h[k][135] = (TH1D*)(file[k]->Get("h_mET"));
      
      h[k][136] = (TH1D*)(file[k]->Get("h_extra_muons"));
      h[k][137] = (TH1D*)(file[k]->Get("h_extra_muons"));
      h[k][138] = (TH1D*)(file[k]->Get("h_extra_muons"));
      h[k][139] = (TH1D*)(file[k]->Get("h_extra_muons"));
      h[k][140] = (TH1D*)(file[k]->Get("h_extra_muons"));
      h[k][141] = (TH1D*)(file[k]->Get("h_extra_muons"));
      h[k][142] = (TH1D*)(file[k]->Get("h_extra_muons"));
      
      h[k][143] = (TH1D*)(file[k]->Get("h_InvMassmminner"));
      h[k][144] = (TH1D*)(file[k]->Get("h_InvMassmmouter"));
      h[k][145] = (TH1D*)(file[k]->Get("h_InvMasseeinner"));
      h[k][146] = (TH1D*)(file[k]->Get("h_InvMasseeouter"));
      h[k][147] = (TH1D*)(file[k]->Get("h_InvMasseminner"));
      h[k][148] = (TH1D*)(file[k]->Get("h_InvMassemouter"));
      h[k][149] = (TH1D*)(file[k]->Get("h_extra_electrons"));
      // cout<<"hello";
      //  h[k][151] = (TH1D*)(file[k]->Get("vertex countmm"));
      //placeholder below
        h[k][151] = (TH1D*)(file[k]->Get("vertex countmm"));
      h[k][152] = (TH1D*)(file[k]->Get("h_InvMasseminner2"));
      h[k][153] = (TH1D*)(file[k]->Get("vertex countee"));
    
      h[k][154] = (TH1D*)(file[k]->Get("h_InvMassnone"));
    
      h[k][155] = (TH1D*)(file[k]->Get("h_InvMassnone"));
      h[k][156] = (TH1D*)(file[k]->Get("h_InvMassnone"));
        h[k][157] = (TH1D*)(file[k]->Get("h_InvMassnone"));
	h[k][158] = (TH1D*)(file[k]->Get("h_InvMassnone"));
     h[k][159] = (TH1D*)(file[k]->Get("h_numJetseminner2"));
     h[k][160] = (TH1D*)(file[k]->Get("h_numJetsemouter"));
     h[k][161] = (TH1D*)(file[k]->Get("h_numJetseeinner"));
     h[k][162] = (TH1D*)(file[k]->Get("h_numJetseeouter"));
     h[k][163] = (TH1D*)(file[k]->Get("h_numJetsmminner")); 
     h[k][164] = (TH1D*)(file[k]->Get("h_numJetsmmouter"));
     h[k][165] = (TH1D*)(file[k]->Get("h_InvMassnone"));
     h[k][166] = (TH1D*)(file[k]->Get("h_InvMassnone"));
     h[k][167] = (TH1D*)(file[k]->Get("h_InvMassnone"));
     h[k][168] = (TH1D*)(file[k]->Get("h_numJetsee"));
     h[k][169] = (TH1D*)(file[k]->Get("h_InvMassnone"));
     h[k][170] = (TH1D*)(file[k]->Get("h_METSignificance2"));
     h[k][171] = (TH1D*)(file[k]->Get("h_X"));
      h[k][172] = (TH1D*)(file[k]->Get("h_X"));
      h[k][173] = (TH1D*)(file[k]->Get("h_ZpTCutmm"));
      h[k][174] = (TH1D*)(file[k]->Get("h_ZpTCutee"));
      h[k][175] = (TH1D*)(file[k]->Get("h_ZpTCutem"));
      h[k][176] = (TH1D*)(file[k]->Get("h_X"));
      h[k][177] = (TH1D*)(file[k]->Get("h_lepton1pt"));
      h[k][178] = (TH1D*)(file[k]->Get("h_lepton2pt"));
      h[k][179] = (TH1D*)(file[k]->Get("h_leptonbothpt"));
      h[k][180] = (TH1D*)(file[k]->Get("h_lepton1smalleta"));
      h[k][181] = (TH1D*)(file[k]->Get("h_lepton2smalleta"));
      h[k][182] = (TH1D*)(file[k]->Get("h_leptonsmalleta"));
      h[k][183] = (TH1D*)(file[k]->Get("h_Cutem"));
      h[k][184] = (TH1D*)(file[k]->Get("h_METSignificanceem"));
      h[k][185] = (TH1D*)(file[k]->Get("h_Sumptem"));
      h[k][186] = (TH1D*)(file[k]->Get("h_numJetsem"));
      h[k][187] = (TH1D*)(file[k]->Get("h_vertexProbem"));
      h[k][188] = (TH1D*)(file[k]->Get("h_acolem"));
      h[k][189] = (TH1D*)(file[k]->Get("h_InvMassbelow"));
      h[k][190] = (TH1D*)(file[k]->Get("h_InvMassabove"));
      h[k][191] = (TH1D*)(file[k]->Get("h_Cutmm"));
      h[k][192] = (TH1D*)(file[k]->Get("h_Cutee"));
      h[k][193] = (TH1D*)(file[k]->Get("h_InvMassem"));
      h[k][194] = (TH1D*)(file[k]->Get("h_InvMassmm"));
      h[k][195] = (TH1D*)(file[k]->Get("h_InvMassee"));
      h[k][196] = (TH1D*)(file[k]->Get("h_leadptmm"));
      h[k][197] = (TH1D*)(file[k]->Get("h_subptmm"));
      h[k][198] = (TH1D*)(file[k]->Get("h_METSignificancemm"));
      h[k][199] = (TH1D*)(file[k]->Get("h_acolmm"));
      h[k][200] = (TH1D*)(file[k]->Get("h_Sumptmm"));
      h[k][201] = (TH1D*)(file[k]->Get("h_vertexProbmm"));
      h[k][202] = (TH1D*)(file[k]->Get("h_ZpTmm"));
      h[k][203] = (TH1D*)(file[k]->Get("h_ZpTee"));
      h[k][204] = (TH1D*)(file[k]->Get("h_liso2mm"));
      h[k][205] = (TH1D*)(file[k]->Get("h_leta1mm"));
      h[k][206] = (TH1D*)(file[k]->Get("h_leta2mm"));
      h[k][207] = (TH1D*)(file[k]->Get("h_lphi1mm"));
      h[k][208] = (TH1D*)(file[k]->Get("h_lphi2mm"));
      h[k][209] = (TH1D*)(file[k]->Get("h_numJetsmm"));
      h[k][210] = (TH1D*)(file[k]->Get("h_leadptee"));
      h[k][211] = (TH1D*)(file[k]->Get("h_subptee"));
      h[k][212] = (TH1D*)(file[k]->Get("h_leta1ee"));
      h[k][213] = (TH1D*)(file[k]->Get("h_leta2ee"));
      h[k][214] = (TH1D*)(file[k]->Get("h_lphi1ee"));
      h[k][215] = (TH1D*)(file[k]->Get("h_lphi2ee"));
            h[k][216] = (TH1D*)(file[k]->Get("h_InvMassnone"));
	    h[k][217] = (TH1D*)(file[k]->Get("h_InvMassee2"));
  //  h[k][200] = (TH2D*)(file[k]->Get("h_InvmassvsMET"));
   cout<<"we got this far"<<endl;
      for(int i=100; i<218; i++)
	{ h[k][i]->Sumw2();
	    cout<<"we got this far"<<i<<endl;
	}
    }
  //     has = (TH1D*)(file[0]->Get("h_phiStarMC"));
  //     has2 = (TH1D*)(file[0]->Get("h_phiStarMC2"));
  //     has3 = (TH1D*)(file[0]->Get("h_phiStarMC2"));

      //  h[0][177]->Sumw2();
      // h[0][178]->Sumw2();

    
      MET = (TH2D*)(file[0]->Get("h_InvmassvsMET"));
      MET1 = (TH2D*)(file[1]->Get("h_InvmassvsMET"));
      MET->Sumw2();
      MET1->Sumw2();
  cout<<"asdfasdf"<<endl;
  //cout<<h[1][149]->Integral(0,76)<<endl;
}

//
//
//


void NickZemu::scale_histogram(int iMC, int i_h)
{
  
  
}
//
//
//
double NickZemu::get_bin(TH1D* _h, int ibin)
{
  return(_h->GetBinContent(ibin));    //ibin = 1 by default
}

//
//
//

void NickZemu::print_summary()        //ugly but it should work...
{
  int k = 0, jindex = 183, bindex = 184;
  double 
    num_events = get_bin(h[k][0]),     num_trilep = get_bin(h[k][1]),
    num_no_lep = get_bin(h[k][2]),     num_pass_ee = get_bin(h[k][3]),    
    num_pass_mm = get_bin(h[k][4]),    num_pass_em = get_bin(h[k][5]),    
    num_pass_ej = get_bin(h[k][6]),    num_pass_mj = get_bin(h[k][7]),
    num_good_ee = get_bin(h[k][8]),    num_good_mm = get_bin(h[k][9]),
    num_good_em = get_bin(h[k][10]),   
    num_good_ej = get_bin(h[k][11]),   num_good_mj = get_bin(h[k][12]),
    num_uds_g   = get_bin(h[k][13]),   num_unfound = get_bin(h[k][14]),
    num_0flavor = get_bin(h[k][15]),   
    num_h1vetos = get_bin(h[k][16]),   num_h2vetos = get_bin(h[k][17]),
    num_Zvetos = get_bin(h[k][18]),    num_3jvetos = get_bin(h[k][19]),
    num_2jvetos = get_bin(h[k][20]),   num_bjvetos = get_bin(h[k][21]),
    num_ivetos = get_bin(h[k][22]),    num_dzvetos = get_bin(h[k][23]),
    num_m1vetos = get_bin(h[k][24]),   num_m2vetos = get_bin(h[k][25]),
    num_ss_ee = get_bin(h[k][26]),     num_ss_mm = get_bin(h[k][27]),
    num_ss_em = get_bin(h[k][28]),     num_3b_xx = get_bin(h[k][29]),
    num_2b_ee = get_bin(h[k][30]),     num_2b_mm = get_bin(h[k][31]),
    num_2b_em = get_bin(h[k][32]),
    num_2b_ej = get_bin(h[k][33]),     num_2b_mj = get_bin(h[k][34]),
    num_miss_b = get_bin(h[k][35]),
    num_unkn_b = get_bin(h[k][36]),    num_true_b = get_bin(h[k][37]),
    num_falseb = get_bin(h[k][38]),    num_unkn_j = get_bin(h[k][39]);

  //13-17 are temp/buffer histograms

  double 
    b_udsg = get_bin(h[k][bindex], 3) + get_bin(h[k][bindex], 4),
    j_udsg = get_bin(h[k][jindex], 3) + get_bin(h[k][jindex], 4),
    jtotal = h[k][jindex]->Integral(),
    b_eff = get_bin(h[k][bindex], 1)/get_bin(h[k][jindex], 1),
    b_mis = b_udsg / j_udsg;
  
  double  
    num_pass1 = num_events - (num_no_lep + num_trilep),
    num_dilep = num_pass_ee + num_pass_mm + num_pass_em,
    num_l_jet = num_pass_ej + num_pass_mj,
    num_pass2 = num_dilep + num_l_jet,
    good_dilep = num_good_ee + num_good_mm + num_good_em,
    good_l_jet = num_good_ej + num_good_mj,
    num_2b_ll = num_2b_ee + num_2b_em + num_2b_mm,
    num_2b_lj = num_2b_ej + num_2b_mj,
    yield_ee = num_2b_ee/num_pass_ee,
    yield_mm = num_2b_mm/num_pass_mm,
    yield_em = num_2b_em/num_pass_em,
    yield_ej = num_2b_ej/num_pass_ej,
    yield_mj = num_2b_mj/num_pass_mj,
    yield_ll = num_2b_ll/num_dilep,
    yield_lj = num_2b_lj/num_l_jet,
    num_os_ee = num_good_ee - num_ss_ee,
    num_os_mm = num_good_mm - num_ss_mm,
    num_os_em = num_good_em - num_ss_em,
    efficiency = num_true_b / (num_true_b + num_miss_b),
    missedrate = num_falseb / num_uds_g;
    //missedrate = num_falseb / (num_falseb + num_true_b);

  cout 
  << endl << "============================================"
  << endl << "    total # of events = " << num_events
  << endl << "# of no lepton events = " << num_no_lep
  << endl << "# of trilepton events = " << num_trilep
  << endl << "   # of passed events = " << num_pass1    
  //<< " [" << num_pass2 << "]" 
  << endl << "--------------------------------------------"
  << endl << "# of dileps dz vetoed = " << num_dzvetos
  << endl << " # of < 2 jets vetoed = " << num_2jvetos
  << endl << " # of lep-jets vetoed = " << num_ivetos
  << endl << " # of no b-jet vetoed = " << num_bjvetos
  << endl << "   # of HT l+j vetoed = " << num_h1vetos
  << endl << "    # of HT ll vetoed = " << num_h2vetos
  << endl << "  # of MET l+j vetoed = " << num_m1vetos
  << endl << "   # of MET ll vetoed = " << num_m2vetos
  << endl << "# of ll Z mass vetoed = " << num_Zvetos
  << endl << "  # of < 4 l+j vetoed = " << num_3jvetos
  << endl << "  # of 3 b-jet vetoed = " << num_3b_xx;
  
  if(num_pass1 != num_pass2)  cout << endl << "!!!" << endl;
  
  cout 
    << endl << "============================================"
  << endl << "dilep\t OS \t SS \t ratio"
  << endl << "--------------------------------------------"
  << endl << " ee  \t" << num_os_ee << "\t" << num_ss_ee 
  << "\t" << num_ss_ee / num_os_ee
  << endl << " mm  \t" << num_os_mm << "\t" << num_ss_mm 
  << "\t" << num_ss_mm / num_os_mm
  << endl << " em  \t" << num_os_em << "\t" << num_ss_em 
  << "\t" << num_ss_em / num_os_em;

  cout 
  << endl << "============================================"
  << endl << "type \t pass \t good \t [2b] \t  yield"
  << endl << "--------------------------------------------"
  << endl << " ee  \t" << num_pass_ee << "\t" << num_good_ee 
	<< "\t" << num_2b_ee << "\t" << yield_ee << endl
  << " mm  \t" << num_pass_mm << "\t" << num_good_mm 
	<< "\t" << num_2b_mm << "\t" << yield_mm << endl
  << " em  \t" << num_pass_em << "\t" << num_good_em 
	<< "\t" << num_2b_em << "\t" << yield_em << endl
  << "--------------------------------------------" << endl
  << "dilep\t" << num_dilep   << "\t" << good_dilep  
	<< "\t" << num_2b_ll << "\t" << yield_ll << endl
  << "--------------------------------------------" << endl
  << " e+j \t" << num_pass_ej << "\t" << num_good_ej 
	<< "\t" << num_2b_ej << "\t" << yield_ej << endl
  << " m+j \t" << num_pass_mj << "\t" << num_good_mj 
	<< "\t" << num_2b_mj << "\t" << yield_mj << endl
  << "--------------------------------------------" << endl
  << " l+j \t" << num_l_jet   << "\t" << good_l_jet  
	<< "\t" << num_2b_lj << "\t" << yield_lj;
  
  cout
  << endl << "============================================"
  << endl << "# of unknowns for non-b-jet = " << num_unkn_j
  << endl << "    # of unknowns for b-jet = " << num_unkn_b
  << endl << "# of non-b-jets w/ missed b = " << num_miss_b
  << endl << "     # of b-jets w/ false b = " << num_falseb
  << endl << "      # of b-jets w/ true b = " << num_true_b
  << endl << "        # of tasteless jets = " << num_0flavor
  << endl << "        # of unmatched jets = " << num_unfound
  << endl << " # of jets for flavor match = " << jtotal
  << endl << "--------------------------------------------"
  << endl << "    pseudo b-jet efficiency = " << efficiency
  << endl << "       alt b-jet efficiency = " << b_eff
  << endl << "   pseduo missed b-tag rate = " << missedrate
  << endl << "      alt missed b-tag rate = " << b_mis
  << endl << "============================================"

  << endl << endl;
  
    
  k = 1, jindex = 183, bindex = 184; 
    num_events = get_bin(h[k][0]),     num_trilep = get_bin(h[k][1]),
    num_no_lep = get_bin(h[k][2]),     num_pass_ee = get_bin(h[k][3]),    
    num_pass_mm = get_bin(h[k][4]),    num_pass_em = get_bin(h[k][5]),    
    num_pass_ej = get_bin(h[k][6]),    num_pass_mj = get_bin(h[k][7]),
    num_good_ee = get_bin(h[k][8]),    num_good_mm = get_bin(h[k][9]),
    num_good_em = get_bin(h[k][10]),   
    num_good_ej = get_bin(h[k][11]),   num_good_mj = get_bin(h[k][12]),
    num_uds_g   = get_bin(h[k][13]),   num_unfound = get_bin(h[k][14]),
    num_0flavor = get_bin(h[k][15]),   
    num_h1vetos = get_bin(h[k][16]),   num_h2vetos = get_bin(h[k][17]),
    num_Zvetos = get_bin(h[k][18]),    num_3jvetos = get_bin(h[k][19]),
    num_2jvetos = get_bin(h[k][20]),   num_bjvetos = get_bin(h[k][21]),
    num_ivetos = get_bin(h[k][22]),    num_dzvetos = get_bin(h[k][23]),
    num_m1vetos = get_bin(h[k][24]),   num_m2vetos = get_bin(h[k][25]),
    num_ss_ee = get_bin(h[k][26]),     num_ss_mm = get_bin(h[k][27]),
    num_ss_em = get_bin(h[k][28]),     num_3b_xx = get_bin(h[k][29]),
    num_2b_ee = get_bin(h[k][30]),     num_2b_mm = get_bin(h[k][31]),
      num_2b_em = get_bin(h[k][32]),
    num_2b_ej = get_bin(h[k][33]),     num_2b_mj = get_bin(h[k][34]),
    num_miss_b = get_bin(h[k][35]),
    num_unkn_b = get_bin(h[k][36]),    num_true_b = get_bin(h[k][37]),
    num_falseb = get_bin(h[k][38]),    num_unkn_j = get_bin(h[k][39]);

  //13-17 are temp/buffer histograms

  
    b_udsg = get_bin(h[k][bindex], 3) + get_bin(h[k][bindex], 4),
      j_udsg = get_bin(h[k][jindex], 3) + get_bin(h[k][jindex], 4),
    jtotal = h[k][jindex]->Integral(),
    b_eff = get_bin(h[k][bindex], 1)/get_bin(h[k][jindex], 1),
    b_mis = b_udsg / j_udsg;
    
    num_pass1 = num_events - (num_no_lep + num_trilep),
    num_dilep = num_pass_ee + num_pass_mm + num_pass_em,
    num_l_jet = num_pass_ej + num_pass_mj,
    num_pass2 = num_dilep + num_l_jet,
    good_dilep = num_good_ee + num_good_mm + num_good_em,
    good_l_jet = num_good_ej + num_good_mj,
    num_2b_ll = num_2b_ee + num_2b_em + num_2b_mm,
    num_2b_lj = num_2b_ej + num_2b_mj,
    yield_ee = num_2b_ee/num_pass_ee,
    yield_mm = num_2b_mm/num_pass_mm,
    yield_em = num_2b_em/num_pass_em,
    yield_ej = num_2b_ej/num_pass_ej,
    yield_mj = num_2b_mj/num_pass_mj,
    yield_ll = num_2b_ll/num_dilep,
    yield_lj = num_2b_lj/num_l_jet,
    num_os_ee = num_good_ee - num_ss_ee,
    num_os_mm = num_good_mm - num_ss_mm,
    num_os_em = num_good_em - num_ss_em,
    efficiency = num_true_b / (num_true_b + num_miss_b),
    missedrate = num_falseb / num_uds_g;
    //missedrate = num_falseb / (num_falseb + num_true_b);

  cout 
  << endl << "============================================"
  << endl << "    total # of events = " << num_events
  << endl << "# of no lepton events = " << num_no_lep
  << endl << "# of trilepton events = " << num_trilep
  << endl << "   # of passed events = " << num_pass1    
  //<< " [" << num_pass2 << "]" 
  << endl << "--------------------------------------------"
  << endl << "# of dileps dz vetoed = " << num_dzvetos
  << endl << " # of < 2 jets vetoed = " << num_2jvetos
  << endl << " # of lep-jets vetoed = " << num_ivetos
  << endl << " # of no b-jet vetoed = " << num_bjvetos
  << endl << "   # of HT l+j vetoed = " << num_h1vetos
  << endl << "    # of HT ll vetoed = " << num_h2vetos
  << endl << "  # of MET l+j vetoed = " << num_m1vetos
  << endl << "   # of MET ll vetoed = " << num_m2vetos
  << endl << "# of ll Z mass vetoed = " << num_Zvetos
  << endl << "  # of < 4 l+j vetoed = " << num_3jvetos
  << endl << "  # of 3 b-jet vetoed = " << num_3b_xx;
  
  if(num_pass1 != num_pass2)  cout << endl << "!!!" << endl;
  
  cout 
    << endl << "============================================"
  << endl << "dilep\t OS \t SS \t ratio"
  << endl << "--------------------------------------------"
  << endl << " ee  \t" << num_os_ee << "\t" << num_ss_ee 
  << "\t" << num_ss_ee / num_os_ee
  << endl << " mm  \t" << num_os_mm << "\t" << num_ss_mm 
  << "\t" << num_ss_mm / num_os_mm
  << endl << " em  \t" << num_os_em << "\t" << num_ss_em 
  << "\t" << num_ss_em / num_os_em;

  cout 
  << endl << "============================================"
  << endl << "type \t pass \t good \t [2b] \t  yield"
  << endl << "--------------------------------------------"
  << endl << " ee  \t" << num_pass_ee << "\t" << num_good_ee 
	<< "\t" << num_2b_ee << "\t" << yield_ee << endl
  << " mm  \t" << num_pass_mm << "\t" << num_good_mm 
	<< "\t" << num_2b_mm << "\t" << yield_mm << endl
  << " em  \t" << num_pass_em << "\t" << num_good_em 
	<< "\t" << num_2b_em << "\t" << yield_em << endl
  << "--------------------------------------------" << endl
  << "dilep\t" << num_dilep   << "\t" << good_dilep  
	<< "\t" << num_2b_ll << "\t" << yield_ll << endl
  << "--------------------------------------------" << endl
  << " e+j \t" << num_pass_ej << "\t" << num_good_ej 
	<< "\t" << num_2b_ej << "\t" << yield_ej << endl
  << " m+j \t" << num_pass_mj << "\t" << num_good_mj 
	<< "\t" << num_2b_mj << "\t" << yield_mj << endl
  << "--------------------------------------------" << endl
  << " l+j \t" << num_l_jet   << "\t" << good_l_jet  
	<< "\t" << num_2b_lj << "\t" << yield_lj;
  
  cout
  << endl << "============================================"
  << endl << "# of unknowns for non-b-jet = " << num_unkn_j
  << endl << "    # of unknowns for b-jet = " << num_unkn_b
  << endl << "# of non-b-jets w/ missed b = " << num_miss_b
  << endl << "     # of b-jets w/ false b = " << num_falseb
  << endl << "      # of b-jets w/ true b = " << num_true_b
  << endl << "        # of tasteless jets = " << num_0flavor
  << endl << "        # of unmatched jets = " << num_unfound
  << endl << " # of jets for flavor match = " << jtotal
  << endl << "--------------------------------------------"
  << endl << "    pseudo b-jet efficiency = " << efficiency
  << endl << "       alt b-jet efficiency = " << b_eff
  << endl << "   pseduo missed b-tag rate = " << missedrate
  << endl << "      alt missed b-tag rate = " << b_mis
  << endl << "============================================"

  << endl << endl;
  
    
}





//
//
//
void NickZemu::make_plots()
  
{
  TString name_c;
  
  for(int i=0;i< N_HIST; i++)
    {
	
      
      name_c = "c";
      name_c += i;
      
      c[i] = new TCanvas(name_c, name_c, 25, 25, 800, 600);
      
      plot(i);
    
    }
  }

//
//
//




void NickZemu::plot(int i)
{
  
  int k=0;
  
  c[i]->cd(1);
  c[i]->SetFillStyle(0);
  c[i]->SetFillColor(kWhite);
  
  h[k][i]->Draw("hist");
}

//
//
//

void NickZemu::n_one_plots()
{
  int FirstValue, SecondValue, ThirdValue, Yield, k=0;

   double bpt1=0, bpt2=0, cpt1=0, cpt2=0, ratio=0, answer=0;
  // for(int j=0; j<N_MC; j++)
  //  {
  //    cout<<"Cut yields. Jet Pt>200; lepton pT>150; MET>100"<<endl;
  //    cout<<"Results for "<<FILE_NAME[j]<<":"<<endl<<endl;
  //    cout<<"Failed   Passed   Total b/cuts  Yield"<<endl;
  //    for(int i=320; i<323; i++)
  //	{
  //	  FirstValue=0, SecondValue=0, ThirdValue=0, Yield=0;
  //	  FirstValue=h[j][i]->GetBinContent(2);
  //	  SecondValue=h[j][i]->GetBinContent(4);
  //	  ThirdValue=h[j][i]->GetBinContent(6);
  //	  
  //	  Yield=(SecondValue*100)/ThirdValue;
  //	  
  //	  cout<<FirstValue<<"     "<<SecondValue<<"      "<<ThirdValue<<"      "<<Yield<<"%"<<endl;
  //	}
  //   cout<<endl; 
  //  }
  // for(int j=323;j<326; j++)
  //  {answer=0;
  //    for(int i=0;i<140; i++)
  //	{
  //	  bpt1=0, bpt2=0, cpt1=0, cpt2=0;
  //	  bpt1=  h[0][j]->Integral(i,140);
  //	  bpt2= h[1][j]->Integral(i,140);
  //	  

	  // 5000000 events in sample, assume 20 fb^-1 of data and ttbar crosssection is 158 picobarns give scaling factor of .632
  //	  cpt1= bpt1*.632;

	  
	  // 107962 tprime events with 20 fb^-1 and cross section of 1 picobarn gives scaling factor of .185
  //	  cpt2=bpt2*.185;

   
	  //ratio
  //	  ratio=cpt2/sqrt(cpt1);
	  
  //	  cout<<endl<<(i*5)<<"       "<<cpt1<<"         "<<cpt2<<"       "<<ratio<<endl;
	  
  //	  if (ratio>answer && ratio<20)
  //	    {answer = ratio;
  //	      k=i;
  //	    }
  //	  
  //	}
      //      cout<<"The highest value obtained was "<<answer<<"at a cut of "<<k*5<<endl;
	
  //  }
  //
  //for(int i=323; i<=325; i++)
  //    {
  //	h[0][i]->Rebin(4);
  //	h[1][i]->Rebin(4);
  //	h[0][i]->Scale(.632);
  //	h[1][i]->Scale(.185);
  //   }
   
  //  THStack *hs1 = new THStack("hs","test");

  //  h[0][323]->GetXaxis()->SetTitle("b-Jet pt given; lpt>150, mET>100");
  // h[1][323]->SetFillColor(2);
  // h[0][323]->Draw("hist");
   //  hs1->Add(h[1][323]);

  cout<<"we got this far"<<endl;
   //  hs1->Add(h[0][323]);
  //hs1->GetHistogram()->GetXaxis()->SetTitle("b-Jet pt give; lpt>150, mET>100");
  
  //hs1->Draw("hist");
  //A .87622
  //    double x=18.798*.945;
  //  data sample C = 7.017
  //  data sample B=4.412
  //dy iso sample = 8.154
  //dy normal sample =8.562
  //CE = 7.055
  //PHMM=1.711
cout<<h[7][194]->Integral()<<" , "<<h[0][194]->Integral()<<endl;
 double y=19.7;
//        double y=0.87622;
//double y=6.98;
  double x=y;
  // 8.306
  double a=x/8.307, b=x/685.0, c= x/28.02, d= x/575., e=x/301.00, f=x/167.0, g=x/12.55*.9753, l=x/24.78;
  //      double a=x/8.562, b=1005.0/48410000000000.0, c= x/28.1, d= x/542.76, e=x/284.00, f=x/180.7, g=x/1.511;
  //   double a=x/7.252, b=x/8000000000.154, c= x/28.1, d= x/542.76, e=x/284.00, f=x/180.7, g=x/1.674;
  //double a=x/12.01, b=x/8000000000.154, c= x/28.167, d= x/544.0, e=x/288.60, f=x/181.44, g=x/1.678;
       double  z=1.0;
    for(int i=100; i<218;i++)
      {  
	h[0][i]->Scale(a);
	h[1][i]->Scale(b);
	h[2][i]->Scale(c);
	h[3][i]->Scale(d);
	h[4][i]->Scale(e);
	h[5][i]->Scale(f);
	h[6][i]->Scale(g);
	h[7][i]->Scale(z);
	h[8][i]->Scale(l);
	
      }  
    // has->Scale(a);
	// has2->Scale(a);
     //     has3->Scale(a);
    //    ratio=  h[7][197]->Integral(0,1040)/(h[0][197]->Integral(0,1040)+h[2][197]->Integral(0,1040)+h[3][197]->Integral(0,1040)+h[4][197]->Integral(0,1040)+h[5][197]->Integral(0,1040)+h[6][197]->Integral(0,1040)); //MM
 
    
   ratio=  h[7][194]->Integral(0,1040)/(h[0][194]->Integral(0,1040)+h[2][194]->Integral(0,1040)+h[3][194]->Integral(0,1040)+h[4][194]->Integral(0,1040)+h[5][194]->Integral(0,1040)); //EE
    cout<<h[7][194]->Integral()<<" , "<<(h[0][194]->Integral()+h[2][194]->Integral(0,1040)+h[3][194]->Integral(0,1040)+h[4][194]->Integral(0,1040)+h[5][194]->Integral(0,1040))<<endl;
   cout<<ratio<<endl<<endl<<endl<<endl;
   cout<<h[7][195]->Integral(0,1040)<<" , "<<(h[0][195]->Integral(0,1040)+h[2][195]->Integral(0,1040)+h[3][195]->Integral(0,1040)+h[4][195]->Integral(0,1040)+h[5][195]->Integral(0,1040))<<endl; //EE
      ratio=  h[7][195]->Integral(0,1040)/(h[0][195]->Integral(0,1040)+h[2][195]->Integral(0,1040)+h[3][195]->Integral(0,1040)+h[4][195]->Integral(0,1040)+h[5][195]->Integral(0,1040)); //EE
      //            cout<<h[7][183]->Integral(0,1040)<<endl;

   cout<<ratio<<endl<<endl<<endl<<endl;
    
   //      ratio=  has2->Integral(0,70)/has->Integral(0,70);
      //    has->Scale(ratio);
  cout<<"we got this far"<<endl; 
  ratio=  h[1][150]->Integral(0,140)/h[0][150]->Integral(0,140);
  //   cout<<ratio;

    for(int i=100; i<164;i++)
      { // h[0][i]->Scale(ratio);
      }       
    //           ratio=  h[7][151]->Integral(0,70)/h[0][151]->Integral(0,70);
    //   cout<<ratio;

     
    //   ratio=  h[1][151]->Integral(0,15)/h[2][151]->Integral(0,50);
   //  cout<<ratio;

    for(int i=151; i<153;i++)
      { // h[2][i]->Scale(ratio);
      }       
    for (int i=100; i<153;i++)
      {//h[2][i]->Scale(.1);
      }
    //h[0][317]->Scale(ratio);
 //h[0][323]->Scale(ratio);
   //ratio=  h[1][324]->Integral(0,20)/h[0][324]->Integral(0,120);
   //  cout<<ratio;
    
    //   for(int i=324; i<329;i++)
    // {h[0][i]->Scale(ratio);
    // }

    // for(int i=318; i<323;i++)
    // {h[0][i]->Scale(ratio);
    // }

  yields();

  
    cout<<"hello"<<endl;
    /*
c93 = new TCanvas("inv mass em ; ZpT>200","X; ZpT>0", 25, 25, 800, 600);

  THStack *hs94 = new THStack("","");
      c93->SetLogy(1);
    hs94->SetMinimum(1);
    h[3][193]->SetFillColor(5);
    h[4][193]->SetFillColor(5);
    h[2][193]->SetFillColor(kGreen+2);
    h[5][193]->SetFillColor(2);
    h[6][193]->SetFillColor(4);
  h[4][193]->Add(h[3][193]);
  hs94->Add(h[4][193]);
  // hs94->Add(h[6][193]);

  hs94->Add(h[5][193]);
  hs94->Add(h[2][193]);
  hs94->Add(h[6][193]);
  hs94->Add(h[0][193]);
  hs94->Draw("hist"); 
  h[1][193]->Draw("same");
hs94->GetXaxis()->SetTitle("Invariant Mass of electron muon pair");
hs94->GetYaxis()->SetTitle("Number Events");
  h[1][193]->Draw("same e");

  //  c93->Print("plots/Invmassem.eps");
  
c94 = new TCanvas("inv mass muon ; ZpT>200","X; ZpT>0", 25, 25, 800, 600);

  THStack *hs95 = new THStack("","");
      c94->SetLogy(1);
    hs95->SetMinimum(1);
    h[3][194]->SetFillColor(5);
    h[4][194]->SetFillColor(5);
    h[2][194]->SetFillColor(kGreen+2);
    h[5][194]->SetFillColor(2);
    h[6][194]->SetFillColor(4);
  h[4][194]->Add(h[3][194]);

  hs95->Add(h[5][194]);
  hs95->Add(h[2][194]);
  hs95->Add(h[6][194]);
  hs95->Add(h[4][194]);
  hs95->Add(h[0][194]);
  hs95->Draw("hist"); 
hs95->GetXaxis()->SetTitle("Invariant Mass of muon pair");
hs95->GetYaxis()->SetTitle("Number Events");
  h[7][194]->Draw("same e");

  c94->Print("plots/Invmassmm.eps");
  
c95 = new TCanvas("inv mass elect ; ZpT>200","X; ZpT>0", 25, 25, 800, 600);

  THStack *hs96 = new THStack("","");
      c95->SetLogy(1);
    hs96->SetMinimum(1);
    h[3][195]->SetFillColor(5);
    h[4][195]->SetFillColor(5);
    h[2][195]->SetFillColor(kGreen+2);
    h[5][195]->SetFillColor(2);
    h[6][195]->SetFillColor(4);
    h[4][195]->Add(h[3][195]);

  hs96->Add(h[5][195]);
  hs96->Add(h[2][195]);
  hs96->Add(h[6][195]);
  hs96->Add(h[4][195]);
  hs96->Add(h[0][195]);
  hs96->Draw("hist"); 
  h[7][195]->Draw("same e");
hs96->GetXaxis()->SetTitle("Invariant Mass of electron pair");
hs96->GetYaxis()->SetTitle("Number Events");
// h[1][195]->Draw("same e");

//  c95->Print("plots/Invmassee.eps");



      c77 = new TCanvas("jet pT ; ZpT>200","X; ZpT>0", 25, 25, 800, 600);

  //   c69->SetLogy(1);

  THStack *hs75 = new THStack("lepton 1 pt","lepton 1 pt");
  //  c77->SetLogy(1);
  h[5][177]->SetFillColor(5);
  h[4][177]->SetFillColor(4);
  h[2][177]->SetFillColor(2);
  h[3][177]->SetFillColor(3);
  hs75->Add(h[4][177]);
  hs75->Add(h[3][177]);
  hs75->Add(h[2][177]);
  hs75->Add(h[0][177]);
  hs75->Draw("hist"); 
  //c77->Print("plots/lepton1pt.eps");

        cout<<"hello"<<endl;

  c79 = new TCanvas("met significance ; ZpT>200","X; ZpT>0", 25, 25, 800, 600);

  THStack *hs80 = new THStack("","");
  //    c79->SetLogy(1);
  //  hs80->SetMinimum(1);

   h[5][184]->SetFillColor(kGreen+2);
   //   h[6][184]->SetFillColor(kGreen+2);
 // h[2][178]->SetFillColor(2);
  //  h[3][178]->SetFillColor(3);
  h[4][184]->Add(h[6][184]);
  //   h[4][184]->Add(h[5][184]);
  h[4][184]->Add(h[3][184]);
   h[4][184]->Add(h[2][184]);
  h[4][184]->Add(h[0][184]);
  //hs80->Add(h[1][184]);
  //  hs80->Add(h[4][184]);
  hs80->Add(h[5][184]);
  //hs80->Add(h[6][184]);

     	hs80->SetMaximum(60);

  hs80->Draw("hist"); 
     	hs80->SetMaximum(60);

hs80->GetXaxis()->SetTitle("mET Significance");
hs80->GetYaxis()->SetTitle("Number Events");
  h[1][184]->Draw("same e");

  //  c79->Print("plots/MetSignificance11.eps");
  c80 = new TCanvas("scalar sum pt ; ZpT>200","X; ZpT>0", 25, 25, 800, 600);

  THStack *hs81 = new THStack("Scalar Sum of pt","Scalar Sum of pt");
  //   c80->SetLogy(1);
  //hs81->SetMinimum(1);

    h[1][185]->SetFillColor(5);
//  h[4][178]->SetFillColor(4);
 // h[2][178]->SetFillColor(2);
  //  h[3][178]->SetFillColor(3);
  h[5][185]->Add(h[6][185]);
  h[5][185]->Add(h[4][185]);
  h[5][185]->Add(h[3][185]);
  h[5][185]->Add(h[2][185]);
  h[5][185]->Add(h[0][185]);
  hs81->Add(h[5][185]);

  
  hs81->Draw("hist"); 
  //  c80->Print("plots/Sumpt.eps");
  h[1][185]->Draw("same e");

  
  c81 = new TCanvas("number jets ; ZpT>200","X; ZpT>0", 25, 25, 800, 600);

  THStack *hs82 = new THStack("","");
  //    c81->SetLogy(1);
  //  hs82->SetMinimum(1);

    h[1][186]->SetFillColor(2);
    h[2][186]->SetFillColor(2);
 // h[2][178]->SetFillColor(2);
  //  h[3][178]->SetFillColor(3);
  h[5][186]->Add(h[6][186]);
  h[5][186]->Add(h[4][186]);
  h[5][186]->Add(h[3][186]);
  //h[5][186]->Add(h[2][186]);
  h[5][186]->Add(h[0][186]);

//  h[1][186]->Draw("same e");
  hs82->Add(h[2][186]);
//  h[5][186]->Add(h[0][186]);
  // hs82->Add(h[5][186]);
  hs82->Draw("hist"); 
hs82->GetXaxis()->SetTitle("Number of Jets");
hs82->GetYaxis()->SetTitle("Number Events");
  h[1][186]->Draw("same e");

  //  c81->Print("plots/NumJetsem.eps");


  c82 = new TCanvas("vertex prob ; ZpT>200","X; ZpT>0", 25, 25, 800, 600);

  THStack *hs83 = new THStack("vertex prob","vertex prob");
  //     c82->SetLogy(1);
  //  hs83->SetMinimum(1);

    h[1][187]->SetFillColor(2);
  h[2][187]->SetFillColor(3);
 // h[2][178]->SetFillColor(2);
  //  h[3][178]->SetFillColor(3);
  h[5][187]->Add(h[6][187]);
  h[5][187]->Add(h[4][187]);
  h[5][187]->Add(h[3][187]);
  hs83->Add(h[2][187]);
  h[5][187]->Add(h[0][187]);
  hs83->Add(h[5][187]);
  hs83->Draw("hist");
  h[1][187]->Draw("same e");
 
  // c82->Print("plots/vertexProbem.eps");


  c83 = new TCanvas("acol; ZpT>200","X; ZpT>0", 25, 25, 800, 600);
  THStack *hs84 = new THStack("acol","acol");
  //    c83->SetLogy(1);
  // hs84->SetMinimum(1);

    h[1][188]->SetFillColor(kGreen+2);
  h[2][188]->SetFillColor(3);
 // h[2][178]->SetFillColor(2);
  //  h[3][178]->SetFillColor(3);
  h[5][188]->Add(h[6][188]);
  h[5][188]->Add(h[4][188]);
  h[5][188]->Add(h[3][188]);
  h[5][188]->Add(h[2][188]);
  h[5][188]->Add(h[0][188]);
  hs84->Add(h[5][188]);
  hs84->Draw("hist");
  h[1][188]->Draw("same e");
 
  //  c83->Print("plots/acolem.eps");

  c84 = new TCanvas("inv mass em; ZpT>200","X; ZpT>0", 25, 25, 800, 600);
  THStack *hs85 = new THStack("inv mass em","inv mass em");
  //    c84->SetLogy(1);
  // hs85->SetMinimum(1);

    h[1][189]->SetFillColor(2);
    h[2][189]->SetFillColor(3);
 // h[2][178]->SetFillColor(2);
  //  h[3][178]->SetFillColor(3);
  h[5][189]->Add(h[6][189]);
  h[5][189]->Add(h[4][189]);
  h[5][189]->Add(h[3][189]);
  h[5][189]->Add(h[2][189]);
  h[5][189]->Add(h[0][189]);
  hs85->Add(h[5][189]);
  hs85->Draw("hist");
  h[1][189]->Draw("same e");
 
  //  c84->Print("plots/InvMassem.eps");

  
  for(double m=100; m>0; m--)
    {
      //      cout<<m<<" , "<<h[1][188]->Integral(0,m)<<" , "<<h[5][188]->Integral(0,m)<<" , "<<sqrt(h[1][188]->Integral(0,m))/sqrt(h[5][188]->Integral(0,m))<<endl;
    }

  cout<<"vertex prob"<<endl<<endl;
  for(double m=100; m>0; m--)
    {
      //      cout<<m<<" , "<<h[1][187]->Integral(0,m)<<" , "<<h[5][187]->Integral(0,m)<<" , "<<sqrt(h[1][187]->Integral(0,m))/sqrt(h[5][187]->Integral(0,m))<<endl;
    }
  //Sum pt stuff
  cout<<"SUM PT:"<<endl<<endl;
 for(int k=100; k>0; k--)
    {
      //n95  cout<<k<<" , "<<n95(m,185,true)<<endl;


    }
  
c78 = new TCanvas("inv mass peak ; ZpT>200","X; ZpT>0", 25, 25, 800, 600);

  //   c69->SetLogy(1);

  THStack *hs79 = new THStack("Invariant mass of the electron and muon","");
  
  // c78->SetLogy(1);
 h[1][183]->SetFillColor(kOrange+8);

 h[3][183]->SetFillColor(1);
    h[4][183]->SetFillColor(1);
    h[2][183]->SetFillColor(kGreen+2);
    h[5][183]->SetFillColor(2);
    h[6][183]->SetFillColor(4);
//  h[4][178]->SetFillColor(4);
 // h[2][178]->SetFillColor(2);
  //  h[3][178]->SetFillColor(3);
  h[3][183]->Add(h[4][183]);
  // hs79->Add(h[6][183]);
  hs79->Add(h[3][183]);
   hs79->Add(h[2][183]);

    hs79->Add(h[5][183]);
  hs79->Add(h[0][183]);

  hs79->Add(h[1][183]);
  hs79->Draw("hist"); 

  hs79->GetXaxis()->SetTitle("Invariant Mass of electron muon pair with BR = 1.5*10^-6");
  hs79->GetYaxis()->SetTitle("Number of events");

  //  c78->Print("plots/SIMpeak.eps");

    */  
 //sone actual n95 programing
 cout<<"n95 practice"<<endl<<endl;
 // for( int k=6; k>-1; k--)
 string tole="helloworld";
 int lll=186;

 
 dataComparisons(194, "MMInvmass", "M_mm");
 dataComparisons(196, "MMleadpt", "lead muon pt");
 dataComparisons(197, "MMsubpt", "subleading muon pt");
 dataComparisons(205, "MMleta1", "leading muon eta");
 dataComparisons(206, "MMleta2", "subleading moun eta");
 dataComparisons(207, "MMlphi1", "");
 dataComparisons(208, "MMlphi2", "");
 dataComparisons(209, "MMNJets", "Number of Jets");
 dataComparisons(151, "MMVtxCount", "Number of Primary Vertices");
 dataComparisons(202, "ZpTmm", "Z_pT");
 dataComparisons(122, "mETmm", "mET");

 dataComparisons(195, "EEInvmass", "M_ee");
 dataComparisons(210, "EEleadpt", "lead electron pt");
 dataComparisons(211, "EEsubpt", "subleading electron pt");
 dataComparisons(212, "EEleta1", "leading electron eta");
 dataComparisons(213, "EEleta2", "subleading electron eta");
 dataComparisons(214, "EElphi1", "");
 dataComparisons(215, "EElphi2", "");
 dataComparisons(123, "mETee", "mET");
 dataComparisons(168, "numJetsee", "number of Jets");
 dataComparisons(203, "ZpTee", "Z_pT");
 dataComparisons(153, "EEVtxCount", "Number of Primary Vertices");


  dataComparisons(114,"IMemnocuts", "M_em");
 dataComparisons(115, "mETemouter", "mET");
  dataComparisons(148,"InvMassemouter", "m_em");
 dataComparisons(116, "mETemouter2", "mET");


 dataComparisons(117, "mETeminner", "mET");
  dataComparisons(147,"InvMasseminner", "M_em");

 dataComparisons(118, "mETeminner2", "mET");
  dataComparisons(152,"InvMasseminner2", "M_em");

 dataComparisons(183, "Cutem", "M_em");
 dataComparisons(170, "mETemcut", "mET");


 dataComparisons(193, "EMInvmass", "");
 
 // dataComparisons(198, "MMMetSig", "");
  // dataComparisons(199, "MMacol", "");
  // dataComparisons(200, "MMsuMpt", "");
  // dataComparisons(201, "MMvertexProb", "");
 dataComparisons(202, "MMZpT", "");
  dataComparisons(203, "EEliso1", "");
  dataComparisons(204, "EEliso2", "");

 //dataComparisons(216,"BarreltestNick", "");
  dataComparisons(159,"NumJetseminner", "Number of Jets");
 dataComparisons(160,"NumJetsemouter", "Number of Jets");
 dataComparisons(186,"NumJetsem", "Number of Jets");
 dataComparisons(161,"NumJetseeinner", "");
 dataComparisons(162,"NumJetseeouter", "");
 dataComparisons(163,"NumJetsmminner", "");
 dataComparisons(164,"NumJetsmmouter", "");
  dataComparisons(173,"CutsZpTmm", "");
  dataComparisons(174,"CutsZpTee", "");
   dataComparisons(175,"CutsZpTem", "");
  dataComparisons(143,"InvMassmminner", "");
  dataComparisons(144,"InvMassmmouter", "");
  dataComparisons(145,"InvMasseeinner", "");
  dataComparisons(146,"InvMasseeouter", "");
  dataComparisons(198,"METSignificanceMM", "");
  dataComparisons(100,"METSignificanceEE", "");
  dataComparisons(184,"METSignificanceEM", "");


  


 dataComparisons(217, "EEInvmass2", "");
 //dataComparisons(214, "EElphi1", "");
 // dataComparisons(215, "EElphi2", "");

 dataComparisons(101, "METresmm", "");
 dataComparisons(102, "NJresmm", "");
 dataComparisons(103, "SPresmm", "");
 dataComparisons(104, "VPresmm", "");
 dataComparisons(105, "METresee", "");
 dataComparisons(106, "NJresee", "");
 dataComparisons(108, "SPresee", "");
 dataComparisons(109, "VPresee", "");
 dataComparisons(110, "METresem", "");
 dataComparisons(111, "NJresem", "");
 dataComparisons(112, "SPresem", "");
 dataComparisons(113, "VPresem", "");

 dataComparisons(124, "mETem", "");



 Ratio(0, 191, 194, "Cutmm", "Invmassmm");

 n95(lll,true, "number_jets");
  n95(184 ,true, "mETSignigicance1");
 n95(189 ,false, "invMassAbove");
 n95(190 ,true, "invMassBelow");
 n95(188 ,true, "acol");
 n95(187 ,true, "-log10vertexprob");
  n95(185 ,true, "sumptBelow");
 n95(185 ,false, "sumptAbove");
 n95(175, true, "n95ZpT"); 

 PhiStar();
 jet_reweighting();
 /*
 c201 = new TCanvas("madgraph vs powheg","this",25,25,800,600);
 c201->SetGrid();

   h[1][299] = (TH1D*)(h[7][194]->Clone());
    cout<<"hello"<<endl;
        h[1][194]->Add(h[2][194]);
       h[1][194]->Add(h[3][194]);
       h[1][194]->Add(h[4][194]);
       h[1][194]->Add(h[5][194]);
         h[1][299]->Divide(h[1][194]);
        h[1][299]->Draw("hist");

        h[1][299]->SetMaximum(1.2);
        h[1][299]->SetMinimum(0.8);

    h[1][201] = (TH1D*)(h[7][194]->Clone());
    h[1][201]->SetMarkerStyle(1);
    h[1][201]->SetFillColor(2);

        h[0][194]->Add(h[2][194]);
       h[0][194]->Add(h[3][194]);
       h[0][194]->Add(h[4][194]);
       h[0][194]->Add(h[5][194]);
       h[0][194]->Add(h[6][194]);

    h[1][201]->Divide(h[0][194]);
    h[1][201]->Draw("same e");
    h[1][201]->SetFillColor(2);
    h[1][201]->SetMaximum(1.2);
    h[1][201]->SetMinimum(0.8);
    leg = new TLegend(0.8,0.8,0.99,0.99);
 leg->AddEntry( h[1][201], "Powheg", "l");
 leg->AddEntry( h[1][299], "MadGraph", "l" );
 leg->Draw();
h[1][201]->GetXaxis()->SetTitle("Invariant Mass");
h[1][201]->GetYaxis()->SetTitle("data/MC");
 c201->Update();
*/
 //  c201->Print("plots/MMcomparison.eps");

 c200 = new TCanvas("invmassmet ; ZpT>200","X; ZpT>0", 25, 25, 800, 600);

  //   c69->SetLogy(1);

 //MET->SetFillColor();
 MET->Draw("box");
 MET1->Draw("same box");
 
 //c200->Print("plots/InvmassvsMET.eps");

 //  
 //  
  //  
  //THStack *hs2 = new THStack("hs","test");
  //h[1][324]->SetFillColor(2);
  //hs2->Add(h[1][324]);
  //  h[0][324]->SetFillColor(1);
  //hs2->Add(h[0][324]);
  // hs2->Draw("hist"); 
  // h[0][324]->GetXaxis()->SetTitle("lepton-pt given bpt>200, mET>100");
  // c2 = new TCanvas("lpt","lpt",25,25,800,600);
  //h[0][324]->Draw("hist");
  //h[1][324]->SetMarkerStyle(20);
  //h[1][324]->Draw("same");
  //c2->Print("plots/lpt.eps");
  
  //  THStack *hs3 = new THStack("hs","test");
  
  // h[0][325]->GetXaxis()->SetTitle("InvMass of the Z");
  //  h[1][325]->SetFillColor(2);
  //hs3->Add(h[1][325]);
  //hs3->Add(h[0][325]);
  //  hs3->Draw("hist");  
  // c3 = new TCanvas("mET","mET",25,25,800,600);
  // h[0][323]->Draw("hist");
  // h[1][323]->SetMarkerStyle(20);
  // h[1][323]->Draw("same");  
  // c3->Print("plots/mET.eps");
  
 
}

//
//

void NickZemu::yields()
{

  cout<<"For the mm case we have"<<endl;
  cout<<"Simulated Signal events: "<<h[1][191]->Integral()<<endl;
  cout<<"Simulated DY events: "<<h[0][191]->Integral()<<endl;
  cout<<"Simulated ttbar events: "<<h[2][191]->Integral()<<endl;
  cout<<"Simulated WW events: "<<h[5][191]->Integral()<<endl;
  cout<<"Simulated WZ events: "<<h[4][191]->Integral()<<endl;
  cout<<"Simulated ZZ events: "<<h[3][191]->Integral()<<endl;
  cout<<"Simulated WJets events: "<<h[6][191]->Integral()<<endl;
  cout<<"Data events: "<<h[7][191]->Integral()<<endl;
  cout<<"tau tau events: "<<h[8][191]->Integral()<<endl;

  cout<<"For the ee case we have"<<endl;

  cout<<"Simulated Signal events: "<<h[1][192]->Integral()<<endl;
 cout<<"Simulated DY events: "<<h[0][192]->Integral()<<endl;
  cout<<"Simulated ttbar events: "<<h[2][192]->Integral()<<endl;
  cout<<"Simulated WW events: "<<h[5][192]->Integral()<<endl;
  cout<<"Simulated WZ events: "<<h[4][192]->Integral()<<endl;
  cout<<"Simulated ZZ events: "<<h[3][192]->Integral()<<endl;
  cout<<"Simulated WJets events: "<<h[6][192]->Integral()<<endl;
  cout<<"Data events: "<<h[7][192]->Integral()<<endl;
  cout<<"tau tau events: "<<h[8][192]->Integral()<<endl;

  cout<<"For the em case we have"<<endl;

  cout<<"Simulated Signal events: "<<h[1][183]->Integral()<<endl;
 cout<<"Simulated DY events: "<<h[0][183]->Integral()<<endl;
  cout<<"Simulated ttbar events: "<<h[2][183]->Integral()<<endl;
  cout<<"Simulated WW events: "<<h[5][183]->Integral()<<endl;
  cout<<"Simulated WZ events: "<<h[4][183]->Integral()<<endl;
  cout<<"Simulated ZZ events: "<<h[3][183]->Integral()<<endl;
  cout<<"Simulated WJets events: "<<h[6][183]->Integral()<<endl;
  cout<<"tau tau events: "<<h[8][183]->Integral()<<endl;

 cout<<"For the em inner case we have"<<endl;

 cout<<"Simulated Signal events: "<<h[1][159]->Integral(0,1)<<endl;
 cout<<"Simulated DY events: "<<h[0][159]->Integral(0,1)<<endl;
 cout<<"Simulated ttbar events: "<<h[2][159]->Integral(0,1)<<endl;
 cout<<"Simulated WW events: "<<h[5][159]->Integral(0,1)<<endl;
 cout<<"Simulated WZ events: "<<h[4][159]->Integral(0,1)<<endl;
 cout<<"Simulated ZZ events: "<<h[3][159]->Integral(0,1)<<endl;
 cout<<"Simulated WJets events: "<<h[6][159]->Integral(0,1)<<endl;
 cout<<"Data events: "<<h[7][159]->Integral(0,1)<<endl;

  cout<<"tau tau events: "<<h[8][159]->Integral()<<endl;

cout<<"For the em outer case we have"<<endl;

 cout<<"Simulated Signal events: "<<h[1][148]->Integral()<<endl;
 cout<<"Simulated DY events: "<<h[0][148]->Integral()<<endl;
 cout<<"Simulated ttbar events: "<<h[2][148]->Integral()<<endl;
 cout<<"Simulated WW events: "<<h[5][148]->Integral()<<endl;
 cout<<"Simulated WZ events: "<<h[4][148]->Integral()<<endl;
 cout<<"Simulated ZZ events: "<<h[3][148]->Integral()<<endl;
 cout<<"Simulated WJets events: "<<h[6][148]->Integral()<<endl;
 cout<<"Data events: "<<h[7][148]->Integral()<<endl;
  cout<<"tau tau events: "<<h[8][148]->Integral()<<endl;

cout<<"For the ee inner case we have"<<endl;

 cout<<"Simulated Signal events: "<<h[1][145]->Integral()<<endl;
 cout<<"Simulated DY events: "<<h[0][145]->Integral()<<endl;
 cout<<"Simulated ttbar events: "<<h[2][145]->Integral()<<endl;
 cout<<"Simulated WW events: "<<h[5][145]->Integral()<<endl;
 cout<<"Simulated WZ events: "<<h[4][145]->Integral()<<endl;
 cout<<"Simulated ZZ events: "<<h[3][145]->Integral()<<endl;
 cout<<"Simulated WJets events: "<<h[6][145]->Integral()<<endl;
cout<<"Data events: "<<h[7][145]->Integral()<<endl;
  cout<<"tau tau events: "<<h[8][145]->Integral()<<endl;

cout<<"For the ee outer case we have"<<endl;

 cout<<"Simulated Signal events: "<<h[1][146]->Integral()<<endl;
 cout<<"Simulated DY events: "<<h[0][146]->Integral()<<endl;
 cout<<"Simulated ttbar events: "<<h[2][146]->Integral()<<endl;
 cout<<"Simulated WW events: "<<h[5][146]->Integral()<<endl;
 cout<<"Simulated WZ events: "<<h[4][146]->Integral()<<endl;
 cout<<"Simulated ZZ events: "<<h[3][146]->Integral()<<endl;
 cout<<"Simulated WJets events: "<<h[6][146]->Integral()<<endl;
 cout<<"Data events: "<<h[7][146]->Integral()<<endl; 
 cout<<"tau tau events: "<<h[8][146]->Integral()<<endl;

cout<<"For the mm inner case we have"<<endl;

 cout<<"Simulated Signal events: "<<h[1][143]->Integral()<<endl;
 cout<<"Simulated DY events: "<<h[0][143]->Integral()<<endl;
 cout<<"Simulated ttbar events: "<<h[2][143]->Integral()<<endl;
 cout<<"Simulated WW events: "<<h[5][143]->Integral()<<endl;
 cout<<"Simulated WZ events: "<<h[4][143]->Integral()<<endl;
 cout<<"Simulated ZZ events: "<<h[3][143]->Integral()<<endl;
 cout<<"Simulated WJets events: "<<h[6][143]->Integral()<<endl;
cout<<"Data events: "<<h[7][143]->Integral()<<endl;
  cout<<"tau tau events: "<<h[8][143]->Integral()<<endl;

cout<<"For the mm outer case we have"<<endl;

 cout<<"Simulated Signal events: "<<h[1][144]->Integral()<<endl;
 cout<<"Simulated DY events: "<<h[0][144]->Integral()<<endl;
 cout<<"Simulated ttbar events: "<<h[2][144]->Integral()<<endl;
 cout<<"Simulated WW events: "<<h[5][144]->Integral()<<endl;
 cout<<"Simulated WZ events: "<<h[4][144]->Integral()<<endl;
 cout<<"Simulated ZZ events: "<<h[3][144]->Integral()<<endl;
 cout<<"Simulated WJets events: "<<h[6][144]->Integral()<<endl;
 cout<<"Data events: "<<h[7][144]->Integral()<<endl;
 cout<<"tau tau events: "<<h[8][144]->Integral()<<endl;

}

//========================================================================//
void NickZemu::dataComparisons(int hist_num, TString name, TString xaxis)
{
  bool pass=true;
  if(hist_num==117 ||hist_num == 186 || hist_num==159 ||  hist_num==147 ||hist_num==118 ||hist_num==152 ||hist_num==170)
    {
            pass=false;
    }

  c[hist_num] = new TCanvas(name,"X; ZpT>0", 0,0,900,750);
  c[hist_num]->Range(0,0,1,1);
  //  c[hist_num]->Divide(1,2);
  // c[hist_num]->cd(1);
       gPad->SetLogy();
gStyle->SetLabelSize(0.060,"Y");
 gStyle->SetLabelSize(0.060,"X");
 gStyle->SetOptStat(0);
 

  ///////new stuff
cc_2 = new TPad("c[histnum]_2", "newpad",0.01,0.35,0.99,0.99);
cc_2->Draw();
cc_2->cd();
cc_2->SetTopMargin(0.1);
cc_2->SetBottomMargin(0.08);
cc_2->SetRightMargin(0.1);
cc_2->SetFillStyle(0);
  ////
   hs[hist_num] = new THStack(name,"");
   // hs[hist_num]->SetMinimum(10);

        h[2][hist_num]->SetFillColor(2);
        h[3][hist_num]->SetFillColor(3);
        h[4][hist_num]->SetFillColor(4);
        h[5][hist_num]->SetFillColor(5);
        h[6][hist_num]->SetFillColor(6);
	//        h[7][hist_num]->SetFillColor(7);
	/*      h[6][hist_num]->Rebin(2);
        h[0][hist_num]->Rebin(2);
        h[2][hist_num]->Rebin(2);
        h[3][hist_num]->Rebin(2);
        h[4][hist_num]->Rebin(2);
        h[5][hist_num]->Rebin(2);
        h[7][hist_num]->Rebin(2);
	*/
	  

	//   hs[hist_num]->Add(h[7][hist_num]);
	   hs[hist_num]->Add(h[2][hist_num]);
	  hs[hist_num]->Add(h[5][hist_num]);
	  hs[hist_num]->Add(h[3][hist_num]);
	  hs[hist_num]->Add(h[4][hist_num]);
	  //	   hs[hist_num]->Add(h[6][hist_num]);
    hs[hist_num]->Add(h[0][hist_num]);

    //    h[6][hist_num]->SetFillColor(3);
    //   g[hist_num]->SetLineColor(2);
    // hs[hist_num]->Add(h[6][hist_num]);
   // g[hist_num]->Draw("hist");
   if(hist_num==185)
     {
       //g[hist_num]->SetMaximum(50);
       //	hs[hist_num]->SetMaximum(50);
	//	h[1][hist_num]->SetMaximum(50);
	//	h[2][hist_num]->SetMaximum(50);
	// 	h[3][hist_num]->SetMaximum(50);
	//	h[4][hist_num]->SetMaximum(50);

     }
   //hs[hist_num]->GetXaxis()->SetTitleOffset(0.30);

   hs[hist_num]->Draw("hist");
   for(int i=0; i<9;i++)
     {
       h[i][hist_num]->Rebin(2);
h[i][hist_num]->SetLabelSize(0.0);
h[i][hist_num]->GetXaxis()->SetTitleSize(0.00);
h[i][hist_num]->GetYaxis()->SetLabelSize(0.06);
h[i][hist_num]->GetYaxis()->SetTitleSize(0.00);
//h[i][hist_num]->GetYaxis()->SetTitleOffset(0.76);
h[i][hist_num]->GetYaxis()->SetTitleOffset(0.40);

     }
   //    hs[hist_num]->SetMinimum(10);
   //  h[0][hist_num]->SetMinimum(10);
   // h[2][hist_num]->SetMinimum(10);
   // h[3][hist_num]->SetMinimum(10);
   // h[4][hist_num]->SetMinimum(10);
   // h[5][hist_num]->SetMinimum(10);
   // h[6][hist_num]->SetMinimum(10);
   // h[7][hist_num]->SetMinimum(10);
  cout<<"madeit"<<endl;

  //hs[hist_num]->GetXaxis()->SetTitle(name);
hs[hist_num]->GetYaxis()->SetTitle("Number Events");
//  g[hist_num]->Draw("same");

//to compare to data
//change false to pass
 if(pass)
   {
      h[7][hist_num]->Draw("same e");
   }
   //to compare to signal
          h[1][hist_num]->Draw("same");

  name+=".eps";
  TString filename = "plots/EE/";
  filename+=name;
    leg = new TLegend(0.70,0.50,0.95,0.95);
    // leg->AddEntry( h[1][201], "Powheg", "l");
    //  leg->AddEntry( h[7][hist_num], "data", "lep" );
 
    leg->AddEntry( h[0][hist_num], "MadGraph DY to LL", "f" );
    leg->AddEntry( h[2][hist_num], "ttbar", "f" );
    leg->AddEntry( h[3][hist_num], "ZZ", "f" );
    leg->AddEntry( h[4][hist_num], "WZ", "f" );
    leg->AddEntry( h[5][hist_num], "WW", "f" );
    //   leg->AddEntry( h[6][hist_num], "MadGraph", "l" );

 leg->Draw();
 // if(hist_num==216)
 //   {
 //     f11->WriteTObject(h[0][hist_num], "D_data");
 //    f11->WriteTObject(h[7][hist_num], "D_PH");
 //
 //   }
  
 // c[hist_num]->cd(2);
 c[hist_num]->cd();
cc_1 = new TPad("c[hist_num]_1", "newpad",0.01,0.0.01,0.99,0.30);
cc_1->Draw();
cc_1->cd();
 gStyle->SetOptStat(0);

cc_1->SetTopMargin(0.04);
cc_1->SetBottomMargin(0.3);
cc_1->SetRightMargin(0.1);
cc_1->SetFillStyle(0);
 cc_1->Update();
gStyle->SetLabelSize(0.080,"Y");
 gStyle->SetLabelSize(0.080,"X");
 gStyle->SetOptStat(0);
 gPad->SetGrid();
 gStyle->SetLabelSize(0.080,"Y");
 gStyle->SetLabelSize(0.080,"X");
 gStyle->SetOptStat(0);

 int clones=hist_num+200;
 if(true)
   {
    h[1][clones] = (TH1D*)(h[7][hist_num]->Clone());
    h[1][clones]->SetMarkerStyle(1);
    h[1][clones]->Draw("e");
    h[1][clones]->GetXaxis()->SetLabelSize(0.08);
    h[1][clones]->GetYaxis()->SetLabelSize(0.08);
    h[0][clones] = (TH1D*)(h[0][hist_num]->Clone());
        h[0][clones]->Add(h[2][hist_num]);
       h[0][clones]->Add(h[3][hist_num]);
       h[0][clones]->Add(h[4][hist_num]);
       h[0][clones]->Add(h[5][hist_num]);
       //       h[0][hist_num]->Add(h[6][hist_num]);
    h[1][clones]->Divide(h[0][clones]);

    h[1][clones]->SetMaximum(1.2);
    h[1][clones]->SetMinimum(0.8);
 h[1][clones]->SetLineWidth(1);

h[1][clones]->GetYaxis()->SetNdivisions(5);

h[1][clones]->GetXaxis()->SetTitleSize(0.12);
h[1][clones]->GetXaxis()->SetLabelSize(0.12);
h[1][clones]->GetYaxis()->SetLabelSize(0.10);
h[1][clones]->GetYaxis()->SetTitleSize(0.12);
h[1][clones]->GetYaxis()->SetTitleOffset(0.12);
h[1][clones]->GetYaxis()->SetRangeUser(0.8,1.2);
h[1][clones]->SetStats(0);
h[1][clones]->GetYaxis()->SetTitle("Data/MC");
h[1][clones]->GetXaxis()->SetTitle(xaxis);
h[1][clones]->GetYaxis()->SetTitleOffset(0.40);


h[1][clones]->SetTitle("");

   }
gStyle->SetLabelSize(0.080,"Y");
 gStyle->SetLabelSize(0.080,"X");
 gStyle->SetOptStat(0);

 c[hist_num]->Update();
 gStyle->SetOptStat(0);

gStyle->SetLabelSize(0.080,"Y");
 gStyle->SetLabelSize(0.080,"X");
  
    f11->WriteTObject(c[hist_num]);
  c[hist_num]->Print(filename);

  
}


void NickZemu::n95(int hist_num, bool below, TString name)
{
  // cout<<"hello";
  double mini;
  double background, efficiency;
  double mini2=100000000.0;
  double blast, minnum;
    double nbins =  h[0][hist_num]->GetNbinsX();
    double lowedge = h[0][hist_num]->GetXaxis()->GetXmin();
        double upedge = h[0][hist_num]->GetXaxis()->GetXmax();
      g[hist_num] = new TH1D(name,"Invariant Mass of the Z em", nbins, lowedge, upedge); 
      // cout<<"monkeu";
   for(int i=0;i<nbins; i++)
   {
   
  if(below)
    {
      background = h[0][hist_num]->Integral(0,i)+h[2][hist_num]->Integral(0,i)+h[3][hist_num]->Integral(0,i)+h[4][hist_num]->Integral(0,i)+h[5][hist_num]->Integral(0,i);
       if(h[1][hist_num]->Integral(0, i)==0)
	 efficiency = 0.001;
       else 
	 {if(h[1][hist_num]->Integral(0,i)>h[1][hist_num]->Integral())
	       efficiency=1;
	    else
	     efficiency =h[1][hist_num]->Integral(0,i)/h[1][hist_num]->Integral();
	 }
    }
 else
   {
     background = h[0][hist_num]->Integral(i,nbins)+h[2][hist_num]->Integral(i,nbins)+h[3][hist_num]->Integral(i,nbins)+h[4][hist_num]->Integral(i,nbins)+h[5][hist_num]->Integral(i,nbins);
     if(h[1][hist_num]->Integral(i,nbins)==0)
	 efficiency = 0.001;
     else 
       {if(h[1][hist_num]->Integral(i,nbins)>h[1][hist_num]->Integral())
	   efficiency=1;
     else 
	 efficiency =h[1][hist_num]->Integral(i,0)/h[1][hist_num]->Integral();
       }
   }
  if(efficiency==0)
    efficiency=0.01;
  blast = 1.64*sqrt(background)/efficiency;
  //  cout<<i<<" , "<<sqrt(background)<<" , "<<efficiency<<" , "<<blast<<endl;

  if(blast<mini2 && blast>0.0)
    { mini2=blast;
      mini=i;
    } 
    g[hist_num]->Fill(i*(upedge-lowedge)/nbins+.5*(upedge-lowedge)/nbins +lowedge, blast);
  

   }
   minnum=mini*(upedge-lowedge)/nbins +lowedge;
   cout<<"the file labled "<<name<<" has a minimum of "<<mini2<<" at "<<minnum<<endl;

   c[hist_num] = new TCanvas(name,"X; ZpT>0", 25, 25, 800, 600);
   hs[hist_num] = new THStack(name,"");
   h[0][hist_num]->Add(h[2][hist_num]);
   h[0][hist_num]->Add(h[3][hist_num]);
   h[0][hist_num]->Add(h[4][hist_num]);
   h[0][hist_num]->Add(h[5][hist_num]);
    hs[hist_num]->Add(h[0][hist_num]);
    //    h[6][hist_num]->SetFillColor(3);
   g[hist_num]->SetLineColor(2);
   // hs[hist_num]->Add(h[6][hist_num]);
   // g[hist_num]->Draw("hist");
   if(hist_num==185)
     {
       g[hist_num]->SetMaximum(50);
     	hs[hist_num]->SetMaximum(50);
     	h[1][hist_num]->SetMaximum(50);
     	h[2][hist_num]->SetMaximum(50);
     	h[3][hist_num]->SetMaximum(50);
     	h[4][hist_num]->SetMaximum(50);

     }
   hs[hist_num]->Draw("hist");
hs[hist_num]->GetXaxis()->SetTitle(name);
hs[hist_num]->GetYaxis()->SetTitle("Number Events");
 c[hist_num]->Update();
//  g[hist_num]->Draw("same");
   h[1][hist_num]->Draw("same e");
  name+=".eps";
  TString filename = "plots/";
  filename+=name;
  c[hist_num]->Print(filename);
    f11->WriteTObject(c[hist_num]);

  // return 1.64*sqrt(background)/efficiency;

}
void NickZemu::jet_reweighting()
{
  double bincontent[12];
  double denom[12]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double total[12];
  double bincontentee[12];
  double denomee[12]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double totalee[12];
  for(Int_t q=1;q<12;q++)
    {
      denom[q]=h[0][209]->GetBinContent(q)+h[2][209]->GetBinContent(q)+h[3][209]->GetBinContent(q)+h[4][209]->GetBinContent(q)+h[5][209]->GetBinContent(q);
      bincontent[q]= h[7][209]->GetBinContent(q);
      total[q]=bincontent[q]/denom[q];

      denomee[q]=h[0][168]->GetBinContent(q)+h[2][168]->GetBinContent(q)+h[3][168]->GetBinContent(q)+h[4][168]->GetBinContent(q)+h[5][168]->GetBinContent(q);
      bincontentee[q]= h[7][168]->GetBinContent(q);
      totalee[q]=bincontentee[q]/denomee[q];
      

      cout<<total[q]<<" , "<<bincontent[q]<<" , "<<totalee[q]<<" , "<<bincontentee[q]<<endl;
      
      //cout<<h[0][209]->Integral(0,12)<<endl;
    }
}
void NickZemu::PhiStar()
{
  
}
void NickZemu::Ratio(int file, int hist1, int hist2, TString name1, TString name2)
{
  double nbins =  h[0][hist1]->GetNbinsX();
 double num= h[file][hist1]->Integral(0,nbins+1);
 double den =  h[file][hist2]->Integral(0,nbins+1);

 double ratio=num/den;
 double error= ratio*sqrt((1/num)+(1/den));
  cout<<"The ratio of "<<name1<<" over "<<name2<<" is "<<ratio<<" +/- "<<error<<endl;
}

void NickZemu::n_one_plots2()
{
  int FirstValue, SecondValue, ThirdValue, Yield, k=0;

  double bpt1=0, bpt2=0, cpt1=0, cpt2=0, ratio=0, answer=0;
  for(int j=0; j<N_MC; j++)
    {
      //  cout<<"Cut yields. Jet Pt>200; lepton pT>60; MET>60"<<endl;
      //  cout<<"Results for "<<FILE_NAME[j]<<":"<<endl<<endl;
      //  cout<<"Failed   Passed   Total b/cuts  Yield"<<endl;
      for(int i=320; i<323; i++)
	{
	  FirstValue=0, SecondValue=0, ThirdValue=0, Yield=0;
	  FirstValue=h[j][i]->GetBinContent(7);
	  SecondValue=h[j][i]->GetBinContent(8);
	  ThirdValue=h[j][i]->GetBinContent(9);
	  
	  Yield=(SecondValue*100)/ThirdValue;
	  
	  //	  cout<<FirstValue<<"     "<<SecondValue<<"      "<<ThirdValue<<"      "<<Yield<<"%"<<endl;
	}
      //   cout<<endl; 
    }
  for(int j=326;j<329; j++)
    {answer=0;
      for(int i=0;i<140; i++)
	{
	  bpt1=0, bpt2=0, cpt1=0, cpt2=0;
	  bpt1=  h[0][j]->Integral(i,140);
	  bpt2=  h[1][j]->Integral(i,140);
	  

	  // 5000000 events in sample, assume 20 fb^-1 of dataand ttbar crosssection is 158 picobarns give scaling factor of .632
	  cpt1= bpt1*.632;

	  
	  // 107962 tprime events with 20 fb^-1 and cross section of 1 picobarn gives scaling factor of .185
	  cpt2=bpt2*.185;

	  
	  //ratio
	  ratio=cpt2/sqrt(cpt1);
	  
	  //  	    cout<<endl<<(i*5)<<"IS this a valid number?: "<<bpt1<<"    "<<cpt1<<"      "<<bpt2<<"      "<<cpt2<<"    "<<ratio<<endl;
	  
	  if (ratio>answer && ratio<20)
	    {answer = ratio;
	      k=i;
	    }
	  
	}
      cout<<"The highest value obtained was "<<answer<<"at a cut of "<<k*5<<endl;
	
    }

}
//
//
//

NickZemu::~NickZemu()
{
  cout <<endl << "..Ending"<<endl;
  f11->Write();
  f11->close(); 
}

//
//
//

void NickZemu()
{
  
  NickZemu s;

  s.open_files();
  s.load_histograms();
  //s.print_summary();
  s.n_one_plots();
  //s.n_one_plots2();
}

//==========================================
//end
//==========================================
