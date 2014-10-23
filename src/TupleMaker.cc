//Spots.cc -- see Spots.h for details

#include "TupleMaker.h"

//===================================================================//
//global functions
//===================================================================//

TString get_true_false(bool ibool)
{
  if(ibool)   
    return ("true");
  
  return ("false");
}
//============================================
double get_rapid(TLorentzVector pll)
{
  return pll.Rapidity();
}
//===================================================================//

double DeltaPhiX(double phi1, double phi2)
{
  double dphi = (phi1 - phi2);

  while(dphi > M_PI)      dphi -= 2.*M_PI;
  while(dphi <= -M_PI)    dphi += 2.*M_PI;
  
  return fabs(dphi);
}

//===================================================================//

double DeltaPhiX(math::XYZVector v1, math::XYZVector v2)
{
  double phi1 = v1.phi(), phi2 = v2.phi();
  
  return DeltaPhiX(phi1, phi2);
}

//===================================================================//

double DeltaRX(double eta1, double eta2, double phi1, double phi2)
{
  double deta = (eta1 - eta2), dphi = DeltaPhiX(phi1, phi2);
  
  return (sqrt(deta*deta + dphi*dphi));
}

//===================================================================//

double DeltaRX(double deta, double phi)
{
  double dphi = DeltaPhiX(phi, 0.);    //safety check
  
  return (sqrt(deta*deta + dphi*dphi));
}

//===================================================================//

double DeltaRX(math::XYZVector v1, math::XYZVector v2)
{
  double deta = (v1.eta() - v2.eta()), dphi = DeltaPhiX(v1, v2);
  
  return (sqrt(deta*deta + dphi*dphi));
}

//===================================================================//

double get_mass(math::XYZTLorentzVector pll)
{
  return pll.mass();
}

//==================================================================//
double acolinearity(math::XYZTLorentzVector lep1, math::XYZTLorentzVector lep2)
{
  double acol=0;
  TVector3 a, b;
  a.SetXYZ(lep1.X(), lep1.Y(), lep1.Z());
  b.SetXYZ(lep2.X(), lep2.Y(), lep2.Z()); 
  
  acol=M_PI-acos((a*b)/(a.Mag()*b.Mag()));
  return acol;
}

//===================================================================//

double get_pT(math::XYZTLorentzVector v)
{
  return sqrt(v.x()*v.x() + v.y()*v.y());
}

//===================================================================//

double get_pt(math::XYZVector v)
{
  return sqrt(v.x()*v.x() + v.y()*v.y());
}

//===================================================================//

double get_abs_p(math::XYZTLorentzVector p)
{
  return p.P();
}

//===================================================================//

double get_MT(double pt, double met, double dphi)
{
  return ( 2. * sqrt(pt*met) * sin(dphi/2.) );
}

//===================================================================//

double get_MT(math::XYZVector vlep, math::XYZVector vmet)
{  
  double lpt = get_pt(vlep), met = get_pt(vmet);
  
  //direct (vlep+vmet).mt() ?
  
  return get_MT(lpt, met, DeltaPhiX(vlep, vmet));
}

//===================================================================//
//begin Spots
//===================================================================//

TupleMaker::TupleMaker(const ParameterSet& config_) 
{
  cout<<"It begins"<<endl;

  setup_stopwatch();
  setup_names();
  setup_counters();
  setup_histograms(config_);
  setup_flags();

  // Nick UPDATE

  // end
  


}

//===================================================================//

void TupleMaker::setup_stopwatch()
{
  watch = new TStopwatch();
  
  cout << endl << "[start]>\t";
  
  time_stamp();
}

//===================================================================//

void TupleMaker::time_stamp()
{
  watch->Stop();  
  watch->Print();
  watch->Continue();  
}

//===================================================================//

void TupleMaker::beginJob()     //this is all done in the constructor
{

}
//===================================================================//

void TupleMaker::setup_names()
{
  for(int i = 0; i < N_HIST; ++i)
  {
    Name_h[i] = "h_";
    
    if(i < 10)          Name_h[i] += "00";
    else if(i < 100)    Name_h[i] += "0";
    
    Name_h[i] += i;
    
    //Title_h[i] = "[no title yet]";
  }
  
  Title_h[  0] = "# of Events Loaded";
  Title_h[  1] = "# of #geq 3 leptons events";
  Title_h[  2] = "# of no lepton events";
  Title_h[  3] = "# of passed ee events w/o jet & btagging";
  Title_h[  4] = "# of passed #mu#mu events w/o jet & btagging";
  Title_h[  5] = "# of passed e#mu events w/o jet & btagging";
  Title_h[  6] = "# of passed e(+jets) events w/o jet & btagging";
  Title_h[  7] = "# of passed #mu(+jets) events w/o jet & btagging";
  Title_h[  8] = "# of good dilepton (ee) events";
  Title_h[  9] = "# of good dilepton (#mu#mu) events";
  Title_h[ 10] = "# of good dilepton (e#mu) events";
  Title_h[ 11] = "# of good lepton+jets (e) events";
  Title_h[ 12] = "# of good lepton+jets (#mu) events";
  


//275->320
}

//===================================================================//

void TupleMaker::setup_histograms(const ParameterSet& config_)
{  

  //     partFlowTag_ = config_.getUntrackedParameter<edm::InputTag>("partFlowTag");


  hlTriggerResults_ = config_.getUntrackedParameter<string>("HLTriggerResults","TriggerResults");
  inputTagIsoDepElectrons_ = config_.getParameter< std::vector<edm::InputTag> >("IsoDepElectron");
   inputTagIsoDepPhotons_ = config_.getParameter< std::vector<edm::InputTag> >("IsoDepPhoton");
  // No longer needed. e/g recommendation (04/04/12)
  //  inputTagIsoValElectronsNoPFId_ = iConfig.getParameter< std::vector<edm::InputTag> >("IsoValElectronNoPF");
  inputTagIsoValElectronsPFId_   = config_.getParameter< std::vector<edm::InputTag> >("IsoValElectronPF");   
  inputTagIsoValPhotonsPFId_   = config_.getParameter< std::vector<edm::InputTag> >("IsoValPhoton");   

 
  Service<TFileService> fs;
  
  h_num_events = fs->make<TH1D>(Name_h[0], Title_h[0], 1, 0.5, 1.5);
  h_base = fs->make<TH1D>("h_base", Title_h[1], 5, 0.5, 5.5);
  h_passacc = fs->make<TH1D>("h_passacc", Title_h[2], 5, 0.5, 5.5);
  h_passrecoeff = fs->make<TH1D>("h_passrecoeff", Title_h[3], 5, 0.5, 5.5);
  h_passseleceff = fs->make<TH1D>("h_passeleceff", Title_h[4], 5, 0.5, 5.5);
  h_num_pass_em = fs->make<TH1D>(Name_h[5], Title_h[5], 6, 0.5, 6.5);
  h_num_pass_ej = fs->make<TH1D>(Name_h[6], Title_h[6], 1, 0.5, 1.5);
  h_num_pass_mj = fs->make<TH1D>(Name_h[7], Title_h[7], 1, 0.5, 1.5);
  h_num_good_ee = fs->make<TH1D>(Name_h[8], Title_h[8], 1, 0.5, 1.5);
  h_num_good_mm = fs->make<TH1D>(Name_h[9], Title_h[9], 1, 0.5, 1.5);
  h_num_good_em = fs->make<TH1D>(Name_h[10], Title_h[10], 1, 0.5, 1.5);
  h_num_good_ej = fs->make<TH1D>(Name_h[11], Title_h[11], 1, 0.5, 1.5);
  h_num_good_mj = fs->make<TH1D>(Name_h[12], Title_h[12], 200, 0, 2);
  h_uppz = fs->make<TH1D>("h_uppz","pz of up quarks", 40, 0, 4000);
  h_uppx = fs->make<TH1D>("h_uppx","Z->?", 101, -50, 50);
  h_uppy = fs->make<TH1D>("h_uppy","Z->q->?", 101, -50, 50);

  h_downpz = fs->make<TH1D>("h_downpz","pz of down quarks", 40, 0, 4000);
  h_downpx = fs->make<TH1D>("h_downpx","px of down quarks", 100, 0, 100);
  h_downpy = fs->make<TH1D>("h_downpy","Z->q->q->?", 101, -50, 50);
  
  h_MCM1pt = fs->make<TH1D>("h_MCM1pt","MC pt 1", 100, 0, 100);
  h_MCM2pt = fs->make<TH1D>("h_MCM2pt","MC pt 2", 100, 0, 100);
  h_MCM1eta = fs->make<TH1D>("h_MCM1eta","MC eta 1", 120, -3, 3);
  h_MCM2eta = fs->make<TH1D>("h_MCM2eta","MC eta 2", 120, -3, 3);
  h_MCMIM = fs->make<TH1D>("h_MCMIM","MC IM", 180, 20, 200);
  h_MCMrapid = fs->make<TH1D>("h_MCMrapid", "MC rapid", 120,-6.,6.);

  h_MCE1pt = fs->make<TH1D>("h_MCE1pt","MC pt 1", 100, 0, 100);
  h_MCE2pt = fs->make<TH1D>("h_MCE2pt","MC pt 2", 100, 0, 100);
  h_MCE1eta = fs->make<TH1D>("h_MCE1eta","MC eta 1", 120, -3, 3);
  h_MCE2eta = fs->make<TH1D>("h_MCE2eta","MC eta 2", 120, -3, 3);
  h_MCEIM = fs->make<TH1D>("h_MCEIM","MC IM", 180, 20, 200);
  h_MCErapid = fs->make<TH1D>("h_MCErapid", "MC rapid", 120,-6.,6.);

  h_MCS1pt = fs->make<TH1D>("h_MCS1pt","MC pt 1", 100, 0, 100);
  h_MCS2pt = fs->make<TH1D>("h_MCS2pt","MC pt 2", 100, 0, 100);
  h_MCS1eta = fs->make<TH1D>("h_MCS1eta","MC eta 1", 120, -3, 3);
  h_MCS2eta = fs->make<TH1D>("h_MCS2eta","MC eta 2", 120, -3, 3);
  h_MCSIM = fs->make<TH1D>("h_MCSIM","MC IM", 180, 20, 200);
  h_MCSrapid = fs->make<TH1D>("h_MCSrapid", "MC rapid", 120,-6.,6.);

  h_MCmu1dxy = fs->make<TH1D>("h_MCmu1dxy","MC mu1 dxy", 20, 0, 1.);
  h_MCmu1dz = fs->make<TH1D>("h_MCmu1dz","MC mu1 dz", 40, 0, 2.);
  h_MCmu1Global = fs->make<TH1D>("h_MCmu1Global","MC mu1 global", 2, -0.5, 1.5);
  h_MCmu1PF = fs->make<TH1D>("h_MCmu1PF","MC mu1 PF", 2, -0.5, 1.5);
  h_MCmu1Chi2 = fs->make<TH1D>("h_MCmu1Chi2","MC mu1 Chi2", 21, -0.5, 20.5);
  h_MCmu1hits = fs->make<TH1D>("h_MCmu1hits","MC mu1 hits", 21, -0.5, 20.5);
  h_MCmu1Stations = fs->make<TH1D>("h_MCmu1Stations","MC mu1 Stations", 21, -0.5, 20.5);

  h_MCmu1Phits = fs->make<TH1D>("h_MCmu1Phits","MC mu1 Phits", 21, -0.5, 20.5);
  h_MCmu1Thits = fs->make<TH1D>("h_MCmu1Thits","MC mu1 Thits", 21, -0.5, 20.5);
  h_MCmu1iso = fs->make<TH1D>("h_MCmu1iso","MC mu1 iso",100 , 0.0, 2.0);

  
  h_MCmu2dxy = fs->make<TH1D>("h_MCmu2dxy","MC mu2 dxy", 20, 0, 1.);
  h_MCmu2dz = fs->make<TH1D>("h_MCmu2dz","MC mu2 dz", 40, 0, 2.);
  h_MCmu2Global = fs->make<TH1D>("h_MCmu2Global","MC mu2 global", 2, -0.5, 1.5);
  h_MCmu2PF = fs->make<TH1D>("h_MCmu2PF","MC mu2 PF", 2, -0.5, 1.5);
  h_MCmu2Chi2 = fs->make<TH1D>("h_MCmu2Chi2","MC mu2 Chi2", 21, -0.5, 20.5);
  h_MCmu2hits = fs->make<TH1D>("h_MCmu2hits","MC mu2 hits", 21, -0.5, 20.5);
  h_MCmu2Stations = fs->make<TH1D>("h_MCmu2Stations","MC mu2 Stations", 21, -0.5, 20.5);

  h_MCmu2Phits = fs->make<TH1D>("h_MCmu2Phits","MC mu2 Phits", 21, -0.5, 20.5);
  h_MCmu2Thits = fs->make<TH1D>("h_MCmu2Thits","MC mu2 Thits", 21, -0.5, 20.5);
  h_MCmu2iso = fs->make<TH1D>("h_MCmu2iso","MC mu2 iso",100 , 0.0, 2.0);


  //



 
  //#275-->320

// ntuple = new TNtuple("ntuple","datafromasciifile","muonpt1:muonpt2:muoneta1:muoneta2:numberjets:leadjetpt:2ndjetpt:3rdjetpt:4thjetpt:5thjetpt:met:qt:InvMass:Vertices:extraElectrons:extraMuons:nothing:somethingelse");
  //tree = new TTree("T", "First attempt at tree");
  tree = fs->make<TTree>("T","T");
  l14vector = new TLorentzVector();
  l24vector = new TLorentzVector();
  MC14vector = new TLorentzVector();
  MC24vector = new TLorentzVector();
  MCN14vector = new TLorentzVector();
  MCN24vector = new TLorentzVector();
  MCD14vector = new TLorentzVector();
  MCD24vector = new TLorentzVector();
  MCB14vector = new TLorentzVector();
  MCB24vector = new TLorentzVector();
  metP4 = new TLorentzVector();

     tree->Branch("lep1pt", &lep1pt, "lep1pt/D");
    tree->Branch("lep2pt", &lep2pt, "lep2pt/D");
   /*  tree->Branch("muoneta1", &meta1, "muoneta1/D");
   tree->Branch("muoneta2", &meta2, "muoneta2/D");
 */
   tree->Branch("numberjets", &num_jets, "numberjets/I");
   tree->Branch("leadjetpt", &jet_pt[0], "leadjetpt/D");
   tree->Branch("2ndjetpt", &jet_pt[1], "2ndjetpt/D");
   tree->Branch("3rdjetpt", &jet_pt[2], "3rdjetpt/D");
   tree->Branch("4thjetpt", &jet_pt[3], "4thjetpt/D");
   tree->Branch("5thjetpt", &jet_pt[4], "5thjetpt/D");
   tree->Branch("met", &mET, "met/D");
   tree->Branch("ZpT", &ZpT, "ZpT/D");
   tree->Branch("InvMass", &InvMass, "InvMass/D");
   tree->Branch("Vertices", &vtxCount, "Vertices/I");
     tree->Branch("extraElectrons", &extraElectrons, "extraElectrons/I");
      tree->Branch("extraMuons", &extraMuons, "extraMuons/I");
   tree->Branch("weight", &weight, "weight/D");
   tree->Branch("nPUVerticesTrue", &nPUVerticesTrue, "nPUVerticesTrue/I");
   // tree->Branch("normChi2", &normChi2, "normChi2/D");
   //  tree->Branch("d_xy",&d_xy, "d_xy/D");
   // tree->Branch("numVH",&numVH, "numVH/I");
   //  tree->Branch("numPH",&numPH, "numPH/I");
   //  tree->Branch("m_iso",&m_iso, "m_iso/D");
      tree->Branch("1stEmuonpt", &Emuon_pt[0], "1stEmuonpt/D");
      tree->Branch("2ndEmuonpt", &Emuon_pt[1], "2ndEmuonpt/D");
     tree->Branch("3rdEmuonpt", &Emuon_pt[2], "3rdEmuonpt/D");
     tree->Branch("1stEmuoneta", &Emuon_eta[0], "1stEmuoneta/D");
     tree->Branch("2ndEmuoneta", &Emuon_eta[1], "2ndEmuoneta/D");
     tree->Branch("3rdEmuoneta", &Emuon_eta[2], "3rdEmuoneta/D");
     tree->Branch("1stEelectronpt", &Eelectron_pt[0], "1stEelectronpt/D");
     tree->Branch("2ndEelectronpt", &Eelectron_pt[1], "2ndEelectronpt/D");
     tree->Branch("3rdEelectronpt", &Eelectron_pt[2], "3rdEelectronpt/D");
     tree->Branch("1stEelectroneta", &Eelectron_eta[0], "1stEelectroneta/D");
      tree->Branch("2ndEelectroneta", &Eelectron_eta[1], "2ndEelectroneta/D");
       tree->Branch("3rdEelectroneta", &Eelectron_eta[2], "3rdEelectroneta/D");
   //   tree->Branch("1stdReEm1", &dReEm1[0], "1stdReEm1/D");
   //   tree->Branch("1stdReEm2", &dReEm2[0], "1stdReEm2/D");
   //   tree->Branch("2nddReEm1", &dReEm1[1], "2nddReEm1/D");
   //   tree->Branch("2nddReEm2", &dReEm2[1], "2nddReEm2/D");
   //   tree->Branch("3rddReEm1", &dReEm1[2], "3rddReEm1/D");
   //   tree->Branch("3rddReEm2", &dReEm2[2], "3rddReEm2/D");
      tree->Branch("rMETparallel", &rMETparallel, "rMETparallel/D");
      tree->Branch("rMETperp", &rMETperp, "rMERperp/D");
      tree->Branch("redMETtotal", &redMETtotal, "redMETtotal/D");
      //    tree->Branch("muon1iso", &muon1iso, "muon1iso/D");
      // tree->Branch("muon2iso", &muon2iso, "muon2iso/D");
      tree->Branch("lepton1phi", &lepton1phi, "lepton1phi/D");
      tree->Branch("lepton2phi", &lepton2phi, "lepton2phi/D");
      tree->Branch("lepton1eta", &lepton1eta, "lepton1eta/D");
      tree->Branch("lepton2eta", &lepton2eta, "lepton2eta/D");      
      tree->Branch("DiLeptonType", &DiLeptonType, "DiLeptonType/I");
      tree->Branch("METSignificance", &METSignificance,"METSignificance/D");
      tree->Branch("acol", &acol, "acol/D");
      tree->Branch("lep1dxy", &lep1dxy, "lep1dxy/D");
      tree->Branch("lep2dxy", &lep2dxy, "lep2dxy/D");
      tree->Branch("lep1Q", &lep1Q, "lep1Q/D");
      tree->Branch("lep2Q", &lep2Q, "lep2Q/D");
      tree->Branch("vertexProb", &vertexProb, "vertexProb/F");
      tree->Branch("MCinvmass", &MCinvmass, "MCinvmass/D");
      tree->Branch("MC1pt", &MC1pt, "MC1pt/D");
      tree->Branch("MC1eta", &MC1eta, "MC1eta/D");
      tree->Branch("MC1phi", &MC1phi, "MC1phi/D");
      tree->Branch("MC2pt", &MC2pt, "MC2pt/D");
      tree->Branch("MC2eta", &MC2eta, "MC2eta/D");
      tree->Branch("MC2phi", &MC2phi, "MC2phi/D");
      tree->Branch("MC1Q", &MC1Q, "MC1Q/D");
      tree->Branch("MC2Q", &MC2Q, "MC2Q/D");
      tree->Branch("MCN1Q", &MCN1Q, "MCN1Q/D");
      tree->Branch("MCN2Q", &MCN2Q, "MCN2Q/D");
      tree->Branch("MCD1Q", &MCD1Q, "MCD1Q/D");
      tree->Branch("MCD2Q", &MCD2Q, "MCD2Q/D");
      tree->Branch("MCB1Q", &MCB1Q, "MCB1Q/D");
      tree->Branch("MCB2Q", &MCB2Q, "MCB2Q/D");
      tree->Branch("MCpTdiff", &MCpTdiff, "MCpTdiff/D");
      tree->Branch("DoubleMuonTrigger", &DoubleMuonTrigger, "DoubleMuonTrigger/O");
      tree->Branch("DoubleElectronTrigger", &DoubleElectronTrigger, "DoubleElectronTrigger/O");
      tree->Branch("MuonElectronTrigger", &MuonElectronTrigger, "MuonElectronTrigger/O");
      tree->Branch("SingleMuonTrigger40",&SingleMuonTrigger40,"SingleMuonTrigger40/O");
 tree->Branch("SingleMuonTriggeriso",&SingleMuonTriggeriso,"SingleMuonTriggeriso/O");
 tree->Branch("SingleMuonTriggeriso2",&SingleMuonTriggeriso2,"SingleMuonTriggeriso2/O");
 tree->Branch("SingleElectronTrigger",&SingleElectronTrigger,"SingleElectronTrigger/O");

 tree->Branch("aftertrigger", &aftertrigger, "aftertrigger/I");
 tree->Branch("liso1", &liso1, "liso1/D");
 tree->Branch("liso2", &liso2, "liso2/D");
 tree->Branch("l14vector", &l14vector, 6400, 0);
 tree->Branch("l24vector", &l24vector, 6400, 0);
 tree->Branch("lep1ptE", &lep1ptE, "lep1ptE/D");
 tree->Branch("lep2ptE", &lep2ptE, "lep2ptE/D");
 tree->Branch("elecReg1", &elecReg1, "elecReg1/D");
 tree->Branch("elecReg2", &elecReg2, "elecReg2/D");
 tree->Branch("eventid",&eventid, "eventid/I");
 tree->Branch("runid",&runid, "runid/I");
 tree->Branch("istautau",&istautau,"istautau/I");
 tree->Branch("MC14vector", &MC14vector, 6400, 0);
 tree->Branch("MC24vector", &MC24vector, 6400, 0);
 tree->Branch("MCN14vector", &MCN14vector, 6400, 0);
 tree->Branch("MCN24vector", &MCN24vector, 6400, 0);
 tree->Branch("MCD14vector", &MCD14vector, 6400, 0);
 tree->Branch("MCD24vector", &MCD24vector, 6400, 0);
 tree->Branch("MCB14vector", &MCB14vector, 6400, 0);
 tree->Branch("MCB24vector", &MCB24vector, 6400, 0);
 tree->Branch("metP4", &metP4, 6400,0);
 tree->Branch("mujetd1", &mujetd1, "mujetd1/D");
 tree->Branch("mujetd2", &mujetd2, "mujetd2/D");
 tree->Branch("isquark", &isquark , "isquark/D");
 tree->Branch("motherId1", &motherId1 , "motherId1/D");
 tree->Branch("motherId2", &motherId2 , "motherId2/D");
 tree->Branch("eff1", &eff1 , "eff1/I");
 
 // these are variables for the W study with isolation 8/2/14
 /*
tree->Branch("Npt", &Npt , "Npt/D");
tree->Branch("Neta", &Neta , "Neta/D");
tree->Branch("Nphi", &Nphi , "Nphi/D");
tree->Branch("Mpt", &Mpt , "Mpt/D");
tree->Branch("Meta", &Meta , "Meta/D");
tree->Branch("Mphi", &Mphi , "Mphi/D");
tree->Branch("Nisolation", &Nisolation , "Nisolation/D");
tree->Branch("Njetd", &Njetd , "Njetd/D");
 tree->Branch("Mq", &Mq , "Mq/D");
tree->Branch("WqT", &WqT , "WqT/D");
 */
}

//===================================================================//

void TupleMaker::fill(TH1D* h, int value)
{
  fill(h, (double)(value));
}

//===================================================================//

void TupleMaker::fill(TH1D* h, double value)
{
  //  cout<<weight<<endl;
  h->Fill(value, weight);
}

//===================================================================//

double TupleMaker::get_weight()
{
  return (DEFAULT_WEIGHT);     //unity by default
}

//===================================================================//

void TupleMaker::setup_counters()
{
  num_events = num_trilep = num_no_lep = 0;
    
  num_2jet_veto = num_bjet_veto = num_isoj_veto = num_lldz_veto = 0;
  num_HT1_veto = num_HT2_veto = num_mET1_veto = num_mET2_veto = 0;
  num_3jet_veto = num_Z_ll_veto = num_0flavor = num_unfound = 0;
  
  num_pass_ee = num_pass_em = num_pass_mm = 0;
  num_pass_ej = num_pass_mj = 0;
  num_good_ee = num_good_em = num_good_mm = 0;
  num_good_ej = num_good_mj = 0;
  
  num_3b_xx = num_ss_ee = num_ss_mm = num_ss_em = 0;
  num_2b_ee = num_2b_mm = num_2b_em = 0;
  num_2b_ej = num_2b_mj = 0;
  
  num_miss_b = num_unkn_b = num_true_b = 0;
  num_falseb = num_unkn_j = num_uds_g = 0;
}

//===================================================================//

void TupleMaker::setup_flags()   //last updated Sept. 21 2012
{
  cout 
  << endl << LINE << endl << " Using these collections:"
      << endl << "\t MET: " << MET_COLLECTION
    //<< endl << "\t Beam Spot: " << BEAMSPOT_COLLECTION
    //<< endl << "\t Vertex: " << VERTEX_COLLECTION
    //<< endl << "\t Electron: " << ELEC_COLLECTION
  << endl << "\t Muon: " << MUON_COLLECTION
    //  << endl << "\t Flavor: " << FLAVOUR_COLLECTION
  << endl << "\t Jet: " << JET_COLLECTION;
    // << endl << "\t btag: " << BTAG_COLLECTION
    // << endl << "\t JEC: " << CORR_COLLECTION << " w/"
    //<< endl << "\t > " << JET_CORR_L1 << "+" << JET_CORR_L2 
    //<< "+" << JET_CORR_L3 << "+" << JET_CORR_MC;
  
  cout 
  << endl << LINE << endl << " Using these (bool) flags:"
  << endl << "\t Running over data = " << get_true_false(IS_DATA)
  << endl << "\t Skip same-sign cut = " << get_true_false(SKIP_SS_CUT)
  << endl << "\t Skip Z mass cut (dileptons) = " << get_true_false(SKIP_Z_MASS_CUT)
  << endl << "\t Skip MET cut = " << get_true_false(SKIP_MET_CUT)
  << endl << "\t Skip HT cut = " << get_true_false(SKIP_HT_CUT);
  
  cout 
  << endl << LINE << endl << " Using these cuts:"
  << endl << "\t vaild tracker hits for muons = " << CUT_TRACK_VHITS
  << endl << "\t tight lepton pT = " << CUT_LEPTON_TRIG_PT
  << endl << "\t loose lepton pT = " << CUT_LEPTON_LOOSE_PT
  << endl << "\t electron eta = " << CUT_ELEC_ETA
  << endl << "\t muon eta = " << CUT_MUON_ETA
  << endl << "\t REI = " << CUT_ELEC_ISO
  << endl << "\t RMI = " << CUT_MUON_ISO
  << endl << "\t muon chi^2 = " << CUT_MUON_CHI2
  << endl << "\t electron |d0| = " << CUT_ELEC_DXY
  << endl << "\t muon |d0| = " << CUT_MUON_DXY
  << endl << "\t DeltaZ(lep1, lep2) = " << CUT_DELTAZ;
  
  
}

//===================================================================//

void TupleMaker::analyze(const Event& event_, const EventSetup& setup_)
{  


  eff1=eff2=eff3=eff4=0;
  aftertrigger=0;
  eleceff3=eleceff4=muoneff3=muoneff4=0;
  ZpT=-1;
  InvMass=-1;
  extraMuons=lep2pt=lep1pt=meta1=meta2=0;
  eiso1=eiso2=liso1=liso2=0;
  extraElectrons=0;
  nPUVertices = 0;
  nPUVerticesTrue = 0;
  weight=1.0;
  acol=0;
  normChi2=100.0;
  d_xy=1.0;
  numVH=0;
  numPH=0;
  lep1dxy=lep2dxy=.3;
  m_iso=1.0;
  lep1pt=lep2pt=0;
  MC1phi=MC1eta=MC1pt=MC2phi=MC2pt=MC2eta=MCinvmass=0;
  MC1Q=MC2Q=0;
  MCrapid=0.0;
  rMETperp=rMETparallel=0;
  rMET.SetPxPyPzE(0., 0., 0., 0.);
  lep1.SetPxPyPzE(0., 0., 0., 0.);
  lep2.SetPxPyPzE(0., 0., 0., 0.);
  metP4->SetPxPyPzE(0., 0., 0., 0.);
  sumJet.SetPxPyPzE(0., 0., 0., 0.);
  elec1P4.SetPxPyPzE(0.,0.,0.,0.);
  elec2P4.SetPxPyPzE(0.,0.,0.,0.);
  l14vector->SetPxPyPzE(0., 0., 0., 0.);
  l24vector->SetPxPyPzE(0., 0., 0., 0.);
  MC14vector->SetPxPyPzE(0., 0., 0., 0.);
  MC24vector->SetPxPyPzE(0., 0., 0., 0.);
  MCN14vector->SetPxPyPzE(0., 0., 0., 0.);
  MCN24vector->SetPxPyPzE(0., 0., 0., 0.);
  MCD14vector->SetPxPyPzE(0., 0., 0., 0.);
  MCD24vector->SetPxPyPzE(0., 0., 0., 0.);
  MCB14vector->SetPxPyPzE(0., 0., 0., 0.);
  MCB24vector->SetPxPyPzE(0., 0., 0., 0.);
  MCB1Q=MCB2Q=MCN1Q=MCN2Q=MCD1Q=MCD2Q=0;
  mujetd1=mujetd2=9.0;

  redMET.SetPxPyPzE(-999., -999., -999., -999.);
  redMETtotal=-1;
  muon1iso=-1.0;
  muon2iso=-1.0;
  elecReg1=elecReg2=0.0;
  num_e=0;
   num_mm=0, num_ee=0, num_me=0;
   lep1ptE=lep2ptE=0.0;
   DiLeptonType=-1;
   istautau=0;
  for(int i=0;i<3;i++)
    {dReEm1[i]=dReEm2[i]=7.0;
    }  
  //comment the next few lines for running over data
  
  //     load_pileup(event_);
  
  // cout<<weight<<endl;
  //  load_beam(event_);
 
     //this should now produce naked, bare and dressed muons for phistar 10/01/14
     // load_PhiStar(event_);


  //this is for the W isolation study 8/02/14
  //  load_Wacceptance(event_, setup_);
	eff1=0;  

	//	   load_acceptance(event_);
 



      bool passpt=false;
      passacc=false;
      //for SIM

 if(eff1>0&& nummuon>1)
      {
	MCrapid=get_rapid(*MC14vector+*MC24vector);
	//	cout<<"mass "<<MCrapid<<endl;
	//	double MCaccInvMass=get_mass(*MC14vector+*MC24vector);
	//	cout<<MC24vector->Pt()<<endl;
	//for SIM only
	weight=1.0;


	//	if(fabs(MCrapid)<2.0)
	  {
	    /*	    	
	if(MCinvmass>20 && MCinvmass<180)
	  {
	int lkl = (int)((MCinvmass-20)/2);
	weight=weights[lkl];
	  }
	    */
	    //	cout<<"weight "<<weight<<" , "<<MCinvmass<<endl;

      h_base->Fill(eff1, weight);

      if(fabs(MC1truth)==13 && fabs(MC2truth)==13)
	eff2=1;
      if(fabs(MC1truth)==11 && fabs(MC2truth)==11)
	eff2=3;
      
      if((fabs(MC1truth)==11 && fabs(MC2truth)==13) || (fabs(MC2truth)==11 && fabs(MC1truth)==13))
	eff2=5;


      if(eff1==1)
	{
      h_MCM1pt->Fill(MC14vector->Pt(), weight);
      h_MCM2pt->Fill(MC24vector->Pt(), weight);
      h_MCM1eta->Fill(MC14vector->Eta(), weight);
      h_MCM2eta->Fill(MC24vector->Eta(), weight);
      h_MCMIM->Fill(MCinvmass, weight);
      h_MCMrapid->Fill(MCrapid, weight);

	}
      if(eff1==3)
	{
      h_MCE1pt->Fill(MC14vector->Pt(), weight);
      h_MCE2pt->Fill(MC24vector->Pt(), weight);
      h_MCE1eta->Fill(MC14vector->Eta(), weight);
      h_MCE2eta->Fill(MC24vector->Eta(), weight);
      h_MCEIM->Fill(MCinvmass, weight);
      h_MCErapid->Fill(MCrapid, weight);

	}
      if(eff1==5)
	{
      h_MCS1pt->Fill(MC14vector->Pt(), weight);
      h_MCS2pt->Fill(MC24vector->Pt(), weight);
      h_MCS1eta->Fill(MC14vector->Eta(), weight);
      h_MCS2eta->Fill(MC24vector->Eta(), weight);
      h_MCSIM->Fill(MCinvmass, weight);
      h_MCSrapid->Fill(MCrapid, weight);

	}
      if(MC14vector->Pt()>MC24vector->Pt())
	{
	  if(MC14vector->Pt()>30 && fabs(MC14vector->Eta())<2.1)
	    passpt=true;
	}
      else
	{
	  if(MC24vector->Pt()>30 && fabs(MC24vector->Eta())<2.1)
	    {	    
	      passpt=true;
	    }
	}
      if(fabs(MC1truth)==11 || fabs(MC1truth)==13)
	{
       
	  if(!(fabs(MC14vector->Eta())<1.44 || fabs(MC14vector->Eta())>1.57))
	    passpt=false;
	}
      if(fabs(MC2truth)==11|| fabs(MC2truth)==13)
	{
	  //	  cout<<MC2truth<<" , "<<eff1<<endl;
	
	  if(!(fabs(MC24vector->Eta())<1.44 || fabs(MC24vector->Eta())>1.57))
	    passpt=false;
	}
      if(passpt && MC14vector->Pt()>20 && fabs(MC14vector->Eta())<2.4 && MC24vector->Pt()>20 && fabs(MC24vector->Eta())<2.4 && MCinvmass>60 && MCinvmass<120)
	{
	  passacc=true;
	  h_passacc->Fill(eff1, weight);

	}
	  }
      }








  load_prep();
  

  load_rMET(event_);

   num_m = load_muon(event_);

  
  load_triggers(event_);

  // cout<<"hello"<<endl;
  //===========================  
  //get MET significance
edm::Handle< edm::View<reco::PFMET> > pfMEThandle;
event_.getByLabel("pfMet", pfMEThandle);
METSignificance = (pfMEThandle->front() ).significance();
//double sigmaX2= (pfMEThandle->front() ).getSignificanceMatrix()(0,0);

  //===========================
 
  num_e= load_elec(event_, setup_, num_m);
  load_vtxs(event_);
  // cout<<"hello"<<endl;
  
  InputTag  vertexLabel(string("offlinePrimaryVertices"));
  Handle<reco::VertexCollection> thePrimaryVertexColl;
  event_.getByLabel(vertexLabel,thePrimaryVertexColl);
  // if (!(primaryVtcs->size()>0))
  Vertex dummy;
  const Vertex *pv = &dummy;
  if (thePrimaryVertexColl->size() != 0) {
    pv = &*thePrimaryVertexColl->begin();
  } else { // create a dummy PV
    Vertex::Error e;
    e(0, 0) = 0.0015 * 0.0015;
    e(1, 1) = 0.0015 * 0.0015;
    e(2, 2) = 15. * 15.;
    Vertex::Point p(0, 0, 0);
    dummy = Vertex(p, e, 0, 0, 0);
  }
  int dileptontype=0;
  // cout<<num_e<<" , "<<num_m<<endl;
  if(num_m>1 && num_e==0)
    {dileptontype=1;
    }
 if(num_m==0 && num_e>1)
    {dileptontype=3;
    }
 if(num_m==1 && num_e==1)
    {dileptontype=5;
    }
 if(num_m>1 && num_e==1)
   {  if(muon_1->pt()>elec1P4.Pt()&&muon_2->pt()>elec1P4.Pt())
    {dileptontype=1;
    }
  else
    dileptontype=5;
   }
 if(num_m==1 && num_e>1)
   {  if(muon_1->pt()<elec1P4.Pt()&&muon_1->pt()<elec2P4.Pt())
    {dileptontype=3;
    }
  else
    dileptontype=5;
   }
  if(num_m>2 && num_e>2)
   { if(muon_1->pt()>elec1P4.Pt()&&muon_2->pt()>elec1P4.Pt()&&muon_2->pt()>elec2P4.Pt())
    {dileptontype=1;
    }
      else if(muon_1->pt()<elec1P4.Pt()&&muon_1->pt()<elec2P4.Pt())
    {dileptontype=3;
    }
      else 
    {dileptontype=5;
    }
    }
  // cout<<dileptontype<<endl;
  //  load_corr(event_, setup_);
      InvMass=0;
  if(dileptontype==1)
    {
        num_mm++;
	liso1=get_m_iso(muon_1);
	liso2=get_m_iso(muon_2);
	//cout<<muon_1->track()->ptError()<<" , "<<muon_1->pt()<<" , "<<muon_1->track()->pt()<<endl;
	lep1ptE=muon_1->track()->ptError();
	lep2ptE=muon_2->track()->ptError();

     InvMass = get_mass( muon_1->p4() + muon_2->p4() );              
      ZpT =  get_pT( muon_1->p4() + muon_2->p4());
      DiLeptonType=1;
   lepton1phi=muon_1->phi();
    lepton2phi=muon_2->phi();
    lepton1eta=muon_1->eta();
    lepton2eta=muon_2->eta();
      lep1pt=muon_1->pt();
      lep2pt=muon_2->pt();
      lep1Q=muon_1->charge();
      lep2Q=muon_2->charge();
      lep1dxy=fabs(muon_1->muonBestTrack()->dxy(pv->position()));
      lep2dxy=fabs(muon_2->muonBestTrack()->dxy(pv->position()));
      /*    
   normChi21 = muon_1->globalTrack()->normalizedChi2();
   numPH21 = muon_1->innerTrack()->hitPattern().numberOfValidPixelHits(); 

   hitPattern1 = muon_1->globalTrack()->hitPattern().numberOfValidMuonHits();
   numberStations1 =  muon_1->numberOfMatchedStations();
   hits1 = muon_1->track()->hitPattern().trackerLayersWithMeasurement(); 
      */

      l14vector->SetPxPyPzE(muon_1->px(), muon_1->py(), muon_1->pz(), muon_1->p4().t());  
    l24vector->SetPxPyPzE(muon_2->px(), muon_2->py(), muon_2->pz(), muon_2->p4().t());  


      //     cout<<lep1dxy<<endl;

          acol = acolinearity(muon_1->p4(), muon_2->p4());

	  	  vertexProb=leptonVertex(muon_1->track(), muon_2->track(), setup_);


		  mujetd1=mujetd2=9.0;  
  event_.getByLabel(JET_COLLECTION, jets);
   
 
  PFJet cjet;
  //  cout<<endl;
  PFJetCollection::const_iterator ijet;
  //cout<<endl;
  for(ijet = jets->begin(); ijet != jets->end(); ++ijet)
  {

    double jpt=0;
    jpt=ijet->pt();
   cjet = correct_jet(ijet, event_, setup_);    


     //const JetCorrector* corrector = JetCorrector::getJetCorrector (JetCorrectionService,setup_);   //Get the jet corrector from the event setup
     //PFJet  correctedJet = *cjet;                                 //copy original jet
     //double jec = corrector->correction(*i_jet,jetRef,event_,setup_);
     //correctedJet.scaleEnergy(jec);        
     //cout<<correctedJet.pt()<<endl;

                // apply the correction
  
	 
	 // if(cut_jet_vtx(cjet))         continue;
   if((DeltaRX(cjet.eta(), muon_1->eta(), cjet.phi(), muon_1->phi())>.1) && (DeltaRX(cjet.eta(), muon_2->eta(), cjet.phi(), muon_2->phi())>.1))      
	   {
	     if((cjet.neutralEmEnergyFraction() < CUT_NEUTR_EM_EN_FRAC) && (cjet.chargedEmEnergyFraction() < CUT_CHARG_EM_EN_FRAC) && (cjet.neutralHadronEnergyFraction() < CUT_NEUTR_HAD_EN_FRAC) && (cjet.chargedHadronEnergyFraction() > CUT_CHARG_HAD_EN_FRAC))
	       {
		 if(DeltaRX(cjet.eta(), muon_1->eta(), cjet.phi(), muon_1->phi())<mujetd1)
		   {
		     mujetd1=DeltaRX(cjet.eta(), muon_1->eta(), cjet.phi(), muon_1->phi());
		     // cout<<mujetd1<<endl;
		     // cout<<muon_1->pt()<<" , "<<cjet.pt()<<endl;
		   }
		 if(DeltaRX(cjet.eta(), muon_2->eta(), cjet.phi(), muon_2->phi())<mujetd2)
		   {
		     mujetd2=DeltaRX(cjet.eta(), muon_2->eta(), cjet.phi(), muon_2->phi());
		     //cout<<mujetd2<<endl;
		   }
		 
	       }
	     
	   }

  }
  //  cout<<"hiho"<<mujetd1<<endl;
    
    
    }
  //  cout<<"hello"<<endl;
    if(dileptontype==3)
    {
      liso1=eiso1;
      liso2=eiso2;
    num_ee++;
  
    l14vector->SetPxPyPzE(elec1P4.Px(), elec1P4.Py(), elec1P4.Pz(), elec1P4.E()); 
    l24vector->SetPxPyPzE(elec2P4.Px(), elec2P4.Py(), elec2P4.Pz(), elec2P4.E());  

     InvMass = get_mass(elec1P4 + elec2P4);              
      ZpT =  get_pT(elec1P4 + elec2P4);
      DiLeptonType=3;
   lepton1phi=elec1phi;
    lepton2phi=elec2phi;
    lepton1eta=elec1eta;
    lepton2eta=elec2eta;
     lep1pt=elec1P4.Pt();
    lep2pt=elec2P4.Pt();
       acol = acolinearity(elec1P4, elec2P4);
       lep1dxy=elec1dxy;
       lep2dxy=elec2dxy;
       lep1Q=elec1Q;
       lep2Q=elec2Q;
       vertexProb=leptonVertex(elec1track, elec2track, setup_);
    }
  //  cout<<"hello"<<endl;
  if(dileptontype==5)
    {
      num_me++;

      l14vector->SetPxPyPzE(muon_1->px(), muon_1->py(), muon_1->pz(), muon_1->p4().t());  
      // cout<<l14vector->Px()<<endl;
    l24vector->SetPxPyPzE(elec1P4.Px(), elec1P4.Py(), elec1P4.Pz(), elec1P4.E());  

      liso1=get_m_iso(muon_1);
      liso2=eiso1;
     InvMass = get_mass(muon_1->p4() + elec1P4 );              
      ZpT =  get_pT( muon_1->p4() + elec1P4);
      DiLeptonType=5;
   lepton1phi=muon_1->phi();
    lepton2phi=elec1phi;
    lepton1eta=muon_1->eta();
    lepton2eta=elec1eta;
    lep1pt=muon_1->pt();
    lep2pt=elec1P4.Pt();
    acol = acolinearity(muon_1->p4(), elec1P4);
    lep1dxy=fabs(muon_1->muonBestTrack()->dxy(pv->position()));
    lep2dxy=elec1dxy;
    lep1Q=muon_1->charge();
    lep2Q=elec1Q;
    vertexProb=leptonVertex(muon_1->track(), elec1track, setup_);

        if((lep1Q==-lep2Q) && InvMass>60 && InvMass<120 && (METSignificance/mET)<4 &&(lep1pt+lep2pt)>50 && acol<2.5)
    	    {
	      //      cout<<"New Event :"<<InvMass<<endl<<endl<<endl;
	      //  cout<<"Muon: "<<lep1pt<<" , "<<lepton1eta<<" , "<<lepton1phi<<endl;
	      //    cout<<"Electron: "<<lep2pt<<" , "<<lepton2eta<<" , "<<lepton2phi<<endl;
    
	      //MC_tau(event_, lepton2phi, lepton2eta);
    	    }
	  //	  cout<<vertexProb<<endl;
    }
  
  // cout<<"hello"<<endl;
  //   if ((num_m+num_e)==2)
    // for phistar i changed this to just muons
  // if (num_m==2 || (MCinvmass>60 && MCinvmass<130))
  //     {//process_=PROCESS_di;
  //    }
 
     //  else
    //   process_=PROCESS_no;


      
 
      
      /*    
      normChi2 = muon_1->globalTrack()->normalizedChi2();
      d_xy = fabs(muon_1->globalTrack()->dxy(beam_xyz));
      
      numVH = muon_1->globalTrack()->hitPattern().numberOfValidTrackerHits();
      numPH = muon_1->globalTrack()->hitPattern().numberOfValidPixelHits();
      m_iso= get_m_iso(muon_1);
*/

     
 

  //get cone size between muon and nearest jet for phistar eff
    //below line should nominally exclude region from 71 to 111
  //  if((num_m+num_e)>1||eff1>0||Wpass>0)
 
 if((num_m+num_e)>1 || nummuon>1)
    {
      process_=PROCESS_di;
    }
  else
    {process_=PROCESS_no;
		veto_event();
      }
    //  cout<<InvMass<<endl;
    /*  
      mpt1=muon_1->pt();
    mpt2=muon_2->pt();
    meta1=muon_1->eta();
    meta2=muon_2->eta();
    lepton1phi=muon_1->phi();
    lepton2phi=muon_2->phi();
    lepton1eta=muon_1->eta();
    lepton2eta=muon_2->eta();
    muon1iso = get_m_iso(muon_1);
    muon2iso = get_m_iso(muon_2);

    lep1.SetPxPyPzE(muon_1->px(), muon_1->py(), muon_1->pz(), muon_1->p4().t());  
    lep2.SetPxPyPzE(muon_2->px(), muon_2->py(), muon_2->pz(), muon_2->p4().t());
    */

    load_jets(event_, setup_);
  
      lepton_search(event_, setup_);

 // load_flvr(event_);
  // load_btag(event_);
  
  //we should have >=2 jets & >=1 b-jets
  // w/ one (or two) high pT lepton
  //  load_dileptons();
  //load_leptonjet();

  // hello=muon_1->pt();

    //    cout<<mET<<" , "<<mETphi<<endl;

    //acceptance and efficiency numbers (eff1, efff2, eff3, eff4)




   if(eff1>0&& nummuon>1)
      {
      if((fabs(MC1truth)==11 && fabs(MC2truth)==13))
	{
	  //	  cout<<MC14vector->Pt()<<" , "<<MC24vector->Pt()<<endl;
	}
    if((fabs(MC1truth)==13 && fabs(MC2truth)==11))
	{
	  //	  cout<<MC24vector->Pt()<<" , "<<MC14vector->Pt()<<endl;
	}
      
 //lepton reco eff
      if(eff2>0&&passacc)
	{
	  // h_passacc->Fill(eff2, weight);

	  if(eff2==1)
	    {
	    h_passrecoeff->Fill(1, weight);
	    h_passrecoeff->Fill(1, weight);
	    if(muoneff3==1)
	      {
		//	cout<<"hello1"<<endl;
	    h_passrecoeff->Fill(2, weight);

	      }
	    if(muoneff3>1)
	      {
		//	cout<<"this should work"<<endl;
	    h_passrecoeff->Fill(2, weight);
	    h_passrecoeff->Fill(2, weight);

	      }
	    }
	  if(eff2==3)
	    {
	    h_passrecoeff->Fill(3, weight);
	    h_passrecoeff->Fill(3, weight);

	    if(eleceff3==1)
	      {
	    h_passrecoeff->Fill(4, weight);

	      }
	    if(eleceff3>1)
	      {
	    h_passrecoeff->Fill(4, weight);
	    h_passrecoeff->Fill(4, weight);

	      }
	    }
	  if(eff2==5)
	    {
	    h_passrecoeff->Fill(1, weight);
	    h_passrecoeff->Fill(3, weight);
	    if(eleceff3>0)
		      {
			h_passrecoeff->Fill(4, weight);
			
		      }
	    if(muoneff3>0)
		      {
			h_passrecoeff->Fill(2, weight);
			
		      }
	    }
	  
	  //	  h_passrecoeff->Fill(eleceff3, weight);
	   if(eleceff4==1)
	    {
	      h_passseleceff->Fill(4, weight);

	    }
	   if(eleceff4>1)
	     {
	      h_num_pass_em->Fill(4, weight);

	       h_passseleceff->Fill(4, weight);
	       h_passseleceff->Fill(4, weight);
	       if(SingleMuonTriggeriso || SingleElectronTrigger)
		 {aftertrigger=3;
		   h_num_pass_em->Fill(3, weight);
		 }
	     }
	   if(muoneff4==1)
	    {
	      h_passseleceff->Fill(2, weight);

	    }
	  if(muoneff4>1)
	    {
	      h_passseleceff->Fill(2, weight);
	      h_passseleceff->Fill(2, weight);
	      h_num_pass_em->Fill(2, weight);
	      if(SingleMuonTriggeriso || SingleElectronTrigger)
		{aftertrigger=1;
		h_num_pass_em->Fill(1, weight);
		}
	    }
	  if(muoneff4==1 && eleceff4==1)
	    {
	      h_num_pass_em->Fill(6, weight);

	      if(SingleMuonTriggeriso || SingleElectronTrigger)
		{h_num_pass_em->Fill(5, weight);
		  aftertrigger=5;
		}
	    }
	}
      }
    //lepton selection eff
    //    cout<<event_.id().event()<<" , "<<event_.id().run()<<endl;
    eventid=event_.id().event();
    runid=event_.id().run();

    //this is for acceptance study. We only trigger if there are 2 true muons which pass certain cuts
    // cout<<"eff1="<<eff1<<endl;
   
	load_NickAnalysis();
      
   
      
}

//===================================================================//
void TupleMaker::load_NickAnalysis()
{
     if(process_ == PROCESS_no)      return;
 
   if(num_jets>10)
      {num_jets=10;
      }
   if(mET>300)
     {mET=300;
     }
   //if(extraElectrons<=0)
   //  { extraElectrons=-1;
   //  }
   // if(extraMuons<=0)
   //  { extraMuons=-1;
   //  }
   int numberleptons = extraElectrons + extraMuons;
    //cout<<numberjets<<endl;
   //   ntuple->Fill(mpt1, mpt2, meta1, meta2, num_jets, jet_pt[0], jet_pt[1], jet_pt[2], jet_pt[3], jet_pt[4], mET, ZpT, InvMass, vtxCount, extraElectrons, extraElectrons);
      
      
   //    redMET = GetReducedMET(sumJet,l14vector,l24vector,metP4,2);
    redMETtotal=redMET.E();
    rMETperp=redMET.Py();
    rMETparallel=redMET.Px();
    //    cout<<"a  "<<DiLeptonType<<endl;
    // cout<<muon1iso<<" "<<muon2iso<<endl;
   // cout<<weight<<endl;
    tree->Fill();
    l14vector->Clear();
    l24vector->Clear();
    MC14vector->Clear();
    MC24vector->Clear();


}

//===================================================================//
void TupleMaker::lepton_search(const Event& event_, const EventSetup& setup_ )
{
  //currently there is incompatability between this module since I added running over electrons and muons. I think it is right below here where we specify muon_2 when there may not be one.
  double phil, phi1, phi2, etal, eta1, eta2;
  double Aeff=0.0;
 if(process_ == PROCESS_no)      return;
 
 //print_debug("load_muon", N_DEBUG);  
 phi1 = lepton1phi;
 eta1 = lepton1eta;
 phi2= lepton2phi; 
 eta2=lepton2eta;
 extraElectrons=extraMuons=0;
 
 
  event_.getByLabel(MUON_COLLECTION, muons);

   Handle<reco::VertexCollection> primaryVtcs;
  event_.getByLabel(VERTEX_COLLECTION, primaryVtcs);
   for(int l=0;l<3;l++)
    {Emuon_pt[l]=Emuon_eta[l]=Eelectron_pt[l]=Eelectron_eta[l]=0;
    } 
  if ((primaryVtcs->size()>0))
    {
  vector<Muon>::const_iterator imuon;
   
  for(imuon = muons->begin(); imuon != muons->end(); ++imuon) 
    {
      
      phil = imuon->phi();
      etal = imuon->eta();
      
             
      
      if(imuon->isGlobalMuon())
	{
	  if(imuon->isTrackerMuon() && imuon->pt()>8 && phi1<2.4 && get_m_iso(imuon) < .2)    
	    {
	      if(DeltaRX(etal, eta1, phil, phi1) > .03 && DeltaRX(etal, eta2, phil, phi2) > .03)
		{
		      if( muon::isSoftMuon(*imuon, *primaryVtcs->begin())) 
			{
		   Emuon_pt[extraMuons]=imuon->pt();
		   Emuon_eta[extraMuons]=eta1;
		  ++extraMuons;
		  //cout<<extraMuons<<endl;
			}
		}
	      
	    }
	}
    }

    }
 



  Handle<edm::ValueMap<float>> mvaTrigV0_handle;
    event_.getByLabel("mvaTrigV0", mvaTrigV0_handle);
    const edm::ValueMap<float> ele_mvaTrigV0 = (*mvaTrigV0_handle.product());

  
    //  evaluate_mvas(event_, setup_);

  InputTag gsfEleLabel(string("gsfElectrons"));
  Handle<GsfElectronCollection> theEGammaCollection;
  event_.getByLabel(gsfEleLabel,theEGammaCollection);
  const GsfElectronCollection theEGamma = *(theEGammaCollection.product());

  edm::Handle<reco::PhotonCollection> photonH;
  event_.getByLabel(inputTagPhotons_,photonH);
   

  InputTag genLable(string("genParticles"));
  Handle<GenParticleCollection> genParticles;
  event_.getByLabel(genLable,genParticles);
  //InputTag  mcTruthLabel(string("generator"));
  //edm::Handle<edm::HepMCProduct> pMCTruth;
  //iEvent.getByLabel(mcTruthLabel,pMCTruth);
  //const HepMC::GenEvent* genEvent = pMCTruth->GetEvent();

  InputTag  vertexLabel(string("offlinePrimaryVertices"));
  Handle<reco::VertexCollection> thePrimaryVertexColl;
  event_.getByLabel(vertexLabel,thePrimaryVertexColl);

  _Rho=0;
  edm::Handle<double> rhoPtr;
 
  const edm::InputTag eventrho("kt6PFJets", "rho");
  event_.getByLabel(eventrho,rhoPtr);
  _Rho=*rhoPtr;
 

  
  edm::Handle<reco::ConversionCollection> hConversions;
  event_.getByLabel("allConversions", hConversions);
  
  Vertex dummy;
  const Vertex *pv = &dummy;
  if (thePrimaryVertexColl->size() != 0) {
    pv = &*thePrimaryVertexColl->begin();
  } else { // create a dummy PV
    Vertex::Error e;
    e(0, 0) = 0.0015 * 0.0015;
    e(1, 1) = 0.0015 * 0.0015;
    e(2, 2) = 15. * 15.;
    Vertex::Point p(0, 0, 0);
    dummy = Vertex(p, e, 0, 0, 0);
  }
  
  InputTag  reducedEBRecHitCollection(string("reducedEcalRecHitsEB"));
  InputTag  reducedEERecHitCollection(string("reducedEcalRecHitsEE"));

  EcalClusterLazyTools lazyTools(event_, setup_, reducedEBRecHitCollection, reducedEERecHitCollection);
  
  edm::ESHandle<TransientTrackBuilder> builder;
  setup_.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);
  TransientTrackBuilder thebuilder = *(builder.product());
  

  bool debug = true;
  bool debugMVAclass = false;
  bool debugMyVar = false;


   event_.getByLabel("offlineBeamSpot", bsHandle);
   const reco::BeamSpot &beamspot = *bsHandle.product();

   typedef std::vector< edm::Handle< edm::ValueMap<reco::IsoDeposit> > > IsoDepositMaps;
   typedef std::vector< edm::Handle< edm::ValueMap<double> > > IsoDepositVals;
  unsigned nTypes=3;
  IsoDepositMaps electronIsoDep(nTypes);

  for (size_t j = 0; j<inputTagIsoDepElectrons_.size(); ++j) {
    event_.getByLabel(inputTagIsoDepElectrons_[j], electronIsoDep[j]);
  }


//  IsoDepositMaps photonIsoDep(nTypes);
//    for (size_t j = 0; j<inputTagIsoDepPhotons_.size(); ++j) {
//      event_.getByLabel(inputTagIsoDepPhotons_[j], photonIsoDep[j]);
//    }
  IsoDepositVals electronIsoValPFId(nTypes);
//  IsoDepositVals photonIsoValPFId(nTypes);
     const IsoDepositVals * electronIsoVals =  &electronIsoValPFId  ;

  for (size_t j = 0; j<inputTagIsoValElectronsPFId_.size(); ++j) {
    event_.getByLabel(inputTagIsoValElectronsPFId_[j], electronIsoValPFId[j]);
  }

//  for (size_t j = 0; j<inputTagIsoValPhotonsPFId_.size(); ++j) {
//    event_.getByLabel(inputTagIsoValPhotonsPFId_[j], photonIsoValPFId[j]);
//  }

 
      for (uint j=0; j<theEGamma.size();j++)
  {

    //cout<<ele_mvaTrigV0.get(j)<<endl;
	  bool elePresel = trainTrigPresel(theEGamma[j]);

	  double mvaTrigMthd1 = -11.;
	  double mvaTrigMthd2 = -11.;

	  double mvaTrigNonIp = -11.;

	  double mvaNonTrigMthd1 = -11;
	  double mvaNonTrigMthd2 = -11;

	  //	  	  mvaNonTrigMthd1 = myMVANonTrigV0->mvaValue((theEGamma[j]),*pv,thebuilder,lazyTools,debugMVAclass);

		  //	cout<<mvaNonTrigMthd1<<endl;
    //
        // get particle flow isolation
	

      // cout<<ie.id()<<endl;
           phil = theEGamma[j].phi();
	   etal = theEGamma[j].eta();

  
	  

	         reco::GsfElectronRef myElectronRef(theEGammaCollection,j);



  
		 if(theEGamma[j].pt()>10  && fabs(eta1)<2.4)
{ // cout<<extraElectrons<<endl;
  
   if(DeltaRX(etal, eta1, phil, phi1) > .03 && DeltaRX(etal, eta2, phil, phi2) > .03)
     {
  if(elePresel)
    {
      
      //cout<<charged<<" ,"<<photon<<" ,"<<neutral<<" the iso"<<iso<<endl;
      // mvaTrigMthd1 = myMVATrigV0->mvaValue((theEGamma[j]),*pv,thebuilder,lazyTools,debugMVAclass);
      
      if((theEGamma[j].pt()<20&&((ele_mvaTrigV0.get(j)>0.00 && fabs(eta1)<.8) || (ele_mvaTrigV0.get(j)>.10 && (fabs(eta1)>.8 && fabs(eta1)<1.479)) || (ele_mvaTrigV0.get(j)>.62 && fabs(eta1)>1.479)))||((ele_mvaTrigV0.get(j)>.94 && fabs(eta1)<.8) || (ele_mvaTrigV0.get(j)>.85 && (fabs(eta1)>.8 && fabs(eta1)<1.479)) || (ele_mvaTrigV0.get(j)>.92 && fabs(eta1)>1.479)))
	{
	  double charged =  (*(*electronIsoVals)[0])[myElectronRef];
	  double photon = (*(*electronIsoVals)[1])[myElectronRef];
	  double neutral = (*(*electronIsoVals)[2])[myElectronRef];
	  if(abs(etal)>2.4)
	    Aeff=0.261;
	  if(abs(etal)>2.3 && abs(etal)<2.4)
	    Aeff=0.194;
	  if(abs(etal)>2.2 && abs(etal)<2.3)
	    Aeff=0.183;
	  if(abs(etal)>2.0 && abs(etal)<2.2)
	    Aeff=0.143;
	  if(abs(etal)>1.479 && abs(etal)<2.0)
	    Aeff=0.115;
	  if(abs(etal)>1.0 && abs(etal)<1.479)
	    Aeff=0.209;
	  if(abs(etal)<1.0)
	    Aeff=.208;
	  
	  
	  
	  
	  
	  
	  
	  
	  double iso =(charged+max(photon+neutral-_Rho*Aeff,0.0))/theEGamma[j].pt();
	  
	  if(iso<.15)
	    {
	      //calculate sip
	      float ip3d    = -999.0;
	      float ip3derr = 1.0;
	      float ip3dSig = 0.0;
	      
	       if (theEGamma[j].gsfTrack().isNonnull()) {
		 const double gsfsign = ( (-theEGamma[j].gsfTrack()->dxy(pv->position())) >=0 ) ? 1. : -1.;
		 
		 
		 const reco::TransientTrack &tt = thebuilder.build(theEGamma[j].gsfTrack());
		 
		 const std::pair<bool,Measurement1D> &ip3dpv = IPTools::absoluteImpactParameter3D(tt,*pv);
		 if (ip3dpv.first) {
		   ip3d = gsfsign*ip3dpv.second.value();
		   ip3derr = ip3dpv.second.error();
		   ip3dSig = ip3d/ip3derr;
		 }
	       }
	       int	misshits = theEGamma[j].gsfTrack()->trackerExpectedHitsInner().numberOfHits();
	       //cout<<misshits<<endl;
	       if(misshits<=1)
		 {
		   // cout<<ip3d<<" ,"<<ip3derr<<" , "<<ip3dSig<<endl;
		   if(ip3dSig<4)
		     {
		       //expected inner hits should be <=1
		       
		       
		       
		       
		       
		       
		       dReEm1[extraElectrons] = DeltaRX(etal, eta1, phil, phi1);
		       dReEm2[extraElectrons] = DeltaRX(etal, eta2, phil, phi2);
		       //		   if(dReEm1<.01 || dReEm2<.01)
		       //  {
		       Eelectron_pt[extraElectrons]=theEGamma[j].pt();
		       Eelectron_eta[extraElectrons]=theEGamma[j].eta();
		       //	      MC_info(etal, phil,Eelectron_pt[extraElectrons], event_);
		       ++extraElectrons;
		       //  }
		     }
		 }  
	    }  
	}
      
      
    }  
     }
		 
}
  }
}
//=========================================================================//
TLorentzVector TupleMaker::GetReducedMET(TLorentzVector sumJet, TLorentzVector lep1, TLorentzVector lep2, TLorentzVector metP4, int version) 
{
    TLorentzVector Q  = lep1 + lep2;
    float bisectorPhi = min(lep1.Phi(), lep2.Phi()) + lep1.DeltaPhi(lep2)/2;

    TVector2 projDilepton(Q.Px(), Q.Py());
    TVector2 projSumJet(sumJet.Px(), sumJet.Py());
    TVector2 projMET(metP4.Px(), metP4.Py());

    //TVector2 delta;
    //TLorentzVector Thrust = lep1 - lep2; 
    //if (fabs(lep1.DeltaPhiX(lep2)) > TMath::Pi()/2) {
    //	delta.Set(0., projDilepton.Rotate(-Thrust.Phi()).Py() + projDilepton.Py());
    //} else {
    //	delta.Set(0, projDilepton.Rotate(-Q.Phi()).Py() + projDilepton.Py());
    //}
    // cout<<Q.Phi()<<" , "<<Q.Px()<<" , "<<Q.Py()<<endl;
    if (version == 1) {
        projDilepton = projDilepton.Rotate(-bisectorPhi);
        projSumJet   = projSumJet.Rotate(-bisectorPhi);
        projMET      = projMET.Rotate(-bisectorPhi);
    } else if (version == 2) {
        projDilepton = projDilepton.Rotate(-Q.Phi());
        projSumJet   = projSumJet.Rotate(-Q.Phi());
        projMET      = projMET.Rotate(-Q.Phi());
    }
    // cout<<" sumJet "<<projSumJet.Px()<<" , "<<projSumJet.Py()<<endl;
    TVector2 unclustered = projSumJet;                     // projDilepton - 1.*(projMET + projDilepton) + delta;
    TVector2 clustered   = -1.*(projDilepton + 1.*projMET);   // + delta;

    TVector2 reducedMET = TVector2(projDilepton.Px()+(fabs(unclustered.Px()) < fabs(clustered.Px()) ? unclustered.Px() : clustered.Px()),
				   projDilepton.Py()+(fabs(unclustered.Py()) < fabs(clustered.Py()) ? unclustered.Py() : clustered.Py()));

    //    TVector2 reducedMET = TVector2(projDilepton.Px()-reducedMET1.Px(),projDilepton.Py()-reducedMET1.Py());
    //cout<<"hello"<<projSumJet.Py()<<" , "<<projMET.Py()<<" , "<<reducedMET.Py()<<endl;
    return TLorentzVector(reducedMET.Px(), reducedMET.Py(), 0, reducedMET.Mod());
}

//=========================================================================//
/*
void TupleMaker:: MC_info(double& etae, double& phie,double& pte,const Event& event_)
{

  //event_.getByLabel( "generator", EvtHandle );
  //const HepMC::GenEvent* Evt = EvtHandle->GetEvent() ;

  cout<<"ID number , pt, phi, eta, |electron pt, phi, eta"<<endl;


   event_.getByLabel("genParticles", genParticles);
   for(size_t i = 0; i < genParticles->size(); i++ ) 
     {
     const GenParticle & p = (*genParticles)[i];
     if(fabs(p.pdgId())==15)
     int id = p.pdgId();
     int st = p.status();  
     const Candidate * mom = p.mother();
     double pt = p.pt(), eta = p.eta(), phi = p.phi(), mass = p.mass();
     double vx = p.vx(), vy = p.vy(), vz = p.vz();
     int charge = p.charge();
     int n = p.numberOfDaughters();
     cout<<"New tau"<<endl;
     for(size_t j = 0; j < n; ++ j) {
       const Candidate * d = p.daughter( j );
       int dauId = d->pdgId();
       cout<<dauId<<" , "<<endl;
       //     double Del = DeltaRX(etae, eta, phie, phi);
   

       //if(Del<.05)
       //  cout<<id<<" , "<<st<<" , "<<pt<<" , "<<eta<<" , "<<phi<<" | "<<pte<<" , "<<Del<<endl;
 
     }


}
*/
 //=========================================================================//
   void TupleMaker::load_triggers(const Event& event_)
{
 

  DoubleMuonTrigger=false;
  DoubleElectronTrigger=false;
  MuonElectronTrigger=false;
  SingleMuonTrigger40=false;
  SingleMuonTriggeriso=false;
  SingleMuonTriggeriso2=false;
  SingleElectronTrigger=false;
  triggerResultsTag_ = InputTag(hlTriggerResults_,"","HLT");
 Handle<TriggerResults> hltResults;
  event_.getByLabel(triggerResultsTag_,hltResults);

  triggerEventTag_ = InputTag("hltTriggerSummaryAOD","","HLT");
  Handle<trigger::TriggerEvent> hltEvent;
  event_.getByLabel(triggerEventTag_,hltEvent);

  const TriggerNames & triggerNames = event_.triggerNames(*hltResults);
  hlNames = triggerNames.triggerNames();
  
  for (int i=0; i < (int)hlNames.size(); ++i) {
    if (!triggerDecision(hltResults, i)) continue;
   
      if (hlNames[i].compare(0, TriggerDoubleMu.length(),TriggerDoubleMu) == 0) {
	    DoubleMuonTrigger=true;
	    //	    cout<<hlNames[i]<<" , "<<endl<<endl;
      
	  }
      if (hlNames[i].compare(0, TriggerSingleMuon.length(),TriggerSingleMuon) == 0) {
	    SingleMuonTrigger40=true;
	    //	    cout<<hlNames[i]<<" , "<<endl<<endl;
      
	  }
      if (hlNames[i].compare(0, TriggerMuoniso2.length(),TriggerMuoniso2) == 0) {
	    SingleMuonTriggeriso2=true;
	    //	    cout<<hlNames[i]<<" , "<<endl<<endl;
      
	  }
          if (hlNames[i].compare(0, TriggerMuoniso.length(),TriggerMuoniso) == 0) {
	    SingleMuonTriggeriso=true;
	    //	    cout<<hlNames[i]<<" , "<<endl<<endl;
      
	  }
     if (hlNames[i].compare(0, TriggerDoubleEl.length(),TriggerDoubleEl) == 0) {
	    DoubleElectronTrigger=true;
	    //	    cout<<hlNames[i]<<" , "<<endl<<endl;
      
	  }
  if (hlNames[i].compare(0, TriggerSingleElectron.length(),TriggerSingleElectron) == 0) {
	    SingleElectronTrigger=true;
	    //	    cout<<hlNames[i]<<" , "<<endl<<endl;
	    
  }
  
  if (hlNames[i].compare(0, TriggerMuonElectron.length(),TriggerMuonElectron) == 0 || (hlNames[i].compare(0, TriggerMuonElectron2.length(),TriggerMuonElectron2) == 0)) {
	    MuonElectronTrigger=true;
	    //	    cout<<hlNames[i]<<" , "<<endl<<endl;
	    
  }
  }
  // cout<<DoubleMuonTrigger<<endl;
}
   //=========================================================================//
   /*
   void TupleMaker:: MC_tau(const Event& event_, double elecPhi, double elecEta)
   {
   
  //event_.getByLabel( "generator", EvtHandle );
  //const HepMC::GenEvent* Evt = EvtHandle->GetEvent() ;
  
  
  InputTag genLable(string("genParticles"));
  Handle<GenParticleCollection> genParticles;
  event_.getByLabel("genParticles", genParticles);
  for(size_t i = 0; i < genParticles->size(); i++ ) 
  {
  const GenParticle & p = (*genParticles)[i];
  
  if(fabs(p.pdgId())==23  &&p.status()==3)
  {
  cout<<"Here is info about the Z decay"<<endl;
  int n = p.numberOfDaughters();
  for(int j = 0; j < n; j++)
  {
  const Candidate * d = p.daughter( j );
  int dauId = d->pdgId();
  cout<<dauId<<" , "<<d->pt()<<" , "<<d->eta()<<" , "<<d->phi()<<endl;
  }
  }
  if(p.pt()>1 && p.mass()<1 && DeltaRX(elecEta, p.eta(), elecPhi, p.phi())<.05)
       {     cout<<"Particles close to the electron"<<endl;
       
       cout<<p.pdgId()<<" , "<<p.pt()<<" , "<<p.eta()<<" , "<<p.phi()<<" , "<<DeltaRX(elecEta, p.eta(), elecPhi, p.phi())<<endl;
       }
       if(fabs(p.pdgId())==13)
       {
       int id = p.pdgId();
       int st = p.status();  
       const Candidate * mom = p.mother();
       double pt = p.pt(), eta = p.eta(), phi = p.phi(), mass = p.mass();
       double vx = p.vx(), vy = p.vy(), vz = p.vz();
       int charge = p.charge();
       int n = p.numberOfDaughters();
       cout<<"New Muon"<<endl<<endl;
       for(int j = 0; j < n; j++)
       {
       const Candidate * d = p.daughter( j );
       int dauId = d->pdgId();
       cout<<dauId<<" , "<<d->pt()<<" , "<<d->eta()<<" , "<<d->phi()<<endl;
       //     double Del = DeltaRX(etae, eta, phie, phi);
       
       if(fabs(dauId) == 12 ||fabs(dauId) == 14 ||fabs(dauId) == 16)
       {
       //cout<<"Neutrino pt "<<d->pt();
       
       //if(Del<.05)
       //  cout<<id<<" , "<<st<<" , "<<pt<<" , "<<eta<<" , "<<phi<<" | "<<pte<<" , "<<Del<<endl;
       }
       
       }
       }
       if(fabs(p.pdgId())==15)
       {
       int id = p.pdgId();
       int st = p.status();  
       const Candidate * mom = p.mother();
       double pt = p.pt(), eta = p.eta(), phi = p.phi(), mass = p.mass();
       double vx = p.vx(), vy = p.vy(), vz = p.vz();
       int charge = p.charge();
       int n = p.numberOfDaughters();
       cout<<"New Tau"<<endl<<endl;
       for(int j = 0; j < n; j++)
       {
       const Candidate * d = p.daughter( j );
       int dauId = d->pdgId();
       cout<<dauId<<" , "<<d->pt()<<" , "<<d->eta()<<" , "<<d->phi()<<endl;
       //     double Del = DeltaRX(etae, eta, phie, phi);
       
       if(fabs(dauId) == 12 ||fabs(dauId) == 14 ||fabs(dauId) == 16)
       {
       //cout<<"Neutrino pt "<<d->pt();
       
       //if(Del<.05)
       //  cout<<id<<" , "<<st<<" , "<<pt<<" , "<<eta<<" , "<<phi<<" | "<<pte<<" , "<<Del<<endl;
       }
       
       }
       }
       }
       }
   */
   //=========================================================//
/*
   void TupleMaker::load_Wacceptance(const Event& event_,  const EventSetup& setup_)
   {

 InputTag  vertexLabel(string("offlinePrimaryVertices"));
  Handle<reco::VertexCollection> thePrimaryVertexColl;
 
 event_.getByLabel(vertexLabel,thePrimaryVertexColl);
  // if (!(primaryVtcs->size()>0))
  Vertex dummy;
  const Vertex *pv = &dummy;
  if (thePrimaryVertexColl->size() != 0) {
    pv = &*thePrimaryVertexColl->begin();
  } else { // create a dummy PV
    Vertex::Error e;
    e(0, 0) = 0.0015 * 0.0015;
    e(1, 1) = 0.0015 * 0.0015;
    e(2, 2) = 15. * 15.;
    Vertex::Point p(0, 0, 0);
    dummy = Vertex(p, e, 0, 0, 0);
  }
  int index, iVertex;
  Handle<PFCandidateCollection> pfCands;
  event_.getByLabel("particleFlow",pfCands);
  const PFCandidateCollection thePfColl = *(pfCands.product());

  Handle<GenParticleCollection> genParticles;
  event_.getByLabel("genParticles", genParticles);
  const Candidate * mom[50];

  int q;
  int Id=0;
  double MC1=-1;
  double MC1pthold=-1;
  double MC2=-1;
  double MC2pthold=-1;
  int nummuon=0;
  int numneutrino=0;
  double neta=0.0;
  double nphi=0.0;
  double pfiso_ch=0.0;
  double pfiso_em=0.0;
  double pfiso_nh=0.0;
  Meta=Mphi=Mpt=0.0;
  Nisolation=9.0;
  Wpass=0;
  WqT=0;
  Mq=0;

  for(size_t i = 0; i < genParticles->size(); i++ )
     {

     const GenParticle & p = (*genParticles)[i];
  
 if(fabs(p.pdgId())==13 && p.status()==1)
       {
	 //	 cout<<"hello"<<endl;
     mom[1] = p.mother();
      q=1;
      Id= fabs(mom[q]->pdgId());
     while(Id==13)
       {
	 q++;
     mom[q] = mom[q-1]->mother();
     Id=fabs(mom[q]->pdgId());
     if(fabs(Id)==24)
       {WqT=mom[q]->mass();
       }
       }
     if(p.pt()>2 && fabs(p.eta())<7 && fabs(Id)==24)
       {
	 nummuon++;

	 if(nummuon==1)
	   {
	     MC1=i;
	     MC1pthold=p.pt();
	     Meta=p.eta();
	     Mphi=p.phi();
	     Mpt=p.pt();
	     Mq=p.charge();
	   }
       }
       }
 if(fabs(p.pdgId())==14 && p.status()==1)
   {
     mom[1] = p.mother();
     q=1;
     Id= fabs(mom[q]->pdgId());
     while(Id==14)
       {
	 q++;
	 mom[q] = mom[q-1]->mother();
	 Id=fabs(mom[q]->pdgId());
	 
       }
     if(p.pt()>2 && fabs(p.eta())<2.4&& fabs(Id)==24)
       {
	 numneutrino++;
	 
	 if(numneutrino==1)
	   {
	     MC2=i;
	     MC2pthold=p.pt();
	     neta=p.eta();
	     nphi=p.phi();
	   }
       }
   }
     }
  //  cout<<nummuon<<" , "<<MC1pthold<<" , "<<numneutrino<<" , "<<MC2pthold<<endl;
  for(PFCandidateCollection::const_iterator it = pfCands->begin(); it<pfCands->end(); it++)
    {
      double teta=it->eta();
      double tphi=it->phi();
      //charged particles
      if(it->particleId() == reco::PFCandidate::h)
	{
	  double dzmin = 10000;
	  double ztrack = it->vertex().z();
	  bool foundVertex = false;
	  index = 0;
	  for(auto iv=thePrimaryVertexColl->begin(); iv!=thePrimaryVertexColl->end(); ++iv, ++index) {
	    
	    double dz = fabs(ztrack - iv->z());
	    if(dz<dzmin) {
	      dzmin = dz; 
	      iVertex = index;
	      foundVertex = true;
	    }
	   
	  }
	  //	  cout<<iVertex<<" , "<<it->et()<<endl;
	  //cout<<it->et()<<" , "<<it->particleId()<<endl;
	  if( (DeltaRX(neta, teta, nphi, tphi))<0.4 && iVertex==0)
	    {
	      pfiso_ch+= it->pt();
	    }
	}
      //neutral particles
      if(it->particleId() == reco::PFCandidate::h0)
	{
	  if((  DeltaRX(neta, teta, nphi, tphi))<0.4)
	    {
	      pfiso_nh+= it->pt();

	    }	  
	  //cout<<it->et()<<" , "<<it->particleId()<<endl;
	}
      //em particles
      if(it->particleId() == reco::PFCandidate::gamma)
	{
	  if((DeltaRX(neta, teta, nphi, tphi))<0.4)
	    {
	      pfiso_em+= it->pt();
	    }
	  // cout<<it->et()<<" , "<<it->particleId()<<endl;
	}
    }
  // cout<<endl<<endl;
  // cout<<pfiso_em<<" , "<<pfiso_ch<<" , "<<pfiso_nh<<endl;
  Nisolation=(pfiso_em+pfiso_ch+pfiso_nh)/MC2pthold;
  Npt=MC2pthold;
  Neta=neta;
  Nphi=nphi;
  

  event_.getByLabel(JET_COLLECTION, jets);
   
  Njetd=9.0;
  PFJet cjet;
  //  cout<<endl;
  PFJetCollection::const_iterator ijet;
  //cout<<endl;
  for(ijet = jets->begin(); ijet != jets->end(); ++ijet)
  {

    double jpt=0;
    jpt=ijet->pt();
   cjet = correct_jet(ijet, event_, setup_);    

   
   if(DeltaRX(cjet.eta(), Neta, cjet.phi(), Nphi)<Njetd)
     {
     Njetd=DeltaRX(cjet.eta(), Neta, cjet.phi(), Nphi);
     }


  }
  if(numneutrino>0)
    {Wpass=1;
    }

   }

*/
//==================================================================================
/*
void TupleMaker::load_PhiStar(const Event& event_)
{
  InputTag  vertexLabel(string("offlinePrimaryVertices"));
  Handle<reco::VertexCollection> thePrimaryVertexColl;
 
  event_.getByLabel(vertexLabel,thePrimaryVertexColl);
  // if (!(primaryVtcs->size()>0))
  Vertex dummy;
  const Vertex *pv = &dummy;
  if (thePrimaryVertexColl->size() != 0) {
    pv = &*thePrimaryVertexColl->begin();
  } else { // create a dummy PV
    Vertex::Error e;
    e(0, 0) = 0.0015 * 0.0015;
    e(1, 1) = 0.0015 * 0.0015;
    e(2, 2) = 15. * 15.;
    Vertex::Point p(0, 0, 0);
    dummy = Vertex(p, e, 0, 0, 0);
  }
  double MC1;
  double MC1pthold;
  double MC2;
  double MC2pthold;
  int q=0, Id=0;
  isquark=0;
  motherId1=motherId2=0;
  int up, upb, down, downb;
  up=upb=down=downb=0;
  nummuon=0;
  muCounter=0;
  const Candidate * mom[50];
  // if (!(primaryVtcs->size()>0))
  //  continue;
  Handle<GenParticleCollection> genParticles;
  event_.getByLabel("genParticles", genParticles);
  //istautau=0;
  //    cout<<"New Event!"<<endl;

  for(size_t i = 0; i < genParticles->size(); i++ )
    {


      const GenParticle & p = (*genParticles)[i];

    if(fabs(p.pdgId())==13 && p.status()==1)
	{

 
	  mom[1] = p.mother();
	  q=1;
	  Id= fabs(mom[q]->pdgId());
	  while(Id==13)
	    {
	      q++;
	      mom[q] = mom[q-1]->mother();
	      Id=fabs(mom[q]->pdgId());
	      
	    }
	  if(p.pt()>5 && Id==23)
	    {
	      nummuon++;
	      
	      if(nummuon==1)
		{
		  MC1=i;
		  MC1pthold=p.pt();
		}
	      if(nummuon==2)
		{
		  if(p.pt()>MC1pthold)
		    {
		      MC2=MC1;
		      MC1=i;
		      MC2pthold=MC1pthold;
		      MC1pthold=p.pt();
		    }
		  else
		    {
		      MC2=i;
		      MC2pthold=p.pt();
		      
		    }
		}
	    }
	}
    }
  if(nummuon>1)
	       {
		 //cout<<"a"<<endl;
     const GenParticle & MCmuon1b = (*genParticles)[MC1];
     const GenParticle & MCmuon2b = (*genParticles)[MC2];

    MCN1Q=MCmuon1b.charge();
    MCN2Q=MCmuon2b.charge();
	       

    MCD1Q=MCmuon1b.charge();
    MCD2Q=MCmuon2b.charge();
     //     cout<<MCmuon1a.px()<<endl;
    // this produces the naked particles
    MCN14vector->SetPxPyPzE(MCmuon1b.px(), MCmuon1b.py(), MCmuon1b.pz(), MCmuon1b.p4().t());  
    MCN24vector->SetPxPyPzE(MCmuon2b.px(), MCmuon2b.py(), MCmuon2b.pz(), MCmuon2b.p4().t());
    MCD14vector->SetPxPyPzE(MCmuon1b.px(), MCmuon1b.py(), MCmuon1b.pz(), MCmuon1b.p4().t());  
    MCD24vector->SetPxPyPzE(MCmuon2b.px(), MCmuon2b.py(), MCmuon2b.pz(), MCmuon2b.p4().t());
	       }


  for(size_t i = 0; i < genParticles->size(); i++ )
    {


      const GenParticle & p = (*genParticles)[i];
      
  if(fabs(p.pdgId())==15 && p.mother()->pdgId()==23)
	   {
	     	     istautau=1;
	     // cout<<"istrue"<<endl;
	   }      
      //this produces the two bare muons
      if(fabs(p.pdgId())==13 && p.mother()->pdgId()==23)
	{
	  muCounter++;
  
	  
	  if(muCounter==1)
	    {MC1a=i;
	      // cout<<"hi"<<MC1a<<endl;
	    }
	  
	  if(muCounter==2)
	    {MC2a=i;  
	      //	      cout<<"Dimuon"<<endl;
	    }
	}
  

      if(p.pdgId()==22&& (fabs(p.mother()->pdgId())==13 ||fabs(p.mother()->pdgId())==22) )
	{
	  if(MCN14vector->Pt()>0.0000001 && MCN24vector->Pt()>0.00000001)
	    {
	      if(DeltaRX(MCN14vector->Eta(), p.eta(), MCN14vector->Phi(), p.phi())<.1)
		{
		  // cout<<MCD14vector->Pt()<<" , "<<p.pt()<<endl;
		  MCD14vector->SetPxPyPzE((MCD14vector->Px()+p.px()),(MCD14vector->Py()+p.py()),(MCD14vector->Pz()+p.pz()),(MCD14vector->E()+p.energy()));
		  // cout<<p.pdgId()<<" , "<<p.pt()<<" , "<<p.status()<<" , "<<p.mother()->pdgId()<<endl;
		  // cout<<MCD14vector->Pt()<<endl;
		}
	      
	      if(DeltaRX(MCN24vector->Eta(), p.eta(), MCN24vector->Phi(), p.phi())<.1)
		{
		    MCD24vector->SetPxPyPzE((MCD24vector->Px()+p.px()),(MCD24vector->Py()+p.py()),(MCD24vector->Pz()+p.pz()),(MCD24vector->E()+p.energy()));
		}
	      
	    }
	}
	     
    }


  if(muCounter>1)
    {
    
      //cout<<"a"<<endl;
      const GenParticle & MCmuon1c = (*genParticles)[MC1a];
      const GenParticle & MCmuon2c = (*genParticles)[MC2a];
      
        MCB1Q=MCmuon1c.charge();
        MCB2Q=MCmuon2c.charge();
      
      
      //     cout<<MCmuon1a.px()<<endl;
      // this produces the naked particles
      MCB14vector->SetPxPyPzE(MCmuon1c.px(), MCmuon1c.py(), MCmuon1c.pz(), MCmuon1c.p4().t());  
      MCB24vector->SetPxPyPzE(MCmuon2c.px(), MCmuon2c.py(), MCmuon2c.pz(), MCmuon2c.p4().t());
	       
    }
	  if(MCN14vector->Pt()>0.0000001 && MCN24vector->Pt()>0.00000001 && MCB14vector->Pt()>0.0000001 && MCB24vector->Pt()>0.00000001)
	    {
	      if(MCD14vector->Pt()-MCB14vector->Pt()>0.1 && MCD14vector->Pt()-MCB24vector->Pt()>0.1)
		{
  cout<<MCB14vector->Eta()<<" , "<<MCD14vector->Eta()<<" , "<<MCN14vector->Eta()<<endl;
  cout<<MCB24vector->Eta()<<" , "<<MCD24vector->Eta()<<" , "<<MCN24vector->Eta()<<endl;
		}
	    }


		 //	     cout<<MCmuon1.pt()<<endl;

     //  cout<<p.pt()<<endl;
     
}    
*/

  //======================================================================//
 /*  
  void TupleMaker::load_acceptance(const Event& event_)
  {
    InputTag  vertexLabel(string("offlinePrimaryVertices"));
    Handle<reco::VertexCollection> thePrimaryVertexColl;
    
    event_.getByLabel(vertexLabel,thePrimaryVertexColl);
    // if (!(primaryVtcs->size()>0))
    Vertex dummy;
    const Vertex *pv = &dummy;
    if (thePrimaryVertexColl->size() != 0) {
      pv = &*thePrimaryVertexColl->begin();
    } 
    else { // create a dummy PV
      Vertex::Error e;
      e(0, 0) = 0.0015 * 0.0015;
      e(1, 1) = 0.0015 * 0.0015;
      e(2, 2) = 15. * 15.;
      Vertex::Point p(0, 0, 0);
      dummy = Vertex(p, e, 0, 0, 0);
    }
    int genmuon=0;
    int genelec=0;
    double MC1;
    double MC1pthold;
    double MC2;
    double MC2pthold;
    double Zmass=0;
    int q=0, Id=0;
    nummuon=0;
    const Candidate * mom[50];
    // if (!(primaryVtcs->size()>0))
    //  continue;
    Handle<GenParticleCollection> genParticles;
    event_.getByLabel("genParticles", genParticles);
    
    for(size_t i = 0; i < genParticles->size(); i++ )
      {
	const GenParticle & p = (*genParticles)[i];
	//	cout<<p.pdgId()<<" , "<<p.mother()->pdgId()<<endl;
	if(fabs(p.pdgId())==15 && p.mother()->pdgId()==23)
	  {
	    	cout<<p.pdgId()<<" , "<<p.mother()->pdgId()<<endl;

		    istautau=1;
		    cout<<istautau<<endl;
	    // cout<<"istrue"<<endl;
	  }  
	if(p.pdgId()==23 && p.status()==3)
	  {
	   if(p.mass()>Zmass)
	     Zmass=p.mass();
	   int n = p.numberOfDaughters();
	   //cout<<"New tau"<<endl;
	   for(int j = 0; j < n; j++) 
	     {
	       const Candidate * d = p.daughter( j );
	       int dauId = d->pdgId();
	       // cout<<dauId<<endl;
	       if(fabs(dauId)==13)
		 genmuon++;
	       if(fabs(dauId)==11)
		 genelec++;
	     }
	 }
     if((fabs(p.pdgId())==13 || fabs(p.pdgId())==11) && p.status()==1)
       {
	 mom[1] = p.mother();
	 q=1;
	 Id= fabs(mom[q]->pdgId());
	 while(Id==13 || Id==11)
	   {
	     q++;
	     mom[q] = mom[q-1]->mother();
	     Id=fabs(mom[q]->pdgId());
	     
	   }
	 if(p.pt()>5 && Id==23)
	   {
	     nummuon++;
	     
	 if(nummuon==1)
	   {
	     MC1=i;
	     MC1pthold=p.pt();
	   }
	   if(nummuon==2)
	     {
	       if(p.pt()>MC1pthold)
		 {
		   MC2=MC1;
		   MC1=i;
		   MC2pthold=MC1pthold;
		   MC1pthold=p.pt();
	       }
	       else
		 {
		   MC2=i;
		   MC2pthold=p.pt();
	       
		   }
	     }
	   
	     //   if(nummuon>2)
	     //  {  if(p.pt()>MC2pthold)
	     //		 {
	     //		  if(p.pt()>MC1pthold)
	     //		     {
	     //		       MC2pthold=MC1pthold;
	     //		       MC2=MC1;
	     //		       MC1=i;
	     //		       MC1pthold=p.pt();
	     //		     }
	     //		   else
	     //		     {
	     //		       MC2=i;
	     //		       MC2pthold=p.pt();
	     //		     }
	     //		 }

	   
       }
       }
     
      }
     
   
 
  
   //     cout<<genmuon<<" , "<<genelec<<endl<<endl;
   // if(Zmass>60 && Zmass<120)
   //	   {
   //this sets eff1 according to the lepton type that comes directly from the Z boson
    
	     if(genmuon==2)
	       eff1=1;
	     if(genelec==2)
	       eff1=3;
	     
	     if(genelec==1 && genmuon==1)
	       eff1=5;
	     //   }
     if(nummuon>1)
	       {
		 const GenParticle & MCmuon1 = (*genParticles)[MC1];
		 const GenParticle & MCmuon2 = (*genParticles)[MC2];
		 TLorentzVector v1(MCmuon1.px(), MCmuon1.py(), MCmuon1.pz(), MCmuon1.p4().E()), v2(MCmuon2.px(), MCmuon2.py(), MCmuon2.pz(), MCmuon2.p4().E());
		 
		 MC14vector->SetPxPyPzE(MCmuon1.px(), MCmuon1.py(), MCmuon1.pz(), MCmuon1.p4().t());  
		 MC24vector->SetPxPyPzE(MCmuon2.px(), MCmuon2.py(), MCmuon2.pz(), MCmuon2.p4().t());	 
		 MCinvmass = get_mass( MCmuon1.p4() + MCmuon2.p4() );
		 MC1eta=MCmuon1.eta();
		 MC1pt=MCmuon1.pt();
		 MC1phi=MCmuon1.phi();
		 MC1truth=MCmuon1.pdgId();
		 MC2truth=MCmuon2.pdgId();
		 MC2eta=MCmuon2.eta();
		 MC2pt=MCmuon2.pt();
		 MC2phi=MCmuon2.phi();
		 MC1Q=MCmuon1.charge();
		 MC2Q=MCmuon2.charge();
		 //	     cou<<MCmuon1.pt()<<endl;
		 
		 //     	 if(Zmass>60 && Zmass<120)
		 //	   {
		 //
		 //
		 //     cout<<genmuon<<" , "<<genelec<<endl;
		 //	       if((TMath::Pi()-v1.Angle(v2.Vect())>0.05))
		 //		   {
		 //		 cout<<"Here is one"<<endl;
		 //		 cout<<MCinvmass<<endl;
		 //		 cout<<Zmass<<endl;
		 //   }
		 //    if(fabs(MCmuon1.pdgId())==13 && fabs(MCmuon2.pdgId())==13)
		 //      eff2=1;
		 //    if(fabs(MCmuon1.pdgId())==11 && fabs(MCmuon2.pdgId())==11)
		 //      eff2=3;
		 //
		 //    if((fabs(MCmuon1.pdgId())==11 && fabs(MCmuon2.pdgId())==13) || (fabs(MCmuon2.pdgId())==11 && fabs(MCmuon1.pdgId())==13))
		 //      eff2=5;
		 //		   
		 //		   }
		 //		 
		 //     }
		 
		 //  cout<<p.pt()<<endl;
		 
    		 
	       

	       }
  }
 */
//=========================================================================//

void TupleMaker::load_prep()
{
  process_ = PROCESS_no;
  dilep_ = DILEP_no;
  lepjet_ = LEPJET_no;
  sign_ = SIGN_na;
  
  for(int j = 0; j < N_JETS; j++)
  {
    is_bjet[j] = false;
    jflavor[j] = ID_H;
    jdiscrim[j] = -100.;
    jet_vz[j] = -999.;
    jflvr_p[j] = ZERO_P4;
  }
  
  num_jets = num_bjets = 0;
  
  mET_x = mET_y = 0.;
  
  ++num_events;
  
  if(!(num_events%100))
  {
    cout << "[check]>\t " << "# " << num_events/1000 << "k ";
  
    time_stamp();
  }
}

//===================================================================//
/*
void TupleMaker::load_pileup(const Event& event_)
{
  int arrayNumber;
  arrayNumber=49;
    Handle<std::vector< PileupSummaryInfo > > PUInfo;
    event_.getByLabel(edm::InputTag("addPileupInfo"), PUInfo);
    std::vector<PileupSummaryInfo>::const_iterator iPV;

    for(iPV = PUInfo->begin(); iPV != PUInfo->end(); ++iPV){
      if (iPV->getBunchCrossing() == 0){

        nPUVertices = iPV->getPU_NumInteractions();
        nPUVerticesTrue = iPV->getTrueNumInteractions();
	arrayNumber= nPUVerticesTrue;

		if(arrayNumber>70)
		  weight =0;
		else
		  //weight = weights[arrayNumber];
			weight=1;
	// WJetsCorrection
	//		weight=weight*WJetsCorrection[arrayNumber];
      }
    }
    //	fill(h_nPUVertices, nPUVertices);
    //	fill(h_nPUVerticesTrue, nPUVerticesTrue);
	//	cout<<iPV<<endl;
      

    
}
*/

//===================================================================//
void TupleMaker::veto_event()
{
  process_ = PROCESS_no;
  dilep_ = DILEP_no;
  lepjet_ = LEPJET_no;
  
  num_jets = num_bjets = 0;
}

//===================================================================//
/*
void TupleMaker::load_beam(const Event& event_)
{
  event_.getByLabel(BEAMSPOT_COLLECTION, beam_spots);
  
  BeamSpot beam_spot = *beam_spots;
  
  beam_xyz = beam_spot.position();   //we'll need this below
}
*/
//===================================================================//
 
void TupleMaker::load_rMET(const Event& event_)
{
  //get raw MET
  
  // //GenMET::InvisibleEt()
  // //Handle<GenMETCollection> metHandle;
	//event_.getByLabel("genMet", metHandle);
  // //metHandle->at(0).px()
  // //vector<GenMET>
  
  //pf_mets
  const edm::InputTag eventMET("pfMetT0pcT1Txy", "", "Demo");

  event_.getByLabel(eventMET, MET);
 
      PFMETCollection::const_iterator imet = MET->begin(); 
  mET_x = imet->px();     //raw x-met
  mET_y = imet->py();     //raw y-met
  
    metP4->SetPxPyPzE(imet->px(), imet->py(), imet->pz(), imet->p4().t());
    //  cout<<metP4.Px()<<" , "<<metP4.Py()<<endl;
  update_MET();
 
  // fill(h_raw_MET, mET);
  // fill(h_raw_METphi, mETphi);
}
 
//===================================================================//
  
void TupleMaker::update_MET()
{
  vec_mET = math::XYZVector(mET_x, mET_y, 0.);
  
  mET = get_MET();
  mETphi = get_METphi();
}
  
//===================================================================//
   
double TupleMaker::get_MET()
{
  return sqrt(mET_x*mET_x + mET_y*mET_y);
}
   
//===================================================================//
    
double TupleMaker::get_METphi()       //-pi to pi
{
  return atan2(mET_y, mET_x);    //y, x
}
    
//===================================================================//


//===================================================================//

int TupleMaker::load_elec(const Event& event_, const EventSetup& setup_, int num_m)
{
  int num_e=0;
  double phil, etal, eta1;
  double elec1pt;
  double Aeff=0;
  // evaluate_mvas(event_, setup_);

  Handle<edm::ValueMap<float>> mvaTrigV0_handle;
    event_.getByLabel("mvaTrigV0", mvaTrigV0_handle);
    const edm::ValueMap<float> ele_mvaTrigV0 = (*mvaTrigV0_handle.product());

  InputTag gsfEleLabel(string("gsfElectrons"));
  Handle<GsfElectronCollection> theEGammaCollection;
  event_.getByLabel(gsfEleLabel,theEGammaCollection);
  const GsfElectronCollection theEGamma = *(theEGammaCollection.product());

  edm::Handle<reco::PhotonCollection> photonH;
  event_.getByLabel(inputTagPhotons_,photonH);
   
  Int_t eee=0;

  //event_.getByLabel(genLable,genParticles);
  //InputTag  mcTruthLabel(string("generator"));
  //edm::Handle<edm::HepMCProduct> pMCTruth;
  //iEvent.getByLabel(mcTruthLabel,pMCTruth);
  //const HepMC::GenEvent* genEvent = pMCTruth->GetEvent();

  //electron regression was removed temporarily for use with the singal file 8/27 some more lines need to be undeleted below as well

  /*
  edm::Handle<edm::ValueMap<double>> regEne_handle;
  event_.getByLabel(edm::InputTag("eleRegressionEnergy","eneRegForGsfEle", "HLT"), regEne_handle);
  const edm::ValueMap<double> ele_regEne = (*regEne_handle.product());

edm::Handle<edm::ValueMap<double>> regErr_handle;
    event_.getByLabel(edm::InputTag("eleRegressionEnergy","eneErrorRegForGsfEle"), regErr_handle);
    const edm::ValueMap<double> ele_regErr = (*regErr_handle.product());
  */


  InputTag  vertexLabel(string("offlinePrimaryVertices"));
  Handle<reco::VertexCollection> thePrimaryVertexColl;
  event_.getByLabel(vertexLabel,thePrimaryVertexColl);



  _Rho=0;
  edm::Handle<double> rhoPtr;
  const edm::InputTag eventrho("kt6PFJets", "rho");
  event_.getByLabel(eventrho,rhoPtr);
  _Rho=*rhoPtr;
  

  
  edm::Handle<reco::ConversionCollection> hConversions;
  event_.getByLabel("allConversions", hConversions);
  
  Vertex dummy;
  const Vertex *pv = &dummy;
  if (thePrimaryVertexColl->size() != 0) {
    pv = &*thePrimaryVertexColl->begin();
  } else { // create a dummy PV
    Vertex::Error e;
    e(0, 0) = 0.0015 * 0.0015;
    e(1, 1) = 0.0015 * 0.0015;
    e(2, 2) = 15. * 15.;
    Vertex::Point p(0, 0, 0);
    dummy = Vertex(p, e, 0, 0, 0);
  }
  
  InputTag  reducedEBRecHitCollection(string("reducedEcalRecHitsEB"));
  InputTag  reducedEERecHitCollection(string("reducedEcalRecHitsEE"));

  EcalClusterLazyTools lazyTools(event_, setup_, reducedEBRecHitCollection, reducedEERecHitCollection);
  
  edm::ESHandle<TransientTrackBuilder> builder;
  setup_.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);
  TransientTrackBuilder thebuilder = *(builder.product());
  

  bool debug = true;
  bool debugMVAclass = false;
  bool debugMyVar = false;


   event_.getByLabel("offlineBeamSpot", bsHandle);
   const reco::BeamSpot &beamspot = *bsHandle.product();

   typedef std::vector< edm::Handle< edm::ValueMap<reco::IsoDeposit> > > IsoDepositMaps;
   typedef std::vector< edm::Handle< edm::ValueMap<double> > > IsoDepositVals;
  unsigned nTypes=3;
  IsoDepositMaps electronIsoDep(nTypes);

  for (size_t j = 0; j<inputTagIsoDepElectrons_.size(); ++j) {
    event_.getByLabel(inputTagIsoDepElectrons_[j], electronIsoDep[j]);
  }


   IsoDepositMaps photonIsoDep(nTypes);
    for (size_t j = 0; j<inputTagIsoDepPhotons_.size(); ++j) {
      event_.getByLabel(inputTagIsoDepPhotons_[j], photonIsoDep[j]);
    }
  IsoDepositVals electronIsoValPFId(nTypes);
  IsoDepositVals photonIsoValPFId(nTypes);
     const IsoDepositVals * electronIsoVals =  &electronIsoValPFId  ;

  for (size_t j = 0; j<inputTagIsoValElectronsPFId_.size(); ++j) {
    event_.getByLabel(inputTagIsoValElectronsPFId_[j], electronIsoValPFId[j]);
  }

  for (size_t j = 0; j<inputTagIsoValPhotonsPFId_.size(); ++j) {
    event_.getByLabel(inputTagIsoValPhotonsPFId_[j], photonIsoValPFId[j]);
  }

 
      for (uint j=0; j<theEGamma.size();j++)
  {

    eee++;
	  bool elePresel = trainTrigPresel(theEGamma[j]);

	  double mvaTrigMthd1 = -11.;
	  double mvaTrigMthd2 = -11.;

	  double mvaTrigNonIp = -11.;
	  
	  double mvaNonTrigMthd1 = -11;
	  double mvaNonTrigMthd2 = -11;

	  //	  mvaNonTrigMthd1 = myMVANonTrigV0->mvaValue((theEGamma[j]),*pv,thebuilder,lazyTools,debugMVAclass);


		  //	cout<<mvaNonTrigMthd1<<endl;
    //
        // get particle flow isolation
	

      // cout<<ie.id()<<endl;
           phil = theEGamma[j].phi();
	   etal = theEGamma[j].eta();

  
	  

	         reco::GsfElectronRef myElectronRef(theEGammaCollection,j);


  if((DeltaRX(theEGamma[j].eta(), MC1eta, theEGamma[j].phi(), MC1phi)<.1) || (DeltaRX(theEGamma[j].eta(), MC2eta, theEGamma[j].phi(), MC2phi)<.1) && passacc)
		       {

			eleceff3++;
		       }

		 if(theEGamma[j].pt()>10 && fabs(eta1)<2.4)
		   {
		     // cout<<extraElectrons<<endl;
		     		   
		     if(elePresel)
		       {

		     //	     mvaNonTrigMthd1 = myMVANonTrigV0->mvaValue((theEGamma[j]),*pv,thebuilder,lazyTools,debugMVAclass);
			 //  mvaTrigMthd1 = myMVATrigV0->mvaValue((theEGamma[j]),*pv,thebuilder,lazyTools,debugMVAclass);
			 double trigValue=ele_mvaTrigV0.get(j);
			 //	 cout<<trigValue<<endl;
			 if((theEGamma[j].pt()<20&&((trigValue>0.00 && fabs(eta1)<.8) || (trigValue>.10 && (fabs(eta1)>.8 && fabs(eta1)<1.479)) || (trigValue>.62 && fabs(eta1)>1.479)))||((trigValue>.94 && fabs(eta1)<.8) || (trigValue>.85 && (fabs(eta1)>.8 && fabs(eta1)<1.479)) || (trigValue>.92 && fabs(eta1)>1.479)))
		    {
		      // cout<<"hello"<<endl;
		     double charged =  (*(*electronIsoVals)[0])[myElectronRef];
		     double photon = (*(*electronIsoVals)[1])[myElectronRef];
		     double neutral = (*(*electronIsoVals)[2])[myElectronRef];
		     if(abs(etal)>2.4)
		       Aeff=0.261;
		     if(abs(etal)>2.3 && abs(etal)<2.4)
		       Aeff=0.194;
		     if(abs(etal)>2.2 && abs(etal)<2.3)
		       Aeff=0.183;
		     if(abs(etal)>2.0 && abs(etal)<2.2)
		       Aeff=0.143;
		     if(abs(etal)>1.479 && abs(etal)<2.0)
		       Aeff=0.115;
		     if(abs(etal)>1.0 && abs(etal)<1.479)
		       Aeff=0.209;
		     if(abs(etal)<1.0)
		       Aeff=.208;





		   	     


		     double iso =(charged+max(photon+neutral-_Rho*Aeff,0.0))/theEGamma[j].pt();
	   //cout<<charged<<" ,"<<photon<<" ,"<<neutral<<" the iso"<<iso<<endl;
		     fill(h_num_good_mj, iso);

		     // cout<<_Rho<<" , "<<Aeff<<" , "<<iso<<endl;
	   if(iso<.15)
	     {
	       //calculate sip
	       float ip3d    = -999.0;
	       float ip3derr = 1.0;
	       float ip3dSig = 0.0;
	       
	       if (theEGamma[j].gsfTrack().isNonnull()) {
		 const double gsfsign = ( (-theEGamma[j].gsfTrack()->dxy(pv->position())) >=0 ) ? 1. : -1.;
		 
		 
		 const reco::TransientTrack &tt = thebuilder.build(theEGamma[j].gsfTrack());
		 
		 const std::pair<bool,Measurement1D> &ip3dpv = IPTools::absoluteImpactParameter3D(tt,*pv);
		 if (ip3dpv.first) {
		   ip3d = gsfsign*ip3dpv.second.value();
		   ip3derr = ip3dpv.second.error();
          ip3dSig = ip3d/ip3derr;
		 }
	       }
		   int	misshits = theEGamma[j].gsfTrack()->trackerExpectedHitsInner().numberOfHits();
		   //cout<<misshits<<endl;
		   if(misshits==0)
		     {

		       bool passconversionveto = !ConversionTools::hasMatchedConversion(theEGamma[j],hConversions,beamspot.position());
		       if(passconversionveto)
			 {
		       // cout<<ip3d<<" ,"<<ip3derr<<" , "<<ip3dSig<<endl;
		       //  if(ip3dSig<4)
		       //	 {
			     if(!(cut_mu_electron(theEGamma[j].eta(), theEGamma[j].phi())))
			    {
		   //expected inner hits should be <=1

		   // elec1P4.SetPxPyPzE(theEGamma[j].px(), theEGamma[j].py(), theEGamma[j].pz(), theEGamma[j].p4().t());

			   //   cout<<elec1P4.Px()<<" , "<<elec1P4.Py()<<" , "<<theEGamma[j].pt()<<endl;
			      
			      if((DeltaRX(theEGamma[j].eta(), MC1eta, theEGamma[j].phi(), MC1phi)<.1) || (DeltaRX(theEGamma[j].eta(), MC2eta, theEGamma[j].phi(), MC2phi)<.1) && passacc)
		       {
			eleceff4++;
		       }
		    num_e++;
		    double ene=0.0;
		    double err=0.0;
		    //		    double ene=ele_regEne.get(eee-1);
		    // double err =ele_regErr.get(eee-1);
		    //	    cout<<theEGamma[j].p4().t()<<" , "<<ene<<endl;
		    // cout<<num_e<<endl;
			     if(num_e==1)
			       {elecReg1=ene;
			       eiso1=iso;
			       elec1pt=theEGamma[j].pt();
			       elec1P4.SetPxPyPzE(theEGamma[j].px(), theEGamma[j].py(), theEGamma[j].pz(), theEGamma[j].p4().t());
			       elec1eta=theEGamma[j].eta();
			       elec1Q=theEGamma[j].charge();
			       elec1phi=theEGamma[j].phi();
			       elec1dxy=fabs(theEGamma[j].gsfTrack()->dxy(pv->position()));
			       elec1track=theEGamma[j].gsfTrack();
			       //       cout<<elec1P4.Px()<<endl;
				 }
			   if(num_e==2)
			     {
			       if(elec1pt<theEGamma[j].pt())
				 {
				   elecReg2=elecReg1;
				   elecReg1=ene;
				   //   elec2pt=elec1pt;
				   // elec1pt=theEGamma[j].pt();
				   eiso2=eiso1;
				   eiso1=iso;
				 elec2P4.SetPxPyPzE(elec1P4.Px(), elec1P4.Py(), elec1P4.Pz(), elec1P4.E());
				 elec1P4.SetPxPyPzE(theEGamma[j].px(), theEGamma[j].py(), theEGamma[j].pz(), theEGamma[j].p4().t());
				 elec2eta=elec1eta;
				 elec2phi=elec1phi;
				 elec1eta=theEGamma[j].eta();
				 elec1phi=theEGamma[j].phi();
				 elec2dxy=elec1dxy;
				 elec1dxy=fabs(theEGamma[j].gsfTrack()->dxy(pv->position()));
				 elec2Q=elec1Q;
			       elec1Q=theEGamma[j].charge();
			       elec2track=elec1track;
			       elec1track=theEGamma[j].gsfTrack();

				 //	 cout<<elec1P4.Px()<<" , "<<elec2P4.Px()<<endl;
				 }
				 else
				 {
				   elecReg2=ene;
				   eiso2=iso;
				   elec2P4.SetPxPyPzE(theEGamma[j].px(), theEGamma[j].py(), theEGamma[j].pz(), theEGamma[j].p4().t());
		       elec2eta=theEGamma[j].eta();
		       elec2phi=theEGamma[j].phi();
		        elec2dxy=fabs(theEGamma[j].gsfTrack()->dxy(pv->position()));
			elec2Q=theEGamma[j].charge();
			elec2track=theEGamma[j].gsfTrack();

		       //   elec2pt=theEGamma[j].pt();
		
		  
				   //   cout<<"else statement  "<<elec1P4.Px()<<" , "<<elec2P4.Px()<<endl;

				 }
			     }
			     }
		     }
			     //	   dReEm1[extraElectrons] = DeltaRX(etal, eta1, phil, phi1);
			     //	   dReEm2[extraElectrons] = DeltaRX(etal, eta2, phil, phi2);
		   //		   if(dReEm1<.01 || dReEm2<.01)
		   //  {
			     //		   Eelectron_pt[extraElectrons]=theEGamma[j].pt();
			     //	   Eelectron_eta[extraElectrons]=theEGamma[j].eta();
		   //	      MC_info(etal, phil,Eelectron_pt[extraElectrons], event_);
			     //		   ++extraElectrons;
		   //  }
			 }
		     }  
	     }      
		   }
		  
		   }
		 //}
  }

     
  return num_e;  
}

//===================================================================//
/*
bool TupleMaker::cut_electron(vector<GsfElectron>::const_iterator ele)
{
    bool isEB           = ele.isEB() ? true : false;
    float pt            = ele.pt();
    float eta           = ele.superCluster()->eta();

    // id variables
    float dEtaIn        = ele.deltaEtaSuperClusterTrackAtVtx();
    float dPhiIn        = ele.deltaPhiSuperClusterTrackAtVtx();
    float sigmaIEtaIEta = ele.sigmaIetaIeta();
    float hoe           = ele.hadronicOverEm();
    float ooemoop       = (1.0/ele.ecalEnergy() - ele.eSuperClusterOverP()/ele.ecalEnergy());

    // impact parameter variables
    float d0vtx         = 0.0;
    float dzvtx         = 0.0;
    if (vtxs->size() > 0) {
        reco::VertexRef vtx(vtxs, 0);    
        d0vtx = ele.gsfTrack()->dxy(vtx->position());
        dzvtx = ele.gsfTrack()->dz(vtx->position());
    } else {
        d0vtx = ele.gsfTrack()->dxy();
        dzvtx = ele.gsfTrack()->dz();
    }

 float eleISOendcap = (thisElec->TrkIso() + thisElec->EmIso() + thisElec->HadIso() - rhoFactor*TMath::Pi()*0.09)/thisElec->Pt(); 
        float eleISObarrel = (thisElec->TrkIso() + thisElec->EmIso() + thisElec->HadIso() - rhoFactor*TMath::Pi()*0.09)/thisElec->Pt(); 
    unsigned int mask = 0;
    float cut_dEtaIn[2]         = {999.9, 999.9};
    float cut_dPhiIn[2]         = {999.9, 999.9};
    float cut_sigmaIEtaIEta[2]  = {999.9, 999.9};
    float cut_hoe[2]            = {999.9, 999.9};
    float cut_ooemoop[2]        = {999.9, 999.9};
    float cut_d0vtx[2]          = {999.9, 999.9};
    float cut_dzvtx[2]          = {999.9, 999.9};
    float cut_iso[2]            = {999.9, 999.9};
    bool cut_vtxFit[2]          = {false, false};
    unsigned int cut_mHits[2]   = {999, 999};



    else if (workingPoint == EgammaCutBasedEleId::LOOSE) {
        cut_dEtaIn[0]        = 0.007; cut_dEtaIn[1]        = 0.009;
        cut_dPhiIn[0]        = 0.150; cut_dPhiIn[1]        = 0.100;
        cut_sigmaIEtaIEta[0] = 0.010; cut_sigmaIEtaIEta[1] = 0.030;
        cut_hoe[0]           = 0.120; cut_hoe[1]           = 0.100;
        cut_ooemoop[0]       = 0.050; cut_ooemoop[1]       = 0.050;
        cut_d0vtx[0]         = 0.020; cut_d0vtx[1]         = 0.020;
        cut_dzvtx[0]         = 0.200; cut_dzvtx[1]         = 0.200;
        cut_vtxFit[0]        = true ; cut_vtxFit[1]        = true;
        cut_mHits[0]         = 1    ; cut_mHits[1]         = 1;
        if (pt >= 20.0) {
            cut_iso[0] = 0.150; cut_iso[1] = 0.150;
        }
        else {
            cut_iso[0] = 0.150; cut_iso[1] = 0.100;
        }
    } 

    unsigned int idx = isEB ? 0 : 1;



  bool dowecut = (get_e_iso(ie) > CUT_ELEC_ISO);  

  double d_xy = fabs(ie->gsfTrack()->dxy(beam_xyz));
  
  dowecut = (dowecut || (fabs(ie->eta()) > CUT_ELEC_ETA) );
  dowecut = (dowecut || (ie->pt() < CUT_LEPTON_LOOSE_PT) );
  
  if(!dowecut)          
  else                  return (true);
    
  dowecut = ( (d_xy > CUT_ELEC_DXY) || cut_mu_electron(ie) );  
  
  return (dowecut);
}
*/
//===================================================================//
/*
bool TupleMaker::cut_tight_electron(int num_electrons)
{
  if(num_electrons < 1)
    return (true);
  
  double pT = elec_1->pt();
  
  return (pT < CUT_LEPTON_TRIG_PT);
}
*/
//===================================================================//
/*
double TupleMaker::get_e_iso(vector<GsfElectron>::const_iterator ie)
{
  //not ready for prime time:
  //  double iso = ie->isolationR03().sumPt, pT = ie->pt();
  //  iso += ie->isolationR03().emEt + ie->isolationR03().hadEt;
  
  double pT = ie->pt(), iso = ie->dr03TkSumPt();
  
  iso += ie->dr03EcalRecHitSumEt() + ie->dr03HcalTowerSumEt();
  
  //IsolationVariables()
  
  if(pT == 0.)
    return 999.;
  
  return ( iso / pT );
}
*/
//===================================================================//

bool TupleMaker::cut_mu_electron(double eeta,double ephi)
{
  int numVH;
  
  double dR, mphi, meta;
  
  vector<Muon>::const_iterator im;
  
  for(im = muons->begin(); im != muons->end(); ++im) 
  {
    if(!im->isGlobalMuon())           continue;
    if(!im->isTrackerMuon())          continue;
    
    numVH = im->globalTrack()->hitPattern().numberOfValidTrackerHits();
    
    if(numVH < CUT_TRACK_VHITS)       continue;
    
    mphi = im->phi();                 meta = im->eta();
    
    dR = DeltaRX(eeta, meta, ephi, mphi);
    
    if(dR < MATCH_DELTA_R_E_MU)       return (true);
  }
  
  return (false);
}

//===================================================================//

int TupleMaker::load_muon(const Event& event_)
{
  //print_debug("load_muon", N_DEBUG);  

  int num_m = 0;

  event_.getByLabel(MUON_COLLECTION, muons);
   Handle<reco::VertexCollection> primaryVtcs;
  event_.getByLabel(VERTEX_COLLECTION, primaryVtcs);
  
  if (!(primaryVtcs->size()>0))
  return 0;

  vector<Muon>::const_iterator imuon;
  //  cout<<"New Event"<<endl;
  //  cout<<"New event"<<" "<<eff1<<endl;
  for(imuon = muons->begin(); imuon != muons->end(); ++imuon) 
    {
      if(imuon->isGlobalMuon())
	{
      //  cout<<DeltaRX(imuon->eta(), MC1eta, imuon->phi(), MC1phi)<<" , "<<DeltaRX(imuon->eta(), MC2eta, imuon->phi(), MC2phi)<<endl;
      if((DeltaRX(imuon->eta(), MC1eta, imuon->phi(), MC1phi)<.1)&& fabs(MC1truth)==13 && passacc)
	    {
	      //    cout<<muoneff3<<endl;
			muoneff3++;

			h_MCmu1Global->Fill(imuon->isGlobalMuon(), weight);
			h_MCmu1PF->Fill(imuon->isPFMuon(), weight);
	  if(imuon->isGlobalMuon())
	    {
				h_MCmu1dxy->Fill(fabs(imuon->muonBestTrack()->dxy(primaryVtcs->begin()->position())),weight);
			h_MCmu1dz->Fill(fabs(imuon->muonBestTrack()->dz(primaryVtcs->begin()->position())),weight);		
			h_MCmu1Chi2->Fill(imuon->globalTrack()->normalizedChi2(), weight);
			h_MCmu1hits->Fill(imuon->globalTrack()->hitPattern().numberOfValidMuonHits(), weight);
			h_MCmu1Stations->Fill(imuon->numberOfMatchedStations(),weight);
			h_MCmu1Phits->Fill(imuon->innerTrack()->hitPattern().numberOfValidPixelHits(),weight);
			h_MCmu1Thits->Fill(imuon->innerTrack()->hitPattern().trackerLayersWithMeasurement(),weight);
			h_MCmu1iso->Fill(get_m_iso(imuon), weight);
	    }


		       }

      if(DeltaRX(imuon->eta(), MC2eta, imuon->phi(), MC2phi)<.1 && fabs(MC2truth)==13 && passacc)
	{
	  muoneff3++;

	  h_MCmu2Global->Fill(imuon->isGlobalMuon(), weight);
	  h_MCmu2PF->Fill(imuon->isPFMuon(), weight);
	  if(imuon->isGlobalMuon())
	    {
	  h_MCmu2dxy->Fill(fabs(imuon->muonBestTrack()->dxy(primaryVtcs->begin()->position())),weight);
	  h_MCmu2dz->Fill(fabs(imuon->muonBestTrack()->dz(primaryVtcs->begin()->position())),weight);
	  h_MCmu2Chi2->Fill(imuon->globalTrack()->normalizedChi2(), weight);
	  h_MCmu2hits->Fill(imuon->globalTrack()->hitPattern().numberOfValidMuonHits(), weight);
	  h_MCmu2Stations->Fill(imuon->numberOfMatchedStations(),weight);
	  h_MCmu2Phits->Fill(imuon->innerTrack()->hitPattern().numberOfValidPixelHits(),weight);
	  h_MCmu2Thits->Fill(imuon->innerTrack()->hitPattern().trackerLayersWithMeasurement(),weight);
	  h_MCmu2iso->Fill(get_m_iso(imuon), weight);
	    }
	}
      if(cut_muon(imuon))
	{


	    	 if(muon::isTightMuon(*imuon, *primaryVtcs->begin())&& (get_m_iso(imuon) < .20))       
	    //	    if(cut_muon(imuon)&& fabs(imuon->muonBestTrack()->dxy(primaryVtcs->begin()->position()))<0.2 && fabs(imuon->muonBestTrack()->dz(primaryVtcs->begin()->position()))<0.5)
		   {
	  //if((DeltaRX(imuon->eta(), MC1eta, imuon->phi(), MC1phi)<.1) || (DeltaRX(imuon->eta(), MC2eta, imuon->phi(), MC2phi)<.1))
	  //	    {
	  //	  if((DeltaRX(MC1eta, imuon->eta(), MC1phi, imuon->phi())<.1) ||(DeltaRX(MC2eta, imuon->eta(), MC2phi, imuon->phi())<.1))
	    
	       if((DeltaRX(imuon->eta(), MC1eta, imuon->phi(), MC1phi)<.1) || (DeltaRX(imuon->eta(), MC2eta, imuon->phi(), MC2phi)<.1) && passacc)
		       {
			muoneff4++;
		       }
	
	       ++num_m;
	      
	      order_muons(num_m, imuon);    
		   }
	}
	}
    }
	

  return num_m;
}
//===================================================================//

bool TupleMaker::cut_muon(vector<Muon>::const_iterator im)
{
   if(!im->isGlobalMuon())     return (false);
   if(!im->isPFMuon())    return (false);

  double  
    normChi2 = im->globalTrack()->normalizedChi2();
  int numPH2 = im->innerTrack()->hitPattern().numberOfValidPixelHits(); 

  int hitPattern = im->globalTrack()->hitPattern().numberOfValidMuonHits();
  int numberStations =  im->numberOfMatchedStations();
  int hits = im->track()->hitPattern().trackerLayersWithMeasurement(); 


  ///  int numVH = im->globalTrack()->hitPattern().numberOfValidTrackerHits();
  //  int numPH = im->globalTrack()->hitPattern().numberOfValidPixelHits();
  //change to 2,4 and 10 for loose cuts
  bool dowecut = ( (fabs(im->eta()) < 2.4) );
  dowecut = ( dowecut && (im->pt() >= 7) );
  //   dowecut = ( dowecut && (get_m_iso(im) < .2) );
   dowecut = ( dowecut && (normChi2 < 10) );
      dowecut = ( dowecut && (hitPattern > 0) );
       dowecut = ( dowecut && (numberStations > 1) );
   //the d_xy cut is bad, it needs to be lower but i think i need to fix things first  
   //    dowecut = ( dowecut && (d_xy < 1) );  
       dowecut = ( dowecut && (numPH2 > 0) ); 
    dowecut = ( dowecut && (hits > 5) ); 
  // //   if(!dowecut)                fill(h_raw_muon_dxy, d_xy);
  //    //  dowecut = ( dowecut && (d_xy < CUT_MUON_DXY) );
  return (dowecut);
}

//===================================================================//

      bool TupleMaker::cut_tight_muon(vector<Muon>::const_iterator im)
{

  if(!im->isGlobalMuon())     return (false);
  cout<<"Still going"<<endl;
  // if(!im->isTrackerMuon())    return (false);
 bool dowecut = ( (fabs(im->eta()) < 2.4) );
   dowecut = ( dowecut && (im->pt() >= 10) );
   //add isolation cut

  dowecut = ( dowecut && (get_m_iso(im) < .2) );
   // dowecut = ( dowecut && (normChi2 <10) );
  //dowecut = ( dowecut && (numVH >= CUT_TRACK_VHITS) );
  //dowecut = ( dowecut && (numPH > 0) );  
  // if(!dowecut)                fill(h_raw_muon_dxy, d_xy);
  //dowecut = ( dowecut && (d_xy < CUT_MUON_DXY) );
  return (dowecut);
}

//===================================================================//

double TupleMaker::get_m_iso(vector<Muon>::const_iterator im)
{
  double pt = im->pt(), iso = im->pfIsolationR04().sumChargedHadronPt;
    
    iso += max(im->pfIsolationR04().sumPhotonEt + im->pfIsolationR04().sumNeutralHadronEt-0.5*im->pfIsolationR04().sumPUPt, 0.0);
  //      iso += max(im->pfIsolationR04().sumNeutralHadronEt-0.5*im->pfIsolationR04().sumPUPt, 0.0);
  
  if(pt == 0.)   
    {   return (100); 
    }
  return (iso/pt);      //relative iso
}

//===================================================================//

void TupleMaker::order_muons(int num_m, vector<Muon>::const_iterator im)
{
  if(num_m < 1)
    cout << "[error]>\t num_muons = " << num_m << endl;
  else if(num_m == 1)          
    muon_1 = im;
  else if(num_m == 2)
  {
    if(im->pt() > muon_1->pt())
    {
      muon_2 = muon_1;
      muon_1 = im;
    }
    else
      muon_2 = im;
  }
  else if(num_m > 2)
  {
    if(im->pt() > muon_1->pt())
    {
      muon_2 = muon_1;
      muon_1 = im;
    }
    else if(im->pt() > muon_2->pt())
      muon_2 = im;
  }
}

//===================================================================//

void TupleMaker::load_vtxs(const Event& event_)
{ 
  //  if(process_ == PROESS_no)    return;
  vtxCount=0;
  //event_.getByLabel(VERTEX_COLLECTION, vertices);
   Handle<reco::VertexCollection> primaryVtcs;
  event_.getByLabel(VERTEX_COLLECTION, primaryVtcs);

  for(VertexCollection::const_iterator iVtx = primaryVtcs->begin(); iVtx!= primaryVtcs->end(); ++iVtx){
    reco::Vertex myVtx = reco::Vertex(*iVtx);
    if(myVtx.isValid() && !(myVtx.isFake())
                && myVtx.ndof()         > 4.  // Number of degrees of freedom
       && fabs(myVtx.z())      <= 24.) // Longitudinal distance from IP
       //     && fabs(myVtx.perp())   <= 2.) // Transverse distance from IP continue;

      {
    //TCPrimaryVtx* vtxCon = new ((*primaryVtx)[vtxCount]) TCPrimaryVtx;
    // vtxCon->SetXYZ(myVtx.x(), myVtx.y(), myVtx.z());
    //vtxCon->SetNDof(myVtx.ndof());
    //vtxCon->SetChi2(myVtx.chi2());
    //vtxCon->SetNtracks(myVtx.nTracks());
    //vtxCon->SetSumPt2Trks(sumPtSquared(myVtx));
    //vtxCon->SetIsFake(myVtx.isFake());
    ++vtxCount;

      }
  }

  //see Andy's code to get lepton vtx
  //VertexCollection::const_iterator ivtx, vtx_SumPt;
  //Vertex::trackRef_iterator it;
 //TrackRefVector jetTRKs = j->getTrackRefs();
  //TrackRefVector::const_iterator jetTRK;
  
  //VertexCollection::const_iterator vtxF = vertices->begin();
  //VertexCollection::const_iterator vtxL = vertices->end();
  //for(VertexCollection::const_iterator vtx = vtxF; vtx != vtxL; ++vtx){
  //  for(Vertex::trackRef_iterator it = vtx->tracks_begin(); it != vtx->tracks_end(); ++it) {
  //    const Track &muTrack = *(mu->track());
  //    const edm::RefToBase<Track> &myTrackRef = *it;
  //    if(myTrackRef.isAvailable()){
  //      const Track &myVertexTrack = *myTrackRef.get();
  //      if (&muTrack == &myVertexTrack) { lepCon->SetpVtx(vtx->x(),vtx->y(),vtx->z()); }
  //    }
  //  }
  //}
  
  //dumb way for now:
  
  //double dz = 999.;
  
  //  if(process_ != PROCESS_di)    return;

  // switch(dilep_)
  // {
  //  case DILEP_ee:
  //    dz = fabs( elec_1->vz() - elec_2->vz() );
  //    break;
  //  case DILEP_mm:
  //    dz = fabs( muon_1->vz() - muon_2->vz() );
  //    break;     
  //  case DILEP_em:
  //    dz = fabs( elec_1->vz() - muon_1->vz() );
  //    break;      
  //  default:
  //    cout << "[error]>\t dilep_ = " << dilep_ << endl;
  //    break;
  //  }
    
  //  if(dz > CUT_DELTAZ)    {  veto_event();    ++num_lldz_veto;  }
}

//===================================================================//
 /*
double TupleMaker::get_lep_vz()
{
  if(process_ == PROCESS_di)    
    switch(dilep_) 
    {
      case DILEP_ee:    return (elec_1->vz());
      case DILEP_mm:    return (muon_1->vz());
      case DILEP_em:    return ( (elec_1->vz() + muon_1->vz())/2. );
      default:          break;
    }
  else if(process_ == PROCESS_lj)
    switch (lepjet_) 
    {
      case LEPJET_ej:   return (elec_1->vz());
      case LEPJET_mj:   return (muon_1->vz());
      default:          break;
    }
  
  return (999.);
}
 */
//===================================================================//

/*void TupleMaker::load_corr(const Event& event_, const EventSetup& setup_)
{
  if(process_ == PROCESS_no)      return;
  
  setup_.get<JetCorrectionsRecord>().get(CORR_COLLECTION, jet_corrections);

  corr_L1 = JetCorrector::getJetCorrector(JET_CORR_L1, setup_);
  corr_L2 = JetCorrector::getJetCorrector(JET_CORR_L2, setup_);
  corr_L3 = JetCorrector::getJetCorrector(JET_CORR_L3, setup_);
  
  if(IS_DATA)
    corr_MC = JetCorrector::getJetCorrector(JET_CORR_MC, setup_);
} 
*/
//===================================================================//

void TupleMaker::load_jets(const Event& event_, const EventSetup& setup_)
{  
  if(process_ == PROCESS_no)      return;
  
  event_.getByLabel(JET_COLLECTION, jets);
   
   for(int i =0; i<10; i++)
       {jet_pt[i]=0.0;
       } 
  PFJet cjet;
  //  cout<<endl;
  PFJetCollection::const_iterator ijet;
  //cout<<endl;
  for(ijet = jets->begin(); ijet != jets->end(); ++ijet)
  {

    double jpt=0;
    jpt=ijet->pt();

    //if(jpt>10)
      //cout<<jpt<<endl;

    // if(lepton_jet(ijet))          continue;
    
     cjet = correct_jet(ijet, event_, setup_);    


     //const JetCorrector* corrector = JetCorrector::getJetCorrector (JetCorrectionService,setup_);   //Get the jet corrector from the event setup
     //PFJet  correctedJet = *cjet;                                 //copy original jet
     //double jec = corrector->correction(*i_jet,jetRef,event_,setup_);
     //correctedJet.scaleEnergy(jec);        
     //cout<<correctedJet.pt()<<endl;

                // apply the correction
     if(cut_jet(cjet))           
       {
	 
	 // if(cut_jet_vtx(cjet))         continue;
	 if(noniso_jet(cjet))      
	   {
	     jet_pt[num_jets] = cjet.pt();
	     ++num_jets;

	     //cout<<sumJet.Px()<<" , "<<sumJet.Py()<<" , "<<sumJet.Pz()<<" , "<<sumJet.E()<<endl;
	     sumJet.SetPxPyPzE(sumJet.Px()+cjet.px(), sumJet.Py()+cjet.py(),sumJet.Pz()+cjet.pz(), sumJet.E()+cjet.p4().t());
	  
	   }

       }
     else{
       
      
     }
    //      order_jets(num_jets, cjet);
    
    //jet_addons(ijet - jets->begin());
    // jet_vz[num_jets-1] = jet_vtx->z();
  }
  
  // fill_pre_jets();
  //cout<<numberjets<<" , "<<num_jets<<endl;
  // if(num_jets < 2)    {  veto_event();    ++num_2jet_veto;  }
  
  // if(num_jets >= N_JETS)
  //  cout << "[error]>\t num_jets = " << num_jets 
  //       << " which is greater than N_JETS!!!" << endl;
}

//===================================================================//

  bool TupleMaker::cut_jet(PFJet j)
{//cout<<j.et()<<endl;
 
 bool 
  jetcut = ((j.pt() > 30) && (fabs(j.eta()) < CUT_JET_ETA)),
  check1 = (j.neutralEmEnergyFraction() < CUT_NEUTR_EM_EN_FRAC),
  check2 = (j.chargedEmEnergyFraction() < CUT_CHARG_EM_EN_FRAC),
  check3 = (j.neutralHadronEnergyFraction() < CUT_NEUTR_HAD_EN_FRAC),
    check4 = (j.chargedHadronEnergyFraction() > CUT_CHARG_HAD_EN_FRAC),
  // cout<<j.chargedMultiplicity()<<" , "<<j.neutralMultiplicity()<<endl;
  // cout<<j.et()<<" , "<<j.pt()<<endl;
    check5= (j.chargedMultiplicity() >0),
    check6= ((j.chargedMultiplicity()+j.neutralMultiplicity()) >0);
    return ( jetcut && check1 && check2 && check3 && check4&&check5&&check6);
  }

 //===================================================================//
 
bool TupleMaker::noniso_jet(PFJet j)
{
  double dR = deltaR_ljet(j);
  
  bool doweveto = (dR > CUT_DELTA_R_LEPJET);
  
  // if(doweveto)   {  veto_event();    ++num_isoj_veto;  }
  
  //if(dR > 0.7)   return doweveto;    //don't fill histograms
  
  //  bool both_em = (dilep_ == DILEP_em);
  
  // if( (dilep_ == DILEP_ee) || both_em || (lepjet_ == LEPJET_ej) )
    //  fill(h_pre_dR_ej, deltaR_ejet(j, 1));
  
  //if( (dilep_ == DILEP_mm) || both_em || (lepjet_ == LEPJET_mj) )
  //  fill(h_pre_dR_mj, deltaR_mjet(j, 1));
  
  // if(dilep_ == DILEP_ee)
    //   fill(h_pre_dR_ej, deltaR_ejet(j, 2));
  
  // if(dilep_ == DILEP_mm)
    
  return doweveto;
  }

 //===================================================================//
 /*
   bool TupleMaker::cut_bjets()
   {
   int n_bjets_max2 = num_bjets;
  
   if(n_bjets_max2 > 2)      n_bjets_max2 = 2;
   
   bool doweveto = false;
   
   switch(n_bjets_max2) 
  {
    case 0:
      return (true);
      break;
    case 2:  
      doweveto = (bjet_2.pt() < CUT_2ND_BJET_PT);      
      //no break
    case 1:
      if(bjet_1.pt() < CUT_1ST_BJET_PT)   
        return (true);
      break;
    default:
      cout << "[error]>\t num_bjets = " << num_bjets << endl;
      break;
  }
  
  return doweveto;
}
*/
//===================================================================//
 /*
bool TupleMaker::lepton_jet(PFJetCollection::const_iterator ijet)
{
  PFJet j = *ijet;
  
  double dR = deltaR_ljet(j);
  
  if(dR < R40)
    fill(h_pre_dR_lj, dR);
  
  return (dR < MATCH_DELTA_R_LEPJET); 
}
 */
//===================================================================//
  
double TupleMaker::deltaR_ljet(PFJet j)
{
  double dR1 = 999., dR2 = 999.; 
  
  if(DiLeptonType==3)
    {
	  dR1 = deltaR_ejet(j, 1);
	  dR2 = deltaR_ejet(j, 2);
    }  
  if(DiLeptonType==1)
    { 
        dR1 = deltaR_mjet(j, 1);
        dR2 = deltaR_mjet(j, 2);
    }
  if(DiLeptonType==5)
    {
          dR1 = deltaR_ejet(j, 1);
	  dR2 = deltaR_mjet(j, 1);
    }
    
    if(dR1 < dR2)      return dR1;
    else               return dR2;
  
   if(process_ == PROCESS_lj)
    switch(lepjet_) 
  {
    case (LEPJET_ej):  
      //   return deltaR_ejet(j, 1);
      break;
    case (LEPJET_mj):
      return deltaR_mjet(j, 1);
      break;
    case (LEPJET_no):
    default:
      break;
  }  
  
  return (999.);
}
  
//===================================================================//
   
double TupleMaker::deltaR_ejet(PFJet j, int which_e)
{
  double etal, etaj, phil, phij;
  
  vector<GsfElectron>::const_iterator ielec;
  
  switch(which_e) 
  {
    case 1: 
 phil = elec1P4.Phi();
  etal = elec1P4.Eta();
      break;
    case 2: 
 phil = elec2P4.Phi();
  etal = elec2P4.Eta();
      break;
    default:
      cout << "[error]>\t which_e = " << which_e << endl;
      return 999.;
      break;
  }
  
  phij = j.phi();
  etaj = j.eta();
  
  return DeltaRX(etal, etaj, phil, phij);
}
   
//===================================================================//
    
double TupleMaker::deltaR_mjet(PFJet j, int which_m)
{
  double etal, etaj, phil, phij;
  
  vector<Muon>::const_iterator imuon;
  
  switch(which_m) 
  {
    case 1:      imuon = muon_1;      break;
    case 2:      imuon = muon_2;      break;
    default:
      cout << "[error]>\t which_m = " << which_m << endl;
      return 999.;
      break;
  }
  
  phil = imuon->phi();        phij = j.phi();
  etal = imuon->eta();        etaj = j.eta();
  
  return DeltaRX(etal, etaj, phil, phij);
}

//===================================================================//
     
PFJet TupleMaker::correct_jet(PFJetCollection::const_iterator org_j, 
			      const Event& event_, const EventSetup& setup_)
{
  PFJet j = *org_j;  
  
  //  double jec_pt, org_pt = j.pt();

  // j.scaleEnergy(corr_L1->correction(j, event_, setup_) );
  //j.scaleEnergy(corr_L2->correction(j, event_, setup_) );
  //j.scaleEnergy(corr_L3->correction(j, event_, setup_) );
  
  // if(IS_DATA)
    // j.scaleEnergy(corr_MC->correction(j, event_, setup_) );
  
    // jec_pt = j.pt();
  
  //  if((jec_pt > CUT_JET_PT) && (org_pt > 0.))
  //{
  //  fill(h_raw_jet_rpt, jec_pt / org_pt);
  //  fill(h_raw_jet_dpt, jec_pt - org_pt);
  // }
  
  return (j);
}

//===================================================================//
      
void TupleMaker::order_jets(int nj, PFJet j)
{  
  int i, k;
  
  PFJet temp;
  
  jet_[nj-1] = j;      //add to the end of the queue
  
  for(k = 1; k < nj; k++)
  {
    temp = jet_[k];
    
    for(i = k - 1; (i >= 0) && (jet_[i].pt() < temp.pt()); i--)
      jet_[i+1] = jet_[i];

    jet_[i+1] = temp;
  }
}
      
//===================================================================//

//===================================================================//

//===================================================================//

bool TupleMaker::cut_same_sign()
{
  if(process_ != PROCESS_di)      return (false);
  
  set_sign();
  
  if(sign_ == SIGN_ss)
  {
    if(dilep_ == DILEP_ee)        ++num_ss_ee;
    if(dilep_ == DILEP_mm)        ++num_ss_mm;
    if(dilep_ == DILEP_em)        ++num_ss_em;
  }
  
  if(SKIP_SS_CUT)                 return (false);
  
  return (sign_ != SIGN_os);
}

//===================================================================//

void TupleMaker::set_sign()
{
  if(process_ != PROCESS_di)      return;
  
  int q1 = 0, q2 = 0;
  
  switch(dilep_)
  {
    case DILEP_ee:    
      q1 = elec_1->charge();      q2 = elec_2->charge();    
      break;
    case DILEP_mm:    
      q1 = muon_1->charge();      q2 = muon_2->charge();    
      break;
    case DILEP_em:    
      q1 = elec_1->charge();      q2 = muon_1->charge();    
      break;
    default:                                              
      break;
  }
  
  switch(q1*q2)
  {
    case -1:    sign_ = SIGN_os;    break;
    case  1:    sign_ = SIGN_ss;    break;
    case  0:    default:
      cout << "[error]>\t q1, q2 = " << q1 << ", " << q2 << endl;
      break;
  }
}

//===================================================================//

//===================================================================//
 double TupleMaker::leptonVertex( reco::TrackRef lep1track, reco::TrackRef lep2track, const EventSetup& setup_)
{

             ESHandle<MagneticField> bField;
             setup_.get<IdealMagneticFieldRecord>().get(bField);
             reco::TransientTrack mu1T(lep1track,bField.product());
             reco::TransientTrack mu2T(lep2track,bField.product());
             reco::Vertex vtx;
             vector<reco::TransientTrack > dmT;
             dmT.push_back(mu1T);
             dmT.push_back(mu2T);
             KalmanVertexFitter KalmanFitter;
             TransientVertex vertex;
             bool isVertex = true;
             try { vertex = KalmanFitter.vertex(dmT); }
             catch ( exception & err) { isVertex = false; }
             if (isVertex){
                 vtx = vertex;
                 float vtxProb1 = (float)TMath::Prob(vtx.chi2(),int(vtx.ndof()));
		 // if (vtxProb1 > vtxProb) vtxProb = vtxProb1;
		 return vtxProb1;
	     }
		 else
		 return 0;

}


//===================================================================//
 double TupleMaker::leptonVertex( reco::GsfTrackRef lep1track, reco::GsfTrackRef lep2track, const EventSetup& setup_)
{
  edm::ESHandle<TransientTrackBuilder> theB;
  setup_.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
  
             ESHandle<MagneticField> bField;
             setup_.get<IdealMagneticFieldRecord>().get(bField);
 
      reco::TransientTrack mu1T = (*theB).build(*lep1track);
      reco::TransientTrack mu2T = (*theB).build(*lep2track);
      //            reco::TransientTrack mu1T(lep1track,bField.product());
      //             reco::TransientTrack mu2T(lep2track,bField.product());
             reco::Vertex vtx;
             vector<reco::TransientTrack > dmT;
             dmT.push_back(mu1T);
             dmT.push_back(mu2T);
             KalmanVertexFitter KalmanFitter;
             TransientVertex vertex;
             bool isVertex = true;
             try { vertex = KalmanFitter.vertex(dmT); }
             catch ( exception & err) { isVertex = false; }
             if (isVertex){
                 vtx = vertex;
                 float vtxProb1 = (float)TMath::Prob(vtx.chi2(),int(vtx.ndof()));
		 // if (vtxProb1 > vtxProb) vtxProb = vtxProb1;
		 return vtxProb1;
	     }
		 else
		 return 0;

}


//===================================================================//
 double TupleMaker::leptonVertex( reco::TrackRef lep1track, reco::GsfTrackRef lep2track, const EventSetup& setup_)
{
 edm::ESHandle<TransientTrackBuilder> theB;
  setup_.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
  //cout<<"1"<<endl;
             ESHandle<MagneticField> bField;
             setup_.get<IdealMagneticFieldRecord>().get(bField);
             reco::TransientTrack mu1T = (*theB).build(*lep1track);
	     // reco::TransientTrack mu2T(lep2track,bField.product());
      reco::TransientTrack mu2T = (*theB).build(*lep2track);
             reco::Vertex vtx;
	     //	     cout<<"hello"<<endl;
             vector<reco::TransientTrack > dmT;
             dmT.push_back(mu1T);
             dmT.push_back(mu2T);
             KalmanVertexFitter KalmanFitter;
             TransientVertex vertex;
             bool isVertex = true;
	     //	     cout<<"still good"<<endl;
             try { vertex = KalmanFitter.vertex(dmT); }
             catch ( exception & err) { isVertex = false; }
             if (isVertex){
                 vtx = vertex;
                 float vtxProb1 = (float)TMath::Prob(vtx.chi2(),int(vtx.ndof()));
		 // if (vtxProb1 > vtxProb) vtxProb = vtxProb1;
		 return vtxProb1;
	     }
		 else
		 return 0;

}


//===================================================================//

//===================================================================//

bool TupleMaker::trainTrigPresel(const reco::GsfElectron& ele) {
  
  bool myTrigPresel = false;
  if(fabs(ele.superCluster()->eta()) < 1.479) {
    if(ele.sigmaIetaIeta() < 0.014 &&
       ele.hadronicOverEm() < 0.15 &&
       ele.dr03TkSumPt()/ele.pt() < 0.2 &&
       ele.dr03EcalRecHitSumEt()/ele.pt() < 0.2 &&
       ele.dr03HcalTowerSumEt()/ele.pt() < 0.2 &&
       ele.gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() == 0)
      myTrigPresel = true;
  }
  else {
    if(ele.sigmaIetaIeta() < 0.035 &&
       ele.hadronicOverEm() < 0.10 &&
       ele.dr03TkSumPt()/ele.pt() < 0.2 &&
       ele.dr03EcalRecHitSumEt()/ele.pt() < 0.2 &&
       ele.dr03HcalTowerSumEt()/ele.pt() < 0.2 &&
       ele.gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() == 0)
      myTrigPresel = true;
  }
  
  
  return myTrigPresel;
}
 //=================================================================//


//==================================================================================//

bool TupleMaker::triggerDecision(edm::Handle<edm::TriggerResults> &hltResults, int iTrigger)
{
  bool triggerPassed = false;
  if(hltResults->wasrun(iTrigger) &&
      hltResults->accept(iTrigger) &&
      !hltResults->error(iTrigger) ){
    triggerPassed = true;
  }
  return triggerPassed;
}


//===================================================================//

void TupleMaker::set_h_bin(TH1D* h, int ibin, int value)
{
  set_h_bin(h, ibin, (double)(value));
}

//===================================================================//

void TupleMaker::set_h_bin(TH1D* h, int ibin, double value)
{
  h->SetBinContent(ibin, value);
}

//===================================================================//


//===================================================================//

void TupleMaker::save_counters()
{
  int BIN_ONE = 1;
  
  set_h_bin(h_num_events, BIN_ONE, num_events);
  
}

//===================================================================//

void TupleMaker::endJob() 
{
  save_counters();
  
  cout << endl << "[*end*]>\t " << "# " << num_events << " ";
  //tree->Write();
  watch->Stop();  
  watch->Print();
  
  cout << endl << LINE << endl;
}

//===================================================================//

void TupleMaker::beginRun(Run const&, EventSetup const&) { }

//===================================================================//

void TupleMaker::endRun(Run const&, EventSetup const&) { }

//===================================================================//

void TupleMaker::beginLuminosityBlock(LuminosityBlock const&, EventSetup const&)
{
  
}

//===================================================================//

void TupleMaker::endLuminosityBlock(LuminosityBlock const&, EventSetup const&)
{
  
}

//===================================================================//

void TupleMaker::fillDescriptions(ConfigurationDescriptions& descriptions) 
{
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//===================================================================//

//void Spots::print_debug(const char* strng, int n_check, int var)
//{
//  print_debug(strng, n_check, (double)(var));  
//}

//===================================================================//

void TupleMaker::print_debug(const char* strng, int n_check, double var)
{
  if(num_events > n_check)
    return;
  //if(abs(num_events - 11281) >= 25)  
    //return;
  
  cout << "[debug]>\t # " << num_events << " @ " << strng;
  
  if(var == JUNK)       cout << endl;
  else                  cout << " w/ x = " << var << endl;
}

//===================================================================//

TupleMaker::~TupleMaker()  { }

//===================================================================//

//define this as a plug-in
DEFINE_FWK_MODULE(TupleMaker);

//===================================================================//
