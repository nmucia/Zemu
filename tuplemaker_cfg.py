import FWCore.ParameterSet.Config as cms
import os
process = cms.Process("Demo")


# real data or MC?
isRealData = True
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.MessageLogger.cerr.threshold = 'ERROR'    #'INFO'
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(-1)  #via Nate
process.MessageLogger.categories.append('Demo')
process.MessageLogger.cerr.INFO = cms.untracked.PSet( limit = cms.untracked.int32(-1) )
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(3000) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(250) )


process.MessageLogger = cms.Service("MessageLogger",
                                    cout = cms.untracked.PSet(
                                                              default = cms.untracked.PSet( limit = cms.untracked.int32(0)),
                                                              FwkJob = cms.untracked.PSet( limit = cms.untracked.int32(0) )
                                                              ) ,
                                    categories = cms.untracked.vstring('FwkJob'),
                                    destinations = cms.untracked.vstring('cout')
                                    )

debugging = False
#debugging = True
if debugging:
  base = os.path.relpath(os.environ.get('CMSSW_BASE'))+'/src'
#process.MessageLogger.cerr.FwkReport.reportEvery = 10
else:
  base = 'src'
#process.MessageLogger.cerr.FwkReport.reportEvery = 1000


### via Nate ###
process.load("Configuration.StandardSequences.Services_cff")
process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
# load recommended met filters
process.load("RecoMET.METFilters.metFilters_cff")

process.load("JetMETCorrections.Type1MET.correctionTermsPfMetType1Type2_cff")
process.load("JetMETCorrections.Type1MET.correctionTermsPfMetType0PFCandidate_cff")
process.load("JetMETCorrections.Type1MET.correctionTermsPfMetShiftXY_cff")
process.load("JetMETCorrections.Type1MET.correctedMet_cff")
#process.printTree = cms.EDAnalyzer("ParticleListDrawer",
#  maxEventsToPrint = cms.untracked.int32(1),
#  printVertex = cms.untracked.bool(False),
#  src = cms.InputTag("genParticles")
#)
### global tag
if (isRealData):
#    process.GlobalTag.globaltag = 'GR_R_44_V14::All'
#    process.GlobalTag.globaltag = 'FT_53_V21_AN3::All'
  process.GlobalTag.globaltag = 'FT53_V21A_AN6::All'
 #    process.GlobalTag.globaltag = 'FT_53_V21_AN4::ALL'
  process.corrPfMetType1.jetCorrLabel = cms.string("ak5PFL1FastL2L3Residual")
  process.corrPfMetShiftXY.parameter = process.pfMEtSysShiftCorrParameters_2012runABCDvsNvtx_data
else:
       process.GlobalTag.globaltag = 'START53_V27::All'
       process.corrPfMetType1.jetCorrLabel = cms.string("ak5PFL1FastL2L3")
       process.corrPfMetShiftXY.parameter = process.pfMEtSysShiftCorrParameters_2012runABCDvsNvtx_mc
#    process.GlobalTag.globaltag = 'START44_V13::All'


process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.load("CondCore.DBCommon.CondDBCommon_cfi")
##Jet Energy Corrections
process.load("RecoJets.Configuration.RecoPFJets_cff")
process.load("RecoJets.Configuration.RecoJPTJets_cff")
process.load("RecoJets.JetAssociationProducers.ak5JTA_cff")
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.load('JetMETCorrections.Configuration.JetCorrectionServices_cff')
process.kt6PFJets.doRhoFastjet = True

#process.ak5JPTL1Offset.algorithm = 'AK5JPT'
process.ak5JetTracksAssociatorAtVertex.useAssigned = cms.bool(True)
#process.ak5JetTracksAssociatorAtVertex.pvSrc = cms.InputTag("offlinePrimaryVertices")


## Ecal Dead Cell Filter
#process.load('RecoMET.METFilters.EcalDeadCellTriggerPrimitiveFilter_cfi')

process.jpt = cms.Sequence(
                       #process.primaryVertexFilter *
                       # process.ak5JTA*process.recoJPTJets *
                       # process.ak5JPTJetsL1L2L3 *
  process.kt6PFJets *
  process.ak5PFJetsL1FastL2L3 
  )
process.load('EgammaAnalysis/ElectronTools/electronIdMVAProducer_cfi')

#isolation
from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso, setupPFMuonIso, setupPFPhotonIso
process.eleIsoSequence = setupPFElectronIso(process, 'gsfElectrons')
process.muIsoSequence = setupPFMuonIso(process, 'muons')
process.phoIsoSequence = setupPFPhotonIso(process, 'photons')


# Electron Regression (post moriond recommendation)
process.load('EgammaAnalysis/ElectronTools/electronRegressionEnergyProducer_cfi')
process.eleRegressionEnergy.inputElectronsTag = cms.InputTag('gsfElectrons')
process.eleRegressionEnergy.inputCollectionType = cms.uint32(0)
process.eleRegressionEnergy.useRecHitCollections = cms.bool(True)
process.eleRegressionEnergy.produceValueMaps = cms.bool(True)
process.eleRegressionEnergy.energyRegressionType = cms.uint32(2)
process.eleRegressionEnergy.regressionInputFile = cms.string("EgammaAnalysis/ElectronTools/data/eleEnergyRegWeights_WithSubClusters_VApr15.root")

# Electron Combination for calibration (post moriond recommendation)
process.load('EgammaAnalysis/ElectronTools/calibratedElectrons_cfi')
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    calibratedElectrons = cms.PSet(
      initialSeed = cms.untracked.uint32(1),
      engineName = cms.untracked.string('TRandom3')
    ),
)
if (isRealData):
  process.calibratedElectrons.isMC = cms.bool(False)
  process.calibratedElectrons.inputDataset = cms.string("22Jan2013ReReco")
else:
  process.calibratedElectrons.isMC = cms.bool(True)
  process.calibratedElectrons.inputDataset = cms.string("Summer12_LegacyPaper")

process.calibratedElectrons.updateEnergyError = cms.bool(True)
process.calibratedElectrons.correctionsType = cms.int32(2)
process.calibratedElectrons.combinationType = cms.int32(3)
process.calibratedElectrons.lumiRatio = cms.double(1.0)
process.calibratedElectrons.verbose = cms.bool(False)
process.calibratedElectrons.synchronization = cms.bool(False)
process.calibratedElectrons.applyLinearityCorrection = cms.bool(True)


JetCorrection = "ak5PFL1FastL2L3"
if (isRealData):
  JetCorrection += "Residual"
process.ak5PFJetsCorr = cms.EDProducer('PFJetCorrectionProducer',
                                        src = cms.InputTag("ak5PFJets"),
                                        correctors = cms.vstring(JetCorrection) # NOTE: use "ak5PFL1FastL2L3" for MC / "ak5PFL1FastL2L3Residual" for Data
                                        )

# load PUJetID

process.load("RecoJets.JetProducers.PileupJetID_cfi")
process.pileupJetIdProducer.jets = "ak5PFJetsCorr"
process.pileupJetIdProducer.residualsTxt = cms.FileInPath("RecoJets/JetProducers/data/mva_JetID_v1.weights.xml")

process.MyAk5PFJetTracksAssociatorAtVertex = cms.EDProducer("JetTracksAssociatorAtVertex",
                                                            process.j2tParametersVX,
                                                            jets = cms.InputTag("ak5PFJetsCorr")
                                                            )
process.MyImpactParameterPFTagInfos = process.impactParameterTagInfos.clone(
  jetTracks = "MyAk5PFJetTracksAssociatorAtVertex"
  )
process.MetSequence = cms.Sequence(process.correctionTermsPfMetType1Type2
                                   * process.correctionTermsPfMetType0PFCandidate
                                   * process.correctionTermsPfMetShiftXY
                                   * process.pfMetT0pcT1Txy
                                                       )

process.JetSequence = cms.Sequence(process.ak5PFJetsCorr
                                   * process.pileupJetIdProducer
                                   * process.MyAk5PFJetTracksAssociatorAtVertex
                                   * process.MyImpactParameterPFTagInfos
                                   #* process.MySecondaryVertexTagInfos
                                   #* process.MyCombinedSecondaryVertexBJetTags
                                   )

if not(isRealData):
 from RecoJets.Configuration.GenJetParticles_cff import genParticlesForJetsNoNu
 process.genParticlesForJetsNoNu = genParticlesForJetsNoNu
 from RecoJets.JetProducers.ak5GenJets_cfi import ak5GenJets
 process.ak5GenJetsNoNu = ak5GenJets.clone( src = cms.InputTag("genParticlesForJetsNoNu") )
 process.JetSequence += process.genParticlesForJetsNoNu
 process.JetSequence += process.ak5GenJetsNoNu


# JEC uncertainty files

#process.ak5PFJets.doAreaFastjet = True


### To get b-tags from ak5PFJets - this works with the 'Reconstruction_cff' above

#process.path = cms.Path(process.btagging)              #djs - for full reconstruction

#this part is from nate
#djs commented out as a test - the above works and should be used when done testing

#djs via #https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookBTagging#DifferentJets
#process.newJetBProbabilityBJetTags = process.jetBProbabilityBJetTags.clone()
#process.newJetBProbabilityBJetTags.tagInfos = cms.VInputTag( cms.InputTag("newImpactParameterTagInfos") )
#djs - keep TCHE and JP for now; note order may be important

#MET corrections --do not work yet
##____________________________________________________________________________||
#process.load('JetMETCorrections.Type1MET.pfMETCorrections_cff')

##____________________________________________________________________________||
#process.load("JetMETCorrections.Type1MET.correctionTermsPfMetType1Type2_cff")

#process.corrPfMetType1.jetCorrLabel = cms.string("ak5PFL1FastL2L3")
# process.corrPfMetType1.jetCorrLabel = cms.string("ak5PFL1FastL2L3Residual")

##____________________________________________________________________________||
#process.load("JetMETCorrections.Type1MET.correctionTermsPfMetType0PFCandidate_cff")

##____________________________________________________________________________||
#process.load("JetMETCorrections.Type1MET.correctionTermsPfMetType0RecoTrack_cff")

##____________________________________________________________________________||
#process.load("JetMETCorrections.Type1MET.correctionTermsPfMetShiftXY_cff")

#process.corrPfMetShiftXY.parameter = process.pfMEtSysShiftCorrParameters_2012runABCDvsNvtx_mc
# process.corrPfMetShiftXY.parameter = process.pfMEtSysShiftCorrParameters_2012runABCDvsNvtx_data

##____________________________________________________________________________||
#process.load("JetMETCorrections.Type1MET.correctedMet_cff")

### MET corrections
#from JetMETCorrections.Type1MET.MetType1Corrections_cff import metJESCorAK5PFJet
#process.metJESCorAK5PF = metJESCorAK5PFJet.clone()
#process.metJESCorAK5PF.inputUncorJetsLabel = "ak5PFJets"
#process.metJESCorAK5PF.metType = "PFMET"
#process.metJESCorAK5PF.inputUncorMetLabel = "pfMet"
#process.metJESCorAK5PF.useTypeII = False
#process.metJESCorAK5PF.jetPTthreshold = cms.double(10.0)
#process.metJESCorAK5PF.corrector = cms.string('ak5PFL1FastL2L3')
### via Nate ###

#electron iso




process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    #'/store/mc/Fall11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0002/A0226FB3-F636-E111-8035-002354EF3BE2.root',
    #'/store/mc/Fall11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0002/A04710B0-ED36-E111-A509-001A92971BB4.root',
    #'/store/mc/Fall11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0002/A0683FE7-F736-E111-BF75-003048679084.root',
    #'/store/mc/Fall11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0002/A070D4DA-F836-E111-8D32-0026189438C2.root',
    #'/store/mc/Fall11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0002/A0731EB6-F836-E111-9F1F-002618943856.root',
    #'/store/mc/Fall11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0002/A07C073E-F736-E111-B296-002618FDA207.root',
    #'/store/mc/Fall11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0002/A085A938-F736-E111-9B9B-001A92971BC8.root',
    #'/store/mc/Fall11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0002/A0967AC6-EB36-E111-BA15-003048FFD740.root',
    #'/store/mc/Fall11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0002/A0981EB2-EB36-E111-9042-002618943836.root',
    #'/store/mc/Fall11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0002/A0A30211-E436-E111-B2F5-001A9281172C.root',
    #'/store/mc/Fall11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0002/A0C9221E-E636-E111-B2B0-003048678E8A.root',
	#'/store/mc/Fall11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0002/A217D591-F736-E111-A274-002618943958.root',
	#'/store/mc/Fall11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0002/A21E3BD3-F836-E111-B140-0018F3D0960A.root',
	#'/store/mc/Fall11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0002/A268397F-EB36-E111-B4DF-00261894397D.root',
	#'/store/mc/Fall11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0002/A282571C-E436-E111-B910-00304867915A.root',
	#'/store/mc/Fall11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0002/A29D5180-E836-E111-ABA1-003048D15DDA.root',
	#'/store/mc/Fall11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0002/A2C29066-FD36-E111-8546-0026189438A5.root',
	#'/store/mc/Fall11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0002/A2CDCE96-ED36-E111-B980-0026189438BC.root',
	#'/store/mc/Fall11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0002/A2DC50C9-ED36-E111-A3B5-0026189438A9.root',
	#'/store/mc/Fall11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0002/A42827D1-F836-E111-A80A-003048FFD752.root',
	#'/store/mc/Fall11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0002/A45478D2-F736-E111-9474-002618943862.root',
	#'/store/mc/Fall11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0002/A4657910-F836-E111-99B2-0018F3D09630.root',
	#'/store/mc/Fall11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0002/A47050A2-ED36-E111-809C-0026189437F2.root',
	#'/store/mc/Fall11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0002/A479E586-E836-E111-ADBE-001A92971BA0.root',
	#'/store/mc/Fall11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0002/A47B95C9-F836-E111-B22A-00261894394D.root',
    #'/store/mc/Fall11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0002/A4C42D2A-F836-E111-83A6-003048678B7C.root',
	#'/store/mc/Fall11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0002/A4C448E6-F836-E111-9CCD-001A92810ADE.root',
	#'/store/mc/Fall11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0002/A4D7D319-E636-E111-8DA7-00248C0BE018.root',
	#'/store/mc/Fall11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0002/A4F20B92-E836-E111-B216-0018F3D096EC.root',
	#'/store/mc/Fall11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0002/A637A303-F836-E111-8C6B-003048D3C010.root'
	#'/store/data/Run2011A/DoubleMu/AOD/21Jun2013-v1/10000/00255071-06DE-E211-8EC3-00261894392D.root'
	#'/store/mc/Fall11/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia/AODSIM/PU_S6_START44_V9B-v1/0000/0013E75A-BB2A-E111-B78E-003048C8EE58.root'
	#'/store/data/Run2011A/DoubleMu/AOD/May10ReReco-v1/0005/E45CB490-F77C-E011-BDEB-002618943959.root'
	#'/store/data/Run2011A/DoubleMu/AOD/05Aug2011-v1/0000/002948D6-D8C0-E011-82C3-002618943862.root'
	#'/store/data/Run2011A/DoubleMu/AOD/08Nov2011-v1/0001/FA913DBA-371B-E111-A78A-00261894380A.root'
	#'/store/mc/Fall11/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/AODSIM/PU_S6_START42_V14B-v1/0000/0007A683-23F6-E011-9703-90E6BA0D09DC.root'
	#'/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START50_V15-v1/0000/003AD0ED-C473-E111-B1F3-F04DA23BCE4C.root'
	#'/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/00037C53-AAD1-E111-B1BE-003048D45F38.root'
	#'root://xrootd.unl.edu//store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph/AODSIM/PU_S10_START53_V7A-v1/0000/001318CD-C5F4-E111-AAD9-001E67398D72.root'
	#'root://xrootd.unl.edu//store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/00037C53-AAD1-E111-B1BE-003048D45F38.root'
	#'root://xrootd.unl.edu//store/data/Run2012A/DoubleMuParked/AOD/22Jan2013-v1/20000/002557DC-16D4-E211-9415-20CF3027A5AB.root'
	#'root://xrootd.unl.edu//store/data/Run2012A/DoubleMuParked/AOD/22Jan2013-v1/20000/006D1F0E-DAD3-E211-B998-20CF3019DEFB.root'
	# '/store/data/Run2012A/DoubleMuParked/AOD/22Jan2013-v1/20000/00A96C5F-0ED4-E211-A18C-485B39800C15.root'
	#'/store/mc/Fall11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0000/001CC6EA-CD36-E111-B13A-0026189438C9.root'
	#'/store/data/Run2011B/DoubleMu/AOD/19Nov2011-v1/0000/0062A788-141C-E111-95B0-003048678BC6.root'
	#'/store/mc/Fall11/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START42_V14B-v1/0000/000D1C75-14F2-E011-A53D-003048678B1A.root'
	#'/store/mc/Fall11/QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6/AODSIM/PU_S6_START44_V9B-v1/0000/00294366-723C-E111-A326-002618943898.root'
	#'/store/mc/Fall11/QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6/AODSIM/PU_S6_START42_V14B-v1/0000/001DB9D2-0325-E111-A955-0018F3D096A6.root'
	#'/store/mc/Summer12/DYJetsToLL_M-10To50filter_8TeV-madgraph/AODSIM/PU_S7_START52_V9-v1/0000/0056FC76-2EA2-E111-B0AC-485B39800B8A.root'
	#'/store/mc/Summer12_DR53X/QCD_Pt_80_170_EMEnriched_TuneZ2star_8TeV_pythia6/AODSIM/PU_S10_START53_V7A-v1/0000/00000C25-E9DE-E111-815E-003048FFD71A.root'
	#'/store/mc/Summer12_DR53X/QCD_Pt_30_80_EMEnriched_TuneZ2star_8TeV_pythia6/AODSIM/PU_S10_START53_V7A-v1/0000/0034FD67-18E8-E111-9E60-0026189438AF.root'
	#'/store/mc/Summer12_DR53X/QCD_Pt_20_MuEnrichedPt_15_TuneZ2star_8TeV_pythia6/AODSIM/PU_S10_START53_V7A-v3/00000/006C5513-8F03-E211-A969-003048D3FC94.root'
	'root://xrootd.unl.edu//store/data/Run2012C/DoubleMuParked/AOD/22Jan2013-v1/10000/0002ACB4-C96C-E211-A96F-20CF3027A628.root'
        #'root://xrootd.unl.edu///store/data/Run2012B/DoubleMuParked/AOD/22Jan2013-v1/10000/1EC938EF-ABEC-E211-94E0-90E6BA442F24.root'
	#'/store/data/Run2012C/DoubleMuParked/AOD/22Jan2013-v1/10000/0002ACB4-C96C-E211-A96F-20CF3027A628.root'
	#'/store/data/Run2012C/SingleElectron/AOD/22Jan2013-v1/10000/00626847-9EAC-E211-ABAA-00259059642E.root'
       #'root://xrootd.unl.edu//store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/00050BBE-D5D2-E111-BB65-001E67398534.root'
       #'root://xrootd.unl.edu//store/mc/Summer12_DR53X/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v2/0000/00144391-81EE-E111-95DA-003048D479F2.root'
       # 'root://xrootd.unl.edu//store/mc/Summer12_DR53X/DYToMuMu_M-20_CT10_TuneZ2star_v2_8TeV-powheg-pythia6/AODSIM/PU_S10_START53_V7A-v1/00000/005307E8-4412-E211-98E2-20CF3027A577.root'
       #'file:/eos/uscms/store/user/nmucia2/SIM6/SIM_100_1_cZr.root'
   )
)



process.demo = cms.EDAnalyzer("TupleMaker",
                                     Electrons = cms.InputTag('calibratedElectrons','calibratedGsfElectrons',''),
                                     Photons = cms.InputTag('photons'),
                                     PFCandidateMap = cms.InputTag('particleFlow:electrons'),
                                     PrintElectrons = cms.bool(True),
                                     PrintPhotons = cms.bool(True),
#                                     IsoDepElectron = cms.VInputTag(cms.InputTag('elPFIsoDepositChargedPFIso'),
#                                                                    cms.InputTag('elPFIsoDepositGammaPFIso'),
#                                                                    cms.InputTag('elPFIsoDepositNeutralPFIso')),
#                                     IsoValElectronPF = cms.VInputTag(cms.InputTag('elPFIsoValueCharged04PFIdPFIso'),
#                                                                     cms.InputTag('elPFIsoValueGamma04PFIdPFIso'),
#                                                                     cms.InputTag('elPFIsoValueNeutral04PFIdPFIso')),
#                                     IsoDepPhoton = cms.VInputTag(cms.InputTag('phPFIsoDepositChargedPFIso'),
#                                                                  cms.InputTag('phPFIsoDepositGammaPFIso'),
#                                                                  cms.InputTag('phPFIsoDepositNeutralPFIso')),
#                                     IsoValPhoton = cms.VInputTag(cms.InputTag('phPFIsoValueCharged03PFIdPFIso'),
#                                                                  cms.InputTag('phPFIsoValueGamma03PFIdPFIso'),
#                                                                  cms.InputTag('phPFIsoValueNeutral03PFIdPFIso'))

)
process.demo.JECuncData = cms.untracked.string(base+'/data/Summer13_V5_DATA_Uncertainty_AK5PF.txt')
process.demo.JECuncMC = cms.untracked.string(base+'/data/Summer13_V5_MC_Uncertainty_AK5PF.txt')




#process.demo = cms.EDAnalyzer('TupleMaker')

JetCorrectionService = cms.string('ak5PFL1FastL2L3')

process.TFileService = cms.Service("TFileService", fileName = cms.string('histos.root') )

#process.p = cms.Path(process.myPartons * process.AK5Flavour * process.ak5PFJetsL1FastL2L3 * process.demo)

process.p =cms.Path(process.jpt * process.mvaTrigV0
                    * process.metFilters
                    * process.JetSequence
                    * process.MetSequence
                    * process.eleRegressionEnergy
                    * process.calibratedElectrons
                    * process.mvaNonTrigV0*(
  process.pfParticleSelectionSequence + 
  process.eleIsoSequence + 
  process.muIsoSequence+
  process.phoIsoSequence+
  process.demo)
                    )

#process.p = cms.Path(process.jpt * process.pfMetT1 * process.demo
#                     )
