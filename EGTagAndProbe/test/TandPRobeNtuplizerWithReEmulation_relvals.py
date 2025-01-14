# RUN THE RE-EMULATION Workflow for MC
# Uses the RECO + RAW Setup

import FWCore.ParameterSet.VarParsing as VarParsing
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
from Configuration.StandardSequences.Eras import eras
from Configuration.Eras.Era_Run3_cff import Run3

isMC = True
isMINIAOD = False

process = cms.Process("TagAndProbe",eras.Run3)
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')


if not isMC: # will use 80X
    from Configuration.AlCa.autoCond import autoCond
    process.GlobalTag.globaltag = ''
    process.load('EGTagAndProbe.EGTagAndProbe.tagAndProbe_cff')
    process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
        'file:/eos/home-a/athachay/workarea/120e6d8c-ec45-40ed-9f5c-986a02b2492b.root'
        )
    )
else:
    print("Processing MC")
    process.GlobalTag.globaltag = '142X_mcRun3_2025_realistic_v4'
    process.load('EGTagAndProbe.EGTagAndProbe.MCanalysis_cff')
    process.source = cms.Source("PoolSource",
     fileNames= cms.untracked.vstring(
    '/store/relval/CMSSW_14_2_0/RelValZEE_14/GEN-SIM-RECO/PU_142X_mcRun3_2025_realistic_v4_Winter25_PU_RV255_6-v1/2590000/00b21376-7244-46a2-9102-305e02d9f3fd.root'
     ),
    secondaryFileNames = cms.untracked.vstring(
        '/store/relval/CMSSW_14_2_0/RelValZEE_14/GEN-SIM-DIGI-RAW/PU_142X_mcRun3_2025_realistic_v4_Winter25_PU_RV255_6-v1/2590000/1167b234-60c4-4b2b-84ba-30ca83db8236.root',
        '/store/relval/CMSSW_14_2_0/RelValZEE_14/GEN-SIM-DIGI-RAW/PU_142X_mcRun3_2025_realistic_v4_Winter25_PU_RV255_6-v1/2590000/611d47c1-e3f3-4bdb-9574-e52b7d066028.root',
        '/store/relval/CMSSW_14_2_0/RelValZEE_14/GEN-SIM-DIGI-RAW/PU_142X_mcRun3_2025_realistic_v4_Winter25_PU_RV255_6-v1/2590000/6f47f4be-33c0-4ac3-9882-11e96e0b4e41.root',
        '/store/relval/CMSSW_14_2_0/RelValZEE_14/GEN-SIM-DIGI-RAW/PU_142X_mcRun3_2025_realistic_v4_Winter25_PU_RV255_6-v1/2590000/b3d121f9-0ff6-48e5-b351-c3d421577963.root',
        '/store/relval/CMSSW_14_2_0/RelValZEE_14/GEN-SIM-DIGI-RAW/PU_142X_mcRun3_2025_realistic_v4_Winter25_PU_RV255_6-v1/2590000/d2dac9d3-a2f6-4b6b-a669-ef785cedc2f3.root',
        '/store/relval/CMSSW_14_2_0/RelValZEE_14/GEN-SIM-DIGI-RAW/PU_142X_mcRun3_2025_realistic_v4_Winter25_PU_RV255_6-v1/2590000/d6387201-2618-4d52-909b-19fc04b4bbb0.root',
        '/store/relval/CMSSW_14_2_0/RelValZEE_14/GEN-SIM-DIGI-RAW/PU_142X_mcRun3_2025_realistic_v4_Winter25_PU_RV255_6-v1/2590000/d9e7dc57-6b66-4d8a-8101-f70bc48130fc.root',
        '/store/relval/CMSSW_14_2_0/RelValZEE_14/GEN-SIM-DIGI-RAW/PU_142X_mcRun3_2025_realistic_v4_Winter25_PU_RV255_6-v1/2590000/f6397966-d225-41fc-bf7e-7758ce55ab28.root',
        '/store/relval/CMSSW_14_2_0/RelValZEE_14/GEN-SIM-DIGI-RAW/PU_142X_mcRun3_2025_realistic_v4_Winter25_PU_RV255_6-v1/2590000/f74c1f2f-e45f-48ad-952b-500f1a463244.root',
        '/store/relval/CMSSW_14_2_0/RelValZEE_14/GEN-SIM-DIGI-RAW/PU_142X_mcRun3_2025_realistic_v4_Winter25_PU_RV255_6-v1/2590001/41cf944a-be68-4b87-8fa3-1319d9c62db9.root'
    )
   )
    
    process.Ntuplizer.useHLTMatch = cms.bool(False) #In case no HLT object in MC sample considered or you're fed up with trying to find the right HLT collections


if isMINIAOD:
    process.egmGsfElectronIDSequence = cms.Sequence()
    process.Ntuplizer.electrons = cms.InputTag("slimmedElectrons")
    process.Ntuplizer.genParticles = cms.InputTag("prunedGenParticles")
    process.Ntuplizer.Vertices = cms.InputTag("offlineSlimmedPrimaryVertices")
else :
    from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
    dataFormat = DataFormat.AOD
    process.load("RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi")
    switchOnVIDElectronIdProducer(process, dataFormat)
    my_id_modules = [
        'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_RunIIIWinter22_iso_V1_cff',
        'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_RunIIIWinter22_noIso_V1_cff',
        'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Winter22_122X_V1_cff'
    ]
    
    process.electronMVAValueMapProducer.src=cms.InputTag("gedGsfElectrons")

    for idmod in my_id_modules:
        setupAllVIDIdsInModule(process, idmod, setupVIDElectronSelection)
    process.Ntuplizer.eleLooseIdMap  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-RunIIIWinter22-V1-loose")
    process.Ntuplizer.eleTightIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-RunIIIWinter22-iso-V1-wp80")
    process.Ntuplizer.eleMediumIdMap  = cms.InputTag("egmGsfElectronIDs:mvaEleID-RunIIIWinter22-iso-V1-wp90")

process.schedule = cms.Schedule()

## L1 emulation stuff

if not isMC:
    from L1Trigger.Configuration.customiseReEmul import L1TReEmulFromRAW 
    process = L1TReEmulFromRAW(process)
else:
    from L1Trigger.Configuration.customiseReEmul import L1TReEmulMCFromRAW
    process = L1TReEmulMCFromRAW(process) 
    from L1Trigger.Configuration.customiseUtils import L1TTurnOffUnpackStage2GtGmtAndCalo 
    process = L1TTurnOffUnpackStage2GtGmtAndCalo(process)


process.load("L1Trigger.L1TCalorimeter.caloParams_2024_v0_2_cfi")

#### handling of cms line options for tier3 submission
#### the following are dummy defaults, so that one can normally use the config changing file list by hand etc.

if options.JSONfile:
    print( "Using JSON: " , options.JSONfile)
    process.source.lumisToProcess = LumiList.LumiList(filename = options.JSONfile).getVLuminosityBlockRange()

if options.inputFiles:
    process.source.fileNames = cms.untracked.vstring(options.inputFiles)

if options.secondaryFilesList:
    listSecondaryFiles = FileUtils.loadListFromFile(options.secondaryFilesList)
    process.source.secondaryFileNames = cms.untracked.vstring(listSecondaryFiles)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

if options.maxEvents >= -1:
    process.maxEvents.input = cms.untracked.int32(options.maxEvents)
if options.skipEvents >= 0:
    process.source.skipEvents = cms.untracked.uint32(options.skipEvents)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.p = cms.Path (
    process.egmGsfElectronIDSequence +
    process.RawToDigi +
    process.L1TReEmul +
    process.NtupleSeq
)


process.schedule = cms.Schedule(process.p) # do my sequence pls

# Silence output
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

# Adding ntuplizer
process.TFileService=cms.Service('TFileService',fileName=cms.string(options.outputFile))
