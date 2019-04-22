import os, sys, re
import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList
from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.AlCa.GlobalTag import GlobalTag

## SET OPTIONS

options = VarParsing('analysis')
options.register ('runOnData', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "runOnData (True/False)")
options.parseArguments()

process = cms.Process("SKIM")

# Load the standard set of configuration modules
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("Geometry.CaloEventSetup.CaloTowerConstituents_cfi")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(500))
process.options  = cms.untracked.PSet(wantSummary = cms.untracked.bool(True), allowUnscheduled = cms.untracked.bool(True))

if (options.runOnData):  process.GlobalTag.globaltag = '94X_dataRun2_v11' 
else:                    process.GlobalTag.globaltag = '94X_mc2017_realistic_v17' 

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(
        'file:///mnt/hadoop/store/mc/RunIIFall17MiniAODv2/GluGluHToMuMu_M-125_13TeV_powheg_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/00000/1E195F26-0284-E811-9812-EC0D9A0B33B0.root'))
    

process.selectedMuons = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("slimmedMuons"),
    cut = cms.string("pt > 10 && abs(eta) < 2.1 && track.isNonnull && isLooseMuon")
)

#Run kinematic fit on first dimuon pair
process.muonKinFit = cms.EDProducer("MuonKinFit",
    srcCands = cms.InputTag("selectedMuons"),
)


process.OutputModule = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string(options.outputFile),
    outputCommands = cms.untracked.vstring('keep *'),
    dataset = cms.untracked.PSet(filterName = cms.untracked.string('')),
    )

process.kinfit = cms.Path(
process.selectedMuons
* process.muonKinFit
)
    
process.e = cms.EndPath(process.OutputModule)

