#proper gt MCHI2_74_V5

hltGetConfiguration /dev/CMSSW_7_4_0/HIon/V26 --full --offline --mc --unprescale --process TEST --globaltag MCHI2_74_V5 --l1-emulator 'stage1,gt' --l1Xml L1Menu_CollisionsHeavyIons2011_v0_nobsc_notau_centrality_q2_singletrack.v1.xml --output none --max-events 50 --input /store/user/mnguyen/Hydjet_Quenched_MinBias_5020GeV/HydjetMB_740pre8_MCHI2_74_V3_53XBS_DIGI-RAW/6da45e4e90741bc03dbd9aec5f36c050/step2_DIGI_L1_DIGI2RAW_HLT_RAW2DIGI_L1Reco_110_1_aB3.root > hlt_MC_stage1.py
echo '' >> hlt_MC_stage1.py

#both data and mc
echo 'process.simCaloStage1Digis.FirmwareVersion = cms.uint32(1)'>> hlt_MC_stage1.py
echo 'process.caloStage1Params.jetSeedThreshold = cms.double(0.)'>> hlt_MC_stage1.py
echo 'process.RCTConfigProducers.eicIsolationThreshold = cms.uint32(7)' >> hlt_MC_stage1.py
echo 'process.RCTConfigProducers.hOeCut = cms.double(999)' >>hlt_MC_stage1.py
echo 'process.RCTConfigProducers.eMinForHoECut = cms.double(999)' >> hlt_MC_stage1.py
echo 'process.RCTConfigProducers.eMaxForHoECut = cms.double(999)' >> hlt_MC_stage1.py
echo 'process.RCTConfigProducers.hMinForHoECut = cms.double(999)' >> hlt_MC_stage1.py
echo 'process.RCTConfigProducers.eMinForFGCut = cms.double(999)' >> hlt_MC_stage1.py

echo '' >> hlt_MC_stage1.py

#openHLT.root
echo 'process.load("HLTrigger.HLTanalyzers.HLTBitAnalyser_cfi")' >> hlt_MC_stage1.py
echo 'process.hltbitanalysis.HLTProcessName = cms.string("TEST")' >> hlt_MC_stage1.py
echo 'process.hltbitanalysis.hltresults = cms.InputTag( "TriggerResults","","TEST" )' >> hlt_MC_stage1.py
echo 'process.hltbitanalysis.l1GtReadoutRecord = cms.InputTag("simGtDigis","","TEST")' >> hlt_MC_stage1.py
echo 'process.hltbitanalysis.UseTFileService = cms.untracked.bool(True)' >> hlt_MC_stage1.py
echo 'process.hltBitAnalysis = cms.EndPath(process.hltbitanalysis)' >> hlt_MC_stage1.py               
echo 'process.TFileService = cms.Service("TFileService",
                                   fileName=cms.string("openHLT.root"))' >> hlt_MC_stage1.py

echo '' >> hlt_MC_stage1.py

echo '
# for HLT
if hasattr(process, "hltCsc2DRecHits"):
    process.hltCsc2DRecHits.wireDigiTag  = cms.InputTag("simMuonCSCDigis","MuonCSCWireDigi")
    process.hltCsc2DRecHits.stripDigiTag = cms.InputTag("simMuonCSCDigis","MuonCSCStripDigi")
# for the L1 emulator
if hasattr(process, "cscReEmulTriggerPrimitiveDigis"):
    process.cscReEmulTriggerPrimitiveDigis.CSCComparatorDigiProducer = cms.InputTag("simMuonCSCDigis","MuonCSCComparatorDigi")
    process.cscReEmulTriggerPrimitiveDigis.CSCWireDigiProducer = cms.InputTag("simMuonCSCDigis","MuonCSCWireDigi")
' >> hlt_MC_stage1.py

echo '' >> hlt_MC_stage1.py

echo '
from L1Trigger.L1TCommon.customsPostLS1 import customiseSimL1EmulatorForPostLS1_HI
customiseSimL1EmulatorForPostLS1_HI(process)
' >> hlt_MC_stage1.py

echo 'process.L1UpgradeAnalyzer = cms.EDAnalyzer("l1t::L1UpgradeAnalyzer",
                                           InputLayer2Collection = cms.InputTag("simCaloStage1FinalDigis"),
                                           InputLayer2TauCollection = cms.InputTag("simCaloStage1FinalDigis:rlxTaus"),
                                           InputLayer2IsoTauCollection = cms.InputTag("simCaloStage1FinalDigis:isoTaus"),
                                           InputLayer2CaloSpareCollection = cms.InputTag("simCaloStage1FinalDigis:HFRingSums"),
                                           InputLayer1Collection = cms.InputTag("simRctUpgradeFormatDigis"),
                                           legacyRCTDigis = cms.InputTag("simRctDigis")
)
process.L1UpgradeAnalyzerPath = cms.EndPath(process.L1UpgradeAnalyzer)' >> hlt_MC_stage1.py