import FWCore.ParameterSet.Config as cms

def get(todo):
    defs = {}

    # L1 trigger configuration - please do not edit this
    defs["L1GTriggerResultsView"] = cms.PSet(
        miniView = cms.string("TriggerResultsView"),
        branchPrefix = cms.untracked.string("trgl1"),
        storePrescales = cms.bool(False),
        triggers = cms.vstring("L1GTAlgo"),
        useGlobalReadoutRecord = cms.bool(False)
        )

    defs["L1GTriggerResultsViewAOD"] = cms.PSet(
        miniView = cms.string("TriggerResultsView"),
        branchPrefix = cms.untracked.string("trgl1"),
        storePrescales = cms.bool(False),
        triggers = cms.vstring("L1GTAlgo"),
        useGlobalReadoutRecord = cms.bool(True)
        )

    defs["AK4CaloJetTriggerResultsView"]  = cms.PSet(
       miniView = cms.string("TriggerResultsView"),
       branchPrefix = cms.untracked.string("trgAK4Calo"),
       #process = cms.string("HLT"),
       storePrescales = cms.bool(False),
       triggers = cms.vstring("Jet30","Jet40","Jet50"),
       Jet30 = cms.vstring("HLT_AK4CaloJet30ForEndOfFill_v1"),
       Jet40 = cms.vstring("HLT_AK4CaloJet40ForEndOfFill_v1"),
       Jet50 = cms.vstring("HLT_AK4CaloJet50ForEndOfFill_v1")
    )

    defs["AK4CaloJetTriggerResultsViewWithPS"]  = cms.PSet(
       miniView = cms.string("TriggerResultsView"),
       branchPrefix = cms.untracked.string("trgAK4Calo"),
       #process = cms.string("HLT"),
       storePrescales = cms.bool(True),
       triggers = cms.vstring("Jet30","Jet40","Jet50"),
       Jet30 = cms.vstring("HLT_AK4CaloJet30ForEndOfFill_v1"),
       Jet40 = cms.vstring("HLT_AK4CaloJet40ForEndOfFill_v1"),
       Jet50 = cms.vstring("HLT_AK4CaloJet50ForEndOfFill_v1")
    )

    defs["FullTrackTriggerResultsView"]  = cms.PSet(
       miniView = cms.string("TriggerResultsView"),
       branchPrefix = cms.untracked.string("trgTracks"),
       #process = cms.string("HLT"),
       storePrescales = cms.bool(False),
       triggers = cms.vstring("FullTrack12"),
       FullTrack12 = cms.vstring("HLT_FullTrack12ForEndOfFill_v1")
    )

    defs["FullTrackTriggerResultsViewWithPS"]  = cms.PSet(
       miniView = cms.string("TriggerResultsView"),
       branchPrefix = cms.untracked.string("trgTracks"),
       #process = cms.string("HLT"),
       storePrescales = cms.bool(True),
       triggers = cms.vstring("FullTrack12"),
       FullTrack12 = cms.vstring("HLT_FullTrack12ForEndOfFill_v1")
    )

    defs["CastorSpecialJetTriggerResultsView"]  = cms.PSet(
        miniView = cms.string("TriggerResultsView"),
        branchPrefix = cms.untracked.string("CasTrg"),
        #process = cms.string("HLT"),
        storePrescales = cms.bool(False),
        triggers = cms.vstring("ZeroBias","MinBias","Random","CastorMedJet","CastorHighJet","CastorDiJet"),
        ZeroBias = cms.vstring("HLT_ZeroBias*"),
        MinBias = cms.vstring("HLT_L1MinimumBias*"),
        Random = cms.vstring("HLT_Random*"),
        CastorMedJet = cms.vstring("HLT_L1CastorMediumJet_v*"),
        CastorHighJet = cms.vstring("HLT_L1CastorHighJet_v*"),
        CastorDiJet = cms.vstring("HLT_L1CastorMediumJet_PFJet15_v*")
    )

    defs["CastorPATriggerResultsView"]  = cms.PSet(
        miniView = cms.string("TriggerResultsView"),
        branchPrefix = cms.untracked.string("CasPATrg"),
        #process = cms.string("HLT"),
        storePrescales = cms.bool(False),
        triggers = cms.vstring("CastorPAMedJet","CastorPAMuon","Random"),
        CastorPAMedJet = cms.vstring("HLT_PAL1CastorMediumJet_BptxAND_v*"),
        CastorPAMuon = cms.vstring("HLT_PAL1CastorHaloMuon_v*"),
        Random = cms.vstring("HLT_Random*")
    )

 
    # main function
    ret = {}
    for t in todo:
        if t not in defs:
            raise Exception("miniView def not known "+t)

        ret[t] = defs[t]
    return ret


