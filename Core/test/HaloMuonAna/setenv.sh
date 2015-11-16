#!/bin/bash
#
#
export SmallXAnaDefFile="/afs/cern.ch/work/m/makbiyik/public/HALOMUON/CMSSW_7_3_5/src/CommonFSQFramework/Core/test/HaloMuonAna/MyAnalysis.py"

# this is the 2013 PA/PP dataset
# export SmallXAnaVersion="CommonFSQFramework.Skim.Samples_MuonAna_0T_20150723"

# this is 2015 LHCf-run dataset(use that one )
export SmallXAnaVersion="CommonFSQFramework.Skim.Samples_MuonAna_0T_20150717"

# this is 2015 LHCf-run dataset (not including all interfill runs)
# export SmallXAnaVersion="CommonFSQFramework.Skim.Samples_MuonAna_0T_20150624"

export HaloMuonOutput="/afs/cern.ch/work/m/makbiyik/public/HALOMUON/CMSSW_7_3_5/src/CommonFSQFramework/Core/test/HaloMuonAna/output"
source /cvmfs/cms.cern.ch/crab3/crab.sh
grid-proxy-init
voms-proxy-init --voms cms 