#!/bin/bash
#
#
export SmallXAnaDefFile="/afs/cern.ch/work/m/makbiyik/public/HALOMUON/CMSSW_7_3_5/src/CommonFSQFramework/Core/test/HaloMuonAna/MyAnalysis.py"
export SmallXAnaVersion="CommonFSQFramework.Skim.Samples_MuonAna_0T_20150717"
# export SmallXAnaVersion="CommonFSQFramework.Skim.Samples_MuonAna_0T_20150624"
export HaloMuonOutput="/afs/cern.ch/work/m/makbiyik/public/HALOMUON/CMSSW_7_3_5/src/CommonFSQFramework/Core/test/HaloMuonAna/output"
source /cvmfs/cms.cern.ch/crab3/crab.sh
grid-proxy-init
voms-proxy-init --voms cms 