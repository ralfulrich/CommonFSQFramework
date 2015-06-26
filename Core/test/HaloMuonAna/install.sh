export SmallXAnaDefFile=./MyAnalysis.py
export SmallXAnaVersion=CommonFSQFramework.Skim.Samples_CastorMuons_201506012
export HaloMuonOutput=`pwd`"/output"
source /cvmfs/cms.cern.ch/crab3/crab.sh
voms-proxy-init --voms cms
cmsenv
