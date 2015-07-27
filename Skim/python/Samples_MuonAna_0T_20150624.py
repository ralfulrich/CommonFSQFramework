anaVersion="MuonAna_0T_20150624"
anaType="MuonAna_0T"

cbSmartCommand="smartCopy"
cbSmartBlackList=""
cbWMS="https://wmscms.cern.ch:7443/glite_wms_wmproxy_server"
skimEfficiencyMethod="getSkimEff"

sam = {}
sam["data_MinimumBias_Run2015A"]={}
sam["data_MinimumBias_Run2015A"]["crabJobs"]=1609
sam["data_MinimumBias_Run2015A"]["GT"]='GR_R_75_V5A'
sam["data_MinimumBias_Run2015A"]["name"]='data_MinimumBias_Run2015A'
sam["data_MinimumBias_Run2015A"]["isData"]=True
sam["data_MinimumBias_Run2015A"]["numEvents"]=80446390
sam["data_MinimumBias_Run2015A"]["pathSE"]='srm://srm-eoscms.cern.ch:8443/srm/v2/server?SFN=/eos/cms/store/group/phys_heavyions/cwohrman/CFF/CastorMuon/MinimumBias/MuonAna_0T_20150624_data_MinimumBias_Run2015A/150624_164551/0000/'
sam["data_MinimumBias_Run2015A"]["pathTrees"]='/XXXTMFTTree/store/group/phys_heavyions/cwohrman/CFF/CastorMuon/MinimumBias/MuonAna_0T_20150624_data_MinimumBias_Run2015A/150624_164551/0000//'
sam["data_MinimumBias_Run2015A"]["json"]='CommonFSQFramework/Skim/lumi/MinBias_CastorMuonRuns_v2.json'
sam["data_MinimumBias_Run2015A"]["lumiMinBias"]=-1
sam["data_MinimumBias_Run2015A"]["XS"]=-1
sam["data_MinimumBias_Run2015A"]["pathPAT"]='/XXXTMFPAT/store/group/phys_heavyions/cwohrman/CFF/CastorMuon/MinimumBias/MuonAna_0T_20150624_data_MinimumBias_Run2015A/150624_164551/0000//'
sam["data_MinimumBias_Run2015A"]["DS"]='/MinimumBias/Run2015A-PromptReco-v1/RECO'

 #sam["data_MinimumBias_Run2015A"]={}
# sam["data_MinimumBias_Run2015A"]["crabJobs"]=1609
# sam["data_MinimumBias_Run2015A"]["GT"]='GR_R_75_V5A'
# sam["data_MinimumBias_Run2015A"]["name"]='data_MinimumBias_Run2015A'
# sam["data_MinimumBias_Run2015A"]["isData"]=True
# sam["data_MinimumBias_Run2015A"]["numEvents"]=80446390
# sam["data_MinimumBias_Run2015A"]["pathSE"]='/tmp/cwohrman/'
# sam["data_MinimumBias_Run2015A"]["pathTrees"]='/tmp/cwohrman/'
# sam["data_MinimumBias_Run2015A"]["json"]='CommonFSQFramework/Skim/lumi/MinBias_CastorMuonRuns_v2.json'
# sam["data_MinimumBias_Run2015A"]["lumiMinBias"]=-1
# sam["data_MinimumBias_Run2015A"]["XS"]=-1
# sam["data_MinimumBias_Run2015A"]["pathPAT"]='/tmp/cwohrman/'
# sam["data_MinimumBias_Run2015A"]["DS"]='/MinimumBias/Run2015A-PromptReco-v1/RECO'

def fixLocalPaths(sam):
        import os,imp
        if "SmallXAnaDefFile" not in os.environ:
            print "Please set SmallXAnaDefFile environment variable:"
            print "export SmallXAnaDefFile=FullPathToFile"
            raise Exception("Whooops! SmallXAnaDefFile env var not defined")

        anaDefFile = os.environ["SmallXAnaDefFile"]
        mod_dir, filename = os.path.split(anaDefFile)
        mod, ext = os.path.splitext(filename)
        f, filename, desc = imp.find_module(mod, [mod_dir])
        mod = imp.load_module(mod, f, filename, desc)

        localBasePathPAT = mod.PATbasePATH
        localBasePathTrees = mod.TTreeBasePATH

        for s in sam:
            if "pathPAT" in sam[s]:
                sam[s]["pathPAT"] = sam[s]["pathPAT"].replace("XXXTMFPAT", localBasePathPAT)
            if "pathTrees" in sam[s]:
                sam[s]["pathTrees"] = sam[s]["pathTrees"].replace("XXXTMFTTree", localBasePathTrees)
            #print sam[s]["pathPAT"]
            #print sam[s]["pathTrees"]
        return sam
sam = fixLocalPaths(sam)
