anaVersion="MuonAna_0T_20150723"
anaType="MuonAna_0T"

cbSmartCommand="smartCopy"
cbSmartBlackList=""
cbWMS="https://wmscms.cern.ch:7443/glite_wms_wmproxy_server"
skimEfficiencyMethod="getSkimEff"

sam = {}

sam["data_PAMinBiasUPC_Run2013"]={}
sam["data_PAMinBiasUPC_Run2013"]["crabJobs"]=5736
sam["data_PAMinBiasUPC_Run2013"]["GT"]='GR_R_75_V5A'
sam["data_PAMinBiasUPC_Run2013"]["name"]='data_PAMinBiasUPC_Run2013'
sam["data_PAMinBiasUPC_Run2013"]["isData"]=True
sam["data_PAMinBiasUPC_Run2013"]["numEvents"]=286814246
sam["data_PAMinBiasUPC_Run2013"]["json"]='CommonFSQFramework/Skim/lumi/MinBias_2013_CastorMuonRuns.json'
sam["data_PAMinBiasUPC_Run2013"]["lumiMinBias"]=-1
sam["data_PAMinBiasUPC_Run2013"]["XS"]=-1
sam["data_PAMinBiasUPC_Run2013"]["DS"]='/PAMinBiasUPC/HIRun2013-PromptReco-v1/RECO'

sam["data_PPMinBias_Run2013"]={}
sam["data_PPMinBias_Run2013"]["crabJobs"]=688
sam["data_PPMinBias_Run2013"]["GT"]='GR_R_75_V5A'
sam["data_PPMinBias_Run2013"]["name"]='data_PPMinBias_Run2013'
sam["data_PPMinBias_Run2013"]["isData"]=True
sam["data_PPMinBias_Run2013"]["numEvents"]=34393913
sam["data_PPMinBias_Run2013"]["json"]='CommonFSQFramework/Skim/lumi/MinBias_2013_CastorMuonRuns_ppInterfill.json'
sam["data_PPMinBias_Run2013"]["lumiMinBias"]=-1
sam["data_PPMinBias_Run2013"]["XS"]=-1
sam["data_PPMinBias_Run2013"]["DS"]='/PPMinBias/Run2013A-PromptReco-v1/RECO'


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
