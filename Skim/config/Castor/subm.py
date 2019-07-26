#!/usr/bin/python

from shutil import copyfile
import os
from os import walk
import sys

files = []
#mypath='/eos/cms/store/t0streamer/Data/HIExpress/000/326/822' # 790 791 815 722
#mypath='/eos/cms/store/hidata/HIRun2018A/HIEmptyBX/RAW/v1/000/326/822/00000/'
#mypath='/eos/cms/store/hidata/HIRun2018A/HIEmptyBX/RAW/v1/000/327/078/00000/'
#mypath='/eos/cms/store/hidata/HIRun2018A/HIForward/AOD/PromptReco-v1/000/326/856/00000/'
#mypath='/eos/cms/store/hidata/HIRun2018A/HIEmptyBX/RAW/v1/000/326/856/00000/'
mypath='/eos/cms/store/hidata/HIRun2018A/HICastor/AOD/PromptReco-v2/000/327/148/00000/' # 90FF9E74-2E83-7845-9EFC-608663B093B4.root
for (dirpath, dirnames, filenames) in walk(mypath):
    for filename in filenames:
        files.append(dirpath.strip() + "/" + filename.strip())

#input = open("files.run326815.missing")
#files = input.readlines()
#input.close()

WORK = os.getcwd()
#CONFIG = '/afs/cern.ch/user/u/ulrich/work/public/CMSSW_10_3_0/src/CommonFSQFramework/Skim/config/Castor/treemaker_CastorFEVT_job.py'
CONFIG = '/afs/cern.ch/user/u/ulrich/work/public/CMSSW_10_3_0/src/CommonFSQFramework/Skim/config/Castor/treemaker_CastorPromptReco_job.py'
OUTPUT = '/afs/cern.ch/user/u/ulrich/work/public/commissioning2018/run327148_prompt_castor'
if not os.path.exists(OUTPUT):
    os.makedirs(OUTPUT)

ijob = 0

for filename in files:
    filename = filename.strip()
    if filename == "":
        continue
    
    ijob += 1
    
    jobdir = OUTPUT + "/job" + str(ijob)
    if not os.path.exists(jobdir):
        os.makedirs(jobdir)
        
    cfg = open("./jobs_tmpl.csh")
    lines = cfg.readlines()
    cfg.close()
    
    job = jobdir + "/job" + str(ijob) + ".csh"
    cfg = open(job, "w")
    for line in lines:
        line = line.replace("@OUTPUT@", jobdir)
        line = line.replace("@FILE@", filename)
        line = line.replace("@WORK@", WORK)
        line = line.replace("@CFG@", CONFIG)
        line = line.replace("@IJOB@", str(ijob) )
        cfg.write(line)
    cfg.close()
    
    #sys.exit(1)
    
    cmd = "bsub -q 1nh -J job"+ str(ijob)+ " < " + job 
    os.system(cmd)
    # 1nh 8nm
    # -R "pool>30000"
    
    # if (ijob>3):
    #sys.exit(1)
