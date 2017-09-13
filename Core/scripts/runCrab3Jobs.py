#!/usr/bin/env python
import sys,os,re,shutil
from optparse import OptionParser

# TODO: voms-proxy-init --voms cms --valid 168:00

import ROOT
ROOT.gROOT.SetBatch(True)

# for fileinpath
from ROOT import *
ROOT.gSystem.Load("libFWCoreFWLite.so")
AutoLibraryLoader.enable()

import CommonFSQFramework.Core.Util

sampleList=CommonFSQFramework.Core.Util.getAnaDefinition("sam")
anaVersion=CommonFSQFramework.Core.Util.getAnaDefinition("anaVersion")
blacklist=CommonFSQFramework.Core.Util.getAnaDefinition("cbSmartBlackList")

def dumpEnvVariable(var):
    ret = "# "+var+"="
    if var in os.environ:
        ret += str(os.environ[var])
    else:
        ret += "<not in env>"
    return ret+"\n"


print "Submitting jobs for", anaVersion


parser = OptionParser(usage="usage: %prog [options] filename",
                        version="%prog 1.0")

parser.add_option("-c", "--cfg", action="store", type="string", dest="cfgInput", default="crabcfg.py" )
parser.add_option("-d", "--dump", action="store_true", dest="dump", default=False )
parser.add_option("-s", "--sample", action="store", type="string", dest="sample" )
parser.add_option("-a", "--allsamples", action="store_true", dest="allsamples", default=False )
#parser.add_option("-d", "--dataOnly", action="store", type="bool", dest="dataOnly" )
(options, args) = parser.parse_args()

if options.dump:
    print ("")
    print ("list of available samples: " + str(sampleList.keys()))
    print ("")
    sys.exit()

elif options.sample:
    sampleListTodo = []
    samplesListFromCLI = options.sample.split(",")
    for s in samplesListFromCLI:
        sampleListTodo.append(s)
elif options.allsamples:
    sampleListTodo = sampleList.keys()
else: 
    parser.print_help()
    sys.exit()

for theSample in sampleListTodo:

  isData=False
  if "isData" in sampleList[theSample]:
    isData=sampleList[theSample]["isData"]

  name=anaVersion+"_"+theSample

  targetPath = anaVersion + "/" + "crab_" + name
  if os.path.exists(targetPath):
    print "Path", targetPath, "already exists. Doing nothing"
    continue    

 


  pycfgextra = []  
  pycfgextra.append("config.General.workArea='"+anaVersion+"'")
  pycfgextra.append("config.General.requestName='"+name+"'")
  pycfgextra.append("config.Data.outputDatasetTag='"+name+"'")
  pycfgextra.append("config.Data.inputDataset='"+sampleList[theSample]["DS"]+"'")
  # customize when running on private datasets
  if "/USER" in sampleList[theSample]["DS"]: 
      print "Submitting jobs with a private USER made input dataset"
      pycfgextra.append("config.Data.inputDBS = 'phys03'")

  
  if isData:
    print isData, sampleList[theSample]["json"]
#    pycfgextra.append("config.Data.splitting='LumiBased'")
#    pycfgextra.append("config.Data.unitsPerJob=100")
    if os.path.isfile(sampleList[theSample]["json"]):
        jsonFile=os.path.abspath(sampleList[theSample]["json"])
    else:
        jsonFile=edm.FileInPath(sampleList[theSample]["json"]).fullPath()
#    else:
#        print ('\nError in json file path/name! fix!\n')
#        sys.exit()
    pycfgextra.append("config.Data.lumiMask='"+jsonFile+"'")
    
  else:
    pycfgextra.append("config.Data.splitting='EventAwareLumiBased'")
    pycfgextra.append("config.Data.unitsPerJob=100000")
  

  # TODO save old value and set it at exit   
  os.environ["TMFSampleName"]=theSample

  if not os.path.isfile(options.cfgInput):
      print ("\nCan't find crabcfg.py. Please specif with -c option!\n")
      sys.exit()

  os.system("cp " + options.cfgInput + " tmp.py")
  with open("tmp.py", "a") as myfile:
    for l in pycfgextra:
        myfile.write(l+"\n")

  os.system("crab submit -c tmp.py")
  shutil.copy("tmp.py", targetPath + "/crabcfg_submitted.py")

  cfgName = None
  with open("tmp.py", "r") as cfg:
    for l in cfg:
        line = l.strip()
        #if "pset=" not in line: continue
        if len(line) > 0 and line[0] == "#": continue
        if "config.JobType.psetName"  not in line: continue
        cfgName = line.split("=")[-1].replace("'","").replace('"',"").strip()


  if not cfgName:
    print "Unable to determine cfg name from crab.cfg!"
  else:
    if not os.path.isfile(cfgName):
        print "Warning: cannot determine the pset. Tried:", cfgName
    else:
        fOut = targetPath + "/" + cfgName
        shutil.copy(cfgName, fOut)

	
sys.exit()  
