setInput("/afs/cern.ch/work/m/makbiyik/public/HALOMUON/CMSSW_7_3_5/src/CommonFSQFramework/Core/test/HalomuonAnalyzer/d17Jun/plotsHalomuon.root")
getAllHistos()

setLegend("data_MinimumBias","Data")
  #SetFillColor(0);
  #SetBorderSize(0);
  #SetFillStyle(0);
 # SetTextFont(23);
#  SetTextSize(22);



draw(["hSector_RMS","hSector_Mean","hist_ch_RMS","hist_ch_Mean"],"int",["data_MinimumBias"])
printCMSPreliminary()
updateCanvas()
saveCanvas()
