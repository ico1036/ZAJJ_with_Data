### Make histogram for **Pileup re-weighting** steps  
1. The purpose of pileup re-weighting is to match MC pileup ( random distribution ) to real pile-up distribution from data ( saved as Json file, you need to convert it as root file histogram )  

2. Prepare Json file  
	- Year 2015 -->  /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/PileUp/pileup_latest.txt
	- Year 2016 -->  /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/pileup_latest.txt
pileup_latest.txt files is the most up-to-date file  

3. Convert Json file to Root file (histogram)  
```python
pileupCalc.py -i PATH_of_Golden_Json_file --inputLumiJSON PATH_of_Input_Lumi_Json_file --calcMode true --minBiasXsec 71300 --maxPileupBin 100 --numPileupBins 100 Output_file_name.root
```

The output file from this code is the standard for MC pile-up re-weighting  

