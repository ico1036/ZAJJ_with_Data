
## Summary code
1. RecoEgamma/EgammaTools/test/runEgammaPostRecoTools.py: Embed ID info to miniAOD and output new miniAOD  
2. plugins/MiniAnalyzer.cc: cmssw code  
3. test/miniaodrun_cfg.py: config code  

## Status  
**CMSSW_9_4_9_cand2** or higher version is needed for analysis 2016+2017 data  
*This info based on following twiki, 2016/2017 Data/MC*  
https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaPostRecoRecipes  
  
Now Im trying this step:  

#### 1. 2016/2017 Data/MC initial setting(See below) -> "RecoEgamma" directory is created  
#### 2. Use **RecoEgamma/EgammaTools/test/runEgammaPostRecoTools.py** and make new MiniAOD whose ElectronIDs are embedded
* Scripts --> *run.sh*  
#### 3. Analyze new MiniAOD files with ElectronIDs -> plugin/MiniAnalyzer.cc and test/miniaodrun_cfg.py
#### 4. Reference code:  
cmssw(c++): https://github.com/aknayak/JetMetTrigger/blob/master/Analyzer/plugins/HLTJetMETNtupleProducer.cc  
config(python): https://github.com/aknayak/JetMetTrigger/blob/master/Analyzer/test/hltJetMETNtuple_cfg.py  
process.egammaPostRecoSeq* must be embedded in process.p = cms.Path(...)



---

## 2016/2017 Data/MC initial setting  
```bash
cmsrel CMSSW_9_4_13
cd CMSSW_9_4_13/src
cmsenv
git cms-init
git cms-merge-topic cms-egamma:EgammaPostRecoTools #just adds in an extra file to have a setup function to make things easier
scram b -j 8
```



