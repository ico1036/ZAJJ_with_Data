## Status  
**CMSSW_9_4_9_cand2** or higher version is needed for analysis 2016+2017 data  
*This info based on following twiki, 2016/2017 Data/MC*  
https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaPostRecoRecipes  
  
Now Im trying this step:  
1. Use *EgammaTools/test/runEgammaPostRecoTools.py* and make new MiniAOD whose ElectronIDs are embedded  
https://github.com/ico1036/cmssw-1/tree/CMSSW_9_4_X/RecoEgamma/EgammaTools/test  
2. Read the new MiniAOD files and find out how to access ElectronIDS  
3. Reference code:  
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



