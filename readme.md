
## Summary 
1. **plugins/MiniAnalyzer.cc**: cmssw code  
2. **test/miniaodrun_cfg.py**: config code  
3. **Analysis**: Analyze Ntuple, make hist and make plots  
4. **JSON**: BrilCalc, Pile up Calc  

  

## 2016/2017 Data/MC initial setting  

**CMSSW_9_4_9_cand2** or higher version is needed for analysis 2016+2017 data  
*This info based on following twiki, 2016/2017 Data/MC*  
https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaPostRecoRecipes  


```bash
cmsrel CMSSW_9_4_9_cand2
cd CMSSW_9_4_9_cand2/src
cmsenv
git cms-init
git cms-merge-topic cms-egamma:EgammaPostRecoTools #just adds in an extra file to have a setup function to make things easier
scram b -j 8
```



