### Partial luminosity can be calculated using brilcalc  

#### 0. Installation must be done in lxplus  

#### 1. Installation  

```bash  
export PATH=$HOME/.local/bin:/afs/cern.ch/cms/lumi/brilconda-1.1.7/bin:$PATH  
pip install --install-option="--prefix=$HOME/.local" brilws  
```  
#### 2. Run example  

```bash
brilcalc lumi -b "STABLE BEAMS" -i test_DoubleEG_GT_Run2016BJSON.txt --normtag /afs/cern.ch/user/l/lumipro/public/Normtags/normtag_DATACERT.json -u /fb 
```  

 - Info of norm tag: https://twiki.cern.ch/twiki/bin/view/CMS/PdmV
 - Info of brilcalws: https://cms-service-lumi.web.cern.ch/cms-service-lumi/brilwsdoc.html#installation
