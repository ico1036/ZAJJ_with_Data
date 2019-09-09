#!/bin/bash

###### --MC

#path_='file:/hcp/data/data02/jwkim2/store/mc/RunIISummer16MiniAODv2/701979BD-51C4-E611-9FC9-C4346BC84780.root'
#out='MC'


##### --Data

path_='file:/hcp/data/data02/jwkim2/store/data/Run2016H/DoubleEG/MINIAOD/03Feb2017_ver3-v1/50000/B064E00E-AAEB-E611-9D8E-0CC47A7E018E.root'
out='Data'




cmsRun runEgammaPostRecoTools.py inputFiles=$path_ outputFile=$out era='2016-Legacy'
