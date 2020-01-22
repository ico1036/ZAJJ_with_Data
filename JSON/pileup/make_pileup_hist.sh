#!/bin/bash

Golden_Json_PATH='/hcp/data/data02/jwkim2/JSON/Final/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt'
Lumi_Json_PATH='pileup_latest.txt'
Outname='MyDataPUHist.root'

pileupCalc.py -i $Golden_Json_PATH --inputLumiJSON $Lumi_Json_PATH --calcMode true --minBiasXsec 69200 --maxPileupBin 100 --numPileupBins 100 $Outname
