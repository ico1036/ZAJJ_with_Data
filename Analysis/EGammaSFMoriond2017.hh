#ifndef EGammaSFMoriond2017_HH
#define EGammaSFMoriond2017_HH

#include "TFile.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TString.h"
#include <iostream>

static TFile* eleReco   ;
static TFile* eleMedium ;


static TH2D* h2Data_eleReco   ; static TH2D* h2MC_eleReco   ; static TH2D* h2SF_eleReco   ;
static TH2D* h2Data_eleMedium ; static TH2D* h2MC_eleMedium ; static TH2D* h2SF_eleMedium ;

void EGammaSFInit() {
	TString dir= "/hcp/data/data02/jwkim2/WORK/CMSSW_9_4_9_cand2/src/MiniAnalyzer/Analysis/scale_factors/";
	eleReco   = TFile::Open( dir + "EGM2D_BtoH_GT20GeV_RecoSF_Legacy2016.root");
	eleMedium = TFile::Open( dir + "2016LegacyReReco_ElectronMedium.root");

	h2Data_eleReco   = (TH2D*)eleReco  ->Get("EGamma_EffData2D") ;
	h2Data_eleMedium = (TH2D*)eleMedium->Get("EGamma_EffData2D") ;

	h2MC_eleReco     = (TH2D*)eleReco  ->Get("EGamma_EffMC2D") ;
	h2MC_eleMedium   = (TH2D*)eleMedium->Get("EGamma_EffMC2D") ;


	h2SF_eleReco   = (TH2D*)eleReco  ->Get("EGamma_SF2D");
	h2SF_eleMedium = (TH2D*)eleMedium->Get("EGamma_SF2D");

}

double GetSFEleReco(double scEta, double pt) {
	double thisPt = (pt > 500.0) ? 499.0 : pt;
	if(pt < 25.0) thisPt = 26.0 ;
	int xBin  = h2SF_eleReco->GetXaxis()->FindBin(scEta ) ;
	int yBin  = h2SF_eleReco->GetYaxis()->FindBin(thisPt) ;
	double sf = h2SF_eleReco->GetBinContent(xBin, yBin) ;
	return sf;
}


double GetSFEleID(double scEta, double pt) {
	TH2D* h2SF_EleID = NULL;
	h2SF_EleID = h2SF_eleMedium;
	//if(h2SF_EleID == NULL) { std::cout << "Error NotFound ElectronID SF " << eleType << std::endl; exit(-99) ;}
	double thisPt = (pt > 500.0) ? 499.0 : pt;
	if(thisPt < 10.0) thisPt = 11.0;
	int xBin  = h2SF_EleID->GetXaxis()->FindBin(scEta ) ;
	int yBin  = h2SF_EleID->GetYaxis()->FindBin(thisPt) ;
	double sf = h2SF_EleID->GetBinContent(xBin, yBin) ;
	return sf;
}





double GetEffEGammaEleReco(double scEta, double pt, bool isData) {
	if(pt > 500.0) pt = 499.0;
	TH2D* h2 = (isData) ? h2Data_eleReco : h2MC_eleReco ;
	int xBin = h2->GetXaxis()->FindBin(scEta) ;
	int yBin = h2->GetYaxis()->FindBin(pt   ) ;
	double eff = h2->GetBinContent(xBin, yBin) ;
	return eff;
}

double GetSFEGammaEle2Reco(double eta1, double pt1, double eta2, double pt2) {
   double effData1 = GetEffEGammaEleReco(eta1, pt1,  true) ;
   double effData2 = GetEffEGammaEleReco(eta2, pt2,  true) ;
   double effMC1   = GetEffEGammaEleReco(eta1, pt1, false) ;
   double effMC2   = GetEffEGammaEleReco(eta2, pt2, false) ;
   double effData = 1.0 - ((1.0 - effData1) * (1.0 - effData2)) ;
   double effMC   = 1.0 - ((1.0 - effMC1  ) * (1.0 - effMC2  )) ;
   return (effData / effMC);
}



#endif
