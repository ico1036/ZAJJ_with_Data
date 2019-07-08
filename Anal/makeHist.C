#include <iostream>
#include <vector>
#include "/hcp/data/data02/jwkim2/WORK/ZAJJ/CMSSW_8_1_0/src/MiniAnalyzer/MiniAnalyzer/src/NpKNU.hh"
using namespace std;


void makeHist(){


	gSystem->Load("/hcp/data/data02/jwkim2/WORK/ZAJJ/CMSSW_8_1_0/src/MiniAnalyzer/MiniAnalyzer/src/libNpKNU.so");
	//TFile *f1 = new TFile("hist_Data.root","recreate"); //for data
	TFile *f1 = new TFile("hist_DYjet.root","recreate"); // for mc
	
	TChain *inChain = new TChain("MiniAnalyzer/NpKNU");
	//inChain->Add("/hcp/data/data02/jwkim2/WORK/ZAJJ/CMSSW_8_1_0/src/MiniAnalyzer/MiniAnalyzer/Ntuple/DoubleEG_GT_Run2016B/Data.root"); //for data
	inChain->Add("/hcp/data/data02/jwkim2/WORK/ZAJJ/CMSSW_8_1_0/src/MiniAnalyzer/MiniAnalyzer/Ntuple/MC/DYjet.root"); // for mc
	

	TClonesArray *eleTCA = new TClonesArray("npknu::Electron");	inChain->SetBranchAddress("electron",&eleTCA);
	TH1D *h1_Mee = new TH1D("h1_Mee","h1_Mee",10000,0,1000);

	int tot_evt = inChain->GetEntries();
	int cnt=0;
	cout <<"Total Event: " <<tot_evt << endl;
	

	
// --EventLoop
	for(int evt_Loop=0; evt_Loop < tot_evt; evt_Loop++){
		inChain->GetEntry(evt_Loop);


		// --Electron Selection
		if(eleTCA->GetEntries() < 2) continue;
		npknu::Electron* elePtr1 = (npknu::Electron*)eleTCA->At(0); 
		npknu::Electron* elePtr2 = (npknu::Electron*)eleTCA->At(1);
		if(elePtr1->charge * elePtr2->charge > 0) continue;

		//Reconstruct Z mass
		TLorentzVector eTVec1 = elePtr1->GetP4();
		TLorentzVector eTVec2 = elePtr2->GetP4();
		TLorentzVector eeTVec = eTVec1+eTVec2;
		double Mee = eeTVec.M();
		
		// Z mass window
		if(Mee < 60 || Mee > 120) continue;
		cnt++;
		h1_Mee->Fill(Mee);
	} // --Event Loop Ended


	cout <<"Electron selection Event: " << cnt << endl;
	f1->Write();

}
