#include <iostream>
#include <vector>
#include "/hcp/data/data02/jwkim2/WORK/CMSSW_9_4_6_patch1/src/MiniAnalyzer/MiniAnalyzer/src/NpKNU.hh"
using namespace std;


void makeHist(){



	gSystem->Load("/hcp/data/data02/jwkim2/WORK/CMSSW_9_4_6_patch1/src/MiniAnalyzer/MiniAnalyzer/src/libNpKNU.so");	
	TChain *inChain = new TChain("MiniAnalyzer/NpKNU");
	
	bool isMC = false;
	//TFile *f1 = new TFile("hist_DYjet.root","recreate"); // for mc
	TFile *f1 = new TFile("hist_Data.root","recreate"); //for data
	
	if(isMC){
		inChain->Add("/hcp/data/data02/jwkim2/WORK/CMSSW_9_4_6_patch1/src/MiniAnalyzer/Ntuple/MC/DYjet.root"); // for mc
	}else{
		inChain->Add("/hcp/data/data02/jwkim2/WORK/CMSSW_9_4_6_patch1/src/MiniAnalyzer/Ntuple/DoubleEG_GT_Run2016B/Data.root"); //for data
	}
	

	TClonesArray *eleTCA	 = new TClonesArray("npknu::Electron");	inChain->SetBranchAddress("electron",&eleTCA);
	TClonesArray *triggerTCA = new TClonesArray("npknu::Trigger");	inChain->SetBranchAddress("trigger",&triggerTCA);
	


	TH1D *h1_Mee = new TH1D("h1_Mee","h1_Mee",10000,0,1000);

	int tot_evt = inChain->GetEntries();
	int cnt=0;
	int tri_cnt=0;
	int ele_tri_idx=309;
	if(isMC){ele_tri_idx=368;}

	cout <<"isMC?: " << isMC << endl;
	cout <<"Electron trigger nameIdx: " << ele_tri_idx << endl;
	cout <<"Total Event: " <<tot_evt << endl;
	

	
// --EventLoop
	for(int evt_Loop=0; evt_Loop < tot_evt; evt_Loop++){
		inChain->GetEntry(evt_Loop);

		int Passing_trigger=0;
	
	// -------START TRIGGER
		// --Trigger Loop
		for(int triLoop=0; triLoop<triggerTCA->GetEntries(); triLoop++){
			npknu::Trigger* triPtr = (npknu::Trigger*)triggerTCA->At(triLoop);
			//cout << "Trigger nameidx: " << triPtr->nameIdx << "	  "  <<"Trigger accept: " << triPtr->accept << endl;
			//if(triPtr->accept != 1){ accept_cnt++; }
			
				if(triPtr->nameIdx == ele_tri_idx){Passing_trigger++;}

		}//End Trigger Loop

		if(Passing_trigger == 0) continue;
		tri_cnt++;



	// -------START Electron Selection
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
		
		// --Fill hist
		h1_Mee->Fill(Mee);
	} // --Event Loop Ended


	cout << "Passing Ele trigger: "<< tri_cnt << endl;
	cout << "Electron selection: " << cnt << endl;
	// --Write file
	//f1->Write();

}
