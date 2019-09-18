#include <iostream>
#include <vector>
#include "/hcp/data/data02/jwkim2/WORK/CMSSW_9_4_6_patch1/src/MiniAnalyzer/MiniAnalyzer/src/NpKNU.hh"
#include "/hcp/data/data02/jwkim2/WORK/CMSSW_9_4_6_patch1/src/MiniAnalyzer/MiniAnalyzer/src/NpKNU.cc"
using namespace std;

void makeHist(){



	gSystem->Load("/hcp/data/data02/jwkim2/WORK/CMSSW_9_4_6_patch1/src/MiniAnalyzer/MiniAnalyzer/src/libNpKNU.so");	
	TChain *inChain = new TChain("MiniAnalyzer/NpKNU");
	
	
	bool isMC = false;
	//bool isMC = true;
	//TFile *f1 = new TFile("hist_DYjet.root","recreate"); // for mc
	TFile *f1 = new TFile("hist_Data.root","recreate"); //for data
	
	if(isMC){
		inChain->Add("/hcp/data/data02/jwkim2/WORK/CMSSW_9_4_9_cand2/src/MiniAnalyzer/Ntuple/MC/DYjet.root"); // for mc
	}else{
		inChain->Add("/hcp/data/data02/jwkim2/WORK/CMSSW_9_4_9_cand2/src/MiniAnalyzer/Ntuple/DoubleEG_GT_Run2016B/Data.root"); //for data
	}
	

	TClonesArray *eleTCA	 = new TClonesArray("npknu::Electron");	inChain->SetBranchAddress("electron",&eleTCA);
	TClonesArray *eleSelTCA	 = new TClonesArray("npknu::Electron");
	
	TClonesArray *triggerTCA = new TClonesArray("npknu::Trigger");	inChain->SetBranchAddress("trigger",&triggerTCA);
	
	// Histograms
	TH1D *h1_Mee = new TH1D("h1_Mee","h1_Mee",10000,0,1000);
	TH1D *h1_e1PT = new TH1D("h1_e1PT","h1_e1PT",10000,0,1000);
	TH1D *h1_e2PT = new TH1D("h1_e2PT","h1_e2PT",10000,0,1000);

	// Electron ID Eff
	 TEfficiency* pEff = new TEfficiency("eff","ElectronID tight",200,0,1000);
	 pEff->SetTitle("ElectronID tight; p_{T} ; Eff");

	int tot_evt = inChain->GetEntries();
	int tri_cnt=0;
	int ele_cnt=0;
	int ele_ID_cnt=0;
	int cnt_two_electrons=0;
	int cnt_Z_window=0;
	int ele_tri_idx=309;
	if(isMC){ele_tri_idx=368;}

	cout <<"isMC?: " << isMC << endl;
	cout <<"Electron trigger nameIdx: " << ele_tri_idx << endl;
	cout <<"Total Event: " <<tot_evt << endl;
	

	
// --EventLoop
	for(int evt_Loop=0; evt_Loop < tot_evt; evt_Loop++){
		inChain->GetEntry(evt_Loop);

		eleSelTCA->Clear("C");


	// -------START TRIGGER
		// --Trigger Loop
		int Passing_trigger=0;
		for(int triLoop=0; triLoop<triggerTCA->GetEntries(); triLoop++){
			npknu::Trigger* triPtr = (npknu::Trigger*)triggerTCA->At(triLoop);
			//cout << "Trigger nameidx: " << triPtr->nameIdx << "	  "  <<"Trigger accept: " << triPtr->accept << endl;
			//if(triPtr->accept != 1){ accept_cnt++; }
			
				if(triPtr->nameIdx == ele_tri_idx){Passing_trigger++;}

		}//End Trigger Loop

		if(Passing_trigger == 0) continue;
		tri_cnt++;

	// -------END TRIGGER


	// -------START Electron Selection
		// ---Electron Loop start
		
		for(int eleLoop=0; eleLoop<eleTCA->GetEntries(); eleLoop++){
            npknu::Electron *elePtr = (npknu::Electron*)eleTCA->At(eleLoop);
			ele_cnt++;
			//if(elePtr->vidIsPassVeto == 0) continue;			
			//if(elePtr->vidIsPassLoose == 0) continue;			
			//if(elePtr->vidIsPassMedium == 0) continue;			
			//if(elePtr->vidIsPassTight == 0) continue;			
			
			//pEff->Fill(elePtr->vidIsPassVeto,elePtr->pt);
			//pEff->Fill(elePtr->vidIsPassLoose,elePtr->pt);
			//pEff->Fill(elePtr->vidIsPassMedium,elePtr->pt);
			pEff->Fill(elePtr->vidIsPassTight,elePtr->pt);


			ele_ID_cnt++;
			new ((*eleSelTCA)[(int)eleSelTCA->GetEntries()]) npknu::Electron(*elePtr);	
	
		}


				

		// ---Grep electron pair
		if(eleSelTCA->GetEntries() < 2) continue;
		npknu::Electron* elePtr1 = (npknu::Electron*)eleSelTCA->At(0); 
		npknu::Electron* elePtr2 = (npknu::Electron*)eleSelTCA->At(1);
		if(elePtr1->charge * elePtr2->charge > 0) continue;
		cnt_two_electrons++;
	

		//Reconstruct Z mass
		TLorentzVector eTVec1 = elePtr1->GetP4();
		TLorentzVector eTVec2 = elePtr2->GetP4();
		TLorentzVector eeTVec = eTVec1+eTVec2;
		double Mee = eeTVec.M();
		
		// Z mass window
		if(Mee < 60 || Mee > 120) continue;
		cnt_Z_window++;
		
		// --Fill hist
		h1_Mee->Fill(Mee);
		h1_e1PT->Fill(elePtr1->pt);
		h1_e2PT->Fill(elePtr2->pt);
	} // --Event Loop Ended


	cout << "Passing Ele trigger: "<< tri_cnt << endl;
	cout << "Passing Electron pair selection: " << cnt_two_electrons << endl;
	cout << "Passing Z mass window: " << cnt_Z_window << endl;
	cout << "	" << endl;
	cout << " Total Electrons: " << ele_cnt << endl;
	cout << " Passing ID Electrons: " << ele_ID_cnt << endl;
	cout << " Electron ID Eff: " << (double)ele_ID_cnt / ele_cnt << endl;

	// --Write file
	f1->Write();


   TCanvas* c1 = new TCanvas("example","",600,400);
   c1->SetFillStyle(1001);
   c1->SetFillColor(kWhite);
	pEff->Draw("AP");
	c1->Print("Eff.png");
}





