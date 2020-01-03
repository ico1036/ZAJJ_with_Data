#include <iostream>
#include <vector>
#include "/hcp/data/data02/jwkim2/WORK/CMSSW_9_4_6_patch1/src/MiniAnalyzer/MiniAnalyzer/src/NpKNU.hh"
#include "/hcp/data/data02/jwkim2/WORK/CMSSW_9_4_6_patch1/src/MiniAnalyzer/MiniAnalyzer/src/NpKNU.cc"
using namespace std;

void makeHist(){



	gSystem->Load("/hcp/data/data02/jwkim2/WORK/CMSSW_9_4_6_patch1/src/MiniAnalyzer/MiniAnalyzer/src/libNpKNU.so");	
	TChain *inChain = new TChain("MiniAnalyzer/NpKNU");



	///  Manually set MC or Data  ///	
	/////////////////////////////////////////////////////////////////////////
	bool isMC = true;
	bool isSig = false;
	TFile *f1 = new TFile("QCDLLAJJ_ele_pho_sel.root","recreate"); // for mc
	//TFile *f1 = new TFile("Data_ele_pho_sel.root","recreate"); //for data
	
	//TFile *f1 = new TFile("LLAJJ_ele_pho_sel.root","recreate"); // for mc(signal)
	/////////////////////////////////////////////////////////////////////////


	if(isMC){
		
		if(isSig){
			inChain->Add("/hcp/data/data02/jwkim2/WORK/CMSSW_9_4_9_cand2/src/MiniAnalyzer/Ntuple/MC/LLAJJ_EWK.root"); // for sig
		}else{
			//inChain->Add("/hcp/data/data02/jwkim2/WORK/CMSSW_9_4_9_cand2/src/MiniAnalyzer/Ntuple/MC/DYjet.root"); // for bkg
			inChain->Add("/hcp/data/data02/jwkim2/WORK/CMSSW_9_4_9_cand2/src/MiniAnalyzer/Ntuple/MC/LLAJJ_QCD.root"); // for bkg
		}

	}else{
			inChain->Add("/hcp/data/data02/jwkim2/WORK/CMSSW_9_4_9_cand2/src/MiniAnalyzer/Ntuple/DoubleEG_GT_Run2016B/Data.root"); //for data
		}	
	

	TClonesArray *eleTCA	 = new TClonesArray("npknu::Electron");	inChain->SetBranchAddress("electron",&eleTCA);
	TClonesArray *eleSelTCA	 = new TClonesArray("npknu::Electron");
	
	TClonesArray *phoTCA	 = new TClonesArray("npknu::Photon");	inChain->SetBranchAddress("photon",&phoTCA);
	TClonesArray *phoSelTCA	 = new TClonesArray("npknu::Photon");
	
	TClonesArray *triggerTCA = new TClonesArray("npknu::Trigger");	inChain->SetBranchAddress("trigger",&triggerTCA);
	
	// Histograms
	TH1D *h1_Mee = new TH1D("h1_Mee","h1_Mee",10000,0,1000);
	TH1D *h1_e1PT = new TH1D("h1_e1PT","h1_e1PT",10000,0,4000);
	TH1D *h1_e2PT = new TH1D("h1_e2PT","h1_e2PT",10000,0,4000);
	TH1D *h1_phoPT = new TH1D("h1_phoPT","h1_phoPT",10000,0,4000);

	
	// Electron ID Eff
	//TEfficiency* pEff = new TEfficiency("eff","ElectronID Tight",1000,0,1000);
	//pEff->SetTitle("ElectronID Tight; p_{T} ; Eff");

	int tot_evt = inChain->GetEntries();
	int tri_cnt=0;
	int ele_tri_idx=309;
	if(isMC){ele_tri_idx=368;}
	
	int ele_cnt=0;
	int ele_ID_cnt=0;
	int ele_Additional=0;
	
	int pho_cnt=0;
	int pho_ID_cnt=0;
	int pho_Additional=0;
	
	int cnt_two_electrons=0;
	int cnt_one_photon=0;
	
	int cnt_Z_window=0;


	cout <<"isMC?: "  << isMC << endl;
	cout <<"isSig?: " << isSig << endl;
	cout <<"Electron trigger nameIdx: " << ele_tri_idx << endl;
	cout <<"Total Event: " <<tot_evt << endl;
	

	
// --EventLoop
	for(int evt_Loop=0; evt_Loop < tot_evt; evt_Loop++){
		inChain->GetEntry(evt_Loop);

		eleSelTCA->Clear("C");
		phoSelTCA->Clear("C");


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


	// ---Electron Loop start
		
		for(int eleLoop=0; eleLoop<eleTCA->GetEntries(); eleLoop++){
            npknu::Electron *elePtr = (npknu::Electron*)eleTCA->At(eleLoop);
			ele_cnt++;
			

			// --Online selection
			//if(elePtr->vidIsPassVeto == 0) continue;			
			//if(elePtr->vidIsPassLoose == 0) continue;			
			if(elePtr->vidIsPassMedium == 0) continue;			
			//if(elePtr->vidIsPassTight == 0) continue;			
			
			//pEff->Fill(elePtr->vidIsPassVeto,elePtr->pt);
			//pEff->Fill(elePtr->vidIsPassLoose,elePtr->pt);
			//pEff->Fill(elePtr->vidIsPassMedium,elePtr->pt);
			//pEff->Fill(elePtr->vidIsPassTight,elePtr->pt);
			
			ele_ID_cnt++;
			// --Offline selection
			if(elePtr->Et() <= 20) continue;
			if(elePtr->Pt() < 20) continue;
			if(fabs(elePtr->eta) >= 2.5) continue;
			if(fabs(elePtr->superCluster_eta) >= 2.5) continue;
			if(fabs(elePtr->superCluster_eta) > 1.4442 &&  fabs(elePtr->superCluster_eta) < 1.566 ) continue;
				

			ele_Additional++;
			new ((*eleSelTCA)[(int)eleSelTCA->GetEntries()]) npknu::Electron(*elePtr);	
			
		} // -- End Electron Loop
	
	// ---Photon Loop start
		
		for(int phoLoop=0; phoLoop<phoTCA->GetEntries(); phoLoop++){
            npknu::Photon *phoPtr = (npknu::Photon*)phoTCA->At(phoLoop);
			pho_cnt++;
			

			// --Online selection
			//if(elePtr->vidIsPassVeto == 0) continue;			
			//if(elePtr->vidIsPassLoose == 0) continue;			
			if(phoPtr->vidIsPassMedium == 0) continue;			
			//if(elePtr->vidIsPassTight == 0) continue;			
			
			//pEff->Fill(elePtr->vidIsPassVeto,elePtr->pt);
			//pEff->Fill(elePtr->vidIsPassLoose,elePtr->pt);
			//pEff->Fill(elePtr->vidIsPassMedium,elePtr->pt);
			//pEff->Fill(elePtr->vidIsPassTight,elePtr->pt);
			
			pho_ID_cnt++;
			// --Offline selection
			
			
			if(fabs(phoPtr->superCluster_eta) > 1.4442 &&  fabs(phoPtr->superCluster_eta) < 1.566 ) continue;
			if(phoPtr->Et() <= 25 ) continue;
			pho_Additional++;

			new ((*phoSelTCA)[(int)phoSelTCA->GetEntries()]) npknu::Photon(*phoPtr);	
	
		} // -- End Photon Loop



//  STEP.1  --- Start Electron pair selection
		
		if(eleSelTCA->GetEntries() < 2) continue;
		npknu::Electron* elePtr1 = (npknu::Electron*)eleSelTCA->At(0); 
		npknu::Electron* elePtr2 = (npknu::Electron*)eleSelTCA->At(1);
		if(elePtr1->charge * elePtr2->charge > 0) continue;
		cnt_two_electrons++;
	

//  STEP.2  --- Start Photon selection

		if(phoSelTCA->GetEntries() < 1) continue;
		npknu::Photon* phoPtr1 = (npknu::Photon*)phoSelTCA->At(0); 
		cnt_one_photon++;		


//  Kinematics --- Z boson window 
		TLorentzVector eTVec1 = elePtr1->GetP4();
		TLorentzVector eTVec2 = elePtr2->GetP4();
		TLorentzVector eeTVec = eTVec1+eTVec2;
		double Mee = eeTVec.M();
		if(Mee < 70 || Mee > 110) continue;
		//if(Mee < 60 || Mee > 120) continue;
		cnt_Z_window++;
		
// --Fill hist
		
	// Electron (STEP 1)
		h1_Mee->Fill(Mee);
		h1_e1PT->Fill(elePtr1->pt);
		h1_e2PT->Fill(elePtr2->pt);
	// Photon   (STEP 2)
		h1_phoPT->Fill(phoPtr1->pt);

	} // --Event Loop Ended


	cout << "Passing Ele trigger: "<< tri_cnt << endl;
	cout << "	" << endl;
	cout << " Electrons ==============================" << endl;
	cout << " Total Electrons: " << ele_cnt << endl;
	cout << " Passing ID Electrons: " << ele_ID_cnt << endl;
	cout << " Passing Ele additional: " << ele_Additional << endl;
	cout << "	" << endl;
	
	cout << " Photons ==============================" << endl;
	cout << " Total Photons: " << pho_cnt << endl;
	cout << " Passing ID Photons: " << pho_ID_cnt << endl;
	cout << " Passing pho additional: " << pho_Additional << endl;
	cout << "	" << endl;

	cout << " Kinematics  ==============================" << endl;

	cout << "Passing Electron pair selection: " << cnt_two_electrons << endl;
	cout << "Passing Photon selection: " << cnt_one_photon << endl;
	cout << "Passing Z mass window: " << cnt_Z_window << endl;

	// --Write file
	f1->Write();

} // END PROGRAM 
