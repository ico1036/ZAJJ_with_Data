#include <iostream>
#include <vector>
#include "/hcp/data/data02/jwkim2/WORK/CMSSW_9_4_6_patch1/src/MiniAnalyzer/MiniAnalyzer/src/NpKNU.hh"
#include "/hcp/data/data02/jwkim2/WORK/CMSSW_9_4_6_patch1/src/MiniAnalyzer/MiniAnalyzer/src/NpKNU.cc"
#include "TClonesArray.h"
#include "TFile.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TSystem.h"
#include "/hcp/data/data02/jwkim2/WORK/CMSSW_9_4_9_cand2/src/MiniAnalyzer/Analysis/scale_factors/EGammaSFMoriond2017.hh"



using namespace std;


double DeltaEta(double eta1, double eta2);
double DeltaPhi(double phi1, double phi2);
double DeltaR(double Deta, double Dphi)	 ;

int main(int argc, char** argv){


	// IO
	gSystem->Load("/hcp/data/data02/jwkim2/WORK/CMSSW_9_4_6_patch1/src/MiniAnalyzer/MiniAnalyzer/src/libNpKNU.so");	
	TChain *inChain = new TChain("MiniAnalyzer/NpKNU");
	
	bool isMC;
	bool isSig = false;
	double lumi=215.149415251; 
	double xsec;
	double gen_evt;
	TString outfilename;


	
	std::string action(argv[1]);
	if (action == "QCDLLAJ"){
		isMC = true;
		outfilename = "QCDLLAJJ_ele_pho_sel.root";
		xsec = 47.46;

	}else if ( action == "DYjets"){
		isMC = true;
		outfilename = "DYjets_SFWeight.root";
		gen_evt = 1310938;
		xsec = 6225.42;

	}else if ( action == "data"){
		isMC = false;
		outfilename = "Data.root";
	}else if ( action == "signal"){
		isMC = true;
		isSig = true;
		outfilename = "LLAJJ_ele_pho_sel.root";		
	}else{
		
		cout << "Invalide input data types!! Select:  mc signal  data" << endl;
	}	

	
	TFile *f1 = new TFile(outfilename,"recreate"); // outfile 
	for(int iFile=2; iFile<argc; iFile++){
		cout << "### InFile " << iFile-1 << " " << argv[iFile] << endl;
		inChain->Add(argv[iFile]);
	
	}	
	int tot_evt = inChain->GetEntries();



	// Copy Class
	TClonesArray *pileupTCA	 = new TClonesArray("npknu::Pileup");	inChain->SetBranchAddress("pileup",&pileupTCA);
	
	TClonesArray *eleTCA	 = new TClonesArray("npknu::Electron");	inChain->SetBranchAddress("electron",&eleTCA);
	TClonesArray *eleSelTCA	 = new TClonesArray("npknu::Electron");
	
	TClonesArray *phoTCA	 = new TClonesArray("npknu::Photon");	inChain->SetBranchAddress("photon",&phoTCA);
	TClonesArray *phoSelTCA	 = new TClonesArray("npknu::Photon");
	
	TClonesArray *jetTCA	 = new TClonesArray("npknu::Jet");  inChain->SetBranchAddress("jet",&jetTCA);
	TClonesArray *jetSelTCA	 = new TClonesArray("npknu::Jet");
	
	TClonesArray *triggerTCA = new TClonesArray("npknu::Trigger");	inChain->SetBranchAddress("trigger",&triggerTCA);
	
	int NumVertex; inChain->SetBranchAddress("ntNumVertex",&NumVertex);
	

	// Histograms
	TH1D *h1_NumVertex			= new TH1D("h1_NumVertex","h1_NumVertex",200,0,200);
	TH1D *h1_Normed_NumVertex   = new TH1D("h1_Normed_NumVertex","h1_Normed_NumVertex",200,0,200); 		
    TH1D *h1_PUWeight_NumVertex	= new TH1D("h1_PUWeight_NumVertex","h1_PUWeight_NumVertex",200,0,200);
    TH1D *h1_tWeight_NumVertex	= new TH1D("h1_tWeight_NumVertex","h1_tWeight_NumVertex",200,0,200);

	TH1D *h1_Mee_step1    = new TH1D("h1_Mee_step1","h1_Mee_step1",10000,0,1000);
	//TH1D *h1_e1PT_step1   = new TH1D("h1_e1PT_step1","h1_e1PT_step1",10000,0,4000);
	//TH1D *h1_e2PT_step1	  = new TH1D("h1_e2PT_step1","h1_e2PT_step1",10000,0,4000);
	//TH1D *h1_Mee_step2    = new TH1D("h1_Mee_step2","h1_Mee_step2",10000,0,1000);
	//TH1D *h1_e1PT_step2   = new TH1D("h1_e1PT_step2","h1_e1PT_step2",10000,0,4000);
	//TH1D *h1_e2PT_step2	  = new TH1D("h1_e2PT_step2","h1_e2PT_step2",10000,0,4000);
	//TH1D *h1_Mee_step3    = new TH1D("h1_Mee_step3","h1_Mee_step3",10000,0,1000);
	//TH1D *h1_e1PT_step3   = new TH1D("h1_e1PT_step3","h1_e1PT_step3",10000,0,4000);
	//TH1D *h1_e2PT_step3	  = new TH1D("h1_e2PT_step3","h1_e2PT_step3",10000,0,4000);
	//
	//TH1D *h1_phoPT_step2  = new TH1D("h1_phoPT_step2","h1_phoPT_step2",10000,0,4000);
	//TH1D *h1_phoPT_step3  = new TH1D("h1_phoPT_step3","h1_phoPT_step3",10000,0,4000);
	//
	//TH1D *h1_jet1PT  = new TH1D("h1_jet1PT","h1_jet1PT",10000,0,5000);
	//TH1D *h1_jet2PT  = new TH1D("h1_jet2PT","h1_jet2PT",10000,0,5000);
	//
	//TH1D *h1_jet1Eta = new TH1D("h1_jet1Eta","h1_jet1Eta",100,-5,5);
	//TH1D *h1_jet2Eta = new TH1D("h1_jet2Eta","h1_jet2Eta",100,-5,5);
	//
	//TH1D *h1_jet1Phi = new TH1D("h1_jet1Phi","h1_jet1Phi",50,-3.15,3.15);
	//TH1D *h1_jet2Phi = new TH1D("h1_jet2Phi","h1_jet2Phi",50,-3.15,3.15);
	//
	//TH1D *h1_Mjj	 = new TH1D("h1_Mjj","h1_Mjj",10000,0,5000);
	//TH1D *h1_dEtajj	 = new TH1D("h1_dEtajj","h1_dEtajj",100,-10,10);
	//TH1D *h1_dPhijj	 = new TH1D("h1_dPhijj","h1_dPhijj",50,-3.15,3.15);
	//TH1D *h1_Zepp	 = new TH1D("h1_Zepp","h1_Zepp",100,0,10);

	
	TString MCFile = "/hcp/data/data02/jwkim2/WORK/CMSSW_9_4_9_cand2/src/MiniAnalyzer/Analysis/MCPUhist.root";
	npknu::PUReweight *puWeight = new npknu::PUReweight("../JSON/pileup/MyDataPUHist.root", MCFile, "pileup");
	



	// Counters
	int per99		= tot_evt/99;
	int per100 =0; 

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
	int cnt_two_jets=0;
	
	int jet_cnt=0;
	int jet_sel_cnt=0;

	int cnt_Z_window=0;


	cout <<"isMC?: "  << isMC << endl;
	cout <<"isSig?: " << isSig << endl;
	cout <<"Electron trigger nameIdx: " << ele_tri_idx << endl;

	cout << "tot Event :" << gen_evt << endl;
	cout << "tot Event normalized :" << xsec * lumi << endl;

	cout <<"Goodrun Event: " <<tot_evt << endl;
	cout <<"Goodrun Event normalised: " <<tot_evt * xsec * lumi / gen_evt  << endl;
	
	EGammaSFInit();

	
// --EventLoop
	for(int evt_Loop=0; evt_Loop < tot_evt; evt_Loop++){
		inChain->GetEntry(evt_Loop);
		if((evt_Loop%per99) ==0) cout << "Run" << per100++ << " %" << endl; // Showing progress 

		eleSelTCA->Clear("C");
		phoSelTCA->Clear("C");
		jetSelTCA->Clear("C");

	// -------START TRIGGER
		// --Trigger Loop
		int Passing_trigger=0;
		for(int triLoop=0; triLoop<triggerTCA->GetEntries(); triLoop++){
			npknu::Trigger* triPtr = (npknu::Trigger*)triggerTCA->At(triLoop);
				if(triPtr->nameIdx == ele_tri_idx){Passing_trigger++;}

		}//End Trigger Loop

		if(Passing_trigger == 0) continue;
		tri_cnt++;

	// -------END TRIGGER



	//  STEP.1  --- Start Electron pair selection
		// ---Electron Loop start
		for(int eleLoop=0; eleLoop<eleTCA->GetEntries(); eleLoop++){
            npknu::Electron *elePtr = (npknu::Electron*)eleTCA->At(eleLoop);
			ele_cnt++;
			

			// --Online selection
			//if(elePtr->vidIsPassVeto == 0) continue;			
			//if(elePtr->vidIsPassLoose == 0) continue;			
			if(elePtr->vidIsPassMedium == 0) continue;			
			//if(elePtr->vidIsPassTight == 0) continue;			
			
			
			ele_ID_cnt++;
			// --Offline selection
			if(elePtr->Et() <= 20) continue;
			if(elePtr->pt < 20) continue;
			if(fabs(elePtr->eta) >= 2.5) continue;
			if(fabs(elePtr->superCluster_eta) >= 2.5) continue;
			if(fabs(elePtr->superCluster_eta) > 1.4442 &&  fabs(elePtr->superCluster_eta) < 1.566 ) continue;
				

			ele_Additional++;
			new ((*eleSelTCA)[(int)eleSelTCA->GetEntries()]) npknu::Electron(*elePtr);	
			
		} // -- End Electron Loop

		if(eleSelTCA->GetEntries() < 2) continue;
		npknu::Electron* elePtr1 = (npknu::Electron*)eleSelTCA->At(0); 
		npknu::Electron* elePtr2 = (npknu::Electron*)eleSelTCA->At(1);
		if(elePtr1->charge * elePtr2->charge > 0) continue;
		cnt_two_electrons++;

	// -------START PU RE_WEIGHTING  
		npknu::Pileup* pileupPtr = (npknu::Pileup*)pileupTCA->At(0);
		double pileup = pileupPtr->TrueNumInteractions;
		double pileupWeight =  (isMC) ?  puWeight->getWeight(pileup) : 1 ;
		double LumiWeight   =  (isMC) ?  lumi * xsec / gen_evt : 1 ;
		double tWeight =  (isMC) ? LumiWeight * pileupWeight : 1 ;

	// -- Electron ID,RECO scale factor
		double SFEleReco1 = (isMC) ? GetSFEleReco(elePtr1->superCluster_eta, elePtr1->pt) : 1  ;
		double SFEleReco2 = (isMC) ? GetSFEleReco(elePtr2->superCluster_eta, elePtr2->pt) : 1  ;
		double SFEleId1   = (isMC) ? GetSFEleID(elePtr1->superCluster_eta, elePtr1->pt)   : 1  ;
		double SFEleId2   = (isMC) ? GetSFEleID(elePtr2->superCluster_eta, elePtr2->pt)   : 1 ;
		double SF = tWeight * SFEleReco1 * SFEleReco2 * SFEleId1 * SFEleId2       ;
		
		h1_NumVertex->Fill(NumVertex);
		h1_Normed_NumVertex		->Fill(NumVertex,LumiWeight);
		h1_PUWeight_NumVertex	->Fill(NumVertex,pileupWeight);
		h1_tWeight_NumVertex	->Fill(NumVertex,tWeight);
			
			
	// -------END PU RE_WEIGHTING

		TLorentzVector eTVec1 = elePtr1->GetP4();
		TLorentzVector eTVec2 = elePtr2->GetP4();
		TLorentzVector eeTVec = eTVec1+eTVec2;
		double Mee = eeTVec.M();
		

		//--Fill hist 
			// Electron 
			//h1_Mee_step1->Fill(Mee,1);
			//h1_Mee_step1->Fill(Mee,tWeight);
		h1_Mee_step1->Fill(Mee,SF);
			//h1_e1PT_step1->Fill(elePtr1->pt);
			//h1_e2PT_step1->Fill(elePtr2->pt);




	/*
	
	//  STEP.2  --- Start Photon selection and Zmass window 
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
		
		if(phoSelTCA->GetEntries() < 1) continue;
		npknu::Photon* phoPtr1 = (npknu::Photon*)phoSelTCA->At(0); 
		cnt_one_photon++;		
		
		// -------- Zmass
		if(Mee < 70 || Mee > 110) continue;
		//if(Mee < 60 || Mee > 120) continue;
		cnt_Z_window++;
 
		// --Fill hist
			// Electron 
			h1_Mee_step2->Fill(Mee);
			h1_e1PT_step2->Fill(elePtr1->pt);
			h1_e2PT_step2->Fill(elePtr2->pt);
			// Photon   
			h1_phoPT_step2->Fill(phoPtr1->pt);
	
		
	//  STEP.3  --- Start Jet selection
		// ---Jet Loop start
		for(int jetLoop=0; jetLoop<jetTCA->GetEntries(); jetLoop++){
            npknu::Jet *jetPtr = (npknu::Jet*)jetTCA->At(jetLoop);
			jet_cnt++;
				
			if(jetPtr->pt < 30) continue;
			if(jetPtr->eta >= 4.7) continue;
			
			double dEta_Aj = DeltaEta(jetPtr->eta,phoPtr1->eta);
			double dPhi_Aj = DeltaPhi(jetPtr->phi,phoPtr1->phi);
			double dR_Aj	 = DeltaR(dEta_Aj,dPhi_Aj);
			
			double dEta_e1j = DeltaEta(jetPtr->eta,elePtr1->eta);
			double dPhi_e1j = DeltaPhi(jetPtr->phi,elePtr1->phi);
			double dR_e1j	 = DeltaR(dEta_e1j,dPhi_e1j);
			
			double dEta_e2j = DeltaEta(jetPtr->eta,elePtr2->eta);
			double dPhi_e2j = DeltaPhi(jetPtr->phi,elePtr2->phi);
			double dR_e2j	 = DeltaR(dEta_e1j,dPhi_e2j);

			
			if(dR_Aj < 0.5) continue;
			if(dR_e1j < 0.5) continue;
			if(dR_e2j < 0.5) continue;
						

			jet_sel_cnt++;		
				
			new ((*jetSelTCA)[(int)jetSelTCA->GetEntries()]) npknu::Jet(*jetPtr);	
	
		} // -- End Jet Loop


		if(jetSelTCA->GetEntries() < 2) continue;
		npknu::Jet* jetPtr1 = (npknu::Jet*)jetSelTCA->At(0); 
		npknu::Jet* jetPtr2 = (npknu::Jet*)jetSelTCA->At(1); 
		
		double dEta_jj	 = DeltaEta(jetPtr1->eta,jetPtr2->eta);
		double dPhi_jj	 = DeltaPhi(jetPtr1->phi,jetPtr2->phi);

		cnt_two_jets++;


		// --Other Kinematics
		TLorentzVector aTVec1  = phoPtr1->GetP4();
		TLorentzVector eeaTVec = eeTVec + aTVec1;

		TLorentzVector jTVec1 = jetPtr1->GetP4();
		TLorentzVector jTVec2 = jetPtr2->GetP4();
		TLorentzVector jjTVec = jTVec1 + jTVec2;
			

		double rapZA = eeaTVec.Rapidity() ;
		double rapJ1 = jTVec1.Rapidity()  ;
		double rapJ2 = jTVec2.Rapidity()  ;
		double zepp  = fabs( rapZA - ( rapJ1 + rapJ2 ) / 2.0 );
		double Mjj	 = jjTVec.M() ;

  	
		 //--Fill hist
		  	
		  // Electron 
		  	h1_Mee_step3->Fill(Mee);
		  	h1_e1PT_step3->Fill(elePtr1->pt);
		  	h1_e2PT_step3->Fill(elePtr2->pt);
		  // Photon   (
		  	h1_phoPT_step3->Fill(phoPtr1->pt);
		  // Jet		
			h1_jet1PT	->Fill(jetPtr1->pt);
		    h1_jet2PT	->Fill(jetPtr2->pt);
		              
		    h1_jet1Eta	->Fill(jetPtr1->eta);
		    h1_jet2Eta	->Fill(jetPtr2->eta);
		              
		    h1_jet1Phi	->Fill(jetPtr1->phi);
		    h1_jet2Phi	->Fill(jetPtr2->phi);
		     
			h1_Mjj      ->Fill(Mjj);
		    h1_dEtajj   ->Fill(dEta_jj);
		    h1_dPhijj	->Fill(dPhi_jj);
		    h1_Zepp		->Fill(zepp);
	*/
	} // --Event Loop Ended



	cout << "Passing Trigger: " << tri_cnt << endl;
	cout << "Passing Trigger exp: " << tri_cnt * xsec * lumi / gen_evt << endl;

	cout << "Passing two or more electrons: " << cnt_two_electrons << endl;
	cout << "Passing two or more electrons exp: " << cnt_two_electrons * xsec * lumi / gen_evt  << endl;

//	cout << "Total events: " << tot_evt << endl;
//	cout << "Passing Ele trigger: "<< tri_cnt << endl;
//	cout << "	" << endl;
//	cout << " Electrons ==============================" << endl;
//	cout << " Total Electrons: " << ele_cnt << endl;
//	cout << " Passing ID Electrons: " << ele_ID_cnt << endl;
//	cout << " Passing Ele additional: " << ele_Additional << endl;
//	cout << "	" << endl;
//	
//	cout << " Photons ==============================" << endl;
//	cout << " Total Photons: " << pho_cnt << endl;
//	cout << " Passing ID Photons: " << pho_ID_cnt << endl;
//	cout << " Passing pho additional: " << pho_Additional << endl;
//	cout << "	" << endl;
//	
//	cout << " Jets ==============================" << endl;
//	cout << " Total Jets: " << jet_cnt << endl;
//	cout << " Passing Selected Jets: " << jet_sel_cnt << endl;
//
//	cout << " Kinematics  ==============================" << endl;
//
//	cout << "Passing Electron pair selection: " << cnt_two_electrons << endl;
//	cout << "Passing Photon selection: " << cnt_one_photon << endl;
//	cout << "Passing Z mass window: " << cnt_Z_window << endl;
//	cout << "Passing Jet pair selection: " << cnt_two_jets << endl;

	// --Write file
	f1->Write();

} // END PROGRAM 




double DeltaEta(double eta1, double eta2){
	return fabs(eta1 - eta2);
}

double DeltaPhi(double phi1, double phi2){
	double Dphi = phi1 - phi2;
	double Deltaphi = ( Dphi > TMath::Pi() ) ? fabs(TMath::TwoPi() - Dphi) : fabs(Dphi) ;

	return Deltaphi;
}

double DeltaR(double Deta, double Dphi){
	return TMath::Sqrt(Deta*Deta + Dphi*Dphi);
}
