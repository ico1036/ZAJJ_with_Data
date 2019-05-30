#include <iostream>
#include <vector>
using namespace std;


void makeHist(){

	TFile *f1 = new TFile("DYjet.root","recreate");
	
	TChain *inChain = new TChain("MiniAnalyzer/outTree");
	inChain->Add("/hcp/data/data02/jwkim2/WORK/ZAJJ/CMSSW_8_1_0/src/MiniAnalyzer/MiniAnalyzer/Ntuple/test_DYNt.root");
	

	   int ele_size;	   inChain->SetBranchAddress("ele_size",&ele_size);
	double ele1_pt;		   inChain->SetBranchAddress("ele1_pt",&ele1_pt);
	double ele1_eta;	   inChain->SetBranchAddress("ele1_eta",&ele1_eta);
	double ele1_phi;	   inChain->SetBranchAddress("ele1_phi",&ele1_phi);
	double ele1_energy;	   inChain->SetBranchAddress("ele1_energy",&ele1_energy);
	double ele1_charge;	   inChain->SetBranchAddress("ele1_charge",&ele1_charge);
		
	double ele2_pt;		   inChain->SetBranchAddress("ele2_pt",&ele2_pt);
	double ele2_eta;	   inChain->SetBranchAddress("ele2_eta",&ele2_eta);
	double ele2_phi;	   inChain->SetBranchAddress("ele2_phi",&ele2_phi);
	double ele2_energy;	   inChain->SetBranchAddress("ele2_energy",&ele2_energy);
	double ele2_charge;	   inChain->SetBranchAddress("ele2_charge",&ele2_charge);

	TH1D *h1_Mee = new TH1D("h1_Mee","h1_Mee",10000,0,1000);

	int tot_evt = inChain->GetEntries();
	int cnt=0;
	cout << tot_evt << endl;
	
// --EventLoop
	for(int evt_Loop=0; evt_Loop < tot_evt; evt_Loop++){
		inChain->GetEntry(evt_Loop);


		if(ele_size < 2) continue;
		if( ele1_charge + ele2_charge != 0 ) continue;
		cnt++;
		
		TLorentzVector eleVec1; eleVec1.SetPtEtaPhiE(ele1_pt,ele1_eta,ele1_phi,ele1_energy);
		TLorentzVector eleVec2; eleVec2.SetPtEtaPhiE(ele2_pt,ele2_eta,ele2_phi,ele2_energy);
		TLorentzVector eeVec = eleVec1 + eleVec2;
		
		h1_Mee->Fill(eeVec.M());
	} // --Event Loop Ended


		cout << cnt << endl;	
		f1->Write();

}
