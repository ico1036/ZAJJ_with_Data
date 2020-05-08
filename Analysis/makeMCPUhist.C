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


using namespace std;



int main(int argc, char** argv){


	// IO
	gSystem->Load("/hcp/data/data02/jwkim2/WORK/CMSSW_9_4_6_patch1/src/MiniAnalyzer/MiniAnalyzer/src/libNpKNU.so");	
	TChain *inChain = new TChain("MiniAnalyzer/NpKNU");

	TString outfilename;
	outfilename = "MCPUhist.root";

	TFile *f = new TFile(outfilename,"recreate");
	for(int iFile=1; iFile<argc; iFile++){
		cout << "### InFile " << iFile << " " << argv[iFile] << endl;
		inChain->Add(argv[iFile]);
	
	}	
	int tot_evt = inChain->GetEntries();
	
	TClonesArray *pileupTCA  = new TClonesArray("npknu::Pileup");   inChain->SetBranchAddress("pileup",&pileupTCA);


	// Histograms
	TH1D *h1_pileup    = new TH1D("pileup","pileup",100,0,100);
	

	// Counters
	int per99		= tot_evt/99;
	int per100 =0; 
	cout <<"Total Event: " <<tot_evt << endl;
	

	
// --EventLoop
	for(int evt_Loop=0; evt_Loop < tot_evt; evt_Loop++){
		inChain->GetEntry(evt_Loop);
		if((evt_Loop%per99) ==0) cout << "Run" << per100++ << " %" << endl; // Showing progress 
	
		

		// START Pileup
		for(int pileupLoop=0; pileupLoop<pileupTCA->GetEntries(); pileupLoop++){
            npknu::Pileup* pileupPtr = (npknu::Pileup*)pileupTCA->At(pileupLoop);
			
			h1_pileup->Fill(pileupPtr->TrueNumInteractions);
        }//End Pileup

	} // --Event Loop Ended


	// --Write file
	f->Write();

} // END PROGRAM 

