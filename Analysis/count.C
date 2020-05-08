#include <iostream>
#include <vector>
#include "TClonesArray.h"
#include "TFile.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TSystem.h"


using namespace std;



int main(int argc, char** argv){


	// IO
	TChain *inChain = new TChain("MiniAnalyzer/NpKNU");

	for(int iFile=1; iFile<argc; iFile++){
		cout << "### InFile " << iFile << " " << argv[iFile] << endl;
		inChain->Add(argv[iFile]);
	
		cout << iFile  << " " << argv[iFile] << " " << inChain->GetEntries() << endl;

	}	
	int tot_evt = inChain->GetEntries();
	cout << "Total evt sum: " << tot_evt << endl;
	

	

} // END PROGRAM 

