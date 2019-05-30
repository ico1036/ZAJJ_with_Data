#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"
#include "TH1D.h"
#include "vector"

class MiniAnalyzer : public edm::EDAnalyzer {
   public:
      explicit MiniAnalyzer (const edm::ParameterSet&);
      ~MiniAnalyzer();

   private:
      
	  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void beginJob() override;
      virtual void endJob() override;

      edm::EDGetTokenT<std::vector<pat::Electron> > electronToken_;
      TTree* outTree;
      
	  std::vector<double> v_ele_pt;
      std::vector<double> v_ele_eta;
      std::vector<double> v_ele_phi;
      std::vector<double> v_ele_energy;
      std::vector<double> v_ele_charge;
	
		 int ele_size; 	
	  double ele1_pt;
	  double ele2_pt;
      double ele1_eta;
      double ele2_eta;
      double ele1_phi;
      double ele2_phi;
      double ele1_energy;
      double ele2_energy;
      double ele1_charge;
      double ele2_charge;






};

MiniAnalyzer::MiniAnalyzer(const edm::ParameterSet& iConfig){
    edm::Service<TFileService> fs;
    outTree = fs->make<TTree> ("outTree","ForTest");
    electronToken_ = consumes<std::vector<pat::Electron> >(iConfig.getParameter<edm::InputTag>("electrons"));
    
	outTree->Branch("ele_size",&ele_size);
	
	outTree->Branch("ele1_pt", &ele1_pt);
	outTree->Branch("ele1_eta", &ele1_eta);
	outTree->Branch("ele1_phi", &ele1_phi);
	outTree->Branch("ele1_energy", &ele1_energy);
	outTree->Branch("ele1_charge", &ele1_charge);

	outTree->Branch("ele2_pt", &ele2_pt);
	outTree->Branch("ele2_eta", &ele2_eta);
	outTree->Branch("ele2_phi", &ele2_phi);
	outTree->Branch("ele2_energy", &ele2_energy);
	outTree->Branch("ele2_charge", &ele2_charge);


}

MiniAnalyzer::~MiniAnalyzer()
{
}

void
MiniAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;
   using namespace reco;

   using std::cout;
   using std::endl;

	v_ele_pt.clear();
	v_ele_eta.clear();
	v_ele_phi.clear();
	v_ele_energy.clear();
	v_ele_charge.clear();




    int numOfEle=0;
    edm::Handle<std::vector<pat::Electron> > electrons;
    iEvent.getByToken(electronToken_, electrons);
    
// --Electron Loop started
	for (const pat::Electron &el : *electrons) {
        if (el.pt() < 5) continue;
		
		v_ele_pt.push_back(el.pt());
    	v_ele_eta.push_back(el.eta());
    	v_ele_phi.push_back(el.phi());
    	v_ele_energy.push_back(el.energy());
    	v_ele_charge.push_back(el.charge());



		numOfEle++;
    
	} // Electron Loop ended

	
	if(v_ele_pt.size() == 1){
		
		ele1_pt		=  v_ele_pt[0];
		ele1_eta	=  v_ele_eta[0];
    	ele1_phi	=  v_ele_phi[0];
    	ele1_energy	=  v_ele_energy[0];
    	ele1_charge =  v_ele_charge[0];

	}else if(v_ele_pt.size() >1 ){
		
        ele1_pt     =  v_ele_pt[0];
        ele1_eta    =  v_ele_eta[0];
        ele1_phi    =  v_ele_phi[0];
        ele1_energy =  v_ele_energy[0];
        ele1_charge =  v_ele_charge[0];

        ele2_pt     =  v_ele_pt[1];
        ele2_eta    =  v_ele_eta[1];
        ele2_phi    =  v_ele_phi[1];
        ele2_energy =  v_ele_energy[1];
        ele2_charge =  v_ele_charge[1];

	}	

	
	ele_size = numOfEle;
	outTree->Fill();



}


void MiniAnalyzer::beginJob() { }  
void MiniAnalyzer::endJob() {

   // outTree->Write();
}



DEFINE_FWK_MODULE(MiniAnalyzer);
