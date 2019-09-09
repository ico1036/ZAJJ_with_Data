#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "MiniAnalyzer/MiniAnalyzer/src/NpKNU.hh"
#include "vector"


class MiniAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit MiniAnalyzer(const edm::ParameterSet&);
      ~MiniAnalyzer();

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

// ---------Electron: 1. EDGet Token 

	edm::EDGetTokenT<edm::View<pat::Electron> > electronToken_;
	edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken;
	edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjectsToken;
	edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescalesToken;
	std::vector<std::string> triggerIdentifiers;
	std::vector<std::string> triggerLabelsName ;
	int TriggerNameToInt(std::string name);
	bool DoPrintTrigger;
     
	TTree* outTree				;
	TClonesArray* evtTCA		;
	TClonesArray* electronTCA	;	  
	TClonesArray* triggerTCA            ;
	TClonesArray* triggerObjectTCA      ;
};


MiniAnalyzer::MiniAnalyzer(const edm::ParameterSet& iConfig){
    edm::Service<TFileService> fs;
    outTree = fs->make<TTree> ("NpKNU","NpKNU");
    
	evtTCA = new TClonesArray("npknu::Evt"); outTree->Branch("evt"      , &evtTCA );
	electronTCA = new TClonesArray("npknu::Electron"); outTree->Branch("electron"      , &electronTCA );
	

// ---------Elecgron: 2. Linking Input tag 
	electronToken_		 = consumes<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("electrons"));


	triggerTCA = new TClonesArray("npknu::Trigger"); outTree->Branch("trigger"      , &triggerTCA );
    triggerObjectTCA = new TClonesArray("npknu::TriggerObject"); outTree->Branch("triggerObject"      , &triggerObjectTCA );
	triggerResultsToken = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"));
	triggerObjectsToken = consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("triggerObjects"));
	triggerPrescalesToken = consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("triggerPrescales"));
    triggerIdentifiers = iConfig.getParameter<std::vector<std::string> >("triggerIdentifiers");
    triggerLabelsName = iConfig.getParameter<std::vector<std::string> >("triggerLabelsName");
    int ThisGlobalTriggerIndex = 0;
   for(std::vector<std::string>::iterator it = triggerLabelsName.begin(); it!= triggerLabelsName.end(); it++) {
      ThisGlobalTriggerIndex++;
      size_t lastSize = it->find("_v");
      if(lastSize != std::string::npos) it->resize(lastSize+2);
      std::cout << "### IndexTriggerLabelName " << ThisGlobalTriggerIndex << " " << *it << std::endl;
  }
    DoPrintTrigger = true;
}

MiniAnalyzer::~MiniAnalyzer(){}

void MiniAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
   using namespace std;
   using namespace reco;

   using std::cout;
   using std::endl;


// START Event
	evtTCA->Clear("C");
	npknu::Evt* evtObjPtr = new ((*evtTCA)[(int)evtTCA->GetEntries()])npknu::Evt() ;
	evtObjPtr->run = (unsigned int) iEvent.id().run();
	evtObjPtr->lumi = (unsigned int) iEvent.id().luminosityBlock();
	evtObjPtr->event = (unsigned long long)iEvent.id().event();
	evtObjPtr->isRealData = iEvent.isRealData();
// END Event


// START Trigger
 
	triggerTCA->Clear("C");
    triggerObjectTCA->Clear("C");
	edm::Handle<edm::TriggerResults> triggerResults ;                     iEvent.getByToken(triggerResultsToken, triggerResults);
	edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;   iEvent.getByToken(triggerObjectsToken, triggerObjects);
	edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;            iEvent.getByToken(triggerPrescalesToken, triggerPrescales);
	const edm::TriggerNames &names = iEvent.triggerNames(*triggerResults);

	if(DoPrintTrigger) {
		std::cout << "\n === TRIGGER PATHS === " << std::endl;
		for (unsigned int i = 0, n = triggerResults->size(); i < n; ++i) {
			int thisTriIdx = TriggerNameToInt(names.triggerName(i));
			std::cout << evtObjPtr << " TriggerNameAll " << names.triggerName(i) << " " << thisTriIdx << std::endl;
        }
        DoPrintTrigger = false;
   }

	std::vector<std::string> checkTriggerPassVec ;
    for (unsigned int i = 0; i < triggerResults->size(); i++) {
		if(triggerResults->accept(i)) {
			npknu::Trigger* triggerPassPtr = new ((*triggerTCA)[(int)triggerTCA->GetEntries()]) npknu::Trigger();
            triggerPassPtr->nameIdx = TriggerNameToInt(names.triggerName(i));
            if((triggerPassPtr->nameIdx) == -1) {
				triggerPassPtr->name = names.triggerName(i);
            }
            triggerPassPtr->accept = triggerResults->accept(i);
            triggerPassPtr->prescaleForIndex = triggerPrescales->getPrescaleForIndex(i);
            checkTriggerPassVec.push_back(names.triggerName(i));
        }
   }

 std::vector<std::string> checkTriggerPathVec ;
    for(pat::TriggerObjectStandAlone triObj : *triggerObjects) { // note: not "const &" since we want to call unpackPathNames
      triObj.unpackPathNames(names);
      triObj.unpackFilterLabels(iEvent, *triggerResults); //2017 update

        npknu::TriggerObject* triggerObjectPtr = new ((*triggerObjectTCA)[(int)triggerObjectTCA->GetEntries()]) npknu::TriggerObject();
        triggerObjectPtr->SetPtEtaPhiE(triObj.pt(), triObj.eta(), triObj.phi(), triObj.energy());


      std::vector<std::string> pathNamesAll  = triObj.pathNames(false);
      std::vector<std::string> pathNamesLast = triObj.pathNames(true);
        size_t numAllCheckToLast = 0;
      for(unsigned i = 0; i < pathNamesAll.size(); i++) {
         bool isNone = triObj.hasPathName(pathNamesAll[i], false, false );
         bool isL3   = triObj.hasPathName(pathNamesAll[i], false, true );
         bool isLF   = triObj.hasPathName(pathNamesAll[i], true, false );
         bool isBoth = triObj.hasPathName(pathNamesAll[i], true, true );
            unsigned int pathNameBits = 0;
         if (isNone) pathNameBits += (1<<0);
         if (isL3  ) pathNameBits += (1<<1);
         if (isLF  ) pathNameBits += (1<<2);
         if (isBoth) pathNameBits += (1<<3);
            if(pathNameBits == 15) numAllCheckToLast++;
            int ThisTriNameIdx = TriggerNameToInt(pathNamesAll[i]);
            triggerObjectPtr->pathNameIdxVec.push_back(ThisTriNameIdx);
            if(ThisTriNameIdx == -1) {
                triggerObjectPtr->pathNameVec.push_back(pathNamesAll[i]);
            }
            triggerObjectPtr->pathBitVec.push_back(pathNameBits);
	}

if(numAllCheckToLast != pathNamesLast.size()) cout << "### Error TriggerPathCheck All " << pathNamesAll.size() << " AllToLast " << numAllCheckToLast << " Last " << pathNamesLast.size() << std::endl;
    } //triggerObject loop over

// END Trigger




// START Electron
	electronTCA->Clear("C");
    



// ---------Elecgron: 3. Handle 
	edm::Handle<edm::View<pat::Electron> > electrons; iEvent.getByToken(electronToken_, electrons);


	// --Electron Loop
	for(edm::View<pat::Electron>::const_iterator ele = electrons->begin(); ele != electrons->end(); ele++){
		if(ele->pt() < 5.0) continue;
		
		const edm::Ptr<pat::Electron> elePtr(electrons, ele - electrons->begin());
		npknu::Electron* electronPtr = new ((*electronTCA)[(int)electronTCA->GetEntries()]) npknu::Electron();
        electronPtr->SetPtEtaPhiE(elePtr->pt(), elePtr->eta(), elePtr->phi(), elePtr->energy());
    



		// ---Electron ID
		electronPtr->vidIsPassVeto	 = ele->electronID("cutBasedElectronID-Summer16-80X-V1-veto")  ;
		electronPtr->vidIsPassLoose  = ele->electronID("cutBasedElectronID-Summer16-80X-V1-loose") ;
		electronPtr->vidIsPassMedium = ele->electronID("cutBasedElectronID-Summer16-80X-V1-medium");
		electronPtr->vidIsPassTight	 = ele->electronID("cutBasedElectronID-Summer16-80X-V1-tight") ;
		
		// ---Other variables
		electronPtr->gsfTrack_Px      = ele->gsfTrack()->px()  ;
		electronPtr->gsfTrack_Py      = ele->gsfTrack()->py()  ;
		electronPtr->gsfTrack_Pz      = ele->gsfTrack()->pz()  ;
		electronPtr->gsfTrack_Pt      = ele->gsfTrack()->pt()  ;
		electronPtr->superCluster_x   = (ele->superCluster().isNonnull()) ?  ele->superCluster()->x()    : std::numeric_limits<float>::max();
		electronPtr->superCluster_y   = (ele->superCluster().isNonnull()) ?  ele->superCluster()->y()    : std::numeric_limits<float>::max();
		electronPtr->superCluster_z   = (ele->superCluster().isNonnull()) ?  ele->superCluster()->z()    : std::numeric_limits<float>::max();
		electronPtr->superCluster_eta = (ele->superCluster().isNonnull()) ?  ele->superCluster()->eta()  : std::numeric_limits<float>::max();
		electronPtr->superCluster_phi = (ele->superCluster().isNonnull()) ?  ele->superCluster()->phi()  : std::numeric_limits<float>::max();
		electronPtr->superCluster_seed_eta = (ele->superCluster().isNonnull() && ele->superCluster()->seed().isNonnull()) ? ele->superCluster()->seed()->eta() : std::numeric_limits<float>::max();
		electronPtr->superCluster_seed_phi = (ele->superCluster().isNonnull() && ele->superCluster()->seed().isNonnull()) ? ele->superCluster()->seed()->phi() : std::numeric_limits<float>::max() ;
		electronPtr->charge                 = ele->charge()                 ;
		electronPtr->caloEnergy             = ele->caloEnergy()             ;
		electronPtr->ecalEnergy             = ele->ecalEnergy()             ;
		electronPtr->superCluster_energy    = ele->superCluster()->energy() ;
		electronPtr->scE1x5                 = ele->scE1x5()         ;
		electronPtr->scE2x5Max              = ele->scE2x5Max()      ;
		electronPtr->scE5x5                 = ele->scE5x5()         ;
		electronPtr->scPixCharge            = ele->scPixCharge()    ;
		electronPtr->scSigmaEtaEta          = ele->scSigmaEtaEta()  ;
		electronPtr->scSigmaIEtaIEta        = ele->scSigmaIEtaIEta();
		electronPtr->e2x5Max                = ele->e2x5Max() ;
		electronPtr->e5x5                   = ele->e5x5()    ;
		electronPtr->e1x5                   = ele->e1x5()    ;
		electronPtr->full5x5_sigmaIetaIeta  = ele->full5x5_sigmaIetaIeta() ;
		electronPtr->sigmaIetaIphi          = ele->sigmaIetaIphi()         ;
		electronPtr->sigmaIphiIphi          = ele->sigmaIphiIphi()         ;
		electronPtr->gsfTrack_dxy           = ele->gsfTrack()->dxy()       ;
		electronPtr->gsfTrack_dz            = ele->gsfTrack()->dz()        ;
		electronPtr->dr03EcalRecHitSumEt              = ele->dr03EcalRecHitSumEt()            ;
		electronPtr->dr03HcalDepth1TowerSumEt         = ele->dr03HcalDepth1TowerSumEt()       ;
		electronPtr->dr03HcalTowerSumEt               = ele->dr03HcalTowerSumEt()             ;
		electronPtr->dr03TkSumPt                      = ele->dr03TkSumPt()                    ;
		electronPtr->ecalDrivenSeed                   = (bool)ele->ecalDrivenSeed()           ;
		electronPtr->deltaEtaSuperClusterTrackAtVtx   = ele->deltaEtaSuperClusterTrackAtVtx() ;
		electronPtr->deltaPhiSuperClusterTrackAtVtx   = ele->deltaPhiSuperClusterTrackAtVtx() ;
		electronPtr->ecalPFClusterIso                 = ele->ecalPFClusterIso()  ;
		electronPtr->hcalPFClusterIso                 = ele->hcalPFClusterIso()  ;
		electronPtr->hcalOverEcal                     = ele->hcalOverEcal()                   ;
		electronPtr->hadronicOverEm                   = ele->hadronicOverEm()                 ;
		//electronPtr->nMissingHits                     = ele->gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS);
		electronPtr->nMissingHits                     = ele->gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS);
		electronPtr->passConversionVeto               = ele->passConversionVeto()   ;
		electronPtr->fbrem                            = ele->fbrem()                ;
		electronPtr->caloIso                          = ele->caloIso()              ;
		electronPtr->ecalIso                          = ele->ecalIso()              ;
		electronPtr->hcalIso                          = ele->hcalIso()              ;
		electronPtr->trackIso                         = ele->trackIso()             ;
		electronPtr->isPF                             = ele->isPF()                 ;
		electronPtr->r9                               = ele->r9()                   ;
		electronPtr->eSuperClusterOverP               = ele->eSuperClusterOverP()   ;
		electronPtr->ooEmooP                          = fabs((1.0 - ele->eSuperClusterOverP()) / ele->ecalEnergy());


	} // EleLoop ended
//END Electron
	
	outTree->Fill();
}


int MiniAnalyzer::TriggerNameToInt(std::string name) {
   unsigned int index=0;
   TString nameStr(name);
   for(std::vector<std::string>::iterator it= triggerLabelsName.begin(); it!= triggerLabelsName.end(); it++) {
      index++;
      if(nameStr.Contains(*it)) return index;
   }
   std::cout << "### NotFound TriggerLabel " << name << std::endl;
   return -1;
}



void MiniAnalyzer::beginJob() { }  
void MiniAnalyzer::endJob() {

   // outTree->Write();
}



DEFINE_FWK_MODULE(MiniAnalyzer);
