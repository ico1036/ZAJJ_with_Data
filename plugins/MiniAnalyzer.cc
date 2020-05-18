#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

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

#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "vector"


class MiniAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit MiniAnalyzer(const edm::ParameterSet&);
      ~MiniAnalyzer();

	

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

// ---------1. EDGet Token 


	// Trigger
	edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken;
	edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjectsToken;
	edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescalesToken;
	std::vector<std::string> triggerIdentifiers;
	std::vector<std::string> triggerLabelsName ;
	
	// Vertex
	edm::EDGetTokenT<reco::VertexCollection> vertexToken ;
	int ntNumVertex ;
    int ntNumGoodVertex ;

	// Pileup
	edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileupToken;
	
	// Particles
	edm::EDGetTokenT<edm::View<pat::Electron> > electronToken_ ;
	edm::EDGetTokenT<edm::View<pat::Photon> > photonToken_ ;
	edm::EDGetTokenT<edm::View<pat::Jet> > jetToken_;


	// Useful definitions....
	int TriggerNameToInt(std::string name);
	bool DoPrintTrigger;
     
	TTree* outTree						;
	TClonesArray* evtTCA				;
	TClonesArray* triggerTCA            ;
	TClonesArray* triggerObjectTCA      ;
	TClonesArray* electronTCA			;	  
	TClonesArray* photonTCA				;	  
	TClonesArray* jetTCA				;	  
	TClonesArray* vertexTCA				;	  
	TClonesArray* pileupTCA				;	  


	int start_num;
	int vertex_num;
	int fill_num;


};   // -- Class End


MiniAnalyzer::MiniAnalyzer(const edm::ParameterSet& iConfig){
    edm::Service<TFileService> fs;
    outTree = fs->make<TTree> ("NpKNU","NpKNU");
    
	evtTCA		= new TClonesArray("npknu::Evt");					outTree->Branch("evt"      , &evtTCA )      ;
	triggerTCA = new TClonesArray("npknu::Trigger");				outTree->Branch("trigger"      , &triggerTCA );
    triggerObjectTCA = new TClonesArray("npknu::TriggerObject");	outTree->Branch("triggerObject"      , &triggerObjectTCA );
    vertexTCA = new TClonesArray("npknu::Vertex");					outTree->Branch("vertex", &vertexTCA );
    pileupTCA = new TClonesArray("npknu::Pileup");					outTree->Branch("pileup", &pileupTCA );
	electronTCA = new TClonesArray("npknu::Electron");				outTree->Branch("electron" , &electronTCA ) ;
	photonTCA	= new TClonesArray("npknu::Photon");				outTree->Branch("photon"   , &photonTCA )   ;
	jetTCA		= new TClonesArray("npknu::Jet");					outTree->Branch("jet"   , &jetTCA )		;
	

// ---------2. Linking Input tag 
	



	// Trigger
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
   //   std::cout << "### IndexTriggerLabelName " << ThisGlobalTriggerIndex << " " << *it << std::endl;
		}
    DoPrintTrigger = true;
	

	// Vertex
	vertexToken          = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertex")) ;
	outTree->Branch("ntNumVertex",&ntNumVertex);
    outTree->Branch("ntNumGoodVertex",&ntNumGoodVertex);

	// Pileup
	pileupToken = (consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("pileup")))   ;

	// Particles
	electronToken_		 = consumes<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("electrons"));
	photonToken_		 = consumes<edm::View<pat::Photon> >(iConfig.getParameter<edm::InputTag>("photons"));
	jetToken_			 = consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jets"));

	start_num  =0 ;
	fill_num   =0 ;
	vertex_num =0 ;

}  // -- Constructer end


MiniAnalyzer::~MiniAnalyzer(){}


// ---------------Event Loop 
void MiniAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;
   using namespace reco;

   using std::cout;
   using std::endl;

	start_num++;
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
		//	std::cout << evtObjPtr << " TriggerNameAll " << names.triggerName(i) << " " << thisTriIdx << std::endl;
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


// START VERTEX
   edm::Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(vertexToken, vertices);
   //if (vertices->empty()) return; //skip the event if no PV found
   reco::VertexCollection::const_iterator firstGoodVertex = vertices->end();
   int firstGoodVertexIdx = 0;
   for(reco::VertexCollection::const_iterator vtx = vertices->begin(); vtx != vertices->end(); ++vtx, ++firstGoodVertexIdx) {
       // Replace isFake() for miniAOD because
       // it requires tracks and miniAOD vertices don't have tracks:
       // Vertex.h: bool isFake() const {return (chi2_==0 && ndof_==0 && tracks_.empty());}
       // bool isFake = vtx->isFake();
       bool isFake =  (vtx->chi2()==0 && vtx->ndof()==0);
       if ( !isFake &&  vtx->ndof()>=4. && vtx->position().Rho()<=2.0 && fabs(vtx->position().Z())<=24.0) {
              firstGoodVertex = vtx;
           break;
           }
    }
   //if ( firstGoodVertex==vertices->end() ) return; // skip event if there are no good PVs
   //std::cout << "good PVs , firstGoodVertexIdx " << firstGoodVertexIdx << std::endl;

    vertexTCA->Clear("C");
    ntNumVertex = 0;
    ntNumGoodVertex = 0;
   for(const reco::Vertex& it : *vertices) {
        ntNumVertex++;
        if( (!it.isFake()) and (it.ndof() >= 4.0) and (it.position().rho() <= 2.0) and (std::fabs(it.z()) <= 24.0) ) ntNumGoodVertex++;
        if(ntNumVertex == 1 || ntNumGoodVertex == 1) {
        npknu::Vertex* vertexPtr = new ((*vertexTCA)[(int)vertexTCA->GetEntries()]) npknu::Vertex();
        vertexPtr->x              = it.x()              ;
        vertexPtr->y              = it.y()              ;
        vertexPtr->z              = it.z()              ;
        vertexPtr->xError         = it.xError()         ;
        vertexPtr->yError         = it.yError()         ;
        vertexPtr->zError         = it.zError()         ;
        vertexPtr->tracksSize     = it.tracksSize()     ;
        vertexPtr->nTracks        = it.nTracks()        ;
        vertexPtr->isFake         = it.isFake()         ;
        vertexPtr->ndof           = it.ndof()           ;
        vertexPtr->position_rho   = it.position().rho() ;
        vertexPtr->chi2           = it.chi2()           ;
        vertexPtr->normalizedChi2 = it.normalizedChi2() ;
        }
   }

	vertex_num++;
// END VERTEX



// START Pileup
    pileupTCA->Clear("C");
    if(!iEvent.isRealData()) {
		edm::Handle<std::vector<PileupSummaryInfo> > PileupInfos; iEvent.getByToken(pileupToken, PileupInfos);
		for (const PileupSummaryInfo& it : *PileupInfos) {
			// In-time bunch crossing
			if(it.getBunchCrossing() == 0) {
				npknu::Pileup* pileupPtr = new ((*pileupTCA)[(int)pileupTCA->GetEntries()]) npknu::Pileup();
				pileupPtr->BunchCrossing       = it.getBunchCrossing();
				pileupPtr->TrueNumInteractions = it.getTrueNumInteractions();
				pileupPtr->PU_NumInteractions  = it.getPU_NumInteractions();
			}
		}
    }
// END Pileup




// START Electron
	electronTCA->Clear("C");
    


// --------- 3. Handle 
	edm::Handle<edm::View<pat::Electron> > electrons; iEvent.getByToken(electronToken_, electrons);


	// --Electron Loop
	for(edm::View<pat::Electron>::const_iterator ele = electrons->begin(); ele != electrons->end(); ele++){
		if(ele->pt() < 5.0) continue;
		
		const edm::Ptr<pat::Electron> elePtr(electrons, ele - electrons->begin());
		npknu::Electron* electronPtr = new ((*electronTCA)[(int)electronTCA->GetEntries()]) npknu::Electron();
        
		// 4 vector
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
	


// START Photon
	photonTCA->Clear("C");
    


// --------- 3. Handle 
	edm::Handle<edm::View<pat::Photon> > photons; iEvent.getByToken(photonToken_, photons);
	
	// --Photon Loop
	for(edm::View<pat::Photon>::const_iterator pho = photons->begin(); pho != photons->end(); ++pho) {
      

	  const edm::Ptr<pat::Photon> phoPtr(photons, pho - photons->begin());
      npknu::Photon* photonPtr =  new ((*photonTCA)[(int)photonTCA->GetEntries()]) npknu::Photon();
		
		// 4-vector
		photonPtr->SetPtEtaPhiE(pho->pt(), pho->eta(), pho->phi(), pho->energy());

		// Photon ID
		photonPtr->vidIsPassLoose    = pho->photonID("cutBasedPhotonID-Spring16-V2p2-loose") ;
		photonPtr->vidIsPassMedium   = pho->photonID("cutBasedPhotonID-Spring16-V2p2-medium");
		photonPtr->vidIsPassTight	 = pho->photonID("cutBasedPhotonID-Spring16-V2p2-tight") ;

		photonPtr->superCluster_eta             = pho->superCluster()->eta()  ;
		photonPtr->superCluster_phi             = pho->superCluster()->phi()  ;

		//std::cout << "Photon PT :" << pho->pt() << std::endl;
		//std::cout << "Photon KNU PT :" << photonPtr->pt << std::endl;

	} // --Photon Loop END

// START JET
	jetTCA->Clear("C");


	// --------3. Handle
	edm::Handle<edm::View<pat::Jet> > jets; iEvent.getByToken(jetToken_, jets);
	
	// --Jet Loop
	for(edm::View<pat::Jet>::const_iterator jet = jets->begin(); jet != jets->end(); ++jet) {
		
		const edm::Ptr<pat::Jet> jtPtr(jets, jet - jets->begin());
		npknu::Jet* jetPtr =  new ((*jetTCA)[(int)jetTCA->GetEntries()]) npknu::Jet();

		jetPtr->SetPtEtaPhiE(jet->pt(), jet->eta(), jet->phi(), jet->energy());
		jetPtr->neutralHadronEnergyFraction   = jet->neutralHadronEnergyFraction() ;
		jetPtr->neutralEmEnergyFraction       = jet->neutralEmEnergyFraction()     ;
		jetPtr->chargedMultiplicity           = jet->chargedMultiplicity()         ;
		jetPtr->neutralMultiplicity           = jet->neutralMultiplicity()         ;
		jetPtr->muonEnergyFraction            = jet->muonEnergyFraction()          ;
		jetPtr->chargedHadronEnergyFraction   = jet->chargedHadronEnergyFraction() ;
		jetPtr->chargedEmEnergyFraction       = jet->chargedEmEnergyFraction()     ;

	}

	outTree->Fill();

	fill_num++;
	
	if(start_num != fill_num){ cout << " ######## Event Number Error " << " " << "run #: " << evtObjPtr->run << " " << "event #: " << evtObjPtr->event << " " << "Start num: "<< start_num << " " << "Vertex num: "<< vertex_num << " " << "Fill num" << fill_num  << endl;
	}
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
