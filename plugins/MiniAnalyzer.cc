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

      //edm::EDGetTokenT<std::vector<pat::Electron> > electronToken_;
      edm::EDGetTokenT<edm::View<pat::Electron> > electronToken_;
      
		TTree* outTree				;
		TClonesArray* evtTCA		;
		TClonesArray* electronTCA	;	  

};


MiniAnalyzer::MiniAnalyzer(const edm::ParameterSet& iConfig){
    edm::Service<TFileService> fs;
    outTree = fs->make<TTree> ("NpKNU","NpKNU");
    //electronToken_ = consumes<std::vector<pat::Electron> >(iConfig.getParameter<edm::InputTag>("electrons"));
    
	evtTCA = new TClonesArray("npknu::Evt"); outTree->Branch("evt"      , &evtTCA );
	electronTCA = new TClonesArray("npknu::Electron"); outTree->Branch("electron"      , &electronTCA );
	electronToken_ = consumes<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("electrons"));
    
	
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

// START Electron
	electronTCA->Clear("C");
    edm::Handle<edm::View<pat::Electron> > electrons; iEvent.getByToken(electronToken_, electrons);

	// EleLoop
	for(edm::View<pat::Electron>::const_iterator ele = electrons->begin(); ele != electrons->end(); ele++){
		if(ele->pt() < 5.0) continue;
		
	const edm::Ptr<pat::Electron> elePtr(electrons, ele - electrons->begin());
	npknu::Electron* electronPtr = new ((*electronTCA)[(int)electronTCA->GetEntries()]) npknu::Electron();
        electronPtr->SetPtEtaPhiE(elePtr->pt(), elePtr->eta(), elePtr->phi(), elePtr->energy());
    
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
      electronPtr->nMissingHits                     = ele->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
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


void MiniAnalyzer::beginJob() { }  
void MiniAnalyzer::endJob() {

   // outTree->Write();
}



DEFINE_FWK_MODULE(MiniAnalyzer);
